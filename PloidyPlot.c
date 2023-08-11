/******************************************************************************************
 *
 *  Hetrez: Analyze and compare het-mers between two k-mer tables to see if they
 *    represent the same sample, different sample but same species, or are from
 *    two different species.
 *
 *  Author:  Gene Myers
 *  Date  :  May, 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <pthread.h>

#define KAMIL

#undef  SOLO_CHECK

#undef  DEBUG_GENERAL
#undef  DEBUG_RECURSION
#undef  DEBUG_THREADS
#undef  DEBUG_BOUNDARY
#undef  DEBUG_SCAN
#undef  DEBUG_BIG_SCAN
#undef  DEBUG_PLOIDY
#undef  ANALYZE_PLOIDY

#include "libfastk.h"
#include "matrix.h"

#include "smu_plot.R.h"

static char *Usage[] = { " [-w<double(6.0)>] [-h<double(4.5)>]",
                         " [-vk] [-lfs] [-pdf] [-T<int(4)>] [-P<dir(/tmp)>]",
                         " [-o<output>] [-e<int(4)>] <source>[.ktab]"
                       };

static int VERBOSE;
static int NTHREADS;   //  At most 64 allowed

#ifdef SOLO_CHECK

static uint8 *CENT;
static int64  CIDX;

#endif

static int ETHRESH;   //  Error threshold

#define SMAX  1000    //  Max. value of CovA+CovB
#define FMAX   500    //  Max. value of min(CovA,CovB)

#define PIXELS  50    //  pixels in area plot (both x & y)
#define PIXEL2 100    //  2*PIXELS

static int BLEVEL;     //  <= 4
static int BWIDTH;     //  = 4^BLEVEL

static int64 MEMORY_LIMIT = 0x100000000ll;  // 4GB 
static int64 Cache_Size;                    // Divided evenly amont the threads

static int KMER;
static int KBYTE;
static int TBYTE;

static int PASS1;

/****************************************************************************************
 *
 *  Print & compare utilities
 *
 *****************************************************************************************/

#define  COUNT_PTR(p)   ((uint16 *) (p+KBYTE))

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
    for (l1 = 0; l1 < 4; l1++)
     for (l0 = 0; l0 < 4; l0++)
       { fmer[i] = t;
         *t++ = dna[l3];
         *t++ = dna[l2];
         *t++ = dna[l1];
         *t++ = dna[l0];
         *t++ = 0;
         i += 1;
       }
}

static void print_hap(uint8 *seq, int len, int half)
{ static int firstime = 1;
  int i, b, h, k;

  if (firstime)
    { firstime = 0;
      setup_fmer_table();
    }

  h = half >> 2;
  b = len >> 2;
  for (i = 0; i < h; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = h << 2; k >= 0; i++)
    { if (i == half)
        printf("%c",dna[(seq[h] >> k) & 0x3]-32);
      else
        printf("%c",dna[(seq[h] >> k) & 0x3]);
      k -= 2;
    }
  for (i = h+1; i < b; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { printf("%c",dna[(seq[b] >> k) & 0x3]);
      k -= 2;
    }
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}


/****************************************************************************************
 *
 *  Find Het Pairs by merging 4 lists:
 *    Out of core: stream indices from fng[a] to end[a] for a in [0,3]
 *    In core:     table pointer from ptr[a] to eptr[a] for a in [0,3]
 *
 *****************************************************************************************/

typedef struct
  { int           level;    //  position of variation
    int64        *bound;    //  17-element array of next level transition points
                          //  Stream scans (big & small):
    Kmer_Stream  *fng[4];   //  4 streams for scan
    int64         end[4];   //  index of end of each scan
                          //  Small scans:
    int           tid;      //  thread assigned to this subtree
    uint8        *cache;    //  in-core cache for subtrees that fit
                          //  In-core scans:
    int64         cidx;     //  Table index of 1st cache entry
    uint8        *ept[4];   //  Each list is in [ptr[a],eptr[a])
    uint8        *ptr[4];   //    in steps of TBYTES
                          //  Plot:
    int64       **plot;     //  Accumulate A+B, B/(A+B) pairs here

#ifdef SOLO_CHECK
    uint8        *cptr;
#endif
  } TP;

static uint8 *Pair;  //  Incidence array (# of table entries)

static uint8 Prefix[4] = { 0x3f, 0x0f, 0x03, 0x00 };
static uint8 Shift[4]  = { 6, 4, 2, 0 };

static void *analysis_thread_1(void *args)
{ TP           *parm  = (TP *) args;
  int64        *end   = parm->end;
  Kmer_Stream **fng   = parm->fng;
  int           level = parm->level;
  int64        *bound = parm->bound;

  int ll = ((level+1)>>2);
  int ls = Shift[(level+1)&0x3];

  int mask  = Prefix[level&0x3];
  int offs  = (level >> 2) + 1;
  int rem   = KBYTE - offs;

  uint8 *ent[4];
  int    lst[4];
  int    in[4], itop;
  int    cnt[4];
  int    mc, hc;
  uint8 *mr, *hr;
  int    a, i, x;
  
  for (a = 0; a < 4; a++)
    { ent[a] = Current_Entry(fng[a],NULL);
      lst[a] = (ent[a][ll] >> ls) & 0x3;
      if (fng[a]->cidx < end[a])
        { x = (ent[a][ll] >> ls) & 0x3;
          while (x > lst[a])
            bound[(a<<2)+(++(lst[a]))] = fng[a]->cidx;
        }
    }

#ifdef DEBUG_SCAN
  for (a = 0; a < 4; a++)
    { printf(" %c %10lld: ",dna[a],fng[a]->cidx);
      print_hap(ent[a],KMER,level);
    }
  printf("\n");
#endif

  while (1)
    { for (a = 0; a < 4; a++)
        if (fng[a]->cidx < end[a])
          break;
      if (a >= 4)
        break;

      mr = ent[a]+offs;
      mc = mr[-1] & mask;
      in[0] = a;
      itop  = 1;
      for (a++; a < 4; a++)
        if (fng[a]->cidx < end[a])
          { hr = ent[a]+offs;
            hc = hr[-1] & mask;
            if (hc == mc)
              { int v = mycmp(hr,mr,rem);
                if (v == 0)
                  in[itop++] = a;
                else if (v < 0)
                  { mc = hc;
                    mr = hr;
                    in[0] = a;
                    itop  = 1;
                  }
              }
            else if (hc < mc)
              { mc = hc;
                mr = hr;
                in[0] = a;
                itop  = 1;
              }
          }

      if (itop > 1)
        { cnt[0] = *((uint16 *) (ent[in[0]]+KBYTE));
          for (i = 1; i < itop; i++)
            { cnt[i] = *((uint16 *) (ent[in[i]]+KBYTE));
              for (a = 0; a < i; a++)
                { x = cnt[a]+cnt[i];
                  if (x <= SMAX)
                    { Pair[fng[in[i]]->cidx] += 1;
                      Pair[fng[in[a]]->cidx] += 1;
                    }
                }
            }
        }

      for (i = 0; i < itop; i++)
        { Kmer_Stream *t;

          a = in[i];
          t = fng[a];
#ifdef DEBUG_SCAN
          if (i == 0) printf("\n");
          printf(" %c %10lld: ",dna[a&0x3],fng[a]->cidx); 
          print_hap(ent[a],KMER,level);
          printf("\n");
#endif
          Next_Kmer_Entry(t);
          if (t->cidx < end[a])
            { Current_Entry(t,ent[a]);
              x = (ent[a][ll] >> ls) & 0x3;
              while (x > lst[a])
                bound[(a<<2)+(++(lst[a]))] = t->cidx;
            }
        }

#ifdef DEBUG_SCAN
      printf("\n");
      for (a = 0; a < 4; a++)
        { printf(" %c %10lld: ",dna[a],fng[a]->cidx);
          print_hap(ent[a],KMER,level);
        }
      printf("\n");
#endif
    }

  for (a = 0; a < 4; a++)
    free(ent[a]);

  return (NULL);
}

static void *analysis_thread_2(void *args)
{ TP           *parm  = (TP *) args;
  int64        *end   = parm->end;
  Kmer_Stream **fng   = parm->fng;
  int           level = parm->level;
  int64        *bound = parm->bound;
  int64       **plot  = parm->plot;

  int ll = ((level+1)>>2);
  int ls = Shift[(level+1)&0x3];

  int mask  = Prefix[level&0x3];
  int offs  = (level >> 2) + 1;
  int rem   = KBYTE - offs;

  uint8 *ent[4];
  int    lst[4];
  int    in[4], itop;
  int    cnt[4];
  int    mc, hc;
  uint8 *mr, *hr;
  int    a, i, x;
  
  for (a = 0; a < 4; a++)
    { ent[a] = Current_Entry(fng[a],NULL);
      lst[a] = (ent[a][ll] >> ls) & 0x3;
      if (fng[a]->cidx < end[a])
        { x = (ent[a][ll] >> ls) & 0x3;
          while (x > lst[a])
            bound[(a<<2)+(++(lst[a]))] = fng[a]->cidx;
        }
    }

#ifdef DEBUG_SCAN
  for (a = 0; a < 4; a++)
    { printf(" %c %10lld: ",dna[a],fng[a]->cidx);
      print_hap(ent[a],KMER,level);
    }
  printf("\n");
#endif

  while (1)
    { for (a = 0; a < 4; a++)
        if (fng[a]->cidx < end[a])
          break;
      if (a >= 4)
        break;

      mr = ent[a]+offs;
      mc = mr[-1] & mask;
      in[0] = a;
      itop  = 1;
      for (a++; a < 4; a++)
        if (fng[a]->cidx < end[a])
          { hr = ent[a]+offs;
            hc = hr[-1] & mask;
            if (hc == mc)
              { int v = mycmp(hr,mr,rem);
                if (v == 0)
                  in[itop++] = a;
                else if (v < 0)
                  { mc = hc;
                    mr = hr;
                    in[0] = a;
                    itop  = 1;
                  }
              }
            else if (hc < mc)
              { mc = hc;
                mr = hr;
                in[0] = a;
                itop  = 1;
              }
          }

      if (itop > 1)
#ifdef SOLO_CHECK
        { for (i = 0; i < itop; i++)
            if (mycmp(ent[in[i]],CENT,KBYTE) == 0)
              for (a = 0; a < itop; a++)
                if (a != i)
                  { printf("  ");
                    print_hap(ent[in[a]],KMER,level);
                    printf(": %d\n",*((uint16 *) (ent[in[a]]+KBYTE)));
                  }
        }
#else
        { cnt[0] = *((uint16 *) (ent[in[0]]+KBYTE));
          for (i = 1; i < itop; i++)
            { cnt[i] = *((uint16 *) (ent[in[i]]+KBYTE));
              if (Pair[fng[in[i]]->cidx] <= 1)
                for (a = 0; a < i; a++)
                  { x = cnt[a]+cnt[i];
                    if (x <= SMAX && Pair[fng[in[a]]->cidx] <= 1)
                      { if (cnt[a] < cnt[i])
                          plot[x][cnt[a]] += 1;
                        else
                          plot[x][cnt[i]] += 1;
                      }
                  }
            }
        }
#endif

      for (i = 0; i < itop; i++)
        { Kmer_Stream *t;

          a = in[i];
          t = fng[a];
#ifdef DEBUG_SCAN
          if (i == 0) printf("\n");
          printf(" %c %10lld: ",dna[a&0x3],fng[a]->cidx); 
          print_hap(ent[a],KMER,level);
          printf("\n");
#endif
          Next_Kmer_Entry(t);
          if (t->cidx < end[a])
            { Current_Entry(t,ent[a]);
              x = (ent[a][ll] >> ls) & 0x3;
              while (x > lst[a])
                bound[(a<<2)+(++(lst[a]))] = t->cidx;
            }
        }

#ifdef DEBUG_SCAN
      printf("\n");
      for (a = 0; a < 4; a++)
        { printf(" %c %10lld: ",dna[a],fng[a]->cidx);
          print_hap(ent[a],KMER,level);
        }
      printf("\n");
#endif
    }

  for (a = 0; a < 4; a++)
    free(ent[a]);

  return (NULL);
}

static void *analysis_in_core_1(void *args)
{ TP          *parm  = (TP *) args;
  uint8      **ptr   = parm->ptr;
  uint8      **ept   = parm->ept;
  int          level = parm-> level;
  uint8      **bound = (uint8 **) (parm->bound);
  uint8       *cache = parm->cache;
  int64        aidx  = parm->cidx;

  int ll = ((level+1)>>2);
  int ls = Shift[(level+1)&0x3];

  int mask  = Prefix[level&0x3];
  int offs  = (level >> 2) + 1;
  int rem   = KBYTE - offs;

  int    lst[4];
  int    in[4], itop;
  int    cnt[4];
  int    mc, hc;
  uint8 *mr, *hr;
  int    a, i, x;

  for (a = 0; a < 4; a++)
    { lst[a] = 0;
      if (ptr[a] < ept[a])
        { x = (ptr[a][ll] >> ls) & 0x3;
          while (x > lst[a])
            bound[(a<<2)+(++(lst[a]))] = ptr[a];
        }
    }

#ifdef DEBUG_SCAN
  for (a = 0; a < 4; a++)
    { printf(" %c %10ld: ",dna[a],(ptr[a]-parm->cache)/TBYTE); 
      print_hap(ptr[a],KMER,level);
    }
  printf("\n");
#endif

  while (1)
    { for (a = 0; a < 4; a++)
        if (ptr[a] < ept[a])
          break;
      if (a >= 4)
        break;

      mr = ptr[a]+offs;
      mc = mr[-1] & mask;
      in[0] = a;
      itop  = 1;
      for (a++; a < 4; a++)
        if (ptr[a] < ept[a])
          { hr = ptr[a]+offs;
            hc = hr[-1] & mask;
            if (hc == mc)
              { int v = mycmp(hr,mr,rem);
                if (v == 0)
                  in[itop++] = a;
                else if (v < 0)
                  { mc = hc;
                    mr = hr;
                    in[0] = a;
                    itop  = 1;
                  }
              }
            else if (hc < mc)
              { mc = hc;
                mr = hr;
                in[0] = a;
                itop  = 1;
              }
          }

      if (itop > 1)
        { cnt[0] = *((uint16 *) (ptr[in[0]]+KBYTE));
          for (i = 1; i < itop; i++)
            { cnt[i] = *((uint16 *) (ptr[in[i]]+KBYTE));
              for (a = 0; a < i; a++)
                { x = cnt[a]+cnt[i];
                  if (x <= SMAX)
                    { Pair[aidx + (ptr[in[i]]-cache)/TBYTE] += 1;
                      Pair[aidx + (ptr[in[a]]-cache)/TBYTE] += 1;
                    }
                }
            }
        }

      for (i = 0; i < itop; i++)
        { a = in[i];
#ifdef DEBUG_SCAN
          if (i == 0) printf("\n");
          printf("%c %10ld: ",dna[a&0x3],(ptr[a]-parm->cache)/TBYTE); 
          print_hap(ptr[a],KMER,level);
          printf("\n");
#endif
          ptr[a] += TBYTE;
          if (ptr[a] < ept[a])
            { x = (ptr[a][ll] >> ls) & 0x3;
              while (x > lst[a])
                bound[(a<<2)+(++(lst[a]))] = ptr[a];
            }
        }

#ifdef DEBUG_SCAN
      for (a = 0; a < 4; a++)
        { printf(" %c %10ld: ",dna[a],(ptr[a]-parm->cache)/TBYTE); 
          print_hap(ptr[a],KMER,level);
        }
      printf("\n");
#endif
    }

  return (NULL);
}

static void *analysis_in_core_2(void *args)
{ TP          *parm  = (TP *) args;
  uint8      **ptr   = parm->ptr;
  uint8      **ept   = parm->ept;
  int          level = parm-> level;
  uint8      **bound = (uint8 **) (parm->bound);
  int64      **plot  = parm->plot;
  uint8       *cache = parm->cache;
  int64        aidx  = parm->cidx;

  int ll = ((level+1)>>2);
  int ls = Shift[(level+1)&0x3];

  int mask  = Prefix[level&0x3];
  int offs  = (level >> 2) + 1;
  int rem   = KBYTE - offs;

  int    lst[4];
  int    in[4], itop;
  int    cnt[4];
  int    mc, hc;
  uint8 *mr, *hr;
  int    a, i, x;

  for (a = 0; a < 4; a++)
    { lst[a] = 0;
      if (ptr[a] < ept[a])
        { x = (ptr[a][ll] >> ls) & 0x3;
          while (x > lst[a])
            bound[(a<<2)+(++(lst[a]))] = ptr[a];
        }
    }

#ifdef DEBUG_SCAN
  for (a = 0; a < 4; a++)
    { printf(" %c %10ld: ",dna[a],(ptr[a]-parm->cache)/TBYTE); 
      print_hap(ptr[a],KMER,level);
    }
  printf("\n");
#endif

  while (1)
    { for (a = 0; a < 4; a++)
        if (ptr[a] < ept[a])
          break;
      if (a >= 4)
        break;

      mr = ptr[a]+offs;
      mc = mr[-1] & mask;
      in[0] = a;
      itop  = 1;
      for (a++; a < 4; a++)
        if (ptr[a] < ept[a])
          { hr = ptr[a]+offs;
            hc = hr[-1] & mask;
            if (hc == mc)
              { int v = mycmp(hr,mr,rem);
                if (v == 0)
                  in[itop++] = a;
                else if (v < 0)
                  { mc = hc;
                    mr = hr;
                    in[0] = a;
                    itop  = 1;
                  }
              }
            else if (hc < mc)
              { mc = hc;
                mr = hr;
                in[0] = a;
                itop  = 1;
              }
          }

      if (itop > 1)
#ifdef SOLO_CHECK
        { for (i = 0; i < itop; i++)
            if (mycmp(ptr[in[i]],CENT,KBYTE) == 0)
              for (a = 0; a < itop; a++)
                if (a != i)
                  { printf("  ");
                    print_hap(ptr[in[a]],KMER,level);
                    printf(": %d\n",*((uint16 *) (ptr[in[a]]+KBYTE)));
                  }
        }
#else
        { cnt[0] = *((uint16 *) (ptr[in[0]]+KBYTE));
          for (i = 1; i < itop; i++)
            { cnt[i] = *((uint16 *) (ptr[in[i]]+KBYTE));
              if (Pair[aidx+(ptr[in[i]]-cache)/TBYTE] <= 1)
                for (a = 0; a < i; a++)
                  { x = cnt[a]+cnt[i];
                    if (x <= SMAX && Pair[aidx+(ptr[in[a]]-cache)/TBYTE] <= 1)
                      { if (cnt[a] < cnt[i])
                          plot[x][cnt[a]] += 1;
                        else
                          plot[x][cnt[i]] += 1;
                      }
                  }
            }
        }
#endif

      for (i = 0; i < itop; i++)
        { a = in[i];
#ifdef DEBUG_SCAN
          if (i == 0) printf("\n");
          printf("%c %10ld: ",dna[a&0x3],(ptr[a]-parm->cache)/TBYTE); 
          print_hap(ptr[a],KMER,level);
          printf("\n");
#endif
          ptr[a] += TBYTE;
          if (ptr[a] < ept[a])
            { x = (ptr[a][ll] >> ls) & 0x3;
              while (x > lst[a])
                bound[(a<<2)+(++(lst[a]))] = ptr[a];
            }
        }

#ifdef DEBUG_SCAN
      for (a = 0; a < 4; a++)
        { printf(" %c %10ld: ",dna[a],(ptr[a]-parm->cache)/TBYTE); 
          print_hap(ptr[a],KMER,level);
        }
      printf("\n");
#endif
    }

  return (NULL);
}


/****************************************************************************************
 *
 *  Find Het Pairs in top level nodes (level < BLEVEL <= 3) by paneling 4 merge intervals
 *    with all threads.
 *
 *****************************************************************************************/

static uint8 *Divpt;

static void big_window(int64 *adiv, int level, TP *parm)
{ int64 bound[17];
  int   a;

#ifdef DEBUG_GENERAL
  printf("Doing big %d: %lld\n",level,adiv[4]-adiv[0]); fflush(stdout);
#endif

  { int64        e;
    uint8        lm;
    int          t, ls;
    Kmer_Stream *T;
#ifndef DEBUG_THREADS
    pthread_t    threads[NTHREADS];
#endif

    T = parm[0].fng[0];

    ls = Shift[level];      //  level < BLEVEL
    lm = 0xff ^ (0x3<<ls);
    for (a = 0; a < 4; a++)
      GoTo_Kmer_Index(parm[0].fng[a],adiv[a]);
    parm[0].level = level;
    parm[0].bound = bound;
    for (t = 1; t < NTHREADS; t++)
      { parm[t].level = level;
        parm[t].bound = bound;
        parm[t-1].end[0] = e = adiv[0] + t*(adiv[1]-adiv[0])/NTHREADS;
        T = parm[t].fng[0];
        GoTo_Kmer_Index(T,e);
        Current_Entry(T,Divpt);
        for (a = 1; a < 4; a++)
          { Divpt[0] = (Divpt[0] & lm) | (a << ls);
            T = parm[t].fng[a];
            GoTo_Kmer_Entry(T,Divpt);
            parm[t-1].end[a] = T->cidx;
          }
      }
    for (a = 0; a < 4; a++)
      parm[NTHREADS-1].end[a] = adiv[a+1];

#ifdef DEBUG_BIG_SCAN
    for (a = 0; a < 4; a++)
      for (t = 0; t < NTHREADS; t++)
        { printf("%c/%d %10lld: ",dna[a],t,parm[t].fng[a]->cidx);
          Current_Entry(parm[t].fng[a],Divpt);
          print_hap(Divpt,KMER,level);
          printf("\n");
        }
#endif

    for (a = 0; a < 16; a += 4)
      { bound[a] = adiv[a>>2];
        for (t = 1; t < 4; t++)
          bound[a+t] = -1;
      }

    if (PASS1)
#ifdef DEBUG_THREADS
      { for (t = 0; t < NTHREADS; t++)
          analysis_thread_1(parm+t);
#else
      { for (t = 1; t < NTHREADS; t++)
          pthread_create(threads+t,NULL,analysis_thread_1,parm+t);
        analysis_thread_1(parm);
        for (t = 1; t < NTHREADS; t++)
          pthread_join(threads[t],NULL);
#endif
      }
    else
#ifdef DEBUG_THREADS
      { for (t = 0; t < NTHREADS; t++)
          analysis_thread_2(parm+t);
#else
      { for (t = 1; t < NTHREADS; t++)
          pthread_create(threads+t,NULL,analysis_thread_2,parm+t);
        analysis_thread_2(parm);
        for (t = 1; t < NTHREADS; t++)
          pthread_join(threads[t],NULL);
#endif
      }

    for (a = 0; a < 16; a += 4)
      for (t = 1; t < 4; t++)
        if (bound[a+t] < 0)
          bound[a+t] = adiv[(a>>2)+1];
    bound[16] = adiv[4];

#ifdef DEBUG_BOUNDARY
    T = parm[0].fng[0];
    for (a = 0; a <= 16; a++)
      { if (a > 0)
          { GoTo_Kmer_Index(T,bound[a]-1);
            Current_Entry(T,Divpt);
            printf("%c %10lld: ",dna[a&0x3],T->cidx);
            print_hap(Divpt,KMER,level+1);
            printf("\n");
          }
        if (a < 16)
          { GoTo_Kmer_Index(T,bound[a]);
            Current_Entry(T,Divpt);
            printf("%c %10lld: ",dna[a&0x3],T->cidx);
            print_hap(Divpt,KMER,level+1);
            printf("\n");
          }
      }
#endif
    }

  level += 1;
  if (level < BLEVEL)
    for (a = 0; a < 16; a += 4)
      big_window(bound+a,level,parm);
}


/****************************************************************************************
 *
 *  Find Het Pairs for a lower level node by list merging
 *
 *****************************************************************************************/

static void in_core_recursion(uint8 **aptr, int level, TP *parm)
{ uint8 *bound[17];
  int    a;

#ifdef SOLO_CHECK
  if (aptr[0] <= parm->cptr && parm->cptr < aptr[4])
    { printf("Inside %ld-%ld (%d %d)\n",
             (aptr[0]-parm->cache)/TBYTE,(aptr[4]-parm->cache)/TBYTE,parm->tid,level);
      printf("     ");
      print_hap(aptr[0],KMER,level);
      printf(" : ");
      print_hap(aptr[4],KMER,level);
      printf("\n");
    }
#endif

  if (aptr[4]-aptr[0] <= TBYTE) return;

#ifdef DEBUG_RECURSION
  printf("Doing in-core %d: %ld\n",level,(aptr[4]-aptr[0])/TBYTE);
#endif

  { int t;

    parm->ptr[0] = aptr[0];
    for (a = 1; a < 4; a++)
      parm->ptr[a] = parm->ept[a-1] = aptr[a];
    parm->ept[3] = aptr[4];

    parm->level = level;
    parm->bound = (int64 *) bound;

    for (a = 0; a < 16; a += 4)
      { bound[a] = aptr[a>>2];
        for (t = 1; t < 4; t++)
          bound[a+t] = NULL;
      }
    bound[16] = aptr[4];

    if (PASS1)
      analysis_in_core_1(parm);
    else
      analysis_in_core_2(parm);

    for (a = 0; a < 16; a += 4)
      for (t = 1; t < 4; t++)
        if (bound[a+t] == NULL)
          bound[a+t] = bound[a+4];

#ifdef DEBUG_BOUNDARY
    { uint8 *cache = parm->cache;

      if (level+1 < KMER)
        for (a = 0; a <= 16; a++)
          { if (a > 0)
              { printf("%c %10ld: ",dna[a&0x3],(bound[a]-cache)/TBYTE); 
                print_hap(bound[a],KMER,level+1);
                printf("\n");
              }
            if (a < 16)
              { printf("%c %c %10ld: ",dna[a>>2],dna[a&0x3],(bound[a]-cache)/TBYTE); 
                print_hap(bound[a],KMER,level+1);
                printf(" - ");
              }
          }
    }
#endif
  }

  level += 1;
  if (level < KMER)
    for (a = 0; a < 16; a += 4)
      in_core_recursion(bound+a,level,parm);
}

static pthread_mutex_t TMUTEX;
static pthread_cond_t  TCOND;

static int *Tstack;
static int  Tavail;

static void small_recursion(int64 *adiv, int level, TP *parm)
{
#ifdef DEBUG_RECURSION
  printf("Doing small %d: %lld [%lld-%lld]\n",level,adiv[4]-adiv[0],adiv[0],adiv[4]);
#endif

  if (adiv[4]-adiv[0] < Cache_Size)
    { uint8       *aptr[5];

      { uint8       *C = parm->cache;
        Kmer_Stream *T;
        int64   i;
        int     a;

#ifdef SOLO_CHECK
        if (adiv[0] <= CIDX && CIDX < adiv[4])
          { printf("Heading in %lld-%lld (%d %d)\n",adiv[0],adiv[4],parm->tid,level);
            parm->cptr = parm->cache + (CIDX - adiv[0])*TBYTE;
          }  
        else
          parm->cptr = NULL;
#endif

        T = parm->fng[0];
        GoTo_Kmer_Index(T,adiv[0]);
        parm->cidx = T->cidx;
        for (i = 0; T->cidx < adiv[4]; i++)
          { Current_Entry(T,C);
            C += TBYTE;
            Next_Kmer_Entry(T);
          } 
        for (a = 1; a <= 4; a++)
          aptr[a] = parm->cache + (adiv[a]-adiv[0])*TBYTE;
        aptr[0] = parm->cache;
      }

      in_core_recursion(aptr,level,parm);

      return;
    }

  { int64 bound[17];
    int   a;

    { int t;

      for (a = 0; a < 4; a++)
        { parm->end[a] = adiv[a+1];
          GoTo_Kmer_Index(parm->fng[a],adiv[a]);
        }
      parm->level = level;
      parm->bound = bound;

      for (a = 0; a < 16; a += 4)
        { bound[a] = adiv[a>>2];
          for (t = 1; t < 4; t++)
            bound[a+t] = -1;
        }

      if (PASS1)
        analysis_thread_1(parm);
      else
        analysis_thread_2(parm);

      for (a = 0; a < 16; a += 4)
        for (t = 1; t < 4; t++)
          if (bound[a+t] < 0)
            bound[a+t] = bound[a+4];
      bound[16] = adiv[4];

#ifdef DEBUG_BOUNDARY
      { Kmer_Stream *T;

        T = parm->fng[0];
        if (level+1 < KMER)
          for (a = 0; a <= 16; a++)
            { if (a > 0)
                { GoTo_Kmer_Index(T,bound[a]-1);
                  Current_Entry(T,Divpt);
                  printf("%c %10lld: ",dna[a&0x3],T->cidx);
                  print_hap(Divpt,KMER,level+1);
                  printf("\n");
                }
              if (a < 16)
                { GoTo_Kmer_Index(T,bound[a]);
                  Current_Entry(T,Divpt);
                  printf("%c %10lld: ",dna[a&0x3],T->cidx); 
                  print_hap(Divpt,KMER,level+1);
                  printf("\n");
                }
            }
      }
#endif
    }

    level += 1;
    if (level < KMER)
      for (a = 0; a < 16; a += 4)
        small_recursion(bound+a,level,parm);
  }
}

static void *small_window(void *args)
{ TP *parm  = (TP *) args;
  int tid   = parm->tid;

  int64 adiv[5];

  { uint8 divpt[TBYTE];
    Kmer_Stream *T;
    int a, x;

    x = parm->level;
    T = parm->fng[0];

    for (a = 0; a < KBYTE; a++)
      divpt[a] = 0;
    for (a = 0; a < 4; a++)
      { if (BLEVEL == 4)
          { divpt[0] = x;
            divpt[1] = (a<<6);
          }
        else
          divpt[0] = (((x<<2) | a) << (6-2*BLEVEL));
        GoTo_Kmer_Entry(T,divpt);
        adiv[a] = T->cidx; 
      }
    if (++x < BWIDTH)
      { divpt[0] = (x << (8-2*BLEVEL));
        divpt[1] = 0;
        GoTo_Kmer_Entry(T,divpt);
        adiv[4] = T->cidx; 
      }
    else
      adiv[4] = T->nels;
  }

  small_recursion(adiv,BLEVEL,parm);

  pthread_mutex_lock(&TMUTEX);
    Tstack[Tavail++] = tid;
  pthread_mutex_unlock(&TMUTEX);

  pthread_cond_signal(&TCOND);

  return (NULL);
}


/****************************************************************************************
 *
 *  Ploidy fitting
 *
 *****************************************************************************************/

char  *Ploidy[5]   = { "diploid", "triploid", "tetraploid", "hexaploid", "octaploid" };
int    Nloidy[6]   = { 2, 3, 4, 6, 8 };

int    Nsmu[5]     = { 1, 1, 2, 3, 4 };
char  *Names[5][4] = { { "AB", NULL, NULL, NULL },
                       { "AAB", NULL, NULL, NULL },
                       { "AAAB", "AABB", NULL, NULL },
                       { "5A1B", "4A2B", "3A3B", NULL },
                       { "7A1B", "6A2B", "5A3B", "4A4B" } };
int    Rail[5][4]  = { { 48, 0, 0, 0 },
                       { 33, 0, 0, 0 },
                       { 25, 48, 0, 0 },
                       { 17, 33, 48, 0 },
                       { 12, 25, 37, 48 } };
int64  Deco[5][4];

double PFrac[5] = { 1., 1.5, 2., 3., 4. };

int determine_ploidy(int cover, int64 *row)
{ double C8[PIXELS];
  double C6[PIXELS];
  double C4[PIXELS];
  double C38[PIXELS];
  double C3[PIXELS];
  double C2[PIXELS];
  double D2, T3, Q4, Q2, H6, H3, H2, O8, O4, O3, O2;
  int64  Df, Tf, Qf, Hf, Of;
  int64  lsqr[5];
  int64  Ds, Ts, Qs, Hs, Os;
  int64  DB1, TB1, QB1, QB2, HB1, HB2, HB3;
  int64  OB1, OB2, OB3, OB4;
  int    ploidy;
  int    j, first, hmax;

  hmax = 0;
  for (j = 0; j < PIXELS; j++)
    if (row[j] > hmax)
      hmax = row[j];

  for (j = 0; j < PIXELS; j++)
    { double del;
      int    k;

#ifdef DEBUG_PLOIDY
      printf(" %2d: %10lld\n",j,row[j]);
#endif
      if (j+1 < PIXELS && row[j] < row[j+1])
        continue;
      if (j > 0 && row[j] < row[j-1])
        continue;
      if (row[j] < .10*hmax)
        continue;
      del = PIXELS;
      for (k = 0; k < 5; k++)
        if (fabs(PIXELS/PFrac[k] - j) < del)
          { del    = fabs(PIXELS/PFrac[k] - j);
            ploidy = k;
          }
#ifdef ANALYZE_PLOIDY
      printf("  Peak at %d (%d)\n",j,ploidy);
#endif
      break;
    }

  for (j = 0; j < PIXELS; j++)
    if (row[j] > 0)
      break;
  first = j;

  //  Compute smudge double Poisson distributions

  { double c8, c6, c4, c38, c3, c2;
    double d8, d6, d4, d38, d3, d2;
    int    a, b, k;

    c8  = cover/8.;
    c6  = cover/6.;
    c4  = cover/4.;
    c38 = cover*.375;
    c3  = cover/3.;
    c2  = cover/2.;

    d8  = exp(-c8)*exp(c8-cover);
    d6  = exp(-c6)*exp(c6-cover);
    d4  = exp(-c4)*exp(c4-cover);
    d38 = exp(-c38)*exp(c38-cover);
    d3  = exp(-c3)*exp(c3-cover);
    d2  = exp(-c2)*exp(c2-cover);
    for (k = 1; k <= cover; k++)
      { d8  *= ((cover-c8)/k);
        d6  *= ((cover-c6)/k);
        d4  *= ((cover-c4)/k);
        d38 *= ((cover-c38)/k);
        d3  *= ((cover-c3)/k);
        d2  *= ((cover-c2)/k);
      }

    b   = 0;
    for (j = first; j < PIXELS; j++)
      { a = (j*cover)/(2*PIXELS);
        for (k = b+1; k <= a; k++)
          { d8  *= (c8/k);
            d6  *= (c6/k);
            d4  *= (c4/k);
            d38 *= (c38/k);
            d3  *= (c3/k);
            d2  *= (c2/k);
          }
        for (k = cover-b; k > cover-a; k--)
          { d8  *= k/(cover-c8);
            d6  *= k/(cover-c6);
            d4  *= k/(cover-c4);
            d38 *= k/(cover-c38);
            d3  *= k/(cover-c3);
            d2  *= k/(cover-c2);
          }
        C8[j]  = d8;
        C6[j]  = d6;
        C4[j]  = d4;
        C38[j] = d38;
        C3[j]  = d3;
        C2[j]  = d2;
        b = a;
      }
  }

  //  Fit DIPOID

  { double a2, o2;
    int64  r, v;

    a2 = o2 = 0.;
    for (j = first; j < PIXELS; j++)
      { a2 += C2[j]*C2[j];
        o2 += C2[j]*row[j];
      }
    D2 = o2/a2;

    Ds = Df = 0;
    DB1 = 0;
    for (j = first; j < PIXELS; j++)
      { v = (int64) (C2[j]*D2);
        r = row[j];
        Df += (v-r)*(v-r);
        Ds += v;
        DB1 += D2*C2[j];
      }
  }

  //  Fit TRIPLOID

  { double a3, o3;
    int64  r, v;

    a3 = o3 = 0.;
    for (j = first; j < PIXELS; j++)
      { a3 += C3[j]*C3[j];
        o3 += C3[j]*row[j];
      }
    T3 = o3/a3;

    Ts = Tf = 0;
    TB1 = 0;
    for (j = first; j < PIXELS; j++)
      { v = (int64) (C3[j]*T3);
        r = row[j];
        Tf += (v-r)*(v-r);
        Ts += v;
        TB1 += T3*C3[j];
      }
  }

  //  Fit QUADRAPLOID

  { double a44, a42, a22, o4, o2;
    Double_Matrix Q, B;
    LU_Factor    *QU;
    double        m[4], q[2];
    int           stable;
    int64         r, v;

    a44 = a42 = a22 = o4 = o2 = 0.;
    for (j = first; j < PIXELS; j++)
      { a44 += C4[j]*C4[j];
        a42 += C4[j]*C2[j];
        a22 += C2[j]*C2[j];
        o4 += C4[j]*row[j];
        o2 += C2[j]*row[j];
      }

    Q.n = B.n = 2;
    Q.m = m;
    B.m = q;
    m[0] = a44; m[1] = a42;
    m[2] = a42; m[3] = a22;
    q[0] = o4;  q[1] = o2;

    QU = LU_Decompose(&Q,&stable);
    LU_Solve(&B,QU);

    Q4 = q[0]; Q2 = q[1];

    Qs = Qf = 0;
    QB2 = 0;
    QB1 = 0;
    for (j = first; j < PIXELS; j++)
      { v = (int64) (C4[j]*Q4+C2[j]*Q2);
        r = row[j];
        Qf += (v-r)*(v-r);
        Qs += v;
        QB1 += Q4*C4[j];
        QB2 += Q2*C2[j];
      }
  }

  //  Fit HEXAPLOID

  { double a66, a63, a62, a33, a32, a22, o6, o3, o2;
    Double_Matrix Q, B;
    LU_Factor    *QU;
    double        m[9], q[3];
    int           stable;
    int64         r, v;

    a66 = a63 = a62 = a33 = a32 = a22 = o6 = o3 = o2 = 0.;
    for (j = first; j < PIXELS; j++)
      { a66 += C6[j]*C6[j];
        a63 += C6[j]*C3[j];
        a62 += C6[j]*C2[j];
        a33 += C3[j]*C3[j];
        a32 += C3[j]*C2[j];
        a22 += C2[j]*C2[j];
        o6 += C6[j]*row[j];
        o3 += C3[j]*row[j];
        o2 += C2[j]*row[j];
      }

    Q.n = B.n = 3;
    Q.m = m;
    B.m = q;
    m[0] = a66; m[1] = a63; m[2] = a62;
    m[3] = a63; m[4] = a33; m[5] = a32;
    m[6] = a62; m[7] = a32; m[8] = a22;
    q[0] = o6;  q[1] = o3;  q[2] = o2;

    QU = LU_Decompose(&Q,&stable);
    LU_Solve(&B,QU);

    H6 = q[0]; H3 = q[1]; H2 = q[2];

    Hs = Hf = 0;
    HB3 = 0;
    HB2 = 0;
    HB1 = 0;
    for (j = first; j < PIXELS; j++)
      { v = (int64) (C2[j]*H2 + C3[j]*H3 + C6[j]*H6);
        r = row[j];
        Hf += (v-r)*(v-r);
        Hs += v;
        HB1 += H6*C6[j];
        HB2 += H3*C3[j];
        HB3 += H2*C2[j];
      }
  }

  //  Fit OCTAPLOID

  { double a88, a84, a83, a82, a44, a43, a42, a33, a32, a22, o8, o4, o3, o2;
    Double_Matrix Q, B;
    LU_Factor    *QU;
    double        m[16], q[4];
    int           stable;
    int64         r, v;

    a88 = a84 = a83 = a82 = a44 = a43 = a42 = a33 = a32 = a22 = o8 = o4 = o3 = o2 = 0.;
    for (j = first; j < PIXELS; j++)
      { a88 += C8[j]*C8[j];
        a84 += C8[j]*C4[j];
        a83 += C8[j]*C38[j];
        a82 += C8[j]*C2[j];
        a44 += C4[j]*C4[j];
        a43 += C4[j]*C38[j];
        a42 += C4[j]*C2[j];
        a33 += C38[j]*C38[j];
        a32 += C38[j]*C2[j];
        a22 += C2[j]*C2[j];

        o8 += C8[j]*row[j];
        o4 += C4[j]*row[j];
        o3 += C38[j]*row[j];
        o2 += C2[j]*row[j];
      }

    Q.n = B.n = 4;
    Q.m = m;
    B.m = q;
    m[ 0] = a88; m[ 1] = a84; m[ 2] = a83; m[ 3] = a82;
    m[ 4] = a84; m[ 5] = a44; m[ 6] = a43; m[ 7] = a42;
    m[ 8] = a83; m[ 9] = a43; m[10] = a33; m[11] = a32;
    m[12] = a82; m[13] = a42; m[14] = a32; m[15] = a22;
    q[ 0] = o8;  q[ 1] = o4;  q[ 2] = o3;  q[ 3] = o2;

    QU = LU_Decompose(&Q,&stable);
    LU_Solve(&B,QU);

    O8 = q[0]; O4 = q[1]; O3 = q[2]; O2 = q[3];

    Os = Of = 0;
    OB4 = 0;
    OB3 = 0;
    OB2 = 0;
    OB1 = 0;
    for (j = first; j < PIXELS; j++)
      { v = (int64) (C2[j]*O2 + C38[j]*O3 + C4[j]*O4 + C8[j]*O8);
        r = row[j];
        Of += (v-r)*(v-r);
        Os += v;
        OB1 += O8*C8[j];
        OB2 += O4*C4[j];
        OB3 += O3*C38[j];
        OB4 += O2*C2[j];
      }
  }

  //  Examine fits
 
#ifdef DEBUG_PLOIDY
  { int    j;
    int64  r, v;

    for (j = first; j < PIXELS; j++)
      { r = row[j];
        printf(" %2d: %10lld",j,row[j]);
        v = (int64) (C2[j]*D2);
        printf(" %10lld (%10lld)",v,v-r);
        v = (int64) (C3[j]*T3);
        printf(" %10lld (%10lld)",v,v-r);
        v = (int64) (C4[j]*Q4+C2[j]*Q2);
        printf(" %10lld (%10lld)",v,v-r);
        v = (int64) (C2[j]*H2 + C3[j]*H3 + C6[j]*H6);
        printf(" %10lld (%10lld)",v,v-r);
        v = (int64) (C2[j]*O2 + C38[j]*O3 + C4[j]*O4 + C8[j]*O8);
        printf(" %10lld (%10lld)",v,v-r);
        printf("\n");
      }
  }
#endif

  Df = lsqr[0] = sqrt((1.*Df)/(PIXELS-first));
  Tf = lsqr[1] = sqrt((1.*Tf)/(PIXELS-first));
  Qf = lsqr[2] = sqrt((1.*Qf)/(PIXELS-first));
  Hf = lsqr[3] = sqrt((1.*Hf)/(PIXELS-first));
  Of = lsqr[4] = sqrt((1.*Of)/(PIXELS-first));

  DB1 = Deco[0][0] = 1000;
  TB1 = Deco[1][0] = 1000;
  QB1 = Deco[2][0] = (1000.*QB1)/Qs;
  QB2 = Deco[2][1] = (1000.*QB2)/Qs;
  HB1 = Deco[3][0] = (1000.*HB1)/Hs;
  HB2 = Deco[3][1] = (1000.*HB2)/Hs;
  HB3 = Deco[3][2] = (1000.*HB3)/Hs;
  OB1 = Deco[4][0] = (1000.*OB1)/Os;
  OB2 = Deco[4][1] = (1000.*OB2)/Os;
  OB3 = Deco[4][2] = (1000.*OB3)/Os;
  OB4 = Deco[4][3] = (1000.*OB4)/Os;

#ifdef ANALYZE_PLOIDY
  printf("\n%10lld : %10lld : AB 100%%\n",Df,Ds);
  printf("%10lld : %10lld : AAB 100%%\n",Tf,Ts);
  printf("%10lld : %10lld : AAAB %.1f%% + AABB %.1f%%\n",Qf,Qs,QB1/10.,QB2/10.);
  printf("%10lld : %10lld : 6A1B %.1f%% + 4A2B %.1f%% + 3A3B %.1f%%\n",
         Hf,Hs,HB1/10.,HB2/10.,HB3/10.);
  printf("%10lld : %10lld : 7A1B %.1f%% + 6A2B %.1f%% + 5A3B %.1f%% + 4A4B %.1f%%\n",
         Of,Os,OB1/10.,OB2/10.,OB3/10.,OB4/10.);
#endif

  { int b, i;

    b = 0;
    for (i = 1; i <= 4; i++)
      if (lsqr[b] > lsqr[i])
        b = i; 

#ifdef ANALYZE_PLOIDY
    if (ploidy != b)
      { printf("Conflict %d vs %d\n",ploidy,b);
        if (ploidy < b)
          { if (Deco[b][0] < .5*Deco[ploidy][0])
              printf("  Peak choice seems OK\n");
            else
              printf("  Switch to greater fit choice\n");
          }
        else
          printf("  Peak choice is greater, so OK\n");
      }
    else
      printf("Fit & Ploidy match %d\n",b);
#endif

    if (ploidy < b)
      { if (Deco[b][0] >= .5*Deco[ploidy][0])
          ploidy = b;
      }
  }

  if (ploidy != 1)
    { int k = 48;
      while (row[k-1] > row[k] && k > 43)
        k -= 1;
      if (k > 43)
        Rail[ploidy][Nsmu[ploidy]-1] = k;
    }

  return (ploidy);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

static char template[15] = "._SPAIR.XXXX";

#ifdef SOLO_CHECK

static uint8 code[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static void compress_norm(char *s, int len, uint8 *t)
{ int    i;
  char   c, d, e;
  char  *s0, *s1, *s2, *s3;

  s0 = s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  c = s0[len];
  d = s1[len];
  e = s2[len];
  s0[len] = s1[len] = s2[len] = 'a';

  for (i = 0; i < len; i += 4)
    *t++ = ((code[(int) s0[i]] << 6) | (code[(int) s1[i]] << 4)
         |  (code[(int) s2[i]] << 2) | code[(int) s3[i]] );

  s0[len] = c;
  s1[len] = d;
  s2[len] = e;
}

#endif

static uint8 comp[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static void compress_comp(char *s, int len, uint8 *t)
{ int    i;
  char   c, d, e;
  char  *s0, *s1, *s2, *s3;

  s0 = s;
  s1 = s0-1;
  s2 = s1-1;
  s3 = s2-1;

  c = s1[0];
  d = s2[0];
  e = s3[0];
  s1[0] = s2[0] = s3[0] = 't';

  for (i = len-1; i >= 0; i -= 4)
    *t++ = ((comp[(int) s0[i]] << 6) | (comp[(int) s1[i]] << 4)
         |  (comp[(int) s2[i]] << 2) | comp[(int) s3[i]] );

  s1[0] = c;
  s2[0] = d;
  s3[0] = e;
}

static void examine_table(Kmer_Stream *T, int *trim, int *sym)
{
  //  Histogram of middle 100M counts and see if trimmed to ETHRESH

  { int64 hist[0x8000];
    int64 frst, last;
    int   hbyte = T->hbyte;
    int   i, nz;

    for (i = 0; i < 0x8000; i++)
      hist[i] = 0;

    if (T->nels+3 < 100000000)
      { frst = 0;
        last = T->nels;
      }
    else
      { frst = T->nels/2 - 50000000;
        last = T->nels/2 + 50000000;
      }

    for (GoTo_Kmer_Index(T,frst); T->cidx < last; Next_Kmer_Entry(T))
      hist[*((int16 *) (T->csuf+hbyte))] += 1;

    for (nz = 1; hist[nz] == 0; nz++)
      ;
    if (nz < ETHRESH)
      *trim = 0;
    else
      *trim = 1;
  }

  //  Walk to a non-palindromic k-mer and see if its complement is in T

  { int64  sidx;
    char  *seq;
    uint8 *cmp;
    int    kmer;

    kmer = T->kmer;

    sidx = 1;
    GoTo_Kmer_Index(T,sidx);
    seq = Current_Kmer(T,NULL);
    cmp = Current_Entry(T,NULL);
    while (1)
      { compress_comp(seq,kmer,cmp);
        if (GoTo_Kmer_Entry(T,cmp))
          { if (T->cidx != sidx)
              { *sym = 1;
                break;
              }
          }
        else
          { *sym = 0;
            break;
          }
        sidx += 1;
        seq = Current_Kmer(T,seq);
      }
    free(cmp);
    free(seq);
  }
}

int main(int argc, char *argv[])
{ Kmer_Stream *T;
  char        *input;
  char        *troot;
  int64      **PLOT;
  int          bypass;

  char  *SORT_PATH;
  int    KEEP;
  int    LINE, FILL, BOTH;
  int    PDF;
  double XDIM, YDIM;
  char  *OUT;
  char  *SRC;

  //  Process command line arguments

  (void) print_hap;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("PloidyPlot");
    XDIM = 6.0;
    YDIM = 4.5;
    PDF  = 0;
    OUT  = NULL;
    ETHRESH  = 4;
    NTHREADS = 4;
    SORT_PATH = "/tmp";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vklfs")
            break;
          case 'e':
            ARG_POSITIVE(ETHRESH,"Error-mer threshold")
            break;
          case 'h':
            ARG_REAL(YDIM);
            break;
          case 'o':
            if (OUT == NULL)
              free(OUT);
            OUT = Strdup(argv[i]+2,"Allocating name");
            if (OUT == NULL)
              exit (1);
            break;
          case 'p':
            if (strcmp("df",argv[i]+2) == 0)
              PDF = 1;
            else
              { fprintf(stderr,"%s: don't recognize option %s\n",Prog_Name,argv[i]);
                exit (1);
              }
            break;
          case 'w':
            ARG_REAL(XDIM);
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            if (NTHREADS > 64)
              { fprintf(stderr,"%s: Warning, only 64 threads will be used\n",Prog_Name);
                NTHREADS = 64;
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    KEEP    = flags['k'];
    LINE    = flags['l'];
    FILL    = flags['f'];
    BOTH    = flags['s'];

    if (LINE+FILL+BOTH == 0)
      LINE = FILL = BOTH = 1;

#ifdef SOLO_CHECK
    if (argc != 3)
#else
    if (argc != 2)
#endif
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -w: width in inches of plots\n");
        fprintf(stderr,"      -h: height in inches of plots\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: draw line plot\n");
        fprintf(stderr,"      -f: draw fill plot\n");
        fprintf(stderr,"      -s: draw stack plot\n");
        fprintf(stderr,"          any combo allowed, none => draw all\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"    -pdf: output .pdf (default is .png)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: root name for output plots\n");
        fprintf(stderr,"          default is root path of <asm> argument\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: count threshold below which k-mers are considered erroneous\n");
        fprintf(stderr,"      -v: verbose mode\n");
        fprintf(stderr,"      -k: keep het-mer table for re-use\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        exit (1);
      }

    SRC = argv[1];
    if (OUT == NULL)
      OUT = Root(argv[1],".ktab");

    troot = mktemp(template);
  }

  //  If appropriately named het-mer table found then ask if reuse

  { FILE *f;
    int   a;

    bypass = 0;
    f = fopen(Catenate(OUT,".smu","",""),"r");
    if (f != NULL)
      { fprintf(stdout,"\n  Found het-table %s.smu, use it? ",OUT);
        fflush(stdout);
        while ((a = getc(stdin)) != '\n')
          if (a == 'y' || a == 'Y')
            bypass = 1;

        if (bypass)
          { PLOT    = Malloc(sizeof(int64 *)*(SMAX+1),"Allocating thread working memory");
            PLOT[0] = Malloc(sizeof(int64)*(SMAX+1)*(FMAX+1),"Allocating plot");
            for (a = 1; a <= SMAX; a++)
              PLOT[a] = PLOT[a-1] + (FMAX+1);
            fread(PLOT[0],sizeof(int64),(SMAX+1)*(FMAX+1),f);
          }

        fclose(f);
      }
  }

  //  Open input table and see if it needs conditioning

  { char *command;
    char *tname;
    int   symm, trim;

    tname   = Malloc(strlen(SRC) + strlen(troot) + 10,"Allocating strings");
    command = Malloc(strlen(SRC) + strlen(troot) + 100,"Allocating strings");
    if (tname == NULL || command == NULL)
      exit (1);

    T = Open_Kmer_Stream(SRC);
    if (T == NULL)
      { fprintf(stderr,"%s: Cannot open k-mer table %s\n",Prog_Name,SRC);
        exit (1);
      }

    KMER  = T->kmer;
    KBYTE = T->kbyte;
    TBYTE = T->tbyte;

    if (bypass)
      goto skip_build;

    examine_table(T,&trim,&symm);

    if (VERBOSE)
      { fprintf(stderr,"\n  The input table is");
        if (trim)
          if (symm)
            fprintf(stderr," trimmed and symmetric\n");
          else
            fprintf(stderr," trimmed but not symmetric\n");
        else
          if (symm)
            fprintf(stderr," untrimmed yet symmetric\n");
          else
            fprintf(stderr," untrimmed and not symmetric\n");
      }

    sprintf(tname,"%s",SRC);
    input = NULL;

    //  Trim source table to k-mers with counts >= ETHRESH if needed

    if (!trim)
      { if (VERBOSE)
          { fprintf(stderr,"\n  Trimming k-mers in table with count < %d\n",ETHRESH);
            fflush(stderr);
          }

        sprintf(command,"Logex -T%d '%s.trim=A[%d-]' %s",NTHREADS,troot,ETHRESH,tname);
        if (system(command) != 0)
          { fprintf(stderr,"%s: Something went wrong with command:\n    %s\n",Prog_Name,command);
            exit (1);
          }

        sprintf(tname,"%s.trim",troot);
      }

    //  Make the (relevant) table symmetric if it is not

    if (!symm)
      { if (VERBOSE)
          { if (trim)
              fprintf(stderr,"\n  Making table symmetric\n");
            else
              fprintf(stderr,"\n  Making trimmed table symmetric\n");
            fflush(stderr);
          }

        sprintf(command,"Symmex -T%d -P%s %s %s.symx",NTHREADS,SORT_PATH,tname,troot);

        if (system(command) != 0)
          { fprintf(stderr,"%s: Something went wrong with command:\n    %s\n",Prog_Name,command);
            exit (1);
          }

        if (!trim)
          { sprintf(command,"Fastrm %s.trim",troot);
            if (system(command) != 0)
              { fprintf(stderr,"%s: Something went wrong with command:\n    %s\n",
                               Prog_Name,command);
                exit (1);
              }
          }

        sprintf(tname,"%s.symx",troot);
      }

    //  input is the name of the relevant conditioned table, unless the original => NULL

    free(command);
    if (!(symm && trim))
      { input = tname;
        Free_Kmer_Stream(T);
        T = Open_Kmer_Stream(input);
      }
    else
      free(tname);
  }

  if (VERBOSE)
    { fprintf(stderr,"\n  Starting to count covariant pairs\n");
      fflush(stderr);
    }

  BLEVEL = 1;
  BWIDTH = 4;
  while (4*NTHREADS > BWIDTH)
    { BWIDTH *= 4;
      BLEVEL += 1;
    }

  Cache_Size = (MEMORY_LIMIT/NTHREADS)/TBYTE;

#ifdef SOLO_CHECK
  if ((int) strlen(argv[3]) != KMER)
    { fprintf(stderr,"%s: string is not of length %d\n",Prog_Name,KMER);
      exit (1);
    }
  CENT = Current_Entry(T,NULL);
  compress_norm(argv[3],KMER,CENT);
  if (GoTo_Kmer_Entry(T,CENT) < 0)
    { fprintf(stderr,"%s: string is not in table\n",Prog_Name);
      exit (1);
    }
  printf("%s: %d\n",argv[3],Current_Count(T));
  CIDX = T->cidx;
#endif

#ifdef DEBUG_GENERAL
  printf("Threads = %d BL = %d(%d) Cache = %lld\n",NTHREADS,BLEVEL,BWIDTH,Cache_Size);
#endif

  { TP      parm[NTHREADS];
    int     a, t;
    int64 **plot;
    uint8  *cache;

    for (t = 0; t < NTHREADS; t++)
      { plot    = Malloc(sizeof(int64 *)*(SMAX+1),"Allocating thread working memory");
        plot[0] = Malloc(sizeof(int64)*(SMAX+1)*(FMAX+1),"Allocating plot");
        for (a = 1; a <= SMAX; a++)
          plot[a] = plot[a-1] + (FMAX+1);
        bzero(plot[0],sizeof(int64)*(SMAX+1)*(FMAX+1));
        parm[t].plot = plot;
      }

    parm[0].fng[0] = T;
    for (t = 0; t < NTHREADS; t++)
      for (a = 0; a < 4; a++)
        if (a+t > 0)
          parm[t].fng[a] = Clone_Kmer_Stream(T);

    Divpt = Current_Entry(T,NULL);
    Pair  = Malloc(sizeof(uint8)*T->nels,"Allocating pair table");
    cache = Malloc(MEMORY_LIMIT,"Allocating cache buffer");
    if (Pair == NULL || cache == NULL)
      exit (1);

    bzero(Pair,T->nels);

    for (PASS1 = 1; PASS1 >= 0; PASS1--)
      {
        //  Analyze the top levels, each threaded, up to level, BLEVEL, where
        //    the number of nodes is greater than the number of threads
        //    by a factor of 4 or more

        { int64 adiv[5];

          for (t = 0; t < KBYTE; t++)
            Divpt[t] = 0;

          adiv[0] = 0;
          for (a = 1; a < 4; a++)
            { Divpt[0] = (a << 6);
              GoTo_Kmer_Entry(T,Divpt);
              adiv[a] = T->cidx;
            }
          adiv[4] = T->nels;

          big_window(adiv,0,parm);
        }

        //  Assign a thread to each subtree at level BLEVEL until all are done

        { pthread_t threads[NTHREADS];
          int       tstack[NTHREADS];
          int       x;

          Tstack = tstack;

          for (t = 0; t < NTHREADS; t++)
            { Tstack[t]   = t;
              parm[t].tid = t;
              parm[t].cache = cache + t*(MEMORY_LIMIT/NTHREADS);
            }
          Tavail = NTHREADS;

          pthread_mutex_init(&TMUTEX,NULL);
          pthread_cond_init(&TCOND,NULL);

          for (x = 0; x < BWIDTH; x++)
            { pthread_mutex_lock(&TMUTEX);

              if (Tavail <= 0)
                pthread_cond_wait(&TCOND,&TMUTEX);

              t = Tstack[--Tavail];

#ifdef DEBUG_GENERAL
              printf("Launching %d on thread %d\n",x,t);
#endif

              pthread_mutex_unlock(&TMUTEX);

              parm[t].level = x;

              pthread_create(threads+t,NULL,small_window,parm+t);
            }

          pthread_mutex_lock(&TMUTEX);
          while (Tavail < NTHREADS)
            pthread_cond_wait(&TCOND,&TMUTEX);
	  pthread_mutex_unlock(&TMUTEX);
        }
      }

    free(cache);
    free(Pair);
    free(Divpt);

    { char  *command;
      int64 *plot0, *plott;
      int    i;

      for (t = NTHREADS-1; t >= 0; t--)
        for (a = 3; a >= 0; a--)
          if (a+t > 0)
            Free_Kmer_Stream(parm[t].fng[a]);
      Free_Kmer_Stream(T);

      for (t = 1; t < NTHREADS; t++)
        for (i = 0; i <= SMAX; i++)
          { plot0 = parm[0].plot[i];
            plott = parm[t].plot[i];
            for (a = 0; a <= FMAX; a++)
              plot0[a] += plott[a];
          }

      for (t = NTHREADS-1; t >= 1; t--)
        { free(parm[t].plot[0]);
          free(parm[t].plot);
        }

      PLOT = parm[0].plot;

      if (input != NULL)
        { command = Malloc(strlen(input)+100,"Allocating strings");
          if (command == NULL)
            exit (1);
          sprintf(command,"Fastrm %s",input);
          if (system(command) != 0)
            { fprintf(stderr,"%s: Something went wrong with command:\n    %s\n",Prog_Name,command);
              exit (1);
            }
          free(command);
          free(input);
        }
    }
  }

#ifdef SOLO_CHECK

  exit (0);

#endif

  if (VERBOSE)
    { fprintf(stderr,"\n  Count complete, plotting\n");
      fflush(stderr);
    }

skip_build:

#ifdef KAMIL

fprintf(stderr,"\n  About to save stuff\n");

FILE  *f;
int    a, i;

f = fopen(Catenate(OUT,"_text.smu","",""),"w");
fprintf(stderr,"\n  Saving stuff\n");

// fprintf(f, "// %dx%d matrix, the i'th number in the j'th row give the number of hetmer pairs (a,b)\n", SMAX,FMAX);
// fprintf(f, "//                     s.t. count(a)+count(b) = j+1 and min(count(a),count(b)) = i+1.\n");
for (a = 0; a <= SMAX; a++)
  {
    for (i = 0; i < FMAX; i++)
      if(PLOT[a][i] > 0)
      {
        fprintf(f,"%i\t%i\t%lld\n",i,a-i,PLOT[a][i]);
      }
  }
fclose(f);

#else

  if (KEEP && !bypass)
    { FILE  *f;

      f = fopen(Catenate(OUT,".smu","",""),"w");
      fwrite(PLOT[0],sizeof(int64),(SMAX+1)*(FMAX+1),f);
      fclose(f);
    }

  //  Create smudge table

  { int    i, j, a;
    int    p, q, r, w;
    int64  low, hgh;
    int    nov, smax, wide;
    int64  inter[SMAX];
    int64 *row, v, vmax;
    double fact;
    FILE  *f;
    char  *command, *capend;
    int64  rsum[SMAX];
    int64  cmax;
    int    cwch, cover, ploidy;

    //  Begin table output

    f = fopen(Catenate(troot,".smu","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".smu","",""));
    else
      printf("Writing %s\n",Catenate(troot,".smu","",""));
    fflush(stdout);
#endif
    fprintf(f,"KF1\tKF2\tCount\n");

    //  Stretch rows to fractional axis in place, scale to preserve row sums!
    //    Row goes from [0,FMAX) to [0,PIXEL)

    for (i = 2*ETHRESH; i < SMAX; i++)
      { int64  old, new;
        double imp;

        row = PLOT[i];

        for (p = 0; p < PIXELS; p++)
          inter[p] = -1;

        old = 0;
        for (a = ETHRESH; a <= i/2; a++)
          { if (i%2 == 0)
              p = (PIXEL2*a)/i;
            else
              p = (PIXEL2*a)/(i-1);
            if (p >= PIXELS)
              p = PIXELS-1;
	    if (inter[p] < 0)
              inter[p] = row[a];
            else
              inter[p] += row[a];
            old += row[a];
          }

        new = 0;
        nov = 1;
        for (p = 0; p < PIXELS; p++)
          if (inter[p] < 0)
            { if (nov)
                row[p] = -1;
              else
                { r = p-1;
                  q = p;
                  while (q < PIXELS && inter[q] < 0)
                    q ++;
                  if (q == PIXELS)
                    hgh = low;
                  else
                    hgh = inter[q];
                  low = inter[r];
                  w = q-r;
                  while (p < q)
                    { row[p] = (low*(q-p) + hgh*(p-r))/w;
                      new += row[p];
                      p += 1;
                    }
                  p -= 1;
                }
            }
          else
            { row[p] = inter[p];
              new += row[p];
              nov = 0;
            }

        rsum[i] = old;
        if (new == 0)
          imp = 1.;
        else
          imp = (1.*old)/new;

        for (p = 0; p < PIXELS; p++)
          if (row[p] > 0)
            row[p] *= imp;
      }

    //  Determine A+B scaling and max in wide/smax

    vmax = 0;
    for (i = 2*ETHRESH; i < SMAX; i++)
      if (vmax < rsum[i])
        { vmax = rsum[i];
          smax = i;
        }
    for (i = SMAX-1; i > 2*smax; i--)
      { if (rsum[i] > .02*vmax)
          break;
      }
    wide = (i-1)/PIXELS+1;
    smax = PIXELS*wide;

    //  reduce A+B resolution by scaling factor and output array

    a = wide;
    while (a <= 2*ETHRESH)
      a += wide;
    for (i = 2*ETHRESH; i < smax; a += wide)
      { fact = wide/(a-i);
        for (p = 0; p < PIXELS; p++)
          { v = -1;
            for (j = i; j < a; j++)
              if (PLOT[j][p] >= 0)
                { if (v < 0)
                    v = PLOT[j][p];
                  else
                    v += PLOT[j][p];
                }
            if (v >= 0)
              { fprintf(f,"%d.5 %d.5 %lld\n",i/wide,p,(int64) (v*fact));
                PLOT[i/wide][p] = v*fact;
              }
            else
              PLOT[i/wide][p] = 0;
          }
        i = a;
      }

    fclose(f);

    //  determine row cwch with max count in reduced space

    cmax = 0;
    cwch = -1;
    for (i = (2*ETHRESH)/wide; i < PIXELS; i++)
      { v = 0;
        for (p = 0; p < PIXELS; p++)
          v += PLOT[i][p];
        if (v > cmax)
          { cmax = v;
            cwch = i;
          }
      }

    //  refine coverage to row in original spcae with max count (among rows comprising cwch)

    cover = cwch*wide;
    for (i = cwch*wide+1; i < cwch*wide + (cwch-1); i++)
      if (rsum[i] > rsum[cover])
        cover = i;

#ifdef ANALYZE_PLOIDY
    printf("\nCoverage is %d\n",cover);
#endif

    //  ploidy analysis on the row cwch

    ploidy = determine_ploidy(cover, PLOT[cwch]);

    //  output analysis for R script

    f = fopen(Catenate(troot,".sma","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".smu","",""));
    else
      printf("Writing %s\n",Catenate(troot,".smu","",""));
    fflush(stdout);
#endif
    fprintf(f,"Name\tx\ty\tAmount\tCover\tPloidy\n");
    for (i = 0; i < Nsmu[ploidy]; i++)
      fprintf(f,"%s\t%d\t%d\t%.1f\t%d\t%d\n",
                Names[ploidy][i],Rail[ploidy][i],cwch,Deco[ploidy][i]/10.,cover,ploidy);
    fclose(f);

    free(PLOT[0]);
    free(PLOT);

    //  Generate the R plot script in another temp file

    f = fopen(Catenate(troot,".R","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".R","",""));
    else
      printf("Generating %s\n",Catenate(troot,".R","",""));
    fflush(stdout);
#endif
    fwrite(smu_plot,strlen(smu_plot),1,f);
    fclose(f);

    //  Call the R plotter with arguments

   command = Malloc(strlen(troot)*3 + strlen(OUT)*2 + 500,"Allocating strings");
   if (command == NULL)
     exit (1);

    sprintf(command,"Rscript %s.R -f %s.smu -a %s.sma -o %s%s -x %g -y %g -s %s",
                    troot,troot,troot,OUT,PDF?" -p":" ",XDIM,YDIM,OUT);
    capend = command+strlen(command);
    if (LINE)
      { sprintf(capend," -t contour 2>/dev/null");
#ifdef DEBUG
        printf("%s\n",command);
        fflush(stdout);
#endif
        system(command);
      }
    if (FILL)
      { sprintf(capend," -t heat 2>/dev/null");
#ifdef DEBUG
        printf("%s\n",command);
        fflush(stdout);
#endif
        system(command);
      }
    if (BOTH)
      { sprintf(capend," -t combo 2>/dev/null");
#ifdef DEBUG
        printf("%s\n",command);
        fflush(stdout);
#endif
        system(command);
      }

    //  Remove the temp files

    sprintf(command,"rm -f %s.smu %s.sma %s.R",troot,troot,troot);
    system(command);

    if ( ! KEEP)
      { sprintf(command,"rm -f %s.smu",OUT);
        system(command);
      }

    free(command);

    if (VERBOSE)
      { fprintf(stderr,"\nAnalysis summary:\n");
        fprintf(stderr,"  k = %d\n",KMER);
        fprintf(stderr,"  p = %d\n",Nloidy[ploidy]);
        fprintf(stderr,"  ploidy = %s\n",Ploidy[ploidy]);
        fprintf(stderr,"  1n = %.1f\n",(1.*cover)/Nloidy[ploidy]);
        fprintf(stderr,"  partition =");
        for (i = 0; i < Nsmu[ploidy]; i++)
          { if (i > 0)
              fprintf(stderr,", ");
            fprintf(stderr," %s:%.1f%%",Names[ploidy][i],Deco[ploidy][i]/10.);
          }
        fprintf(stderr,"\n");
      }
  }

#endif  //  KAMIL

  free(OUT);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
