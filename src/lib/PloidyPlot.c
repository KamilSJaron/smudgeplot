/******************************************************************************************
 *
 * PloidyPlot: A C backend searching quickly for hetmers:
 *                  unique k-mer pairs different by exactly one nucleotide
 *
 *  Author:  Gene Myers
 *  Date  :  May, 2021
 *  Reduced to the k-mer pair search by Kamil Jaron, August 2023
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

#undef  SOLO_CHECK

#undef  DEBUG_GENERAL
#undef  DEBUG_RECURSION
#undef  DEBUG_THREADS
#undef  DEBUG_BOUNDARY
#undef  DEBUG_SCAN
#undef  DEBUG_BIG_SCAN

#include "libfastk.h"
#include "matrix.h"

static char *Usage[] = { " [-v] [-T<int(4)>] [-P<dir(/tmp)>]",
                         " [-o<output>] [-e<int(4)>] <source>[.ktab]"
                       };

static int VERBOSE;
static int NTHREADS;   //  At most 64 allowed
static int ETHRESH;

#ifdef SOLO_CHECK

static uint8 *CENT;
static int64  CIDX;

#endif

#define SMAX  1000    //  Max. value of CovA+CovB
#define FMAX   500    //  Max. value of min(CovA,CovB)

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
    { ent[a] = NULL;
      if (fng[a]->cidx <= 0)
        lst[a] = 0;
      else
        { GoTo_Kmer_Index(fng[a],fng[a]->cidx-1);
          ent[a] = Current_Entry(fng[a],NULL);
          if (((ent[a][level>>2] >> Shift[level&0x3]) & 0x3) == a)
            lst[a] = (ent[a][ll] >> ls) & 0x3;
          else
            lst[a] = 0;
          Next_Kmer_Entry(fng[a]);
        }
      if (fng[a]->cidx < end[a])
        { ent[a] = Current_Entry(fng[a],ent[a]);
          x = (ent[a][ll] >> ls) & 0x3;
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
    { ent[a] = NULL;
      if (fng[a]->cidx <= 0)
        lst[a] = 0;
      else
        { GoTo_Kmer_Index(fng[a],fng[a]->cidx-1);
          ent[a] = Current_Entry(fng[a],NULL);
          if (((ent[a][level>>2] >> Shift[level&0x3]) & 0x3) == a)
            lst[a] = (ent[a][ll] >> ls) & 0x3;
          else
            lst[a] = 0;
          Next_Kmer_Entry(fng[a]);
        }
      if (fng[a]->cidx < end[a])
        { ent[a] = Current_Entry(fng[a],ent[a]);
          x = (ent[a][ll] >> ls) & 0x3;
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
    int          b, t, ls;
    Kmer_Stream *T;
#ifndef DEBUG_THREADS
    pthread_t    threads[NTHREADS];
#endif

    b = 0;
    for (a = 1; a < 4; a++)
      if (adiv[a+1]-adiv[a] > adiv[b+1]-adiv[b])
        b = a;
    if (adiv[b+1]-adiv[b] == 0)
      return;

    ls = Shift[level];      //  level < BLEVEL
    lm = 0xff ^ (0x3<<ls);
    for (a = 0; a < 4; a++)
      GoTo_Kmer_Index(parm[0].fng[a],adiv[a]);
    parm[0].level = level;
    parm[0].bound = bound;
    for (t = 1; t < NTHREADS; t++)
      { parm[t].level = level;
        parm[t].bound = bound;
        parm[t-1].end[b] = e = adiv[b] + t*(adiv[b+1]-adiv[b])/NTHREADS;
        T = parm[t].fng[b];
        GoTo_Kmer_Index(T,e);
        Current_Entry(T,Divpt);
        for (a = 0; a < 4; a++)
          if (a != b)
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
          if (parm[t].fng[a]->cidx < parm[t].fng[a]->nels)
            { Current_Entry(parm[t].fng[a],Divpt);
              print_hap(Divpt,KMER,level);
            }
          else
            printf(" EOT");
          printf("\n");
        }
#endif

    for (a = 0; a < 16; a += 4)
      { bound[a] = adiv[a>>2];
        for (t = 1; t < 4; t++)
          bound[a+t] = -1;
      }
    bound[16] = adiv[4];

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

    for (a = 15; a > 0; a--)
      if (bound[a] < 0)
        bound[a] = bound[a+1];

#ifdef DEBUG_BOUNDARY
    T = parm[0].fng[0];
    for (a = 0; a <= 16; a++)
      { if (a > 0)
          { GoTo_Kmer_Index(T,bound[a]-1);
            if (T->cidx < T->nels)
              { Current_Entry(T,Divpt);
                printf("%c %10lld: ",dna[a&0x3],T->cidx);
                print_hap(Divpt,KMER,level+1);
              }
            else
              printf("%c %10lld: EOT",dna[a&0x3],T->cidx);
            printf("\n");
          }
        if (a < 16)
          { GoTo_Kmer_Index(T,bound[a]);
            if (T->cidx < T->nels)
              { Current_Entry(T,Divpt);
                printf("%c %10lld: ",dna[a&0x3],T->cidx);
                print_hap(Divpt,KMER,level+1);
              }
            else
              printf("%c %10lld: EOT",dna[a&0x3],T->cidx);
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

    for (a = 15; a > 0; a--)
      if (bound[a] == NULL)
        bound[a] = bound[a+1];

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
      bound[16] = adiv[4];

      if (PASS1)
        analysis_thread_1(parm);
      else
        analysis_thread_2(parm);

      for (a = 15; a > 0; a--)
        if (bound[a] < 0)
          bound[a] = bound[a+1];

#ifdef DEBUG_BOUNDARY
      { Kmer_Stream *T;

        T = parm->fng[0];
        if (level+1 < KMER)
          for (a = 0; a <= 16; a++)
            { if (a > 0)
                { GoTo_Kmer_Index(T,bound[a]-1);
                  if (T->cidx < T->nels)
                    { Current_Entry(T,Divpt);
                      printf("%c %10lld: ",dna[a&0x3],T->cidx);
                      print_hap(Divpt,KMER,level+1);
                    }
                  else
                    printf("%c %10lld: EOT",dna[a&0x3],T->cidx);
                  printf("\n");
                }
              if (a < 16)
                { GoTo_Kmer_Index(T,bound[a]);
                  if (T->cidx < T->nels)
                    { Current_Entry(T,Divpt);
                      printf("%c %10lld: ",dna[a&0x3],T->cidx); 
                      print_hap(Divpt,KMER,level+1);
                    }
                  else
                    printf("%c %10lld: EOT",dna[a&0x3],T->cidx);
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

  char  *SORT_PATH;
  char  *OUT;
  char  *SRC;

  //  Process command line arguments

  (void) print_hap;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("hetmers");

    OUT      = NULL;
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
          case 'o':
            if (OUT != NULL)
              free(OUT);
            OUT = Strdup(argv[i]+2,"Allocating output name");
            if (OUT == NULL)
              exit (1);
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

#ifdef SOLO_CHECK
    if (argc != 3)
#else
    if (argc != 2)
#endif
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: root name for output table\n");
        fprintf(stderr,"            default is root of <source> argument\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: count threshold below which k-mers are considered erroneous\n");
        fprintf(stderr,"      -v: verbose mode\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        exit (1);
      }

    SRC = argv[1];
    if (OUT == NULL)
      OUT = PathnRoot(argv[1],".ktab");

    troot = mktemp(template);
  }

  //  If appropriately named het-mer table found then ask if reuse

  { FILE *f;
    int   a;

    f = fopen(Catenate(OUT,".smu","",""),"r");
    if (f != NULL)
      { int bypass;

        bypass = 0;
        fprintf(stdout,"\n  Found het-table %s.smu, use it? ",OUT);
        fflush(stdout);
        while ((a = getc(stdin)) != '\n')
          if (a == 'y' || a == 'Y')
            bypass = 1;
        if (bypass)
          { fprintf(stderr,"\n  Using the found het-table, done\n");
            fclose(f);
            exit (0);
          }
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
        SystemX(command);

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

        SystemX(command);

        if (!trim)
          { sprintf(command,"Fastrm %s.trim",troot);
            SystemX(command);
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
  if ((int) strlen(argv[2]) != KMER)
    { fprintf(stderr,"%s: string is not of length %d\n",Prog_Name,KMER);
      exit (1);
    }
  CENT = Current_Entry(T,NULL);
  compress_norm(argv[2],KMER,CENT);
  if (GoTo_Kmer_Entry(T,CENT) < 0)
    { fprintf(stderr,"%s: string is not in table\n",Prog_Name);
      exit (1);
    }
  printf("%s: %d\n",argv[2],Current_Count(T));
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
          SystemX(command);
          free(command);
          free(input);
        }
    }
  }

#ifndef SOLO_CHECK

  if (VERBOSE)
    { fprintf(stderr,"\n  Count complete, outputting table\n");
      fflush(stderr);
    }

  { int   a, i;
    FILE *f;

    f = fopen(Catenate(OUT,".smu","",""),"w");
    if (f == NULL)
      { fprintf(stderr,"Could not open %s.smu\n",troot);
        exit (1);
      }

    for (a = 0; a <= SMAX; a++)
      for (i = 0; i < FMAX; i++)
        if (PLOT[a][i] > 0)
          fprintf(f,"%i\t%i\t%lld\n",i,a-i,PLOT[a][i]);
    fclose(f);
  }

#endif

  free(PLOT[0]);
  free(PLOT);
  free(OUT);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
