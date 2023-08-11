/*******************************************************************************************
 *
 *  C library routines to access and operate upon FastK histogram, k-mer tables, and profiles
 *
 *  Author:  Gene Myers
 *  Date  :  November 2020
 *
 *******************************************************************************************/

#include "libfastk.h"

#include "gene_core.c"

/*********************************************************************************************\
 *
 *  HISTOGRAM CODE
 *
 *********************************************************************************************/

//  Toggle histogram from unique to instance counts or vice versa

static void toggle_histogram(Histogram *H)
{ int64 *hist = H->hist;
  int    low  = H->low;
  int    high = H->high;
  int64  x;
  int    i;

  if (H->unique)
    { for (i = low+1; i < high; i++)
        hist[i] *= i;
      H->unique = 0;
    }
  else
    { for (i = low+1; i < high; i++)
        hist[i] /= i;
      H->unique = 1;
    }

  x = hist[high+1];
  hist[high+1] = hist[low];
  hist[low] = x;

  x = hist[high+2];
  hist[high+2] = hist[high];
  hist[high] = x;
}

//  Read histogram encoded in 'name' and create and return Histogram instance of it

Histogram *Load_Histogram(char *name)
{ Histogram *H;
  int        kmer, low, high;
  int64      ilowcnt, ihighcnt;
  int64     *hist;
  char      *dir, *root, *full;
  int        f;

  dir  = PathTo(name);
  root = Root(name,".hist");
  full = Malloc(strlen(dir)+strlen(root)+10,"Histogram name allocation");
  if (full == NULL)
    exit (1);
  sprintf(full,"%s/%s.hist",dir,root);
  f = open(full,O_RDONLY);
  if (f < 0)
    return (NULL);
  free(full);
  free(root);
  free(dir);

  read(f,&kmer,sizeof(int));
  read(f,&low,sizeof(int));
  read(f,&high,sizeof(int));
  read(f,&ilowcnt,sizeof(int64));
  read(f,&ihighcnt,sizeof(int64));

  H    = Malloc(sizeof(Histogram),"Allocating histogram");
  hist = Malloc(sizeof(int64)*((high-low)+3),"Allocating histogram");
  if (H == NULL || hist == NULL)
    exit (1);

  read(f,hist,sizeof(int64)*((high-low)+1));
    
  close(f);

  H->kmer = kmer;
  H->low  = low;
  H->high = high;
  H->hist = hist = hist-low;
  H->unique = 1;
  hist[high+1] = ilowcnt;     // boundary counts for opposite mode hidden at top of histogram
  hist[high+2] = ihighcnt;

  return (H);
}

//  Modify histogram so its ragne is 'low'..'hgh' and it displays counts as specified by 'unique'

void Modify_Histogram(Histogram *H, int low, int high, int unique)
{ int64 *hist = H->hist;
  int64  under, over;
  int64  hunder, hover;
  int    i;
 
  if (H->low > low || H->high < high)
    return;

  under  = hist[H->low];
  over   = hist[H->high];
  for (i = H->low+1; i <= low; i++)
    under += hist[i];
  for (i = H->high-1; i >= high; i--)
    over += hist[i];

  hunder = hist[H->high+1];
  hover  = hist[H->high+2];
  if (H->unique)
    { for (i = H->low+1; i <= low; i++)
        hunder += hist[i]*i;
      for (i = H->high-1; i >= high; i--)
        hover += hist[i]*i;
    }
  else
    { for (i = H->low+1; i <= low; i++)
        hunder += hist[i]/i;
      for (i = H->high-1; i >= high; i--)
        hover += hist[i]/i;
    }

  if (low != H->low)
    memmove(H->hist+H->low,H->hist+low,((high-low)+1)*sizeof(int64));

  H->hist += H->low;
  H->hist  = Realloc(H->hist,((high-low)+3)*sizeof(int64),"Reallocating histogram");
  H->hist -= low;
  H->low   = low;
  H->high  = high;

  H->hist[low]    = under;
  H->hist[high]   = over;
  H->hist[high+1] = hunder;
  H->hist[high+2] = hover;

  if ((H->unique == 0) != (unique == 0))
    toggle_histogram(H);
}

//  Write histogram in FastK format to 'name'

int Write_Histogram(char *name, Histogram *H)
{ int64 *hist = H->hist;
  int    low  = H->low;
  int    high = H->high;
  char  *dir, *root, *full;
  int    f;

  if (H->unique == 0)
    toggle_histogram(H);

  dir  = PathTo(name);
  root = Root(name,".hist");
  full = Malloc(strlen(dir)+strlen(root)+10,"Histogram name allocation");
  if (full == NULL)
    exit (1);
  sprintf(full,"%s/%s.hist",dir,root);
  f = open(full,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
  if (f < 0)
    return (1);
  free(full);
  free(root);
  free(dir);

  write(f,&H->kmer,sizeof(int));
  write(f,&low,sizeof(int));
  write(f,&high,sizeof(int));
  write(f,hist+(high+1),sizeof(int64));
  write(f,hist+(high+2),sizeof(int64));
  write(f,hist+low,sizeof(int64)*((high-low)+1));
  close(f);

  if (H->unique == 0)
    toggle_histogram(H);

  return (0);
}

//  Free all memory for histogram

void Free_Histogram(Histogram *H)
{ free(H->hist+H->low);
  free(H);
}

/****************************************************************************************
 *
 *  K-MER TABLE CODE
 *
 *****************************************************************************************/

//  Private view of a Kmer_Table

typedef struct
  { int     kmer;         //  kmer length
    int     minval;       //  the minimum count of a k-mer in the table
    int64   nels;         //  # of unique, sorted k-mers in the table
                       // hidden fields
    int     ibyte;        //  # of prefix bytes
    int     kbyte;        //  kmer encoding in bytes (= ceiling(kmer/4))
    int     tbyte;        //  kmer+count entry in bytes (= kbyte + 2)
    int     hbyte;        //  kmer suffix in bytes (= kbyte - ibyte)
    int     pbyte;        //  kmer,count suffix in bytes (= tbyte - ibyte)
    int     ixlen;        //  length of prefix index (= 4^(4*ibyte))
    uint8  *table;        //  the (huge) table in memory
    int64  *index;        //  prefix compression index
    int    *inver;        //  inverse prefix index
    int     shift;        //  shift for inverse mapping
  } _Kmer_Table;

#define TABLE(T) ((_Kmer_Table *) T)

/****************************************************************************************
 *
 *  Basic compressed sequence utilities
 *
 *****************************************************************************************/

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
  
static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n--)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}

static int *inverse_index(int ixlen, int64 nels, int64 *index, int *pshift)
{ int64 step, pow;
  int   shift, inlen;
  int  *inver;
  int64 i, j, k;

  step = nels/ixlen;
  pow  = 2;
  for (shift = 0; pow <= step; shift++)
    pow <<= 1; 
  pow >>= 1;

  inlen = nels/pow;
  inver = (int *) Malloc(sizeof(int)*(inlen+1),"Allocating inverse prefix array");

  i = k = 0;
  for (j = 0; j < inlen; j++)
    { i = (j << shift); 
      while (index[k] <= i)
        k += 1;
      inver[j] = k;
    }
  inver[inlen] = ixlen-1;

  *pshift = shift;
  return (inver);
}


/****************************************************************************************
 *
 *  Loading and Freeing a Kmer_Table
 *
 *****************************************************************************************/

//  Reads over 2GB don't work on some systems, patch to overcome said

static inline int64 big_read(int f, uint8 *buffer, int64 bytes)
{ int64 v, x;

  v = 0;
  while (bytes > 0x70000000)
    { x = read(f,buffer,0x70000000);
      if (x < 0)
        return (-1);
      v += x;
      bytes  -= 0x70000000;
      buffer += 0x70000000;
    }
  x = read(f,buffer,bytes);
  if (x < 0)
    return (-1);
  return (v+x);
} 

//  Load table encoded in file 'name' and create Kmer_Table object of entries
//    with minimum count 'cut_off'

Kmer_Table *Load_Kmer_Table(char *name, int cut_off)
{ Kmer_Table  *T;
  Kmer_Stream *S;
  int          kmer, tbyte, kbyte, minval, ibyte, pbyte, hbyte;
  int64        nels;
  uint8       *table;
  int64       *index, ixlen;
  int         *inver, shift;

  int    f, flen;
  char  *dir, *root, *full;
  int    smer, nthreads;

  setup_fmer_table();

  //  Open stub file and get # of parts

  dir  = PathTo(name);
  root = Root(name,".ktab");
  full = Malloc(strlen(dir)+strlen(root)+20,"Histogram name allocation");
  if (full == NULL)
    exit (1);
  sprintf(full,"%s/%s.ktab",dir,root);
  f = open(full,O_RDONLY);
  sprintf(full,"%s/.%s.ktab.",dir,root);
  flen = strlen(full);
  free(root);
  free(dir);
  if (f < 0)
    { free(full);
      return (NULL);
    }

  read(f,&smer,sizeof(int));
  read(f,&nthreads,sizeof(int));
  read(f,&minval,sizeof(int));
  read(f,&ibyte,sizeof(int));

  kmer  = smer;
  kbyte = (kmer+3)>>2;
  tbyte = kbyte+2;
  pbyte = tbyte-ibyte;
  hbyte = kbyte-ibyte;
  ixlen = (1 << (8*ibyte));

  index = Malloc(ixlen*sizeof(int64),"Allocating table prefix index\n");
  if (index == NULL)
    exit (1);

  //  Find all parts and accumulate total size

  nels = 0;
  if (cut_off > minval)

    { bzero(index,ixlen*sizeof(int64));
      close(f);

      S = Open_Kmer_Stream(name);
      for (First_Kmer_Entry(S); S->csuf != NULL; Next_Kmer_Entry(S))
        if (Current_Count(S) >= cut_off)
          nels += 1;
      Free_Kmer_Stream(S);
    }

  else

    { int    p;
      int64  n;

      read(f,index,ixlen*sizeof(int64));
      close(f);

      for (p = 1; p <= nthreads; p++)
        { sprintf(full+flen,"%d",p);
          f = open(full,O_RDONLY);
          if (f < 0)
            { fprintf(stderr,"Table part %s is missing ?\n",full);
              exit (1);
            }
          read(f,&kmer,sizeof(int));
          read(f,&n,sizeof(int64));
          nels += n;
          if (kmer != smer)
            { fprintf(stderr,"Table part %s does not have k-mer length matching stub ?\n",
                             full);
              exit (1);
            }
          close(f);
        }
    }

  //  Allocate in-memory table

  T     = Malloc(sizeof(Kmer_Table),"Allocating table record");
  table = Malloc(nels*pbyte,"Allocating k-mer table\n");
  if ( T == NULL || table == NULL)
    exit (1);

  //  Load the table parts into memory

  if (cut_off > minval)

    { uint8 *jptr;
      int64  off;
      int    x;

      S = Open_Kmer_Stream(name);
      jptr = table;
      for (First_Kmer_Entry(S); S->csuf != NULL; Next_Kmer_Entry(S))
        if (Current_Count(S) >= cut_off)
          { mycpy(jptr,S->csuf,pbyte);
            jptr += pbyte;
            index[S->cpre] += 1;
          }
      Free_Kmer_Stream(S);

      off = 0;
      for (x = 0; x <  ixlen; x++)
        { off += index[x];
          index[x] = off;
        }

      minval = cut_off;
    }

  else

    { int    p;
      int64  n;
 
      nels = 0;
      for (p = 1; p <= nthreads; p++)
        { sprintf(full+flen,"%d",p);
          f = open(full,O_RDONLY);
          read(f,&kmer,sizeof(int));
          read(f,&n,sizeof(int64));
          big_read(f,table+nels*pbyte,n*pbyte);
          nels += n;
          close(f);
        }
    }

  free(full);

  inver = inverse_index(ixlen,nels,index,&shift);

  //  Finalize table record

  T->kmer   = kmer;
  T->minval = minval;
  T->nels   = nels;
  TABLE(T)->tbyte = tbyte;
  TABLE(T)->kbyte = kbyte;
  TABLE(T)->ibyte = ibyte;
  TABLE(T)->pbyte = pbyte;
  TABLE(T)->hbyte = hbyte;
  TABLE(T)->ixlen = ixlen;
  TABLE(T)->table = table;
  TABLE(T)->index = index;
  TABLE(T)->inver = inver;
  TABLE(T)->shift = shift;

  return (T);
}

//  Free all memory for table

void Free_Kmer_Table(Kmer_Table *T)
{ free(TABLE(T)->table);
  free(TABLE(T)->index);
  free(TABLE(T)->inver);
  free(T);
}


/****************************************************************************************
 *
 *  Fetch entry info
 *
 *****************************************************************************************/

  //  Asssumes i is in range

char *Fetch_Kmer(Kmer_Table *_T, int64 i, char *seq)
{ _Kmer_Table *T = TABLE(_T);
  int    hbyte = T->hbyte;
  int64 *index = T->index;
  int64  idx;

  if (seq == NULL)
    { seq = (char *) Malloc(T->kmer+3,"Reallocating k-mer buffer");
      if (seq == NULL)
        exit (1);
      if (T->nels == 0)
        return (seq);
    }

  idx = T->inver[i>>T->shift];
  while (index[idx] <= i)
    idx += 1;

  { int    j;
    uint8 *a;
    char  *s;

    s = seq;
    switch (T->ibyte)
    { case 3:
        memcpy(s,fmer[idx>>16],4);
        s += 4;
        memcpy(s,fmer[idx>>8 & 0xff],4);
        s += 4;
        memcpy(s,fmer[idx&0xff],4);
        s += 4;
        break;
      case 2:
        memcpy(s,fmer[idx>>8],4);
        s += 4;
        memcpy(s,fmer[idx&0xff],4);
        s += 4;
        break;
      case 1:
        memcpy(s,fmer[idx],4);
        s += 4;
        break;
    }

    a = T->table + i*T->pbyte;
    for (j = 0; j < hbyte; j++, s += 4)
      memcpy(s,fmer[a[j]],4);
    seq[T->kmer] = '\0';
  }

  return (seq);
}

  //  Asssumes i is in range

inline int Fetch_Count(Kmer_Table *T, int64 i)
{ return (*((uint16 *) (TABLE(T)->table+i*TABLE(T)->pbyte+TABLE(T)->hbyte))); }


/****************************************************************************************
 *
 *  Find k-mer in table
 *
 *****************************************************************************************/

static uint8 code[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static uint8 comp[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static int is_minimal(char *seq, int len)
{ int j, k;
  int x, y;
  
  for (k = 0, j = len-1; k < j; k++, j--)
    { x = code[(int) seq[k]];
      y = comp[(int) seq[j]];
      if (x < y)
        return (1);
      if (x > y)
        return (0);
    }
  if (k <= j)
    { x = code[(int) seq[k]];
      if (x < 2)
        return (1);
      else
        return (0);
    }
  else
    return (1);
}

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

int64 Find_Kmer(Kmer_Table *_T, char *kseq)
{ _Kmer_Table *T = (_Kmer_Table *) _T;
  int    kmer  = T->kmer;
  int    ibyte = T->ibyte;
  int    hbyte = T->hbyte;
  int    pbyte = T->pbyte;
  uint8 *table = T->table;
  int64 *index = T->index;

  uint8  cmp[T->kbyte], *c;
  int64  l, r, m, t;

  //  kseq must be at least kmer bp long

  if (is_minimal(kseq,kmer))
    compress_norm(kseq,kmer,cmp);
  else
    compress_comp(kseq,kmer,cmp);

  c = cmp;
  m = *c++;
  for (l = 1; l < ibyte; l++)
    m = (m << 8) | *c++;
  if (m == 0)
    l = 0;
  else
    l = index[m-1];
  if (l >= T->nels)
    return (-1);
  r = t = index[m];
  if (r <= l)
    return (-1);

  // smallest l s.t. KMER(l) >= (kmer) c  (or nels if does not exist)

  while (l < r)
    { m = ((l+r) >> 1);
      if (mycmp(table+m*pbyte,c,hbyte) < 0)
        l = m+1;
      else
        r = m;
    }

  if (l >= t || mycmp(table+l*pbyte,c,hbyte) != 0)
    return (-1);

  return (l);
}

/****************************************************************************************
 *
 *  K-MER STREAM CODE
 *
 *****************************************************************************************/

typedef struct
  { int    kmer;       //  Kmer length
    int    minval;     //  The minimum count of a k-mer in the stream
    int64  nels;       //  # of elements in entire table
                   //  Current position (visible part)
    int64  cidx;       //  current element index
    uint8 *csuf;       //  current element suffix
    int    cpre;       //  current element prefix
                   //  Other useful parameters
    int    ibyte;      //  # of bytes in prefix
    int    kbyte;      //  Kmer encoding in bytes
    int    tbyte;      //  Kmer+count entry in bytes
    int    hbyte;      //  Kmer suffix in bytes (= kbyte - ibyte)
    int    pbyte;      //  Kmer,count suffix in bytes (= tbyte - ibyte)
                   //  Hidden parts
    int    ixlen;      //  length of prefix index (= 4^(4*ibyte))
    int    shift;      //  shift for inverse mapping
    uint8 *table;      //  The (huge) table in memory
    int64 *index;      //  Prefix compression index
    int   *inver;      //  inverse prefix index
    int    copn;       //  File currently open
    int    part;       //  Thread # of file currently open
    int    nthr;       //  # of thread parts
    int    nlen;       //  length of path name
    char  *name;       //  Path name for table parts (only # missing)
    uint8 *ctop;       //  Ptr top of current table block in buffer
    int64 *neps;       //  Size of each thread part in elements
    int    clone;      //  Is this a clone?
  } _Kmer_Stream;

#define STREAM(S) ((_Kmer_Stream *) S)

#define STREAM_BLOCK 1024

/****************************************************************************************
 *
 *  Open a table and return as a Kmer_Stream object
 *
 *****************************************************************************************/

//  Load up the table buffer with the next STREAM_BLOCK suffixes (if possible)

static void More_Kmer_Stream(_Kmer_Stream *S)
{ int    pbyte = S->pbyte;
  uint8 *table = S->table;
  int    copn  = S->copn;
  uint8 *ctop;

  if (S->part > S->nthr)
    return;
  while (1)
    { ctop = table + read(copn,table,STREAM_BLOCK*pbyte);
      if (ctop > table)
        break;
      close(copn);
      S->part += 1;
      if (S->part > S->nthr)
        { S->csuf = NULL;
          return;
        }
      sprintf(S->name+S->nlen,"%d",S->part);
      copn = open(S->name,O_RDONLY);
      lseek(copn,sizeof(int)+sizeof(int64),SEEK_SET);
    }
  S->csuf = table;
  S->ctop = ctop;
  S->copn = copn;
}

Kmer_Stream *Open_Kmer_Stream(char *name)
{ _Kmer_Stream *S;
  int           kmer, tbyte, kbyte, minval, ibyte, pbyte, hbyte, ixlen;
  int64         nels;
  int           copn;
  int           shift;

  int    f, p;
  char  *dir, *root, *full;
  int    smer, nthreads;
  int64  n;

  setup_fmer_table();

  //  Open stub file and read header values

  dir  = PathTo(name);
  root = Root(name,".ktab");
  full = Malloc(strlen(dir)+strlen(root)+20,"Histogram name allocation");
  if (full == NULL)
    exit (1);
  sprintf(full,"%s/%s.ktab",dir,root);
  f = open(full,O_RDONLY);
  sprintf(full,"%s/.%s.ktab.",dir,root);
  free(root);
  free(dir);
  if (f < 0)
    { free(full);
      return (NULL);
    }
  read(f,&smer,sizeof(int));
  read(f,&nthreads,sizeof(int));
  read(f,&minval,sizeof(int));
  read(f,&ibyte,sizeof(int));

  //  Set size variables and allocate space for components

  kmer  = smer;
  kbyte = (kmer+3)>>2;
  tbyte = kbyte+2;
  pbyte = tbyte-ibyte;
  hbyte = kbyte-ibyte;
  ixlen = (1 << (8*ibyte));

  S        = Malloc(sizeof(_Kmer_Stream),"Allocating table record");
  S->name  = full;
  S->nlen  = strlen(full);
  S->table = Malloc(STREAM_BLOCK*pbyte,"Allocating k-mer buffer\n");
  S->neps  = Malloc(nthreads*sizeof(int64),"Allocating parts table of Kmer_Stream");
  S->index = Malloc(ixlen*sizeof(int64),"Allocating table prefix index\n");
  if (S == NULL || S->table == NULL || S->neps == NULL || S->index == NULL)
    exit (1);

  //  Read in index from stub and then close it

  read(f,S->index,ixlen*sizeof(int64));
  close(f);

  //  Read header of each part aaccumulating # of elements

  nels = 0;
  for (p = 1; p <= nthreads; p++)
    { sprintf(S->name+S->nlen,"%d",p);
      copn = open(S->name,O_RDONLY);
      if (copn < 0)
        { fprintf(stderr,"%s: Table part %s is missing ?\n",Prog_Name,S->name);
          exit (1);
        }
      read(copn,&kmer,sizeof(int));
      read(copn,&n,sizeof(int64));
      nels += n;
      S->neps[p-1] = nels;
      if (kmer != smer)
        { fprintf(stderr,"%s: Table part %s does not have k-mer length matching stub ?\n",
                         Prog_Name,S->name);
          exit (1);
        }
      close(copn);
    }

  //  Create inverse index and set all object parameters

  S->inver = inverse_index(ixlen,nels,S->index,&shift);

  S->kmer   = kmer;
  S->minval = minval;
  S->tbyte  = tbyte;
  S->kbyte  = kbyte;
  S->nels   = nels;
  S->ibyte  = ibyte;
  S->pbyte  = pbyte;
  S->ixlen  = ixlen;
  S->shift  = shift;
  S->hbyte  = hbyte;
  S->nthr   = nthreads;
  S->clone  = 0;

  //  Set position to beginning

  sprintf(S->name+S->nlen,"%d",1);
  copn = open(S->name,O_RDONLY);
  lseek(copn,sizeof(int)+sizeof(int64),SEEK_SET);

  S->copn  = copn;
  S->part  = 1;

  More_Kmer_Stream(S);

  S->cidx  = 0;

  if (S->cidx >= S->nels)
    { S->csuf = NULL;
      S->cpre = S->ixlen;
      S->part = S->nthr+1;
    }
  else
    { S->cpre  = 0;
      while (S->index[S->cpre] <= 0)
        S->cpre += 1;
    }

  return ((Kmer_Stream *) S);
}

Kmer_Stream *Clone_Kmer_Stream(Kmer_Stream *O)
{ _Kmer_Stream *S;
  int copn;

  S = Malloc(sizeof(_Kmer_Stream),"Allocating table record");
  if (S == NULL)
    exit (1);

  *S = *STREAM(O);
  S->clone = 1;

  S->table = Malloc(STREAM_BLOCK*STREAM(O)->pbyte,"Allocating k-mer buffer\n");
  S->name  = Malloc(S->nlen+20,"Allocating k-mer buffer\n");
  if (S->table == NULL || S->name == NULL)
    exit (1);
  strncpy(S->name,STREAM(O)->name,S->nlen);

  //  Set position to beginning

  sprintf(S->name+S->nlen,"%d",1);
  copn = open(S->name,O_RDONLY);
  lseek(copn,sizeof(int)+sizeof(int64),SEEK_SET);

  S->copn  = copn;
  S->part  = 1;

  More_Kmer_Stream(S);
  S->cidx  = 0;

  if (S->cidx >= S->nels)
    { S->csuf = NULL;
      S->cpre = S->ixlen;
      S->part = S->nthr+1;
    }
  else
    { S->cpre  = 0;
      while (S->index[S->cpre] <= 0)
        S->cpre += 1;
    }

  return ((Kmer_Stream *) S);
}

void Free_Kmer_Stream(Kmer_Stream *_S)
{ _Kmer_Stream *S = STREAM(_S);

  if (!S->clone)
    { free(S->neps);
      free(S->index);
      free(S->inver);
    }
  free(S->name);
  free(S->table);
  if (S->copn >= 0)
    close(S->copn);
  free(S);
}


/****************************************************************************************
 *
 *  Hidden, for 1-code production useful to now the size of the largest prefix bucket
 *
 *****************************************************************************************/

int Max_Bucket(Kmer_Stream *_S)
{ _Kmer_Stream *S = STREAM(_S);
  int64  pidx, *curr;
  int    i, s, max;

  pidx = 0;
  curr = S->index;
  max  = 0;
  for (i = 0; i < S->ixlen; i++)
    { s = curr[i] - pidx;
      if (s > max)
        max = s;
      pidx = curr[i];
    }
  return (max);
}


/****************************************************************************************
 *
 *  Hidden, specialized version of Open, Clone, & Free used by Fastmerge
 *
 *****************************************************************************************/

int Open_Kmer_Cache(Kmer_Stream *_S, int64 bidx, int64 eidx, int bpre, int epre,
                    char *where, uint8 *buffer, int buflen)
{ _Kmer_Stream *S = STREAM(_S);
  int   f, g;
  int   p, len;
  int64 beps, leps;
  int64 beg, end;
  char *name, *r;

  memmove(S->index,S->index+bpre,((epre-bpre)+1)*sizeof(int64));
  S->index = Realloc(S->index,((epre-bpre)+1)*sizeof(int64),"Reallocating prefix index");
  S->index -= bpre;

  if (S->part <= S->nthr)
    close(S->copn);

  name = Malloc(S->nlen+strlen(where)+10,"Reallocing table name");
  if (name == NULL)
    exit (1);

  r = rindex(S->name,'/');
  S->name[S->nlen] = '\0';
  sprintf(name,"%s/%scache",where,r+2);
  S->name[S->nlen] = '.';

  g = open(name,O_CREAT|O_TRUNC|O_RDWR,S_IRWXU);
  if (g < 0)
    return (1);

  p = 0;
  while (bidx >= S->neps[p])
    p += 1;

  beps = bidx;
  if (p == 0)
    leps = 0;
  else
    leps = S->neps[p-1];

  while (eidx > leps)
    { sprintf(S->name+S->nlen,"%d",p+1);
      f = open(S->name,O_RDONLY);
      lseek(f,sizeof(int)+sizeof(int64),SEEK_SET);

      if (eidx > S->neps[p])
        end = S->neps[p]-leps;
      else
        end = eidx-leps;
      beg = beps - leps;

      end *= S->pbyte;
      beg *= S->pbyte;
      if (beg > 0)
        lseek(f,beg,SEEK_CUR);

      while (beg < end)
        { if (beg+buflen > end)
            len = end-beg;
          else
            len = buflen;
          read(f,buffer,len);
          if (write(g,buffer,len) < 0)
            return (1);
          beg += len;
        }

      close(f);

      beps = leps = S->neps[p++];
    }

  free(S->name);
  S->name = name;

  lseek(g,0,SEEK_SET);
  S->copn = g;
  S->part = 0;
  S->cidx = bidx;
  S->cpre = bpre;
  S->nels = eidx;
  More_Kmer_Stream(S);

  return (0);
}

Kmer_Stream *Clone_Kmer_Cache(Kmer_Stream *O, int64 fidx, int64 bidx, int bpre)
{ _Kmer_Stream *S;
  int copn;

  S = Malloc(sizeof(_Kmer_Stream),"Allocating table record");
  if (S == NULL)
    exit (1);

  *S = *STREAM(O);
  S->clone = 1;

  S->table = Malloc(STREAM_BLOCK*S->pbyte,"Allocating k-mer buffer\n");

  //  Set position to beginning

  copn = open(S->name,O_RDONLY);
  lseek(copn,(bidx-fidx)*S->pbyte,SEEK_SET);

  S->copn  = copn;
  S->cidx  = bidx;
  S->cpre  = bpre;

  More_Kmer_Stream(S);

  return ((Kmer_Stream *) S);
}

void Free_Kmer_Cache(Kmer_Stream *_S, int bpre)
{ _Kmer_Stream *S = STREAM(_S);

  if (!S->clone)
    { free(S->neps);
      free(S->index + bpre);
      free(S->inver);
    }
  free(S->table);
  if (S->copn >= 0)
    close(S->copn);
  if (!S->clone)
    { unlink(S->name);
      free(S->name);
    }
  free(S);
}


/****************************************************************************************
 *
 *  Free, Iterate, and Get Stream Entries
 *
 *****************************************************************************************/

inline void First_Kmer_Entry(Kmer_Stream *_S)
{ _Kmer_Stream *S = STREAM(_S);
  int64 *index = S->index;

  if (S->cidx != 0)
    { if (S->part != 1)
        { if (S->part <= S->nthr)
            close(S->copn);
          sprintf(S->name+S->nlen,"%d",1);
          S->copn = open(S->name,O_RDONLY);
          S->part = 1;
        }

      lseek(S->copn,sizeof(int)+sizeof(int64),SEEK_SET);

      More_Kmer_Stream(S);
      S->cidx = 0;
      S->cpre = 0;
      while (index[S->cpre] <= 0)
        S->cpre += 1;
    }
}

inline void Next_Kmer_Entry(Kmer_Stream *_S)
{ _Kmer_Stream *S = STREAM(_S);

  S->csuf += S->pbyte;
  S->cidx += 1;
  if (S->csuf >= S->ctop)
    { if (S->cidx >= S->nels)
        { S->csuf = NULL;
          S->cpre = S->ixlen;
          S->part = S->nthr+1;
          return;
        }
      More_Kmer_Stream(S);
    }

  while (S->index[S->cpre] <= S->cidx)
    S->cpre += 1;
}

char *Current_Kmer(Kmer_Stream *_S, char *seq)
{ _Kmer_Stream *S = STREAM(_S);
  int    cpre  = S->cpre;
  int    hbyte = S->hbyte;

  if (seq == NULL)
    { seq = (char *) Malloc(S->kmer+3,"Reallocating k-mer buffer");
      if (seq == NULL)
        exit (1);
      if (S->csuf == NULL)
        return (seq);
    }

  { int    j;
    uint8 *a;
    char  *s;

    s = seq;
    switch (S->ibyte)
    { case 3:
        memcpy(s,fmer[cpre>>16],4);
        s += 4;
        memcpy(s,fmer[cpre>>8 & 0xff],4);
        s += 4;
        memcpy(s,fmer[cpre&0xff],4);
        s += 4;
        break;
      case 2:
        memcpy(s,fmer[cpre>>8],4);
        s += 4;
        memcpy(s,fmer[cpre&0xff],4);
        s += 4;
        break;
      case 1:
        memcpy(s,fmer[cpre],4);
        s += 4;
        break;
    }

    a = S->csuf;
    for (j = 0; j < hbyte; j++, s += 4)
      memcpy(s,fmer[a[j]],4);
    seq[S->kmer] = '\0';
  }

  return (seq);
}

inline int Current_Count(Kmer_Stream *S)
{ return (*((uint16 *) (S->csuf+STREAM(S)->hbyte))); }


uint8 *Current_Entry(Kmer_Stream *_S, uint8 *ent)
{ _Kmer_Stream *S = STREAM(_S);
  int    cpre  = S->cpre;
  int    pbyte = S->pbyte;

  if (ent == NULL)
    { ent = (uint8 *) Malloc(S->pbyte,"Reallocating k-mer buffer");
      if (ent == NULL)
        exit (1);
      if (S->csuf == NULL)
        return (ent);
    }

  { int    j;
    uint8 *a;
    uint8 *e;

    e = ent;
    switch (S->ibyte)
    { case 3:
        *e++ = (cpre>>16);
        *e++ = (cpre>>8 & 0xff);
        *e++ = (cpre&0xff);
        break;
      case 2:
        *e++ = (cpre>>8);
        *e++ = (cpre&0xff);
        break;
      case 1:
        *e++ = cpre;
        break;
    }

    a = S->csuf;
    for (j = 0; j < pbyte; j++)
      *e++ = a[j];
  }

  return (ent);
}

  //  Asssumes i is in range

inline void GoTo_Kmer_Index(Kmer_Stream *_S, int64 i)
{ _Kmer_Stream *S = STREAM(_S);
  int64 *index = S->index;
  int    p;

  if (S->cidx == i)
    return;

  S->cidx = i;

  p = S->inver[i>>S->shift];
  while (index[p] <= i)
    p += 1;
  S->cpre = p;

  p = 0;
  while (i >= S->neps[p])
    p += 1;
  
  if (p > 0)
    i -= S->neps[p-1];
  p += 1;

  if (S->part != p)
    { if (S->part <= S->nthr)
        close(S->copn);
      sprintf(S->name+S->nlen,"%d",p);
      S->copn = open(S->name,O_RDONLY);
      S->part = p;
    }

  lseek(S->copn,sizeof(int) + sizeof(int64) + i*S->pbyte,SEEK_SET);

  More_Kmer_Stream(S);
}

int GoTo_Kmer_String(Kmer_Stream *S, char *seq)
{ uint8 entry[S->kbyte];

  if (is_minimal(seq,S->kmer))
    compress_norm(seq,S->kmer,entry);
  else
    compress_comp(seq,S->kmer,entry);

  return (GoTo_Kmer_Entry(S,entry));
}

int GoTo_Kmer_Entry(Kmer_Stream *_S, uint8 *entry)
{ _Kmer_Stream *S = STREAM(_S);

  int64 *index = S->index;
  int    ibyte = S->ibyte;
  int    hbyte = S->hbyte;
  int    pbyte = S->pbyte;
  int64  proff = sizeof(int) + sizeof(int64);

  uint8  kbuf[hbyte];
  int    p, f;
  int64  l, r, m, lo, hi;

  m = *entry++;
  for (l = 1; l < ibyte; l++)
    m = (m << 8) | *entry++;

  if (m == 0)
    l = 0;
  else
    l = index[m-1];
  if (l >= S->nels)
    { if (S->part <= S->nthr)
        close(S->copn);
      S->csuf = NULL;
      S->cidx = S->nels;
      S->cpre = S->ixlen;
      S->part = S->nthr+1;
      return (0);
    }
  r = index[m];
  if (r <= l)
    { GoTo_Kmer_Index(_S,l);
      return (0);
    }
  S->cpre = m;

  if (S->part <= S->nthr)
    close(S->copn);

  hi = r;
  lo = 0;
  for (p = 1; p <= S->nthr; p++)
    { if (l < S->neps[p-1])
        break;
      lo = S->neps[p-1];
    }

  l -= lo;
  r -= lo;

  sprintf(S->name+S->nlen,"%d",p);
  f = open(S->name,O_RDONLY);
  S->part = p;
  S->copn = f;

  // smallest l s.t. KMER(l) >= entry  (or S->neps[p] if does not exist)

  while (r-l > STREAM_BLOCK)
    { m = ((l+r) >> 1);
      lseek(f,proff+m*pbyte,SEEK_SET);
      read(f,kbuf,hbyte);
      if (mycmp(kbuf,entry,hbyte) < 0)
        l = m+1;
      else
        r = m;
    }

  if (l >= S->nels)
    { close(S->copn);
      S->csuf = NULL;
      S->cidx = S->nels;
      S->cpre = S->ixlen;
      S->part = S->nthr+1;
      return (0);
    }

  lseek(f,proff+l*pbyte,SEEK_SET);

  More_Kmer_Stream(S);
  S->cidx = l + lo;

  while (S->cidx < hi)
    { m = mycmp(S->csuf,entry,hbyte);
      if (m >= 0)
        return (m == 0);
      Next_Kmer_Entry(_S);
    }
  return (0);
}


/*********************************************************************************************\
 *
 *  PROFILE CODE
 *
 *****************************************************************************************/

typedef struct
  { int    kmer;     //  Kmer length
    int    nparts;   //  # of threads/parts for the profiles
    int    nreads;   //  total # of reads in data set
    int64 *nbase;    //  nbase[i] for i in [0,nparts) = id of last read in part i + 1
    int64 *index;    //  index[i] for i in [0,nreads) = offset in relevant part of
                     //    compressed profile for read i
                  //  hidden parts
    int    clone;    //  set if a clone
    int    cfile;    //  current open part file (-1 if none)
    int    cpart;    //  index of current open part (-1 if none)
    int    nlen;     //  length of part prefix
    char  *name;     //  part file name prefix
    uint8 *count;    //  decompression buffer
  } _Profile_Index;

#define PROF_BUF0 4096
#define PROF_BUF1 4095

#define PROFILE(P) ((_Profile_Index *) P)


/****************************************************************************************
 *
 *  Open a profile as a Profile_Index.  Index to compressed profiles is in memory,
 *    but compressed profiles are left on disk and reaad only when requested.
 *
 *****************************************************************************************/

Profile_Index *Open_Profiles(char *name)
{ _Profile_Index *P;
  int             kmer, nparts;
  int64           nreads, *nbase, *index;
  uint8          *count;

  int    f, x;
  char  *dir, *root, *full;
  int    smer, nthreads;
  int64  n;

  //  Open stub file and get # of parts

  dir    = PathTo(name);
  root   = Root(name,".prof");
  full   = Malloc(strlen(dir)+strlen(root)+20,"Allocating hidden file names\n");
  sprintf(full,"%s/%s.prof",dir,root);
  f = open(full,O_RDONLY);
  sprintf(full,"%s/.%s.",dir,root);
  x = strlen(full);
  free(root);
  free(dir);
  if (f < 0)
    return (NULL);
  read(f,&smer,sizeof(int));
  read(f,&nthreads,sizeof(int));
  close(f);

  //  Find all parts and accumulate total size

  nreads = 0;
  for (nparts = 0; nparts < nthreads; nparts++)
    { sprintf(full+x,"pidx.%d",nparts+1);
      f = open(full,O_RDONLY);
      if (f < 0)
        { fprintf(stderr,"Profile part %s is misssing ?\n",full);
          exit (1);
        }
      read(f,&kmer,sizeof(int));
      read(f,&n,sizeof(int64));
      read(f,&n,sizeof(int64));
      nreads += n;
      if (kmer != smer)
        { fprintf(stderr,"Profile part %s does not have k-mer length matching stub ?\n",full);
          exit (1);
        }
      close(f);
    }

  //  Allocate in-memory table

  P     = Malloc(sizeof(_Profile_Index),"Allocating profile record");
  index = Malloc((nreads+1)*sizeof(int64),"Allocating profile index");
  nbase = Malloc(nparts*sizeof(int64),"Allocating profile index");
  count = Malloc(PROF_BUF0,"Allocating profile index");
  if (P == NULL || index == NULL || nbase == NULL || count == NULL)
    exit (1);

  nreads = 0;
  index[0] = 0;
  for (nparts = 0; nparts < nthreads; nparts++)
    { sprintf(full+x,"pidx.%d",nparts+1);
      f = open(full,O_RDONLY);
      read(f,&kmer,sizeof(int));
      read(f,&n,sizeof(int64));
      read(f,&n,sizeof(int64));
      read(f,index+(nreads+1),n*sizeof(int64));
      nreads += n;
      nbase[nparts] = nreads;
      close(f);

      sprintf(full+x,"prof.%d",nparts+1);
      f = open(full,O_RDONLY);
      if (f < 0)
        { fprintf(stderr,"Profile part %s is misssing ?\n",full);
          exit (1);
        }
      close(f);
    }

  P->kmer   = kmer;
  P->nparts = nparts;
  P->nreads = nreads;
  P->index  = index;
  P->nbase  = nbase;

  P->clone  = 0;
  P->name   = full;
  P->nlen   = x;
  P->cpart  = -1;
  P->cfile  = -1;
  P->count  = count;

  return ((Profile_Index *) P);
}

Profile_Index *Clone_Profiles(Profile_Index *P)
{ _Profile_Index *Q;
  uint8          *count;
  char           *name;

  //  Allocate in-memory table

  Q     = Malloc(sizeof(_Profile_Index),"Allocating profile record");
  count = Malloc(PROF_BUF0,"Allocating profile index");
  name  = Malloc(PROFILE(P)->nlen + 20,"Allocating profile index");
  if (Q == NULL || count == NULL || name == NULL)
    exit (1);

  *Q = *PROFILE(P);
  strncpy(name,PROFILE(P)->name,Q->nlen);

  Q->clone = 1;
  Q->cpart = -1;
  Q->cfile = -1;
  Q->count = count;
  Q->name  = name;

  return ((Profile_Index *) Q);
}

/****************************************************************************************
 *
 *  Free a Profile_Index and fetch a profile
 *
 *****************************************************************************************/

#undef SHOW_RUN

void Free_Profiles(Profile_Index *_P)
{ _Profile_Index *P = PROFILE(_P);

  if (!P->clone)
    { free(P->index);
      free(P->nbase);
    }
  if (P->cfile >= 0)
    close(P->cfile);
  free(P->name);
  free(P->count);
  free(P);
}

  //  Places uncompressed profile for read id (0-based) in profile of length plen.
  //    Returns the length of the uncompressed profile.  If the plen is less than
  //    this then only the first plen counts are uncompressed into profile

int Fetch_Profile(Profile_Index *_P, int64 id, int plen, uint16 *profile)
{ _Profile_Index *P = PROFILE(_P);

  uint8 *count, *cend;
  int    f;
  int    w, len;
  uint8 *p, *q;
  uint16 x, d, i;
  int    n;

  for (w = 0; w < P->nparts; w++)
    if (id < P->nbase[w])
      break;
  if (w >= P->nparts)
    { fprintf(stderr,"Id %lld is out of range [1,%lld]\n",id,P->nbase[P->nparts-1]);
      exit (1);
    }

  if (w != P->cpart)
    { if (P->cfile >= 0)
        close(P->cfile);
      sprintf(P->name+P->nlen,"prof.%d",w+1);
      f = open(P->name,O_RDONLY);
      if (f < 0)
        { fprintf(stderr,"Profile part %s is misssing ?\n",P->name);
          exit (1);
        }
      P->cfile = f;
      P->cpart = w;
    }
  f = P->cfile;

  if (id == 0 || (w > 0 && id == P->nbase[w-1]))
    { lseek(f,0,SEEK_SET);
      len = P->index[id+1];
    }
  else
    { int64 off = P->index[id];
      lseek(f,off,SEEK_SET);
      len = P->index[id+1] - off;
    }

  if (len == 0)
    return (len);

  count = P->count;
  cend  = count + PROF_BUF1;

  read(f,count,PROF_BUF0);

  p = count;
  q = count + len;

  x = *p++;
  if ((x & 0x80) != 0)
    d = ((x & 0x7f) << 8) | *p++;
  else
    d = x;
  n = 1;

  if (plen > 0)
    { profile[0] = d;
#ifdef SHOW_RUN
  printf(" %d\n",d);
#endif

      while (p < q)
        { if (p >= cend)
            { if (p == cend)
                { *count = *p; 
                  read(f,count+1,PROF_BUF1);
                  q -= PROF_BUF1;
                }
              else
                { read(f,count,PROF_BUF0);
                  q -= PROF_BUF0;
                }
              p = count;
            }
          x = *p++;
          if ((x & 0xc0) == 0)
            { if (n+x > plen)
                { n += x;
                  break;
                }
              for (i = 0; i < x; i++)
                profile[n++] = d;
#ifdef SHOW_RUN
              printf(" [%hu]\n",x);
#endif
            }
          else
            { if ((x & 0x80) != 0)
                { if ((x & 0x40) != 0)
                    x <<= 8;
                  else
                    x = (x << 8) & 0x7fff;
                  x |= *p++;
                  d = (d+x) & 0x7fff;
#ifdef SHOW_RUN
                  printf(" %hd+(%d)\n",x,d);
#endif
                }
              else
                { if ((x & 0x20) != 0)
                    d += (x & 0x1fu) | 0xffe0u;
                  else
                    d += (x & 0x1fu);
#ifdef SHOW_RUN
                  if ((x & 0x20) != 0)
                    printf(" -%d(%d)\n",32-(x&0x1fu),d);
                  else
                    printf(" +%d(%d)\n",x&0x1fu,d);
#endif
                }
              if (n >= plen)
                { n += 1;
                  break;
                }
              profile[n++] = d;
            }
        }
    }

  while (p < q)
    { if (p >= cend)
        { if (p == cend)
            { *count = *p; 
              read(f,count+1,PROF_BUF1);
              q -= PROF_BUF1;
            }
          else
            { read(f,count,PROF_BUF0);
              q -= PROF_BUF0;
            }
          p = count;
        } 
      x = *p++;
      if ((x & 0xc0) == 0)
        n += x;
      else
        { if ((x & 0x80) != 0)
            p += 1;
          n += 1;
        }
    }

  return (n);
}
