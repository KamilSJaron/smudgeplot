#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <dirent.h>
#include <zlib.h>

#include "gene_core.h"

/*******************************************************************************************
 *
 *  GENERAL UTILITIES
 *
 ********************************************************************************************/

char *Prog_Name;

void *Malloc(int64 size, char *mesg)
{ void *p;

  if ((p = malloc(size)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

void *Realloc(void *p, int64 size, char *mesg)
{ if (size <= 0)
    size = 1;
  if ((p = realloc(p,size)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

char *Strdup(char *name, char *mesg)
{ char *s;

  if (name == NULL)
    return (NULL);
  if ((s = strdup(name)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (s);
}

char *Strndup(char *name, int len, char *mesg)
{ char *s;

  if (name == NULL)
    return (NULL);
  if ((s = strndup(name,len)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (s);
}

char *PathTo(char *name)
{ char *path, *find;

  if (name == NULL)
    return (NULL);
  if ((find = rindex(name,'/')) != NULL)
    path = Strndup(name,find-name,"Extracting path from");
  else
    path = Strdup(".","Allocating default path");
  return (path);
}

char *Root(char *name, char *suffix)
{ char *path, *find, *dot;
  int   epos;

  if (name == NULL)
    return (NULL);
  find = rindex(name,'/');
  if (find == NULL)
    find = name;
  else
    find += 1;
  if (suffix == NULL)
    { dot = strrchr(find,'.');
      path = Strndup(find,dot-find,"Extracting root from");
    }
  else
    { epos  = strlen(find);
      epos -= strlen(suffix);
      if (epos > 0 && strcasecmp(find+epos,suffix) == 0)
        path = Strndup(find,epos,"Extracting root from");
      else
        path = Strdup(find,"Allocating root");
    }
  return (path);
}

char *Catenate(char *path, char *sep, char *root, char *suffix)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  if (path == NULL || root == NULL || sep == NULL || suffix == NULL)
    { free(cat);
      max = -1;
      return (NULL);
    }
  len =  strlen(path);
  len += strlen(sep);
  len += strlen(root);
  len += strlen(suffix);
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name for %s)\n",Prog_Name,root);
          return (NULL);
        }
    }
  sprintf(cat,"%s%s%s%s",path,sep,root,suffix);
  return (cat);
}

char *Numbered_Suffix(char *left, int num, char *right)
{ static char *suffix = NULL;
  static int   max = -1;
  int len;

  if (left == NULL || right == NULL)
    { free(suffix);
      max = -1;
      return (NULL);
    }
  len =  strlen(left);
  len += strlen(right) + 40;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((suffix = (char *) realloc(suffix,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making number suffix for %d)\n",Prog_Name,num);
          return (NULL);
        }
    }
  sprintf(suffix,"%s%d%s",left,num,right);
  return (suffix);
}


#define  COMMA  ','

//  Print big integers with commas/periods for better readability

void Print_Number(int64 num, int width, FILE *out)
{ if (width == 0)
    { if (num < 1000ll)
        fprintf(out,"%lld",num);
      else if (num < 1000000ll)
        fprintf(out,"%lld%c%03lld",num/1000ll,COMMA,num%1000ll);
      else if (num < 1000000000ll)
        fprintf(out,"%lld%c%03lld%c%03lld",num/1000000ll,
                                           COMMA,(num%1000000ll)/1000ll,COMMA,num%1000ll);
      else
        fprintf(out,"%lld%c%03lld%c%03lld%c%03lld",num/1000000000ll,
                                                   COMMA,(num%1000000000ll)/1000000ll,
                                                   COMMA,(num%1000000ll)/1000ll,COMMA,num%1000ll);
    }
  else
    { if (num < 1000ll)
        fprintf(out,"%*lld",width,num);
      else if (num < 1000000ll)
        { if (width <= 4)
            fprintf(out,"%lld%c%03lld",num/1000ll,COMMA,num%1000ll);
          else
            fprintf(out,"%*lld%c%03lld",width-4,num/1000ll,COMMA,num%1000ll);
        }
      else if (num < 1000000000ll)
        { if (width <= 8)
            fprintf(out,"%lld%c%03lld%c%03lld",num/1000000ll,COMMA,(num%1000000ll)/1000ll,
                                               COMMA,num%1000ll);
          else
            fprintf(out,"%*lld%c%03lld%c%03lld",width-8,num/1000000ll,COMMA,(num%1000000ll)/1000ll,
                                                COMMA,num%1000ll);
        }
      else
        { if (width <= 12)
            fprintf(out,"%lld%c%03lld%c%03lld%c%03lld",num/1000000000ll,COMMA,
                                                       (num%1000000000ll)/1000000ll,COMMA,
                                                       (num%1000000ll)/1000ll,COMMA,num%1000ll);
          else
            fprintf(out,"%*lld%c%03lld%c%03lld%c%03lld",width-12,num/1000000000ll,COMMA,
                                                        (num%1000000000ll)/1000000ll,COMMA,
                                                        (num%1000000ll)/1000ll,COMMA,num%1000ll);
        }
    }
}

//  Return the number of symbols to print num, base 10 (without commas as above)

int  Number_Digits(int64 num)
{ int digit;

  if (num == 0)
    return (1);
  if (num < 0)
    { num = -num;
      digit = 1;
    }
  else
    digit = 0;
  while (num >= 1)
    { num /= 10;
      digit += 1;
    }
  return (digit);
}


/*******************************************************************************************
 *
 *  READ AND ARROW COMPRESSION/DECOMPRESSION UTILITIES
 *
 ********************************************************************************************/

//  Compress read into 2-bits per base (from [0-3] per byte representation

void Compress_Read(int len, char *s)
{ int   i; 
  char  c, d;
  char *s0, *s1, *s2, *s3;
  
  s0 = s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;
  
  c = s1[len];
  d = s2[len];
  s0[len] = s1[len] = s2[len] = 0;
  
  for (i = 0; i < len; i += 4)
    *s++ = (char ) ((s0[i] << 6) | (s1[i] << 4) | (s2[i] << 2) | s3[i]);
  
  s1[len] = c;
  s2[len] = d;
}

//  Uncompress read form 2-bits per base into [0-3] per byte representation

void Uncompress_Read(int len, char *s)
{ int   i, tlen, byte;
  char *s0, *s1, *s2, *s3;
  char *t;

  s0 = s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  tlen = (len-1)/4;

  t = s+tlen;
  for (i = tlen*4; i >= 0; i -= 4)
    { byte = *t--;
      s0[i] = (char) ((byte >> 6) & 0x3);
      s1[i] = (char) ((byte >> 4) & 0x3);
      s2[i] = (char) ((byte >> 2) & 0x3);
      s3[i] = (char) (byte & 0x3);
    }
  s[len] = 4;
}

//  Convert read in [0-3] representation to ascii representation (end with '\n')

void Lower_Read(char *s)
{ static char letter[4] = { 'a', 'c', 'g', 't' };

  for ( ; *s != 4; s++)
    *s = letter[(int) *s];
  *s = '\0';
}

void Upper_Read(char *s)
{ static char letter[4] = { 'A', 'C', 'G', 'T' };

  for ( ; *s != 4; s++)
    *s = letter[(int) *s];
  *s = '\0';
}

void Letter_Arrow(char *s)
{ static char letter[4] = { '1', '2', '3', '4' };

  for ( ; *s != 4; s++)
    *s = letter[(int) *s];
  *s = '\0';
}

//  Convert read in ascii representation to [0-3] representation (end with 4)

void Number_Read(char *s)
{ static char number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };

  for ( ; *s != '\0'; s++)
    *s = number[(int) *s];
  *s = 4;
}

void Number_Arrow(char *s)
{ static char arrow[128] =
    { 3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 0, 1, 2, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 2,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
    };

  for ( ; *s != '\0'; s++)
    *s = arrow[(int) *s];
  *s = 4;
}

void Change_Read(char *s)
{ static char change[128] =
    {   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0, 'a',   0, 'c',   0,   0,   0, 'g',
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0, 't',   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
        0, 'A',   0, 'C',   0,   0,   0, 'G',
        0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0, 'T',   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,
    };

  for ( ; *s != '\0'; s++)
    *s = change[(int) *s];
}
