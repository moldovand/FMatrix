#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "FileScan.h"
#include "defVM.h"

int FeedComment(FILE *fp)
{
  int c;

  while((c = fgetc(fp)) != '\n') ;

  return true;
}

int ScanForScan(FILE *fp)
{
  int c;

  c = fgetc(fp);
  while((c != EOF) 
	&& (c != '#') && !isdigit(c) && (c != '-') && (c != 'p')) {
    if (c == '%') FeedComment(fp);
    c = fgetc(fp);
  }

  ungetc(c, fp);
  
  if (c == EOF)      return EOF;
  else if (c == '#') return IDETA;
  else if (c == 'p') return P_NUM;
  else               return DIGIT;
}
