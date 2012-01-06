//   file FREADER.CPP

#include "freader.h"
#include "string.h"            /*  for function isNext()  */
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

// some lousy compilers (e.g., AT&T for Sun) do not recognize standard stuff:
#ifndef SEEK_SET
#define SEEK_SET 0
#endif


FileReader :: FileReader (char* s)
   // Constructor. Creates a file reader on the file which name is s.
{
   int i ;

   fileName = new char [strlen(s)+1] ;
   strcpy (fileName,s) ;

   file = fopen(fileName,"r") ;
   if (! file) {
      printf ("cannot open file %s \n",fileName) ;
      exit(0) ;}

   atEnd            = FALSE ;
   numberOfKeywords = 15 ;       // modified for XFEM

   keywords         = new long [numberOfKeywords] ;
   for (i=0 ; i<numberOfKeywords ; i++)
      keywords[i] = 0L ;
   strcpy (lastKeyword,"") ;
   strcpy (lastFirstWord,"") ;
   lastPosition   = 0L ;
}


FileReader :: ~FileReader ()
   // Destructor.
{
   fclose(file) ;
   delete fileName ;
   delete keywords ;
}


int  FileReader :: carriageReturn ()
   // Positions the receiver at the beginning of the next line. Returns
   // Failure if end of file is reached, else returns Success.
{
   char nextLine[LINE] ;

   if (fgets(nextLine,LINE,file))
      return SUCCESS ;
   else {
      atEnd = TRUE ;
      return FAILURE ;}
}


int  FileReader :: get (int i, int mode, char* s)
   // Initilizes 's' to the i-th following word on the current line, and
   // returns Success. If reaches the end of line, then :
   //  .if mode is FIND, issues a warning, returns Failure ;
   //  .if mode is LOOK, returns Failure.
{
   char* p = s ;

   while (i-- > 0) {
      if (! this->nextWord(s) || ! strcmp(s,"*"))     // reached end-of-line
	 goto test ;
      else if (! i)                                   // reached i-th value
	 return SUCCESS;}

   test :
   *p = '\0' ;
   if (mode == FIND)
      printf ("Warning : cannot find the %d-th value for this data\n",i) ;
   return FAILURE ;
}


char  FileReader :: next ()
   // Returns the next character in the file. Positions the cursor after
   // that character.
{
   char c ;

   if ((c=fgetc(file)) == EOF) {
      atEnd = TRUE ;
      printf ("Warning : EOF found while reading\n") ;}
   return c ;
}


int  FileReader :: nextWord (char* s)
   // Initializes 's' to the next word on the current line of the file.
   // s is assumed to be large enough to contain such word. Positions the
   // cursor after the first blank following that word. Returns Failure if
   // there is no word between the cursor and the end of line.
{
   char  c ;
   char* p = s ;                          // a pointer to s

   c = this->next() ;                     // initialize the search
   while (isspace(c)) {                   // skip blanks
      if (atEnd || c=='\n') 
         goto failure ;
      c = this->next() ;}
   while (! isspace(c)) {                 // until next blank is reached...
      if (atEnd) 
         goto failure ;
      *p++ = c ;                          // ... store the characters
      c = this->next() ;}
   if (ungetc(c,file)==EOF) {
      printf ("Error in ungetc()\n");
      exit(0) ;}
   *p = '\0' ;
   return SUCCESS ;

   failure :
   *p = '\0' ;
   return FAILURE ;
}


int  FileReader :: nextWordSkipCr (char* s)
   // Initializes 's' to the next word in the file. By contrast with method
   // 'nextWord', this method goes to the next line if end-of-line is
   // reached.
   // s is assumed to be large enough to contain such word. Positions the
   // cursor after the first blank following that word. Returns Failure if
   // there is no word between the cursor and the end of line.
{
   char  c ;
   char* p = s ;                          // a pointer to s

   c = this->next() ;                     // initialize the search
   while (isspace(c)) {                   // skip blanks
      if (atEnd) 
         goto failure ;
      c = this->next() ;}
   while (! isspace(c)) {                 // until next blank is reached...
      if (atEnd) 
         goto failure ;
      *p++ = c ;                          // ... store the characters
      c = this->next() ;}
   if (! ungetc(c,file)) {
      printf ("Error in ungetc\n") ;
      exit(0) ;}
   *p = '\0' ;
   return SUCCESS ;

   failure :
   *p = '\0' ;
   return FAILURE ;
}


void  FileReader :: read (char* k, char* s)
   // Initializes 's' to the word following the keyword k.
{
   if (this -> upToKeyword(k)) {
      if (this -> get(1,FIND,s))
	 return ;
      else
	 goto error ;}
   else
      goto error ;

   error :
   printf ("Error :   cannot read the word following keyword '%s'\n",k) ;
   exit(0) ;
}


void  FileReader :: read (char* k, char* f, char* d, int i, char* s)
   // Reads, in the section of the keyword k, at the line starting with word
   // f, the i-th word following the word d. Initializes s to that string.
{
   if (this -> upToKeywordAndFirstWord(k,f,FIND)) {
      if (this -> upToData(d,FIND)) {
	 if (this -> get(i,FIND,s))
	    return ;
	 else
	    goto error ;}
      else
	 goto error ;}
   else
      goto error ;

   error :
   printf ("Error :   cannot read the value No %d of data '%s' ",i,d) ;
   printf ("in line '%s'\n          of section '%s'\n",f,k) ;
   exit(0) ;
}


int  FileReader :: readIfHas (char* k, char* f, char* d, char* s)
   // If, in the section of the keyword k, the line starting with word f
   // has the word d, initializes s to the word following d, and returns True,
   // else returns False.
{
   if (this -> upToKeywordAndFirstWord(k,f,FIND)) {
      if (this -> upToData(d,LOOK)) {
	 if (this -> nextWord(s))
	    return TRUE ;
	 else
	    goto error ;}
      else
	 return FALSE ;}
   else
      goto error ;

   error :
   printf ("Error :   cannot read if has data '%s' in line '%s' ",d,f) ;
   printf ("of section '%s'\n",k) ;
   exit(0) ;
   return 0 ; // SC

}


int  FileReader :: searchKeyword (char* k)
   // Searches word after word in the file for the string k. A keyword
   // is a string of characters lying at the beginning of a line. Returns
   // Success or Failure.
{
   char s[SIZE] ;

   while (! atEnd) {
      if (! this -> carriageReturn())
	 return FAILURE ;
      if (! this -> nextWordSkipCr(s))
	 return FAILURE ;
      if (! strcmp(s,k))                       // if words match
	 return SUCCESS ;}

   return FAILURE ;
}


int  FileReader :: upToData (char* d, int mode)
   // Positions the receiver past the next appearence of data 'd' in
   // the current line. Returns Success if finds d, else returns Failure
   // (if in "look" mode) or an error messsage (if in "find" mode).
{
   char s[SIZE] ;

   strcpy(s,"") ;

   while (strcmp(s,"*")) {
      if (! strcmp(s,d))                            // found d
	 return SUCCESS ;
      else                                          // search on
	 if (! this->nextWord(s))                   // reached end of line
	    goto test ;}

   test  :
   if (mode == FIND)
      printf("Warning : cannot find data '%s'\n",d);
   return FAILURE ;
}


int  FileReader :: upToFirstWord (char* f, int mode)
   // Positions the receiver past the next appearence of the word f as the
   // first word of a line in the current section. Returns Success if finds
   // f, else returns Failure (with a warning message if in "find" mode).
{
   char s[SIZE] ;

   strcpy(s,"") ;

   while (strcmp(s,"**")) {
      if (! strcmp(s,f)) {                 // found f
	 strcpy(lastFirstWord,f) ;         // keep track of f
	 lastPosition = ftell(file) ;      // keep track of current position
	 return SUCCESS ;}
      else if (! this->carriageReturn() || ! this->nextWordSkipCr(s))
	 goto test ;}

   test :
   if (mode == FIND)
      printf("Warning : cannot find a line starting with '%s'\n",f) ;
   return FAILURE ;
}


int  FileReader :: upToKeyword (char* k)
   // Positions the receiver after the keyword k (a keyword is a string of
   // characters at the beginning of a line). If that position is still un-
   // known, searches the keyword in the file. Returns Success or Failure.
   // Note : if the list of possible keywords is to be modified, do not
   //        forget to update 'numberOfKeywords' in method FileReader(s).
{
   int n ;

   strcpy (lastKeyword,k) ;          // keep track of k
   strcpy (lastFirstWord,"") ;       // no first word yet to remember

   if (! strcmp(k,"Element"))
      n = 0 ;
   else if (! strcmp(k,"Node"))
      n = 1 ;
   else if (! strcmp(k,"Material"))
      n = 2 ;
   else if (! strcmp(k,"Load"))
      n = 3 ;
   else if (! strcmp(k,"LoadTimeFunction"))
      n = 4 ;
   else if (! strcmp(k,"TimeStep"))
      n = 5 ;
   else if (! strcmp(k,"TimeIntegrationScheme"))
      n = 6 ;
   else if (! strcmp(k,"NLSolver"))
      n = 7 ;
   else if (! strcmp(k,"EnrichmentItem"))
      n = 8 ;
   else if (! strcmp(k,"GeometryEntity"))
      n = 9 ;
   else if (! strcmp(k,"EnrichmentFunction"))
      n = 10 ;
	else if (! strcmp(k,"CrackGrowthDirectionLaw"))
      n = 11 ;
	else if (! strcmp(k,"CrackGrowthIncrementLaw"))
      n = 12 ;
	else if (! strcmp(k,"ExportStress"))
      n = 13 ;
	else if (! strcmp(k,"ExportDisplacement"))
      n = 14 ;
   else {
      printf ("%s : unknown keyword for reading data \n",k) ;
      exit(0) ;}

   if (keywords[n] == 0L) {            // the position of k is still unknown
      if (this -> searchKeyword(k))
	 keywords[n] = ftell(file) ;
      else {
	 printf ("warning : cannot find keyword '%s'\n",k) ;
	 return FAILURE ;}}

   else                                // the position of k is already known
      fseek(file,keywords[n],SEEK_SET) ;

   return SUCCESS ;
}


int  FileReader :: upToKeywordAndFirstWord (char* k, char* f, int mode)
   // Positions the receiver after f, where f if the first word of a line
   // in section k. Returns Success or Failure.
   // This method executes fast if k=lastKeyword and f=lastFirstWord.
{
   if (!strcmp(k,lastKeyword) && !strcmp(f,lastFirstWord)) {      // fast !
      fseek(file,lastPosition,SEEK_SET) ;
      return SUCCESS ;}

   else {
      if (this->upToKeyword(k))                                   // slow !
	 return this->upToFirstWord(f,mode) ;
      else
	 return FAILURE ;}
}
