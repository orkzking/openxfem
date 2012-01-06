//   *************************
//   *** CLASS FILE READER ***
//   *************************
 

#include <stdio.h>
#define SIZE 32             /* the maximum length of any word in the file */
#define LINE 128            /* the maximum length of any line in the file */


class FileReader
/*!
   This class provides access to data for objects of the finite element
   method.
 DESCRIPTION :
   The main attribute of a FileReader is a FILE structure (attribute 'file')
   on the file 'fileName'. Attribute 'atEnd' is True if 'file' is positioned
   at end of file.
   In order to speed up information reading, additional attributes are
   introduced :
   .the file reader keeps track (in 'keywords') of the position in the
    file of all keywords. It can therefore position itself instantaneously
    at any known keyword. The keyword list can be easily modified (see
    method 'upToKeyword') ;
   .the file reader keeps track of the last object (eg, material 10) it
    has read data for. Here 'lastKeyword' is set to "Material", 'last-
    FirstWord' to "10", and 'lastPosition' to the position of the begin-
    ning of the line of material 10 in the file.
    This is useful if material 10 reads its data sequentially.
 TASKS :
   The main methods of a FileReader are the 3 'read' methods. For example,
   the number of the 4-th node of the 10-th element is obtained by the
   message aFileReader->read("Element","10","node",4,answer.
 REMARKS :
   .The objects read their input in the data file by means of a FileReader.
    They write their output in the same data file by means of a simple FILE.
    Both these streams are attributes of the domain.
 */
{
   enum { FALSE , TRUE } ;
   enum { FAILURE , SUCCESS } ;
   enum { LOOK , FIND } ;

   private :
      char*  fileName ;
      FILE*  file ;
      int    atEnd ;
      int    numberOfKeywords ;
      long*  keywords ;
      char   lastKeyword [SIZE] ;
      char   lastFirstWord [SIZE] ;
      long   lastPosition ;

   public :
      FileReader (char*) ;                         // constructor
      ~FileReader () ;                             // destructor

      int          carriageReturn () ;
      int          get (int,int,char*) ;
      char         next () ;
      int          nextWord (char*) ;
      int          nextWordSkipCr (char*) ;
      void         read (char*,char*) ;
      void         read (char*,char*,char*,int,char*) ;
      int          readIfHas (char*,char*,char*,char*) ;
      FileReader*  reset()                 { rewind(file) ; return this ;}
      int          searchKeyword (char*) ;
      int          upToData (char*,int) ;
      int          upToFirstWord (char*,int) ;
      int          upToKeyword (char*) ;
      int          upToKeywordAndFirstWord (char*,char*,int) ;
} ;









