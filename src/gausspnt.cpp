//   file GAUSSPNT.CPP
 
#include "gausspnt.h"
#include "element.h"
#include "domain.h"
#include "string.h"
#include "material.h"
#include <math.h>



GaussPoint :: GaussPoint (Element* e, int n, Mu::Point* a, double w, int aSize)
   // Constructor. Creates a Gauss point belonging to element e, with number
   // n, with coordinates a, with weight w.
{
   element      = e ;
   number       = n ;
   coordinates  = a ;
   weight       = w ;
   strainVector = NULL ;
   stressVector = NULL ;
   backStressVector = NULL ;
   previousStressVector = new FloatArray(aSize) ;
   previousBackStressVector = new FloatArray(aSize) ;
   plasticCode = 0;
   stressLevel = 0;
   deltaGamma = 0;
   yieldRadius = 0;
   previousYieldRadius = 0;
   sizeOfVector = aSize;
}


GaussPoint :: ~GaussPoint ()
   // Destructor.
{
   delete coordinates ;
   delete strainVector ;
   delete stressVector ;
   delete previousStressVector ;
   delete backStressVector ;
   delete previousBackStressVector ;
}


void  GaussPoint :: printOutput (FILE* strFile)
   // Prints the strains and stresses on the data file.
{
   int  i,n ;

   fprintf (strFile,"  GP %d :  strains ",number) ;
   n = strainVector->giveSize() ;
   for (i=1 ; i<=n ; i++)
	   fprintf (strFile," % .4e",strainVector->at(i)) ;

   fprintf (strFile,"\n          stresses",number) ;
   n = stressVector -> giveSize() ;
   for (i=1 ; i<=n ; i++)
	   fprintf (strFile," % .4e",stressVector->at(i)) ;

   fprintf (strFile,"\n          plastic code : %d", plasticCode);
   fprintf (strFile,"\n          stress level : %.4e", stressLevel);

   fprintf (strFile,"\n") ;
}

void  GaussPoint :: printBinaryResults (FILE* s01File)
   // Prints the strains and stresses on the data file.
{
	float a[11];
	if (s01File != NULL) {
		//to put it at the end of the file
		int pos = fseek(s01File, 0L, SEEK_END);
		//writing
		this->getBinaryRecord(a);
		fwrite(a, sizeof(a), 1, s01File);
	}
}

void  GaussPoint :: getBinaryRecord (float* a)
{
  // SC - 12.98
  // check to integrate trusses results for postpro!!!
  int aNumberOfGP = element->giveNumberOfGaussPoints();

  a[0] = (float)aNumberOfGP;
  for (int i = 1 ; i <= 4; i++) 
  {
	 a[i]   = (float)this->stressVector->at(i);
	 a[i+4] = (float)this->strainVector->at(i);
  }

  a[9] = (float)this->plasticCode;
  a[10] = (float)this->stressLevel;
}


void  GaussPoint :: updateYourself ()
   // Performs end-of-step updates.
{  //on delete previousStressVector (il ne pointe sur plus rien du tout)
   //on le fait pointer sur stressVector (sn-1 = sn)
   //on met s = NULL MAIS ON NE DELETE PAS s !!! (sinon, par la meme
   //occasion on deleterait aussi sn-1!!!)
   delete previousStressVector;
   previousStressVector = stressVector;
   delete previousBackStressVector;
   previousBackStressVector = backStressVector;

   delete strainVector ;
   strainVector = NULL ;

   stressVector = NULL ;
   backStressVector = NULL ;

   previousYieldRadius = yieldRadius;

   plasticCode = 0;
   stressLevel = 0;
   deltaGamma = 0;
}

void GaussPoint::updateYourselfForXFEM()
// *************************************
// 2005-09-13
{
   delete previousStressVector ;
   previousStressVector = NULL ;

   delete previousBackStressVector;
   previousBackStressVector = NULL ;

   delete strainVector ;
   strainVector = NULL ;

   stressVector = NULL ;
   backStressVector = NULL ;

   plasticCode = 0;
   stressLevel = 0;
   deltaGamma = 0;
}

void  GaussPoint :: computeStressLevel ()  
{        
	Material* mat;
	mat    = this->element->giveMaterial();
	stressLevel = mat->computeStressLevelFor(this);
}

double GaussPoint :: computeVonMisesStress()
//******************************************
// sigma_e = \sqrt{3}(II_s)^{1/2}
{
  double II = stressVector->computeInvariantJ2();
  double equiStress = sqrt(3.0 * II) ;

  return equiStress ;
}

void GaussPoint::letStrainVectorBe (FloatArray* v)
{ 
	if (strainVector)
		delete strainVector ;
	strainVector = v;
}


void GaussPoint::letStressVectorBe (FloatArray* v) 
{ 
	if (stressVector)
		delete stressVector ;
	stressVector = v;
}


void GaussPoint::letBackStressVectorBe (FloatArray* v) 
{ 
	if (backStressVector)
		delete backStressVector ;
	backStressVector = v;
}
