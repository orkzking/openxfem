//   file MATERIAL.CPP
 
#include "material.h"
#include "elasticmaterial.h"
#include "vonmisesmaterial.h"
#include "vonmisesmaterial_h.h"
#include "dictionr.h"
#include "domain.h"
#include "gausspnt.h"
#include "element.h"
#include <stdlib.h>

Material::Material(int n,Domain* d) : FEMComponent(n,d)
{	propertyDictionary = new Dictionary() ;
}

Material::~Material ()                 
{	delete propertyDictionary ;
}

char* Material :: giveClassName (char* s)
{	return strcpy(s,"Material") ;
}

void  Material :: giveKeyword (char* keyword)
{
	strcpy(keyword,"Material") ;
}

double  Material :: give (char aProperty)
   // Returns the value of the property aProperty (e.g. the Young's modulus
   // 'E') of the receiver.
{
   double  value ;
 
   if (propertyDictionary -> includes(aProperty))
      value = propertyDictionary -> at(aProperty) ;
   else {                                         // read and store the value
      value = this -> read(aProperty) ;
      propertyDictionary -> add(aProperty,value) ;}
   return value ;
}

void  Material :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("Material with properties : \n") ;
   propertyDictionary -> printYourself() ;
}

Material*  Material :: ofType (char* aClass)
   // Returns a new material, which has the same number than the receiver,
   // but belongs to aClass (ElasticMaterial, or VonMisesMaterial)
{
   Material* newMaterial ;

   if (! strcmp(aClass,"ElasticMaterial"))
      newMaterial = new ElasticMaterial(number,domain) ;
   else if (! strcmp(aClass,"VonMisesMaterial"))
      newMaterial = new VonMisesMaterial(number,domain) ;
   else if (! strcmp(aClass,"VonMisesMaterial_H"))
      newMaterial = new VonMisesMaterial_H(number,domain) ;
   else {
      printf ("%s : unknown material type \n",aClass) ;
      exit(0) ;}

   return newMaterial ;
}

Material*  Material :: typed ()
   // Returns a new material, which has the same number than the receiver,
   // but belongs to aClass (ElasticMaterial, or VonMisesMaterial)

{
   Material* newMaterial ;
   char     type[32] ;

   this -> readString("class",type) ;
   newMaterial = this -> ofType(type) ;

   return newMaterial ;
}






