#include<typeinfo>

class EnrichmentItem;
using namespace std;

//! A functor
/*! Implementes a functor which will be used in STL algorithms
    This functor check if the typeid() of two objects are identical or not
 */
template<class U,class V>
class IsType{
public:
  bool operator ()(V* t){
	 return (typeid(*t)==typeid(U))? true:false;
  }
};