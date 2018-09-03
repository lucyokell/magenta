// Ensure the full strain template definition can be seen
#include "strain.h"
#include "strain.cpp"

template class Strain<unsigned char>;
template class Strain<unsigned short>; 
template class Strain<unsigned long>; 
template class Strain<unsigned long long>; 
  
  // instantiate other templates here