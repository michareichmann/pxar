#ifndef SHIST600_H
#define SHIST600_H

#include "pxardllexport.h"

class DLLEXPORT shist600{
public:
  shist600();
  ~shist600();

  void  fill(int x, float w = 1.); 
  void  clear();
  float get(int i); 
  float get(float i); 
  float getSumOfWeights(); 
  
private:
  static const int NBINS = 556;

  float fX[NBINS + 2];
  // fX[0]   = underflow_low
  // fX[1]   =   0 ..   1 low
  // fX[256] = 255 .. 256 low
  // fX[257] = overflow_low
  // fX[258 .. 299] = unused
  // fX[300] = underflow_high
  // fX[301] = 0 .. 1 high
  // fX[556] = 255 .. 256 high
  // fX[557] = overflow_high
  // fX[558 .. 599] = unused
};

#endif
