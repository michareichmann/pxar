#include "shist600.hh"

// ----------------------------------------------------------------------
shist600::shist600() {
  clear();
}

// ----------------------------------------------------------------------
shist600::~shist600() {
}

// ----------------------------------------------------------------------
void shist600::clear() {
  for (int i = 0; i < NBINS+2; ++i) fX[i] = 0.;
}

// ----------------------------------------------------------------------
void shist600::fill(int x, float w) {
  if (x < 0) {
    fX[0] += w; 
  } else if (256 < i && i < 299) {
      fX[257] += w;
  } else if (299 <= i && i < 300) {
      fX[300] += w;
  } else if (x > 556) {
    fX[NBINS+1] += w; 
  } else {
    fX[x+1] += w; 
  }
}

// ----------------------------------------------------------------------
float shist600::get(int i) {
  if (i < 0) {
    return fX[0];
  } else if(i > 256 && i < 299){
      return fX[257];
  } else if(299 <= i && i < 300){
      return fX[300];
  } else if (i > 556) {
    return fX[NBINS+1];
  } else {
    return fX[i+1];
  }
}

// ----------------------------------------------------------------------
float shist600::get(float i) {
  int ii(i); 
  if (i < 0) ii = i - 1.;
  return get(ii);
}

// ----------------------------------------------------------------------
float shist600::getSumOfWeights() {
  float sum(0.); 
  for (int i = 0; i < NBINS+2; ++i) sum += fX[i];
  return sum;
}
