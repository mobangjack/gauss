#ifndef __GAUSS_H__
#define __GAUSS_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct Gauss
{
  /* Critical Section */
  uint16_t N; // buffer length
  /* Critical Section */
  
  float* buf; // buffer
  uint16_t n; // initialization counter
  uint16_t i; // current enqueuing buffer index
  float sum;  // sum of buffer
  float diff; // diff of sum
  float res0; // residual of the number 0 parameter
  float resn; // residual of the number (n - 1) parameter
  float mean; // mean of buffer
  float last_mean; // last mean
  float delta_mean; // delta of mean
  float sse; // sum of square error
  float delta_sse; // delta of sse
  float mse; // mean square error
  float last_mse; // last mse
  float delta_mse; // delta mse
}Gauss;

struct Gauss* GaussCreate(uint16_t N);
void GaussReset(struct Gauss* gauss);
float GaussFilter(struct Gauss* gauss, float x);
void GaussDestroy(struct Gauss* gauss);

#ifdef __cplusplus
}
#endif

#endif

