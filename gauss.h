/**
 * Copyright (c) 2011-2016, Jack Mo (mobangjack@foxmail.com).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __GAUSS_H__
#define __GAUSS_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <float.h>

typedef struct
{
  /* Critical Section */
  uint16_t N; // buffer length
  /* Critical Section */
  
  float* buf; // buffer
  uint16_t n; // initialization counter
  uint16_t i; // current enqueuing buffer index
  float sum;  // sum of buffer
  float delta_sum;  // delta of sum
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
  float err; // gauss err
}Gauss_t;

Gauss_t* GaussCreate(uint16_t N);
void GaussReset(Gauss_t* gauss);
void GaussProc(Gauss_t* gauss, float x);
void GaussDestroy(Gauss_t* gauss);

#ifdef __cplusplus
}
#endif

#endif

