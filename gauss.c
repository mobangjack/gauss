#include "gauss.h"

#define SQR(x) (x*x)

struct Gauss* GaussCreate(uint16_t N)
{
  Gauss* gauss = (Gauss*)malloc(sizeof(Gauss));
  if (gauss == NULL) return NULL;
  gauss->buf = (float*)malloc(N * sizeof(float));
  if (gauss->buf == NULL) {
    free(gauss);
    gauss = NULL;
    return NULL;
  }
  gauss->N = N;
  GaussReset(gauss);
  return gauss;
}

void GaussReset(struct Gauss* gauss)
{
  memset(gauss->buf, 0, gauss->N * sizeof(float));
  gauss->n = 0;
  gauss->i = 0;
  gauss->sum = 0;
  gauss->diff = 0;
  gauss->res0 = 0;
  gauss->resn = 0;
  gauss->mean = 0;
  gauss->last_mean = 0;
  gauss->delta_mean = 0;
  gauss->sse = 0;
  gauss->delta_sse = 0;
  gauss->mse = 0;
  gauss->last_mse = 0;
  gauss->delta_mse = 0;
}

float GaussFilter(struct Gauss* gauss, float x) {
  gauss->last_mean = gauss->mean;
  gauss->last_mse = gauss->mse;
  if (gauss->n < gauss->N) {
    gauss->buf[gauss->n++] = x;
    gauss->sum += x;
    gauss->mean = gauss->sum / gauss->n;
    gauss->delta_mean = gauss->mean - gauss->last_mean;
    gauss->resn = x - gauss->mean;
    gauss->sse += SQR(gauss->resn);
    gauss->mse = gauss->sse / gauss->n;
    gauss->delta_mse = gauss->mse - gauss->last_mse;
  } else {
    if (gauss->i == gauss->N) gauss->i = 0;
    gauss->diff = x - gauss->buf[gauss->i];
    gauss->sum += gauss->diff;
    gauss->mean = gauss->sum / gauss->n;
    gauss->delta_mean = gauss->mean - gauss->last_mean;
    gauss->res0 = gauss->buf[gauss->i] - gauss->last_mean;
    gauss->resn = x - gauss->mean;
    gauss->delta_sse = (gauss->diff) * (gauss->resn + gauss->res0);
    gauss->sse += gauss->delta_sse;
    gauss->mse = gauss->sse / gauss->n;
    gauss->delta_mse = gauss->mse - gauss->last_mse;
    gauss->buf[gauss->i++] = x;
  }
  return gauss->mean;
}

void GaussDestroy(struct Gauss* gauss)
{
  if (gauss != NULL) {
    free(gauss->buf);
    gauss->buf = NULL;
    free(gauss);
    gauss = NULL;
  }
}
