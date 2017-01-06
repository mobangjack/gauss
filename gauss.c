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

#include "gauss.h"

#define SQR(x) (x*x)

Gauss_t* GaussCreate(uint16_t N)
{
	Gauss_t* gauss = (Gauss_t*)malloc(sizeof(Gauss_t));
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

void GaussReset(Gauss_t* gauss)
{
	memset(gauss->buf, 0, gauss->N * sizeof(float));
	gauss->n = 0;
	gauss->i = 0;
	gauss->sum = 0;
	gauss->delta_sum = 0;
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
	gauss->err = FLT_MAX;
}

void GaussProc(Gauss_t* gauss, float x) {
	gauss->last_mean = gauss->mean;
	gauss->last_mse = gauss->mse;
	if (gauss->n < gauss->N) {
		gauss->buf[gauss->n++] = x;
		gauss->delta_sum = x;
		gauss->sum += gauss->delta_sum;
		gauss->mean = gauss->sum / gauss->n;
		gauss->delta_mean = gauss->mean - gauss->last_mean;
		gauss->resn = x - gauss->mean;
		gauss->sse += SQR(gauss->resn);
		gauss->mse = gauss->sse / gauss->n;
		gauss->delta_mse = gauss->mse - gauss->last_mse;
		gauss->err = SQR(gauss->delta_mean) + SQR(gauss->delta_mse);
	} else {
		if (gauss->i == gauss->N) gauss->i = 0;
		gauss->delta_sum = x - gauss->buf[gauss->i];
		gauss->sum += gauss->delta_sum;
		gauss->mean = gauss->sum / gauss->n;
		gauss->delta_mean = gauss->mean - gauss->last_mean;
		gauss->res0 = gauss->buf[gauss->i] - gauss->last_mean;
		gauss->resn = x - gauss->mean;
		gauss->delta_sse = (gauss->delta_sum) * (gauss->resn + gauss->res0);
		gauss->sse += gauss->delta_sse;
		gauss->mse = gauss->sse / gauss->n;
		gauss->delta_mse = gauss->mse - gauss->last_mse;
		gauss->err = SQR(gauss->delta_mean) + SQR(gauss->delta_mse);
		gauss->buf[gauss->i++] = x;
	}
}

void GaussDestroy(Gauss_t* gauss)
{
	if (gauss != NULL) {
		if (gauss->buf != NULL) {
			free(gauss->buf);
			gauss->buf = NULL;
		}
		free(gauss);
		gauss = NULL;
	}
}
