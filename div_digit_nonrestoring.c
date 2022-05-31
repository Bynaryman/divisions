#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

static inline
int32_t max (int32_t a ,int32_t b)
{
	return (a > b) ? a : b;
}

static inline
uint32_t maxu (uint32_t a ,uint32_t b)
{
	return (a > b) ? a : b;
}

static inline
int32_t min (int32_t a ,int32_t b)
{
	return (a < b) ? a : b;
}

static inline
uint32_t minu (uint32_t a ,uint32_t b)
{
	return (a < b) ? a : b;
}

static inline
uint32_t mul (uint32_t a ,uint32_t b)
{
	uint64_t t0 = a;
	uint64_t t1 = b;
	uint64_t t2 = (t0 * t1) >> 32;
	return t2;
}

static inline
int32_t mul64h (int32_t a ,int32_t b)
{
	int64_t t0 = a;
	int64_t t1 = b;
	int64_t t2 = (t0 * t1) >> 32;
	return t2;
}

// number of leading zeros
static inline
uint32_t nlz(uint32_t x)
{
	uint32_t z = 0;
	if (x == 0) return(32);
	if (x <= 0x0000FFFF) {z = z + 16; x = x << 16;}
	if (x <= 0x00FFFFFF) {z = z + 8; x = x << 8;}
	if (x <= 0x0FFFFFFF) {z = z + 4; x = x << 4;}
	if (x <= 0x3FFFFFFF) {z = z + 2; x = x << 2;}
	if (x <= 0x7FFFFFFF) {z = z + 1;}
	return z;
}

uint32_t div_digit_recurrence_nonrestoring(float x, float y) {
	//casting to integers
	int32_t X = *(uint32_t*)&x;
	int32_t Y = *(uint32_t*)&y;

	// Special value handling for binary32
	uint32_t absX, absY, Sr, absXm1, absYm1, Max, Inf;
	absX = X & 0x7FFFFFFF; absY = Y & 0x7FFFFFFF; Sr = (X ^ Y) & 0x80000000;
	absXm1 = absX - 1;
	absYm1 = absY - 1;
	if (maxu(absXm1, absYm1) >= 0x7F7FFFFF)
	{
		Max = maxu(absX, absY); Inf = Sr | 0x7F800000;
		if (Max > 0x7F800000 || absX == absY)
			return Inf | 0x00400000 | Max; // qNaN with payload encoded in
						       // the last 22 bits of X or Y
		if (absX < absY) return Sr;
		return Inf;
	}

	// Computing normalized significands mpx mpy and shift value c
	uint32_t MX, MY, mpX, mpY, c;
	MX = maxu(nlz(absX), 8);
	MY = maxu(nlz(absY), 8);
	mpX = (X << MX) | 0x80000000;
	mpY = (Y << MY) | 0x80000000;
	//MX = max(nlz(absX), 8);
	//MY = max(nlz(absY), 8);
	c = mpX >= mpY;

	// Computing D-1 (D=d+emax, d=epx-epy-1+c)
	uint32_t Ex, Ey, nx, ny;
	uint32_t Dm1;
	Ex = absX >> 23;
	Ey = absY >> 23;

	nx = absX >= 0x00800000;
	ny = absY >= 0x00800000;
	Dm1 = (Ex - nx) - (Ey - ny) - (MX - MY) + (125 + c);

	// getting the correctly rounded result
	// many algorithms exits for this task
	// 3 families:
	//  - digit recurrence algorithms (SRT, etc..)
	//  - functional iteration algorithms (Newton-Raphson iteration, Goldschmidt, etc..)
	//  - polynomial approximation algorithms

	// HERE we do digit recurrence non restoring, highly sequential algorithm
	// p iterations for p mantissa bits (here 24 in float single)


	uint32_t N, M, L, j;
	int32_t Wj;
	N = mpX >> (7 + c);
	M = mpY >> 8;
	L = 1;
	Wj = N - M;
	Wj = (Wj << 1) - M;
	for (j = 1; j < 24; j++)
	{
		// Non restoring iteration j
		if (Wj >= 0) {
			L = (L << 1) | 1;
			Wj = (Wj << 1) - M;
		} else {
			Wj = (Wj << 1) + M;
			L = L << 1;
		}

	}
	// Correction iteration j
	if (Wj >= 0) {
		L = (L << 1) | 1;
	} else {
		Wj = (Wj << 1) + M;
		L = L << 1;
	}

	// final packing of encoded float
	return (Sr | (Dm1 << 23)) + ((L >> 1) + (L & 1));
}


int main(int argc, char *argv[]) {

	float X=0.0f;
	float Y=0.0f;
	X=atof(argv[1]);
	Y=atof(argv[2]);
	uint32_t res=0;
	res = div_digit_recurrence_nonrestoring(X,Y);
	printf("result digit recurrence non restoring: %f\n", *(float*)&res);
	return 0;
}
