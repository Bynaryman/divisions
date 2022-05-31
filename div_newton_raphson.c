#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

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

float div_newton_raphson(float x, float y) {
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

	// are normals
	nx = absX >= 0x00800000;
	ny = absY >= 0x00800000;
	Dm1 = (Ex - nx) - (Ey - ny) - (MX - MY) + (125 + c);

	// getting the correctly rounded result
	// many algorithms exits for this task
	// 3 families:
	//  - digit recurrence algorithms (SRT, etc..)
	//  - functional iteration algorithms (Newton-Raphson iteration, Goldschmidt, etc..)
	//  - polynomial approximation algorithms

	// Here we try Newton Raphson.

	// uint32_t N, D;
	//float N = ldexpf(x,-(Ex-127+1));
	// Here I use divsion but it is witch a power of 2, easy to perform with integer bit manipulation
	float N = fabs(x/(1<<(Ex-127)+1));
	float D = fabs(y/(1<<(Ey-127)+1));
	//float D = ldexpf(y,-(Ey-127+1));
	//D = mpY >> 8;
	float XX = (48.0f/17.0f) - (32.0f/17.0f) * D;
	// 3 steps for single precision
	// ceil(log2((P+1)/log2(17))) where P is precision
	for (int k=0 ; k < 3 ; ++k) {
		XX = XX+XX*(1-D*XX);
	}
	//printf("res %f:\n", N*XX);
	float res = N*XX*(1<<(Ex-Ey));
	float sign = 1.0f;
	if (Sr>0) sign=-1.0f;
	return res*sign;
}


int main(int argc, char *argv[]) {

	float X=0.0f;
	float Y=0.0f;
	X=atof(argv[1]);
	Y=atof(argv[2]);
	float res = div_newton_raphson(X,Y);
	printf("result Newton Raphson: %f\n", res);
	return 0;
}
