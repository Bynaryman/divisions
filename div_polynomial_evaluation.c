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

uint32_t div_digit_recurrence_restoring(float x, float y) {
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

	// Here we try polyninial approx.


	// Computing the one-sided approximation v
	// Although a functional iteration approach (the Newton-Raphson iteration or
	// one of its variants; see Section 4.7 and Section 7.6.2) could be used to compute
	// v, more ILP can be exposed by evaluating a suitable polynomial that approxi-
	// mates the exact quotient l.

	uint32_t S, T;
	S = mpX >> c ;
	T = (X<<1) << MX;
	printf("s=%x, t=%x\n", S, T);

	uint32_t r0 = mul(T, 0xffffe7d7); // 0.32
	uint32_t r1 = 0xffffffe8 - r0; // 0.32
	uint32_t r2 = mul(S, r1); // 2.30
	uint32_t r3 = 0x00000020 + r2; // 2.30
	uint32_t r4 = mul(T, T); // 0.32
	uint32_t r5 = mul(S, r4); // 2.30
	uint32_t r6 = mul(T, 0xffbad86f); // 0.32
	uint32_t r7 = 0xfffbece7 - r6; // 0.32
	uint32_t r8 = mul(r5, r7); // 2.30
	uint32_t r9 = r3 + r8; // 2.30
	uint32_t r10 = mul(r4, r5); // 2.30
	uint32_t r11 = mul(T, 0xf3672b51); // 0.32
	uint32_t r12 = 0xfd9d3a3e - r11; // 0.32
	uint32_t r13 = mul(T, 0x9a3c4390); // 0.32
	uint32_t r14 = 0xd4d2ce9b - r13; // 0.32
	uint32_t r15 = mul(r4, r14); // 0.32
	uint32_t r16 = r12 + r15; // 0.32
	uint32_t r17 = mul(r10, r16); // 2.30
	uint32_t r18 = r9 + r17; // 2.30
	uint32_t r19 = mul(r4, r4); // 0.32
	uint32_t r20 = mul(T, 0x1bba92b3); // 0.32
	uint32_t r21 = 0x525a1a8b - r20; // 0.32
	uint32_t r22 = mul(r4, 0x0452b1bf); // 0.32
	uint32_t r23 = r21 + r22; // 0.32
	uint32_t r24 = mul(r19, r23); // 0.32
	uint32_t r25 = mul(r10, r24); // 2.30
	uint32_t V = r18 + r25; // 2.30

	printf("r22 is : %x\n", r22);
	printf("r23 is : %x\n", r23);
	printf("r24 is : %x\n", r24);
	printf("r25 is : %x\n", r25);
	printf("V is : %x\n", V);


	// final packing of encoded float
	uint32_t tmp = Ex + Ey -nx -ny -MX -MY + c -110;
	uint32_t Lr = max( 0 , 0x0 -tmp);
	uint32_t U = ((V >> minu(Lr,25)) & 0xFFFFFFC0);//(Sr | (Dm1 << 23));
	uint32_t N = U << minu(Lr, 25);
	uint32_t P = mul(N, mpY);
	return P;
	//return (Sr | (Dm1 <<23)) + ((V >> minu(Lr,25)) & 0xFFFFFFC0);//(Sr | (Dm1 << 23));
}


int main(int argc, char *argv[]) {

	float X=0.0f;
	float Y=0.0f;
	X=atof(argv[1]);
	Y=atof(argv[2]);
	uint32_t res=0;
	res = div_digit_recurrence_restoring(X,Y);
	printf("result polynomial evaluation: %f\n", *(float*)&res);
	return 0;
}
