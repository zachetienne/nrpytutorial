
// If compiled with AVX512F SIMD instructions enabled:
#ifdef __AVX512F__
#include <immintrin.h>
#define REAL_SIMD_ARRAY __m512d
#define SIMD_width 8 // 8 doubles per loop iteration
#define ReadSIMD(a) _mm512_loadu_pd(a)
#define WriteSIMD(a,b) _mm512_storeu_pd(a,(b))
#define ConstSIMD(a) _mm512_set1_pd(a)
#define AddSIMD(a,b) _mm512_add_pd((a),(b))
#define SubSIMD(a,b) _mm512_sub_pd((a),(b))
#define MulSIMD(a,b) _mm512_mul_pd((a),(b))
#define DivSIMD(a,b) _mm512_div_pd((a),(b))
#define ExpSIMD(a) _mm512_exp_pd((a))
#define SinSIMD(a) _mm512_sin_pd((a))
#define CosSIMD(a) _mm512_cos_pd((a))
// All AVX512 chips have FMA enabled
#define FusedMulAddSIMD(a,b,c) _mm512_fmadd_pd((a),(b),(c))
#define FusedMulSubSIMD(a,b,c) _mm512_fmsub_pd((a),(b),(c))

// If compiled with AVX SIMD instructions enabled:
#elif __AVX__
#include <immintrin.h>
#define REAL_SIMD_ARRAY __m256d
#define SIMD_width 4 // 4 doubles per loop iteration
#define ReadSIMD(a) _mm256_loadu_pd(a)
#define WriteSIMD(a,b) _mm256_storeu_pd(a,(b))
#define ConstSIMD(a) _mm256_set1_pd(a)
#define AddSIMD(a,b) _mm256_add_pd((a),(b))
#define SubSIMD(a,b) _mm256_sub_pd((a),(b))
#define MulSIMD(a,b) _mm256_mul_pd((a),(b))
#define DivSIMD(a,b) _mm256_div_pd((a),(b))
#define ExpSIMD(a) _mm256_exp_pd((a))
#define SinSIMD(a) _mm256_sin_pd((a))
#define CosSIMD(a) _mm256_cos_pd((a))
#ifdef __FMA__
#define FusedMulAddSIMD(a,b,c) _mm256_fmadd_pd((a),(b),(c))
#define FusedMulSubSIMD(a,b,c) _mm256_fmsub_pd((a),(b),(c))
#else
#define FusedMulAddSIMD(a,b,c) _mm256_add_pd(_mm256_mul_pd((a),(b)), (c)) // a*b+c
#define FusedMulSubSIMD(a,b,c) _mm256_sub_pd(_mm256_mul_pd((a),(b)), (c)) // a*b-c
#endif


// If compiled with SSE2 SIMD instructions enabled:
#elif __SSE2__
#include <emmintrin.h>
#define REAL_SIMD_ARRAY __m128d
#define SIMD_width 2 // 2 doubles per loop iteration
#define ReadSIMD(a) _mm_loadu_pd(a)
#define WriteSIMD(a,b) _mm_storeu_pd(a,(b))
#define ConstSIMD(a) _mm_set1_pd(a)
#define AddSIMD(a,b) _mm_add_pd((a),(b))
#define SubSIMD(a,b) _mm_sub_pd((a),(b))
#define MulSIMD(a,b) _mm_mul_pd((a),(b))
#define DivSIMD(a,b) _mm_div_pd((a),(b))
#define ExpSIMD(a) _mm_exp_pd((a))
#define SinSIMD(a) _mm_sin_pd((a))
#define CosSIMD(a) _mm_cos_pd((a))
#ifdef __FMA__ // Unlikely that any SSE2 chip has FMA, but this will work fine.
#define FusedMulAddSIMD(a,b,c) _mm_fmadd_pd((a),(b),(c))
#define FusedMulSubSIMD(a,b,c) _mm_fmsub_pd((a),(b),(c))
#else
#define FusedMulAddSIMD(a,b,c) _mm_add_pd(_mm_mul_pd((a),(b)), (c)) // a*b+c
#define FusedMulSubSIMD(a,b,c) _mm_sub_pd(_mm_mul_pd((a),(b)), (c)) // a*b-c
#endif

#else
// If SIMD instructions unavailable:
#define REAL_SIMD_ARRAY REAL
#define SIMD_width 1 // 1 double per loop iteration
#define ConstSIMD(a) (a)
#define AddSIMD(a,b) ((a)+(b))
#define SubSIMD(a,b) ((a)-(b))
#define MulSIMD(a,b) ((a)*(b))
#define DivSIMD(a,b) ((a)/(b))
#define FusedMulAddSIMD(a,b,c) ((a)*(b) + (c))
#define FusedMulSubSIMD(a,b,c) ((a)*(b) - (c))
#define SqrtSIMD(a) (sqrt(a))
#define ExpSIMD(a) (exp(a))
#define SinSIMD(a) (sin(a))
#define CosSIMD(a) (cos(a))
#define WriteSIMD(a,b) *(a)=(b)
#define ReadSIMD(a) *(a)
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0
#endif
