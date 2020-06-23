
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
#define NegFusedMulAddSIMD(a,b,c) _mm512_fnmadd_pd((a),(b),(c))
#define NegFusedMulSubSIMD(a,b,c) _mm512_fnmsub_pd((a),(b),(c))
// In the case of 512-bit SIMD:
//    The result from this comparison is: result[i] = (a OP b) ? 1 : 0, stored in an 8-bit mask array.
//    Then if result==1 we set upwind = 0+1, and if result==0 we set upwind = 0
#define UPWIND_ALG(a) _mm512_mask_add_pd(upwind_Integer_0,  _mm512_cmp_pd_mask( (a) , upwind_Integer_0, _CMP_GT_OQ), upwind_Integer_0 ,upwind_Integer_1)

// If compiled with AVX SIMD instructions enabled:
#elif __AVX__
#include <immintrin.h>
#define REAL_SIMD_ARRAY __m256d
#define SIMD_width 4 // 4 doubles per loop iteration
// Upwind algorithm notes, in the case of 256-bit SIMD:
// Sources: https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_mm256_cmp_pd&expand=736
//  ...and: https://stackoverflow.com/questions/37099874/is-avx-intrinsic-mm256-cmp-ps-supposed-to-return-nan-when-true
//     The result from _mm256_cmp_pd is 0 if a>0 and NaN otherwise.
//     Thus if OP is >, then: if a > b then the result is NaN, and if a <= b then the result is 0.
//     We want the result to be 1 if a>b and 0 otherwise, so we simply perform a logical AND operation
//     on the result, against the number 1, because AND(NaN,1)=1, and AND(0,1)=0,
//     where NaN=0xffffff... in double precision.
#define UPWIND_ALG(a) _mm256_and_pd(_mm256_cmp_pd( (a), upwind_Integer_0, _CMP_GT_OQ ), upwind_Integer_1)
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
#define NegFusedMulAddSIMD(a,b,c) _mm256_fnmadd_pd((a),(b),(c))
#define NegFusedMulSubSIMD(a,b,c) _mm256_fnmsub_pd((a),(b),(c))
#else
#define FusedMulAddSIMD(a,b,c) _mm256_add_pd(_mm256_mul_pd((a),(b)), (c)) // a*b+c
#define FusedMulSubSIMD(a,b,c) _mm256_sub_pd(_mm256_mul_pd((a),(b)), (c)) // a*b-c
#define NegFusedMulAddSIMD(a,b,c) _mm256_sub_pd( (c), _mm256_mul_pd((a),(b)) ) // c-a*b
// NegFusedMulSubSIMD(a,b,c) = -a*b - c
//                           = c - (c+a*b+c)
//                           = SubSIMD(c, AddSIMD(c, AddSIMD(MulSIMD(a,b), c)))
#define NegFusedMulSubSIMD(a,b,c) _mm256_sub_pd( (c), _mm256_add_pd( (c), _mm256_add_pd( _mm256_mul_pd((a),(b)), (c) )))
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
// See description above UPWIND_ALG for __AVX__:
#define UPWIND_ALG(a) _mm_and_pd(_mm_cmpgt_pd( (a), upwind_Integer_0 ), upwind_Integer_1)

#ifdef __FMA__ // There are no mainstream non-AVX+ chips that have FMA, but we include the following for completeness.
#define FusedMulAddSIMD(a,b,c) _mm_fmadd_pd((a),(b),(c))
#define FusedMulSubSIMD(a,b,c) _mm_fmsub_pd((a),(b),(c))
#define NegFusedMulAddSIMD(a,b,c) _mm_sub_pd( (c), _mm_mul_pd((a),(b)) ) // c-a*b
#define NegFusedMulSubSIMD(a,b,c) _mm_sub_pd( (c), _mm_sub_pd( (c), _mm_sub_pd( _mm_mul_pd((a),(b)), (c) ))) // -a*b-c = c-c-a*b-c = SubSIMD(c,SubSIMD(c,SubSIMD(MulSIMD(a,b),c)
#else
#define FusedMulAddSIMD(a,b,c) _mm_add_pd(_mm_mul_pd((a),(b)), (c)) // a*b+c
#define FusedMulSubSIMD(a,b,c) _mm_sub_pd(_mm_mul_pd((a),(b)), (c)) // a*b-c
#define NegFusedMulAddSIMD(a,b,c) _mm_sub_pd( (c), _mm_mul_pd((a),(b)) ) // c-a*b
// NegFusedMulSubSIMD(a,b,c) = -a*b - c
//                           = c - (c+a*b+c)
//                           = SubSIMD(c, AddSIMD(c, AddSIMD(MulSIMD(a,b), c)))
#define NegFusedMulSubSIMD(a,b,c) _mm_sub_pd( (c), _mm_add_pd( (c), _mm_add_pd( _mm_mul_pd((a),(b)), (c) )))
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
#define NegFusedMulAddSIMD(a,b,c) ((c) - (a)*(b))
#define NegFusedMulSubSIMD(a,b,c) (-((a)*(b) + (c))) // -a*b-c = -(a*b+c)
#define SqrtSIMD(a) (sqrt(a))
#define ExpSIMD(a) (exp(a))
#define SinSIMD(a) (sin(a))
#define CosSIMD(a) (cos(a))
#define WriteSIMD(a,b) *(a)=(b)
#define ReadSIMD(a) *(a)
// Algorithm for upwinding, SIMD-disabled version.
// *NOTE*: This upwinding is backwards from
//  usual upwinding algorithms, because the
//  upwinding control vector in BSSN (the shift)
//  acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0
#endif
