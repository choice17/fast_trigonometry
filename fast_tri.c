#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// https://zh.wikipedia.org/wiki/%E6%B3%B0%E5%8B%92%E7%BA%A7%E6%95%B0
// https://stackoverflow.com/posts/9073675/revisions

#define sign(x) ((x >= 0)? 1 : -1)
#define abs(x) ((x > 0)? (x) : -(x))

#ifdef DEB
#define logf(x) printf(#x ":%.6f ", x)
#define logfn(x) printf(#x ":%.6f\n", x)
#define logd(x) printf(#x ":%d ", x)
#define logdn(x) printf(#x ":%d\n", x)
#define log(fmt, x...) printf("%s %d " x, __func__,__LINE__, ##x) 
#else
#define logf(x)
#define logfn(x)
#define logd(x)
#define logdn(x)
#define log(fmt, x...)
#endif

#define fastsmall (0.05f)
#define fastpi (3.14159265358f)
#define fastpi_2 (1.5707963267948966f)
#define fastpi_4 (0.7853981633974483f)
#define fast2pi (6.283185307179586)

//const static float fastpi2 = 3.14159265358f * 3.14159265358f;

static inline int mod(int a, int b)
{
    return a - ((a/b)*b);
}

/* should be power of 2 */
static inline int mod2(int a, int b)
{
    return a & (b - 1);
}

static inline int fast_fmod2(float x, float m, int q, float *ax)
{
    float _x = x / m;
    int _ix = _x;
    *ax = (_x - _ix) * m;
    return mod2(_ix, q);
}

#define cubic(ax, b, a0, a1, a2) ((b) + (ax) * ((a0) + (ax) * ((a1) + (ax) * (a2)))) 
#define quad(ax, b, a0, a1) ((b) + (ax) * ((a0) + (ax) * (a1)))

#define cubicArr(ax, a) ((a[0]) + (ax) * ((a[1]) + (ax) * ((a[2]) + (ax) * (a[3])))) 
#define quadArr(ax, a) ((a[0]) + (ax) * ((a[1]) + (ax) * (a[2])))

static inline float lerp(float a, float b, float x)
{
    return a * x + b;   
}

float fast_atan(float x)
{
    const static float a0[4] = {
        -0.0037002278175466974,
        1.04710252,
        -0.18627361,
        -0.07857707
    };
    const static float a1[4] = { 
        -0.03835873083805175,
        1.22213379,
        -0.47333034,
        0.0749329
    };
    const static float a2[4] = { 
        0.23748654401677938,
        0.74395186,
        -0.19253416,
        0.01897149
    };
    const static float a3[4] = {
        0.6655143527424704,
        0.31987618,
        -0.05084753,
        0.00303765
    };
    const static float a4[4] = {
        1.029536240433053,
        0.10961277,
        -0.00974906,
        0.00032041
    };

    const static float atan_high_thr     = 10.0f;
    const static float atan_higl_thr     = 5.0f;
    const static float atan_midh_thr     = fastpi;
    const static float atan_mid_thr      = fastpi_2;
    const static float atan_mids_thr     = fastpi_4;
    const static float atan_smal_thr     = 0.1;

    int signX = sign(x);
    float ax = abs(x);

    // https://math.stackexchange.com/questions/982838/asymptotic-approximation-of-the-arctangent/982859
    if (ax > atan_high_thr)       return (signX*(fastpi_2) - 1/x);
    // for small x
    else if (ax < atan_smal_thr)  return x;
    // by regression
    else if (ax < atan_mids_thr)  return signX * cubicArr(ax, a0);
    // by regression
    else if (ax < atan_mid_thr)   return signX * cubicArr(ax, a1);
    // by regression
    else if (ax < atan_midh_thr)  return signX * cubicArr(ax, a2);
    // by regression
    else if (ax < atan_higl_thr)  return signX * cubicArr(ax, a3);
    else                          return signX * cubicArr(ax, a4);
} 

float fast_cos_quad2(float x)
{
    const static float c0=0.9994880721226505;
    const static float c1=-0.49574695;
    const static float c2=0.03677138;
    float x2 = x * x;
    return quad(x2, c0, c1, c2);
    //return 1 + x2 * (-0.5 + 0.041666666666666664 * x2);
}

float fast_cos(float x)
{
    float ax = abs(x);
    //if (ax < fastsmall) return x;
    int quad = fast_fmod2(ax, fastpi_2, 4, &ax);

    if      (quad == 0) return fast_cos_quad2(ax);
    else if (quad == 1) return -fast_cos_quad2(fastpi_2 - ax);
    else if (quad == 2) return -fast_cos_quad2(ax);
    else if (quad == 3) return fast_cos_quad2(fastpi_2 - ax);
    else assert(0);
}

float fast_sin(float x)
{
    log("%.6f\n", x);
    return fast_cos(fastpi_2 - x);
}

float fast_tan(float x)
{
    /* https://mae.ufl.edu/~uhk/ACCURATE-TANGENT.pdf */
    const static float t0[] = { // 0.05 to 0.8
        -0.005306908868911164,
        1.07958615,
        -0.31961611,
        0.72791356
    };
    const static float t1[] = { // 0.8 to 1.2
        -6.798523672037293,
        24.20035856,
        -26.75328616,
        10.9062448
    };
    const static float t2[] = { // 1.2 to 1.37
        -288.476741666932,
        715.84338056,
        -593.63206871,
        166.00749918
    };
    const static float t3[] = { // 1.37 to 1.45
        -4043.7026279735396,
        8896.89927632,
        -6536.19587117,
        1605.2422637
    };
    const static float t4[] = { // 1.45 to 1.472
        1513.2746348922383,
        -2142.93149953,
        762.05357123
    };
    const static float t5[] = { // 1.472 to 1.571
       8.23809275e+00, 8.32297410e+00, 8.40959713e+00, 8.49801630e+00,
       8.58828831e+00, 8.68047230e+00, 8.77462995e+00, 8.87082562e+00,
       8.96912650e+00, 9.06960275e+00, 9.17232772e+00, 9.27737805e+00,
       9.38483395e+00, 9.49477931e+00, 9.60730201e+00, 9.72249410e+00,
       9.84045207e+00, 9.96127709e+00, 1.00850754e+01, 1.02119584e+01,
       1.03420433e+01, 1.04754533e+01, 1.06123178e+01, 1.07527733e+01,
       1.08969634e+01, 1.10450395e+01, 1.11971613e+01, 1.13534974e+01,
       1.15142258e+01, 1.16795348e+01, 1.18496236e+01, 1.20247029e+01,
       1.22049964e+01, 1.23907409e+01, 1.25821881e+01, 1.27796052e+01,
       1.29832767e+01, 1.31935050e+01, 1.34106127e+01, 1.36349437e+01,
       1.38668654e+01, 1.41067703e+01, 1.43550786e+01, 1.46122405e+01,
       1.48787389e+01, 1.51550926e+01, 1.54418593e+01, 1.57396399e+01,
       1.60490821e+01, 1.63708858e+01, 1.67058076e+01, 1.70546676e+01,
       1.74183553e+01, 1.77978378e+01, 1.81941682e+01, 1.86084950e+01,
       1.90420739e+01, 1.94962800e+01, 1.99726227e+01, 2.04727622e+01,
       2.09985289e+01, 2.15519462e+01, 2.21352559e+01, 2.27509490e+01,
       2.34018014e+01, 2.40909153e+01, 2.48217690e+01, 2.55982756e+01,
       2.64248532e+01, 2.73065087e+01, 2.82489401e+01, 2.92586589e+01,
       3.03431416e+01, 3.15110147e+01, 3.27722851e+01, 3.41386275e+01,
       3.56237468e+01, 3.72438401e+01, 3.90181892e+01, 4.09699325e+01,
       4.31270796e+01, 4.55238646e+01, 4.82025764e+01, 5.12160763e+01,
       5.46313220e+01, 5.85344018e+01, 6.30378929e+01, 6.82918980e+01,
       7.45011028e+01, 8.19520672e+01, 9.10587111e+01, 1.02441914e+02,
       1.17077345e+02, 1.36591117e+02, 1.63910235e+02, 2.04888709e+02,
       2.73185895e+02, 4.09779859e+02, 8.19560938e+02, 1.63312394e+16
   };
    const static float tan0_th = 0.8;
    const static float tan1_th = 1.2;
    const static float tan2_th = 1.37;
    const static float tan3_th = 1.45;
    const static float tan4_th = 1.472;

    float _x;
    int _ix, _quad, x2;

    int signX = sign(x);
    float ax = abs(x);

    if (ax < fastsmall)
        return x;
    if (ax > fastpi_2) {
        _quad = fast_fmod2(ax, fastpi_2, 2, &ax);
        if (signX >= 0) { // > 0
            if (_quad == 1) {
                ax = fastpi_2 - ax;
                signX = -1;
            } else signX = 1;
        } else {  // < 0
            if (_quad == 0) {}
            else {
                ax = fastpi_2 - ax;
                signX = 1;
            }
                        printf("MOD>> %.4f => %.4f Q[%d] => %s%.4f\n",
                x,ax,_quad,signX>=0?"+":"-",ax);
        }
        if (ax < fastsmall) {
            return signX * ax;
        }
    }
    if (ax < tan0_th)
        return signX * cubicArr(ax, t0);
    else if (ax < tan1_th)
        return signX * cubicArr(ax, t1);
    else if (ax < tan2_th)
        return signX * cubicArr(ax, t2);
    else if (ax < tan3_th)
        return signX * cubicArr(ax, t3);
    else if (ax < tan4_th)
        return signX * quadArr(ax, t4);
    else {
        // 1.471 to 1.571
        float ax_up = ax * 1000;
        int iax_up = ax_up;
        int idx = iax_up - 1472;
        float r = ax_up - iax_up;
        //printf(">>>%.3f=>[%.3f]Q[%d] %.4f %.4f rate:%.4f => %.4f\n",
        //    x, ax, _quad, t5[idx], t5[idx+1], r, signX * lerp(t5[idx], t5[idx+1], r));
        return signX * lerp(t5[idx], t5[idx+1], r);
    }
}

float fast_asin(float x)
{
  float negate = x < 0;
  float ax = abs(x);
  float ret = -0.0187293;
  ret *= ax;
  ret += 0.0742610;
  ret *= ax;
  ret -= 0.2121144;
  ret *= ax;
  ret += 1.5707288;
  ret = fastpi * 0.5 - sqrt(1.0 - ax) * ret;
  return ret - 2 * negate * ret;
}

//https://stackoverflow.com/questions/3380628/fast-arc-cos-algorithm
float fast_acos(float x)
{
  #if 0
  float negate = x < 0;
  float ax = abs(x);
  float ret = -0.0187293;
  ret = ret * ax;
  ret = ret + 0.0742610;
  ret = ret * ax;
  ret = ret - 0.2121144;
  ret = ret * ax;
  ret = ret + 1.5707288;
  ret = ret * sqrt(1.0-ax);
  ret = ret - 2 * negate * ret;
  return negate * 3.14159265358979 + ret;
  
    static const float C  = 0.10501094f;
    float r, s, t, u;
    t = (x < 0) ? (-x) : x;  // handle negative arguments
    u = 1.0f - t;
    s = sqrtf (u + u);
    r = C * u * s + s;  // or fmaf (C * u, s, s) if FMA support in hardware
    if (x < 0) r = fastpi - r;  // handle negative arguments
    return r;
  #endif
  return fastpi_2 - fast_asin(x);
}


#define MATH_OPS(x) (double (*)(double))(x)
#define OPS_NUM (6)
int app(int argc, char **argv)
{
    float b = atof(argv[2]);
    int range = atoi(argv[3]);
    double (*math_ops[OPS_NUM])(double) = { MATH_OPS(atan), MATH_OPS(sin), MATH_OPS(cos), MATH_OPS(tan), MATH_OPS(acos), MATH_OPS(asin)};
    float (*fast_math[OPS_NUM])(float) = {fast_atan, fast_sin, fast_cos, fast_tan, fast_acos, fast_asin};
    const char *metstr[] = {"atan", "sin", "cos", "tan", "acos", "asin"};
    unsigned int met = atoi(argv[4]);
    assert(met < OPS_NUM+1);

    float a = 0.0;
    if (strcmp(argv[1], "-0")==0) {
        for (int i = 0; i < range; i++) {
            a = math_ops[met](b);
        }
    } else {
        for (int i = 0; i < range; i++) {
            a = fast_math[met](b);
        }
    }
    printf("selected %s\n", metstr[met]);
    printf("retval(%2.2f) = %.5f\n\n", b, a);
    printf("atan(%2.2f) = %.5f\n", b, atan(b));
    printf("sin(%2.2f) = %.5f\n", b, sin(b));
    printf("cos(%2.2f) = %.5f\n", b, cos(b));
    printf("tan(%2.2f) = %.5f\n", b, tan(b));
    printf("acos(%2.2f) = %.5f\n", b, acos(b));
    printf("asin(%2.2f) = %.5f\n", b, asin(b));
    printf("fast_atan(%2.2f) = %.5f\n", b, fast_atan(b));
    printf("fast_sin(%2.2f) = %.5f\n", b, fast_sin(b));
    printf("fast_cos(%2.2f) = %.5f\n", b, fast_cos(b));
    printf("fast_tan(%2.2f) = %.5f\n", b, fast_tan(b));
    printf("fast_acos(%2.2f) = %.5f\n", b, fast_acos(b));
    printf("fast_asin(%2.2f) = %.5f\n", b, fast_asin(b));
    return 0;
}

static inline float max(float x, float y)
{
    return (x > y)? x: y;
}
static inline float min(float x, float y)
{
    return (x < y)? x: y;
}
static inline float absf(float x)
{
    return (x < 0)? -x: x;
}

int ut_test(int argc, char **argv)
{
    float b = atof(argv[2]);
    int range = atoi(argv[3]);
    double (*math_ops[OPS_NUM])(double) = { MATH_OPS(atan), MATH_OPS(sin), MATH_OPS(cos), MATH_OPS(tan), MATH_OPS(acos), MATH_OPS(asin)};
    float (*fast_math[OPS_NUM])(float) = {fast_atan, fast_sin, fast_cos, fast_tan, fast_acos, fast_asin};
    const char *metstr[] = {"atan", " sin", " cos", " tan", "acos", "asin"};
    //unsigned int met = atoi(argv[4]);
    //assert(met < 5);

    float _max[OPS_NUM] = { 0 };
    float _min[OPS_NUM] = { 0 };
    float _max_val[OPS_NUM] = { 0 };
    float _min_val[OPS_NUM] = { 0 };
    float _sum[OPS_NUM] = { 0 };
    float _mean[OPS_NUM] = { 0 };
    float *_diff = (float *)malloc(sizeof(float)*range*4);
    float inc = b/range;
    float in = -b/2;
    for (int i = 0; i < range; i++) {
        for (int j = 0; j < OPS_NUM; j++) {
            float val0 = math_ops[j](in);
            float val1 = fast_math[j](in);
            float diff = absf(val0-val1);
            _diff[i*4 + j] = diff;
            if (diff > _max[j]) {
                _max[j] = diff;
                _max_val[j] = in;
            }
            if (diff < _min[j]) {
                _min[j] = diff;
                _min_val[j] = in;
            }
            _sum[j] += diff;
        }
        in += inc;
    }
    for (int i = 0; i < OPS_NUM; i++) {
        _mean[i] = _sum[i] / range;
    }
    printf("========= unit-test acc summary ========\n\n");
    printf(" test range [%.6f --- %.6f]\n", -b/2 , b/2 );
    printf("                 math vs fast_tri\n");
    for (int i = 0; i < OPS_NUM; i++) {
       printf("%8s:: max:%.6f(%.4f) min:%.6f(%.4f) mean:%.6f\n",
        metstr[i], _max[i], _max_val[i], _min[i], _min_val[i], _mean[i]);
    }
    FILE *fp = fopen("diff.bin", "wb");
    fwrite(_diff, 1, sizeof(float)*range*4, fp);
    fclose(fp);
    free(_diff);
    return 0;
}

void help()
{
    printf(
        "===============\n"
        "\targc:5 => [-0(math)/else(fast_math)] [input float] [iters] [fast_ops_method]\n"
        "\targc:4 => [dummy] [range * -pi/2 to pi/2] [arr size]\n"
        );
}
int main(int argc, char **argv)
{
 
    double (*math)(double) = atan;    
    if (argc == 5)
        app(argc, argv);
    else if (argc == 4)
        ut_test(argc, argv);
    else
        help();
    return 0;
}

/* time ./app {-0} 7.2 100000000 {0}
| | math | fast math |
| atan | 2.501 |  0.526 |
| sin | 5.367 | 1.944 |
| cos | 4.626 | 1.700 |
| tan | | |
*/