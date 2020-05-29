#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define sign(x) ((x >= 0)? 1 : -1)
#define abs(x) ((x >= 0)? (x) : -(x))

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

const static float fastsmall = 0.05f;
const static float fastpi = 3.14159265358f;
const static float fastpi_2 = 3.14159265358f/2;
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

float fast_atan(float x)
{
    const static float a0 = 1.10382257f;
    const static float a1 = -0.36105955f;
    const static float a2 = 0.04400221f;
    const static float b = -0.00492334997743038f;
    int signX = sign(x);
    float ax = abs(x);
    if (ax > fastpi)
        // https://math.stackexchange.com/questions/982838/asymptotic-approximation-of-the-arctangent/982859
        return (signX*(fastpi_2) - 1/x);
    else if (ax < fastsmall)
        // for small x
        return x;
    else {
        // by regression
        return signX * (b + ax * ( a0 + ax * ( a1 + ax * a2)));
    }
} 

float fast_cos(float x)
{
    const float c0 = 0.9970072879558844f;
    const float c1 = 0.0385896;
    const float c2 = -0.61311836;
    const float c3 = 0.11737913;

    float _x; 
    int _ix, _quad;

    int signX = sign(x);
    float ax = abs(x);
    log("%.6f->%.6f < %.6f\n", x, ax, fastsmall);
    if (ax < fastsmall)
        return 1;
    else if (ax > fastpi_2) {
#if 1   /* for fixed point */
        _x = ax / fastpi_2;
        _ix = (int)_x;
        _quad = mod2(_ix, 4);
        ax = (_x - _ix) * fastpi_2;
#else   
        _quad = mod(ax / fastpi_2, 4);
#endif
        signX = 1;
        if (_quad == 1) {
            ax = fastpi_2 - ax;
            signX = -1;
        } else if (_quad == 2) {
            signX = -1;
        } else if (_quad == 3) {
            ax = fastpi_2 - ax;
        }
        if (ax < fastsmall)
            return signX;
        logd(signX);logd(_quad);logdn(_ix);
        logf(_x); logfn(ax);
    } else
        signX = 1;
    return signX * (c0 + ax * (c1 + ax * (c2 + ax * c3)));
}

float fast_sin(float x)
{
    log("%.6f\n", x);
    return fast_cos(fastpi_2 - x);
}

float fast_tan(float x)
{
    /* https://mae.ufl.edu/~uhk/ACCURATE-TANGENT.pdf */
    const float t0 = 0.0031084695854767252f;
    const float t1 = 0.9140496f;
    const float t2 = 0.53540681;
    const float t3 = -0.85751509;
    const float t4 = 0.9570163;

    float _x;
    int _ix, _quad, x2;

    int signX = sign(x);
    float ax = abs(x);

    if (ax < fastsmall)
        return x;
    if (ax > 1) {
        _x = ax / fastpi_2;
        _ix = (int)_x;
        _quad = mod2(_ix, 2);
        if (_quad == 1) {
            ax = fastpi_2 - ax;
            signX = -1;
        } else
            signX = 1;
        if (ax > 1) {
            ax = (_x - _ix) * fastpi_2;
            x2 = ax * ax;
            return signX * (ax * (945 - x2 * (105 - x2))) / (945 - x2 * (420 - 15 * x2));
        }
    }
    return signX * (t0 + ax * (t1 + ax * (t2 + ax * t3 + ax * t4)));
}


int app(int argc, char **argv)
{
    float b = atof(argv[2]);
    int range = atoi(argv[3]);
    double (*math_ops[4])(double) = { (double *)atan, (double *)sin, (double *)cos, (double *)tan};
    float (*fast_math[4])(float) = {fast_atan, fast_sin, fast_cos, fast_tan};
    const char *metstr[] = {"atan", "sin", "cos", "tan"};
    unsigned int met = atoi(argv[4]);
    assert(met < 5);

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
    printf("fast_atan(%2.2f) = %.5f\n", b, fast_atan(b));
    printf("fast_sin(%2.2f) = %.5f\n", b, fast_sin(b));
    printf("fast_cos(%2.2f) = %.5f\n", b, fast_cos(b));
    printf("fast_tan(%2.2f) = %.5f\n", b, fast_tan(b));
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
    double (*math_ops[4])(double) = { (double *)atan, (double *)sin, (double *)cos, (double *)tan};
    float (*fast_math[4])(float) = {fast_atan, fast_sin, fast_cos, fast_tan};
    const char *metstr[] = {"atan", " sin", " cos", " tan"};
    //unsigned int met = atoi(argv[4]);
    //assert(met < 5);

    float _max[4] = { 0 };
    float _min[4] = { 0 };
    float _max_val[4] = { 0 };
    float _min_val[4] = { 0 };
    float _sum[4] = { 0 };
    float _mean[4] = { 0 };
    float *_diff = (float *)malloc(sizeof(float)*range*4);
    float inc = b*fastpi/range;
    float in = -b/2 * fastpi;
    for (int i = 0; i < range; i++) {
        for (int j = 0; j < 4; j++) {
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
    for (int i = 0; i < 4; i++) {
        _mean[i] = _sum[i] / range;
    }
    printf("========= unit-test acc summary ========\n\n");
    printf(" test range [%.6f --- %.6f]\n", -b/2 * fastpi, b/2 * fastpi);
    printf("                 math vs fast_tri\n");
    for (int i = 0; i < 4; i++) {
       printf("    %s:: max:%.6f(%.4f) min:%.6f(%.4f) mean:%.6f\n",
        metstr[i], _max[i], _max_val[i], _min[i], _min_val[i], _mean[i]);
    }
    FILE *fp = fopen("diff.bin", "wb");
    fwrite(_diff, 1, sizeof(float)*range*4, fp);
    fclose(fp);
    free(_diff);
    return 0;
}

int main(int argc, char **argv)
{
    
    if (argc == 5)
        app(argc, argv);
    else
        ut_test(argc, argv);
    return 0;
}

/* time ./app {-0} 7.2 100000000 {0}
| | math | fast math |
| atan | 2.501 |  0.526 |
| sin | 5.367 | 1.944 |
| cos | 4.626 | 1.700 |
| tan | | |
*/