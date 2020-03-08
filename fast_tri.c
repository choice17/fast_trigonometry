#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define sign(x) ((x >= 0)? 1 : -1)
#define abs(x) ((x >= 0)? (x) : -(x))

const static float fastpi = 3.14159265358f;
const static float fastpi_2 = 3.14159265358f/2;

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
    else if (ax < 0.05)
        // for small x
        return x;
    else {
        // by regression
        float x2 = ax * ax;
        float x3 = x2 * ax;
        return signX * (b + a0 * ax + a1 * x2 + a2 * x3);
    }
} 

int main(int argc, char **argv)
{
    float b = atof(argv[2]);
    int range = atoi(argv[3]);
    float a = 0.0;
    if (strcmp(argv[1], "-0")==0) {
        for (int i = 0; i < range; i++) {
            a = atan(b);
        }
    } else {
        for (int i = 0; i < range; i++) {
                a = fast_atan(b);
        }
    }
    printf("atan(%2.2f) = %.5f\n", b, a);
}