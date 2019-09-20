#include <stdio.h>
#include <math.h>
#include <cassert>
#include <cstdlib>

const int INFINITE_ROOTS = -1;
const double EPS = 1e-7;
const int MAX_TRIES = 100;


/**
 * @brief solves a quadratic equation Ax^2 + Bx + C = 0.
 * @param a x^2 coefficient of equation.
 * @param b x coefficient of equation.
 * @param c free coefficient of equation.
 * @param x1 pointer to the root.
 * @param x2 pointer to the root.
 * @return returns number of roots. In case of any x, returns INFINITE_ROOTS
 */

int solve_equation(double a, double b, double c, double* x1, double* x2);

/**
 * @param text Invite string for user.
 * @param size Number of values in input.
 * @return pointer to calloc() generated array, containing input values. Free after using.
 */

double* input(const char text[], size_t size);

/**
 * @brief compare double and integer value
 * @param d double value
 * @param i integer value
 * @return 1 if equal, 0 if not equal
 */

int is_equal(double d, int i);

/**
 * @brief read coefficients of quadratic equation and print roots.
 * @return exit code 1: impossible number of roots.
 */

int main() {
    double* coeff = input("Enter coefficients",  3);
    double x1 = NAN, x2 = NAN;
    int nRoots = solve_equation(coeff[0], coeff[1], coeff[2], &x1, &x2);
    switch(nRoots) {
        case 0:
            printf("No roots");
            break;
        case 1:
            printf("Number of roots: 1\nx = %g", (is_equal(x1, 0)) ? 0 : x1);
            break;
        case 2:
            printf("Number of roots: 2\nx1 = %g\nx2 = %g", (is_equal(x1, 0)) ? 0 : x1, (is_equal(x2, 0)) ? 0 : x2);
            break;
        case INFINITE_ROOTS:
            printf("Infinite number of solutions");
            break;
        default:
            fprintf(stderr, "main(): invalid number of roots");
            return 1;
    }
    free(coeff);
    return 0;
}

int solve_equation(double a, double b, double c, double* x1, double* x2) {
    assert(isfinite(a));
    assert(isfinite(b));
    assert(isfinite(c));
    assert(x1 != NULL);
    assert(x2 != NULL);
    assert(x1 != x2);
    if(is_equal(a, 0)) {
        if(is_equal(a,0)) {
            if(is_equal(a,0)) {
                return INFINITE_ROOTS;
            }
            return 0;
        }
        *x1 = -c / b;
        return 1;
    }
    double discriminant = b * b - 4 * a * c;
    if(discriminant < 0 ) {
        return 0;
    }
    if(is_equal(discriminant, 0)) {
        *x1 = -b / 2 / a;
        return 1;
    }
    //if(discriminant > 0)
    *x1 = (-b + sqrt(discriminant)) / 2 / a;
    *x2 = (-b - sqrt(discriminant)) / 2 / a;
    return 2;
}

double* input(const char text[], size_t size){
    printf("%s\n", text);
    double* data = (double*)calloc(size, sizeof(data[0]));
    for (int i = 0; i < size; ++i) {
        printf("[%d / %d]\n", i + 1, size);
        for (int j = 0; j < MAX_TRIES; ++j) {
            if (scanf("%lf", &data[i]) == 1) {
                break;
            }
            printf("Try enter [%d / %d] again\n", i + 1, size);
            fflush(stdin);
        }
    }
    return data;
}

int is_equal(double d, int i) {
    return (abs(d - i) < EPS) ? 1 : 0;
}