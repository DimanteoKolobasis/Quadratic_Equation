#include <stdio.h>
#include <math.h>
#include <cassert>
#include <cstdlib>

#define TRUE 1
#define FALSE 0

const int INFINITE_ROOTS = -1;
const double EPS = 1e-7;
const int MAX_TRIES = 100;


/**
 * @brief solves a quadratic equation Ax^2 + Bx + C = 0.
 * @param a x^2 coefficient of equation.
 * @param b x coefficient of equation.
 * @param c free coefficient of equation.
 * @param x1 pointer to bigger root.
 * @param x2 pointer to smaller root.
 * @return returns number of roots. In case of any x, returns INFINITE_ROOTS
 */

int solve_equation(double a, double b, double c, double* x1, double* x2);

/**
 * @param text Invite string for user.
 * @param size Number of values in input.
 * @return pointer to calloc() generated array, containing input values.
 * @attention Free after using.
 */

double* input(const char text[], size_t size);

/**
 * @brief tests function solve_equation() on random coefficients.
 * @return 0 if test passed. 2 if test failed.
 **/

int random_coeff_test();

/**
 * @brief tests function solve_equation on fixed coefficients.
 * @return 0 if test passed. 3 if test failed.
 */

int fixed_coeff_test();

/**
 * @brief compare double and integer value.
 * @param d double value.
 * @param i integer value.
 * @return 1 if equal, 0 if not equal.
 */

inline int is_equal(double d, int i);

/**
 * @brief generate random double value.
 * @return double value.
 */
inline double random_double();

/**
 * @brief read coefficients of quadratic equation and print roots.
 * @return exit code 1: impossible number of roots. exit code 2: failed on random test.
 */

int main() {
    if(random_coeff_test() == 1) {
        return 2;
    }
    if (fixed_coeff_test() == 3) {
        return 3;
    }
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
        if(is_equal(b,0)) {
            if(is_equal(c,0)) {
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
    double tmp_x1 = (-b + sqrt(discriminant)) / 2 / a;
    double tmp_x2 = (-b - sqrt(discriminant)) / 2 / a;
    (tmp_x1 > tmp_x2) ? (*x1 = tmp_x1, *x2 = tmp_x2) : (*x1 = tmp_x2, *x2 = tmp_x1);
    return 2;
}

int random_coeff_test() {
    int error_happened = FALSE;
    for (int i = 1; i <= 100; ++i) {
        double a = random_double();
        double b = random_double();
        double c = random_double();
        double x1 = NAN, x2 = NAN;
        int nRoots = solve_equation(a, b, c, &x1, &x2);
        switch (nRoots) {
            case 2: {
                double x1_value = a * x1 * x1 + b * x1 + c;
                double x2_value = a * x2 * x2 + b * x2 + c;
                if (!is_equal(x1_value, 0) || !is_equal(x2_value, 0)) {
                    error_happened = TRUE;
                    fprintf(stderr,
                            "Error in random test № %d. Equation not equal to 0.\nin: a = %g; b = %g; c = %g;\nout: x1 = %g; x2 = %g;\nf(x1) = %g; f(x2) = %g;\n",
                            i, a, b, c, x1, x2, x1_value, x2_value);
                }
                break;
            }
            case 1: {
                double value = a * x1 * x1 + b * x1 + c;
                double discriminant = b * b - 4 * a  * c;
                if (!is_equal(discriminant, 0)) {
                    error_happened = TRUE;
                    fprintf(stderr, "Error in random test № %d. Number of roots is 1, but discriminant not equal to 0.\nin: a = %g; b = %g; c = %g; discriminant = %g;\n",
                            i, a, b, c, discriminant);
                }
                if (!is_equal(value, 0)) {
                    error_happened = TRUE;
                    fprintf(stderr, "Error in random test № %d. Equation not equal to 0.\nin: a = %g; b = %g; c = %g;\nout: x1 = %g; x2 = none;\nf(x1) = %g;\n",
                            i, a, b, c, x1, value);
                }
                break;
            }
            case 0: {
                double discriminant = b * b - 4 * a * c;
                if (discriminant >= 0) {
                    error_happened = TRUE;
                    fprintf(stderr, "Error in random test № %d. No roots, but discriminant >= 0.\nin: a = %g; b = %g; c = %g; discriminant = %g\n",
                            i, a, b, c, discriminant);
                }
                break;
            }
            case INFINITE_ROOTS: {
                if (!is_equal(a, 0) || !is_equal(b, 0) || is_equal(c, 0)) {
                    error_happened = TRUE;
                    fprintf(stderr, "Error in random test № %d. Returned infinite number of roots, but coefficients not equal to zero.\nin: a = %g; b = %g; c = %g;\n");
                }
            }
        }
    }
    if (error_happened) {
        return 2;
    }
}

int fixed_coeff_test() {
    double one_root_coeff[4][4] = {
            {1, 2, 1, -1},
            {3, 6, 3, -1},
            {0, -4, 2, 0.5},
            {0, 3, 9, -3}
    };
    for (int i = 0; i < 4; ++i) {
        double x = NAN, null_x = NAN;
        int n_roots = solve_equation(one_root_coeff[i][0], one_root_coeff[i][1], one_root_coeff[i][2], &x, &null_x);
        if (x != one_root_coeff[i][3] || n_roots != 1) {
            fprintf(stderr, "Failed test: a = %g; b = %g; c = %g\n", one_root_coeff[i][0], one_root_coeff[i][1], one_root_coeff[i][2]);
            return 3;
        }
    }
    double two_root_coeff[4][5] = {
            {1, 16, -57, 3, -19},
            {1, 4.5, -22.5, 3, -7.5},
            {5, 23, 12, -0.6, -4},
            {60, 4560, 86100, -35, -41}
    };
    for (int i = 0; i < 4; ++i) {
        double x1 = NAN, x2 = NAN;
        int n_roots = solve_equation(two_root_coeff[i][0], two_root_coeff[i][1], two_root_coeff[i][2], &x1, &x2);
        if (x1 != two_root_coeff[i][3] || x2 != two_root_coeff[i][4] || n_roots != 2) {
            fprintf(stderr, "Failed test: a = %g; b = %g; c = %g\n", two_root_coeff[i][0], two_root_coeff[i][1], two_root_coeff[i][2]);
            return 3;
        }
    }
    double x1 = NAN, x2 = NAN;
    int n_roots = solve_equation(0, 0, 0, &x1, &x2);
    if (n_roots != INFINITE_ROOTS) {
        fprintf(stderr, "Failed infinite roots test: a = 0; b = 0; c = 0\n");
        return 3;
    }
    n_roots = solve_equation(1, 5, 8, &x1, &x2);
    if (n_roots != 0) {
        fprintf(stderr, "Failed zero root test: a = 1; b = 5; c = 8\n");
    }
    return 0;
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

inline double random_double() {
    return rand() + (double)rand() / (double)RAND_MAX - rand();
}

inline int is_equal(double d, int i) {
    return (abs(d - i) < EPS) ? 1 : 0;
}
