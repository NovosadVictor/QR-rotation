#ifndef QR_FUNCTIONS_H
#define QR_FUNCTIONS_H

#include <math.h>


// Set vector of sins
void set_vector(int n, int k, double *matrix, double *d1, double *d2) {
    matrix[k * n + k - 1] *= matrix[k * n + k - 1];

    int is_zeros = 1;
    for (int i = k + 1; i < n; i++)
       if (fabs(matrix[i * n + k - 1]) > exp(-15)) {
           is_zeros = 0;
           break;
    }
    if (is_zeros == 1) {
        for (int i = 0; i < n - 1; i++) {
            d2[i] = 1.0;
            d1[i] = 0.0;
        }
        return ;
    }

    for (int i = k + 1; i < n; i++) {
        double x = matrix[k * n + k - 1];
        matrix[k * n + k - 1] += matrix[i * n + k - 1] * matrix[i * n + k - 1];
        double sin = matrix[i * n + k - 1] / sqrt(matrix[k * n + k - 1]);
        double cos = sqrt(x) / sqrt(matrix[k * n + k - 1]);
//        printf("\nk = %d, cos = %lf, sin = %lf\n", k, cos, sin);
        matrix[i * n + k - 1] = 0.0;
        d1[i - k - 1] = -sin;
        d2[i - k - 1] = cos;
    }

    matrix[k * n + k - 1] = sqrt(matrix[k * n + k - 1]);
}


// Matrix to almost triangle matrix
void to_almost_triangle(int n, double *matrix, double*d1, double *d2) {
    if (n == 2)
        return;
    for (int k = 1; k < n - 1; k++) {
        set_vector(n, k, matrix, d1, d2);
        for (int i = k + 1; i < n; i++) {
            for (int j = k; j < n; j++) {
                double matrix_k = matrix[k * n + j];
                double matrix_i = matrix[i * n + j];
                matrix[k * n + j] = matrix_k * d2[i - k - 1] - matrix_i * d1[i - k - 1];
                matrix[i * n + j] = matrix_k * d1[i - k - 1] + matrix_i * d2[i - k - 1];
            }
/*            printf("\nmatrix after mult left, k = %d\n", k);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++)
                    printf("%lf ", matrix[i * n + j]);
                printf("\n");
            }
            printf("\n");
*/
        }
        for (int l = k + 1; l < n; l++) {
            for (int i = k; i < n; i++) {
                double matrix_k = matrix[i * n + k];
                double matrix_i = matrix[i * n + l];
                matrix[i * n + k] = matrix_k * d2[l - k - 1] - matrix_i * d1[l - k - 1];
                matrix[i * n + l] = matrix_i * d2[l - k - 1] + matrix_k * d1[l - k - 1];
            }
/*            printf("\nmatrix after mult right, k = %d\n", k);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++)
                    printf("%lf ", matrix[i * n + j]);
                printf("\n");
            }
            printf("\n");
*/
        }
    }
}


//QR for almost triangle matrix
void qr(int n, double *matrix, double *d1, double *d2) {
    for (int k = 1; k < n; k++) {
        double x = matrix[(k - 1) * n + k - 1];
        if (fabs(x) < exp(-15) && fabs(matrix[k * n + k - 1]) < exp(-15)) {
            d1[k - 1] = 0.0;
            d2[k - 1] = 1.0;
            continue;
        }
        matrix[(k - 1) * n + k - 1] *= matrix[(k - 1) * n + k - 1];
        matrix[(k - 1) * n + k - 1] += matrix[k * n + k - 1] * matrix[k * n + k - 1];
        double sin = matrix[k * n + k - 1] / sqrt(matrix[(k - 1) * n + k - 1]);
        double cos = x / sqrt(matrix[(k - 1) * n + k - 1]);
        d1[k - 1] = -sin;
        d2[k - 1] = cos;
//        printf("qr k = %d, cos = %lf, sin = %lf\n", k, cos, sin);
        matrix[(k - 1) * n + k - 1] = sqrt(matrix[(k - 1) * n + k - 1]);
        matrix[k * n + k - 1] = 0.0;
    }

/*    printf("QR for almost triangle\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%lf ", matrix[i * n + j]);
        printf("\n");
    }
    printf("\n");
*/
}


// From QR to RQ
void set_rq(int n, double *matrix, double *d1, double *d2) {
    for (int k = 0; k < n - 1; k++) {
        for (int i = 0; i < n; i++) {
            double matrix_k = matrix[i * n + k];
            double matrix_i = matrix[i * n + k + 1];
            matrix[i * n + k] = matrix_k * d2[k] - matrix_i * d1[k];
            matrix[i * n + k + 1] = matrix_i * d2[k] + matrix_k * d1[k];
        }
/*        printf("RQ k = %d\n", k);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                printf("%lf ", matrix[i * n + j]);
            printf("\n");
        }
        printf("\n");
*/
    }
}


// Finding eigenvalues
void eigenvalues(int n, double *matrix, double *d1, double *d2) {
    // !!!!NEED FIX!!!!
    for (int s = 0; s < 20; s++) {
        to_almost_triangle(n, matrix, d1, d2);
        qr(n, matrix, d1, d2);
        set_rq(n, matrix, d1, d2);
        printf("RQ iteration %d\n", s);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                printf("%lf ", matrix[i * n + j]);
            printf("\n");
        }
        printf("\n");
    }

    for (int i = 0; i < n - 1; i++)
        d1[i] = matrix[i * n + i];
    d2[0] = matrix[(n - 1) * n + (n - 1)];
}
















#endif // QR_FUNCTIONS_H
