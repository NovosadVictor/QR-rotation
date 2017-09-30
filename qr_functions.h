#ifndef QR_FUNCTIONS_H
#define QR_FUNCTIONS_H

#include <math.h>


// Set vector of sins
void set_vector(int n, int k, double *matrix) {
    matrix[k * n + k - 1] *= matrix[k * n + k - 1];
    for (int i = k + 1; i < n; i++) {
        matrix[k * n + k - 1] += matrix[i * n + k - 1] * matrix[i * n + k - 1];
        double angle = matrix[i * n + k - 1] / sqrt(matrix[k * n + k - 1]);
        matrix[i * n + k - 1] = (-1) * angle;
    }

    matrix[k * n + k - 1] = sqrt(matrix[k * n + k - 1]);
}


// Matrix to almost triangle matrix
void to_almost_triangle(int n, double *matrix) {
    if (n == 2)
        return;
    for (int k = 1; k < n - 1; k++) {
        set_vector(n, k, matrix);
        for (int i = k + 1; i < n; i++) {
            for (int j = k; j < n; j++) {
                double matrix_k = matrix[k * n + j];
                double matrix_i = matrix[i * n + j];
                double cos = sqrt(1 - (matrix[i * n + (k - 1)] * matrix[i * n + (k - 1)]));
                matrix[k * n + j] = matrix_k * cos - matrix_i * matrix[i * n + (k - 1)];
                matrix[i * n + j] = matrix_k * matrix[i * n + (k - 1)] + matrix_i * cos;
            }
/*            for (int i = 0; i < n; i++) {
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
                double cos = sqrt(1 - (matrix[l * n + (k - 1)] * matrix[l * n + (k - 1)]));
                matrix[i * n + k] = matrix_k * cos + matrix_i * matrix[l * n + (k - 1)];
                matrix[i * n + l] = matrix_i * cos - matrix_k * matrix[l * n + (k - 1)];
            }
/*            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++)
                    printf("%lf ", matrix[i * n + j]);
                printf("\n");
            }
            printf("\n");
*/
        }
    }

    for (int j = 0; j < n - 2; j++)
        for (int i = j + 2; i < n; i ++)
            matrix[i * n + j] = 0.0;

}


//QR for almost triangle matrix
void qr(int n, double *matrix) {
    for (int k = 1; k < n; k++) {
        matrix[(k - 1) * n + k - 1] *= matrix[(k - 1) * n + k - 1];
        matrix[(k - 1) * n + k - 1] += matrix[k * n + k - 1] * matrix[k * n + k - 1];
        double angle = matrix[k * n + k - 1] / sqrt(matrix[(k - 1) * n + k - 1]);
        matrix[k * n + k - 1] = (-1) * angle;
    }

    matrix[(k - 1) * n + k - 1] = sqrt(matrix[(k - 1) * n + k - 1]);
}


// From QR to RQ
void set_rq(int n, double *matrix) {
/*!!!!    for (int l = k + 1; l < n; l++) {
        for (int i = k; i < n; i++) {
            double matrix_k = matrix[i * n + k];
            double matrix_i = matrix[i * n + l];
            double cos = sqrt(1 - (matrix[l * n + (k - 1)] * matrix[l * n + (k - 1)]));
            matrix[i * n + k] = matrix_k * cos + matrix_i * matrix[l * n + (k - 1)];
            matrix[i * n + l] = matrix_i * cos - matrix_k * matrix[l * n + (k - 1)];
        }
*/ // Need FIX!!!!!
}
















#endif // QR_FUNCTIONS_H
