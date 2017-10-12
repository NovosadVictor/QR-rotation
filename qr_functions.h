#ifndef QR_FUNCTIONS_H
#define QR_FUNCTIONS_H

#include <math.h>


// Set vector of sins
void set_vector(int n, int k, double *matrix, double *d1, double *d2) {
    matrix[k * n + k - 1] *= matrix[k * n + k - 1];

    for (int i = k + 1; i < n; i++) {
        double x = matrix[k * n + k - 1];
        matrix[k * n + k - 1] += matrix[i * n + k - 1] * matrix[i * n + k - 1];
        if (matrix[k * n + k - 1] < exp(-15)) {
            d1[i - k - 1] = 0.0;
            d2[i - k - 1] = 1.0;
            continue;
        }
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
//                printf("%lf , %lf, cos = %lf, sin = %lf\n",
//                       matrix[k * n + j], matrix[i * n + j], d2[i - k - 1], d1[i - k - 1]);
                matrix[k * n + j] = matrix_k * d2[i - k - 1] - matrix_i * d1[i - k - 1];
                matrix[i * n + j] = matrix_k * d1[i - k - 1] + matrix_i * d2[i - k - 1];
//                printf("matrix k, i %lf , %lf\n", matrix[k * n + j], matrix[i * n + j]);
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
            for (int i = 0; i < n; i++) {
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
void qr(int n, double *matrix, double *d1, double *d2, double *s_k) {
/*    printf("before QR matrix\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%lf ", matrix[i * n + j]);
        printf("\n");
    }
    printf("\n");
*/

    s_k[0] = matrix[(n - 1) * n + (n - 1)];

    if (fabs(s_k[0]) < exp(-20))
        s_k[0] = 1.0;

    printf("s_k = %lf\n", s_k[0]);

    for (int i = 0; i < n; i++)
        matrix[i * n + i] -= s_k[0];

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

        for (int j = k; j < n; j++) {
            double matrix_k = matrix[(k - 1) * n + j];
            double matrix_i = matrix[k * n + j];
            matrix[(k - 1) * n + j] = matrix_k * d2[k - 1] - matrix_i * d1[k - 1];
            matrix[k * n + j] = matrix_k * d1[k - 1] + matrix_i * d2[k - 1];
        }

    }
}


// From QR to RQ
void set_rq(int n, double *matrix, double *d1, double *d2, double *s_k) {
    printf("QR matrix\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%lf ", matrix[i * n + j]);
        printf("\n");
    }
    printf("\n");

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

    for (int i = 0; i < n; i++)
        matrix[i * n + i] += s_k[0];
}


// Finding eigenvalues
void eigenvalues(int n, double *matrix, double *d1, double *d2, double *s_k) {
/*    double sum = 0;
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++)
            sum += matrix[i * n + j] * matrix[i * n + j];
            */
//    int actual_size = n;
    int s = 0;
    while(n > 2) {
        s++;
        to_almost_triangle(n, matrix, d1, d2);
        qr(n, matrix, d1, d2, s_k);
        set_rq(n, matrix, d1, d2, s_k);

        printf("RQ iteration %d\n", s);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                printf("%lf ", matrix[i * n + j]);
            printf("\n");
        }
        printf("\n");

        if (fabs(matrix[(n - 1) * n + (n - 2)]) < exp(-20)) {
            d1[n - 2] = matrix[(n - 1) * n + (n - 1)];
            double *new_matr = (double *) malloc((n - 1) * (n - 1) * sizeof(double));
            for (int i = 0; i < n - 1; i++)
                for (int j = 0; j < n - 1; j++)
                    new_matr[i * (n - 1) + j] = matrix[i * n + j];
            free(matrix);
            matrix = (double *) malloc((n - 1) * (n - 1) * sizeof(double));
            for (int i = 0; i < n - 1; i++)
                for (int j = 0; j < n - 1; j++)
                    matrix[i * (n - 1) + j] = new_matr[i * (n - 1) + j];
            free(new_matr);
            n--;
            continue;
        }

 /*       sum = 0;
        for (int i = 1; i < n; i++)
            for (int j = 0; j < i; j++)
                sum += matrix[i * n + j] * matrix[i * n + j];
                */
    }

    double trace = matrix[0] + matrix[3];
    double det = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    double D = trace * trace - 4 * det;
//    printf("tr = %lf, det = %lf, D = %lf\n", trace, det, D);

    d1[0] = (trace + sqrt(D)) / 2;
    d2[0] = (trace - sqrt(D)) / 2;

//    for (int i = 0; i < n - 1; i++)
//        d1[i] = matrix[i * n + i];
//    d2[0] = matrix[(n - 1) * n + (n - 1)];
    free(matrix);
}


#endif // QR_FUNCTIONS_H
