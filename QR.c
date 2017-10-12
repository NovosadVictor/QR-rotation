#include <malloc.h>
#include <time.h>

#include "qr_functions.h"
#include "formulas.h"



int main(int argc, char **argv) {
    if (argc != 3) {
        printf("You should choose mode(1 - formula, 2 - file)\n");
        printf("And write matrix size(for mode 1) of filename(for mode 2)\n");
        return -1;
    }

    if (atoi(argv[1]) == 2) {
        FILE *fi = fopen(argv[2], "r");
        if (fi == NULL) {
            printf("incorrect filename\n");
            return -1;
        }

        int n;
        if (fscanf(fi, "%d", &n) != 1) {
            printf("Error in matrix size\n");
            return -1;
        }

        double *matrix = (double *) malloc(n * n * sizeof(double));
        double *d1 = (double *) malloc((n - 1) * sizeof(double));
        double *d2 = (double *) malloc((n - 1) * sizeof(double));
        double *s_k = (double *) malloc(sizeof(double));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double elem;
                if (fscanf(fi, "%lf", &elem) != 1) {
                    printf("Incorrect matrix in file\n");
                    return -1;
                }
                matrix[i * n + j] = elem;
            }
        char c;
        if (fscanf(fi, "%c", &c) > 0) {
            printf("NOT EOF\n");
            return -1;
        }

        printf("Input matrix\n");
        print(matrix, n);

        time_t start = clock();

        eigenvalues(n, matrix, d1, d2, s_k);

        time_t end = clock();

        printf("Result:\n\tEigenvalues:\n");
        for (int i = 0; i < n - 1; i++)
            printf("%lf, ", d1[i]);
        printf("%lf\n", d2[0]);

        printf("execution time: %f sec\n", (double) (end - start) / CLOCKS_PER_SEC);
 //       free(matrix);
        free(d1);
        free(d2);
        free(s_k);
        fclose(fi);
    }
    if (atoi(argv[1]) == 1) {
        int n = atoi(argv[2]);

        double *matrix = (double *) malloc(n * n * sizeof(double));
        double *d1 = (double *) malloc((n - 1) * sizeof(double));
        double *d2 = (double *) malloc((n - 1) * sizeof(double));
        double *s_k = (double *) malloc(sizeof(double));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                matrix[i * n + j] = formula(i, j, n);

        printf("Input matrix\n");
        print(matrix, n);

        time_t start = clock();

        eigenvalues(n, matrix, d1, d2, s_k);

        time_t end = clock();

        printf("Result:\n\tEigenvalues:\n");
        for (int i = 0; i < n - 1; i++)
            printf("%lf, ", d1[i]);
        printf("%lf\n", d2[0]);

        printf("execution time: %f sec\n", (double) (end - start) / CLOCKS_PER_SEC);
//        free(matrix);
        free(d1);
        free(d2);
        free(s_k);
    }
    return 0;
}