#include <sys/time.h>
#include </usr/local/opt/openblas/include/cblas.h>
#include "Matrix.h"

#define TIME_START gettimeofday(&start, NULL);
#define TIME_END                                                                     \
    gettimeofday(&end, NULL);                                                        \
    duration = (end.tv_sec - start.tv_sec) + (1e-6) * (end.tv_usec - start.tv_usec); \
    printf("duration = %lfs\n", duration);

int main(int args, char *argv[])
{
    struct timeval start, end;
    double duration;
    pMatrix A, B, C, D, E;
    FILE *file;
    file = fopen("Matrix_A_8K.txt", "r");
    if (file == NULL)
    {
        printf("File does not exist!\n");
        return 1;
    }
    A = readFile(file);
    fclose(file);
    // printMatrix(A);
    file = fopen("Matrix_B_8K.txt", "r");
    if (file == NULL)
    {
        printf("File does not exist!\n");
        return 1;
    }
    B = readFile(file);
    fclose(file);
    // printMatrix(B);
    // E = CreateMatrix(A->row, B->col, NULL);
    // C = matmul_plain(A, B);
    // D = matmul_improved(A, B);
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, D->data, E->col);

    // TIME_START
    // C = matmul_plain(A, B);
    // TIME_END

    TIME_START
    D = matmul_improved(A, B);
    TIME_END

    E = CreateMatrix(A->row, B->col, NULL);
    TIME_START
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    TIME_END
    printf("%lf\n", similar(C, D));
    printf("%lf\n", similar(C, E));
    printf("%lf\n", similar(D, E));
    return 0;
}
