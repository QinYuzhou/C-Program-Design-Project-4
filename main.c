#include <sys/time.h>
#include <cblas.h>
#include "Matrix.h"

#define TIME_START gettimeofday(&start, NULL);
#define TIME_END(NAME)                                                               \
    gettimeofday(&end, NULL);                                                        \
    duration = (end.tv_sec - start.tv_sec) + (1e-6) * (end.tv_usec - start.tv_usec); \
    printf("%s's duration = %lfs\n", (NAME), duration);

inline double rand_num()
{
    return ((double)rand() / RAND_MAX) * 10;
}

pMatrix create_random_Matrix(size_t size)
{
    pMatrix ret = CreateMatrix(size, size, NULL);
    size *= size;
    for (size_t i = 0; i < size; i++)
        ret->data[i] = rand_num();
    return ret;
}

int main(int args, char *argv[])
{
    struct timeval start, end;
    double duration, cnt1 = 0, cnt2 = 0, cnt3 = 0;
    int size = 16;
    pMatrix A, B, C, D, E;
    printf("Create A\n");
    A = create_random_Matrix(1ll << size);
    printf("Create B\n");
    B = create_random_Matrix(1ll << size);
    printf("Create Over\n");

    C = matmul_plain(A, B);
    TIME_START
    C = matmul_plain(A, B);
    TIME_END("plain")
    cnt1 += duration;
    TIME_START
    C = matmul_plain(A, B);
    TIME_END("plain")
    cnt1 += duration;
    TIME_START
    C = matmul_plain(A, B);
    TIME_END("plain")
    cnt1 += duration;
    TIME_START
    C = matmul_plain(A, B);
    TIME_END("plain")
    cnt1 += duration;
    TIME_START
    C = matmul_plain(A, B);
    TIME_END("plain")
    cnt1 += duration;

    D = matmul_improved(A, B);
    TIME_START
    D = matmul_improved(A, B);
    TIME_END("improved")
    cnt2 += duration;
    TIME_START
    D = matmul_improved(A, B);
    TIME_END("improved")
    cnt2 += duration;
    TIME_START
    D = matmul_improved(A, B);
    TIME_END("improved")
    cnt2 += duration;
    TIME_START
    D = matmul_improved(A, B);
    TIME_END("improved")
    cnt2 += duration;
    TIME_START
    D = matmul_improved(A, B);
    TIME_END("improved")
    cnt2 += duration;

    E = CreateMatrix(A->row, B->col, NULL);
    if (sizeof(data_Type) == 8)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    if (sizeof(data_Type) == 4)
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    TIME_START
    E = CreateMatrix(A->row, B->col, NULL);
    if (sizeof(data_Type) == 8)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    if (sizeof(data_Type) == 4)
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    TIME_END("OpenBLAS")
    cnt3 += duration;
    TIME_START
    E = CreateMatrix(A->row, B->col, NULL);
    if (sizeof(data_Type) == 8)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    if (sizeof(data_Type) == 4)
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    TIME_END("OpenBLAS")
    cnt3 += duration;
    TIME_START
    E = CreateMatrix(A->row, B->col, NULL);
    if (sizeof(data_Type) == 8)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    if (sizeof(data_Type) == 4)
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    TIME_END("OpenBLAS")
    cnt3 += duration;
    TIME_START
    E = CreateMatrix(A->row, B->col, NULL);
    if (sizeof(data_Type) == 8)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    if (sizeof(data_Type) == 4)
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    TIME_END("OpenBLAS")
    cnt3 += duration;
    TIME_START
    E = CreateMatrix(A->row, B->col, NULL);
    if (sizeof(data_Type) == 8)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    if (sizeof(data_Type) == 4)
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->row, B->col, A->col, 1, A->data, A->col, B->data, B->col, 0, E->data, E->col);
    TIME_END("OpenBLAS")
    cnt3 += duration;

    printf("%.6lf\n", similar(C, D));
    printf("%.6lf\n", similar(C, E));
    printf("%.6lf\n", similar(D, E));
    printf("%.7lfs\n", cnt1 / 5);
    printf("%.7lfs\n", cnt2 / 5);
    printf("%.7lfs\n", cnt3 / 5);

    return 0;
}
