#include "Matrix.h"

pMatrix CreateMatrix(size_t _row, size_t _col, data_Type *_data) //创建一个_row行_col列的矩阵，并更新其中元素，如果_data为NULL则将矩阵中的元素全部定为0
{
    pMatrix ret = malloc(sizeof(Matrix));
    ret->row = _row;
    ret->col = _col;
    if (_data == NULL)
        ret->data = aligned_alloc(256, _row * _col * sizeof(data_Type));
    else
        ret->data = _data;
    return ret;
}

void deleteMatrix(pMatrix mat) //删除矩阵,释放空间
{
    if (mat == NULL)
        return;
    if (mat->data != NULL)
        free(mat->data);
    free(mat);
    return;
}

pMatrix readFile(FILE *file) //从文件读入矩阵
{
    pMatrix mat = malloc(sizeof(Matrix));
    fscanf(file, "%zu %zu", &mat->row, &mat->col);
    mat->data = aligned_alloc(256, sizeof(data_Type) * mat->row * mat->col);
    for (int i = 0; i < mat->row * mat->col; i++)
        fscanf(file, "%lf", &mat->data[i]);
    return mat;
}

pMatrix matmul_plain(pMatrix mat_1, pMatrix mat_2) //使用三重循环计算矩阵乘法
{
    if (mat_1 == NULL || mat_2 == NULL)
        return NULL;
    if (mat_1->data == NULL || mat_2->data == NULL)
        return NULL;
    if (mat_1->col != mat_2->row)
        return NULL;
    pMatrix ret = CreateMatrix(mat_1->row, mat_2->col, NULL);
    for (size_t i = 0; i < mat_1->row; i++)
        for (size_t j = 0; j < mat_1->col; j++)
            for (size_t k = 0; k < mat_2->col; k++)
                ret->data[i * ret->col + k] += mat_1->data[i * mat_1->col + j] * mat_2->data[j * mat_2->col + k];
    return ret;
}

void plus_for_Strassen(pMatrix mat_1, size_t startIndex_1, pMatrix mat_2, size_t startIndex_2, pMatrix mat_3, size_t startIndex_3, size_t size) // Strassen算法中的加法
{
#pragma omp parallel for
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            mat_3->data[startIndex_3 + i * mat_3->col + j] = mat_1->data[startIndex_1 + i * mat_1->col + j] + mat_2->data[startIndex_2 + i * mat_2->col + j];
}

void minus_for_Strassen(pMatrix mat_1, size_t startIndex_1, pMatrix mat_2, size_t startIndex_2, pMatrix mat_3, size_t startIndex_3, size_t size) // Strassen算法中的减法
{
#pragma omp parallel for
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            mat_3->data[startIndex_3 + i * mat_3->col + j] = mat_1->data[startIndex_1 + i * mat_1->col + j] - mat_2->data[startIndex_2 + i * mat_2->col + j];
}

pMatrix Strassen(pMatrix A, size_t index_A, pMatrix B, size_t index_B, size_t size) // Strassen算法，返回A*B的结果
{
    if (size == 64)
    {
        pMatrix ret = CreateMatrix(size, size, NULL);
        memset(ret->data, 0, ret->col * ret->row * sizeof(data_Type));
#pragma omp parallel for
        for (size_t i = 0; i < size; i++)
            for (size_t j = 0; j < size; j++)
                for (size_t k = 0; k < size; k++)
                    ret->data[i * ret->col + k] += A->data[index_A + i * A->col + j] * B->data[index_B + j * B->col + k];
        return ret;
    }
    size_t new_size = (size >> 1);
    size_t a11 = index_A, a12 = a11 + new_size, a21 = a11 + new_size * A->col, a22 = a21 + new_size;
    size_t b11 = index_B, b12 = b11 + new_size, b21 = b11 + new_size * B->col, b22 = b21 + new_size;

    pMatrix S1 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S2 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S3 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S4 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S5 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S6 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S7 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S8 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S9 = CreateMatrix(new_size, new_size, NULL);
    pMatrix S10 = CreateMatrix(new_size, new_size, NULL);

    minus_for_Strassen(B, b12, B, b22, S1, 0, new_size);
    plus_for_Strassen(A, a11, A, a12, S2, 0, new_size);
    plus_for_Strassen(A, a21, A, a22, S3, 0, new_size);
    minus_for_Strassen(B, b21, B, b11, S4, 0, new_size);
    plus_for_Strassen(A, a11, A, a22, S5, 0, new_size);
    plus_for_Strassen(B, b11, B, b22, S6, 0, new_size);
    minus_for_Strassen(A, a12, A, a22, S7, 0, new_size);
    plus_for_Strassen(B, b21, B, b22, S8, 0, new_size);
    minus_for_Strassen(A, a11, A, a21, S9, 0, new_size);
    plus_for_Strassen(B, b11, B, b12, S10, 0, new_size);

    pMatrix P1 = Strassen(A, a11, S1, 0, new_size);
    deleteMat(S1);
    pMatrix P2 = Strassen(S2, 0, B, b22, new_size);
    deleteMat(S2);
    pMatrix P3 = Strassen(S3, 0, B, b11, new_size);
    deleteMat(S3);
    pMatrix P4 = Strassen(A, a22, S4, 0, new_size);
    deleteMat(S4);
    pMatrix P5 = Strassen(S5, 0, S6, 0, new_size);
    deleteMat(S5);
    deleteMat(S6);
    pMatrix P6 = Strassen(S7, 0, S8, 0, new_size);
    deleteMat(S7);
    deleteMat(S8);
    pMatrix P7 = Strassen(S9, 0, S10, 0, new_size);
    deleteMat(S9);
    deleteMat(S10);

    pMatrix C = CreateMatrix(size, size, NULL);
    size_t c11 = 0, c12 = c11 + new_size, c21 = c11 + new_size * size, c22 = c21 + new_size;
    plus_for_Strassen(P5, 0, P4, 0, C, c11, new_size);
    minus_for_Strassen(C, c11, P2, 0, C, c11, new_size);
    plus_for_Strassen(C, c11, P6, 0, C, c11, new_size);
    deleteMat(P6);
    plus_for_Strassen(P1, 0, P2, 0, C, c12, new_size);
    deleteMat(P2);
    plus_for_Strassen(P3, 0, P4, 0, C, c21, new_size);
    deleteMat(P4);
    plus_for_Strassen(P5, 0, P1, 0, C, c22, new_size);
    deleteMat(P5);
    deleteMat(P1);
    minus_for_Strassen(C, c22, P3, 0, C, c22, new_size);
    deleteMat(P3);
    minus_for_Strassen(C, c22, P7, 0, C, c22, new_size);
    deleteMat(P7);

    return C;
}

pMatrix matmul_improved(pMatrix mat_1, pMatrix mat_2) //优化矩阵乘法
{
    if (mat_1 == NULL || mat_2 == NULL)
        return NULL;
    if (mat_1->data == NULL || mat_2->data == NULL)
        return NULL;
    if (mat_1->col != mat_2->row)
        return NULL;
    if (mat_1->col != mat_1->row || mat_2->col != mat_2->row)
        return NULL;
    return Strassen(mat_1, 0, mat_2, 0, mat_1->col);
}

void printMatrix(pMatrix mat) //输出矩阵
{
    if (mat == NULL)
        return;
    if (mat->data == NULL)
        return;
    for (int i = 0; i < mat->row * mat->col; i++)
    {
        if (i % mat->col == 0)
            printf("\n");
        printf("%.1f ", mat->data[i]);
    }
    printf("\n");
}

double similar(pMatrix mat_1, pMatrix mat_2) //比较两个矩阵是否相似
{
    if (mat_1 == NULL || mat_2 == NULL)
        return -1;
    if (mat_1->col != mat_2->col || mat_1->row != mat_2->row)
        return -1;
    double differ = 0;
    for (int i = 0; i < mat_1->col * mat_2->row; i++)
        differ = fmax(differ, fabs(mat_1->data[i] - mat_2->data[i]));
    return differ;
}