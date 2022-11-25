#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <math.h>
#include <omp.h>

#define _OMP_THREAD_ 4

#ifndef deleteMat
#define deleteMat(Mat) (deleteMatrix(Mat), Mat = NULL)
#endif

typedef double data_Type;

typedef struct Matrix
{
    size_t row;
    size_t col;
    data_Type *data;
} Matrix;
typedef Matrix *pMatrix;

pMatrix CreateMatrix(size_t, size_t, data_Type *); //创建一个_row行_col列的矩阵，并更新其中元素，如果_data为NULL则将矩阵中的元素全部定为0

void deleteMatrix(pMatrix); //删除矩阵,释放空间

pMatrix readFile(FILE *); //从文件读入矩阵

pMatrix matmul_plain(pMatrix, pMatrix); //使用三重循环计算矩阵乘法

void plus_for_Strassen(pMatrix, size_t, pMatrix, size_t, pMatrix, size_t, size_t); // Strassen算法中的加法

void minus_for_Strassen(pMatrix, size_t, pMatrix, size_t, pMatrix, size_t, size_t); // Strassen算法中的减法

pMatrix Strassen(pMatrix, size_t, pMatrix, size_t, size_t); // Strassen算法，返回A*B的结果

pMatrix matmul_improved(pMatrix, pMatrix); //优化矩阵乘法

void printMatrix(pMatrix); //输出矩阵

double similar(pMatrix, pMatrix); //比较两个矩阵是否相似
