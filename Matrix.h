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

typedef float data_Type;

typedef struct Matrix
{
    size_t row;
    size_t col;
    data_Type *data;
} Matrix;

Matrix *CreateMatrix(size_t, size_t, data_Type *); //创建一个_row行_col列的矩阵，并更新其中元素，如果_data为NULL则将矩阵中的元素全部定为0

void deleteMatrix(Matrix *); //删除矩阵,释放空间

Matrix *readFile(FILE *); //从文件读入矩阵

Matrix *matmul_plain(Matrix *, Matrix *); //使用三重循环计算矩阵乘法

void plus_for_Strassen(Matrix *, size_t, Matrix *, size_t, Matrix *, size_t, size_t); // Strassen算法中的加法

void minus_for_Strassen(Matrix *, size_t, Matrix *, size_t, Matrix *, size_t, size_t); // Strassen算法中的减法

Matrix *Strassen(Matrix *, size_t, Matrix *, size_t, size_t); // Strassen算法，返回A*B的结果

Matrix *matmul_improved(Matrix *, Matrix *); //优化矩阵乘法

void printMatrix(Matrix *); //输出矩阵

double similar(Matrix *, Matrix *); //比较两个矩阵是否相似
