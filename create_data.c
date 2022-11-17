#include <stdio.h>
#include <stdlib.h>
#include <string.h>
double rand_num()
{
    return ((double)rand()/RAND_MAX)*10;
}
int main()
{
    srand(20030408);
    char *file_name_A = malloc(sizeof(char) * 20);
    char *file_name_B = malloc(sizeof(char) * 20);
    size_t size;
    FILE *file_A;
    FILE *file_B;
    //创造两个16*16的矩阵
    size=(1<<4);
    strcpy(file_name_A, "Matrix_A_16.txt");
    strcpy(file_name_B, "Matrix_B_16.txt");
    file_A=fopen(file_name_A,"w");
    file_B=fopen(file_name_B,"w");
    fprintf(file_A,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_A,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_A,"%.1lf ",rand_num());
    fclose(file_A);
    fprintf(file_B,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_B,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_B,"%.1lf ",rand_num());
    fclose(file_B);
    //创造两个128*128的矩阵
    size=(1<<7);
    strcpy(file_name_A, "Matrix_A_128.txt");
    strcpy(file_name_B, "Matrix_B_128.txt");
    file_A=fopen(file_name_A,"w");
    file_B=fopen(file_name_B,"w");
    fprintf(file_A,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_A,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_A,"%.1lf ",rand_num());
    fclose(file_A);
    fprintf(file_B,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_B,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_B,"%.1lf ",rand_num());
    fclose(file_B);
    //创造两个1K*1K的矩阵
    size=(1<<10);
    strcpy(file_name_A, "Matrix_A_1K.txt");
    strcpy(file_name_B, "Matrix_B_1K.txt");
    file_A=fopen(file_name_A,"w");
    file_B=fopen(file_name_B,"w");
    fprintf(file_A,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_A,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_A,"%.1lf ",rand_num());
    fclose(file_A);
    fprintf(file_B,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_B,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_B,"%.1lf ",rand_num());
    fclose(file_B);
    //创造两个8K*8K的矩阵
    size=(1<<13);
    strcpy(file_name_A, "Matrix_A_8K.txt");
    strcpy(file_name_B, "Matrix_B_8K.txt");
    file_A=fopen(file_name_A,"w");
    file_B=fopen(file_name_B,"w");
    fprintf(file_A,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_A,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_A,"%.1lf ",rand_num());
    fclose(file_A);
    fprintf(file_B,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_B,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_B,"%.1lf ",rand_num());
    fclose(file_B);
    //创造两个64K*64K的矩阵
    size=(1<<16);
    strcpy(file_name_A, "Matrix_A_64K.txt");
    strcpy(file_name_B, "Matrix_B_64K.txt");
    file_A=fopen(file_name_A,"w");
    file_B=fopen(file_name_B,"w");
    fprintf(file_A,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_A,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_A,"%.1lf ",rand_num());
    fclose(file_A);
    fprintf(file_B,"%lu %lu\n",size,size);
    for(int i=1;i<=size;i++,fprintf(file_B,"\n"))
        for(int j=1;j<=size;j++)
            fprintf(file_B,"%.1lf ",rand_num());
    fclose(file_B);
    return 0;
}