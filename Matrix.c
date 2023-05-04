#include "Matrix.h"
#include <stdio.h>
#include <stdlib.h>

/*Implementation*/

Matrix createMatrix(const int a,const int b)    //tested
{
    Matrix M = NULL; 

    if(a > 0 && b > 0)
    {
        M = (Matrix) malloc(sizeof(MatrixDef));
        if(M != NULL)
        {
            int i, size;    
            size = a * b;
            M->rows = a;
            M->columns = b;
            M->matrix = (double*) malloc(size * sizeof(double));

            for(i = 0; i < size; i++)
                M->matrix[i] = 0;
        }
        else
            puts("Failed to allocate Memory");
    }
    else
        puts("ERROR: negative number passed as argument\n");

    return M;
}

Matrix createIdentityMatrix(const int n)   //tested
{
    Matrix M = createMatrix(n,n);
    if(M != NULL)
    {
        int i;
        for ( i = 0; i < n; i++)
            M->matrix[i * (n + 1)] = 1.0;
    }
        
    return M;
}

Matrix sumMatrices(Matrix M1,Matrix M2)  //tested
{
    Matrix M = NULL;
    
    if(M1 != NULL && M2 != NULL)
    {
        if(M1->columns == M2->columns && M1->rows == M2->rows)
        {
            M = createMatrix(M1->rows, M2->columns);
            if(M != NULL)
            {
                int i, size;
                size = M1->columns * M1->rows;

                for(i = 0; i < size; i++)
                    M->matrix[i] = M1->matrix[i] + M2->matrix[i];
            }
        }
        else
            puts("ERROR: Matrices don't have the same size");
    }
    else
        puts("ERROR: Empty Matrix passed as argument");
    
    return M;
}

Matrix multiplyScalarMatrix(const double lambda,Matrix M)   //tested
{
    if(M != NULL)
    {
        int i, size;
        size = M->columns * M->rows;

        for(i = 0; i < size; i++)
            M->matrix[i] *= lambda;
    }
    else
        printf("ERROR: Empty Matrix Passed as Argument\n\n");
    
    return M;
}

Matrix multiplyMatrices(Matrix A,Matrix B)
{
    Matrix C = NULL;
    if(A != NULL && B != NULL && A->columns == B->rows)
    {
        C = createMatrix(A->rows,B->columns);
        if(C != NULL)
        {
            int i, j, k;
            int column, row;
            row = C->rows;
            column = C->columns;
            for ( i = 0; i < row; i++)
                for ( j = 0; j < column; j++)
                    for( k = 0; k < column; k++)
                        C->matrix[i * column + j] += A->matrix[i * column + k] * B->matrix[k * column + j];
        }
    }

    return C;
}

Matrix transposeMatrix(Matrix M)    //tested
{
    if(M != NULL)
    {
        double aux;
        int i, j;

        for (i = 0; i < M->rows; i++)
            for(j = i; j < M->columns; j++)
            {
                aux = M->matrix[i * M->columns + j];
                M->matrix[i * M->columns + j] = M->matrix[j * M->rows + i];
                M->matrix[j * M->rows + i] = aux;
            }
        i = M->columns;
        M->columns = M->rows;
        M->rows = i; 
    }
    else
        printf("ERROR: Empty Matrix Passed as Argument\n\n");

    return M;
}

Matrix inverseMatrix(Matrix M) //find a better way
{
    Matrix M_inv = NULL, M_copy;
    if(M != NULL && M->rows == M->columns)     
    {
        M_copy = copyMatrix(M);
        M_inv = createIdentityMatrix(M->rows);
        double aux;

        /*Actual trash O(n^4), find another way*/
        for (int i = 0; i < M->rows; i++)
        {
            if(M_copy->matrix[i * (M->rows + 1)] == 0.0)
            {
                int k = i + 1;
                while (k < M->rows && M_copy->matrix[k * M->rows + i] == 0.0) 
                    k++;
                if(k < M->rows)
                    M_copy = interchangeRow(M_copy, i + 1, k + 1);
                else if(i + 1 < M->rows)
                    i++;
            }

            aux = M_copy->matrix[i * (M->rows + 1)];
            for (int j = i; j < M->rows; j++)
            {
                M_copy->matrix[i * M->rows + j] /= aux;
                M_inv->matrix[i * M->rows + j] /= aux;

                for (int k = i + 1; k < M->rows; k++) 
                {
                    for (int l = 0; l < M->rows; l++) 
                        M_inv->matrix[k * M->rows + l] -= M_inv->matrix[i * M->rows + l] * M_copy->matrix[k * M->rows + i];
                    M_copy->matrix[k * M->rows + j] -= M_copy->matrix[i * M->rows + j] * M_copy->matrix[k * M->rows + i];
                }
            }
        }
        emptyMatrix(M_copy);
    }

    return M_inv;
}

Matrix adjugateMatrix(Matrix M)
{
    int i, j, n;
    Matrix M_adj = NULL;
    if(M != NULL && M->rows == M->columns)
    {
        n = M->rows;
        M_adj = createMatrix(n,n);
        if(M_adj != NULL)
        {
            for ( i = 0; i < n; i++)
                for( j = 0; j < n; j++)
                    M_adj->matrix[i * n + j] = cofactor(M, i + 1, j + 1);
            M_adj = transposeMatrix(M_adj);
        }
    }

    return M_adj;
}

double cofactor(Matrix M, int a, int b)
{
    int n;
    Matrix M_aux;
    double cofactorM = 0.0;
    if(M != NULL)
    {
        n = M->columns;
        if(n == M->rows && a > 0 && a <= n && b > 0 && b <= n)
        {
            M_aux = minorOfMatrix(M,a,b);
            cofactorM = determinant(M_aux);
            if((a + b) % 2 == 1)
                cofactorM *= -1;
            emptyMatrix(M_aux);
        }
    }

    return cofactorM;
}

double determinant(Matrix M)
{
    double det = 0.0;
    if(M != NULL && M->rows == M->columns)
    {
            Matrix M_aux = copyMatrix(M);
            double max; 
            int M_size = M->rows;
            int i;
            det = 1.0;

            for (i = 0; i < M_size; i++) 
            {
                int search = i;
                if(M_aux->matrix[i * (M_size + 1)] == 0.0)
                {
                    while (M_aux->matrix[i * M_size + search] == 0.0 && search < M_size) 
                        search++;
                    if(search < M_size)
                        M_aux = interchangeCollumn(M_aux, search + 1, i + 1);
                    else
                        return 0.0;
                }

                max = M_aux->matrix[i * (M_size + 1)];
                det *= max;
                for (int j = M_size - 1; j >= i; j--) 
                {
                    M_aux->matrix[j * M_size + i] /= max; 
                    for (int k = i + 1; k < M_size; k++)
                        M_aux->matrix[j * M_size + k] -= M_aux->matrix[j * M_size + i] * M_aux->matrix[i * M_size + k]; 
                }
            }

            emptyMatrix(M_aux);
    }
    else {
        puts("ERROR, either a NULL pointer was passed or the Matrix is not Squared");
    }

    return det;
}

Matrix copyMatrix(Matrix M) //tested
{
    Matrix M_copy = NULL;
    if(M != NULL)
    {
        M_copy = createMatrix(M->rows, M->columns);
        if(M_copy != NULL)
        {
            int i, size;
            size = M->rows * M->columns;
            M_copy->rows = M->rows;
            M_copy->columns = M->columns; 

            for (i = 0; i < size; i++)
                M_copy->matrix[i] = M->matrix[i];
        }
    }
    
    return M_copy;
}

Matrix minorOfMatrix(Matrix M, int a, int b)    //tested
{
    a--,b--;
    Matrix M_minor;
    
    if(M != NULL && M->rows == M->columns && M->rows > 1)
    {
        M_minor = createMatrix(M->rows - 1, M->rows - 1);
        if(M_minor != NULL)
        {
            int i, j, k, size;

            size = M->rows * M->rows;
            k = 0;
            for(i = 0;  i < size; i++)
                if(i % M->columns != b && !( (j = i - a*M->columns) >= 0 && j < M->columns ))
                    M_minor->matrix[k] = M->matrix[i], k++;
        }
        else {
            puts("ERROR: Could not allocate Memory");
        }
        /* 
            i = a * column + j, if i is a coordinate in the row a, also note that j is in [0,column[
            
            i = q * column + b, if i is a coordinate in the collumn b, note that i ≡ b (mod column)

            The if above is a negation of the conditions necessary for either one of the conditions
        */
    }

    return M_minor;
}

Matrix interchangeCollumn(Matrix M, int a, int b)   //tested
{
    a--,b--;
    int i;
    int column = M->columns;
    double aux;
    if(M != NULL)
        for(i = 0; i < M->rows; i++)
        {
            aux = M->matrix[i*column + a];
            M->matrix[i*column + a] = M->matrix[i*column + b];
            M->matrix[i*column + b] = aux;
        }

    return M;
}

Matrix interchangeRow(Matrix M,int a,int b)
{
    a--,b--;
    int i;
    int column = M->columns;
    double aux;
    if(M != NULL)
        for(i = 0; i < column; i++)
        {
            aux = M->matrix[a * column + i];
            M->matrix[a * column + i] = M->matrix[b * column + i];
            M->matrix[b * column + i] = aux;
        }
    
    return M;
}

void printMatrix(Matrix M)  //tested
{
    int i, size;

    if(M != NULL)
    {
        size = M->columns * M->rows;

        for(i = 0; i < size; i++)
        {
            if(i % M->columns == 0)
                printf("\n");
            printf("%7.2lf  ",M->matrix[i]);
        }
    }
    else 
        printf("Empty Matrix");
    printf("\n\n");
}

Matrix changeInCoordinate(Matrix M,int x,int y,double b) //tested
{
    x--, y--;
    if(M != NULL)
    {
        if(x >= 0 && y >= 0 && x < M->rows && y < M->columns)
            M->matrix[x * M->columns + y] = b;
    }
    else
        printf("ERROR: Empty Matrix Passed as Argument\n\n");
    
    return M;
}

boolean equalMatrix(Matrix A,Matrix B)
{
    int i, size;
    boolean equal = FALSE;
    if(A != NULL && B != NULL)
    {
        if(A->rows == B->rows && A->columns == B->columns)
        {
            i = 0;
            size = A->rows;
            while (i < size && A->matrix[i] == B->matrix[i])
            {i++;}
            if(i == size)
                equal = TRUE;
        }
    }
    else if(A == B)
        equal = TRUE;

    return equal;
}

boolean symmetricMatrix(Matrix M)
{
    boolean symmetric;
    Matrix M_trsp;
    M_trsp = copyMatrix(M);
    M_trsp = transposeMatrix(M_trsp);
    symmetric = equalMatrix(M_trsp,M);
    emptyMatrix(M_trsp);
    return symmetric;
}

void emptyMatrix(Matrix M)  //tested
{
    if(M != NULL)
    {
        free(M->matrix);
        free(M);
    }   
}

Vector createVector(int dimension)
{
    Vector V = createMatrix(1,dimension);
    return V;
}

double dotProduct(Vector V, Vector U)
{
    int i;
    double dot = 0.0;
    if(V != NULL && U != NULL && V->rows == 1 && U->rows == 1 && V->columns == U->columns)
    {   
        for(i = 0; i < V->columns; i++)
            dot += V->matrix[i] * U->matrix[i];
    }

    return dot;
}

Vector crossProduct(Vector V, Vector U)
{
    double det;
    Vector W;
    if(V != NULL && U != NULL && V->rows == 1 && U->rows == 1 && V->columns == 3 && U->columns == 3)
    {
        W = createVector(3);
        det = U->matrix[1] * V->matrix[2] - U->matrix[2] * V->matrix[1];
        W->matrix[0] = det;
        det = U->matrix[0] * V->matrix[2] - U->matrix[2] * V->matrix[0];
        W->matrix[1] = -det;
        det = U->matrix[0] * V->matrix[1] - U->matrix[1] * V->matrix[0];
        W->matrix[2] = det;
    }
    else
        W = NULL;

    return W;
}

double angleBetween(Vector U,Vector V)
{
    double angle = 0.0;
    double cosine;

    if(U != NULL && V != NULL && U->rows == 1 && V->rows == 1 && V->columns == U->columns)
    {
        cosine = dotProduct(U,V)/(normOf(U) * normOf(V));
        angle = acos(cosine);
    }

    return angle;
}

double normOf(Vector V)
{
    int i;
    double norm = 0.0;
    if(V != NULL && V->rows == 1)
    {
        for(i = 0; i < V->columns; i++)
            norm += V->matrix[i] * V->matrix[i];
        norm = sqrt(norm);
    }

    return norm;
}
