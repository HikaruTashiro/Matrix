#include "Matrix.c"

int main()
{
    Matrix M;
    M = createIdentityMatrix(5);
    puts("Matrix M:\n");
    printMatrix(M);
    multiplyScalarMatrix(3.14,M);
    puts("Matrix pi * M:\n");
    printMatrix(M);

    Matrix A;
    A = createIdentityMatrix(5);
    puts("Matrix A:\n");
    printMatrix(A);
    A = multiplyScalarMatrix(-1.0, A);
    puts("Matrix -A:\n");
    printMatrix(A);

    Matrix B = sumMatrices(A,M);
    puts("Matrix pi * M - A:\n");
    printMatrix(B);

    B = changeInCoordinate(B,5,4,5.4);
    B = changeInCoordinate(B,3,2,9.8);
    B = changeInCoordinate(B,1,3,10.1);
    B = changeInCoordinate(B,1,2,11.1);
    B = changeInCoordinate(B,2,5,2.3);
    puts("Matrix B after changes:\n");
    printMatrix(B);
    puts("Transpose of Matrix B:\n");
    B = transposeMatrix(B);
    printMatrix(B);
    puts("Interchange Collumns 2 and 5 of Matrix B:\n");
    B = interchangeCollumn(B,1,4);
    printMatrix(B);

    Matrix Aux;
    puts("Matrix M :\n");
    printMatrix(M);
    Aux = M;
    puts("Minor (2,1) of Matrix M :\n");
    M = minorOfMatrix(M,2,1);
    printMatrix(M);
    emptyMatrix(Aux);

    Matrix C;
    C = createIdentityMatrix(3);
    C = changeInCoordinate(C,1,1,1.0);
    C = changeInCoordinate(C,1,2,0.0);
    C = changeInCoordinate(C,1,3,0.0);
    //C = changeInCoordinate(C,1,4,4.0);
    C = changeInCoordinate(C,2,1,3.0);
    C = changeInCoordinate(C,2,2,1.0);
    C = changeInCoordinate(C,2,3,0.0);
    //C = changeInCoordinate(C,2,4,0.0);
    C = changeInCoordinate(C,3,1,5.0);
    C = changeInCoordinate(C,3,2,7.0);
    C = changeInCoordinate(C,3,3,1.0);
    //C = changeInCoordinate(C,3,4,3.0);
    //C = changeInCoordinate(C,4,1,3.0);
    //C = changeInCoordinate(C,4,2,0.0);
    //C = changeInCoordinate(C,4,3,1.0);
    //C = changeInCoordinate(C,4,4,0.0);
    puts("Matrix C:\n");
    printMatrix(C);
    printf("Cofactor C = %lf\n\n",cofactor(C,1,1));
    Matrix C_adj, C_inv;
    
    C_adj = adjugateMatrix(C);
    C_inv = inverseMatrix(C);
    puts("Adjugate of Matrix C:\n");
    printMatrix(C_adj);
    puts("Inverse of Matrix C:\n");  
    printMatrix(C_inv);
    printf("det C = %lf\n",determinant(C));
    puts("C Adjugate * Matrix C:\n");
    Matrix MulC;
    MulC = multiplyMatrices(C,C_adj);
    printMatrix(MulC);

    emptyMatrix(MulC);
    emptyMatrix(C_adj);
    emptyMatrix(C_inv);
    emptyMatrix(C);

    Matrix D;
    D = createIdentityMatrix(4);
    D = changeInCoordinate(D,1,1,1.0);
    D = changeInCoordinate(D,1,2,2.0);
    D = changeInCoordinate(D,1,3,4.0);
    D = changeInCoordinate(D,1,4,8.0);
    D = changeInCoordinate(D,2,1,1.0);
    D = changeInCoordinate(D,2,2,3.0);
    D = changeInCoordinate(D,2,3,9.0);
    D = changeInCoordinate(D,2,4,27.0);
    D = changeInCoordinate(D,3,1,1.0);
    D = changeInCoordinate(D,3,2,5.0);
    D = changeInCoordinate(D,3,3,25.0);
    D = changeInCoordinate(D,3,4,125.0);
    D = changeInCoordinate(D,4,1,1.0);
    D = changeInCoordinate(D,4,2,7.0);
    D = changeInCoordinate(D,4,3,49.0);
    D = changeInCoordinate(D,4,4,343.0);
    printf("Matrix D:\n");
    printMatrix(D);

    printf("det D = %.2lf\n",determinant(D));

    Matrix E;
    E = createIdentityMatrix(3);
    E = changeInCoordinate(E,1,1,-3.0);
    E = changeInCoordinate(E,1,2,6.0);
    E = changeInCoordinate(E,1,3,12.0);
    E = changeInCoordinate(E,2,1,-1.0);
    E = changeInCoordinate(E,2,2,3.0);
    E = changeInCoordinate(E,2,3,5.0);
    E = changeInCoordinate(E,3,1,-1.0);
    E = changeInCoordinate(E,3,2,9.0);
    E = changeInCoordinate(E,3,3,25.0);
    printf("Matrix E:\n");
    printMatrix(E);

    printf("det E = %.2lf\n",determinant(E));

    Matrix F;
    F = createIdentityMatrix(1);
    F = changeInCoordinate(F,1,1,25.0);
    printf("Matrix F:\n");
    printMatrix(F);

    printf("det F = %.2lf\n",determinant(F));
    
    emptyMatrix(A);
    emptyMatrix(B);
    emptyMatrix(D);
    emptyMatrix(E);
    emptyMatrix(F);
    emptyMatrix(M);
}
