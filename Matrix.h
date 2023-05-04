/*
    File:   Matrix.h
    Name:   Jonas Edward Tashiro
    Date:   01/05/2022
    Description:    Implements the concept of a Matrix
                    Unary operations alter the contents and return the address of the argument passed
                    Binary operations create new Matrices end return their address 
                        // Therefore when doing binary operations create a new Matrix
*/

/*Directives*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Type definitions*/
typedef struct{
    int rows;
    int columns;
    double *matrix;
} MatrixDef;

typedef MatrixDef* Matrix;
typedef enum{FALSE,TRUE} boolean;
typedef Matrix Vector;

/*Prototypes*/
    /*Constructor*/
    Matrix createMatrix(const int,const int);
    Matrix createIdentityMatrix(const int);

    /*Matrix Manipulation*/
    Matrix sumMatrices(Matrix,Matrix);                  //creates new Matrix and returns the address
    Matrix multiplyScalarMatrix(const double,Matrix);   //alters the contents of Matrix without creating new Matrix
    Matrix multiplyMatrices(Matrix,Matrix);             //creates new Matrix and returns the address
    Matrix transposeMatrix(Matrix);                     //alters the contents of Matrix without creating new Matrix
    Matrix inverseMatrix(Matrix);
    Matrix adjugateMatrix(Matrix);
    double determinant(Matrix);                         //returns only the result
    double cofactor(Matrix,int,int);                    //returns only the result
    Matrix copyMatrix(Matrix);                          //creates new Matrix and returns the address
    Matrix minorOfMatrix(Matrix,int,int);               //creates new Matrix and returns the address
    Matrix interchangeCollumn(Matrix,int,int);          //alters the contents of Matrix without creating new Matrix
    Matrix interchangeRow(Matrix,int,int);              //alters the contents of Matrix without creating new Matrix

    /*Access*/
    void printMatrix(Matrix);
    Matrix changeInCoordinate(Matrix,int,int,double);
    boolean equalMatrix(Matrix,Matrix);
    boolean symmetricMatrix(Matrix);
    void emptyMatrix(Matrix);

    /*Analytic Geometry Operations*/
        /*Constructor*/
        Vector createVector(int);

        /*Vector Manipulation*/
        double dotProduct(Vector,Vector);
        Vector crossProduct(Vector,Vector);
        double angleBetween(Vector,Vector);
        double normOf(Vector);

