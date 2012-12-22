/******************************************************************************

    Copyright (C) 2009 - 2012 Sebastian Pancratz

******************************************************************************/

load "DFMatrix.m";

Q := Rationals();

M := ZeroMatrix(Q, 2, 2);
M[1,1] := 2;
M[1,2] := 1;
M[2,1] := 3;
M[2,2] := 7;

print "Magma: Characteristic polynomial";
time CharacteristicPolynomial(M);
time Adjoint(M);
time Determinant(M);

print "DFMatrix: Characteristic polynomial";
time DFCharacteristicPolynomial(M);
time DFAdjoint(M);
time DFDeterminant(M);

