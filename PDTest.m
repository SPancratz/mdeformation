/******************************************************************************

    Copyright (C) 2009 - 2012 Sebastian Pancratz

******************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// This file contains code to test the polynomial decomposition and 
// reduction process.
// 

load "GMConnection.m";

R<t> := FieldOfFractions(PolynomialRing(Rationals()));
S<X,Y,Z> := PolynomialRing(R, 3);

P := (1-t) * (X^4 + Y^4 + Z^4) + t * (X^3 * Y + Y^3 * Z + Z^3 * X);
Q := X^5 - X^4 *Y - X^2 * Z^3 + X * Y^4 - X * Y^3 * Z + X * Z^4;

// S<X,Y,Z> := PolynomialRing(Rationals(), 3);
// P := X^4 + Y^4 + Z^4;
// Q := X^5 - X^4 * Y - X^2 * Z^3 + X * Y^4 - X * Y^3 * Z + X * Z^4;
// P := X^3 + Y^3 + Z^3;
// Q := X^3 - X^2 * Y + X * Z^2 + Y^2 * Z; 

A, Det := PolynomialDecomposition(P, Q);

for i:=1 to #A do
	A[i] := A[i] / Det;
end for;

print "P =", P;
print "Q =", Q;
print "A =", A;

B := BasisSets(S, TotalDegree(P));

print "Basis sets:", B;

M, ROWS, COLS := AuxMatrix(P, 2);

LHS := ToPolynomialCols(ToVectorRows(Q, ROWS) * Transpose(M^(-1)), COLS);
for i:=1 to 3 do
	LHS := LHS - A[i];
end for;

print "M^(-1) Q - ( A(0) - A(1) - A(2) ) =", LHS;

print "Monomials in the difference:", Monomials(LHS);

print "Is the difference in B(k)?", IsInSpan(LHS, B[2]);
