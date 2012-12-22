/******************************************************************************

    Copyright (C) 2009 - 2012 Sebastian Pancratz

******************************************************************************/

// Given a matrix A, this method returns the minimum valuation of its 
// entries at p.
// 
// Assumes that the call "Valuation(a, p)" is supported for each entry
// a of A.
// 
function MatrixValuation(A, p)
	Q := Rationals();
	v := Valuation(Q ! A[1, 1], p);
	for i:=1 to Nrows(A) do
		for j:=1 to Ncols(A) do
			v := Min(v, Valuation(Q ! A[i, j], p));
		end for;
	end for;
	return v;
end function;

function PowerSeriesValuation(C, p)
	L := [* *];
	for i:=1 to #C do
		Append(~L, MatrixValuation(C[i], p));
	end for;
	return L;
end function;

load "GMConnection.m";
load "PowerSeries.m";

///////////////////////////////////////////////////////////////////////////////
// Test setup
// 

print "============================";
print "Compute GM Connection Matrix";
print "============================";

R<t>      := PolynomialRing(Rationals());
F         := FieldOfFractions(R);
S<X,Y,Z>  := PolynomialRing(R, 3);

// 
// R<t>     := PolynomialRing(pAdicField(5, 26));
// S<X,Y,Z> := PolynomialRing(R, 3);
// 

P := (1-t) * (X^3 + Y^3 + Z^3) + t * (X^3 + Y^3 + Z^3 + X*Y*Z);
// P = (1-t) * (X^3 + Y^3 + Z^3) + t * (X^3 + Y^3 + Z^3 + (X+2*Y+Z)*(X+3*Y+4*Z)^2);
// P := (1-t) * (X^4 + Y^4 + Z^4) + t * (X^3 * Y + Y^3 * Z + Z^3 * X);
// P := (1-t) * (X^4 + Y^4 + Z^4) + t * (X^4+Y^4+Z^4+(X+2*Y+Z)^2*(X+3*Y+4*Z)^2);
// P := (1-t) * (X^5 + Y^5 + Z^5) + t * (X^5+Y^5+Z^5+(X+2*Y+Z)^2*(X+3*Y+4*Z)^3);

M, r := GMConnection(P);

print M;
print r;

print "=====";
print "Done.";
print "=====";

/******************************************************************************

// 
// LCM of denominators
// 

b := ZeroMatrix(R, Nrows(M), Ncols(M));
r := R ! 1;

for i:=1 to Nrows(M) do
	for j:=1 to Ncols(M) do
		r := LeastCommonMultiple(r, R ! Denominator(M[i, j]));
	end for;
end for;

for i:=1 to Nrows(b) do
	for j:=1 to Ncols(b) do
		b[i, j] := R ! (r * M[i, j]);
	end for;
end for;

// Output
// 
print "Matrix after clearing denominators b =", b;
print "LCM of denominators r =", r;

// 
// Check that r is square-free
//

dr := Derivative(r);
g := GreatestCommonDivisor(r, dr);

// Output
// 
print "Check that r is square-free:";
print "GCD(r, dr/dt) =", g;

///////////////////////////////////////////////////////////////////////////////
// Roots
// 

FACT := Factorization(r);

for k:=1 to #FACT do
	f := FACT[k][1];
	print "Consider the factor f =", f;
	K<c> := ext<Rationals() | f>;
	RM := ZeroMatrix(K, Nrows(b), Ncols(b));
	for i:=1 to Nrows(RM) do
		for j:=1 to Ncols(RM) do
			RM[i, j] := Evaluate(b[i, j], c) / Evaluate(dr, c);
		end for;
	end for;
	
	// Output
	// 
	// print "Residue matrix RM =", RM;
	print "Eigenvalues of RM =", Eigenvalues(RM);
end for;

///////////////////////////////////////////////////////////////////////////////
// Infinity
// 

// Let d be the degree of r, and D the degree of the matrix b.
// 

d := Degree(r);
D := Degree(b[1, 1]);
for i:=1 to Nrows(b) do
	for j:=1 to Ncols(b) do
		D := Max(D, Degree(b[i, j]));
	end for;
end for;

// Output
// 
print "Degree(r) =", d;
print "Degree(b) =", D;

// Let N be the matrix of (d-1)th coefficients of b.
// 

N := ZeroMatrix(CoefficientRing(R), Nrows(b), Ncols(b));

if D ge d-1 then
	for i:=1 to Nrows(b) do
		for j:=1 to Ncols(b) do
			N[i, j] := Coefficient(b[i, j], d-1);
		end for;
	end for;
end if;

// Output
// 
print "N is the matrix of (d-1)th coefficients in b.";
print "N =", N;

E := ZeroMatrix(Integers(), Nrows(b), Ncols(b));

for i:=1 to Nrows(b) do
	for j:=1 to Ncols(b) do
		E[i, j] := Degree(b[i, j]) - (d-1);
	end for;
end for;

// Output
// 
print "E is the matrix defined by E[i,j] = deg(b[i,j]) - (d-1).";
print "E =", E;

// Compute power series
//
 C := PowerSeriesP(M, 1000, 7, 10);
 V := PowerSeriesValuation(C, 7);

******************************************************************************/

