/******************************************************************************

    Copyright (C) 2009 - 2012 Sebastian Pancratz

******************************************************************************/

// This file contains implementations of division-free algorithms for the 
// computation of the characteristic polynomial and adjoint of a matrix over 
// a ring.
// 

// Computes the charateristic polynomial.
// 
function DFCharacteristicPolynomial(M)
	
	// Extract parameters
	//
	m := Nrows(M);
	n := Ncols(M);
	R := CoefficientRing(M); 
	
	// Validate assertions
	//
	assert m eq n;
	
	// F  A matrix over R such that F[p, t] is the coefficient of X^{t-p} 
	//    in the polynomial $F_t$.
	// a  A list of lists such that a[p, t] is a vector of length t.
	// A  A matrix over R.
	// 
	F := ZeroMatrix(R, n, n);
	a := [* *];
	for i:=1 to n do
		Append(~a, [* *]);
		for j:=1 to n do
			Append(~(a[i]), [* *]);
		end for;
	end for;
	A := ZeroMatrix(R, n, n);
	
	F[1, 1] := - M[1, 1];
	for t:=2 to n do
		// Set a(1, t) to be M(<=t, t)
		// 
		for i:=1 to t do
			Append(~(a[1, t]), M[i, t]);
		end for;
		// Set A[1, t] to be the (t)th entry in a[1, t]
		// 
		A[1, t] := M[t, t];
		for p:=2 to t-1 do
			// Set a(p, t) to the product of M[<=t, <=t] * a(p-1, t)
			// 
			for i:=1 to t do
				s := R ! 0;
				for j:=1 to t do
					s := s + M[i, j] * a[p-1, t][j];
				end for;
				Append(~(a[p, t]), s);
			end for;
			// Set A[p, t] to be the (t)th entry in a[p, t]
			A[p, t] := a[p, t][t];
		end for;
		// Set A[t, t] to be M[t, <=t] * a(p-1, t)
		//
		for j:=1 to t do
			A[t, t] := A[t, t] + M[t, j] * a[t-1, t][j];
		end for;
		for p:=1 to t-1 do
			F[p, t] := F[p, t-1];
			for k:=1 to p-1 do
				F[p, t] := F[p, t] - A[k, t] * F[p-k, t];
			end for;
			F[p, t] := F[p, t] - A[p, t];
		end for;
		for k:=1 to t-1 do
			F[t, t] := F[t, t] - A[k, t] * F[t-k, t];
		end for;
		F[t, t] := F[t, t] - A[t, t];
	end for;
	
	S<X> := PolynomialRing(R);
	f := S ! (X^n);
	for p:=1 to n do
		f := f + F[p, n] * X^(n-p);
	end for;
	
	return f;
end function;

// Computes the adjoint.
// 
// Assumes that f is the characteristic polynomial of M.
// 
function _DFAdjoint(M, f)
	
	// Extract parameters
	//
	m := Nrows(M);
	n := Ncols(M);
	R := CoefficientRing(M); 
	
	// Validate assertions
	//
	assert m eq n;
	assert Degree(f) eq n;
	
	// Let A denote the adjoint of M, which we want to compute, and 
	// N denote a copy of M used to store powers of M.
	// 
	A := ScalarMatrix(n, Coefficient(f, 1));
	N := ScalarMatrix(n, R ! 1);
	for i:=1 to n-1 do
		// Set N to be M^i
		// 
		N := N * M;
		A := A + Coefficient(f, i+1) * N;
	end for;
	if IsEven(n) then
		A := - A;
	end if;
	
	return A;
end function;

function DFAdjoint(M)

	// Extract parameters
	//
	m := Nrows(M);
	n := Ncols(M);
	R := CoefficientRing(M); 
	
	// Validate assertions
	//
	assert m eq n;
	
	f := DFCharacteristicPolynomial(M);
	
	return _DFAdjoint(M, f);
end function;

function _DFDeterminant(M, f)
		
	// Extract parameters
	//
	m := Nrows(M);
	n := Ncols(M);
	R := CoefficientRing(M); 
	
	// Validate assertions
	//
	assert m eq n;
	assert Degree(f) eq n;
	
	d := Coefficient(f, 0);
	if IsOdd(n) then
		d := -d;
	end if;
	
	return d;
end function;

function DFDeterminant(M)
	
	// Extract parameters
	//
	m := Nrows(M);
	n := Ncols(M);
	R := CoefficientRing(M); 
	
	// Validate assertions
	//
	assert m eq n;
	
	f := DFCharacteristicPolynomial(M);
	
	return _DFDeterminant(M, f);
end function;

// Computes the determinant, the adjoint and the characteristic polynomial of 
// the matrix M.
// 
function DFComputeAll(M)
	
	// Extract parameters
	//
	m := Nrows(M);
	n := Ncols(M);
	R := CoefficientRing(M); 
	
	// Validate assertions
	//
	assert m eq n;
	
	f := DFCharacteristicPolynomial(M);
	A := _DFAdjoint(M, f);
	d := _DFDeterminant(M, f);
	
	return d, A, f;
end function;

