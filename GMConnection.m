/******************************************************************************

    Copyright (C) 2009 - 2012 Sebastian Pancratz

******************************************************************************/

// This file contains code to compute the Gauss-Manin connection matrix of 
// a homogeneous polynomial.
// 

///////////////////////////////////////////////////////////////////////////////
// External files
// 

load "DFMatrix.m";
load "Basis.m";

// Given a vector v of coefficients and a list M of monomials, returns the
// corresponding polynomial.
//  
function ToPolynomialRows(v, M)
	S := Parent(M[1]);
	Q := S!0;
	for i:=1 to #M do
		Q := Q + v[i] * M[i];
	end for;
	return Q;
end function;

// Given a polynomial Q and a list M of monomials, returns a vector with
// the corresponding monomial coefficients.
// 
function ToVectorRows(Q, M)
	S := Parent(Q);
	R := CoefficientRing(S);
	v := Vector([R | MonomialCoefficient(Q, f) : f in M]);
	return v;
end function;

// Given a vector v of coefficients and a list M of lists of monomials, 
// for example, the list COLS in the code further below, returns the 
// corresponding polynomial.
// 
function ToPolynomialCols(v, M)
	// 
	// Find the first existing entry so we can set the parent.
	// 
	k := 1;
	while IsEmpty(M[k]) do
		k := k + 1;
	end while;
	S := Parent(M[k][1]);
	Q := S!0;
	i := 1;
	for k:=1 to #M do
		for f in M[k] do
			Q := Q + v[i] * f;
			i := i + 1;
		end for;
	end for;
	return Q;
end function;

// Given a polynomial Q and a list M of lists of monomials, for example, the
// list COLS in the code further below, returns the corresponding vector 
// of monomial coefficients.
// 
function ToVectorCols(Q, M)
	S := Parent(Q);
	R := CoefficientRing(S);
	L := [ ];
	for k:=1 to #M do
		for f in M[k] do
			Append(~L, MonomialCoefficient(Q, f));
		end for;
	end for;
	v := Vector(R, L);
	return v;
end function;

// Suppose f is a partial monomial with only entries in the first i-1
// exponents.  Suppose the current total degree of f is sum.  This function
// recursively calls itself to consider all extensions of f to the ith
// exponent.  If i exceeds the rank of S, that is, if the monomial is already
// completely specified, and the total degree of the monomial is d then it is
// added to the list L.
// 
procedure _GenerateMonomials(S, d, ~L, f, sum, i)
	if (i lt Rank(S)) then
		for exp:=0 to d-sum do
			_GenerateMonomials(S, d, ~L, f*(S.i)^exp, sum+exp, i+1);
		end for;
	else
		Append(~L, f*(S.i)^(d-sum));
	end if;
end procedure;

// Suppose S is a univariate or multivariate polynomial ring.  This function 
// generates all monomials of degree precisely d.
// 
function GenerateMonomials(S, d)
	L := [* *];
	_GenerateMonomials(S, d, ~L, S!1, 0, 1);
	return Reverse(L);
end function;

// Determines whether the monomial f is divisible by the monomial g.
// 
function MonomialIsDivisibleBy(f, g)
	S := Parent(f);
	n := Rank(S)-1;
	for i:=1 to n+1 do
		if Degree(f, i) lt Degree(g, i) then
			return false;
		end if;
	end for;
	return true;
end function;

// Returns the quotient of the monomials f and g, assuming that g divisides f.
// 
function MonomialDivision(f, g)
	assert MonomialIsDivisibleBy(f, g);
	S := Parent(f);
	n := Rank(S)-1;
	h := S!1;
	for i:=1 to n+1 do
		h := h * (S.i)^(Degree(f, i) - Degree(g, i));
	end for;
	return h;
end function;

// Consider the following set up, in which all polynomials are in the multivariate 
// polynomial ring S = R[X0,...,Xn] and they are assumed to be homogeneous.  We have
// Q of degree kd - (n+1), P of degree d, and look for polynomials A0, ..., An such 
// that Q = A0 P0 + ... + An Pn, where Pi is the partial derivative of P with 
// respect to Xi.
// 
// This method returns a matrix M in order to write this as Q = M (A0,...,An) modulo
// the space spanned by Bk, the set of monomials of total degree kd-(n+1) and 
// individual degree between 0 and d-2, inclusively.
// 
// We make the following assumptions.
//   o The (implicit) parameters n, d are natural numbers.
//   o The parameter P is a homogeneous polynomial of degree d in n+1 variables.
//   o The integer k is in the range 1, 2, ..., n+1.
//   o kd - (n+1) >= 0, and (k-1)d - n >= 0.  The second of these is stronger.
// 
function AuxMatrix(P, k)
	// 
	// Extract basic parameters.
	// 
	S := Parent(P);
	R := CoefficientRing(S);
	n := Rank(S)-1;
	d := TotalDegree(P);
	// 
	// Evaluate assertions!
	// 
	assert (k-1)*d ge n;
	assert k le n+1;
	// 
	// Compute derivatives DP[i].
	// 
	DP := [* *];
	for j:=1 to n+1 do
		Append(~DP, Derivative(P, j));
	end for;
	// 
	// Compute lists of all monomials of degree kd-(n+1) and (k-1)d-n, respectively, 
	// call these ROWSt and COLSt.
	// 
	ROWSt := GenerateMonomials(S, k*d - (n+1));
	COLSt := GenerateMonomials(S, (k-1)*d - n);
	// 
	// Build the row list, by extracting all monomials which are divisible by 
	// X_i^{d-1} for some parameter i.
	// 
	ROWS := [* *];
	for f in ROWSt do
		for j:=1 to n+1 do
			if Degree(f, S.j) ge d-1 then
				Append(~ROWS, f);
				break;
			end if;
		end for;
	end for;
	// 
	// Build the column list
	// 
	COLS := [* *];
	Append(~COLS, COLSt);
	for j:=1 to n do
		Append(~COLS, [* *]);
		for g in COLS[j] do
			if Degree(g, j) lt d-1 then
				Append(~(COLS[j+1]), g);
			end if;
		end for;
	end for;
	// 
	// Remove the two now unnecessary lists
	// 
	delete ROWSt;
	delete COLSt;
	// 
	// Compute the size of the desired matrix
	// 
	rows := #ROWS;
	cols := 0;
	for j:=1 to n+1 do
		cols := cols + #(COLS[j]);
	end for;
	// 
	// Compute the desired matrix
	// 
	M := ZeroMatrix(R, rows, cols);
	for row:=1 to rows do
		f := ROWS[row];
		col := 1;
		for j:=1 to n+1 do
			for g in COLS[j] do
				if MonomialIsDivisibleBy(f, g) then
					M[row,col] := MonomialCoefficient(DP[j], MonomialDivision(f, g));
				end if;
				col := col+1;
			end for;
		end for;
	end for;
	return M, ROWS, COLS;
end function;

// Given suitable polynomials P and Q, we look for polynomials A0, ..., An
// such that Q = P0 A0 + ... + Pn An.  But note that this is not quite what 
// we provide in the end.  Let us write this in matrix form, Q = M A.  Thus
// A = (1/det M) Adjoint(M) Q.
// 
// This function returns a rescaled list of homogeneous polynomials 
// A = (det M) [*A0 A1 ... An*], to ensure A0, ..., An lie in the same
// polynomial ring and avoid passing to the quotient field, and the 
// determinant det M.
// 
// N.B.  Not maintained since Week 5, Trinity Term 2009.
// 
function PolynomialDecomposition(P, Q)
	// 
	// Extract basic parameters
	// 
	S := Parent(P);
	R := CoefficientRing(S);
	n := Rank(S)-1;
	d := TotalDegree(P);
	k := (TotalDegree(Q)+(n+1)) div d;
	// 
	// Evaluate assertions!
	// 
	// The first assertion says that the degree of Q is of the right form.
	// The second assertion says that the problem might(?) be soluble, 
	// although we maybe also need the polynomials P to define a non-
	// singular variety for this to be the case.
	// 
	assert TotalDegree(Q) eq k*d - (n+1);
	assert (k-1)*d ge n;
	// 
	// Obtain the matrix M
	// 
 	M, ROWS, COLS := AuxMatrix(P, k);
	MAdj := Adjoint(M);
	Det := Determinant(M);
	// 
	// Extract Q in vector form
	// Compute the LHS of the matrix equation
	// 
	Qv := ToVectorRows(Q, ROWS);
	// 
	// DEBUG Output
	// 
	print "As a vector, Q =", Qv;
	// 
	// Compute the LHS of the matrix equation, i.e., 
	// Adj(M) Q = det(M) [A0 ... An].
	// 
	Av := Qv * Transpose(MAdj);
	// 
	// DEBUG Output
	// 
	print "As a vector, Adj(M) * Q =", Av;
	// 
	// Extract the desired polynomials
	// 
	A := [* *];
	col := 1;
	for j:=1 to n+1 do
		f := S!0;
		for g in COLS[j] do
			f := f + Av[col] * g;
			col := col+1;
		end for;
		Append(~A, f);
	end for;
	return A, Det;
end function;

// The method uses the same algorithm as in the method 
// "PolynomialDecomposition", but includes further parameters to avoid repeated
// computations when computing connection matrices.
// 
// Given the auxiliary matrix M = M(P,k) with determinant Delta = Delta(P, k), 
// let MAdj be the adjoint of M.  The lists ROWS and COLS denote the respective
// index sets for M.  Note that COLS is given as a list of lists.
// 
function _PolynomialDecomposition(MAdj, ROWS, COLS, P, Q)
	// 
	// Extract basic parameters
	// 
	S := Parent(P);
	R := CoefficientRing(S);
	n := Rank(S)-1;
	d := TotalDegree(P);
	k := (TotalDegree(Q)+(n+1)) div d;
	// 
	// Evaluate assertions!
	// 
	// The first assertion says that the degree of Q is of the right form.
	// The second assertion says that the problem might(?) be soluble, 
	// although we maybe also need the polynomials P to define a non-
	// singular variety for this to be the case.
	// 
	assert TotalDegree(Q) eq k*d - (n+1);
	assert (k-1)*d ge n;
	assert k le n+1;
	// 
	// Extract Q in vector form
	// Compute the LHS of the matrix equation
	// 
	Qv := ToVectorRows(Q, ROWS);
	Av := Qv * Transpose(MAdj);
	A := [* *];
	col := 1;
	for j:=1 to n+1 do
		f := S!0;
		for g in COLS[j] do
			f := f + Av[col] * g;
			col := col + 1;
		end for;
		Append(~A, f);
	end for;
	return A;
end function;

// Suppose that P is a multivariate polynomial over the ring R, and that 
// R is a ring with one free variable supporting the call Derivative(-), 
// e.g. R = Q[t] or Q(t).  This method returns the derivative of P with 
// respect to t.
// 
function TDerivative(P)
	// 
	// Extract basic parameters.
	// 
	S := Parent(P);
	// 
	// Extract coefficients and monomials from P.
	//
	C := Coefficients(P);
	M := Monomials(P);
	// 
	// Compute the derivative of P with respect to t.
	// 
	Pd := S!0;
	for j:=1 to #C do
		Pd := Pd + Derivative(C[j]) * M[j];
	end for;
	return Pd;
end function;

// Given a multivariate polynomial P over a polynomial ring of rank n+1,
// returns a list of the n+1 derivatives.
// 
function Derivatives(P)
	// 
	// Extract basis parameters
	// 
	S := Parent(P);
	n := Rank(S)-1;
	// 
	// Compute derivatives DP[i].
	// 
	DP := [* *];
	for j:=1 to n+1 do
		Append(~DP, Derivative(P, j));
	end for;
	return DP;
end function;

// Given a polynomial f in a multivariate polynomial ring over a ring R and a
// list L of monomials, returns whether f lies in the R-span of M.
// 
function IsInSpan(f, L)
	M := Monomials(f);
	S := { Parent(f) | };
	for m in L do
		Include(~S, m);
	end for;
	for m in M do
		if not(m in S) then
			return false;
		end if;
	end for;
	return true;
end function;

// Assumes that Q is a polynomial of degree kd - (n+1) that is to be reduced,
// where k = 2, ..., n, n+1 is sufficiently large for this to make sense.
// 
// Returns a list [* gamma(1) ... gamma(k) *] of elements gamma(i) in the span 
// of Bi such that Q is given by their sum.
// 
function Reduce(P, dP, M, ROWS, COLS, MAdj, MDet, B, Q, k)
	// 
	// Extract basic parameters.
	// 
	S := Parent(P);
	n := Rank(S)-1;
	// 
	// Set up the output variable.
	// 
	L := [* *];
	// 
	// Note that gamma(n+1) is necessarily zero.  Thus we first reduce to
	// the case where k is at most n.
	// 
	if k eq n+1 then
		A := _PolynomialDecomposition(MAdj[k], ROWS[k], COLS[k], P, Q);
		for i:=1 to n+1 do
			A[i] := A[i] / MDet[k];
		end for;
		/**************************************************************
		// DEBUG Output
		LHS := Q;
		for i:=1 to n+1 do
			LHS := LHS - A[i] * dP[i];
		end for;
		print "  Reduce Q =", Q;
		print "  k =", k;
		print "  Q - A(0) dP(0) - ... - A(n) dP(0) =", LHS;
		print "  As B(n+1) is empty, this should be 0.";
		**************************************************************/
		Q := S!0;
		for i:=1 to n+1 do
			Q := Q + Derivative(A[i], i) / (k-1);
		end for;
		k := k - 1;
	end if;
	// 
	// Now we may assume that k is between 2 and n.
	// 
	while not IsInSpan(Q, B[k]) do
		A := _PolynomialDecomposition(MAdj[k], ROWS[k], COLS[k], P, Q);
		for i:=1 to n+1 do
			A[i] := A[i] / MDet[k];
		end for;
		/**************************************************************
		// DEBUG Output
		LHS := Q;
		for i:=1 to n+1 do
			LHS := LHS - A[i] * dP[i];
		end for;
		print "  Reduce Q =", Q;
		print "  k =", k;
		print "  Q - A(0) dP(0) - ... - A(n) dP(0) =", LHS;
		print "  Does the difference lie in B(k)?", IsInSpan(LHS, B[k]);
		**************************************************************/
		// 
		// Use Q now for what will be gamma(k).
		// 
		for i:=1 to n+1 do
			Q := Q - A[i] * dP[i];
		end for;
		Append(~L, Q);
		Q := S!0;
		for i:=1 to n+1 do
			Q := Q + Derivative(A[i], i) / (k-1);
		end for;
		k := k - 1;
	end while;
	if k ge 1 then
		Append(~L, Q);
		k := k - 1;
	end if;
	while k ge 1 do
		Append(~L, S!0);
		k := k - 1;
	end while;
	/**********************************************************************
	// DEBUG Output
	print "  Return L =", Reverse(L);
	**********************************************************************/
	return Reverse(L);
end function;

// Assumes that Q is a polynomial of degree kd - (n+1) that is to be reduced,
// where k = 2, ..., n, n+1 is sufficiently large for this to make sense.
// 
// Let K be Max(k, n).  This function returns two lists 
// [* gamma'(1) ... gamma'(K) *] and [* c(1) ... c(K) *] where $\gamma'(j)$ is
// in the span of $b_j$ such that, setting gamma(j) = gamma'(j) / c(j), we have
// $Q \Omega / P^k = \sum_{j=1}^K \gamma(j) \Omega / P^j$.
// 
function _Reduce(P, dP, M, ROWS, COLS, MAdj, MDet, B, Q, k)
	// 
	// Extract basic parameters.
	// 
	S := Parent(P);
	R := CoefficientRing(S);
	n := Rank(S)-1;
	// 
	// Set up the output variable.
	// 
	L := [* *];
	C := [* *];
	// 
	// The variable to contain the denominator.
	// 
	c := R ! 1;
	// 
	// Note that gamma(n+1) is necessarily zero.  Thus we first reduce to
	// the case where k is at most n.
	// 
	if k eq n+1 then
		A := _PolynomialDecomposition(MAdj[k], ROWS[k], COLS[k], P, Q);
		c := c * MDet[k];
		Q := S ! 0;
		for i:=1 to n+1 do
			Q := Q + Derivative(A[i], i);
		end for;
		c := c * (k - 1);
		k := k - 1;
	end if;
	// 
	// Now we may assume that k is between 2 and n.
	// 
	while not IsInSpan(Q, B[k]) do
		A := _PolynomialDecomposition(MAdj[k], ROWS[k], COLS[k], P, Q);
		// 
		// Use Q now for what will be gamma(k).
		// 
		Q := MDet[k] * Q;
		for i:=1 to n+1 do
			Q := Q - A[i] * dP[i];
		end for;
		c := c * MDet[k];
		Append(~L, Q);
		Append(~C, c);
		Q := S ! 0;
		for i:=1 to n+1 do
			Q := Q + Derivative(A[i], i);
		end for;
		c := c * (k - 1);
		k := k - 1;
	end while;
	if k ge 1 then
		Append(~L, Q);
		Append(~C, c);
		k := k - 1;
	end if;
	while k ge 1 do
		Append(~L, S ! 0);
		Append(~C, R ! 1);
		k := k - 1;
	end while;
	return Reverse(L), Reverse(C);
end function;

function GMConnection(P)
	// 
	// Extract basic parameters.
	// 
	S := Parent(P);
	R := CoefficientRing(S);
	n := Rank(S)-1;
	d := TotalDegree(P);
	
	// 
	// Output:
	//   Print input
	// 
	print "Base ring R =", R;
	print "Polynomial ring S =", S;
	print "n =", n;
	print "d =", d;
	print "Defining polynomial P =", P;
	
	// 
	// Compute the derivative of P with respect to t, the variable of the 
	// coefficient ring.
	//
	dP   := Derivatives(P);
	dPdT := TDerivative(P);
	
	// 
	// Output:
	//   Print the derivative
	// 
	print "dP/dX =", dP;
	print "dP/dt =", dPdT;
	
	// 
	// Construct the index set I and its "size" lI.
	// 
	I  := BasisSets(S, d);
	lI := 0;
	for B in I do
		lI := lI + #B;
	end for;
	
	// 
	// Output:
	//   Print the sets Bk and their combined size
	// 
	print "Basis sets =", I;
	print "Total size =", lI;
	print "Constructing auxiliary matrices..";
	
	// 
	// Construct the auxiliary matrices
	// Note that a priori we might need them for k = 1, ..., n, n+1. 
	// The rows and columns sets are not non-empty unless 
	// (k-1)*d >= n.
	//
	M    := [* *];
	ROWS := [* *];
	COLS := [* *];
	MAdj := [* *];
	MDet := [* *];
	k := 1;
	while (k-1)*d lt n do
		// 
		// Output
		// 
		print "  Considering k =", k, "..";
		
		Append(~M, [* *]);
		Append(~ROWS, [* *]);
		Append(~COLS, [* *]);
		Append(~MAdj, [* *]);
		Append(~MDet, [* *]);
		
		// Output
		//
		print "    M(P,k) is a 0 x 0 matrix";
		
		k := k + 1;
	end while;
	while k le n+1 do
		// 
		// Output
		// 
		print "  Considering k =", k, "..";
		
		tempM, tempROWS, tempCOLS := AuxMatrix(P, k);
		
		tempd, tempA := DFComputeAll(tempM);
		
		Append(~M, tempM);
		Append(~ROWS, tempROWS);
		Append(~COLS, tempCOLS);
		Append(~MAdj, tempA);
		Append(~MDet, tempd);
		
		// Output
		//
		print "    M(P,k) is a", Nrows(tempM), "x", Ncols(tempM), "matrix";
		
		k := k + 1;
	end while;
	delete k;
	
	//
	// Output
	// 
	print "Done.";
	print "Construct the connection matrix..";
	
	// 
	// Construct the matrix.
	//   The Gauss-Manin connection matrix the has entries 
	//   GM[i, j] / D[i, j].  Later, we arrange this in the form 
	//   GM[i, j] / r, where r is an element of R, the polynomial
	//   least common multiple of the denominators.
	// 
	GM  := ZeroMatrix(R, lI, lI);
	D   := ScalarMatrix(R, lI, R ! 1);
	col := 1;
	for colk := 1 to n do
		for g in I[colk] do
			// 
			// Output
			// 
			print "  Consider column:", col;
			print "  Reduce the element corresponding to g =", g;
			// 
			// Set and reduce Q
			//
			Q := - colk * g * dPdT;
			L, C := _Reduce(P, dP, M, ROWS, COLS, MAdj, MDet, I, Q, colk+1);

			// Output
			// 
			print "  L = ", L;
			print "  C = ", C;

 			// 
			// Extract the column vector
			// 
			// Iterate over all rows..
			row := 1;
			for rowk := 1 to #L do
				for f in I[rowk] do
					GM[row, col] := MonomialCoefficient(L[rowk], f);
					D[row, col]  := C[rowk];
					
					// For some unknown reason, the GCD is 
					// computed to be 0 by Magma, leading 
					// to a division by 0.  Thus we take 
					// out this code.
					// 
					// Cancel common factors
					gcd := GCD(GM[row, col], D[row, col]);
					if gcd ne R ! 1 then
						GM[row, col] := GM[row, col] / gcd;
						D[row, col] := D[row, col] / gcd;
					end if;
					
					// Increment the row index
					row := row + 1;
				end for;
			end for;
			// Increment the column index
			col := col + 1;
		end for;
	end for;
	
	// 
	// Output
	//
	print "Clearing denomiators..";
	
	//
	// Write the Gauss-Manin connection matrix in the form GM / r, 
	// where r is the least common multiple of all denominators.
	//
	r := R ! 1;
	for i:=1 to lI do
		for j:=1 to lI do
			print "  i =", i, ", j =", j;
			print "  r =", r;
			print "  D[i,j] =", D[i,j];
			r := LeastCommonMultiple(r, D[i, j]);
		end for;
	end for;
	for i:=1 to lI do
		for j:=1 to lI do
			print "  i =", i, "j =", j;
			GM[i, j] := GM[i, j] * (r div D[i, j]);
		end for;
	end for;

	return GM, r;
end function;

