/******************************************************************************

    Copyright (C) 2009 - 2012 Sebastian Pancratz

******************************************************************************/

// This files contains code to numerically compute a power series solution to 
// the matrix differential equation (d/dt + M) C = 0.  The set up of this is 
// as follows:
// 
// Let R be a polynomial ring in the indeterminate t over a ring A.  Let F 
// be the field of fractions of R, i.e. the field of rational functions in t.
// M is a square matrix with entries in F, and if we clear denominators by 
// writing M = b/r with b a matrix with entries in R and r a polynomial in R
// then we assume that the constant term of r is non-zero.
// 

// Works out a power series solution to the matrix differential equation.
// Considers the coefficients up to t^(N-1), that is, returns a list with 
// N entries.
// 
function PowerSeries(M, N)
    // 
    // Extract basis parameters.
    // 
    F := Parent(M[1, 1]);
    A := CoefficientRing(F);
    R := PolynomialRing(A);
    
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
    
    // 
    // Let blist[k+1] be the matrix of kth coefficients of b.
    // 
    blist := [* *];
    for k:=0 to D do
        T := ZeroMatrix(A, Nrows(b), Ncols(b));
        for i:=1 to Nrows(b) do
            for j:=1 to Ncols(b) do
                T[i, j] := Coefficient(b[i, j], k);
            end for;
        end for;
        Append(~blist, T);
    end for;
    
    // 
    // Compute power series
    // 
    L := [* *];
    Append(~L, ScalarMatrix(Nrows(M), A ! 1));
    for i:=0 to N-2 do
        C := ZeroMatrix(A, Nrows(M), Ncols(M));
        for j:=Max(0,i-D) to i do
            C := C - blist[i-j+1] * L[j+1];
        end for;
        for j:= Max(0, i-d) to i-1 do
            C := C - Coefficient(r, i-j) * (j+1) * L[j+1+1];
        end for;
        C := ( 1 / (Coefficient(r, 0) * (i+1)) ) * C;
        Append(~L, C);
    end for;
    return L;
end function;

// Works out a power series solution to the matrix differential equation.
// Considers the coefficients up to t^(N-1), that is, returns a list with 
// N entries.
// 
// Here we assume that A is the field of rational numbers.  But instead of
// performing exact arithmetic over the rationals, we coerce the data to 
// $\Q_p$ with fixed precision to speed up computations.  The result is coerced
// to the rationals again.
// 
function PowerSeriesP(M, N, p, precision)
    // 
    // Extract basis parameters.
    // 
    F := Parent(M[1, 1]);
    A := CoefficientRing(F);
    R := PolynomialRing(A);
    
    // 
    // Validate some assertions!
    // 
    assert A eq Rationals();
    
    Q   := Rationals();
    Qp  := pAdicField(p, precision);
    MQ  := MatrixAlgebra(Q, Nrows(M));
    MQp := MatrixAlgebra(Qp, Nrows(M));
    
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
    
    // 
    // Let blist[k+1] be the matrix of kth coefficients of b.
    // 
    blist := [* *];
    for k:=0 to D do
        T := ZeroMatrix(Qp, Nrows(b), Ncols(b));
        for i:=1 to Nrows(b) do
            for j:=1 to Ncols(b) do
                T[i, j] := Qp ! Coefficient(b[i, j], k);
            end for;
        end for;
        Append(~blist, T);
    end for;
    
    // 
    // Compute power series
    // 
    L := [* *];
    Append(~L, ScalarMatrix(Nrows(M), Qp ! 1));
    for i:=0 to N-2 do
        C := ZeroMatrix(Qp, Nrows(M), Ncols(M));
        for j:=Max(0,i-D) to i do
            C := C - blist[i-j+1] * L[j+1];
        end for;
        for j:= Max(0, i-d) to i-1 do
            C := C - Qp ! Coefficient(r, i-j) * (j+1) * L[j+1+1];
        end for;
        C := ( 1 / (Qp ! Coefficient(r, 0) * (i+1)) ) * C;
        Append(~L, MQ ! C);
    end for;
    return L;
end function;

