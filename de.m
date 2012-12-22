/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

load "padic_mat.m";

/*
    Given a matrix over a polynomial ring, returns a list of matrices 
    over the base ring.
 */
function MatrixToList(A)

    R := CoefficientRing(A);
    m := NumberOfRows(A);
    n := NumberOfColumns(A);

    d := -1;
    for i := 1 to b do
        for j := 1 to b do
            d := Max(d, Degree(A[i,j]));
        end for;
    end for;

    L := [];
    k := 0;
    while k le d do
        B := ZeroMatrix(R, m, n);
        for i := 1 to m do
            for j := 1 to n do
                B[i,j] := Coefficient(A[i,j], k);
            end for;
        end for;
        Append(~L, B);
        k := k + 1;
    end while;

    return L;
end function;

/*
    Returns the power series expansion around the origin
    of the solution to the matrix differential equation
    $(d/dt + M/r) C = 0$ with $C_0 = I$.

    Returns the solution in the form $p^{-c} C$, where
    $C$ is a matrix over integer polynomials and $c$ is
    an integer.

    Assumes that $M$ is a matrix over integer polynomials
    and $r$ is an integer polynomial.

    TODO:

    o   In general, we should allow $M$ to be a matrix 
        of polynomials over the rationals, and $r$ and 
        integer polynomial.
 */
procedure de_solve(~C, ~c, M, r, K, p, N, Nw)

    Z := Integers();
    b := NumberOfRows(C);

    r0   := Coefficient(r, 0);
    lenR := Degree(r) + 1;

    // TODO:  Support rational M
    M    := MatrixToList(M);
    lenB := #M;
    for k := 1 to lenB do
        M[k] := M[k] mod p^Nw;

    lD := [];
    ld := [];
    D  := IdentityMatrix(Z, b, b);
    d  := Z!1;
    Append(~lD, D);
    Append(~ld, d);

    for i := 1 to K-1 do
        D := ZeroMatrix(Z, b, b);
        d := Z!1;

        // TODO:  Support rational M
        for j := Max(0, i - lenB) to i-1 do
            // C[i] := C[i] + M[i-j] * C[j]
            v := Min(d, ld[j]);
            D := p^(d - v) * D + p^(ld[j] - v) * M[i-j] * lD[j];
            d := v;
        end for;

        for j := Max(0, i - lenR) + 1 to i-1 do
            // C[i] := C[i] + r[i-j] * j * C[j]
            v := Min(d, ld[j]);
            D := p^(d - v) * D + Coeff(r, i-j) * j * p^(ld[j] - v) * lD[j];
            d := v;
        end for;

        // D := D scalardiv t mod p^Nw
        t := - (i + 1) * r0;
        padic_mat_canonicalise(D, d, p);
        v := Z!0;
        _remove(~t, ~v, t, p);
        d := d + v;
        if d lt Nw then
            t := Modinv(t, p^(Nw - d));
            D := t * D;
        else 
            d := Z!0;
            D := ZeroMatrix(Z, b, b);
        end if;

        Append(~lD, D);
        Append(~ld, d);
    end for;

    for k := 1 to K do
        padic_mat_reduce(~lD[k], ~ld[k], p, N);
    end for;

    // Convert the pairs of lists (lD, ld) to (C, c)

    // Find valuation, abort if all lD[i] are zero
    isZero := true;
    c := Z!0;
    i := Z!1;
    while i lt #lD do
        if lD[i] ne 0 then
            if isZero then
                isZero := false;
                c := MatrixValuation(lD[i], ld[i]);
            else
                c := Min(c, MatrixValuation(lD[i], ld[i]));
            end if;
        end if;
    end while;
    if isZero then
        print "Exception (de_solve).  lD == 0.";
        exit;
    end if;

    S<X> := PolynomialRing(Z);
    C    := ZeroMatrix(S, b, b);

    k := K;
    while k ge 1 do
        for i := 1 to b do
            for j := 1 to b do
                // C[i,j][k-1] := p^(ld[k] - c) * lD[k][i,j];
                C[i,j] := C[i,j] + (p^(ld[k] - c) * lD[k][i,j]) * X^(k-1);
            end for;
        end for;
        k := k - 1;
    end while;

end procedure;
