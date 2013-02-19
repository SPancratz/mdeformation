/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

    Exports the following functions:

        diagfrob(~numF, ~denF, A, n, d, p, N)

******************************************************************************/

load "basis.m";
load "padic_mat.m";

/*
    Returns the rising factorial $x (x + 1) \dotsm (x + k - 1)$ 
    for integers $x$ and $k$ with $k >= 0$.
 */
function _rfac(x, k)
    Z := Integers();
    y := Z!1;
    for i := 0 to k-1 do
        y := y * (x + i);
    end for;
    return y;
end function;

/*
    Computes the expression 
    \[
    \mu_m = 
        \sum_{k=0}^{\floor{m/2}} 2^{\floor{3m/4}-\nu_m-k} \frac{1}{(m-2k)! k!}
    \]
    modulo $p^N$ in the case that $p = 2$.

    Assumes that $m \geq 0$ and $N \geq 1$.
 */
function mu_2(m, N)
    Z := Integers();

    if (m eq 3) then  /* 4/3 */
        if (N le 2) then
            return Z!0;
        else
            return 4 * Modinv(3, 2^(N-2));
        end if;
    elif (m eq 7) then  /* 232/315 */
        if (N le 3) then
            return Z!0;
        else
            return (232 * Modinv(315, 2^(N - 3))) mod 2^N;
        end if;
    else
        f := Factorial(m);
        s := Z!1;
        t := Z!1;
        for k := 0 to (m div 2) - 1 do
            h := ((m - 2 * k - 1) * (m - 2 * k)) div 2;
            t := (t * h) div (k + 1);
            s := s + t;
        end for;

        // Now we have that 
        // f = m!
        // s = \sum_{k=0}^{\floor{m/2}} 2^{-k} \frac{m!}{(m-pk)! k!}

        v := Z!0;
        w := Z!0;
        _remove(~s, ~v, s, 2);
        _remove(~f, ~w, f, 2);
        v := v - w + ((3*m) div 4);

        if (v ge N) then
            return Z!0;
        else
            t := 2^(N - v);
            f := (Modinv(f, t) * s) mod t;
            f := f * 2^v;
            return f;
        end if;
    end if;
end function;

/*
    Computes the expression 
    \[
    \mu_m = \sum_{k=0}^{\floor{m/p}} p^{\floor{m/p} - k} \frac{1}{(m-pk)! k!}
    \]
    modulo $p^N$ in the case that $p > 2$ is an odd prime.

    The restriction that $p > 2$ is crucial as otherwise 
    the expression is not $p$-adically integral.

    Assumes that $m \geq 0$ and $N \geq 1$.
 */
function mu_p(m, p, N)
    Z := Integers();

    s := Z!1;
    t := Z!1;
    for k := 0 to (m div p) - 1 do
        t := t * _rfac(m - p*k - (p - 1), p);
        t := t div (p * (k + 1));
        s := s + t;
    end for;

    f := Factorial(m);

    // Now we have that 
    // f = m!
    // s = \sum_{k=0}^{\floor{m/p}} p^{-k} \frac{m!}{(m-pk)! k!}

    v := Z!0;
    w := Z!0;
    _remove(~s, ~v, s, p);
    _remove(~f, ~w, f, p);
    v := v - w + (m div p);

    if (v ge N) then
        return Z!0;
    else
        t := p^(N-v);
        f := (Modinv(f, t) * s) mod t;
        f := f * p^v;
        return f;
    end if;
end function;

function mu(m, p, N)
    if p eq 2 then
        return mu_2(m, N);
    else
        return mu_p(m, p, N);
    end if;
end function;

/*
    Computes the sequence $[\mu_1, \dotsc, \mu_M]$.
 */
function precompute_mu(M, p, N)
    return [ mu(m, p, N) : m in [1..M] ];
end function;

/*
    Let $R = \floor{M/p}$.  This functions computes the list of 
    powers of inverses $(1/d)^r$ modulo $p^N$ for $r = 1, \dotsc, R$.

    Assumptions:

        * $p$ prime
        * $p \nmid d$
        * $N \geq 1$
        * $M \geq 0$
 */
function precompute_dinv(M, d, p, N)
    R := M div p;
    L := [ 1..R ];
    if (R ge 1) then
        pN := p^N;
        L[1] := Modinv(d, pN);
        for r := 2 to R do
            L[r] := (L[r-1] * L[1]) mod pN;
        end for;
    end if;
    return L;
end function;

/*
    Returns a sequence of Teichmueller lifts of the elements of $A$
    to precision $p^N$.

    Assumes that the elements of $A$ are integers and returns the 
    list of Teichmueller lifts as integers.
 */
function precompute_lifts(A, n, p, N)
    Z := Integers();
    F := GF(p);
    R := pAdicRing(p, N);
    Q := quo< R | p^N >;

    L := [];

    for i := 1 to n+1 do
        t := TeichmuellerLift(F!A[i], Q);
        Append(~L, Z!t);
    end for;

    return L;
end function;

/*
    Computes the double sum involved in the expression for 
    $\alpha_{u+1, v+1}$ when $p = 2$, namely
    \[
    \sum_{m,r} 
        \Bigl(\frac{u_i + 1}{d}\Bigr)_r 2^{- \floor{(m+1)/4} + \nu} \mu_m.
    \]
 */
function dsum_2(DINV, MU, ui, vi, M, n, d, N)
    Z := Integers();

    m0 := (2 * (ui + 1) - (vi + 1)) div d;
    u  := ui + 1;
    pN := 2^N;

    x  := Z!0;
    f0 := Z!0;
    f1 := Z!0;
    f2 := Z!0;

    m := m0;
    while m le M do

        r := m div 2;

        // f0, f1, f2
        if (r eq 0) then
            f2 := Z!1;
        elif (r eq 1) then
            f1 := Z!1;
            f2 := u;
        elif (r eq 5) then
            f1 := f2;
            f2 := (u * (u + d) * (u + 2*d) * (u + 3*d) * (u + 4*d));
            if (m0 eq 0) then 
                f2 := f2 div 4;
            else
                f2 := f2 div 8;
            end if;
        else
            f0 := f1;
            f1 := f2;
            f2 := f0 * (((u + (r-2)*d) * (u + (r-1)*d)) div 2);
            f2 := f2 mod pN;
        end if;

        // g = f_{r} d^{-r} \mu_m
        //   = f2 * dinv[r] * mu[m]
        g := f2;
        if r gt 0 then
            g := (g * DINV[r]) mod pN;
        end if;
        if m gt 0 then
            g := (g * MU[m]) mod pN;
        end if;
        x := x + g;

        m := m + 2;
    end while;

    return x mod pN;
end function;

/*
    Computes the double sum involved in the expression for 
    $\alpha_{u+1, v+1}$ when $p$ is an odd prime, namely
    \[
    \sum_{m,r} \Bigl(\frac{u_i + 1}{d}\Bigr)_r p^{r - \floor{m/p}} \mu_m.
    \]
 */
function dsum_p(DINV, MU, ui, vi, M, n, d, p, N)
    Z := Integers();

    pN := p^N;
    x  := Z!0;

    r  := Z!0;
    m  := (p * (ui + 1) - (vi + 1)) div d;

    // Unrolled first iteration
    if m le M then
        h := p^(r - (m div p));
        if m eq 0 then
            g := h mod pN;
        else
            g := (h * MU[m]) mod pN;
        end if;
        x := x + g;
    end if;
    r := r + 1;
    m := m + p;

    f := Z!1;
    while m le M do

        f := f * (ui + 1 + (r - 1) * d);
        g := f * DINV[r];
        h := p^(r - (m div p));
        g := (g * h * MU[m]) mod pN;
        x := x + g;

        r := r + 1;
        m := m + p;
    end while;

    return x mod pN;
end function;

function dsum(DINV, MU, ui, vi, M, n, d, p, N)
    if p eq 2 then
        return dsum_2(DINV, MU, ui, vi, M, n, d, N);
    else
        return dsum_p(DINV, MU, ui, vi, M, n, d, p, N);
    end if;
end function;

/*
    Computes $\alpha_{u+1,v+1}$ modulo $p^N$ in the case when $p = 2$.
 */
function alpha_2(U, V, DINV, MU, M, n, d, N)

    ud := (n + 1 + &+U) div d;

    pN := 2^N;
    x  := 2^ud;

    for i := 1 to n+1 do
        g := dsum(DINV, MU, U[i], V[i], M, n, d, 2, N);
        x := (x * g) mod pN;
    end for;

    if (ud mod 2 ne 0) and (x ne 0) then
        x := pN - x;
    end if;

    return x;
end function;

/*
    Computes $\alpha_{u+1,v+1}$ modulo $p^N$ when $p > 2$ is an odd prime.
 */
function alpha_p(U, V, A, DINV, MU, M, n, d, p, N)

    ud := (n + 1 + &+U) div d;

    pN := p^N;
    x  := p^ud;

    for i := 1 to n+1 do
        e := (p * (U[i] + 1) - (V[i] + 1)) div d;
        f := A[i]^e mod pN;
        g := dsum(DINV, MU, U[i], V[i], M, n, d, p, N);
        x := (x * f * g) mod pN;
    end for;

    if (ud mod 2 ne 0) and (x ne 0) then
        x := pN - x;
    end if;

    return x;
end function;

function alpha(U, V, A, DINV, MU, M, n, d, p, N)
    if p eq 2 then
        return alpha_2(U, V, DINV, MU, M, n, d, N);
    else
        return alpha_p(U, V, A, DINV, MU, M, n, d, p, N);
    end if;
end function;

/*
    Computes the entry at position $(U,V)$ in the Frobenius matrix.

    Expects two sequences (that is, not lists) $U$ and $V$ so that 
    we can apply the '&' operator.
 */
procedure entry(~u, ~v, U, V, A, DINV, MU, M, C, n, d, p, N)
    Z := Integers();

    // Compute $f := (-1)^{u'+v'} (v'-1)! p^n$ exactly.
    ud := (n + 1 + &+U);
    vd := (n + 1 + &+V);

    f := Factorial(vd - 1);
    h := p^n;
    f := f * h;
    if (ud + vd) mod 2 ne 0 then
        f := -f;
    end if;

    // Compute $g := (u'-1)! \alpha_{u+1,v+1}$ to precision $N - n + 2 C$.
    g := Factorial(ud - 1);
    h := alpha(U, V, A, DINV, MU, M, n, d, p, N - n + 2 * C);
    g := g * h;

    // Set $(u, v)$ to the product of $f$ and $g^{-1}$.
    v := Z!0;
    w := Z!0;
    _remove(~f, ~v, f, p);
    _remove(~g, ~w, g, p);
    v := v - w;

    if v ge N then
        u := 0;
        v := 0;
    else
        t := p^(N - v);
        g := Modinv(g, t);
        u := (f * g) mod t;
    end if;
end procedure;

/*
    Computes the matrix for the action of Frobenius on the 
    smooth diagonal hypersurface of degree $d$ given by the 
    list of coefficients $A$ (expected as integers) modulo $p$.
 */
procedure diagfrob(~numF, ~denF, A, n, d, p, N)
    Z := Integers();
    S := PolynomialRing(Z, n + 1);

    r  := Valuation(Factorial(n - 1), p);
    s  := (n + 1) * Floor(Log(p, n - 1));
    C  := n + 2 * r + s;
    N2 := N - n + 2 * C;
    M  :=  ((p*p*N2) div (p-1)) + p*p*Floor(Log(p, (N2 div (p-1)) + 2)) + p*p * 4;

    ALIFT := precompute_lifts(A, n, p, N2);
    DINV  := precompute_dinv(M, d, p, N2);
    MU    := precompute_mu(M, p, N2);

    I := BasisSets(S, d);
    B := [];
    b := 0;
    for T in I do
        b := b + #T;
        for t in T do
            Append(~B, t);
        end for;
    end for;

    numF := ZeroMatrix(Z, b, b);
    denF := Z!1;

    for i := 1 to b do
        for j := 1 to b do

            // Extract exponents U and V from B[i] and B[j]
            U := [];
            V := [];
            k := 1;
            c := true;
            while k le n+1 and c do
                Append(~U, Degree(B[i], k));
                Append(~V, Degree(B[j], k));
                if (p * (U[k] + 1) - (V[k] + 1)) mod d ne 0 then
                    c := false;
                end if;
                k := k + 1;
            end while;

            if c eq false then
                numF[i, j] := Z!0;
            else
                u := Z!0;
                v := Z!0;
                entry(~u, ~v, U, V, ALIFT, DINV, MU, M, C, n, d, p, N);
                numF[i, j] := u * p^(v + (r + s));
            end if;
        end for;
    end for;

    denF := -(r + s);

    padic_mat_canonicalise(~numF, ~denF, p);
end procedure;

