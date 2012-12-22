/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

    This file contains auxiliary functionality for treating 
    pairs $(U, v)$ $p$-adic matrices $p^v U$.

******************************************************************************/

/*
    Let $f = p^v u$.
 */
procedure _remove(~u, ~v, f, p)
    if f eq 0 then
        v := 0;
        u := 0;
    else
        v := Valuation(f, p);
        u := f div (p^v);
    end if;
end procedure;

/*
    Returns the rising factorial $x (x + 1) \dotsm (x + k - 1)$.
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
    Returns the $p$-adic valuation of the integer matrix $A$.

    By our convention, the valuation of the zero matrix is $0$.
 */
function MatrixValuation(A, p)
    m := NumberOfRows(A);
    n := NumberOfColumns(A);

    if m eq 0 or n eq 0 then
        return 0;
    else
        c := true;  // Zero?
        v := 0;
        for i := 1 to m do
            for j := 1 to n do
                if A[i,j] ne 0 then
                    if c then
                        v := Valuation(A[i,j], p);
                        c := false;
                    else
                        v := Min(v, Valuation(A[i,j], p));
                    end if;
                end if;
            end for;
        end for;

        if c then
            return 0;
        else
            return v;
        end if;
    end if;
end function;

/*
    Given an integer matrix $U$ and an integral exponent $v$,
    representing the matrix $p^{-v} U$, brings $(U,v)$ into
    canonical form.
 */
procedure padic_mat_canonicalise(~U, ~v, p)
    if U eq 0 then
        v := 0;
    else
        w := MatrixValuation(U, p);
        if w ne 0 then
            U := U div p^w;
            v := v + w;
        end if;
    end if;
end procedure;

/*
    Reduces the matrix $p^v U$ modulo $p^N$.
 */
procedure padic_mat_reduce(~U, ~v, p, N)

    padic_mat_canonicalise(~U, ~v, p);

    if v ge N then
        P := Parent(U);
        U := P!0;
    else
        U := U mod p^(N-v);
    end if;
end procedure;

