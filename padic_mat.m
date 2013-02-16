/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

    This file contains auxiliary functionality for treating pairs $(U, v)$ 
    of integer matrices and integers  as $p$-adic matrices $p^v U$.

******************************************************************************/

/*
    Writes the integer $f$ in the form $f = p^v u$ where $u$ is a $p$-adic 
    unit whenever $f$ is non-zero.  When $f = 0$, sets $u = v = 0$.
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
    Returns the $p$-adic valuation of the integer matrix $A$ 
    whenever $A$ is non-zero and $0$ otherwise.
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
    Given an integer matrix $U$ and an integral exponent $v$ 
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

