/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz

    Exports the following functions:

        frobBound(n, p);
        prec_zeta_function(n, d, p, a)
        prec_frobq(n, d, p, a, N0)
        prec_frobp(n, p, a, N1)

******************************************************************************/

function frobBound(n, p)
    delta := Valuation(Factorial(n-1), p);
    for i := 2 to n-1 do
        delta := delta + Floor(Log(p,i));
    end for;
    return delta;
end function;

/*
    Returns a $p$-adic precision bound for the computation of 
    the zeta-function of a smooth projective hypersurface in 
    projective space $\mathbf{P}^n$ of degree $d$, defined over 
    $\mathbf{F}_q$ where $q = p^a$.

    If the reverse characteristic polynomial of $q^{-1} F_q$ is 
    computed to $p$-adic precision at least the return value, 
    then this uniquely determines the correct integer polynomial.

    Let $b$ denote the dimension of the cohomology vector space, 
    or equivalently the degree of the characteristic polynomial. 
    Whenever $b$ is even, this function returns a bound such 
    that one can compute the bottom half of the coefficients 
    correctly and then use the functional equation to complete 
    the computation.  Whenever $b$ is odd, the sign in the 
    functional equation is unknown in general and this function 
    returns a bound sufficient to determine all coefficients 
    directly.
 */

function prec_zeta_function(n, d, p, a)

    b := BasisSize(n, d);

    Z := Integers();
    N := Z!0;

    if (n eq 3 and p ne 2) then 
        f := Binomial(d - 1, 3);
        g := Binomial(b, b div 2);
        g := g * 2;
        N := a * f * Ceiling(Log(p, g));
    elif (n mod 2 eq 0) then 
        f := Binomial(b, b div 2);
        g := p^( (a * (b div 2) * (n - 1) + 1) div 2 );
        f := f * g;
        f := f * 2;
        N := Floor(Log(p, f)) + 1;
    else
        max := Z!0;
        for i := b div 2 to b do
            f := Binomial(b, i);
            g := p^( (a * i * (n - 1) + 1) div 2 );
            f := f * g;
            f := f * 2;
            max := Max(max, f);
        end for;
        N := Floor(Log(p, max)) + 1;
    end if;

    return N;
end function;

/*
    Returns $N_1s$ such that, in order to determine the coefficients 
    of the characteristic polynomial of $p^{-1} F_p$ to $p$-adic 
    precision $N_0$ it suffices to compute the matrix to precision $N_1$.

    When $X$ is a smooth projective surface and $p > 2$, we can do better 
    by instead computing $q^{h_{0,2}} f(T/q)$, in which case the extraction 
    of the polynomial requires an extra $a$ digits of precision.
 */

function prec_frobq(n, d, p, a, n, N0)
    if (n eq 3 and p ne 2) then
        return N0 + a;
    else
        return N0 + frobBound(n, p);
    end if;
end function;

function prec_frobp(n, p, a, N1)
    return N1 + (a - 1) * frobBound(n, p);
end function

