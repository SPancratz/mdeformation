/******************************************************************************

    (C) 2013 Sebastian Pancratz

******************************************************************************/

load "precisions.m";

function _revcharpoly(F, n, N0)

    Qq := BaseRing(Parent(F));
    p  := Prime(Qq);
    a  := Degree(Qq);
    b  := Nrows(F);
    if (n mod 2 eq 0) then
        hi := b div 2;
    else
        hi := b;
    end if;

    // Compute the reverse characteristic polynomial and 
    // reduce modulo $p^{N_0}$ in a signed fashion.
    f  := CharacteristicPolynomial(F);
    L  := [Coefficient(Coefficient(f, b - i), 0) : i in [0..b]];
    Z  := Integers();
    LZ := [Z!ChangePrecision(L[i+1], N0) : i in [1..b+1]];
    t  := p^N0 - (p^N0 div 2);
    for i := 0 to hi do
        if (LZ[i+1] ge t) then 
            LZ[i+1] := LZ[i+1] - p^N0;
        end if;
    end for;

    // Use the functional equation.
    if (n mod 2 eq 0) then 
        q := p^a;
        for i := hi + 1 to b do
            LZ[i+1] := LZ[b-i+1] * q^( (n - 1) * i - ((n - 1) * b) div 2 );
        end for;
    end if;

    return PolynomialRing(Z)!LZ;
end function;

function _revcharpoly_surfaces(F, d, N0)

    Qq := BaseRing(Parent(F));
    p  := Prime(Qq);
    a  := Degree(Qq);
    b  := Nrows(F);

    f  := CharacteristicPolynomial(F);
    Z  := Integers();
    L  := [Z!Coefficient(Coefficient(f, b-i), 0) : i in [0..b]];

    t   := p^N0 - (p^N0 div 2);
    h02 := Binomial(d - 1, 3);
    q   := p^a;
    for i := 0 to b do
        s := q^Abs(h02 - i);
        if (h02 gt i) then
            L[i+1] := L[i+1] * s;
        elif (h02 lt i) then
            L[i+1] := L[i+1] div s;
        end if;
        
        L[i+1] := L[i+1] mod p^N0;
        if (L[i+1] ge t) then
            L[i+1] := L[i+1] - p^N0;
        end if;

        if (h02 gt i) then
            L[i+1] := L[i+1] div s;
        elif (h02 lt i) then 
            L[i+1] := L[i+1] * s;
        end if;
    end for;

    return PolynomialRing(Z)!L;
end function;

function revcharpoly(F, n, d, N0)

    Qq := BaseRing(Parent(F));
    p  := Prime(Qq);

    if (n eq 3 and p gt 2) then
        return _revcharpoly_surfaces(F, d, N0);
    else
        return _revcharpoly(F, n, N0);
    end if;
end function;
