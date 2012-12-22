/******************************************************************************

    Copyright (C) 2009-2012 Sebastian Pancratz

******************************************************************************/

/*
    Computes a list [*b1 b2 ... bn*], where bk is defined to be the set 
    (given here as a list) of monomials X^i such that |i| = kd-(n+1) 
    and 0 <= i < d-1.

    We make the following assumptions:
        o S is a polynomial ring in n+1 variables.
        o d is a natural number.
 */
function BasisSets(S, d)
	n := Rank(S)-1;
	B := [* *];
	for k:=1 to n do
		// 
		// Pick up the first few cases were bk is the empty set.
		// 
		if k*d lt n+1 then
			Append(~B, [* *]);
			continue;
		end if;
		// 
		// Now in the general case..
		// 
		f := S!0;
		for i:=1 to n+1 do
			f := f + S.i;
		end for;
		Lk := Monomials(f^(k*d - (n+1)));
		Bk := [* *];
		for f in Lk do
			for i:=1 to n+1 do
				if Degree(f, i) ge d-1 then
					continue f;
				end if;
			end for;
			Append(~Bk, f);
		end for;
		Append(~B, Bk);
	end for;
	return B;
end function;

/* 
    Returns the dimension of $H_{dR}^n(U/S)$.
 */
function BasisSize(n, d)
    return ((d - 1) * ((d - 1)^n - (-1)^n)) div d;
end function;

