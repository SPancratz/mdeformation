/******************************************************************************

    Copyright (C) 2013 Sebastian Pancratz, Jan Tuitman

******************************************************************************/

load "diagfrob.m";

Z := Integers();

nRuns := 10;
nFail := 0;

Lp := [2,3,5,7,11,13,17];
Ln := [2,3,4,5,6,7];
Ld := [3,4,5];

for run := 1 to nRuns do

    p := Random(Lp);
    n := Random(Ln);
    d := Random(Ld);
    while d mod p eq 0 do
        d := Random(Ld);
    end while;

    Zmodp := ResidueClassRing(p);
    A:=[];
    for i := 1 to n+1 do
      A[i] := Zmodp!Random(1,p-1);
    end for;

    delta := frobBound(n, p);

    b := BasisSize(n, d);
    numF := ZeroMatrix(Z, b, b);
    denF := Z!0;
    diagfrob(~numF, ~denF, A, n, d, p, 50);

    if (denF lt -delta) then
        print "Error.  ord_p(F) is not within bounds.";
        print p, n, d;
        nFail := nFail + 1;
    end if;
end for;

print "Failures:", nFail;

