/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

load "Diagfrob.m";

Z := Integers();
n := 4;             
d := 5;
A := [1,1,1,1,1];
b := BasisSize(n, d);
numF := ZeroMatrix(Z, b, b);
denF := Z!0;
diagfrob(~numF, ~denF, A, n, d, 2, 50);


