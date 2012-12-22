/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

// This file contains code to test the general functionality of the code to 
// compute the Gauss-Manin connection matrix.  There should be no assumption 
// made on the base rings etc.
// 

load "GMConnection.m";

///////////////////////////////////////////////////////////////////////////////
// Test setup
// 

print "============================";
print "Compute GM Connection Matrix";
print "============================";

R<t>     := PolynomialRing(Rationals());
//R<t>     := PolynomialRing(pAdicField(7, 4));
S<X,Y,Z> := PolynomialRing(R, 3);

// P := (1-t) * (X^4 + Y^4 + Z^4) + t * (X^3 * Y + Y^3 * Z + Z^3 * X);
// P := (1-t) * (X^4 + Y^4 + Z^4) + t * (X^4+Y^4+Z^4+(X+2*Y+Z)^2*(X+3*Y+4*Z)^2);
// P := (1-t) * (X^5 + Y^5 + Z^5) + t * (X^5+Y^5+Z^5+(X+2*Y+Z)^2*(X+3*Y+4*Z)^3);
P := (1-t) * (X^3 + Y^3 + Z^3) + t * (X^3 + Y^3 + Z^3 + (X+2*Y+Z)*(X+3*Y+4*Z)^2);

time M, r := GMConnection(P);

print "=====";
print "Done.";
print "=====";

// Output
// 
print "Matrix M =", M;
print "Polynomial r =", r;

// 
// Check that r is square-free
//

dr := Derivative(r);
g := GreatestCommonDivisor(r, dr);

// Output
// 
print "Check that r is square-free:";
print "GCD(r, dr/dt) =", g;

