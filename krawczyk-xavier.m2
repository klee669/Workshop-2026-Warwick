restart

needsPackage "MonodromySolver"

krawczykOperator = (r,p0,F) -> (
    b := point{{interval(-1,1), interval(-1,1)}};
    p1:= r*matrix b + matrix p0;
    evalJFp:= evaluateJacobian(F,p0);
    evalF:= evaluate(F,p0);
    evalJFp1:= evaluateJacobian(F, p1);
    invevalJF:= inverse evaluateJacobian(F, midpoint matrix p0);
    term1:= transpose matrix p0 - (invevalJF* transpose evalF);
    term2:= id_(RR^2)-invevalJF*evalJFp1;
    final:= term1 + (term2)*(r*transpose matrix b)
    )

krawczykTest = (rho, r,p0,F) -> (
    b := point{{interval(-1,1), interval(-1,1)}};
    K:= krawczykOperator(r,p0,F);
    intervalsK:= flatten entries K;
    intervalsrb:= flatten entries (matrix p0 + rho*r* matrix b); 
    all apply(2, i -> isSubset(intervalsK#i, intervalsrb#i))
    )


end


restart
load("krawczyk-xavier.m2")

declareVariable \ {x, y, z} 
varMatrix = gateMatrix{{x,y}}
f = gateMatrix{{x^2 + y^2  - 1},{x-y^2}}
F = gateSystem(varMatrix, f)
p0 = point{{interval(0.65),interval(0.75) }}

krawczykTest (7/8, 0.1, p0, F)
