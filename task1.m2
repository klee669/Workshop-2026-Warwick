restart
needsPackage "MonodromySolver"
krawczykOperator = (F,p,r) -> (
    EJ:=evaluateJacobian(F, transpose midpoint p);
    n:=numcols EJ;
    EF:=transpose evaluate(F,transpose p);
    IEJ:=inverse EJ;
   

    B:=transpose matrix{{interval(-1, 1),interval(-1,1)}};
    EJrB:=evaluateJacobian(F, transpose(p+ r*B));
    Id:=id_(RR^n);

    p - IEJ*EF + (Id - IEJ* EJrB) * r*B
)


krawczykTest = (F,rho,p,r) ->  (
    K := krawczykOperator(F,p, r);
    fK := flatten entries K;
    B:=transpose matrix{{interval(-1, 1),interval(-1,1)}};
    fKO := flatten entries(r*rho*B+p);
    all apply(numrows B, i -> isSubset(fK#i, fKO#i))


)



tangentNormalFrame = (FJ,m,n) -> (
    M := transpose((SVD FJ)#2);
    (M_{0..m-1}, M_{m..n-1})
    )


krawczykSurfaceOperator = (F,p,Rt,Rn) -> (
    FJ:=evaluateJacobian(F, transpose p);
    m:=numrows FJ;
    n:=numcols FJ;
    (Vn, Vt):=tangentNormalFrame(FJ, m, n);
    
    Bn:=transpose matrix{apply(m, i -> interval(-1,1))};
    Bt:=transpose matrix{apply(n-m, i -> interval(-1,1))}; 

    g1:=transpose evaluate(F ,transpose(Rt*Vt*Bt +p));
    g2:= evaluateJacobian(F, transpose (Rt*Vt*Bt + p +Rn*Vn*Bn))*Vn;

    Y:= inverse(evaluateJacobian(F, transpose p)* Vn);
    Id :=id_(RR^m);
    K := -Y*g1+ (Id - Y*g2)*Rn*Bn;
    (K, maxNorm K)
    )

maxNorm = I -> (
    max apply(flatten entries I, i -> right abs i)
    )

krawczykSurfaceTest = (F,p,Rt,Rn,rho) -> (
    (K, normK) := krawczykSurfaceOperator(F, p, Rt, Rn);
    normK < Rn*rho
    )

krawczykSurfacePatch = (F,p,Rt,Rn) -> (
    FJ:=evaluateJacobian(F, transpose p);
    m:=numrows FJ;
    n:=numcols FJ;
    (Vn, Vt):=tangentNormalFrame(FJ, m, n);

    Bt:=transpose matrix{apply(n-m, i -> interval(-1,1))};
    tangentPart := Rt*Vt*Bt;

    (Knormal, normK) := krawczykSurfaceOperator(F, p, Rt, Rn);
    Kambient := p + tangentPart + Vn*Knormal;

    (Knormal, Kambient, normK)
    )

krawczykSurfacePatchTest = (F,p,Rt,Rn,rho) -> (
    (Knormal, Kambient, normK) := krawczykSurfacePatch(F, p, Rt, Rn);
    (normK < Rn*rho, Knormal, Kambient, normK)
    )




newtonIterate = (F, JF, p0, Vn) -> (
    Y := inverse (JF*Vn);
    pt := p0 -  Vn*Y*transpose(evaluate(F, transpose p0));
    pt
    )

newtonRefinement = (F,p,episilon) -> (
--    episilon := 10^-5;
    JF := evaluateJacobian(F, transpose p);
    (m,n) := (numRows JF, numColumns JF);
    (Vn, Vt) = tangentNormalFrame(JF, m, n);
    while norm evaluate(F, transpose p) > episilon do (
        p = newtonIterate(F, JF, p, Vn);
        );
    p
    )


growRadius = (F,p,Rt,Rn,rho, increasingFactor) -> (
    while krawczykSurfaceTest(F, p, Rt, Rn, 7/8) do (
        Rt = increasingFactor*Rt;
        Rn = increasingFactor*Rn;
    );
    (Rt, Rn)
)


decreaseRadius = (F,p,Rt,Rn,rho, decreasedFactor) -> (
    while not krawczykSurfaceTest(F,p,Rt,Rn,rho) do (
        Rt = Rt*decreasedFactor;
        Rn = Rn*decreasedFactor;
    );
    (Rt,Rn)
)


refinePoint = method(Options => {
        epsilon => 10^-5,
        growthFactor => 1.25,
        decreaseFactor => .5,
        rho => 7/8,
    })
refinePoint(GateSystem, Matrix, Number, Number) := o -> (F, p, Rt, Rn) -> (
    p = newtonRefinement(F, p, o.epsilon);
    (Rt, Rn) = growRadius(F, p, Rt, Rn, Vn, Vt, o.rho, o.growthFactor);
    (Rt, Rn) = decreaseRadius(F, p, Rt, Rn, o.rho, o.decreaseFactor);
    
    (p, Rt, Rn)
    )
refinePoint(GateSystem, Point, Number, Number) := o -> (F, p, Rt, Rn) -> refinePoint(F, transpose matrix p, Rt, Rn)

end

restart
load("task1.m2")
declareVariable \ {x, y, z}
varMatrix = gateMatrix{{x,y}}

rho=7/8
r=0.1
f = gateMatrix{{x^2+y^2-1},{x-y^2}}
F = gateSystem(varMatrix, f)
p = midpoint matrix{{interval(.5,.8)},{interval(.6,.9)}}
ans1=krawczykTest(F, rho,p,r)


Rt = .1
Rn = .1


varMatrix = gateMatrix{{x,y}}
j = gateMatrix{{x^2+y^2-1}}
J = gateSystem(varMatrix, j)
p1= matrix {{0_RR},{2_RR}}
refinePoint(J,p1,Rt,Rn)

K = krawczykSurfaceOperator(J,p1,Rt,Rn)
krawczykSurfaceTest(J, p1, Rt, Rn, 7/8)

varMatrix = gateMatrix{{x,y,z}}
h = gateMatrix{{x^2+y^2-1}, {z}}
H = gateSystem(varMatrix, h)
p2= matrix {{0_RR},{1_RR}, {0}}
p2= matrix {{1/sqrt(2)},{1/sqrt(2)+1}, {0}}
K = krawczykSurfaceOperator(H,p2,Rt,Rn)
krawczykSurfaceTest(H, p2, .05, Rn, 7/8)
refinePoint(H,p2,Rt,Rn)


varMatrix = gateMatrix{{x,y,z}}
g = gateMatrix{{x^2+y^2+z^2-1}}
G = gateSystem(varMatrix, g)
p2= matrix {{0_RR},{0_RR}, {1}}
K = krawczykSurfaceOperator(G,p2,Rt,Rn)
krawczykSurfaceTest(G, p2, .05, Rn, 7/8)
