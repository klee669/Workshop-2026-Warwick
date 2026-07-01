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
    M := (SVD FJ)#2;
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
p1= matrix {{0_RR},{1_RR}}
K = krawczykSurfaceOperator(J,p1,Rt,Rn)
krawczykSurfaceTest(J, p1, Rt, Rn, 7/8)

varMatrix = gateMatrix{{x,y,z}}
h = gateMatrix{{x^2+y^2-1}, {z}}
H = gateSystem(varMatrix, h)
p2= matrix {{0_RR},{1_RR}, {0}}
p2= matrix {{1/sqrt(2)},{1/sqrt(2)}, {0}}
K = krawczykSurfaceOperator(H,p2,Rt,Rn)
krawczykSurfaceTest(H, p2, .05, Rn, 7/8)


varMatrix = gateMatrix{{x,y,z}}
g = gateMatrix{{x^2+y^2+z^2-1}}
G = gateSystem(varMatrix, g)
p2= matrix {{0_RR},{0_RR}, {1}}
K = krawczykSurfaceOperator(G,p2,Rt,Rn)
krawczykSurfaceTest(G, p2, .05, Rn, 7/8)



restart
load("taskM.m2")

declareVariable \ {x,y,z}

varMatrix = gateMatrix{{x,y}}
j = gateMatrix{{x^2+y^2-1}}
J = gateSystem(varMatrix, j)

p1 = matrix{{0_RR},{1_RR}}

Rt = 0.1
Rn = 0.1
rho = 7/8

krawczykSurfaceTest(J,p1,Rt,Rn,rho)

(RtNew,RnNew) = decreaseRadius(J,p1,Rt,Rn,rho)

RtNew
RnNew

krawczykSurfaceTest(J,p1,RtNew,RnNew,rho)



varMatrix = gateMatrix{{x,y,z}}
h = gateMatrix{{x^2+y^2-1},{z}}
H = gateSystem(varMatrix, h)

p2 = matrix{{0_RR},{1_RR},{0_RR}}

krawczykSurfaceTest(H,p2,Rt,Rn,rho)

(RtNew,RnNew) = decreaseRadius(H,p2,Rt,Rn,rho)

RtNew
RnNew

krawczykSurfaceTest(H,p2,RtNew,RnNew,rho)
K = krawczykSurfaceOperator(H,p2,RtNew,RnNew)




restart
load("taskM.m2")

decreaseRadius = (F,p,Rt,Rn,rho, decreasedFactor) -> (
    while not krawczykSurfaceTest(F,p,Rt,Rn,rho) do (
        Rt = Rt*decreasedFactor;
        Rn = Rn*decreasedFactor;
    );
    (Rt,Rn)
)

declareVariable \ {x,y,z}

rho = 7/8
Rt = 0.1
Rn = 0.1
decreaseFactor = 1/2

varMatrix = gateMatrix{{x,y,z}}
h = gateMatrix{{x^2+y^2-1}, {z}}
H = gateSystem(varMatrix, h)

p2 = matrix{{1_RR/sqrt(2_RR)}, {1_RR/sqrt(2_RR)}, {0_RR}}

krawczykSurfaceTest(H, p2, Rt, Rn, rho)

(RtH,RnH) = decreaseRadius(H, p2, Rt, Rn, rho)

RtH
RnH

krawczykSurfaceTest(H, p2, RtH, RnH, rho)

KH = krawczykSurfaceOperator(H, p2, RtH, RnH)
KH


Rt = 10
Rn = 10
decreasedFactor = 1/3

varMatrix = gateMatrix{{x,y,z}}
g = gateMatrix{{x^2+y^2+z^2-1}}
G = gateSystem(varMatrix, g)

p3 = matrix{{0_RR}, {0_RR}, {1_RR}}

krawczykSurfaceTest(G, p3, Rt, Rn, rho)

(RtG,RnG) = decreaseRadius(G, p3, Rt, Rn, rho)

RtG
RnG

krawczykSurfaceTest(G, p3, RtG, RnG, rho)

KG = krawczykSurfaceOperator(G, p3, RtG, RnG)
KG