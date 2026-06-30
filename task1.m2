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

krawczykSurfaceOperator = (F,p,Rt,Rn) -> (
    FJ:=evaluateJacobian(F, transpose p);
    m:=numrows FJ;
    n:=numcols FJ;
    (s,U,V):=SVD(FJ);
    Vn:=V_{0..m-1};
    Vt:=V_{m..n-1};
    Bn:=transpose matrix{apply(m, i -> interval(-1,1))};
    Bt:=transpose matrix{apply(n-m, i -> interval(-1,1))}; 

    g1:=transpose evaluate(F ,transpose(Rt*Vt*Bt +p));
    g2:= evaluateJacobian(F, transpose (Rt*Vt*Bt + p +Rn*Vn*Bn))*Vn;

    Y:= inverse(evaluateJacobian(F, transpose p)* Vn);
    Id :=id_(RR^m);
    ans:= -Y*g1+ (Id - Y*g2)*Rn*Bn
    
)


end

restart
load("task1.m2")
declareVariable \ {x, y, z}
varMatrix = gateMatrix{{x,y}}
n=numcols varMatrix

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

varMatrix = gateMatrix{{x,y,z}}
h = gateMatrix{{x^2+y^2-1}, {z}}
H = gateSystem(varMatrix, h)
p2= matrix {{0_RR},{1_RR}, {0}}
p2= matrix {{1/sqrt(3)},{1/sqrt(3)}, {1/sqrt(3)}}
K = krawczykSurfaceOperator(H,p2,Rt,Rn)


varMatrix = gateMatrix{{x,y,z}}
g = gateMatrix{{x^2+y^2+z^2-1}}
G = gateSystem(varMatrix, g)
p2= matrix {{0_RR},{0_RR}, {1}}
K = krawczykSurfaceOperator(G,p2,Rt,Rn)


JJ=evaluateJacobian(J, transpose p1)
(s,U,V)=SVD(JJ)
rn=V_{0..m-1}
rt=V_{m..n-1}

varMatrix = gateMatrix{{x,y}}
Rt=0.1
Rn=0.1
Bn=matrix{apply(m, i -> interval(-1,1))}

K = krawczykSurfaceOperator(H,p2,Rt,Rn)

krawczykSurfaceTest = (J,p1,Rt, Rn) -> (
    
    K = krawczykSurfaceOperator(J,p1,Rt,Rn)
    fK = flatten entries K
    fKO = flatten entries(Rn*rho);
    all apply(numrows Bn, i -> isSubset(fK#i, fKO#i))
)
    