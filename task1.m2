restart
needsPackage "MonodromySolver"

krawczykoperator = (F,p,r) -> (
    EJ:=evaluateJacobian(F, transpose midpoint p);
    n:=numcols EJ;
    EF:=transpose evaluate(F,transpose p);
    IEJ:=inverse EJ;
   

    B:=transpose matrix{{interval(-1, 1),interval(-1,1)}};
    EJrB:=evaluateJacobian(F, transpose(p+ r*B));
    Id:=id_(RR^n);

    p - IEJ*EF + (Id - IEJ* EJrB) * r*B
)


krawczyktest = (F,rho,p,r) ->  (
    K := krawczykoperator(F,p, r);
    fK := flatten entries K;
    B:=transpose matrix{{interval(-1, 1),interval(-1,1)}};
    fKO := flatten entries(r*rho*B+p);
    all apply(numrows B, i -> isSubset(fK#i, fKO#i))

)

end

restart
load("task1.m2")
declareVariable \ {x, y}
varMatrix = gateMatrix{{x,y}}


rho=7/8
r=0.1
f = gateMatrix{{x^2+y^2-1},{x-y^2}}
F = gateSystem(varMatrix, f)
p = midpoint matrix{{interval(.5,.8)},{interval(.6,.9)}}

ans=krawczyktest(F, rho,p,r)


