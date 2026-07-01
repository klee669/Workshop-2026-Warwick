restart 

needsPackage "MonodromySolver"

growRadius = (F,p,Rt,Rn,rho) -> (
    while krawczykSurfaceTest(F, p, Rt, Rn, 7/8) do (
        Rt = 2*Rt;
        print Rt;
        Rn = 2*Rn;
    );
    (Rt, Rn)
)
end
restart
load("task1.m2")
load("Radius.m2")
declareVariable \ {x, y, z}
Rt = .05
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
krawczykSurfacePatch(H,p2,0.05,0.01)
K = krawczykSurfaceOperator(H,p2,0.05,0.01)
krawczykSurfaceTest(H, p2, 0.05, 0.1, 7/8)

varMatrix = gateMatrix{{x,y,z}}
g = gateMatrix{{x^2+y^2+z^2-1}}
G = gateSystem(varMatrix, g)
p2= matrix {{0_RR},{0_RR}, {1}}
K = krawczykSurfaceOperator(G,p2,Rt,Rn)
krawczykSurfaceTest(G, p2, .05, Rn, 7/8)

rmax= 10^8

Rt = .05
Rn = .1

growRadius(G,p2, Rt, Rn, 7/8)
