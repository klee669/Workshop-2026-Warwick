restart
needsPackage "MonodromySolver"

tangentNormalFrame = (F,p0) -> (
    J := evaluateJacobian(F, p0);
    (m, n) := (numrows J, numcols J);
    M := (SVD J)#2;
    (M_{0..m-1}, M_{m..n-1})
    )

end
restart
load("svddecomp.m2")

declareVariable \ {x, y, z}

---circle in 2D---
varMatrix1 = gateMatrix{{x,y}}
f1 = gateMatrix{{x^2 + y^2  - 1}}
F1 = gateSystem(varMatrix1, f1)
n = 2

p0 = point{{0,1_RR}}

tangentNormalFrame(F1,p0)

---sphere in 3D---
varMatrix2 = gateMatrix{{x,y,z}}
f2 = gateMatrix{{x^2 + y^2 + z^2  - 1}}
F2 = gateSystem(varMatrix2, f2)
n = 3
p01 = point{{0,0,1_RR}}

tangentNormalFrame(F2,p01)

---circle in 3D---
varMatrix3 = gateMatrix{{x,y,z}}
f3 = gateMatrix{{x^2 + y^2  - 1},{z}}
F3 = gateSystem(varMatrix3, f3)
p02 = point{{0,1_RR,0}}

tangentNormalFrame(F3,p02)


    

