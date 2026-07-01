restart
newtonIterate = (F, JF, p0, Vn) -> (
    Y := inverse (JF*Vn);
    pt := transpose matrix p0 -  Vn*Y*transpose(evaluate(F, p0));
    point(pt)
    )

newtonRefinement = (F,p) -> (
    episilon := 10^-5;
    JF := evaluateJacobian(F, p);
    (m,n) := (numRows JF, numColumns JF);
    (Vn, Vt) = tangentNormalFrame(JF, m, n);
    while norm evaluate(F, p) > episilon do (
        p = newtonIterate(F, JF, p, Vn);
        );
    p
    )

end

restart
load "task1.m2"
declareVariable \ {x,y,z}

--circle in 2D

varMatrix1 = gateMatrix{{x,y}}
f1 = gateMatrix{{x^2 + y^2  - 1}}
F1 = gateSystem(varMatrix1, f1)

p0 = point{{0,2_RR}} --point not on circle
JF1 = evaluateJacobian(F1, p0)
(Vn, Vt) = tangentNormalFrame(JF1, 1, 2)

newtonRefinement(F1, p0)

--prototype newtonRefinement, iterating by hand (for circle in 2D)
Y = inverse (JF1*Vn)
p1 = transpose matrix p0 -  Vn*Y*evaluate(F1, p0)
p0 = point p1
JF1 = evaluateJacobian(F1, p0)

---sphere in 3D---
varMatrix2 = gateMatrix{{x,y,z}}
f2 = gateMatrix{{x^2 + y^2 + z^2  - 1}}
F2 = gateSystem(varMatrix2, f2)

p01 = point{{1/sqrt(3)+1,1/sqrt(3)+1,1/sqrt(3)+2}}
JF2 = evaluateJacobian(F2, p01)
(Vn, Vt) = tangentNormalFrame(JF2, 1, 3)

)

---prototype newtonRefinement
Y = inverse (JF2*Vn)
p11 = transpose matrix p01 -  Vn*Y*transpose(evaluate(F2, p01))
p01 = point p11
JF2 = evaluateJacobian(F2, p01)

---circle in 3D---

varMatrix3 = gateMatrix{{x,y,z}}
f3 = gateMatrix{{x^2 + y^2  - 1},{z}}
F3 = gateSystem(varMatrix3, f3)
p03 = point{{0,2_RR,0}}

JF3 = evaluateJacobian(F3, p03)
(Vn, Vt) = tangentNormalFrame(JF3, 2, 3)

Y = inverse (JF3*Vn)
p13 = transpose matrix p03 -  Vn*Y*transpose(evaluate(F3, p03))
p03 = point p13
JF3 = evaluateJacobian(F3, p03)


