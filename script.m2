load "task1.m2"


declareVariable \ {x,y,z}

--sphere in 3D
varMatrix = gateMatrix{{x,y,z}}
f = gateMatrix{{x^2 + y^2 + z^2  - 1}}
F = gateSystem(varMatrix, f)

p0 = matrix {{1_RR},{0},{0}} --point on variety
p1 = matrix {{1/sqrt(3)},{1/sqrt(3)},{1/sqrt(3)}} --point on variety
p2 = matrix {{2_RR},{0},{0}} --point not on variety



(Vn0, Vt0) = tangentNormalFrame(evaluateJacobian(F, point p0),1,3)
krawczykSurfaceTest (F, p0, 0.1, 0.1, 7/8, Vn0, Vt0)
--Krawczyk test at (1,0,0) for 0.1^3 box, with rho = 7/8
    
(Vn1, Vt1) = tangentNormalFrame(evaluateJacobian(F, point p1),1,3)
krawczykSurfaceTest (F, p1, 1, 1, 7/8, Vn1, Vt1)
--Krawczyk test at (1/sqrt{3},1/sqrt{3},1/sqrt{3}) for 0.1^3 box, with rho = 7/8

(Vn2, Vt2) = tangentNormalFrame(evaluateJacobian(F, point p2),1,3)
krawczykSurfaceTest (F, p2, 0.1, 0.1, 7/8, Vn2, Vt2)
--Krawczyk test at (2,0,0) for 0.1^3 box, with rho = 7/8


refinePoint(F,p2,0.1,0.1)
-- moves a point onto variety via Newton method
-- resizes box containing point to optimise portion of variety it contain



