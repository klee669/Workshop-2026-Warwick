restart
needsPackage "MonodromySolver"
load("task1.m2")
declareVariable \ {x, y, z}
Rt = .1
Rn = .1


varMatrix = gateMatrix{{x,y}}
j = gateMatrix{{x^2+y^2-1}}
J = gateSystem(varMatrix, j)
p1= matrix {{0_RR},{2_RR}}
p, Rt, Rn, Vt=refinePoint(J,p1,Rt,Rn)

FJ = evaluateJacobian(J, transpose p)
m,n = (numRows FJ, numColumns FJ)

queue={}
d= numcols Vt
for i from 0 to d-1 do (
    v = Vt_(i);
    queue = append(queue, p + Rt*v);
    queue = append(queue, p - Rt*v);
)

