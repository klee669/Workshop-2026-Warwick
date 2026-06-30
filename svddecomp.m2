restart
needsPackage "MonodromySolver"


declareVariable \ {x, y, z} 
varMatrix = gateMatrix{{x,y}}
f = gateMatrix{{x^2 + y^2  - 1}}
F = gateSystem(varMatrix, f)
p0 =point {{0,1_RR}}

J =  evaluateJacobian(F, p0)
svdJ = SVD J
svdJ#2
