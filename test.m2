restart

needsPackage "MonodromySolver"
declareVariable \ {x, y, z} 
varMatrix = gateMatrix{{x,y,z}}
f = gateMatrix{{x^2 + y^2 + z^2 - 1}}
F = gateSystem(varMatrix, f)

p = point{{1,1,1}}
evaluate(F,p)

jf=diff(varMatrix,f)
jF = gateSystem(varMatrix, transpose jf)
evaluate(jF,p)

interval(1,2)
I =point{{interval(1,2), interval(3,4), interval(1,2)}}
evaluate(jF, I)
interval(1+ii, 2+3*ii)
