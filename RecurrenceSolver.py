x, n = var('x n')
#Solves Equations of the form a_k*x_(n+k) + a_(k-1)*x_(n+k-1) + ... + a_0*x_n = p_1*e_1 + p_2*e_2 + ... + p_z*e_z
#Where a_i are complex coeifficients, x_i is the unknown sequence, p_i are polynomial function in n, and e_i are exponential functions
coeifficients = [] #Array of coeifficients in order i.e. a_k, a_(k-1), ... , a_0.
polynomialFunc = [0] #Array of polynomials in order i.e. p_1, p_2, ... p_z
#NOTE: To input a constant polynomial in n one must write c + 0*n in Sage
exponentialFunc = [0] #Array of exponentials in order i.e. e_1, e_2, ... e_z
initialConditions = [1, 1] #Intials conditions for the sequence i.e. x_0, x_1, ..., x_k
#NOTE: There MUST be k initial conditions for the sequence to be defined!
homogeneousAdvancementEquation = 0;
particularFunction(n) = 0;
for i in range(len(polynomialFunc)):
    particularFunction += polynomialFunc[i]*exponentialFunc[i]^n
for i in range(len(coeifficients)-1,-1,-1):
    homogeneousAdvancementEquation += coeifficients[len(coeifficients) - i - 1]*x^i
for i in range(len(exponentialFunc)):
    if polynomialFunc[i] != 0:
        homogeneousAdvancementEquation *= (x-exponentialFunc[i])^(polynomialFunc[i].degree(n)+1)
roots = homogeneousAdvancementEquation.roots(x)
generalSolution = [];
for i in range(len(roots)):
    for j in range(roots[i][1]):
        f(n) = n^j*roots[i][0]^n
        generalSolution.append(f(n))
while len(initialConditions) < homogeneousAdvancementEquation.degree(x):
    next = 0
    for i in range(1,len(coeifficients)):
        next -= coeifficients[i] * initialConditions[len(initialConditions) - i]
    next += particularFunction(len(initialConditions) - len(coeifficients) + 1)
    next /= coeifficients[0]
    initialConditions.append(next)
M = MatrixSpace(ComplexField(),homogeneousAdvancementEquation.degree(x),len(generalSolution) + 1)
system= []
for i in range(homogeneousAdvancementEquation.degree(x)):
    for j in range(len(generalSolution)):
        f(n) = generalSolution[j]
        system.append(f(i))
    system.append(initialConditions[i])
system = M(system)
constants = system.echelon_form().matrix_from_columns([len(generalSolution)])
finalAnswer = 0
for i in range(len(generalSolution)):
    finalAnswer += constants[i]*generalSolution[i]
print(finalAnswer)
for i in range(10):
    print(finalAnswer(i).n())
