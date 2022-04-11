userInput = input("> ")
row = 0
matrix = []


def readIn(userInput):
    arr = userInput.split(' ')
    user_file = arr[len(arr) - 1]

    transport = open(user_file)
    content = transport.readlines()
    global row
    row = int(content[0])
    i = 1
    while i <= row:
        matrix.append(content[i].strip().split())
        i += 1

    global output
    output = content[row + 1].strip().split()

    for i in range(0, row):
        for j in range(0, row):
            matrix[i][j] = float(matrix[i][j])

    for i in range(0, row):
        output[i] = float(output[i])

    if len(arr) == 3:
        SPPGaussian(matrix, output)
    else:
        naiveGaussian(matrix, output)




#NGE
def fwdElimination(coeff, const):
    for k in range(0, row - 1):
        for i in range(k + 1, row):
            mult = coeff[i][k] / coeff[k][k]
            for j in range(k, row):
                coeff[i][j] = coeff[i][j] - mult * coeff[k][j]
            const[i] = const[i] - mult * const[k]

def backSub(coeff, const, sol):
    sol[row - 1] = const[row - 1] / coeff[row - 1][row - 1]
    for i in range(row - 1, -1, -1):
        sum = const[i]
        for j in range(i + 1, row):
            sum = sum - coeff[i][j] * sol[j]
        sol[i] = sum / coeff[i][i]

def naiveGaussian(coeff, const):
    sol = []
    sol = [0 for i in range(row)]
    fwdElimination(coeff, const)
    backSub(coeff, const, sol)
    f = open("sys1.sol", "w")
    for x in sol:
        f.write(str(x) + " ")
    f.close()




#SPP
def SPPFwdElim(coeff, const, ind):
    scaling = [0] * row
    for i in range(0, row):
        smax = 0

        for j in range(0, row):
            smax = max(smax, abs(coeff[i][j]))

        scaling[i] = smax

    for k in range(0, row - 1): #check
        rmax = 0
        maxInd = k

        for i in range(k, row):
            r = abs(coeff[ind[i]][k] / scaling[ind[i]])
            if r > rmax:
                rmax = r
                maxInd = i

        ind[maxInd], ind[k] = ind[k], ind[maxInd] #swap

        for i in range(k + 1, row): #check
            mult = coeff[ind[i]][k] / coeff[ind[k]][k]

            for j in range(k, row): #check
                coeff[ind[i]][j] = coeff[ind[i]][j] - mult * coeff[ind[k]][j]

            const[ind[i]] = const[ind[i]] - mult * const[ind[k]]

def SPPBackSub(coeff, const, sol, ind):
    sol[row - 1] = const[ind[row - 1]] / coeff[ind[row - 1]][row - 1]
    for i in range(row - 1, - 1, - 1): #check
        sum = const[ind[i]]
        for j in range(i + 1, row): #check
            sum = sum - coeff[ind[i]][j] * sol[j]

        sol[i] = sum / coeff[ind[i]][i]

def SPPGaussian(coeff, const):
    sol = [0] * row
    ind = [0] * row

    for i in range(0, row):
        ind[i] = i

    SPPFwdElim(coeff, const, ind)
    SPPBackSub(coeff, const, sol, ind)
    f = open("sys1.sol", "w")
    for x in sol:
        f.write(str(x) + " ")
    f.close()




#gaussian sys1.lin
#gaussian --spp sys1.lin
readIn(userInput)


