def gaussJordanElimination(mat, N):
    forwardElim(mat, N)
    backSubstitute(mat, N)

def swap_row(mat, i, j):
    for k in range(N + 1):
        temp = mat[i][k]
        mat[i][k] = mat[j][k]
        mat[j][k] = temp

def forwardElim(mat, N):
    for k in range(N):    
        i_max = k
        v_max = mat[i_max][k]
        
        for i in range(k + 1, N):
            if (abs(mat[i][k]) > v_max):
                v_max = mat[i][k]
                i_max = i
        
        if (i_max != k):
            swap_row(mat, k, i_max)
        
        for i in range(N):
            if i != k:
                f = mat[i][k] / mat[k][k]
                for j in range(k, N + 1):
                    mat[i][j] -= f * mat[k][j]
                mat[i][k] = 0

def backSubstitute(mat, N):
    for i in range(N):
        if mat[i][i] == 0:
            print("Many solutions or no solution")
            return
        
        mat[i][N] /= mat[i][i]

def answer(mat, N):
    for i in range(N):
        print("x{} = {}".format(i+1, mat[i][N]))

N = int(input("Enter the size of N: "))
mat = []
print("Enter the matrix values row-wise:")
for _ in range(N):k
    row = list(map(float, input().split()))
    mat.append(row)

gaussJordanElimination(mat, N)
answer(mat, N)

