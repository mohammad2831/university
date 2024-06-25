def gaussianElimination(mat, N):  
    forwardElim(mat, N)
    solve(mat, N)

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
        
        for i in range(k + 1, N):
            f = mat[i][k]/mat[k][k]
            for j in range(k + 1, N + 1):
                mat[i][j] -= mat[k][j]*f
            mat[i][k] = 0

    for i in range(N):
        x=0
        for j in range(N):
            if mat[i][j] !=0:
                x = x+1
        if x==0:
            print("many solution")
            return 0           
        
        


def solve(mat, N):
    x = [0] * N
    for i in range(N-1, -1, -1):
        x[i] = mat[i][N]
        for j in range(i + 1, N):
            x[i] -= mat[i][j] * x[j]
        x[i] /= mat[i][i]
    
    print("\nSolution for the system:")
    for val in x:
        print("{:.8f}".format(val))

N = int(input("Enter the size of N: "))
mat = []
print("Enter the matrix values row-wise:")
for _ in range(N):
    row = list(map(float, input().split()))
    mat.append(row)

gaussianElimination(mat, N)
