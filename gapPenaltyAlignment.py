
def preencher_primeira_coluna(matriz,valor):
    for i in range(len(matriz)):
            matriz[i][0]=valor
    return matriz


def preencher_primeira_linha(matriz,valor):
    for j in range(len(matriz[0])):
            matriz[0][j]=valor
    return matriz


def criar_matriz(dim1,dim2):
    matriz=[[0 for x in range(dim1)] for y in range(dim2)]
    return matriz


def score(caracter1, caracter2, custo_acertar, custo_errar ):
    if caracter1==caracter2:
        return custo_acertar
    if caracter1!=caracter2:
        return custo_errar


def find_first_move(seq1,seq2,M,I_x,I_y):
    tam_seq1=len(seq1)
    tam_seq2=len(seq2)
    first_move=[]
    score_max= max(M[tam_seq1][tam_seq2],I_x[tam_seq1][tam_seq2],I_y[tam_seq1][tam_seq2])
    if score_max==M[tam_seq1][tam_seq2]:
        first_move.append(['DIAGONAL',tam_seq1,tam_seq2])
    if score_max==I_x[tam_seq1][tam_seq2]:
        first_move.append(['UP',tam_seq1,tam_seq2])
    if score_max==I_y[tam_seq1][tam_seq2]:
        first_move.append(['LEFT',tam_seq1,tam_seq2])
    return first_move
        

def find_next_move(seq1,seq2,first_move,M,I_x,I_y,custo_acertar,custo_errar,h,g):
    for m in first_move:
        if m[1] == 0 and m[2] == 0:
            return None
        elif m[1] == 0 and m[2] != 0:
                return [['LEFT',0,m[2]-1]]
        elif m[2] == 0 and m[1] != 0:
                return [['UP',m[1]-1,0]]
    
    next_move=[]
    
    for m in first_move:
        if m[0]=='DIAGONAL':
            if M[m[1]][m[2]]==M[m[1]-1][m[2]-1]+score(seq1[m[1]-1],seq2[m[2]-1],custo_acertar,custo_errar):
                next_move.append(['DIAGONAL',m[1]-1,m[2]-1])
            if M[m[1]][m[2]]==I_x[m[1]-1][m[2]-1]+score(seq1[m[1]-1],seq2[m[2]-1],custo_acertar,custo_errar):
                next_move.append(['UP',m[1]-1,m[2]-1])
            if M[m[1]][m[2]]==I_y[m[1]-1][m[2]-1]+score(seq1[m[1]-1],seq2[m[2]-1],custo_acertar,custo_errar):
                next_move.append(['LEFT',m[1]-1,m[2]-1])
        if m[0]=='LEFT':
            if I_y[m[1]][m[2]]==M[m[1]][m[2]-1]+h+g:
                next_move.append(['DIAGONAL',m[1],m[2]-1])
            if I_y[m[1]][m[2]]==I_y[m[1]][m[2]-1]+g:
                next_move.append(['LEFT',m[1],m[2]-1])
        if m[0]=='UP':
            if I_x[m[1]][m[2]]==M[m[1]-1][m[2]]+h+g:
                next_move.append(['DIAGONAL',m[1]-1,m[2]])
            if I_x[m[1]][m[2]]==I_x[m[1]-1][m[2]]+g:
                next_move.append(['UP',m[1]-1,m[2]])

    return next_move


def genes(seq1,seq2,path):
    s1=[]
    s2=[]
  
    for i in range(len(path)-1):
        if path[i][0]=='DIAGONAL':
            s1.append(seq1[path[i][1]-1])
            s2.append(seq2[path[i][2]-1])
        elif path[i][0]=='UP':
            s1.append(seq1[path[i][1]-1])
            s2.append('-')
        elif path[i][0]=='LEFT':
            s1.append('-')
            s2.append(seq2[path[i][2]-1])
        
    print('seq1:', s1[::-1])
    print('seq2:', s2[::-1])


def path(current_path,first_moves,seq1,seq2,M,I_x,I_y,custo_acertar,custo_errar,h,g):
    for f in first_moves:
        if f[1]==0 and f[2]==0:
            p=list(current_path)
            p.append(f)
            print(genes(seq1,seq2,p))
            return
        next_moves=find_next_move(seq1,seq2,[f],M,I_x,I_y,custo_acertar,custo_errar,h,g)
        copy=list(current_path)
        copy.append(f)
        for n in next_moves:
            path(copy,[n],seq1,seq2,M,I_x,I_y,custo_acertar,custo_errar,h,g)


def alinhamento_afim(seq1, seq2, h, g, custo_acertar, custo_errar):
    tam_seq1=len(seq1)
    tam_seq2=len(seq2)
    menos_infinito=float("-inf")

    M=criar_matriz(tam_seq2+1,tam_seq1+1)
    I_x=criar_matriz(tam_seq2+1,tam_seq1+1)
    I_y=criar_matriz(tam_seq2+1,tam_seq1+1)

    M=preencher_primeira_coluna(M,menos_infinito)
    M=preencher_primeira_linha(M,menos_infinito)
    M[0][0]=0

    I_x=preencher_primeira_linha(I_x,menos_infinito)
    for i in range(1,tam_seq1+1):
        I_x[i][0]=h+g*i

    I_y=preencher_primeira_coluna(I_y,menos_infinito)
    for j in range(1,tam_seq2+1):
        I_y[0][j]=h+g*j

    for i in range(1, tam_seq1+1):
        for j in range(1, tam_seq2+1):
            M[i][j] = max(M[i-1][j-1], I_x[i-1][j-1], I_y[i-1][j-1])+score(seq1[i-1],seq2[j-1],custo_acertar,custo_errar)
            I_x[i][j] = max(M[i-1][j]+h+g,I_x[i-1][j]+g)
            I_y[i][j] = max(M[i][j-1]+h+g,I_y[i][j-1]+g)
    
    
    first_move=find_first_move(seq1,seq2,M,I_x,I_y)
  
    print(path([],first_move,seq1,seq2,M,I_x,I_y,custo_acertar,custo_errar,h,g))


print(alinhamento_afim('AAT','ACACT',-3,-1,1,-1))
print(alinhamento_afim('GCGCGTTAGACTAGCACCG','GGGTTGCACCG',-5,-1,3,-2))
