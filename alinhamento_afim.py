def count_open_gaps(seq1,seq2):
    count=0
    for i in range(len(seq1)-1):
        if seq1[i]!='-' and seq1[i+1]=='-' and seq2[i]!='-' and seq2[i+1]!='-':
            count+=1
        elif seq2[i]!='-' and seq2[i+1]=='-' and seq1[i]!='-' and seq1[i+1]!='-':
            count+=1
        
    return count


def is_left(pos1,pos2):
    if pos2[1]== pos1[1]-1 and pos2[0]==pos1[0]:
        return True
    return False


def is_up(pos1,pos2):
    if pos2[0]==pos1[0]-1 and pos2[1]==pos1[1]:
        return True
    return False


def is_diag(pos1,pos2):
    if pos2[1]==pos1[1]-1 and pos2[0]==pos1[0]-1:
        return True
    return False


def print_genes(paths,seq2,seq1,h,g,custo_acertar,custo_errar,score_max):
    for p in paths:
        l1=[]
        l2=[]
        score=[]
        
        for i in range(0,len(p)-1):
            
            if is_left(p[i],p[i+1]):
                l1.insert(0,'-')
                l2.insert(0,seq1[p[i][1]-1])
                score.insert(0,g)

                
            elif is_up(p[i],p[i+1]):
                l1.insert(0,seq2[p[i][0]-1])
                l2.insert(0,'-')
                score.insert(0,g)
                
                
            elif is_diag(p[i],p[i+1]):
                l1.insert(0,seq2[p[i][0]-1])
                l2.insert(0,seq1[p[i][1]-1])
                if l1[0]==l2[0] and l1[0]!='-' and l2[0]!='-':
                    score.insert(0,custo_acertar)
                elif l1[0]!=l2[0]:
                    score.insert(0,custo_errar)
                
        gaps=count_open_gaps(l1,l2)
        if sum(score[:])+gaps*h==score_max:
            print "seq1->",l2
            print "seq2->",l1
            print "score max=", score_max
            print('-'*50)


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


def find_next_move(seq1, seq2, M,I_x,I_y, custo_acertar, custo_errar, i, j):
    if i == 0 and j == 0:
        return None
    elif i == 0 and j != 0:
        return ['LEFT']
    elif j == 0 and i != 0:
        return ['UP']
    
    moves=[]
    if M[i][j]==max(M[i][j],I_x[i][j],I_y[i][j]):
         moves.append('DIAGONAL')
    if I_x[i][j]==max(M[i][j],I_x[i][j],I_y[i][j]):
        moves.append('UP')
    if I_y[i][j]==max(M[i][j],I_x[i][j],I_y[i][j]):
        moves.append('LEFT')
    
    return moves


def find_path(seq1, seq2, M, I_x, I_y, custo_acertar, custo_errar, path, all_paths, i, j):
    next_moves = find_next_move(seq1, seq2, M, I_x, I_y, custo_acertar, custo_errar, i, j)

    if next_moves is None:
        all_paths.append(path + [[i, j]])
        return

    for move in next_moves:
        if move == 'UP':
            find_path(seq1, seq2, M, I_x, I_y, custo_acertar, custo_errar, path + [[i, j]], all_paths,i - 1,j)
        elif  move == 'LEFT':
            find_path(seq1, seq2, M, I_x, I_y, custo_acertar, custo_errar, path + [[i,j]], all_paths, i, j - 1)
        else:
            find_path(seq1, seq2, M, I_x, I_y, custo_acertar, custo_errar, path + [[i,j]], all_paths, i - 1, j - 1)
    return


def traceback(seq1, seq2, M, I_x, I_y, tam_seq1, tam_seq2, custo_acertar, custo_errar):
    i, j = tam_seq2, tam_seq1
    all_paths = [[]]
    find_path(seq1, seq2, M, I_x, I_y, custo_acertar, custo_errar, [], all_paths, i, j)
    return all_paths


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
            I_x[i][j] = max(M[i-1][j]+h+g,I_x[i-1][j]+g,I_y[i-1][j]+h+g)
            I_y[i][j] = max(M[i][j-1]+h+g,I_y[i][j-1]+g,I_x[i][j-1]+h+g)
    
    score_max= max(M[tam_seq1][tam_seq2],I_x[tam_seq1][tam_seq2],I_y[tam_seq1][tam_seq2])
    
    all_paths = traceback(seq1, seq2, M, I_x, I_y, tam_seq2, tam_seq1, custo_acertar, custo_errar)
    
    print_genes(all_paths,seq1,seq2,h,g,custo_acertar,custo_errar,score_max)


print(alinhamento_afim('GGGTTGCACCG','GCGCGTTAGACTAGCACCG',-5,-1,3,-2))
