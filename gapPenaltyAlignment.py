from BLOSUM62 import *
from DNAfull import *


def fill_first_column(matrix, value):
    for i in range(len(matrix)):
        matrix[i][0] = value
    return matrix


def fill_first_row(matrix, value):
    for j in range(len(matrix[0])):
        matrix[0][j] = value
    return matrix


def create_matrix(dim1, dim2):
    matrix = [[0 for x in range(dim1)] for y in range(dim2)]
    return matrix


def score(character1, character2, is_dna):
    if is_dna:
        return DNAfull[letters_map_DNA[character1]][letters_map_DNA[character2]]
    else:
        return blossum62[letters_map_proteins[character1]][letters_map_proteins[character2]]


def find_first_move(seq1, seq2, m, i_x, i_y):
    tam_seq1 = len(seq1)
    tam_seq2 = len(seq2)
    first_move = []
    score_max = max(m[tam_seq1][tam_seq2], i_x[tam_seq1][tam_seq2], i_y[tam_seq1][tam_seq2])
    if score_max == m[tam_seq1][tam_seq2]:
        first_move.append(['DIAGONAL', tam_seq1, tam_seq2])
    if score_max == i_x[tam_seq1][tam_seq2]:
        first_move.append(['UP', tam_seq1, tam_seq2])
    if score_max == i_y[tam_seq1][tam_seq2]:
        first_move.append(['LEFT', tam_seq1, tam_seq2])
    return first_move
        

def find_next_move(seq1, seq2, first_move, m, i_x, i_y, h, g, is_dna):
    for move in first_move:
        if move[1] == 0 and move[2] == 0:
            return None
        elif move[1] == 0 and move[2] != 0:
                return [['LEFT', 0, move[2]-1]]
        elif move[2] == 0 and move[1] != 0:
                return [['UP', move[1]-1, 0]]
    
    next_move = []
    
    for move in first_move:
        if move[0] == 'DIAGONAL':
            if m[move[1]][move[2]] == m[move[1]-1][move[2]-1]+score(seq1[move[1]-1], seq2[move[2]-1], is_dna):
                next_move.append(['DIAGONAL', move[1]-1, move[2]-1])
            if m[move[1]][move[2]] == i_x[move[1]-1][move[2]-1]+score(seq1[move[1]-1], seq2[move[2]-1], is_dna):
                next_move.append(['UP', move[1]-1, move[2]-1])
            if m[move[1]][move[2]] == i_y[move[1]-1][move[2]-1]+score(seq1[move[1]-1], seq2[move[2]-1], is_dna):
                next_move.append(['LEFT', move[1]-1, move[2]-1])
        if move[0] == 'LEFT':
            if i_y[move[1]][move[2]] == m[move[1]][move[2]-1]+h+g:
                next_move.append(['DIAGONAL', move[1], move[2]-1])
            if i_y[move[1]][move[2]] == i_y[move[1]][move[2]-1]+g:
                next_move.append(['LEFT', move[1], move[2]-1])
        if move[0] == 'UP':
            if i_x[move[1]][move[2]] == m[move[1]-1][move[2]]+h+g:
                next_move.append(['DIAGONAL', move[1]-1, move[2]])
            if i_x[move[1]][move[2]] == i_x[move[1]-1][move[2]]+g:
                next_move.append(['UP', move[1]-1, move[2]])

    return next_move


def build_result(all_paths, seq1, seq2, is_dna):
    result = [[0] * 3 for i in range(len(all_paths))]
    counter = 0

    for p in all_paths:
        s1 = []
        s2 = []
        alignment = []

        for i in range(len(p)-1):
            if p[i][0] == 'DIAGONAL':
                s1.append(seq1[p[i][1]-1])
                s2.append(seq2[p[i][2]-1])
                if seq2[p[i][2] - 1] == seq1[p[i][1] - 1]:
                    alignment.insert(0, '|')
                else:
                    if score(seq1[p[i][1] - 1], seq2[p[i][2] - 1], is_dna) >= 1:
                        alignment.insert(0, ':')
                    else:
                        alignment.insert(0, '.')
            elif p[i][0] == 'UP':
                s1.append(seq1[p[i][1]-1])
                s2.append('-')
                alignment.insert(0, ' ')
            elif p[i][0] == 'LEFT':
                s1.append('-')
                s2.append(seq2[p[i][2]-1])
                alignment.insert(0, ' ')

        result[counter][0] = ''.join(s1[::-1])
        result[counter][1] = ''.join(alignment)
        result[counter][2] = ''.join(s2[::-1])
        counter += 1

    return result


def path(all_paths, current_path, first_moves, seq1, seq2, m, i_x, i_y, h, g, is_dna):
    for f in first_moves:
        if f[1] == 0 and f[2] == 0:
            p = list(current_path)
            p.append(f)
            all_paths.append(current_path + [f])
            return
        next_moves = find_next_move(seq1, seq2, [f], m, i_x, i_y, h, g, is_dna)
        copy = list(current_path)
        copy.append(f)
        for n in next_moves:
            path(all_paths, copy, [n], seq1, seq2, m, i_x, i_y, h, g, is_dna)


def gap_penalty_align(seq1, seq2, h, g, is_dna):
    seq1 = seq1.split('\n')[1::]
    seq2 = seq2.split('\n')[1::]

    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)

    tam_seq1 = len(seq1)
    tam_seq2 = len(seq2)
    minus_infinity = float("-inf")

    m = create_matrix(tam_seq2+1, tam_seq1+1)
    i_x = create_matrix(tam_seq2+1, tam_seq1+1)
    i_y = create_matrix(tam_seq2+1, tam_seq1+1)

    m = fill_first_column(m, minus_infinity)
    m = fill_first_row(m, minus_infinity)
    m[0][0] = 0

    i_x = fill_first_row(i_x, minus_infinity)
    for i in range(1, tam_seq1+1):
        i_x[i][0] = h+g*i

    i_y = fill_first_column(i_y, minus_infinity)
    for j in range(1, tam_seq2+1):
        i_y[0][j] = h+g*j

    for i in range(1, tam_seq1+1):
        for j in range(1, tam_seq2+1):
            m[i][j] = max(m[i-1][j-1], i_x[i-1][j-1], i_y[i-1][j-1])+score(seq1[i-1], seq2[j-1], is_dna)
            i_x[i][j] = max(m[i-1][j]+h+g, i_x[i-1][j]+g)
            i_y[i][j] = max(m[i][j-1]+h+g, i_y[i][j-1]+g)

    first_move = find_first_move(seq1, seq2, m, i_x, i_y)
    all_paths = []
  
    path(all_paths, [], first_move, seq1, seq2, m, i_x, i_y, h, g, is_dna)
    return build_result(all_paths, seq1, seq2, is_dna)
