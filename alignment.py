from BLOSUM62 import *
from DNAfull import *


def is_left(pos1, pos2):
    if pos2[1] == pos1[1] - 1 and pos2[0] == pos1[0]:
        return True
    return False


def is_up(pos1, pos2):
    if pos2[0] == pos1[0] - 1 and pos2[1] == pos1[1]:
        return True
    return False


def is_diagonal(pos1, pos2):
    if pos2[1] == pos1[1] - 1 and pos2[0] == pos1[0] - 1:
        return True
    return False


def print_result(paths, seq1, seq2):
    for p in paths:
        l1 = []
        l2 = []
        alignment = []

        for i in range(0, len(p) - 1):

            if is_left(p[i], p[i + 1]):
                l1.insert(0, '-')
                l2.insert(0, seq1[p[i][1] - 1])
                alignment.insert(0, ' ')

            elif is_up(p[i], p[i + 1]):
                l1.insert(0, seq2[p[i][0] - 1])
                l2.insert(0, '-')
                alignment.insert(0, ' ')

            elif is_diagonal(p[i], p[i + 1]):
                l1.insert(0, seq2[p[i][0] - 1])
                l2.insert(0, seq1[p[i][1] - 1])
                if seq2[p[i][0] - 1] == seq1[p[i][1] - 1]:
                    alignment.insert(0, '|')
                else:
                    alignment.insert(0, '.')

        print(''.join(l2))
        print(''.join(alignment))
        print(''.join(l1))


def score(character1, character2, is_dna):
    if is_dna:
        return DNAfull[letters_map_DNA[character1]][letters_map_DNA[character2]]
    else:
        return blossum62[letters_map_proteins[character1]][letters_map_proteins[character2]]


def find_next_move(seq1, seq2, matrix_scores, gap, i, j, is_dna):
    if i == 0 and j == 0:
        return None
    elif i == 0 and j != 0:
        return ['LEFT']
    elif j == 0 and i != 0:
        return ['UP']

    actual_score = matrix_scores[i][j]
    diagonal_score = matrix_scores[i - 1][j - 1]
    left_score = matrix_scores[i][j - 1]
    up_score = matrix_scores[i - 1][j]

    moves = []
    if actual_score == up_score + gap:
        moves.append('UP')
    if actual_score == left_score + gap:
        moves.append('LEFT')
    if actual_score == diagonal_score + score(seq1[j - 1], seq2[i - 1], is_dna):
        moves.append('DIAGONAL')
    return moves


def find_path(seq1, seq2, matrix_scores, gap, path, all_paths, i, j, is_dna):
    next_moves = find_next_move(seq1, seq2, matrix_scores, gap, i, j, is_dna)

    if next_moves is None:
        all_paths.append(path + [[i, j]])
        return

    for move in next_moves:
        if move == 'UP':
            find_path(seq1, seq2, matrix_scores, gap, path + [[i, j]], all_paths, i - 1, j, is_dna)
        elif move == 'LEFT':
            find_path(seq1, seq2, matrix_scores, gap, path + [[i, j]], all_paths, i, j - 1, is_dna)
        else:
            find_path(seq1, seq2, matrix_scores, gap, path + [[i, j]], all_paths, i - 1, j - 1, is_dna)
    return


def traceback(seq1, seq2, matrix_scores, tam_seq1, tam_seq2, gap, is_dna):
    i, j = tam_seq2, tam_seq1
    all_paths = []
    find_path(seq1, seq2, matrix_scores, gap, [], all_paths, i, j, is_dna)
    return all_paths


def needleman_wunsch(seq1, seq2, gap, is_dna):
    tam_seq1, tam_seq2 = len(seq1), len(seq2)
    matrix_scores = [[0 for x in range(tam_seq1 + 1)] for y in range(tam_seq2 + 1)]

    for i in range(tam_seq1 + 1):
        matrix_scores[0][i] = gap * i
    for j in range(tam_seq2 + 1):
        matrix_scores[j][0] = gap * j
    for i in range(1, tam_seq2 + 1):
        for j in range(1, tam_seq1 + 1):
            matrix_scores[i][j] = max(
                matrix_scores[i - 1][j - 1] + score(seq1[j - 1], seq2[i - 1], is_dna),
                matrix_scores[i - 1][j] + gap, matrix_scores[i][j - 1] + gap)

    all_paths = traceback(seq1, seq2, matrix_scores, tam_seq1, tam_seq2, gap, is_dna)
    print_result(all_paths, seq1, seq2)