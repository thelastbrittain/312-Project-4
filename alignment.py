def align(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        banded_width=-1,
        gap='-'
) -> tuple[float, str | None, str | None]:
    """
        Align seq1 against seq2 using Needleman-Wunsch
        Put seq1 on left (j) and seq2 on top (i)
        => matrix[i][j]
        :param seq1: the first sequence to align; should be on the "left" of the matrix
        :param seq2: the second sequence to align; should be on the "top" of the matrix
        :param match_award: how many points to award a match
        :param indel_penalty: how many points to award a gap in either sequence
        :param sub_penalty: how many points to award a substitution
        :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
        :param gap: the character to use to represent gaps in the alignment strings
        :return: alignment cost, alignment 1, alignment 2
    """
    if banded_width == -1:
        return alignUnrestricted(seq1, seq2, match_award, indel_penalty, sub_penalty, gap)
    else:
        return alignBanded(seq1, seq2, match_award, indel_penalty, sub_penalty, banded_width, gap)



def alignUnrestricted(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        gap='-'
) -> tuple[float, str | None, str | None]:
    table: dict[tuple[int, int], tuple[float, tuple[int,int] | None]] = dict() # key = location, value = (cost, location). Location = (i, j)
    list_seq1 = list(seq1)
    list_seq2 = list(seq2)

    # fill columns
    table[(0, 0)] = (0, None)
    for i in range(1, len(seq1) + 1):
        table[(i, 0)] = (i * indel_penalty, (i - 1, 0))
    # fill rows
    for j in range(1, len(seq2) + 1):
        table[(0,j)] = (j * indel_penalty, (0, j -1))

    # fill in whole table
    for i in range(1,len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            table[(i, j)] = calcCostAndPrev(i,j, table, match_award, indel_penalty, sub_penalty, list_seq1, list_seq2)
    
    final_cell: tuple[float, tuple[int,int] | None] = table[(len(seq1), len(seq2))]
    final_seq1, final_seq2 = findFinalStrings(list_seq1, list_seq2, table)
    final_cost: float = final_cell[0]
    return (final_cost, final_seq1, final_seq2)


def alignBanded(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        banded_width=-1,
        gap='-'
) -> tuple[float, str | None, str | None]:
    # make sure banded width is less than length of sequences
    
    table: dict[tuple[int, int], tuple[float, tuple[int,int] | None]] = dict() # key = location, value = (cost, location). Location = (i, j)
    list_seq1 = list(seq1)
    list_seq2 = list(seq2)

    # fill columns ## banded width 
    table[(0, 0)] = (0, None)
    for i in range(1, banded_width + 1):
        table[(i, 0)] = (i * indel_penalty, (i - 1, 0))

    # fill rows  ## banded width
    for j in range(1, banded_width + 1):
        table[(0,j)] = (j * indel_penalty, (0, j -1))

    # fill in whole table
    for i in range(1,len(seq1) + 1):
        for j in range(bandedMinimum(i, banded_width), bandedMaximum(i, banded_width, len(seq2))):  # change for banded width
            table[(i, j)] = calcCostAndPrev(i,j, table, match_award, indel_penalty, sub_penalty, list_seq1, list_seq2)
    
    final_cell: tuple[float, tuple[int,int] | None] = table[(len(seq1), len(seq2))]
    final_seq1, final_seq2 = findFinalStrings(list_seq1, list_seq2, table)
    final_cost: float = final_cell[0]
    return (final_cost, final_seq1, final_seq2)

def bandedMinimum(i: int, banded_width: int) -> int:
    min = i - banded_width
    if min < 1:
        return 1
    else:
        return min

def bandedMaximum(i: int, banded_width: int, length_seq_2: int) -> int: # finish this
    max = i + banded_width + 1
    if max > length_seq_2 + 1:
        return length_seq_2 + 1
    else: 
        return max


def findFinalStrings(list_seq1: list[str], list_seq2: list[str], table: dict[tuple[int, int], tuple[float, tuple[int,int] | None]]) -> tuple[str, str]:
    # get final cell
    final_sequence1: list[str] = []
    final_sequence2: list[str] = []
    # final_cell: tuple[float, tuple[int,int] | None] = table[(len(list_seq1), len(list_seq2))]

    current_location: tuple[int,int] = (len(list_seq1), len(list_seq2))
    prev_location: tuple[int,int] | None = table[current_location][1]
    while prev_location is not None:

        # if up, keep left (seq1), add dash to top (seq2)
        if prev_location[1] == current_location[1]: # if we didn't move columns (we must have moved rows, so we came from up)
            final_sequence1.append(list_seq1[current_location[0] - 1])
            final_sequence2.append("-")
        elif prev_location[0] == current_location[0]: # if we didn't move rows, we must have moved columns, so we came from left
             # if left, add dash to left (seq1), keep top (seq2)
            final_sequence2.append(list_seq2[current_location[1] - 1])
            final_sequence1.append("-")
        else:   # if diaganol, then no changes
            final_sequence1.append(list_seq1[current_location[0] - 1])
            final_sequence2.append(list_seq2[current_location[1] - 1])

        # go to next one
        current_location = prev_location
        prev_location = table[current_location][1] 
    
    final1 = ''.join(final_sequence1[::-1])
    final2 = ''.join(final_sequence2[::-1])
    print(final1)
    print(final2)

    return final1, final2


def calcCostAndPrev(i: int, j: int, table: dict[tuple[int, int], tuple[float, tuple[int,int] | None]], match_award: int,
        indel_penalty: int,
        sub_penalty: int, list_seq1: list[str], list_seq2: list[str]) -> tuple[float, tuple[int, int]]:
    
    costLeft = float('inf')
    costUp = float('inf')
    costDiagonal =  float('inf')
    left = (i, j-1)
    up = (i-1, j)
    diaganol = (i-1, j-1)

    # calculate costs
    if left in table.keys(): # O(1)
        costLeft = indel_penalty + table[left][0] # O(1)
    if up in table.keys():
        costUp = indel_penalty + table[up][0]    
    costDiagonal = diff(i, j, match_award, sub_penalty, table, list_seq1, list_seq2) +  table[diaganol][0] # O(1)
    
    # return lowest cost
    if costDiagonal <= costUp and costDiagonal <= costLeft:
        return (costDiagonal, diaganol)
    elif costLeft <= costUp and costLeft <= costLeft:
        return (costLeft, left)
    else:
        return (costUp, up)


def diff(i: int, j: int, match_award: int, sub_penalty: int, table: dict[tuple[int, int], tuple[float, tuple[int,int] | None]], list_seq1: list[str], list_seq2: list[str] ) -> int:
    # if values are the same, return match, else return subsitute
    if list_seq1[i - 1] == list_seq2[j - 1]:
        return match_award
    else:
        return sub_penalty