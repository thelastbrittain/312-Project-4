from alignment import alignUnrestricted

seq1 = 'ATGCATGC'
seq2 = 'ATGGTGC'

seq1Hard = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatctctttataaacggcacttcc'
seq2Hard = 'ataagagtgattggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'



def test_table_setup_works():
    alignUnrestricted(seq1=seq1, seq2=seq2)
    assert(True)

def test_simple_correct_cost():
    cost = alignUnrestricted(seq1, seq2)[0]
    assert cost == -12

def test_harder_cost_correct():
    score = alignUnrestricted(seq1Hard, seq2Hard)[0]
    assert score == -116