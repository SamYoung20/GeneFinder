# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Samantha Young

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        Addion to double check if reverse choixes works

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    """
    # TODO: implement this
    # pass
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    else:
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        No need for additions because it only checks one thing, there are no special cases.

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    dna_reversed = dna[:: -1]
    dna_revLis = list(dna_reversed)
    for i in range(len(dna_revLis)):
        dna_revLis[i] = get_complement(dna_revLis[i])
    return "".join(dna_revLis)

    # TODO: implement this
    # pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        Addition because there was no "no in frame stop codon" case

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGA")
    'ATGAGA'
    """
    dna_str = ""
#    print(index)
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if (codon == 'TAG') or (codon == 'TAA') or (codon == 'TGA'):
            index = i
            break
        else:
            dna_str = dna_str + codon
    return(dna_str)

    # Stop codons are TAG TAA TGA
    # TODO: implement this
    # pass


def find_all_ORFs_oneframe(dna):
    """Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        Addition: There is no nested DNA checks
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATATGCATTGATTTATGCCCTAG")
    ['ATGCATATGCAT', 'ATGCCC']
    """
    orfs = []
    i = 0
    while (i+2)< len(dna):
        codon = dna[i:i+3]
        if (codon == 'ATG'):
            orfs.append(rest_of_ORF(dna[i:]))
            addStrand = len(rest_of_ORF(dna[i:]))
            i = i + addStrand
        i= i+3
    return(orfs)

    # TODO: implement this
    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        additions: Checks in case there are no start codon at the beginning

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("TTTTAGATGCCCTAG")
    ['ATGCCC']
    """
    all_ORF = []
    frame1 = find_all_ORFs_oneframe(dna)
    frame2 = find_all_ORFs_oneframe(dna[1:])
    frame3 = find_all_ORFs_oneframe(dna[2:])
    all_ORF.extend(frame1)
    all_ORF.extend(frame2)
    all_ORF.extend(frame3)
    return(all_ORF)

    # TODO: implement this
    # pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        Additions: No need because this just puts together all the other tests
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    bothStrands = []
    bothStrands.extend(find_all_ORFs(dna))
    bothStrands.extend(find_all_ORFs(get_reverse_complement(dna)))
    return(bothStrands)
    # TODO: implement this
    # pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
#   get_reverse_complement("ATGC")
    doctest.testmod()
#   doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)

#   doctest.run_docstring_examples(get_reverse_complement, globals(), verbose=True)
