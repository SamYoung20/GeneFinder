# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Samantha Young
"""

import random
from amino_acids import aa, codons, aa_table   #you may find these useful
from load import load_seq
dna1 = load_seq("./data/X73525.fa")


class DNAfinder:
    def __init__(self, dnadata):
        self.dna = dnadata

    def shuffle_string(self, s):
        """Shuffles the characters in the input string
            NOTE: this is a helper function, you do not
            have to modify this in any way """
        return ''.join(random.sample(s, len(s)))

    def get_complement(self, nucleotide):
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

        """ old code
        if nucleotide == 'A':
            return 'T'
        elif nucleotide == 'T':
            return 'A'
        elif nucleotide == 'C':
            return 'G'
        else:
            return 'C'
    """
        cods = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'} # changed to a dictionary for readability
        return cods[nucleotide]

    def get_reverse_complement(self, dna):
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
            dna_revLis[i] = self.get_complement(dna_revLis[i])
        return "".join(dna_revLis)

    def rest_of_ORF(self, dna):
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
        for i in range(0, len(dna), 3):
            codon = dna[i:i+3]
            if (codon == 'TAG') or (codon == 'TAA') or (codon == 'TGA'):
                index = i
                break
            else:
                dna_str = dna_str + codon
        return(dna_str)

        # Stop codons are TAG TAA TGA

    def find_all_ORFs_oneframe(self, dna):
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
        while (i+2) < len(dna):
            codon = dna[i:i+3]
            if (codon == 'ATG'):
                orfs.append(self.rest_of_ORF(dna[i:]))
                addStrand = len(self.rest_of_ORF(dna[i:]))
                i = i + addStrand
            i = i + 3
        return(orfs)

    def find_all_ORFs(self, dna):
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
        frame1 = self.find_all_ORFs_oneframe(dna)
        frame2 = self.find_all_ORFs_oneframe(dna[1:])
        frame3 = self.find_all_ORFs_oneframe(dna[2:])
        all_ORF.extend(frame1)
        all_ORF.extend(frame2)
        all_ORF.extend(frame3)
        return(all_ORF)

    def find_all_ORFs_both_strands(self, dna):
        """ Finds all non-nested open reading frames in the given DNA sequence on both
            strands.
            Additions: No need because this just puts together all the other tests
            dna: a DNA sequence
            returns: a list of non-nested ORFs
        >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
        ['ATGCGAATG', 'ATGCTACATTCGCAT']
        """
        bothStrands = []
        bothStrands.extend(self.find_all_ORFs(dna))
        bothStrands.extend(self.find_all_ORFs(self.get_reverse_complement(dna)))
        return(bothStrands)

    def longest_ORF(self, dna):
        """ Finds the longest ORF on both strands of the specified DNA and returns it
            as a string

        >>> longest_ORF("ATGCGAATGTAGCATCAAA")
        'ATGCTACATTCGCAT'
        """
        all_orfs = self.find_all_ORFs_both_strands(dna)
        return max(all_orfs, key=len)

    def longest_ORF_noncoding(self,dna, num_trials):
        """ Computes the maximum length of the longest ORF over num_trials shuffles
            of the specfied DNA sequence

            dna: a DNA sequence
            num_trials: the number of random shuffles
            returns: the maximum length longest ORF """
        dnaLs = list(dna)
        the_longest = []
        for i in range(0, num_trials):
            random.shuffle(dnaLs)
            dnaStr = "".join(dnaLs)
            the_longest.append(self.longest_ORF(dnaStr))
            # print(the_longest)
        # print(len(max(the_longest, key=len)))
        return max(the_longest, key=len)

    def coding_strand_to_AA(self, dna):
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
        amino_acidLs = []
        dif = len(dna) % 3
        for i in range(0, len(dna) - dif, 3):
            codonNew = dna[i:i+3]
            amino_acidLs.append(aa_table[codonNew])
        aminoStr = "".join(amino_acidLs)
        return(aminoStr)

    def gene_finder(self, dna):
        """ Returns the amino acid sequences that are likely coded by the specified dna

            dna: a DNA sequence
            returns: a list of all amino acid sequences coded by the sequence dna.
        """
        threshold = len(self.longest_ORF_noncoding(dna, 1500))
        list_o_Orfs = self.find_all_ORFs_both_strands(dna)
        final_orf = [n for n in list_o_Orfs if len(n) >= threshold]
        aminoAcidCodes = [self.coding_strand_to_AA(n) for n in final_orf]
        print(aminoAcidCodes)
        return aminoAcidCodes


if __name__ == "__main__":
    mygenefinder = DNAfinder(dna1)
    mygenefinder.gene_finder(mygenefinder.dna)
