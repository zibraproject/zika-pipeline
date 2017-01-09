#!/usr/bin/env python
import ctypes as ct
import ssw_lib
from Bio import SeqIO
import argparse
import os,sys


def get_parser():
    parser = argparse.ArgumentParser(
        description="""A simple read demultiplexer for Oxford Nanopore data.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input", action=FileExist,
        help="Path to fasta file.")
    parser.add_argument("--barcodes",default="NB_barcodes.fasta", type=str,
        help="Relative path to fasta file describing barcodes.")
    parser.add_argument("--threshold", default=90, type=int,
        help="Minimum match score to accept called barcodes.")


    return parser

def align_seq(seq,args):
    resultdict=dict()
    for bc_name in barcode_dict:
        match,score = nucl_align(seq,barcode_dict[bc_name],"query",bc_name)
        resultdict[match]=dict()
        resultdict[match]["score"]=score


    results = sorted([(resultdict[x]["score"],x,resultdict[x]) for x in resultdict.keys()])[::-1]
    #for result in results:
    #    print result
    result = results[0]
    score,ide,details=result
    #print ide.split("_")[0],score,details

    if score >= args.threshold:
        next
    else:
        ide = "unclassified"

    return ide.split("_")[0],score

def nucl_align(sQSeq,sRSeq,query,target):
    #pathtolibssw=pkg_resources.resource_filename('nanonet', 'libssw.so')
    #ospathtolibssw=os.path.dirname(pathtolibssw)
    sQId=query
    sRId=target
    lEle = []
    dRc = {}
    dEle2Int = {}
    dInt2Ele = {}
    lEle = ['A', 'C', 'G', 'T', 'N']
    dRc = {'A':'C', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A'}
    for i,ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    nEleNum = len(lEle)
    lScore = [0 for i in xrange(nEleNum**2)]
    for i in xrange(nEleNum-1):
        for j in xrange(nEleNum-1):
            if lEle[i] == lEle[j]:
                lScore[i*nEleNum+j] = 3
            else:
                lScore[i*nEleNum+j] = -1
# translate score matrix to ctypes
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore
# set flag
    nFlag = 0
    # This line should be the path to libssw.so but I can't get it to work.
    ssw = ssw_lib.CSsw(".")


# build query profile
    qNum = to_int(sQSeq, lEle, dEle2Int)
    qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)
# set mask len
    if len(sQSeq) > 30:
        nMaskLen = len(sQSeq) / 2
    else:
        nMaskLen = 15

# iter target sequence
    rNum = to_int(sRSeq, lEle, dEle2Int)

# format ofres: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
    res = align_one(ssw, qProfile, rNum, len(sRSeq), 3, 1, nFlag, nMaskLen)
    resRc = None
# build cigar and trace back path
    strand = 0
    if resRc == None or res[0] > resRc[0]:
        resPrint = res
        strand = 0
        sCigar, sQ, sA, sR = buildPath(sQSeq, sRSeq, res[4], res[2], res[8])
    else:
        resPrint = resRc
        strand = 1
        sCigar, sQ, sA, sR = buildPath(sQRcSeq, sRSeq, resRc[4], resRc[2], resRc[8])
    #print 'target_name: {}\nquery_name: {}\noptimal_alignment_score: {}\t'.format(sRId, sQId, resPrint[0])
    #print 'suboptimal_alignment_score: {}\t'.format(resPrint[1])
    #print res
    ssw.init_destroy(qProfile)
    return sRId,resPrint[0]

def to_int(seq, lEle, dEle2Int):
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i,ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n

    return num

def align_one(ssw, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    """
    align one pair of sequences
    @param  qProfile   query profile
    @param  rNum   number array for reference
    @param  nRLen   length of reference sequence
    @param  nFlag   alignment flag
    @param  nMaskLen   mask length
    """
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, nMaskLen)

    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)

    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)

def buildPath(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = 'MIDNSHP=X'
    sCigar = ''
    sQ = ''
    sA = ''
    sR = ''
    nQOff = nQryBeg
    nROff = nRefBeg
    for x in lCigar:
        n = x >> 4
        m = x & 15
        if m > 8:
            c = 'M'
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == 'M':
            sQ += q[nQOff : nQOff+n]
            sA += ''.join(['|' if q[nQOff+j] == r[nROff+j] else '*' for j in xrange(n)])
            sR += r[nROff : nROff+n]
            nQOff += n
            nROff += n
        elif c == 'I':
            sQ += q[nQOff : nQOff+n]
            sA += ' ' * n
            sR += '-' * n
            nQOff += n
        elif c == 'D':
            sQ += '-' * n
            sA += ' ' * n
            sR += r[nROff : nROff+n]
            nROff += n

    return sCigar, sQ, sA, sR

class FileExist(argparse.Action):
    """Check if the input file exist."""
    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.exists(values):
             raise RuntimeError("File/path for '{}' does not exist, {}".format(self.dest, values))
        setattr(namespace, self.dest, values)

def parse_barcodes(barcode_file):
    #print "parsing barcodes"
    barcode_list = list()
    barcode_list.append("uncalssified")
    barcode_dict = dict()
    barcode_sequences = SeqIO.parse(open(barcode_file),'fasta')
    for barcode in barcode_sequences:
        name, sequence = barcode.id, str(barcode.seq)
        barcode_dict[name]=sequence
        barcode_list.append(name)
        barcode_dict[name+"_rev"]=str(barcode.reverse_complement().seq)
    #print barcode_list
    #for barcode in barcode_dict:
    #    print barcode, barcode_dict[barcode]

    #sys.exit()
    return barcode_dict,barcode_list

def main():

    args = get_parser().parse_args()
    global barcode_dict
    barcode_dict,barcode_list=parse_barcodes(args.barcodes)

    """barcode_dict = {
        'NB01': 'GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTT',
        'NB01_rev': 'AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCAGGTTA',
        'NB02': 'GGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCTTT',
        'NB02_rev': 'AAGGTTAAACAGACGACTACAAACGGAATCGACAGCACCAGGTTA',
        'NB03': 'GGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCTTT',
        'NB03_rev': 'AAGGTTAACCTGGTAACTGGGACACAAGACTCCAGCACCAGGTTA',
        'NB04': 'GGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCTTT',
        'NB04_rev': 'AAGGTTAATAGGGAAACACGATAGAATCCGAACAGCACCAGGTTA',
        'NB05': 'GGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTT',
        'NB05_rev': 'AAGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACCAGGTTA',
        'NB06': 'GGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCTTT',
        'NB06_rev': 'AAGGTTAAGACTACTTTCTGCCTTTGCGAGAACAGCACCAGGTTA',
        'NB07': 'GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCTTT',
        'NB07_rev': 'AAGGTTAAAAGGATTCATTCCCACGGTAACACCAGCACCAGGTTA',
        'NB08': 'GGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCTTT',
        'NB08_rev': 'AAGGTTAAACGTAACTTGGTTTGTTCCCTGAACAGCACCAGGTTA',
        'NB09': 'GGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCTTT',
        'NB09_rev': 'AAGGTTAAAACCAAGACTCGCTGTGCCTAGTTCAGCACCAGGTTA',
        'NB10': 'GGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCTTT',
        'NB10_rev': 'AAGGTTAAGAGAGGACAAAGGTTTCAACGCTTCAGCACCAGGTTA',
        'NB11': 'GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCTTT',
        'NB11_rev': 'AAGGTTAATCCATTCCCTCCGATAGATGAAACCAGCACCAGGTTA',
        'NB12': 'GGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCTTT',
        'NB12_rev': 'AAGGTTAATCCGATTCTGCTTCTTTCTACCTGCAGCACCAGGTTA'
    }
    barcode_list = ('NB01','NB02','NB03','NB04','NB05','NB06','NB07','NB08','NB09','NB10','NB11','NB12','unclassified')
    """
    resultdict=dict()
    input_file = args.input

    fasta_sequences = SeqIO.parse(open(input_file),'fasta')

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        #new_sequence = some_function(sequence)
        #print ">"+str(name)
        #print sequence

        id_,score=align_seq(sequence,args)
        print str(name),id_,score
        if id_ not in resultdict:
            resultdict[id_]=dict()
            resultdict[id_]["counter"]=0
            resultdict[id_]["score"]=list()
            resultdict[id_]["sequences"]=list()
        resultdict[id_]["counter"]+=1
        resultdict[id_]["score"].append(score)
        resultdict[id_]["sequences"].append(fasta)

    ##print resultdict
    print "Score Threshold:",args.threshold
    for ids in barcode_list:
        if ids in resultdict.keys():
            print ids,
            print resultdict[ids]["counter"],
            print "Mean:", (sum(resultdict[ids]["score"])/resultdict[ids]["counter"])
            output_handle=open(os.path.join(os.path.dirname(input_file),ids+"_"+os.path.basename(input_file)),"w")
            SeqIO.write(resultdict[ids]["sequences"], output_handle, "fasta")
            output_handle.close()
        else:
            print ids,"0","Mean:N/A"




if __name__ == "__main__":
    main()
