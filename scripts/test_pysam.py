
"""
Test pysam
"""
import pysam
import numpy as np

def Test0():
    import pysam
    bf = pysam.AlignmentFile("in.bam", "rb" )
    for pileupcolumn in bf.pileup("chrM", 300, 301):

        for read in [al for al in pileupcolumn.pileups if al.alignment.mapq>30]:

            if not read.is_del and not read.is_refskip:
                if pileupcolumn.pos + 1 == 301:
                    print read.alignment.query_name, read.alignment.query_sequence[read.query_position]

    bf.close()


def Test1():
    samfile = pysam.Samfile('1006-01.bam', 'rb')

    print samfile.header['SQ'][0]
    for pileupColumn in samfile.pileup('1', 10000, 10001 ) :

        a = [ reads.alignment.mapq for reads in pileupColumn.pileups if reads.alignment.mapq >= 0 ]
        print pileupColumn.tid, pileupColumn.pos, pileupColumn.n, np.array(a).mean(), a

        for reads in pileupColumn.pileups :
            reads.alignment.tags += [('AS', 1000)]
            print reads.alignment.opt('RG'), \
                  reads.alignment.opt('AS'), dict(reads.alignment.tags)['AS'],  \
                  reads.alignment.flag, reads.alignment.isize, '%', reads.alignment.pos, reads.alignment.qstart, \
                  reads.alignment.rlen, '#', reads.alignment

    samfile.close()

def Test2 () :
    samfile = pysam.Samfile('1006-01.bam', 'rb')

    for pileupColumn in samfile.pileup('1', 10000, 10005 ) :

        for reads in [ al for al in pileupColumn.pileups if al.alignment.mapq >= 0 ]:
            print pileupColumn.pos+1, len(reads.alignment.positions), reads.alignment.aend, reads.alignment.alen, \
                  reads.alignment.inferred_length,  reads.alignment.is_duplicate, samfile.header['SQ'][reads.alignment.rname], \
                  reads.alignment.is_paired, '#',len(reads.alignment.query),len(reads.alignment.seq), reads.alignment.qstart, \
                  reads.alignment.qend, reads.alignment.rlen,len(reads.alignment.qual)

    samfile.close()

def Test3 () :
    samfile = pysam.Samfile('1006-01.bam', 'rb')

    for read in samfile.fetch( '1', 10000, 10005 ) :
        print read.opt('AS'), read.opt('RG'), read.tags

#####################################################################################################
if __name__ == '__main__' :

    Test1()

