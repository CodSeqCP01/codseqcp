#!/usr/bin/env python3

import sys, os, shutil, argparse, subprocess
from collections import defaultdict

_doc_ = '''
            PROGRAM
                    codseqcp - clusters homologous genomic or coding sequences from
                    fastq files of multiple individuals or species for population,
                    phylogenetics, etc, analyses

            USAGE

                    python codseqcp.py -in samples.txt -e 1e-5 -sp Danio_rerio -o cds
                    to cluster coding sequences

                    or

                    python codseqcp.py -in samples.txt -e 1e-5 -o genomic
                    to cluster genomic sequences


            DESCRIPTION
                    samplefilenames is a list of fastq names without ".fastq" extension
                    followed by tab delimitted annotations:

                    for Illumina-based single read library this would be

                    samplename\tsingle

                    for IonTorrent-based unpaired library this would be

                    samplename\tiontorrent

                    for Illumina-based paired-reads in single file this would be

                    samplename\tpaired\tinterlaced

                    for Illumina-based paired-reads in separate files this would be

                    samplename\tpaired\tseparate
      
                    The program assumes that the fastq files are in the same working directory
                    where the program is located. If paired reads are in separate files,
                    the file with forward reads is suffixed with _R1
                    e.g. samplename_R1.fastq
                    and the file with reverse reads is suffixed with _R2
                    e.g. samplename_R2.fastq

       
                    The final result is a cluster of sequences in FASTA formatted file named

                    clustered_cds.fasta

                        or

                    clustered_genomic.fasta

                    '''.format(sys.argv[0])
            


__author__ = "Edson Ishengoma"
__version__ = "0.1"
__email__ = "edson.ishengoma@muce.ac.tz"

try:
    shutil.rmtree(os.path.join(os.getcwd(), 'store'))
except FileNotFoundError:
    pass

store = os.path.join(os.getcwd(), 'store')

try:
    os.mkdir(store)
except OSError:
    sys.exit()  

def CDSdb():
    # initially parse augustus gff files to fasta files
    sample_names = [f for f in os.listdir(store) if f.endswith("gff")]
    for f in sample_names:
        print("**********************\n\nConverting sample {} into fasta ".format(f))
        parsegff(f)
    

    with open(store + "/" + "pooled_cds.fasta","a") as out_file:
        for f in [fn for fn in os.listdir(store) if fn.endswith("codingseq")]:
            with open(os.path.join(store, f),"r") as infile:
                shutil.copyfileobj(infile, out_file)

    # make a nucleotide blast database of cds
    pooled_cds = store + "/" + "pooled_cds.fasta"
    if os.path.isfile(pooled_cds):
        makeblastdb_cmd = "makeblastdb -dbtype nucl -in {} {}".format(pooled_cds, ">> blast.out")
        print('...Making local blast sequence database with the following arguments...' + '\n' + makeblastdb_cmd)
        exitcode = subprocess.call(makeblastdb_cmd,shell=True)

        if exitcode:
            sys.stderr.write("Something went wrong making of blast nucleotide database")
       

def genomicdb():
    with open(args.samples_names_file) as f:
        spades_workdir = '/'.join(f.readline().split('\t')[0].split('/')[:-1])
        Dirs = [SPAdir for SPAdir in os.listdir(spades_workdir) if os.path.isdir(os.path.join(spades_workdir, SPAdir))]
        paths = map(lambda Dir : Dir + "/scaffolds.fasta", Dirs)
        for path in paths:
            scaffold_file = spades_workdir + "/" + path
            try:
                print("making database for " + scaffold_file)
                parsescaffolds(scaffold_file)
            except IOError as ex:
                print(ex.strerror)
                pass

    with open(store + "/" + "pooled_scaffolds.fasta","a") as out_file:
        for f in [fn for fn in os.listdir(store) if fn.endswith("genomic")]:
            with open(os.path.join(store, f),"r") as infile:
                shutil.copyfileobj(infile, out_file)
               

    # make genomic blast database
    pooled_scaffolds = store + "/" + "pooled_scaffolds.fasta"
    if os.path.isfile(pooled_scaffolds):
        makeblastdb_cmd = "makeblastdb -dbtype nucl -in {} {}".format(pooled_scaffolds, ">> blast.out")
        print('...Making local blast sequence database with the following arguments...' + '\n' + makeblastdb_cmd)
        exitcode = subprocess.call(makeblastdb_cmd,shell=True)

        if exitcode:
            sys.stderr.write("Something went wrong making of blast nucleotide database")

    
def blastn(seq_file, evalue):
        db = seq_file
        blastn_cmd_list = ["blastn","-outfmt 6"]
        e_value = "{} {}".format("-evalue", evalue)
        blastdb = "{} {}".format("-db", db)
        query = "{} {}".format("-query", seq_file)
        blast_out = "{} {}".format(">", os.path.join(store, "blasthits.txt"))
        blastn_cmd_list.extend([blastdb,query,e_value,blast_out])
        blastn_command = " ".join(blastn_cmd_list)
        print( "...Performing blast search on sequences with arguments... " + "\n" +  blastn_command)
        exitcode = subprocess.call(blastn_command,shell=True)
        if exitcode:
            sys.stderr.write("ERROR: blastn was unsuccessful")


def getEntry(name,seq,dict):
    descr = name
    c = 1
    while descr in dict:
        descr = "%s_%i" % (name,c)
        c += 1
    dict[descr] = (seq)


def parsescaffolds(scaffolds_path):
    basename = scaffolds_path.split("/")[-2]
    with open(store + '/' + basename + ".genomic", "a") as outfile:
        for key, seq in parsefasta(open(scaffolds_path)).items():
            node = key.split("_")[1]
            outfile.write(">{}.genomic{} \n".format(basename, node))

            while len(seq) > 0:
                outfile.write(seq[0:60]+"\n")
                seq = seq[60:]

    
def parsegff(file):
    fh = open(os.path.join(store, file))
    result = {}
    name, seq = "",[]
    copy = False
    copy_start = "# coding sequence"
    copy_end = "# protein sequence"
    buffer = ""

    for line in fh.readlines():
        if line.startswith("# start"):
            if name:
                getEntry(name,"".join(seq),result)

            name = line.strip().split(".")[-1]
            seq = []

        else:
            if line.startswith(copy_start):
                copy = True

            elif line.startswith(copy_end):
                copy = False
                buffer += "\n"

            if copy:
                buffer += line
                seqline = line.strip().replace("# coding sequence = [", "").replace("# ", "").replace("]", "")
                seq.append(seqline)

    if name:
        getEntry(name,"".join(seq),result)

    basename = file.split(".")[0]
    with open(store + '/' + basename + ".codingseq", "a") as outfile:
        for key, seq in result.items():
            node = key.replace(key, key[1:]) #remove the g
            outfile.write(">{}.codingseq{} \n".format(basename, node))

            while len(seq) > 0:
                outfile.write(seq[0:60]+"\n")
                seq = seq[60:]


def parsefasta(file):
    result = {}
    name, seq = "",[]
    for line in file.readlines():
        if line.startswith(">"):
            if name:
                getEntry(name,"".join(seq),result)

            name = line[1:].strip()
            seq = []
        else:
            seqline = line.strip().replace(" ", "")
            seq.append(seqline)
    if name:
        getEntry(name,"".join(seq),result)

        return result

def parseblasthits(blasthitsfile):
    qdefs = set()
    for line in open(blasthitsfile):
        if line:
            qdef = line.split("\t")[0]
            qdefs.add(qdef)

    blstxtrct = defaultdict(list) # extracted key blast results 

    for line in open(blasthitsfile):
        qdef = line.split("\t")[0]
        hdef = line.split("\t")[1]  
        pid = line.split("\t")[2]
        qstart = line.split("\t")[6]
        qend = line.split("\t")[7]
        hitstart = line.split("\t")[8]
        hitend = line.split("\t")[9]

        for seqid in qdefs:
            if hdef == seqid and hdef == line.split("\t")[1]:
                blstxtrct[seqid].append((qdef,pid,qstart,qend,hitstart,hitend))

    # first filter results to retain top hits per sample based on percent identity
    tophits = {}
    for sbjct, hits in blstxtrct.items():
        sortedhits = sorted(hits, key=lambda x: float(x[1]), reverse=True) #sort by percent identity

        #print(hdef)
        flg = set()
        filtered = []
        for elem in sortedhits:
            if elem[0].split(".")[0] not in flg:
                filtered.append(elem)
                flg.add(elem[0].split(".")[0])
        tophits[sbjct] = filtered

    # remove redundant sample clusters
    # iteration 1
    seen = set()
    keep = {}

    for elem, tp_elms in tophits.items():
        if all(tup in seen for tup in tp_elms):
            continue
        keep[elem] = tp_elms
        seen.update(tp_elms)

    # iteration 2
    final = {}
    s = set()

    for k, v in keep.items():
        if k not in s:
            final[k] = v
            for item in v:
                s.add(item[0])
    return final

def clusterfasta(blasthitfile, seqfile, clusterfile):
    blsthits = parseblasthits(blasthitfile)
    seqs = parsefasta(open(seqfile))
    keys = [i+1 for i in range(len(blsthits))]
    numclust = dict(zip(keys, list(blsthits.values()))) # clusters numbered numerically

    with open(clusterfile, "a") as outfile:
        for clust, hits in numclust.items():
            for rec in hits:
                qdef = rec[0]
                qstart = int(rec[2])
                qend = int(rec[3])
                if qdef in seqs:
                    qbasename = qdef.split(".")[0]
                    seq = seqs[qdef].upper()
                    alignedseq = seq[qstart-1:qend] # extract only aligned blastn hsps
                    outfile.write(">"+ qdef.replace(qdef, "Cluster_" + str(clust) + "_" + qbasename) + "\n")

                    while len(alignedseq) > 0:
                        #write sixty characters per line
                        outfile.write(str(alignedseq[0:60]) +"\n")
                        alignedseq = str(alignedseq [60:])

def clusterGenomic(samplenamesfile, evalue, seq_type):
    #check the intergrity of sample configuration file
    if not os.path.isfile(samplenamesfile):
        print("\n" + "------Problem with assembly of sample reads------" + "\n")
        print("Something may be wrong with your sample reads or configuration of your sample file" + "\n")
        print(_doc_)
        sys.exit()
    else:
        print("\n" + ".....CodSeqCP pipeline has started..." + "\n")
        spades(samplenamesfile)
        genomicdb()
        seqfile = store + "/" + "pooled_scaffolds.fasta"
        blastn(seqfile, evalue)

    # Get pooled blastn outfmt 6 file
    blasthitfile = os.path.join(store, "blasthits.txt")
    clusterfile = os.path.join(os.getcwd(), "clustered_genomic.fasta")
    return clusterfasta(blasthitfile, seqfile, clusterfile)
    


def clusterCDS(samplenamesfile, evalue, seq_type, species_name):
    if not os.path.isfile(samplenamesfile):
        print("------Problem with assembly of sample reads------" + "\n")
        print("**Something may be wrong with your sample reads or configuration of your sample file**" + "\n")
        print(_doc_)
        sys.exit()
    else:
        print("\n" + ".....CodSeqCP pipeline has started..." + "\n")
        spades(samplenamesfile)
        augustus(samplenamesfile, augustus_species=species_name)
        CDSdb()
        seqfile = os.path.join(store, "pooled_cds.fasta")
        blastn(seqfile, evalue)
        # Get pooled blastn outfmt 6 file
        blasthitfile = os.path.join(store, "blasthits.txt")
        clusterfile = os.path.join(os.getcwd(), "clustered_cds.fasta")
        return clusterfasta(blasthitfile, seqfile, clusterfile)
    
def spades(samplenamesfile):
    spades_scrpt_path = os.path.join(os.getcwd(), "run_spades.py")
    spades_command = "python3 {} {} {}".format(spades_scrpt_path, samplenamesfile, '>> spades.out')
    exitcode = subprocess.call(spades_command,shell=True)

def augustus(samplenamesfile, augustus_species=None):
    # e.g. python augustus_runner.py Danio_rerio
    augustus_scrpt_path = os.path.join(os.getcwd(), "run_augustus.py")
    augustus_command = "python3 {} {} {}".format(augustus_scrpt_path, samplenamesfile, augustus_species)
    exitcode = subprocess.call(augustus_command,shell=True)


parser = argparse.ArgumentParser()
parser.add_argument('-in', '--input', action='store', dest='samples_names_file',
                    help='Enter full path of file containing your samples names', required=True)
parser.add_argument('-o','--output', action='store', dest='output_seq_type',
                    help='Enter cds to output coding sequences or enter genomic to output genomic sequences', required=True)
parser.add_argument('-e','--evalue', action='store', dest='e_value',
                    help='Enter blast search evalue threshold', required=True)
parser.add_argument('-sp','--species', action='store', dest='species_name',
                    help='Enter the species name preconfigured within Augustus as a model for gene prediction', required=False)
args = parser.parse_args()
    
   
def main():
    if args.output_seq_type == 'genomic':
        # Dont need Augustus predicted cds, only genomic dna
        clusterGenomic(args.samples_names_file, args.e_value, args.output_seq_type)  

    elif args.output_seq_type == 'cds':
        clusterCDS(args.samples_names_file, args.e_value, args.output_seq_type, args.species_name)

    else:
        print(_doc_)
        sys.exit()
        
if __name__ == "__main__":
    main()
