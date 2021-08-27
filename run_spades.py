#!/usr/bin/env python3

import argparse, os, sys, shutil, subprocess

help_text = "Usage: python run_spades.py samplenamesfile"

def single_lib_paired_end_cmds(filename, interlaced=True):
    outdir = "{} {}".format("-o", filename)
    cmd = ["spades.py", "-k 21,33,55,77", outdir]
    if not interlaced:
        cmd.append("{} {} {} {}".format("-1", filename + "_R1.fastq", "-2", filename + "_R2.fastq"))
    else:
        cmd.insert(1, "{} {}".format("--12", filename + ".fastq")) 

    cmds = " ".join(cmd)
    return cmds
        

def single_lib_single_read_cmds(filename):
    outdir = "{} {}".format("-o", filename)
    cmd = ["spades.py","-s", filename + ".fastq", "-k 21,33,55,77", outdir]
    cmds = " ".join(cmd)
    return cmds


def iontorrent_unpaired_cmds(filename):
    outdir = "{} {}".format("-o", filename)
    cmd = ["spades.py","--iontorrent","-s", filename + ".fastq", "-k 21,33,55,77", outdir]
    cmds = " ".join(cmd)
    return cmds


def spades(filename, paired=True, interlaced=False, single=False, IonTorrent=False):
    if paired:
        if interlaced:
            cmd = single_lib_paired_end_cmds(filename, interlaced=True)
            sys.stderr.write("Paired-end SPAdes assembly has started on sample {}\n".format(filename))
            sys.stderr.write("Called with arguments " + cmd + "\n")
            execution = subprocess.call(cmd,shell=True)
            if execution:
                sys.stderr.write("SPAdes assembly unsucessful for sample " + filename + "\n")

        else:
            cmd = single_lib_paired_end_cmds(filename, interlaced=False)
            sys.stderr.write("Paired-end SPAdes assembly has started on sample {}\n".format(filename))
            sys.stderr.write("Called with arguments " + cmd + "\n")
            execution = subprocess.call(cmd,shell=True)
            if execution:
                sys.stderr.write("SPAdes assembly unsucessful for sample " + filename + "\n")

    elif single:
        cmd = single_lib_single_read_cmds(filename)
        sys.stderr.write("Single reads SPAdes assembly has started on sample {}\n".format(filename))
        sys.stderr.write("Called with arguments " + cmd + "\n")
        execution = subprocess.call(cmd,shell=True)
        if execution:
            sys.stderr.write("SPAdes assembly unsucessful for sample " + filename + "\n")

    elif IonTorrent:
        cmd = iontorrent_unpaired_cmds(filename)
        sys.stderr.write("IonTorrent reads SPAdes assembly has started on sample {}\n".format(filename))
        sys.stderr.write("Called with arguments " + cmd + "\n")
        execution = subprocess.call(cmd,shell=True)
        if execution:
            sys.stderr.write("SPAdes assembly unsucessful for sample " + filename + "\n")

    else:
        print(help_text)
        
        
def main():
    if len(sys.argv) != 2:
        print (help_text)
        sys.exit()

    samplenamesfile = sys.argv[1]
    if not samplenamesfile:
        print ("You have not specified samples file")
        print (help_text)
        sys.exit()

    fp = open(samplenamesfile)
    for filename in fp:
        try:
            if len(filename.strip().split('\t')) == 2:
                fn = filename.strip().split('\t')[0]
                readformat = filename.strip().split('\t')[1]
                if readformat == "single":
                    spades(fn, paired=False, interlaced=False, single=True, IonTorrent=False)
                elif readformat == "iontorrent":
                    spades(fn, paired=False, interlaced=False, single=False, IonTorrent=True)
                else:
                    print ("There is something wrong with single end library specifications")
                    sys.exit()
            else:
                fn = filename.strip().split('\t')[0]
                readformat = filename.strip().split('\t')[1]
                fileformat = filename.strip().split('\t')[2]
                if readformat == "paired" and fileformat == "separate":
                    spades(fn, paired=True, interlaced=False, single=False, IonTorrent=False)
                elif readformat == "paired" and fileformat == "interlaced":
                    spades(fn, paired=True, interlaced=True, single=False, IonTorrent=False)
                else:
                    print ("There is something wrong with paired end library specifications")
                    sys.exit()
        except Exception as e:
            pass

if __name__ == "__main__":
    main()

        

