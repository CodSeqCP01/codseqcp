# CodSeqCP
Homologous Coding or Genomic Sequence Clustering Pipeline
	PROGRAM
		codseqcp - clusters homologous genomic or coding sequences from
		fastq files of multiple individuals or species for population,
		phylogenetics, etc, analyses

	USAGE
		python codeqcp.py samplenamesfile evalue cds species_name
		to cluster coding sequences

		or

		python codseqcp.py samplenamesfiles evalue genomic
		to cluster genomic sequences

	DESCRIPTION
		samplefilenames is a list of fastq names without ".fastq" extension
		followed by tab delimitted annotations:
		
		for Illumina-based single read library

		samplename	single

		for IonTorrent-based unpaired library

		samplename	iontorrent

		for Illumina-based paired-reads in single file

		samplename	paired	interlaced

		for Illumina-based paired-reads in separate files

		samplename	paired	separate
      
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

	REQUIREMENTS
		CodSeqCP runs on linux and assumes that below programs are in your path
			- SPAdes (for reads assembly)

			- Augustus (for gene prediction)

			- blast+ specifically blastn binary (for sequence alignment)


	SUPPORT
		Contact an author directly: edson.ishengoma@muce.ac.tz
    
