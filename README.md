# CodSeqCP
Homologous Coding or Genomic Sequence Clustering Pipeline
	
	PROGRAM
		codseqcp - clusters homologous genomic or coding sequences from
		fastq files of multiple individuals or species for population,
		phylogenetics, etc, analyses

	USAGE
		python codeqcp.py samplenamesfile evalue cds species_name
		example: python codseqcp.py samples.txt 1e-5 cds Danio_rerio
		run the above command to cluster coding sequences using the 
		data set in the test folder

		or

		python codseqcp.py samplenamesfiles evalue genomic
		example: python codseqcp.py samples.txt 1e-5 genomic
		run the above command to cluster genomic sequences
		using the data set in the test folder

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
		
		evalue is an integer value that describes the number of hits one can "expect" 
		to return just by chance during a blast search

		species_name is a species pre-trained in Augustus for gene prediction.
		Convention used in CodSeqCP is to have a species name connected by an underscore "_".
		e.g. Danio_rerio

		List of supported Augustus species is as follows
		
		identifier	species	CodSeqCP 
		human	Homo sapiens	Homo_sapiens
		fly	Drosophila melanogaster	Drosophila_melanogaster
		arabidopsis	Arabidopsis thaliana	Arabidopsis_thaliana
		brugia	Brugia malayi	Brugia_malayi
		aedes	Aedes aegypti	Aedes_aegypti
		tribolium	Tribolium castaneum	Tribolium_castaneum
		schistosoma	Schistosoma mansoni	Schistosoma_mansoni
		tetrahymena	Tetrahymena thermophila	Tetrahymena_thermophila
		galdieria	Galdieria sulphuraria	Galdieria_sulphuraria
		maize	Zea mays	Zea_mays
		toxoplasma	Toxoplasma gondii	Toxoplasma_gondii
		caenorhabditis	Caenorhabditis elegans	Caenorhabditis_elegans
		aspergillus_fumigatus	Aspergillus fumigatus	Aspergillus_fumigatus
		aspergillus_nidulans	Aspergillus nidulans	Aspergillus_nidulans
		aspergillus_oryzae	Aspergillus oryzae	Aspergillus_oryzae
		aspergillus_terreus	Aspergillus terreus	Aspergillus_terreus
		botrytis_cinerea	Botrytis cinerea	Botrytis_cinerea
		candida_albicans	Candida albicans	Candida_albicans
		candida_guilliermondii	Candida guilliermondii	Candida_guilliermondii
		candida_tropicalis	Candida tropicalis	Candida_tropicalis
		chaetomium_globosum	Chaetomium globosum	Chaetomium_globosum
		coccidioides_immitis	Coccidioides immitis	Coccidioides_immitis
		coprinus	Coprinus cinereus	Coprinus_cinereus
		coyote_tobacco	Nicotiana attenuata	Nicotiana_attenuata
		cryptococcus_neoformans_gattii	Cryptococcus neoformans gattii	Cryptococcus_neoformans_gattii
		cryptococcus_neoformans_neoformans_B	Cryptococcus neoformans neoformans	Cryptococcus_neoformans_neoformans
		debaryomyces_hansenii	Debaryomyces hansenii	Debaryomyces_hansenii
		encephalitozoon_cuniculi_GB	Encephalitozoon cuniculi	Encephalitozoon_cuniculi
		eremothecium_gossypii	Eremothecium gossypii	Eremothecium_gossypii
		fusarium_graminearum	Fusarium graminearum	Fusarium_graminearum
		histoplasma_capsulatum	Histoplasma capsulatum	Histoplasma_capsulatum
		kluyveromyces_lactis	Kluyveromyces lactis	Kluyveromyces_lactis
		laccaria_bicolor	Laccaria bicolor	Laccaria_bicolor
		lamprey	Petromyzon marinus	Petromyzon_marinus
		leishmania_tarentolae	Leishmania tarentolae	Leishmania_tarentolae
		lodderomyces_elongisporus	Lodderomyces elongisporus	Lodderomyces_elongisporus
		magnaporthe_grisea	Magnaporthe grisea	Magnaporthe_grisea
		neurospora_crassa	Neurospora crassa	Neurospora_crassa
		phanerochaete_chrysosporium	Phanerochaete chrysosporium	Phanerochaete_chrysosporium
		pichia_stipitis	Pichia stipitis	Pichia_stipitis
		rhizopus_oryzae	Rhizopus oryzae	Rhizopus_oryzae
		saccharomyces_cerevisiae_S288C	Saccharomyces cerevisiae	Saccharomyces_cerevisiae
		schizosaccharomyces_pombe	Schizosaccharomyces pombe	Schizosaccharomyces_pombe
		thermoanaerobacter_tengcongensis	Thermoanaerobacter tengcongensis	Thermoanaerobacter_tengcongensis
		trichinella	Trichinella spiralis	Trichinella_spiralis
		ustilago_maydis	Ustilago maydis	Ustilago_maydis
		yarrowia_lipolytica	Yarrowia lipolytica	Yarrowia_lipolytica
		nasonia	Nasonia vitripennis	Nasonia_vitripennis
		tomato	Solanum lycopersicum	Solanum_lycopersicum
		chlamydomonas	Chlamydomonas reinhardtii	Chlamydomonas_reinhardtii
		amphimedon	Amphimedon queenslandica	Amphimedon_queenslandica
		pneumocystis	Pneumocystis jirovecii	Pneumocystis_jirovecii
		wheat	Triticum aestivum	Triticum_aestivum
		chicken	Gallus gallus	Gallus_gallus
		zebrafish	Danio rerio	Danio_rerio
		E_coli_K12	Escherichia coli	Escherichia_coli
		s_aureus	Staphylococcus aureus	Staphylococcus_aureus
		volvox	Volvox carteri	Volvox_carteri
		
		The final result is a cluster of sequences in FASTA formatted file named

		clustered_cds.fasta

			or

		clustered_genomic.fasta

	REQUIREMENTS
		CodSeqCP runs on linux and assumes that below programs are in your path
			- SPAdes (for reads assembly)

			- Augustus (for gene prediction)

			- blast+ (for sequence alignment)


	SUPPORT
		Contact an author directly: edson.ishengoma@muce.ac.tz


    
