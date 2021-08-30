# CodSeqCP
Homologous Coding or Genomic Sequence Clustering Pipeline
Clusters genomic or coding sequences from raw fastq file of multiple individuals or species 
for population, phylogenetics, etc analyses

**Usage**

Command line arguments for clustering genomic sequences
		
	python codseqcp.py -in SAMPLES_NAMES_FILE, -o OUTPUT_SEQ_TYPE -e E_VALUE
	
	#example
	python codseqcp.py -in samples.txt -o genomic -e 1e-5 
	
For clustering coding sequences
				
				
				
		python codseqcp.py -in SAMPLES_NAMES_FILE, -o OUTPUT_SEQ_TYPE -e E_VALUE -sp SPECIES_NAME

		#example
		python codseqcp.py -in samples.txt -e 1e-5 -o cds -sp Danio_rerio
		

**Description**

_SAMPLES_NAMES_FILE_  
A text file which specify the location, the name of your sample file and whether the library is single or paired (e.g.)

for Illumina-based single libraries
				
		samplename	single

for IonTorrent-based unpaired library

		samplename	iontorrent

for Illumina-based paired-reads in single file

		samplename	paired	interlaced

for Illumina-based paired-reads in separate files

		samplename	paired	separate
      
If paired reads are in separate files,
		the file with forward reads is suffixed with _R1
		e.g. samplename_R1.fastq
		and the file with reverse reads is suffixed with _R2
		e.g. samplename_R2.fastq
		
_OUTPUT_SEQ_TYPE_
 
 Can either be cds in which case the pipeline will run augustus on assembled sequences to extract to output sequences with
 coding potential, or genomic in which case the pipeline will cluster genomic sequences right after the assembly step.
 
 _E_VALUE_

An integer value that describes the number of hits one can "expect" to return just by chance during a blast search

_SPECIES_NAME_

Is the name of species whose gene models are already pre-configured in Augustus for gene prediction. CODSEQCP co-opt these
species as proxies for cds prediction

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
		

**Output results**

The final result is a cluster of sequences in FASTA formatted file named

		clustered_cds.fasta

or

		clustered_genomic.fasta

	
**Prerequisites**
		

			- Python
			- Linux 64 bits
			- SPAdes (for reads assembly)
			- Augustus (for gene prediction)
			- blast+ (for sequence alignment)
			
**Installation**

You can directly clone the pipeline directly as so,

	# get CodSeqCP set of scripts
		git clone https://github.com/CodSeqCP01/codseqcp.git 
		
		# cd codseqcp
		# modify the path of your samples in samples.txt to indicate the actual location of your samples; 
		# be sure BLAST, SPAdes and Augustus are installed and available in your path then run;
		
		python codseqcp.py -in samples.txt -e 1e-5 -o cds -sp Danio_rerio
		
		to test if the program is running
		
or you can use the CODSEQCP conda installation to escape the hurdles of installing individual requirements separately
	
	# be sure you have anaconda or its ligter version "miniconda" installed;
	
	# download the conda environment file here;
	
	wget https://raw.githubusercontent.com/CodSeqCP01/codseqcp/master/codseqcp-0.0.1-linux64.yml
	
	# install CODSEQCP and all its dependencies in a conda environment named codseqcp-0.0.1
	conda env create -f codseqcp-0.0.1-linux64.yml
	
	# activate your codseqcp-0.0.1 environment
	conda activate codseqcp-0.0.1
	
	# CODSEQCP set of scripts will be located at something like  ~/miniconda3/envs/codeqcp-0.0.1/lib/python3.8/site-packages/CODSEQCP
	# you can move these scripts in some convenient working directory and run the main codseqcp.py script while you are still in the 
		active codseqcp conda environment
	
	
**Support**
		
		Contact edson.ishengoma@muce.ac.tz or clintr@sun.ac.za for technical and/or scientific support


    
