#!/usr/bin/env python3

import argparse, os, sys, shutil, subprocess

Help_text= "Usage: python run_augustus.py species e.g python run_augustus.py samplenamesfile Danio_rerio"

# Currently supported augustus species
d = {}
#A
d["Acyrthosiphon_pisum"] = "pea_aphid"
d["Aedes_aegypti"] = "aedes"
d["Amphimedon_queenslandica"] = "amphimedon"
d["Ancylostoma_ceylanicum"] = "ancylostoma_ceylanicum"
d["Apis_dorsata"] = "adorsata"
d["Apis_mellifera"] = "honeybee1"
d["Arabidopsis_thaliana"] = "arabidopsis"
d["Aspergillus_fumigatus"] = "aspergillus_fumigatus"
d["Aspergillus_nidulans"] = "aspergillus_nidulans"
d["Aspergillus_oryzae"] = "aspergillus_oryzae"
d["Aspergillus_terreus"] = "aspergillus_terreus"

#B
d["Bombus_impatiens"] = "bombus_impatiens1"
d["Bombus_terrestris"] = "bombus_terrestris2"
d["Botrytis_cinerea"] = "botrytis_cinerea"
d["Brugia_malayi"] = "brugia"
d["Burkholderia_pseudomallei"] = "b_pseudomallei"

#C
d["Caenorhabditis_elegans"] = "caenorhabditis"
d["Callorhinchus_milii"] = "elephant_shark"
d["Camponotus_floridanus"] = "camponotus_floridanus"
d["Candida_albicans"] = "candida_albicans"
d["Candida_guilliermondii"] = "candida_guilliermondii"
d["Candida_tropicalis"] = "candida_tropicalis"
d["Chaetomium_globosum"] = "chaetomium_globosum"
d["Chiloscyllium_punctatum"] = "chiloscyllium"
d["Chlamydomonas_reinhardtii"] = "chlamy2011"
d["Chlorella_sp."] = "chlorella"
d["Ciona_intestinalis"] = "ciona"
d["Coccidioides_immitis"] = "coccidioides_immitis"
d["Conidiobolus_coronatus"] = "Conidiobolus_coronatus"
d["Coprinus_cinereus"] = "coprinus_cinereus"
d["Cryptococcus_neoformans gattii"] = "cryptococcus_neoformans_gattii"
d["Cryptococcus_neoformans neoformans"] = "cryptococcus_neoformans_neoformans_B"
d["Culex_pipiens"] = "culex"

#D
d["Danio_rerio"] = "zebrafish"
d["Debaryomyces_hansenii"] = "debaryomyces_hansenii"
d["Drosophila_melanogaster"] = "fly"

#E
d["Encephalitozoon_cuniculi"] = "encephalitozoon_cuniculi_GB"
d["Eremothecium_gossypii"] = "eremothecium_gossypii"
d["Escherichia_coli K12"] = "E_coli_K12"

#F
d["Fusarium_graminearium"] = "fusarium_graminearum"

#G
d["Galdieria_sulphuraria"] = "galdieria"
d["Gallus_gallus_domesticus"] = "chicken"

#H
d["Heliconius_melpomene"] = "heliconius_melpomene1"
d["Histoplasma_capsulatum"] = "histoplasma_capsulatum"
d["Homo_sapiens"] = "human"

#K
d["Kluyveromyces_lactis"] = "kluyveromyces_lactis"

#L
d["Laccaria_bicolor"] = "laccaria_bicolor"
d["Leishmania_tarentolae"] = "leishmania_tarentolae"
d["Lethenteron_camtschaticum"] = "japaneselamprey"
d["Lodderomyces_elongisporus"] = "lodderomyces_elongisporus"

#M
d["Magnaporthe_grisea"] = "magnaporthe_grisea"
d["Mnemiopsis_leidyi"] = "mnemiopsis_leidyi"

#N
d["Nasonia_vitripennis"] = "nasonia"
d["Nematostella_vectensis"] = "Neurospora crassa"
d["Neurospora_crassa"] = "neurospora_crassa"
d["Nicotiana_attenuata"] = "coyote_tobacco"

#O
d["Oryza_sativa"] = "rice"

#P
d["Parasteatoda_sp."] = "parasteatoda"
d["Petromyzon_marinus"] = "sealamprey"
d["Phanerochaete_chrysosporium"] = "phanerochaete_chrysosporium"
d["Pichia_stipitis"] = "pichia_stipitis"
d["Pisaster_ochraceus"] = "pisaster"
d["Plasmodium_falciparum"] = "pfalciparum"
d["Pneumocystis_jirovecii"] = "pneumocystis"

#R
d["Rhincodon_typus"] = "rhincodon"
d["Rhizopus_oryzae"] = "rhizopus_oryzae"
d["Rhodnius_prolixus"] = "rhodnius"
d["Saccharomyces_cerevisiae"] = "saccharomyces_cerevisiae_rm11-1a_1"
d["Schistosoma_mansoni"] = "schistosoma"
d["Schizosaccharomyces_pombe"] = "schizosaccharomyces_pombe"
d["Scyliorhinus_torazame"] = "scyliorhinus"
d["Solanum_lycopersicum"] = "tomato"
d["Staphylococcus_aureus"] = "s_aureus"
d["Streptococcus_pneumoniae"] = "s_pneumoniae"
d["Strongylocentrotus_purpuratus"] = "strongylocentrotus_purpuratus"
d["Sulfolobus_solfataricus"] = "Sulfolobus solfataricus"

#T
d["Tetrahymena_thermophila"] = "tetrahymena"
d["Theobroma_cacao"] = "cacao"
d["Thermoanaerobacter_tengcongensis"] = "thermoanaerobacter_tengcongensis"
d["Toxoplasma_gondii"] = "toxoplasma"
d["Tribolium_castaneum"] = "tribolium2012"
d["Trichinella_sp."] = "trichinella"
d["Triticum_sp."] = "wheat"

#U
d["Ustilago_maydis"] = "ustilago_maydis"

#V
d["Verticillium_albo_atrum"] = "verticillium_albo_atrum1"
d["Verticillium_longisporum"] = "verticillium_longisporum1"
d["Volvox_sp."] = "volvox"

#X
d["Xipophorus_maculatus"] = "Xipophorus_maculatus"

#Y
d["Yarrowia_lipolytica"] = "yarrowia_lipolytica"

#Z
d["Zea_mays"] = "maize"

#set species code fuction, scf
def scf(sps):
    if sps in d:
        return d[sps]

 
def augustus(scaffolds, outfilename, species=None):
    #basic augustus commands
    augustus_cmd_list = ["augustus --gff3=on --stopCodonExcludedFromCDS=false --uniqueGeneId=true --codingseq=on"]
    if species:
        spcode = "{}{}".format("--species=", scf(species))
        augustus_cmd_list.append(spcode)

    #input scaffold file
    if scaffolds:
        augustus_cmd_list.append(scaffolds)
    
    #output file in gff format
    augustus_cmd_list.append("{} {}.gff".format(">", os.path.join(os.getcwd(), 'store', outfilename)))
    cmds = " ".join(augustus_cmd_list)
    sys.stderr.write("Running coding sequence prediction on {}\n\n".format(scaffolds))
    sys.stderr.write(cmds + "\n")
    exitcode = subprocess.call(cmds,shell=True)
    if exitcode:
        sys.stderr.write("ERROR: augustus gene prediction was unsuccessful for " + scaffolds + "\n\n")
              
def main():
    if len(sys.argv) != 3:
        print (Help_text)
        sys.exit()

    sps = any(elem for elem in d.keys())
    samples_names_file = sys.argv[1]
    with open(samples_names_file) as f:
        spades_workdir = '/'.join(f.readline().split('\t')[0].split('/')[:-1])
        Dirs = [SPAdir for SPAdir in os.listdir(spades_workdir) if os.path.isdir(os.path.join(spades_workdir, SPAdir))]
        paths = map(lambda Dir : Dir + "/scaffolds.fasta", Dirs)
        for path in paths:
            scaffold_file = spades_workdir + "/" + path
            if scaffold_file:
                if sps:
                    sps = sys.argv[2]
                    augustus(scaffold_file, path.split("/")[0], species=sps)
                else:
                    print("Oops! Seems invalid arguments given for calling genes with augustus")
                    sys.exit()

if __name__ == "__main__":
    main()

        


