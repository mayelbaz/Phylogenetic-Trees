import pylab
from Bio import Entrez
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Align.Applications import MuscleCommandline

# keys - the genes that are tool can work on
genes = {
    "P53":{
        # keys - the organism, value - the NCBI barcode for the gene and organism
        "human":"NM_000546.6",
        "cat": "NM_001009294.1",
        "mouse":"NM_011640.3",
        "lizard": "XM_035135727.1",
        "zebra_fish": "NM_001271820.1",
        "sheep": "NM_001009403.1",
        "rabbit": "XM_008270660.2",
        "horse": "NM_001202405.1",
        "chimp": "XM_016931470.2",
        "swamp eel":"XM_020585891.1",
    },
    "BRCA1":{
        # keys - the organism, value - the NCBI barcode for the gene and organism
        "human": "NM_007294.4",
        "mouse": "NM_009764.3",
        "chicken": "XM_015299561.2",
        "alligator": "XM_014609224.2",
        "duck": "XM_035347820.1",
        "pigeon": "XM_021281038.1",
        "turkey": "XM_010724647.3",
        "penguin": "XM_019473014.1",
        "hedgehog": "XM_016194170.1",
        "fox": "XM_025987104.1",
    },
    "APOE":{
        # keys - the organism, value - the NCBI barcode for the gene and organism
        "human": "NM_000041.4",
        "rat": "NM_138828.3",
        "zebra_fish": "NM_001020565.1",
        "chimpanzee": "NM_001009007.1",
        "fox": "XM_026013924.1",
        "seal": "XM_032388067.1",
        "killer whale": "XM_004271187.2",
        "goat": "XM_018062561.1",
        "polar bear": "XM_008684597.1",
        "koala": "XM_021003616.1",
    },
    "VEGFA":{
        # keys - the organism, value - the NCBI barcode for the gene and organism
        "human": "NM_003376.6",
        "rat": "NM_031836.3",
        "rabbit": "XM_017345155.1",
        "dog": "NM_001003175.2",
        "ferret": "XM_004740321.2",
        "crocodile": "XM_019549669.1",
        "lizard": "XM_035108188.1",
        "pigeon": "XM_021289394.1",
        "eagle": "XM_009919507.1",
        "puma": "XM_025926021.1",
    }

}

job_done_precentage = 0
photo_ready = False
saved_file_name = "tree.png"
alignment_file_name = ""

# input - the gene required , output - the list of organisms that we have for this gene
def get_animals_for_gene(gene):
    if gene in genes.keys():
        return genes[gene].keys()
    else:
        return []

# output - the optional genes
def get_optional_genes():
    return genes.keys()

# draws the input tree and saves in file - saved_file_name
def display_tree(tree):
    global photo_ready
    tree.ladderize()  # Flip branches so deeper clades are displayed at top
    Phylo.draw(tree,show_confidence=False,do_show=False)
    pylab.savefig(saved_file_name)
    photo_ready = True

# opens up the input alignment_file and creates a tree object
def alignment_to_tree(alignment_file):
    aln = AlignIO.read(alignment_file, 'fasta')
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    return tree

# opens the fasta input file and creates an alignment file
def create_alignment(input_file):
    in_file = input_file
    out_file = "aligned.fasta"
    # the location of the muscle EXE on the computer
    # installation guide can be found here
    # https://2018-03-06-ibioic.readthedocs.io/en/latest/install_muscle.html
    muscle_exe = r"/Users/elizabethfedorovsky/miniconda3/bin/muscle"

    muscle_cline = MuscleCommandline(muscle_exe,input=in_file,
                                     out=out_file,
                                     diags=True,
                                     maxiters=1,
                                     log="align_log.txt")
    muscle_cline()
    return out_file

# input - the gene and organisms required
# output - the created FASTA file
def create_fasta_file(gene,animals):

    # for live updates on the process -
    global job_done_precentage
    if len(animals) > 0:
        percentage_to_add = 50/len(animals)
    else:
        percentage_to_add = 0

    Entrez.email = "fedorlizzie@campus.technion.com"  # Always tell NCBI who you are
    filename = gene+"_for_alignment.fasta"
    out_handle = open(filename, "w")
    full_fasta = ""
    for key in genes[gene]:
        if key in animals:
            print("KEY:" ,key)
            # Downloading...
            net_handle = Entrez.efetch(
                db="nucleotide", id=genes[gene][key], rettype="fasta", retmode="text"
            )
            fasta = net_handle.read()
            if full_fasta != "":
                full_fasta = full_fasta +"\n" +fasta
            else:
                full_fasta = fasta
            net_handle.close()
            print("Saved")
            job_done_precentage = job_done_precentage+percentage_to_add

    out_handle.write(full_fasta)
    out_handle.close()
    return filename

# changes the names on the tree's clade using the dict provided
def change_names_in_clade(my_dict, myclade):
    if len(myclade.clades) != 0:
        myclade.name = ""
    for clade in myclade.clades:
        if clade.name in my_dict.keys():
            clade.name = my_dict[clade.name]
        else:
            clade.name = ""
        change_names_in_clade(my_dict, clade)

# changes the names on the tree's clade using the gene and animals
def change_names_in_tree(gene,animals,tree):
    #creating the dict of NCBI barcode -> animal name
    dict_gene_id_to_animal = {}
    for animal in animals:
        dict_gene_id_to_animal[genes[gene][animal]] = animal

    tree.name = ""
    change_names_in_clade(dict_gene_id_to_animal,tree.clade)

# running the process for a gene and a list of organisms
def start_program(gene,animals):
    global alignment_file_name, job_done_precentage
    print(gene)
    print(animals)
    if len(animals) >= 3:
        job_done_precentage = 0
        fasta_filename = create_fasta_file(gene,animals)
        alignment_file_name = create_alignment(fasta_filename)
        job_done_precentage = 75
        tree = alignment_to_tree(alignment_file_name)
        job_done_precentage = 85
        change_names_in_tree(gene,animals,tree)
        display_tree(tree)
        job_done_precentage = 100
