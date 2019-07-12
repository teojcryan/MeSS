from Bio import Entrez
import random as rd
import pandas as pd

def ref_or_not(taxid, search_term):
    '''
    function that searches for reference genomes in the NCBI assembly database
    its output is a dictionary with assembly info (or empty if no reference genome found)
    taxid=taxonomy id searched from the list of species
    search term= search term defined by the user
    '''
    id_list = Entrez.read(Entrez.esearch(db='assembly', term=search_term % taxid))['IdList']
    print('\nSearch term:', search_term)
    ref = True
    if len(id_list) > 1:
        print('\nfound %i reference genomes' % len(id_list))
        random_pick = rd.choice(id_list)
        print('\nRandom choice of %s (assembly ID) as the reference' % random_pick)
        genome_summary = Entrez.read(Entrez.esummary(db='assembly', id=random_pick))
    if len(id_list) == 1:
        print('\nfound 1 reference genome')
        genome_summary = Entrez.read(Entrez.esummary(db='assembly', id=id_list))
    if len(id_list) == 0:
        ref = False
        print('\nDid not find reference genome ! (Restricting search term)')

    if ref:
        acc = genome_summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
        name = genome_summary['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
        status = genome_summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
        reference_status = genome_summary['DocumentSummarySet']['DocumentSummary'][0]['RefSeq_category']
        ftppath_refseq = genome_summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        dic = {'acc': acc, 'name': name, 'status': status, 'taxid': taxid, 'reference': reference_status,
               'path': ftppath_refseq}
    if not ref:
        dic = None
    return dic

def generate_even_prop(species_names, pv, pb, ph):
    '''
    function that calculates the number of viral human and bacterial genomes from a list
    and assigns even proportions for each genome

    ph=proportion of human reads
    pb=proportion of bacterial reads
    pv=proportion of viral reads
    species_names=list of names inputted by the user
    '''
    read_prop = [0] * len(species_names)
    nh = 0
    nv = 0
    nb = 0
    for i in species_names:
        if i == 'Homo sapiens':
            nh += 1
        if 'virus' in i:
            nv += 1
    nb = len(species_names) - nh - nv  # we assume that the rest of the list has only bacterial genomes
    print('nb of genomes \nhuman:', nh, '\nbacterial:', nb, '\nviral', nv)
    # print('number of human genomes%i\n number of bacterial genomes%i\n number of viral genomes%i'%nh,nb,nv)
    for i in range(len(species_names)):
        read_prop[i] = pb / nb
        if species_names[i] == 'Homo sapiens':
            read_prop[i] = ph / nh
        if 'virus' in species_names[i]:
            read_prop[i] = pv / nv
    even_reads = {'name': species_names, 'percent_reads': read_prop}
    return even_reads

def genome_list_to_csv(input_list, read_percentage):
    '''
    function that iterates ref_or_not() and restricts its search terms if the output is None
    input_list=input of species names (from summary_table.tsv)
    read_percentage=array generated by generate_even_prop()
    '''
    acc_list = []
    prop_list = read_percentage
    name_list = []
    status_list = []
    taxid_list = []
    ref_list = []
    path_list = []
    for i in input_list:
        print('\n searching for', i, '...')
        taxid = Entrez.read(Entrez.esearch(db='taxonomy', term=i))['IdList'][0]
        search_term = 'txid%s[Organism:exp] AND ("complete genome"[filter]) AND ("latest refseq"[filter]) AND ("representative genome"[filter] OR "reference genome"[filter])'
        genome_of_choice = ref_or_not(taxid, search_term)

        if taxid == '9606':#For the human id the previous search term does not work because the assembly status == chromosome
            search_term = 'txid%s[Organism:exp] AND ("latest refseq"[filter])'
            print(search_term)
            genome_of_choice = ref_or_not(taxid, search_term)

        if genome_of_choice == None:
            search_term = 'txid%s[Organism:exp] AND ("complete genome"[filter])'#restrict the search term to complete genomes
            print(search_term)
            genome_of_choice = ref_or_not(taxid, search_term)

        acc_list.append(genome_of_choice['acc'])
        name_list.append(genome_of_choice['name'])
        status_list.append(genome_of_choice['status'])
        ref_list.append(genome_of_choice['reference'])
        taxid_list.append(genome_of_choice['taxid'])
        path_list.append(genome_of_choice['path'])

    data = {'AssemblyAccession': acc_list, 'PercentReads': prop_list, 'Name': name_list, 'TaxID': taxid_list,
            'AssemblyStatus': status_list, 'ReferenceStatus': ref_list, 'RefseqPath': path_list}
    df = pd.DataFrame(data)
    return df


'''
Main
'''
Entrez.email = snakemake.params['NCBI_email']

Entrez.api_key = snakemake.params['NCBI_key']

query=list(pd.read_csv(snakemake.input[0],delimiter='\t')['UserInputName'])
print(query)
vrp=snakemake.params['proportion_reads']['virus']
hrp=snakemake.params['proportion_reads']['human']
brp=snakemake.params['proportion_reads']['bacteria']


even_prop_list=generate_even_prop(query,vrp,brp,hrp)['percent_reads']

table=genome_list_to_csv(query,even_prop_list)

'''
Writing the summary table
'''
table.to_csv(snakemake.output[0],sep='\t',header=True)