import sys
import pssm
import MySQLdb
import re, gzip, cPickle
db = MySQLdb.connect(host='localhost', user='tcell', passwd='tcell', db='tcell')

doInserts = True

#################################################################
## Load synonym thesaurus to get UCSC ID to Entrez ID          ##
#################################################################
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def is_symbol(s):
    re1 = re.compile('[A-Z][a-z]+(\d)+')
    if re1.match(s):
        return True
    else:
        return False

ucscConverter = {}
ucsc2symbol = {}
symbol2ucsc = {}
ucsc2entrez = {}
entrez2ucsc = {}
symbol2entrez = {}
entrez2symbol = {}
inFile = open('ucsc2entrez2symbol.csv','r')
#inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    ucscId = splitUp[0]
    entrez = splitUp[1]
    symbol = splitUp[2]
    if not ucscId in ucscConverter:
        ucscConverter[ucscId] = {'symbol':symbol,'entrez':entrez}
    ucsc2symbol[ucscId] = symbol
    if not symbol in symbol2ucsc:
        symbol2ucsc[symbol] = [ucscId]
    else:
        symbol2ucsc[symbol].append(ucscId)
    if not entrez=='NA':
        ucsc2entrez[ucscId] = entrez
        if not symbol=='NA':
            symbol2entrez[symbol] = entrez
            entrez2symbol[entrez] = symbol
            if not entrez in entrez2ucsc:
                entrez2ucsc[entrez] = [ucscId]
            else:
                entrez2ucsc[entrez].append(ucscId)

inFile.close()
print len(ucsc2entrez)

"""
inFile = gzip.open('synonymThesaurus.csv.gz','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    ucscId = splitUp[0]
    symbol = [i for i in splitUp[1].split(';') if is_symbol(i)]
    entrez = [i for i in splitUp[1].split(';') if is_number(i)]
    if not ucscId in ucscConverter:
        ucscConverter[ucscId] = {'symbol':'NA','entrez':'NA'}
    if len(symbol)==1 and ucscConverter[ucscId]['symbol']=='NA':
        ucscConverter[ucscId]['symbol'] = symbol[0]
    if len(entrez)==1 and ucscConverter[ucscId]['entrez']=='NA':
        ucscConverter[ucscId]['entrez'] = entrez[0]
    if len(entrez)==1:
        if not ucscId in ucsc2entrez:
            ucsc2entrez[ucscId] = entrez[0]
            if not entrez[0] in entrez2ucsc:
                entrez2ucsc[entrez[0]] = [ucscId]
            else:
                entrez2ucsc[entrez[0]].append(ucscId)
        if len(symbol)==1:
            if not symbol[0] in symbol2entrez:
                symbol2entrez[symbol[0]] = entrez[0]
                entrez2symbol[entrez[0]] = symbol[0]
    if len(symbol)==1:
        if not ucscId in ucsc2symbol:
            ucsc2symbol[ucscId] = symbol[0]
            if not symbol[0] in entrez2ucsc:
                symbol2ucsc[symbol[0]] = [ucscId]
            else:
                symbol2ucsc[symbol[0]].append(ucscId)

inFile.close()
print len(ucsc2entrez)

### GENE ###
# Load in gene information
genesMMU = {}
inFile = open('Mus_musculus.gene_info','r')
inFile.readline() # Skip header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split('\t')
    genesMMU[splitUp[1]] = splitUp[2]

inFile.close()

# Load in gene information
inFile = open('knownToLocusLink.txt','r')
inFile.readline() # Skip header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split('\t')
    if not splitUp[0] in ucscConverter:
        ucscConverter[splitUp[0]] = {'symbol':splitUp[1],'entrez':'NA'}
    if splitUp[0] in genesMMU and ucscConverter[splitUp[0]]['entrez']=='NA':
        ucscConverter[splitUp[0]]['entrez'] = genesMMU[splitUp[0]]

inFile.close()
"""

def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a)<len(b):
        if a[-3:]=='-3p':
            re1 = re.compile(a+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(a+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(b):
            return 1
    else:
        if b[-3:]=='-3p':
            re1 = re.compile(b+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(b+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(a):
            return 1
    return 0

def miRNAInDict(miRNA, dict1):
    retMe = []
    for i in dict1.keys():
        if compareMiRNANames(miRNA, i):
            retMe.append(miRNAIDs[i])
    return retMe

### SHOW TABLES ###
c = db.cursor()
c.execute("""SHOW TABLES""")
print c.fetchall()

### GENE ###

# Enter in gene information into the database
#(('id', 'int(10) unsigned', 'NO', 'PRI', None, 'auto_increment'), ('symbol', 'varchar(100)', 'YES', 'MUL', None, ''), ('entrez', 'int(11)', 'YES', '', None, ''))
c.execute("""SHOW COLUMNS FROM gene""")
print c.fetchall()

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genesInDb = [i[1] for i in tmp]

done = []
if doInserts:
    for gene in ucscConverter:
        if (not gene in done) and (not gene in genesInDb):
            done.append(gene)
            if not (ucscConverter[gene]['symbol']=='NA' or ucscConverter[gene]['entrez']=='NA'):
                c.execute("""INSERT INTO gene (ucsc, symbol, entrez) VALUES (%s, %s, %s)""", [gene, ucscConverter[gene]['symbol'], ucscConverter[gene]['entrez']])
            elif not ucscConverter[gene]['symbol']=='NA':
                c.execute("""INSERT INTO gene (ucsc, symbol) VALUES (%s,%s)""", [gene, ucscConverter[gene]['symbol']])
            elif not ucscConverter[gene]['entrez']=='NA':
                c.execute("""INSERT INTO gene (ucsc, entrez) VALUES (%s, %s)""", [gene, ucscConverter[gene]['entrez']])

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
print len(done), len(tmp)

# Load up TF families
tfFamilies = {}
inFile = open('tfFamilies_musculus.csv','r')
inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    tfFamilies[splitUp[0]] = splitUp[2].split(' ')

inFile.close()

# Enter in TF family information into the database
c.execute("""SHOW COLUMNS FROM tf_family""")
print c.fetchall()

c.execute("""SELECT * FROM tf_family""")
tmp = c.fetchall()
tfFamInDb = [i[1] for i in tmp]

if doInserts:
    for tfFam in tfFamilies.keys():
        if not tfFam in tfFamInDb:
            c.execute("""INSERT INTO tf_family (family_name) VALUES (%s)""", [tfFam])
    c.connection.commit()

c.execute("""SELECT * FROM tf_family""")
tmp = c.fetchall()
print len(tfFamilies), len(tmp)

# Enter in TF family gene information into the database
c.execute("""SHOW COLUMNS FROM tf_fam_gene""")
print c.fetchall()

c.execute("""SELECT * FROM tf_fam_gene""")
tmp = c.fetchall()
tfFamGeneInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

c.execute("""SELECT * FROM tf_family""")
tmp = c.fetchall()
tfFams = dict(zip([str(i[1]) for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = {}
for i in tmp:
    if not str(i[3]) in genes:
        genes[str(i[3])] = []
    genes[str(i[3])].append(i[0])

if doInserts:
    for tfFam in tfFamilies:
        for gene1 in tfFamilies[tfFam]:
            if gene1 in genes:
                for gene2 in genes[gene1]:
                    if not str(tfFams[tfFam])+'_'+str(gene2) in tfFamGeneInDb:
                        c.execute("""INSERT INTO tf_fam_gene (tf_family_id, gene_id) VALUES (%s,%s)""", [tfFams[tfFam], gene2])
    c.connection.commit()

c.execute("""SELECT * FROM tf_fam_gene""")
tmp = c.fetchall()
print len(tmp)


### exp_condS ###
# Load in exp_cond information
exp_conds = []
inFile = open('conditions.txt','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    exp_conds.append(line.strip())

inFile.close()

# Enter in exp_cond information into the database
c.execute("""SHOW COLUMNS FROM exp_cond""")
print c.fetchall()

c.execute("""SELECT * FROM exp_cond""")
tmp = c.fetchall()
patsInDb = [i[1] for i in tmp]

if doInserts:
    for exp_cond in exp_conds:
        if not exp_cond in patsInDb:
            c.execute("""INSERT INTO exp_cond (name) VALUES (%s)""", [exp_cond])
    c.connection.commit()

c.execute("""SELECT * FROM exp_cond""")
tmp = c.fetchall()
print len(exp_conds), len(tmp)


### GO_BP ###
geneOntology = []
inFile = open('gene_ontology_ext.obo','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    if line.strip()=='[Term]':
        goIds = []
        id1 = inFile.readline().strip().split('id: ')[1]
        name = inFile.readline().strip().split('name: ')[1]
        category = inFile.readline().strip().split('namespace: ')[1]
        if category=='biological_process':
            geneOntology.append([id1, name])
    if line.strip()[0:6]=='alt_id' and category=='biological_process':
        geneOntology.append([line.strip().split('alt_id: ')[1], name])

inFile.close()

geneOntology.append(['GO:0036017','response to erythropoietin'])
geneOntology.append(['GO:0036018','cellular response to erythropoietin'])
geneOntology.append(['GO:0038084','vascular endothelial growth factor signaling pathway'])
geneOntology.append(['GO:0044380','protein localization to cytoskeleton'])
geneOntology.append(['GO:1900046','regulation of hemostasis'])
geneOntology.append(['GO:2001141','regulation of RNA biosynthetic process'])
geneOntology.append(['GO:2001198','regulation of dendritic cell differentiation'])

# Enter in GO BP information into the database
c.execute("""SHOW COLUMNS FROM go_bp""")
print c.fetchall()

c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
gobpInDb = [i[1] for i in tmp]

if doInserts:
    for go in geneOntology:
        if not go[0] in gobpInDb:
            c.execute("""INSERT INTO go_bp (go_id, name, description) VALUES (%s, %s, %s)""", [go[0],go[1],'NA'])
    c.connection.commit()

c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
print len(geneOntology), len(tmp)


### GO_BP TO GENE ###
c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
goTerms = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([str(i[3]) for i in tmp],[i[0] for i in tmp]))

"""
go2gene = []
inFile = open('gene2go.mmu','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split('\t')
    if splitUp[1] in genes and splitUp[2] in goTerms:
        go2gene.append([goTerms[splitUp[2]], genes[splitUp[1]]])

inFile.close()
"""

go2gene = []
inFile = open('entrez_go.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    tmpGO = splitUp[1].split(';')
    for go1 in tmpGO:
        if splitUp[0] in genes and go1 in goTerms:
            go2gene.append([goTerms[go1], genes[splitUp[0]]])

inFile.close()

# Enter in GO BP information into the database
c.execute("""SHOW COLUMNS FROM go_gene""")
print c.fetchall()

c.execute("""SELECT * FROM go_gene""")
tmp = c.fetchall()
goGeneInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

if doInserts:
    for goGene in go2gene:
        if not str(goGene[0])+'_'+str(goGene[1]) in goGeneInDb:
            c.execute("""INSERT INTO go_gene (go_bp_id, gene_id) VALUES (%s, %s)""", goGene)
    c.connection.commit()

c.execute("""SELECT * FROM go_gene""")
tmp = c.fetchall()
print len(go2gene), len(tmp)


### miRNAs ###
# Load in miRNA information
miRNAs = []
miRNAIDs = {}
rev_miRNAIDs = {}
inFile = open('mmu.mature.fa','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    tmp = line.strip().split(' ')
    tmp[0] = re.sub('-5p','',tmp[0])
    miRNAs.append(tmp)
    miRNAIDs[tmp[0].lower()] = tmp[1]
    rev_miRNAIDs[tmp[1]] = tmp[0].lower()

inFile.close()

# Enter in exp_cond information into the database
c.execute("""SHOW COLUMNS FROM mirna""")
print c.fetchall()

c.execute("""SELECT * FROM mirna""")
tmp = c.fetchall()
miRNAInDb = [i[1] for i in tmp]

if doInserts:
    for miRNA in miRNAs:
        if not miRNA[1] in miRNAInDb:
            c.execute("""INSERT INTO mirna (mature_seq_id, name) VALUES (%s, %s)""", [miRNA[1], miRNA[0].lower()])
    c.connection.commit()

c.execute("""SELECT * FROM mirna""")
tmp = c.fetchall()
print len(miRNAs), len(tmp)


### KNOWN MOTIFS ###
knownMotifs = {}
sources = {'jaspar':'jasparCoreVertebrata_redundant.pkl', 'transfac':'transfac_2012.1_PSSMs_vertabrate.pkl', 'uniprobe':'uniprobePSSMsNonRedundant.pkl', 'selex':'selexPSSMsNonRedundant.pkl'}
for i in sources:
    pklFile = open(sources[i],'rb')
    knownMotifs[i] = cPickle.load(pklFile)
    pklFile.close()
    for pssm1 in knownMotifs[i]:
        knownMotifs[i][pssm1].setMethod('meme')

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([i[3] for i in tmp],[i[0] for i in tmp]))

# Load up gene ids for
motif2gene = {}
inFile = open('mouseTFs_All.csv','r')
inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    motif2gene[splitUp[0]] = { 'symbol':splitUp[1].lower().capitalize(), 'entrez':int(splitUp[2]) }
inFile.close()

# Enter in exp_cond information into the database
c.execute("""SHOW COLUMNS FROM known_motif""")
print c.fetchall()

c.execute("""SELECT * FROM known_motif""")
tmp = c.fetchall()
kms = [i[2] for i in tmp]

missingMots2gene = []
missingGene = []
if doInserts:
    for source in knownMotifs:
        for pssm1 in knownMotifs[source]:
            if not knownMotifs[source][pssm1].getName() in kms:
                # Insert motif
                c.execute("""INSERT INTO motif VALUES ()""")
                c.execute("""SELECT LAST_INSERT_ID()""")
                motifId = c.fetchall()[0][0]
                # Insert motif columns
                if knownMotifs[source][pssm1].getName() in motif2gene:
                    if motif2gene[knownMotifs[source][pssm1].getName()]['entrez'] in genes.keys():
                        for i in range(len(knownMotifs[source][pssm1].getMatrix())):
                            c.execute("""INSERT INTO motif_column (motif_id, col_index, a, c, g, t) VALUES (%s, %s, %s, %s, %s, %s)""", [motifId, i, knownMotifs[source][pssm1].getMatrix()[i][0], knownMotifs[source][pssm1].getMatrix()[i][1], knownMotifs[source][pssm1].getMatrix()[i][2], knownMotifs[source][pssm1].getMatrix()[i][3]])
                        c.execute("""INSERT INTO known_motif (source_database, motif_name, motif_id, gene_id) VALUES (%s, %s, %s, %s)""", [source, pssm1, motifId, genes[motif2gene[knownMotifs[source][pssm1].getName()]['entrez']]])
                    else:
                        missingGene.append(pssm1+','+str(motif2gene[knownMotifs[source][pssm1].getName()]['entrez']))
                else:
                    missingMots2gene.append(knownMotifs[source][pssm1].getName())
    #c.connection.commit()

outFile = open('missingMotifs.csv','w')
outFile.write('\n'.join(missingMots2gene))
outFile.close()

outFile = open('missingGenes.csv','w')
outFile.write('\n'.join([str(i) for i in missingGene]))
outFile.close()

c.execute("""SELECT * FROM known_motif""")
tmp = c.fetchall()
print len(tmp)


### BICLUSTER ###
# Bicluster,Number_Of_Genes,Genes,Number_of_Tumors_in_Bicluster,Tumors_in_Bicluster,TCGA_Var_Exp_First_PC,TCGA_Var_Exp_First_PC_Perm_P_Value,TCGA_Survival,TCGA_Survival_P_Value,French_Var_Exp_First_PC,French_Average_Random_Var_Exp_First_PC,French_Var_Exp_First_PC_Perm_P_Value,French_Survival,French_Survival_P_Value,REMBRANDT_Var_Exp_First_PC,REMBRANDT_Average_Random_Var_Exp_First_PC,REMBRANDT_Var_Exp_First_PC_Perm_P_Value,REMBRANDT_Survival,REMBRANDT_Survival_P_Value,GSE7696_Var_Exp_First_PC,GSE7696_Average_Random_Var_Exp_First_PC,GSE7696_Var_Exp_First_PC_Perm_P_Value,GSE7696_Survival,GSE7696_Survival_P_Value,Enriched_GO_BP_Terms,Sustained_Angiogenesis,Insensitivity_To_Antigrowth_Signals,Evading_Apoptosis,Limitless_Replicative_Potential,Evading_Immune_Detection,Tissue_Invasion_And_Metastasis,Self_Sufficiency_In_Growth_Signals,Tumor_Promoting_Inflammation,Reprogramming_Energy_Metabolism,Genome_Instability_And_Mutation
# Load in bicluster information
biclusters = []
inFile = open('biclusters.csv','r')
header = inFile.readline().strip().split(',')
while 1:
    line = inFile.readline()
    if not line:
        break
    biclusters.append(dict(zip(header,line.strip().split(','))))

inFile.close()

# Enter in exp_cond information into the database
c.execute("""SHOW COLUMNS FROM bicluster""")
print c.fetchall()

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
biclustersInDb = [i[1] for i in tmp]

if doInserts:
    for bic1 in biclusters:
        if not bic1['Bicluster'] in biclustersInDb:
            c.execute("""INSERT INTO bicluster (name, var_exp_fpc, var_exp_fpc_p_value) VALUES (%s, %s, %s)""", [bic1['Bicluster'], bic1['Var. Exp. First PC'], bic1['Var. Exp. First PC Perm. P-Value']])
    c.connection.commit()

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
print len(biclusters), len(tmp)


### BICLUSTER GENES ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([str(i[1]) for i in tmp],[i[0] for i in tmp]))

# Enter in exp_cond information into the database
c.execute("""SHOW COLUMNS FROM bic_gene""")
print c.fetchall()

c.execute("""SELECT * FROM bic_gene""")
tmp = c.fetchall()
bicGeneInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
bic_genes = []
missing1 = []
missing2 = []
for bic1 in biclusters:
    for gene in bic1['Genes'].split(' '):
        if gene in genes:
            bic_genes.append([bics[bic1['Bicluster']],genes[gene]])
    missing2 += [i for i in bic1['Genes'].split(' ') if not i in ucscConverter]

#outFile = open('missingGenes.csv','w')
#outFile.write('\n'.join(missing))
#outFile.close()

if doInserts:
    for bg1 in bic_genes:
        if not str(bg1[0])+'_'+str(bg1[1]) in bicGeneInDb:
            c.execute("""INSERT INTO bic_gene (bicluster_id, gene_id) VALUES (%s, %s)""", bg1)
    c.connection.commit()

c.execute("""SELECT * FROM bic_gene""")
tmp = c.fetchall()
print len(bic_genes), len(tmp)


### BICLUSTER exp_cond ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM exp_cond""")
tmp = c.fetchall()
exp_conds = dict(zip([str(i[1]) for i in tmp],[i[0] for i in tmp]))

# Enter in exp_cond information into the database
c.execute("""SHOW COLUMNS FROM bic_con""")
print c.fetchall()

c.execute("""SELECT * FROM bic_con""")
tmp = c.fetchall()
bicPatInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
bic_con = []
for bic1 in biclusters:
    for exp_cond in bic1['Conditions'].split(' '):
        if exp_cond in exp_conds:
            bic_con.append([bics[bic1['Bicluster']],exp_conds[exp_cond]])

if doInserts:
    for bp1 in bic_con:
        if not str(bp1[0])+'_'+str(bp1[1]) in bicPatInDb:
            c.execute("""INSERT INTO bic_con (bicluster_id, exp_cond_id) VALUES (%s, %s)""", bp1)
    c.connection.commit()

c.execute("""SELECT * FROM bic_con""")
tmp = c.fetchall()
print len(bic_con), len(tmp)


### BICLUSTER GO_BP ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
terms = dict(zip([str(i[1]) for i in tmp],[i[0] for i in tmp]))

# Enter in exp_cond information into the database
c.execute("""SHOW COLUMNS FROM bic_go""")
print c.fetchall()

c.execute("""SELECT * FROM bic_go""")
tmp = c.fetchall()
bicGoInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
bic_go = []
missing = []
for bic1 in biclusters:
    if not bic1['GO_Term_BP']=='NA':
        for go in bic1['GO_Term_BP'].split(';'):
            if go in terms:
                bic_go.append([bics[bic1['Bicluster']], terms[go]])
            else:
                #print bic1['Bicluster'], go
                missing.append(go)

#outFile = open('missingTerms.csv','w')
#outFile.write('\n'.join(missing))
#outFile.close()

if doInserts:
    for bg1 in bic_go:
        if not str(bg1[0])+'_'+str(bg1[1]) in bicGoInDb:
            c.execute("""INSERT INTO bic_go (bicluster_id, go_bp_id) VALUES (%s, %s)""", bg1)
    c.connection.commit()

c.execute("""SELECT * FROM bic_go""")
tmp = c.fetchall()
print len(bic_go), len(tmp)


### TF REGULATORS ###
c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

"""
tfRegulators = []
inFile = open('data/tfRegulators_V2.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    tfs = splitUp[1].split(' ')
    for tf in tfs:
        if tf in genes:
            tfRegulators.append([bics[splitUp[0]], genes[tf], splitUp[2]])
        else:
            print tf, splitUp[0]

inFile.close()
"""

c.execute("""SELECT * FROM discovered_motif""")
tmp = c.fetchall()
dms = [i[2] for i in tmp]

# Need to shove motifs from MEME and WEEDER into database
for run1 in ['pita','targetscan','tfbs_db']:
    # MEME
    pklFile = open(run1+'/meme_upstream.pkl','rb')
    memeMots = cPickle.load(pklFile)
    pklFile.close()
    if doInserts:
        for bic1 in memeMots:
            for pssm1 in memeMots[bic1]:
                if not run1+'_'+pssm1.name in dms:
                    # Insert motif
                    c.execute("""INSERT INTO motif VALUES ()""")
                    c.execute("""SELECT LAST_INSERT_ID()""")
                    motifId = c.fetchall()[0]
                    # Insert motif columns
                    for i in range(len(pssm1.matrix)):
                        c.execute("""INSERT INTO motif_column (motif_id, col_index, a, c, g, t) VALUES (%s, %s, %s, %s, %s, %s)""", [motifId, i, pssm1.matrix[i][0], pssm1.matrix[i][1], pssm1.matrix[i][2], pssm1.matrix[i][3]])
                    bicId = bics[run1+'_'+pssm1.name.split('_')[0]]
                    c.execute("""INSERT INTO discovered_motif (method, motif_id, motif_name, bicluster_id, score) VALUES (%s, %s, %s, %s, %s)""", ['MEME', motifId, run1+'_'+pssm1.name, bicId, pssm1.evalue])
        c.connection.commit()

    # WEEDER
    pklFile = open(run1+'/weeder_upstream.pkl','rb')
    weederMots = cPickle.load(pklFile)
    pklFile.close()
    if doInserts:
        for bic1 in weederMots:
            for pssm1 in weederMots[bic1]:
                if not pssm1.name in dms:
                    # Insert motif
                    c.execute("""INSERT INTO motif VALUES ()""")
                    c.execute("""SELECT LAST_INSERT_ID()""")
                    motifId = c.fetchall()[0]
                    # Insert motif columns
                    for i in range(len(pssm1.matrix)):
                        c.execute("""INSERT INTO motif_column (motif_id, col_index, a, c, g, t) VALUES (%s, %s, %s, %s, %s, %s)""", [motifId, i, pssm1.matrix[i][0], pssm1.matrix[i][1], pssm1.matrix[i][2], pssm1.matrix[i][3]])
                    bicId = bics[run1+'_'+pssm1.name.split('_')[0]]
                    c.execute("""INSERT INTO discovered_motif (method, motif_id, motif_name, bicluster_id, score) VALUES (%s, %s, %s, %s, %s)""", ['WEEDER', motifId, run1+'_'+pssm1.name, bicId, pssm1.evalue])
        c.connection.commit()

c.execute("""SELECT * FROM tf_regulator""")
tmp = c.fetchall()
tfregs = [str(i[1])+'_'+str(i[2]) for i in tmp]

# bicluster, TF symbol, method, correlation, p_value
for bic1 in biclusters:
    # MEME
    # Motif 1
    #print bic1['Up.MEME Motif1 Correlated Matches']
    direct = [i.split(':')[0] for i in bic1['Up.MEME Motif1 Expanded Matches'].split(' ')]
    for tfReg in bic1['Up.MEME Motif1 Correlated Matches'].split(' '):
        tmp = tfReg.split(':')
        dirExp = 'Expanded'
        if tmp[0] in genes:
            if tmp[0] in direct:
                dirExp = 'Direct'
            if doInserts:
                bicId = bics[bic1['Bicluster']]
                if not str(bicId)+'_'+str(genes[tmp[0]]) in tfregs:
                   c.execute("""INSERT INTO tf_regulator (bicluster_id, gene_id, cor, p_value, ordinal) VALUES (%s, %s, %s, %s, %s)""", [bicId, genes[tmp[0]], tmp[1], tmp[2], dirExp])

    # Motif 2
    #print bic1['Up.MEME Motif2 Correlated Matches']
    direct = [i.split(':')[0] for i in bic1['Up.MEME Motif2 Expanded Matches'].split(' ')]
    for tfReg in bic1['Up.MEME Motif2 Correlated Matches'].split(' '):
        tmp = tfReg.split(':')
        dirExp = 'Expanded'
        if tmp[0] in genes:
            if tmp[0] in direct:
                dirExp = 'Direct'
            if doInserts:
                bicId = bics[bic1['Bicluster']]
                if not str(bicId)+'_'+str(genes[tmp[0]]) in tfregs:
                    c.execute("""INSERT INTO tf_regulator (bicluster_id, gene_id, cor, p_value, ordinal) VALUES (%s, %s, %s, %s, %s)""", [bicId, genes[tmp[0]], tmp[1], tmp[2], dirExp])

    # WEEDER
    # Motif 1
    #print bic1['Up.WEEDER Motif1 Correlated Matches']
    direct = [i.split(':')[0] for i in bic1['Up.WEEDER Motif1 Expanded Matches'].split(' ')]
    for tfReg in bic1['Up.WEEDER Motif1 Correlated Matches'].split(' '):
        tmp = tfReg.split(':')
        dirExp = 'Expanded'
        if tmp[0] in genes:
            if tmp[0] in direct:
                dirExp = 'Direct'
            if doInserts:
                bicId = bics[bic1['Bicluster']]
                if not str(bicId)+'_'+str(genes[tmp[0]]) in tfregs:
                    c.execute("""INSERT INTO tf_regulator (bicluster_id, gene_id, cor, p_value, ordinal) VALUES (%s, %s, %s, %s, %s)""", [bicId, genes[tmp[0]], tmp[1], tmp[2], dirExp])

    # Motif 2
    #print bic1['Up.WEEDER Motif2 Correlated Matches']
    direct = [i.split(':')[0] for i in bic1['Up.WEEDER Motif2 Expanded Matches'].split(' ')]
    for tfReg in bic1['Up.WEEDER Motif2 Correlated Matches'].split(' '):
        tmp = tfReg.split(':')
        dirExp = 'Expanded'
        if tmp[0] in genes:
            if tmp[0] in direct:
                dirExp = 'Direct'
            if doInserts:
                bicId = bics[bic1['Bicluster']]
                if not str(bicId)+'_'+str(genes[tmp[0]]) in tfregs:
                    c.execute("""INSERT INTO tf_regulator (bicluster_id, gene_id, cor, p_value, ordinal) VALUES (%s, %s, %s, %s, %s)""", [bicId, genes[tmp[0]], tmp[1], tmp[2], dirExp])

    # TFBS_DB
    #print bic1['TFBS_DB.Correlated Matches']
    direct = [i.split(':')[0] for i in bic1['TFBS_DB.Exapnded Matches'].split(' ')]
    for tfReg in bic1['TFBS_DB.Correlated Matches'].split(' '):
        tmp = tfReg.split(':')
        dirExp = 'Expanded'
        if tmp[0] in genes:
            if tmp[0] in direct:
                dirExp = 'Direct'
            if doInserts:
                bicId = bics[bic1['Bicluster']]
                if not str(bicId)+'_'+str(genes[tmp[0]]) in tfregs:
                    c.execute("""INSERT INTO tf_regulator (bicluster_id, gene_id, cor, p_value, ordinal) VALUES (%s, %s, %s, %s, %s)""", [bicId, genes[tmp[0]], tmp[1], tmp[2], dirExp])

c.connection.commit()

c.execute("""SELECT * FROM tf_regulator""")
tmp = c.fetchall()
print len(tmp)
#print len(tfRegulators)

# Put TF target gene interactions into database
tftgDB = {}
inFile = open('tfbsDb_plus_and_minus_5000_entrez.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    if not splitUp[0] in tftgDB:
        tftgDB[splitUp[0]] = []
    tftgDB[splitUp[0]].append(splitUp[1])

inFile.close()

# Load into database
c.execute("""SELECT id, motif_name FROM known_motif""")
tmp = c.fetchall()
knownMotifs = dict(zip([str(i[1]).replace('$','_') for i in tmp],[str(i[0]) for i in tmp]))

c.execute("""SELECT id, entrez FROM gene""")
tmp = c.fetchall()
genes = dict(zip([str(i[1]) for i in tmp],[str(i[0]) for i in tmp]))

c.execute("""SELECT known_motif_id, gene_id FROM tf_targets""")
tmp = c.fetchall()
tfTargets = [str(i[0])+'_'+str(i[1]) for i in tmp]

if doInserts:
    for mot1 in tftgDB:
        if mot1 in knownMotifs:
            for g1 in tftgDB[mot1]:
                if g1 in genes:
                    if not str(knownMotifs[mot1])+'_'+str(genes[g1]) in tfTargets:
                        c.execute("""INSERT INTO tf_targets (known_motif_id, gene_id) VALUES (%s, %s)""", [knownMotifs[mot1], genes[g1]])

c.connection.commit()

c.execute("""SELECT * FROM tf_targets""")
tmp = c.fetchall()
print len(tmp)

### miRNA REGULATORS ###
#c.execute("""SELECT * FROM mirna""")
#tmp = c.fetchall()
#mirnas = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

#c.execute("""SELECT * FROM bicluster""")
#tmp = c.fetchall()
#bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

#mirnaRegulators = []
#inFile = open('data/miRNARegulators.csv','r')
#while 1:
#    line = inFile.readline()
#    if not line:
#        break
#    splitUp = line.strip().split(',')
#    miRs = splitUp[1].lower().split(' ')
#    for miR in miRs:
#        miRs2 = miRNAInDict(miR, miRNAIDs)
#        if not len(miRs2)==0:
#            for miR2 in miRs2:
#                mirnaRegulators.append([bics[splitUp[0]], mirnas[miR2]])
#        else:
#            print miR
#inFile.close()


# Enter in miRNA reglator information into the database
#c.execute("""SHOW COLUMNS FROM mirna_regulator""")
#print c.fetchall()

#c.execute("""SELECT * FROM mirna_regulator""")
#tmp = c.fetchall()
#mirnaRegInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

#if doInserts:
#    for mirna in mirnaRegulators:
#        if not str(miRNA[0])+'_'+str(miRNA[1]) in mirnaRegInDb:
#            c.execute("""INSERT INTO mirna_regulator (bicluster_id, mirna_id) VALUES (%s, %s)""", mirna)
#    c.connection.commit()

#c.execute("""SELECT * FROM mirna_regulator""")
#tmp = c.fetchall()
#print len(mirnaRegulators), len(tmp)

