#!/usr/bin/env python
import sys
import traceback as tb
import logging
import json
import os
import csv
from functools import wraps

import MySQLdb
from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask
from functools import wraps

import math
import itertools
import gzip
import pandas
import numpy as np
import rpy2.robjects as robjects


NUM_PARTS = 5

#convert = {'Evading apoptosis':'cellDeath.gif', 'Evading immune detection':'avoidImmuneDestruction.gif', 'Genome instability and mutation':'genomicInstability.gif', 'Insensitivity to antigrowth signals':'evadeGrowthSuppressors.gif', 'Limitless replicative potential':'immortality.gif', 'Reprogramming energy metabolism':'cellularEnergetics.gif', 'Self sufficiency in growth signals':'sustainedProliferativeSignalling.gif', 'Sustained angiogenesis':'angiogenesis.gif', 'Tissue invasion and metastasis':'invasion.gif', 'Tumor promoting inflammation':'promotingInflammation.gif'}

app = Flask(__name__)
app.config.from_envvar('TCELL_SETTINGS')

######################################################################
#### General helpers
######################################################################

def check_auth(username, password):
    return username == app.config['BASIC_USER'] and password == app.config['BASIC_PASS']

def authenticate():
    return Response('Could not verify your access level for that URL.\n'
            'You have to login with proper credentials', 401,
            {'WWW-Authenticate': 'Basic realm="Login Required"'})

def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password):
            return authenticate()
        return f(*args, **kwargs)
    return decorated

def dbconn():
    return MySQLdb.connect(host=app.config['HOST'], user=app.config['USER'],
                           passwd=app.config['PASS'], db=app.config['DB'])


def read_exps():
    with gzip.open(app.config['GENE_EXPR_FILE'], 'rb') as f:
        return pandas.read_csv(f, sep=',', index_col=0, header=0)

######################################################################
#### Graph/Visualization functionality
######################################################################

GRAPH_COLOR_MAP = {
    'Ag85B_0H':'#a50026',
    'Ag85B_2H':'#d73027',
    'Ag85B_4H':'#f46d43',
    'Ag85B_8H':'#fdae61',
    'ESAT6_0H':'#313695',
    'ESAT6_2H':'#4575b4',
    'ESAT6_4H':'#74add1',
    'ESAT6_8H':'#abd9e9',
}

def phyper(q, m, n, k, lower_tail=False):
    """calls the R function phyper"""
    r_phyper = robjects.r['phyper']
    kwargs = {'lower.tail': lower_tail}
    return float(r_phyper(float(q), float(m), float(n), float(k), **kwargs)[0])


def submat_data(submat, col_indexes):
    """given a sub matrix and a list of column indexes
    that specify the columns, of the matrix, return a list
    of (col_idx, median, min, max, lower_quartile, upper_quartile)
    tuples
    """
    col_medians = np.median(submat, axis=0)
    col_mins = np.min(submat, axis=0)
    col_maxs = np.max(submat, axis=0)
    col_upper_quarts = np.percentile(submat, q=75.0, axis=0)
    col_lower_quarts = np.percentile(submat, q=25.0, axis=0)
    data = [[idx,
             col_mins[i],
             col_lower_quarts[i],
             col_medians[i],
             col_upper_quarts[i],
             col_maxs[i]
             ]
            for i, idx in enumerate(col_indexes)]
    return sorted(data, key=lambda x: x[3])

def subarray_data(subarray, col_indexes):
    """given a sub array and a list of column indexes
    that specify the columns, of the matrix, return a list
    of (col_idx, median, min, max, lower_quartile, upper_quartile)
    tuples
    """
    col_medians = np.median(submat, axis=0)
    col_mins = np.min(submat, axis=0)
    col_maxs = np.max(submat, axis=0)
    col_upper_quarts = np.percentile(submat, q=75.0, axis=0)
    col_lower_quarts = np.percentile(submat, q=25.0, axis=0)
    data = [[idx,
             col_mins[i],
             col_lower_quarts[i],
             col_medians[i],
             col_upper_quarts[i],
             col_maxs[i]
             ]
            for i, idx in enumerate(col_indexes)]
    return sorted(data, key=lambda x: x[3])

def cluster_data_OLD(cursor, cluster_id, df):
    cond_map = {name: index
                   for index, name in enumerate(df.columns.values)}
    gene_map = {name.upper(): index for index, name in enumerate(df.index)}
    cursor.execute("""select g.ucsc from bic_gene bg
join gene g on bg.gene_id=g.id where bicluster_id=%s""", [cluster_id])
    genes = [row[0] for row in cursor.fetchall()]

    gene_indexes = sorted([gene_map[g.upper()] for g in genes])

    Ag85B_indexes = sorted([cond_map[c] for c in ['Ag85B_0H_1','Ag85B_0H_2','Ag85B_0H_3','Ag85B_0H_4','Ag85B_2H_2','Ag85B_2H_3','Ag85B_2H_4','Ag85B_4H_2','Ag85B_4H_3','Ag85B_4H_4','Ag85B_8H_1','Ag85B_8H_3','Ag85B_8H_4']])
    submat1 = df.values[np.ix_(gene_indexes, Ag85B_indexes)]
    Ag85B_data = submat_data(submat1, Ag85B_indexes)
    ESAT6_indexes = sorted([cond_map[c] for c in ['ESAT6_0H_1','ESAT6_0H_2','ESAT6_0H_3','ESAT6_0H_4','ESAT6_2H_1','ESAT6_2H_2','ESAT6_2H_3','ESAT6_2H_4','ESAT6_4H_1','ESAT6_4H_2','ESAT6_4H_3','ESAT6_4H_4','ESAT6_8H_2','ESAT6_8H_3','ESAT6_8H_4']])
    submat2 = df.values[np.ix_(gene_indexes, ESAT6_indexes)]
    ESAT6_data = submat_data(submat2, ESAT6_indexes)
    return Ag85B_data, ESAT6_data

def cluster_data(cursor, cluster_name, df):
    cond_map = {name: index
                   for index, name in enumerate(df.columns.values)}

    bic_map = {name.upper(): index for index, name in enumerate(df.index)}
    bic_index = bic_map[cluster_name.upper()]
    col_indexes = [[cond_map[i] for i in c] for c in [['ESAT6_0H_1','ESAT6_0H_2','ESAT6_0H_3','ESAT6_0H_4'],['ESAT6_2H_1','ESAT6_2H_2','ESAT6_2H_3','ESAT6_2H_4'],['ESAT6_4H_1','ESAT6_4H_2','ESAT6_4H_3','ESAT6_4H_4'],['ESAT6_8H_2','ESAT6_8H_3','ESAT6_8H_4'],['Ag85B_0H_1','Ag85B_0H_2','Ag85B_0H_3','Ag85B_0H_4'],['Ag85B_2H_2','Ag85B_2H_3','Ag85B_2H_4'],['Ag85B_4H_2','Ag85B_4H_3','Ag85B_4H_4'],['Ag85B_8H_1','Ag85B_8H_3','Ag85B_8H_4']]]

    print cond_map,col_indexes
    data = []
    idx = ['ESAT6_0H','ESAT6_2H','ESAT6_4H','ESAT6_8H','Ag85B_0H','Ag85B_2H','Ag85B_4H','Ag85B_8H']
    for i in range(len(col_indexes)):
        ci = col_indexes[i]
        curDF = df.values[np.ix_([bic_index], ci)][0]
        print curDF
        col_median = np.median(curDF)
        col_min = np.min(curDF)
        col_max = np.max(curDF)
        col_upper_quart = np.percentile(curDF, q=75.0)
        col_lower_quart = np.percentile(curDF, q=25.0)
        data.append([idx[i], col_min, col_lower_quart, col_median, col_upper_quart, col_max])
    return data

def subtype_enrichment(cursor, cluster_id, df):
    patient_map = {name: index
                   for index, name in enumerate(df.columns.values)}
    gene_map = {name.upper(): index for index, name in enumerate(df.index)}

    cursor.execute("""select g.symbol, g.ucsc from bic_gene bg
join gene g on bg.gene_id=g.id where bicluster_id=%s""", [cluster_id])
    genes = [row[0] for row in cursor.fetchall()]


    cursor.execute("""select name from bic_pat bp
join patient p on bp.patient_id=p.id where bicluster_id=%s""",
                   [cluster_id])
    included_patients = [row[0] for row in cursor.fetchall()]

    cursor.execute("""select name from patient where id not in
(select patient_id from bic_pat where bicluster_id=%s)""",
                   [cluster_id])
    excluded_patients = [row[0] for row in cursor.fetchall()]

    gene_indexes = sorted([gene_map[g.upper()] for g in genes])

    # above is exactly like boxplot
    # now make a phenotype map
    cursor.execute("""select p.name, pt.name from patient p join phenotypes pt on p.phenotype_id=pt.id
where pt.name <> 'NA'""")
    ptmap = {patient: phenotype for patient, phenotype in cursor.fetchall()}
    all_patients = {patient for patient in ptmap.keys()}
    phenotypes = {phenotype for phenotype in ptmap.values()}

    in_patient_indexes = sorted([patient_map[p] for p in included_patients if p in all_patients])
    ex_patient_indexes = sorted([patient_map[p] for p in excluded_patients if p in all_patients])

    # we use the submat_data function to sort our patients
    in_submat = df.values[np.ix_(gene_indexes, in_patient_indexes)]
    in_data = submat_data(in_submat, in_patient_indexes)
    sorted_in_indexes = [row[0] for row in in_data]
    ex_submat = df.values[np.ix_(gene_indexes, ex_patient_indexes)]
    ex_data = submat_data(ex_submat, ex_patient_indexes)
    sorted_ex_indexes = [row[0] for row in ex_data]

    # sorted by median pValue
    sorted_patient_indexes = sorted_in_indexes + sorted_ex_indexes

    # group patients into phenotype groups.
    # NOTE: the ptmap items need to be sorted, otherwise groupby fails to group correctly
    pt_patients = itertools.groupby(sorted(ptmap.items(), key=lambda pair: pair[1]),
                                    key=lambda pair: pair[1])
    pt_patients = {phenotype: set(map(lambda p: p[0], patients))
                   for phenotype, patients in pt_patients}

    num_columns = len(all_patients)
    cols_per_part = int(math.floor(num_columns / NUM_PARTS))
    #print "# columns: %d # cols/part: %d" % (num_columns, cols_per_part)

    pvalues = []
    min_pvalue = 100.0
    max_pvalue = -100.0

    for i in range(NUM_PARTS):
        part_pvalues = {}
        start = cols_per_part * i
        end = (cols_per_part * (i + 1)) - 1

        # adjust end for the last part
        if i == (NUM_PARTS - 1) and end != (num_columns - 1):
            end = num_columns - 1
        cur_patients = [df.columns.values[p_i] for p_i in sorted_patient_indexes[start:end + 1]]
        #print "Part %d, %d-%d, # current patients: %d" % (i, start, end, len(cur_patients))
        for phenotype in phenotypes:
            q = len([p for p in cur_patients if p in pt_patients[phenotype]])
            k = len(cur_patients)
            m = len(pt_patients[phenotype])
            n = num_columns - m
            pvalue = phyper(q, m, n, k)
            #print "part %d, phenotype: %s, q=%d, k=%d, m=%d, n=%d" % (i, phenotype, q, k, m, n)
            if pvalue == 0.0:
                pvalue = 10e-10

            if pvalue <= 0.5:
                pvalue = -math.log10(2 * pvalue)
            else:
                pvalue = math.log10(2 * (1.0 - pvalue))

            if math.isinf(pvalue):
                signum = -1 if pvalue < 0 else 1
                pvalue = signum * -math.log10(10e-10)

            if pvalue < min_pvalue:
                min_pvalue = pvalue
            if pvalue > max_pvalue:
                max_pvalue = pvalue

            part_pvalues[phenotype] = pvalue
        pvalues.append(part_pvalues)

    return pvalues, min_pvalue, max_pvalue


######################################################################
#### Available application paths
######################################################################

@app.errorhandler(Exception)
def unhandled_exception(e):
    app.logger.exception(e)
    return render_template('unknown_error.html')

@app.route('/')
@requires_auth
def index():
    return render_template('index.html')


def __heatmap_data(conn, bicluster):
    cur = conn.cursor()
    cur.execute("select g.ucsc, g.symbol from gene g join bic_gene bg on bg.gene_id=g.id join bicluster b on bg.bicluster_id=b.id where b.name=%s", [bicluster])
    tmp = cur.fetchall()
    ucsc_ids = [row[0] for row in tmp]
    geneSymbols = []
    for i in tmp:
        if not i[1]=='NA':
            geneSymbols.append(i[1])
        else:
            geneSymbols.append(i[0])
    df = pandas.read_csv(app.config['CMONKEY_EXPR_FILE'], sep=',', index_col=0)
    col_prefixes = sorted(set([name[:-2] for name in df.columns.values]))
    sel = df[df.index.isin(ucsc_ids)]
    target_df = None
    for col_prefix in col_prefixes:
        mysel = sel.filter(regex=col_prefix + "_.*")
        median_df = mysel.median(axis=1).to_frame(name=col_prefix)
        if target_df is None:
            target_df = median_df
        else:
            target_df = target_df.join(median_df)
    heatmap_min = df.values.min()
    heatmap_max = df.values.max()
    absmax = max(abs(heatmap_max), abs(heatmap_min))
    return target_df, geneSymbols, -absmax, absmax

@app.route('/bicluster/<bicluster>')
@requires_auth
def bicluster(bicluster=None):
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id,name,var_exp_fpc,var_exp_fpc_p_value FROM bicluster WHERE name=%s""", [bicluster])
    bc_pk, bc_name, bc_varexp_fpc, bc_varexp_fpc_pval = c.fetchone()
    bic_info = {
        'pk': bc_pk,
        'name': bc_name,
        'varexp_fpc': bc_varexp_fpc,
        'varexp_fpc_pval': bc_varexp_fpc_pval,
        'varexp_flag': bc_varexp_fpc_pval <= 0.05
        }

    c.execute("""SELECT g.id, g.symbol, g.ucsc FROM bic_gene bg join gene g on bg.gene_id=g.id where bg.bicluster_id=%s order by g.symbol""", [bc_pk])
    tmp = list(c.fetchall())
    genes = []
    for i in range(len(tmp)):
        if tmp[i][1]==None:
            genes.append([tmp[i][0],tmp[i][2],tmp[i][2]])
        else:
            genes.append([tmp[i][0],tmp[i][1],tmp[i][2]])
    print genes
    #c.execute("""SELECT * FROM exp_cond ec join bic_con bc on ec.id=bc.exp_cond_id WHERE bc.bicluster_id=%s ORDER BY ec.name""", [bc_pk])
    #conds = list(c.fetchall())

    # Regulators
    elements = []
    elements.append({'data': { 'id': 'bc%d' % bc_pk, 'name': bc_name}, 'classes': 'bicluster' })

    regulators = []
    c.execute("""SELECT g.id, g.symbol, tfr.cor, tfr.p_value, tfr.ordinal FROM tf_regulator tfr join gene g on tfr.gene_id=g.id WHERE tfr.bicluster_id=%s""", [bc_pk])
    tfs = list(c.fetchall())
    tfList = []
    targets = []
    for tf in tfs:
        if not tf[1] in tfList:
            tfList.append(tf[1])
            action = 'Repressor'
            if float(tf[2]) > 0:
                action = 'Activator'
            regulators.append(['TF', tf[0], tf[1], action, tf[2], tf[3], tf[4]])
            elements.append({'data': { 'id': 'reg%d' % tf[0], 'name': tf[1] }, 'classes': 'tf' })
            elements.append({'data': { 'id': 'tfbc%d' % tf[0], 'source': 'reg%d' % tf[0], 'target': 'bc%d' % bc_pk }, 'classes': action.lower() })
        elif tf[4]=='Direct':
            for i in range(len(regulators)):
                if regulators[i][1]==tf[1]:
                    regulators[i][6] = 'Direct'

        # For each tf regulator grab all target genes
        c.execute("""SELECT * FROM known_motif WHERE gene_id=%s""", [tf[0]])
        knownMotifs = list(c.fetchall())
        for mot1 in knownMotifs:
            c.execute("""SELECT * FROM tf_targets, bic_gene WHERE bic_gene.bicluster_id=%s AND tf_targets.known_motif_id=%s AND bic_gene.gene_id=tf_targets.gene_id""", [bc_pk, mot1[0]])
            tmp = c.fetchall()
            if len(tmp)>0:
                targGenes = []
                for i in tmp:
                    c.execute("""SELECT * FROM gene WHERE gene.id=%s""", [i[2]])
                    tmp2 = c.fetchall()[0]
                    if tmp2[1]==None:
                        targGenes.append([tmp2[0],tmp2[2],tmp2[2]])
                    else:
                        targGenes.append([tmp2[0],tmp2[1],tmp2[2]])
                targGenes = sorted(targGenes, key=lambda x: x[2])
                if not [tf[1], mot1[3], ', '.join([i[2] for i in targGenes])] in targets:
                    targets.append([tf[1], mot1[3], len(targGenes), [[i[2],i[1]] for i in targGenes]])
        if tf[4]=='Expanded':
            #
            c.execute("""SELECT tf_family_id FROM tf_fam_gene WHERE gene_id=%s""", [tf[0]])
            tmp = c.fetchall()
            if len(tmp)>0:
                print 'tmp',tmp
                c.execute("""SELECT tfg.gene_id, g.symbol FROM tf_fam_gene as tfg, gene as g WHERE tfg.tf_family_id=%s AND tfg.gene_id!=%s AND tfg.gene_id=g.id""", [tmp[0][0], tf[0]])
                for tmpGene in c.fetchall():
                    c.execute("""SELECT * FROM known_motif WHERE gene_id=%s""", [tmpGene[0]])
                    knownMotifs = list(c.fetchall())
                    for mot1 in knownMotifs:
                        c.execute("""SELECT * FROM tf_targets, bic_gene WHERE bic_gene.bicluster_id=%s AND tf_targets.known_motif_id=%s AND bic_gene.gene_id=tf_targets.gene_id""", [bc_pk, mot1[0]])
                        tmp = c.fetchall()
                        if len(tmp)>0:
                            targGenes = []
                            for i in tmp:
                                c.execute("""SELECT * FROM gene WHERE gene.id=%s""", [i[2]])
                                tmp2 = c.fetchall()[0]
                                if tmp2[1]==None:
                                    targGenes.append([tmp2[0],tmp2[2],tmp2[2]])
                                else:
                                    targGenes.append([tmp2[0],tmp2[1],tmp2[2]])
                            targGenes = sorted(targGenes, key=lambda x: x[2])
                            if not [tf[1]+' ('+tmpGene[1]+')', mot1[3], ', '.join([i[2] for i in targGenes])] in targets:
                                targets.append([tf[1]+' ('+tmpGene[1]+')', mot1[3], len(targGenes), [[i[2],i[1]] for i in targGenes]])

    c.execute("""SELECT mirna.id, mirna.name FROM mirna_regulator mr join mirna on mirna.id=mr.mirna_id WHERE mr.bicluster_id=%s""", [bc_pk])
    mirnas = list(c.fetchall())

    mirnaList = []
    for mirna in mirnas:
        if not mirna[0] in mirnaList:
            regulators.append(['miRNA', mirna[0], mirna[1], 'Repressor', 'Direct'])
            mirnaList.append(mirna[1])
            elements.append({'data': { 'id': 'reg%d' % mirna[0], 'name': mirna[1]}, 'classes': 'mirna' })
            elements.append({'data': { 'id': 'mirnabc%d' % mirna[0], 'source': 'reg%d' % mirna[0], 'target': 'bc%d' % bc_pk }, 'classes': 'repressor' })

    regulators = sorted(regulators, key=lambda x: x[2])

    # GO
    c.execute("""SELECT go_bp.id, go_bp.go_id, go_bp.name FROM bic_go, go_bp WHERE go_bp.id=bic_go.go_bp_id and bic_go.bicluster_id=%s""", [bc_pk])
    tmps = list(c.fetchall())
    gobps = []
    for gobp in tmps:
        c.execute("""SELECT distinct gene.symbol FROM go_gene, gene, bic_gene WHERE go_gene.go_bp_id=%s AND bic_gene.bicluster_id=%s AND go_gene.gene_id=gene.id AND gene.id=bic_gene.gene_id order by gene.symbol""", [gobp[0], bc_pk])
        tmp = c.fetchall()
        if len(tmp)>0:
            gobps.append(list(gobp) + [[row[0] for row in tmp]])

    # Prepare graph plotting data
    exp_data = read_exps()
    js_boxplot_data = cluster_data(c, bic_info['name'], exp_data)
    ratios_mean = np.mean(exp_data.values)
    conditions = ['ESAT6_0H','ESAT6_2H','ESAT6_4H','ESAT6_8H','Ag85B_0H','Ag85B_2H','Ag85B_4H','Ag85B_8H']
    boxplot_colors = [GRAPH_COLOR_MAP[c] for c in conditions]

    # Heatmap
    heatmap_df, geneSymbols, heatmap_min, heatmap_max = __heatmap_data(db, bicluster)
    #heatmap_genes = [g for g in heatmap_df.index]
    heatmap_genes = geneSymbols
    heatmap_values = []
    num_rows = len(heatmap_genes)
    for y in range(num_rows):
        i_y = num_rows - y
        for x in range(heatmap_df.values.shape[1]):
            heatmap_values.append({'x': x, 'y': i_y, 'value': heatmap_df.values[y][x], 'name': heatmap_genes[y]})

    db.close()
    return render_template('bicluster.html', **locals())


@app.route('/search')
@requires_auth
def search():
    gene = request.args.get('gene')
    db = dbconn()
    c = db.cursor()
    type ='gene'
    if not gene:
        return render_template('index.html')
    if gene.find('mmu-')==-1:
        c.execute("""SELECT symbol FROM gene WHERE symbol=%s""", [gene])
        geneData = c.fetchall()
    else:
        c.execute("""SELECT name FROM mirna WHERE name=%s""", [gene])
        geneData = c.fetchall()
        type = 'mirna'
    db.close()

    if len(geneData)==0:
        return render_template('index.html')

    symbol = geneData[0][0]
    if type == 'gene':
        return redirect(url_for('gene', symbol=symbol, regType=type))
    else:
        return redirect(url_for('mirna', symbol=symbol, regType=type))


def __get_muts(c, gene_pk, symbol):
    # Get causal flows downstream of mutation in gene
    muts = {}
    c.execute("""SELECT * FROM somatic_mutation WHERE mutation_type='gene' AND ext_id=%s""", [gene_pk])
    tmp_muts = c.fetchall()
    if len(tmp_muts)==1:
        muts['name'] = symbol
        c.execute("""SELECT * FROM causal_flow WHERE somatic_mutation_id=%s""", [tmp_muts[0][0]])
        tmp_cf = c.fetchall()
        muts['flows'] = 0
        muts['regs'] = []
        muts['tfs'] = []
        muts['miRNAs'] = []
        muts['biclusters'] = []
        muts['data'] = []
        for cf1 in tmp_cf:
            g1 = ''
            if cf1[3]=='tf':
                c.execute("""SELECT * FROM tf_regulator WHERE gene_id=%s AND bicluster_id=%s""", [cf1[2], cf1[4]])
                if len(c.fetchall())>0:
                    c.execute("""SELECT symbol FROM gene WHERE id=%s""", [cf1[2]])
                    g1 = c.fetchall()[0][0]
                    if not g1 in muts['regs']:
                        muts['regs'].append(g1)
                        muts['tfs'].append(g1)
            else:
                c.execute("""SELECT * FROM mirna_regulator WHERE mirna_id=%s AND bicluster_id=%s""", [cf1[2], cf1[4]])
                if len(c.fetchall())>0:
                    c.execute("""SELECT name FROM mirna WHERE id=%s""", [cf1[2]])
                    g1 = c.fetchall()[0][0]
                    if not g1 in muts['regs']:
                        muts['regs'].append(g1)
                        muts['miRNAs'].append(g1)
            if not g1=='':
                c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [cf1[4]])
                b1 = c.fetchall()[0]
                if not b1 in muts['biclusters']:
                    muts['biclusters'].append(b1[0])
                c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id  WHERE bh.bicluster_id=%s",
                          [cf1[4]])

                tmp1 = c.fetchall()
                h1 = list(set([convert[i[0]] for i in tmp1]))
                h2 = [[i[0],convert[i[0]]] for i in tmp1]
                muts['data'].append([symbol, g1, b1[0], b1[1], b1[2], h2])
    return muts


def __get_regulators(c, symbol, regType):
    regs = {}
    tmp_regs = c.fetchall()
    if len(tmp_regs)>0:
        regs['name'] = symbol
        regs['biclusters'] = len(set([i[1] for i in tmp_regs]))
        regs['data'] = []
        # Collect all biclusters downstream regulated by TF or miRNA
        done = []
        for reg in tmp_regs:
            if not reg[1] in done:
                done.append(reg[1])
                action = 'Repressor'
                if regType=='gene' and float(reg[3]) > 0:
                    action = 'Activator'
                c.execute("""SELECT name FROM bicluster WHERE id=%s""", [reg[1]])
                b1 = c.fetchall()[0]
                c.execute("""SELECT go_bp_id FROM bic_go WHERE bic_go.bicluster_id=%s""", [reg[1]])
                bg1 = c.fetchall()
                #c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id WHERE bh.bicluster_id=%s", [reg[1]])
                #tmp1 = c.fetchall()
                #h1 = list(set([convert[i[0]] for i in tmp1]))
                #h2 = [[i[0],convert[i[0]]] for i in tmp1]
                regs['data'].append([symbol, action, b1[0], len(bg1)]) #, b1[1], b1[2]) , h2])
    return regs


@app.route('/mirna')
@app.route('/mirna/<symbol>')
@requires_auth
def mirna(symbol=None, regType=None, defaults={'symbol': None, 'regType': None}):
    # Get biclusters regulated by mirna
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id FROM mirna WHERE name=%s""", [symbol])
    mirna_pk = c.fetchone()[0]
    #muts = __get_muts(c, mirna_pk, symbol)
    c.execute("""SELECT * FROM mirna_regulator WHERE mirna_id=%s""", [mirna_pk])
    regs = __get_regulators(c, symbol, 'mirna')
    db.close()
    return render_template('search.html', gene=symbol, regs=regs, bics={})


@app.route('/gene')
@app.route('/gene/<symbol>')
@requires_auth
def gene(symbol=None, regType=None, defaults={'symbol': None, 'regType': None}):
    db = dbconn()
    c = db.cursor()
    #c.execute("""SELECT id FROM gene WHERE symbol=%s""", [symbol])
    #gene_pk = c.fetchone()[0]
    #muts = __get_muts(c, gene_pk, symbol)
    c.execute("""SELECT * FROM tf_regulator, gene WHERE tf_regulator.gene_id=gene.id AND gene.symbol=%s""", [symbol])
    regs = __get_regulators(c, symbol,'gene')

    # Get biclusters that gene resides
    bics = {}
    #c.execute("SELECT * FROM gene g join bic_gene bg on g.id=bg.gene_id join bicluster b on bg.bicluster_id=b.id where g.symbol=%s", [symbol])
    c.execute("SELECT * FROM bic_gene bg, bicluster b, gene g where g.symbol=%s AND bg.gene_id=g.id AND b.id=bg.bicluster_id", [symbol])
    tmp_bics = c.fetchall()
    if len(tmp_bics) > 0:
        bics['name'] = gene
        bics['biclusters'] = len(tmp_bics)
        bics['data'] = []
        for bic1 in tmp_bics:
            #c.execute("SELECT hm.name FROM bic_hal bh join hallmark hm on bh.hallmark_id=hm.id WHERE bh.bicluster_id=%s", [bic1[3]])
            #tmp1 = c.fetchall()
            #h1 = list(set([convert[i[0]] for i in tmp1]))
            #h2 = [[i[0],convert[i[0]]] for i in tmp1]
            c.execute("SELECT * FROM bic_go WHERE bicluster_id=%s", [bic1[3]])
            tmp1 = c.fetchall()
            bics['data'].append([bic1[4], bic1[5], bic1[6], len(tmp1)]) #, bic1[7], bic1[8]]) #, h2])
    db.close()
    return render_template('search.html', gene=symbol, regs=regs, bics=bics)


@app.route('/network')
@requires_auth
def network():
    return render_template('network.html')


@app.route('/about')
@requires_auth
def about():
    return render_template('about.html')


@app.route('/download')
@requires_auth
def download():
    return render_template('download.html')


@app.route('/citation')
@requires_auth
def citation():
    return render_template('citation.html')


@app.route('/genecompletions')
@requires_auth
def genecompletions():
    term = request.args.get('term')
    db = dbconn()
    try:
        c = db.cursor()
        c.execute("""SELECT symbol FROM gene WHERE symbol LIKE %s""", [str(term)+'%'])
        tmpGene = [i[0] for i in c.fetchall()]
        c.execute("""SELECT name FROM mirna WHERE name LIKE %s""", [str(term)+'%'])
        tmpMiRNA = [i[0] for i in c.fetchall()]
        json1 = json.dumps(list(set(tmpGene))+list(set(tmpMiRNA)))
    finally:
        db.close()
    return Response(response=json1, status=200, mimetype='application/json')


@app.route('/combinatorial_network')
@requires_auth
def combinatorial_network():
    with open(app.config['NODES_FILE'], 'r') as infile:
        csvreader = csv.reader(infile, delimiter=',')
        csvreader.next()
        nodes = {node_id: {'id': node_id, 'tf_ko': tf_ko, 'in_gbm': in_gbm}
                 for node_id, tf_ko, in_gbm in csvreader}

    with open(app.config['EDGES_FILE'], 'r') as infile:
        csvreader = csv.reader(infile, delimiter=',')
        csvreader.next()
        edges = []
        for edge, sig_coocc in csvreader:
            source, edge_type, target = edge.split()
            edges.append({'source': source, 'target': target, 'type': edge_type,
                          'sig_coocc': sig_coocc})

    graph_data = []
    for node_id, node_data in nodes.items():
        classes = []
        if node_id.startswith('hsa-miR'):
            classes.append('mirna')
        else:
            classes.append('gene')

        if node_data['tf_ko'] == 'Yes':
            classes.append('crispr')
        if node_data['in_gbm'] == 'Yes':
            classes.append('in_gbm')
        if 'in_gbm' in classes and 'crispr' in classes:
            classes.append('crispr_gbm')
        graph_data.append({ 'data': { 'id': node_id }, 'classes': ' '.join(classes) })

    for i, edge in enumerate(edges):
        if edge['sig_coocc'] == 'Yes':
            graph_data.append({ 'data': { 'id': 'e%d' % i, 'source': edge['source'], 'target': edge['target'] }, 'classes': 'sigcoocc' })
        else:
            graph_data.append({ 'data': { 'id': 'e%d' % i, 'source': edge['source'], 'target': edge['target'] } })

    return render_template('combinatorial_network.html', **locals())


if __name__ == '__main__':
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    app.debug = True
    app.secret_key = 'supercalifragilistic'
    app.logger.addHandler(handler)
    app.run(host='0.0.0.0', debug=True)

