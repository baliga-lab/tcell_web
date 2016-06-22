drop table bic_con, bic_gene, bic_go, discovered_motif;
drop table go_gene, tftg_enrichment, tf_targets;
drop table tf_regulator, tf_fam_gene, tf_family, motif_matches;
drop table motif_column, mirna_regulator, go_bp, exp_cond, bicluster;
drop table known_motif, mirna_targets, mirna, motif;
drop table gene;

CREATE TABLE bicluster (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    name varchar(20),
    var_exp_fpc float,
    var_exp_fpc_p_value float,
    PRIMARY KEY (id),
    KEY idx_name (name)
);

CREATE TABLE gene (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    ucsc varchar(100),
    symbol varchar(100),
    entrez integer,
    PRIMARY KEY (id),
    KEY idx_name (symbol)
);

CREATE TABLE bic_gene (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    bicluster_id integer unsigned,
    gene_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (bicluster_id) REFERENCES bicluster (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE exp_cond (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    name varchar(50),
    PRIMARY KEY (id),
    KEY idx_name (name)
);

CREATE TABLE bic_con (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    bicluster_id integer unsigned,
    exp_cond_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (bicluster_id) REFERENCES bicluster (id),
    FOREIGN KEY (exp_cond_id) REFERENCES exp_cond (id)
);

CREATE TABLE tf_family (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    family_name varchar(100),
    PRIMARY KEY (id)
);

CREATE TABLE tf_fam_gene (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    tf_family_id integer unsigned,
    gene_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (tf_family_id) REFERENCES tf_family (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE motif (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    PRIMARY KEY (id)
);

CREATE TABLE motif_column (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    motif_id integer unsigned,
    col_index integer,
    a float,
    c float,
    g float,
    t float,
    PRIMARY KEY (id),
    FOREIGN KEY (motif_id) REFERENCES motif (id)
);

CREATE TABLE known_motif (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    source_database varchar(50),
    motif_id integer unsigned,
    motif_name varchar(150),
    gene_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (motif_id) REFERENCES motif (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE discovered_motif (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    method varchar(50),
    motif_id integer unsigned,
    motif_name varchar(150),
    bicluster_id integer unsigned,
    score float,
    PRIMARY KEY (id),
    FOREIGN KEY (motif_id) REFERENCES motif (id),
    FOREIGN KEY (bicluster_id) REFERENCES bicluster (id)
);

CREATE TABLE motif_matches (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    discovered_motif_id integer unsigned,
    known_motif_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (known_motif_id) REFERENCES known_motif (id)
);

CREATE TABLE tftg_enrichment (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    bicluster_id integer unsigned,
    known_motif_id integer unsigned,
    p_value float,
    pecTargets float,
    PRIMARY KEY (id),
    FOREIGN KEY (bicluster_id) REFERENCES bicluster (id),
    FOREIGN KEY (known_motif_id) REFERENCES known_motif (id)
);

CREATE TABLE tf_regulator (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    bicluster_id integer unsigned,
    gene_id integer unsigned,
    cor float,
    p_value float,
    ordinal varchar(15),
    PRIMARY KEY (id),
    FOREIGN KEY (bicluster_id) REFERENCES bicluster (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE tf_targets (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    tf_id integer unsigned,
    gene_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (tf_id) REFERENCES gene (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE mirna (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    mature_seq_id varchar(20),
    name varchar(50),
    PRIMARY KEY (id),
    KEY idx_name (name)
);

CREATE TABLE mirna_regulator (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    bicluster_id integer unsigned,
    mirna_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (bicluster_id) REFERENCES bicluster (id),
    FOREIGN KEY (mirna_id) REFERENCES mirna (id)
);

CREATE TABLE mirna_targets (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    mirna_id integer unsigned,
    gene_id integer unsigned,
    source varchar(11),
    PRIMARY KEY (id),
    FOREIGN KEY (mirna_id) REFERENCES mirna (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE go_bp (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    go_id varchar(15),
    name text,
    description text,
    PRIMARY KEY (id),
    KEY idx_go_id (go_id)
);

CREATE TABLE bic_go (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    bicluster_id integer unsigned,
    go_bp_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (bicluster_id) REFERENCES bicluster (id),
    FOREIGN KEY (go_bp_id) REFERENCES go_bp (id)
);

CREATE TABLE go_gene (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    go_bp_id integer unsigned,
    gene_id integer unsigned,
    PRIMARY KEY (id),
    FOREIGN KEY (go_bp_id) REFERENCES go_bp (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);

