CREATE TABLE bicluster (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    name varchar(20),
    var_exp_fpc float,
    var_exp_fpc_p_value float,
    survival float,
    survival_p_value float,
    PRIMARY KEY (id),
    KEY idx_name (name)
);

CREATE TABLE gene (
    id integer unsigned NOT NULL AUTO_INCREMENT,
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

CREATE TABLE tf_regulator (
    id integer unsigned NOT NULL AUTO_INCREMENT,
    bicluster_id integer unsigned,
    gene_id integer unsigned,
    action varchar(15),
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
    mir2disease varchar(50),
    hmdd binary,
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

