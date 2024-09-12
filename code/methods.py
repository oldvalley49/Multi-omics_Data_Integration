import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
from scipy.io import mmread
import csv

from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams

scglue.plot.set_publication_params()

def load_rna(tissue):
    if tissue == "pbmc":
        rna_counts = mmread("/users/tfurutan/data/glue/rna.mtx")
        rna_counts = rna_counts.transpose().tocsr()
        rna = ad.AnnData(rna_counts)

        #cell names
        file = open("/users/tfurutan/data/glue/cellnames.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        rna.obs_names = [data[i+1][1] for i in range(rna.n_obs)]

        #gene names
        file = open("/users/tfurutan/data/glue/genenames.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        rna.var_names = [data[i+1][1] for i in range(rna.n_vars)]

        #set counts layer
        rna.layers["counts"] = rna.X.copy()
        #set origin
        rna.obs["domain"] = "scRNA-seq"
    else:
        rna_counts = mmread("/fastscratch/myscratch/tfurutan/"+tissue+"/encode_rna.mtx")
        rna_counts = rna_counts.transpose().tocsr()
        rna = ad.AnnData(rna_counts)

        #cell names
        file = open("/fastscratch/myscratch/tfurutan/"+tissue+"/encode_cellnames_rna.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        rna.obs_names = [data[i+1][1] for i in range(rna.n_obs)]

        #gene names
        file = open("/fastscratch/myscratch/tfurutan/"+tissue+"/encode_genenames.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        rna.var_names = [data[i+1][1] for i in range(rna.n_vars)]

        #set counts layer
        rna.layers["counts"] = rna.X.copy()
        #set origin
        rna.obs["domain"] = "scRNA-seq"
    return rna

def rna_preprocess(rna, var_num, gtf_file):
    sc.pp.highly_variable_genes(rna, n_top_genes=var_num, flavor="seurat_v3")
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")
    sc.pp.neighbors(rna, metric="cosine")
    scglue.data.get_gene_annotation(
        rna, gtf=gtf_file,
        gtf_by="gene_name"
    )
    rna = rna[:, rna.var['chrom'].dropna().index]
    return rna

def load_atac(tissue):
    if tissue == "pbmc":
        atac_counts = mmread("/users/tfurutan/data/glue/atac.mtx")
        atac_counts = atac_counts.transpose().tocsr()
        atac = ad.AnnData(atac_counts)

        #cell names
        file = open("/users/tfurutan/data/glue/cellnames_atac.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        atac.obs_names = [data[i+1][1] for i in range(atac.n_obs)]

        #granges names
        file = open("/users/tfurutan/data/glue/granges.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        atac.var_names = [data[i+1][1] for i in range(atac.n_vars)]

        #set origin
        atac.obs["domain"] = "scATAC-seq"
    else:
        atac_counts = mmread("/fastscratch/myscratch/tfurutan/"+tissue+"/encode_atac.mtx")
        atac_counts = atac_counts.transpose().tocsr()
        atac = ad.AnnData(atac_counts)

        #cell names
        file = open("/fastscratch/myscratch/tfurutan/"+tissue+"/encode_cellnames_atac.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        atac.obs_names = [data[i+1][1] for i in range(atac.n_obs)]

        #granges names
        file = open("/fastscratch/myscratch/tfurutan/"+tissue+"/encode_granges.csv", "r")
        data = list(csv.reader(file, delimiter=","))
        file.close
        atac.var_names = [data[i+1][1] for i in range(atac.n_vars)]

        #set origin
        atac.obs["domain"] = "scATAC-seq"
    return atac

def atac_preprocess(atac):
    scglue.data.lsi(atac, n_components=100, n_iter=15)
    sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
    split = atac.var_names.str.split(r"[-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    return atac

def training(rna, atac, var_num, tissue):
    #scglue.config.BEDTOOLS_PATH = "/opt/homebrew/bin"
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    scglue.graph.check_graph(guidance, [rna, atac])
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi"
    )
    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        fit_kws={"directory": "glue"}
    )
    dx = scglue.models.integration_consistency(
        glue, {"rna": rna, "atac": atac}, guidance_hvf
    )
    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    combined = ad.concat([rna, atac])
    feature_embeddings = glue.encode_graph(guidance_hvf)
    feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
    feature_embeddings.iloc[:5, :]
    rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
    atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()
    coembedding = pd.DataFrame(combined.obsm["X_glue"])
    coembedding.columns = ["PC" + str(i) for i in range(1, 51)]
    coembedding.index = combined.obs_names
    filename = "/users/tfurutan/integration/var/coembed/"+tissue+"/glue_v"+str(var_num)+".csv"
    coembedding.to_csv(filename)

def runGlue(tissue, var_num, gtf_file):
    rna = load_rna(tissue)
    atac = load_atac(tissue)
    rna = rna_preprocess(rna, var_num, gtf_file)
    atac = atac_preprocess(atac)
    training(rna, atac, var_num, tissue)
    return None

for var in [2000, 2500, 3000, 4000]:
    runGlue("pancreas", var, "/users/tfurutan/integration_benchmark/annotations.gtf")


#/users/tfurutan/integration_benchmark/annotations.gtf
#/users/tfurutan/data/10x_annotation.gtf