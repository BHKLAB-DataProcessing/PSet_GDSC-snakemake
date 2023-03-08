from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(
#     access_key_id=config["key"],
#     secret_access_key=config["secret"],
#     host=config["host"],
#     stay_on_remote=False
# )

prefix = config["prefix"]
rna_tool = 'Kallisto-0.46.1'
rna_ref = 'Gencode_v33'
version = "v1"
sens_version = "8.0"
microarray_ver = "u133a"

basePath = "https://orcestradata.blob.core.windows.net/gdsc/GDSC/2019"
rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

data_ver = "gdscv1" if version == "v1" else "gdscv2"

rule get_pset:
    input:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + "download/" + rna_tool_dir + ".tar.gz",
        prefix + "download/GDSC_molecular.zip",
        prefix + "download/" + rna_ref_file,
        prefix + "download/Ensemblv99annotation.RData",
        prefix + "processed/GDSC_" + microarray_ver + "_ENSG.RData",
        prefix + "download/gdsc_mutation_w5.csv",
        prefix + "download/GDSC_rnaseq_meta.txt",
        prefix + "download/mutations_latest.csv",
        prefix + "processed/" + data_ver + "_sens_info_" + sens_version + ".rds",
        prefix + "processed/" + data_ver + "_sens_raw_" + sens_version + ".rds",
        prefix + "processed/profiles_" + sens_version + ".RData",
        prefix + "processed/drugInfo_" + sens_version + ".RData",
        prefix + "processed/cellInfo_" + sens_version + ".RData"
    output:
        prefix + "gdsc.txt"
    shell:
        """
        touch {prefix}gdsc.txt
        """


rule recalculate_and_assemble:
    input:
        prefix + "processed/raw_sense_slices_" + sens_version + ".zip",
    output:
        prefix + "processed/profiles_" + sens_version + ".RData",
    shell:
        """
        Rscript {prefix}scripts/recalculateAndAssembleSlice.R {prefix} {sens_version}
        """

rule process_raw_sens_data:
    input:
        prefix + "download/GDSC1_public_raw_data_17Jul19.csv",
        prefix + "download/GDSC2_public_raw_data_17Jul19.csv",
        prefix + "download/GDSC1_public_raw_data_25Feb20.csv",
        prefix + "download/GDSC2_public_raw_data_25Feb20.csv",
        prefix + "processed/drugInfo_8.0.RData",
        prefix + "processed/drugInfo_8.2.RData",
        prefix + "processed/cellInfo_8.0.RData",
        prefix + "processed/cellInfo_8.2.RData"
    output:
        prefix + "processed/gdscv1_sens_info_8.0.rds",
        prefix + "processed/gdscv1_sens_raw_8.0.rds",
        prefix + "processed/gdscv2_sens_info_8.2.rds",
        prefix + "processed/gdscv2_sens_raw_8.2.rds",
        prefix + "processed/raw_sense_slices_8.0.zip",
        prefix + "processed/raw_sense_slices_8.2.zip"
    shell:
        """
        Rscript {prefix}scripts/processRawSensData.R {prefix} {version}
        """

rule process_gdsc_drugs:
    input:
        prefix + "download/screened_compounds_rel_8.0.csv",
        prefix + "download/screened_compunds_rel_8.2.csv",
        prefix + "download/drugs_with_ids.csv",
    output:
        prefix + "processed/drugInfo_8.0.RData",
        prefix + "processed/drugInfo_8.2.RData"
    shell:
        """
        Rscript {prefix}scripts/processGDSCDrugs.R {prefix}
        """

rule download_gdsc_drugs:
    output:
        prefix + "download/screened_compounds_rel_8.0.csv",
        prefix + "download/screened_compunds_rel_8.2.csv"
    shell:
        """
        Rscript {prefix}scripts/downloadGDSCDrugs.R {prefix}
        """

rule process_gdsc_cells:
    input:
        prefix + "download/gdsc_cellinfo.xlsx",
        prefix + "download/gdsc_cellinfo_8.2.xlsx",
        prefix + "download/cell_annotation_all.csv"
    output:
        prefix + "processed/cellInfo_8.0.RData",
        prefix + "processed/cellInfo_8.2.RData"
    shell:
        """
        Rscript {prefix}scripts/processGDSCCells.R {prefix}
        """

rule download_gdsc_cells:
    output:
        prefix + "download/gdsc_cellinfo.xlsx",
        prefix + "download/gdsc_cellinfo_8.2.xlsx"
    shell:
        """
        Rscript {prefix}scripts/downloadGDSCCells.R {prefix}
        """

rule normalize_microarray:
    input:
        prefix + "brain_array/hgu219hsensgcdf_20.0.0.tar.gz",
        prefix + "brain_array/hgu219hsensgprobe_20.0.0.tar.gz",
        prefix + "brain_array/pd.hgu219.hs.ensg_20.0.0.tar.gz",
        prefix + "brain_array/hthgu133ahsensgcdf_20.0.0.tar.gz",
        prefix + "brain_array/hthgu133ahsensgprobe_20.0.0.tar.gz",
        prefix + "brain_array/pd.hthgu133a.hs.ensg_20.0.0.tar.gz",
        prefix + "brain_array/hgu133ahsensgcdf_20.0.0.tar.gz",
        prefix + "brain_array/hgu133ahsensgprobe_20.0.0.tar.gz",
        prefix + "microarray/gdsc_ge_sampleinfo_" + microarray_ver + ".txt",
        prefix + "microarray/celfile_timestamp_" + microarray_ver + ".csv",
        prefix + "microarray/gdsc_array_" + microarray_ver + ".zip",
        prefix + "microarray/gdsc_ge_sampleinfo_" + microarray_ver + ".txt"
    output:
        prefix + "processed/GDSC_" + microarray_ver + "_ENSG_RAW.RData",
        prefix + "processed/GDSC_" + microarray_ver + "_ENSG.RData"
    shell:
        """
        Rscript {prefix}scripts/normalize_{microarray_ver}.R {prefix}
        """

rule download_microarray:
    output:
        prefix + "microarray/gdsc_ge_sampleinfo_" + microarray_ver + ".txt",
        prefix + "microarray/celfile_timestamp_" + microarray_ver + ".csv",
        prefix + "microarray/gdsc_array_" + microarray_ver + ".zip"
    shell:
        """
        Rscript {prefix}scripts/download_{microarray_ver}.R {prefix}
        """

rule download_brain_array:
    output:
        prefix + "brain_array/hgu219hsensgcdf_20.0.0.tar.gz",
        prefix + "brain_array/hgu219hsensgprobe_20.0.0.tar.gz",
        prefix + "brain_array/pd.hgu219.hs.ensg_20.0.0.tar.gz",
        prefix + "brain_array/hthgu133ahsensgcdf_20.0.0.tar.gz",
        prefix + "brain_array/hthgu133ahsensgprobe_20.0.0.tar.gz",
        prefix + "brain_array/pd.hthgu133a.hs.ensg_20.0.0.tar.gz",
        prefix + "brain_array/hgu133ahsensgcdf_20.0.0.tar.gz",
        prefix + "brain_array/hgu133ahsensgprobe_20.0.0.tar.gz"
    shell:
        """
        Rscript {prefix}scripts/downloadBrainArray.R {prefix}
        """

rule download_gdsc_array:
    output:
        prefix + "download/gdsc_array.zip"
    shell:
        """
        Rscript {prefix}scripts/downloadGDSCArray.R {prefix}
        """

rule download_sensitivity:
    output:
        prefix + "download/GDSC1_public_raw_data_17Jul19.csv",
        prefix + "download/GDSC2_public_raw_data_17Jul19.csv",
        prefix + "download/GDSC1_public_raw_data_25Feb20.csv",
        prefix + "download/GDSC2_public_raw_data_25Feb20.csv"
    shell:
        """
        wget ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.0/GDSC1_public_raw_data_17Jul19.csv \
            -O {prefix}download/GDSC1_public_raw_data_17Jul19.csv
        wget ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.0/GDSC2_public_raw_data_17Jul19.csv \
            -O {prefix}download/GDSC2_public_raw_data_17Jul19.csv
        wget ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.2/GDSC1_public_raw_data_25Feb20.csv \
            -O {prefix}download/GDSC1_public_raw_data_25Feb20.csv
        wget ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.2/GDSC2_public_raw_data_25Feb20.csv \
            -O {prefix}download/GDSC2_public_raw_data_25Feb20.csv
        """

rule download_annotation:
    output:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + "download/" + rna_ref_file,
        prefix + "download/GDSC_rnaseq_meta.txt",
        prefix + "download/Ensemblv99annotation.RData"
    shell:
        """
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv' \
            -O {prefix}download/drugs_with_ids.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/cell_annotation_all.csv' \
            -O {prefix}download/cell_annotation_all.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/{rna_ref_file}' \
            -O {prefix}download/{rna_ref_file}
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/GDSC_rnaseq_meta.txt' \
            -O {prefix}download/GDSC_rnaseq_meta.txt
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/Ensembl.v99.annotation.RData' \
            -O {prefix}download/Ensemblv99annotation.RData
        """

rule download_data:
    output:
        prefix + "download/" + rna_tool_dir + ".tar.gz",
        prefix + "download/CCLE_mutations.csv",
        prefix + "download/sample_info.csv",
        prefix + "download/GDSC_molecular.zip",
        prefix + "download/gdsc_mutation_w5.csv",
        prefix + "download/celline.gdsc.RData",
        prefix + "download/mutations_latest.csv"
    shell:
        """
        wget '{basePath}/RNA-seq/{rna_tool_dir}.tar.gz' -O {prefix}download/{rna_tool_dir}.tar.gz
        wget '{basePath}/GDSC_molecular.zip' -O {prefix}download/GDSC_molecular.zip
        wget '{basePath}/celline.gdsc.RData' -O {prefix}download/celline.gdsc.RData
        wget '{basePath}/mutations_latest.csv' -O {prefix}download/mutations_latest.csv
        wget 'https://ndownloader.figshare.com/files/24613355' -O {prefix}download/CCLE_mutations.csv
        wget 'https://ndownloader.figshare.com/files/24613394' -O {prefix}download/sample_info.csv
        wget 'ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-5.0/gdsc_mutation_w5.csv' -O {prefix}download/gdsc_mutation_w5.csv
        """
