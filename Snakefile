from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
rna_tool = 'Kallisto-0.46.1'
rna_ref = 'Gencode_v33'
data_version = 'v1'
sens_version = '8.0'
standardize = True
microarray_ver = "u133a"

basePath = "https://orcestradata.blob.core.windows.net/gdsc/GDSC/2019"
rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

data_ver = "gdscv1" if data_version == "v1" else "gdscv2"
year = "2019" if sens_version == '8.0' else "2020"

rule get_pset:
    input:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + "download/" + rna_tool_dir + ".tar.gz"),
        S3.remote(prefix + "download/GDSC_molecular.zip"),
        S3.remote(prefix + "download/" + rna_ref_file),
        S3.remote(prefix + "download/Ensemblv99annotation.RData"),
        S3.remote(prefix + "processed/GDSC_" + microarray_ver + "_ENSG.RData"),
        S3.remote(prefix + "download/gdsc_mutation_w5.csv"),
        S3.remote(prefix + "download/GDSC_rnaseq_meta.txt"),
        S3.remote(prefix + "download/mutations_latest.csv"),
        S3.remote(prefix + "processed/" + data_ver +
                  "_sens_info_" + sens_version + ".rds"),
        S3.remote(prefix + "processed/" + data_ver +
                  "_sens_raw_" + sens_version + ".rds"),
        S3.remote(prefix + "processed/profiles_" + sens_version + ".RData"),
        S3.remote(prefix + "processed/drugInfo_" + sens_version + ".RData"),
        S3.remote(prefix + "processed/cellInfo_" + sens_version + ".RData")
    output:
        prefix + "GDSC_" + year + \
            "(" + data_version + "-" + sens_version + ").rds"
    shell:
        """
        Rscript scripts/getGDSC.R \
            {prefix} {rna_tool} {rna_ref} {data_version} {sens_version} {standardize} {microarray_ver} \
            cnv mutation microarray fusion
        """

rule recalculate_and_assemble:
    input:
        S3.remote(prefix + "processed/raw_sense_slices_" +
                  sens_version + ".zip"),
    output:
        S3.remote(prefix + "processed/profiles_" + sens_version + ".RData"),
    shell:
        """
        Rscript scripts/recalculateAndAssembleSlice.R {prefix} {sens_version}
        """

rule process_raw_sens_data:
    input:
        S3.remote(prefix + "download/GDSC1_public_raw_data_17Jul19.csv"),
        S3.remote(prefix + "download/GDSC2_public_raw_data_17Jul19.csv"),
        S3.remote(prefix + "download/GDSC1_public_raw_data_25Feb20.csv"),
        S3.remote(prefix + "download/GDSC2_public_raw_data_25Feb20.csv"),
        S3.remote(prefix + "processed/drugInfo_8.0.RData"),
        S3.remote(prefix + "processed/drugInfo_8.2.RData"),
        S3.remote(prefix + "processed/cellInfo_8.0.RData"),
        S3.remote(prefix + "processed/cellInfo_8.2.RData")
    output:
        S3.remote(prefix + "processed/" + data_ver + "_sens_info_8.0.rds"),
        S3.remote(prefix + "processed/" + data_ver + "_sens_raw_8.0.rds"),
        S3.remote(prefix + "processed/" + data_ver + "_sens_info_8.2.rds"),
        S3.remote(prefix + "processed/" + data_ver + "_sens_raw_8.2.rds"),
        S3.remote(prefix + "processed/raw_sense_slices_8.0.zip"),
        S3.remote(prefix + "processed/raw_sense_slices_8.2.zip")
    shell:
        """
        Rscript scripts/processRawSensData.R {prefix} {data_version}
        """

rule process_gdsc_drugs:
    input:
        S3.remote(prefix + "download/screened_compounds_rel_8.0.csv"),
        S3.remote(prefix + "download/screened_compunds_rel_8.2.csv"),
        S3.remote(prefix + "download/drugs_with_ids.csv"),
    output:
        S3.remote(prefix + "processed/drugInfo_8.0.RData"),
        S3.remote(prefix + "processed/drugInfo_8.2.RData")
    shell:
        """
        Rscript scripts/processGDSCDrugs.R {prefix}
        """

rule download_gdsc_drugs:
    output:
        S3.remote(prefix + "download/screened_compounds_rel_8.0.csv"),
        S3.remote(prefix + "download/screened_compunds_rel_8.2.csv")
    shell:
        """
        Rscript scripts/downloadGDSCDrugs.R {prefix}
        """

rule process_gdsc_cells:
    input:
        S3.remote(prefix + "download/gdsc_cellinfo.xlsx"),
        S3.remote(prefix + "download/gdsc_cellinfo_8.2.xlsx"),
        S3.remote(prefix + "download/cell_annotation_all.csv")
    output:
        S3.remote(prefix + "processed/cellInfo_8.0.RData"),
        S3.remote(prefix + "processed/cellInfo_8.2.RData")
    shell:
        """
        Rscript scripts/processGDSCCells.R {prefix}
        """

rule download_gdsc_cells:
    output:
        S3.remote(prefix + "download/gdsc_cellinfo.xlsx"),
        S3.remote(prefix + "download/gdsc_cellinfo_8.2.xlsx")
    shell:
        """
        Rscript scripts/downloadGDSCCells.R {prefix}
        """

rule normalize_microarray:
    input:
        S3.remote(prefix + "download/celline.gdsc.RData"),
        S3.remote(prefix + "download/Ensemblv99annotation.RData"),
        S3.remote(prefix + "brain_array/hgu219hsensgcdf_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hgu219hsensgprobe_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/pd.hgu219.hs.ensg_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hthgu133ahsensgcdf_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hthgu133ahsensgprobe_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/pd.hthgu133a.hs.ensg_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hgu133ahsensgcdf_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hgu133ahsensgprobe_20.0.0.tar.gz"),
        S3.remote(prefix + "microarray/gdsc_ge_sampleinfo_" +
                  microarray_ver + ".txt"),
        S3.remote(prefix + "microarray/celfile_timestamp_" +
                  microarray_ver + ".csv"),
        S3.remote(prefix + "microarray/gdsc_array_" + microarray_ver + ".zip"),
        S3.remote(prefix + "microarray/gdsc_ge_sampleinfo_" +
                  microarray_ver + ".txt")
    output:
        S3.remote(prefix + "processed/GDSC_" +
                  microarray_ver + "_ENSG_RAW.RData"),
        S3.remote(prefix + "processed/GDSC_" + microarray_ver + "_ENSG.RData")
    shell:
        """
        Rscript scripts/normalize_{microarray_ver}.R {prefix}
        """

rule download_microarray:
    output:
        S3.remote(prefix + "microarray/gdsc_ge_sampleinfo_" +
                  microarray_ver + ".txt"),
        S3.remote(prefix + "microarray/celfile_timestamp_" +
                  microarray_ver + ".csv"),
        S3.remote(prefix + "microarray/gdsc_array_" + microarray_ver + ".zip")
    shell:
        """
        Rscript scripts/download_{microarray_ver}.R {prefix}
        """

rule download_brain_array:
    output:
        S3.remote(prefix + "brain_array/hgu219hsensgcdf_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hgu219hsensgprobe_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/pd.hgu219.hs.ensg_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hthgu133ahsensgcdf_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hthgu133ahsensgprobe_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/pd.hthgu133a.hs.ensg_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hgu133ahsensgcdf_20.0.0.tar.gz"),
        S3.remote(prefix + "brain_array/hgu133ahsensgprobe_20.0.0.tar.gz")
    shell:
        """
        Rscript scripts/downloadBrainArray.R {prefix}
        """

rule download_gdsc_array:
    output:
        S3.remote(prefix + "download/gdsc_array.zip")
    shell:
        """
        Rscript scripts/downloadGDSCArray.R {prefix}
        """

rule download_sensitivity:
    output:
        S3.remote(prefix + "download/GDSC1_public_raw_data_17Jul19.csv"),
        S3.remote(prefix + "download/GDSC2_public_raw_data_17Jul19.csv"),
        S3.remote(prefix + "download/GDSC1_public_raw_data_25Feb20.csv"),
        S3.remote(prefix + "download/GDSC2_public_raw_data_25Feb20.csv")
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
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + "download/" + rna_ref_file),
        S3.remote(prefix + "download/GDSC_rnaseq_meta.txt"),
        S3.remote(prefix + "download/Ensemblv99annotation.RData")
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
        S3.remote(prefix + "download/" + rna_tool_dir + ".tar.gz"),
        S3.remote(prefix + "download/CCLE_mutations.csv"),
        S3.remote(prefix + "download/sample_info.csv"),
        S3.remote(prefix + "download/GDSC_molecular.zip"),
        S3.remote(prefix + "download/gdsc_mutation_w5.csv"),
        S3.remote(prefix + "download/celline.gdsc.RData"),
        S3.remote(prefix + "download/mutations_latest.csv")
    shell:
        """
        wget '{basePath}/RNA-seq/{rna_tool_dir}.tar.gz' -O {prefix}download/{rna_tool_dir}.tar.gz
        wget '{basePath}/GDSC_molecular.zip' -O {prefix}download/GDSC_molecular.zip
        wget '{basePath}/celline.gdsc.RData' -O {prefix}download/celline.gdsc.RData
        wget '{basePath}/Mutation/mutations_latest.csv' -O {prefix}download/mutations_latest.csv
        wget 'https://ndownloader.figshare.com/files/24613355' -O {prefix}download/CCLE_mutations.csv
        wget 'https://ndownloader.figshare.com/files/24613394' -O {prefix}download/sample_info.csv
        wget 'ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-5.0/gdsc_mutation_w5.csv' -O {prefix}download/gdsc_mutation_w5.csv
        """
