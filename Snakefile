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
basePath = "https://orcestradata.blob.core.windows.net/gdsc/GDSC/2019"

rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

rule get_pset:
    input:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + "download/" + rna_tool_dir + '.tar.gz',
        prefix + 'download/' + rna_ref_file
    output:
        prefix + "gdsc.txt"
    shell:
        """
        touch {prefix}gdsc.txt
        """

rule download_annotation:
    output:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + 'download/' + rna_ref_file,
        prefix + "download/GDSC_rnaseq_meta.txt",
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
        """

rule download_data:
    output:
        prefix + "download/" + rna_tool_dir + '.tar.gz',
        prefix + "download/CCLE_mutations.csv",
        prefix + "download/sample_info.csv",
        prefix + "download/GDSC_molecular.zip"
    shell:
        """
        wget '{basePath}/RNA-seq/{rna_tool_dir}.tar.gz' -O {prefix}download/{rna_tool_dir}.tar.gz
        wget '{basePath}/GDSC_molecular.zip' -O {prefix}download/GDSC_molecular.zip
        wget 'https://ndownloader.figshare.com/files/24613355' -O {prefix}download/CCLE_mutations.csv
        wget 'https://ndownloader.figshare.com/files/24613394' -O {prefix}download/sample_info.csv
        """
