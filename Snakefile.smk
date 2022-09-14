from datetime import datetime
from pathlib import Path

configfile: "[path]/repos/design_specific_effector_primers/config/config_snakefile.yaml"
scripts_directory = config["working_directory"] + "scripts/"
log_directory = config["save_general_results"] + "logs/"
log_versions = config["save_general_results"] + "logs/versions.txt"

domain_gff3_output_path = config["save_general_results"] + "rpsblast_effector_domain/"
domain_gff3_files = ["conserved_sites.gff3", "featdata.txt", "hitdata.txt",
                     "Non-specific.gff3", "Specific.gff3", "Superfamily.gff3"]

pipeline_output_path = config["save_general_results"] + "Pipeline_results/"
pipeline_output_names = ["Effector_primers_general_results", "Effector_primers_mfe_primer",
                         "Effector_primers_primer3", "Effector_primers_amplicons", "Effector_primers_genome_names"]
pipeline_output_extensions = ["tsv", "html"]

run_start = datetime.now()
webserver_dir_name = run_start.strftime("%d_%m_%Y_%H_%M_%S")

rule all:
    input:
        config["save_general_results"] + "Mash_tree.png",
        expand("""{path}{file}""", path=domain_gff3_output_path, file=domain_gff3_files),
        expand("""{path}{name}.{ext}""", path=pipeline_output_path + "html_sub/", name=pipeline_output_names, ext=pipeline_output_extensions),
        config["save_general_results"] + "Amplicon_MSA/" + "regular",
        config["save_general_results"] + "Amplicon_MSA/" + "extended"


# Check if the user provided a data.tsv file.
if config["genome_metadata"]:
    genome_metadata = " --metadata " + config["genome_metadata"]
else:
    genome_metadata = config["genome_metadata"]


rule mash:
    input:
        config["genomes_ingroup"],
        config["genomes_outgroup"]
    conda:
        config["effector_env_packages"]
    params:
        config["save_general_results"],
        scripts_directory + "mash_dendrogram.py",
        config["save_general_results"] + "Input_genome_list.txt",
        config["save_general_results"] + "Mash_distance.mashdist",
        config["save_general_results"] + "Mash_sketch.msh",
        config["save_general_results"] + "Labels.csv",
        config["save_general_results"] + "Labels_ingroup.csv",
        config["genus"],
        scripts_directory + "mash_newick_dendrogram.py"
    output:
        config["save_general_results"] + "Mash_tree.png",
        config["save_general_results"] + "Mash_tree_ete3.png"
    shell:
        """
        mkdir -p {params[0]}
        python3 {params[1]} --ingroup {input[0]} --outgroup {input[1]} --output {params[0]} --genus {params[7]} {genome_metadata}
        python3 {params[8]} --mashdistance {params[3]} --labels {params[5]} --labels_ingroup {params[6]} --output {params[0]}
        rm -f {params[2]} {params[3]} {params[4]} {params[5]} {params[6]}
        version=`mash --version`
        echo mash version $version >> {log_versions}
        """


rule prokka:
    input:
        config["genome"]
    conda:
        config["prokka_env_packages"]
    params:
        config["save_general_results"],
        config["save_general_results"] + "Prokka/",
        config["save_general_results"] + "Prokka/Prokka_out.log"
    log:
        config["log_file_output"] + "Prokka.log"
    output:
        config["save_general_results"] + "Prokka/Prokka_out.faa",
        config["save_general_results"] + "Prokka/Prokka_out.gff",
        config["save_general_results"] + "Prokka/Prokka_out.ffn"
    shell:
        '''
        mkdir -p {params[0]}
        mkdir -p {log_directory}
        [[ -d {log_directory} ]] || mkdir {log_directory} 
        prokka --force --prefix Prokka_out --outdir {params[1]} {input}
        mv {params[2]} {log}
        version=`head -n 1 {log}`
        echo "${{version#*This is }}" >> {log_versions}
        '''


rule t3sepp:
    input:
        config["save_general_results"] + "Prokka/Prokka_out.faa"
    conda:
        config["t3sepp_env_packages"]
    params:
        scripts_directory + "run_T3SEpp.py",
        config["signalp_executable"],
        config["tmhmm_executable"],
        config["psortb_executable"],
        config["t3sepp_executable"],
        config["t3sepp_cutoff"],
        config["save_general_results"]
    log:
        config["log_file_output"] + "T3SEpp.log"
    output:
        config["save_general_results"] + "T3SEpp.out.txt"
    shell:
        # Set the working directory for T3SEpp to the directory that contains the T3SEpp executable.
        # Remove the input file & PSORTb output that is copied to the T3SEpp directory to prevent interweaving of input.
        # T3SEpp version is hardcoded.
        """
        cd {params[4]}
        python3 {params[0]} --fasta {input} --signalp {params[1]} --tmhmm {params[2]} --psortb {params[3]} --t3sepp {params[4]} --cutoff {params[5]} --output {params[6]} | tee {log} 
        FILE={input}
        basename "$FILE"
        f="$(basename -- $FILE)"
        rm {params[4]}"$f"
        sed -i {log} -e 's/\x1b\[[0-9;]*m//g'
        {params[1]} -V >> {log_versions}
        version=`grep 'version' /5_workspace/tools/TMHMM/tmhmm-2.0c/bin/tmhmm | head -n 1`
        echo "${{version##*is }}" >> {log_versions}
        {params[3]} --version >> {log_versions}
        echo T3SEpp version 1.0 >> {log_versions}
        """


rule get_effectors:
    input:
        config["save_general_results"] + "T3SEpp.out.txt",
        config["save_general_results"] + "Prokka/Prokka_out.ffn",
        config["save_general_results"] + "Prokka/Prokka_out.faa"
    conda:
        config["effector_env_packages"]
    params:
        scripts_directory + "T3SEpp_effector_identification.py",
        config["save_general_results"] + "Effector_genes.fasta",
        config["save_general_results"] + "Effector_proteins.fasta"
    log:
        config["log_file_output"] + "effector_genes.log"
    output:
        config["save_general_results"] + "Effector_genes.fasta",
        config["save_general_results"] + "Effector_proteins.fasta"
    shell:
        """
        python3 {params[0]} --t3sepp {input[0]} -ffn {input[1]} -faa {input[2]} --output_genes {params[1]} --output_proteins {params[2]} | tee {log}
        sed -i {log} -e 's/\x1b\[[0-9;]*m//g'
        """


rule rpsblast:
    input:
        config["save_general_results"] + "Effector_proteins.fasta"
    conda:
        config["rpsblast"]
    params:
        config["save_general_results"] + "rpsblast_effector_domain/",
        config["rpsblast_db"],
        config["threads"],
        config["save_general_results"] + "rpsblast_effector_domain/effectors.cdd.rpsblast"
    output:
        config["save_general_results"] + "rpsblast_effector_domain/effectors.cdd.rpsblast"
    shell:
        """
        mkdir -p {params[0]}
        rpsblast -query {input} -db {params[1]} -outfmt 11 -num_threads {params[2]} -out {params[3]}
        rpsblast -version | head -n 1 >> {log_versions}
        """


rule rpsbproc:
    input:
        config["save_general_results"] + "rpsblast_effector_domain/effectors.cdd.rpsblast"
    conda:
        config["rpsblast"]
    params:
        config["rpsbproc_executable"],
        config["save_general_results"] + "rpsblast_effector_domain/effectors.cdd.rpsbproc"
    log:
        config["log_file_output"] + "rpsbproc.log"
    output:
        config["save_general_results"] + "rpsblast_effector_domain/effectors.cdd.rpsbproc"
    shell:
        """
        {params[0]} -i {input} -o {params[1]} | tee {log}
        {params[0]} -version | head -n 1 >> {log_versions}
        """


rule domain_gff3:
    input:
        config["save_general_results"] + "rpsblast_effector_domain/effectors.cdd.rpsbproc",
        config["save_general_results"] + "Prokka/Prokka_out.gff",
        config["save_general_results"] + "Effector_proteins.fasta",
    conda:
        config["rpsblast"]
    params:
        config["domain_gff3"],
        config["save_general_results"] + "rpsblast_effector_domain/",
        config["cddid"]
    log:
        config["log_file_output"] + "domain_gff3.log"
    output:
        expand("""{path}{file}""", path=domain_gff3_output_path, file=domain_gff3_files)
    shell:
        """
        python3 {params[0]} --rpbsproc_file {input[0]} --inputGFF3 {input[1]} --protein {input[2]} --outpath {params[1]} --database {params[2]} | tee {log}
        """


rule primer3:
    input:
        config["save_general_results"] + "Effector_genes.fasta"
    conda:
        config["effector_env_packages"]
    params:
        scripts_directory + "Primer_design_script.py",
        config["working_directory"] + "config/Primer_design_config.tsv",
        config["cut_off_right"],
        config["save_general_results"] + "Primers.tsv",
        config["save_general_results"]
    log:
        config["log_file_output"] + "Primer3.log"
    output:
        config["save_general_results"] + "Primers.tsv"
    shell:
        """
        python3 {params[0]} --input {input} --parameters {params[1]} --cut_off_right {params[2]} --output {params[3]} | tee {log}
        sed -i {log} -e 's/\x1b\[[0-9;]*m//g'
        cp {params[1]} {params[4]}
        """


rule mfeprimer:
    input:
        config["save_general_results"] + "Primers.tsv",
        config["genomes_ingroup"],
        config["genomes_outgroup"]
    conda:
        config["effector_env_packages"]
    params:
        scripts_directory + "run_mfeprimer.py",
        scripts_directory + "Prepare_input_files.py",
        scripts_directory + "Analyze_MFEprimer_output.py",
        config["database_save"] + config["genus"],
        config["mfeprimer_executable"],
        config["save_general_results"] + "MFEprimer/",
        config["threads"],
        config["save_jsons"],
        config["working_directory"] + "config/Acat_control.fasta",
        config["control_genomes"]
    log:
        config["log_file_output"] + "MFEprimer.log"
    output:
        directory(config["save_general_results"] + "MFEprimer/" + "Ingroup_MFEprimer_dataframe_results"),
        directory(config["save_general_results"] + "MFEprimer/" + "Outgroup_MFEprimer_dataframe_results")
    shell:
        """
        [ -d "test" ] && tar -zxf {params[3]}
        mkdir -p {params[3]}
        python3 {params[0]} --input {input[0]} --genome {input[1]} --python_prepare {params[1]} --mfe_analysis {params[2]} --database_save {params[3]} --mfeprimer {params[4]} --output_fa {params[5]} --output_json {params[5]} --threads {params[6]} --save_json {params[7]} --group Ingroup --control {params[8]} --control_genomes {params[9]} | tee {log}
        python3 {params[0]} --input {input[0]} --genome {input[2]} --python_prepare {params[1]} --mfe_analysis {params[2]} --database_save {params[3]} --mfeprimer {params[4]} --output_fa {params[5]} --output_json {params[5]} --threads {params[6]} --save_json {params[7]} --group Outgroup --control {params[8]} --control_genomes {params[9]} | tee -a {log}
        tar -zcf {params[3]}.tar.gz {params[3]} 
        rm -r {params[3]}
        sed -i {log} -e 's/\x1b\[[0-9;]*m//g'
        {params[4]} version >> {log_versions}
        """

# Add the labels flag if the data.tsv file was given and processed by the mash rule.
if Path(config["save_general_results"] + "Labels.json").is_file():
    labels = " --labels " + config["save_general_results"] + "Labels.json"
else:
    labels = ""


rule combine_results:
    input:
        config["save_general_results"] + "T3SEpp.out.txt",
        config["save_general_results"] + "Primers.tsv",
        config["save_general_results"] + "MFEprimer/" + "Ingroup_MFEprimer_dataframe_results",
        config["save_general_results"] + "MFEprimer/" + "Outgroup_MFEprimer_dataframe_results",
        config["save_general_results"] + "Mash_tree_ete3.png"
    conda:
        config["effector_env_packages"]
    params:
        scripts_directory + "Combine_pipeline_results.py",
        pipeline_output_path,
        config["working_directory"] + "src/img/",
        config["webserver"] + config["genus"] + "_" + webserver_dir_name + "/",
        scripts_directory + "Search_for_multiplex_assay_combinations.py",
        labels
    output:
        expand("""{path}{name}.{ext}""", path=pipeline_output_path + "html_sub/", name=pipeline_output_names, ext=pipeline_output_extensions),
        pipeline_output_path + "html_sub/Amplicons_raw_dataframe.tsv"
    log:
        config["log_file_output"] + "Pipeline_results.log"
    shell:
        """
        mkdir -p {params[1]}
        python3 {params[0]} --t3sepp {input[0]} --primer3 {input[1]} --mfe_ingroup {input[2]} --mfe_outgroup {input[3]} {params[5]} --output ./ --multiplex_script {params[4]} | tee {log}
        mv Main_menu.html {params[1]}
        mv html_sub/ {params[1]}
        mkdir {params[1]}html_sub/img
        cp {params[2]}* {params[1]}html_sub/img ; rm {params[1]}html_sub/img/design_specific_effector_primers.drawio.png
        cp -r {params[1]} {params[3]}
        cp {input[4]} {params[3]}html_sub/img
        rm {params[3]}html_sub/*.tsv
        sed -i {log} -e 's/\x1b\[[0-9;]*m//g'
        """


rule clustalw_amplicons:
    input:
        pipeline_output_path + "html_sub/Amplicons_raw_dataframe.tsv",
        pipeline_output_path + "html_sub/Effector_primers_genome_names.tsv"
    conda:
        config["effector_env_packages"]
    params:
        scripts_directory + "run_MSA_for_produced_amplicons.py",
        config["clustalw"],
        config["save_general_results"] + "Amplicon_MSA/"
    output:
        directory(config["save_general_results"] + "Amplicon_MSA/" + "regular"),
        directory(config["save_general_results"] + "Amplicon_MSA/" + "extended")
    log:
        config["log_file_output"] + "MSA_for_amplicons.log"
    shell:
        # ClustalW version is hardcoded because the interface does not allow for easy extraction.
        """
        mkdir -p {params[2]}
        python3 {params[0]} --input {input[0]} -id {input[1]} --clustalw {params[1]} --output {params[2]} | tee {log} || true
        sed -i {log} -e 's/\x1b\[[0-9;]*m//g'
        echo ClustalW version 2.1 >> {log_versions}
        """
