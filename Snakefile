configfile: "config.yaml"
SEED = config['seed']
SAMPLES = sorted(glob_wildcards("data/input/{sample}.rds").sample)
MERGING = len(SAMPLES) > 1

# Define the final outputs for individual sample processing
sample_outputs = expand("results/datasets/{sample}/{sample}_report.tsv", sample = SAMPLES) + \
                 expand("results/datasets/{sample}/{sample}_final_seqtab.rds", sample = SAMPLES) + \
                 expand("results/datasets/{sample}/{sample}_removed_seqs.txt", sample = SAMPLES)

# If merging is required, include merging outputs
if MERGING:
    merging_outputs = [
        "results/final/clusters.tsv",
        "results/final/final_seqtab.rds",
        "results/final/overall_report.tsv"
    ]
else:
    merging_outputs = []

rule all:
    input:
        sample_outputs + merging_outputs

rule filter_seqtab:
    input:
        seqtab = "data/input/{sample}.rds",
        script = "scripts/filter_seqtab.R"
    output:
        filtered_seqtab="results/datasets/{sample}/1-filtering/{sample}_filtering_seqtab.rds",
        removed_seqs="results/datasets/{sample}/1-filtering/{sample}_filtering_removed_seqs.txt"
    params:
        min_abundance = config['filtering']['min_abundance'],
        min_occurrence = config['filtering']['min_occurrence'],
        min_bp = config['filtering']['min_bp']
    shell:
        "Rscript {input.script} "
            "--seqtab {input.seqtab} "
            "--min_abundance {params.min_abundance} "
            "--min_occurrence {params.min_occurrence} "
            "--min_bp {params.min_bp} "
            "--out_seqtab {output.filtered_seqtab} "
            "--out_removed_seqs {output.removed_seqs}"

rule cluster_seqtab:
    input:
        filtered_seqtab = "results/datasets/{sample}/1-filtering/{sample}_filtering_seqtab.rds",
        script = "scripts/cluster_seqtab.R"
    output:
        clustered_seqtab = "results/datasets/{sample}/2-clustering/{sample}_clustering_seqtab.rds", 
        out_clusters = "results/datasets/{sample}/2-clustering/{sample}_clusters.tsv",  # Correspondence between ASVs and cluster representatives
        removed_seqs="results/datasets/{sample}/2-clustering/{sample}_clustering_removed_seqs.txt"
    params:
        perc_identity = config['clustering']['perc_identity'],
        min_coverage = config['clustering']['min_coverage'],
        representative_method = config['clustering']['representative_method']
    shell:
        "Rscript {input.script} "
            "--seqtab {input.filtered_seqtab} "
            "--perc_identity {params.perc_identity} "
            "--min_coverage {params.min_coverage} "
            "--representative_method {params.representative_method} "
            "--out_seqtab {output.clustered_seqtab} "
            "--out_clusters {output.out_clusters} "
            "--out_removed_seqs {output.removed_seqs}"

rule align_asvs:
    input:
        seqtab = "results/datasets/{sample}/2-clustering/{sample}_clustering_seqtab.rds",  
        script = "scripts/seqtab_to_fasta.R"  # The script that generates a FASTA file from seqtab
    output:
        fasta_file = temp("results/datasets/{sample}/3-hmmsearch/{sample}.fasta"),
        aligned_fasta = temp("results/datasets/{sample}/3-hmmsearch/{sample}.pir")
    shell:
        # Step 1: Generate FASTA from seqtab
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--out_fasta {output.fasta_file} && "
        
        # Step 2: Align FASTA using mafft
        "mafft "
        "--auto "
        "--reorder "
        "--thread {threads} "
        "{output.fasta_file} > {output.aligned_fasta}"

rule hmmbuild:
    input:
        aligned_fasta = "results/datasets/{sample}/3-hmmsearch/{sample}.pir"  # The aligned FASTA from the align_asvs rule
    output:
        hmm_profile = temp("results/datasets/{sample}/3-hmmsearch/{sample}.hmm")  # The generated HMM profile
    shell:
        # Run hmmbuild to generate an HMM profile
        "hmmbuild --seed {SEED} {output.hmm_profile} {input.aligned_fasta} > /dev/null"

rule hmmsearch:
    input:
        hmm_profile = "results/datasets/{sample}/3-hmmsearch/{sample}.hmm",  # The HMM profile from hmmbuild
        fasta = "results/datasets/{sample}/3-hmmsearch/{sample}.fasta"
    output:
        hmm_search_results = temp("results/datasets/{sample}/3-hmmsearch/{sample}.hmmsearch")  # The HMM search results
    shell:
        # Run hmmsearch using the generated HMM profile against the FASTA file
        "hmmsearch --seed {SEED} --tblout {output.hmm_search_results} {input.hmm_profile} {input.fasta} > /dev/null"

rule process_hmmer:
    input:
        seqtab ="results/datasets/{sample}/2-clustering/{sample}_clustering_seqtab.rds",
        hmmer = "results/datasets/{sample}/3-hmmsearch/{sample}.hmmsearch",
        script = "scripts/process_hmmer.R"
    output:
        filtered_seqtab = "results/datasets/{sample}/3-hmmsearch/{sample}_hmmsearch_seqtab.rds",
        removed_seqs = "results/datasets/{sample}/3-hmmsearch/{sample}_hmmsearch_removed_seqs.txt"
    params:
        min_evalue = config['hmmer_search']['min_evalue']
    shell:
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--hmmer {input.hmmer} "
        "--min_evalue {params.min_evalue} "
        "--out_seqtab {output.filtered_seqtab} "
        "--out_removed_seqs {output.removed_seqs}"

rule make_blast_db:
    input:
        seqtab = "results/datasets/{sample}/3-hmmsearch/{sample}_hmmsearch_seqtab.rds",
        script = "scripts/seqtab_to_fasta.R"
    output:
        fasta = temp("results/datasets/{sample}/4-internal_gaps/{sample}.fasta"),
        blast_db = temp(directory("results/datasets/{sample}/4-internal_gaps/{sample}_blastdb"))
    shell:
        # Generate FASTA from the filtered seqtab
        # Add new names to avoid problems with BLAST and order sequences by decreasing abundance (needed to detect internal gaps)
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--add_names "
        "--out_fasta {output.fasta} && "
        
        # Create BLAST database
        "makeblastdb -in {output.fasta} -out {output.blast_db}/{wildcards.sample} -dbtype nucl > /dev/null"

rule blast_internal_gaps:
    input:
        blast_db = "results/datasets/{sample}/4-internal_gaps/{sample}_blastdb",  # BLAST database prefix
        fasta = "results/datasets/{sample}/4-internal_gaps/{sample}.fasta"
    output:
        blast_gaps = temp("results/datasets/{sample}/4-internal_gaps/{sample}_blast_gaps.txt")
    params:
        perc_identity = config['internal_gaps']['blast_perc_identity']
    shell:
        # Run BLAST to identify internal gaps
        "blastn "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' "
        "-perc_identity {params.perc_identity} "
        "-db {input.blast_db}/{wildcards.sample} "
        "-query {input.fasta} "
        "-num_threads {threads} "
        "> {output.blast_gaps}"

rule process_blast_gaps:
    input:
        seqtab = "results/datasets/{sample}/3-hmmsearch/{sample}_hmmsearch_seqtab.rds",
        blast = "results/datasets/{sample}/4-internal_gaps/{sample}_blast_gaps.txt",
        fasta = "results/datasets/{sample}/4-internal_gaps/{sample}.fasta",
        script = "scripts/process_blast_internal_gaps.R"
    output:
        filtered_seqtab = "results/datasets/{sample}/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        removed_seqs = "results/datasets/{sample}/4-internal_gaps/{sample}_internal-gaps_removed_seqs.txt"
    params:
        min_query_coverage = config["internal_gaps"]["min_query_coverage"], # Minimum fraction of the query sequence that must be aligned for removal consideration (default: 99%).
        min_gaps_subject = config["internal_gaps"]["min_gaps_subject"], # Minimum length of internal gaps in the subject sequence for removal (default: 15 bases).
        diff_ranks_log = config["internal_gaps"]["diff_ranks_log"] # Log-transformed rank difference threshold to remove the subject instead of the query (default: -0.4).
    shell:
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--blast {input.blast} "
        "--fasta {input.fasta} "
        "--min_query_coverage {params.min_query_coverage} "
        "--min_gaps_subject {params.min_gaps_subject} "
        "--diff_ranks_log {params.diff_ranks_log} "
        "--out_seqtab {output.filtered_seqtab} "
        "--out_removed_seqs {output.removed_seqs}"

rule trim_fasta_tails:
    input:
        seqtab = "results/datasets/{sample}/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        script = "scripts/seqtab_to_fasta.R"  # The script to generate the FASTA file from seqtab
    output:
        fasta = temp("results/datasets/{sample}/5-chimeras/{sample}.fasta"),
        fasta_trimmed = temp("results/datasets/{sample}/5-chimeras/{sample}_trimmed_1-180.fasta")
    shell:
        # Step 1: Generate FASTA file from seqtab with the -add_names flag
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--out_fasta {output.fasta} "
        "--add_names && "
        
        # Step 2: Trim the FASTA to the first 180 bp using seqkit
        "seqkit subseq -r 1:180 {output.fasta} > {output.fasta_trimmed} 2> /dev/null"

rule make_blast_db_2: # we deleted seqs with gaps after the first db, so we need to redo it with the new fasta (and new names) 
    input:
        fasta = "results/datasets/{sample}/5-chimeras/{sample}.fasta",
    output:
        blast_db = temp(directory("results/datasets/{sample}/5-chimeras/{sample}_blastdb"))
    shell:        
        "makeblastdb -in {input.fasta} -out {output.blast_db}/{wildcards.sample} -dbtype nucl > /dev/null"

rule blast_chimeras_1:
    input:
        fasta_trimmed = "results/datasets/{sample}/5-chimeras/{sample}_trimmed_1-180.fasta",
        blast_db = "results/datasets/{sample}/5-chimeras/{sample}_blastdb"
    output:
        blast_chimeras_1 = temp("results/datasets/{sample}/5-chimeras/{sample}_blast_chimeras_1.txt")
    shell:
        # Perform the BLAST search on the trimmed FASTA file
        "blastn "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' "
        "-perc_identity 100 "
        "-db {input.blast_db}/{wildcards.sample} "
        "-query {input.fasta_trimmed} "
        "-num_threads {threads} "
        "-evalue 1e-50 "
        "> {output.blast_chimeras_1}"

rule process_chimeras_1:
    input:
        seqtab = "results/datasets/{sample}/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        blast = "results/datasets/{sample}/5-chimeras/{sample}_blast_chimeras_1.txt",
        fasta = "results/datasets/{sample}/5-chimeras/{sample}.fasta",
        script = "scripts/process_chimeras_1.R"
    output:
        out_list = temp("results/datasets/{sample}/5-chimeras/{sample}_putative_chimeras.txt")
    shell:
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--blast {input.blast} "
        "--fasta {input.fasta} "
        "--out_list {output.out_list}"

rule trim_fasta_heads:
    input:
        fasta = "results/datasets/{sample}/5-chimeras/{sample}.fasta",
        putative_chimeras = "results/datasets/{sample}/5-chimeras/{sample}_putative_chimeras.txt"
    output:
        fasta_trimmed = temp("results/datasets/{sample}/5-chimeras/{sample}_trimmed_241-end.fasta"),
    shell:
        "seqkit grep -w 0 -f {input.putative_chimeras} {input.fasta} 2> /dev/null | "
        "seqkit subseq -r 241:-1 > {output.fasta_trimmed} 2> /dev/null && "
        "rm {input.fasta}.seqkit.fai" # remove fasta index created by seqkit

rule blast_chimeras_2:
    input:
        fasta_trimmed = "results/datasets/{sample}/5-chimeras/{sample}_trimmed_241-end.fasta",
        blast_db = "results/datasets/{sample}/5-chimeras/{sample}_blastdb"
    output:
        blast_chimeras_2 = temp("results/datasets/{sample}/5-chimeras/{sample}_blast_chimeras_2.txt")
    params:
        perc_identity = 95,
        evalue = 1e-30
    shell:
        "blastn "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' "
        "-perc_identity {params.perc_identity} "
        "-db {input.blast_db}/{wildcards.sample} "
        "-query {input.fasta_trimmed} "
        "-num_threads {threads} "
        "-evalue {params.evalue} "
        "> {output.blast_chimeras_2} 2> /dev/null"

rule process_chimeras_2:
    input:
        seqtab = "results/datasets/{sample}/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        blast1 = "results/datasets/{sample}/5-chimeras/{sample}_blast_chimeras_1.txt",
        blast2 = "results/datasets/{sample}/5-chimeras/{sample}_blast_chimeras_2.txt",
        fasta = "results/datasets/{sample}/5-chimeras/{sample}.fasta",
        script = "scripts/process_chimeras_2.R"
    output:
        final_seqtab = "results/datasets/{sample}/5-chimeras/{sample}_chimeras_seqtab.rds",
        removed_seqs = "results/datasets/{sample}/5-chimeras/{sample}_removed_seqs.txt",
        out_chimeras = "results/datasets/{sample}/5-chimeras/{sample}_chimeras.txt"
    params:
        max_pident_chimera_parent = config["chimeras"]["max_pident_chimera_parent"],
        max_pident_parents = config["chimeras"]["max_pident_parents"],
        max_chimera_occurrence = config["chimeras"]["max_chimera_occurrence"]
    shell:
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--blast1 {input.blast1} "
        "--blast2 {input.blast2} "
        "--fasta {input.fasta} "
        "--out_seqtab {output.final_seqtab} "
        "--out_removed_seqs {output.removed_seqs} "
        "--out_chimeras {output.out_chimeras} "
        "--max_pident_chimera_parent {params.max_pident_chimera_parent} "
        "--max_pident_parents {params.max_pident_parents} "
        "--max_chimera_occurrence {params.max_chimera_occurrence}"

rule create_symlink_final_seqtab: # create symlink of final seqtab to the root of the sample directory
    input:
       seqtab = "results/datasets/{sample}/5-chimeras/{sample}_chimeras_seqtab.rds"
    output:
       final_seqtab = "results/datasets/{sample}/{sample}_final_seqtab.rds"
    shell:
        "ln -s ${{PWD}}/{input.seqtab} {output.final_seqtab}"

rule generate_report:
    input:
        original_seqtab = "data/input/{sample}.rds",  # Original sequence table
        final_seqtab = "results/datasets/{sample}/{sample}_final_seqtab.rds",  # Final cleaned seqtab
        removed_seqs = [
            "results/datasets/{sample}/1-filtering/{sample}_filtering_removed_seqs.txt",
            "results/datasets/{sample}/2-clustering/{sample}_clustering_removed_seqs.txt",
            "results/datasets/{sample}/3-hmmsearch/{sample}_hmmsearch_removed_seqs.txt", 
            "results/datasets/{sample}/4-internal_gaps/{sample}_internal-gaps_removed_seqs.txt",
            "results/datasets/{sample}/5-chimeras/{sample}_removed_seqs.txt" 
        ],
        script = "scripts/final_report.R"  # Path to the final report script
    output:
        out_report = "results/datasets/{sample}/{sample}_report.tsv",  # Final report
        out_removed_seqs = "results/datasets/{sample}/{sample}_removed_seqs.txt"  # File of all removed sequences
    shell:
        "Rscript {input.script} "
        "--original_seqtab {input.original_seqtab} "
        "--final_seqtab {input.final_seqtab} "
        "--removed_seqs {input.removed_seqs} "
        "--out_report {output.out_report} "
        "--out_removed_seqs {output.out_removed_seqs}"

rule merge_seqtabs:
    input:
        seqtabs = expand("results/datasets/{sample}/{sample}_final_seqtab.rds", sample=SAMPLES),
        script = "scripts/merge_seqtabs.R"
    output:
        merged_seqtab = temp("results/final/merged_seqtab.rds")
    shell:
        "Rscript {input.script} "
        "--seqtabs {input.seqtabs} "
        "--out_seqtab {output.merged_seqtab}"

rule cluster_merged_seqtab:
    input:
        merged_seqtab = "results/final/merged_seqtab.rds",
        script = "scripts/cluster_seqtab.R"
    params:
        perc_identity = config['clustering_merged']['perc_identity'],
        min_coverage = config['clustering_merged']['min_coverage'],
        representative_method = config['clustering_merged']['representative_method']
    output:
        clustered_seqtab = "results/final/final_seqtab.rds", 
        out_clusters = "results/final/clusters.tsv"  # Correspondence between ASVs and cluster representatives
    shell:
        "Rscript {input.script} "
            "--seqtab {input.merged_seqtab} "
            "--perc_identity {params.perc_identity} "
            "--min_coverage {params.min_coverage} "
            "--representative_method {params.representative_method} "
            "--out_seqtab {output.clustered_seqtab} "
            "--out_clusters {output.out_clusters} "
            "--out_removed_seqs /dev/null"

rule generate_overall_report:
    input:
        merged_seqtab = "results/final/merged_seqtab.rds",  # Merged seqtab file
        clust_seqtab = "results/final/final_seqtab.rds", # Clustered merged seqtab file
        reports = expand("results/datasets/{sample}/{sample}_report.tsv", sample=SAMPLES)  # Reports from each sample
    output:
        out_report = "results/final/overall_report.tsv",  # Final overall report
    shell:
        "Rscript scripts/overall_report.R "
        "--merged_seqtab {input.merged_seqtab} "
        "--clust_seqtab {input.clust_seqtab} "
        "--reports {input.reports} "
        "--out_report {output.out_report}"