rule filter_seqtab:
    input:
        seqtab = "data/input/{sample}.rds",
        script = "scripts/filter_seqtab.R"
    output:
        filtered_seqtab="results/1-filtering/{sample}_filtering_seqtab.rds",
        removed_seqs="results/1-filtering/{sample}_filtering_removed_seqs.txt"
    params:
        min_abundance = 2,
        min_occurrence = 1,
        min_bp = 300
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
        filtered_seqtab = "results/1-filtering/{sample}_filtering_seqtab.rds",
        script = "scripts/cluster_seqtab.R"
    output:
        clustered_seqtab = "results/2-clustering/{sample}_clustering_seqtab.rds", 
        out_clusters = "results/2-clustering/{sample}_clusters.tsv",  # Correspondence between ASVs and cluster representatives
        removed_seqs="results/2-clustering/{sample}_clustering_removed_seqs.txt"
    params:
        perc_identity = 100,
        min_coverage = 0.9,
        representative_method = "abundance"
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
        clustered_seqtab = "results/2-clustering/{sample}_clustering_seqtab.rds",  
        script = "scripts/seqtab_to_fasta.R"  # The script that generates a FASTA file from seqtab
    output:
        fasta_file = "results/3-hmmsearch/{sample}.fasta",
        aligned_fasta = "results/3-hmmsearch/{sample}.pir"
    shell:
        # Step 1: Generate FASTA from seqtab
        "Rscript {input.script} "
        "--seqtab {input.clustered_seqtab} "
        "--out_fasta {output.fasta_file} && "
        
        # Step 2: Align FASTA using mafft
        "mafft "
        "--auto "
        "--reorder "
        "--thread {threads} "
        "{output.fasta_file} > {output.aligned_fasta}"

rule hmmbuild:
    input:
        aligned_fasta = "results/3-hmmsearch/{sample}.pir"  # The aligned FASTA from the align_asvs rule
    output:
        hmm_profile = "results/3-hmmsearch/{sample}.hmm"  # The generated HMM profile
    shell:
        # Run hmmbuild to generate an HMM profile
        "hmmbuild {output.hmm_profile} {input.aligned_fasta} > /dev/null"

rule hmmsearch:
    input:
        hmm_profile = "results/3-hmmsearch/{sample}.hmm",  # The HMM profile from hmmbuild
        fasta = "results/3-hmmsearch/{sample}.fasta"
    output:
        hmm_search_results = "results/3-hmmsearch/{sample}.hmmsearch"  # The HMM search results
    shell:
        # Run hmmsearch using the generated HMM profile against the FASTA file
        "hmmsearch --tblout {output.hmm_search_results} {input.hmm_profile} {input.fasta} > /dev/null"

rule process_hmmer:
    input:
        seqtab = "results/2-clustering/{sample}_clustering_seqtab.rds",
        hmmer = "results/3-hmmsearch/{sample}.hmmsearch",
        script = "scripts/process_hmmer.R"
    output:
        filtered_seqtab = "results/3-hmmsearch/{sample}_hmmsearch_seqtab.rds",
        removed_seqs = "results/3-hmmsearch/{sample}_hmmsearch_removed_seqs.txt"
    params:
        min_evalue = 1e-5
    shell:
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--hmmer {input.hmmer} "
        "--min_evalue {params.min_evalue} "
        "--out_seqtab {output.filtered_seqtab} "
        "--out_removed_seqs {output.removed_seqs}"

rule make_blast_db:
    input:
        seqtab = "results/3-hmmsearch/{sample}_hmmsearch_seqtab.rds",
        script = "scripts/seqtab_to_fasta.R"
    output:
        fasta = "results/4-internal_gaps/{sample}.fasta",
        blast_db = directory("results/4-internal_gaps/{sample}_blastdb")
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
        blast_db = "results/4-internal_gaps/{sample}_blastdb",  # BLAST database prefix
        fasta = "results/4-internal_gaps/{sample}.fasta"
    output:
        blast_gaps = "results/4-internal_gaps/{sample}_blast_gaps.txt"
    params:
        perc_identity = 95
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
        seqtab = "results/3-hmmsearch/{sample}_hmmsearch_seqtab.rds",
        blast = "results/4-internal_gaps/{sample}_blast_gaps.txt",
        fasta = "results/4-internal_gaps/{sample}.fasta",
        script = "scripts/process_blast_internal_gaps.R"
    output:
        filtered_seqtab = "results/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        removed_seqs = "results/4-internal_gaps/{sample}_internal-gaps_removed_seqs.txt"
    params:
        min_query_coverage = 0.99,  # Minimum fraction of the query sequence that must be aligned for removal consideration (default: 99%).
        min_gaps_subject = 15,  # Minimum length of internal gaps in the subject sequence for removal (default: 15 bases).
        diff_ranks_log = -0.4  # Log-transformed rank difference threshold to remove the subject instead of the query (default: -0.4).
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
        seqtab = "results/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        script = "scripts/seqtab_to_fasta.R"  # The script to generate the FASTA file from seqtab
    output:
        fasta = "results/5-chimeras/{sample}.fasta",
        fasta_trimmed = "results/5-chimeras/{sample}_trimmed_1-180.fasta"
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
        fasta = "results/5-chimeras/{sample}.fasta",
    output:
        blast_db = directory("results/5-chimeras/{sample}_blastdb")
    shell:        
        "makeblastdb -in {input.fasta} -out {output.blast_db}/{wildcards.sample} -dbtype nucl > /dev/null"

rule blast_chimeras_1:
    input:
        fasta_trimmed = "results/5-chimeras/{sample}_trimmed_1-180.fasta",
        blast_db = "results/5-chimeras/{sample}_blastdb"
    output:
        blast_chimeras_1 = "results/5-chimeras/{sample}_blast_chimeras_1.txt"
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
        seqtab = "results/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        blast = "results/5-chimeras/{sample}_blast_chimeras_1.txt",
        fasta = "results/5-chimeras/{sample}.fasta",
        script = "scripts/process_chimeras_1.R"
    output:
        out_list = "results/5-chimeras/{sample}_putative_chimeras.txt"
    shell:
        "Rscript {input.script} "
        "--seqtab {input.seqtab} "
        "--blast {input.blast} "
        "--fasta {input.fasta} "
        "--out_list {output.out_list}"

rule trim_fasta_heads:
    input:
        fasta = "results/5-chimeras/{sample}.fasta",
        putative_chimeras = "results/5-chimeras/{sample}_putative_chimeras.txt"
    output:
        fasta_trimmed = "results/5-chimeras/{sample}_trimmed_241-end.fasta"
    shell:
        "seqkit grep -w 0 -f {input.putative_chimeras} {input.fasta} 2> /dev/null | "
        "seqkit subseq -r 241:-1 > {output.fasta_trimmed} 2> /dev/null"

rule blast_chimeras_2:
    input:
        fasta_trimmed = "results/5-chimeras/{sample}_trimmed_241-end.fasta",
        blast_db = "results/5-chimeras/{sample}_blastdb"
    output:
        blast_chimeras_2 = "results/5-chimeras/{sample}_blast_chimeras_2.txt"
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
        seqtab = "results/4-internal_gaps/{sample}_internal-gaps_seqtab.rds",
        blast1 = "results/5-chimeras/{sample}_blast_chimeras_1.txt",
        blast2 = "results/5-chimeras/{sample}_blast_chimeras_2.txt",
        fasta = "results/5-chimeras/{sample}.fasta",
        script = "scripts/process_chimeras_2.R"
    output:
        final_seqtab = "results/5-chimeras/{sample}_chimeras_seqtab.rds",
        removed_seqs = "results/5-chimeras/{sample}_removed_seqs.txt",
        out_chimeras = "results/5-chimeras/{sample}_chimeras.txt"
    params:
        max_pident_chimera_parent = 99,
        max_pident_parents = 95,
        max_chimera_occurrence = 5
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

rule create_final_dir: # create dir and copy the final seqtab there
    input:
       seqtab = "results/5-chimeras/{sample}_chimeras_seqtab.rds"
    output:
       final_seqtab = "results/final/{sample}_final_seqtab.rds"
    shell:
        "cp {input.seqtab} {output.final_seqtab}"

rule generate_report:
    input:
        original_seqtab = "data/input/{sample}.rds",  # Original sequence table
        final_seqtab = "results/final/{sample}_final_seqtab.rds",  # Final cleaned seqtab
        removed_seqs = [
            "results/1-filtering/{sample}_filtering_removed_seqs.txt",  # Step 1 - Filtering
            "results/2-clustering/{sample}_clustering_removed_seqs.txt",  # Step 2 - Clustering
            "results/3-hmmsearch/{sample}_hmmsearch_removed_seqs.txt",  # Step 3 - HMMER
            "results/4-internal_gaps/{sample}_internal-gaps_removed_seqs.txt",  # Step 4 - Internal Gaps
            "results/5-chimeras/{sample}_removed_seqs.txt"  # Step 5 - Chimeras
        ],
        script = "scripts/final_report.R"  # Path to the final report script
    output:
        out_report = "results/final/{sample}_report.tsv",  # Final report
        out_removed_seqs = "results/final/{sample}_removed_seqs.txt"  # File of all removed sequences
    shell:
        "Rscript {input.script} "
        "--original_seqtab {input.original_seqtab} "
        "--final_seqtab {input.final_seqtab} "
        "--removed_seqs {input.removed_seqs} "
        "--out_report {output.out_report} "
        "--out_removed_seqs {output.out_removed_seqs}"
