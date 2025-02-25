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
        # Add new names to avoid problems with BLAST
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
        min_query_coverage = 0.99,
        min_gaps_subject = 15,
        diff_ranks_log = -0.4
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

rule trim_fasta:
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
        "seqkit subseq -r 1:180 {output.fasta} > {output.fasta_trimmed}"

rule blast_chimeras_1:
    input:
        fasta_trimmed = "results/5-chimeras/{sample}_trimmed_1-180.fasta",
        blast_db = "results/4-internal_gaps/{sample}_blastdb"  
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