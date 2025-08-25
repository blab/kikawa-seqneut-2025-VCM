MIN_DATE = "2024-08-01"
MAX_DATE = "2025-08-22"
MAX_SEQUENCES = 500
LINEAGES = [
    "h1n1pdm",
    "h3n2",
]
GENES = [
    "SigPep",
    "HA1",
    "HA2",
]

NEXTCLADE_DATASET_BY_LINEAGE = {
    "h1n1pdm": "nextstrain/flu/h1n1pdm/ha",
    "h3n2": "flu_h3n2_ha_broad",
}
SUBCLADE_URL_BY_LINEAGE_AND_SEGMENT = {
    "h1n1pdm": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H1N1pdm_HA/main/.auto-generated/subclades.tsv",
    },
    "h3n2": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/main/.auto-generated/subclades.tsv",
    },
}

wildcard_constraints:
    lineage = r'h1n1pdm|h3n2',

rule all:
    input:
        tree_annotations=expand("tables/{lineage}_annotations.tsv", lineage=LINEAGES),
        trees=expand("auspice/kikawa-seqneut-2025-VCM_{lineage}.json", lineage=LINEAGES),
        frequencies=expand("auspice/kikawa-seqneut-2025-VCM_{lineage}_tip-frequencies.json", lineage=LINEAGES),

rule download_kikawa_2025_SH_VCM_strain_sequences:
    output:
        data="data/strains-kikawa-2025-SH-VCM.fasta",
    shell:
        r"""
        curl -L \
            -o {output.data:q} \
            https://raw.githubusercontent.com/jbloomlab/flu-seqneut-2025/refs/heads/main/results/viral_strain_seqs/circulating_2025_HA_ectodomain_nts.fa
        """

rule convert_strain_sequences_to_table:
    input:
        strains="data/strains-kikawa-2025-SH-VCM.fasta",
    output:
        strains="data/strains-kikawa-2025-SH-VCM.tsv",
    shell:
        r"""
        seqkit fx2tab -H --no-qual {input.strains} \
            | csvtk rename -t -C "" -f "#name" -n "bloom_strain" > {output.strains}
        """

rule get_kikawa_2025_SH_VCM_gisaid_strains_by_lineage:
    input:
        data="data/strains-kikawa-2025-SH-VCM.tsv",
    output:
        data="data/gisaid-strains-kikawa-2025-SH-VCM-{lineage}.tsv",
    params:
        lineage=lambda wildcards: {"h3n2": "H3N2", "h1n1pdm": "H1N1"}.get(wildcards.lineage)
    shell:
        r"""
        tsv-select -H -f bloom_strain {input.data} \
            | tsv-filter -H --str-in-fld "bloom_strain:{params.lineage}" \
            | csvtk mutate -T -f bloom_strain -n strain -p "^(.+)_{params.lineage}$" > {output.data}
        """

rule download_parsed_sequences:
    output:
        sequences="data/{lineage}/ha.fasta"
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/trials/ingest/{lineage}/ha/sequences.fasta.xz",
    shell:
        r"""
        aws s3 cp {params.s3_path} - | xz -c -d > {output.sequences}
        """

rule download_parsed_metadata:
    output:
        metadata="data/{lineage}/metadata.tsv",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/trials/ingest/{lineage}/metadata.tsv.xz",
    shell:
        r"""
        aws s3 cp {params.s3_path} - | xz -c -d > {output.metadata}
        """

rule map_nextstrain_to_gisaid_strains_by_lineage:
    input:
        metadata="data/{lineage}/metadata.tsv",
    output:
        strains="data/strain_to_gisaid_strain_{lineage}.tsv",
    shell:
        r"""
        tsv-filter -H --str-ne "passage_category:egg" {input.metadata} \
            | tsv-select -H -f strain,gisaid_strain \
            | csvtk replace -t -f gisaid_strain -p " " -r "_" > {output.strains}
        """

rule join_kikawa_2025_SH_VCM_strains_and_nextstrain_strains:
    input:
        kikawa_strains="data/gisaid-strains-kikawa-2025-SH-VCM-{lineage}.tsv",
        all_strains="data/strain_to_gisaid_strain_{lineage}.tsv",
    output:
        strains="data/strains-kikawa-2025-SH-VCM-{lineage}.tsv",
    shell:
        r"""
        tsv-join -H --filter-file {input.kikawa_strains} --key-fields strain --data-fields gisaid_strain --append-fields bloom_strain {input.all_strains} > {output.strains}
        """

rule get_kikawa_2025_SH_VCM_nextstrain_strains:
    input:
        strains="data/strains-kikawa-2025-SH-VCM-{lineage}.tsv",
    output:
        strains="data/strains-kikawa-2025-SH-VCM-{lineage}.txt",
    shell:
        r"""
        tsv-select -H -f strain {input.strains} \
            | sed 1d \
            | sort -k 1,1 > {output.strains}
        """

rule select_kikawa_sequences_and_metadata:
    input:
        sequences="data/{lineage}/ha.fasta",
        metadata="data/{lineage}/metadata.tsv",
        strains="data/strains-kikawa-2025-SH-VCM-{lineage}.txt",
    output:
        sequences="builds/{lineage}/kikawa_sequences.fasta",
        metadata="builds/{lineage}/kikawa_metadata.tsv",
    shell:
        r"""
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.strains} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """

rule annotate_ids_to_kikawa_metadata:
    input:
        metadata="builds/{lineage}/kikawa_metadata.tsv",
        ids="data/strains-kikawa-2025-SH-VCM-{lineage}.tsv",
    output:
        metadata="builds/{lineage}/kikawa_metadata_with_ids.tsv",
    shell:
        r"""
        augur merge \
            --metadata NEXTSTRAIN={input.metadata} IDS={input.ids} \
            --output-metadata {output.metadata}
        """

rule select_contextual_sequences_and_metadata:
    input:
        sequences="data/{lineage}/ha.fasta",
        metadata="data/{lineage}/metadata.tsv",
        strains="data/strains-kikawa-2025-SH-VCM-{lineage}.txt",
        references="config/{lineage}/references.txt",
        outliers="config/{lineage}/outliers.txt",
    output:
        sequences="builds/{lineage}/contextual_sequences.fasta",
        metadata="builds/{lineage}/contextual_metadata.tsv",
    params:
        min_date=MIN_DATE,
        max_date=MAX_DATE,
        subsample_max_sequences=MAX_SEQUENCES,
    shell:
        r"""
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.strains} {input.outliers} \
            --include {input.references} \
            --exclude-ambiguous-dates-by any \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --group-by region year month \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """

rule merge_metadata_and_sequences:
    input:
        kikawa_metadata="builds/{lineage}/kikawa_metadata_with_ids.tsv",
        kikawa_sequences="builds/{lineage}/kikawa_sequences.fasta",
        contextual_metadata="builds/{lineage}/contextual_metadata.tsv",
        contextual_sequences="builds/{lineage}/contextual_sequences.fasta",
    output:
        sequences="builds/{lineage}/sequences.fasta",
        metadata="builds/{lineage}/metadata.tsv",
    shell:
        r"""
        augur merge \
            --sequences {input.kikawa_sequences} {input.contextual_sequences} \
            --metadata kikawa={input.kikawa_metadata} contextual={input.contextual_metadata} \
            --source-columns "{{NAME}}" \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """

rule get_nextclade_dataset:
    output:
        nextclade_dir=directory("nextclade_dataset/{lineage}_ha/"),
        reference="nextclade_dataset/{lineage}_ha/reference.fasta",
        annotation="nextclade_dataset/{lineage}_ha/genome_annotation.gff3",
    params:
        name=lambda wildcards: NEXTCLADE_DATASET_BY_LINEAGE.get(wildcards.lineage, f"nextstrain/flu/{wildcards.lineage}/ha"),
        nextclade_server_arg=lambda wildcards: f"--server={shquotewords(config['nextclade_server'])}" if config.get("nextclade_server") else "",
    shell:
        r"""
        nextclade dataset get \
            -n '{params.name}' \
            {params.nextclade_server_arg} \
            --output-dir {output.nextclade_dir}
        """

rule align:
    input:
        sequences="builds/{lineage}/sequences.fasta",
        nextclade_dataset="nextclade_dataset/{lineage}_ha/",
    output:
        alignment="builds/{lineage}/aligned.fasta",
        translations=expand("builds/{{lineage}}/translations/{gene}.fasta", gene=GENES),
        annotations="builds/{lineage}/nextclade.tsv",
    params:
        translations_dir="builds/{lineage}/translations",
    threads: 8
    shell:
        r"""
        nextclade run \
            {input.sequences} \
            --jobs {threads} \
            --input-dataset {input.nextclade_dataset} \
            --gap-alignment-side right \
            --include-reference \
            --output-fasta {output.alignment} \
            --output-translations "{params.translations_dir}/{{cds}}.fasta" \
            --output-tsv {output.annotations}
        """

rule merge_nextclade_with_metadata:
    input:
        metadata="builds/{lineage}/metadata.tsv",
        nextclade="builds/{lineage}/nextclade.tsv",
    output:
        merged="builds/{lineage}/metadata_with_nextclade.tsv",
    params:
        metadata_id="strain",
        nextclade_id="seqName",
    shell:
        r"""
        augur merge \
           --metadata \
             metadata={input.metadata} \
             nextclade={input.nextclade} \
           --metadata-id-columns \
             metadata={params.metadata_id} \
             nextclade={params.nextclade_id} \
           --output-metadata {output.merged}
        """

rule tree:
    input:
        alignment="builds/{lineage}/aligned.fasta",
    output:
        tree="builds/{lineage}/tree_raw.nwk",
    threads: 8
    shell:
        r"""
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads} \
            --output {output.tree}
        """

rule get_reference_name:
    input:
        reference="nextclade_dataset/{lineage}_ha/reference.fasta",
    output:
        reference="builds/{lineage}/reference_name.txt",
    shell:
        r"""
        seqkit fx2tab --name --only-id {input.reference} > {output.reference}
        """

def clock_rate(wildcards):
    # these rates are from 12y runs on 2019-10-18
    rate = {
        ('h1n1pdm', 'ha'): 0.00329,
        ('h3n2', 'ha'): 0.00382,
    }
    return rate.get((wildcards.lineage, "ha"), 0.001)

def clock_std_dev(wildcards):
    return clock_rate(wildcards) / 5

rule refine:
    input:
        tree="builds/{lineage}/tree_raw.nwk",
        alignment="builds/{lineage}/aligned.fasta",
        metadata="builds/{lineage}/metadata_with_nextclade.tsv",
        reference="builds/{lineage}/reference_name.txt",
    output:
        tree="builds/{lineage}/tree.nwk",
        node_data="builds/{lineage}/branch_lengths.json",
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = clock_rate,
        clock_std_dev = clock_std_dev,
    shell:
        r"""
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root $(cat {input.reference}) \
            --remove-outgroup \
            --stochastic-resolve \
            --timetree \
            --use-fft \
            --no-covariance \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    input:
        tree="builds/{lineage}/tree.nwk",
        alignment="builds/{lineage}/aligned.fasta",
        translations=expand("builds/{{lineage}}/translations/{gene}.fasta", gene=GENES),
        reference="nextclade_dataset/{lineage}_ha/reference.fasta",
        annotation="nextclade_dataset/{lineage}_ha/genome_annotation.gff3",
    output:
        node_data="builds/{lineage}/muts.json",
        translations=expand("builds/{{lineage}}/translations/{gene}_withInternalNodes.fasta", gene=GENES),
    params:
        inference="joint",
        genes=GENES,
        input_translations="builds/{lineage}/translations/%GENE.fasta",
        output_translations="builds/{lineage}/translations/%GENE_withInternalNodes.fasta",
    shell:
        r"""
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root-sequence {input.reference} \
            --annotation {input.annotation} \
            --genes {params.genes} \
            --translations "{params.input_translations}" \
            --output-node-data {output.node_data} \
            --output-translations "{params.output_translations}" \
            --inference {params.inference}
        """

rule download_subclades:
    output:
        subclades="config/{lineage}/ha/subclades.tsv",
    params:
        url=lambda wildcards: SUBCLADE_URL_BY_LINEAGE_AND_SEGMENT.get(wildcards.lineage, {}).get("ha"),
    shell:
        r"""
        curl -o {output.subclades} "{params.url}"
        """

rule subclades:
    input:
        tree="builds/{lineage}/tree.nwk",
        muts="builds/{lineage}/muts.json",
        subclades="config/{lineage}/ha/subclades.tsv",
    output:
        node_data="builds/{lineage}/subclades.json",
    params:
        membership_name = "subclade",
        label_name = "Subclade",
    shell:
        r"""
        augur clades \
            --tree {input.tree} \
            --mutations {input.muts} \
            --clades {input.subclades} \
            --membership-name {params.membership_name} \
            --label-name {params.label_name} \
            --output {output.node_data}
        """

rule emerging_haplotypes:
    input:
        nextclade="builds/{lineage}/metadata_with_nextclade.tsv",
        haplotypes="config/{lineage}/emerging_haplotypes.tsv",
    output:
        haplotypes_table="builds/{lineage}/emerging_haplotypes.tsv",
        node_data="builds/{lineage}/emerging_haplotypes.json",
    params:
        clade_column="subclade",
        membership_name="emerging_haplotype",
    shell:
        r"""
        python scripts/assign_haplotypes.py \
            --substitutions {input.nextclade:q} \
            --haplotypes {input.haplotypes:q} \
            --clade-column {params.clade_column:q} \
            --haplotype-column-name {params.membership_name:q} \
            --output-table {output.haplotypes_table:q} \
            --output-node-data {output.node_data:q}
        """

rule export:
    input:
        tree="builds/{lineage}/tree.nwk",
        metadata="builds/{lineage}/metadata_with_nextclade.tsv",
        node_data=[
            "builds/{lineage}/branch_lengths.json",
            "builds/{lineage}/muts.json",
            "builds/{lineage}/subclades.json",
            "builds/{lineage}/emerging_haplotypes.json",
            "config/{lineage}/vaccine.json",
        ],
        auspice_config="config/{lineage}/auspice_config.json",
    output:
        auspice_json="auspice/kikawa-seqneut-2025-VCM_{lineage}.json",
    shell:
        r"""
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

rule tip_frequencies:
    input:
        tree="builds/{lineage}/tree.nwk",
        metadata="builds/{lineage}/metadata_with_nextclade.tsv",
    output:
        frequencies="auspice/kikawa-seqneut-2025-VCM_{lineage}_tip-frequencies.json",
    params:
        narrow_bandwidth=1 / 12.0,
        proportion_wide=0.0,
        min_date=MIN_DATE,
        max_date=MAX_DATE,
        pivot_interval=4,
        pivot_interval_units="weeks",
    shell:
        r"""
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --pivot-interval {params.pivot_interval} \
            --pivot-interval-units {params.pivot_interval_units} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output.frequencies}
        """

rule export_tree_annotations:
    input:
        tree="auspice/kikawa-seqneut-2025-VCM_{lineage}.json",
    output:
        tree_annotations="tables/{lineage}_annotations.tsv",
    params:
        attributes=[
            "div",
            "num_date",
            "kikawa",
            "subclade",
        ],
    shell:
        r"""
        python scripts/auspice_tree_to_table.py \
            {input.tree} \
            {output.tree_annotations} \
            --include-internal-nodes \
            --attributes {params.attributes}
        """
