# kikawa-seqneut-2025-VCM
Nextstrain analysis of Kikawa et al. neutralization assay experiments for 2025 Southern Hemisphere VCM

## Configure the workflow

The workflow uses sequences that the Bloom lab used in their neutralization experiments which do not have any metadata stored with them.
This workflow uses minimal metadata for each sequence which requires a mapping of Bloom lab strain names to Nextstrain strain names to link sequences and metadata.
This map must be manually curated, but the workflow can generate an initial map to work from.

Build the initial map for each lineage like so.

``` bash
nextstrain build --docker . \
    data/h1n1pdm/strains-kikawa-2025-SH-VCM.tsv \
    data/h3n2/strains-kikawa-2025-SH-VCM.tsv
```

Then, copy those files to the corresponding lineage directory in the `config/` directory.
Manually update the files, as needed to create the map between Nextstrain, Bloom lab, and GISAID strain names.

## Run the workflow

Run the workflow with the default Nextstrain Docker image.

``` bash
nextstrain build --docker .
```

## View the trees

View the trees locally with Auspice.

``` bash
nextstrain view --docker auspice
```

Or view [the H1N1pdm HA tree](https://nextstrain.org/groups/blab/kikawa-seqneut-2025-VCM/h1n1pdm?d=tree,map,frequencies&f_kikawa=present_1&p=grid) or [the H3N2 HA tree](https://nextstrain.org/groups/blab/kikawa-seqneut-2025-VCM/h3n2?d=tree,map,frequencies&f_kikawa=present_1&p=grid) in the "blab" Nextstrain Group.
