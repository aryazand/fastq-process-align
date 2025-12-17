# FASTQ-PROCESS-ALIGN

This is a clone and modification of [[https://github.com/MPUSP/snakemake-simple-mapping]]. The primary modifcations are

1. Added support for **trim_galore**
2. Added support for **UMIs** 
3. Removed support for variant analysis since I prefer this to be a separate workflow 

## TODO 

- Add in usage information
- Add in DAG information
- Add support for CRAM files
- Make certain intermediate files at temporary

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-simple-mapping
```

Adjust options in the default config file `config/config.yml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```

To run the workflow with test files using **apptainer**:

```bash
snakemake --cores 2 --sdm conda apptainer --directory .test
```