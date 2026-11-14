# All of Us Controlled Tier workflow

The notebooks in `notebooks/` implement the external disease-association analysis:

1. prepare and harmonize VAE/OLS targets;
2. build disease flags and ancestry resources;
3. prepare proteomics;
4. prepare clinical chemistry features;
5. fit locked disease-control and spline-interaction models;
6. summarize the external disease-association and directional-concordance results;
7. summarize demographics and model-specific sample sizes.

## Environment

Inside an authorized All of Us workspace, set:

```bash
export GOOGLE_PROJECT="your-workbench-project"
export WORKSPACE_CDR="your-project.your_cdr_dataset"
export WORKSPACE_ROOT="$PWD/aou_workspace"
export AOU_COMBINED_REPLICATION_ROOT="$PWD/scripts/aou_workbench"
```

Optional file lists use the platform path separator:

```bash
export AOU_PROTEOMICS_FILES="/controlled/path/olink.tsv.gz"
export AOU_PROTEOMICS_MANIFESTS="/controlled/path/olink_manifest.tsv"
export AOU_ANCESTRY_FILES="/controlled/path/ancestry_pcs.tsv"
```

Install the Python dependencies and configure the workspace paths. The
historical filename `src/aou_replication_utils.py` is retained for compatibility;
the analysis is not a same-estimand replication.

The aggregate file `data/aou/replication_summary.tsv` is sufficient
to regenerate Main Figure 6 and Extended Data Figure 4 without Controlled Tier
access. Its filename and fields are retained for compatibility with the frozen
workspace outputs.
