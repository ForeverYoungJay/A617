# Alloy 617 tomography pipeline

Local source of truth is `code/`; the current Slice 001 files remain at the
project root for backward compatibility with the existing analysis scripts.

On `viper`, the pipeline uses:

```text
~/a617/code/       Python source
~/a617/raw/        downloaded NRDS instrument files (temporary)
~/a617/processed/  lossless ANG + EDS feature tables
~/a617/train/      model inputs/checkpoints (to be added with a defined model)
```

Run from the local project:

```bash
rsync -av code/ viper:~/a617/code/
ssh viper 'bash ~/a617/code/../download_dataset.sh'
ssh viper 'bash ~/a617/code/../process_and_cleanup.sh'
```

The training command is intentionally absent until the target (grain labels,
edge classification, or enrichment regression) and model are specified.
