#!/usr/bin/env bash
set -euo pipefail
ssh viper 'mkdir -p ~/a617/code ~/a617/processed ~/a617/train'
rsync -av code/ viper:~/a617/code/
rsync -av server/download_dataset.sh server/process_and_cleanup.sh viper:~/a617/
