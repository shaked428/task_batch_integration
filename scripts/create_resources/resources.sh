#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"
RAW_DATA=resources/common
DATASET_DIR=resources/task_batch_integration

# cat > /tmp/params.yaml << 'HERE'
# input_states: s3://openproblems-data/resources/datasets/**/state.yaml
# rename_keys: 'input:output_dataset'
# output_state: '$id/state.yaml'
# settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad"}'
# publish_dir: s3://openproblems-data/resources/task_batch_integration/datasets/
# HERE



# tw launch https://github.com/openproblems-bio/task_batch_integration.git \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/process_datasets/main.nf \
#   --workspace 53907369739130 \
#   --compute-env 5DwwhQoBi0knMSGcwThnlF \
#   --params-file /tmp/params.yaml \
#   --entry-name auto \
#   --config common/nextflow_helpers/labels_tw.config \
#   --labels task_batch_integration,process_datasets


for dir in "$RAW_DATA"/*; do
  if [ -d "$dir" ] && [ -f "$dir/state.yaml" ] && [ -f "$dir/dataset.h5ad" ]; then
    id=$(basename "$dir")
    
    cat > /tmp/params.yaml << HERE
input_states: $dir/state.yaml
rename_keys: 'input:output_dataset'
output_state: 'state.yaml'
settings: '{"output_dataset": "dataset.h5ad", "output_solution": "solution.h5ad"}'
publish_dir: $DATASET_DIR/$id
HERE

    echo "Processing dataset: $id"

    nextflow run . \
      -main-script target/nextflow/workflows/process_datasets/main.nf \
      -profile singularity \
      -params-file /tmp/params.yaml \
      --input "$dir/dataset.h5ad" \
      --entry-name auto \
      --config common/nextflow_helpers/labels_ci.config \
      --labels task_batch_integration,process_datasets

  else
    echo "Skipping: $dir (missing state.yaml or dataset.h5ad)"
  fi
done


# cat > /tmp/params.yaml << 'HERE'
# input_states: $RAW_DATA/**/state.yaml
# rename_keys: 'input:output_dataset'
# output_state: 'state.yaml'
# settings: '{"output_dataset": "dataset.h5ad", "output_solution": "solution.h5ad"}'
# publish_dir: $DATASET_DIR/$id
# HERE

# nextflow run . \
#   -main-script target/nextflow/workflows/process_datasets/main.nf \
#   -profile singularity \
#   -params-file /tmp/params.yaml \
#   --input $RAW_DATA/dkd/dataset.h5ad \
#   --entry-name auto \
#   --config common/nextflow_helpers/labels_ci.config \
#   --labels task_batch_integration,process_datasets