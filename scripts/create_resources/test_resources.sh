#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# remove this when you have implemented the script
echo "TODO: replace the commands in this script with the sequence of components that you need to run to generate test_resources."
echo "  Inside this script, you will need to place commands to generate example files for each of the 'src/api/file_*.yaml' files."
exit 1

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/task_batch_integration

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  --publish_dir "$DATASET_DIR" \
  --id "pancreas" \
  --input "$RAW_DATA/cxg_mouse_pancreas_atlas/dataset.h5ad" \
  --output_train '$id/train.h5ad' \
  --output_test '$id/test.h5ad' \
  --output_solution '$id/solution.h5ad' \
  --output_state '$id/state.yaml'

# run one method
viash run src/methods/knn/config.vsh.yaml -- \
    --input_train $DATASET_DIR/cxg_mouse_pancreas_atlas/train.h5ad \
    --input_test $DATASET_DIR/cxg_mouse_pancreas_atlas/test.h5ad \
    --output $DATASET_DIR/cxg_mouse_pancreas_atlas/prediction.h5ad

# run one metric
viash run src/metrics/accuracy/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/cxg_mouse_pancreas_atlas/prediction.h5ad \
    --input_solution $DATASET_DIR/cxg_mouse_pancreas_atlas/solution.h5ad \
    --output $DATASET_DIR/cxg_mouse_pancreas_atlas/score.h5ad

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  "$DATASET_DIR" s3://openproblems-data/resources_test/task_batch_integration \
  --delete --dryrun
