# f2dace-qe-source

This repo contains source files for `pw.x` of Quantum Espresso (`develop` branch as of April 2026).

All preprocessed files to be ingested are in `src_f90_to_parse/`.

The script `run.sh` contains the execution script used in the previous fparser-based frontend. This script attempts pruning for the selected top subroutine `vexx_b_k_gpu`.

The original fparser AST constructed is `./ckpt/ast_v0.f90`.

The resulting checkpoint (before errors are thrown) is in `./ckpt/ast_v1_vexx_bp_k_gpu.f90`.

## Examples
Code examples of specific cases are in `./examples/`.

Example of EXTERNAL: The function `local_kpoint_index` (found in `divide_et_impera.f90`) is referenced in multiple subroutines in `exx_bp.f90`.

