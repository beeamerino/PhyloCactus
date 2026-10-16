# Generate the SLURM Job Scripts of a Chunked Run of Step 7

Writes the scripts that run
[`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md)
as a SLURM job array, one chunk per task on one core, and
[`merge_barcoding_chunks()`](https://beeamerino.github.io/PhyloCactus/reference/merge_barcoding_chunks.md)
once every task has finished without error. It follows the conventions
of
[`generate_ml_search_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_ml_search_script.md)
and
[`generate_bootstrap_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_bootstrap_script.md):
`cluster_*` arguments, `load_module`, one log per array task, and the
versions the job ran with printed at its start.

## Usage

``` r
generate_barcoding_job_scripts(
  library_dir,
  folds_dir,
  classifier_dir,
  method = c("idtaxa", "nn"),
  n_chunks,
  job_dir = file.path(classifier_dir, "job"),
  seed = 1L,
  threshold = 60,
  schemes = c("species", "genus"),
  loci = NULL,
  max_folds = NULL,
  cluster_job_name = "cactus_idtaxa",
  cluster_partition = "main",
  cluster_mem = "4G",
  cluster_time = "04:00:00",
  cluster_queue = NULL,
  cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
  load_module = NULL,
  r_lib = NULL,
  rscript_exec = "Rscript"
)
```

## Arguments

- library_dir, folds_dir:

  Character. As in
  [`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md),
  on the cluster.

- classifier_dir:

  Character. `output_dir` of the run, on the cluster.

- method, seed, threshold, schemes, loci, max_folds:

  As in
  [`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md).

- n_chunks:

  Integer. Number of chunks, one array task each.

- job_dir:

  Character. Where the scripts are written. Defaults to `job/` inside
  `classifier_dir`.

- cluster_job_name:

  Character. SLURM job name; the merge job adds `_merge`. Defaults to
  `"cactus_idtaxa"`.

- cluster_partition:

  Character. SLURM partition. Defaults to `"main"`.

- cluster_mem:

  Character. Memory per task (`--mem`). Defaults to `"4G"`.

- cluster_time:

  Character. Time limit per task (`-t`). Defaults to `"04:00:00"`.

- cluster_queue:

  Character. Optional SLURM QoS (`-q`). Defaults to `NULL`.

- cluster_mail_user:

  Character. Address for SLURM mail. Defaults to
  `Sys.getenv("MY_EMAIL", "")`; an empty value writes no mail line. The
  array mails on failure only; the merge mails at its end.

- load_module:

  Character vector of environment modules to load, or `NULL`.

- r_lib:

  Character. R library to put first (`R_LIBS_USER`), or `NULL`.

- rscript_exec:

  Character. The `Rscript` command. Defaults to `"Rscript"`.

## Value

Invisibly, the paths of the five files.

## Details

Five files are written in `job_dir`: `run_chunk.R` and `run_merge.R`,
which hold the arguments of the run; `job_array.sh`, the array
(`--array=1-n_chunks`, one core per task); `job_merge.sh`, the merge;
and `submit.sh`, which submits the array and then the merge with
`--dependency=afterok:<array>`. Submit with `bash submit.sh` from any
directory.

The paths are written as given and must be those of the cluster: pass
absolute paths, or paths relative to the directory the job is submitted
from. Nothing about a particular cluster is written unless it is passed:
with `load_module = NULL` no module is loaded, with
`cluster_queue = NULL` no QoS is requested, and with `r_lib = NULL` the
default R library is used.

## See also

[`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md),
[`merge_barcoding_chunks()`](https://beeamerino.github.io/PhyloCactus/reference/merge_barcoding_chunks.md).

## Examples

``` r
if (FALSE) { # \dontrun{
generate_barcoding_job_scripts(
  library_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/4_library",
  folds_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/5_folds",
  classifier_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/7_classifier",
  method = "idtaxa", n_chunks = 60,
  load_module = "R/4.4.0", r_lib = "~/R/phylocactus/4.4"
)
} # }
```
