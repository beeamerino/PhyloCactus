# -------------------------------------------------------------
# PhyloCactus: Tutorial 6 - Identifying a Sample (11_barcoding/10_identify/)
# -------------------------------------------------------------
# Worked examples of the identification of PhyloCactus 0.5.0, one for each kind of input the package
# accepts, all of them public records that are not in the library: a Sanger record, a specimen of
# unknown species, two specimens with several loci, two GenBank plastomes, a sample sheet, a genome
# assembly and paired Illumina reads. Each example is downloaded in a block of its own and then
# identified with the path of its file, as with a sample of your own.
#
# It needs the library of Tutorial 5 in this folder: 11_barcoding/4_library/ (step 4),
# 11_barcoding/9_threshold/ (step 9, the species thresholds of rule R-c) and 11_barcoding/11_metrics/
# (step 10). The downloads need a network connection and the package rentrez. The assembly needs
# minimap2 and samtools, the reads need fastp and GetOrganelle (with its seed databases embplant_pt and
# embplant_nr); in a conda environment, activate it before starting R from the same terminal.
#
# Each identification writes its table, its run record and its HTML report in
# 11_barcoding/10_identify/ (the reads in 11_barcoding/10_identify/<run_name>/).
# -------------------------------------------------------------
library(PhyloCactus)

tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)
dir.create("examples", showWarnings = FALSE)

# A GenBank record, or several, saved as a flat file: the package reads the accession of each record
# as the name of its query.
fetch_genbank <- function(accessions, file) {
  if (!file.exists(file)) {
    writeLines(rentrez::entrez_fetch("nuccore", id = accessions, rettype = "gb", retmode = "text"), file)
  }
  file
}

# -------------------------------------------------------------
# Example 1. One Sanger record: LT745702.1, matK of Pseudorhipsalis amazonica subsp. panamensis
# (GenBank, 1 record, about 1 kb).
# -------------------------------------------------------------
fetch_genbank("LT745702.1", "examples/LT745702.1.gb")

identify_barcoding_query("examples/LT745702.1.gb", run_name = "LT745702.1")

# -------------------------------------------------------------
# Example 2. A specimen of unknown species: PZ714536.1, matK of Opuntia sp. (GenBank, 1 record).
# -------------------------------------------------------------
fetch_genbank("PZ714536.1", "examples/PZ714536.1.gb")

identify_barcoding_query("examples/PZ714536.1.gb", run_name = "PZ714536.1")

# -------------------------------------------------------------
# Example 3. One specimen, four loci: Coryphantha elephantidens subsp. bumamma, voucher Vazquez 2628
# (GenBank OK340410.1, OK340349.1, OK340287.1 and OK340224.1). The four records are one sample, read
# together in one report.
# -------------------------------------------------------------
fetch_genbank(c("OK340410.1", "OK340349.1", "OK340287.1", "OK340224.1"), "examples/OK340410.1_specimen.gb")

identify_barcoding_query("examples/OK340410.1_specimen.gb", sample = "Coryphantha_elephantidens_bumamma",
                         run_name = "OK340410.1")

# -------------------------------------------------------------
# Example 4. One specimen, four loci, and a wrong species: Xiquexique gounellei
# (GenBank OP189074.1, OP188885.1, OP188764.1 and OP187233.1).
# -------------------------------------------------------------
fetch_genbank(c("OP189074.1", "OP188885.1", "OP188764.1", "OP187233.1"), "examples/OP189074.1_specimen.gb")

identify_barcoding_query("examples/OP189074.1_specimen.gb", sample = "Xiquexique_gounellei",
                         run_name = "OP189074.1")

# -------------------------------------------------------------
# Example 5. A GenBank plastome named to species: MN517611.1, Mammillaria zephyranthoides
# (GenBank, 107 kb).
# -------------------------------------------------------------
fetch_genbank("MN517611.1", "examples/MN517611.1.gb")

identify_barcoding_query("examples/MN517611.1.gb", run_name = "MN517611.1")

# -------------------------------------------------------------
# Example 6. A GenBank plastome named as another species: MT387452.1, Rhipsalis teres (GenBank, 120 kb).
# -------------------------------------------------------------
fetch_genbank("MT387452.1", "examples/MT387452.1.gb")

identify_barcoding_query("examples/MT387452.1.gb", run_name = "MT387452.1")

# -------------------------------------------------------------
# Example 7. Several samples at once, from a sample sheet (one row per sample; the input types are
# those of the matrix of accepted inputs, inst/extdata/barcoding_accepted_inputs.csv).
# -------------------------------------------------------------
sheet <- data.frame(
  sample_id = c("LT745702.1", "PZ714536.1", "MN517611.1"),
  input_type = c("genbank", "genbank", "genbank"),
  file_1 = c("examples/LT745702.1.gb", "examples/PZ714536.1.gb", "examples/MN517611.1.gb"),
  declared_species = c("Pseudorhipsalis amazonica", "", "Mammillaria zephyranthoides")
)
utils::write.csv(sheet, "examples/samples.csv", row.names = FALSE)

identify_barcoding_samples("examples/samples.csv")

# -------------------------------------------------------------
# Example 8. A genome assembly: GCA_024363205.1, Cereus fernambucensis (NCBI Assembly, about 0.9 Gb,
# a compressed FASTA of a few hundred MB). The FTP folder of an assembly is read from its NCBI summary:
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/363/205/<accession>_<assembly name>/
# -------------------------------------------------------------
asm_summary <- rentrez::entrez_summary("assembly",
                                       rentrez::entrez_search("assembly", "GCA_024363205.1[Assembly Accession]")$ids[1])
asm_ftp <- sub("^ftp://", "https://", asm_summary$ftppath_genbank)
asm_file <- file.path("examples", paste0(basename(asm_ftp), "_genomic.fna.gz"))
if (!file.exists(asm_file)) {
  utils::download.file(paste0(asm_ftp, "/", basename(asm_ftp), "_genomic.fna.gz"), asm_file, mode = "wb")
}

identify_barcoding_assembly(asm_file, run_name = "GCA_024363205.1")

# -------------------------------------------------------------
# Example 9. Paired Illumina reads: SRR26398994, a genome skim of Opuntia comonduensis (ENA, two
# FASTQ files, about 0.56 GB). The links and the md5 of the two files are read from the file report
# of ENA; a file that does not match its md5 is downloaded again.
# -------------------------------------------------------------
ena <- utils::read.delim(paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR26398994",
                                "&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"),
                         colClasses = "character")
fq_urls <- paste0("https://", strsplit(ena$fastq_ftp[1], ";")[[1]])
fq_md5 <- strsplit(ena$fastq_md5[1], ";")[[1]]
fq <- file.path("examples", basename(fq_urls))
for (k in seq_along(fq)) {
  if (!file.exists(fq[k]) || unname(tools::md5sum(fq[k])) != fq_md5[k]) {
    utils::download.file(fq_urls[k], fq[k], mode = "wb")
  }
  stopifnot(unname(tools::md5sum(fq[k])) == fq_md5[k])
}

identify_barcoding_reads(fq[1], fq[2], run_name = "SRR26398994", threads = 4L)

# -------------------------------------------------------------
# Figures: the answer of one identification, with the sampling of the named genus, and the answers of
# all the examples side by side. Both are drawn from the tables of the reports.
# -------------------------------------------------------------
id_dir <- file.path("11_barcoding", "10_identify")
plot_barcoding_identification(file.path(id_dir, "REPORT_LT745702.1_answer.csv"),
                              sampling = file.path(id_dir, "REPORT_LT745702.1_sampling.csv"))

answers <- list(
  "LT745702.1 (Sanger)" = file.path(id_dir, "REPORT_LT745702.1_answer.csv"),
  "PZ714536.1 (Sanger, Opuntia sp.)" = file.path(id_dir, "REPORT_PZ714536.1_answer.csv"),
  "Coryphantha elephantidens (4 loci)" = file.path(id_dir, "REPORT_OK340410.1_answer.csv"),
  "Xiquexique gounellei (4 loci)" = file.path(id_dir, "REPORT_OP189074.1_answer.csv"),
  "MN517611.1 (plastome)" = file.path(id_dir, "REPORT_MN517611.1_answer.csv"),
  "MT387452.1 (plastome)" = file.path(id_dir, "REPORT_MT387452.1_answer.csv"),
  "GCA_024363205.1 (assembly)" = file.path(id_dir, "REPORT_GCA_024363205.1_answer.csv"),
  "SRR26398994 (reads)" = file.path(id_dir, "SRR26398994", "REPORT_SRR26398994_answer.csv")
)
plot_barcoding_identification_summary(answers)
