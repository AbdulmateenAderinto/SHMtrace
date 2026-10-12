# Export a SHazaM 5-mer targeting model in the CSV format shmtrace reads.
#
#   Rscript scripts/export_shazam_model.R HH_S5F  HH_S5F.csv
#   Rscript scripts/export_shazam_model.R MK_RS5NF MK_RS5NF.csv
#
# shmtrace does not ship these models: SHazaM is AGPL-3 and shmtrace is MIT, so
# the model is written from your own SHazaM installation. Cite Yaari et al.,
# Front Immunol 2013 (S5F) and the SHazaM package when you use the result.
args <- commandArgs(trailingOnly = TRUE)
name <- if (length(args) >= 1) args[1] else "HH_S5F"
out  <- if (length(args) >= 2) args[2] else paste0(name, ".csv")
suppressMessages(library(shazam))
m <- get(name)
mut <- m@mutability
sub <- m@substitution
keep <- grepl("^[ACGT]{5}$", names(mut))
mut <- mut[keep]
df <- data.frame(fivemer = names(mut), mutability = as.numeric(mut),
                 to_A = sub["A", names(mut)], to_C = sub["C", names(mut)],
                 to_G = sub["G", names(mut)], to_T = sub["T", names(mut)])
write.csv(df, out, row.names = FALSE)
cat(sprintf("wrote %s: %d 5-mers from %s (%s), shazam %s\n", out, nrow(df), m@name,
            m@species, as.character(packageVersion("shazam"))))
