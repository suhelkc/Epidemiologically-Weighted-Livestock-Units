# 02_wahis_attack_rates.R -------------------------------------------------
# Reads the WAHIS early-warning extract, keeps domestic FMD outbreak-species rows,
# and estimates pooled within-outbreak attack rates by species with an
# event-clustered bootstrap.
source("R/00_setup.R")
set.seed(SEED)

# Reading the 55 MB workbook with readxl needs several GB of RAM. If data/processed/wahis_subset.csv
# exists (the 14 columns used here, exported once), it is read instead; otherwise the xlsx is parsed.
keep <- c("epi_event_id","disease_eng","iso_code","country","region","Outbreak_id","Outbreak_start_date",
          "Longitude","Latitude","is_wild","Species","susceptible","cases","IN_FUR","Is_cluster","N_outbreaks_cluster")
csv  <- file.path(P$proc, "wahis_subset.csv")
if (file.exists(csv)) {
  raw <- fread(csv)
} else {
  f     <- file.path(P$raw, "WAHIS_data.xlsx")
  sheet <- grep("^infur", excel_sheets(f), value = TRUE)[1]
  hdr   <- names(read_excel(f, sheet = sheet, n_max = 1))
  raw   <- as.data.table(read_excel(f, sheet = sheet, col_types = ifelse(hdr %in% keep, "guess", "skip")))
  fwrite(raw[, ..keep], csv)
}
raw[, Outbreak_start_date := as.Date(Outbreak_start_date)]
raw[, is_wild := as.logical(is_wild)]

fmd <- raw[grepl("Foot and mouth", disease_eng, ignore.case = TRUE) & (is_wild %in% FALSE)]
cat("FMD domestic rows:", nrow(fmd), "| outbreaks:", uniqueN(fmd$Outbreak_id), "| countries:", uniqueN(fmd$country), "\n")
reg <- unique(fmd[, .(Outbreak_id, region)])[, .N, by = region]
print(reg[order(-N)])

grp <- c("Cattle" = "Cattle", "Buffaloes" = "Buffalo", "Swine" = "Pigs", "Sheep" = "Sheep",
         "Goats" = "Goats", "Sheep/goats (mixed herd)" = "Sheep/goats mixed")
ar <- copy(fmd)[, sp := grp[Species]][!is.na(sp) & susceptible > 0 & !is.na(cases) & cases <= susceptible]
ar <- dedupe_species_rows(ar)   # one row per outbreak x species (IN preferred over FUR)
cat("Events used:", uniqueN(ar$epi_event_id), "\n")

pooled <- function(d) d[, .(cases = sum(cases), susc = sum(susceptible)), by = sp][, ar := cases / susc][]
est <- pooled(ar)
ev  <- split(ar, ar$epi_event_id)
B   <- 1000
bs  <- rbindlist(lapply(seq_len(B), function(b) {
  pick <- sample(names(ev), length(ev), replace = TRUE)
  pooled(rbindlist(ev[pick]))[, b := b]
}))
ci  <- bs[, .(lo = quantile(ar, 0.025), hi = quantile(ar, 0.975)), by = sp]
out <- merge(est, ci, by = "sp")
out <- merge(out, ar[, .(n_records = .N), by = sp], by = "sp")
setorder(out, -n_records)
fwrite(out, file.path(P$out, "T4_attack_rates.csv"))
print(out[, .(sp, n_records, attack_rate = round(100 * ar, 1), lo = round(100 * lo, 1), hi = round(100 * hi, 1))])
