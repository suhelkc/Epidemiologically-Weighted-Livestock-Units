# 03_build_blocks.R -------------------------------------------------------
# Harmonise the GLW4 rasters, assign cells to countries, aggregate to 0.5-degree blocks nested in
# countries, and attach WAHIS FMD outbreaks to blocks.
source("R/00_setup.R")
R_BLOCK <- 6L                     # 6 x 6 cells of 5 arcmin = 0.5 degree

## 1. Rasters -> counts per cell ------------------------------------------
tmpl <- rast(file.path(P$raw, GLW_FILES[["cattle"]]))
area <- cellSize(tmpl, unit = "km")                       # km2 per cell (geodesic)
lay  <- lapply(GLW_FILES, function(f) {
  r <- rast(file.path(P$raw, f)); stopifnot(compareGeom(r, tmpl)); r * area
})
names(lay) <- names(GLW_FILES)
glob <- data.table(species = c(SPECIES, "LU"),
                   total_millions = sapply(lay[c(SPECIES, "LU")], function(r) global(r, "sum", na.rm = TRUE)[1, 1]) / 1e6)
fwrite(glob, file.path(P$proc, "global_totals.csv")); print(glob)

## 2. Country raster --------------------------------------------------------
adm  <- vect(file.path(P$raw, "world-administrative-boundaries.shp"))
isos <- unique(adm$iso3)
adm$iso_id <- match(adm$iso3, isos)
iso_r <- rasterize(adm, tmpl, field = "iso_id")

## 3. Cell table -> block-by-country units ---------------------------------
dt <- data.table(cell = seq_len(ncell(tmpl)), iso_id = values(iso_r)[, 1])
for (nm in names(lay)) set(dt, j = nm, value = values(lay[[nm]])[, 1])
dt <- dt[!is.na(iso_id)]
dt[, c("row", "col") := .(rowFromCell(tmpl, cell), colFromCell(tmpl, cell))]
dt[, c("br", "bc") := .((row - 1L) %/% R_BLOCK, (col - 1L) %/% R_BLOCK)]
cols <- c(SPECIES, "LU")
units <- dt[, c(lapply(.SD, sum, na.rm = TRUE), list(n_cells = .N)), by = .(br, bc, iso_id), .SDcols = cols]
units[, iso3 := isos[iso_id]]
units[, unit_id := .I]
rm(dt); gc()
fwrite(units, file.path(P$proc, "units.csv"))
cat("Units (block x country):", nrow(units), "| countries:", uniqueN(units$iso3), "\n")

## 4. Outbreaks -> units ----------------------------------------------------
wah <- fread(file.path(P$proc, "wahis_subset.csv"))
wah[, is_wild := as.logical(is_wild)]
ob <- wah[grepl("Foot and mouth", disease_eng, ignore.case = TRUE) & (is_wild %in% FALSE)]
ob <- unique(ob[, .(Outbreak_id, epi_event_id, iso_code, region, Longitude, Latitude,
                    year = as.integer(substr(Outbreak_start_date, 1, 4)))], by = "Outbreak_id")
n_all <- nrow(ob)
ob <- ob[!is.na(Longitude) & !is.na(Latitude)]
cell <- cellFromXY(tmpl, as.matrix(ob[, .(Longitude, Latitude)]))
ob <- ob[!is.na(cell)][, cell := cell[!is.na(cell)]]
ob[, c("row", "col") := .(rowFromCell(tmpl, cell), colFromCell(tmpl, cell))]
ob[, c("br", "bc") := .((row - 1L) %/% R_BLOCK, (col - 1L) %/% R_BLOCK)]
ob[, iso_poly := isos[values(iso_r)[cell, 1]]]
ob[, iso3 := fifelse(!is.na(iso_poly), iso_poly, iso_code)]       # fall back to WAHIS ISO code
ob <- merge(ob, units[, .(br, bc, iso3, unit_id)], by = c("br", "bc", "iso3"), all.x = TRUE)
cat(sprintf("Outbreaks: %d domestic FMD; %d with coordinates; %d matched to a unit (%.1f%%)\n",
            n_all, nrow(ob), sum(!is.na(ob$unit_id)), 100 * mean(!is.na(ob$unit_id))))
cat(sprintf("Polygon ISO equals WAHIS ISO: %.1f%%\n", 100 * mean(ob$iso_poly == ob$iso_code, na.rm = TRUE)))
fwrite(ob[, .(Outbreak_id, epi_event_id, iso_code, region, year, iso3, unit_id)], file.path(P$proc, "fmd_outbreaks_units.csv"))
