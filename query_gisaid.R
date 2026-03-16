## ============================================================
## Query GISAID for avian influenza sequences from Madagascar
## and Africa using GISAIDR
##
## Usage:
##   module load R/4.4.0 libxml2/2.9.10
##   export GISAIDR_USERNAME="your_username"
##   export GISAIDR_PASSWORD="your_password"
##   LD_LIBRARY_PATH=/opt/gensoft/lib/libxml2/2.9.10/lib:/pasteur/appa/homes/cduitama/R/lib:$LD_LIBRARY_PATH \
##     Rscript query_gisaid.R
## ============================================================

lib <- "/pasteur/appa/homes/cduitama/R/library/R-4.4.0"
.libPaths(c(lib, .libPaths()))

library(GISAIDR)

# ── Credentials (set via environment variables, never hardcode) ──────────────
username <- Sys.getenv("GISAIDR_USERNAME")
password <- Sys.getenv("GISAIDR_PASSWORD")

if (username == "" || password == "") {
  stop("Please set GISAIDR_USERNAME and GISAIDR_PASSWORD environment variables.")
}

# ── Option 1: Try EpiFlu (avian influenza database) ─────────────────────────
cat("Attempting login to GISAID EpiFlu...\n")
credentials <- tryCatch(
  login(username = username, password = password, database = "EpiFlu"),
  error = function(e) {
    cat("EpiFlu login failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(credentials)) {
  cat("Login successful. Querying for avian influenza sequences...\n")

  # African countries in GISAID location format
  africa_locations <- c(
    "Africa / Madagascar",
    "Africa / Kenya",
    "Africa / Tanzania",
    "Africa / South Africa",
    "Africa / Egypt",
    "Africa / Nigeria",
    "Africa / Cameroon",
    "Africa / Senegal",
    "Africa / Ghana",
    "Africa / Ethiopia",
    "Africa / Uganda",
    "Africa / Burkina Faso",
    "Africa / Ivory Coast",
    "Africa / Niger",
    "Africa / Mali",
    "Africa / Guinea",
    "Africa / Mozambique",
    "Africa / Zimbabwe",
    "Africa / Zambia",
    "Africa / Democratic Republic of the Congo"
  )

  all_results <- data.frame()

  for (loc in africa_locations) {
    cat("  Querying:", loc, "...\n")
    df <- tryCatch(
      query(credentials = credentials, location = loc, load_all = TRUE),
      error = function(e) {
        cat("    Error:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(df) && nrow(df) > 0) {
      cat("    Found", nrow(df), "sequences\n")
      all_results <- rbind(all_results, df)
    } else {
      cat("    No sequences found\n")
    }
  }

  if (nrow(all_results) > 0) {
    cat("\nTotal sequences found:", nrow(all_results), "\n")

    # Select key metadata columns
    cols_of_interest <- intersect(
      c("accession_id", "virus_name", "location", "collection_date",
        "submission_date", "length", "host", "subtype", "segment",
        "originating_lab", "submitting_lab"),
      colnames(all_results)
    )
    metadata <- all_results[, cols_of_interest]

    # Add completeness flag (avian flu genome ~13,500 nt for full 8 segments)
    if ("length" %in% colnames(metadata)) {
      metadata$complete <- metadata$length >= 13000
    }

    # Save results
    outfile <- "avian_flu_Africa_metadata.tsv"
    write.table(metadata, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Results saved to:", outfile, "\n")

    # Summary
    cat("\n── Summary ──────────────────────────────\n")
    cat("Sequences by location:\n")
    print(sort(table(metadata$location), decreasing = TRUE))
    if ("complete" %in% colnames(metadata)) {
      cat("\nComplete genomes:", sum(metadata$complete, na.rm = TRUE), "\n")
      cat("Partial genomes:", sum(!metadata$complete, na.rm = TRUE), "\n")
    }
  } else {
    cat("No sequences found across all queried locations.\n")
    cat("NOTE: GISAIDR may not fully support EpiFlu — see fallback option below.\n")
  }

} else {
  # ── Option 2: Fallback — guide for manual GISAID EpiFlu download ───────────
  cat("\n── FALLBACK: Manual GISAID EpiFlu download ──────────────────────────────\n")
  cat("GISAIDR does not officially support EpiFlu (avian influenza).\n")
  cat("Please download metadata manually from GISAID:\n\n")
  cat("  1. Go to https://www.gisaid.org and log in\n")
  cat("  2. Click on 'EpiFlu' database\n")
  cat("  3. Search with these filters:\n")
  cat("     - Type: A\n")
  cat("     - Location: Africa\n")
  cat("     - Host: avian / bird\n")
  cat("  4. Download metadata as CSV/TSV\n")
  cat("  5. Place the file as 'avian_flu_Africa_metadata.tsv' in this directory\n")
  cat("  6. Run analyse_metadata.R to process it\n")
}
