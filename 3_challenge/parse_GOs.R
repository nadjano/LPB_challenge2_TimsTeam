# Check if both input files were provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript parse_GOs.R <goa.gaf> <ids_input_file>", call. = FALSE)
}
# Check for requiered packages and install them if not already installed
required_packages <- c("tidyverse")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
# Import arguments
file_gaf <- args[1]
file_ids <- args[2]

# Import input files
uniprot_ids <- readLines(file_ids)
goa_names = c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier","GO_ID","DB_Reference", "Evidence_Code",
              "With_From", "Aspect" , "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon",
              "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID")

goa <- read_tsv(file_gaf, comment = "!", col_names = goa_names)

# Subset to uniprot IDs provided, remove IEA evidence elements and deduplicate
goa_uniprot <- goa %>%
  filter(DB_Object_ID %in% uniprot_ids) %>%
  filter(Evidence_Code != "IEA") %>%
  distinct(DB_Object_ID, GO_ID, Aspect)

# Rank GO terms by their occurences in our subset
go_rank <- goa_uniprot %>% count(GO_ID, name = "Rank")

# Join dfs
goa_uniprot_ranked <- goa_uniprot %>%
  left_join(go_rank)

# Keep only max 10 top GO terms per Uniprot protein
results <- map(uniprot_ids, \(id)
    goa_uniprot_ranked %>%
      group_by(Aspect, DB_Object_ID) %>%
      filter(DB_Object_ID == id) %>%
      slice_max(order_by = Rank, n = 10, with_ties = FALSE) %>%
      summarise(GO_terms = paste(GO_ID, collapse = ", "), .groups = "drop")
    ) %>%
  bind_rows()

# Reformat df to expected output format
formatted_results <- results %>% pivot_wider(names_from = Aspect, values_from = GO_terms) %>%
  rowwise() %>%
  mutate(
    output = paste0(
      "UniProtID: ", DB_Object_ID, "\n",
      "Biological Process: ", P, "\n",
      "Molecular Function: ", F, "\n",
      "Cellular Component: ", C, "\n"
    )
  ) %>%
  pull(output)

# Write formatted output to file
writeLines(formatted_results, "output.txt")

print("Script finished. Output was saved to output.txt.")