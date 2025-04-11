
rm(list = ls())

source("project_support.R")

spp <- read.csv("data/corvid_species.csv")

spp$mimicry_category <- case_when(
  spp$mimicry_primary == 1 & spp$mimicry_secondary == 1 ~ "1° & 2° evidence",
  spp$mimicry_primary == 1 & spp$mimicry_secondary == 0 ~ "1° evidence only",
  (is.na(spp$mimicry_primary) | spp$mimicry_primary == 0) & spp$mimicry_secondary == 1 ~ "2° evidence only",
  (is.na(spp$mimicry_primary) | spp$mimicry_primary == 0) & spp$mimicry_secondary == 0 ~ "no mimicry seen"
)

stopifnot(!any(is.na(spp$mimicry_category)))

spp$mimicry_category <- factor(spp$mimicry_category, levels = c("no mimicry seen", "2° evidence only", "1° evidence only", "1° & 2° evidence"))

spp |>
  group_by(mimicry_category) |>
  summarize(
    n_species = n(),
    `2° sources per species` = round(sum(corvid_database_entries) / n_species),
    `1° sources per species` = round(sum(primary_count) / n_species),
  ) -> cat_counts

write_pandoc_table(cat_counts, "figures/cat_counts.md", style = "rmarkdown")

# tinyplot orders things in the same order as the factor order

dot_cols = c(
  "no mimicry seen" = "gray",
  "2° evidence only" = "orange",
  "1° evidence only" = "blue",
  "1° & 2° evidence" = "red"
)

png("figures/double_evidence_detection.png", res = 300, units = "in", height = 7, width = 10)

tinyplot(spp$secondary_count, spp$primary_count, by = spp$mimicry_category, log = "xy", xlab = "number of secondary sources", ylab = "number of primary sources", col = col_alpha(dot_cols, 0.5), pch = 16, xlim = c(1, 800), ylim = c(1, 3000), legend = list(title = "mimicry category"), main = "Corvid mimicry status (n = 128 spp.)")

tar <- which(spp$scientific_name_short == "C. corax")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$scientific_name_short == "G. glandarius")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$scientific_name_short == "P. pica")
text(spp$secondary_count[tar], spp$primary_count[tar],
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 4)

tar <- which(spp$scientific_name_short == "C. stelleri")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$scientific_name_short == "C. yncas")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$scientific_name_short == "C. moneduloides")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$scientific_name_short == "A. coerulescens")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$scientific_name_short == "P. infaustus")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$scientific_name_short == "C. brachyrhynchos")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 4)

tar <- which(spp$scientific_name_short == "C. cornix")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 2)

tar <- which(spp$scientific_name_short == "C. monedula")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$scientific_name_short == "C. cristata")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$scientific_name_short == "C. corone")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 4)

tar <- which(spp$scientific_name_short == "P. hudsonia")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$scientific_name_short == "N. caryocatactes")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$scientific_name_short == "C. frugilegus")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$scientific_name_short == "C. chrysops")
text(spp$secondary_count[tar], spp$primary_count[tar], 
  labels = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

dev.off()
