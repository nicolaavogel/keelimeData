# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
setwd("/Users/bfj994/Documents/keelime/")
# Define a custom color palette
tool_colors <- c("keelime_reckless" = "#1b9e77", "keelime_normal" = "#33a02c", "keelime_strict" =  "#00441b","endoCaller" = "#e41a1c", "SPAdes" = "#377eb8")

### First Plot

# Read data from the Excel file
data <- read_excel("/Users/bfj994/Documents/keelime/keelimeResults.xlsx", sheet = "edDEndoLow")

# Create a new column combining Simulated Age and Simulated Substitutions for labels
data <- data %>%
  mutate(Simulated_Label = paste(Simulated.Age, "(", Simulated.Substitutions, ")", sep = ""))

# Sort data by Simulated Substitutions
data <- data %>%
  arrange(Simulated.Substitutions)

# Convert Simulated_Label to a factor with ordered levels based on Simulated Substitutions
data$Simulated_Label <- factor(data$Simulated_Label, levels = unique(data$Simulated_Label))

# Custom labeller function
custom_labeller <- as_labeller(c(
  "500" = "~1.3X coverage (500 simulated reads)",
  "5000" = "~13.2X coverage (5000 simulated reads)",
  "25000" = "~66,2X coverage (25000 simulated reads)"
))

# Plot
plot1 <- ggplot(data, aes(x = Simulated_Label, y = No.Correct.Bases, color = Comparison.Sequence.Tool)) +
  geom_point(size = 3) +
  geom_line(aes(group = Comparison.Sequence.Tool), size = 1) +
  labs(x = "Simulated Age (Simulated Substitutions)",
       y = "Number of correct bases ",
       color = "Tool") +
  theme_bw(base_size = 15) +
  facet_wrap(~ data$`Number of reads`, labeller = custom_labeller)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = tool_colors)

plot1
ggsave("diver.png", plot1, dpi = 600)
### Second Plot

# Read data
data2 <- read_excel("/Users/bfj994/Documents/keelime/keelimeResults.xlsx", sheet = "edDNovoLow")

# Clean coverage: replace commas with dots for consistency
data2$Simulated.Coverage <- str_replace_all(data2$Simulated.Coverage, ",", ".")

# Reorder Simulated.Coverage by Simulated.Number.of.reads
data2 <- data2 %>%
  mutate(Simulated.Number.of.reads = as.numeric(Simulated.Number.of.reads)) %>%
  arrange(desc(Simulated.Number.of.reads)) %>%
  mutate(Simulated.Coverage = factor(Simulated.Coverage, levels = unique(Simulated.Coverage)))

# Convert damage to factor
data2$Simulated.ancient.damage <- factor(data2$Simulated.ancient.damage, levels = c("No Damage", "High Damage"))

# Plot
plot2 <- ggplot(data2, aes(x = Simulated.Coverage, y = No.Correct.Bases, color = Comparison.Sequence.Tool)) +
  geom_point(size = 3) +
  geom_line(aes(group = Comparison.Sequence.Tool), size = 1) +
  labs(
    x = "Simulated Coverage",
    y = "Number of Correct Bases",
    color = "Tool"
  ) +
  facet_wrap(~ Simulated.ancient.damage) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = tool_colors)
plot2
# Save the plot
ggsave("downsample_kee.png", plot2, dpi = 600)
# Print the plots
print(plot1)
print(plot2)


# Read data from the Excel file
data3 <- read_excel("/Users/bfj994/Documents/keelime/keelimeResults.xlsx", sheet = "time_kee")
# Convert data to long format for faceting
data_long <- data3 %>%
  pivot_longer(cols = c("User time (min)", "Memory (GB)"),
               names_to = "Metric",
               values_to = "Value")

# Convert 'No of reads' to a factor to preserve the order and reverse the levels
data_long$`No of reads` <- factor(data_long$`No of reads`, levels = rev(unique(data_long$`No of reads`)))

# Plot
ggplot(data_long, aes(x = `No of reads`, y = Value, color = Metric)) +
  geom_point(size = 3) +
  geom_line(aes(group = Metric), size = 1) +
  labs(title = "No of Reads vs User Time and Memory Usage",
       x = "Number of Reads",
       y = "Value",
       color = "Metric") +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = rev(levels(data_long$`No of reads`)))

library(ggplot2)
library(dplyr)
library(readxl)
library(ggtree)

# Read the data
data4 <- read_excel("/Users/bfj994/Documents/keelime/keelimeResults.xlsx", sheet = "MissingSeq")

# Convert Coverage to factor based on NoFrag order
data4 <- data4 %>%
  mutate(Coverage = factor(Coverage, levels = unique(Coverage[order(NoFrag)]))) %>%
  mutate(Organelle = factor(Organelle, levels = c("Mitochondrion", "Chloroplast")))

# Define colors for organelles and Miss labels
organelle_colors <- c("Mitochondrion" = "#FFB3B3", 
                      "Chloroplast" = "#B3FFB3") 

miss_color <- "lightgrey"

# Create the plot
plot4 <- ggplot(data4, aes(x = Coverage, y = Correct.Bases, color = Tool)) +
  geom_point(size = 3) +
  geom_line(aes(group = Tool), size = 1) +
  labs(
    x = "Simulated Coverage",
    y = "Number of Correct Bases",
    color = "Tool"
  ) +
  facet_wrap(~Organelle + Miss, scales = "free_y", nrow = 2, ncol = 2) +  # Facet by Organelle and Miss
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 13, color = "black"),
    strip.background = element_rect(fill = miss_color),  # Default grey for all Miss
    panel.spacing = unit(1, "lines")  # Space between facets
  ) +
  scale_color_manual(values = tool_colors)

plot4
ggsave("Missing_seq.all.png", plot4, width = 15, height = 10, dpi = 600)

################################################################################
## EMP plots

# Load necessary libraries
library(ape)
library(ggtree)
library(rentrez)
library(phangorn)
library(stringr)

# Function to fetch scientific name from GenBank or return input for special labels
get_scientific_name <- function(input_string) {
  # If the label is "Consensus" or starts with "Kap", return it directly
  if (grepl("^Consensus$|^Kap", input_string)) {
    return(input_string)
  }
  
  # Extract the accession number before the first '.'
  accession_number <- substr(input_string, 1, regexpr("\\.", input_string) - 1)
  
  # Fetch the GenBank record
  record <- tryCatch({
    entrez_fetch(db = "nucleotide", id = accession_number, rettype = "gb", retmode = "text")
  }, error = function(e) {
    return(NA)
  })
  
  # Extract scientific name
  if (!is.na(record)) {
    organism_pattern <- "ORGANISM\\s+(.*?)\\n"
    matches <- regmatches(record, regexpr(organism_pattern, record, perl = TRUE))
    if (length(matches) > 0) {
      scientific_name <- trimws(sub("ORGANISM\\s+", "", matches[[1]]))
    } else {
      scientific_name <- "Unknown"
    }
  } else {
    scientific_name <- "FetchError"
  }
  
  return(scientific_name)
}


# Set path to your RAxML tree file
empBs <- "/Users/bfj994/Documents/keelime/empTrees/RAxML_bestTree.empBstrict"
empBn <- "/Users/bfj994/Documents/keelime/empTrees/RAxML_bestTree.empBnorm"
empBr <- "/Users/bfj994/Documents/keelime/empTrees/RAxML_bestTree.empBreck"

empAs <- "/Users/bfj994/Documents/keelime/empTrees/RAxML_bestTree.empAstrict"
empAn <- "/Users/bfj994/Documents/keelime/empTrees/RAxML_bestTree.empAnorm"
empAr <- "/Users/bfj994/Documents/keelime/empTrees/RAxML_bestTree.empAreck"

empKs <- "/Users/bfj994/Documents/keelime/RAxML_bestTree.keeKapstrict_nonu_cleaned_rotated_raxml"
empKn <- "/Users/bfj994/Documents/keelime/RAxML_bestTree.kapK_nonu_cleaned_normal_rotated"
empKr <- "/Users/bfj994/Documents/keelime/RAxML_bestTree.kapK_nonu_cleaned_reckless_rotated"


midpoint.root <- function(tree) {
  d <- cophenetic(tree)
  tips <- which(d == max(d), arr.ind = TRUE)[1, ]
  root(tree, node = getMRCA(tree, tips))
}

# Extract accession number
extract_accession <- function(label) {
  sub("([A-Z]{1,2}_?[0-9]+\\.[0-9]+).*", "\\1", label)
}
# Read the tree
tBs <- read.tree(empBs)
label_map <- setNames(sapply(tBs$tip.label, get_scientific_name), tBs$tip.label)
tBs$tip.label <- label_map[tBs$tip.label]

ggtree(tBs) +                             
  geom_tiplab(aes(label = label, color = label == "Consensus"), size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  geom_nodelab(geom = 'label', size = 3) + 
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with strict mode in keelime") +
  theme(
    plot.margin = margin(0.1, 8, 0.1, 0.1, "cm"),
    legend.position = "none")+
  coord_cartesian(clip = "off")
ggsave("empBs.png", width = 13, height = 8, dpi = 600)

# Read the tree
tBs <- read.tree(empBn)
label_map <- setNames(sapply(tBs$tip.label, get_scientific_name), tBs$tip.label)
tBs$tip.label <- label_map[tBs$tip.label]

ggtree(tBs) +                             
  geom_tiplab(aes(label = label, color = label == "Consensus"), size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  geom_nodelab(geom = 'label', size = 3) + 
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with normal mode in keelime") +
  theme(
    plot.margin = margin(0.1, 8, 0.1, 0.1, "cm"),
    legend.position = "none")+
  coord_cartesian(clip = "off")
ggsave("empBn.png", width = 13, height = 8, dpi = 600)

# Read the tree
tBs <- read.tree(empBr)
label_map <- setNames(sapply(tBs$tip.label, get_scientific_name), tBs$tip.label)
tBs$tip.label <- label_map[tBs$tip.label]

ggtree(tBs) +                             
  geom_tiplab(aes(label = label, color = label == "Consensus"), size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  geom_nodelab(geom = 'label', size = 3) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with reckless mode in keelime") +
  theme(
    plot.margin = margin(0.1, 8, 0.1, 0.1, "cm"),
    legend.position = "none")+
  coord_cartesian(clip = "off")
ggsave("empBr.png", width = 13, height = 8, dpi = 600)

# Read the tree
tBs <- read.tree(empAs)
tBs_m <- midpoint.root(tBs)
label_map <- setNames(sapply(tBs$tip.label, get_scientific_name), tBs$tip.label)
tBs_m$tip.label <- label_map[tBs$tip.label]

ggtree(tBs_m) +                             
  geom_tiplab(aes(label = label, color = label == "Consensus"), size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  geom_nodelab(geom = 'label', size = 3) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with strict mode in keelime") +
  theme(
    plot.margin = margin(0.1, 8, 0.1, 0.1, "cm"),
    legend.position = "none")+
  coord_cartesian(clip = "off")
ggsave("empAs.png", width = 15, height = 12, dpi = 600)

# Read the tree
tBs <- read.tree(empAn)
tBs_m <- midpoint.root(tBs)
label_map <- setNames(sapply(tBs$tip.label, get_scientific_name), tBs$tip.label)
tBs_m$tip.label <- label_map[tBs$tip.label]
ggtree(tBs_m) +                             
  geom_tiplab(aes(label = label, color = label == "Consensus"), size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  geom_nodelab(geom = 'label', size = 3) + 
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with normal mode in keelime") +
  theme(
    plot.margin = margin(0.1, 15, 0.1, 0.1, "cm"),
    legend.position = "none")+
  coord_cartesian(clip = "off")
ggsave("empAn.png", width = 15, height = 12, dpi = 600)

# Read the tree
tBs <- read.tree(empAr)
tBs_m <- midpoint.root(tBs)
label_map <- setNames(sapply(tBs$tip.label, get_scientific_name), tBs$tip.label)
tBs_m$tip.label <- label_map[tBs$tip.label]

ggtree(tBs_m) +                             
  geom_tiplab(aes(label = label, color = label == "Consensus"), size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  geom_nodelab(geom = 'label', size = 3) + 
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with reckless mode in keelime") +
  theme(
    plot.margin = margin(0.1, 15, 0.1, 0.1, "cm"),
    legend.position = "none")+ 
  coord_cartesian(clip = "off")
ggsave("empAr.png", width = 15, height = 12, dpi = 600)



########### BETULA #################
# Read the tree
tAs <- read.tree(empKs)
label_map <- setNames(sapply(tAs$tip.label, get_scientific_name), tAs$tip.label)
tAs$tip.label <- label_map[tAs$tip.label]

ggtree(tAs) +                             
  geom_tiplab(aes(label = label, 
                  color = grepl("^Kap", label)), 
              size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  theme_tree2(base_size = 16) +
  geom_nodelab(geom = 'label', size = 3) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with strict mode in keelime") +
  theme(
    plot.margin = margin(0.1, 3, 0.1, 0.1, "cm"),
    legend.position = "none") +
  coord_cartesian(clip = "off")

ggsave("empKs.png", width = 15, height = 12, dpi = 600)

tAn <- read.tree(empKn)
label_map <- setNames(sapply(tAn$tip.label, get_scientific_name), tAn$tip.label)
tAn$tip.label <- label_map[tAn$tip.label]
ggtree(tAn) +                             
  geom_tiplab(aes(label = label, 
                  color = grepl("^Kap", label)), 
              size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  theme_tree2(base_size = 16) +
  geom_nodelab(geom = 'label', size = 3) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with strict mode in keelime") +
  theme(
    plot.margin = margin(0.1, 3, 0.1, 0.1, "cm"),
    legend.position = "none") +
  coord_cartesian(clip = "off")

ggsave("empKn.png", width = 15, height = 12, dpi = 600)

tAn <- read.tree(empKr)
label_map <- setNames(sapply(tAn$tip.label, get_scientific_name), tAn$tip.label)
tAn$tip.label <- label_map[tAn$tip.label]
ggtree(tAn) +                             
  geom_tiplab(aes(label = label, 
                  color = grepl("^Kap", label)), 
              size = 4, hjust = 0.0005) + 
  geom_tippoint() +
  theme_tree2(base_size = 16) +
  geom_nodelab(geom = 'label', size = 3) + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  ggtitle("Consensus call with strict mode in keelime") +
  theme(
    plot.margin = margin(0.1, 3, 0.1, 0.1, "cm"),
    legend.position = "none") +
  coord_cartesian(clip = "off")

ggsave("empKr.png", width = 15, height = 12, dpi = 600)


#############################################################
library(ape)
library(ggtree)
library(ggplot2)
library(rentrez)  # needed if using get_scientific_name

# Load trees
E <- "/Users/bfj994/Documents/keelime/empTrees/Elephantidae.new.dnd"
E_w <- "/Users/bfj994/Documents/keelime/empTrees/Elephantidae_without.new.dnd"
Et <- read.tree(E)
Et_m <- midpoint.root(Et)
Ewt <- read.tree(E_w)

# Extract accession number
extract_accession <- function(label) {
  sub("([A-Z]{1,2}_?[0-9]+\\.[0-9]+).*", "\\1", label)
}

Et_accessions <- sapply(Et$tip.label, extract_accession)
Ewt_accessions <- sapply(Ewt$tip.label, extract_accession)
missing_accessions <- setdiff(Et_accessions, Ewt_accessions)
missing_tips <- Et$tip.label[Et_accessions %in% missing_accessions]

# Fetch scientific names
label_map <- setNames(sapply(Et$tip.label, get_scientific_name), Et$tip.label)
Et$tip.label <- label_map[Et$tip.label]
missing_tips_named <- label_map[missing_tips]
# List of labels to highlight in bold
bold_tips <- c("Mammut americanum", "Elephas maximus")  # Replace these with the actual labels of the tips you want to bold

# Build base plot
ggtree(Et) +
  geom_tiplab(aes(color = label %in% missing_tips_named, 
                  fontface = ifelse(label %in% bold_tips, "bold", "plain")), 
              size = 5, hjust = 0) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  theme(legend.position = "none",
        plot.margin = margin(0.5, 7, 0.5, 0.5, "cm")) +
  coord_cartesian(clip = "off")

ggsave("supp_tree_Elephantidae.png", height = 8, width = 9, dpi = 600)

E <- "/Users/bfj994/Documents/keelime/empTrees/Ursidae.new.dnd"
E_w <- "/Users/bfj994/Documents/keelime/empTrees/Ursidae_without.new.dnd"
Et <- read.tree(E)
Et_m <- midpoint.root(Et)
Ewt <- read.tree(E_w)

# Extract accession number
extract_accession <- function(label) {
  sub("([A-Z]{1,2}_?[0-9]+\\.[0-9]+).*", "\\1", label)
}

Et_accessions <- sapply(Et$tip.label, extract_accession)
Ewt_accessions <- sapply(Ewt$tip.label, extract_accession)
missing_accessions <- setdiff(Et_accessions, Ewt_accessions)
missing_tips <- Et$tip.label[Et_accessions %in% missing_accessions]

# Fetch scientific names
label_map <- setNames(sapply(Et$tip.label, get_scientific_name), Et$tip.label)
Et_m$tip.label <- label_map[Et$tip.label]
missing_tips_named <- label_map[missing_tips]
# List of labels to highlight in bold
bold_tips <- c("Ursus arctos", "Ursus spelaeus")  # Replace these with the actual labels of the tips you want to bold

# Build base plot
ggtree(Et_m) +
  geom_tiplab(aes(color = label %in% missing_tips_named, 
                  fontface = ifelse(label %in% bold_tips, "bold", "plain")), 
              size = 5, hjust = 0) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  theme(legend.position = "none",
        plot.margin = margin(0.5, 7, 0.5, 0.5, "cm")) +
  coord_cartesian(clip = "off")

ggsave("supp_tree_Ursidae.png", height = 8, )


E <- "/Users/bfj994/Documents/keelime/empTrees/Sciuridae.new.dnd"
E_w <- "/Users/bfj994/Documents/keelime/empTrees/Sciuridae_without.new.dnd"
Et <- read.tree(E)
Et_m <- midpoint.root(Et)
Ewt <- read.tree(E_w)

# Extract accession number
extract_accession <- function(label) {
  sub("([A-Z]{1,2}_?[0-9]+\\.[0-9]+).*", "\\1", label)
}

Et_accessions <- sapply(Et$tip.label, extract_accession)
Ewt_accessions <- sapply(Ewt$tip.label, extract_accession)
missing_accessions <- setdiff(Et_accessions, Ewt_accessions)
missing_tips <- Et$tip.label[Et_accessions %in% missing_accessions]

# Fetch scientific names
label_map <- setNames(sapply(Et$tip.label, get_scientific_name), Et$tip.label)
Et_m$tip.label <- label_map[Et$tip.label]
missing_tips_named <- label_map[missing_tips]
# List of labels to highlight in bold
bold_tips <- c("Marmota baibacina", "Marmota himalayana")  # Replace these with the actual labels of the tips you want to bold

# Build base plot
ggtree(Et_m) +
  geom_tiplab(aes(color = label %in% missing_tips_named, 
                  fontface = ifelse(label %in% bold_tips, "bold", "plain")), 
              size = 5, hjust = 0) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  geom_treescale(x = 0, y = -1, width = 0.03) +
  theme_tree2(base_size = 16) +
  theme(legend.position = "none",
        plot.margin = margin(0.5, 7, 0.5, 0.5, "cm")) +
  coord_cartesian(clip = "off")

ggsave("supp_tree_Sc.png", height = 20, width = 15, dpi = 600)

E <- "/Users/bfj994/Documents/keelime/Betula.new.dnd"
E_w <- "/Users/bfj994/Documents/keelime/Betula_wo.new.dnd"
Et <- read.tree(E)
Et_m <- midpoint.root(Et)
Ewt <- read.tree(E_w)

# Extract accession number
extract_accession <- function(label) {
  sub("([A-Z]{1,2}_?[0-9]+\\.[0-9]+).*", "\\1", label)
}

Et_accessions <- sapply(Et$tip.label, extract_accession)
Ewt_accessions <- sapply(Ewt$tip.label, extract_accession)
missing_accessions <- setdiff(Et_accessions, Ewt_accessions)
missing_tips <- Et$tip.label[Et_accessions %in% missing_accessions]

# Fetch scientific names
label_map <- setNames(sapply(Et$tip.label, get_scientific_name), Et$tip.label)
Et_m$tip.label <- label_map[Et$tip.label]
label_map[label_map == "FetchError"] <- "Betula platyphylla"
missing_tips_named <- label_map[missing_tips]
# List of labels to highlight in bold
bold_tips <- c("Betula nana", "Betula pubescens")  # Replace these with the actual labels of the tips you want to bold

# Build base plot
ggtree(Et_m) +
  geom_tiplab(aes(color = label %in% missing_tips_named, 
                  fontface = ifelse(label %in% bold_tips, "bold", "plain")), 
              size = 5, hjust = 0) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  geom_treescale(x = 0, y = -1, width = 0.003) +
  theme_tree2(base_size = 16) +
  theme(legend.position = "none",
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")) +
  coord_cartesian(clip = "off")

ggsave("supp_tree_Betula.png", height = 20, width = 15, dpi = 600)



