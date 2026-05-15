## Final figure ----
### Data ----

require(dplyr)
url <- "https://raw.githubusercontent.com/breesie/growthadvantagecompare/main/sequences_static/metadata.tsv"
metadata <- read.delim(url, sep = "\t", header = TRUE)

#https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid/nextstrain_clades/global.tsv.gz
url <- "https://raw.githubusercontent.com/breesie/growthadvantagecompare/refs/heads/main/sequences_static/global-2.tsv"
clades <- read.delim(url, sep = "\t", header = TRUE)

#https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid/pango_lineages/global.tsv.gz
url <- "https://raw.githubusercontent.com/breesie/growthadvantagecompare/refs/heads/main/sequences_static/global.tsv"
pangos <- read.delim(url, sep = "\t", header = TRUE)

cladenum <- clades %>% rename(name = clade)

pangonum <- pangos %>% 
  mutate(variant = "pango") %>% 
  rename(name = clade)

metadata_clade_pango <- metadata %>%
  group_by(clade_membership) %>%
  count(pango_lineage) %>% 
  mutate(clade_membership = stringr::str_extract(clade_membership, "^[0-9]+[A-Z]"))

#old clade/ pango pairings
# names <- tibble(
#   clade = c("21K", "21L", "22B", "23A", "24E", "24A", "25C"),
#   pango = c("BA.1", "BA.2", "BA.5", "XBB.1.5", "KP.3.1.1", "JN.1", "XFG")
# )

cladetotn_usa <- cladenum %>% filter(location == "USA") %>% 
  group_by(date) %>% summarise(totn = sum(sequences)) %>% 
  mutate("variant" = "clade")

pangototn_usa <- pangonum %>% filter(location == "USA") %>% 
  group_by(date) %>% summarise(totn = sum(sequences)) %>% 
  mutate("variant" = "pango")

pangoprops_usa <- left_join(pangototn_usa, pangonum %>% filter(location == "USA"), by = "date") %>% 
  mutate(pangoprop = sequences/totn) 

cladeprops_usa <- left_join(cladetotn_usa, cladenum %>% filter(location == "USA"), by = "date") %>%  
  mutate(cladeprop = sequences/totn) 


# Figure 1 latest 
library(tidyverse)
library(scales)
require(zoo)


# Base clade colors
clade_colors <- c(
  "21K" = "#1b9e77",
  "21L" = "#d95f02",
  "22B" = "#7570b3",
  "23A" = "#e7298a",
  "24E" = "#66a61e"
)

clades <- names(clade_colors)

# Build data (your logic)
get_clade_data <- function(clade) {
  
  mindate <- cladeprops_usa %>% 
    filter(name == clade) %>% 
    filter(cladeprop == max(cladeprop, na.rm = TRUE)) %>% 
    summarise(mindate = min(as.Date(date)) - 150) %>% 
    pull(mindate)
  
  p_in_c <- metadata_clade_pango %>% 
    filter(clade_membership == clade) %>% 
    pull(pango_lineage)
  
  k = 7 # 7 day rolling mean
  
  pangos <- pangoprops_usa %>% 
    group_by(name) %>% 
    mutate(pangoprop = zoo::rollmean(pangoprop, k = k, fill = NA, align = "center")) %>% 
    ungroup() %>% 
    filter(name %in% p_in_c,
           date >= mindate,
           date <= mindate + 300) %>%
    mutate(clade = clade)
  
  clade_line <- cladeprops_usa %>% 
    group_by(name) %>% 
    mutate(cladeprop = zoo::rollmean(cladeprop, k = k, fill = NA, align = "center")) %>% 
    ungroup() %>% 
    filter(name == clade,
           date >= mindate,
           date <= mindate + 300) %>%
    mutate(clade = clade)
  
  list(pangos = pangos, clade_line = clade_line)
}

data_list <- map(clades, get_clade_data)

pangos_all <- map_dfr(data_list, "pangos") %>% 
  mutate(date=as.Date(date))
clade_all  <- map_dfr(data_list, "clade_line") %>% 
  mutate(date=as.Date(date))

# Rank lineages within clade
lineage_rank <- pangos_all %>%
  group_by(clade, name) %>%
  summarise(max_prop = max(pangoprop, na.rm = TRUE), .groups = "drop") %>%
  group_by(clade) %>%
  mutate(rank = rank(-max_prop, ties.method = "first")) %>%
  ungroup()

pangos_all <- pangos_all %>%
  left_join(lineage_rank, by = c("clade", "name"))



# Build gradient colors (FIXED + STRONG TOP 3 CONTRAST)
pango_color_vector <- c()

for (clade in clades) {
  
  base_color <- clade_colors[clade]
  
  df <- lineage_rank %>% 
    filter(clade == !!clade) %>%
    arrange(rank)
  
  n <- nrow(df)
  
  cols <- character(n)
  
  for (i in seq_len(n)) {
    
    r <- df$rank[i]
    
    if (r == 1) {
      
      # MOST DOMINANT → darkest + most saturated
      cols[i] <- colorspace::lighten(base_color, 0.1)
      
    } else if (r == 2) {
      
      # SECOND → slightly lighter
      cols[i] <- colorspace::lighten(base_color, 0.40)
      
    } else if (r == 3) {
      
      # THIRD → neutral mid tone
      cols[i] <- colorspace::lighten(base_color, 0.70)
      
    } else {
      
      # REST → fade out quickly
      x <- (r - 3) / max(1, (n - 3))
      
      cols[i] <- colorspace::lighten(
        colorspace::desaturate(base_color, amount = 0.4 * x),
        amount = 0.2 + 0.7 * x
      )
    }
  }
  
  names(cols) <- df$name
  
  pango_color_vector <- c(pango_color_vector, cols)
}

# safety cleanup
pango_color_vector <- pango_color_vector[!duplicated(names(pango_color_vector))]




# Line width mapping
pangos_all <- pangos_all %>%
  mutate(linewidth = scales::rescale(rank,
                                     to = c(0.4, 1.2),
                                     from = range(rank, na.rm = TRUE))) %>% 
  mutate(date=as.Date(date))

# Line width mapping
# invert so rank 1 (largest) = thickest
pangos_all$linewidth <- max(pangos_all$linewidth) - pangos_all$linewidth + min(pangos_all$linewidth)

# Top lineages for legend
top_n <- 3

top_lineages <- lineage_rank %>%
  group_by(clade) %>%
  slice_min(rank, n = top_n) %>%
  pull(name)



pangos_all_top <- 
  pangos_all %>%
  filter(name %in% top_lineages) %>% 
  mutate(date=as.Date(date))


### Plot ----

clade_labels <- c(
  "21K" = "21K (BA.1)",
  "21L" = "21L (BA.2)",
  "22B" = "22B (BA.5)",
  "23A" = "23A (XBB.1.5)",
  "24E" = "24E (KP.3.1.1)",
  "24A" = "24A (JN.1)",
  "25C" = "25C (XFG)"
)

# Top pangos per clade (already ranked)
top_by_clade <- lineage_rank %>%
  group_by(clade) %>%
  arrange(rank, .by_group = TRUE) %>%
  slice_min(rank, n = top_n) %>%
  summarise(pangos = list(name), .groups = "drop")

# Interleave: each clade followed by its own top pangos
legend_breaks <- purrr::flatten_chr(purrr::map(clades, function(cl) {
  pangos_for_clade <- top_by_clade$pangos[top_by_clade$clade == cl]
  c(cl, if (length(pangos_for_clade)) pangos_for_clade[[1]] else character(0))
}))

# Labels: clades use friendly names; pangos as-is
label_map <- c(setNames(top_lineages, top_lineages), clade_labels)
legend_labels <- label_map[legend_breaks]

figure_1_real <- 
  
  ggplot() +
  
  geom_line(
    data = pangos_all_top,
    aes(x = date,
        y = pangoprop,
        group = name,
        color = name,
        linewidth = .4),
    alpha = 0.7
  ) +
  
  scale_color_manual(
    values = c(pango_color_vector, clade_colors),
    breaks = legend_breaks,
    labels = legend_labels
  ) +
  
  scale_linewidth_identity() +
  
  # CLade line (strongest)
  geom_line(
    data = clade_all,
    aes(x = date, y = cladeprop,
        group = name, color = clade),
    linewidth = 1
  ) +
  
  facet_wrap(~dplyr::recode(clade, !!!clade_labels), scales = "free_x") +
  
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "3 months",
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  
  labs(
    title = "Clades and Constituent Pango Lineages in the United States",
    x = "Date",
    y = "Proportion of sequences",
    color = "Clades and \nmost common \nPango lineages"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


# Size and save figure ----
figure_1_real +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  )

ggsave("figure_1.jpg",
       width = 6,
       height = 3,
       units = c("in"))