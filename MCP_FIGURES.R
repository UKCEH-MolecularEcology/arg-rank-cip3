required_packages <- c(
  "tidyverse",
  "vegan",
  "cowplot",
  "patchwork",
  "rstatix",
  "readxl",
  "readr",
  "ggsignif",
  "ggpubr",
  "broom",
  "broom.mixed",
  "ggrain",
  "lme4",
  "lmerTest",
  "ggtext",
  "grid",
  "scales"
)

for(pkg in required_packages){

  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }

  library(pkg, character.only = TRUE)
}

# Visual settings

base_font <- "Helvetica"

met_cols <- list(

  stage = c(
    "Influent" = "#C44E52",
    "Effluent" = "#4C72B0"
  ),

  stage_lower = c(
    "influent" = "#C44E52",
    "final effluent" = "#4C72B0",
    "effluent" = "#4C72B0"
  ),

  delta = c(
    "#d73027",
    "#fdae61",
    "#abd9e9",
    "#4575b4"
  ),

  neutral = c(
    "dark" = "#2B2B2B",
    "mid"  = "#7A7A7A",
    "light" = "#D9D9D9"
  )
)

# Global theme

theme_met <- function(){

  theme_bw(base_family = base_font) +

    theme(

      text = element_text(
        color = "black",
        size = 10
      ),

      axis.title = element_text(
        size = 12,
        face = "bold"
      ),

      axis.text = element_text(
        size = 10
      ),

      strip.text = element_text(
        size = 12,
        face = "bold"
      ),

      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.8
      ),

      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      legend.title = element_text(
        size = 11,
        face = "bold"
      ),

      legend.text = element_text(
        size = 10
      ),

      legend.background = element_blank(),

      plot.title = element_text(
        size = 14,
        face = "bold"
      ),

      plot.tag = element_text(
        size = 16,
        face = "bold"
      ),

      plot.margin = margin(
        8, 8, 8, 8
      )
    )
}

theme_set(theme_met())

# Helper functions

save_fig <- function(plot,
                     filename,
                     width,
                     height,
                     dpi = 600){

  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

clean_numeric <- function(x){

  x <- str_trim(as.character(x))

  is_censored <- str_detect(x, "^<")
  is_censored[is.na(is_censored)] <- FALSE

  num <- suppressWarnings(
    as.numeric(
      str_remove(x, "^<\\s*")
    )
  )

  num[is_censored] <- num[is_censored] / 2

  return(num)
}

theme_met_classic <- function(){

  theme_classic(base_family = base_font) +

    theme(

      text = element_text(
        color = "black"
      ),

      axis.title = element_text(
        size = 12,
        face = "bold"
      ),

      axis.text = element_text(
        size = 10
      ),

      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),

      legend.title = element_blank(),

      plot.tag = element_text(
        size = 16,
        face = "bold"
      )
    )
}

# Figure 1 -- alpha diversity and PCoA


arg_long <- read.csv(
  "all_ARG_counts.csv",
  stringsAsFactors = FALSE
)

arg_wide <- arg_long %>%
  pivot_wider(
    names_from = ARG,
    values_from = Count,
    values_fill = 0
  )

arg <- as.data.frame(arg_wide)

rownames(arg) <- arg$Sample
arg$Sample <- NULL

arg <- arg[rowSums(arg) > 0, ]

arg_rel <- decostand(
  arg,
  method = "total"
)

prevalence <- colSums(arg_rel > 0) / nrow(arg_rel)

arg_rel <- arg_rel[
  ,
  prevalence >= 0.05
]

meta2 <- data.frame(
  sample_id = c(
    "E101","E105","E111","E115","E121","E127","E133","E139","E144","E149",
    "E157","E161","E167","E172","E178","E184","E190","E196","E201","E205",
    "E317","E321","E326","E332","E338","E344","E350","E356","E361","E365",
    "E372","E376","E382","E387","E393","E399","E405","E411","E416","E420",
    "E424","E426","E433","E439","E445","E451","E457","E463","E468","E472",
    "E480","E484","E496","E502","E509","E520",
    "E108","E114","E117","E124","E129","E136","E141","E145","E153","E164",
    "E170","E174","E181","E186","E193","E198","E202","E209","E324","E329",
    "E334","E341","E346","E353","E358","E362","E369","E379","E385","E389",
    "E396","E401","E408","E413","E417","E425","E430","E436","E441","E448",
    "E453","E460","E465","E469","E476","E487","E498","E506","E511","E517","E522"
  ),
  Group.y = c(
    rep("influent", 56),
    rep("final effluent", 51)
  )
)

meta2$Group.y <- factor(
  meta2$Group.y,
  levels = c(
    "influent",
    "final effluent"
  )
)

common_ids <- intersect(
  rownames(arg_rel),
  meta2$sample_id
)

arg_rel <- arg_rel[common_ids, ]

meta2 <- meta2[
  match(common_ids, meta2$sample_id),
]

# Alpha diversity

shannon <- diversity(
  arg_rel,
  index = "shannon"
)

diversity_df <- data.frame(
  sample_id = rownames(arg_rel),
  shannon = shannon
) %>%
  left_join(meta2, by = "sample_id")

# Beta diversity

arg_hel <- decostand(
  arg_rel,
  method = "hellinger"
)

bray_dist <- vegdist(
  arg_hel,
  method = "bray"
)

pcoa <- cmdscale(
  bray_dist,
  eig = TRUE,
  k = 2
)

eig <- pcoa$eig

var_exp <- round(
  100 * eig / sum(eig),
  2
)

pcoa_points <- as.data.frame(
  pcoa$points
)

colnames(pcoa_points) <- c(
  "PCoA1",
  "PCoA2"
)

pcoa_points$PCoA1 <- -pcoa_points$PCoA1
pcoa_points$PCoA2 <- -pcoa_points$PCoA2

pcoa_points$sample_id <- rownames(pcoa_points)

plot_data <- left_join(
  pcoa_points,
  meta2,
  by = "sample_id"
)

# Panel A

p_alpha <- ggplot(
  diversity_df,
  aes(
    x = Group.y,
    y = shannon,
    fill = Group.y
  )
) +

  geom_violin(
    width = 0.7,
    alpha = 0.35,
    color = NA,
    trim = FALSE
  ) +

  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.8,
    color = "black"
  ) +

  geom_jitter(
    width = 0.08,
    size = 1.5,
    alpha = 0.5,
    color = "black"
  ) +

  scale_fill_manual(
    values = met_cols$stage_lower
  ) +

  labs(
    x = NULL,
    y = "Shannon diversity"
  ) +

  theme_met_classic() +

  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(
      color = "grey92",
      linewidth = 0.4
    )
  )

# Panel B

p_pcoa <- ggplot(
  plot_data,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Group.y
  )
) +

  geom_point(
    size = 3,
    alpha = 0.8
  ) +

  stat_ellipse(
    type = "t",
    linewidth = 1,
    linetype = 2,
    alpha = 0.3
  ) +

  scale_color_manual(
    values = met_cols$stage_lower
  ) +

  labs(
    x = paste0(
      "PCoA1 (",
      var_exp[1],
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      var_exp[2],
      "%)"
    )
  ) +

  theme_met_classic() +

  theme(
    legend.position = "top",
    panel.grid.major = element_line(
      color = "grey92",
      linewidth = 0.4
    )
  )

fig1 <- p_alpha + p_pcoa +
  plot_annotation(tag_levels = "A")

save_fig(
  fig1,
  "METFIG1.png",
  width = 14,
  height = 6
)

# Figure 2 -- MC2 Risk Scores


df <- read_excel("full_merged_dataset.xlsx") %>%
  rename(
    PlantID = Site,
    MonthID = `Sample Round`,
    Sample = sample_id,
    Group = Group.x
  ) %>%
  mutate(
    MonthID = factor(MonthID),

    Group = case_when(
      str_to_lower(Group) == "influent" ~ "Influent",
      TRUE ~ "Effluent"
    ),

    ScoreType = str_to_title(str_trim(ScoreType)),

    ScoreTypeLabel = recode(
      ScoreType,
      "Ecological" = "ERR",
      "Human Health" = "HHRR"
    )
  ) %>%
  filter(
    Group %in% c("Influent", "Effluent"),
    ScoreTypeLabel %in% c("ERR", "HHRR")
  )

df$Group <- factor(
  df$Group,
  levels = c("Influent", "Effluent")
)

run_mixed_model <- function(score_label){

  sub <- df %>%
    filter(ScoreTypeLabel == score_label)

  model <- lmer(
    RiskScore ~ Group + (1 | PlantID / MonthID),
    data = sub
  )

  tidy(model, effects = "fixed") %>%
    filter(term == "GroupEffluent") %>%
    mutate(ScoreTypeLabel = score_label)
}

mixed_ERR  <- run_mixed_model("ERR")
mixed_HHRR <- run_mixed_model("HHRR")

mixed_results <- bind_rows(
  mixed_ERR,
  mixed_HHRR
)

print("Mixed-effects model results")
print(mixed_results)

rain_plot <- function(score){

  ggplot(
    df %>% filter(ScoreTypeLabel == score),

    aes(
      x = Group,
      y = RiskScore,
      fill = Group
    )
  ) +

    geom_rain(
      alpha = 0.5,
      rain.side = "r",

      violin.args = list(
        alpha = 0.4
      ),

      boxplot.args = list(
        color = "black",
        linewidth = 0.6,
        outlier.shape = NA
      ),

      point.args = list(
        size = 0.7,
        alpha = 0.35
      )
    ) +

    scale_fill_manual(
      values = met_cols$stage
    ) +

    labs(
      x = NULL,
      y = "Risk score"
    ) +

    theme(
      legend.position = "none"
    )
}

rain_ERR <- rain_plot("ERR")
rain_HHRR <- rain_plot("HHRR")

paired <- df %>%
  select(
    PlantID,
    MonthID,
    ScoreTypeLabel,
    Group,
    RiskScore
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = RiskScore
  ) %>%
  filter(
    !is.na(Influent),
    !is.na(Effluent)
  ) %>%
  mutate(
    Delta = Effluent - Influent
  )

paired_stats <- function(score_label){

  sub <- paired %>%
    filter(ScoreTypeLabel == score_label)

  wil <- wilcox.test(
    sub$Effluent,
    sub$Influent,
    paired = TRUE
  )

  tibble(
    ScoreTypeLabel = score_label,
    Median_Influent = median(sub$Influent),
    Median_Effluent = median(sub$Effluent),
    Median_Delta = median(sub$Delta),
    Wilcoxon_p = wil$p.value
  )
}

paired_ERR  <- paired_stats("ERR")
paired_HHRR <- paired_stats("HHRR")

paired_results <- bind_rows(
  paired_ERR,
  paired_HHRR
)

print("Paired influent–effluent results")
print(paired_results)

delta_max <- max(abs(paired$Delta))

plot_score <- function(score_type){

  sub <- paired %>%
    filter(ScoreTypeLabel == score_type)

  slopes <- ggplot(sub) +

    geom_segment(
      aes(
        x = 0,
        xend = 1,
        y = Influent,
        yend = Effluent,
        color = Delta
      ),

      size = 0.9,
      alpha = 0.9,

      arrow = arrow(
        length = unit(0.1, "cm")
      )
    ) +

    geom_point(
      aes(
        x = 0,
        y = Influent,
        color = Delta
      ),
      size = 2
    ) +

    geom_point(
      aes(
        x = 1,
        y = Effluent,
        color = Delta
      ),
      size = 2
    ) +

    scale_x_continuous(
      breaks = c(0, 1),
      labels = c(
        "Influent",
        "Effluent"
      )
    ) +

    scale_color_gradientn(
      colors = met_cols$delta,
      limits = c(
        -delta_max,
        delta_max
      )
    ) +

    labs(
      x = NULL,
      y = "Risk score"
    )

  box <- ggplot(
    sub,
    aes(
      y = Delta,
      x = 1
    )
  ) +

    geom_violin(
      fill = met_cols$neutral["light"],
      alpha = 0.7,
      color = NA,
      trim = FALSE
    ) +

    geom_boxplot(
      width = 0.2,
      fill = "white",
      color = "black",
      linewidth = 0.6,
      outlier.shape = NA
    ) +

    geom_jitter(
      width = 0.08,
      size = 1.3,
      alpha = 0.6
    ) +

    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = met_cols$neutral["mid"]
    ) +

    labs(
      x = NULL,
      y = expression(
        Delta ~
          "Risk score (Effluent - Influent)"
      )
    ) +

    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  list(
    slopes_plot = slopes,
    box_plot = box
  )
}

plots <- list(
  plot_score("ERR"),
  plot_score("HHRR")
)

fig2 <-

  (rain_ERR +
     plots[[1]]$slopes_plot +
     plots[[1]]$box_plot) /

  (rain_HHRR +
     plots[[2]]$slopes_plot +
     plots[[2]]$box_plot) +

  plot_layout(
    widths = c(
      1.1,
      1.2,
      0.8
    )
  ) +

  plot_annotation(
    tag_levels = "A"
  )

save_fig(
  fig2,
  "METFIG2.png",
  width = 13,
  height = 6.5
)


for(row in 1:nrow(mixed_results)){

  r <- mixed_results[row, ]

  cat(
    r$ScoreTypeLabel,
    ": β = ",
    round(r$estimate, 2),
    ", p = ",
    format(r$p.value, scientific = TRUE),
    "\n",
    sep = ""
  )
}

for(row in 1:nrow(paired_results)){

  r <- paired_results[row, ]

  cat(
    r$ScoreTypeLabel,
    ": median Δ = ",
    round(r$Median_Delta, 2),
    ", Wilcoxon p = ",
    format(r$Wilcoxon_p, scientific = TRUE),
    "\n",
    sep = ""
  )
}


# Figure 3 -- Site-Level Risk Scores

df_raw <- read_csv(
  "arg_metacompare_scores.csv",
  show_col_types = FALSE
)

df <- df_raw %>%
  select(
    sample_id,
    ScoreType,
    Group,
    RiskScore
  ) %>%
  mutate(

    Group = tolower(Group),

    Group = case_when(
      Group == "influent" ~ "Influent",
      Group == "final effluent" ~ "Effluent",
      TRUE ~ Group
    ),

    ScoreType = tools::toTitleCase(
      trimws(ScoreType)
    ),

    ScoreType = recode(
      ScoreType,
      "Ecological" = "ERR",
      "Human Health" = "HHRR"
    ),

    Group = factor(
      Group,
      levels = c(
        "Influent",
        "Effluent"
      )
    )
  )



site_mapping <- c(
  "E105"="Driffield","E108"="Driffield","E111"="Exmouth",
  "E114"="Exmouth","E115"="Finham","E117"="Finham",
  "E121"="Pembury","E124"="Pembury","E127"="Pen-y-bont",
  "E129"="Pen-y-bont","E133"="Saltford","E136"="Saltford",
  "E139"="Seaton Carew","E141"="Seaton Carew",
  "E144"="Stockport","E145"="Stockport",
  "E149"="Chalton","E153"="Chalton",
  "E161"="Driffield","E164"="Driffield",
  "E167"="Exmouth","E170"="Exmouth",
  "E172"="Finham","E174"="Finham",
  "E178"="Pembury","E181"="Pembury",
  "E184"="Pen-y-bont","E186"="Pen-y-bont",
  "E190"="Saltford","E193"="Saltford",
  "E196"="Seaton Carew","E198"="Seaton Carew",
  "E201"="Stockport","E202"="Stockport",
  "E205"="Chalton","E209"="Chalton",
  "E321"="Driffield","E324"="Driffield",
  "E326"="Exmouth","E329"="Exmouth",
  "E332"="Finham","E334"="Finham",
  "E338"="Pembury","E341"="Pembury",
  "E344"="Pen-y-bont","E346"="Pen-y-bont",
  "E350"="Saltford","E353"="Saltford",
  "E356"="Seaton Carew","E358"="Seaton Carew",
  "E361"="Stockport","E362"="Stockport",
  "E365"="Chalton","E369"="Chalton",
  "E376"="Driffield","E379"="Driffield",
  "E382"="Exmouth","E385"="Exmouth",
  "E387"="Finham","E389"="Finham",
  "E393"="Pembury","E396"="Pembury",
  "E399"="Pen-y-bont","E401"="Pen-y-bont",
  "E405"="Saltford","E408"="Saltford",
  "E411"="Seaton Carew","E413"="Seaton Carew",
  "E416"="Stockport","E417"="Stockport",
  "E420"="Chalton","E425"="Chalton",
  "E426"="Driffield","E430"="Driffield",
  "E433"="Exmouth","E436"="Exmouth",
  "E439"="Finham","E441"="Finham",
  "E445"="Pembury","E448"="Pembury",
  "E451"="Pen-y-bont","E453"="Pen-y-bont",
  "E457"="Saltford","E460"="Saltford",
  "E463"="Seaton Carew","E465"="Seaton Carew",
  "E468"="Stockport","E469"="Stockport",
  "E472"="Chalton","E476"="Chalton",
  "E484"="Driffield","E487"="Driffield",
  "E496"="Finham","E498"="Finham",
  "E502"="Pembury","E506"="Pembury",
  "E509"="Pen-y-bont","E511"="Pen-y-bont",
  "E517"="Saltford","E520"="Seaton Carew",
  "E522"="Seaton Carew"
)

df$Site <- site_mapping[df$sample_id]

df <- df %>%
  filter(!is.na(Site))

run_mixed <- function(score){

  sub <- df %>%
    filter(ScoreType == score)

  model <- lmer(
    RiskScore ~ Group + (1 | Site),
    data = sub
  )

  tidy(model, effects = "fixed") %>%
    filter(term == "GroupEffluent")
}

mixed_ERR  <- run_mixed("ERR")
mixed_HHRR <- run_mixed("HHRR")

mixed_results <- bind_rows(
  mixed_ERR,
  mixed_HHRR
)

print("Mixed-effects model results")
print(mixed_results)

wilcox_by_site <- df %>%
  group_by(Site, ScoreType) %>%
  summarise(
    p = tryCatch(
      wilcox.test(RiskScore ~ Group)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  )

print("=== Per-site Wilcoxon tests ===")
print(wilcox_by_site)

df_for_signif <- df %>%
  group_by(Site, ScoreType) %>%
  summarise(
    p = tryCatch(
      wilcox.test(RiskScore ~ Group)$p.value,
      error = function(e) NA
    ),

    y = max(RiskScore, na.rm = TRUE) * 1.05,

    .groups = "drop"
  ) %>%
  mutate(
    group1 = "Influent",
    group2 = "Effluent",
    label = sprintf("p = %.3g", p)
  )

p3 <- ggplot(
  df,
  aes(
    x = Group,
    y = RiskScore,
    fill = Group
  )
) +

  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    color = "black"
  ) +

  geom_jitter(
    position = position_jitter(width = 0.2),
    size = 2,
    alpha = 0.8
  ) +

  facet_grid(
    ScoreType ~ Site,
    scales = "free_y",
    switch = "y"
  ) +

  scale_fill_manual(
    values = met_cols$stage
  ) +

  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = met_cols$neutral["mid"]
  ) +

  stat_pvalue_manual(
    df_for_signif,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y",
    tip.length = 0
  ) +

  labs(
    x = NULL,
    y = "Risk score"
  ) +

  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    ),
    legend.position = "bottom",
    legend.title = element_blank()
  )

save_fig(
  p3,
  "METFIG3.png",
  width = 16,
  height = 10
)

cat("\n================ SUMMARY ================\n")

print(mixed_results)
print(wilcox_by_site)


# Figure 4 -- Physicochemical vs RiskScore


full_merged_dataset <- read_excel(
  "full_merged_dataset.xlsx"
)

log_transform_metrics <- c(
  "Alkalinity_as_CaCO3_(mg/l)",
  "Ammoniacal_Nitrogen_as_N_(mg/l)",
  "Biochemical_Oxygen_Demand_Total_(mg/l)",
  "Chemical_Oxygen_Demand_Total_(mg/l)",
  "Dissolved_ammonium_(NH4)_(mg/l)",
  "Dissolved_chloride_(mg_Cl/L)",
  "Dissolved_fluoride_(mg_F/L)",
  "Dissolved_nitrate_(NO3)",
  "Dissolved_nitrite_(mg_NO2/L)",
  "Dissolved_organic_carbon_(mg/L)",
  "Dissolved_sulphate_(mg_SO4/L)",
  "Soluble_reactive_phosphorus_(?g/L)",
  "Suspended_Solids_(mg/l)",
  "Suspended_solids_mg/L",
  "Total_dissolved_nitrogen_(mg_N/L)",
  "Total_dissolved_phosphorus_(?g/L)",
  "Total_phosphorus_(?g/L)"
)

physchem_cols <- c(
  "pH",
  log_transform_metrics
)

existing_cols <- physchem_cols[
  physchem_cols %in%
    colnames(full_merged_dataset)
]

full_merged_dataset <- full_merged_dataset %>%
  mutate(
    across(
      all_of(existing_cols),
      clean_numeric
    )
  )

full_merged_dataset <- full_merged_dataset %>%
  mutate(
    Group.x = case_when(
      Group.x == "influent" ~ "influent",
      Group.x == "final effluent" ~ "effluent",
      TRUE ~ Group.x
    ),

    Group.x = factor(
      Group.x,
      levels = c(
        "influent",
        "effluent"
      )
    )
  )

metric_labels <- c(
  "pH" = "pH"
)

full_merged_dataset <- full_merged_dataset %>%
  mutate(
    ScoreType = recode(
      ScoreType,
      "Ecological" = "ERR",
      "Human health" = "HHRR"
    )
  )

make_sig_physchem_plot <- function(df,
                                   score_type,
                                   y_limits = NULL){

  long_df <- df %>%
    filter(ScoreType == score_type) %>%
    pivot_longer(
      cols = all_of(existing_cols),
      names_to = "PhysChem_Metric",
      values_to = "PhysChem_Value"
    ) %>%
    filter(!is.na(PhysChem_Value)) %>%
    mutate(
      Value_log = case_when(
        PhysChem_Metric %in% log_transform_metrics ~
          log10(pmax(PhysChem_Value, 1e-3)),
        TRUE ~ PhysChem_Value
      )
    )

  cor_results <- long_df %>%
    group_by(
      PhysChem_Metric,
      Group.x
    ) %>%
    summarise(
      test = list(
        cor.test(
          Value_log,
          RiskScore,
          method = "pearson"
        )
      ),
      .groups = "drop"
    ) %>%
    mutate(
      tidy = map(test, broom::tidy)
    ) %>%
    unnest(tidy)

  sig_metrics <- cor_results %>%
    group_by(PhysChem_Metric) %>%
    summarise(
      any_sig = any(
        p.value < 0.05,
        na.rm = TRUE
      )
    ) %>%
    filter(any_sig) %>%
    pull(PhysChem_Metric)

  long_df <- long_df %>%
    filter(
      PhysChem_Metric %in% sig_metrics
    )

  ggplot(
    long_df,
    aes(
      x = Value_log,
      y = RiskScore,
      color = Group.x
    )
  ) +

    geom_point(
      alpha = 0.7,
      size = 3
    ) +

    geom_smooth(
      method = "lm",
      se = TRUE,
      linetype = "dashed",
      linewidth = 0.8
    ) +

    stat_cor(
      aes(
        label = paste(
          ..r.label..,
          ..p.label..,
          sep = "~`,`~"
        ),
        group = Group.x
      ),

      method = "pearson",

      size = 4,

      color = "black",

      label.x.npc = "left",
      label.y.npc = "top"
    ) +

    facet_wrap(
      ~PhysChem_Metric,
      scales = "free_x"
    ) +

    scale_color_manual(
      values = met_cols$stage_lower
    ) +

    scale_y_continuous(
      limits = y_limits
    ) +

    labs(
      x = NULL,
      y = score_type
    ) +

    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

p_err <- make_sig_physchem_plot(
  full_merged_dataset,
  "ERR",
  c(0, 60)
)

p_hhrr <- make_sig_physchem_plot(
  full_merged_dataset,
  "HHRR",
  c(0, 5)
)

fig4 <- p_err / p_hhrr +
  plot_annotation(tag_levels = "A")

save_fig(
  fig4,
  "METFIG4.png",
  width = 14,
  height = 18
)

# Figure 5 -- ESKAPEE Burden

esk <- read_tsv(
  "ESKAPEE_combined.tsv",
  col_names = FALSE,
  show_col_types = FALSE
)

meta <- read_excel(
  "full_merged_dataset.xlsx"
) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )

esk_clean <- esk %>%
  filter(X4 == "species") %>%
  mutate(
    species_group = case_when(

      X5 == "Escherichia coli" ~
        "Escherichia coli",

      X5 == "Enterococcus faecium" ~
        "Enterococcus faecium",

      X5 == "Staphylococcus aureus" ~
        "Staphylococcus aureus",

      X5 == "Klebsiella pneumoniae" ~
        "Klebsiella pneumoniae",

      X5 == "Acinetobacter baumannii" ~
        "Acinetobacter baumannii",

      X5 == "Pseudomonas aeruginosa" ~
        "Pseudomonas aeruginosa",

      str_detect(X5, "^Enterobacter") ~
        "Enterobacter spp.",

      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(species_group))

esk_species <- esk_clean %>%
  group_by(X1, species_group) %>%
  summarise(
    n_ESKAPEE = n(),
    .groups = "drop"
  ) %>%
  rename(sample_id = X1)

plot_data5 <- esk_species %>%
  left_join(
    meta %>%
      select(
        sample_id,
        Site,
        `Sample Round`,
        Group.x,
        nContigs
      ),
    by = "sample_id"
  ) %>%
  mutate(
    ESKAPEE_norm = n_ESKAPEE / nContigs,

    Group = case_when(
      tolower(Group.x) == "influent" ~
        "Influent",
      TRUE ~ "Effluent"
    )
  )

paired_data <- plot_data5 %>%
  group_by(
    Site,
    `Sample Round`,
    species_group
  ) %>%
  filter(
    n_distinct(Group) == 2
  ) %>%
  ungroup()

stat_results <- paired_data %>%
  group_by(species_group) %>%
  wilcox_test(
    ESKAPEE_norm ~ Group,
    paired = TRUE
  ) %>%
  add_significance()

stat_results <- stat_results %>%
  left_join(
    paired_data %>%
      group_by(species_group) %>%
      summarise(
        y.position =
          max(ESKAPEE_norm) * 1.2
      ),
    by = "species_group"
  )

p5 <- ggplot(
  paired_data,
  aes(
    x = Group,
    y = ESKAPEE_norm + 1e-8,
    group = interaction(
      Site,
      `Sample Round`
    )
  )
) +

  geom_line(
    color = met_cols$neutral["mid"],
    linewidth = 0.6,
    alpha = 0.7
  ) +

  geom_point(
    aes(color = Group),
    size = 2.8,
    alpha = 0.9,
    position = position_jitter(
      width = 0.05
    )
  ) +

  scale_color_manual(
    values = met_cols$stage
  ) +

  scale_y_log10() +

  facet_wrap(
    ~ species_group,
    nrow = 1,
    scales = "free_y"
  ) +

  stat_pvalue_manual(
    stat_results,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.4
  ) +

  labs(
    x = NULL,
    y = "Normalised ESKAPEE burden (log10)"
  ) +

  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

save_fig(
  p5,
  "METFIG5.png",
  width = 14,
  height = 3.8
)

# COMPLETE

cat("  METFIG1.png\n")
cat("  METFIG2.png\n")
cat("  METFIG3.png\n")
cat("  METFIG4.png\n")
cat("  METFIG5.png\n")

