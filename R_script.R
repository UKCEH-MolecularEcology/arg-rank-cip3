library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(performance)

theme_set(
  theme_bw(base_family = "Helvetica") +
    theme(
      strip.background = element_rect(
        colour = "black",
        fill = "white"
      ),
      strip.text = element_text(
        size = 14,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 90,
        size = 9,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = 11,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 13,
        face = "bold"
      ),
      axis.title.x = element_text(
        size = 13,
        face = "bold"
      ),
      legend.text = element_text(
        size = 11
      ),
      legend.title = element_text(
        size = 11,
        face = "bold"
      ),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      )
    )
)

group_cols <- c(
  "Influent" = "#D55E00",
  "Effluent" = "#0072B2"
)

secondary_cols <- c(
  "AS" = "#009E73",
  "TF" = "#E69F00",
  "BAFF" = "#56B4E9",
  "AS_TF" = "#CC79A7"
)

df_raw <- read_excel(
  "full_merged_dataset.xlsx"
)

colnames(df_raw) <- make.names(
  trimws(colnames(df_raw))
)

df <- df_raw %>%
  select(
    Site,
    SampleRound = Sample.Round,
    ScoreType,
    Group = Group.x,
    RiskScore
  ) %>%
  mutate(

    Group = tolower(Group),

    Group = ifelse(
      grepl("influent", Group),
      "Influent",
      "Effluent"
    ),

    ScoreType = recode(
      tools::toTitleCase(trimws(ScoreType)),
      "Ecological" = "ERR",
      "Human Health" = "HHRR"
    )

  ) %>%
  filter(
    ScoreType %in% c(
      "ERR",
      "HHRR"
    )
  )

secondary_mapping <- c(

  "Driffield" = "AS",
  "Exmouth" = "BAFF",
  "Finham" = "AS",
  "Pembury" = "TF",
  "Pen-y-bont" = "AS",
  "Saltford" = "TF",
  "Seaton Carew" = "AS",
  "Stockport" = "AS",
  "Chalton" = "AS_TF"

)

df$Secondary <- secondary_mapping[df$Site]

df <- df %>%
  filter(
    !is.na(Secondary)
  )

df$Group <- factor(
  df$Group,
  levels = c(
    "Influent",
    "Effluent"
  )
)

df$Secondary <- factor(
  df$Secondary,
  levels = c(
    "AS",
    "TF",
    "BAFF",
    "AS_TF"
  )
)

df$ScoreType <- factor(
  df$ScoreType,
  levels = c(
    "ERR",
    "HHRR"
  )
)

df_wide <- df %>%
  pivot_wider(
    id_cols = c(
      Site,
      Secondary,
      ScoreType,
      SampleRound
    ),
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

model <- lmer(
  RiskScore ~ 
    Group * Secondary +
    ScoreType +
    (1 | Site) +
    (1 | Site:SampleRound),
  data = df
)

print(summary(model))

model_anova <- anova(
  model
)

print(model_anova)

write.csv(
  as.data.frame(model_anova),
  "MET_Model_ANOVA.csv",
  row.names = FALSE
)

check_model(model)

model_performance <- performance::model_performance(
  model
)

write.csv(
  as.data.frame(model_performance),
  "MET_Model_Performance.csv",
  row.names = FALSE
)

emm <- emmeans(
  model,
  ~ Group | Secondary * ScoreType
)

emm_df <- as.data.frame(
  emm
)

write.csv(
  emm_df,
  "MET_Emmeans_Group_By_System.csv",
  row.names = FALSE
)

contrast_df <- as.data.frame(
  contrast(
    emm,
    method = "revpairwise"
  )
)

contrast_df$p.adj <- p.adjust(
  contrast_df$p.value,
  method = "BH"
)

write.csv(
  contrast_df,
  "MET_Treatment_Contrasts.csv",
  row.names = FALSE
)

order_levels <- df_wide %>%
  group_by(Secondary) %>%
  summarise(
    med = median(
      Delta,
      na.rm = TRUE
    )
  ) %>%
  arrange(med) %>%
  pull(Secondary)

df$Secondary <- factor(
  df$Secondary,
  levels = order_levels
)

df_wide$Secondary <- factor(
  df_wide$Secondary,
  levels = order_levels
)

emm_df$Secondary <- factor(
  emm_df$Secondary,
  levels = order_levels
)

p1 <- ggplot(
  df,
  aes(
    x = Secondary,
    y = RiskScore,
    fill = Group
  )
) +

  geom_boxplot(
    position = position_dodge(0.75),
    width = 0.55,
    linewidth = 0.4,
    outlier.shape = NA,
    alpha = 0.75
  ) +

  geom_jitter(
    aes(
      colour = Group
    ),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 1.3,
    alpha = 0.55
  ) +

  facet_wrap(
    ~ScoreType,
    scales = "free_y"
  ) +

  scale_fill_manual(
    values = group_cols,
    name = "Sample type"
  ) +

  scale_colour_manual(
    values = group_cols,
    name = "Sample type"
  ) +

  labs(
    x = NULL,
    y = "MetaCompare risk score"
  )

p2 <- ggplot(
  emm_df,
  aes(
    x = Secondary,
    y = emmean,
    colour = Group,
    group = Group
  )
) +

  geom_point(
    position = position_dodge(0.35),
    size = 2.8
  ) +

  geom_line(
    position = position_dodge(0.35),
    linewidth = 0.6
  ) +

  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.15,
    linewidth = 0.5,
    position = position_dodge(0.35)
  ) +

  facet_wrap(
    ~ScoreType,
    scales = "free_y"
  ) +

  scale_colour_manual(
    values = group_cols,
    name = "Sample type"
  ) +

  labs(
    x = NULL,
    y = "Estimated risk score"
  )

p3 <- ggplot(
  df_wide,
  aes(
    x = Secondary,
    y = Delta,
    fill = Secondary
  )
) +

  geom_violin(
    trim = FALSE,
    alpha = 0.5
  ) +

  geom_boxplot(
    width = 0.18,
    linewidth = 0.4,
    fill = "white",
    outlier.shape = NA
  ) +

  geom_jitter(
    width = 0.1,
    size = 1.2,
    alpha = 0.7
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey40"
  ) +

  facet_wrap(
    ~ScoreType,
    scales = "free_y"
  ) +

  scale_fill_manual(
    values = secondary_cols,
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = "Change in risk score\n(Effluent - Influent)"
  )

final_fig <- 
  (p1 | p2) /
  p3 +
  plot_annotation(
    tag_levels = "A",
    title =
      "Treatment-associated reduction in ARG risk across wastewater treatment systems"
  ) &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold"
    ),
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  "MET-FIG-Treatment-Risk-Reduction.png",
  final_fig,
  width = 11,
  height = 8,
  dpi = 600
)

paired_summary <- df_wide %>%
  group_by(
    ScoreType,
    Secondary
  ) %>%
  summarise(
    n = n(),
    median_delta = median(
      Delta,
      na.rm = TRUE
    ),
    IQR_delta = IQR(
      Delta,
      na.rm = TRUE
    ),
    mean_delta = mean(
      Delta,
      na.rm = TRUE
    ),
    sd_delta = sd(
      Delta,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  paired_summary,
  "MET_Paired_Risk_Reduction_By_System.csv",
  row.names = FALSE
)
