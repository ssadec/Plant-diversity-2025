


# =========================
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggsci)
library(stringr)

# =========================

df <- read.csv("/Users/yusha/Fig1c_number_of_species.csv", stringsAsFactors = FALSE)

# =========================

df_clean <- df %>%
  mutate(
    Zone = str_to_lower(trimws(Zone)),
    Zone = case_when(
      Zone %in% c("tropic","tropical") ~ "Tropical zone",
      Zone == "temperate"              ~ "Temperate zone",
      TRUE                             ~ NA_character_
    ),
    CP_Ct = case_when(
      CP_Ct == "Food_crop" ~ "Food crop",
      CP_Ct == "Cash_crop" ~ "Cash crop",
      TRUE ~ NA_character_
    ),
    CP_C3 = case_when(
      CP_C3 == "Herbaceous_plant" ~ "Herbaceous crop",
      CP_C3 == "Woody_plant"      ~ "Woody crop",
      TRUE ~ NA_character_
    )
  )

# =========================


# All studies：
df_all <- df_clean %>% mutate(Group = "All studies")

# Zone：
df_zone <- df_clean %>% filter(!is.na(Zone)) %>% rename(Group = Zone)

# CP_Ct：
df_ct <- df_clean %>% filter(!is.na(CP_Ct)) %>% rename(Group = CP_Ct)

# CP_C3：
df_c3 <- df_clean %>% filter(!is.na(CP_C3)) %>% rename(Group = CP_C3)

# =========================
df_plot <- bind_rows(df_all, df_zone, df_ct, df_c3)


# =========================
df_plot <- bind_rows(df_all, df_zone, df_ct, df_c3)

# =========================
n_all  <- nrow(df_all)  # All studies 
n_zone <- df_zone %>% group_by(Group) %>% summarise(N = n_distinct(S_ID), .groups="drop")
n_ct   <- df_ct   %>% group_by(Group) %>% summarise(N = n_distinct(S_ID), .groups="drop")
n_c3   <- df_c3   %>% group_by(Group) %>% summarise(N = n_distinct(S_ID), .groups="drop")

# =========================
n_table <- bind_rows(
  tibble(Group = "All studies", N = n_all),
  n_zone, n_ct, n_c3
)

# =========================
df_plot <- df_plot %>%
  left_join(n_table, by="Group") %>%
  mutate(x_label = paste0(Group, "\n(N=", N, ")"))

# =========================
df_plot$x_label <- factor(
  df_plot$x_label,
  levels = c(
    paste0("All studies\n(N=", n_table$N[n_table$Group=="All studies"], ")"),
    paste0("Herbaceous crop\n(N=", n_table$N[n_table$Group=="Herbaceous crop"], ")"),
    paste0("Woody crop\n(N=", n_table$N[n_table$Group=="Woody crop"], ")"),
    paste0("Food crop\n(N=", n_table$N[n_table$Group=="Food crop"], ")"),
    paste0("Cash crop\n(N=", n_table$N[n_table$Group=="Cash crop"], ")"),
    paste0("Temperate zone\n(N=", n_table$N[n_table$Group=="Temperate zone"], ")"),
    paste0("Tropical zone\n(N=", n_table$N[n_table$Group=="Tropical zone"], ")")
  )
)



# =========================

summary_df <- df_plot %>%
  group_by(x_label) %>%
  summarise(mean_val = mean(NCP_N, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(df_plot, aes(x = x_label, y = NCP_N, fill = Group, color = Group)) +

  geom_half_violin(
    color = NA,
    side = "r",
    trim = FALSE,
    adjust = 1.2,
    width = 0.8,
    position = position_nudge(x = 0.12)
  ) +
  
  
  geom_boxplot(
    color = "black",
    width = 0.15,
    size = 0.6,
    outlier.size = 1,
    position = position_nudge(x = -0.05)
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "yellow",
    color = "black",
    position = position_nudge(x = -0.05)
  ) 
  

  geom_point(
    data = df_plot %>% distinct(S_ID, x_label, NCP_N, Group),
    aes(x = as.numeric(x_label) - 0.15),
    size = 2,
    position = position_jitter(width = 0.05)
  ) +
  

  scale_fill_jama() +
  scale_color_jama() +

  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 2, 5, 10, 20, 40),
    limits = c(0, 40),
    expand = expansion(mult = c(0, 0.001))
  ) +
  
  labs(y = "Number of added plant species") +
  
  theme_test(base_family = "Arial") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, family = "Arial", color = "black"),
    axis.text.x  = element_text(size = 24, hjust = 0.5, vjust = 0.5, family = "Arial", color = "black"),
    axis.text.y  = element_text(size = 24, family = "Arial", color = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10),
    aspect.ratio = 1/3
  ) +
  coord_cartesian(clip = "off")

p2

