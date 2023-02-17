## Summary plots of the data


### load packages
source("source_packages.R")


###########################
### Plot sequence reads ###
###########################

##----------------------------------------##
## histograms of number of sequence reads ##
##----------------------------------------##


## load insect COI data
coi_flies <- read.csv("./data/flies_COI.csv")
coi_mosq <- read.csv("./data/nbfmosq_COI.csv")

head(coi_flies)
head(coi_mosq)

unique(coi_flies$plot)
coi_flies[coi_flies$locality == "Spandau_1",]

## load and prepare sample and locality names
sample_data <- read.csv("./data/sample_summary.csv")

sample_data

sample_flies <- sample_data %>% 
  filter(type == "Flies") %>% 
  dplyr::select(-c(id_BIBS, date))

sample_nbf <- sample_data %>% 
  filter(type == "NBF") 

glimpse(sample_flies)
glimpse(sample_nbf)


#######################
## Fly COI histogram 

## add sample data to insect COI tables

unique(coi_flies$plot)
unique(sample_flies$plot)

coi_flies2 <- coi_flies %>% 
  mutate(plot = case_when(
    plot == "OH-7" ~ "OH-07",
    TRUE ~ plot
  )) %>% 
  mutate(fly_det = case_when(
    fly_sp == "No_match" ~ "No",
    TRUE ~ "Yes")) %>% 
  dplyr::select(-c(library, database, count, percentage))

flies_tmp <- coi_flies2 %>% 
  mutate(sample = as.character(sample)) %>% 
  full_join(sample_flies, by = c("sample", "plot", "locality")) %>% 
  arrange(sample) %>% 
  mutate(type2 = "FlyPool") %>%
  mutate(sample_num = as.integer(as.factor(sample))) %>% 
  mutate(sample_name = paste(type2, sample_num, locality, sep = "_")) %>% 
  mutate(short_loc = case_when(
    locality == "Grunewald_1" ~ "_Gw1",
    locality == "Grunewald_2" ~ "_Gw2",
    locality == "Mueggelsee" ~ "_Ms",
    locality == "Spandau_1" ~ "_Sp1"
  )) %>% 
  mutate(short_name = paste0("FlyP", sample_num, short_loc)) %>% 
  group_by(sample_name, short_name, sample, fly_det, sample_num, COI_results) %>% 
  summarise(n = length(unique(fly_sp))) %>%
  ungroup() %>% 
  arrange(sample_num) %>% 
  mutate(primer = "COI_Fly") %>% 
  mutate(n2 = case_when(
    fly_det == "No" ~ as.numeric(0), 
    fly_det == "Yes" ~ as.numeric(n)
  ))


# reordering to plot based on location order
my_order <- unique(flies_tmp$sample_name)
my_order2 <- unique(flies_tmp$short_name)

flies_toplot <- flies_tmp %>% 
  group_by(sample_name, short_name, sample, sample_num, primer, COI_results) %>% 
  summarise(n_taxa = sum(n2)) %>% 
  mutate(sample_name2 = as.factor(sample_name)) %>% 
  mutate(sample_name2 = fct_relevel(sample_name2, my_order)) %>% 
  mutate(short_name = as.factor(short_name)) %>% 
  mutate(short_name = fct_relevel(short_name, my_order2)) 


# create vector to colour x axis labels 
flies_toplot <- arrange(flies_toplot, flies_toplot$short_name)
flies_toplot <- flies_toplot %>% 
  mutate(colours = case_when(
    is.na(n_taxa) ~ "grey50",
    TRUE ~ "black")) 

(fly_hist <- ggplot(flies_toplot, aes(x = as.factor(sample_name2), y = n_taxa)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
           alpha = 0.7, fill = "forestgreen") +
  xlab("Sample id") +
  ylab("N OTUs detected") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, colour = flies_toplot$colours),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid = element_blank()))



#######################
## NBF COI histogram 

## summary data
mosq_summary <- coi_mosq %>% 
  group_by(sample) %>% 
  summarise(n_total = n(),
            n_taxat = length(unique(mosq_sp)),
            n_nomatch = sum(mosq_sp == "No_match")) %>% 
  mutate(positive_ident = n_total - n_nomatch)
mosq_summary

# add sample data to insect COI tables
unique(coi_mosq$plot)
unique(sample_nbf$plot)

mosq_samples_filter <- unique(sample_nbf$sample)

coi_mosq2 <- coi_mosq %>% 
  filter(sample %in% mosq_samples_filter) %>% 
  mutate(mosq_det = case_when(
    mosq_sp == "No_match" ~ "No",
    TRUE ~ "Yes")) %>% 
  dplyr::select(-c(library, plot, count, percentage))

mosq_tmp <- coi_mosq2 %>% 
  mutate(sample = as.character(sample)) %>% 
  full_join(sample_nbf, by = c("sample")) %>% 
  ## change name of sample for plotting
  mutate(sample2 = case_when(
    sample == "M160" ~ "M147.5",
    TRUE ~ sample
  )) %>% 
  arrange(sample2) %>% 
  mutate(type2 = "NBFPool") %>%
  mutate(sample_num = as.integer(as.factor(sample2))) %>% 
  mutate(sample_name = paste(type2, sample_num, locality, sep = "_")) %>% 
  mutate(short_loc = case_when(
    locality == "Grunewald_1" ~ "_Gw1",
    locality == "Grunewald_2" ~ "_Gw2",
    locality == "Mueggelsee" ~ "_Ms",
    locality == "Spandau_1" ~ "_Sp1"
  )) %>% 
  mutate(short_name = paste0("NBFP", sample_num, short_loc)) %>% 
  group_by(sample_name, short_name, sample2, mosq_det, sample_num, COI_results) %>% 
  summarise(n = length(unique(mosq_sp))) %>%
  ungroup() %>% 
  arrange(sample_num) %>% 
  mutate(primer = "COI_NBF") %>% 
  mutate(n2 = case_when(
    mosq_det == "No" ~ as.numeric(0), 
    mosq_det == "Yes" ~ as.numeric(n)
  ))

# reordering to plot based on location order
my_order <- unique(mosq_tmp$sample_name)
my_order2 <- unique(mosq_tmp$short_name)

mosq_toplot <- mosq_tmp %>% 
  group_by(sample_name, short_name, sample2, sample_num, primer, COI_results) %>% 
  summarise(n_taxa = sum(n2)) %>% 
  mutate(sample_name2 = as.factor(sample_name)) %>% 
  mutate(sample_name2 = fct_relevel(sample_name2, my_order)) %>% 
  mutate(short_name = as.factor(short_name)) %>% 
  mutate(short_name = fct_relevel(short_name, my_order2)) 

# create vector to colour x axis labels 
mosq_toplot <- arrange(mosq_toplot, mosq_toplot$short_name)
mosq_toplot <- mosq_toplot %>% 
  mutate(colours = case_when(
    is.na(n_taxa) ~ "grey50",
    TRUE ~ "black")) 
  
                        
mosq_hist <- ggplot(mosq_toplot, aes(x = as.factor(sample_name2), y = n_taxa)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
           alpha = 0.7, fill = "darkorange") +
  xlab("Sample id") +
  ylab("N OTUs detected") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, colour = mosq_toplot$colours),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid = element_blank())

mosq_hist



##---------------------------------------##
## Histogram Mammal identification reads ##
##---------------------------------------##


## plot mammal identification reads

mammal_reads <- read.csv("./data/mammal_reads.csv")
sample_mammals <- read.csv("./data/samples_mammals_data.csv")


head(mammal_reads)
sample_mammals

## remove not mammals
unique(mammal_reads$species)

## add sample data
mammals2 <- mammal_reads %>% 
  full_join(sample_mammals, by = "sample") %>% 
  left_join(sample_flies %>% 
              dplyr::select(sample, locality), by = "sample") %>% 
  left_join(sample_nbf %>% 
              dplyr::select(sample, locality), by = "sample") %>% 
  mutate(locality = case_when(
    is.na(locality.x) ~ locality.y,
    is.na(locality.y) ~ locality.x)) %>% 
  dplyr::select(-c(locality.x, locality.y)) %>% 
  filter(!is.na(locality)) %>%  # remore sample with error
  arrange(sample) %>% 
  mutate(mam_det = case_when(
    is.na(species) ~ "No",
    TRUE ~ "Yes")) 

## prepare data based on insect, location and marker
mammals_tmp <- mammals2 %>% 
  mutate(type2 = case_when(
    grepl("M", sample) == TRUE ~ "MosqPool",
    TRUE ~ "FlyPool"
  )) %>% 
  mutate(sample2 = case_when(
    sample == "M160" ~ "M147.5", # this is needed because this sample appears in the wrong order otherwise
    TRUE ~ sample
  )) %>% 
  arrange(sample2) %>% 
  mutate(sample_num = as.integer(as.factor(sample2))) %>% 
  mutate(sample_name = paste(type2, sample_num, locality, sep = "_")) %>% 
  mutate(short_loc = case_when(
    locality == "Grunewald_1" ~ "_Gw1",
    locality == "Grunewald_2" ~ "_Gw2",
    locality == "Mueggelsee" ~ "_Ms",
    locality == "Spandau_1" ~ "_Sp1"
  )) %>% 
  mutate(short_type = case_when(
    type2 == "FlyPool" ~ "FlyP",
    type2 == "MosqPool" ~ "NBFP",
  )) %>% 
  mutate(short_name = paste0(short_type, sample_num, short_loc)) %>% 
  group_by(sample_name, short_name, sample, mam_det, sample_num, primer) %>% 
  summarise(n = length(unique(species))) %>%
  ungroup() %>% 
  arrange(sample_num) %>% 
  mutate(n2 = case_when(
    mam_det == "No" ~ as.numeric(0), 
    mam_det == "Yes" ~ as.numeric(n)
  ))

# reordering to plot based on location order
my_order <- unique(mammals_tmp$sample_name)
my_order2 <- unique(mammals_tmp$short_name)

mammals_toplot <- mammals_tmp %>% 
  group_by(sample_name, short_name, sample, sample_num, primer) %>% 
  summarise(n_taxa = sum(n2)) %>% 
  mutate(sample_name2 = as.factor(sample_name)) %>% 
  mutate(sample_name2 = fct_relevel(sample_name2, my_order)) %>% 
  mutate(short_name = as.factor(short_name)) %>% 
  mutate(short_name = fct_relevel(short_name, my_order2)) 

mammals_hist <- ggplot(mammals_toplot, aes(x = as.factor(sample_name2), y = n_taxa, fill = primer)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
           alpha = 0.7) +
  xlab("Sample id") +
  ylab("N OTUs detected") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid = element_blank())

mammals_hist


##---------------##
## Combined plot ##
##---------------##

flies_toplot <- flies_toplot %>% mutate(type = "COI_Flies")
mosq_toplot <- mosq_toplot %>% mutate(type = "COI_NBF")
mammals_toplot <- mammals_toplot %>% mutate(type = "Mammal identification")

all_taxa <- rbind(mammals_toplot, flies_toplot, mosq_toplot)


all_taxa <- all_taxa %>% 
  mutate(legend = as.factor(as.character(primer))) %>% 
  mutate(legend = fct_relevel(legend, "COI_Fly", "COI_NBF", "12S", "16S", "NA"))

# create all plots individually
(fly_hist <- ggplot(flies_toplot, aes(x = as.factor(short_name), y = n_taxa)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
             alpha = 0.7, fill = "blue") +
    xlab("") +
    ylab("N OTUs detected") +
    ylim(0,8) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, colour = flies_toplot$colours),
          axis.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.grid = element_blank())
)


(mosq_hist <- ggplot(mosq_toplot, aes(x = as.factor(short_name), y = n_taxa)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
             alpha = 0.7, fill = "darkred") +
    xlab("") +
    ylab("N OTUs detected") +
    ylim(0,8) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, colour = mosq_toplot$colours),
          axis.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.grid = element_blank())
)


mammals_toplot <- mammals_toplot %>% 
  arrange(short_name) %>% 
  group_by(short_name) %>% 
  slice_head() %>% 
  ungroup() %>% 
  arrange(short_name) %>% 
  mutate(colour = case_when(
    is.na(primer) ~ "grey50",
    TRUE ~ "black"))

(mammals_hist <- ggplot(mammals_toplot, aes(x = as.factor(short_name), y = n_taxa, fill = primer)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
             alpha = 0.7, show.legend = FALSE) +
    xlab("Sample id") +
    ylab("N OTUs detected") +
    scale_fill_manual(values = c("forestgreen", "darkorange"), na.value = "white", 
                      name = "Primer", labels = c("12S", "16S", "")) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, colour = mammals_toplot$colour),
          axis.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.grid = element_blank())
)


# create a plot with all data to extract the legend for later
## plot histograms with only fly and mosquito pools divided by primer
(both_hist <- ggplot(all_taxa, aes(x = sample_name2, y = n_taxa, fill = legend)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
             alpha = 0.7, width = 1
    ) +
    scale_fill_manual(values = c("blue", "darkred", "forestgreen", "darkorange"),
                      na.value = "white", 
                      name = "Marker", 
                      labels = c( "COI-Fly", "COI-NBF", "12S", "16S", "")) +
    xlab("Sample id") +
    ylab("Unique Taxa (OTUs)") +
    facet_wrap(~type, 
               ncol = 2,
               nrow = 2
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_text(size = 14, colour = "black", face = "bold"), 
      strip.text = element_text(size = 12, colour = "black", face = "bold"),
      legend.title = element_text(size = 14, colour = "black", face = "bold"),
      legend.text = element_text(size = 14, colour = "black")
    )
)

#extract legend
my_legend <- get_legend(both_hist)


# combine the plots together
plots <- align_plots(mammals_hist, fly_hist, align = "v", axis = "l")
top_row <- plot_grid(plots[[2]], mosq_hist, labels = c("A", "B"), align = "h", label_size = 12)
allplots <- plot_grid(top_row, plots[[1]], labels = c("", "C"), label_size = 12, ncol = 1)

final_plot <- plot_grid(allplots, my_legend, ncol = 2, 
                        rel_widths = c(1.5,0.2)) +
  theme(plot.background = element_rect(fill = "white"))

ggsave(plot = final_plot, 
       "./output/Fig2_hist_uniqueOTUs_COI_marker.png",  
       width = 12, height = 10, dpi = 600)




################################################################################
################################################################################


###################################
## Plot heatmap mammal detection ##
###################################


## fly data
fly_dat <- read.csv("./data/fly_mammal_diversity.csv")
head(fly_dat)
str(fly_dat)

unique(fly_dat$Locality)

# prepare clean data in long format for heatmap
fly_dat_2 <- fly_dat %>% 
  dplyr::select(Locality, Marker, Canis_lupus:Rattus_norvegicus) %>% 
  mutate(Locality_short = case_when(
    Locality == "Grunewald-1" ~ "Grunewald_1",
    Locality == "Grunewald-2" ~ "Grunewald_2",
    Locality == "Muggelsee" ~ "Mueggelsee",
    Locality == "Spandau" ~ "Spandau_1")) %>% 
  mutate(Marker_new = case_when(
    Marker == "16S" ~ "Flies-16S",
    Marker == "12S" ~ "Flies-12S")) %>% 
  mutate(LocMarker = paste0(Marker_new, "-", Locality_short), 
         iDNA = "flies") %>% 
  group_by(LocMarker) %>% 
  summarise(across(where(is.numeric), ~ sum(.x))) %>% 
  mutate(Canis_lupus_domesticus = Canis_lupus) %>% 
  dplyr::select(-Canis_lupus)

fly_dat_2  


## varibles: locality, marker 
data <- fly_dat_2 %>% 
  column_to_rownames("LocMarker") %>%
  dplyr::select_if(., is.numeric) %>% 
  mutate_all(~replace(., . == 0, NA)) %>% 
  as.matrix()
summary(data)


## BF mosquito data
bf_tmp <- read.csv("./data/bfmosq_mammalID.csv")

bf_tmp <- bf_tmp %>% 
  filter(Kiez != "BF") # remove extra row

### prepare 12s data
tmp_12s <- bf_tmp %>% 
  dplyr::select(DNA_no, Locality, Kiez, X12S_340_sanger) %>% 
  mutate(n_12s = case_when(
    !is.na(X12S_340_sanger) ~ 1, 
    TRUE ~ 0)) %>% 
  pivot_wider(names_from = X12S_340_sanger, values_from = n_12s) %>% 
  dplyr::select(-c(`NA`, `Turdus philomelos`)) %>% 
  clean_names(case = "none") %>% 
  mutate(Marker = "12S")

tmp_12s

### prepare 16s data
tmp_16s <- bf_tmp %>% 
  dplyr::select(DNA_no, Locality, Kiez, X16S) %>% 
  mutate(n_16s = case_when(
    !is.na(X16S) ~ 1, 
    TRUE ~ 0)) %>% 
  pivot_wider(names_from = X16S, values_from = n_16s) %>% 
  dplyr::select(-c(`NA`)) %>% 
  clean_names(case = "none") %>% 
  mutate(Vulpes_vulpes = NA, # add empty column to merge with previous table
         Dama_dama = NA,
         Marker = "16S")

tmp_16s

### put both markers together
bf_dat <- rbind(tmp_12s, tmp_16s)

unique(bf_dat$Marker)

# prepare data for models, group by locality
bf_dat_2 <- bf_dat %>% 
  dplyr::select(Kiez, Marker, Capreolus_capreolus:Dama_dama) %>% 
  group_by(Kiez, Marker) %>% 
  summarise(Capreolus_capreolus = sum(Capreolus_capreolus, na.rm = TRUE),
            Dama_dama = sum(Dama_dama, na.rm = TRUE), 
            Sus_scrofa = sum(Sus_scrofa, na.rm = TRUE), 
            Vulpes_vulpes = sum(Vulpes_vulpes, na.rm = TRUE)) %>% 
  as.data.frame() %>% 
  mutate(Locality_short = case_when(
    Kiez == "Grunewald_1" ~ "Grunewald_1",
    Kiez == "Grunewald_2" ~ "Grunewald_2",
    Kiez == "Muggelsee" ~ "Mueggelsee",
    Kiez == "Spandau_1" ~ "Spandau_1",
    Kiez == "Spandau_2" ~ "Spandau_2")) %>% 
  mutate(Marker_new = case_when(
    Marker == "16S" ~ "BF-16S",
    Marker == "12S" ~ "BF-12S")) %>% 
  mutate(LocMarker = paste0(Marker_new, "-", Locality_short), 
         iDNA = "BF mosquitos") %>% 
  group_by(LocMarker) %>% 
  summarise(across(where(is.numeric), ~ sum(.x))) 

bf_dat_2

## columns to mosquitoes that were in flies
toaddcols <- colnames(fly_dat_2)[!colnames(fly_dat_2) %in% colnames(bf_dat_2) ]
toaddcols <- toaddcols %>% 
  janitor::make_clean_names(case = "none")

bf_dat_longer <- bf_dat_2 %>% 
  add_column(!!!toaddcols) %>% 
  janitor::clean_names(case = "none") %>% 
  mutate(across(all_of(toaddcols), ~ 0)) %>% 
  pivot_longer(cols = -LocMarker, names_to = "Species")


## add columns to flies that were in moquitoes
toaddcols_2 <- colnames(bf_dat_2)[!colnames(bf_dat_2) %in% colnames(fly_dat_2) ]

fly_dat_longer <- fly_dat_2 %>% 
  add_column(!!!toaddcols_2) %>% 
  janitor::clean_names(case = "none") %>% 
  mutate(across(all_of(toaddcols_2), ~ 0)) %>% 
  pivot_longer(cols = -LocMarker, names_to = "Species")

## put together flies and mosquitoes
all_data <- rbind(fly_dat_longer, bf_dat_longer)


##################
## Plot all with common names
unique(all_data$Species)

all_data_names <- all_data %>% 
  mutate(Species_english = case_when(
    Species == "Myodes_glareolus" ~ "Bank vole",
    Species == "Sus_scrofa" ~ "Wild boar",
    Species == "Vulpes_vulpes" ~ "Red fox",
    Species == "Apodemus_sylvaticus" ~ "Wood mouse",
    Species == "Meles_meles" ~ "European badger",
    Species == "Siciurus_vulgaris" ~ "Red squirrel",
    Species == "Capreolus_capreolus" ~ "Roe deer",
    Species == "Oryctolagus_cuniculus" ~ "European rabbit",
    Species == "Procyon_lotor" ~ "Raccoon",
    Species == "Rattus_norvegicus" ~ "Brown rat",
    Species == "Canis_lupus_domesticus" ~ "Domestic dog",
    Species == "Dama_dama" ~ "Fallow deer"))


(heat_flies_2 <- ggplot(all_data_names, aes(Species_english, LocMarker, fill= as.factor(value))) + 
    geom_tile() +
    scale_fill_scico_d(palette = "bilbao", name = "N detections") +
    xlab("Mammal species") +
    ylab("Sample type (primer-insect-locality") +
    theme(panel.border = element_rect(colour = "black", fill = "transparent"), 
          axis.text.x = element_text(size = 12, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black", face = "bold")) +
    geom_hline(yintercept = 14.5, colour = "grey30", lwd = 0.2) +
    geom_hline(yintercept = 10.5, colour = "grey30",lwd = 0.5) +
    geom_hline(yintercept = 5.5, colour = "grey30", lwd = 0.2))

ggsave(plot = heat_flies_2, 
       filename = "./output/Fig3_heatmap_primer_insect_locality_english.png", 
       dpi = 600, width = 8, height = 8)




################################################################################
################################################################################


#########################
## Accumulation curves ##
#########################


fly_dat <- read.csv("./Data/fly_mammal_diversity.csv")
head(fly_dat)
str(fly_dat)

fly_accum <- fly_dat %>% 
  dplyr::select(Locality, Marker, Canis_lupus:Rattus_norvegicus) %>% 
  group_by(Locality) %>% 
  mutate(id = 1:n()) %>% 
  ungroup() %>% 
  mutate(for_rows = paste0(Locality, id)) %>% 
  as.data.frame()

mammals_matrix <-fly_accum %>% 
  column_to_rownames("for_rows") %>% 
  dplyr::select(-Locality, -Marker, -id)

env_matrix <- fly_accum %>% 
  column_to_rownames("for_rows") %>% 
  dplyr::select(Marker, Locality)

Accum.1 <- BiodiversityR::accumcomp(mammals_matrix, y = env_matrix, factor = "Marker", method = "exact", conditioned = FALSE, plotit = FALSE)
accum.long1 <- BiodiversityR::accumcomp.long(Accum.1, ci=NA, label.freq=5)
head(accum.long1)

plotgg1 <- ggplot(data = accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  geom_line(aes(colour = Grouping), size=2) +
  geom_point(data = subset(accum.long1, labelit == TRUE), 
             aes(colour = Grouping, shape = Grouping), size=5) +
  geom_ribbon(aes(fill = Grouping), colour = "transparent", alpha=0.2, show.legend=FALSE) + 
  scale_colour_manual(values = c("grey20", "orange")) +
  scale_fill_manual(values = c("grey20", "orange")) +
  labs(x = "# Samples", y = "Mammal species", colour = "Marker", shape = "Marker") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey80", size = 0.1),
    panel.grid.major = element_line(colour = "grey80", size = 0.2),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 14, colour = "black", face = "bold"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black", face = "bold"),
    legend.text = element_text(size = 12, colour = "black"))

plotgg1

ggsave(plot = plotgg1, 
       filename = "./output/FigS1_accum_curve_flies_markers.png", 
       dpi = 600, width = 6, height = 5)

