
## GLMs and PCOA of mammal community detected by flies and BF mosquitoes


### load packages
source("source_packages.R")



###########################################################################
################################ Fly pools ################################
###########################################################################

# importing and preparing the data
fly_dat<- read.csv("./data/fly_mammal_diversity.csv")
head(fly_dat)
str(fly_dat)

# prepare data for models
fly_dat_2 <- fly_dat %>% 
  mutate(tot_diversity = fly_dat %>% 
           dplyr::select(Canis_lupus:Rattus_norvegicus) %>% 
           rowSums(),
         potential_diversity = fly_dat %>% 
           dplyr::select(Canis_lupus:Rattus_norvegicus) %>% 
           ncol(), 
         missing_diversity = potential_diversity - tot_diversity, 
         rl_diversity = 42, # 42 = mammals in red list without bats
         missing_rl_diversity = rl_diversity - tot_diversity)  
fly_dat_2  

mean(fly_dat_2$tot_diversity)
min(fly_dat_2$tot_diversity)
max(fly_dat_2$tot_diversity)



##---------------------------------------##
## GLM models for total mammal diversity ##  
##---------------------------------------##


## Null Model 
mod_null <- glm(cbind(tot_diversity, missing_rl_diversity) ~ 1,
                data = fly_dat_2, 
                family = binomial())
summary(mod_null)


## Marker Model
mod_marker <- glm(cbind(tot_diversity, missing_rl_diversity) ~ Marker,
                  data = fly_dat_2, 
                  family = binomial())
summary(mod_marker)
plot(DHARMa::simulateResiduals(mod_marker))

## Locality Model 
mod_loc <- glm(cbind(tot_diversity, missing_rl_diversity) ~ Locality,
               data = fly_dat_2, 
               family = binomial())
summary(mod_loc)
plot(DHARMa::simulateResiduals(mod_loc))


## Marker + Locality Model 
mod_markerloc <- glm(cbind(tot_diversity, missing_rl_diversity) ~ Marker * Locality,
                     data = fly_dat_2, 
                     family = binomial())
summary(mod_markerloc)
plot(DHARMa::simulateResiduals(mod_markerloc))



## LRT to test explanatory variable effects

# marker
lrtest(mod_marker, mod_null) # p > 0.05 no effect

# locality
lrtest(mod_loc, mod_null) # p > 0.05 no effect

# marker and locality
lrtest(mod_markerloc, mod_null) # p > 0.05 no effect


## Plot effects locality x marker 
pred_markerloc <- ggpredict(mod_markerloc, terms = c("Locality", "Marker"))
names(pred_markerloc)

(myplot <- ggplot(pred_markerloc, aes(x = x, y = predicted, colour = group)) +
  geom_point(position=position_dodge(width=0.5), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.5, size = 1, 
                position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0, 0.5)) +
    scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8, 
                           name = "Marker") +
    ylim(0, 0.18) +
    theme_bw() +
    labs(x = "Locality",
         y = "Proportion detected") +
    theme(axis.title = element_text(size = 14, colour = "black", face = "bold"),
          axis.text =  element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          legend.text = element_text(size = 12, colour = "black")))

ggsave("./output/fig4_GLM_flyes_markerLoc_effects.png", plot = myplot, 
       dpi = 600, height = 5, width = 6)



##---------------------##
## PCoA BETA DIVERSITY ##
##---------------------##

# select only mammal species in data frame
my_species <- fly_dat %>% 
  dplyr::select(Canis_lupus:Rattus_norvegicus)
colSums(my_species) # double check there are no zeros

# distance matrix
beta_dist <- vegdist(my_species, index = "jaccard")

# Ordination: PCoA
pcoa_jac <- cmdscale(beta_dist, k = 2, eig = T)

pcoa_jac_plotting <- as.data.frame(pcoa_jac$points)
colnames(pcoa_jac_plotting) <- c("axis_1", "axis_2")
pcoa_jac_plotting$Locality <- fly_dat$Locality
pcoa_jac_plotting$Marker <- fly_dat$Marker

pcoa_jac$eig[1]/(sum(pcoa_jac$eig))
## [1] 00.5312973
pcoa_jac$eig[2]/(sum(pcoa_jac$eig))
## [1] 0.3038629


## Permanova for locality

# Test dispersion homogeneity between groups

permutest(betadisper(beta_dist, pcoa_jac_plotting$Marker))
permutest(betadisper(beta_dist, pcoa_jac_plotting$Locality))

# Test effect of locality 
adonis(my_species ~ Locality, data = pcoa_jac_plotting, permutations = 999, method = "jaccard")
# p > 0.05: no significant differences in community composition of localities


# Test effect of marker
adonis(my_species ~ Marker, data = pcoa_jac_plotting, permutations = 999, method = "jaccard")
# p = 0.002: Significant differences!


## Species to the plot 
bio.fit <- envfit(pcoa_jac, my_species, perm = 999)
bio.fit

# Get the vectors 
species_plotting <- scores(bio.fit,display=c("vectors"))
species_plotting <- as.data.frame(species_plotting)
species_plotting <- cbind.data.frame(species_plotting, pvals = bio.fit$vectors$pvals)

# assign colour based on significance
color_species <- ifelse(species_plotting$pvals < 0.05, "black", "#808080")
alpha_species <- ifelse(species_plotting$pvals < 0.05, 0.7, 0.5)

# get species common names
species_plotting2 <- species_plotting %>% 
  rownames_to_column("Species") %>% 
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
    Species == "Canis_lupus" ~ "Domestic dog",
    Species == "Dama_dama" ~ "Fallow deer"))

## FINAL PCoA
(plot_marker_final_english <- ggplot(pcoa_jac_plotting, aes(x = axis_1, y = axis_2)) +
    # geom_point(size = 3) +
    stat_ellipse(data = pcoa_jac_plotting, aes(x = axis_1, y = axis_2, fill = Marker),
                 geom="polygon",level=0.95,alpha=0.2) +
    geom_segment(data=species_plotting, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(0.2, "cm")), color = color_species, alpha=alpha_species) +
    geom_text(
      data = species_plotting,
      aes(Dim1 * 1.2, Dim2 * 1.2, label = species_plotting2$Species_english),
      color = color_species,
      alpha = alpha_species,
      hjust = c(0.5, 0.3, 0.5, 1, 0.5, 0.8, 1, 1, 0, 0.4, 0.1),
      vjust = c(0.5, 0.5, -0.5, 1, 2, -1, 0.5, 0.5, 0.5, 0.5, 0.5)
    ) +
    scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8, direction = 1) +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8, direction = -1) +
    geom_point(aes(x = axis_1, y = axis_2, colour = Locality, shape = Marker), size = 3) +
    theme_bw() + 
    xlab("PCoA 1 (53.1%)") +
    ylab("PCoA 2 (30.4%)") +
    theme(
      axis.title = element_text(colour = "black", face = "bold", size = 14),
      legend.title = element_text(colour = "black", face = "bold", size = 14),
      legend.text = element_text(colour = "black", size = 12)
    ))

ggsave(filename = "./output/Fig5_PCoA_flyes_markerLoc_english.png", plot = plot_marker_final_english, 
       dpi = 600, height = 6, width = 7)




############################################################
################### Blood Fed Mosquitoes ###################
############################################################

# importing and preparing the data
bf_dat <- read.csv("./data/bf_mammal_diversity.csv")
head(bf_dat)
str(bf_dat)

# prepare data for models, group by locality
bf_dat_2 <- bf_dat %>% 
  group_by(Locality, Date) %>% 
  summarise(Capreolus_capreolus = sum(Capreolus_capreolus),
            Dama_dama = sum(Dama_dama), 
            Sus_scrofa = sum(Sus_scrofa), 
            Vulpes_vulpes = sum(Vulpes_vulpes)) %>% 
  as.data.frame() %>% 
  mutate(tot_abund = rowSums(across(Capreolus_capreolus:Vulpes_vulpes)), 
         rl_diversity = 42) # 42 = mammals in red list without bats

bf_dat_2$tot_diversity <- rowSums(bf_dat_2[,3:6]>0)
bf_dat_2$missing_rl_diversity <- bf_dat_2$rl_diversity - bf_dat_2$tot_diversity


##---------------------------------------##
## GLM models for total mammal diversity ##  
##---------------------------------------##

## Null model
mod_null <- glm(cbind(tot_diversity, missing_rl_diversity) ~ 1,
                data = bf_dat_2, 
                family = binomial())
summary(mod_null)

## Locality model
mod_locality <- glm(cbind(tot_diversity, missing_rl_diversity) ~ Locality,
                  data = bf_dat_2, 
                  family = binomial())
summary(mod_locality)


## LRT to test explanatory variable effects

# locality
lrtest(mod_locality, mod_null) # p > 0.05 no effect

## Predict locality    
pred_locality <- ggpredict(mod_locality, terms = "Locality")
names(pred_locality)


##---------------------##
## PCoA BETA DIVERSITY ##
##---------------------##

## JACCARD ##

# select only species
my_species_bf <- bf_dat %>% 
  dplyr::select(Capreolus_capreolus:Vulpes_vulpes)
colSums(my_species_bf) # double check there are no zeros
rowSums(my_species_bf)

# remove rows with zeros
my_species_bf2 <- my_species_bf[rowSums(my_species_bf)!=0,]

bf_dat_pcoa <- bf_dat[rowSums(my_species_bf)!=0,]

# distance matrix
beta_dist_bf <- vegdist(my_species_bf2, index = "jaccard")

# Ordination: PCoA
pcoa_jac_bf <- cmdscale(beta_dist_bf, k = 2, eig = T)

pcoa_jac_plotting <- as.data.frame(pcoa_jac_bf$points)
colnames(pcoa_jac_plotting) <- c("axis_1", "axis_2")
pcoa_jac_plotting$Locality <- bf_dat_pcoa$Locality
pcoa_jac_plotting$Mosquito <- bf_dat_pcoa$Mosquito

pcoa_jac_bf$eig[1]/(sum(pcoa_jac_bf$eig))
## [1] 0.5791113
pcoa_jac_bf$eig[2]/(sum(pcoa_jac_bf$eig))
## [1] 0.3820594


## Permanova for locality
permutest(betadisper(beta_dist_bf, pcoa_jac_plotting$Locality))
plot(betadisper(beta_dist_bf, pcoa_jac_plotting$Locality))
adonis(my_species_bf2 ~ Locality, data = pcoa_jac_plotting, permutations = 999, method = "jaccard")
# p = 0.015: significant differences!

## Permanova for mosquito sp
permutest(betadisper(beta_dist_bf, pcoa_jac_plotting$Mosquito))
plot(betadisper(beta_dist_bf, pcoa_jac_plotting$Mosquito))
pcoa_jac_mosq <- pcoa_jac_plotting[!is.na(pcoa_jac_plotting$Mosquito),]
my_species_bf3 <- my_species_bf2[!is.na(pcoa_jac_plotting$Mosquito),]
adonis(my_species_bf3 ~ Mosquito, data = pcoa_jac_mosq, permutations = 999, method = "jaccard")
# p > 0.05: no significant effect




