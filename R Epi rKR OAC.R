####################          Epidemiology of rKR          ####################

pacman::p_load(pacman, data.table, rio, tidyverse)

# Set working directory
directory <- "A:\DATA"
setwd(directory)



####################          Import and prepare NJR dataset          ####################

# Load data
load(file = "njr.RData")

# Convert to dataframes
rk <- as.data.frame(rkslim)
pk <- as.data.frame(pkslim)

# Remove duplicates (same day procedures)
rk <- 
  rk %>% 
  group_by(nn_nid, side, op_date) %>% 
  mutate(dups = row_number()) %>% 
  filter(dups==1) %>%
  select(-dups) %>% 
  ungroup()

no_uniq_rkr_supp <- nrow(rk)

# Create year fields & filter down to years 2006-2020
pk <- pk %>% mutate(ryear = format(primary_op_date, format = "%Y")) %>% filter(ryear >=2006, ryear<2021)
rk <- rk %>% mutate(ryear = format(op_date, format = "%Y")) %>% filter(ryear >=2006, ryear<2021)

# Count remaining records for data cleaning flowchart
good_dates <- nrow(rk)

# Create age fields
pk <- pk %>% rename(age = age_at_primary)
rk <- rk %>% rename(age = age_at_revision)

# Filter out NJR procedures in patients <18y
rk <- rk %>% filter(age>17 & age<105)
pk <- pk %>% filter(age>17 & age<105)

good_age <- nrow(rk)

# Create field to 'number' the revision procedure (revno)
# Staged revisions are considered to be one procedure if: (a) ITT as staged and (b) second stage <365 days
# However, we will ignore this initially and consider all procedures as independent

# Group by nn_nid, side -> Sort by op_date -> Number sequence -> Then link back to primary if possible
rk <- rk %>% arrange(nn_nid, side, op_date) %>% group_by(nn_nid, side) %>% mutate(seq = row_number())

# Get primary keys (primary_njr_index_no, primary_procedure_id) AND fields to identify ops
pk_id <- pk %>% select(nn_nid, side, primary_op_date, primary_njr_index_no, primary_procedure_id, patient_gender, age)

# The first rKR has seq==1, so create this as a merge field
pk_id <- pk_id %>% mutate(seq = 1)

# Merge
rk <- merge(rk, pk_id, by=c("nn_nid", "side", "seq"), all.x = TRUE)

# The first revision with a linked primary should be labelled revno==1
rk <- rk %>% mutate(dummy_revno = case_when(!is.na(primary_njr_index_no) ~ 1, TRUE ~ 0))

# For a given nn_nid, side - Did the first rKR link to a primary?
rk <- rk %>% group_by(nn_nid, side) %>% mutate(linked = max(dummy_revno))

# If it did link, then revno == seq
# Else, code it as 0
rk <- rk %>% mutate(revno = case_when(linked ==1 ~ seq, TRUE ~ 0L))

# Re-code as 0,1,2,3 and create factor
rk <- rk %>% mutate(revno = case_when(revno >3 ~ 3L, TRUE ~ revno))

rk$revno <- factor(rk$revno, levels=c(1,2,3,0), ordered=TRUE, labels = c("First linked rKR", "Second linked rKR", "Third or more linked rKR", "No linked primary"))

## Now consider staged revisions to be one procedure if: (a) ITT as staged and (b) second stage <365 days

rk <-
  rk %>% 
  group_by(nn_nid, side) %>% 
  arrange(nn_nid, side, op_date) %>% 
  mutate(prev_proc = lag(procedure_type, n=1, default = NA)) %>% 
  mutate(date_prev = lag(op_date, n=1, default = NA)) %>% 
  mutate(date_diff = op_date - date_prev) %>% 
  mutate(toseq = case_when(
    procedure_type == "Knee Stage 2 of 2 Stage Revision" & date_diff <365 & prev_proc == "Knee Stage 1 of 2 Stage Revision" ~ 0,
    TRUE ~ 1)) %>%
  filter(toseq==1) %>%
  arrange(nn_nid, side, op_date) %>% 
  group_by(nn_nid, side) %>% 
  mutate(seq1 = case_when(
    revno == "No linked primary" ~ 0L,
    TRUE ~ row_number())) %>% 
  select(-revno) %>% 
  rename(revno = seq1) %>% 
  mutate(revno = case_when(revno >3 ~ 3L, TRUE ~ revno))

rk$revno <- factor(rk$revno, levels=c(1,2,3,0), ordered=TRUE, labels = c("First linked rKR", "Second linked rKR", "Third or more linked rKR", "No linked primary"))

indep_rk <- nrow(rk)

save.image("epi.RData")



####################          Import and prepare ONS dataset          ####################

# Import ONS
ons <- import("MYEB1_detailed_population_estimates_series_UK_(2020_geog21).csv")

# Filter out Scotland
ons <- ons %>% filter(country!="S")

# Total population in 2006
pop2006 <-
  ons %>%
  summarise(pop = sum(population_2006))

# Total adult population in 2006
adultpop2006 <-
  ons %>%
  filter(age>17) %>% 
  summarise(pop = sum(population_2006))

# Total population in 2019
pop2019 <-
  ons %>%
  summarise(pop = sum(population_2019))

# Total adult population in 2019
adultpop2019 <-
  ons %>%
  filter(age>17) %>% 
  summarise(pop = sum(population_2019))

# Total population in 2020
pop2020 <-
  ons %>%
  summarise(pop = sum(population_2020))

# Total adult population in 2020
adultpop2020 <-
  ons %>%
  filter(age>17) %>% 
  summarise(pop = sum(population_2020))

# Sex: M=1, F=2 (i.e. opposite to NJR)
ons <- ons %>% mutate(patient_gender = recode(sex, `1` = "M", `2` = "F"))

# Remove those in the general population <18y
ons <- 
  ons %>% 
  filter(!age<18) 

# Reshape
ons$id <- 1:nrow(ons)

# Note that 'reshape' only works with balanced datasets
onsl <- 
  reshape(ons, 
          direction='long',
          idvar = "id",
          varying= 6:25, sep = "_")

# Collapse by year
ons_denom <- onsl %>% select(time, population) %>% group_by(time) %>% summarise_all(sum)
ons_denom <- ons_denom %>% rename(ryear = time) %>% filter(ryear >=2006, ryear<2021)



####################          Counts & Crude Incidence Rates          ####################

# Counts
pkr_num <- pk %>% count(ryear)

rk <- rk %>% ungroup()
rkr_num <- rk %>% count(ryear)

# Put numerator and denominator together
pkrpop <- merge(pkr_num, ons_denom, by="ryear", all.x = TRUE)
rkrpop <- merge(rkr_num, ons_denom, by="ryear", all.x = TRUE)

# Crude rates
p_load(PHEindicatormethods)
pkr_cr <-  pkrpop %>% group_by(ryear) %>% phe_rate(n, population, type = "full", confidence = 0.95, multiplier = 100000)
rkr_cr <-  rkrpop %>% group_by(ryear) %>% phe_rate(n, population, type = "full", confidence = 0.95, multiplier = 100000)

# Figures 1a and 1b
# Plot pKR and rKR counts with crude rates overlaid
scaleFactor_pk <- max(pkr_cr$n, na.rm=TRUE) / max(pkr_cr$value, na.rm=TRUE)

# pKR
plot.count.rate.pk <-
  pkr_cr %>% 
  ggplot()+
  geom_line(aes(x=ryear, y= value), col="black", group="method") +
  geom_point(aes(x=ryear, y= value), col="black") +
  geom_col(aes(x=ryear, y=n/scaleFactor_pk), alpha=0.4) +
  scale_y_continuous(name="Rate per 100,000 adult population", sec.axis=sec_axis(~.*scaleFactor_rk, name="Number of procedures")) +
  labs(x= "Year of surgery") +
  theme(
    legend.title = element_blank(),
    axis.title.y.left=element_text(colour="black"),
    axis.text.y.left=element_text(colour="black"),
    axis.title.y.right=element_text(color="grey"),
    axis.text.y.right=element_text(color="grey"),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Primary KR") + theme(plot.title = element_text(face = "bold", hjust=0.5)) +
  theme(plot.margin=unit(c(0.5,0.1,0.5,0.1), "cm"))

# rKR
scaleFactor_rk <- max(rkr_cr$n, na.rm=TRUE) / max(rkr_cr$value, na.rm=TRUE)

plot.count.rate.rk <-
  rkr_cr %>% 
  ggplot()+
  geom_line(aes(x=ryear, y= value), col="black", group="method") +
  geom_point(aes(x=ryear, y= value), col="black") +
  geom_col(aes(x=ryear, y=n/scaleFactor_rk), alpha=0.4) +
  scale_y_continuous(name="Rate per 100,000 adult population", sec.axis=sec_axis(~.*scaleFactor_rk, name="Number of procedures")) +
  labs(x= "Year of surgery") +
  theme(
    legend.title = element_blank(),
    axis.title.y.left=element_text(colour="black"),
    axis.text.y.left=element_text(colour="black"),
    axis.title.y.right=element_text(color="grey"),
    axis.text.y.right=element_text(color="grey"),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Revision KR") + theme(plot.title = element_text(face = "bold", hjust=0.5)) +
  theme(plot.margin=unit(c(0.5,0.1,0.5,0.1), "cm"))


cr_pr <- gridExtra::grid.arrange(plot.count.rate.pk, plot.count.rate.rk, nrow = 1)



####################          Crude rates by age group          ####################

# Cut age +/- sex in a few different ways to explore categorisation
rk <- rk %>% rename(age = age.x)
rk <- rk %>% rename(patient_gender = patient_gender.x)

df_list <- list(ons, pk, rk)
names(df_list) <- c("ons", "pk", "rk")

df_list <-
  df_list %>%
  map(~mutate(., 
              agegrp = cut(age, include.lowest = TRUE, breaks = c(18,50,55,60,65,70,75,80,85,110), labels = c("18-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+")),
              agegrp1 = cut(age, include.lowest = TRUE, breaks = c(18,50,60,70,80,110), labels = c("18-49", "50-59", "60-69", "70-79", "80+")),
              agesex = case_when(
                patient_gender == "F" & agegrp1 == "18-49" ~ "Female 18-49y",
                patient_gender == "F" & agegrp1 == "50-59" ~ "Female 50-59y",
                patient_gender == "F" & agegrp1 == "60-69" ~ "Female 60-69y",    
                patient_gender == "F" & agegrp1 == "70-79" ~ "Female 70-79y",
                patient_gender == "F" & agegrp1 == "80+" ~ "Female 80y+",
                patient_gender == "M" & agegrp1 == "18-49" ~ "Male 18-49y",
                patient_gender == "M" & agegrp1 == "50-59" ~ "Male 50-59y",
                patient_gender == "M" & agegrp1 == "60-69" ~ "Male 60-69y",    
                patient_gender == "M" & agegrp1 == "70-79" ~ "Male 70-79y",
                patient_gender == "M" & agegrp1 == "80+" ~ "Male 80y+",
                TRUE ~ NA_character_)
  ))

list2env(df_list, .GlobalEnv)

# ONS
ons$id <- 1:nrow(ons)
onsl <- 
  reshape(ons, 
          direction='long',
          idvar = "id",
          varying= 6:25, sep = "_")

ons_denom_age <- onsl %>% select(agegrp1, time, population) %>% group_by(agegrp1, time) %>% summarise_all(sum)

ons_denom_age <- ons_denom_age %>% rename(ryear = time) %>% filter(ryear >=2006, ryear<2021)

# Calculate crude rates
pkr_num <- pk %>% count(agegrp1,ryear) %>% complete(agegrp1,ryear) %>% filter(!is.na(agegrp1))
rkr_num <- rk %>% count(agegrp1,ryear) %>% complete(agegrp1,ryear) %>% filter(!is.na(agegrp1))

# Put numerator and denominator together
pkrpop <- merge(pkr_num, ons_denom_age, by=c("agegrp1","ryear"))
rkrpop <- merge(rkr_num, ons_denom_age, by=c("agegrp1","ryear"))

# Age-specific rates
pkr_asr <-  pkrpop %>% group_by(agegrp1, ryear) %>% phe_rate(n, population, type = "full", confidence = 0.95, multiplier = 100000)
rkr_asr <-  rkrpop %>% group_by(agegrp1, ryear) %>% phe_rate(n, population, type = "full", confidence = 0.95, multiplier = 100000)

# Figure 3
# Plot age-specific rates
df_list <- list(pkr_asr,rkr_asr)

plot.asr.agegrp1 <- 
  lapply(df_list, function(x) {
    x_plot <- x %>% 
      group_by(ryear,agegrp1) %>%
      filter(!is.na(value), !is.na(lowercl), !is.na(uppercl)) %>% 
      summarise_at(vars(value, lowercl, uppercl),
                   list(name = sum)) %>% 
      ggplot(aes(x=ryear, y = value_name, group=factor(agegrp1), colour=factor(agegrp1))) +
      geom_errorbar(aes(ymin=lowercl_name, ymax = uppercl_name), width=.1) +
      geom_point() +
      geom_line() +
      labs(x= "Year of surgery", y= "Rate per 100,000 population", colour="Age group (years)") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16))
  })

# Plot for revisions
plot.asr.agegrp1[[2]]



####################          Directly Standardised Rates          ####################

# Appendix B
# Create a truncated European Standard Population with combined groups
# "18-49", "50-59", "60-69", "70-79", "80+"
esp2013age <- c((sum(esp2013[5:10])+2200), sum(esp2013[11:12]), sum(esp2013[13:14]), sum(esp2013[15:16]), sum(esp2013[17:19]))

# Make this new ESP add back up to 100,000
esp2013age <- esp2013age * 100000/sum(esp2013age)

# Prep work above for agegrp1 valid
# Crude rates
pkr_cr <-  pkrpop %>% group_by(ryear) %>% phe_rate(n, population, type = "full", confidence = 0.95, multiplier = 100000) %>% mutate(joint=1)
rkr_cr <-  rkrpop %>% group_by(ryear) %>% phe_rate(n, population, type = "full", confidence = 0.95, multiplier = 100000) %>% mutate(joint=3)

# DSRs
pkr_dsr <-  pkrpop %>% group_by(ryear) %>% phe_dsr(n, population, stdpop = esp2013age, stdpoptype = "vector") %>% mutate(joint=2)
rkr_dsr <-  rkrpop %>% group_by(ryear) %>% phe_dsr(n, population, stdpop = esp2013age, stdpoptype = "vector") %>% mutate(joint=4)

pr.k <- rbind(pkr_cr, pkr_dsr, rkr_cr, rkr_dsr) %>%  select(ryear, value, lowercl, uppercl, joint)
pr.k$joint <- factor(pr.k$joint, levels = c(1,2,3,4), labels = c("Primary KR (Crude)", "Primary KR (DSR)","Revision KR (Crude)", "Revision KR (DSR)"))

dsr.pr.k <-
  ggplot(pr.k, aes(x=ryear, y= value, group=joint, colour=joint)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lowercl, ymax = uppercl), width=.1) +
  labs(x= "Year of surgery", y= "Rate per 100,000 adult population") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16)
    )

dsr.pr.k



####################          Linked and unlinked rKR          ####################

# Calculate proportion linked
prop_linked <-
  rk %>%
  select(linked, ryear) %>%
  group_by(ryear) %>%
  summarise_at(vars(linked),
               list(n_patients = length, linked = sum)) %>%
  mutate(unlinked = n_patients - linked) %>% 
  mutate(perc = round(linked * 100 / n_patients, 3)) %>% 
  arrange(ryear)

# Stacked bar of count linked/unlinked rKR (Figure 2)
plot.prop.linked <-
  prop_linked %>%
  rename(Linked = linked) %>% 
  rename(`Not linked` = unlinked) %>% 
  gather(key, value, -c(ryear,n_patients, perc)) %>% 
  ggplot() +
  geom_col(aes(ryear, value, fill = key), alpha = 0.5) +
  labs(x = "Year of surgery", y= "Number of procedures") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))



####################          Trends in main diagnosis over time (Crude rates)          ####################

# Create hierarchy to classify main indication for revision (ifrhier) from binary indicators in NJR dataset
rk <-
  rk %>% mutate(ifrhier = case_when(
    ind_rev_other==1 ~ 10,
    ind_rev_pain ==1 ~ 9,
    ind_rev_stiffness ==1 ~ 8,
    ind_rev_progressive_arthritis_remaining ==1 ~ 7,
    ind_rev_peri_prosthetic_fracture ==1 ~ 6,
    ind_rev_mds1patella_maltracking==1|ind_rev_component_dissociation ==1|ind_rev_dislocation_subluxation ==1|ind_rev_instability==1 ~ 5,
    ind_rev_mds1wear_of_tibia==1 | ind_rev_mds1wear_of_patella==1|ind_rev_implant_fracture_tibial==1|ind_rev_implant_fracture_patella ==1|ind_rev_implant_fracture_femoral ==1|ind_rev_implant_fracture ==1|ind_rev_wear_of_polyethylene_component==1 ~ 4,
    ind_rev_mds1lysis ==1|ind_rev_lysis_tibia ==1|ind_rev_lysis_femur ==1|ind_rev_mds1aseptic_loosening ==1|ind_rev_aseptic_loosening_patella ==1|ind_rev_aseptic_loosening_tibia==1|ind_rev_aseptic_loosening_femur  ==1 ~ 3,
    ind_rev_mds1incorrect_sizing==1|ind_rev_malalignment ==1 ~ 2,
    ind_rev_infection ==1 ~1,
    TRUE ~ 11))

# Label hierarchy
rk$ifrhier <- factor(rk$ifrhier, levels=c(1:11), ordered=TRUE, labels = c("Infection", "Malalignment", "Loosening/Lysis", "Component Wear", "Instability", "Fracture", "Progressive Arthritis", "Stiffness", "Unexplained pain", "Other", "Not specified")) 

# Create counts for each indication in each year
rkr_num <- rk %>% count(ryear,ifrhier) %>% complete(ryear,ifrhier)

# Put numerator and denominator together
rkrpop <- merge(rkr_num, ons_denom, by="ryear")

# Calculate crude rates
rkr_ifr_cr <-  rkrpop %>% group_by(ryear, ifrhier) %>% phe_rate(n, population, type = "full", confidence = 0.95, multiplier = 100000)

# Plot as grouped bar plot (Figure 4)
ind.time <-
  rkr_ifr_cr %>%
  select(ryear, ifrhier, value) %>% 
  rename(variable = ryear) %>% 
  ggplot(aes(x=ifrhier, y=(value), fill=forcats::fct_rev(variable))) +
  labs(x = "Diagnosis", y= "Rate per 100,000 adult population", fill = "Year of surgery") +
  scale_x_discrete(limits = rev(levels(rkr_ifr_cr$ifrhier))) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))



####################          Trends in main diagnosis over time (Proportions)          ####################

# Calculate proportional share of each diagnosis 
p_load(janitor)
ifr_year <-
  as.data.frame(rk %>% tabyl(ifrhier, ryear) %>% adorn_percentages("col"))

# Create a horizontal grouped barplot (Figure 8)
p_load(reshape2)
iyl_prop_time <- melt(ifr_year, id = "ifrhier")

# Proportion of each diagnosis by year for all rKR
ind.prop.time <-
  ggplot(data=iyl_prop_time, aes(x=ifrhier, y=(value*100), fill=forcats::fct_rev(variable))) +
  labs(x = "Diagnosis", y= "Annual percentage of revision diagnoses (%)", fill = "Year of surgery") +
  scale_x_discrete(limits = rev(levels(iyl_prop_time$ifrhier))) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))

ind.prop.time



####################          Trends in main diagnosis of re-revision procedures (2014-2019 data only)          ####################

# Filter down to years of interest and calculate percentage share
ifr1419 <- 
  as.data.frame(rk %>% filter(ryear >=2014, ryear<2020) %>% tabyl(ifrhier, revno) %>% adorn_percentages("col")) 

ifr_tot1419 <- 
  rk %>% filter(ryear >=2014, ryear<2020) %>% tabyl(ifrhier) %>% 
  select(ifrhier, percent) %>% 
  rename(`All rKR` = percent)

ifr1419 <-
  merge(ifr1419, ifr_tot1419, by="ifrhier")

ifr1419 <-
  ifr1419 %>% 
  select(ifrhier, `All rKR`, everything())

# Create a horizontal grouped barplot (Figure 5)
p_load(reshape2)
iyl1419 <- melt(ifr1419, id = "ifrhier")
iyl1419 <- iyl1419 %>% filter(ifrhier != "Not specified")

# As standard barplot
iyl1419$variable2 <- factor(iyl1419$variable, levels=rev(levels(iyl1419$variable)))
iyl1419$ifrhier <- droplevels(iyl1419$ifrhier)

ind.rr.bar1419 <-
  ggplot(data=iyl1419, aes(x=ifrhier, y=(value*100), colour=factor(variable2), fill=factor(variable2))) +
  labs(x = "Diagnosis", y= "Percentage of revision diagnoses") +
  scale_x_discrete(limits = rev(levels(iyl1419$ifrhier))) +
  scale_fill_discrete(guide=guide_legend(reverse=T)) +
  scale_colour_discrete(guide=guide_legend(reverse=T)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() +
  theme(legend.title = element_blank())

ind.rr.bar1419



####################          Table 1: Annual totals, Crude rates, Age-specific rates          ####################

# Set up headers
tot_head <- c("Annual totals of rKR procedures","","","","","","","","","","","","","","","")
inc_head <- c("Incidence rates of rKR","","","","","","","","","","","","","","","")
ar_head <- c("Age-specific incidence rate per 100,000 persons","","","","","","","","","","","","","","","")

# Count data
p_load(janitor)
rk <- as.data.frame(rk)
revno <- rk %>% tabyl(revno, ryear) %>% adorn_totals() %>% unname() %>% as.data.frame()

# Incidence rates
row_ir <- c("Crude incidence rate per 100,000 adults", 
            paste(round(unname(rkr_cr[1,4]),1), "(",round(unname(rkr_cr[1,5]),1),"-",round(unname(rkr_cr[1,6]),1),")"),
            paste(round(unname(rkr_cr[2,4]),1), "(",round(unname(rkr_cr[2,5]),1),"-",round(unname(rkr_cr[2,6]),1),")"),
            paste(round(unname(rkr_cr[3,4]),1), "(",round(unname(rkr_cr[3,5]),1),"-",round(unname(rkr_cr[3,6]),1),")"),
            paste(round(unname(rkr_cr[4,4]),1), "(",round(unname(rkr_cr[4,5]),1),"-",round(unname(rkr_cr[4,6]),1),")"),
            paste(round(unname(rkr_cr[5,4]),1), "(",round(unname(rkr_cr[5,5]),1),"-",round(unname(rkr_cr[5,6]),1),")"),
            paste(round(unname(rkr_cr[6,4]),1), "(",round(unname(rkr_cr[6,5]),1),"-",round(unname(rkr_cr[6,6]),1),")"),
            paste(round(unname(rkr_cr[7,4]),1), "(",round(unname(rkr_cr[7,5]),1),"-",round(unname(rkr_cr[7,6]),1),")"),
            paste(round(unname(rkr_cr[8,4]),1), "(",round(unname(rkr_cr[8,5]),1),"-",round(unname(rkr_cr[8,6]),1),")"),
            paste(round(unname(rkr_cr[9,4]),1), "(",round(unname(rkr_cr[9,5]),1),"-",round(unname(rkr_cr[9,6]),1),")"),
            paste(round(unname(rkr_cr[10,4]),1), "(",round(unname(rkr_cr[10,5]),1),"-",round(unname(rkr_cr[10,6]),1),")"),
            paste(round(unname(rkr_cr[11,4]),1), "(",round(unname(rkr_cr[11,5]),1),"-",round(unname(rkr_cr[11,6]),1),")"),
            paste(round(unname(rkr_cr[12,4]),1), "(",round(unname(rkr_cr[12,5]),1),"-",round(unname(rkr_cr[12,6]),1),")"),
            paste(round(unname(rkr_cr[13,4]),1), "(",round(unname(rkr_cr[13,5]),1),"-",round(unname(rkr_cr[13,6]),1),")"),
            paste(round(unname(rkr_cr[14,4]),1), "(",round(unname(rkr_cr[14,5]),1),"-",round(unname(rkr_cr[14,6]),1),")"),
            paste(round(unname(rkr_cr[15,4]),1), "(",round(unname(rkr_cr[15,5]),1),"-",round(unname(rkr_cr[15,6]),1),")"))

# Age-specific rates
ages <- rkr_asr %>% filter(agegrp1 == "18-49")

a1849 <- c("18-49 years",
           paste(round(unname(ages[1,5]),1), "(",round(unname(ages[1,6]),1),"-",round(unname(ages[1,7]),1),")"),
           paste(round(unname(ages[2,5]),1), "(",round(unname(ages[2,6]),1),"-",round(unname(ages[2,7]),1),")"),
           paste(round(unname(ages[3,5]),1), "(",round(unname(ages[3,6]),1),"-",round(unname(ages[3,7]),1),")"),
           paste(round(unname(ages[4,5]),1), "(",round(unname(ages[4,6]),1),"-",round(unname(ages[4,7]),1),")"),
           paste(round(unname(ages[5,5]),1), "(",round(unname(ages[5,6]),1),"-",round(unname(ages[5,7]),1),")"),
           paste(round(unname(ages[6,5]),1), "(",round(unname(ages[6,6]),1),"-",round(unname(ages[6,7]),1),")"),
           paste(round(unname(ages[7,5]),1), "(",round(unname(ages[7,6]),1),"-",round(unname(ages[7,7]),1),")"),
           paste(round(unname(ages[8,5]),1), "(",round(unname(ages[8,6]),1),"-",round(unname(ages[8,7]),1),")"),
           paste(round(unname(ages[9,5]),1), "(",round(unname(ages[9,6]),1),"-",round(unname(ages[9,7]),1),")"),
           paste(round(unname(ages[10,5]),1), "(",round(unname(ages[10,6]),1),"-",round(unname(ages[10,7]),1),")"),
           paste(round(unname(ages[11,5]),1), "(",round(unname(ages[11,6]),1),"-",round(unname(ages[11,7]),1),")"),
           paste(round(unname(ages[12,5]),1), "(",round(unname(ages[12,6]),1),"-",round(unname(ages[12,7]),1),")"),
           paste(round(unname(ages[13,5]),1), "(",round(unname(ages[13,6]),1),"-",round(unname(ages[13,7]),1),")"),
           paste(round(unname(ages[14,5]),1), "(",round(unname(ages[14,6]),1),"-",round(unname(ages[14,7]),1),")"),
           paste(round(unname(ages[15,5]),1), "(",round(unname(ages[15,6]),1),"-",round(unname(ages[15,7]),1),")"))

ages <- rkr_asr %>% filter(agegrp1 == "50-59")

a5059 <- c("50-59 years",
           paste(round(unname(ages[1,5]),1), "(",round(unname(ages[1,6]),1),"-",round(unname(ages[1,7]),1),")"),
           paste(round(unname(ages[2,5]),1), "(",round(unname(ages[2,6]),1),"-",round(unname(ages[2,7]),1),")"),
           paste(round(unname(ages[3,5]),1), "(",round(unname(ages[3,6]),1),"-",round(unname(ages[3,7]),1),")"),
           paste(round(unname(ages[4,5]),1), "(",round(unname(ages[4,6]),1),"-",round(unname(ages[4,7]),1),")"),
           paste(round(unname(ages[5,5]),1), "(",round(unname(ages[5,6]),1),"-",round(unname(ages[5,7]),1),")"),
           paste(round(unname(ages[6,5]),1), "(",round(unname(ages[6,6]),1),"-",round(unname(ages[6,7]),1),")"),
           paste(round(unname(ages[7,5]),1), "(",round(unname(ages[7,6]),1),"-",round(unname(ages[7,7]),1),")"),
           paste(round(unname(ages[8,5]),1), "(",round(unname(ages[8,6]),1),"-",round(unname(ages[8,7]),1),")"),
           paste(round(unname(ages[9,5]),1), "(",round(unname(ages[9,6]),1),"-",round(unname(ages[9,7]),1),")"),
           paste(round(unname(ages[10,5]),1), "(",round(unname(ages[10,6]),1),"-",round(unname(ages[10,7]),1),")"),
           paste(round(unname(ages[11,5]),1), "(",round(unname(ages[11,6]),1),"-",round(unname(ages[11,7]),1),")"),
           paste(round(unname(ages[12,5]),1), "(",round(unname(ages[12,6]),1),"-",round(unname(ages[12,7]),1),")"),
           paste(round(unname(ages[13,5]),1), "(",round(unname(ages[13,6]),1),"-",round(unname(ages[13,7]),1),")"),
           paste(round(unname(ages[14,5]),1), "(",round(unname(ages[14,6]),1),"-",round(unname(ages[14,7]),1),")"),
           paste(round(unname(ages[15,5]),1), "(",round(unname(ages[15,6]),1),"-",round(unname(ages[15,7]),1),")"))

ages <- rkr_asr %>% filter(agegrp1 == "60-69")

a6069 <- c("60-69 years",
           paste(round(unname(ages[1,5]),1), "(",round(unname(ages[1,6]),1),"-",round(unname(ages[1,7]),1),")"),
           paste(round(unname(ages[2,5]),1), "(",round(unname(ages[2,6]),1),"-",round(unname(ages[2,7]),1),")"),
           paste(round(unname(ages[3,5]),1), "(",round(unname(ages[3,6]),1),"-",round(unname(ages[3,7]),1),")"),
           paste(round(unname(ages[4,5]),1), "(",round(unname(ages[4,6]),1),"-",round(unname(ages[4,7]),1),")"),
           paste(round(unname(ages[5,5]),1), "(",round(unname(ages[5,6]),1),"-",round(unname(ages[5,7]),1),")"),
           paste(round(unname(ages[6,5]),1), "(",round(unname(ages[6,6]),1),"-",round(unname(ages[6,7]),1),")"),
           paste(round(unname(ages[7,5]),1), "(",round(unname(ages[7,6]),1),"-",round(unname(ages[7,7]),1),")"),
           paste(round(unname(ages[8,5]),1), "(",round(unname(ages[8,6]),1),"-",round(unname(ages[8,7]),1),")"),
           paste(round(unname(ages[9,5]),1), "(",round(unname(ages[9,6]),1),"-",round(unname(ages[9,7]),1),")"),
           paste(round(unname(ages[10,5]),1), "(",round(unname(ages[10,6]),1),"-",round(unname(ages[10,7]),1),")"),
           paste(round(unname(ages[11,5]),1), "(",round(unname(ages[11,6]),1),"-",round(unname(ages[11,7]),1),")"),
           paste(round(unname(ages[12,5]),1), "(",round(unname(ages[12,6]),1),"-",round(unname(ages[12,7]),1),")"),
           paste(round(unname(ages[13,5]),1), "(",round(unname(ages[13,6]),1),"-",round(unname(ages[13,7]),1),")"),
           paste(round(unname(ages[14,5]),1), "(",round(unname(ages[14,6]),1),"-",round(unname(ages[14,7]),1),")"),
           paste(round(unname(ages[15,5]),1), "(",round(unname(ages[15,6]),1),"-",round(unname(ages[15,7]),1),")"))

ages <- rkr_asr %>% filter(agegrp1 == "70-79")

a7079 <- c("70-79 years",
           paste(round(unname(ages[1,5]),1), "(",round(unname(ages[1,6]),1),"-",round(unname(ages[1,7]),1),")"),
           paste(round(unname(ages[2,5]),1), "(",round(unname(ages[2,6]),1),"-",round(unname(ages[2,7]),1),")"),
           paste(round(unname(ages[3,5]),1), "(",round(unname(ages[3,6]),1),"-",round(unname(ages[3,7]),1),")"),
           paste(round(unname(ages[4,5]),1), "(",round(unname(ages[4,6]),1),"-",round(unname(ages[4,7]),1),")"),
           paste(round(unname(ages[5,5]),1), "(",round(unname(ages[5,6]),1),"-",round(unname(ages[5,7]),1),")"),
           paste(round(unname(ages[6,5]),1), "(",round(unname(ages[6,6]),1),"-",round(unname(ages[6,7]),1),")"),
           paste(round(unname(ages[7,5]),1), "(",round(unname(ages[7,6]),1),"-",round(unname(ages[7,7]),1),")"),
           paste(round(unname(ages[8,5]),1), "(",round(unname(ages[8,6]),1),"-",round(unname(ages[8,7]),1),")"),
           paste(round(unname(ages[9,5]),1), "(",round(unname(ages[9,6]),1),"-",round(unname(ages[9,7]),1),")"),
           paste(round(unname(ages[10,5]),1), "(",round(unname(ages[10,6]),1),"-",round(unname(ages[10,7]),1),")"),
           paste(round(unname(ages[11,5]),1), "(",round(unname(ages[11,6]),1),"-",round(unname(ages[11,7]),1),")"),
           paste(round(unname(ages[12,5]),1), "(",round(unname(ages[12,6]),1),"-",round(unname(ages[12,7]),1),")"),
           paste(round(unname(ages[13,5]),1), "(",round(unname(ages[13,6]),1),"-",round(unname(ages[13,7]),1),")"),
           paste(round(unname(ages[14,5]),1), "(",round(unname(ages[14,6]),1),"-",round(unname(ages[14,7]),1),")"),
           paste(round(unname(ages[15,5]),1), "(",round(unname(ages[15,6]),1),"-",round(unname(ages[15,7]),1),")"))

ages <- rkr_asr %>% filter(agegrp1 == "80+")

a80 <- c("80+ years",
           paste(round(unname(ages[1,5]),1), "(",round(unname(ages[1,6]),1),"-",round(unname(ages[1,7]),1),")"),
           paste(round(unname(ages[2,5]),1), "(",round(unname(ages[2,6]),1),"-",round(unname(ages[2,7]),1),")"),
           paste(round(unname(ages[3,5]),1), "(",round(unname(ages[3,6]),1),"-",round(unname(ages[3,7]),1),")"),
           paste(round(unname(ages[4,5]),1), "(",round(unname(ages[4,6]),1),"-",round(unname(ages[4,7]),1),")"),
           paste(round(unname(ages[5,5]),1), "(",round(unname(ages[5,6]),1),"-",round(unname(ages[5,7]),1),")"),
           paste(round(unname(ages[6,5]),1), "(",round(unname(ages[6,6]),1),"-",round(unname(ages[6,7]),1),")"),
           paste(round(unname(ages[7,5]),1), "(",round(unname(ages[7,6]),1),"-",round(unname(ages[7,7]),1),")"),
           paste(round(unname(ages[8,5]),1), "(",round(unname(ages[8,6]),1),"-",round(unname(ages[8,7]),1),")"),
           paste(round(unname(ages[9,5]),1), "(",round(unname(ages[9,6]),1),"-",round(unname(ages[9,7]),1),")"),
           paste(round(unname(ages[10,5]),1), "(",round(unname(ages[10,6]),1),"-",round(unname(ages[10,7]),1),")"),
           paste(round(unname(ages[11,5]),1), "(",round(unname(ages[11,6]),1),"-",round(unname(ages[11,7]),1),")"),
           paste(round(unname(ages[12,5]),1), "(",round(unname(ages[12,6]),1),"-",round(unname(ages[12,7]),1),")"),
           paste(round(unname(ages[13,5]),1), "(",round(unname(ages[13,6]),1),"-",round(unname(ages[13,7]),1),")"),
           paste(round(unname(ages[14,5]),1), "(",round(unname(ages[14,6]),1),"-",round(unname(ages[14,7]),1),")"),
           paste(round(unname(ages[15,5]),1), "(",round(unname(ages[15,6]),1),"-",round(unname(ages[15,7]),1),")"))

# Bind table together
t1 <- rbind(tot_head,revno, inc_head, row_ir, ar_head, a1849,a5059,a6069,a7079,a80)

# Set up rownames/colnames
colnames(t1) <- c(" ","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020")

# Change order and rename total to 'All rKR'
t1 <- t1[c(1,6,5,2,3,4,7:nrow(t1)), ]
t1[2,1] <- "All rKR" 
rownames(t1) <- c()

write.csv(t1,"A:\DATA\epit1.csv", row.names = FALSE)



####################          Export images          ####################

setwd("A:\DATA\figs")

plotlist = list()
plotlist[[1]] = plot.count.rate.pk
plotlist[[2]] = plot.count.rate.rk
plotlist[[3]] = plot.asr.agegrp1[[2]]
plotlist[[4]] = plot_rrkr_cr
plotlist[[5]] = ind.time
plotlist[[6]] = ind.rr.bar1419
plotlist[[7]] = dsr.pr.k
plotlist[[8]] = ind.prop.time
plotlist[[9]] = plot.prop.linked
plotlist[[10]] = cr_pr

for(i in 1:10) {
    ggsave(plot = plotlist[[i]], 
    width = 200,
    height = 200,
    units = "mm",
    dpi = 1000,
    file = paste("file",i,".png",sep=""))
  }

setwd(directory)

save.image("epi.RData")