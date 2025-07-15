### apply masks to cortmem data ###
## grab the significant edges of each matrices per condition per block
sigedges_pos_tc_cort = matrix(ncol = nrow(sigedges_pos_tc), nrow = length(cor_matrices_cort))
sigedges_pos_tc_plac = matrix(ncol = nrow(sigedges_pos_tc), nrow = length(cor_matrices_plac))
sigedges_neg_tc_cort = matrix(ncol = nrow(sigedges_neg_tc), nrow = length(cor_matrices_cort))
sigedges_neg_tc_plac = matrix(ncol = nrow(sigedges_neg_tc), nrow = length(cor_matrices_plac))

# Loop through each subject
for (sub in seq_along(cor_matrices_cort)) {
  cor_matrix_cort = cor_matrices_cort[[sub]]
  cor_matrix_plac = cor_matrices_plac[[sub]]
  # Loop through each row of sig_edges
  for (i in 1:nrow(sigedges_pos_tc)) {
    # Extract the node pairs from sig_edges
    node1_pos = sigedges_pos_tc[i,1]
    node2_pos = sigedges_pos_tc[i,2]
    # Extract the corresponding value from the correlation matrix and store it
    sigedges_pos_tc_cort[sub,i] = cor_matrix_cort[node1_pos, node2_pos]
    sigedges_pos_tc_plac[sub,i] = cor_matrix_plac[node1_pos, node2_pos]} 
  for (i in 1:nrow(sigedges_neg_tc)) {
    node1_neg = sigedges_neg_tc[i,1]
    node2_neg = sigedges_neg_tc[i,2]
    sigedges_neg_tc_cort[sub,i] = cor_matrix_cort[node1_neg, node2_neg]
    sigedges_neg_tc_plac[sub,i] = cor_matrix_plac[node1_neg, node2_neg]} }

# convert the result matrix to a data.frame
cort_pos = as.data.frame(sigedges_pos_tc_cort)
plac_pos = as.data.frame(sigedges_pos_tc_plac)
cort_neg = as.data.frame(sigedges_neg_tc_cort)
plac_neg = as.data.frame(sigedges_neg_tc_plac)

# set the row names to be the node pairs from sig_edges
colnames(cort_pos) = apply(sigedges_pos_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(plac_pos) = apply(sigedges_pos_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(cort_neg) = apply(sigedges_neg_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(plac_neg) = apply(sigedges_neg_tc, 1, function(row) paste0(row[1],"-",row[2]))

### compute avg strength
cortmem = read.csv("./CortMem_cortisol_data_allsub.csv")
cortmem = cortmem %>% select(Subject, week, Condition)
colnames(cortmem) = c("sub", "week", "Condition")
cortmem = distinct(cortmem)

sub = unique(cortmem$sub)
cort_pos$avg_nodes = rowMeans(cort_pos)
cort_pos$cond = "cort"
cort_pos$sub = sub
plac_pos$avg_nodes = rowMeans(plac_pos)
plac_pos$cond = "plac"
plac_pos$sub = sub
cort_neg$avg_nodes = rowMeans(cort_neg)
cort_neg$cond = "cort"
plac_neg$avg_nodes = rowMeans(plac_neg)
plac_neg$cond = "plac"

### cortmem motion ###
motion_cort = list.files(path = "./cortmem_motion", pattern = "\\.txt$", full.names = TRUE)
motion_cortmem = lapply(motion_cort, function(file) {
  read.table(file, header = FALSE, sep = "", fill = TRUE, stringsAsFactors = FALSE)})
# extract only rows with "rest1" from each participant's data
rest1_df = tibble(file = motion_cort) %>%
  mutate(sub = tools::file_path_sans_ext(basename(file))) %>%
  mutate(sub = str_sub(sub, 1, 4)) %>%    # keep only first 4 chars
  mutate(data = map(file, ~ read_table(.x, col_names = FALSE))) %>%
  mutate(data = map(data, ~ filter(.x, X2 == "rest1"))) %>%
  unnest(data) %>%
  rename(week = X1, run = X2, motion = X3)
# combine all into a single data frame
rest1_df = data.frame(rest1_df)[2:5]
rest1_df = unique(rest1_df)
rest1_df$sub = as.integer(rest1_df$sub)
cmem_motion = rest1_df 

cortmem = read.csv("./CortMem_cortisol_data_allsub.csv")
cortmem = cortmem %>% select(Subject, week, Condition)
colnames(cortmem) = c("sub", "week", "Group")
cortmem = cortmem %>%
  mutate(week = dplyr::recode(week, "Week1" = "week_a", "Week2" = "week_b"),
         Group = dplyr::recode(Group, "Cortisol" = "Stress", "Placebo" = "Control"))
cortmem = distinct(cortmem) %>% left_join(rest1_df)


##### permutation test ######

### setup data/sample ###
# mask
pos_mask_tc = read.csv("./cpm/net_masks/pos_mask_modage.csv", header=F)
neg_mask_tc = read.csv("./cpm/net_masks/neg_mask_modage.csv", header=F)

######### load labels for network ######### 

# create matrix of just the upper triangle
pos_edges_tc = data.frame(which(pos_mask_tc == 1 & upper.tri(pos_mask_tc, diag = FALSE), arr.ind = TRUE))
neg_edges_tc = data.frame(which(neg_mask_tc == 1 & upper.tri(neg_mask_tc, diag = FALSE), arr.ind = TRUE))
colnames(pos_edges_tc) = c("col1","col2")
colnames(neg_edges_tc) = c("col1","col2")


# load labels 
labels = read.csv("./misc/xilin_liz_combined.csv")
# put labels with sig edges
net_map = setNames(labels$net_names, labels$node)
roi_map = setNames(labels$BA_othername, labels$node)

# apply to sigedges files
sigedges_pos_tc_tc = pos_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
sigedges_neg_tc_tc = neg_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])

### set the networks
sigedges_df = rbind(sigedges_pos_tc_tc[1:2], sigedges_neg_tc_tc[1:2])

no_edges = nrow(sigedges_pos_tc) # set number of edges for pos net

### set the data working with 
# define the cor matrices data (this is for cortmem)
cor_matrices_Stress = cor_matrices_cort
cor_matrices_Control = cor_matrices_plac

pos_cortmem = rbind(cort_pos, plac_pos)
pos_cortmem$Group = ifelse(pos_cortmem$cond=="cort", "Stress", "Control")
pos_cortmem = left_join(pos_cortmem, cortmem %>% select(sub, motion, Group))
pos_cortmem$avg_nodes_resid = resid(lm(avg_nodes ~ motion, pos_cortmem)) # create motion residualized data

stress = pos_cortmem %>% filter(Group=="Stress")
control = pos_cortmem %>% filter(Group=="Control")

### draw random edges ###
### same network size as positive network (that we are testing) ###

set.seed(100)
# define the number of nodes and draws
num_nodes = 377
no_draws = 1000

# create a list to store results of each run
results_stress = vector("list", no_draws)
results_control = vector("list", no_draws)
# generate all possible pairs of nodes
all_pairs = t(combn(1:num_nodes, 2))
all_pairs_df = data.frame(all_pairs)
colnames(all_pairs_df) = c("node1", "node2")
# find pairs not in trauma network
colnames(sigedges_df) = c("node1", "node2")
non_sig_pairs = setdiff(all_pairs_df, sigedges_df)

### create summary files to store values from random draws ### 
summary_list = list(Stress = vector("list", no_draws), Control = vector("list", no_draws))
group_matrices = list(Stress = cor_matrices_Stress, Control = cor_matrices_Control)

### draw random network ###
# loop through number of iterations
for (run in 1:no_draws) {
  # sample random node pairs (edges) of the same number/size as cpm network that do not belong in the cpm net
  sampled_pairs = data.frame(non_sig_pairs[sample(nrow(non_sig_pairs), no_edges), ], row.names = NULL)
  # loop through each group (stress/control)
  for (group in names(group_matrices)) { 
    group_data = group_matrices[[group]]
    random_edges = matrix(ncol = no_edges, nrow = length(group_data))
    # extract edge values for each random edge
    for (sub in seq_along(group_data)) {
      cor_mat = group_data[[sub]]
      
      for (i in 1:nrow(sampled_pairs)) {
        node1 = sampled_pairs[[i, 1]]
        node2 = sampled_pairs[[i, 2]]
        random_edges[sub, i] = cor_mat[node1, node2]}}
    
    # create mean connectivity values
    df_random = as.data.frame(random_edges)
    subject_means = rowMeans(df_random, na.rm = TRUE)
    summary_df = data.frame(subject = 1:nrow(df_random),
                            mean_connectivity = subject_means,
                            df_random[1:no_edges])
    summary_list[[group]][[run]] = summary_df}}

#### do group comparisons (count neg for each edge) ####
# create summary files -- these contains subject no, mean conn, and each edge value of each random draw
stress_list = summary_list$Stress
control_list = summary_list$Control

# create results lists of group comparisons for each edge
perm_test_res = matrix(nrow = no_draws, ncol = 2)
colnames(perm_test_res) = c("count_neg", "coef_diff")

# loop through each iteration
for (i in 1:no_draws) {
  # extract the i-th data frame from each list
  stress_df = stress_list[[i]][, -c(1, 2)] # excl first two vars (subject & mean_conn)
  control_df = control_list[[i]][, -c(1, 2)]
  stress_df = cbind(stress_df, Group=stress$Group, motion=stress$motion, mean_conn=stress_list[[i]][, 2])
  control_df = cbind(control_df, Group=control$Group, motion=control$motion, mean_conn=control_list[[i]][, 2])
  df = rbind(stress_df, control_df)
  
  # residualize each edge for motion
  resid_stress = apply(stress_df[, 1:no_edges], 2, function(x) resid(lm(x ~ stress_df$motion)))
  resid_control = apply(control_df[, 1:no_edges], 2, function(x) resid(lm(x ~ control_df$motion)))
  
  # perform regression on each edge and grab the resulting coefficients
  beta = numeric(no_edges)
  for (j in 1:no_edges) {
    beta[j] = t.test(resid_stress[, j], resid_control[, j], paired = TRUE)$statistic}
  
  # count the number of times the mean difference is negative
  count_neg = sum(beta < 0)
  perm_test_res[i, 1] = count_neg
  
  # compute the mean difference
  model = lm(mean_conn ~ Group + motion, df)
  coef = model$coefficients[2]
  perm_test_res[i, 2] = coef}

# convert perm_test_res to a data frame
# perm_test_res is a list of no_draw x 2 that stores the number of edges w neg beta at each draw (count_neg) and mean diff (coef_diff)
perm_test_res = as.data.frame(perm_test_res) 
perm_test_res$count_neg_prop = perm_test_res$count_neg / no_edges * 100 # prop of random net edge w lower FC (stress-control) 

# find p-val 
hist(perm_test_res$count_neg_prop) # distribution of proportion of edges with negative betas (across x no_draws)
mean(perm_test_res$count_neg) # 49.732 -- mean value (from x random draws)
mean(perm_test_res$count_neg_prop) # 55.25778 -- mean value (prop)
# find the true value
mean_diff = numeric(no_edges)
predictor_df = pos_cortmem[ , 1:no_edges] 
pos_cortmem_resid = map_dfc(predictor_df, ~ resid(lm(.x ~ pos_cortmem$motion)))
names(pos_cortmem_resid) = paste0(names(predictor_df), "_resid")

stress_resid = data.frame(pos_cortmem_resid[1:27,])
cont_resid = data.frame(pos_cortmem_resid[28:54,])
for (j in 1:no_edges) {
  test_result = t.test(stress_resid[, j], cont_resid[, j], paired = TRUE)
  mean_diff[j] = test_result$estimate}
# create summary file
result = data.frame(edge = 1:no_edges, mean_diff = mean_diff)
# Count number of negative group effects (Stress < Control)
true_count_neg = sum(result$mean_diff < 0) # 58
true_count_neg_prop = true_count_neg / no_edges * 100 # 64.44

# compute p-val - probability of true count_neg to be more negative than random count_neg
1 - (length(which(true_count_neg > perm_test_res$count_neg)))/no_draws # 0.002

# look at the mean difference
mean(perm_test_res$coef_diff)
true_b = coef(lm(avg_nodes_resid ~ Group, data = pos_cortmem))["GroupStress"]
1 - (length(which(true_b < perm_test_res$coef_diff)))/no_draws # 0.497

## visualize 
perm_test_res %>% 
  ggplot(aes(x = count_neg_prop)) +
  geom_histogram(fill = "#5F5091", color = "#5F5091", alpha = 0.3, binwidth = 2) +
  geom_vline(xintercept = true_count_neg_prop, color = "black", linetype = "dashed", linewidth=1) +
  labs(x = "% Network Edge with Lower FC (Cortisol - Placebo)", y = "Frequency") + theme_classic()

# create bar plot
ggplot(result, aes(x = reorder(edge, mean_diff), y = mean_diff)) +
  geom_bar(stat = "identity", fill = "#5F5091") +
  coord_flip() + theme_classic() + theme(axis.text.y = element_blank()) +
  labs(x = "Positive network edge",
       y = "Mean Difference (Cort - Placebo)")

#### cort and PANAS levels between conditions
cortmem_panas = read.csv("./cortmem_panas_longformat.csv")
cortmem_panas$panas_change_neg = cortmem_panas$PANAS_neg_postscan - cortmem_panas$PANAS_neg_prescan
cortmem_cort = read.csv("./CortMem_cortisol_data_allsub.csv")
cortmem_cort = cortmem_cort %>% filter(desc==" Post_drug" | desc==" baseline") %>% select(Subject, Condition, desc, Avg) %>% 
  pivot_wider(names_from = desc, values_from = Avg)
cortmem_cort$cort_change = (cortmem_cort$` Post_drug`) - (cortmem_cort$` baseline`)
cortmem_cort %>% group_by(Condition) %>% summarize(round(mean(cort_change, na.rm=T),2), round(sd(cort_change, na.rm=T),2))
t.test((cortmem_cort %>% filter(Condition=="Cortisol"))$cort_change, (cortmem_cort %>% filter(Condition=="Placebo"))$cort_change, paired=T)

panas_cortmem = 
  cortmem_panas %>% 
  mutate(Pill = if_else(Pill == "Placebo", "Placebo", Pill)) %>%
  ggplot(aes(x = Pill, y = PANAS_neg_prescan, fill=Pill, color=Pill)) + 
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.7, 
               size = 0.75, outlier.shape = NA) +  # Boxplot with thicker outline and no outliers
  scale_fill_manual(values=c("Placebo" = "white", "Cortisol" = "#5F5091")) + # White fill for "Placebo"
  scale_color_manual(values=c("Placebo" = "#5F5091", "Cortisol" = "#5F5091")) + # Purple outline for both
  geom_point(position = position_jitter(width = 0.15), alpha = 0.2, size=1) + 
  labs(x="Condition", y="Negative Affect") + theme_classic() + theme(legend.position="none")
cort_cortmem =
  cortmem_cort %>% 
  mutate(Condition = if_else(Condition == "Placebo", "Placebo", Condition)) %>%  # No leading space here
  ggplot(aes(x = Condition, y = cort_change, fill=Condition, color=Condition)) + 
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.7, 
               size = 0.75, outlier.shape = NA) +  # Boxplot with thicker outline and no outliers
  scale_fill_manual(values=c("Placebo" = "white", "Cortisol" = "#5F5091")) + # White fill for "Placebo"
  scale_color_manual(values=c("Placebo" = "#5F5091", "Cortisol" = "#5F5091")) + # Purple outline for both
  geom_point(position = position_jitter(width = 0.15), alpha = 0.15, size=1) + 
  labs(x="Condition", y=expression(Delta~"Cortisol Levels (" * mu * "g/dL)")) + theme_classic() + theme(legend.position="none")
ggarrange(panas_cortmem, cort_cortmem, ncol=2, nrow=1)
