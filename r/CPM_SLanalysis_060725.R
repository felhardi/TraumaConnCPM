### apply the TC pos mask to the SL data ###
# mask
pos_mask_tc = read.csv("./cpm/net_masks/pos_mask_modage.csv", header=F)
neg_mask_tc = read.csv("./cpm/net_masks/neg_mask_modage.csv", header=F)

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
sigedges_pos_tc = pos_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
sigedges_neg_tc = neg_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
## grab the significant edges of each matrices per condition per block
sigedges_pos_tc_block1 = matrix(ncol = nrow(sigedges_pos_tc), nrow = length(cor_matrices_sl_block1))
sigedges_pos_tc_block2 = matrix(ncol = nrow(sigedges_pos_tc), nrow = length(cor_matrices_sl_block2))
sigedges_pos_tc_block3 = matrix(ncol = nrow(sigedges_pos_tc), nrow = length(cor_matrices_sl_block3))
sigedges_pos_tc_block4 = matrix(ncol = nrow(sigedges_pos_tc), nrow = length(cor_matrices_sl_block4))

sigedges_neg_tc_block1 = matrix(ncol = nrow(sigedges_neg_tc), nrow = length(cor_matrices_sl_block1))
sigedges_neg_tc_block2 = matrix(ncol = nrow(sigedges_neg_tc), nrow = length(cor_matrices_sl_block2))
sigedges_neg_tc_block3 = matrix(ncol = nrow(sigedges_neg_tc), nrow = length(cor_matrices_sl_block3))
sigedges_neg_tc_block4 = matrix(ncol = nrow(sigedges_neg_tc), nrow = length(cor_matrices_sl_block4))

# loop through each subject
for (sub in seq_along(cor_matrices_sl_block1)) {
  cor_matrix1 = cor_matrices_sl_block1[[sub]]
  cor_matrix2 = cor_matrices_sl_block2[[sub]]
  cor_matrix3 = cor_matrices_sl_block3[[sub]]
  cor_matrix4 = cor_matrices_sl_block4[[sub]]
  # Loop through each row of sig_edges
  for (i in 1:nrow(sigedges_pos_tc)) {
    # Extract the node pairs from sig_edges
    node1_pos = sigedges_pos_tc[i,1]
    node2_pos = sigedges_pos_tc[i,2]
    # Extract the corresponding value from the correlation matrix and store it
    sigedges_pos_tc_block1[sub,i] = cor_matrix1[node1_pos, node2_pos]
    sigedges_pos_tc_block2[sub,i] = cor_matrix2[node1_pos, node2_pos]
    sigedges_pos_tc_block3[sub,i] = cor_matrix3[node1_pos, node2_pos]
    sigedges_pos_tc_block4[sub,i] = cor_matrix4[node1_pos, node2_pos]} 
  for (i in 1:nrow(sigedges_neg_tc)) {
    node1_neg = sigedges_neg_tc[i,1]
    node2_neg = sigedges_neg_tc[i,2]
    sigedges_neg_tc_block1[sub,i] = cor_matrix1[node1_neg, node2_neg]
    sigedges_neg_tc_block2[sub,i] = cor_matrix2[node1_neg, node2_neg]
    sigedges_neg_tc_block3[sub,i] = cor_matrix3[node1_neg, node2_neg]
    sigedges_neg_tc_block4[sub,i] = cor_matrix4[node1_neg, node2_neg]} }

# convert the result matrix to a data.frame
block1_pos = as.data.frame(sigedges_pos_tc_block1)
block2_pos = as.data.frame(sigedges_pos_tc_block2)
block3_pos = as.data.frame(sigedges_pos_tc_block3)
block4_pos = as.data.frame(sigedges_pos_tc_block4)

block1_neg = as.data.frame(sigedges_neg_tc_block1)
block2_neg = as.data.frame(sigedges_neg_tc_block2)
block3_neg = as.data.frame(sigedges_neg_tc_block3)
block4_neg = as.data.frame(sigedges_neg_tc_block4)

# set the row names to be the node pairs from sig_edges
colnames(block1_pos) = apply(sigedges_pos_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block2_pos) = apply(sigedges_pos_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block3_pos) = apply(sigedges_pos_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block4_pos) = apply(sigedges_pos_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block1_neg) = apply(sigedges_neg_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block2_neg) = apply(sigedges_neg_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block3_neg) = apply(sigedges_neg_tc, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block4_neg) = apply(sigedges_neg_tc, 1, function(row) paste0(row[1],"-",row[2]))

# create average score (all nodes)
block1_pos$pos_edgestr_block1 = rowMeans(block1_pos)
block2_pos$pos_edgestr_block2 = rowMeans(block2_pos)
block3_pos$pos_edgestr_block3 = rowMeans(block3_pos)
block4_pos$pos_edgestr_block4 = rowMeans(block4_pos)

block1_neg$neg_edgestr_block1 = rowMeans(block1_neg)
block2_neg$neg_edgestr_block2 = rowMeans(block2_neg)
block3_neg$neg_edgestr_block3 = rowMeans(block3_neg)
block4_neg$neg_edgestr_block4 = rowMeans(block4_neg)

# create a combined data frame with all the summary edge
avgedge_allblocks_tc =
  as.data.frame(cbind(subid_sl,
                      block1_pos$pos_edgestr_block1, block2_pos$pos_edgestr_block2, block3_pos$pos_edgestr_block3, block4_pos$pos_edgestr_block4,
                      block1_neg$neg_edgestr_block1, block2_neg$neg_edgestr_block2, block3_neg$neg_edgestr_block3, block4_neg$neg_edgestr_block4), 
                stringsAsFactors = FALSE)
colnames(avgedge_allblocks_tc) = c("sub","pos_b1","pos_b2","pos_b3","pos_b4",
                                   "neg_b1","neg_b2","neg_b3","neg_b4")
avgedge_allblocks_tc = avgedge_allblocks_tc %>% mutate_at(c(1:9), as.numeric)

### look by condition ### 
# extract the subject IDs from sub_data_sl and format them to match SL_cond
SL_cond = read.csv("./misc/SL_stresscon.csv")
SL_cond$sub = as.numeric(substring(SL_cond$StressLearn_ID, 2,4))
trauma = read.csv("./trauma.csv")
motion = read.csv("./motion.csv")

# merge with avg edge data frame
avgedge_allblocks_tc = left_join(avgedge_allblocks_tc, SL_cond %>% select(sub, Group))
avgedge_allblocks_tc = left_join(avgedge_allblocks_tc, trauma)
avgedge_allblocks_tc = left_join(avgedge_allblocks_tc, motion)
avgedge_allblocks_tc$scan = ifelse(avgedge_allblocks_tc$sub > 55, 1, 0) # scanner type


##### examine group differences at each run ####
### main sl analysis ###
summary(lm(pos_b1 ~ Group, avgedge_allblocks_tc))
summary(lm(pos_b2 ~ Group, avgedge_allblocks_tc))
summary(lm(pos_b3 ~ Group, avgedge_allblocks_tc))
summary(lm(pos_b4 ~ Group, avgedge_allblocks_tc))
p.adjust(c(0.291, 0.00513, 0.165, 0.0369), method = "fdr")

# control for covariates
summary(lm(pos_b1 ~ Group + sex + motion + scan, avgedge_allblocks_tc))
summary(lm(pos_b2 ~ Group + sex + motion + scan, avgedge_allblocks_tc))
summary(lm(pos_b3 ~ Group + sex + motion + scan, avgedge_allblocks_tc))
summary(lm(pos_b4 ~ Group + sex + motion + scan, avgedge_allblocks_tc))
p.adjust(c(0.263, 0.00518, 0.168, 0.0389), method = "fdr")

### visualize
avgedge_allblocks_tc %>%
  select(pos_b1, pos_b2, pos_b3, pos_b4, neg_b1, neg_b2, neg_b3, neg_b4, Group) %>%
  pivot_longer(cols = -Group, values_to = "fc", 
               names_to = c("network","block"), names_pattern="(.+)_(b[1234])") %>% 
  mutate(network = dplyr::recode(network, pos = "Positive", neg = "Negative")) %>%
  mutate(block = dplyr::recode(block, b1 = "Run 1\n(Baseline)", b2 = "Run 2", b3 = "Run 3", b4 = "Run 4")) %>%
  mutate(network = factor(network, levels = c("Positive", "Negative"))) %>%  
  ggplot(aes(x = block, y = fc, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width = 0.9), width = 0.3) + 
  scale_fill_manual(values=c("#588b8b", "#f28f3b")) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.15, size=1, aes(y=fc)) + 
  scale_color_manual(values=c("#588b8b", "#f28f3b")) + 
  labs(x="Time", y="Network Mean FC") + facet_grid(network ~ ., scales = "free_y") + theme_classic() + 
  theme(legend.position="none") + geom_vline(xintercept = 3, linetype = "dashed", color = "grey70", size = 0.5) +
  scale_x_discrete(
    limits = c("", "Run 1\n(Baseline)", "", "Run 2", "Run 3", "Run 4", ""), 
    expand = expansion(add = c(0.05, 0.05)))


##### permutation test: see if more than chance / random network #####
### set parameters ###
sigedges_df = rbind(sigedges_pos_tc[1:2], sigedges_neg_tc[1:2])
no_edges = nrow(sigedges_pos_tc) # set number of edges for pos net

avgedge_allblocks = avgedge_allblocks_tc

# define the cor matrices data (this is for SL, block 2)
cor_matrices_sl_block2_with_groups = mapply(
  function(matrix, group) list(matrix = matrix, group = group),
  cor_matrices_sl_block2, avgedge_allblocks$Group, SIMPLIFY = FALSE)
cor_matrices_Stress = cor_matrices_sl_block2_with_groups %>% keep(~ .x$group == "Stress") %>% map("matrix")
cor_matrices_Control = cor_matrices_sl_block2_with_groups %>% keep(~ .x$group == "Control") %>% map("matrix")
stress = avgedge_allblocks %>% filter(Group =="Stress")
control = avgedge_allblocks %>% filter(Group =="Control")

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
  # loop through each group (stress/control)
  for (group in names(group_matrices)) { 
    group_data = group_matrices[[group]]
    random_edges = matrix(ncol = no_edges, nrow = length(group_data))
    # extract edge values for each random edge
    for (sub in seq_along(group_data)) {
      cor_mat = group_data[[sub]]
      # sample random node pairs (edges) of the same number/size as cpm network that do not belong in the cpm net
      sampled_pairs = data.frame(non_sig_pairs[sample(nrow(non_sig_pairs), no_edges), ], row.names = NULL)
      for (i in 1:nrow(sampled_pairs)) {
        node1 = sampled_pairs[[i, 1]]
        node2 = sampled_pairs[[i, 2]]
        random_edges[sub, i] = cor_mat[node1, node2]}}
    # create mean connectivity values
    df_random = as.data.frame(random_edges)
    subject_means = rowMeans(df_random)
    summary_df = data.frame(subject = 1:nrow(df_random),
                            mean_connectivity = subject_means,
                            df_random[1:no_edges])
    summary_list[[group]][[run]] = summary_df}}

# create summary files -- these contains subject no, mean conn, and each edge value of each random draw
stress_list = summary_list$Stress
control_list = summary_list$Control

### do group comparisons (count neg for each edge) ####
# create results lists of group comparisons for each edge
perm_test_res = matrix(nrow = no_draws, ncol = 2)
colnames(perm_test_res) = c("count_neg", "coef_diff")

# loop through each iteration
for (i in 1:no_draws) {
  # extract the i-th data frame from each list
  stress_df = stress_list[[i]][, -c(1, 2)] # excl first two vars (subject & mean_conn)
  control_df = control_list[[i]][, -c(1, 2)]
  stress_df = cbind(stress_df, Group=stress$Group, motion=stress$motion, sex=stress$sex, scan=stress$scan, mean_conn=stress_list[[i]][, 2])
  control_df = cbind(control_df, Group=control$Group, motion=control$motion, sex=control$sex, scan=control$scan, mean_conn=control_list[[i]][, 2])
  df = rbind(stress_df, control_df)
  
  # perform regression on each edge and grab the resulting coefficients
  beta = numeric(no_edges)
  for (j in 1:no_edges) {
    reg_result = lm(df[, j] ~ Group + motion + sex + scan, df)
    beta[j] = coef(reg_result)["GroupStress"]}
  # count the number of times the mean difference is negative
  count_neg = sum(beta < 0)
  perm_test_res[i, 1] = count_neg
  
  # compute the mean difference
  model = lm(mean_conn ~ Group + motion + sex + scan, df)
  coef = coef(model)["GroupStress"]
  perm_test_res[i, 2] = coef}

# convert perm_test_res to a data frame
# perm_test_res is a list of no_draw x 2 that stores the number of edges w neg beta at each draw (count_neg) and mean diff (coef_diff)
perm_test_res = as.data.frame(perm_test_res) 
perm_test_res$count_neg_prop = perm_test_res$count_neg / no_edges * 100 # prop of random net edge w lower FC (stress-control) 

# find p-val 
hist(perm_test_res$count_neg_prop) # distribution of proportion of edges with negative betas (across x no_draws)
mean(perm_test_res$count_neg) # mean value (from x random draws)
mean(perm_test_res$count_neg_prop) # mean value (proportion of neg edges)
# find the true value
block2_pos_alledge = cbind(block2_pos[1:no_edges], Group=avgedge_allblocks$Group, motion=avgedge_allblocks$motion, sex=avgedge_allblocks$sex)
# find beta value and proportion
# apply linear regression to each edge (adjusting for motion and sex)
lm_results = apply(block2_pos_alledge[ , 1:no_edges], 2, function(edge_values) {
  lm(edge_values ~ Group + motion + sex, data = block2_pos_alledge)})
# extract p-values, coef, and fdr
p_values = sapply(lm_results, function(model) summary(model)$coefficients["GroupStress", "Pr(>|t|)"])
group_coefs = sapply(lm_results, function(model) coef(model)["GroupStress"])
fdr_p = p.adjust(p_values, method = "fdr")
# combine into a result dataframe
result = data.frame(Variable = names(p_values), Group_Effect = group_coefs, P_Value = p_values, fdr_p = fdr_p)
# count number of negative group effects (Stress < Control)
true_count_neg = sum(result$Group_Effect < 0) # 87 (true count)
true_count_neg_prop = true_count_neg / no_edges * 100 # 96.66667 (true proportion)

# compute p-val - probability of true count_neg to be more negative than random count_neg
1 - (length(which(true_count_neg > perm_test_res$count_neg)))/no_draws # 0.038

##### visualize permutation test #####
## visualize (tc = #f28f3b (bin: 1.5))
perm = perm_test_res %>% 
  ggplot(aes(x = count_neg_prop)) +
  geom_histogram(fill = "#f28f3b", color = "#f28f3b", alpha = 0.3, binwidth = 1.5) +
  geom_vline(xintercept = true_count_neg_prop, color = "black", linetype = "dashed", linewidth=1) +
  labs(x = "% Network Edge with Lower FC (Stress - Control)", y = "Frequency") + theme_classic()
edgeplot = ggplot(result, aes(x = reorder(Variable, Group_Effect), y = Group_Effect)) +
  geom_bar(stat = "identity", fill = "#fbd3aa") +
  coord_flip() + theme_classic() + theme(axis.text.y = element_blank()) +
  labs(x = "Positive Network Edge",
       y = "Mean Difference (Stress - Control)")
ggarrange(edgeplot, perm)


##### test association with depressive symptoms #####
# PHQ9
PHQ = read.csv("./qualtrics/StressMem fMRI Questionnaires_October 13, 2024_17.33.csv", stringsAsFactors=F)
PHQ = PHQ %>% select(src_subject_id, starts_with("phq")) %>% slice(-1:-3) %>%
  filter(if_all(everything(), ~ . != "")) %>%
  mutate(across(starts_with("phq"), ~ dplyr::recode(.x, 
                                                    "Not at all" = 0,
                                                    "Several days" = 1,
                                                    "More than half the days" = 2,
                                                    "Nearly every day" = 3,
                                                    .default = NA_real_))) %>% 
  mutate(sub = as.numeric(gsub("s(\\d+)", "\\1", src_subject_id))) 
PHQ$PHQ = rowSums(PHQ[2:10])
PHQ[95,]$sub = 95 # duplicate

# merge with main df
avgedge_allblocks_tc = left_join(avgedge_allblocks_tc, PHQ %>% select(sub, PHQ))

# test linear association
summary(lm(PHQ ~ pos_b1, avgedge_allblocks_tc %>% filter(Group=="Stress")))
summary(lm(PHQ ~ pos_b1, avgedge_allblocks_tc %>% filter(Group=="Control")))
summary(lm(PHQ ~ pos_b2, avgedge_allblocks_tc %>% filter(Group=="Stress")))
summary(lm(PHQ ~ pos_b2, avgedge_allblocks_tc %>% filter(Group=="Control")))

# test nonlinear associations
library(mgcv)
avgedge_allblocks_tc$Group = as.factor(avgedge_allblocks_tc$Group)

summary(gam(PHQ ~ s(pos_b2, by=Group) + Group, data = avgedge_allblocks_tc))
summary(gam(PHQ ~ s(pos_b2, by=Group) + Group + s(pos_b1), data = avgedge_allblocks_tc)) # control for baseline
summary(gam(PHQ ~ s(pos_b2, by=Group) + Group + s(pos_b1) + sex, data = avgedge_allblocks_tc)) # control for sex

# compare linear vs nonlinear
lin = (lm(PHQ ~ pos_b2 * Group + pos_b1, avgedge_allblocks_tc))
nonlin = (gam(PHQ ~ s(pos_b2, by=Group) + Group + s(pos_b1), data = avgedge_allblocks_tc))
AIC(lin, nonlin)

# visualize
# bar plot (median split for viz only)
bar_med = avgedge_allblocks_tc %>%
  mutate(PHQ_high = ifelse(PHQ > median(PHQ, na.rm = TRUE), 1, 0)) %>% 
  select(pos_b2, Group, PHQ_high) %>%
  pivot_longer(cols = -c(Group, PHQ_high), names_to = c("network", "block"), 
               names_pattern = "(.+)_(b[2])", values_to = "fc") %>%
  mutate(Group = as.factor(Group), PHQ = as.factor(PHQ_high)) %>%
  mutate(PHQ = dplyr::recode(PHQ, "0" = "Low", "1" = "High")) %>%
  mutate(block = dplyr::recode(block, b2 = "Run 2")) %>%
  ggplot(aes(x = interaction(Group, block, sep = " "), y = fc, fill = Group, pattern = PHQ)) +
  geom_bar_pattern(stat = "summary", fun = "mean", position = position_dodge(width = 0.7), 
                   pattern_density = 0.1, pattern_fill = "white", pattern_spacing = 0.025, width = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.7), width = 0.2) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), size = 1, alpha = 0.15) +
  scale_fill_manual(values = c("#588b8b", "#f28f3b")) + scale_pattern_manual(values = c("none", "stripe")) +
  labs(x = "Post Acute Stress Positive Network FC", y = "Positive Network Mean FC") +
  theme_classic() + guides(fill = guide_legend(override.aes = list(pattern = "none"))) + 
  geom_hline(aes(yintercept = 0.2425565), color = "grey50", linetype = "dashed") +
  scale_x_discrete(labels = c("Control Run 2", "Stress Run 2"))

# line plot
avgedge_allblocks_tc$pos_b2res = lm(pos_b2 ~ pos_b1, avgedge_allblocks_tc)$resid # create residualized variable
phq_line = avgedge_allblocks_tc %>%
  ggplot(aes(x=pos_b2res, y=PHQ, color=Group)) + 
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs"), aes(fill = Group), alpha = 0.2) + 
  geom_point(alpha = 0.3, size=1.5) + facet_wrap(~ Group) + theme_classic() + 
  theme(legend.position = "top") + 
  labs(x = "Post Acute Stress Positive Network FC (baseline-adjusted)", y = "PHQ", color = "Group", fill = "Group") +
  scale_color_manual(values = c("Control" = "#588b8b", "Stress" = "#f28f3b")) +
  scale_fill_manual(values = c("Stress" = "grey60", "Control" = "grey60"))

ggarrange(bar_med, phq_line)
