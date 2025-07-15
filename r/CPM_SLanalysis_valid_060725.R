###### apply KN posmask to sl data #####
# load CPM sig edges (KN)
pos_mask_kn = read.csv("../net_masks/pos_mask_kn.csv", header=F)
neg_mask_kn = read.csv("../net_masks/neg_mask_kn.csv", header=F)
# create matrix of just the upper triangle
pos_edges_kn = data.frame(which(pos_mask_kn == 1 & upper.tri(pos_mask_kn, diag = FALSE), arr.ind = TRUE))
neg_edges_kn = data.frame(which(neg_mask_kn == 1 & upper.tri(neg_mask_kn, diag = FALSE), arr.ind = TRUE))
colnames(pos_edges_kn) = c("col1","col2")
colnames(neg_edges_kn) = c("col1","col2")

# put labels with sig edges
labels = read.csv("./xilin_liz_combined.csv")
net_map = setNames(labels$net_names, labels$node)
roi_map = setNames(labels$BA_othername, labels$node)
sigedges_pos_kn = pos_edges_kn %>% 
  mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
sigedges_neg_kn = neg_edges_kn %>% 
  mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])

## grab the significant edges of each matrices per condition per block
sigedges_pos_block1 = matrix(ncol = nrow(sigedges_pos_kn), nrow = length(cor_matrices_sl_block1))
sigedges_pos_block2 = matrix(ncol = nrow(sigedges_pos_kn), nrow = length(cor_matrices_sl_block2))
sigedges_pos_block3 = matrix(ncol = nrow(sigedges_pos_kn), nrow = length(cor_matrices_sl_block3))
sigedges_pos_block4 = matrix(ncol = nrow(sigedges_pos_kn), nrow = length(cor_matrices_sl_block4))

sigedges_neg_block1 = matrix(ncol = nrow(sigedges_neg_kn), nrow = length(cor_matrices_sl_block1))
sigedges_neg_block2 = matrix(ncol = nrow(sigedges_neg_kn), nrow = length(cor_matrices_sl_block2))
sigedges_neg_block3 = matrix(ncol = nrow(sigedges_neg_kn), nrow = length(cor_matrices_sl_block3))
sigedges_neg_block4 = matrix(ncol = nrow(sigedges_neg_kn), nrow = length(cor_matrices_sl_block4))

# Loop through each subject
for (sub in seq_along(cor_matrices_sl_block1)) {
  cor_matrix1 = cor_matrices_sl_block1[[sub]]
  cor_matrix2 = cor_matrices_sl_block2[[sub]]
  cor_matrix3 = cor_matrices_sl_block3[[sub]]
  cor_matrix4 = cor_matrices_sl_block4[[sub]]
  # Loop through each row of sig_edges
  for (i in 1:nrow(sigedges_pos_kn)) {
    # Extract the node pairs from sig_edges
    node1_pos = sigedges_pos_kn[i,1]
    node2_pos = sigedges_pos_kn[i,2]
    # Extract the corresponding value from the correlation matrix and store it
    sigedges_pos_block1[sub,i] = cor_matrix1[node1_pos, node2_pos]
    sigedges_pos_block2[sub,i] = cor_matrix2[node1_pos, node2_pos]
    sigedges_pos_block3[sub,i] = cor_matrix3[node1_pos, node2_pos]
    sigedges_pos_block4[sub,i] = cor_matrix4[node1_pos, node2_pos]} 
  for (i in 1:nrow(sigedges_neg_kn)) {
    node1_neg = sigedges_neg_kn[i,1]
    node2_neg = sigedges_neg_kn[i,2]
    sigedges_neg_block1[sub,i] = cor_matrix1[node1_neg, node2_neg]
    sigedges_neg_block2[sub,i] = cor_matrix2[node1_neg, node2_neg]
    sigedges_neg_block3[sub,i] = cor_matrix3[node1_neg, node2_neg]
    sigedges_neg_block4[sub,i] = cor_matrix4[node1_neg, node2_neg]} }

# convert the result matrix to a data.frame
block1_pos = as.data.frame(sigedges_pos_block1)
block2_pos = as.data.frame(sigedges_pos_block2)
block3_pos = as.data.frame(sigedges_pos_block3)
block4_pos = as.data.frame(sigedges_pos_block4)
block1_neg = as.data.frame(sigedges_neg_block1)
block2_neg = as.data.frame(sigedges_neg_block2)
block3_neg = as.data.frame(sigedges_neg_block3)
block4_neg = as.data.frame(sigedges_neg_block4)

# set the row names to be the node pairs from sig_edges
colnames(block1_pos) = apply(sigedges_pos_kn, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block2_pos) = apply(sigedges_pos_kn, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block3_pos) = apply(sigedges_pos_kn, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block4_pos) = apply(sigedges_pos_kn, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block1_neg) = apply(sigedges_neg_kn, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block2_neg) = apply(sigedges_neg_kn, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block3_neg) = apply(sigedges_neg_kn, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block4_neg) = apply(sigedges_neg_kn, 1, function(row) paste0(row[1],"-",row[2]))

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
subid_sl = read.csv("./subid_sl.csv")
avgedge_allblocks_kn =
  as.data.frame(cbind(subid_sl,
                      block1_pos$pos_edgestr_block1, block2_pos$pos_edgestr_block2, block3_pos$pos_edgestr_block3, block4_pos$pos_edgestr_block4,
                      block1_neg$neg_edgestr_block1, block2_neg$neg_edgestr_block2, block3_neg$neg_edgestr_block3, block4_neg$neg_edgestr_block4), 
                stringsAsFactors = FALSE)
colnames(avgedge_allblocks_kn) = c("sub","pos_b1","pos_b2","pos_b3","pos_b4", "neg_b1","neg_b2","neg_b3","neg_b4")
avgedge_allblocks_kn = avgedge_allblocks_kn %>% mutate_at(c(1:9), as.numeric)

# extract the subject IDs from sub_data_sl and format them to match SL_cond
SL_cond = read.csv("./SL_stresscon.csv")
SL_cond$sub = as.numeric(substring(SL_cond$StressLearn_ID, 2,4))
trauma = read.csv("./trauma.csv")
motion = read.csv("./motion.csv")

# merge with avg edge data frame
avgedge_allblocks_kn = left_join(avgedge_allblocks_kn, SL_cond %>% select(sub, Group))
avgedge_allblocks_kn = left_join(avgedge_allblocks_kn, trauma)
avgedge_allblocks_kn = left_join(avgedge_allblocks_kn, motion)
avgedge_allblocks_kn$scan = ifelse(avgedge_allblocks_kn$sub > 55, 1, 0) # scanner type


##### apply SHAPES posmask to sl data #####
# load CPM sig edges (shapes)
pos_mask_shapes = read.csv("../net_masks/pos_mask_shapes_fi.csv", header=F)
neg_mask_shapes = read.csv("../net_masks/neg_mask_shapes_fi.csv", header=F)
# create matrix of just the upper triangle (symmetric matrix)
pos_edges_shapes = data.frame(which(pos_mask_shapes == 1 & upper.tri(pos_mask_shapes, diag = FALSE), arr.ind = TRUE))
neg_edges_shapes = data.frame(which(neg_mask_shapes == 1 & upper.tri(neg_mask_shapes, diag = FALSE), arr.ind = TRUE))
colnames(pos_edges_shapes) = c("col1","col2")
colnames(neg_edges_shapes) = c("col1","col2")

sigedges_pos_shapes = pos_edges_shapes %>% 
  mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
sigedges_neg_shapes = neg_edges_shapes %>% 
  mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])

## grab the significant edges of each matrices per condition per block
sigedges_pos_shapes_block1 = matrix(ncol = nrow(sigedges_pos_shapes), nrow = length(cor_matrices_sl_block1))
sigedges_pos_shapes_block2 = matrix(ncol = nrow(sigedges_pos_shapes), nrow = length(cor_matrices_sl_block2))
sigedges_pos_shapes_block3 = matrix(ncol = nrow(sigedges_pos_shapes), nrow = length(cor_matrices_sl_block3))
sigedges_pos_shapes_block4 = matrix(ncol = nrow(sigedges_pos_shapes), nrow = length(cor_matrices_sl_block4))

sigedges_neg_shapes_block1 = matrix(ncol = nrow(sigedges_neg_shapes), nrow = length(cor_matrices_sl_block1))
sigedges_neg_shapes_block2 = matrix(ncol = nrow(sigedges_neg_shapes), nrow = length(cor_matrices_sl_block2))
sigedges_neg_shapes_block3 = matrix(ncol = nrow(sigedges_neg_shapes), nrow = length(cor_matrices_sl_block3))
sigedges_neg_shapes_block4 = matrix(ncol = nrow(sigedges_neg_shapes), nrow = length(cor_matrices_sl_block4))

# loop through each subject
for (sub in seq_along(cor_matrices_sl_block1)) {
  cor_matrix1 = cor_matrices_sl_block1[[sub]]
  cor_matrix2 = cor_matrices_sl_block2[[sub]]
  cor_matrix3 = cor_matrices_sl_block3[[sub]]
  cor_matrix4 = cor_matrices_sl_block4[[sub]]
  # loop through each row of sig_edges
  for (i in 1:nrow(sigedges_pos_shapes)) {
    # extract the node pairs from sig_edges
    node1_pos = sigedges_pos_shapes[i,1]
    node2_pos = sigedges_pos_shapes[i,2]
    # extract the corresponding value from the correlation matrix and store it
    sigedges_pos_shapes_block1[sub,i] = cor_matrix1[node1_pos, node2_pos]
    sigedges_pos_shapes_block2[sub,i] = cor_matrix2[node1_pos, node2_pos]
    sigedges_pos_shapes_block3[sub,i] = cor_matrix3[node1_pos, node2_pos]
    sigedges_pos_shapes_block4[sub,i] = cor_matrix4[node1_pos, node2_pos]} 
  for (i in 1:nrow(sigedges_neg_shapes)) {
    node1_neg = sigedges_neg_shapes[i,1]
    node2_neg = sigedges_neg_shapes[i,2]
    sigedges_neg_shapes_block1[sub,i] = cor_matrix1[node1_neg, node2_neg]
    sigedges_neg_shapes_block2[sub,i] = cor_matrix2[node1_neg, node2_neg]
    sigedges_neg_shapes_block3[sub,i] = cor_matrix3[node1_neg, node2_neg]
    sigedges_neg_shapes_block4[sub,i] = cor_matrix4[node1_neg, node2_neg]} }

# convert the result matrix to a data.frame
block1_pos = as.data.frame(sigedges_pos_shapes_block1)
block2_pos = as.data.frame(sigedges_pos_shapes_block2)
block3_pos = as.data.frame(sigedges_pos_shapes_block3)
block4_pos = as.data.frame(sigedges_pos_shapes_block4)

block1_neg = as.data.frame(sigedges_neg_shapes_block1)
block2_neg = as.data.frame(sigedges_neg_shapes_block2)
block3_neg = as.data.frame(sigedges_neg_shapes_block3)
block4_neg = as.data.frame(sigedges_neg_shapes_block4)

# set the row names to be the node pairs from sig_edges
colnames(block1_pos) = apply(sigedges_pos_shapes, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block2_pos) = apply(sigedges_pos_shapes, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block3_pos) = apply(sigedges_pos_shapes, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block4_pos) = apply(sigedges_pos_shapes, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block1_neg) = apply(sigedges_neg_shapes, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block2_neg) = apply(sigedges_neg_shapes, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block3_neg) = apply(sigedges_neg_shapes, 1, function(row) paste0(row[1],"-",row[2]))
colnames(block4_neg) = apply(sigedges_neg_shapes, 1, function(row) paste0(row[1],"-",row[2]))

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
avgedge_allblocks_shapes =
  as.data.frame(cbind(subid_sl, block1_pos$pos_edgestr_block1, block2_pos$pos_edgestr_block2, block3_pos$pos_edgestr_block3, block4_pos$pos_edgestr_block4,
                      block1_neg$neg_edgestr_block1, block2_neg$neg_edgestr_block2, block3_neg$neg_edgestr_block3, block4_neg$neg_edgestr_block4), 
                stringsAsFactors = FALSE)
colnames(avgedge_allblocks_shapes) = c("sub", "pos_b1","pos_b2","pos_b3","pos_b4",
                                       "neg_b1","neg_b2","neg_b3","neg_b4")
avgedge_allblocks_shapes = avgedge_allblocks_shapes %>% mutate_at(c(1:9), as.numeric)

# merge with avg edge data frame
avgedge_allblocks_shapes = left_join(avgedge_allblocks_shapes, SL_cond %>% select(sub, Group))
avgedge_allblocks_shapes = left_join(avgedge_allblocks_shapes, trauma)
avgedge_allblocks_shapes = left_join(avgedge_allblocks_shapes, motion)
avgedge_allblocks_shapes$scan = ifelse(avgedge_allblocks_shapes$sub > 55, 1, 0) # scanner type


####################################
####################################
##### main sl analysis #####
t.test((avgedge_allblocks_tc %>% filter(Group=="Stress"))$pos_b2, (avgedge_allblocks_tc %>% filter(Group=="Control"))$pos_b2, var.equal=T)
t.test((avgedge_allblocks_kn %>% filter(Group=="Stress"))$pos_b2, (avgedge_allblocks_kn %>% filter(Group=="Control"))$pos_b2, var.equal=T)
t.test((avgedge_allblocks_shapes %>% filter(Group=="Stress"))$pos_b2, (avgedge_allblocks_shapes %>% filter(Group=="Control"))$pos_b2, var.equal=T)

summary(lm(pos_b2 ~ Group + motion + as.factor(sex) + as.factor(scan), avgedge_allblocks_tc))
summary(lm(pos_b2 ~ Group + motion + as.factor(sex) + as.factor(scan), avgedge_allblocks_kn))
summary(lm(pos_b2 ~ Group + motion + as.factor(sex) + as.factor(scan), avgedge_allblocks_shapes))

# create combined df
avgedge_allblocks_comb = bind_rows(
  avgedge_allblocks_tc %>% select(starts_with("pos_b"), Group) %>% mutate(sample = "tc"),
  avgedge_allblocks_kn %>% select(starts_with("pos_b"), Group) %>% mutate(sample = "kn"), 
  avgedge_allblocks_shapes %>% select(starts_with("pos_b"), Group) %>% mutate(sample = "shapes"))

summary(aov(pos_b1 ~ sample, avgedge_allblocks_comb))
summary(aov(pos_b2 ~ sample, avgedge_allblocks_comb))
summary(aov(pos_b3 ~ sample, avgedge_allblocks_comb))
summary(aov(pos_b4 ~ sample, avgedge_allblocks_comb))
# no sig differences between samples across all timepoints

# visualize side by side
avgedge_allblocks_comb %>%
  pivot_longer(cols = starts_with("pos_b"), names_to = "block", values_to = "value") %>%
  mutate(block = factor(block), sample = factor(sample, levels = c("tc", "kn", "shapes"))) %>%
  mutate(block = factor(block, levels = c("pos_b1", "pos_b2", "pos_b3", "pos_b4"),
                        labels = c("Run 1", "Run 2", "Run 3", "Run 4"))) %>%
  ggplot(aes(x = block, y = value, fill = sample)) + theme(legend.position=none) + 
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8, preserve = "single"), width = 0.8) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
    position = position_dodge(width = 0.8, preserve = "single"), width = 0.2) +
  scale_fill_manual(values = c("tc" = "#f28f3b", "kn" = "#fbd3aa", "shapes" = "#d1543f"), 
                    labels = c("tc" = "Original", "kn" = "Subsample", "shapes" = "Independent")) +
  labs(x = "Time", y = "Positive Network Mean FC", fill = "Sample") + theme_classic()

# plot them side by side #
avgedge_allblocks_comb %>%
  pivot_longer(cols = c(pos_b1, pos_b2, pos_b3, pos_b4),
    names_to = "run", values_to = "value") %>% 
  mutate(run = factor(run, levels = c("pos_b1", "pos_b2", "pos_b3", "pos_b4"),
                        labels = c("Run 1", "Run 2", "Run 3", "Run 4"))) %>%
  group_by(sample, Group, run) %>%
  summarise(mean_pos = mean(value, na.rm = TRUE),
    se_pos = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  mutate(sample_facet = factor(sample, levels = c("tc", "kn", "shapes"),
                          labels = c("Original", "Subsample", "Independent"))) %>%
  ggplot(aes(x = run, y = mean_pos, fill = interaction(Group, sample))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pos - se_pos, ymax = mean_pos + se_pos),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ sample_facet, scales = "free_y") +
  scale_fill_manual(
    values = c("Stress.tc" = "#f28f3b", "Control.tc" = "#588b8b",
      "Stress.kn" = "#fbd3aa", "Control.kn" = "#cbdad5",
      "Stress.shapes" = "#d1543f", "Control.shapes" = "#3d6464")) +
  labs(x = "Run", y = "Mean Network FC", fill = "Group") + theme_classic() + theme(legend.position="none")

# plot predictive power
behav_pred_tc = read.csv("./behav_pred_modage.csv", header=F)
behav_pred_tc = cbind(behav_pred_tc, sub = trauma$sub, trauma = trauma$total_5)
behav_pred_kn = read.csv("./cpm/net_masks/behav_pred_kn.csv", header=F)
trauma = read.csv("./trauma.csv")
trauma_kn = trauma %>% filter(!study=="SL")
behav_pred_kn = cbind(behav_pred_kn, sub = trauma_kn$sub, trauma = trauma_kn$total_5)
behav_pred_shapes = read.csv("./cpm/cpm_shapes/behav_pred_shapes_fi.csv", header=F)
distal = read.csv("./cpm/cpm_shapes/shapes_data_fi.csv")
behav_pred_shapes = cbind(behav_pred_shapes, sub = as.integer(distal$ucla_a_id), trauma = distal$distal_all)

pred = bind_rows(behav_pred_tc %>% mutate(sample = "tc"), behav_pred_kn %>% mutate(sample = "kn"), 
          behav_pred_shapes %>% mutate(sample = "shapes")) %>% 
  mutate(sample_facet = factor(sample, levels = c("tc", "kn", "shapes"),
                               labels = c("Original", "Subsample", "Independent"))) %>%
  ggplot(aes(x=trauma, y=V1)) + geom_point(alpha=0.4) + 
  geom_smooth(se=T, method="lm", color="black") + facet_wrap(~ sample_facet, scales = "free") +
  theme_classic() + labs(x="Predicted value", y="Observed value")
pred


### plot association with trauma 

# compute meanconn for TC
cor_matrices = readRDS("/gpfs/milgram/scratch60/gee_dylan/fah29/cor_matrices_tc_fin.rds")
weight_sigedges_pos_tc = data.frame(subject = 1:length(cor_matrices))
weight_sigedges_neg_tc = data.frame(subject = 1:length(cor_matrices))

# load CPM sig edges (TC)
pos_mask_tc = read.csv("./cpm/net_masks/pos_mask_modage.csv", header=F)
neg_mask_tc = read.csv("./cpm/net_masks/neg_mask_modage.csv", header=F)
# create matrix of just the upper triangle
pos_edges_tc = data.frame(which(pos_mask_tc == 1 & upper.tri(pos_mask_tc, diag = FALSE), arr.ind = TRUE))
neg_edges_tc = data.frame(which(neg_mask_tc == 1 & upper.tri(neg_mask_tc, diag = FALSE), arr.ind = TRUE))
colnames(pos_edges_tc) = c("col1","col2")
colnames(neg_edges_tc) = c("col1","col2")
sigedges_pos_tc = pos_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
sigedges_neg_tc = neg_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])

# loop over each subject and extract connectivity values for specified node pairs
for (sub in 1:length(cor_matrices)) {
  mat = cor_matrices[[sub]]
  for (i in 1:nrow(sigedges_pos_tc)) {
    node1 = sigedges_pos_tc$col1[i]
    node2 = sigedges_pos_tc$col2[i]
    column_name = paste0(sigedges_pos_tc$net1[i], "-", sigedges_pos_tc$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_pos_tc[sub, column_name] = connectivity_strength}
  for (i in 1:nrow(sigedges_neg_tc)) {
    node1 = sigedges_neg_tc$col1[i]
    node2 = sigedges_neg_tc$col2[i]
    column_name = paste0(sigedges_neg_tc$net1[i], "-", sigedges_neg_tc$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_neg_tc[sub, column_name] = connectivity_strength}}

# create summary sum and mean scores
weight_sigedges_pos_tc$meanconn_pos = rowMeans(weight_sigedges_pos_tc[, -1])
weight_sigedges_neg_tc$meanconn_neg = rowMeans(weight_sigedges_neg_tc[, -1])

# merge with trauma data 
trauma = cbind(trauma, weight_sigedges_pos_tc, weight_sigedges_neg_tc)
# sanity check
cor.test(trauma$total_5, trauma$meanconn_pos, method="spearman") # .47
cor.test(trauma$total_5, trauma$meanconn_neg, method="spearman") # -.55

### compute meanconn for KN
# initialize a data frame to store the results
cor_matrices_kn = readRDS("/gpfs/milgram/scratch60/gee_dylan/fah29/cor_matrices_kn_fin.rds")
weight_sigedges_pos_kn = data.frame(subject = 1:length(cor_matrices_kn))
weight_sigedges_neg_kn = data.frame(subject = 1:length(cor_matrices_kn))

# loop over each subject and extract connectivity values for specified node pairs
for (sub in 1:length(cor_matrices_kn)) {
  mat = cor_matrices_kn[[sub]]
  for (i in 1:nrow(sigedges_pos_kn)) {
    node1 = sigedges_pos_kn$col1[i]
    node2 = sigedges_pos_kn$col2[i]
    column_name = paste0(sigedges_pos_kn$net1[i], "-", sigedges_pos_kn$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_pos_kn[sub, column_name] = connectivity_strength}
  for (i in 1:nrow(sigedges_neg_kn)) {
    node1 = sigedges_neg_kn$col1[i]
    node2 = sigedges_neg_kn$col2[i]
    column_name = paste0(sigedges_neg_kn$net1[i], "-", sigedges_neg_kn$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_neg_kn[sub, column_name] = connectivity_strength}}

# create summary sum and mean scores
weight_sigedges_pos_kn$meanconn_pos = rowMeans(weight_sigedges_pos_kn[, -1])
weight_sigedges_neg_kn$meanconn_neg = rowMeans(weight_sigedges_neg_kn[, -1])

# merge with trauma data 
trauma_kn = cbind(trauma_kn, weight_sigedges_pos_kn, weight_sigedges_neg_kn)
names(trauma_kn) = make.names(names(trauma_kn), unique = TRUE)

# check correlations
cor.test(trauma_kn$total_5, trauma_kn$meanconn_pos, method="spearman")
cor.test(trauma_kn$total_5, trauma_kn$meanconn_neg, method="spearman")


### compute meanconn for shapes
cor_matrices_shapes = readRDS("/gpfs/milgram/scratch60/gee_dylan/fah29/cor_matrices_shapes.rds")

weight_sigedges_pos_shapes = data.frame(subject = 1:length(cor_matrices_shapes))
weight_sigedges_neg_shapes = data.frame(subject = 1:length(cor_matrices_shapes))

# loop over each subject and extract connectivity values for specified node pairs
for (sub in 1:length(cor_matrices_shapes)) {
  mat = cor_matrices_shapes[[sub]]
  for (i in 1:nrow(sigedges_pos_shapes)) {
    node1 = sigedges_pos_shapes$col1[i]
    node2 = sigedges_pos_shapes$col2[i]
    column_name = paste0(sigedges_pos_shapes$net1[i], "-", sigedges_pos_shapes$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_pos_shapes[sub, column_name] = connectivity_strength}
  for (i in 1:nrow(sigedges_neg_shapes)) {
    node1 = sigedges_neg_shapes$col1[i]
    node2 = sigedges_neg_shapes$col2[i]
    column_name = paste0(sigedges_neg_shapes$net1[i], "-", sigedges_neg_shapes$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_neg_shapes[sub, column_name] = connectivity_strength}}

# create summary sum and mean scores
weight_sigedges_pos_shapes$meanconn_pos = rowMeans(weight_sigedges_pos_shapes[, -1])
weight_sigedges_neg_shapes$meanconn_neg = rowMeans(weight_sigedges_neg_shapes[, -1])

# load distal data for shapes
distal = read.csv("./cpm/cpm_shapes/shapes_data_fi.csv")
distal = cbind(distal, weight_sigedges_pos_shapes, weight_sigedges_neg_shapes)
# check correlations
cor.test(distal$distal_all, distal$meanconn_pos, method="spearman") # .53
cor.test(distal$distal_all, distal$meanconn_neg, method="spearman") # -.52

# visualize plot side by side 
cortrauma = bind_rows(
  data.frame(meanconn_pos = trauma$meanconn_pos, trauma = trauma$total_5, sample = "tc"),
  data.frame(meanconn_pos = trauma_kn$meanconn_pos, trauma = trauma_kn$total_5, sample = "kn"),
  data.frame(meanconn_pos = distal$meanconn_pos, trauma = distal$distal_all, sample = "shapes")) %>%
  mutate(sample_facet = factor(sample, levels = c("tc", "kn", "shapes"),
                          labels = c("Original", "Subsample", "Independent"))) %>%
  ggplot(aes(y = meanconn_pos, x = trauma)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#0064dd", fill = "#0064dd") + theme_classic() + 
  facet_wrap(~ sample_facet, scales = "free") + labs(y = "Positive Network Mean FC", x = "Total traumatic life events")

ggarrange(pred, cortrauma, nrow=2)
