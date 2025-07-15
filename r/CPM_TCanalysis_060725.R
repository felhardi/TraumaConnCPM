#### load cpm masks ####
pos_mask_tc = read.csv("../net_masks/pos_mask_modage.csv", header=F)
neg_mask_tc = read.csv("../net_masks/neg_mask_modage.csv", header=F)

# create matrix of just the upper triangle
pos_edges_tc = data.frame(which(pos_mask_tc == 1 & upper.tri(pos_mask_tc, diag = FALSE), arr.ind = TRUE))
neg_edges_tc = data.frame(which(neg_mask_tc == 1 & upper.tri(neg_mask_tc, diag = FALSE), arr.ind = TRUE))
colnames(pos_edges_tc) = c("col1","col2")
colnames(neg_edges_tc) = c("col1","col2")

# load labels 
labels = read.csv("./xilin_liz_combined.csv")
# put labels with sig edges
net_map = setNames(labels$net_names, labels$node)
roi_map = setNames(labels$BA_othername, labels$node)

# apply to sigedges files
sigedges_pos_tc = pos_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
sigedges_neg_tc = neg_edges_tc %>% mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])

# initialize a data frame to store the results
weight_sigedges_pos = data.frame(subject = 1:length(cor_matrices_tc_fin))
weight_sigedges_neg = data.frame(subject = 1:length(cor_matrices_tc_fin))

# loop over each subject and extract connectivity values for specified node pairs
for (sub in 1:length(cor_matrices_tc_fin)) {
  mat = cor_matrices_tc_fin[[sub]]
  for (i in 1:nrow(sigedges_pos_tc)) {
    node1 = sigedges_pos_tc$col1[i]
    node2 = sigedges_pos_tc$col2[i]
    column_name = paste0(sigedges_pos_tc$net1[i], "-", sigedges_pos_tc$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_pos[sub, column_name] = connectivity_strength}
  for (i in 1:nrow(sigedges_neg_tc)) {
    node1 = sigedges_neg_tc$col1[i]
    node2 = sigedges_neg_tc$col2[i]
    column_name = paste0(sigedges_neg_tc$net1[i], "-", sigedges_neg_tc$net2[i], "-", i)
    connectivity_strength = mat[node1, node2]
    weight_sigedges_neg[sub, column_name] = connectivity_strength}}

# create summary mean scores
weight_sigedges_pos$meanconn_pos = rowMeans(weight_sigedges_pos[, -1])
weight_sigedges_neg$meanconn_neg = rowMeans(weight_sigedges_neg[, -1])

# merge with trauma data 
trauma = read.csv("./trauma.csv")
trauma = cbind(trauma, weight_sigedges_pos, weight_sigedges_neg)


#### predictive power ####

# plot predicted value and R observed value
behav_pred = read.csv("../net_masks/behav_pred_modage.csv", header=F)
behav_pred = cbind(behav_pred, trauma %>% select(sub, total_5))

# sanity check: test predictive power
cor.test(behav_pred$V1, behav_pred$total_5, method="spearman") # .41

# plot predictive power
behav_pred %>% 
  ggplot(aes(x=total_5, y=V1)) + geom_point(alpha=0.4) + 
  geom_smooth(se=T, method="lm", color="black") + 
  theme_classic() + labs(x="Predicted value", y="Observed value")


##### check correlations with total trauma #####
cor.test(trauma$total_5, trauma$meanconn_pos, method="spearman") # .47
cor.test(trauma$total_5, trauma$meanconn_neg, method="spearman") # -.55

# plot association with traumatic life events
ggarrange((trauma %>% select(meanconn_pos, total_5) %>%
             ggplot(aes(x=total_5, y=meanconn_pos)) + geom_point(alpha=0.4) + 
             geom_smooth(se=T, method="lm", color="#0064dd", fill="#0064dd") + 
             theme_classic() + labs(x="Total traumatic life events", y="Positive Network Mean FC")),
          (trauma %>% select(meanconn_neg, total_5) %>%
             ggplot(aes(x=total_5, y=meanconn_neg)) + geom_point(alpha=0.4) + 
             geom_smooth(se=T, method="lm", color="#bf1212", fill="#bf1212") + 
             theme_classic() + labs(x="Total traumatic life events", y="Negative Network Mean FC")), ncol=1, nrow=2)


##### check correlations with other measures ####

## PSS ##
pss = read.csv("./pss.csv")
pss = cbind(pss, weight_sigedges_pos %>% select(meanconn_pos), weight_sigedges_neg %>% select(meanconn_neg))

# check correlation
cor.test(pss$pss, pss$meanconn_pos) # .06; n=119
cor.test(pss$pss, pss$meanconn_neg) # .07; n=119

# grab subid of people with complete pss data 
pss_subid = pss %>% filter(is.na(pss)==F & is.na(meanconn_pos)==F) %>% select(sub)

# create plots
pss_pos = pss %>% 
  ggplot(aes(x = pss, y = meanconn_pos)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#0064dd", fill = "#0064dd") + theme_classic() + 
  labs(y = "Positive Network Mean FC", x = "PSS")
pss_neg = pss %>% 
  ggplot(aes(x = pss, y = meanconn_neg)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#bf1212", fill = "#bf1212") + theme_classic() + 
  labs(y = "Negative Network Mean FC", x = "PSS")

### get AUDIT data ###
audit = read.csv("./audit.csv")
audit = cbind(audit, weight_sigedges_pos %>% select(meanconn_pos), weight_sigedges_neg %>% select(meanconn_neg))
audit = left_join(pss_subid, audit)

cor.test(audit$audit_mean, audit$meanconn_pos) # .01; n=119
cor.test(audit$audit_mean, audit$meanconn_neg) # -0.07; n=119

aud_pos = audit %>% 
  ggplot(aes(x = audit_mean, y = meanconn_pos)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#0064dd", fill = "#0064dd") + theme_classic() + 
  labs(y = "Positive Network Mean FC", x = "AUDIT")
aud_neg = audit %>% 
  ggplot(aes(x = audit_mean, y = meanconn_neg)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#bf1212", fill = "#bf1212") + theme_classic() + 
  labs(y = "Negative Network Mean FC", x = "AUDIT")

# plot trauma scores as a comparison (n = 119 with available data)
trauma_pss = left_join(pss_subid, trauma %>% select(sub, meanconn_pos, meanconn_neg, total_5))
cor.test(trauma_pss$total_5, trauma_pss$meanconn_pos) # .32
cor.test(trauma_pss$total_5, trauma_pss$meanconn_neg) # -.33

trauma_pos = trauma_pss %>% 
  ggplot(aes(x = total_5, y = meanconn_pos)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#0064dd", fill = "#0064dd") + theme_classic() + 
  labs(y = "Positive Network Mean FC", x = "Total traumatic life events")
trauma_neg = trauma_pss %>% 
  ggplot(aes(x = total_5, y = meanconn_neg)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#bf1212", fill = "#bf1212") + theme_classic() + 
  labs(y = "Negative Network Mean FC", x = "Total traumatic life events")

# edu
trauma_pss = trauma_pss %>% left_join(trauma %>% select(sub, edu_num))
edu_pos = trauma_pss %>% 
  ggplot(aes(x = edu_num, y = meanconn_pos)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#0064dd", fill = "#0064dd") + theme_classic() + 
  labs(y = "Positive Network Mean FC", x = "Educational attainment")
edu_neg = trauma_pss %>% 
  ggplot(aes(x = edu_num, y = meanconn_neg)) + geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "#bf1212", fill = "#bf1212") + theme_classic() + 
  labs(y = "Negative Network Mean FC", x = "Educational attainment")
cor.test(trauma_pss$meanconn_pos, trauma_pss$edu_num) # -0.004
cor.test(trauma_pss$meanconn_neg, trauma_pss$edu_num) # 0.05

# create a combined plot for the supp
ggarrange(trauma_pos, edu_pos, aud_pos, pss_pos, trauma_neg, edu_neg, aud_neg, pss_neg, nrow=2, ncol=4)


### node overlap across positive and negative networks ###
pos_nodeg = bind_rows(select(sigedges_pos_tc, node = col1), select(sigedges_pos_tc, node = col2)) %>% 
  count(node, name = "degree_tc")
neg_nodeg = bind_rows(select(sigedges_neg_tc, node = col1), select(sigedges_neg_tc, node = col2)) %>% 
  count(node, name = "degree_tc")
node_posneg = intersect(pos_nodeg[1], neg_nodeg[1])
node_posneg = left_join(node_posneg, labels)

node_posdeg = node_posneg %>% left_join(select(sigedges_pos_tc, node = col1))
node_posdeg = node_posdeg %>% left_join(select(sigedges_pos_tc, node = col2))
node_posdeg = node_posdeg %>% group_by(node) %>% count()

node_negdeg = node_posneg %>% left_join(select(sigedges_neg_tc, node = col1))
node_negdeg = node_negdeg %>% left_join(select(sigedges_neg_tc, node = col2))
node_negdeg = node_negdeg %>% group_by(node) %>% count()
