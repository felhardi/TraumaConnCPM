######## load pos_mask and neg_mask
pos_mask_all = readMat('./cpm/virtual_lesion/pos_mask_all.mat')$pos.mask.all
neg_mask_all = readMat('./cpm/virtual_lesion/neg_mask_all.mat')$neg.mask.all

pos_mask_list = vector("list", 170)
neg_mask_list = vector("list", 170)
for (i in 1:170) {
  pos_mask_list[[i]] = pos_mask_all[,,i]
  neg_mask_list[[i]] = neg_mask_all[,,i]}

# load labels
labels = read.csv("./misc/xilin_liz_combined.csv")
networks = sort(unique(labels$net_names))

# create combinations of networks (including self-combinations)
combs = combn(networks, 2, simplify = FALSE)
combs = c(combs, lapply(networks, function(x) c(x, x)))

# put labels with sig edges
net_map = setNames(labels$net_names, labels$node)
roi_map = setNames(labels$BA_othername, labels$node)

# create matrix of just the upper triangle
pos_edges = which(pos_mask == 1 & upper.tri(pos_mask, diag = FALSE), arr.ind = TRUE)
neg_edges = which(neg_mask == 1 & upper.tri(neg_mask, diag = FALSE), arr.ind = TRUE)

sigedges_pos = data.frame(pos_edges) %>% mutate(roi1 = roi_map[row], roi2 = roi_map[col],
                                                net1 = net_map[row], net2 = net_map[col])
sigedges_neg = data.frame(neg_edges) %>% mutate(roi1 = roi_map[row], roi2 = roi_map[col],
                                                net1 = net_map[row], net2 = net_map[col])

# create a matrix for each unique net-net pair
for(i in 1:length(combs)) {
  comb <- combs[[i]]
  net1 <- comb[1]
  net2 <- comb[2]
  # Initialize a matrix of zeros
  pos_mask <- matrix(0, nrow = 377, ncol = 377)
  # Filter the dataframe for the current combination of net1 and net2
  selected_rows <- sigedges_pos[(sigedges_pos$net1 == net1 & sigedges_pos$net2 == net2) | 
                                  (sigedges_pos$net1 == net2 & sigedges_pos$net2 == net1), ]
  # Populate the appropriate cells with 1s
  for(j in 1:nrow(selected_rows)) {
    pos_mask[selected_rows$row[j], selected_rows$col[j]] <- 1}
  # Assign the matrix to a variable in the format pos_mask_## (e.g., pos_mask_01)
  assign(paste0('pos_mask_', sprintf('%02d', i)), pos_mask)}
for(i in 1:length(combs)) {
  comb <- combs[[i]]
  net1 <- comb[1]
  net2 <- comb[2]
  # Initialize a matrix of zeros
  neg_mask <- matrix(0, nrow = 377, ncol = 377)
  # Filter the dataframe for the current combination of net1 and net2
  selected_rows <- sigedges_neg[(sigedges_neg$net1 == net1 & sigedges_neg$net2 == net2) | 
                                  (sigedges_neg$net1 == net2 & sigedges_neg$net2 == net1), ]
  # Populate the appropriate cells with 1s
  for(j in 1:nrow(selected_rows)) {
    neg_mask[selected_rows$row[j], selected_rows$col[j]] <- 1}
  # Assign the matrix to a variable in the format pos_mask_## (e.g., pos_mask_01)
  assign(paste0('neg_mask_', sprintf('%02d', i)), neg_mask)}

# check which one contains all zero and exclude
# List to store the names of matrices that contain all zeros, save out the rest
zero_matrices_pos <- list()
library(R.matlab)
# Check each matrix
for(i in 1:length(combs)) {
  # Get the matrix name
  matrix_name <- paste0('pos_mask_', sprintf('%02d', i))
  # Retrieve the matrix
  matrix <- get(matrix_name)
  # Check if all elements are zero
  if(all(matrix == 0)) {
    zero_matrices_pos <- c(zero_matrices_pos, matrix_name)
  } else {
    # Save the non-zero matrix to a .mat file
    writeMat(con = paste0(matrix_name, ".mat"), matrix = matrix)}}
# check which one contains all zero and exclude
# List to store the names of matrices that contain all zeros, save out the rest
zero_matrices_neg <- list()
# Check each matrix
for(i in 1:length(combs)) {
  # Get the matrix name
  matrix_name <- paste0('neg_mask_', sprintf('%02d', i))
  # Retrieve the matrix
  matrix <- get(matrix_name)
  # Check if all elements are zero
  if(all(matrix == 0)) {
    zero_matrices_neg <- c(zero_matrices_neg, matrix_name)
  } else {
    # Save the non-zero matrix to a .mat file
    writeMat(con = paste0(matrix_name, ".mat"), matrix = matrix)}}

######## load up results
combs_list <- do.call(rbind, combs)
combs_df <- data.frame(combs_list, stringsAsFactors = FALSE)
colnames(combs_df) <- c("net1", "net2")

# load all the matrices
emptymat_pos = data.frame(Mask_Index=do.call(rbind, zero_matrices_pos))
emptymat_pos$Mask = as.numeric(gsub("pos_mask_", "", emptymat_pos$Mask_Index))
emptymat_pos$empty = 1
emptymat_neg = data.frame(Mask_Index=do.call(rbind, zero_matrices_neg))
emptymat_neg$Mask = as.numeric(gsub("neg_mask_", "", emptymat_neg$Mask_Index))
emptymat_neg$empty = 1

posmask_results = read.csv("./cpm/virtual_lesion/posmask_results_020325.csv")
colnames(posmask_results) = c("Mask", "R_pos", "P_pos")
posmask_results %>% filter(R_pos < 0) %>% count() 
posmask_results = left_join(posmask_results, emptymat_pos)
posmask_results$R_pos = ifelse(is.na(posmask_results$empty)==F, NA, posmask_results$R_pos)
posmask_results$R_pos = ifelse(posmask_results$R_pos < 0, 0, posmask_results$R_pos)
posmask_results = cbind(posmask_results, combs_df)

negmask_results = read.csv("./cpm/virtual_lesion/negmask_results_020425.csv")
negmask_results %>% filter(R_neg < 0) %>% count() 
negmask_results = left_join(negmask_results, emptymat_neg)
negmask_results$R_neg = ifelse(is.na(negmask_results$empty)==F, NA, negmask_results$R_neg)
negmask_results$R_neg = ifelse(negmask_results$R_neg < 0, NA, negmask_results$R_neg)
negmask_results = cbind(negmask_results, combs_df)

mask_results = cbind(posmask_results, negmask_results %>% select(R_neg, P_neg))

vl_results = mask_results
library(psych)
vl_results$r.pos_z = fisherz(vl_results$R_pos)
vl_results$r.neg_z = fisherz(vl_results$R_neg)

vl_results = vl_results %>% mutate(net1 = if_else(net1 == "Cerebellum", "Cer", net1), 
                                   net2 = if_else(net2 == "Cerebellum", "Cer", net2))

# plot
net1 = vl_results$net1
net2 = vl_results$net2
vl_pos = 
  data.frame(rbind(cbind(net1, net2), cbind(net2, net1)), r.pos_z=vl_results$r.pos_z) %>% rename("Network"=net1) %>% rename("Network2"=net2) %>%
  filter(!is.na(r.pos_z)) %>% mutate(Network = if_else(Network == "Cerebellum", "Cer", Network)) %>%
  ggplot(aes(x = fct_reorder(Network, r.pos_z, .fun = mean, .desc = T), y = r.pos_z, fill = Network)) + 
  geom_boxplot() + scale_fill_viridis_d() + theme_classic() +
  geom_point(aes(color=Network2), alpha = 0.5, position = position_jitter(width = 0.2, height = 0)) + 
  scale_color_viridis_d() + geom_hline(yintercept = 0.18, linetype = "dashed", color = "grey30") +
  labs(y = "Predicted x Observed R (Pos)", x = "Included Network Connections") + theme(legend.position = "none")

vl_neg =   
  data.frame(rbind(cbind(net1, net2), cbind(net2, net1)), r.neg_z=vl_results$r.neg_z) %>% rename("Network"=net1) %>% rename("Network2"=net2) %>%
  filter(!is.na(r.neg_z)) %>% mutate(Network = if_else(Network == "Cerebellum", "Cer", Network)) %>%
  ggplot(aes(x = fct_reorder(Network, r.neg_z, .fun = mean, .desc = TRUE), y = r.neg_z, fill = Network)) + 
  geom_boxplot() + scale_fill_viridis_d() + theme_classic() +
  geom_point(aes(color=Network2), alpha = 0.5, position = position_jitter(width = 0.2, height = 0)) + 
  scale_color_viridis_d() + geom_hline(yintercept = 0.19, linetype = "dashed", color = "grey30") +
  labs(y = "Predicted x Observed R (Neg)", x = "Included Network Connections") + theme(legend.position = "none")

library(ggpubr)
ggarrange(vl_pos, vl_neg, ncol=2)


### test if these are statistically diff from null (actual)
null_pos = 0.1813
null_neg = 0.1874
vl_pos_long = vl_results %>%
  select(net1, net2, r.pos_z) %>% 
  pivot_longer(cols = c(starts_with("net")), names_to = "net_type", values_to = "net") %>%
  filter(!is.na(r.pos_z))
vl_neg_long = vl_results %>%
  select(net1, net2, r.neg_z) %>% 
  pivot_longer(cols = c(starts_with("net")), names_to = "net_type", values_to = "net") %>%
  filter(!is.na(r.neg_z))

library(broom)
vl_pos_long %>%
  group_by(net) %>%
  summarise(t_test = list(t.test(r.pos_z, mu = null_pos) %>% tidy())) %>%
  unnest(t_test) # MF, Motor, SN, Cer, DMN do not differ from actual
vl_neg_long %>%
  group_by(net) %>%
  summarise(t_test = list(t.test(r.neg_z, mu = null_neg) %>% tidy())) %>%
  unnest(t_test) # none of them differs from actual 
