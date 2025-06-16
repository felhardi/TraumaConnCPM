# mask
pos_mask_tc = read.csv("../net_masks/pos_mask_modage.csv", header=F)
neg_mask_tc = read.csv("../net_masks/neg_mask_modage.csv", header=F)

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
sigedges_pos_tc = pos_edges_tc %>% 
  mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])
sigedges_neg_tc = neg_edges_tc %>% 
  mutate(roi1 = roi_map[col1], roi2 = roi_map[col2], net1 = net_map[col1], net2 = net_map[col2])

# create node-to-node edge count df for each network
sigedges_pos_tc$pair = apply(sigedges_pos_tc, 1, function(x) {
  paste(sort(x[5:6]), collapse = "-")})
sigedges_neg_tc$pair = apply(sigedges_neg_tc, 1, function(x) {
  paste(sort(x[5:6]), collapse = "-")})

pos_edge_count = as.data.frame(table(sigedges_pos_tc$pair))
pos_split_pairs = strsplit(as.character(pos_edge_count$Var1), "-")
pos_edge_pairs = data.frame(do.call(rbind, pos_split_pairs))
colnames(pos_edge_pairs) = c("Node1", "Node2")
pos_edge_pairs$Count = pos_edge_count$Freq # add counts to edge pairs
pos_edge_pairs = rbind(c("Cerebellum", "Cerebellum", NA), c("Cerebellum", "DMN", NA), 
                       c("Cerebellum", "Cerebellum", NA), pos_edge_pairs) # add empty rows/cols
pos_edge_pairs = rbind(pos_edge_pairs, c("VI", "VI", NA), c("VII", "VII", NA))

net_count = labels %>% count(net_names) # number of nodes within each network
pos_edge_pairs = pos_edge_pairs %>%
  left_join(net_count, by = c("Node1" = "net_names"), relationship="many-to-many") %>% rename(node1_count = n)
pos_edge_pairs = pos_edge_pairs %>%
  left_join(net_count, by = c("Node2" = "net_names"), relationship="many-to-many") %>% rename(node2_count = n)
pos_edge_pairs$totedge = pos_edge_pairs$node1_count * pos_edge_pairs$node2_count
pos_edge_pairs$Count = as.numeric(pos_edge_pairs$Count)
pos_edge_pairs$Count_prop = as.numeric(pos_edge_pairs$Count) / pos_edge_pairs$totedge * 100
pos_edge_pairs$Count_prop_adj = pos_edge_pairs$Count_prop / 90 * 100
pos_edge_pairs$Count_adj = pos_edge_pairs$Count / 90 * 100

neg_edge_count = as.data.frame(table(sigedges_neg_tc$pair))
neg_split_pairs = strsplit(as.character(neg_edge_count$Var1), "-")
neg_edge_pairs = data.frame(do.call(rbind, neg_split_pairs))
colnames(neg_edge_pairs) = c("Node1", "Node2")
neg_edge_pairs$Count = neg_edge_count$Freq # add counts to edge pairs
neg_edge_pairs = rbind(neg_edge_pairs, c("VII", "VII", NA))
neg_edge_pairs = neg_edge_pairs %>%
  left_join(net_count, by = c("Node1" = "net_names"), relationship="many-to-many") %>% rename(node1_count = n)
neg_edge_pairs = neg_edge_pairs %>%
  left_join(net_count, by = c("Node2" = "net_names"), relationship="many-to-many") %>% rename(node2_count = n)
neg_edge_pairs$totedge = neg_edge_pairs$node1_count * neg_edge_pairs$node2_count
neg_edge_pairs$Count = as.numeric(neg_edge_pairs$Count)
neg_edge_pairs$Count_prop = as.numeric(neg_edge_pairs$Count) / neg_edge_pairs$totedge * 100
neg_edge_pairs$Count_prop_adj = neg_edge_pairs$Count_prop / 414 * 100
neg_edge_pairs$Count_adj = neg_edge_pairs$Count / 414 * 100

# create summary edge count df
count_sigedges_roi_pos = sigedges_pos_tc %>% count(col1) %>% rename(node = col1, n1_pos = n) %>%
  full_join(sigedges_pos_tc %>% count(col2) %>% rename(node = col2, n2_pos = n))
count_sigedges_roi_pos$roi_sum_pos = rowSums(count_sigedges_roi_pos[2:3], na.rm=T)
count_sigedges_roi_pos = full_join(labels %>% select(node), count_sigedges_roi_pos %>% select(node, roi_sum_pos))
count_sigedges_roi_pos$roi_sum_pos = ifelse(is.na(count_sigedges_roi_pos$roi_sum_pos)==T, 0, count_sigedges_roi_pos$roi_sum_pos)
count_sigedges_roi_neg = sigedges_neg_tc %>% count(col1) %>% rename(node = col1, n1_neg = n) %>%
  full_join(sigedges_neg_tc %>% count(col2) %>% rename(node = col2, n2_neg = n))
count_sigedges_roi_neg$roi_sum_neg = rowSums(count_sigedges_roi_neg[2:3], na.rm=T)
count_sigedges_roi_neg = full_join(labels %>% select(node), count_sigedges_roi_neg %>% select(node, roi_sum_neg))
count_sigedges_roi_neg$roi_sum_neg = ifelse(is.na(count_sigedges_roi_neg$roi_sum_neg)==T, 0, count_sigedges_roi_neg$roi_sum_neg)

count_sigedges_pos = sigedges_pos_tc %>% count(net1) %>% rename(net = net1, n1_pos = n) %>%
  full_join(sigedges_pos_tc %>% count(net2) %>% rename(net = net2, n2_pos = n))
count_sigedges_pos$net_sum_pos = rowSums(count_sigedges_pos[2:3], na.rm=T)
count_sigedges_neg = sigedges_neg_tc %>% count(net1) %>% rename(net = net1, n1_neg = n) %>%
  full_join(sigedges_neg_tc %>% count(net2) %>% rename(net = net2, n2_neg = n))
count_sigedges_neg$net_sum_neg = rowSums(count_sigedges_neg[2:3], na.rm=T)

count_sigedges = left_join(count_sigedges_pos[c(1, 4)], count_sigedges_neg[c(1, 4)])

# compute the number of edges per network
# Cerebellum 17; DMN 32; FP 74; MF 55; Motor 60; SC 12; SN 52; VI 30; VII 13; Vas 32
count_sigedges$edge_tot = c(((17 * (17-1)) / 2) + (17 * (377-17)), 
                            ((32 * (32-1)) / 2) + (32 * (377-32)),
                            ((74 * (74-1)) / 2) + (74 * (377-74)),
                            ((55 * (55-1)) / 2) + (55 * (377-55)),
                            ((60 * (60-1)) / 2) + (60 * (377-60)),
                            ((12 * (12-1)) / 2) + (12 * (377-12)),
                            ((52 * (52-1)) / 2) + (52 * (377-52)),
                            ((30 * (30-1)) / 2) + (30 * (377-30)),
                            ((13 * (13-1)) / 2) + (13 * (377-13)),
                            ((32 * (32-1)) / 2) + (32 * (377-32)))

# compute the proportion (density - edge sum / all possible edge)
count_sigedges$net_prop_pos = count_sigedges$net_sum_pos / count_sigedges$edge_tot * 100
count_sigedges$net_prop_neg = count_sigedges$net_sum_neg / count_sigedges$edge_tot * 100

# adjust for number of network edge (network edge proportionate to network density)
count_sigedges$net_prop_pos_adj = count_sigedges$net_prop_pos / 90 * 100
count_sigedges$net_prop_neg_adj = count_sigedges$net_prop_neg / 414 * 100

# visualize
count_sigedges %>% ### this is the final one that we settled with
  select(net, net_prop_pos_adj, net_prop_neg_adj) %>% 
  pivot_longer(
    cols = c(net_prop_pos_adj, net_prop_neg_adj),
    names_to = c(".value", "Network"),
    names_pattern = "(.+)_(.+)_adj") %>%
  mutate(Network = dplyr::recode(Network, "neg" = "Negative", "pos" = "Positive")) %>% 
  mutate(net = dplyr::recode(net, "Cerebellum" = "Cer")) %>% 
  ggplot(aes(x = net, y = net_prop, fill = net)) + 
  geom_bar(aes(alpha = Network), stat = "identity", position = "dodge") + 
  theme_classic() + 
  ylab("Network edge count/degree (adjusted)") + xlab("Network") +
  scale_fill_viridis_d(option = "D", name = "Subnetwork") + 
  scale_alpha_manual(values = c(0.5, 1))
