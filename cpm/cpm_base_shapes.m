clear;
clc;

% ------------ INPUTS -------------------

load('mat_shapes_sh.mat');
load('shapes_data.mat');

matrix = mat_shapes;
tc = shapes_select;

% fisher z-transform the matrices
tolerance = 1e-12; % to avoid log(0)
matrix(matrix >= 1 - tolerance) = 1 - tolerance;
matrix(matrix <= -1 + tolerance) = -1 + tolerance;

z_mats = 0.5 * log((1 + matrix) ./ (1 - matrix));

all_mats = z_mats;
all_behav = tc.distal_all;
age = tc.ucla_a_age;

% threshold for feature selection
thresh = 0.01;

% ---------------------------------------

no_sub = size(all_mats,3);
no_node = size(all_mats, 1);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);
behav_pred = zeros(no_sub,1);

r_mat_all = zeros(no_node, no_node, no_sub);
p_mat_all = zeros(no_node, no_node, no_sub);

pos_mask_all = zeros(no_node, no_node, no_sub);
neg_mask_all = zeros(no_node, no_node, no_sub);

% to store correlations at each iteration
R_pos_iter = zeros(no_sub, 1);
P_pos_iter = zeros(no_sub, 1);
R_neg_iter = zeros(no_sub, 1);
P_neg_iter = zeros(no_sub, 1);

% leave one out train/test

for leftout = 1:no_sub

    fprintf('\n Leaving out subj # %6.3f',leftout);
    
    % leave out subject from matrices and behavior
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats,[],size(train_mats,3));

    train_behav = all_behav;
    train_behav(leftout) = [];
    train_age = double(age);
    train_age(leftout) = [];

    % partial correlation
    [r_mat, p_mat] = partialcorr(train_vcts', train_behav, train_age, "Type", "Spearman");

    r_mat = reshape(r_mat,no_node,no_node);
    p_mat = reshape(p_mat,no_node,no_node);
    
    % set threshold and define masks into pos/neg networks
    pos_mask = zeros(no_node, no_node);
    neg_mask = zeros(no_node, no_node);
    
    % find the significant edges
    pos_edge = find(r_mat > 0 & p_mat < thresh);
    neg_edge = find(r_mat < 0 & p_mat < thresh);
    
    pos_mask(pos_edge) = 1;
    neg_mask(neg_edge) = 1;

    % save out pos_mask and neg_mask for this iteration
    pos_mask_all(:, :, leftout) = pos_mask;
    neg_mask_all(:, :, leftout) = neg_mask;

    % get sum of all edges in TRAIN subs 
    % (divide by 2 to control for the fact that matrices are symmetric)
    train_sumpos = zeros(no_sub-1,1);
    train_sumneg = zeros(no_sub-1,1);
    
    for ss = 1:size(train_sumpos)
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
    end
    
    % build model on TRAIN subs
    fit_pos = regress(train_behav, [train_sumpos, ones(no_sub-1,1)]);
    fit_neg = regress(train_behav, [train_sumneg, ones(no_sub-1,1)]);
    % combining both positive and negative features
    b = regress(train_behav, [train_sumpos, train_sumneg, ones(no_sub-1,1)]);

    % run model on TEST sub    
    test_mat = all_mats(:,:,leftout);
    test_sumpos = sum(sum(test_mat.*pos_mask))/2;
    test_sumneg = sum(sum(test_mat.*neg_mask))/2;
    
    behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
    behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
    behav_pred(leftout) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);
    
end

% compare predicted and observed scores
include_idx = setdiff(1:no_sub, 13); # NA pred

[R_pos, P_pos] = corr(behav_pred_pos(include_idx), all_behav(include_idx), "type", "Spearman")
[R_neg, P_neg] = corr(behav_pred_neg(include_idx), all_behav(include_idx), "type", "Spearman")
[R_comb, P_comb] = corr(behav_pred(include_idx), all_behav(include_idx), "type", "Spearman")

% find edges that were consistently significant across all iterations
same_pos_mask = all(pos_mask_all == 1, 3);
same_neg_mask = all(neg_mask_all == 1, 3);

% find indices for the upper triangular part and combine
[p_row, p_col] = find(triu(same_pos_mask, 1));  % exclude diagonal
[n_row, n_col] = find(triu(same_neg_mask, 1));  % exclude diagonal
pos_mask_indices = [p_row, p_col];
neg_mask_indices = [n_row, n_col];

% save metrics
writematrix(same_pos_mask, 'pos_mask_shapes.csv'); % pos mask
writematrix(same_neg_mask, 'neg_mask_shapes.csv'); % neg mask
%save('pos_mask_all_shapes_alldist.mat', 'pos_mask_all'); % pos mask at every iter
%save('neg_mask_all_shapes_alldist.mat', 'neg_mask_all'); % neg mask at every iter
%writematrix(pos_mask_indices, 'pos_mask_indices_shapes_alldist.csv'); % pos list
%writematrix(neg_mask_indices, 'neg_mask_indices_shapes_alldist.csv'); % neg list
%writematrix(behav_pred, 'behav_pred_shapes_old19.csv'); % predicted comb
