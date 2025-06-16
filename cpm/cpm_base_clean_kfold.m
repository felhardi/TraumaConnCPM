clear;
clc;

% ------------ INPUTS -------------------

load(fullfile('roi_matrices', 'corr_matrices.mat'));
load(fullfile('trauma_data', 'trauma_select.mat'));

matrix = matrix_array;
tc = trauma_select;

% fisher z-transform the matrices
tolerance = 1e-12; % to avoid log(0)
matrix(matrix >= 1 - tolerance) = 1 - tolerance;
matrix(matrix <= -1 + tolerance) = -1 + tolerance;

z_mats = 0.5 * log((1 + matrix) ./ (1 - matrix));

all_mats = z_mats;
all_behav = tc.total_5;
age = tc.age;

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

% Create 10 folds for cross-validation
cv = cvpartition(no_sub, 'KFold', 10);

% leave one out train/test

for fold = 1:cv.NumTestSets

    fprintf('\nProcessing Fold %d/%d', fold, cv.NumTestSets);
    
    test_idx = test(cv, fold); % Logical index for test subjects
    train_idx = training(cv, fold); % Logical index for train subjects

    % Extract train and test data
    train_mats = all_mats(:,:,train_idx);
    train_vcts = reshape(train_mats,[],size(train_mats,3));
    train_behav = all_behav(train_idx);
    train_age = age(train_idx);

    test_mats = all_mats(:,:,test_idx);
    test_behav = all_behav(test_idx);

    % partial correlation
    [r_mat, p_mat] = partialcorr(train_vcts', train_behav, train_age, "Type", "Spearman");

    r_mat = reshape(r_mat,no_node,no_node);
    p_mat = reshape(p_mat,no_node,no_node);
    
    % Identify significant positive and negative edges
    pos_mask = (r_mat > 0 & p_mat < thresh);
    neg_mask = (r_mat < 0 & p_mat < thresh);

    % Store masks for this fold
    pos_mask_all(:,:,test_idx) = repmat(pos_mask, [1, 1, sum(test_idx)]);
    neg_mask_all(:,:,test_idx) = repmat(neg_mask, [1, 1, sum(test_idx)]);

    % Compute summed network strengths for training data
    train_sumpos = squeeze(sum(sum(train_mats .* pos_mask, 1), 2)) / 2;
    train_sumneg = squeeze(sum(sum(train_mats .* neg_mask, 1), 2)) / 2;

    % Train regression models
    fit_pos = regress(train_behav, [train_sumpos, ones(size(train_sumpos))]);
    fit_neg = regress(train_behav, [train_sumneg, ones(size(train_sumneg))]);
    b = regress(train_behav, [train_sumpos, train_sumneg, ones(size(train_sumpos))]);

    % Apply model to test subjects
    for i = find(test_idx)'
        test_mat = all_mats(:,:,i);
        test_sumpos = sum(sum(test_mat .* pos_mask)) / 2;
        test_sumneg = sum(sum(test_mat .* neg_mask)) / 2;

        behav_pred_pos(i) = fit_pos(1) * test_sumpos + fit_pos(2);
        behav_pred_neg(i) = fit_neg(1) * test_sumneg + fit_neg(2);
        behav_pred(i) = b(1) * test_sumpos + b(2) * test_sumneg + b(3);
    end
end

% compare predicted and observed scores
[R_pos, P_pos] = corr(behav_pred_pos,all_behav,"type","Spearman")
[R_neg, P_neg] = corr(behav_pred_neg,all_behav,"type","Spearman")
[R_comb, P_comb] = corr(behav_pred, all_behav,"type","Spearman")

% find edges that were consistently significant across all iterations
%same_pos_mask = all(pos_mask_all == 1, 3);
%same_neg_mask = all(neg_mask_all == 1, 3);

% find indices for the upper triangular part and combine
%[p_row, p_col] = find(triu(same_pos_mask, 1));  % exclude diagonal
%[n_row, n_col] = find(triu(same_neg_mask, 1));  % exclude diagonal
%pos_mask_indices = [p_row, p_col];
%neg_mask_indices = [n_row, n_col];

% save metrics
%writematrix(same_pos_mask, 'pos_mask_modage.csv'); % pos mask
%writematrix(same_neg_mask, 'neg_mask_modage.csv'); % neg mask
%save('pos_mask_all.mat', 'pos_mask_all'); % pos mask at every iter
%save('neg_mask_all.mat', 'neg_mask_all'); % neg mask at every iter
%writematrix(pos_mask_indices, 'pos_mask_indices_age.csv'); % pos list
%writematrix(neg_mask_indices, 'neg_mask_indices_age.csv'); % neg list
%writematrix(behav_pred, 'behav_pred_modage_10k.csv'); % predicted neg
