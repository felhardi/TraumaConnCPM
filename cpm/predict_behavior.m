function [R_pos, R_neg, R_comb] = predict_behavior(all_mats, all_behav, age)

% threshold for feature selection
thresh = 0.01;

% Initialization
no_sub = size(all_mats, 3);
no_node = size(all_mats, 1);

behav_pred_pos = zeros(no_sub, 1);
behav_pred_neg = zeros(no_sub, 1);
behav_pred =  zeros(no_sub, 1);

% leave one out train/test

    parfor leftout = 1:no_sub

        fprintf('\n Leaving out subj # %6.3f',leftout);
    
        % leave out subject from matrices and behavior
    
        train_mats = all_mats;
        train_mats(:,:,leftout) = [];
        train_vcts = reshape(train_mats,[],size(train_mats,3));

        train_behav = all_behav;
        train_behav(leftout) = [];
        train_age = age;
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
    [R_pos, ~] = corr(behav_pred_pos,all_behav,"type","Spearman");
    [R_neg, ~] = corr(behav_pred_neg,all_behav,"type","Spearman");
    [R_comb, ~] = corr(behav_pred, all_behav,"type","Spearman");

end