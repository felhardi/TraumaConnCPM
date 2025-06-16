clear;
clc;

% ------------ INPUTS -------------------

load('corr_matrices_kn.mat');  
load('trauma_kn.mat');

matrix = matrix_array_kn;
tc = trauma_kn;

% Fisher z-transform the matrices
matrix(matrix >= 1 - (1e-12)) = 1 - (1e-12);
matrix(matrix <= -1 + (1e-12)) = -1 + (1e-12);

z_mats = 0.5 * log((1 + matrix) ./ (1 - matrix));

all_mats = z_mats;
all_behav = tc.total_5;
age = tc.age;

% Calculate the true prediction correlation
[R_pos, R_neg, R_comb] = predict_behavior(all_mats, all_behav, age);

% Number of iterations for permutation testing
no_iterations = 1000;
prediction_r = zeros(no_iterations, 3);

% Store the true prediction result
prediction_r(1, :) = [R_pos, R_neg, R_comb];

pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

% Create an estimate distribution of the test statistic via random shuffles of data labels
parfor it = 2:no_iterations
    fprintf('\n Performing iteration %d out of %d', it, no_iterations);
    new_behav = all_behav(randperm(length(all_behav)));
    [temp_r_pos, temp_r_neg, temp_r_comb] = predict_behavior(all_mats, new_behav, age);
    
    % Store results in the temporary array, then move to main array outside the loop
    prediction_r(it, :) = [temp_r_pos, temp_r_neg, temp_r_comb];
end

% Sorting and calculating p-values
sorted_prediction_r_pos = sort(prediction_r(:, 1), 'descend');
position_pos = find(sorted_prediction_r_pos == R_pos);
pval_pos = position_pos(1) / no_iterations;

sorted_prediction_r_neg = sort(prediction_r(:, 2), 'descend');
position_neg = find(sorted_prediction_r_neg == R_neg);
pval_neg = position_neg(1) / no_iterations;

sorted_prediction_r_comb = sort(prediction_r(:, 3), 'descend');
position_comb = find(sorted_prediction_r_comb == R_comb);
pval_comb = position_comb(1) / no_iterations;

display(pval_pos)
display(pval_neg)
display(pval_comb)

%save('permtest_wage_wsex.mat', 'pval_pos', 'pval_neg', 'pval_comb');