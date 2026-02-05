% run_param_tables

%% Print parameter estimates 
clear all; close all; 
addpath('matlabinclude');
addpath('autotrade');

% set up parameters and indexes consistent with loaded data
load('results/estimation/estimation_data.mat');

% load parameter estimates
load('results/estimation/mle_converged.mat');

%% Print parameter esimates
mp=mp_mle;

DOLATEX=1; 

% 1st stage driving estimates 
dktax.print_driving_parameters(mp,est,DOLATEX); 
dktax.print_driving_parameters(mp,est,DOLATEX, 'results/tables/tab_est_driving.tex'); 

se_mle=getfields(estim.update_mp(sqrt(diag(Avar_mle)), mp_mle), mp.pnames); 

fid = 1; % file identifier - this should be pointing to different files for each table 

% ntypes * ncartypes 
vv = {'u_0', 'u_a', 'u_a_sq'}; 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [mp.ntypes, mp.ncartypes]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    fprintf('\n--- %s ---\n\n', v); 
    tab = cell2table(mp_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', mp_mle.lbl_types); 
    tab_se = cell2table(se_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', mp_mle.lbl_types); 
    % dktax.table2latex(tab, fid); 
    dktax.table2latex_se(tab, tab_se, fid); 
end


% ncartypes
vv = {'acc_0','acc_a'};
tab = cell(1, numel(vv)); 
tab_se = cell(1, numel(vv)); 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [1, mp.ncartypes]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    fprintf('\n--- %s ---\n\n', v); 
    tab{i} = cell2table(mp_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', {v}); 
	tab_se{i} = cell2table(se_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', {v});     
end
dktax.table2latex_se(vertcat(tab{:}),vertcat(tab_se{:}), fid); 



% ntypes 
vv = {'mum'};
tab = cell(numel(vv),1); 
tab_se = cell(numel(vv),1); 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [mp.ntypes, 1]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    fprintf('\n--- %s ---\n\n', v); 
    tab{i} = cell2table(mp_mle.(v), 'VariableNames', {v}, 'RowNames', mp_mle.lbl_types); 
    tab_se{i} = cell2table(se_mle.(v), 'VariableNames', {v}, 'RowNames', mp_mle.lbl_types); 
end
dktax.table2latex_se([tab{:}], [tab_se{:}], fid); 

% ntypes 
vv = {'psych_transcost', 'psych_transcost_nocar'};
tab = cell(numel(vv),1); 
tab_se = cell(numel(vv),1); 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [mp.ntypes, 1]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    fprintf('\n--- %s ---\n\n', v); 
    tab{i} = cell2table(mp_mle.(v), 'VariableNames', {v}, 'RowNames', mp_mle.lbl_types); 
    tab_se{i} = cell2table(se_mle.(v), 'VariableNames', {v}, 'RowNames', mp_mle.lbl_types); 
end
dktax.table2latex_se([tab{:}], [tab_se{:}], fid); 


% scalars
fprintf('\n--- scalars ---\n\n'); 
vv = {'sigma_s', 'tc_sale', 'tc_sale_even'}; 
x = []; 
x_se = []; 
for i=1:numel(vv) 
    if ismember(vv{i},mp_mle.pnames)
	    x_ = mp_mle.(vv{i}); 
	    x_se_ = se_mle.(vv{i}); 
	    if ~isscalar(x_)
	        warning('Parameter %s appears to be non-scalar: move it elsewhere', vv{i}); 
	        x_ = nan; 
	    end
	    x = [x; x_];
	    x_se = [x; x_se_];
	end
end
tab = array2table(x, 'rownames', vv, 'variablenames', {'Estimate'}); 
tab_se = array2table(x, 'rownames', vv, 'variablenames', {'Estimate'}); 
dktax.table2latex_se(tab, tab_se, fid); 

return
