clear;
clc;
%% load data
load('kdata_full_0727_post_recon.mat');
[nx, ny, nslices, ~] = size(iField);
Recons = abs(iField);

%% generating dictionary
alpha0 =  - pi/22.5;  % - pi/22.5
TEs = 0.002436;
TR = 0.007288; 
nframes = 128;
IR_eff = -1;
R = [];

% Increments of 5% for both parameters
T1 =  .1; while T1(end)<5, T1 = [T1, T1(end)*1.05]; end
T2 = .01; while T2(end)<0.5, T2 = [T2, T2(end)*1.05]; end
T2s = 0;
alpha = [alpha0];

% calculate the dictionary:
D = MRF_dictionary_3contrasts(T1, T2, T2s, alpha, IR_eff, R, TR, TEs, nframes);

%% plot some dictionary items
indices = 2760:10:3080;
idx_legend = 1;
legendText = {};
figure;
plot(1:3, D.magnetization(:, indices(1)));
xticks([1 2 3]);
xlabel('Contrasts');
ylabel('Signal Intensity');
for i = 2 : length(indices)
    hold all;
    plot(1:3, D.magnetization(:, indices(i)));
    if mod(i, 5) == 2
        legendText{idx_legend} = sprintf(['T1 = ', num2str(floor(D.lookup_table(indices(i), 1)*1000)), ' ms, T2 = ', ...
                         num2str(floor(D.lookup_table(indices(i), 2)*1000)), ' ms']);
        idx_legend = idx_legend + 1;
    end
end
legend(legendText);

%% T1T2 matching
qMaps_all = zeros(nx, ny, nslices, 3);
for idx_slice = 1:nslices
    disp(idx_slice);
    x = abs(squeeze(cat(4, Recons(:, :, idx_slice, end-1), Recons(:, :, idx_slice, 1), Recons(:, :, idx_slice, end))));
    x = reshape(x, [size(x,1)*size(x,2) size(x,3)]);

    % Dictionary matching:
    clear c idx
    for q=size(x,1):-1:1
        [c(q,1),idx(q,1)] = max(x(q,:) * abs(D.magnetization), [], 2);
    end
    PD    = c ./ D.normalization(idx).';
    PD    = reshape(PD, [nx ny]);
    qMaps = D.lookup_table(idx,:);
    qMaps = reshape(qMaps, [nx, ny, size(D.lookup_table,2)]);
    qMaps = cat(3, qMaps(:, :, 1:2), PD);
    qMaps_all(:, :, idx_slice, :) = qMaps;
end
save('qMaps_all.mat', 'qMaps_all');