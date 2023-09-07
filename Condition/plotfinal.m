candidate_file = 'candidates_siginj_allfiles_sorted.txt';
injectedsigfile = 'injectedsigs_sorted.txt';
%Get all candidates
candidates = textread(candidate_file, '%s', 'delimiter', '\n');
snrs = []; tau0 = []; tau1p5 = [];
for i = 1:length(candidates)
chirpvals = str2num(candidates{i});
tau0 = [tau0, chirpvals(2)];
tau1p5 = [tau1p5, chirpvals(3)];
snrs = [snrs, chirpvals(end)];
end


% Get injected signals
injectedtau1p5 = []; injectedtau0 = [];
candidates = textread(injectedsigfile, '%s', 'delimiter', '\n');
for i = 1:length(candidates)
chirpvals = str2num(candidates{i});
injectedtau0 = [injectedtau0, chirpvals(2)];
injectedtau1p5 = [injectedtau1p5, chirpvals(3)];
end


%Make Plot
figure; hold on;
colormap winter;
imagesc(x,y, Fmatrix'); axis xy;
hold on;
boundary_plot;
hold on;
scatter(tau0, tau1p5, snrs,'black','DisplayName','PSO Estimated Locations');
hold on;
scatter(injectedtau0, injectedtau1p5, 50,'x','red','DisplayName','Injected Signal Locations');
xlim([0,90]);
ylim([0,2]);
hold off;
legend;