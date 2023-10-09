candidate_file = 'GLITCH_VETO_STUDY/GVS_initialrun_sorted.txt';
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

%{
% Get injected signals
injectedtau1p5 = []; injectedtau0 = [];
candidates = textread(injectedsigfile, '%s', 'delimiter', '\n');
for i = 1:length(candidates)
chirpvals = str2num(candidates{i});
injectedtau0 = [injectedtau0, chirpvals(2)];
injectedtau1p5 = [injectedtau1p5, chirpvals(3)];
end
%}

%Get image matrix
S = load('imagvals.mat');
Fmatrix = S.fitvals;
Fmatrix(Fmatrix > 0) = 1;

%Make Plot
sigcolor = 'black';
psocolor = 'red';
x = linspace(0, 90, 10000);
y = linspace(0, 2, 10000);
figure; hold on;
xlim([0,90]);
ylim([0,2]);
map = [153 204 255;
          255 255 133]./255;
colormap(map);
imagesc(x,y, Fmatrix'); axis xy;
hold on;
boundary_plot;
hold on;
scatter(tau0, tau1p5, snrs + 20,psocolor,'DisplayName','PSO Estimated Locations');
hold on;
% scatter(injectedtau0, injectedtau1p5, 70,'x',sigcolor,'DisplayName','Injected Signal Locations');

hold off;
legend;