
filepos = 'candidates_siginj_allfiles_sorted.txt';
fileneg = 'candidates_siginj_negchirptimes_allfiles_v2_sorted.txt';


%Get positive candidates
poscandidates = textread(filepos, '%s', 'delimiter', '\n');
numcandidates = length(poscandidates);

positivesnrs = []; positivetoas = [];
for i = 1:numcandidates
posvals = str2num(poscandidates{i});
positivesnrs = [positivesnrs, posvals(end)];
positivetoas = [positivetoas, posvals(end-1)];
end

%Get negative candidates
negcandidates = textread(fileneg, '%s', 'delimiter', '\n');
negativesnrs = []; negativetoas = []; injectedsnrs = [];
for i = 1:numcandidates
negvals = str2num(negcandidates{i});
injectedsnrs = [injectedsnrs, negvals(end)];
negativesnrs = [negativesnrs, negvals(end-1)];
negativetoas = [negativetoas, negvals(end-2)];
end



%Get relative snr and toa diff

reltoadiff = positivetoas - negativetoas;
abstoadiff = abs(reltoadiff);

relsnrdiff = (positivesnrs - negativesnrs)./positivesnrs;
absrelsnrdiff = abs(relsnrdiff);

%Segregate into GLITCH, INJECTED SIGNALS and NOISE ONLY
glitchidxs = []; glitchrelsnrdiff = []; glitchreltoadiff = []; 
noiseidxs = []; noiserelsnrdiff = []; noisereltoadiff =[]; 
injsigidxs = []; injsigrelsnrdiff = []; injsigreltoadiff = []; 
for i = 1:numcandidates
    if injectedsnrs(i) ~=0
        injsigidxs = [injsigidxs, i];
        injsigrelsnrdiff = [injsigrelsnrdiff, relsnrdiff(i)];
        injsigreltoadiff = [injsigreltoadiff, reltoadiff(i)];
    else
        if positivesnrs(i) >= 8.5
            glitchidxs = [glitchidxs,i];
            glitchrelsnrdiff = [glitchrelsnrdiff, relsnrdiff(i)];
            glitchreltoadiff = [glitchreltoadiff, reltoadiff(i)];
        else
            noiseidxs = [noiseidxs,i];
            noiserelsnrdiff = [noiserelsnrdiff, relsnrdiff(i)];
            noisereltoadiff =[noisereltoadiff, reltoadiff(i)];
        end
    end
end


absglitchrelsnrdiff = abs(glitchrelsnrdiff);
absglitchreltoadiff = abs(glitchreltoadiff);
absnoiserelsnrdiff = abs(noiserelsnrdiff);
absnoisereltoadiff = abs(noisereltoadiff);
absinjsigrelsnrdiff = abs(injsigrelsnrdiff);
absinjsigreltoadiff = abs(injsigreltoadiff);


%Make figures
figure;
sz= 100;
hold on;
scatter(glitchrelsnrdiff, glitchreltoadiff,sz,'red','filled','DisplayName','GLITCHES');
hold on;
scatter(noiserelsnrdiff, noisereltoadiff,sz, 'blue','filled','DisplayName','NOISE ONLY');
hold on;
scatter(injsigrelsnrdiff, injsigreltoadiff, sz,'black','filled','DisplayName','INJECTED SIGNALS');
hold off;
legend;
xlabel('\Delta SNR/SNR');
ylabel('\Delta TOA');
title('Relative SNR Difference vs TOA Difference');

figure;
hold on;
scatter(absglitchrelsnrdiff, absglitchreltoadiff,sz,'red','filled','MarkerEdgeColor','black','DisplayName','GLITCHES');
hold on;
scatter(absnoiserelsnrdiff, absnoisereltoadiff,sz, 'blue','filled','MarkerEdgeColor','black','DisplayName','NOISE ONLY');
hold on;
scatter(absinjsigrelsnrdiff, absinjsigreltoadiff, sz,'green','filled','MarkerEdgeColor','black','DisplayName','INJECTED SIGNALS');
hold off;
legend;

xlabel('|\Delta SNR/SNR|');
ylabel('|\Delta TOA|');
ax = gca;
ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
legend('FontSize',20);
title('Absolute Relative SNR Difference vs Absolute TOA Difference', 'FontSize',40);