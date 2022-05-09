clear all
close all

root = 'C:\Users\olijm\Google Drive\CDI_modelling\fire3';
branch = 'simulationResults_ContactData_F';
twigs = {'0pt5.mat','0pt625.mat','0pt75.mat','0pt875.mat','1.mat','1pt125.mat','1pt25.mat','1pt375.mat','1pt5.mat','1pt625.mat','1pt75.mat','1pt875.mat','2.mat'};

figure(1)
ax1 = gca;
hold(ax1,'on')

for b = 1:size(twigs,2)
    load(fullfile(root,[branch,twigs{b}]))
    maxInd = size(contactSet{1},1);

    rmsd = applyRmsd(maxInd,data,fieldSettings.dt);
    plot(ax1,log((0:(maxInd-1))*fieldSettings.dt),log(rmsd))
end