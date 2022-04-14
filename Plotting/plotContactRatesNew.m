% This script relies on you having generated the contactSet variable from
% the SPR simulations; can be made by calling the simulation function as
% [PCs,endField,contactSet] = simulateWensinkField(intermediateField,fieldSettings,dispSettings);

clear all
close all

root = 'C:\Users\olijm\Google Drive\CDI_modelling\fire3';
branch = 'simulationResults_ContactData_F';
twigs = {'0pt5.mat','0pt625.mat','0pt75.mat','0pt875.mat','1.mat','1pt125.mat','1pt25.mat','1pt375.mat','1pt5.mat','1pt625.mat','1pt75.mat','1pt875.mat','2.mat'};

figure(1)
ax1 = gca;
hold(ax1,'on')

meanContNoStore = zeros(size(twigs));
contRateStore = zeros(size(twigs));
vStore = zeros(size(twigs));

minContTime = 5; %Discounts extremely transient contacts

for b = 1:size(twigs,2)
    load(fullfile(root,[branch,twigs{b}]))
    maxInd = size(contactSet{1},1);
    
    meanContNoTmp = zeros(size(contactSet{1},1),1);
    contRateTmp = zeros(size(contactSet{1},1),1);

    %Outer loop - through each rod
    for r = 1:size(contactSet{1},1)
        expandedConts = zeros(maxInd,size(contactSet,2)); %Will contain 1s when a particular rod is contacted by this one, 0s otherwise
        
        %Inner loop - through each timepoint
        for t = 1:size(contactSet,2)
            expandedConts(contactSet{t}{r},t) = 1;
        end
        
        %Store the mean number of rods this rod is touching at any given
        %time
        meanContNoTmp(r) = mean(sum(expandedConts,1));
        
        %Set any initial contacts to zero (contact rate is based on new
        %contacts only)
        expandedConts(expandedConts(:,1)==1,:) = 0;
        contRateTmp(r) = sum(sum(expandedConts,2)>minContTime)/(size(contactSet,2)*fieldSettings.dt);
    end

    vStore(b) = mean(arrayfun(@(x)mean(x.vmag),data));
    meanContNoStore(b) = mean(meanContNoTmp);
    contRateStore(b) = mean(contRateTmp);

    currCol = [(b-1)/size(twigs,2),1-((b-1)/size(twigs,2)),1];
    plot(ax1,vStore(b),contRateStore(b),'.','MarkerSize',20,'Color',currCol)
end

alphK = vStore'\contRateStore';
plot([0,0.3],[0,0.3*alphK],'--','LineWidth',1.5,'Color','k')

ax1.LineWidth = 1.5;
ax1.Box = 'on';

xlabel(ax1,'Average rod velocity')
ylabel(ax1,'Average lost contact rate')

disp(['alphaK is ', num2str(alphK)])