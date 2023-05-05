% This script relies on you having generated the contactSet variable from
% the SPR simulations; can be made by calling the simulation function as
% [PCs,endField,contactSet] = simulateWensinkField(intermediateField,fieldSettings,dispSettings);

clear all
close all

root = '/home/omeacock/Documents/SPRruns/ContactRateRuns';
branch = 'SimulationResults_F=%f.mat';
% twigs = {'0pt5.mat','0pt625.mat','0pt75.mat','0pt875.mat','1.mat','1pt125.mat','1pt25.mat','1pt375.mat','1pt5.mat','1pt625.mat','1pt75.mat','1pt875.mat','2.mat'};
fList = 0.5:0.125:2;

figure(1)
ax1 = gca;
hold(ax1,'on')
ax2 = axes('Position',[0.2,0.6,0.25,0.3],'Units','Normalized');
hold(ax2,'on')

meanContNoStore = zeros(size(fList));
contRateStore = zeros(size(fList));
vStore = zeros(size(fList));

minContTime = 5; %Discounts extremely transient contacts

for fInd = 1:size(fList,2)
    f = fList(fInd);
    load(fullfile(root,sprintf(branch,f)))
    maxInd = size(contactSet{1},1);
    
    vDt = fieldSettings.FrameSkip * fieldSettings.motiledt; %The time between sampled points for trajectories
    cDt = fieldSettings.FireSkip * fieldSettings.motiledt; %The time between sampled points for contact detection
    
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
        contRateTmp(r) = sum(sum(expandedConts,2)>minContTime)/(size(contactSet,2)*cDt);
    end

    vStore(fInd) = mean(arrayfun(@(x)mean(x.vmag),data));
    meanContNoStore(fInd) = mean(meanContNoTmp);
    contRateStore(fInd) = mean(contRateTmp);

    currCol = [(fInd-1)/size(fList,2),1-((fInd-1)/size(fList,2)),1];
    plot(ax1,vStore(fInd),contRateStore(fInd),'.','MarkerSize',20,'Color',currCol)    
    plot(ax2,vStore(fInd),meanContNoStore(fInd),'.','MarkerSize',20,'Color',currCol)
end

alphK = vStore'\contRateStore';
plot(ax1,[0,0.3],[0,0.3*alphK],'--','LineWidth',1.5,'Color','k')
plot(ax2,[0,0.3],[5,5],'--','LineWidth',1.5,'Color','k')

ax1.LineWidth = 1.5;
ax1.Box = 'on';

ax2.LineWidth = 1.5;
ax2.Box = 'on';
axis(ax2,[0,0.3,0,7])

xlabel(ax1,'Average rod velocity')
ylabel(ax1,'Average lost contact rate')

disp(['alphaK is ', num2str(alphK)])
