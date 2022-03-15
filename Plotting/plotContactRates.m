% This script relies on you having generated the contactSet variable from
% the SPR simulations; can be made by calling the simulation function as
% [PCs,endField,contactSet] = simulateWensinkField(intermediateField,fieldSettings,dispSettings);

clear all
close all

root = 'C:\Users\olijm\Google Drive\CDI_modelling\fire3';
branch = 'simulationResults_ContactData_F';
twigs = {'0pt5.mat','0pt625.mat','0pt75.mat','0pt875.mat','1.mat','1pt125.mat','1pt25.mat','1pt375.mat','1pt5.mat','1pt625.mat','1pt75.mat','1pt875.mat','2.mat'};

lenThresh = 91; % Threshold below which the tail of very transient contacts gets cut off in the curve fitting part of the algorithm. These don't fit the exponential model well

figure(1)
ax1 = gca;
hold(ax1,'on')

rateStore = zeros(size(twigs));
vStore = zeros(size(twigs));

for b = 1:size(twigs,2)
    load(fullfile(root,[branch,twigs{b}]))

    contRateStore = [];

    %Outer loop - through each rod
    for r = 1:size(contactSet{1},1)
        currConts = contactSet{1}{r};
        contLens = ones(size(currConts));
        contNoList = [];

        %Inner loop - through each timepoint
        for t = 2:size(contactSet,2)
            nextConts = contactSet{t}{r};

            %Find the contacts between target object and field of other objects
            %that are created or lost between previous and current timepoint
            [newConts,newContsInds] = setdiff(nextConts,currConts);
            [lostConts,lostContsInds] = setdiff(currConts,nextConts);

            %Store data for number of lost non-transitory contacts this timepoint
            longConts = contLens(lostContsInds) >= (lenThresh/fieldSettings.dt);
            contNoList = [contNoList;sum(longConts)];

            %Store data for lost contacts and wipe from temporary length storage
            contLens(lostContsInds) = [];
            
            %Increment remaining contact lengths
            contLens = contLens + 1;

            %Create new length storage for new contacts
            for i = 1:size(newContsInds,1)
                insertInd = newContsInds(i);
                contLens = [contLens(1:(insertInd-1));1;contLens(insertInd:end)];
            end

            %Update currConts
            currConts = nextConts;
        end
        contRateStore = [contRateStore;sum(contNoList)/(fieldSettings.dt*size(contactSet,2))];
    end

    currCol = [(b-1)/size(twigs,2),1-((b-1)/size(twigs,2)),1];

    rateStore(b) = mean(contRateStore);
    vStore(b) = mean(arrayfun(@(x)mean(x.vmag),data));
    
    plot(ax1,vStore(b),rateStore(b),'.','MarkerSize',20,'Color',currCol)
end

alphK = vStore'\rateStore';
plot([0,0.3],[0,0.3*alphK],'--','LineWidth',1.5,'Color','k')

%axis(ax1,[0,100,0,0.075])
ax1.LineWidth = 1.5;
ax1.Box = 'on';

xlabel(ax1,'Average rod velocity')
ylabel(ax1,'Average lost contact rate')

disp(['alphaK is ', num2str(alphK)])