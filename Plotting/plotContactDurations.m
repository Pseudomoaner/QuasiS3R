% This script relies on you having generated the contactSet variable from
% the SPR simulations; can be made by calling the simulation function as
% [PCs,endField,contactSet] = simulateWensinkField(intermediateField,fieldSettings,dispSettings);

clear all
close all

root = 'C:\Users\olijm\Google Drive\CDI_modelling\fire3';
branch = 'simulationResults_ContactData_F';
twigs = {'0pt5.mat','0pt625.mat','0pt75.mat','0pt875.mat','1.mat','1pt125.mat','1pt25.mat','1pt375.mat','1pt5.mat','1pt625.mat','1pt75.mat','1pt875.mat','2.mat'};
showPlot = [1,0,1,0,1,0,1,0,1,0,1,0,1];

lenThresh = 51; % Threshold below which the tail of very transient contacts gets cut off in the curve fitting part of the algorithm. These don't fit the exponential model well

figure(1)
ax1 = gca;
ax2 = axes('Position',[0.5,0.5,0.35,0.4],'Units','Normalized');
hold(ax1,'on')
hold(ax2,'on')

vStore = zeros(size(twigs));
gamStore = zeros(size(twigs));

for b = 1:size(twigs,2)
    load(fullfile(root,[branch,twigs{b}]))

    lenStore = [];

    %Outer loop - through each rod
    for r = 1:size(contactSet{1},1)
        currConts = contactSet{1}{r};
        contLens = ones(size(currConts));

        %Inner loop - through each timepoint
        for t = 2:size(contactSet,2)
            nextConts = contactSet{t}{r};

            %Find the contacts between target object and field of other objects
            %that are created or lost between previous and current timepoint
            [newConts,newContsInds] = setdiff(nextConts,currConts);
            [lostConts,lostContsInds] = setdiff(currConts,nextConts);

            %Store data for lost contacts and wipe from temporary length storage
            lenStore = [lenStore;contLens(lostContsInds)*fieldSettings.dt];
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
    end
    
    disp(['Average number of contacts =', num2str(sum(lenStore(lenStore > lenThresh))/(size(contactSet,2)*size(contactSet{1},1)))])

    [N,Edges] = histcounts(lenStore,'BinWidth',fieldSettings.dt,'Normalization','pdf');

    edgeLim = find(diff(Edges > lenThresh) == 1);
    ft = fit(Edges(edgeLim:end-1)'+0.5,N(edgeLim:end)','exp1');
    
    currCol = [(b-1)/size(twigs,2),1-((b-1)/size(twigs,2)),1];
    
    vStore(b) = mean(arrayfun(@(x)mean(x.vmag),data));
    gamStore(b) = -ft.b;

    if showPlot(b)
        plot(ax1,Edges(1:end-1)+0.5,N,'LineWidth',1.5,'Color',currCol)
        plot(ax1,Edges(1:end-1)+0.5,ft.a*exp(ft.b*(Edges(1:end-1)+0.5)),':','LineWidth',1.5,'Color',currCol)

        plot(ax2,vStore(b),gamStore(b),'.','MarkerSize',20,'Color',currCol)
    else
        plot(ax2,vStore(b),gamStore(b),'.','MarkerSize',14,'Color',currCol)
    end

end

axis(ax1,[0,500,0,0.015])
ax1.LineWidth = 1.5;
ax1.Box = 'on';

xlabel(ax1,'Rod-rod contact time')
ylabel(ax1,'pdf')

ell = gamStore'\vStore';

plot(ax2,[0,0.3],[0,0.3/ell],'k--','LineWidth',1.5)
axis(ax2,[0,0.3,0,0.02])

ax2.LineWidth = 1.5;
ax2.Box = 'on';

xlabel(ax2,'Mean rod velocity')
ylabel(ax2,'Fitted contact length scale')

disp(['ell is ', num2str(ell)])