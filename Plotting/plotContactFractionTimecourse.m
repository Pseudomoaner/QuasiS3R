clear all
close all

Root = 'C:\Users\olijm\Desktop\SeanAna\ContactFracs';

setA = 'SimulationResults_Rho%s_2_ContactFracs.mat';
setB = 'SimulationResults_Rho%s_ContactFracs.mat';

twigs = {'0.01','0.0065','0.003','0.0018','0.001','0.00065','0.0003','0.00018','0.0001'};
densVals = cellfun(@(x)str2num(x),twigs);

showPlot = logical([1,0,1,0,1,0,1,0,0]); %Example traces to be shown in ax1

vAvg = 0.1; %Average velocity of rods in these simulations

figure
ax1 = gca;
hold(ax1,'on')
ax2 = axes('Position',[0.6,0.2,0.25,0.3],'Units','Normalized');
hold(ax2,'on')

fitStore = zeros(size(twigs,2),2);

for i = 1:size(twigs,2)
    currColInd = (i-1)/(size(twigs,2)-1);
    currCol = [currColInd,0.6,1-currColInd];

    load(fullfile(Root,sprintf(setA,twigs{i})))
    loglog(ax2,str2num(twigs{i}),yfit.c/vAvg,'.','MarkerSize',10,'Color',currCol)
    fitStore(i,1) = yfit.c/vAvg;

    load(fullfile(Root,sprintf(setB,twigs{i})))
    loglog(ax2, str2num(twigs{i}),yfit.c/vAvg,'.','MarkerSize',10,'Color',currCol)
    fitStore(i,2) = yfit.c/vAvg;
    
    if showPlot(i)
        x = 0:5:1000;
        yDat = yfit.a*(1+yfit.b*exp(-yfit.c*x));
        plot(ax1,x,ttosContactFrac,'Color',currCol)
        plot(ax1,x,yDat,'--','Color',currCol,'LineWidth',1.5)

        loglog(ax2, str2num(twigs{i}),yfit.c/vAvg,'.','MarkerSize',15,'Color',currCol)
    end
end

logx = [log(densVals)',log(densVals)'];
logy = log(fitStore);

ft = fittype('0.5*x + a'); %Fit of type alphM * x^1/2, log transformed
alphMFit = fit(logx(:),logy(:),ft,'StartPoint',1);
alphM = exp(alphMFit.a);

plot(ax2, densVals, sqrt(densVals)*alphM,'k--','LineWidth',1.5)

ax1.Box = 'on';
ax1.LineWidth = 1.5;
axis(ax1,[0,1000,0,0.9])

ax2.XScale = 'log';
ax2.YScale = 'log';
ax2.Box = 'on';
ax2.LineWidth = 1.5;
axis(ax2,[0.00005,0.02,0.0065,0.1])

fprintf('Value of alphM is %d.\n',alphM)