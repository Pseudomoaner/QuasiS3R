clear all
close all

root = '/home/omeacock/Documents/SPRruns/DiffusionConstRuns';
branch = 'SimulationResults';
twigs = {'0.5.mat','0.75.mat','1.mat','1.25.mat','1.5.mat','1.75.mat','2.mat'};

figure(1)
ax1 = gca;
hold(ax1,'on')
ax2 = axes('Position',[0.6,0.2,0.25,0.3],'Units','Normalized');
hold(ax2,'on')

Ds = zeros(size(twigs));
Vs = zeros(size(twigs));

for b = 1:size(twigs,2)
    load(fullfile(root,[branch,twigs{b}]))
    maxInd = fieldSettings.maxF + 1;

    rmsd = applyRmsd(maxInd,data,fieldSettings.dt);
    times = (0:(maxInd-1))*fieldSettings.dt;
    locGrad = diff(log10(rmsd))./diff(log10(times'));
    
    halfGradXs = find(diff(locGrad>0.5));
    sampPt = halfGradXs(2);
    y0 = log10(rmsd(sampPt)) - 0.5*log10(times(sampPt));
    
    Ds(b) = 0.25*(10^(2*y0));
    Vs(b) = mean(arrayfun(@(x)mean(x.vmag),data));
    
    currCol = [(b-1)/(size(twigs,2)-1),1-(b-1)/(size(twigs,2)-1),1];
    
    loglog(ax1,times,rmsd,'LineWidth',0.8,'Color',currCol)
    loglog(ax1,[times(2),times(end)],2*sqrt(Ds(b)*[times(2),times(end)]),'--','Color',currCol,'LineWidth',1.5)

    plot(ax2,Vs(b),Ds(b),'.','MarkerSize',15,'Color',currCol)
end

% axis(ax1,'equal')
ax1.XScale = 'log';
ax1.YScale = 'log';
axis(ax1,[times(1),times(end),0.1,200])
ax1.Box = 'on';
ax1.LineWidth = 1.5;

alphD = Vs'\Ds';
plot(ax2,[0,0.3],[0,0.3]*alphD,'k:','LineWidth',1.5)
ax2.Box = 'on';
ax2.LineWidth = 1.5;
fprintf('alphD is %f\n',alphD)