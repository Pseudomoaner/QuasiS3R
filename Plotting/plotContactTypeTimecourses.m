clear all
close all

RootSim = '~/Documents/SPRruns/Fig1Runs';
outputDirNames = {'Patchy','FrozenPatchy','Nonpatchy','FrozenNonpatchy'};
outputSimName = 'SimulationResults.mat';

for i = 1:size(outputDirNames,2)
    load(fullfile(RootSim,outputDirNames{i},outputSimName))
    
    atatConts = zeros(size(contactSet));
    atseConts = zeros(size(contactSet));
    athitConts = zeros(size(contactSet));
    
    for t = 1:size(contactSet,2)
        for j = 1:size(contactSet{t},1)
            if trackableData.FireRate{t}(j) > 0
                sensConts = and(trackableData.FireRate{t}(contactSet{t}{j}) == 0, trackableData.Hit{t}(contactSet{t}{j})==0);
                atConts = trackableData.FireRate{t}(contactSet{t}{j})>0;
                hitConts = ~or(sensConts,atConts);
                
                atatConts(t) = atatConts(t) + sum(atConts);
                atseConts(t) = atseConts(t) + sum(sensConts);
                athitConts(t) = athitConts(t) + sum(hitConts);
            end
        end
    end
    
    totConts = atatConts + atseConts + athitConts;
    
    subplot(2,2,i)
    hold on
    
    %Attacker-attacker contacts
    plot(atatConts./totConts,'Color',[32,128,196]/255,'LineWidth',1.5)
    %Attacker-sensitive contacts
    plot(atseConts./totConts,'Color',[255, 190, 11]/255,'LineWidth',1.5)
    %Attacker-hit contacts
    plot(athitConts./totConts,'Color',[0.9,0,0.4],'LineWidth',1.5)
    
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
end