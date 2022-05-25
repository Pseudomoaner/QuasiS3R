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
    
    atFracs = cellfun(@(x)sum(x == 't'),trackableData.Population)/size(trackableData.Population{1},1);
    seFracs = cellfun(@(x)sum(x == 0),trackableData.Hit)/size(trackableData.Hit{1},1) - atFracs;
    hitFracs = cellfun(@(x)sum(x > 0),trackableData.Hit)/size(trackableData.Hit{1},1);
    
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
    plot(atatConts./totConts,'--','Color',[32,128,196]/255,'LineWidth',1.5)
    %Attacker-sensitive contacts
    plot(atseConts./totConts,'--','Color',[255, 120, 10]/255,'LineWidth',1.5)
    %Attacker-hit contacts
    plot(athitConts./totConts,'--','Color',[1,0.85,0.75],'LineWidth',1.5)
    
    %Attacker fractions
    plot(atFracs,'Color',[32,128,196]/255,'LineWidth',1.5)
    %Sensitives fractions
    plot(seFracs,'Color',[255, 120, 10]/255,'LineWidth',1.5)
    %Hit fractions
    plot(hitFracs,'Color',[1,0.85,0.75],'LineWidth',1.5)
    
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    axis([0,1000,0,1])
end