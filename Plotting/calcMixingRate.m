clear all
close all

root = '/home/omeacock/Documents/SPRruns/CDI_Mixing_Fitting/';
branch = 'SimulationResults_Rho';
twigs = {'0.01_2','0.0065_2','0.003_2','0.0018_2','0.001_2','0.00065_2','0.0003_2','0.00018_2','0.0001_2'};
extension = '.mat';

rateStore = zeros(size(twigs));

%Loop through all your simulations
for i = 1:size(twigs,2)
    load([root,filesep,branch,twigs{i},extension])

    %Extract the contact fraction timecourses
    [stotContactFrac,ttosContactFrac] = calculateContactFracs(trackableData,fieldSettings); %stotContactFrac means sensitive to toxic contact fraction, and vica versa for ttosContactFrac

    wellMixedFrac = sum(trackableData.Population{1} == 's')/size(trackableData.Population{1},1); %Expected final value of ttosContactFrac
    expCorr = (2*sqrt(patchSettings.seedDensity)*fieldSettings.lam-1); %Correction factor to account that the initial contact fraction won't be exactly 0.
    
    %Find best fit exponential mixing rate
    time = fieldSettings.dt*(0:size(trackableData.Population,2)-1);

    ft = fittype('a*(1+b*exp(-c*x))','problem',{'a','b'});
    [yfit, gof] = fit(time',ttosContactFrac',ft,'problem',{wellMixedFrac,expCorr},'StartPoint',0.1);

    rateStore(i) = yfit.c;
    
    %Plot the results (can comment out for batch analysis)
    plot(yfit)
    hold on
    plot(time,ttosContactFrac)
    xlabel('Time')
    ylabel('Mixed toxic/sensitive contacts:total toxic contacts')
    
    save([root,filesep,branch,twigs{i},'_ContactFracs',extension],'wellMixedFrac','yfit','stotContactFrac','ttosContactFrac')
    
    yfit.c
    pause(0.1)
end