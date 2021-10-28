clear all
close all

root = 'C:\Users\olijm\Desktop\CDI_modelling\fire3';
branches = {'SimulationResults'};
extension = '.mat';

rateStore = zeros(size(branches));

%Loop through all your simulations
for i = 1:size(branches,1)
    load([root,filesep,branches{i},extension])

    %Extract the contact fraction timecourses
    [stotContactFrac,ttosContactFrac] = calculateContactFracs(trackableData,fieldSettings); %stotContactFrac means sensitive to toxic contact fraction, and vica versa for ttosContactFrac

    wellMixedFrac = sum(trackableData.Population{1} == 's')/size(trackableData.Population{1},1); %Expected final value of ttosContactFrac
    
    %Find best fit exponential mixing rate
    time = fieldSettings.dt*(0:size(trackableData.Population,2)-1);

    ft = fittype('a*(1-exp(-b*x))','problem','a');
    [yfit, gof] = fit(time',ttosContactFrac',ft,'problem',wellMixedFrac,'StartPoint',0.1);

    rateStore(i) = yfit.b;
    
    %Plot the results (can comment out for batch analysis)
    plot(yfit)
    hold on
    plot(time,ttosContactFrac)
    xlabel('Time')
    ylabel('Mixed toxic/sensitive contacts:total toxic contacts')
end