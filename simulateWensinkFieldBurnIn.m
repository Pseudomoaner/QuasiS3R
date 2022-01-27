function [field] = simulateWensinkFieldBurnIn(startField,fS,dS)
%Separate simulation to allow simulation to settle into passive configuration without cell motility, growth or tilting taking place.
field = startField;

for i = 1:fS.burnInSteps
    fprintf('Burn-in frame is %i of %i.\n',i,fS.burnInSteps)
    
    field = field.stepModel(fS.burnIndt,0,fS.divThresh,fS.postDivMovement,dS.colourCells,'burnIn');
end