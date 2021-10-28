function contactLengths = findContactLengthDistribution(contactSet)
%FINDCONTACTLENGTHDISTRIBUTION finds the length of contacts between rods,
%as calculated instantaneously using the WensinkField.calculateContacts
%function.
%
%   INPUTS:
%       -contactSet: A Tx1 cell array, with each cell containing the output
%       of WensinkField.calculateContacts at a given timepoint.
%
%   OUTPUTS:
%       -contactLengths: 1D vector indicating the number of contacts that 
%       last for one timeunit, two timeunits etc.
%
%   Author: Oliver J. Meacock, (c) 2021

contactLengths = zeros(size(contactSet));

for i = 1:size(contactSet,2)
    for j = 1:size(contactSet{i},1)
        for k = 1:size(contactSet{i}{j},1)
            currLen = 1;
            ID = contactSet{i}{j}(k);
            contactSet{i}{j}(contactSet{i}{j} == ID) = NaN;
            
            %Go through each subsequent timepoint, either removing the object
            %from that object's storage and incrementing the length counter and
            %continuing, or terminating the loop if the ID can't be found at
            %this timepoint
            if 1 + i < size(contactSet,2)
                while sum(contactSet{i + currLen}{j} == ID) == 1
                    contactSet{i + currLen}{j}(contactSet{i + currLen}{j} == ID) = NaN;
                    currLen = currLen + 1;
                    
                    if currLen + i > size(contactSet,2)
                        break
                    end
                end
            end
            
            if ~isnan(ID)
                contactLengths(currLen) = contactLengths(currLen) + 1;
            end
        end
    end
end