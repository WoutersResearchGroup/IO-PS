function [OpporValue] = OpporValue(Prox,RCA,Products,ProductCompInd,ProxSums)

%DistanceAndOpporGain:  %1) HsCode; 2) Distance; 3) Distance if opportunity; 4) OpporGain; 5) OpporGain if Opportunity 6) Density; 7) Density if Oppor
% This code calculates the "OpenForest" / Opportunity value for a country

Densities = zeros(size(Products,1),1);


%% Calculate densities; Distance and Open Forest


%Calculate Total Opp value for country:

OpporValue = 0;

for i=1:size(Products,1) % Run through all potential products
    
    if RCA(i,2) <= 1 % If opportunity not exploited yet
        
        DensityNumerator = 0;

        for j=1:size(Products,1) %Run through all columns of the proximity matrix

            if RCA(j,2) > 1 % If you are already competitive; calculate Density from competitive to opportunity

                if i ~= j % Make sure product is not not itself

                    DensityNumerator = DensityNumerator + Prox(i,j);

                end

            end
        end

        Densities(i) = DensityNumerator / (ProxSums(i) - 1);


        %if RCA(i,2) <= 1 % If opportunity not exploited yet, add to unexploited column

        OpporValue = OpporValue + Densities(i) * ProductCompInd(i);
        
    end
    
end


end