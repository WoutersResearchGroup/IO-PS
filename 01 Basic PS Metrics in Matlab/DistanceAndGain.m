function [DistanceAndOpporGain,Densities,Distance,ProxSums] = DistanceAndGain(Prox,RCA,Products,ProductCompInd)

%Progress = 'DistanceAndGain_Start'

%DistanceAndOpporGain:  %1) HsCode; 2) Distance; 3) Distance if opportunity; 4) OpporGain; 5) OpporGain if Opportunity 6) Density; 7) Density if Oppor


% This code creates calculates the "OpenForest" for a country and the
% contribution to open forest of each of the development opportunities at
% the raw data level

% Oppors  Format = (1 HSCodes; 2 Distances; 3 Distance if Unexploited; 4 RCA of Product; 5 Opportunity gain; 6 Adapted Opportunity gain; 7 Densities; 8 Densities if unexploited)

% It creates the intermediate matrixes:
%1) SAPS (1 hscode; 2 RCA)
%2) BelPS (1 hscode; 2 RCA)

%It requires the loading of:
%1) Prox(Products,Products)
%2) SARCAMat ( 1 hs92code; 2 Export of Country; 3 World Export of prduct; 4 SA Good Percentage of world export of good; 5 % SA Good Percentage of SA exports; 6  Good Percentage of world exports; 7 SA RCA)
%3) Products
%4) ProductCompInd

Products;
RCA;  % Format (1 hs92code; SA RCA of good)

Oppors = zeros(size(Products,1),8);
Oppors(:,1) = Products;
Oppors(:,4) = RCA(:,2);

Densities = zeros(size(Products,1),1);
Distance = zeros(size(Products,1),1);

CheckSumSA = zeros(size(Products,1),1); 
OpenForestContrSum = zeros(size(Products,1),1); 
ContributionToOpenForest = zeros(size(Products,1),1); 
ContributionToOpenForestSAAdapted = zeros(size(Products,1),1); 

%% Calculate Sums Of Proximities per product

ProxSums = zeros(size(Products,1),1); 

for i=1:size(Products,1)
   
    ProxSums(i) = sum(Prox(i,:));
    
    i;
    
end


%% Calculate densities; Distance and Open Forest


%Calculate Total Opp value for country:

OpporValue = 0;

for i=1:size(Products,1) % Run through all potential products
    
    %Calculate SA Densities and Open Forest
    
    DensityNumerator = 0;
    DistanceNumerator = 0;
    
    
    
    for j=1:size(Products,1) %Run through all columns of the proximity matrix
        
        if RCA(j,2) > 1 % If you are already competitive; calculate Density from competitive to opportunity otherwise Distance and open forest
            
            if i ~= j % Make sure product is not not itself
        
                DensityNumerator = DensityNumerator + Prox(i,j);
                
            end
            
        else % Consider other products that could be unlocked from opportunity
            
            if i ~= j % Make sure product is not not itself
                
                 DistanceNumerator = DistanceNumerator + Prox(i,j); 
                
                
                if RCA(i,2) < 1 % Check if original activity is indeed an opportunity
                
                    OpenForestContrSum(i) = OpenForestContrSum(i) + Prox(i,j)*ProductCompInd(j)/(ProxSums(j) - 1); 
                    
                
                end
                
            end
            
        end
    end
    
%     CheckSumSA(i) = CheckSumSA(i) / (sum(Prox6D(i,:)) - 1);
    
    Densities(i) = DensityNumerator / (ProxSums(i) - 1);
    Distance(i) = DistanceNumerator / (ProxSums(i) - 1);
    
    
    
    ContributionToOpenForest(i) = OpenForestContrSum(i) - Densities(i)* ProductCompInd(i);
    ContributionToOpenForestSAAdapted(i) = OpenForestContrSum(i);
    
    
    Oppors(i,2) = Distance(i,1);
  
    Oppors(i,5) = ContributionToOpenForest(i);
    Oppors(i,6) = ContributionToOpenForestSAAdapted(i);
    Oppors(i,7) = Densities(i);
    
    if RCA(i,2) <= 1 % If opportunity not exploited yet, add to unexploited column
        
        Oppors(i,8) = Oppors(i,7);
        Oppors(i,3) = Oppors(i,2);
        
        OpporValue = OpporValue + Densities(i) * ProductCompInd(i);
        
    end
    
        
        
    
    i; %Keep track of execution
    
end

% UnitySA = CheckSumSA + SADensities6D; 
% UnityBel = CheckSumBel + BelDensities6D; 

DistanceAndOpporGain = Oppors(:,[1,2,3,5,6,7,8]);  %1) HsCode; 2) Distance; 3) Distance if opportunity; 4) OpporGain; 5) OpporGain if Opportunity 6) Density; 7) Density if Oppor

% dlmwrite('DensitiesAndOpenForestSA6D.txt',DensitiesAndOpenForestSA6D,'precision',10)

%Progress = 'DistanceAndGain_Finish'


end

