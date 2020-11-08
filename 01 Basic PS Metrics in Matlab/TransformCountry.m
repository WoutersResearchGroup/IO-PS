function [SARCAMat,ProductSum] = TransformCountry(AllMat,CountryRaw,Products)

%This macro transforms the raw data to calculate the country RCA
%It also Sums The world production for a specific good and calculates the
%RCA for the involved countries

% The resulting Arrays are as follows:
SARCAMat = zeros(size(Products,1),7); % Final format = ( 1 hs92code; 2 Export of Country; 3 World Export of product; 4 SA Good Percentage of world export of good; 5 % SA Good Percentage of SA exports; 6  Good Percentage of world exports; 7 SA RCA)
ProductSum = zeros(size(Products,1),2); %Final format = 1 hs92code; Sum of world production

%AllMat; % Format = ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA)
%CountryRaw; % Format =  ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA)

Progress = 'TransformCountry_Start'

SARCAMat(:,1) = Products;
ProductSum(:,1) = Products;


for i = 1:size(Products,1) % run through all HS92s
    
    for j = 1:size(CountryRaw,1) %Run through all Country Exports
    
        if  CountryRaw(j,3) == SARCAMat(i,1)
            
            SARCAMat(i,2) = SARCAMat(i,2) + CountryRaw(j,4); %Add country export to column 2 
            
        end
        
        
    end
    
    
    for j = 1:size(AllMat,1) %Run through all raw world Exports entries
        
        if  AllMat(j,3) == ProductSum(i,1)
            
            ProductSum(i,2) = ProductSum(i,2) + AllMat(j,4); % Add up world production to column 2 of Product Sum Matrix
            
        end
        
    end
    
end

SARCAMat(:,3) = ProductSum(:,2); 

TotalSA = sum(SARCAMat(:,2))
TotalWorld = sum(ProductSum(:,2))

SARCAMat(:,4) = SARCAMat(:,2) ./ SARCAMat(:,3); % SA Good Percentage of world export of good
SARCAMat(:,5) = SARCAMat(:,2) ./ TotalSA;  % SA Good Percentage of SA exports
SARCAMat(:,6) = SARCAMat(:,3) ./ TotalWorld;  % Good Percentage of world exports
SARCAMat(:,7) = SARCAMat(:,5) ./ SARCAMat(:,6); %SA RCA

% xlswrite('ProductSum.xlsx',ProductSum)
% xlswrite('SARCA20146D.xlsx',CountryRCAMat)

Progress = 'TransformCountry_Finish'


end