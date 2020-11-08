function [M,MAbs,Countries] = generateM(Products,AllMat,AllCA)

Progress = 'generateM_Start'

%This code generates the M and Mabs (with the actual RCA)

%IT requires the loading of:
%1) AllCA; Format = % Format = ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA) (It requires that the countries are in alphabetical Order)
%2) AllMat; % Format = ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA) (It requires that the countries are in alphabetical Order)
%3) hs92codes; 

CountriesDupl = AllCA(:,2); %Reads all country names into array (Sorted that names of countries are in alphabetical order)

Countries = unique(CountriesDupl); %Removes duplicates

ProductRCA = AllMat(:,[3 6]); %Filter only necessary parts of AllMat (Sorted that names of countries are in alphabetical order)

MAbs = zeros(size(Countries,1),size(Products,1)); %Instantiate MAbs
M = zeros(size(Countries,1),size(Products,1)); %Instantiate M 
  
c = 1; %country number

for k = 1:size(ProductRCA,1) %Run through all rows of All2014 array
    
    if k > 1
        if strcmp(CountriesDupl(k),CountriesDupl(k-1))
        else
            c=c+1;
        end
    end
   
        for j = 1:size(Products,1) %Run through the M array columns

            if Products(j) == ProductRCA(k,1)

                MAbs(c,j) = ProductRCA(k,2);

                if MAbs(c,j) > 1

                    M(c,j) = 1;

                end                 

            end
            
        end
    
end

Progress = 'generateM_Finish'
    
end