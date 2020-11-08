function [NumIters,Products,M,MAbs,Countries,CountryCompInd,ProductCompInd,CP,Prox,Centrality,DistanceAndOpporGain,Densities,Distance] = CallAllBasicPSMetrics(CountryLetters,CSVFileName)

MainTime = tic;

ReadInData = tic;

Mat = readmatrix(CSVFileName);
CA = readcell(CSVFileName,'Range',[2 1]);

Mat(isnan(Mat))=0;
indices = find(Mat(:,4)==0);
Mat(indices,:) = [];
CA(indices,:) = [];   

DataReadComplete = toc(ReadInData)

%Mat; % Format = ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA)
%CA;  % Format = ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA)
%CountryRaw;  % Format = ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA)

%% Declarations (Changeable inputs)
NumIters = 18; % Number of iterations required for the Method of Reflections

%% UpFront calculations:
Products = unique(Mat(:,3));

%% Tranform country data to create country RCA matrix and sum world production

%Outputs:
%SARCAMat   % Format = ( 1 hs92code; 2 Export of Country; 3 World Export of product; 4 SA Good Percentage of world export of good; 5 % SA Good Percentage of SA exports; 6  Good Percentage of world exports; 7 SA RCA)
%ProductSum % Format = 1 hs92 code; 2 World exports of product

%CountryRaw;  % Format = ( 1 Year; 2 Country; 3 hs92code; 4 Export Value; 5 Import value: 6 Export RCA; 7 Import RCA)

CountryRows = strcmp(CA(:,2), CountryLetters);
CountryEntries = sum(CountryRows);
CountryRaw = zeros(CountryEntries,7);
CountryRaw(:,1) = Mat(CountryRows,1); 
CountryRaw(:,[3 4]) = Mat(CountryRows,[3 4]); 

[CountryRCAMat,~] = TransformCountry(Mat,CountryRaw,Products); 


%% Generate M matrixes

[M,MAbs,Countries] = generateM(Products,Mat,CA);

%% Calculate complexity

%CountryCompInd 1) complexities)
%ProductCompInd 1) complexities
[CountryCompInd,ProductCompInd] = CalcComplexity(M,Products,Countries,NumIters);

%% Calculate conditional probability and proximity matrix

%This code generates:
%1)CP(i,j) -> the conditional probability matrix. 
%       It indicates the probability of an RCA in j given an RCA in i. 
%       Its size is CP(Products,Products)
%2)Prox(i,j) -> the matrix of product proximities
%3)Centrality(i,2) a Matrix of ProductCentrality per product (1)
%HsCode; 2) Centrality )

[CP,Prox,Centrality] = GenerateCPandProxMat(M,Products,Countries);

%% Calculate Opportunity Gain and Distance

%RCA   % Format = ( 1 hs92code; 2 SA RCA)

RCA = CountryRCAMat(:,[1 7]); 

%DistanceAndOpporGain 1) HsCode; 2) Distance; 3) Distance if opportunity; 4) OpporGain; 5) OpporGain if Opportunity 6) Density; 7) Density if Oppor
[DistanceAndOpporGain,Densities,Distance,ProxSums] = DistanceAndGain(Prox,RCA,Products,ProductCompInd);

writematrix(Products,'Products.xlsx');
writematrix(M,'M.xlsx');
writematrix(MAbs,'MAbs.xlsx');
writecell(Countries,'Countries.xlsx');
writematrix(CountryCompInd,'CountryCompInd.xlsx');
writematrix(ProductCompInd,'ProductCompInd.xlsx');
writematrix(CP,'CP.xlsx');
writematrix(Prox,'Prox.xlsx');
writematrix(Centrality,'Centrality.xlsx');
writematrix(DistanceAndOpporGain,'DistanceAndOpporGain.xlsx');
writematrix(Densities,'Densities.xlsx');
writematrix(Distance,'Distance.xlsx');
writematrix(ProxSums,'ProxSums.xlsx');


TotalTime = toc(MainTime)


end