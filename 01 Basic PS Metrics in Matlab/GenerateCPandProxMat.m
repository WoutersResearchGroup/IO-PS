function [CP,Prox,Centrality] = GenerateCPandProxMat(M,Products,Countries)

Progress = 'GenerateCPandProxMat_Start'

%This code generates:
%1)CP(i,j) -> the conditional probability matrix, . 
%       It indicates the probability of an RCA in j given an RCA in i. 
%       Its size is CP(Products,Products)

%2)Prox(i,j) -> the matrix of product proximities

%3)Centrality(i,2) a Matrix of ProductCentrality per product (1)
%HsCode; 2) Centrality )

%To Run the code it requires:
%1) M matrix (with 1's and 0's) 
%2) Products 
%3) Countries

CP = zeros(size(Products,1),size(Products,1)); %Instantiate CP
  
c = 1; %country number
p1 = 1; %product 1 (i) number
p2 = 1; %product 2 (j) number
CountP1 = 0;
CountP2GivenP1 = 0;

for p1 = 1:size(Products,1) %product 1 (i) number
    
    for p2 = 1:size(Products,1) %product 2 (j) number
        
        CountP1 = 0;
        CountP2GivenP1 = 0;
        
        for c = 1:size(Countries,1) %run through production of all countries for product i and j
    
            if M(c,p1) == 1
               
                CountP1 = CountP1 + 1;
                
                if M(c,p2) == 1
                    
                    CountP2GivenP1 = CountP2GivenP1 + 1;
                    
                end
                
            end
            
        end
        
        CP(p1,p2) = CountP2GivenP1/CountP1;
        
    end
    
    p1;
    
end

Prox = CP; %Instantiate Prox equal to CP before minimum mirroring


for p1 = 1:size(Products,1) %product 1 (i) number
    
   for p2 = 1:size(Products,1) %product 2 (j) number
        
       if p1 ~= p2

            if Prox(p1,p2) > Prox(p2,p1)

                Prox(p1,p2) = Prox(p2,p1);

            else

                Prox(p2,p1) = Prox(p1,p2);

            end

        end

        
   end
    
   p1;
    
end

Centrality = zeros(size(Products,1),2); %Format = 1) Product; 2) Centrality

Centrality(:,1) = Products(:,1);

for i = 1:size(Products,1) %Run through all products

        Centrality(i,2) = (sum(Prox(i,:)) - 1)  / (size(Products,1) - 1) ; 

end
    
%dlmwrite('ProductCentrality6D.txt',ProductCentrality6D,'precision',10)

Progress = 'GenerateCPandProxMat_Finish'

end