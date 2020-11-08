function [CountryCompInd,ProductCompInd] = CalcComplexity(M,Products,Countries,NumIters)

Progress = 'CalcComplexity_Start'

%This code generates the complexity of countries and products

%Output:
%CountryComplexityIndex6D(countries,complexities)
%ProductComplexityIndex6D(Products,complexities)
%Countries
%Products6D

%To Run the code it requires:
%1) M2D matrix (with 1's and 0's) 
%2) Products2D 
%3) Countries

M = M; 

Kc0 = zeros(size(Countries,1),1); %Indicates the number of products each country produces (array of diversity)

for i = 1:size(Countries,1)
    
    Kc0(i) = sum(M(i,:)) ;
    
end

Kp0 = zeros(size(Products,1),1); %Indicates the number of countries that produce each product (array of ubiquity)

for i = 1:size(Products,1)
    
    Kp0(i) = sum(M(:,i)) ;
    
end


Kc1 = zeros(size(Countries,1),1); %Indicates the next iteration of complexity of each country (array of country complexity 1)

for i = 1:size(Countries,1) %Run through countries
    
    ComplexitySum = 0;
    
    for j = 1:size(Products,1) %Run through products
    
        ComplexitySum = ComplexitySum + M(i,j)*Kp0(j);
        
    end
    
    Kc1(i) = ( 1/Kc0(i) ) * ComplexitySum; 
    
    
end


Kp1 = zeros(size(Products,1),1); %Indicates the next iteration of complexity products (array of product complexity 1)

for i = 1:size(Products,1) %Run through Products
    
    ComplexitySum = 0;
    
    for j = 1:size(Countries,1) %Run through Countries
    
        ComplexitySum = ComplexitySum + M(j,i)*Kc0(j);
        
    end
    
    Kp1(i) = ( 1 / Kp0(i) ) * ComplexitySum; 
    
end


Kc2 = zeros(size(Countries,1),1); %Indicates the next iteration of complexity of each country (array of country complexity 1)
Kp2 = zeros(size(Products,1),1); %Indicates the next iteration of complexity products (array of product complexity 1)


for n = 2:NumIters  %Run through number of iterations as required
    
    KcLast = Kc1;
    KpLast = Kp1;
    
    for i = 1:size(Countries,1) %Run through countries

        ComplexitySum = 0;

        for j = 1:size(Products,1) %Run through products

            ComplexitySum = ComplexitySum + M(i,j)*Kp1(j);

        end

        Kc2(i) = ( 1 / Kc0(i) ) * ComplexitySum; 


    end
    
    Kc1 = Kc2;

    for i = 1:size(Products,1) %Run through Products

        ComplexitySum = 0;

        for j = 1:size(Countries,1) %Run through Countries

            ComplexitySum = ComplexitySum + M(j,i)*Kc1(j);

        end

        Kp2(i) = ( 1 / Kp0(i) ) * ComplexitySum; 

    end
    
    Kp1 = Kp2;
    
    n;

end
KcN = Kc2;
KpN = Kp2;

CountryCompInd = zeros(size(Countries,1),1);
ProductCompInd = zeros(size(Products,1),1);

for i = 1:size(Countries,1)

    CountryCompInd(i) = ( KcN(i) - mean(KcN(:)) ) / std(KcN(:));

end

for i = 1:size(Products,1)

    ProductCompInd(i) = ( KpN(i) - mean(KpN(:)) ) / std(KpN(:));

end


% dlmwrite('CountryComplexity6D.txt',CountryComplexityIndex6D,'precision',10)
% dlmwrite('ProductComplexity6D.txt',ProductComplexityIndex6D,'precision',10)
% dlmwrite('Products6D.txt',Products6D,'precision',6)
% dlmwrite('Countries.txt',Countries)

Progress = 'CalcComplexity_Finish'

end