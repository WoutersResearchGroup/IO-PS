function [CompInPS,DistToPS,NumberWithRCA,SumCompInPS] = PSCompAndDist(DistanceAndOpporGain,ProductCompInd,SARCAMat)

%1) DistanceAndOpporGain %1) HsCode; 2) Distance; 3) Distance if opportunity; 4) OpporGain; 5) OpporGain if Opportunity 6) Density; 7) Density if Oppor
%2) ProductCompInd
%6) SARCAMat 1) hs92code; 2) Export of Country; 3) World Export of product; 4) SA Good Percentage of world export of good; 5) SA Good Percentage of SA exports; 6)  Good Percentage of world exports; 7) SA RCA

SumCompInPS = 0; 
NumberWithRCA = 0;

SumDistToPS = 0;
NumberCurOutPS = 0;

for i = 1:size(ProductCompInd,1)
    
    if SARCAMat(i,7) > 1 %Calculate avg complexity of products within the PS with RCA > 1
        
        SumCompInPS = SumCompInPS + ProductCompInd(i);
        NumberWithRCA = NumberWithRCA + 1;
         
    else % Calculate avg Distance to products in the PS with RCA < 1
        
        SumDistToPS = SumDistToPS + DistanceAndOpporGain(i,2);
        NumberCurOutPS = NumberCurOutPS + 1;
        
    end        
        
end

CompInPS = SumCompInPS / NumberWithRCA;
DistToPS = SumDistToPS / NumberCurOutPS;

end