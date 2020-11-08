function IOPS_Figure()

Round1Metrics = readmatrix('Round1Metrics.txt'); % [CompContr OpporGainScenarios1 DistReq1];

BaselineMetrics = readmatrix('BaselineMetrics.txt'); % 1) Opportunity Value; 2) Average complexity of products within the value chain with RCA > 1; 3) Avg distance to Products in the VC with RCA < 1; 4) Average complexity of products within the PS with RCA > 1; 5) Average distance to products in the product space with RCA <1; 6) Num in GVC; 7) Sum of complexity in GVC; 8) Num in PS; 9) Sum of complexity in PS

AvgComp = Round1Metrics(:,1);
OpporGain = Round1Metrics(:,2);
Dist = Round1Metrics(:,3);

DeltaComp = AvgComp - BaselineMetrics(2);

FigTry1 = figure('Name','Distance versus complexity change and opportunity gain.','NumberTitle','off');
hold on;

Category = zeros(size(Dist));
PositiveOpporGain = 0.00001*ones(size(Dist));

Dist_NorthWest = NaN(size(Category));
DeltaComp_NorthWest = NaN(size(Category));

Dist_NorthEast = NaN(size(Category));
DeltaComp_NorthEast = NaN(size(Category));

Dist_West = NaN(size(Category));
DeltaComp_West = NaN(size(Category));

Dist_SouthEast = NaN(size(Category));
DeltaComp_SouthEast = NaN(size(Category));

Dist_East = NaN(size(Category));
DeltaComp_East = NaN(size(Category));

Dist_SouthWest = NaN(size(Category));
DeltaComp_SouthWest = NaN(size(Category));



DataForCompPareto = [Dist,-1*DeltaComp];
DataForOpporPareto = [Dist,-1*OpporGain];


for i = 1:size(Category,1)
    
%     if DataForCompPareto(i,2) > 0
%         DataForCompPareto(i,:) = NaN;
%     end
%     
%     if DataForOpporPareto(i,2) > 0
%         DataForOpporPareto(i,:) = NaN;
%     end
    
    
    
    if i == 1 || i == 2 || i == 14 || i == 35
       Dist(i) =  NaN;
       DeltaComp(i) = NaN;
%       PositiveOpporGain(i)  = NaN;
        
    end
    
    Category(i) = i; 
    if OpporGain(i) > 0
        PositiveOpporGain(i) = OpporGain(i);
    end
    
    
    if  i == 22 || i == 23 || i == 38 || i == 11 
       Dist_NorthWest(i) =  Dist(i);
       DeltaComp_NorthWest(i) = DeltaComp(i);
       
    elseif  i == 5 || i == 8
       Dist_West(i) =  Dist(i) - 0.05;
       DeltaComp_West(i) = DeltaComp(i);
       
    elseif i == 32 || i == 40
       Dist_SouthEast(i) =  Dist(i);
       DeltaComp_SouthEast(i) = DeltaComp(i);
       
    elseif i == 3 || i == 7  || i == 28 || i == 44
       Dist_East(i) =  Dist(i) + 0.05;
       DeltaComp_East(i) = DeltaComp(i); 
       
    elseif i == 41
       Dist_SouthWest(i) =  Dist(i);
       DeltaComp_SouthWest(i) = DeltaComp(i);
       
    else
       Dist_NorthEast(i) =  Dist(i);
       DeltaComp_NorthEast(i) = DeltaComp(i);
        
    end
    
    
    
end

% CompParetoSet = paretoQS(DataForCompPareto);
% OpporParetoSet = paretoQS(DataForOpporPareto);

CompParetoSet = [4 6 12 13 22 24 26 30 34];
OpporParetoSet = [4 6 10 12 13 30 31 32 37 45]; 

DotSizes = PositiveOpporGain * 100;

scatter(Dist,DeltaComp,DotSizes,'x','b');
scatter(Dist(CompParetoSet),DeltaComp(CompParetoSet),DotSizes(CompParetoSet),'s','k');
scatter(Dist(OpporParetoSet),DeltaComp(OpporParetoSet),DotSizes(OpporParetoSet),'o','m');


Cat2 = num2cell(Category);

text(Dist_NorthEast,DeltaComp_NorthEast,Cat2,'VerticalAlignment','bottom','HorizontalAlignment','left');%'VerticalAlignment','bottom'
text(Dist_NorthWest,DeltaComp_NorthWest,Cat2,'VerticalAlignment','bottom','HorizontalAlignment','right'); %,'VerticalAlignment','top'
text(Dist_West,DeltaComp_West,Cat2,'VerticalAlignment','middle','HorizontalAlignment','right'); %,'VerticalAlignment','top'
text(Dist_SouthEast,DeltaComp_SouthEast,Cat2,'VerticalAlignment','top','HorizontalAlignment','left'); %,'VerticalAlignment','top'
text(Dist_East,DeltaComp_East,Cat2,'VerticalAlignment','middle','HorizontalAlignment','left'); %,'VerticalAlignment','top'
text(Dist_SouthWest,DeltaComp_SouthWest,Cat2,'VerticalAlignment','top','HorizontalAlignment','right'); %,'VerticalAlignment','top'

%h = labelpoints(Dist,DeltaComp,Cat2,'E',0.15,1);

xlabel('Sum of distances to products in category');
ylabel('Change in average value chain complexity');

legend({'Categories','Complexity pareto set','Opportunity gain pareto set'},'Location','southeast')