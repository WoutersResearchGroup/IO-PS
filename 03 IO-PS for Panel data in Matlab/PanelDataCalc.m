function PanelDataCalc()

%% Set file name

MatCSVFileName = 'InputData.csv';

%% Extract Mat and CA

Mat = readmatrix(MatCSVFileName);
CA = readcell(MatCSVFileName,'Range',[2 1]);

Mat(isnan(Mat))=0;
indices = find(Mat(:,4)==0);
Mat(indices,:) = [];
CA(indices,:) = [];   

%% Set names for variable per year

%     Round1MetricsYear = sprintf('Round1Metrics%d',k);
%     BaseLineMetricsYear = sprintf('BaseLineMetrics%d',k);
%     GVCActYear = sprintf('GVCActYear%d',k);
%     GVCFullYear = sprintf('GVCFullYear%d',k);
%     GVCTierYear = sprintf('GVCTierYear%d',k);
%     
%% Calculate all variables for relevant year

[Round1Metrics,BaseLineMetrics,GVCAct,GVCFull,GVCTier] = IOPSwithCountrySelectPSBasedExSimulation(Mat,CA,'zaf');


%% Write out Round1MetricsYear

Round1MetricsOut = [Round1Headings;num2cell(Round1Metrics)];
writecell(Round1MetricsOut,sprintf('Round1Metrics%d.xlsx',k));
dlmwrite(sprintf('Round1Metrics%d.txt',k),Round1Metrics,'precision',10);

%% Write out BaseLineMetrics
BaseLineMetricsOut = [BaseLineMetricsHeadings;num2cell(transpose(BaseLineMetrics))];
dlmwrite(sprintf('BaseLineMetrics%d.txt',k),BaseLineMetrics,'precision',10);
writecell(BaseLineMetricsOut,sprintf('BaseLineMetrics%d.xlsx',k));

%% Write out GVC Tier Act and Full

GVCTierOut = [GVCTierHeadings;num2cell(GVCTier)];
GVCActOut = [GVCActHeadings;num2cell(GVCAct)];
GVCFullOut = [GVCFullHeadings;num2cell(GVCFull)];

writecell(GVCTierOut,sprintf('GVCTierResults%d.xlsx',k));
writecell(GVCActOut,sprintf('GVCActResults%d.xlsx',k));
writecell(GVCFullOut,sprintf('GVCFullResults%d.xlsx',k));
    
end




