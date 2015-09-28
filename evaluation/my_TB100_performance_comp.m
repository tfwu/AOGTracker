function my_performance_comp
%% modified from tarcker_benchmark_v1.0 (http://cvlab.hanyang.ac.kr/tracker_benchmark/benchmark_v10.html)

close all;

% save dir
saveDir = '/home/matt/Data/Results/TB100/Plots/';
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% tracker result dir

% where resultDir_AOGTracker/baseline/*.mat files (for AOGBaseline) are generated by AOGTracker_genPerfMat.m
%       resultDir_AOGTracker/full/*.mat files (for AOGTracker) are generated by AOGTracker_genPerfMat.m
resultDir_AOGTracker = '/home/matt/Data/Results/TB100/AOGTracker/'; 

% where .mat files are generated by genPerfMat.m in tarcker_benchmark_v1.0 (current version)
resultDir_Others = '/home/matt/Data/Results/TB100/OtherTrackers/';  

% other trackers to be compared
numOtherTrackers = 29;

% attribute annotations
[seqs, attrNames, attrAnnotation] = my_seq_attr;

evalTypes = {'OPE', 'TRE', 'SRE'};
metricTypes = {'error', 'overlap'};
rankTypes = {'AUC', 'threshold'};
for t = 1 : 3
    evalType = evalTypes{t};
    for subsetIdx = 0 : 2  % 0: TB100, 1: TB50, 2: TB-CVPR2013
        [seqNames, note] = getUsedSeqNames(subsetIdx);
        curSaveDir = [saveDir note filesep];
        if ~exist(curSaveDir,'dir')
            mkdir(curSaveDir);
        end
        
        % metric
        for m = 1 : 2
            metricType = metricTypes{m};
            % AOGTracker results
            resultFile_AOG_Baseline = [resultDir_AOGTracker 'baseline/aveSuccessRatePlot_AOGTracker_' metricType '_' evalType  '_' note '.mat'];
            resultFile_AOG = [resultDir_AOGTracker 'full/aveSuccessRatePlot_AOGTracker_' metricType '_' evalType  '_' note '.mat'];
            if ~exist(resultFile_AOG, 'file') 
                fprintf('can not find result files %s\n', resultFile_AOG);
                continue;
            end
            
            % Results of other trackers
            resultFile_Other = [resultDir_Others filesep note filesep 'aveSuccessRatePlot_' num2str(numOtherTrackers) 'alg_' metricType '_' evalType  '.mat'];
            
            % ranking
            for r = 1 : 2
                rankType = rankTypes{r};
                if strcmp(metricType,'error') && strcmp(rankType,'AUC')
                    continue;
                end               
                switch subsetIdx
                    case 0
                        figTitle = [evalType '-100(sequence average)'];
                    case 1
                        figTitle = [evalType '-50(sequence average)'];
                    case 2
                        figTitle = [evalType '-CVPR2013(51, sequence average)'];
                end
                figSaveName= [curSaveDir  evalType '_' metricType '_' rankType '_All_quality_plot' ];
                idxSeqs = 1 : length(seqNames);
                showComp(resultFile_AOG, resultFile_AOG_Baseline, resultFile_Other, seqNames, idxSeqs, ...
                    metricType, rankType, figTitle, figSaveName);
                
                % different attr
                for a = 1 : length(attrNames)
                    idxSeqsTmp = find(attrAnnotation(:, a));
                    if  length(idxSeqsTmp) < 2
                        continue;
                    end
                    
                    selSeqs = seqs(idxSeqsTmp);
                    idxSeqs = [];
                    for i = 1 : length(selSeqs)
                        idxSeqs = [idxSeqs, find(strcmp(selSeqs{i}, seqNames))];
                    end
                    
                    figSaveName= [curSaveDir  evalType '_' metricType '_' rankType '_' attrNames{a} ];
                    switch metricType
                        case 'overlap'
                            figTitle = ['Success plots of ' evalType ' - ' attrNames{a} ' (' num2str(length(idxSeqs)) ')'];
                        case 'error'
                            figTitle = ['Precision plots of ' evalType ' - ' attrNames{a} ' (' num2str(length(idxSeqs)) ')'];
                    end
                    showComp(resultFile_AOG, resultFile_AOG_Baseline, resultFile_Other, seqNames, idxSeqs, ...
                        metricType, rankType, figTitle, figSaveName);
                end % a
            end % r
        end % m
    end % subsetIdx
end % t


function showComp(resultFile_AOG, resultFile_AOG_Baseline, resultFile_Other, seqNames, idxSeqs, ...
    metricType, rankType, figTitle, figSaveName)

thresholdSetOverlap = 0:0.05:1;
thresholdSetError = 0:50;

switch metricType
    case 'overlap'
        thresholdSet = thresholdSetOverlap;
        rankIdx = 11;
        xLabelName = 'thresholds';
        yLabelName = 'Success rate';
    case 'error'
        thresholdSet = thresholdSetError;
        rankIdx = 21;
        xLabelName = 'thresholds';
        yLabelName = 'Precision';
end

% get results
trackerNames = {'AOGTracker'};

tmp = load(resultFile_AOG);
seqOrderAOGTracker = lower(tmp.usedSeqNames);
for i = 1 : length(idxSeqs)
    j = find(strcmp(seqNames{idxSeqs(i)}, seqOrderAOGTracker));
    try
        aveSuccessRatePlot(1, i, :) = tmp.aveSuccessRatePlot(j, :);
    catch
        found = 1;
    end
end

rowIdx = 2;
if exist(resultFile_AOG_Baseline,'file')
    tmp = load(resultFile_AOG_Baseline);
    seqOrderAOGTracker = lower(tmp.usedSeqNames);
    for i = 1 : length(idxSeqs)
        j = find(strcmp(seqNames{idxSeqs(i)}, seqOrderAOGTracker));
        try
            aveSuccessRatePlot(2, i, :) = tmp.aveSuccessRatePlot(j, :);
        catch
            found = 1;
        end
    end
    rowIdx = 3;
    trackerNames =[trackerNames; 'AOGBaseline'];
end

if exist(resultFile_Other, 'file')
    tmp = load(resultFile_Other);
    result = tmp.aveSuccessRatePlot(:, idxSeqs, :);
    aveSuccessRatePlot(rowIdx:size(result,1)+rowIdx-1, :, :) = result;
    trackerNames = [trackerNames; tmp.nameTrkAll];
end

for idxTrk = 1 : length(trackerNames)
    %each row is the sr plot of one sequence
    tmp=aveSuccessRatePlot(idxTrk, :, :);
    aa=reshape(tmp,[length(idxSeqs),size(aveSuccessRatePlot,3)]);
    aa=aa(sum(aa,2)>eps,:);
    bb=mean(aa);
    switch rankType
        case 'AUC'
            perf(idxTrk) = mean(bb);
        case 'threshold'
            perf(idxTrk) = bb(rankIdx);
    end
end

[~,indexSort]=sort(perf,'descend');

% plot
rankNum = 10; % show top 10 tracker in colors
fontSize = 16;
fontSizeLegend = 10;

% plot settings
plotDrawStyleAll = {...
    struct('color',[1,0,0],'lineStyle','--'),...
    struct('color',[0,1,0],'lineStyle','-'),...
    struct('color',[0,0,1],'lineStyle','-'),...
    struct('color',[0,0,0],'lineStyle','--'),...%    struct('color',[1,1,0],'lineStyle','-'),...%yellow
    struct('color',[1,0,1],'lineStyle','-'),...%pink
    struct('color',[0,1,1],'lineStyle','--'),...
    struct('color',[136,0,21]/255,'lineStyle','-'),...%dark red
    struct('color',[255,127,39]/255,'lineStyle','--'),...%orange
    struct('color',[0,162,232]/255,'lineStyle','-'),...%Turquoise
    struct('color',[163,73,164]/255,'lineStyle','--'),...%purple    %%%%%%%%%%%%%%%%%%%%
    struct('color',[0.5,0.5,0.5],'lineStyle','-.'),...%gray-25%
    };

fig = figure;
axes1 = axes('Parent',fig,'FontSize',14);

i=1;
for idxTrk=indexSort
    tmp=aveSuccessRatePlot(idxTrk,:,:);
    aa=reshape(tmp,[length(idxSeqs),size(aveSuccessRatePlot,3)]);
    aa=aa(sum(aa,2)>eps,:);
    bb=mean(aa) .* 100;
    
    switch rankType
        case 'AUC'
            score = mean(bb);
            tmp=sprintf('%.2f', score);
            loc = 'NorthEast';
        case 'threshold'
            score = bb(rankIdx);
            tmp=sprintf('%.2f', score);
            loc = 'SouthEast';
    end
    
    
    if i <= rankNum
        h(i) = plot(thresholdSet,bb,'color',plotDrawStyleAll{i}.color,...
            'lineStyle', plotDrawStyleAll{i}.lineStyle,'lineWidth', 4,'Parent',axes1);
        tmpName{i} = [trackerNames{idxTrk} ' [' tmp ']'];
    else
        h(i) = plot(thresholdSet,bb,'color',[0.5,0.5,0.5],...
            'lineStyle', '-.','lineWidth', 1,'Parent',axes1);
    end
    hold on
    i=i+1;
end

legend1=legend(tmpName,'Interpreter', 'none','fontsize',fontSizeLegend, ...
    'Location',loc);
title(figTitle,'fontsize',fontSize);
xlabel(xLabelName,'fontsize',fontSize);
% ylabel(yLabelName,'fontsize',fontSize);

grid on;
hold off;
set(gcf, 'PaperPosition', [0 0 7.78,5.83]);
set(gcf, 'PaperSize', [7.78,5.83]);
% saveas(gcf,figSaveName,'png');
saveas(gcf,figSaveName,'pdf');

close(fig);

function [seqNames, note] = getUsedSeqNames(subsetIdx)

% names in seqNames must be in the order other trackers in TB100 used to
% generate results.

if subsetIdx == 0
    seqNames = {'Basketball','Biker','Bird1','Bird2','BlurBody',...
        'BlurCar1','BlurCar2','BlurCar3','BlurCar4','BlurFace','BlurOwl',...
        'Board','Bolt','Bolt2','Box','Boy','Car1','Car2','Car24','Car4',...
        'CarDark','CarScale','ClifBar','Coke','Couple','Coupon','Crossing',...
        'Crowds','Dancer','Dancer2','David','David2','David3','Deer',...
        'Diving','Dog','Dog1','Doll','DragonBaby','Dudek','FaceOcc1',...
        'FaceOcc2','Fish','FleetFace','Football','Football1','Freeman1','Freeman3',...
        'Freeman4','Girl','Girl2','Gym','Human2','Human3','Human4.2','Human5',...
        'Human6','Human7','Human8','Human9','Ironman','Jogging.1','Jogging.2','Jump','Jumping',...
        'KiteSurf','Lemming','Liquor','Man','Matrix','Mhyang','MotorRolling',...
        'MountainBike','Panda','RedTeam','Rubik','Shaking','Singer1','Singer2',...
        'Skater','Skater2','Skating1','Skating2.1','Skating2.2','Skiing','Soccer','Subway',...
        'Surfer','Suv','Sylvester','Tiger1','Tiger2','Toy','Trans','Trellis',...
        'Twinnings','Vase','Walking','Walking2','Woman'};
    note = 'TB-100';
elseif subsetIdx == 1
    seqNames = {'Basketball', 'Biker', 'Bird1', 'BlurBody', 'BlurCar2', 'BlurFace', ...
        'BlurOwl', 'Bolt', 'Box','Car1','Car4','CarDark','CarScale','ClifBar','Couple', ...
        'Crowds','David','Deer','Diving','DragonBaby','Dudek','Football',...
        'Freeman4','Girl','Human3','Human4.2','Human6','Human9','Ironman','Jump','Jumping',...
        'Liquor','Matrix','MotorRolling','Panda','RedTeam','Shaking','Singer2','Skating1',...
        'Skating2.1','Skating2.2','Skiing','Soccer','Surfer','Sylvester','Tiger2','Trellis','Walking',...
        'Walking2','Woman' };
    note = 'TB-50';
elseif subsetIdx == 2
    seqNames = {'cardark', 'car4', 'david', 'david2', 'sylvester',...
        'trellis', 'fish', 'mhyang', 'soccer', 'matrix',...
        'ironman', 'deer', 'skating1', 'shaking', 'singer1', ...
        'singer2', 'coke', 'bolt', 'boy', 'dudek', 'crossing',...
        'couple', 'football1', 'jogging.1', 'jogging.2', 'doll',...
        'girl', 'walking2', 'walking', 'fleetface', 'freeman1',...
        'freeman3', 'freeman4', 'david3', 'jumping', 'carscale',...
        'skiing', 'dog1', 'suv', 'motorrolling', 'mountainbike',...
        'lemming', 'liquor', 'woman', 'faceocc1', 'faceocc2','basketball'...
        'football', 'subway', 'tiger1', 'tiger2'};
    note = 'TB-CVPR13';
end

seqNames = lower(seqNames);