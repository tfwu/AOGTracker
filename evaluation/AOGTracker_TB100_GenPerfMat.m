function AOGTracker_genPerfMat(resultFileFolder, perfMatPath)
%% modified from tarcker_benchmark_v1.0 (http://cvlab.hanyang.ac.kr/tracker_benchmark/benchmark_v10.html)

if nargin < 1
    resultFileFolder = '/home/matt/Data/TB100/';  % change them accordingly
    perfMatPath = '/home/matt/Data/TB100/'; 
end

evalTypes = {'OPE', 'TRE', 'SRE'};
for t = 1 : 3
    evalType = evalTypes{t};
    % "full" corresponds to the "note" in tracker_config_release.xml, e.g. "full" for AOGTracker, "baseline" for AOGBaseline
    resultFile = [resultFileFolder 'AOGTracker_Result_' evalType '_full.txt']; 
    for subsetIdx = 0 : 2  % 0: TB100, 1: TB50, 2: TB-CVPR2013
        if subsetIdx == 0
            note = '_TB-100';
        elseif subsetIdx == 1
            seqNames = {'Basketball', 'Biker', 'Bird1', 'BlurBody', 'BlurCar2', 'BlurFace', ...
                'BlurOwl', 'Bolt', 'Box','Car1','Car4','CarDark','CarScale','ClifBar','Couple', ...
                'Crowds','David','Deer','Diving','DragonBaby','Dudek', 'Football',...
                'Freeman4','Girl','Human3','Human4.2','Human6','Human9','Ironman','Jump','Jumping',...
                'Liquor','Matrix','MotorRolling','Panda','RedTeam','Shaking','Singer2','Skating1',...
                'Skating2.1','Skating2.2','Skiing','Soccer','Surfer','Sylvester','Tiger2','Trellis','Walking',...
                'Walking2','Woman' };
            seqNames = lower(seqNames);
            note = '_TB-50';
        elseif subsetIdx == 2
            seqNames = {'cardark', 'car4', 'david', 'david2', 'sylvester',...
                'trellis', 'fish', 'mhyang', 'soccer', 'matrix',...
                'ironman', 'deer', 'skating1', 'shaking', 'singer1', ...
                'singer2', 'coke', 'bolt', 'boy', 'dudek', 'crossing',...
                'couple', 'football1', 'jogging.1', 'jogging.2', 'doll',...
                'girl', 'walking2', 'walking', 'fleetface', 'freeman1',...
                'freeman3', 'freeman4', 'david3', 'jumping', 'carscale',...
                'skiing', 'dog1', 'suv', 'motorrolling', 'mountainbike',...
                'lemming', 'liquor', 'woman', 'faceocc1', 'faceocc2','basketball',...
                'football', 'subway', 'tiger1', 'tiger2'};            
            seqNames = lower(seqNames);
            note = '_TB-CVPR13';
        end
        
        % each row in the resultFile is
        % seqName \t seqResultFile \t gtFile \t startFrameIdx0 \t startFrameIdx \t endFrameIdx \t timeConsumed \t speed \t shiftType(TRE)
        if ~exist(resultFile, 'file')
            return;
        end
        fid = fopen(resultFile);
        allRes = textscan(fid, '%s %s %s %d %d %d %f %f %s', 'Delimiter', '\t');
        fclose(fid);
        numResult = length(allRes{1});          
        
        [~, sIdx] = sort(allRes{3});        
        
        thresholdSetOverlap = 0:0.05:1;
        thresholdSetError = 0:50;
        
        idxSeq = 1;
        idxRes = 1;        
        aveSuccessRatePlot = [];
        aveSuccessRatePlotErr = [];
        usedSeqNames = [];
        while idxRes <= numResult
            seqName = allRes{1}{sIdx(idxRes)};
            gtFile = allRes{3}{sIdx(idxRes)};
            idxResNext = idxRes + 1;
            while idxResNext <= numResult && strcmp(seqName, allRes{1}{sIdx(idxResNext)}) && strcmp(gtFile, allRes{3}{sIdx(idxResNext)})
                idxResNext = idxResNext + 1;
            end
            idxNum = idxResNext - idxRes;
            s = struct;
            for i = 1 : idxNum
                s(i).name = allRes{1}{sIdx(idxRes)};
                s(i).resultFile = allRes{2}{sIdx(idxRes)};
                s(i).gtFile = allRes{3}{sIdx(idxRes)};
                s(i).startFrame0 = allRes{4}(sIdx(idxRes));
                s(i).startFrame = allRes{5}(sIdx(idxRes));
                s(i).endFrame = allRes{6}(sIdx(idxRes)) - 1;
                s(i).timeConsumed = allRes{7}(sIdx(idxRes));
                s(i).speed = allRes{8}(sIdx(idxRes));
                s(i).shift = allRes{9}{sIdx(idxRes)};
                idxRes = idxRes + 1;
            end
            
            [~, gtBasename, ~] = fileparts(gtFile);
            idxTmp = strfind(gtBasename, '.');
            if ~isempty(idxTmp)
                seqName1 = [seqName '.' gtBasename(idxTmp(1)+1:end)];
            else
                seqName1 = seqName;
            end
            
            if subsetIdx > 0
                if isempty(strmatch(lower(seqName1), seqNames, 'exact'))
                    continue;
                end
            end
            
%             disp([num2str(idxSeq) ' ' seqName1 ' ' num2str(idxNum)]);            
            usedSeqNames{idxSeq} = seqName1;
            
            rect_anno = dlmread(s(1).gtFile);
            
            lenALL = 0;
            
            successNumOverlap = zeros(idxNum,length(thresholdSetOverlap));
            successNumErr = zeros(idxNum,length(thresholdSetError));
            
            for idx = 1:idxNum
                res = dlmread(s(idx).resultFile);
                res(:, 1) = res(:, 1) + 1;
                res(:, 2) = res(:, 2) + 1;
                istart = s(idx).startFrame - s(idx).startFrame0 + 1;
                iend = min(size(rect_anno, 1), s(idx).endFrame - s(idx).startFrame0 + 1);
                anno = rect_anno(istart : iend, :);
                len = size(anno,1);
                
                res = res(1:iend-istart+1, :);
                [aveCoverage, aveErrCenter, errCoverage, errCenter] = calcSeqErr(res, anno);
                
                for tIdx=1:length(thresholdSetOverlap)
                    successNumOverlap(idx,tIdx) = sum(errCoverage >thresholdSetOverlap(tIdx));
                end
                
                for tIdx=1:length(thresholdSetError)
                    successNumErr(idx,tIdx) = sum(errCenter <= thresholdSetError(tIdx));
                end
                
                lenALL = lenALL + len;
            end
            
            if strcmp(evalType, 'OPE')
                aveSuccessRatePlot(idxSeq,:) = successNumOverlap/(lenALL+eps);
                aveSuccessRatePlotErr(idxSeq,:) = successNumErr/(lenALL+eps);
            else
                aveSuccessRatePlot(idxSeq,:) = sum(successNumOverlap)/(lenALL+eps);
                aveSuccessRatePlotErr(idxSeq,:) = sum(successNumErr)/(lenALL+eps);
            end
            idxSeq = idxSeq + 1;
        end
        
        %
        dataName1=[perfMatPath 'aveSuccessRatePlot_AOGTracker_overlap_' evalType  note '.mat'];
        save(dataName1,'aveSuccessRatePlot', 'usedSeqNames');
        
        aa=aveSuccessRatePlot;
        aa=aa(sum(aa,2)>eps,:);
        bb=mean(aa);
        disp(['****** Overlap: numSeq ' num2str(idxSeq-1) ' ' evalType note ': ' num2str(mean(bb))]);
               
        
        dataName2=[perfMatPath 'aveSuccessRatePlot_AOGTracker_error_' evalType note '.mat'];
        aveSuccessRatePlot = aveSuccessRatePlotErr;
        save(dataName2,'aveSuccessRatePlot', 'usedSeqNames');
                
        aa=aveSuccessRatePlot;
        aa=aa(sum(aa,2)>eps,:);
        bb=mean(aa);
        disp(['       Error: numSeq ' num2str(idxSeq-1) ' ' evalType note ': ' num2str(mean(bb))]);      
        
    end
end



function [aveErrCoverage, aveErrCenter,errCoverage, errCenter] = calcSeqErr(results, rect_anno)

seq_length = size(results, 1);

for i = 2:seq_length
    r = results(i,:);
    r_anno = rect_anno(i,:);
    if (isnan(r) | r(3)<=0 | r(4)<=0) & (~isnan(r_anno))
        results.res(i,:)=results.res(i-1,:);
    end
end

centerGT = [rect_anno(:,1)+(rect_anno(:,3)-1)/2 rect_anno(:,2)+(rect_anno(:,4)-1)/2];

rectMat = results;

rectMat(1,:) = rect_anno(1,:);

center = [rectMat(:,1)+(rectMat(:,3)-1)/2 rectMat(:,2)+(rectMat(:,4)-1)/2];

errCenter = sqrt(sum(((center(1:seq_length,:) - centerGT(1:seq_length,:)).^2),2));

index = rect_anno > 0;
idx=(sum(index,2)==4);
tmp = calcRectInt(rectMat(idx,:),rect_anno(idx,:));

errCoverage=-ones(length(idx),1);
errCoverage(idx) = tmp;
errCenter(~idx)=-1;

aveErrCoverage = sum(errCoverage(idx))/length(idx);

aveErrCenter = sum(errCenter(idx))/length(idx);

function overlap = calcRectInt(A,B)
%
%each row is a rectangle.
% A(i,:) = [x y w h]
% B(j,:) = [x y w h]
% overlap(i,j) = area of intersection
% normoverlap(i,j) = overlap(i,j) / (area(i)+area(j)-overlap)
%
% Same as built-in rectint, but faster and uses less memory (since avoids repmat).


leftA = A(:,1);
bottomA = A(:,2);
rightA = leftA + A(:,3) - 1;
topA = bottomA + A(:,4) - 1;

leftB = B(:,1);
bottomB = B(:,2);
rightB = leftB + B(:,3) - 1;
topB = bottomB + B(:,4) - 1;

tmp = (max(0, min(rightA, rightB) - max(leftA, leftB)+1 )) .* (max(0, min(topA, topB) - max(bottomA, bottomB)+1 ));
areaA = A(:,3) .* A(:,4);
areaB = B(:,3) .* B(:,4);
overlap = tmp./(areaA+areaB-tmp);

