classdef prtDataSetTuf < prtDataSetClass
    % prtDataSetTuf < prtDataSetClass
    %  An extention of prtDataSetClass that includes tools for dealing with
    %  landmine detection.
    %
    % Creation:
    %  If you have a tufresults structure in your workspace, 
    %       ds = prtDataSetTuf(tufresults)
    %  will generate the data set.  Alternatively:
    %       tufResultsMat = ...\TUFResults\...\201203_Yuma_STMR_Big__r2_ at 2012-04-13 13.27.59.mat';
    %       s = load(tufResultsMat);
    %       ds = prtDataSetTuf(s)
    %  or: 
    %       ds = prtDataSetTuf(tufResultsMat);
    %  work.
    %
    % Scoring:
    %   [far,pd,thresh,names] = ds.scoreRoc;
    %   or:
    %   ds.scoreRoc; %makes plots
    %
    % Removing clutter:
    %   ds = ds.removeClutter;
    %
    % Removing or retaining target types:
    %   ds = ds.applyObjectTypeExclusions({'VS1.6','VS2.2','155MM'}); 
    %   ds = ds.applyObjectTypeRetentions('*155*'); %note - wildcards work
    %
    % Processing individual lanes:
    %   uLanes = ds.getUniqueLanes;
    %   for i = 1:length(uLanes);
    %      scoreRoc(ds.retainLanes(uLanes{i}));
    %      title(uLanes{i});
    %      pause;
    %   end
    % 
    % Retain shallow targets:
    %   dsShallow = ds.retainObjectDepths(@(x)x <= 3);
    %   dsShallow.scoreRoc;
    %
    % Change labeling:
    %   ds = ds.maryizeTargets;
    %
    % Merge target names (which can be weirdly precise):
    %   ds = ds.mergeObjectTypes('*155*','155');
    %   ds = ds.maryizeTargets;
    %
    % Get data (requires TUF, tuf.locate, and MDR files):
    %   [gprData,ch] = ds.readGprData(10);
    
    properties
        excursionArray
    end
    
    properties (Hidden)
        %don't get cute with this if you don't know regexp.  WHen this is
        %true, *HOLE* matches what 99% of people think it should. 
        ezRegexp = true;
    end
    
    methods (Access = protected)
        function self = retainExcursions(self,retain)
            %self = retainExcursions(self,retain)
            % This is protected; you should get to this only through
            % retainLanes or retainLsdKvs
            self.excursionArray = self.excursionArray(retain);
            self.userData.scorableArea = sum([self.excursionArray.area]);
        end 
    end
    
    methods
        
        function self = eventicizeObservationInfo(self)
            
            for i = 1:self.nObservations
                if ~isempty(self.observationInfo(i).event_id)
                    baseTime = self.observationInfo(i).event_timestamp;
                    break;
                end
            end
            
            obsInfo = self.observationInfo;
            for i = 1:self.nObservations;
                if isempty(obsInfo(i).event_id)
                    x = round(obsInfo(i).xUTM*100);
                    y = round(obsInfo(i).yUTM*100);
                    
                    id = sprintf('%dE_%dN',x,y);
                    obsInfo(i).event_id = id;
                    obsInfo(i).event_type = 'object';
                    obsInfo(i).event_timestamp = baseTime;     
                end
            end
            self.observationInfo = obsInfo;
            
        end
        
        function featureStruct = locateFeatures(self,featureName)
            
            obsInfo = self.observationInfo;
            featureStruct(1) = tuf.get_features(featureName,obsInfo(1));
            featureStruct = repmat(featureStruct,self.nObservations,1);
            
            h = prtUtilProgressBar(0,'Extracting features','autoClose',true);
            
            if any(cellfun(@(c)isempty(c),{self.observationInfo.event_timestamp}))
                error('Some events appear to be empty... use ds = eventicizeObservationInfo(ds)');
            end
            step = 100;
            for i = 1:step:self.nObservations
                if ~mod(i-1,step)
                    h.update(i./self.nObservations);
                end
                inds = i:min([i+step,self.nObservations]);
                featureStruct(inds) = tuf.get_features(featureName,obsInfo(inds));
            end
            h.update(1);
            
        end
    
        function [n,p,summary] = plotDetectionSummary(self,targetFars,varargin)
            % [n,p,summary] = plotDetectionSummary(self,targetFars)
            %
            % load tempDs.mat ds
            % ds = ds.removeClutter;
            % ds = ds.hmdsStandardMergeObjectTypes;
            % ds = ds.maryizeTargets;
            % ds.plotDetectionSummary(.01);
            
            if self.nFeatures > 1
                fprintf('Using first feature\n');
                self = self.retainFeatures(1);
            end
            [far,pd,thresh,names] = scoreRoc(self.binarizeTargets); %#ok<NASGU,ASGLU>
            
            %convert FAR's to thresholds
            threshold = nan(size(targetFars));
            for farInd = 1:length(targetFars)
                rocInd(farInd) = find(far(:) < targetFars(farInd),1,'last');
                threshold(farInd) = thresh(rocInd(farInd));
            end
            
            summary = self.getDetectionSummary(threshold,varargin{:});
            
            [n,p] = prtDataSetTuf.summaryToMatrices(summary);
            
            close all;
            for i = 1:length(summary); 
                
                subplot(2,4,[1:3,5:7]); cla;
                localPie(1:size(p{1},2),1:size(p{1},1),p{i}'./n{i}',n{i}'); 
                set(gca,'ytick',1:length(summary(1).types));
                set(gca,'yticklabel',cellfun(@(s)sprintf('%s',s(1:min([7,length(s)]))),summary(1).types,'UniformOutput',false));
                set(gca,'xtick',(0:length(summary(1).depthBins))+.5);
                xlabels = cellfun(@(s)sprintf('%d',s),num2cell(summary(1).depthBins),'UniformOutput',false);
                set(gca,'xticklabel',xlabels);
                hold off;
                
                subplot(2,4,4); cla;
                h = plot(far,pd); xlim([0 .01]); ylim([0 1]);
                set(h,'linewidth',3);
                hold on; plot(far(rocInd(i)),pd(rocInd(i)),'k*');
                hold off;
                subplot(2,4,8); cla;
                h = plot(thresh,far./max(far),thresh,pd);
                labs = get(gca,'ytick');
                trueLab = labs.*max(far);
                set(gca,'yticklabel',arrayfun(@(x)sprintf('%.2d',x),trueLab,'UniformOutput',false));
                hold on; plot(thresh(rocInd(i)),far(rocInd(i))./max(far),'kx',thresh(rocInd(i)),pd(rocInd(i)),'kx');
                axis tight;
                tickOff;
                hold off;
                set(h,'linewidth',3);
                
                
                if length(summary) > 1
                    pause;
                    
                    if i < length(summary)
                        delete(gca);
                    end
                    %                     close all;
                end
            end
        end
        
        function summary = getDetectionSummary(self,varargin)
            p = inputparser;
            p.addRequired('threshold');
            p.addParamValue('depthBins',[-inf,1,8,12,inf]);
            p.parse(varargin{:});
            inputs = p.Results;
            
            tempSelf = self.removeFalseAlarms;
            if self.nFeatures > 1
                fprintf('Using first feature');
                tempSelf = tempSelf.retainFeatures(1);
            end
            tempSelf = tempSelf.removeFalseAlarms;
            tempSelf = tempSelf.removeClutter;
            
            uClasses = tempSelf.uniqueClasses;
            classNames = tempSelf.classNames;
                    
            depths = tempSelf.getAlarmObjectDepths;
            binnedDepths = sum(bsxfun(@lt,inputs.depthBins,depths(:)),2);
            globalDepths = binnedDepths;
            
            summary = struct;
            csummary = struct;
            for iThreshold = 1:length(inputs.threshold)
                cThreshold = inputs.threshold(iThreshold);
                
                for iType = 1:tempSelf.nClasses;
                    indices = tempSelf.targets == uClasses(iType);
                    dsType = tempSelf.retainObservations(indices);
                    
                    depths = dsType.getAlarmObjectDepths;
                    binnedDepths = sum(bsxfun(@lt,inputs.depthBins,depths(:)),2);
                    depths = binnedDepths;
                    
                    uTypeDepths = unique(depths);
                    found = nan(1,length(uTypeDepths));
                    total = found;
                    for depthInd = 1:length(uTypeDepths)
                        retain = depths == uTypeDepths(depthInd);
                        current = dsType.retainObservations(retain);
                        
                        found(depthInd) = sum(current.getX >= cThreshold);
                        total(depthInd) = current.nObservations;
                    end
                    csummary(iType).depths = uTypeDepths;
                    csummary(iType).depthBins = inputs.depthBins;
                    csummary(iType).foundPerDepth = found;
                    csummary(iType).totalPerDepth = total;
                    csummary(iType).found = sum(found);
                    csummary(iType).total = sum(total);
                    csummary(iType).class = classNames{iType};
                end
                summary(iThreshold).localSummary = csummary;
                summary(iThreshold).threshold = cThreshold;
                summary(iThreshold).globalDepths = unique(globalDepths);
                summary(iThreshold).depthBins = inputs.depthBins;
                summary(iThreshold).types = classNames;
            end
            
        end
        
        function self = sortAlarmsByConfidence(self)
            x = self.getX;
            [~,sortInd] = sort(x(:,1),'descend');
            self = self.retainObservations(sortInd);
        end
        
        function [data,ch,obsStruct] = readGprData(self,index,dtLength)
            if nargin == 2
                dtLength = 30;
            end
            dt = dukeTufAlarmToDtCh(self.observationInfo(index));
            
            dataFile = tuf.locate(self.observationInfo(index).lsd_kv,'data');
            start = max([dt-dtLength,1]);
            data = readMdrMex(dataFile,start,dt+dtLength);
            obsStruct = self.observationInfo(index);
            ch = self.observationInfo(index).ch;
        end
        
        
        function self = binarizeTargets(self)
            %self.targets = 1 for threats, 0 for all others
            isClutter = logical([self.observationInfo.isClutter]);
            isHit = logical([self.observationInfo.hit]);
            
            binaryTargets = ~isClutter(:) & isHit(:);
            self.targets = binaryTargets;
            self.classNames = {'Non-Threat','Threat'};
        end
        
        function self = trinarizeTargets(self)
            %self.targets = 1 for threats, -1 for clutter, 0 for all others
            isClutter = logical([self.observationInfo.isClutter]);
            isHit = logical([self.observationInfo.hit]);
            
            targets = zeros(self.nObservations,1);
            targets(isHit) = 1;
            targets(isClutter) = -1;
            self.targets = targets;
            self.classNames = {'Clutter','False Alarm','Threat'};
        end
        
        function self = maryizeTargets(self) 
            %self = maryizeTargets(self) 
            %   Use the target types to make class labels.  TYpically you
            %   may want to remove clutter, and/or mergeObjectTypes before
            %   this.
            
            types = self.getAlarmObjectTypes;
            [y,uTypes] = prtUtilStringsToClassNumbers(types);
            uTypes = uTypes(:);
            
            faInd = find(strcmpi(uTypes,'False Alarm'));
            nonFaInd = setdiff(1:length(uTypes),faInd);
            if ~isempty(faInd)
                y(y == faInd) = 0;
                y(y > faInd) = y(y > faInd)-1;
                uTypes = uTypes([faInd;nonFaInd(:)]);
            end
            
            
            
            self.targets = y;
            self.classNames = uTypes;
        end
        
        function self = hmdsStandardMergeObjectTypes(self)
            self = mergeObjectTypes(self,'VS2.2*','VS2.2');
            self = mergeObjectTypes(self,'*155*','155');
            self = mergeObjectTypes(self,'M15*','M15');
            self = mergeObjectTypes(self,'PP*','PP');
            self = mergeObjectTypes(self,'*JUG*','JUG');
            self = mergeObjectTypes(self,'*PAIL*','PAIL');
            self = mergeObjectTypes(self,'*PAIL*','PAIL');
        end
        
        
        function self = mergeObjectTypes(self,typeList,newType)
            %self = mergeObjectTypes(self,typeList,newType)
            % Combine all the object types matching typeList into a single
            % type - newType.
            %
            % For example, 
            %
            % dsEasy = ds.mergeObjectTypes('155*','155')
            %
            % Turns all 155*'s into 155s.  And
            %
            % dsEasy = ds.mergeObjectTypes({'155*','VS*'},'155')
            %
            % Turns all VS's and 155*s into 155's.  I'm not sure why you'd
            % ever want to do that, but you *could*.
            
            if ~isa(typeList,'cell')
                typeList = {typeList};
            end
            pat = prtDataSetTuf.strCellToRegexpPattern(typeList,self.ezRegexp);
            
            types = self.getAlarmObjectTypes;
            types = lower(types);
            matches = regexp(types,pat);
            matches = cellfun(@(c)~isempty(c),matches);
            foundMatches = find(matches);
            
            obsInfo = self.observationInfo;
            for i = 1:length(foundMatches)
                obsInfo(foundMatches(i)).targetId = newType;
            end
            self.observationInfo = obsInfo;
        end
        
        function types = getAlarmObjectTypes(self)
            types = {self.observationInfo.targetId};
            types = types(:);
        end
        
        function [self,retain] = retainObjectDepths(self,shallow,deep)
            % self = removeObjectDepths(self,shallow = -inf,deep = inf)
            % self = removeObjectDepths(self,fnHandle)
            depths = getAlarmObjectDepths(self);
            
            if nargin == 2 && isa(shallow,'function_handle')
                fnHandle = shallow;
                retain = fnHandle(depths);
            elseif nargin == 3 || nargin == 1
                if nargin == 1
                    shallow = -inf;
                    deep = inf;
                else
                    if isempty(shallow)
                        shallow = -inf;
                    end
                    if isempty(deep)
                        deep = inf;
                    end
                end
                retain = depths >= shallow & depths <= deep;
            end
            retain = retain | isnan(depths);
            self = self.retainObservations(retain);
        end
                
        
        function depths = getAlarmObjectDepths(self)
            %depths = getAlarmObjectDepths(self) 
            targets = cat(1,{self.observationInfo.objects});
            depths = nan(length(targets),1);
            for i = 1:length(targets)
                if ~isempty(targets{i})
                    depths(i) = targets{i}.depth;
                end
            end
        end
        
        function [self,removed] = hmdsApplyObjectTypeExclusions(self)
            %self = hmdsApplyObjectTypeExclusions(self)
            if self.ezRegexp
                exclusions = {'*WIRE*','*SIM*','*SAND*','*1GP_HME*'};
            else
                %just in case someone set this...
                exclusions = {'.*WIRE.*','.*SIM.*','.*SAND.*','.*1GP_HME.*'};
            end
            [self,removed] = applyObjectTypeExclusions(self,exclusions);
        end
        
        function [self,removed] = removeClutter(self)
            %[self,removed] = removeClutter(self)
            removed = [self.observationInfo.isClutter];
            removed = logical(removed);
            self = self.removeObservations(removed(:));
        end
        
        function [self,removed] = removeFalseAlarms(self)
            types = self.getAlarmObjectTypes;
            removed = strcmpi(types,'False Alarm');
            self = self.removeObservations(removed);
        end
        
        function [self,retained] = retainFalseAlarms(self)
            types = self.getAlarmObjectTypes;
            retained = strcmpi(types,'False Alarm');
            self = self.retainObservations(retained);
        end
        
        function [self,retained] = applyObjectTypeRetentions(self,retentions)
            %self = applyObjectTypeRetentions(self,exclusions)
            % See applyObjectTypeExclusions, but, the opposite.  
            
            if ~isa(retentions,'cell')
                retentions = {retentions};
            end
            
            %turn a cell array into a regexp pattern of or's
            pat = prtDataSetTuf.strCellToRegexpPattern(retentions,self.ezRegexp);
            
            types = self.getAlarmObjectTypes;
            if self.ezRegexp
                types = lower(types);
            end
            
            matches = regexp(types,pat);
            retained = cellfun(@(c)~isempty(c),matches);
            %make sure to retain these:
            retained = retained | strcmpi(types,'False Alarm');
            
            self = self.retainObservations(retained); %doesn't remove area... 
        end
        
        function [self,removed] = applyObjectTypeExclusions(self,exclusions)
            %self = applyObjectTypeExclusions(self,exclusions)
            %   Apply type exclusions to the data set.  exclusions should
            %   be a cell array of regexp-friendly strings.  If you don't
            %   know what that means, just use strings like normal - 
            %
            %   {'155','155M'} matches any object with the type 155 or
            %   155M, but not objects of the type 155MM, or 155MA.
            %
            %   {'155*'} matches 155MM and 155MMA, but not x155
            %
            %   {'*155*'} matches anything with 155 anywhere in it.
            %
            %  Note: applyObjectTypeExclusions will never exclude false
            %  alarms - i.e., objects with type matching "False Alarm"
            
            if ~isa(exclusions,'cell')
                exclusions = {exclusions};
            end
            
            %turn a cell array into a regexp pattern of or's
            pat = prtDataSetTuf.strCellToRegexpPattern(exclusions,self.ezRegexp);
            
            types = self.getAlarmObjectTypes;
            if self.ezRegexp
                types = lower(types);
            end
            
            matches = regexp(types,pat);
            removed = cellfun(@(c)~isempty(c),matches);
            %make sure to NOT remove these:  Maybe warn.
            if any(removed & strcmpi(types,'False Alarm'))
                warning('prtDataSetTuf:triedToExcludeFalseAlarms','Some exclusion patterns seem to match "False Alarm"; False Alarms are being retained');
            end
            removed = removed & ~strcmpi(types,'False Alarm');
            
            self = self.removeObservations(removed); %doesn't remove area... 
        end
        
        function self = insertMissingAlarmLocations(self)
            %self = insertMissingAlarmLocations(self)
            
            lsdKvs = self.getLsdKvs;
            uLsdKvs = unique(lsdKvs);
            obsInfo = self.observationInfo;
            for i = 1:length(uLsdKvs)
                
                dt = [obsInfo.dt];
                alarmsInNeed = find(dt(:) == 0 & strcmpi(lsdKvs(:),uLsdKvs{i}));
                
                disp(i./length(uLsdKvs));
                %this is slow, so only do it if we actually miss anything
                if length(alarmsInNeed) > 0
                    dataFile = tuf.locate(uLsdKvs{i},'data');
                    pos = readMdrScanHeadersCache(dataFile);
                    tempD = readMdrMex(dataFile,1,1);
                    nChannels = size(tempD,2);
                    %                     pos = as.niitekPos2PosMat(pos,nChannels);
                
                    
                    for alarmIndex = 1:length(alarmsInNeed)
                        n = obsInfo(alarmsInNeed(alarmIndex)).yUTM;
                        e = obsInfo(alarmsInNeed(alarmIndex)).xUTM;
                        [dtxt,inside,distance] = ne2dtxt([e,n],pos,nChannels);
                        
                        obsInfo(alarmsInNeed(alarmIndex)).dt = dtxt(1);
                        obsInfo(alarmsInNeed(alarmIndex)).ch = dtxt(2);
                    end
                end
            end
            self.observationInfo = obsInfo;
        end
    
        
        function self = prtDataSetTuf(varargin)
            %self = prtDataSetTuf(tufresults)
            %self = prtDataSetTuf(tufResultsMatFile)
            % Generate a prtDataSetTuf using either tufresults (what ends
            % up in your workspace after running TUF), or a tufResults MAT
            % file.  
            %
            % I prefer saving tufresults structures and using those, since
            % that lets you get a data set with a bunch of algorithms all
            % at once.
            
            self = self@prtDataSetClass; %super-class
            
            if isa(varargin{1},'Scoring.SimpleScore')
                tufresults.scoreObj = varargin{1}; %try and be cool
            elseif isa(varargin{1},'char') 
                tufresults = load(varargin{1});
                self = prtDataSetTuf(tufresults.scoreobj);
                return;
            elseif isa(varargin{1},'struct')
                tufresults = varargin{1};
            else
                error('First input argument must be a file name, or a tuf results structure');
            end
            
            if nargin == 1 || isempty(varargin{2})
                algorithmNames = fieldnames(tufresults);
            else
                algorithmNames = varargin{2};
                if isa(algorithmNames,'char')
                    algorithmNames = {algorithmNames};
                end
            end
            algorithmNames = algorithmNames(~strcmpi(algorithmNames,'all_score_args'));
            
            ds = cell(length(algorithmNames),1);
            for nameIndex = 1:length(algorithmNames)
                algorithmName = algorithmNames{nameIndex};
                fprintf('Loading data set for %s\n',algorithmName);
                ds{nameIndex} = parseTufResultsStructure(self,tufresults.(algorithmName));
                ds{nameIndex}.name = algorithmName;
            end
            
            %given the ds's from all the tufresults scoring objects,
            %generate the prtDataSet; we should really check that all the
            %observationInfo and userData structures match here...
            self.userData = ds{1}.userData;
            self.observationInfo = ds{1}.observationInfo;
            self.targets = ds{1}.targets; 
            for i = 1:length(ds)
                self = self.catFeatures(ds{i});
            end
            %build in some utility stuff
            self = self.setCrossValidationFolds;
            self = self.generateExcursionArray;
        end
        
        function lanes = getUniqueLanes(self)
            lanes = unique(self.getLanes);
        end
        
        function lanes = getLanes(self)
            lanes = {self.observationInfo.lane};
            lanes = lanes(:);
        end
        
        function lsd_kvs = getUniqueLsdKvs(self)
            lsd_kvs = unique(self.getLsdKvs);
        end
        
        function lsd_kvs = getLsdKvs(self)
            lsd_kvs = {self.observationInfo.lsd_kv};
            lsd_kvs = lsd_kvs(:);
        end
        
        
        function self = retainLanes(self,laneSpec)
            %self = retainLanes(self,laneSpec)
            if ~iscell(laneSpec)
                laneSpec = {laneSpec};
            end
            
            lanes = self.getLanes;
            retain = false(self.nObservations,1);
            for i = 1:length(laneSpec);
                retain = retain | strcmpi(lanes,laneSpec);
            end
            self = self.retainObservations(retain);
            
            laneExcursions = {self.excursionArray.lane};
            laneExcursions = laneExcursions(:);
            retain = false(length(laneExcursions),1);
            for i = 1:length(laneSpec);
                retain = retain | strcmpi(laneExcursions,laneSpec);
            end
            self = self.retainExcursions(retain);
            
        end
        
        function self = removeLanes(self,laneSpec)
            %self = removeLanes(self,laneSpec)
            if ~iscell(laneSpec)
                laneSpec = {laneSpec};
            end
            lanes = self.getUniqueLanes;
            retainLaneSpec = setdiff(lanes,laneSpec);
            self = self.retainLanes(retainLaneSpec);
        end
        
        function self = retainLsdKv(self,lsdKvSpec)
            %self = retainLsdKv(self,laneSpec)
            if ~iscell(lsdKvSpec)
                lsdKvSpec = {lsdKvSpec};
            end
            
            lsd_kv = self.getLsdKvs;
            retain = false(self.nObservations,1);
            for i = 1:length(lsdKvSpec);
                retain = retain | strcmpi(lsd_kv,lsdKvSpec);
            end
            self = self.retainObservations(retain);
            
            lsdExcursions = self.getUniqueLsdKvs;
            retain = false(length(lsdExcursions),1);
            for i = 1:length(lsdKvSpec);
                retain = retain | strcmpi(lsdExcursions,lsdKvSpec);
            end
            self = self.retainExcursions(retain);
            
        end
        
        function self = removeLsdKv(self,lsdKvSpec)
            %self = removeLsdKv(self,laneSpec)
            if ~iscell(lsdKvSpec)
                lsdKvSpec = {lsdKvSpec};
            end
            lsd_kv = {self.observationInfo.lsd_kv};
            retainLsdKvSpec = setdiff(lsd_kv,lsdKvSpec);
            self = self.retainLanes(retainLsdKvSpec);
        end
        
        
        
        function self = generateExcursionArray(self)
            %self = generateExcursionArray(self)
            % Build the lane excursion database; and store it in
            % .excursionArray
            lsds = {self.observationInfo.lsd_kv};
            lanes = {self.observationInfo.lane};
            
            [uniqueLsds,uInds] = unique(lsds);
            uniqueLanes = lanes(uInds);
            
            for i = 1:length(uniqueLsds)
                theLsd = uniqueLsds{i};
                tlane = tuf.locate(theLsd,'tlane');
                dataFile = tuf.locate(theLsd,'data');
                
                sensorExtent = gprScoreMdrToExtent(dataFile);
                
                objectSet = as.asc2ObjectSet(tlane);
                tlaneExtent = extentPolygon('polygon',objectSet.Extent.parameters);
                excursionExtent = intersection(tlaneExtent,sensorExtent);
                
                self.excursionArray(i).lsd_kv = theLsd;
                self.excursionArray(i).lane = uniqueLanes{i};
                self.excursionArray(i).extent = excursionExtent;
                self.excursionArray(i).area = excursionExtent.getArea;
            end
        end
        
        function [far,pd,threshold,algoNames] = scoreRoc(self)
            %[far,pd,threshold,algoNames] = scoreRoc(self)
            
            if self.nFeatures > 10
                error('That seems like a lot of ROCs...');  %to do, make this check with you, the user
            end
            
            [nfa,pd,threshold] = prtScoreRocNfa(self);
            far = nfa;
            if iscell(nfa)
                for i = 1:length(nfa)
                    far{i} = nfa{i}./self.userData.scorableArea;
                end
            else
                far = nfa./self.userData.scorableArea;
            end
            algoNames = self.featureNames;
            
            if nargout == 0
                if iscell(nfa)
                    for i = 1:length(nfa)
                        h = plot(cat(2,far{:}),cat(2,pd{:}));
                    end
                else
                     h = plot(far,pd);
                end
                ylim([0 1]);
                xlim([0 .02]);
                legend(h,algoNames,'interpreter','none','location','southEast');
            end
        end
        
        
        function self = setCrossValidationFolds(self,xValType)
            %self = setCrossValidationFolds(self,xValType)
            
            if nargin < 2 || isempty(xValType)
                xValType = 'stratified_alarm';
            end
            
            switch lower(xValType)
                case 'random'
                    folder = tuf.crossval.Random;
                case 'region_based'
                    folder = tuf.crossval.Region_Based;
                case 'site_based'
                    folder = tuf.crossval.Site_Based;
                case 'stratified_alarm'
                    folder = tuf.crossval.Stratified_Alarm;
                otherwise
                    error('xValType must be one of {''random'', ''region_based'', ''site_based'', ''stratified_alarm''}');
            end
            
            if ~(isfield(self.userData,'nAlarmsPerSsid') && isfield(self.userData,'ssids'))
                % Must add these fields
                
                sids = arrayfun(@(s)tuf.get_sample_entry(s.lsd_kv).sid,self.observationInfo,'uniform',false);
                [uSids, ~, uSidsY] = unique(sids);
                
                % Make sure things are sorted properly
                [~, sortingInds] = sort(uSidsY);
                self = self.retainObservations(sortingInds); % Resort the dataset
                
                nUSids = hist(uSidsY,1:length(uSids));
                
                self.userData.nAlarmsPerSsid = nUSids(:);
                self.userData.ssids = uSids;
            end
            
            alarmSet = mat2cell(self.observationInfo,self.userData.nAlarmsPerSsid,1); % This line assumes your alarms are ordered according to the SIDS be careful!!!
            folder.make_folds(self.userData.ssids, alarmSet);
            
            nFolds = folder.num_folds;
            
            nAlarmsAbove = cat(1,0, cumsum(self.userData.nAlarmsPerSsid(1:end-1)));
            
            crossvalInd = zeros(self.nObservations,1);
            for iFold = 1:nFolds
                [~, fold_map] = folder.test_fold(iFold,alarmSet);
                
                
                for iSsid = 1:length(fold_map)
                    cInds = fold_map{iSsid} + nAlarmsAbove(iSsid);
                    crossvalInd(cInds) = iFold;
                end
            end
            
            % Check and make sure we did this right
            % testAlarms = folder.test_fold(2,alarmSet);
            % ourTestAlarms = self.observationInfo(crossvalInd==2)';
            %
            % isequal(testAlarms,ourTestAlarms)
            %  %ans =
            %  %   1
            self.userData.crossValidationType = xValType;
            self = self.setObservationInfo('crossValidationIndex', crossvalInd);
        end
        
        function ds = parseTufResultsStructure(self,tufresults) %#ok<MANU>
            %ds = parseTufResultsStructure(self,tufresults)
            
            % Interpret SimpleScore Obj
            % tufresults is now a SimpleScore object, tufresults (as written to the
            % workspace by TUF is a structure containing fields with the names of the
            % prescreeners/algorithms that were run. Above we rip out a particulare
            % algorithm and call tufresults just the SimpleScore object.
            
            
            % Re run the object if necessary, forces recalculation of the areas
            % But this recalculation uses tufs cached stuff or stuff saved in the
            % object, I don't know but its' pretty fast only slow because it does all
            % of the linking. Which is a shame because we do it again below.
            tufresults.run(false);
            
            % SimpleScore saves the targets and clutter after linking, but it does it
            % wacky. Very wacky. So we redo it. This is the slow part. I wish
            % simplescore saved this information in a usable way. It saves the targets
            % and clutter but not the false alarms. Since we like to look at our top
            % false alarms we need to do it again.
            
            removeOutsideLane = true;
            runStrs = tufresults.alarms.keys';
            alarmsCell = tufresults.alarms.values;
            linkerCell = tufresults.linkers_.values;
            targets = [];
            targetInds = [];
            clutter = [];
            clutterInds = [];
            fas = [];
            
            %generate the structures we need; these are targets, clutter,
            %and fas
            for iLinker = 1:length(linkerCell)
                cAlarms = alarmsCell{iLinker};
                cLinker = linkerCell{iLinker};
                %for some reason, alarms from algorithms dont have lsd_kv's
                %internally.  Why?  We don't know.  But we can guess at the
                %right lsd_kv
                [cTargetObjects, cClutterObjects, cFalseAlarms] = cLinker.link_objects(cAlarms, tufresults.ignore_clutter);
                %                 X = cat(1,cat(1,cTargetObjects.x),cat(1,cClutterObjects.x),cat(1,cFalseAlarms.xUTM));
                %                 Y = cat(1,cat(1,cTargetObjects.y),cat(1,cClutterObjects.y),cat(1,cFalseAlarms.yUTM));
                %                 plot(X(:),Y(:),'b.');
                %                 pause;
                if (~isempty(cTargetObjects) && ~isfield(cTargetObjects,'lsd_kv')) || (~isempty(cClutterObjects) && ~isfield(cClutterObjects,'lsd_kv')) || (~isempty(cFalseAlarms) && ~isfield(cFalseAlarms,'lsd_kv')) 
                    for i = 1:length(cTargetObjects)
                        cTargetObjects(i).lsd_kv = runStrs{iLinker};
                        cTargetObjects(i).sample = runStrs{iLinker};
                    end
                    for i = 1:length(cClutterObjects)
                        cClutterObjects(i).lsd_kv = runStrs{iLinker};
                        cClutterObjects(i).sample = runStrs{iLinker};
                    end
                    for i = 1:length(cFalseAlarms)
                        cFalseAlarms(i).lsd_kv = runStrs{iLinker};
                        cFalseAlarms(i).sample = runStrs{iLinker};
                    end
                end
                targets = cat(1,targets,cTargetObjects(:));
                targetInds = cat(1,targetInds,ones(length(cTargetObjects),1)*iLinker);
                
                clutter = cat(1,clutter,cClutterObjects(:));
                clutterInds = cat(1,clutterInds,ones(length(cClutterObjects),1)*iLinker);
                
                fas = cat(1,fas,cFalseAlarms(:));
            end
            
            %generate alarm structures for all the targets.
            for iTarget = 1:length(targets)
                cTarget = targets(iTarget);
                cAlarms = cTarget.alarms;
                
                if ~isempty(cAlarms)
                    cConfs = cat(1,cAlarms.conf);
                    [~, bestAlarm] = max(cConfs);
                    cAlarm = cAlarms(bestAlarm);
                    cAlarm.isRealAlarm = true;
                else
                    cAlarm = prtDataSetTuf.targetToFakeAlarm(cTarget);
                end
                cAlarm.targetId = cTarget.info.id;
                cAlarm.lsd_kv = runStrs{targetInds(iTarget)};
                cAlarm.tufTarget = cTarget;
                cAlarm.isClutter = false; % oofta, it's a target
                %                 cAlarm.isTarget = true;
                
                cAlarm = prtDataSetTuf.standardizeAlarm(cAlarm);
                
                if iTarget == 1
                    targetAlarms = repmat(cAlarm, length(targets),1);
                else
                    targetAlarms(iTarget) = cAlarm;
                end
            end
            if removeOutsideLane
                targetAlarms = targetAlarms(cat(1,targets.on_lane)); % Remove targets not in the lane.
            end
            
            %generate alarm structures for all the clutter
            for iClutter = 1:length(clutter)
                cClutter = clutter(iClutter);
                cAlarms = cClutter.alarms;
                
                if ~isempty(cAlarms)
                    cConfs = cat(1,cAlarms.conf);
                    [~, bestAlarm] = max(cConfs);
                    cAlarm = cAlarms(bestAlarm);
                    cAlarm.isRealAlarm = true;
                else
                    cAlarm = prtDataSetTuf.clutterToFakeAlarm(cClutter);
                end
                
                cAlarm.targetId = cClutter.info.id;
                cAlarm.id = cClutter.info.id;
                %Fixed bug, 2012.09.15
                cAlarm.lsd_kv = runStrs{clutterInds(iClutter)};
                cAlarm.isClutter = true; % you betcha
                %                 cAlarm.isTarget = false;
                
                cAlarm = prtDataSetTuf.standardizeAlarm(cAlarm);
                
                if iClutter == 1
                    clutterAlarms = repmat(cAlarm, length(clutter),1);
                else
                    clutterAlarms(iClutter) = cAlarm;
                end
            end
            if removeOutsideLane            
                clutterAlarms = clutterAlarms(cat(1,clutter.on_lane)); % Remove clutter not in the lane.
            end
            
            for iFalseAlarm = 1:length(fas)
                tempStruct = prtDataSetTuf.standardizeAlarm(fas(iFalseAlarm));
                tempStruct.targetId = 'False Alarm';
                tempStruct.isRealAlarm = true;
                tempStruct.isClutter = false; % nope, not as far as we know
                %                 tempStruct.isTarget = false;
                falseAlarms(iFalseAlarm) = tempStruct;
            end
            
            alarms = cat(1,targetAlarms(:),clutterAlarms(:),falseAlarms(:)); %NOTE - FALSE ALARMS MAY CONTAIN CLUTTER
            
            %we do this here instead of up there... because... i don't really know why
            for i = 1:length(alarms)
                dd = tuf.get_sample_entry(alarms(i).lsd_kv);
                alarms(i).lane = dd.region.description;
            end
            
            %ds = prtDataSetClass(cat(1,alarms.conf),cat(1,alarms.isTarget));
            ds = prtDataSetClass(cat(1,alarms.conf),cat(1,alarms.hit));
            ds = ds.setFeatureNames({tufresults.display_name});
            ds = ds.setClassNames({'Non-Mine','Mine'});
            ds = ds.setObservationInfo(alarms(:));
            
            ds.userData.laneArea = tufresults.area_;
            ds.userData.targetArea = tufresults.scored_target_area;
            ds.userData.scorableArea = tufresults.scored_area;
            ds = ds.setObservationInfo('isClutterFalseAlarm',arrayfun(@(s)~isempty(s.objects)&~s.hit,ds.observationInfo));
        end 
        
        function explore(self)
            AdditionalOptions.additionalOnClickFunction = @(ds,clickedObsInd,currentlyDisplayedFeatures)self.exploreClickFunction(ds,clickedObsInd,currentlyDisplayedFeatures);
            explore@prtDataSetClass(self, AdditionalOptions)
        end
        
        function excursionCell = toExcursions(self)
            
            lsds = {self.observationInfo.lsd_kv};
            uLsds = unique(lsds);
            for lsdInd = 1:length(uLsds);
                cData = self.retainObservations(strcmpi(lsds,uLsds{lsdInd}));
                cAlarms = cData.retainObservations([cData.observationInfo.isDetected]==1);
                xy = cat(2,[cAlarms.observationInfo.xUTM]',[cAlarms.observationInfo.yUTM]');
                conf = cAlarms.X;
                
                mdrFile = tuf.locate(cAlarms.observationInfo(1).sample,'data');
                tlane = tuf.locate(cAlarms.observationInfo(1).sample,'tlane');
                
                extent = gprScoreNiitekPosToExtent(mdrFile);
                alarmSet = gprAlarmSet;
                alarmSet = alarmSet.fromGpsLocAndConf(xy,conf,extent);
                objectSet = gprObjectSet;
                objectSet = objectSet.ascRead(tlane);
                
                excursion = gprExcursion('alarmSet',alarmSet,'objectSet',objectSet);
                excursion = excursion.link;
                excursionCell{lsdInd} = excursion;
            end
                
        end
        
        function exploreClickFunction(self, ds,clickedObsInd,currentlyDisplayedFeatures)
            obsStruct = ds.observationInfo(clickedObsInd);
            
            dtxt = [obsStruct.dt obsStruct.ch];
            
            xy = [obsStruct.xUTM obsStruct.yUTM];
            mdrFile = tuf.locate(ds.observationInfo(clickedObsInd).sample,'data');
            
            if any(dtxt==0)
                % you need to look up the dtxt from the pos and utm
                scanHeaders = readMdrScanHeadersCache(mdrFile);
                dTemp = readMdrMex(mdrFile,1,1);
                dtxt = gpsToDtXt(xy,niitekPos2PosMat(scanHeaders,size(dTemp,2)));
            end
            
            alarmSet = gprAlarmSet('dataFile', mdrFile);
            figure;
            displayDtXt(alarmSet, dtxt);
            %             set(gcf,'color',[1 1 1]);
            if ~isempty(obsStruct.tufTarget)
                titleString = sprintf('%s @ %din - %s',obsStruct.tufTarget.name,obsStruct.tufTarget.depth,mat2str(obsStruct.conf,4));
            else
                titleString = sprintf('%s',mat2str(obsStruct.conf,4));
            end
            
            title(titleString,'FontSize',14,'Interpreter','none'); % expand vertically by 0.02 normalized units
                            
            % Position string
            posString = sprintf(' DT: %d, CT: %d; N: %f, E: %f;', dtxt(1), dtxt(2), xy(2), xy(1));
            
            figurePixelPos = getpixelposition(gcf);
            
            hText = uicontrol(gcf,'style','text',...
                'string',cat(2, posString, ' ', mdrFile),...
                'units','pixels',...
                'position',[0 0 figurePixelPos(3) 15],...
                'HorizontalAlignment','left','backgroundcolor',[1 1 1]);
            
        end
    end
    
    methods (Static)
        
        function alarmOut = standardizeAlarm(cAlarm)
            %clean up alarm information so they all act presentable.
            translateFields = {'HIT','hit'; 'TARGET','target'; 'MineIsInsideLane','mineIsInsideLane'; 'DIST','objectDistance'; 'IsInsideLane','isInsideLane'};
            requiredFieldDefaults = {'xUTM',nan; 'yUTM',nan; 'conf',nan; 'dt',0;'ch',0; 'objects',struct; 'hit',nan; 'target',struct; ...
                'targetId',''; 'category',nan; 'isDetected',nan;  'isInsideLane',nan; 'mineIsInsideLane',nan; 'objectDistance',nan; 'objectdist',nan; ...
                'lsd_kv',''; 'sensor_type','unknown'; 'prescreener','unknown'; 'tufTarget',[]; 'isClutter',false; 'isRealAlarm', false; 'sample',''; ...
                'event_id',''; 'event_timestamp',''; 'event_type',''};
            
            %Translate from different versions to uniform camelCase
            for i = 1:size(translateFields,1)
                if isfield(cAlarm,translateFields{i,1})
                    cAlarm.(translateFields{i,2}) = cAlarm.(translateFields{i,1});
                    cAlarm = rmfield(cAlarm,translateFields{i,1});
                end
            end
            cAlarm.sample = cAlarm.lsd_kv;
            
            %make the output with the defaults; this makes sure the alarm
            %fields are always in the same order:
            tempRequiredFieldDefaults = requiredFieldDefaults';
            tempRequiredFieldDefaults = tempRequiredFieldDefaults(:)';
            alarmOut = struct(tempRequiredFieldDefaults{:});
            
            %For all the fields already set in cAlarm, move them into
            %alarmOut.  This implicityly removes all fields not in
            %requiredFieldDefaults(:,1)
            for i = 1:size(requiredFieldDefaults,1)
                if isfield(cAlarm,requiredFieldDefaults{i,1})
                    alarmOut.(requiredFieldDefaults{i,1}) = cAlarm.(requiredFieldDefaults{i,1});
                end
            end
        end
        
        function fakeAlarm = targetToFakeAlarm(target)
            fakeAlarm.xUTM = target.geom.getCentroid.getX;
            fakeAlarm.yUTM = target.geom.getCentroid.getY;
            fakeAlarm.conf = nan;
            fakeAlarm.dt = 0;
            fakeAlarm.ch = 0;
            fakeAlarm.objects = target;
            fakeAlarm.HIT = 1;
            fakeAlarm.TARGET.ID = target.info.id;
            fakeAlarm.TARGET.DEPTH = target.depth;
            fakeAlarm.TARGET.CORNER = []; % Don't know what this is?
            fakeAlarm.TARGET.N = target.y;
            fakeAlarm.TARGET.E = target.x;
            %fakeAlarm.TARGET.LANE = target.lane;
            fakeAlarm.TARGET.LANE = target.extra;
            fakeAlarm.TARGET.RAD = [];
            fakeAlarm.TARGET.NUMID = target.id;
            
            fakeAlarm.MineIsInsideLane = 1;
            fakeAlarm.DIST = 0;
            fakeAlarm.objectdist = 0;
            fakeAlarm.category = nan;
            fakeAlarm.isDetected = 0;
            fakeAlarm.IsInsideLane = 1;
            fakeAlarm.isRealAlarm = false;
        end
        
        function fakeAlarm = clutterToFakeAlarm(clutter)
            fakeAlarm.xUTM = clutter.geom.getCentroid.getX;
            fakeAlarm.yUTM = clutter.geom.getCentroid.getY;
            fakeAlarm.conf = nan;
            fakeAlarm.dt = 0;
            fakeAlarm.ch = 0;
            fakeAlarm.objects = clutter;
            fakeAlarm.HIT = 0;
            fakeAlarm.TARGET.ID = clutter.info.id;
            fakeAlarm.TARGET.DEPTH = clutter.depth;
            fakeAlarm.TARGET.CORNER = []; % Don't know what this is?
            fakeAlarm.TARGET.N = clutter.y;
            fakeAlarm.TARGET.E = clutter.x;
            %fakeAlarm.TARGET.LANE = clutter.lane;
            fakeAlarm.TARGET.LANE = clutter.extra;
            fakeAlarm.TARGET.RAD = [];
            fakeAlarm.TARGET.NUMID = clutter.id;
            
            fakeAlarm.MineIsInsideLane = 1;
            fakeAlarm.DIST = 0;
            fakeAlarm.objectdist = 0;
            fakeAlarm.category = nan;
            fakeAlarm.isDetected = 0;
            fakeAlarm.IsInsideLane = 1;
            fakeAlarm.isRealAlarm = false;
        end
        
        function pat = strCellToRegexpPattern(cellStr,ez)
            %turn a cell array into a regexp pattern of or's
            pat = '';
            if ~isa(cellStr,'cell')
                cellStr = {cellStr};
            end
            for i = 1:length(cellStr)
                currentStr = cellStr{i};
                if ez
                    currentStr = sprintf('^%s$',currentStr); %regexp is over eager without the leadning ^ and trailing $
                    currentStr = strrep(currentStr,'*','.*'); %regexp is over eager without the trailing $
                    currentStr = lower(currentStr); %It's EZ!
                end
                if i == 1
                    pat = sprintf('(%s)',currentStr);
                else
                    pat = sprintf('%s|(%s)',pat,currentStr);
                end
            end
        end
        
        function [n,p] = summaryToMatrices(summary)
            n = cell(1);
            p = cell(1);
            for i = 1:length(summary)
                for j = 1:length(summary(i).localSummary)
                    for k = 1:length(summary(i).localSummary(j).depths);
                        n{i}(j,summary(i).localSummary(j).depths(k) == summary(1).globalDepths) = summary(i).localSummary(j).totalPerDepth(k);
                        p{i}(j,summary(i).localSummary(j).depths(k) == summary(1).globalDepths) = summary(i).localSummary(j).foundPerDepth(k);
                    end
                end
                n{i} = round(n{i});
                p{i} = round(p{i});
            end
        end
    end
end
