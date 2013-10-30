classdef dufEventSet < prtDataSetClass
    
    properties
        excursionSummaries
        excursions
        scoring
        eventSummary
    end
    
    methods
        function self = dufEventSet(varargin)
            
            warning('We have things we need to add here in dufEventSet.dufEventSet().');
            % Add a check for uid and make one if necessary (use integers
            % as strings).
            % Anything else we need to assert here?
            
            
            if nargin == 1
                self.X = cat(1,varargin{1}.confidence);
                self.observationInfo = varargin{1};
            else
                self = prtUtilAssignStringValuePairs(self, varargin{:});
            end
        end
        
        function self = setTargetsFromObservationInfo(self,obsProcessFn)
            % ds = ds.setTargetsFromObservationInfo(@(o)o.isOnObject);
            % 
            oInfo = self.observationInfo;
            class = zeros(length(oInfo),1);
            for i = 1:length(oInfo)
                class(i,1) = obsProcessFn(oInfo(i));
            end
            self.targets = class;
            
            
        end
        
        function self = aggregate(self,aggregationType)
            
            if nargin < 2
                aggregationType = 'max';
            end
            % Handle stuff with truth linking here?
            objectInfo = cat(1,self.observationInfo);
            objectIds = nan(length(objectInfo),1);
            isAggregated = false(length(objectInfo),1);
            
            for eventInd = 1:length(objectInfo)
                if ~isempty(objectInfo(eventInd).objects)
                    objectIds(eventInd) = objectInfo(eventInd).objects(1).info.numID;
                end
            end
            uObjects = unique(objectIds);
            uObjects = uObjects(~isnan(uObjects));
            
            allCollectionInfo = cat(1,self.observationInfo.collectionInfo);
            uidList = {allCollectionInfo.uid}';
            [uids, lastOccurance, collectionIndex] = unique(uidList);
            
            for uidInd = 1:length(uids)
                inUid = strcmpi(uids{uidInd},uidList);
                
                for objectInd = 1:length(uObjects)
                    currentObjectInd = find(inUid & objectIds == uObjects(objectInd));
                    if length(currentObjectInd) > 1
                        switch aggregationType
                            case 'max'
                                [v,maxInd] = max(self.X(currentObjectInd,1)); % Note we only look at this thing
                                isAggregated(currentObjectInd(setdiff(1:length(currentObjectInd),maxInd))) = true;
                            otherwise
                                error('invalid type');
                        end
                    end
                end
            end
            self = self.setObservationInfo('isAggregated',isAggregated);
        end
        
        function self = score(self)
    
            
            if self.nFeatures > 10
                yourAnswer = input('That seems like a lot of ROCs. Do you really want to do this? yes or no?','s');
                if isempty(yourAnswer) || ~strcmpi(yourAnswer(1),'y');
                    return
                end
            end
            
            if isempty(self.excursions)
                self = self.link();
            end
            
            excursionIds = {self.excursionSummaries.collectionInfo.uid};
            
            allAreas = zeros(length(excursionIds),1);
            for iExcursion = 1:length(excursionIds)
                allAreas(iExcursion) = self.excursions{iExcursion}.intersectedExtent.getArea;
            end
            warning('There is a weird issue here in dufEventSet.score(). What if a lane has no alarms? I don''t know what to do.');
            % What if an excursion has no objects? Must take care of that
            % in run.
            
            temp = self.retainObservations(~[self.observationInfo.isAggregated]);
            [nfa,pd] = prtScoreRocNfa(temp);
            plot(nfa,pd);
            self.scoring = struct('dataSetAggregated',temp,'nfa',nfa,'pd',pd,'scorableArea',sum(allAreas));
        end
    end
end