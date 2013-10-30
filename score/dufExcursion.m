classdef dufExcursion
    
    properties
        alarmSet
        objectSet
        linkedAlarmsObjects
        intersectedExtent
        
        scoring
        
        ogAlarmSet
        ogObjectSet
        objectLabels
        
        extractedConfidences
        extractedLabels
        extractedTypes
        extractedDepths
        extractedLocations
        
        extentPlotOptions = struct('edgeColor','g','facecolor','none','linewidth',1,'linestyle','--');
        
        removeObjectFunction = @(s)s.isClutter;
        setObjectLabelFunction = @(s)s.isThreat;
    end
    
    methods
        
        function self = reset(self)
            if ~isempty(self.ogAlarmSet)
                self = dufExcursion('alarmSet',self.ogAlarmSet,'objectSet',self.ogObjectSet);
            else
                return; %nothing to reset to.
            end
        end
        
        function self = dufExcursion(varargin)
            self = prtUtilAssignStringValuePairs(self,varargin{:});
        end
        
        function [hObjects,hAlarms] = plot(self)
            
            hObjects = self.objectSet.plot;
            set(hObjects.extent,'hittest','off');
            set(cat(1,hObjects.objectHandle{:}),'hittest','off');
            set(cat(1,hObjects.haloHandle{:}),'hittest','off');
            
            hold on;
            hAlarms = self.alarmSet.plot;
            set(hAlarms.extent,'hittest','off');
            set(hAlarms.alarmHandle,'hittest','off');
            
            hold on;
            if ~isempty(self.intersectedExtent)
                hExtent = self.intersectedExtent.plot;
                dufScoreUtilApplyPlotOptions(hExtent,self.extentPlotOptions);
                set(hExtent,'hitTest','off');
            end
            
            hold off;
        end
        
        function self = link(self)
            %self = link(self)
            points = [];
            for i = 1:length(self.alarmSet.alarmArray)
                points = cat(1,points,self.alarmSet.alarmArray(i).extent.toPolygon);
            end
            self.linkedAlarmsObjects = sparse(length(self.alarmSet.alarmArray),length(self.objectSet.objectArray));
            for i = 1:length(self.objectSet.objectArray)
                within = self.objectSet.objectArray(i).extent.withinHalo(points,self.objectSet.halo);
                
                %Put objects into alarms:
                if any(within)
                    withinInds = find(within);
                    for ind = withinInds(:)';
                        self.linkedAlarmsObjects(ind,i) = 1;
                        self.alarmSet.alarmArray(ind).onObject = true;
                        self.alarmSet.alarmArray(ind).object = self.objectSet.objectArray(i);
                    end
                end
            end
        end
        
        function [self,alarmMarkedForDeletion] = intersect(self)
            
            self.intersectedExtent = intersection(self.objectSet.extent,self.alarmSet.extent);
           
            objectMarkedForDeletion = false(size(self.objectSet.objectArray));
            alarmMarkedForDeletion = false(size(self.alarmSet.alarmArray));
            
            self.ogAlarmSet = self.alarmSet;
            self.ogObjectSet = self.objectSet;
            %If object center is outside intersected region; throw it out!
            for i = 1:length(self.objectSet.objectArray)
                if ~withinHalo(self.intersectedExtent,self.objectSet.objectArray(i).extent.getCenter,0)
                    objectMarkedForDeletion(i) = true;
                    alarmMarkedForDeletion(logical(self.linkedAlarmsObjects(:,i))) = true;
                end
                %To do: calculate area of removed objects
            end
            
            for i = 1:length(self.alarmSet.alarmArray)
                if ~withinHalo(self.intersectedExtent,self.alarmSet.alarmArray(i).extent.getCenter,0)
                    alarmMarkedForDeletion(i) = true;
                end
            end
            self.alarmSet.alarmArray = self.alarmSet.alarmArray(~alarmMarkedForDeletion);
            self.objectSet.objectArray = self.objectSet.objectArray(~objectMarkedForDeletion);
            self.linkedAlarmsObjects = self.linkedAlarmsObjects(~alarmMarkedForDeletion,~objectMarkedForDeletion);
        end
        
        %formerly asSummarize
        function self = removeObjects(self,fnHandle)
            %self = removeObjects(self,fnHandle)
            %
            
            objectMarkedForDeletion = false(size(self.objectSet.objectArray));
            alarmMarkedForDeletion = false(size(self.alarmSet.alarmArray));
            for i = 1:length(self.objectSet.objectArray)
                objectMarkedForDeletion(i) = fnHandle(self.objectSet.objectArray(i));
            end
            for iAlarm = 1:length(self.alarmSet.alarmArray)
                cObjects = find(logical(self.linkedAlarmsObjects(iAlarm,:)));
                
                if ~isempty(cObjects)
                    % Linked to at least one object
                    
                    if all(objectMarkedForDeletion(cObjects))
                        % All objects that it is linked to are going to be
                        % deleted -> Delete alarm
                        alarmMarkedForDeletion(iAlarm) = true;
                    end
                    
                end
            end
            
            self.alarmSet.alarmArray = self.alarmSet.alarmArray(~alarmMarkedForDeletion);
            self.objectSet.objectArray = self.objectSet.objectArray(~objectMarkedForDeletion);
            self.linkedAlarmsObjects = self.linkedAlarmsObjects(~alarmMarkedForDeletion,~objectMarkedForDeletion);
        end
        
        function self = setObjectLabels(self,fnHandle)
            self.objectLabels = nan(size(self.objectSet.objectArray));
            for i = 1:length(self.objectSet.objectArray)
                self.objectLabels(i) = fnHandle(self.objectSet.objectArray(i));
            end
        end
        
        function self = aggregateAlarms(self)
            
            conf = nan(size(self.alarmSet.alarmArray));
            for i = 1:length(self.alarmSet.alarmArray)
                conf(i) = self.alarmSet.alarmArray(i).confidence;
            end
            
            ds = nan(length(self.objectSet.objectArray),1);
            y = ds;
            for i = 1:length(self.objectSet.objectArray)
                
                if any(self.linkedAlarmsObjects(:,i))
                    [ds(i),maxInd] = max(conf(logical(self.linkedAlarmsObjects(:,i))));
                    
                    localAlarms = self.alarmSet.alarmArray(logical(self.linkedAlarmsObjects(:,i)));
                    theAlarm = localAlarms(maxInd);
                    self.extractedLocations(i,:) = theAlarm.extent.getCenter;
                else
                    self.extractedLocations(i,:) = self.objectSet.objectArray(i).extent.getCenter;
                end
                
                y(i) = self.objectLabels(i);
                self.extractedTypes{i,1} = self.objectSet.objectArray(i).name;
                self.extractedDepths(i,1) = self.objectSet.objectArray(i).depth;
            end
            
            falseAlarmLogical = ~full(any(self.linkedAlarmsObjects == 1,2));
            falseAlarms = self.alarmSet.alarmArray(falseAlarmLogical);
            for i = 1:length(falseAlarms)
                self.extractedLocations(end+1,:) = falseAlarms(i).extent.getCenter;
            end
            confFa = conf(falseAlarmLogical);
            confFa = confFa(:);
            yFa = zeros(size(confFa));
            
            if isempty(self.extractedTypes)
                self.extractedTypes = {};
            end
            self.extractedTypes(end+1:end+length(confFa),1) = {'False Alarm'};
            self.extractedDepths(end+1:end+length(confFa),1) = nan;
            
            self.extractedConfidences = cat(1,ds,confFa);
            self.extractedLabels = cat(1,y,yFa);
            
        end
        
        function self = score(self)
            
            self = self.link;
            self = self.intersect;
            self = self.removeObjects(self.removeObjectFunction);
            self = self.setObjectLabels(self.setObjectLabelFunction);
            self = self.aggregateAlarms;
            
            self.scoring.area = self.intersectedExtent.getArea;
            
            try
                if length(unique(self.extractedLabels)) == 1
                    self.extractedLabels = cat(1,self.extractedLabels(:),0);
                    self.extractedConfidences = cat(1,self.extractedConfidences(:),nan);
                end
                self.scoring.confidence = self.extractedConfidences;
                self.scoring.labels = self.extractedLabels;
                [self.scoring.pf,self.scoring.pd,self.scoring.rocThresholds,self.scoring.auc] = prtScoreRoc(self.extractedConfidences,self.extractedLabels);
                [self.scoring.nfa,self.scoring.pd,self.scoring.rocThresholds,self.scoring.auc] = prtScoreRocNfa(self.extractedConfidences,self.extractedLabels);
            catch ME
                disp('Could not make ROC curve')
                disp(ME);
            end
        end
        
        
        function gui(self)
            h = plot(self);
            set(gca,'buttondownfcn',@(e,v)self.mouseButtonDown(e,v));
        end
        
        function mouseButtonDown(self,e,v)
            disp('unfinished; use gprExcursion');
        end
        
    end
    
end
