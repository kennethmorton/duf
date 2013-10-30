classdef dufObjectSet
    
    properties
        objectArray
        extent
        source
        halo = 0.25;
        
        plotOptions = struct('edgecolor','none','facecolor',[0.9 0.9 0.9]);
    end
    
    methods
        
        function self = dufObjectSet(varargin)
            self = prtUtilAssignStringValuePairs(self,varargin{:});
        end
        
        function handleArray = plot(self)
            handleArray.extent = plot(self.extent);
            dufScoreUtilApplyPlotOptions(handleArray.extent,self.plotOptions);
            for i = 1:length(self.objectArray)
                hold on;
                if i == 1
                    handleArray.objectHandle{1} = plot(self.objectArray(i));
                    handleArray.haloHandle{1} = plotHalo(self.objectArray(i),self.halo);
                else
                    handleArray.objectHandle{i} = plot(self.objectArray(i));
                    handleArray.haloHandle{i} = plotHalo(self.objectArray(i),self.halo);
                end
            end
            hold off;
        end
        
        function toTruthEle(self,eleFile)
            
            for objectInd = 1:length(self.objectArray)
                curObj = self.objectArray(objectInd);
                center = getCenter(curObj.extent);
                
                alarmArray(objectInd) = dufAlarm;
                alarmArray(objectInd).confidence = 1;
                alarmArray(objectInd).extent = extentInstantaneous('point',center);
            end
            
            alarmSet = dufAlarmSet;
            alarmSet.alarmArray = alarmArray; 
            alarmSet.eleWrite(eleFile);
        end
        
        function self = ascRead(self,ascFile)
            %self = ascRead(self,ascFile)
            
            %Check file exists
            self.source = ascFile;
            objectSet = as.asc2ObjectSet(ascFile);
            
            self.extent = extentPolygon('polygon',objectSet.Extent.parameters);
            for i = 1:length(objectSet.Objects)
                if i == 1
                    self.objectArray = dufObject.fromObjectStruct(objectSet.Objects(i));
                else
                    self.objectArray(i) = dufObject.fromObjectStruct(objectSet.Objects(i));
                end
            end
            
        end
        
    end
end