classdef dufAlarmSet
    
    properties
        alarmArray
        extent
        
        alarmFile
        dataFile
        
        extentPlotOptions = struct('edgeColor','k','faceColor','none','linewidth',2);
        
        alarmColorBins = linspace(2,5,256);
        alarmColorMap = flipud(autumn(256));
    end
    
    methods
        
        
        function self = dufAlarmSet(varargin)
            self = prtUtilAssignStringValuePairs(self,varargin{:});
        end
        
        function handleArray = plot(self)
            handleArray.extent = plot(self.extent);
            dufScoreUtilApplyPlotOptions(handleArray.extent,self.extentPlotOptions);
            
            alarmColors = self.getAlarmColors();
            
            hold on;
            for i = 1:length(self.alarmArray)
                handleArray.alarmHandle(i) = plot(self.alarmArray(i));
                set(handleArray.alarmHandle(i),'color',alarmColors(i,:));
            end
            hold off;
        end
        
        function self = fromGpsLocAndConf(self,gpsLoc,conf,extentIn)
            % alarmSet = fromGpsLocAndConf(alarmSet,gpsLoc,conf)
            % alarmSet = fromGpsLocAndConf(alarmSet,gpsLoc,conf,extent)
            for i = 1:size(gpsLoc,1);
                alarm = struct('location',gpsLoc(i,:),'confidence',conf(i));
                if i == 1
                    self.alarmArray = dufAlarm.fromAlarmStruct(alarm);
                else
                    self.alarmArray(i) = dufAlarm.fromAlarmStruct(alarm);
                end
            end
            if nargin < 4
                %Tight:
                hullInds = convhull(gpsLoc(:,1),gpsLoc(:,2));
                self.extent = extentPolygon('polygon',gpsLoc(hullInds,:));
            else
                self.extentIn = extentIn;
            end
        end
        
        function self = eleRead(self,eleFile,extent)
            %self = eleRead(self,eleFile,extent)
            
            self.alarmFile = eleFile;
            alarmSet = as.ele2AlarmSet(eleFile);
            
            if nargin < 3
                self.extent = extentEverywhere; %can't get extent from ELE
            else
                self.extent = extent;
            end
            
            for i = 1:length(alarmSet.Alarms)
                if i == 1
                    self.alarmArray = dufAlarm.fromAlarmStruct(alarmSet.Alarms(i));
                else
                    self.alarmArray(i) = dufAlarm.fromAlarmStruct(alarmSet.Alarms(i));
                end
            end
        end
        
        function eleWrite(self,eleFile)
            assertDir(fileparts(eleFile));
            fid = fopen(eleFile,'w');
            for i = 1:length(self.alarmArray)
                alarm = self.alarmArray(i);
                p = getCenter(alarm.extent);
                
                fprintf(fid,'%d %d %d 1 1 N 1 1 E 0.000000 %f %f %c %f\r\n',i,nan,nan,...
                    p(2),p(1),alarm.sensor,alarm.confidence);
            end
            fclose(fid);
        end
        
        function colors = getAlarmColors(self)
            confs = cat(1, self.alarmArray.confidence);
            [~, colorMapInds] = histc(confs, cat(2,self.alarmColorBins,inf));
            colorMapInds(colorMapInds < 1) = 1;
            colorMapInds(colorMapInds > size(self.alarmColorMap,1)) = size(self.alarmColorMap,1);
            
            colors = self.alarmColorMap(colorMapInds,:);
        end
    end
end