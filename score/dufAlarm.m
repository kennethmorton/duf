classdef dufAlarm
    
    properties
        extent
        confidence
        
        onObject
        object
        
        %Storage for all the other crap
        info
        sensor
        
        plotOptions = struct('markersize',12,'color','r','marker','.');
    end
    
    methods
       
        function self = dufAlarm(varargin)
            self = prtUtilAssignStringValuePairs(self,varargin{:});
        end
        
        function h = plot(self)
            h = plot(self.extent);
            dufScoreUtilApplyPlotOptions(h,self.plotOptions);
        end
        
        function s = toString(self,escapeLatex)
            if nargin < 2
                escapeLatex = true;
            end
            if isempty(self.onObject)
                s = sprintf('Conf: %.2f; @ [%.2f %.2f]; False Alarm',self.confidence,self.extent.point(1),self.extent.point(2));
            else
                s = sprintf('Conf: %.2f; @ [%.2f %.2f]; On %s @ %d',self.confidence,self.extent.point(1),self.extent.point(2),self.object.name,self.object.depth);
            end
            if escapeLatex
                s = strrep(s,'_','\_');
            end
        end
    end
    
    methods (Static)
        function alarm = fromAlarmStruct(alarmStruct)
            %translate from a structure with fields "confidence" and
            % "location" to a dufAlarm.

            alarm = dufAlarm;
            alarm.confidence = alarmStruct.confidence;
            
            alarm.extent = extentInstantaneous('point',alarmStruct.location);
        end
    end
    
end