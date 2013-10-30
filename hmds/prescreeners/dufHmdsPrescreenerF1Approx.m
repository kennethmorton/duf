classdef dufHmdsPrescreenerF1Approx < dufEventCreator
    properties (SetAccess = private)
        name = 'F1 (Block)'
        nameAbbreviation = 'F1'
    end
    
    properties
        gprAlgo = gprF1BlockApproximation;
        threshold = 1.2;
    end
    
    methods
         function self = dufHmdsPrescreenerF1Approx(varargin)
            self = prtUtilAssignStringValuePairs(self,varargin{:});
         end
    end
    methods (Access = protected, Hidden = true)
        function eventDs = runAction(self,fileDs)
            % Do prescreening here
            
            allEvents = cell(fileDs.nObservations,1);
            for iFile = 1:fileDs.nObservations
                cFile = fileDs.X(iFile);
            
                % Actually prescreen
                [~, confidence] = processFile(self.gprAlgo,cFile.mdrFile);
                confidence = squeeze(confidence);
                
                confidence(:,1:30) = 0; 
                confidence(:,end-29:end) = 0;
                
                % From the 2D confidence map, make alarms
                pos = readMdrScanHeadersCache(cFile.mdrFile);
                posMat = as.niitekPos2PosMat(pos,size(confidence,1));
                
                options = as.optionsConfidence2AlarmSet;
                options.extractionFunction = @(C)as.thresholdedConnectedRegionMax(C,self.threshold);
                cAlarmSet = as.confidence2AlarmSetLandmines(confidence',posMat,options);
                
                % Interpret the output of old alarm score into events
                % The last two steps could be combined instead of
                % converting.
                cEvents = repmat(struct('confidence',[],'xUTM',[],'yUTM',[],'dt',[],'xt',[],'collectionInfo',[]),length(cAlarmSet.Alarms),1);
                for iAlarm = 1:length(cAlarmSet.Alarms)
                    cAlarm = cAlarmSet.Alarms(iAlarm);
                    
                    cEvents(iAlarm).confidence = cAlarm.confidence;
                    cEvents(iAlarm).xUTM = cAlarm.Extent.parameters(1);
                    cEvents(iAlarm).yUTM = cAlarm.Extent.parameters(2);
                    cEvents(iAlarm).dt = cAlarm.Info.downTrack;
                    cEvents(iAlarm).xt = cAlarm.Info.crossTrack;
                    cEvents(iAlarm).collectionInfo = cFile;
                end
                
                allEvents{iFile} = cEvents;
            end
            
            eventDs = dufEventSetHmds(cat(1,allEvents{:}));
        end
    end
end