classdef dufEventCreator < dufAction
    properties (SetAccess = protected)
        isSupervised = false;
        isCrossValidateValid = false;
    end
    methods
        function self = dufEventCreator(varargin)
            
            self.classTrain = 'dufCollection';
            self.classRun = 'dufCollection';
            self.classRunRetained = false; % Outputs a dufEventSet
            
            self = prtUtilAssignStringValuePairs(self,varargin{:});
        end
        
        function dsOut = run(self, dsIn)
            warning('We have things we need to add here in dufEventCreator.run().');
            % Add looping so that run is called on each file individually
            % Loop post run processing call to add collectionInfo if it is not already there.
            % any additional fields we might need?
            % We also need to make sure that files that yield no alarms are
            % included somehow.
            
            if ~self.isTrained
                self = train(self,dsIn);
            end
            dsOut = run@prtAction(self, dsIn);
        end
    
    end
end