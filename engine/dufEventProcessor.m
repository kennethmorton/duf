classdef dufEventProcessor < dufAction
    properties (SetAccess = protected)
        isSupervised = false;
        isCrossValidateValid = true;
    end
    
    methods
        function self = dufEventProcessor(varargin)
            
            self.classTrain = 'dufEventSet';
            self.classRun = 'dufEventSet';
            self.classRunRetained = true;
            
            self = prtUtilAssignStringValuePairs(self, varargin{:});
        end
    end
end