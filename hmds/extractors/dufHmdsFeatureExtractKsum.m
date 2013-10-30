classdef dufHmdsFeatureExtractKsum < dufEventProcessor
    properties (SetAccess = private)
        name = 'KSum'
        nameAbbreviation = 'KSUM'
    end
    
    properties
        
    end
    
    methods
         function self = dufHmdsFeatureExtractKsum(varargin)
            self = prtUtilAssignStringValuePairs(self,varargin{:});
         end
    end
    methods (Access = protected, Hidden = true)
        function eventDs = runAction(self,eventDs)
            % Do prescreening here
            
            eventDs.X = randn(eventDs.nObservations,2);
            
            %gprAlignPolyfit + gprHaircut + gprCfar
        end
    end
end