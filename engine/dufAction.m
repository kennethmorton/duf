classdef dufAction < prtAction
    % We implement default train and run methods that do nothing.
    methods (Access = protected, Hidden = true)
        function  self = trainAction(self, ds) %#ok<INUSD>
        end
        function ds = runAction(self,ds) %#ok<INUSL>
        end
    end
    
    methods
        function dsOut = run(self, dsIn)
            if ~self.isTrained
                self = train(self,dsIn);
            end
            dsOut = run@prtAction(self, dsIn);
        end
    end
end