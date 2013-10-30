classdef dufCollection < prtDataSetClass
    % dufCollection
    % 
    %
    % % Example usage:
    % directoryMain = '//collins-vm-00.egr.duke.edu/kraken/data/forwardLooking/ASTO/';
    % runData = xlsToStructure(fullfile(directoryMain,'selexMetaData.csv'));
    % 
    % myRuns = runData(1);
    % dcSmall = dufCollection(myRuns);
    % 
    % dcBig = dufCollection(runData(1:10));
    
    methods
        function self = dufCollection(runData,varargin)
            self = self@prtDataSetClass(varargin{:});
            self.X = runData;
            self.observationInfo = runData; %why not?
        end
    end
    
    % Overloaded from prtDataSetClass
    methods 
        function Summary = summarize(Obj)
            Summary.nObservations = Obj.nObservations;
        end
    end
end