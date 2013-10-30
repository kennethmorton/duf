classdef dufExcursionSet
    
    properties
        excursionArray
        
        scoring
        uniqueTargets
        uniqueDepths
    end
    
    methods
        
        function self = dufExcursionSet(varargin)
            self = prtUtilAssignStringValuePairs(self,varargin{:});
        end
        
        function self = loadFromFileListStruct(self,fileListStruct)
            % Generates a dufExclusionSet from a structure array where each
            % elment has fields: 
            %   .mdrFile
            %   .eleFile
            %   .ascFile
            %
            
            for i = 1:length(fileListStruct)
                as = dufAlarmSet;
                as = as.loadFromEle(fileListStruct(i).eleFile,fileListStruct(i).mdrFile);
                
                os = dufObjectSet;
                os = os.loadFromAsc(fileListStruct(i).ascFile);
                
                if i == 1
                    excursionArr = dufExcursion('alarmSet',as,'objectSet',os);
                else
                    excursionArr(i) = dufExcursion('alarmSet',as,'objectSet',os);
                end
            end
            self.excursionArray = excursionArr;
        end
        
        function self = score(self)
            
            ds = [];
            y = [];
            types = [];
            depths = [];
            areas = [];
            for i = 1:length(self.excursionArray)
                self.excursionArray(i) = self.excursionArray(i).score;
                ds = cat(1,ds,self.excursionArray(i).extractedConfidences);
                y = cat(1,y,self.excursionArray(i).extractedLabels);
                types = cat(1,types,self.excursionArray(i).extractedTypes);
                depths = cat(1,depths,self.excursionArray(i).extractedDepths);
                areas = cat(1,areas,self.excursionArray(i).scoring.area);
            end
            
            self.scoring.ds = ds;
            self.scoring.y = y;
            self.scoring.types = types;
            self.scoring.depths = depths;
            self.scoring.areas = areas;
            self.scoring.area = sum(areas);
        end
        
    end
end