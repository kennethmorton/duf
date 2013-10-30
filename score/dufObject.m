classdef dufObject
    
    properties
        extent
        
        name
        depth
        class 
        
        %these could all theoretically be different classes of things
        isThreat
        isClutter
        isMine
        %%isHole %not implemented, b/c storage of object types is a mess
        
        metalContent
        
        objectPlotOptions = struct('edgecolor',[.3 .3 .3],'facecolor',[.3 .3 .3],'linewidth',2);
        haloPlotOptions = struct('color',[.3 .3 .3],'linewidth',1,'linestyle','--');
        textPlotOptions = struct('fontweight','bold');
        %Storage for all the other crap
        info
    end
    
    methods
        
        function [h,hText] = plot(self)
            h = plot(self.extent);
            dufScoreUtilApplyPlotOptions(h,self.objectPlotOptions);
            c = self.extent.getCenter;
            hText = text(c(1)+.3,c(2),self.toStringEscapeUnderscores);
        end
        
        
        function s = toStringEscapeUnderscores(self)
            s = toString(self);
            s = strrep(s,'_','\_');
        end
        
        function s = toString(self)
            s = sprintf('%s @ %d',self.name,self.depth);
        end
        
        function h = plotHalo(self,halo)
            h = plotHalo(self.extent,halo);
            dufScoreUtilApplyPlotOptions(h,self.haloPlotOptions);
        end
        
    end
    
    methods (Static)
        function dufObj = fromObjectStruct(objectStruct)
            %translate from an old objectStruct to a dufObject
            %
            
            dufObj = dufObject;
            dufObj.name = objectStruct.Info.objectID;
            dufObj.depth = objectStruct.Info.depth;
            
            dufObj.class = objectStruct.Info.objectInfo.purpose;
            dufObj.metalContent = objectStruct.Info.objectInfo.metalContent;
                        
            dufObj.isThreat = ~strcmpi(dufObj.class,'clutter');
            dufObj.isClutter = strcmpi(dufObj.class,'clutter');
            dufObj.isMine = objectStruct.Info.objectInfo.isMine;
            
            dufObj.info.numID = objectStruct.Info.numID;
            dufObj.info.isSquare = objectStruct.Info.isSquare;
            dufObj.info.squareCorners = objectStruct.Info.squareCorners;
            dufObj.info.squareCenter = objectStruct.Info.squareCenter;
            
            dufObj.info.laneNumber = objectStruct.Info.laneNumber;
            
            switch lower(objectStruct.Extent.shape)
                case 'circle'
                    dufObj.extent = extentCircle('radius',objectStruct.Extent.parameters.radius,'center',objectStruct.Extent.parameters.center);
                case 'polygon'
                    %dufObj.extent = extentPolygon('polygon',objectStruct.Extent.parameters);
                    dufObj.extent = extentPolygonConvex('polygon',objectStruct.Extent.parameters);
                otherwise
                    error('unknown type');
            end
            
        end
    end
    
end