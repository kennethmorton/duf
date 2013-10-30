function dufScoreUtilApplyPlotOptions(h,plotOptions)

f = fieldnames(plotOptions);
for i = 1:length(f);
    try
        set(h,f{i},plotOptions.(f{i}));
    catch ME
        
        fprintf('trouble setting field %s of object of type %s\n',f{i},get(h(1),'type'));
        %disp(ME.message);
    end
end