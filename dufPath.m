function dufPath(varargin)

if isdeployed
    return;
end

P = genpath(dufRoot);
addpath(P);

%Remove some paths we don't need (we remove all directories that start with
% a . or a ]
removePath = [];
[string,remString] = strtok(P,pathsep);
while ~isempty(string);
    if ~isempty(strfind(string,[filesep '.'])) || ~isempty(strfind(string,[filesep ']']))
        removePath = cat(2,removePath,pathsep,string);
    end
    [string,remString] = strtok(remString,pathsep); %#ok
end
if ~isempty(removePath)
    rmpath(removePath);
end

% Add in the specified ] directories
for iArg = 1:length(varargin)
    cArg = varargin{iArg};
    cDir = fullfile(dufRoot,cat(2,']',cArg));
    assert(logical(exist(cDir,'file')),']%s is not a directory in %s',cArg,dufRoot);
    P = genpath(cDir);
    addpath(P);
    
    removePath = [];
    [string,remString] = strtok(P,pathsep);
    while ~isempty(string);
        if ~isempty(strfind(string,[filesep '.']))
            removePath = cat(2,removePath,pathsep,string);
        end
        [string,remString] = strtok(remString,pathsep); %#ok
    end
    if ~isempty(removePath)
        rmpath(removePath);
    end
end