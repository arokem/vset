function figHdl = vcNewGraphWin(figHdl, fType)
% Open a window for plotting
%
%    figHdl = vcNewGraphWin([fig handle],[figure type])
%
% A graph window figure handle is returned and stored in the currernt
% vcSESSION.GRAPHWIN entry.
%
% A few figure shapes can be defaulted
%   fType:  Default - Matlab normal figure position
%           upper left
%           tall 
%           wide
%   This list may grow.
%
% Examples
%  vcNewGraphWin;
%  vcNewGraphWin([],'upper left')   
%  vcNewGraphWin([],'tall')
%  vcNewGraphWin([],'wide')
%
% See also:
%
% Copyright ImagEval Consultants, LLC, 2005.

if ieNotDefined('figHdl'), figHdl = figure; end
if ieNotDefined('fType'),  fType = 'upper left'; end

set(figHdl,'Name','ISET GraphWin','NumberTitle','off');
set(figHdl,'CloseRequestFcn','ieCloseRequestFcn');
set(figHdl,'Color',[1 1 1]);

% Position the figure
fType = ieParamFormat(fType);
switch(fType)
    case 'upperleft'
        set(figHdl,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
    case 'tall'
        set(figHdl,'Units','normalized','Position',[0.007 0.055 0.28 0.85]);
    case 'wide'
        set(figHdl,'Units','normalized','Position',[0.007 0.62  0.7  0.3]);
    otherwise % default
end

ieSessionSet('graphwinfigure',figHdl);
ieSessionSet('graphwinhandle',guidata(figHdl));

return;
