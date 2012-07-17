%s_InitISET
%
% It is often convenient to use ISET at the head of a script, even though
% it is not required.  This script 
%   * Closes down previous instantiations of ISET
%   * Clears the workspace
%   * Starts a fresh version
%   * Hides the main window
%
% We might want to clear out the ISET session file to prevent it from
% loading with the new ISET session.
%
% Copyright ImagEval Consultants, LLC, 2011.

if exist('vcSESSION','var')
  ieMainClose;
  clx
end
%% Initialize variables
% If the user already has a vcSESSION, clear it.  We are going to load the
% new one from the isetSession file.
clear vcSESSION;

% Define vc global variables and structures.  vc[UPPERCASE] are globals.
global vcSESSION;

% thisVersion =  1.001;    % Until October 11, 2005
% this Version = 1.01;     % Started Feb. 20, 2006
% thisVersion = 3.0;       % Started September, 2007
thisVersion = 4.0;         % Started August, 2009

ieSessionSet('version',thisVersion);

%% Check Matlab version number
% These are versions where we have run ISET 
expectedMatlabVersion = {'7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9'};  
version = ver('Matlab');
v = version.Version(1:3);
if ~ismember(v, expectedMatlabVersion)    % (matlabVersion ~= expectedMatlabVersion)
    str = sprintf('You are running version %s. ',version.Version);
    str = [str,'Supported versions {'];
    for ii=1:length(expectedMatlabVersion)
        str = [str,' ',expectedMatlabVersion{ii}];
    end
    str = [str,' }'];
    errordlg(str);
end

% disp(['ISET ',num2str(vcSESSION.VERSION),', Matlab ',version.Version]);

clear expectedMatlabVersion version matlabVersion

%% Default session file name.
% We check for a session file named iset-dateTime
%
%   * If one exists, we load it.  
%   * If several exist, we load the latest one
%   * The user can load a sessison file with a different name from the Main
%     Window.

d = dir('iset-*.mat');
if isempty(d)
    sessionFileName = 'isetSession.mat';
else
    sessionFileName = d(end).name;
end

sessionDir = pwd;
year = date; year = year(end-3:end);
% fprintf('------------------\n');
% fprintf('Copyright ImagEval Consulting, LLC, 2003-%s\n',year);
if exist([sessionDir,filesep,sessionFileName],'file')
    % The user has already created a vcSESSION file and we start with those
    % structures.  We adjust the session directory, however, to be consistent
    % with the current location of the file.
    cmd = ['load ',sessionFileName]; eval(cmd); clear cmd
    % disp(iePrintSessionInfo);
    
    ieSessionSet('dir',pwd);
    ieSessionSet('name',sessionFileName);
    
    % This is the version.  We check major release compatibility.
    if round(vcSESSION.VERSION) ~= round(thisVersion)
        error('isetSession file not current (%.03f). Delete the session file and restart.',thisVersion);
    end
    
    if ieSessionGet('initHelp'), iseHelp; end

    % We have the material from the last session.  So we are ready to jump
    % below and start opening windows based on the current state.
    
else
    % fprintf('No ISET session file found in %s.\nInitiating a session.\n',pwd);
    
    % Sets up the local working directories for an ISET session and guide
    % the user through the selection of default values.  Creates the
    % vcSESSION global variable.
    ieInitSession;
    % iseHelp;
    % At this point, we have null ISA, PIXEL and so forth defined.  We need to
    % open a window and start filling in these variables.  Which do we
    % start with? 
end

% wp = ieSessionGet('whitePoint');
% fprintf('White point for Scene/OI rendering set to %s\n',wp);

clear sessionFileName
clear thisVersion
clear sessionDir
clear v

% I think these can go away
% ieMainW
% ieMainW('visible','off')

% ---