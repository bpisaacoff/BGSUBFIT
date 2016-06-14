function BGSUBFit(directoryname,dfrlmsz,avgwin,moloffwin,varargin)
%% BGSUBFit
% BGSUBFit is the wrapper function to do real background substraction of
% single molecule imaging movies so that the molecules can be fit and their
% intensity accurately measured.
%
%
%
%%%% Inputs %%%%
%%% required
%   directoryname   is the name of the directory where the movies will be
% selected, if there is an error finding the directory the program will
% open uigetfile in the current working directory
%
%   dfrlmsz   is the size of a diffraction limited spot in pixels. It's the
% nominal diameter, NOT the FWHM or something similar. Must be an integer!
% For an expected diffraction limited standard deviation, std, using the
% full width at 20% max, dfrlmsz = std*(2*sqrt(2*log(5)))
%
%   avgwin   is the window size (in frames) to be used for the average
% subtraction. Needs to be an odd integer
%
%   moloffwin   is the window size (in frames) to be checked to determine
%   in which frames that molecule was off and are thus safe to subtract.
%   Needs to be an even integer
%
%%% optional
% See the list of optional parameters in parameters section below. They are
% called with a name value pair as is standard for Matlab functions. e.g.,
% 'bpthrsh',90 would set the parameter bpthrsh=90
%
%
%%%% Outputs %%%%

%%%% Dependencies %%%%
% AVGSUB_tiffs
% saveastiff
% TIFFStack
% Guessing
% bpass
% Mol_off_frames
% MLEwG (for MLE fitting)
% gaussfit (for least squares fitting)
% Subtract_then_fit
% Track_3D2
% Tracking
% Track_filter
% hungarian

%%%%
% Written by Benjamin P Isaacoff at the University of Michigan
% last update 6/11/16 BPI

%% parameters & defaults
% You are of course welcome to change the default values, but I would
% strongly urge you to instead set them as inputs to the function using a
% name value pair. The default parameter values are:

% actions
% Do guessing
params.makeGuesses = 1;
% check guesses
params.check_guesses = 0;
% Make the off frames list
params.makeOffFrames = 1;
% Do fitting
params.fitting = 1;
% Do tracking
params.tracking = 1;
% make the ViewFits movie
params.makeViewFits = 1;

%AVGSUB running average boolean
params.runningavg = 1;
%AVGSUB offset
params.offset = 1000;

%Guessing parameters
params.bpthrsh = 90;
%how many pixels to ignore around the edge of the frame
params.egdesz = dfrlmsz;
%compare brightnesses in each frame? if not, use entire movie
params.pctile_frame = 0;

%Fitting parameters
% do MLE fitting? If not least squares fitting will be used
params.MLE_fit = 0;
% Goodfit parameters. See Subtract_mol_off_frames for the details
params.maxdistfrac = 0.75;
params.stdtol = 5;
params.maxerr = 3; % if you change this please also change the statement after the next loop

%Tracking parameters
% minimum merit
params.trackparams(1)=0.01;
% Integration time (ms)
params.trackparams(2)=200;
% gamma
params.trackparams(3)=1;
% maximum step size
params.trackparams(4)=3;
% minimum track length
params.trackparams(5)=3;
% speed estimation window halfsize
params.trackparams(6)=1;
% time delay between consecutive frames (ms)
params.trackparams(7)=0;

%ViewFits parameters
% use the original movie? if not, use the avgsub movie
params.orig_movie = 1;
% diameter of the circles showing the fits
params.circ_D = dfrlmsz;
% linewidth of the circles
params.linewidth = 1;
% write a .avi movie showing the fits. If not, goes to debug mode
params.write_mov=1;
% autoscale frame by frame?
params.autoscale_on = 1;


paramsnames=fieldnames(params);

% if any sim parameters are included as inputs, change the parameter
% mentioned
if nargin>4
    for ii=1:2:nargin-5
        whichField = strcmp(paramsnames,varargin{ii});
        if all(~whichField)
            warning('Check spelling. Parameter change may have not occurred.')
        end
        eval(['params.' paramsnames{whichField} ' = varargin{ii+1};'])
    end
end

%changing the default error if doing MLE fitting
if params.MLE_fit && params.maxerr==3
    params.maxerr=0.1; 
end

%% Select the movies
display('Select the movies.')
try
    [datalist,dataloc,findex]=uigetfile([directoryname filesep '*.tif*'],'multiselect','on');
catch
    curdir=pwd;
    [datalist,dataloc,findex]=uigetfile([curdir filesep '*.tif*'],'multiselect','on');
end
if findex==0
    fprintf('no data selected\n')
    return
end
if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);

%% The Average Subtraction
% check which of the selected movies already have the avgsub tif stack
% made, if not then make it. If you want to remake the avgsub movie, then
% delete it or movie it from the directory

% AVGSUB_tiffs will save an average subtracted .tif stack and write a .txt
% file with the parameters, called moviename_bgsub_info.txt
if params.makeGuesses
    bFiles=dir([dataloc,'*_avgsub.tif']);%list all the avgsub tif stacks
    for ii=1:numel(dlocs);
        %compare the chosen files with the avgsub list and choose ones
        %without an avgsub movie
        if ~any(ismember({bFiles.name},[dnames{ii},'_avgsub.tif']))&&...
                ~any(ismember({bFiles.name},dnames{ii}))
            %do the avgsub
            AVGSUB_tiffs([dlocs{ii},filesep,dnames{ii},'.tif'],params.runningavg,avgwin,params.offset);
        end
    end
end

%% Make Guesses
% loop through each movie and make the guesses, which will be saved in a
% .mat file called moviename_guesses.mat
if params.makeGuesses
    for ii=1:numel(dlocs);
        Guessing([dlocs{ii},filesep,dnames{ii},'_avgsub.tif'],dfrlmsz,...
            params.bpthrsh,params.egdesz,params.pctile_frame,params.check_guesses);
    end
end

%% Make the off frames list
% loop through all of the movies and using the guesses .mat file will write
% the off frames list to a .mat file with _Mol_off_frames appended to the
% name
if params.makeOffFrames
    for ii=1:numel(dlocs);
        Mol_off_frames([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses'],dfrlmsz,moloffwin);
    end
    
end

%% Subtract and fit
% loop through the movies and fit the subtracted images using the off
% frames list. Out puts a .mat file with AccBGSUB_fits appended
if params.fitting
    for ii=1:numel(dlocs);
        Subtract_then_fit([dlocs{ii},filesep,dnames{ii},'.tif'],...
            [dlocs{ii},filesep,dnames{ii},'_avgsub_guesses_Mol_off_frames'],...
            [dlocs{ii},filesep,dnames{ii},'_avgsub_guesses'],...
            params.MLE_fit,params.egdesz,params.maxdistfrac,params.stdtol,params.maxerr);
    end
end

%% Tracking
% loop through each movie and does tracking and appends the tracks to the
% fits .mat file. Also will append a logical vector called trk_filt which
% indicates if the fit passed was successfully tracked and wasn't the first
% or last frame in a track
if params.tracking
    for ii=1:numel(dlocs);
        Track_filter([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits'],1,params.trackparams);
    end
end

%% Make the ViewFits movie
% loop through each movie and make a ViewFits movie, or just go into debug
% mode, to look at the results. Outpits an avi file called
% moviename_ViewFits.avi
if params.makeViewFits
    for ii=1:numel(dlocs);
        if params.orig_movie
            VF_fname=[dlocs{ii},filesep,dnames{ii},'.tif'];
        else
            VF_fname=[dlocs{ii},filesep,dnames{ii},'_avgsub.tif'];
        end
         ViewFits(VF_fname,[dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits'],...
             params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
    end
end

end






























