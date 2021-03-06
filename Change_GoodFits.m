function Change_GoodFits(directoryname,new_maxerr,new_stdtol)
%% Change_GoodFits 
% This is a function to change the goodfits vector in the fits array. If
% you don't want to change either maxerr or stdtol, then enter the
% corresponding input (new_maxerr or new_stdtol) as []. See either
% Subtract_then_fit or BGSUBFIT User Guide for an explanation of what
% maxerr and stdtol are.
%
% directoryname   is the name of the directory where the fits .mat files
% will be selected, if there is an error finding the directory the program
% will open uigetfile in the current working directory
%
%  new_maxerr is the new max_err that you want to use. If you don't want to
%  change this and only change the stdtol, then enter [] for new_maxerr
%
%  new_stdtol is the new stdtol that you want to use. If you don't want to
%  change this and only change the max_err, then enter [] for new_stdtol
%
%% Select the fits .mat files
display('Select the fit .mat files.')
try
    [datalist,dataloc,findex]=uigetfile([directoryname filesep '*.mat*'],'multiselect','on');
catch
    curdir=pwd;
    [datalist,dataloc,findex]=uigetfile([curdir filesep '*.mat*'],'multiselect','on');
end
if findex==0
    fprintf('no data selected\n')
    return
end
if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);

%% Make the new good_fit vector

h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,'Changing GoodFits');

for ii=1:numel(dlocs);
    try; waitbar(ii/numel(dlocs),h1); end
    
    fits_fname=[dlocs{ii},filesep,dnames{ii}];
    load(fits_fname,'fits','stdtol','maxerr','dfrlmsz')%load in the necessary variables
    % change the parameters if they were sent
    if ~isempty(new_maxerr);maxerr=new_maxerr;end
    if ~isempty(new_stdtol);stdtol=new_stdtol;end
    
    %the conversion between dfrlmsz and the STD of the Gaussian, reccomended
    %using the full width at 20% max given by (2*sqrt(2*log(5)))
    dfD2std=(2*sqrt(2*log(5)));
    %the guessed std
    gesss=dfrlmsz/dfD2std;
    
    %%% make the new goodfits vector %%%
    goodfits=(fits(:,4)<=(stdtol*gesss) & fits(:,4)>=(gesss/stdtol)) & ...
        all(fits(:,[2,3,6,8])>=0,2) & fits(:,7)<maxerr;
    
    %update the fits array
    fits(:,9)=goodfits;
    %save it
    save(fits_fname,'fits','stdtol','maxerr','dfrlmsz','-append')%append the new vector and the used parameters
end
end

