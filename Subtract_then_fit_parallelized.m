function  Subtract_then_fit(mov_fname,Mol_off_frames_fname,guessfname,MLE_fit,edgedist,maxdistfrac,stdtol,maxerr)
%% Subtract_mol_off_frames
% subtracts the average intensity of off frames for each guess
% stored in Mol_off_frames_fname.

%%%% Inputs %%%%
% mov_fname the filename of the tiff stack movie

% Mol_off_frames_fname is the filename .mat file output from the function
% Mol_off_frames

% guessfname is the filename for the guesses .mat file

% MLE_fit  a Boolean determining whether or not MLE fitting is used. Set to
% 1 to use MLE and to 0 to use least squares. Default is 0. Note that MLE
% is quite slow, and so its nots recommended for a large number of guesses

% edgedist is the distance in pixels from the edge of the frame to ignore.
% default is 10

% maxdistfrac is max x&y distance of fit from guess,as a fraction of
% dfrlmsz, default is 0.75

% stdtol is tolerance on fit Gaussian STD, to leae filtering options for
% later, default value is 5

% maxerr is the maximum error of the fit for MLE fit, using variance default
% 0.1 (can't be above this) for LSQR fit, using the 95% confidence interval
% on the position, default max is 3

%%%% Output %%%%
% a .mat file, importantly containing the fits array with columns:
% 1. frame number, 2. row pos (px),3. column pos (px), 4. width (px) ,5. offset,
% 6. amplitude, 7. error, 8. sum(:), 9. goodfit boolean

%%%% Dependencies %%%%
% TIFFStack
% MLEwG (for MLE fitting)
% gaussfit (for least squares fitting)

%  Copyright 2016 Benjamin P Isaacoff
%
% Licensed under the Apache License, Version 2.0 (the "License"); you
% may not use this file except in compliance with the License. You may
% obtain a copy of the License at
%
%   http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
% implied. See the License for the specific language governing
% permissions and limitations under the License.


if nargin<4;MLE_fit=0;end
if nargin<5;edgedist=10;end
if nargin<6;maxdistfrac=0.75;end
if nargin<7;stdtol=5;end
if nargin<8;
    if MLE_fit
        maxerr=0.1;
    else
        maxerr=3;
    end
end

tic;%for measuring the time to run the entire program
%% Import the data
%note some of the weird syntax in this section is due to the parfor loop
%later on
% load off frames list and some parameters
dataload=load(Mol_off_frames_fname);
dfrlmsz=dataload.dfrlmsz;off_frames=dataload.off_frames;moloffwin=dataload.moloffwin;

%create A `TIFFStack` object  which behaves like a read-only memory
%mapped TIFF file
tfstk=TIFFStack(mov_fname);
movsz=size(tfstk);%the size of the movie
[pathstr,fname,~] = fileparts(mov_fname);

% load the guesses
dataload=load(guessfname);
guesses=dataload.guesses;

clear dataload

%check number of fits vs length of off frames
if size(guesses,1)~=numel(off_frames);error('Unequal number of fits and number of off frames lists');end

% import the first moloffwin+1 frames
curframes=1:(moloffwin+1);
mov=double(tfstk(:,:,curframes));

%% The Averaging and Subtraction

%the conversion between dfrlmsz and the STD of the Gaussian, reccomended
%using the full width at 20% max given by (2*sqrt(2*log(5)))
dfD2std=(2*sqrt(2*log(5)));

%the guessed std
gesss=dfrlmsz/dfD2std;

%max x&y distance of fit from guess
mxdst=maxdistfrac*dfrlmsz;

%fit info is [frame number,row,col,width,offset,amplitude,variance,sum(:),goodfit boolean]
fits=NaN(size(guesses,1),9);

h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,['Fitting ',fname]);
for jj=1:movsz(3)
    try; waitbar(jj/movsz(3),h1); end
    
    curfrmnum=jj;
    
    %determine the frame list of frames to check for the current frame
    if curfrmnum<=(moloffwin/2)%the first group of frames
        frmlst=curfrmnum+(-(curfrmnum-1):(moloffwin/2));
    elseif curfrmnum>=(movsz(3)-moloffwin/2)%the last group of frames
        frmlst=movsz(3)+(-moloffwin:0);
    else %all the frames in the middle
        frmlst=curfrmnum+((-moloffwin/2):(moloffwin/2));
    end
    
    %import appropriate movie frames
    if frmlst(end)>curframes(end)
        numnewfrmsend=frmlst(end)-curframes(end);
        numnewfrmsbeg=frmlst(1)-curframes(1);
        %new current frames list
        curframes=frmlst;
        %new movie frames
        if numnewfrmsend<(moloffwin+1)
            mov=cat(3,mov(:,:,(numnewfrmsbeg+1):(moloffwin+1)),...
                double(tfstk(:,:,curframes(end)+(-(numnewfrmsend-1):0))));
        else
            mov=double(tfstk(:,:,curframes));
            warning('There was a big jump in the frames without a fit')
        end
        %check to make sure that everything is the right size
        if curfrmnum>(moloffwin) && (length(curframes)~=(moloffwin+1) || size(mov,3)~=(moloffwin+1))
            error('Error determining the correct frames to import')
        end
    end
    
    guessessjj=guesses(guesses(:,1)==curfrmnum,:);
    gessinds=find(guesses(:,1)==curfrmnum);
    %preallocating the fits array for the parfor loop
    curfits=NaN(size(guessessjj,1),9);
    
    
    parfor ii=1:size(guessessjj,1)
        
        %current molecule's position
        molr=guessessjj(ii,2);
        molc=guessessjj(ii,3);
        
        %checking that it's not outside the frame and that off_frames for this
        %guess isn't empty
        if (molc>edgedist && molc<(movsz(2)-edgedist) && molr>edgedist && molr<(movsz(1)-edgedist)) && ...
                (molc>dfrlmsz && molc<(movsz(2)-dfrlmsz) && molr>dfrlmsz && molr<(movsz(1)-dfrlmsz))&& ...
                ~isempty(off_frames{gessinds(ii)})
            %the average frame
            mean_mov=mean(mov(molr+(-dfrlmsz:dfrlmsz),molc+(-dfrlmsz:dfrlmsz),off_frames{gessinds(ii)}-frmlst(1)+1),3);
            %the molecule image
            molim=mov(molr+(-dfrlmsz:dfrlmsz),molc+(-dfrlmsz:dfrlmsz),curfrmnum-frmlst(1)+1);
            %the subtracted image
            data=molim-mean_mov;
            
            %%%% Fitting %%%%
            plot_on=0;%for debugging purposes only! Note doesn't work while in a parfor loop
            %the guessed tail intensity of the gaussian
            gessb=min(data(:));
            %the guessed amplitude, using the formula in MLEwG
            gessN=range(data(:))*(4*pi*gesss^2);
            %fit guess vector
            params0=[dfrlmsz,dfrlmsz,gesss,gessb,gessN];
            %The sum(:) of the the data
            sumsum=sum(data(:));
            
            if MLE_fit
                %fitting with MLE
                [paramsF,varianceF] = MLEwG (data,params0,1,plot_on,1);
                paramsF=[paramsF,varianceF];
                %recalculating the values based on their equations to match
                paramsF(5)=paramsF(5)*(2*pi*paramsF(3)^2);
                paramsF(4)=sqrt(paramsF(4));
                errbad=varianceF>maxerr;%too much error on fit?
            else
                %fitting with least squares
                try
                    [fitPars,conf95,~,~]=gaussFit(data,'searchBool',0,'nPixels',2*dfrlmsz+1,'checkVals',0);
                catch
                    conf95=[inf,inf];
                    fitPars=[0,0,0,0,0,0];
                    warning('gaussfit error')
                end
                %converting the variables to match the output of MLEwG
                paramsF=[fitPars(1),fitPars(2),fitPars(3),fitPars(5),...
                    fitPars(4),mean(conf95([1,2]))];
                errbad=mean(conf95([1,2]))>maxerr;%too much error on fit?
            end
            %Convert back into full frame coordinates, NOTE the -1!
            act_r=paramsF(1)-dfrlmsz-1+molr;
            act_c=paramsF(2)-dfrlmsz-1+molc;
            
            %%%% Fit Checks %%%%
            % fits is [frame number,row pos,col pos,width, offset,amplitude,err,sum(:),goodfit boolean]
            if (paramsF(3)<=(stdtol*params0(3)) && paramsF(3)>=(params0(3)/stdtol)) && ... %Compare width with diffraction limit
                    ~errbad && ... %too much error on fit?
                    (paramsF(1)<=(params0(1)+mxdst) && paramsF(1)>=(params0(1)-mxdst)) && ... %check row position
                    (paramsF(2)<=(params0(2)+mxdst) && paramsF(2)>=(params0(2)-mxdst)) && ...  %check col position
                    ~any([paramsF([1,2,3,5]),sumsum]<0) %none of the fitted parameters should be negative, except the offset!
                
                %Put the results into the array
                curfits(ii,:)=[curfrmnum,act_r,act_c,paramsF(3:6),sumsum,1];
            else
                %track the guesses that don't get fit, for debugging/viewfits
                %purposes
                curfits(ii,:)=[curfrmnum,act_r,act_c,paramsF(3:6),sumsum,0];
            end
             %debugging
            if plot_on
                figure(12)
                subplot(1,3,1)
                imshow(mean_mov,[])
                title('Mean BG')
                subplot(1,3,2)
                imshow(molim,[])
                title('Raw Molecule')
                subplot(1,3,3)
                imshow(data,[])
                title('BGSUB')
                
                keyboard
            end
        end
    end
    fits(gessinds,:)=curfits;    
end

fits_col_headers={'frame num','row pos (px)','column pos (px)','sigma (px)','offset','N','error','sum(:)','goodfit boolean'};
tictoc=toc;%the time to run the entire program
%save the data
save([pathstr,filesep,fname,'_AccBGSUB_fits.mat'],'fits','fits_col_headers','mov_fname','Mol_off_frames_fname','guessfname',...
    'MLE_fit','maxdistfrac','stdtol','maxerr','dfrlmsz','tictoc')

try
    close(h1)
end
end

