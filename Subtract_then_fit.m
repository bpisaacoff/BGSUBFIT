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
% 1. frame number, 2. x (px),3. y (px), 4. width (px) ,5. offset,
% 6. amplitude, 7. error, 8. sum(:), 9. goodfit boolean

%%%% Dependencies %%%%
% TIFFStack
% MLEwG (for MLE fitting)
% gaussfit (for least squares fitting)

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
%% Import the data
% load off frames list and some parameters
load(Mol_off_frames_fname,'off_frames','dfrlmsz','moloffwin')

%create A `TIFFStack` object  which behaves like a read-only memory
%mapped TIFF file
tfstk=TIFFStack(mov_fname);
movsz=size(tfstk);%the size of the movie
[pathstr,fname,~] = fileparts(mov_fname);

% load the guesses
load(guessfname,'guesses');

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

%fit info is [frame number,x,y,width,offset,amplitude,variance,sum(:),goodfit boolean]
fits=NaN(size(guesses,1),9);

h1=waitbar(0);
waitbar(0,h1,['Fitting ',fname]);
set(findall(h1,'type','text'),'Interpreter','none');
for ii=1:size(guesses,1)   
    try
        waitbar(ii/size(guesses,1),h1)
    catch
    end
        
    curfrmnum=guesses(ii,1);
    
    %determine the frame list of frames to check for the current frame
    if curfrmnum<=(moloffwin/2)%the first group of frames
        frmlst=curfrmnum+(-(curfrmnum-1):(moloffwin/2));
    elseif curfrmnum>=(movsz(3)-moloffwin/2)%the last group of frames
        frmlst=curfrmnum+((-moloffwin/2):(movsz(3)-curfrmnum));
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
            mov=cat(3,mov(:,:,(numnewfrmsbeg+1):(moloffwin+1)),double(tfstk(:,:,curframes(end)+(-(numnewfrmsend-1):0))));
        else
            mov=double(tfstk(:,:,curframes));
            warning('There was a big jump in the frames without a fit')
        end
        %check to make sure that everything is the right size
        if curfrmnum>(moloffwin) && (length(curframes)~=(moloffwin+1) || size(mov,3)~=(moloffwin+1))
            error('Error determining the correct frames to import')
        end
    end
    
    %current molecule's position
    molx=guesses(ii,2);
    moly=guesses(ii,3);
    
    %checking that it's not outside the frame and that off_frames for this
    %guess isn't empty
    if (moly>edgedist && moly<(movsz(2)-edgedist) && molx>edgedist && molx<(movsz(1)-edgedist)) && ...
            (moly>dfrlmsz && moly<(movsz(2)-dfrlmsz) && molx>dfrlmsz && molx<(movsz(1)-dfrlmsz))&& ... 
            ~isempty(off_frames{ii})
        %the average frame
        mean_mov=mean(mov(molx+(-dfrlmsz:dfrlmsz),moly+(-dfrlmsz:dfrlmsz),off_frames{ii}-frmlst(1)+1),3);
        %the molecule image
        molim=mov(molx+(-dfrlmsz:dfrlmsz),moly+(-dfrlmsz:dfrlmsz),curfrmnum-frmlst(1)+1);
        %the subtracted image
        data=molim-mean_mov;        
        
        %%%% Fitting %%%%
        plot_on=0;%for debugging purposes only!
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
                if plot_on;figure(11);end
                [fitPars,conf95,~,~]=gaussFit(data,'searchBool',0,'nPixels',2*dfrlmsz+1,'checkVals',plot_on);
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
        act_x=paramsF(1)-dfrlmsz-1+molx;
        act_y=paramsF(2)-dfrlmsz-1+moly;
        
        %%%% Fit Checks %%%%
        % fits is [frame number,x,y,width, offset,amplitude,err,sum(:),goodfit boolean]
        if (paramsF(3)<=(stdtol*params0(3)) && paramsF(3)>=(params0(3)/stdtol)) && ... %Compare width with diffraction limit
                ~errbad && ... %too much error on fit?
                (paramsF(1)<=(params0(1)+mxdst) && paramsF(1)>=(params0(1)-mxdst)) && ... %check x position
                (paramsF(2)<=(params0(2)+mxdst) && paramsF(2)>=(params0(2)-mxdst)) && ...  %check y position
                ~any([paramsF([1,2,3,5]),sumsum]<0) %none of the fitted parameters should be negative, except the offset!
            
            %Put the results into the array
            fits(ii,:)=[curfrmnum,act_x,act_y,paramsF(3:6),sumsum,1];            
        else
            %track the guesses that don't get fit, for debugging/viewfits
            %purposes
            fits(ii,:)=[curfrmnum,act_x,act_y,paramsF(3:6),sumsum,0];
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

fits_col_headers={'frame num','x pos (px)','y pos (px)','sigma (px)','offset','N','error','sum(:)','goodfit boolean'};

%save the data
save([pathstr,filesep,fname,'_AccBGSUB_fits.mat'],'fits','fits_col_headers','mov_fname','Mol_off_frames_fname','guessfname',...
    'MLE_fit','maxdistfrac','stdtol','maxerr','dfrlmsz')

try
    close(h1)
end
end

