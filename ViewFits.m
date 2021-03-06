function  ViewFits(movfname,fits_fname,circ_D,write_mov,autoscale_on,linewidth)
%ViewFits plots or writes a view fits movie using the tiff stack movie
%specified by mov_fname, and the fits from Subtract_mol_off_frames mat file
%specficied by fits_fname.
%
% circ_D is the diameter of the circle to show the fit, default is 7
%
% write_mov is a boolean determining whether the viewfits movie will be
% written to an avi. If set to 0 this function will be used in debug mode
%
% autoscale_on is a boolean determining if the movie grayscale will
% be set frame by frame. If set to 0 a handful of frames throughout the movie
% are used to set the grayscale
% 
% linewidth is the linewidth of the circles in the movie. Default is 1
%
%%%% Color Scheme %%%%
% green circles are good fits
% red circles are bad fits 
% blue circles are identified as being on the nanoparticle (from GNR_enhancement)
% magenta circles are fits which passed the tracking filter

if nargin<3;circ_D=7;end
if nargin<4;write_mov=0;end
if nargin<5;autoscale_on=0;end
if nargin<6;linewidth=1;end
%%

%load in fits
load(fits_fname);

%create A `TIFFStack` object  which behaves like a read-only memory
%mapped TIFF file
tfstk=TIFFStack(movfname);
movsz=size(tfstk);%the size of the movie
[pathstr,name,~] = fileparts(movfname);

if write_mov
    v = VideoWriter([pathstr,filesep,name,'_ViewFits.avi'],'Uncompressed AVI');
    open(v);
    
    disp(['Making ViewFits for ',name]);
end

if ~autoscale_on
    frms4scl=100 ;%number of frames for the scaling
    %pull out some frames to find the percentiles of
    frmscl=double(tfstk(:,:,round(linspace(1,movsz(3),frms4scl))));
    % the intensity bounds for not autoscaling
    int_bounds=prctile(frmscl(frmscl>0),[.1,99.8]);
end

figure
set(gcf,'Position',[1281,1,1280,948])
for ii=1:movsz(3)
    
    curfrm=double(tfstk(:,:,ii));
    
    if autoscale_on
        int_bounds=prctile(curfrm(curfrm>0),[.1,99.8]);
    end
    imshow(curfrm,int_bounds)
    
    %which molecules appear in this frame
    thisfrm_gf=fits(fits(:,1)==ii & fits(:,9)==1,:);%goodfits
    thisfrm_bf=fits(fits(:,1)==ii & fits(:,9)==0,:);%badfits 
    
    if exist('identity_vec','var')
        thisfrm_ongnr=fits(fits(:,1)==ii & identity_vec~=0,:);%on gnr
    else
        thisfrm_ongnr=[];
    end 
    
    if exist('trk_filt','var')
        thisfrm_trk=fits(fits(:,1)==ii & trk_filt,:);%on gnr
    else
        thisfrm_trk=[];
    end 
    
    %note that for plotting row and column are switched...
    if ~isempty(thisfrm_gf)
        vcs=viscircles([thisfrm_gf(:,3),thisfrm_gf(:,2)],repmat(circ_D,[length(thisfrm_gf(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color','green')
    end
    if ~isempty(thisfrm_bf)
        vcs=viscircles([thisfrm_bf(:,3),thisfrm_bf(:,2)],repmat(circ_D,[length(thisfrm_bf(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color','red')
    end
    if ~isempty(thisfrm_ongnr)
        vcs=viscircles([thisfrm_ongnr(:,3),thisfrm_ongnr(:,2)],repmat(circ_D+2,[length(thisfrm_ongnr(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color','blue')
    end
    if ~isempty(thisfrm_trk)
        vcs=viscircles([thisfrm_trk(:,3),thisfrm_trk(:,2)],repmat(circ_D+4,[length(thisfrm_trk(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color','magenta')
    end
    
    if write_mov
        frame = getframe;
        writeVideo(v,frame);
    else
        title([name,' frame #',num2str(ii)],'Interpreter', 'none')
        keyboard
    end    
end

if write_mov
    close(v)
end
end

