function varargout = iSim(varargin)
% ISIM MATLAB code for iSim.fig
%      ISIM, by itself, creates a new ISIM or raises the existing
%      singleton*.
%
%      H = ISIM returns the handle to a new ISIM or the handle to
%      the existing singleton*.
%
%      ISIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ISIM.M with the given input arguments.
%
%      ISIM('Property','Value',...) creates a new ISIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iSim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iSim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iSim

% Last Modified by GUIDE v2.5 19-Feb-2024 08:56:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iSim_OpeningFcn, ...
                   'gui_OutputFcn',  @iSim_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before iSim is made visible.
function iSim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iSim (see VARARGIN)

% Choose default command line output for iSim
handles.output = hObject;

global WAVE_LENGTH
global DESIGNED_RI
global NA
global Angle
global k
global PX_RES
PX_RES = 0.075; % [um]

WAVE_LENGTH = 0.445;
DESIGNED_RI = 1.52;
NA = 1.42;
Angle = asin(NA/DESIGNED_RI);
k = 2.*pi/WAVE_LENGTH;

global tg
global tg0
global ti
global ti0
tg = 170;
tg0 = 170;
ti = 120;
ti0 = 120;

global ng
global ng0
global ni
global ni0
global np
global ns
ni = 1.50;
ni0 = 1.52;
ng = 1.52;
ng0 = 1.52;

np = 1.62;
ns = 1.40;

global ROI_SIZE_MIN
global ROI_SIZE_MAX

ROI_SIZE_MIN = 25;
ROI_SIZE_MAX = 25;

global ROI_SIZE

global radius_mean
global radius_min
global radius_max
radius_mean = 0.01;
radius_min = radius_mean - (radius_mean*0.05);
radius_max = radius_mean + (radius_mean*0.05);

global radius
radius = radius_mean;
set(handles.radius,'String',sprintf('%.3f',radius));

global zf
zf = 2.0;
set(handles.zf,'String',sprintf('%.3f',zf));

global zp
zp = 0.0;
set(handles.zp,'String',sprintf('%.3f',zp));

global zp_low
global zp_high
zp_low = 0.0;
zp_high = 5.0;

global zp_low_limit
global zp_high_limit
zp_low_limit = 0.0;    
zp_high_limit = 5.0;

global VISCOSITY
VISCOSITY = 1;
set(handles.VISCOSITY,'String',sprintf('%d',VISCOSITY));

global NP
NP = 50;
set(handles.NP,'String',sprintf('%d',NP));

global NumFrame
NumFrame = 1;
set(handles.NumFrame,'String',sprintf('%d',NumFrame));

global FOV_HEIGHT
global FOV_WIDTH
FOV_HEIGHT = 128;
FOV_WIDTH = 128;
set(handles.FOV_HEIGHT,'String',sprintf('%d',FOV_HEIGHT));
set(handles.FOV_WIDTH,'String',sprintf('%d',FOV_WIDTH));

global OUT_HEIGHT
global OUT_WIDTH
OUT_HEIGHT  = 128;
OUT_WIDTH = 128;
set(handles.OUT_HEIGHT,'String',sprintf('%d',OUT_HEIGHT));
set(handles.OUT_WIDTH,'String',sprintf('%d',OUT_WIDTH));

global OUT_IMAGE
OUT_IMAGE = zeros(OUT_HEIGHT,OUT_WIDTH);

global CANVAS_HEIGHT
global CANVAS_WIDTH
global MARGIN_SCALE
MARGIN_SCALE = 2;
CANVAS_HEIGHT = OUT_HEIGHT + ROI_SIZE_MAX * MARGIN_SCALE;
CANVAS_WIDTH = OUT_WIDTH + ROI_SIZE_MAX * MARGIN_SCALE;

global CANVAS_IMAGE
CANVAS_IMAGE = zeros(CANVAS_HEIGHT,CANVAS_WIDTH);

global SCAT_AMP
SCAT_AMP = 1*10^9;

global BACK_INT
global BACK_RANDOM_P

BACK_INT = 0.005;
BACK_RANDOM_P = 3;
set(handles.BACK_INT,'String',sprintf('%.3f',BACK_INT));
set(handles.BACK_RANDOM_P,'String',sprintf('%d',BACK_RANDOM_P));

global FLAG_DISPLAY_PSF
FLAG_DISPLAY_PSF = 1;
set(handles.CHECK_DISPLAY_PSF,'Value',FLAG_DISPLAY_PSF);

global FLAG_DISPLAY_CANVAS
FLAG_DISPLAY_CANVAS = 1;
set(handles.CHECK_DISPLAY_CANVAS,'Value',FLAG_DISPLAY_CANVAS);

global FLAG_DISPLAY_NOISE
FLAG_DISPLAY_NOISE = 1;
set(handles.CHECK_DISPLAY_NOISE,'Value',FLAG_DISPLAY_NOISE);

global FLAG_DISPLAY_OUT
FLAG_DISPLAY_OUT = 1;
set(handles.CHECK_DISPLAY_OUT,'Value',FLAG_DISPLAY_OUT);

global TICK_TIME
global EXPOSURE_TIME
TICK_TIME = 600;
EXPOSURE_TIME = 1200;
set(handles.TICK_TIME,'String',sprintf('%d',TICK_TIME));
set(handles.EXPOSURE_TIME,'String',sprintf('%d',EXPOSURE_TIME));

global WITH_DAT_FILE
WITH_DAT_FILE = 1;
set(handles.WITH_DAT_FILE,'Value',WITH_DAT_FILE);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iSim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = iSim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in Random.
function Random_Callback(hObject, eventdata, handles)
% hObject    handle to Random (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global NP
% STR = get(handles.NP,'String'); NP = str2num(STR);
% 
% global zp_low
% global zp_high
% 
% global CANVAS_HEIGHT
% global CANVAS_WIDTH
% global ROI_SIZE_MAX
% global ROI_SIZE_MIN
% 
% xStart_rand_pos = zeros(NP);
% yStart_rand_pos = zeros(NP);
% zStart_rand_pos = zeros(NP);
% 
% xMinLimit= ROI_SIZE_MAX/2 + 1;
% xMaxLimit= CANVAS_WIDTH+ROI_SIZE_MAX - ROI_SIZE_MAX/2 - 1;
% xStart_rand_pos = (xMaxLimit-xMinLimit).*rand(NP,1) + xMinLimit;
% 
% yMinLimit= ROI_SIZE_MAX/2 + 1;
% yMaxLimit= CANVAS_HEIGHT+ROI_SIZE_MAX - ROI_SIZE_MAX/2 - 1;
% yStart_rand_pos = (yMaxLimit-yMinLimit).*rand(NP,1) + yMinLimit;
% 
% global zf
% zMinLimit= zp_low;
% zMaxLimit= zp_high;
% zStart_rand_pos = (zMaxLimit-zMinLimit).*rand(NP,1) + zMinLimit;
% 
% global radius
% global radius_min
% global radius_max
% global radius_mean
% radius_rand = (radius_max-radius_min).*rand(NP,1) + radius_mean;
% 
% global NumFrame
% STR = get(handles.NumFrame,'String'); NumFrame = str2num(STR);
% 
% 
% global TICK_TIME
% global EXPOSURE_TIME
% STR = get(handles.EXPOSURE_TIME,'String'); EXPOSURE_TIME = str2num(STR);
% STR = get(handles.TICK_TIME,'String'); TICK_TIME = str2num(STR);
% NUM_FRAME_ACCUMULATION = uint16(EXPOSURE_TIME/TICK_TIME);
% RealNumTracePoint = NumFrame * NUM_FRAME_ACCUMULATION;
% 
% 
% 
% 
% global zf
% for pid = 1:NP
%     radius = radius_rand(pid);
%     
%     p = double(floor(200 / radius));
%      
%     r = xMinLimit + (xMaxLimit - xMinLimit)*sum(rand(RealNumTracePoint,p),2)/p;
%     r = r - mean(r);
%     x_px(:,pid)  = cumsum(r) + xStart_rand_pos(pid);
% 
%     r = yMinLimit + (yMaxLimit - yMinLimit)*sum(rand(RealNumTracePoint,p),2)/p;
%     r = r - mean(r);
%     y_px(:,pid)  = cumsum(r) + yStart_rand_pos(pid);
%     
%     zMinLimit= zp_low-zf;
%     zMaxLimit= zp_high-zf;
% 
%     r = (zMaxLimit - zMinLimit)*sum(rand(RealNumTracePoint,p),2)/p + zStart_rand_pos(pid); 
%     r = r - mean(r);
%     z_um(:,pid)  = cumsum(r) + zStart_rand_pos(pid);
% end
% 
% for pid = 1:NP
%     row = [];
%     [row col] = find(z_um(:,pid)<0);
%     if(isempty(row) ~= 1)
%         z_um(:,pid) = z_um(:,pid) + abs(min(z_um(:,pid)));
%     end
% end

global NP
STR = get(handles.NP,'String'); NP = str2double(STR);

global zp_low
global zp_high

global CANVAS_HEIGHT
global CANVAS_WIDTH
global ROI_SIZE_MAX
global ROI_SIZE_MIN


global xStart_rand_pos
global yStart_rand_pos
global zStart_rand_pos

% xStart_rand_pos = zeros(NP);
% yStart_rand_pos = zeros(NP);
% zStart_rand_pos = zeros(NP);

% xMinLimit= ROI_SIZE_MAX/2 + 1;
% xMaxLimit= CANVAS_WIDTH+ROI_SIZE_MAX - ROI_SIZE_MAX/2 - 1;
% xStart_rand_pos = (xMaxLimit-xMinLimit).*rand(NP,1) + xMinLimit;
% 
% yMinLimit= ROI_SIZE_MAX/2 + 1;
% yMaxLimit= CANVAS_HEIGHT+ROI_SIZE_MAX - ROI_SIZE_MAX/2 - 1;
% yStart_rand_pos = (yMaxLimit-yMinLimit).*rand(NP,1) + yMinLimit;
% 
% zMinLimit= zp_low;
% zMaxLimit= zp_high;
% zStart_rand_pos = (zMaxLimit-zMinLimit).*rand(NP,1) + zMinLimit;


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

global NumFrame
STR = get(handles.NumFrame,'String'); NumFrame = str2double(STR);

global TICK_TIME
global EXPOSURE_TIME
global TIME_PERIOD

STR = get(handles.TICK_TIME,'String'); TICK_TIME = str2double(STR);
STR = get(handles.EXPOSURE_TIME,'String'); EXPOSURE_TIME = str2double(STR);

NUM_FRAME_ACCUMULATION = uint16(EXPOSURE_TIME/TICK_TIME);
RealNumTracePoint = uint32(NumFrame) * uint32(NUM_FRAME_ACCUMULATION);

global zf
global MARGIN_SCALE

global PX_RES
global VISCOSITY
STR = get(handles.VISCOSITY,'String'); VISCOSITY = str2double(STR);

xMinLimit= ROI_SIZE_MAX/2 + 1;
xMaxLimit= CANVAS_WIDTH+ROI_SIZE_MAX*MARGIN_SCALE - ROI_SIZE_MAX*MARGIN_SCALE/2 - 1;
xStart_rand_pos = (xMaxLimit-xMinLimit).*rand(NP,1) + xMinLimit;
xStart_rand_pos(1,1) = (xMaxLimit-xMinLimit)/2;

yMinLimit= ROI_SIZE_MAX/2 + 1;
yMaxLimit= CANVAS_HEIGHT+ROI_SIZE_MAX*MARGIN_SCALE - ROI_SIZE_MAX*MARGIN_SCALE/2 - 1;
yStart_rand_pos = (yMaxLimit-yMinLimit).*rand(NP,1) + yMinLimit;
yStart_rand_pos(1,1) = (yMaxLimit-yMinLimit)/2;

zMinLimit= zp_low;
zMaxLimit= zp_high;
zStart_rand_pos = (zMaxLimit-zMinLimit).*rand(NP,1) + zMinLimit;
zStart_rand_pos(1,1) = (zMaxLimit-zMinLimit)/2;

global radius
global radius_min
global radius_max
global radius_mean
global radius_rand

radius_rand = (radius_max-radius_min).*randn(NP,1) + radius_mean;
% histogram(handles.dist_radius,radius_rand);
% xlabel(handles.dist_radius,'radius (um)');
% drawnow

Dcoeff = 200./(VISCOSITY.*(1000*radius_rand));
Delta = sqrt(2*Dcoeff.*TICK_TIME*(10^-6))';

% histogram(handles.dist_dcoeff_cal,Dcoeff);
% xlabel(handles.dist_dcoeff_cal,'D (um^{2}/s)');
% drawnow
% x_pos(1,:) = PX_RES.*xStart_rand_pos';
% y_pos(1,:) = PX_RES.*yStart_rand_pos';
% z_pos(1,:) = zStart_rand_pos';
% 
% Localization_error = 0.01; % 10nm
% 
% for idx = 2: RealNumTracePoint
%     x_pos(idx,:)  = x_pos(idx-1,:) + Delta.*randn(1,NP) + Localization_error.*randn(1,NP);
%     y_pos(idx,:)  = y_pos(idx-1,:) + Delta.*randn(1,NP) + Localization_error.*randn(1,NP);
%     z_pos(idx,:)  = z_pos(idx-1,:) + Delta.*randn(1,NP) + Localization_error.*randn(1,NP);
% 
% %     for n=1:NP
% %         if(z_pos(idx,n) < 0)
% %             z_pos(idx,n) = abs(z_pos(idx,n));
% %         end
% %     end
% end

x_pos = zeros(NumFrame,NP);
y_pos = zeros(NumFrame,NP);
z_pos = zeros(NumFrame,NP);

xDelta = zeros(NumFrame-1,NP);
yDelta = zeros(NumFrame-1,NP);
zDelta = zeros(NumFrame-1,NP);


x_pos(1,:) = PX_RES.*xStart_rand_pos';
y_pos(1,:) = PX_RES.*yStart_rand_pos';
z_pos(1,:) = zStart_rand_pos';

for idx = 1:RealNumTracePoint - 1
    xDelta(idx,:) = Delta.*randn(1,NP); %+ Localization_error.*randn(1,NP);
    yDelta(idx,:) = Delta.*randn(1,NP); %+ Localization_error.*randn(1,NP);
    zDelta(idx,:) = Delta.*randn(1,NP); %+ Localization_error.*randn(1,NP); 
end

for idx = 2: RealNumTracePoint    
    for n=1:NP
        x_pos(idx,n)  = x_pos(idx-1,n) + xDelta(idx-1,n);
        y_pos(idx,n)  = y_pos(idx-1,n) + yDelta(idx-1,n);

        z_pos(idx,n)  = z_pos(idx-1,n) + zDelta(idx-1,n);
        if(z_pos(idx,n) < zp_low)
            z_pos(idx,n)  = zp_low + radius_rand(n,1);
        elseif(z_pos(idx,n) > zp_high)
            z_pos(idx,n)  = zp_high - radius_rand(n,1);
        end
    end
end

if(RealNumTracePoint >= 2)
    hFig = figure(2);
    set(hFig,'Position',[100,100,500,300]);

    subplot(2,3,1);
    for n=1:length(x_pos(1,:))
        plot(x_pos(:,n));
        hold on
    end
    hold off
    ylim([0 CANVAS_WIDTH*PX_RES]);
    grid on
    subplot(2,3,2);
    for n=1:length(y_pos(1,:))
        plot(y_pos(:,n));
        hold on
    end
    hold off
    ylim([0 CANVAS_WIDTH*PX_RES]);
    grid on
    subplot(2,3,3);
    for n=1:length(z_pos(1,:))
        plot(z_pos(:,n));
        hold on
    end
    hold off
    ylim([zp_low CANVAS_WIDTH*PX_RES]);
    grid on
    subplot(2,3,4);
    for n=1:length(xDelta(1,:))
        plot(xDelta(:,n));
        hold on
    end
    hold off
    grid on
    subplot(2,3,5);
    for n=1:length(yDelta(1,:))
        plot(yDelta(:,n));
        hold on
    end
    hold off
    grid on
    subplot(2,3,6);
    for n=1:length(zDelta(1,:))
        plot(zDelta(:,n));
        hold on
    end
    hold off
    grid on

    image_name = sprintf('Trace plots.png');
    saveas(hFig,image_name,'png');
    image_name = sprintf('Trace plots.fig');
    saveas(hFig,image_name,'fig');
end


global x_px
global y_px
global z_um

x_px = zeros(RealNumTracePoint,NP);
y_px = zeros(RealNumTracePoint,NP);
z_um = zeros(RealNumTracePoint,NP);

for idx = 1:RealNumTracePoint
    x_px(idx,:) = x_pos(idx,:)./PX_RES;
    y_px(idx,:) = y_pos(idx,:)./PX_RES;
    z_um(idx,:) = z_pos(idx,:);
end

hFig = figure(1);
set(hFig,'Position',[100,100,300,200]);
subplot(2,1,1); histogram(radius_rand);
xlabel ('Radius (um)'); ylabel('Counts'); grid on
subplot(2,1,2); histogram(Dcoeff);
xlabel ('D (um^{2}/s)'); ylabel('Counts'); grid on



% --- Executes on button press in Simulation.
function Simulation_Callback(hObject, eventdata, handles)
% hObject    handle to Simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Simulation, 'Userdata', []);

global zp
global CANVAS_IMAGE

global CANVAS_IMAGE_NOISE_ADD

global PSF_ROI
global ROI_SIZE

global OUT_HEIGHT
global OUT_WIDTH
global FOV_HEIGHT
global FOV_WIDTH

global OUT_IMAGE

YS_OUT_IMAGE = OUT_HEIGHT/2 - round(FOV_HEIGHT/2)+1;
YE_OUT_IMAGE = OUT_HEIGHT/2 + round(FOV_HEIGHT/2);
XS_OUT_IMAGE = OUT_WIDTH/2 - round(FOV_WIDTH/2)+1;
XE_OUT_IMAGE = OUT_WIDTH/2 + round(FOV_WIDTH/2);

global FLAG_DISPLAY_CANVAS


% =========================================================================
TextInfoName = sprintf('Setting information.txt');
fWriteTXT = fopen(TextInfoName,'wt');

global WAVE_LENGTH
global DESIGNED_RI
global NA
global Angle
global k
global PX_RES

fprintf(fWriteTXT,'Wavelength (um) = %.3f\n',WAVE_LENGTH);
fprintf(fWriteTXT,'RI(Objective) = %.2f\n',DESIGNED_RI);
fprintf(fWriteTXT,'NA = %.2f\n',NA);
fprintf(fWriteTXT,'Pixel resolution (um) = %.3f\n',PX_RES);

global tg
global tg0
global ti
global ti0
fprintf(fWriteTXT,'(tg, tg0) = %.3f, %.3f\n',tg, tg0);
fprintf(fWriteTXT,'(ti, ti0) = %.3f, %.3f\n',ti, ti0);

global ng
global ng0
global ni
global ni0
fprintf(fWriteTXT,'(ng, ng0) = %.3f, %.3f\n',ng, ng0);
fprintf(fWriteTXT,'(ni, ni0) = %.3f, %.3f\n',ni, ni0);

global np
global ns
fprintf(fWriteTXT,'(np, ns) = %.3f, %.3f\n',np, ns);

global radius_mean
global radius_min
global radius_max
fprintf(fWriteTXT,'radius(min,mean,max) = %.3f, %.3f, %.3f\n',radius_min, radius_mean, radius_max);

global VISCOSITY
fprintf(fWriteTXT,'VISCOSITY = %d\n',VISCOSITY);

global ROI_SIZE_MIN
global ROI_SIZE_MAX
fprintf(fWriteTXT,'ROI(min,max) = %d. %d\n',ROI_SIZE_MIN,ROI_SIZE_MAX);

global zf
fprintf(fWriteTXT,'zf (um) = %d. %d\n',zf);

global zp
global zp_low
global zp_high
fprintf(fWriteTXT,'zp (min,mean,max) = %.1f, %.1f, %.1f\n',zp_low, zp, zp_high);

global zp_low_limit
global zp_high_limit
fprintf(fWriteTXT,'zp limits(low,high) = %.3f, %.3f\n',zp_low_limit,zp_high_limit);

global NP
fprintf(fWriteTXT,'number of particles = %d\n',NP);

global NumFrame
fprintf(fWriteTXT,'number of franes = %d\n',NumFrame);

global OUT_WIDTH
global OUT_HEIGHT
fprintf(fWriteTXT,'Image size (width, hight) = %d, %d\n',OUT_WIDTH,OUT_HEIGHT);

global TICK_TIME
global EXPOSURE_TIME
fprintf(fWriteTXT,'Time (tick, exposure) = %d, %d\n',TICK_TIME,EXPOSURE_TIME);

global BACK_INT
global BACK_RANDOM_P
global SCAT_AMP
fprintf(fWriteTXT,'Back_int, SCAT Amp. = %.3f, %e\n',BACK_INT,SCAT_AMP);

fclose(fWriteTXT);


% =========================================================================
global CANVAS_HEIGHT
global CANVAS_WIDTH


FPS = round(1/(EXPOSURE_TIME*10^-6));
SIZE = uint16(radius_mean*1000);
Out_FileName = sprintf('SIMULATION-R%d-E%d-NP%d-FPS%d-TD.dat',SIZE,VISCOSITY,NP,FPS);
Out_TD = fopen(Out_FileName,'wb');


NUM_FRAME_ACCUMULATION = uint16(EXPOSURE_TIME/TICK_TIME);
RealNumTracePoint = uint32(NumFrame) * uint32(NUM_FRAME_ACCUMULATION);


% MovieFileName = sprintf('CANVAS-SIMULATION-R%d-E%d-NP%d-FPS%d-TD.avi',SIZE,VISCOSITY,NP,FPS);
% myVideo_a = VideoWriter(MovieFileName);     % Create videowrite handle
% myVideo_a.FrameRate = 30;         % Set a frame per seconds 
% myVideo_a.Quality = 100;                    % Set quality
% open(myVideo_a);

MovieFileName = sprintf('CANVAS-NOISE-SIMULATION-R%d-E%d-NP%d-FPS%d-TD.avi',SIZE,VISCOSITY,NP,FPS);
myVideo = VideoWriter(MovieFileName);     % Create videowrite handle
myVideo.FrameRate = 30;         % Set a frame per seconds 
myVideo.Quality = 100;                    % Set quality
open(myVideo);


global radius_rand

global x_px
global y_px
global z_um

global zp_low_limit
global zp_high_limit


global xc_adj
global yc_adj
global size_mean

global WITH_DAT_FILE

if(NumFrame == 1)
    MovieFileName = sprintf('Canvas generation example.avi');
    myVideo_gen_example = VideoWriter(MovieFileName);     % Create videowrite handle
    myVideo_gen_example.FrameRate = 2;         % Set a frame per seconds 
    myVmyVideo_gen_exampleideo.Quality = 100;                    % Set quality
    open(myVideo_gen_example);
end

global MARGIN_SCALE

global FrameIdx

global OUT_IMAGE
global FLAG_DISPLAY_OUT

for FrameIdx=1:NumFrame
    
    set(handles.FrameIdx,'String',sprintf('%d',FrameIdx));
    drawnow

    TextSaveName = sprintf("Learning data %d.txt",FrameIdx);
    fWrite = fopen(TextSaveName,'wt');
    fprintf(fWrite,'xs\tys\txb\tyb\n');
    
    pCount = 0;

    global data_xy
    data_xy = [];
    
    CANVAS_IMAGE = zeros(CANVAS_HEIGHT,CANVAS_WIDTH);
    for pid = 1:NP
        CANVAS_IMAGE_TEMP = zeros(CANVAS_HEIGHT,CANVAS_WIDTH);
        
        GenImgCount = 0;
        for GenImgIdx = (FrameIdx-1)*NUM_FRAME_ACCUMULATION+1:(FrameIdx-1)*NUM_FRAME_ACCUMULATION+NUM_FRAME_ACCUMULATION
            radius = radius_rand(pid);        
            xc = uint16(x_px(GenImgIdx,pid));
            yc = uint16(y_px(GenImgIdx,pid));
            zp = z_um(GenImgIdx,pid);
            
            ROI_SIZE = round( ROI_SIZE_MIN + abs(abs(zp)*(ROI_SIZE_MAX - ROI_SIZE_MIN)*(5/ROI_SIZE_MAX) ) );
           
            if(rem(ROI_SIZE,2) == 0)
               ROI_SIZE = ROI_SIZE + 1;
            end
            if(ROI_SIZE >= ROI_SIZE_MAX) 
                ROI_SIZE = ROI_SIZE_MAX;
            end
            
            fprintf('%.3f, %d\n',zp, ROI_SIZE);

            GenImgCount = GenImgCount + 1;
            xc_acc(GenImgCount,1) = xc;
            yc_acc(GenImgCount,1) = yc;
            zp_acc(GenImgCount,1) = zp;
            roi_size_acc(GenImgCount,1) = ROI_SIZE;
            
            set(handles.zp,'String',sprintf('%.3f',zp));
            set(handles.radius,'String',sprintf('%.3f',radius));
            
            
            Generate_Callback(hObject, eventdata, handles);

            Amp = PSF_ROI(round(ROI_SIZE/2),round(ROI_SIZE/2));
            fprintf('%d: %d | (%d, %d, %.3f) | size: %d | r:%.3f | A:%.3f\n',FrameIdx, pid, xc, yc, zp, ROI_SIZE, radius, Amp);
            
            
            [CANVAS_IMAGE_TEMP] = ROI_MERGE(yc,xc,PSF_ROI,ROI_SIZE,CANVAS_IMAGE_TEMP,CANVAS_HEIGHT,CANVAS_WIDTH);

            CANVAS_IMAGE = CANVAS_IMAGE + CANVAS_IMAGE_TEMP;
        end
        

        xc_adj = mean(xc_acc) - ROI_SIZE_MAX * MARGIN_SCALE/2;
        yc_adj = mean(yc_acc) - ROI_SIZE_MAX * MARGIN_SCALE/2;
        size_mean = mean(roi_size_acc);
          
        if( (xc_adj - round(size_mean/2)+1 >= XS_OUT_IMAGE) && ...
            (xc_adj + round(size_mean/2)   <= XE_OUT_IMAGE) && ...
            (yc_adj - round(size_mean/2)+1 >= YS_OUT_IMAGE) && ...
            (yc_adj + round(size_mean/2)   <= YE_OUT_IMAGE) && ...
            (zp >= zp_low_limit) && ...
            (zp <= zp_high_limit) )
        
            fprintf(fWrite,'%d\t%d\t%d\t%d\n',...
                uint8(xc_adj - round(size_mean/2))+1,...
                uint8(yc_adj - round(size_mean/2))+1,...
                uint8(size_mean),...
                uint8(size_mean));

            pCount = pCount + 1;
            data_xy(pCount,1) = uint8(xc_adj)+1;
            data_xy(pCount,2) = uint8(yc_adj)+1;
            data_xy(pCount,3) = uint8(size_mean);


        end
        
        %fprintf('%d: %d, (%d, %d, %.3f), %d\n',FrameIdx, pid, xc, yc, zp, size_mean);
        %Generate_Callback(hObject, eventdata, handles);
        %[CANVAS_IMAGE_TEMP] = ROI_MERGE(yc,xc,PSF_ROI,ROI_SIZE,CANVAS_IMAGE_TEMP,CANVAS_HEIGHT,CANVAS_WIDTH);
        %CANVAS_IMAGE = CANVAS_IMAGE + CANVAS_IMAGE_TEMP;
        
        if(FLAG_DISPLAY_CANVAS == 1)
            CANVAS_IMAGE_NORM = uint8(iNORM_2D_ARRAY(CANVAS_IMAGE,1,255));
            CANVAS_IMAGE_NORM_RGB = cat(3, CANVAS_IMAGE_NORM, CANVAS_IMAGE_NORM, CANVAS_IMAGE_NORM);
            imshow((CANVAS_IMAGE_NORM_RGB),'Parent',handles.CANVAS);
            drawnow

            if(NumFrame == 1)
                writeVideo(myVideo_gen_example,uint8(CANVAS_IMAGE_NORM_RGB));
            end

        end
    end
    fclose(fWrite);
    CANVAS_IMAGE = CANVAS_IMAGE./double(NUM_FRAME_ACCUMULATION);
    
    global CROP_CANVAS
    CROP_CANVAS = zeros(FOV_HEIGHT,FOV_WIDTH);
    YS_CANVAS = CANVAS_HEIGHT/2 - round(FOV_HEIGHT/2)+1;
    YE_CANVAS = CANVAS_HEIGHT/2 + round(FOV_HEIGHT/2);
    XS_CANVAS = CANVAS_WIDTH/2 - round(FOV_WIDTH/2)+1;
    XE_CANVAS = CANVAS_WIDTH/2 + round(FOV_WIDTH/2);
    CROP_CANVAS = CANVAS_IMAGE(YS_CANVAS: YE_CANVAS, XS_CANVAS: XE_CANVAS);
    
    global CANVAS_IMAGE_CROP
    CANVAS_IMAGE_CROP = zeros(OUT_HEIGHT,OUT_WIDTH);
    YS_OUT_IMAGE = OUT_HEIGHT/2 - round(FOV_HEIGHT/2)+1;
    YE_OUT_IMAGE = OUT_HEIGHT/2 + round(FOV_HEIGHT/2);
    XS_OUT_IMAGE = OUT_WIDTH/2 - round(FOV_WIDTH/2)+1;
    XE_OUT_IMAGE = OUT_WIDTH/2 + round(FOV_WIDTH/2);
    CANVAS_IMAGE_CROP(YS_OUT_IMAGE: YE_OUT_IMAGE, XS_OUT_IMAGE: XE_OUT_IMAGE) = CROP_CANVAS; 

    CANVAS_IMAGE_NORM = uint8(iNORM_2D_ARRAY(CANVAS_IMAGE_CROP,1,255));
    CANVAS_IMAGE_NORM_RGB = cat(3, CANVAS_IMAGE_NORM, CANVAS_IMAGE_NORM, CANVAS_IMAGE_NORM);
    %writeVideo(myVideo_a,uint8(CANVAS_IMAGE_NORM_RGB));
    
    CANVAS_IMAGE_CROP = imgaussfilt(CANVAS_IMAGE_CROP,1);

    ADD_NOISE_Callback(hObject, eventdata, handles);
    
    global CANVAS_IMAGE_NOISE_CROP
    CANVAS_IMAGE_NOISE_CROP = CANVAS_IMAGE_NOISE_ADD(YS_CANVAS: YE_CANVAS, XS_CANVAS: XE_CANVAS);
    CANVAS_IMAGE_NOISE_ADD_CROP(YS_OUT_IMAGE: YE_OUT_IMAGE, XS_OUT_IMAGE: XE_OUT_IMAGE) = CANVAS_IMAGE_NOISE_CROP; 

    CANVAS_IMAGE_NOISE_ADD_NORM = uint8(iNORM_2D_ARRAY(CANVAS_IMAGE_NOISE_ADD_CROP,1,255));
    CANVAS_IMAGE_NOISE_ADD_NORM_RGB = cat(3, CANVAS_IMAGE_NOISE_ADD_NORM, CANVAS_IMAGE_NOISE_ADD_NORM, CANVAS_IMAGE_NOISE_ADD_NORM);
    CANVAS_IMAGE_NOISE_ADD_NORM_MARK_RGB = CANVAS_IMAGE_NOISE_ADD_NORM_RGB;

    CANVAS_IMAGE_NORM_RGB = insertShape(CANVAS_IMAGE_NORM_RGB,'Rectangle',[1  1  FOV_WIDTH FOV_HEIGHT],'LineWidth', 5,'Color',[68 146 197]);
    CANVAS_IMAGE_NOISE_ADD_NORM_RGB = insertShape(CANVAS_IMAGE_NOISE_ADD_NORM_RGB,'Rectangle',[1  1  FOV_WIDTH FOV_HEIGHT],'LineWidth', 5,'Color',[236 52 36]);
    CANVAS_IMAGE_NOISE_ADD_NORM_MARK_RGB = insertShape(CANVAS_IMAGE_NOISE_ADD_NORM_MARK_RGB,'Rectangle',[1  1  FOV_WIDTH FOV_HEIGHT],'LineWidth', 5,'Color',[255 255 0]);

    if(isempty(data_xy) == 0)
        PEAK_RectArea = zeros(pCount,4);
        PEAK_RectArea(:,4) = data_xy(:,3); %BOX_SIZE;
        PEAK_RectArea(:,3) = data_xy(:,3); %BOX_SIZE;
        PEAK_RectArea(:,2) = floor(data_xy(:,2) - data_xy(:,3)/2);
        PEAK_RectArea(:,1) = floor(data_xy(:,1) - data_xy(:,3)/2);
        
        CANVAS_IMAGE_NOISE_ADD_NORM_MARK_RGB = insertShape(CANVAS_IMAGE_NOISE_ADD_NORM_MARK_RGB,'Rectangle',...
        [PEAK_RectArea(:,1)  PEAK_RectArea(:,2)  PEAK_RectArea(:,3) PEAK_RectArea(:,4)],...
        'LineWidth', 1,'Color','yellow');
    end
    
    

    if(FLAG_DISPLAY_OUT == 1)
        imshow((CANVAS_IMAGE_NOISE_ADD_NORM_MARK_RGB),'Parent',handles.IMAGE_OUT_MARK);
        drawnow
    end

    MOVIE_BUFFER(1:FOV_HEIGHT, FOV_WIDTH*0+1 :FOV_WIDTH*1, : ) = CANVAS_IMAGE_NORM_RGB;
    MOVIE_BUFFER(1:FOV_HEIGHT, FOV_WIDTH*1+1 :FOV_WIDTH*2, : ) = CANVAS_IMAGE_NOISE_ADD_NORM_RGB;
    MOVIE_BUFFER(1:FOV_HEIGHT, FOV_WIDTH*2+1 :FOV_WIDTH*3, : ) = CANVAS_IMAGE_NOISE_ADD_NORM_MARK_RGB;

    writeVideo(myVideo,uint8(MOVIE_BUFFER));
    
    CROP_IMAGE_Callback(hObject, eventdata, handles);
    
    

%     OUT_IMAGE = imgaussfilt(OUT_IMAGE,1);
    OUT_IMAGE_NORM = uint8(iNORM_2D_ARRAY(OUT_IMAGE,1,255));
    OUT_IMAGE_NORM_RGB = cat(3, OUT_IMAGE_NORM, OUT_IMAGE_NORM, OUT_IMAGE_NORM);
    
    out_image_name = sprintf("Learning data %d.png",FrameIdx);
    imwrite(uint8(OUT_IMAGE_NORM_RGB),out_image_name);

    fwrite(Out_TD,OUT_WIDTH,'*uint32');
    fwrite(Out_TD,OUT_HEIGHT,'*uint32');
    BIN = 2;
    fwrite(Out_TD,BIN,'*uint32');
    BYTE_TEMP = 4;
    fwrite(Out_TD,BYTE_TEMP,'*uint32');
    PCOIndex = FrameIdx;
    fwrite(Out_TD,PCOIndex,'*uint32');
    PCO2Index = 0;
    fwrite(Out_TD,PCO2Index,'*uint32');
    ANDORIndex = 0;
    fwrite(Out_TD,ANDORIndex,'*uint32');
    Reserved1 = 0;
    fwrite(Out_TD,Reserved1,'*uint32');
    Reserved2 = 0;
    fwrite(Out_TD,Reserved2,'*uint32');
    Reserved3 = 0;
    fwrite(Out_TD,Reserved3,'*uint32');
    Reserved4 = 0;
    fwrite(Out_TD,Reserved4,'*uint32');
    RemotePos = 0;
    fwrite(Out_TD,RemotePos,'*double');
    StageXPos = 0;
    fwrite(Out_TD,StageXPos,'*double');
    StageYPos = 0;
    fwrite(Out_TD,StageYPos,'*double');
    StageZPos = 0;
    fwrite(Out_TD,StageZPos,'*double');

    ALL_CONTRASTS = reshape(OUT_IMAGE.',1,[])';
    fwrite(Out_TD,ALL_CONTRASTS,'*single');

    ALL_CONTRASTS = reshape(OUT_IMAGE.',1,[])';
    fwrite(Out_TD,ALL_CONTRASTS,'*single');

    ALL_CONTRASTS = reshape(OUT_IMAGE.',1,[])';
    fwrite(Out_TD,ALL_CONTRASTS,'*single');
    

    curval = get(handles.Simulation, 'UserData');
    if ~isempty(curval) && strcmp(curval, 'stop')
        break;
    end

    if(NumFrame == 1)
        close(myVideo_gen_example);
    end
end

fclose(Out_TD);
close(myVideo);




% --- Executes on button press in EXPORT.
function EXPORT_Callback(hObject, eventdata, handles)
% hObject    handle to EXPORT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FrameIdx

global OUT_HEIGHT
global OUT_WIDTH
global PX_RES

global xc_adj
global yc_adj
global size_mean

global data_xy

global CANVAS_IMAGE_NOISE_CROP
global OUT_IMAGE

NUM_OF_PARTICLES = length(data_xy(:,1));

OUT_DATA = [];
Count = 0;
for y = 1:OUT_HEIGHT
    for x=1:OUT_WIDTH
        Count = Count + 1;
        OUT_DATA(Count,1) = x-1;
        OUT_DATA(Count,2) = (x-1).*PX_RES;
        OUT_DATA(Count,3) = OUT_HEIGHT-y;
        OUT_DATA(Count,4) = (OUT_HEIGHT-y).*PX_RES;
    end
end
ValueLinearArray = reshape((OUT_IMAGE).',1,[])';
OUT_DATA(:,5) = ValueLinearArray;

for n = 1:NUM_OF_PARTICLES
    Count = 0;
    for y = 1:OUT_HEIGHT
        for x=1:OUT_WIDTH
            Count = Count + 1;
            
            if(round(data_xy(n,1)) == x && round(data_xy(n,2)) == y)
                OUT_DATA(Count,6) = (x-2);
                OUT_DATA(Count,7) = (OUT_HEIGHT-y+1);
                OUT_DATA(Count,8) = (x-2).*PX_RES;
                OUT_DATA(Count,9) = (OUT_HEIGHT-y+1).*PX_RES;
            end
    
        end
    end
end

ML_SEG_MARK_CANVAS = zeros(OUT_HEIGHT,OUT_WIDTH);

for k=1:NUM_OF_PARTICLES
    BOX_SIZE = data_xy(k,3);
    xS = round(data_xy(k,1) - floor(BOX_SIZE/2))-1;
    xE = round(data_xy(k,1) + floor(BOX_SIZE/2))-1;
    yS = round(data_xy(k,2) - floor(BOX_SIZE/2))-1;
    yE = round(data_xy(k,2) + floor(BOX_SIZE/2))-1;

    if(yS<1) 
        yE = BOX_SIZE+yS+1;
        yS = 1; 
    end
    if(xS<1)
        xS = BOX_SIZE+xS+1;
        xS = 1; 
    end
    if(yE>OUT_HEIGHT)
        yS = OUT_HEIGHT +1 - BOX_SIZE; yE = OUT_HEIGHT;
    end
    if(xE>OUT_WIDTH)
        xS = OUT_WIDTH + 1 - BOX_SIZE; xE = OUT_WIDTH;
    end
    
    ML_SEG_MARK_CANVAS(yS:yE,xS:xS) = 1; % Draw y-line on xS
    ML_SEG_MARK_CANVAS(yS:yE,xE:xE) = 1; % Draw y-line on xE
    ML_SEG_MARK_CANVAS(yS:yS,xS:xE) = 1; % Draw x-line on yS
    ML_SEG_MARK_CANVAS(yE:yE,xS:xE) = 1; % Draw x-line on yE
end
ValueLinearArray = reshape((ML_SEG_MARK_CANVAS).',1,[])';
OUT_DATA(:,10) = ValueLinearArray;


Export_name = sprintf('[DATAGRAPH] Learning DATA %d.txt',FrameIdx);
save (Export_name, 'OUT_DATA', '-ascii');


% --- Executes on button press in STOP.
function STOP_Callback(hObject, eventdata, handles)
% hObject    handle to STOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Simulation, 'Userdata', 'stop');



% --- Executes on button press in ADD_NOISE.
function ADD_NOISE_Callback(hObject, eventdata, handles)
% hObject    handle to ADD_NOISE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global CANVAS_IMAGE

global BACK_INT
global BACK_RANDOM_P
STR = get(handles.BACK_INT,'String'); BACK_INT = str2num(STR);
STR = get(handles.BACK_RANDOM_P,'String'); BACK_RANDOM_P = str2num(STR);

[h v] = size(CANVAS_IMAGE);
CANVAS_RANDOM_BACK = zeros(h,v);

xmin= -BACK_INT;
xmax= BACK_INT;
n=h*v;
r = xmin + (xmax - xmin)*sum(rand(n,BACK_RANDOM_P),2)/BACK_RANDOM_P;
CANVAS_RANDOM_BACK = reshape(r, [h,v]);

global CANVAS_IMAGE_NOISE_ADD
CANVAS_IMAGE_NOISE_ADD = CANVAS_IMAGE + CANVAS_RANDOM_BACK;
    
global FLAG_DISPLAY_NOISE

if(FLAG_DISPLAY_NOISE == 1)
    CANVAS_IMAGE_NOISE_ADD_NORM = uint8(iNORM_2D_ARRAY(CANVAS_IMAGE_NOISE_ADD,1,255));
    CANVAS_IMAGE_NOISE_ADD_NORM_RGB = cat(3, CANVAS_IMAGE_NOISE_ADD_NORM, CANVAS_IMAGE_NOISE_ADD_NORM, CANVAS_IMAGE_NOISE_ADD_NORM);
    imshow((CANVAS_IMAGE_NOISE_ADD_NORM_RGB),'Parent',handles.CANVAS_ADD);
    drawnow
end




% --- Executes on button press in CROP_IMAGE.
function CROP_IMAGE_Callback(hObject, eventdata, handles)
% hObject    handle to CROP_IMAGE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global FOV_HEIGHT
global FOV_WIDTH
global OUT_HEIGHT
global OUT_WIDTH

global CANVAS_HEIGHT
global CANVAS_WIDTH

global BACK_INT
global BACK_RANDOM_P


global CANVAS_IMAGE

global OUT_IMAGE

CROP_CANVAS = zeros(FOV_HEIGHT,FOV_WIDTH);
YS_CANVAS = CANVAS_HEIGHT/2 - round(FOV_HEIGHT/2)+1;
YE_CANVAS = CANVAS_HEIGHT/2 + round(FOV_HEIGHT/2);
XS_CANVAS = CANVAS_WIDTH/2 - round(FOV_WIDTH/2)+1;
XE_CANVAS = CANVAS_WIDTH/2 + round(FOV_WIDTH/2);
CROP_CANVAS = CANVAS_IMAGE(YS_CANVAS: YE_CANVAS, XS_CANVAS: XE_CANVAS);

OUT_IMAGE = zeros(OUT_HEIGHT,OUT_WIDTH);
YS_OUT_IMAGE = OUT_HEIGHT/2 - round(FOV_HEIGHT/2)+1;
YE_OUT_IMAGE = OUT_HEIGHT/2 + round(FOV_HEIGHT/2);
XS_OUT_IMAGE = OUT_WIDTH/2 - round(FOV_WIDTH/2)+1;
XE_OUT_IMAGE = OUT_WIDTH/2 + round(FOV_WIDTH/2);
OUT_IMAGE(YS_OUT_IMAGE: YE_OUT_IMAGE, XS_OUT_IMAGE: XE_OUT_IMAGE) = CROP_CANVAS; 

xmin= -BACK_INT;
xmax= BACK_INT;

STR = get(handles.BACK_INT,'String'); BACK_INT = str2num(STR);
STR = get(handles.BACK_RANDOM_P,'String'); BACK_RANDOM_P = str2num(STR);

n=OUT_HEIGHT*OUT_WIDTH;
r = xmin + (xmax - xmin)*sum(rand(n,BACK_RANDOM_P),2)/BACK_RANDOM_P;
OUT_IMAGE_RANDOM_BACK = reshape(r, [OUT_HEIGHT,OUT_WIDTH]);

OUT_IMAGE = OUT_IMAGE + OUT_IMAGE_RANDOM_BACK;

global FLAG_DISPLAY_OUT

if(FLAG_DISPLAY_OUT == 1)
    OUT_IMAGE_NORM = uint8(iNORM_2D_ARRAY(OUT_IMAGE,1,255));
    OUT_IMAGE_NORM_RGB = cat(3, OUT_IMAGE_NORM, OUT_IMAGE_NORM, OUT_IMAGE_NORM);
    imshow((OUT_IMAGE_NORM_RGB),'Parent',handles.IMAGE_OUT);
    drawnow
end







function [OUT_IMG] = ROI_MERGE(Y,X,ROI_IMG,ROI_SIZE,CANVAS,Height,Width)

yS = floor(Y - ROI_SIZE/2);
yE = yS + ROI_SIZE-1;
xS = floor(X - ROI_SIZE/2);    
xE = xS + ROI_SIZE-1;

if(yS <= 0)
    yS = 1;
    yE = ROI_SIZE;
end
if(xS <= 0)
    xS = 1;
    xE = ROI_SIZE;
end
if(yE >= Height)
    yS = Height - ROI_SIZE+1;
    yE = Height;
end
if(xE >= Width)
    xS = Width - ROI_SIZE+1;
    xE = Width;
end

[H, W] = size(CANVAS);
if(yS<1) 
    yS = 1; yE = ROI_SIZE;
end
if(xS<1)
    yS = 1; yE = ROI_SIZE;
end
if(yE>H)
    yS = H +1 - ROI_SIZE; yE = H;
end
if(xE>W)
    xS = W + 1 - ROI_SIZE; xE = W;
end

%CANVAS(yS:yE,xS:xE) = ROI_IMG;


Radius = floor(ROI_SIZE/2);
x = -Radius:Radius;
y = -Radius:Radius;
[xc, yc] = meshgrid(x,y);
[theta, rho] = cart2pol(xc,yc);

PSF = fspecial('gaussian',Radius,Radius);
ROI_IMG = edgetaper(ROI_IMG,PSF);
for y=1:ROI_SIZE
    for x=1:ROI_SIZE
        rPosition = rho(x,y);
        if(rPosition < Radius)
            CANVAS(yS+y-1,xS+x-1) = ROI_IMG(y,x);
        end
    end
end

OUT_IMG = CANVAS;



% --- Executes on button press in Generate.
function Generate_Callback(hObject, eventdata, handles)
% hObject    handle to Generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% global tg
% global tg0
% global ti
% global ti0
% 
% global ng
% global ng0
% global ni
% global ni0
% global np
% global ns
% 
% global Angle
% 
% global zf
% global zp
% global k
% global radius
% global ROI_SIZE
% global WAVE_LENGTH
% global PX_RES
% 
% global Efficiency
% global SCAT_AMP
% global NA
% Efficiency = (1/pi)*asin(min(NA/ns,1));
% SCAT_AMP = 1;
% 
% WAIT_STEP = 0;
% 
% global PSF_ROI
% number_thetas = 1000;
% max_angle = asin(NA / ni0);
% d_theta = max_angle / number_thetas;
% 
% R = ((ng - ns)/(ng + ns))^2;
% ks = k * ns;
% alpha = 4*pi*radius^3 * (np^2 - ns^2)/(np^2 + 2*ns^2);
% C = ks^4 / (6*pi) * abs(alpha)^2;
% mu = 1/pi * asin(NA/ni);
% 
% M = 170;
% 
% if gpuDeviceCount("available") >= 1
%     gpu_available = 1;
% else
%     gpu_available = 0;
% end
% thetas = linspace(0, max_angle, number_thetas);
% if gpu_available
%     thetas = gpuArray(thetas);
% end
% thetas = reshape(thetas, number_thetas, 1);
% 
% 
% %% Caluculation of OPD
% ti = 0 -zf + ni*(tg0/ng0 - tg/ng + ti0/ni0 - 0/ns);
% OPD_ref_in = ni*ti - ni0*ti0 + ng*tg - ng0*tg0;
% OPD_ref_out = ni*ti - ni0*ti0 + ng*tg - ng0*tg0;
% OPD_ref = OPD_ref_in + OPD_ref_out;
% 
% ti = zp - zf + ni*(tg0/ng0 - tg/ng + ti0/ni0 - zp/ns);
% OPD_scat_in = ns*zp +ni*ti - ni0*ti0 + ng*tg - ng0*tg0;
% ns2_ni2sin2 = ( complex(ns^2 - ni^2 * sin(thetas).^2) ).^0.5;
% ni2_ni2sin2 = ( complex(ni^2 - ni^2 * sin(thetas).^2) ).^0.5;
% ni02_ni2sin2 = ( complex(ni0^2 - ni^2 * sin(thetas).^2) ).^0.5;
% ng2_ni2sin2 = ( complex(ng^2 - ni^2 * sin(thetas).^2) ).^0.5;
% ng02_ni2sin2 = ( complex(ng0^2 - ni^2 * sin(thetas).^2) ).^0.5;
% OPD_scat_out = zp*ns2_ni2sin2 + ti*ni2_ni2sin2 - ti0*ni02_ni2sin2 + tg*ng2_ni2sin2 - tg0*ng02_ni2sin2;
% OPD_scat = OPD_scat_in + OPD_scat_out;
% %%
% 
% %% Calculation of Fresnel reflection and transmission corfficients
% n1 = ng;
% n2 = ns;
% n1cosi = n1 * cos(thetas);
% n2cosi = n2 * cos(thetas);
% n1_n2sini = ( complex(1 - (n1/n2 * sin(thetas)).^2 ) ).^0.5;
% rs = (n1cosi - n2*n1_n2sini) ./ (n1cosi + n2*n1_n2sini);
% rp = (n2cosi - n1*n1_n2sini) ./ (n1*n1_n2sini + n2cosi);
% ts = 2*n2*n1_n2sini ./ (n1cosi + n2*n1_n2sini);
% tp = 2*n2*n1_n2sini ./ (n2cosi + n1*n1_n2sini);
% %%
% 
% %%
% cam_x = linspace(-(ROI_SIZE-1)/4, (ROI_SIZE-1)/4, ROI_SIZE);
% cam_y = linspace(-(ROI_SIZE-1)/4, (ROI_SIZE-1)/4, ROI_SIZE);
% [cam_x, cam_y] = meshgrid(cam_x, cam_y);
% pixel_res = 0.075;
% 
% r_d = pixel_res*M*(cam_x.^2 + cam_y.^2).^0.5;
% % phi_d = angle(pixel_res*M*(cam_x + 1i*cam_y));
% phi_d = atan(cam_y ./ cam_x);
% phi_d(:, (ROI_SIZE+1)/2) = 0;
% if gpu_available
%     r_d = gpuArray(r_d);
%     phi_d = gpuArray(phi_d);
% end
% r_d = reshape(r_d, 1, ROI_SIZE^2);
% phi_d = reshape(phi_d, 1, ROI_SIZE^2);
% sin2pd = sin(2*phi_d);
% cos2pd = cos(2*phi_d);
% sin_td = (ni/M) * sin(thetas);
% cos_td = (complex(1 - sin_td.^2)).^0.5;
% %%
% 
% %%
% bessel0 = besselj(0, k * sin_td * r_d);
% bessel2 = besselj(2, k * sin_td * r_d);
% %%
% 
% %%
% % C_scat = mu * sqrt(C) * sqrt(T) * exp(1i*angle(alpha));
% % C_ref = mu * sqrt(R);
% C_scat = 1;
% C_ref = 1;
% 
% prefactor = k/4 * (ni/M)^2;
% prefactor_scat = prefactor * C_scat;
% prefactor_ref = prefactor * C_ref;
% rps = rp + rs;
% 
% mu = 1;
% mu1 = 1;
% S1 = 0;
% S2 = 0;
% for n = 1 : 3
%     nn = (2*n+1)/(n*(n+1));
%     x = 2 * pi * ns * radius / WAVE_LENGTH;
%     m = np / ns;
%     Vn = legendre(n, cos(thetas)); % n
%     pn = Vn(2,:)';
%     tn = n*cos(thetas).*pn;
%     jn_x = (pi/(2*x))^0.5 * besselj(n+0.5, x);
%     jn1_x = (pi/(2*x))^0.5 * besselj(n+1.5, x);
%     jn_mx = (pi/(2*m*x))^0.5 * besselj(n+0.5, m*x);
%     jn1_mx = (pi/(2*m*x))^0.5 * besselj(n+1.5, m*x);
%     xjn_x_df = jn_x + x*(n/x*jn_x - jn1_x);
%     mxjn_mx_df = jn_mx + m*x*(n/(m*x)*jn_mx - jn1_mx);
%     hn_x = (pi/(2*x))^0.5 * besselh(n+0.5, x);
%     hn1_x = (pi/(2*x))^0.5 * besselh(n+1.5, x);
%     xhn_x_df = hn_x + x*(n/x*hn_x - hn1_x);
% 
%     an = (m^2 * jn_mx * xjn_x_df - jn_x * mxjn_mx_df)...
%         /(m^2 * jn_mx * xhn_x_df - hn_x * mxjn_mx_df);
%     bn = (jn_mx * xjn_x_df - jn_x * mxjn_mx_df)...
%         /(jn_mx * xhn_x_df - hn_x * mxjn_mx_df);
% 
%     S1 = S1 + 1i * nn * (an.*pn + bn.*tn);
%     S2 = S2 + 1i * nn * (an.*tn + bn.*pn);
% end
% 
% E_ref_x = rps .* (1-cos_td) .* bessel2 .* sin2pd .* ...
%     exp(1i*k*(OPD_ref + zf.*cos_td)) .* ...
%     sqrt(ni.*cos_td./cos(thetas)) .* sin(2*thetas);
% E_ref_y = -rps .* ( (1+cos_td).*bessel0 + (1-cos_td).*bessel2.*cos2pd ).*...
%     exp(1i*k*(OPD_ref + zf.*cos_td)) .* ...
%     sqrt(ni.*cos_td./cos(thetas)) .* sin(2*thetas);
% E_scat_x = (-1i) .* ( (S1.*ts + S2.*tp.*cos_td ).*bessel0...
%     - (S2.*tp.*cos_td - S1.*ts  ).*bessel2.*exp(1i.*2.*phi_d) ).*...
%     exp(1i*k*(OPD_scat + zf.*cos_td)) .* ...
%     sqrt(ni.*cos_td./cos(thetas)) .* sin(2.*thetas);
% E_scat_y = ( (S1.*ts + S2.*tp.*cos_td).*bessel0...
%     + (S2.*tp.*cos_td - S1.*ts).*bessel2.*exp(1i.*2.*phi_d) ).*...
%     exp(1i*k*(OPD_scat + zf.*cos_td)) .* ...
%     sqrt(ni.*cos_td./cos(thetas)) .* sin(2*thetas);
% 
% E_scat_x = sum(E_scat_x, 1) * d_theta;
% E_scat_y = sum(E_scat_y, 1) * d_theta;
% E_ref_x = sum(E_ref_x, 1) * d_theta;
% E_ref_y = sum(E_ref_y, 1) * d_theta;
% 
% E_scat_x = reshape(E_scat_x, ROI_SIZE, ROI_SIZE);
% E_scat_y = reshape(E_scat_y, ROI_SIZE, ROI_SIZE);
% E_ref_x = reshape(E_ref_x, ROI_SIZE, ROI_SIZE);
% E_ref_y = reshape(E_ref_y, ROI_SIZE, ROI_SIZE);
% 
% %% Calculation of electric fields
% E_scat = cat(3, E_scat_x, E_scat_y);
% E_ref = cat(3, E_ref_x, E_ref_y);
% 
% E_scat = SCAT_AMP * prefactor_scat * sum(E_scat, 3);
% E_ref = Efficiency*prefactor_ref * sum(E_ref, 3);
% 
% PSF_SCAT = real(sum((E_scat'.*E_scat), 3));
% PSF_iSCAT = real(sum((E_ref.*E_scat' + E_ref'.*E_scat), 3));
% 
% PSF_ROI = gather(PSF_SCAT + PSF_iSCAT)*1e12;
% 
% global FLAG_DISPLAY_PSF
% 
% if(FLAG_DISPLAY_PSF == 1)
%     imagesc(handles.IMAGE_PSF,PSF_ROI);
%     drawnow
% end

% ========================================================================
% ========================================================================

global tg
global tg0
global ti
global ti0

global ng
global ng0
global ni
global ni0
global np
global ns

global Angle

global zf
global zp
global k
global radius
global ROI_SIZE
global WAVE_LENGTH
global PX_RES

global Efficiency
global NA
Efficiency = (1/pi)*asin(min(NA/ns,1));

global SCAT_AMP


STR = get(handles.radius,'String'); radius = str2num(STR);
STR = get(handles.zf,'String'); zf = str2num(STR);
STR = get(handles.zp,'String'); zp = str2num(STR);


M = 250;
WAIT_STEP = 0;

global PSF_ROI
PSF_ROI = zeros(ROI_SIZE,ROI_SIZE);

% for xp = -(ROI_SIZE/2)  + 0.5 : round(ROI_SIZE/2)-1
%     for yp = -(ROI_SIZE/2)  + 0.5 : round(ROI_SIZE/2)-1      % Variable : X, Y Image coordinates
%         WAIT_STEP = WAIT_STEP + 1;
%         WaitPercents = round((WAIT_STEP / (ROI_SIZE)^2)*100);
%         %set(handles.WaitPercents,'String',spZaC rintf('%3d',WaitPercents));
%         %drawnow
%         
%         xCoord = double(xp) * (PX_RES) * M;
%         yCoord = double(yp) * (PX_RES) * M;
% 
%         E_REF_INTEGRATION = 0;
%         E_SCAT_INTEGRATION = 0;
%         E_iSCAT_INTEGRATION = 0;
% 
%         % ==================================================================================
%         % Measurement loop
%         for Theta = 0.0 : 0.01 : Angle   % Variable : Numerical apertuuare: angle [Rad.]  
%         % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%         % START of the integraion loop for the NA of objectives
%         % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%             % =========================================================================================================================================
%             % Observation point: rD, phiD 
%             rD = sqrt(xCoord.^2+yCoord.^2);
%             if(xCoord == 0)
%                 YoverX = 0; % y/x = y/0 = 0;
%             else
%                 YoverX = (yCoord/xCoord);
%             end
%             phiD = atan(YoverX);
%             % =========================================================================================================================================                    
% 
%             % =========================================================================================================================================
%             % Calculate the reference field
%             % =========================================================================================================================================                    
%             n1 = ng;
%             n2 = ns;
%             rs = ( n1*cos(Theta) - n2*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) ) / (n1*cos(Theta) + n2*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) ); 
%             rp = ( n2*cos(Theta) - n1*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) ) / (n2*cos(Theta) + n1*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) );
% 
%             % 'zp' substitutes '0' 
%             ti = 0 - zf + ni*(-0/ns - tg/ng + tg0/ng0 + ti0/ni0);
%             OPD_REF_IN = (ng*tg + ni*ti) - (ng0*tg0 + ni0*ti0) ;
%             OPD_REF_OUT = (ng*tg + ni*ti) - (ng0*tg0 + ni0*ti0);
%             OPD_REF = OPD_REF_IN + OPD_REF_OUT;
% 
%             sin_td = (ni/M)*sin(Theta);
%             cos_td = sqrt(1 - ( (ni/M)*sin(Theta))^2 );
% 
%             E_REF_X = (-1i*k/2)*(-1i*(rp-rs))*...
%                 ( (-1)*(1-cos_td)*besselj(2,k*rD*sin_td)*sin(2*phiD) )*...
%                 exp(1i*k*OPD_REF)*...
%                 exp(1i*k*zf*cos_td)*...
%                 sqrt(ni*cos_td/cos(Theta))*...
%                 (ni/M)^2*(1/2)*sin(2*Theta);
%             E_REF_Y = (-1i*k/2)*(-1i*(rp-rs))*...
%                 ( (1+cos_td)*besselj(0,k*rD*sin_td)+...
%                   (1-cos_td)*besselj(2,k*rD*sin_td)*cos(2*phiD) )*...
%                 exp(1i*k*OPD_REF)*...
%                 exp(1i*k*zf*cos_td)*...
%                 sqrt(ni*cos_td/cos(Theta))*...
%                 (ni/M)^2*(1/2)*sin(2*Theta);
%             
%             E_REF = Efficiency.*(E_REF_X + E_REF_Y);
% 
%             % =========================================================================================================================================
%             % Calculate the scattering field
%             % =========================================================================================================================================
%             n1 = ns;
%             n2 = ng;
%             tp = 2*n1*sqrt(1-((n2/n1)*sin(Theta))^2)/ ( n1*cos(Theta) + n2*sqrt(1-((n2/n1)*sin(Theta))^2) );
%             ts = 2*n1*sqrt(1-((n2/n1)*sin(Theta))^2)/ ( n1*sqrt(1-((n2/n1)*sin(Theta))^2) + n2*cos(Theta) );
% 
%             ti = zp - zf + ni*(-zp/ns - tg/ng + tg0/ng0 + ti0/ni0);
%             OPD_SCAT_IN = (zp*ns + ti*ni + tg*ng ) - (ni0*ti0 + ng0*tg0);
%             OPD_SCAT_OUT = sqrt(ns^2-(ni^2)*sin(Theta)^2)*zp... 
%             + sqrt(ni^2-(ni^2)*sin(Theta)^2)  * ti...            
%             - sqrt(ni0^2-(ni^2)*sin(Theta)^2) * ti0...
%             + sqrt(ng^2-(ni^2)*sin(Theta)^2)  * tg...
%             - sqrt(ng0^2-(ni^2)*sin(Theta)^2) * tg0;
%             OPD_SCAT = OPD_SCAT_IN + OPD_SCAT_OUT;      
% 
%             S1 = 0;
%             S2 = 0;
%             ES1 = 0;     
%             ES2 = 0;     
%             for n = 1:3
%                 % Calculate P_(n-1)
%                 Vn_1 = legendre(n-1,cos(Theta)); % n-1
%                 if(n == 1) 
%                     %Pn_1 = Vn_1(1,:);
%                     Pn_1 = 0;  % P_(1-1)(1) = 0; if |m|>L : P_(L)^(m) = 0.
%                 else
%                     Pn_1 = Vn_1(2,:);
%                 end
%                 % Calculate P_(n)
%                 Vn = legendre(n,cos(Theta)); % n
%                 Pn = Vn(2,:);
% 
%                 x = 2*pi*ns*(radius)./WAVE_LENGTH;   % Expected zp = radius
%                 m = np./ns;
% 
%                 an_numerator = x*(n+1)*(m^2-1)*Bessel(n,m*x)*Bessel(n,x) - m*x^2*(m*Bessel(n,m*x)*Bessel(n+1,x) - Bessel(n,x)*Bessel(n+1,m*x));
%                 an_denominator = an_numerator + 1i*( x*(n+1)*(m^2-1)*Bessel(n,m*x)*Neumann(n,x) - m*x^2*( m*Bessel(n,m*x)*Neumann(n+1,x) - Bessel(n+1,m*x)*Neumann(n,x) ) );
%                 an = an_numerator/an_denominator;
% 
%                 bn_numerator = m*Bessel(n+1,m*x)*Bessel(n,x) - Bessel(n+1,x)*Bessel(n,m*x);
%                 bn_denominator = bn_numerator + 1i*( m*Bessel(n+1,m*x)*Neumann(n,x) - Bessel(n,m*x)*Neumann(n+1,x));
%                 bn = bn_numerator/bn_denominator;
% 
%                 S1 = ( (2*n+1)/(n*(n+1)) ) * ( (an + bn*n*cos(Theta))*Pn  - bn*(n+1)*Pn_1 );
%                 S2 = ( (2*n+1)/(n*(n+1)) ) * ( (an*n*cos(Theta) + bn)*Pn  - an*(n+1)*Pn_1 );
%                 ES1 = ES1 + S1; % Accumulation : Es_phi
%                 ES2 = ES2 + S2; % Accumulation : Es_theta
%             end
%             ES1 = (1/(-1i))*ES1;
%             ES2 = (1/(-1i))*ES2;
% 
%             E_SCAT_X = (-1i*k/2)*...
%                 ( (ES1*ts + ES2*tp*cos_td )*besselj(0,k*rD*sin_td)...
%                 - (ES2*tp*cos_td - ES1*ts  )*besselj(2,k*rD*sin_td)*exp(1i*2*phiD) )*...
%                 exp(1i*k*OPD_SCAT)*...    
%                 exp(1i*k*zf*cos_td)*...                        
%                 sqrt(ni*cos_td/cos(Theta))*...
%                 (ni/M)^2*(1/2)*sin(2*Theta);
% 
%             E_SCAT_Y = (-1i*k/2)*(1i)*...
%                 ( (ES1*ts + ES2*tp*cos_td)*besselj(0,k*sin_td)...
%                 + (ES2*tp*cos_td - ES1*ts)*besselj(2,k*sin_td)*exp(1i*2*phiD) )*...
%                 exp(1i*k*OPD_SCAT)*...    
%                 exp(1i*k*zf*cos_td)*...                        
%                 sqrt(ni*cos_td/cos(Theta))*...
%                 (ni/M)^2*(1/2)*sin(2*Theta);
% 
%             E_SCAT = SCAT_AMP.*(E_SCAT_X + E_SCAT_Y);
%             % =========================================================================================================================================     
%             % Accumulation for reference fields
%             E_REF_INTEGRATION = E_REF_INTEGRATION + E_REF;
%             % Accumulation for scattering fields
%             E_SCAT_INTEGRATION = E_SCAT_INTEGRATION + E_SCAT;
%             % Accumulation for iSCAT fields
%             E_iSCAT_INTEGRATION = E_iSCAT_INTEGRATION + E_REF*E_SCAT' + E_REF'*E_SCAT;
%         end
%         % =============================================================================================================================================
%         PSF_REF(yp+round(ROI_SIZE/2),xp+round(ROI_SIZE/2)) = real(E_REF_INTEGRATION.*E_REF_INTEGRATION');
%         PSF_SCAT(yp+round(ROI_SIZE/2),xp+round(ROI_SIZE/2)) = real(E_SCAT_INTEGRATION.*E_SCAT_INTEGRATION');
%         PSF_iSCAT(yp+round(ROI_SIZE/2),xp+round(ROI_SIZE/2)) = real(E_iSCAT_INTEGRATION);
% 
%         PSF_ROI(yp+round(ROI_SIZE/2),xp+round(ROI_SIZE/2)) = PSF_ROI(yp+round(ROI_SIZE/2),xp+round(ROI_SIZE/2)) + real(E_iSCAT_INTEGRATION);
%     end
% end

PSF_SIZE_X = ROI_SIZE;
PSF_SIZE_Y = ROI_SIZE;
PX_SIZE = PX_RES;
SCAT_COEFF = 0.1;

ixp = -(PSF_SIZE_X/2)  + 0.5 : round(PSF_SIZE_X/2)-1;
iyp = -(PSF_SIZE_Y/2)  + 0.5 : round(PSF_SIZE_Y/2)-1;

[y, x] = meshgrid( -(PSF_SIZE_Y/2)  + 0.5 : round(PSF_SIZE_Y/2)-1,  -(PSF_SIZE_X/2)  + 0.5 : round(PSF_SIZE_X/2)-1);
resultArray = [y(:), x(:), zeros(PSF_SIZE_Y*PSF_SIZE_X,1)];
xp =resultArray(:,2);
yp =resultArray(:,1);
    

xCoord = double(xp) * PX_SIZE *M;
yCoord = double(yp) * PX_SIZE *M;

E_REF_INTEGRATION = 0;
E_SCAT_INTEGRATION = 0;
E_iSCAT_INTEGRATION = 0;
       
dAngle = 0.005;

% ==================================================================================
% Measurement loop
for Theta = 0.0 : dAngle : Angle   % Variable : Numerical apertuuare: angle [Rad.]  
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% START of the integraion loop for the NA of objectives
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    % =========================================================================================================================================
    % Observation point: rD, phiD 
    rD = sqrt(xCoord.^2+yCoord.^2);
    mask_tmp = xCoord == 0; 
    YoverX = (yCoord./xCoord);
    YoverX(mask_tmp) = 0;
    phiD = atan(YoverX);

    zd = 0;
    % =========================================================================================================================================                    
    
    % =========================================================================================================================================
    % Calculate the reference field
    % =========================================================================================================================================                    
    n1 = ng;
    n2 = ns;
    rs = ( n1*cos(Theta) - n2*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) ) / (n1*cos(Theta) + n2*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) ); 
    rp = ( n2*cos(Theta) - n1*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) ) / (n2*cos(Theta) + n1*sqrt(1 - ((n1/n2)*sin(Theta))^2 ) );

    % 'zp' substitutes '0' 
    ti = 0 - zf + ni*(-0/ns - tg/ng + tg0/ng0 + ti0/ni0);
    OPD_REF_IN = (ng*tg + ni*ti) - (ng0*tg0 + ni0*ti0) ;
    OPD_REF_OUT = (ng*tg + ni*ti) - (ng0*tg0 + ni0*ti0);
    OPD_REF = OPD_REF_IN + OPD_REF_OUT;

    
    sin_td = (ni/M)*sin(Theta);
    cos_td = sqrt(1 - ( (ni/M)*sin(Theta))^2 );
     
        E_REF_X = (-1i*k/2)*(-1i*(rp-rs))*...
            ( (-1)*(1-cos_td)*besselj(2,k.*rD*sin_td).*sin(2*phiD) )*...
            exp(1i*k*OPD_REF)*...
            exp(1i*k*zd*cos_td)*...
            sqrt(ni*cos_td/cos(Theta))*...
            (ni/M)^2*(1/2)*sin(2*Theta);

        E_REF_Y = (-1i*k/2)*(-1i*(rp-rs))*...
            ( (1+cos_td).*besselj(0,k.*rD*sin_td)+...
              (1-cos_td).*besselj(2,k*rD*sin_td).*cos(2*phiD) )*...
            exp(1i*k*OPD_REF)*...
            exp(1i*k*zd*cos_td)*...
            sqrt(ni*cos_td/cos(Theta))*...
            (ni/M)^2*(1/2)*sin(2*Theta);
    
    E_REF = E_REF_X + E_REF_Y;
    
    % =========================================================================================================================================
    % Calculate the scattering field
    % =========================================================================================================================================
    n1 = ns;
    n2 = ng;
    tp = 2*n1*sqrt(1-((n2/n1)*sin(Theta))^2)/ ( n1*cos(Theta) + n2*sqrt(1-((n2/n1)*sin(Theta))^2) );
    ts = 2*n1*sqrt(1-((n2/n1)*sin(Theta))^2)/ ( n1*sqrt(1-((n2/n1)*sin(Theta))^2) + n2*cos(Theta) );

    ti = zp - zf + ni*(-zp/ns - tg/ng + tg0/ng0 + ti0/ni0);
    OPD_SCAT_IN = (zp*ns + ti*ni + tg*ng ) - (ni0*ti0 + ng0*tg0);
    OPD_SCAT_OUT = sqrt(ns^2-(ni^2)*sin(Theta)^2)*zp... 
    + sqrt(ni^2-(ni^2)*sin(Theta)^2)  * ti...            
    - sqrt(ni0^2-(ni^2)*sin(Theta)^2) * ti0...
    + sqrt(ng^2-(ni^2)*sin(Theta)^2)  * tg...
    - sqrt(ng0^2-(ni^2)*sin(Theta)^2) * tg0;
    OPD_SCAT = OPD_SCAT_IN + OPD_SCAT_OUT;      

    S1 = 0;
    S2 = 0;
    ES1 = 0;     
    ES2 = 0;     
    for n = 1:3
        % Calculate P_(n-1)
        Vn_1 = legendre(n-1,cos(Theta)); % n-1
        if(n == 1) 
            %Pn_1 = Vn_1(1,:);
            Pn_1 = 0;  % P_(1-1)(1) = 0; if |m|>L : P_(L)^(m) = 0.
        else
            Pn_1 = Vn_1(2,:);
        end
        % Calculate P_(n)
        Vn = legendre(n,cos(Theta)); % n
        Pn = Vn(2,:);

        x = 2*pi*ns*(radius)./WAVE_LENGTH; % Expected zp = radius
        m = np./ns;

        an_numerator = x*(n+1)*(m^2-1)*Bessel(n,m*x)*Bessel(n,x) - m*x^2*(m*Bessel(n,m*x)*Bessel(n+1,x) - Bessel(n,x)*Bessel(n+1,m*x));
        an_denominator = an_numerator + 1i*( x*(n+1)*(m^2-1)*Bessel(n,m*x)*Neumann(n,x) - m*x^2*( m*Bessel(n,m*x)*Neumann(n+1,x) - Bessel(n+1,m*x)*Neumann(n,x) ) );
        an = an_numerator/an_denominator;

        bn_numerator = m*Bessel(n+1,m*x)*Bessel(n,x) - Bessel(n+1,x)*Bessel(n,m*x);
        bn_denominator = bn_numerator + 1i*( m*Bessel(n+1,m*x)*Neumann(n,x) - Bessel(n,m*x)*Neumann(n+1,x));
        bn = bn_numerator/bn_denominator;

        S1 = ( (2*n+1)/(n*(n+1)) ) * ( (an + bn*n*cos(Theta))*Pn  - bn*(n+1)*Pn_1 );
        S2 = ( (2*n+1)/(n*(n+1)) ) * ( (an*n*cos(Theta) + bn)*Pn  - an*(n+1)*Pn_1 );
        ES1 = ES1 + S1; % Accumulation : Es_phi
        ES2 = ES2 + S2; % Accumulation : Es_theta
    end
    ES1 = (1/(-1i))*ES1;
    ES2 = (1/(-1i))*ES2;
    
        E_SCAT_X = (-1i*k/2)*...
            ( (ES1*ts + ES2*tp*cos_td ).*besselj(0,k*rD*sin_td)...
            - (ES2*tp*cos_td - ES1*ts  ).*besselj(2,k*rD*sin_td).*exp(1i*2*phiD) )*...
            exp(1i*k*OPD_SCAT)*...    
            exp(1i*k*zd*cos_td)*...                        
            sqrt(ni*cos_td/cos(Theta))*...
            (ni/M)^2*(1/2)*sin(2*Theta);

        E_SCAT_Y = (-1i*k/2)*(1i)*...
            ( (ES1*ts + ES2*tp*cos_td)*besselj(0,k*sin_td)...
            + (ES2*tp*cos_td - ES1*ts)*besselj(2,k*sin_td).*exp(1i*2*phiD) )*...
            exp(1i*k*OPD_SCAT)*...    
            exp(1i*k*zd*cos_td)*...                        
            sqrt(ni*cos_td/cos(Theta))*...
            (ni/M)^2*(1/2)*sin(2*Theta);
                 
    E_SCAT =  SCAT_COEFF*(E_SCAT_X + E_SCAT_Y);
    
    % =========================================================================================================================================     
    % Accumulation for reference fields
    E_REF_INTEGRATION = E_REF_INTEGRATION + E_REF;
    % Accumulation for scattering fields
    E_SCAT_INTEGRATION = E_SCAT_INTEGRATION + E_SCAT;
    % Accumulation for iSCAT fields
   E_iSCAT_INTEGRATION = E_iSCAT_INTEGRATION + E_REF.*transpose(E_SCAT') + transpose(E_REF').*E_SCAT;
end

PSF_E_REF_TEMP = reshape(E_REF_INTEGRATION, PSF_SIZE_Y, PSF_SIZE_X);
PSF_E_SCAT_TEMP = reshape(E_SCAT_INTEGRATION, PSF_SIZE_Y, PSF_SIZE_X);


PSF_E_REF = PSF_E_REF_TEMP;
PSF_E_SCAT = PSF_E_SCAT_TEMP;


PSF_IMG_REF = real(PSF_E_REF.*PSF_E_REF');
PSF_IMG_SCAT = real(PSF_E_SCAT.*PSF_E_SCAT');
E_iSCAT_INTEGRATION_2d = reshape(E_iSCAT_INTEGRATION, PSF_SIZE_Y, PSF_SIZE_X);
PSF_IMG_iSCAT = real(E_iSCAT_INTEGRATION_2d);

PSF_ROI = (PSF_IMG_SCAT + PSF_IMG_iSCAT)*SCAT_AMP;

PSF_ROI = imgaussfilt(PSF_ROI,1);

global FLAG_DISPLAY_PSF

if(FLAG_DISPLAY_PSF == 1)
    imagesc(handles.IMAGE_PSF,PSF_ROI);
    drawnow
end
% ========================================================================
% ========================================================================


function NP_Callback(hObject, eventdata, handles)
% hObject    handle to NP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NP as text
%        str2double(get(hObject,'String')) returns contents of NP as a double


% --- Executes during object creation, after setting all properties.
function NP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NumFrame_Callback(hObject, eventdata, handles)
% hObject    handle to NumFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumFrame as text
%        str2double(get(hObject,'String')) returns contents of NumFrame as a double


% --- Executes during object creation, after setting all properties.
function NumFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BACK_INT_Callback(hObject, eventdata, handles)
% hObject    handle to BACK_INT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BACK_INT as text
%        str2double(get(hObject,'String')) returns contents of BACK_INT as a double


% --- Executes during object creation, after setting all properties.
function BACK_INT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BACK_INT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BACK_RANDOM_P_Callback(hObject, eventdata, handles)
% hObject    handle to BACK_RANDOM_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BACK_RANDOM_P as text
%        str2double(get(hObject,'String')) returns contents of BACK_RANDOM_P as a double


% --- Executes during object creation, after setting all properties.
function BACK_RANDOM_P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BACK_RANDOM_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OUT_WIDTH_Callback(hObject, eventdata, handles)
% hObject    handle to OUT_WIDTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OUT_WIDTH as text
%        str2double(get(hObject,'String')) returns contents of OUT_WIDTH as a double


% --- Executes during object creation, after setting all properties.
function OUT_WIDTH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OUT_WIDTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OUT_HEIGHT_Callback(hObject, eventdata, handles)
% hObject    handle to OUT_HEIGHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OUT_HEIGHT as text
%        str2double(get(hObject,'String')) returns contents of OUT_HEIGHT as a double


% --- Executes during object creation, after setting all properties.
function OUT_HEIGHT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OUT_HEIGHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FOV_WIDTH_Callback(hObject, eventdata, handles)
% hObject    handle to FOV_WIDTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FOV_WIDTH as text
%        str2double(get(hObject,'String')) returns contents of FOV_WIDTH as a double


% --- Executes during object creation, after setting all properties.
function FOV_WIDTH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FOV_WIDTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FOV_HEIGHT_Callback(hObject, eventdata, handles)
% hObject    handle to FOV_HEIGHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FOV_HEIGHT as text
%        str2double(get(hObject,'String')) returns contents of FOV_HEIGHT as a double


% --- Executes during object creation, after setting all properties.
function FOV_HEIGHT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FOV_HEIGHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CHECK_DISPLAY_PSF.
function CHECK_DISPLAY_PSF_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_DISPLAY_PSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_DISPLAY_PSF
global FLAG_DISPLAY_PSF

Value = get(handles.CHECK_DISPLAY_PSF,'Value');
if(Value == 1)
    FLAG_DISPLAY_PSF = 1;
elseif(Value == 0)
    FLAG_DISPLAY_PSF = 0;
end

% --- Executes on button press in CHECK_DISPLAY_CANVAS.
function CHECK_DISPLAY_CANVAS_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_DISPLAY_CANVAS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_DISPLAY_CANVAS
global FLAG_DISPLAY_CANVAS
Value = get(handles.CHECK_DISPLAY_CANVAS,'Value');
if(Value == 1)
    FLAG_DISPLAY_CANVAS = 1;
elseif(Value == 0)
    FLAG_DISPLAY_CANVAS = 0;
end

% --- Executes on button press in CHECK_DISPLAY_NOISE.
function CHECK_DISPLAY_NOISE_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_DISPLAY_NOISE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_DISPLAY_NOISE
global FLAG_DISPLAY_NOISE
Value = get(handles.CHECK_DISPLAY_NOISE,'Value');
if(Value == 1)
    FLAG_DISPLAY_NOISE = 1;
elseif(Value == 0)
    FLAG_DISPLAY_NOISE = 0;
end

% --- Executes on button press in CHECK_DISPLAY_OUT.
function CHECK_DISPLAY_OUT_Callback(hObject, eventdata, handles)
% hObject    handle to CHECK_DISPLAY_OUT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHECK_DISPLAY_OUT
global FLAG_DISPLAY_OUT
Value = get(handles.CHECK_DISPLAY_OUT,'Value');
if(Value == 1)
    FLAG_DISPLAY_OUT = 1;
elseif(Value == 0)
    FLAG_DISPLAY_OUT = 0;
end



function TICK_TIME_Callback(hObject, eventdata, handles)
% hObject    handle to TICK_TIME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TICK_TIME as text
%        str2double(get(hObject,'String')) returns contents of TICK_TIME as a double


% --- Executes during object creation, after setting all properties.
function TICK_TIME_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TICK_TIME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EXPOSURE_TIME_Callback(hObject, eventdata, handles)
% hObject    handle to EXPOSURE_TIME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EXPOSURE_TIME as text
%        str2double(get(hObject,'String')) returns contents of EXPOSURE_TIME as a double


% --- Executes during object creation, after setting all properties.
function EXPOSURE_TIME_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EXPOSURE_TIME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zp_Callback(hObject, eventdata, handles)
% hObject    handle to zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zp as text
%        str2double(get(hObject,'String')) returns contents of zp as a double


% --- Executes during object creation, after setting all properties.
function zp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zf_Callback(hObject, eventdata, handles)
% hObject    handle to zf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zf as text
%        str2double(get(hObject,'String')) returns contents of zf as a double


% --- Executes during object creation, after setting all properties.
function zf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radius_Callback(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius as text
%        str2double(get(hObject,'String')) returns contents of radius as a double


% --- Executes during object creation, after setting all properties.
function radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VISCOSITY_Callback(hObject, eventdata, handles)
% hObject    handle to VISCOSITY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VISCOSITY as text
%        str2double(get(hObject,'String')) returns contents of VISCOSITY as a double


% --- Executes during object creation, after setting all properties.
function VISCOSITY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VISCOSITY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in WITH_DAT_FILE.
function WITH_DAT_FILE_Callback(hObject, eventdata, handles)
% hObject    handle to WITH_DAT_FILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WITH_DAT_FILE
global WITH_DAT_FILE

Value = get(handles.WITH_DAT_FILE,'Value');
if(Value == 1)
    WITH_DAT_FILE = 1;
elseif(Value == 0)
    WITH_DAT_FILE = 0;
end



function FrameIdx_Callback(hObject, eventdata, handles)
% hObject    handle to FrameIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameIdx as text
%        str2double(get(hObject,'String')) returns contents of FrameIdx as a double


% --- Executes during object creation, after setting all properties.
function FrameIdx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
