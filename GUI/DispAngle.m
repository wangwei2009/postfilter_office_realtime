function varargout = DispAngle(varargin)
% DISPANGLE MATLAB code for DispAngle.fig
%      DISPANGLE, by itself, creates a new DISPANGLE or raises the existing
%      singleton*.
%
%      H = DISPANGLE returns the handle to a new DISPANGLE or the handle to
%      the existing singleton*.
%
%      DISPANGLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPANGLE.M with the given input arguments.
%
%      DISPANGLE('Property','Value',...) creates a new DISPANGLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DispAngle_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DispAngle_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DispAngle

% Last Modified by GUIDE v2.5 05-Dec-2017 19:23:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DispAngle_OpeningFcn, ...
                   'gui_OutputFcn',  @DispAngle_OutputFcn, ...
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

% --- Executes just before DispAngle is made visible.
function DispAngle_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DispAngle (see VARARGIN)

% Choose default command line output for DispAngle
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using DispAngle.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

global NS;
NS = 0;
global start;
start = 0;
global recObj;
global dir;
dir = 1;
global fs;
global buffer_size;
buffer_size = 1024;
fs = 16000;
global deviceReader;
global channels;
channels = 6;
deviceReader = audioDeviceReader('NumChannels',channels,'SampleRate',fs,'SamplesPerFrame',buffer_size);
devices = getAudioDevices(deviceReader)
% SoundCardNum = input('please select XMOS sound card number:');
SoundCardNum = 3;
deviceReader = audioDeviceReader('NumChannels',channels,'SampleRate',fs,'SamplesPerFrame',buffer_size,'Device',devices{SoundCardNum});

%devices = getAudioDevices(deviceReader);
setup(deviceReader);
global deviceWriter;
deviceWriter = audioDeviceWriter('SampleRate',fs);
setup(deviceWriter,zeros(buffer_size,1));

% UIWAIT makes DispAngle wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DispAngle_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% axes(handles.axes1);
% cla;

global deviceReader;
global deviceWriter;
global fs;
global NS;
global channels;

global dir;
global start;

if NS == 0
    set(handles.pushbutton3,'string','NS OFF....');
    set(handles.pushbutton3, 'BackgroundColor',[1 0 0]);
else
    set(handles.pushbutton3,'string','NS ON....');
    set(handles.pushbutton3, 'BackgroundColor',[0 1 0]);
end


if start == 0
    start = 1;
    set(handles.pushbutton1,'string','recording....');
    set(handles.pushbutton1, 'BackgroundColor',[0 1 0]);
else
    start = 0;
    
    set(handles.pushbutton1,'string','stop');
    set(handles.pushbutton1, 'BackgroundColor',[1 0 0]);
    
    NS = 0;
    set(handles.pushbutton3,'string','NS OFF....');
    set(handles.pushbutton3, 'BackgroundColor',[1 0 0]);
end

global angle;
angle_dir = (dir-1)*30;
angle = [angle_dir,30]/180*pi;
M = 4;        %Channels
global r;
r = 0.032;
%%

%% Frequency domain delay-sum,time alignment
% [ DelaySumOut, x] = DelaySumURA(x,fs,N_FFT,N_FFT,N_FFT/2,r,angle);


%%
global N_FFT;
N_FFT = 256;

window = hamming(N_FFT);


P_len = N_FFT/2+1;

Pxii_pre = ones(M,P_len);

Pxij_pre = ones((M*M-M)/2,P_len);
%%

inc = 128;
frameLength = 256;
overlap = frameLength - inc;


last_acquiredAudio = zeros(overlap,channels);
y_last_tail = zeros(overlap,1);
% playData = zeros(chunk_size,1);

L = 1;
% w = randn(L,1);
ang = [90;0];
timeCount = 0;
global buffer_size;
noisePeriod = round(0.1*16000/buffer_size);%
global noise;
noise = zeros(buffer_size*noisePeriod,M);

global Fvv;
Fvv = zeros(N_FFT/2+1,M,M);
while start
    acquiredAudio = deviceReader();
    if(timeCount<noisePeriod)
        noise(timeCount*buffer_size+1:timeCount*buffer_size+buffer_size,:) = acquiredAudio(:,[2,3,4,5]);
        timeCount = timeCount+1;
    else
        x = [last_acquiredAudio(:,[2,3,4,5]);acquiredAudio(:,[2,3,4,5])];
        if NS
            %% Frequency domain delay-sum,time alignment
            [ DelaySumOut, x] = DelaySumURA(x,fs,N_FFT,N_FFT,N_FFT/2,r,angle);
        
            %% postfilter
            [y,Pxii_pre,Pxij_pre]= post_filter_func( x,fs,Fvv,Pxii_pre,Pxij_pre,angle);
            threshold = 0.2;
            %% updata noise coherence
            if(sum(abs(x(:,1)))/size(x,1)<threshold)
                t = 0;
                alpha = 0.9;
                for i = 1:M-1
                    for j = i+1:M
                        t = t+1;
                        Fvv_ij = Pxij_pre(t,:)./sqrt(Pxii_pre(i,:).*Pxii_pre(i,:));
                        Fvv(:,i,j) = alpha*Fvv(:,i,j)+(1-alpha)*real(Fvv_ij)';
                        index = find(Fvv(:,i,j)>0.90);
                        if(size(index,1)>0)
                            Fvv(index,i,j)=0.90;
                        end
                        Fvv(1,i,j) = 0.90;
                    end
                end
            end

            y = y';
            
            %% compensate window
            playData = y(1:end-overlap);
            playData(1:overlap) = playData(1:overlap)+y_last_tail;
        else
            y = x(:,1);
            playData = y(1:end-overlap);
        end

        deviceWriter(real(playData));
        
        y_last_tail = y(end-overlap+1:end);

        %% concatenate
        last_acquiredAudio = acquiredAudio(end-overlap+1:end,:);
    end
    

%     set(handles.text1,'string',beta);
    theta = 0:0.01:pi/4;
    rho = sin(2*theta).*cos(2*theta);
%     polarplot(theta-22.5/180*pi+ang(1)/180*pi,rho)
    polarplot(theta-22.5/180*pi+(dir-1)*30/180*pi,rho)
    drawnow
end



% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dir;
if(dir == 12)
    dir = 1;
else
    dir = dir + 1;
end
set(handles.pushbutton2,'string',dir);
global angle;
angle_dir = (dir-1)*30;
angle = [angle_dir,30]/180*pi;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NS;
if NS == 1
    NS = 0;
    set(handles.pushbutton3,'string','NS OFF....');
    set(handles.pushbutton3, 'BackgroundColor',[1 0 0]);
else
    NS = 1;
        set(handles.pushbutton3,'string','NS ON....');
    set(handles.pushbutton3, 'BackgroundColor',[0 1 0]);
end
