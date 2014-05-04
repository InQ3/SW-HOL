function varargout = qualityComparisonDialog(varargin)
% QUALITYCOMPARISONDIALOG M-file for qualityComparisonDialog.fig
%      QUALITYCOMPARISONDIALOG, by itself, creates a new QUALITYCOMPARISONDIALOG or raises the existing
%      singleton*.
%
%      H = QUALITYCOMPARISONDIALOG returns the handle to a new QUALITYCOMPARISONDIALOG or the handle to
%      the existing singleton*.
%
%      QUALITYCOMPARISONDIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUALITYCOMPARISONDIALOG.M with the given input arguments.
%
%      QUALITYCOMPARISONDIALOG('Property','Value',...) creates a new QUALITYCOMPARISONDIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before qualityComparisonDialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to qualityComparisonDialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qualityComparisonDialog

% Last Modified by GUIDE v2.5 24-May-2012 22:15:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qualityComparisonDialog_OpeningFcn, ...
                   'gui_OutputFcn',  @qualityComparisonDialog_OutputFcn, ...
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


% --- Executes just before qualityComparisonDialog is made visible.
function qualityComparisonDialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to qualityComparisonDialog (see VARARGIN)

% Choose default command line output for qualityComparisonDialog
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes qualityComparisonDialog wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = qualityComparisonDialog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
%Execute selected refinement algorithms and build evolution graph
val = get(handles.popupmenu1,'Value');
switch val
case 1
%Load regular Tetrahedron from mat file
load t_regular;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end


title('Mesh Quality Evolution','FontSize',12); %Title of graph


case 2
%Load Cap Tetrahedron from mat file
load capTet;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph
    
case 3
    
%Load Liu Joe Tetrahedron from mat file    
load liujoeTet;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph
    
case 4
    
%Load Needle Tetrahedron from mat file  
load needleTet;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph
    
case 5
%Load Sliver Tetrahedron from mat file 
load sliverTet;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph
    
case 6
%Load Bey Tetrahedron from mat file 
load tbey;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph
    
case 7
%Load Wedge Tetrahedron from mat file 
load wedgeTet;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph
    
case 8
%Load Cube from mat file   
load cubo;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph
    
case 9
%Load Sample Tetrahedral Mesh from mat file    
load tetmesh;
% Update global data structure from file loaded
TetraDT=Tes;
TetraCoordinates=X;

%Get iteration steps

iterationNumber = str2double(get(handles.edit2,'string'));
if  (iterationNumber<=1)
  errordlg('You must enter a valid iteration number','Bad Input','modal')
  return 
end


%Get Refinement Algorithm
firstAlgorithm = get(handles.popupmenu2,'Value');
secondAlgorithm =get(handles.popupmenu3,'Value');

if (firstAlgorithm==secondAlgorithm)
   errordlg('Must select different algorithms','modal');
   return
end

%Get Refinement Type Option
refineType= get(handles.popupmenu5,'Value');

%Get Quality Measure Type
qMeasure=get(handles.popupmenu4,'Value');
 
switch firstAlgorithm

    case 1
       meanQualityValues=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);            
    case 2
       meanQualityValues=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);  
    case 3 
       meanQualityValues=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 4
       meanQualityValues=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 5
       meanQualityValues=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
       meanQualityValues=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
       meanQualityValues=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
        
end

switch secondAlgorithm

    case 1
       meanQualityValues2=LE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 2 
       meanQualityValues2=ThreeTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure); 
    case 3 
       meanQualityValues2=LETrisection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 4
       meanQualityValues2=LEFourSection(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 5
        meanQualityValues2=FourTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 6
        meanQualityValues2=EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);
    case 7
        meanQualityValues2=FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure);

end

%Build Quality Evolution Graph

figure; %Show figure
hold on; %hold graph on screen
plot(1:iterationNumber,meanQualityValues,1:iterationNumber,meanQualityValues2,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationNumber 0 1 ]); %Axis scaling


switch qMeasure

    case 1
    legend('Quality Formula Etha'); %Legend of graph , quality formula  
    case 2
    legend('Quality Formula Whiteh'); %Legend of graph , quality formula      
    case 3
    legend('Quality Formula Ratio'); %Legend of graph , quality formula        
    case 4
    legend('Quality Formula Solid Angle'); %Legend of graph , quality formula        
end

title('Mesh Quality Evolution','FontSize',12); %Title of graph

end    



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function text1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
