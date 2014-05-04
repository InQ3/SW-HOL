function varargout = main(varargin)
% MAIN M-file for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 20-Oct-2012 23:27:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NewMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to NewMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RefinementMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to RefinementMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LEMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to LEMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function TetrahedralSelectionMenutem_Callback(hObject, eventdata, handles)
% hObject    handle to TetrahedralSelectionMenutem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RefineAllMenuItem_Callback(hObject, eventdata, handles)
% Uniform Longest Edge Refinement Algorithm 
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end    
  
    
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID 
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  TetraDT(1:Tetracount,:)=[];

    
 %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
  
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info





% --------------------------------------------------------------------
function ExitMenuItem_Callback(hObject, eventdata, handles)
% Exit Application
exit;


% --------------------------------------------------------------------
function QualityMenuItem_Callback(hObject, eventdata, handles) 
% hObject    handle to QualityMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function TetraQualityMenuItem_Callback(hObject, eventdata, handles)
%Measure of the quality of tetrahedral mesh
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element

quality =zeros(1,Tetracount); %preallocating for improving performance
for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
       
end

figure; %Show figure
hold on; %hold graph on screen
bar(quality,1:Tetracount,'hist'); %plot bar graph as histogram
grid on; % Turn on grid lines for this plot
xlabel('Tetrahedra Quality');  %label of axis x
ylabel('#Tetrahedra');  %label of axis y
axis([0 1 0 Tetracount ]); %Axis scaling 
title('Tetrahedral Mesh Quality'); %Title of graph





 


% --------------------------------------------------------------------
function ViewMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ViewMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ZoomMenuItem_Callback(hObject, eventdata, handles)
   zoom on; % Enable Zoom on figure

% hObject    handle to ZoomMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PanMenuItem_Callback(hObject, eventdata, handles)
  pan on; %Enable Pan on figure
% hObject    handle to PanMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RotationMenuItem_Callback(hObject, eventdata, handles)
 rotate3d on; % Enable Rotate 3D on figure
% hObject    handle to RotationMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HistogramMenuItem_Callback(hObject, eventdata, handles)
%Classify the classes of tets based on the quality formula Etha
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element

quality =zeros(1,Tetracount); %preallocating for improving performance
for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %quality measure Etha
       
end


%Count number of tetrahedra in quality intervals
count1=0;  %Initialization of counter variables
count2=0;
count3=0;
count4=0;
for i=1:Tetracount
              
        if ((quality(i)>=0.0000) && (quality(i)<=0.2000))
         count1 =count1+1;   
        end
        
         if ((quality(i)>0.2000) && (quality(i)<=0.4000))
         count2 =count2+1;   
        end
        
         if ((quality(i)>0.4000) && (quality(i)<=0.6000))
         count3 =count3+1;   
         end
        
         if ((quality(i)>0.6000) && (quality(i)<=1.0000))
         count4 =count4+1;   
        end
       
end

NumberTets =[count1 count2 count3 count4]; %Number of tetrahedra for each interval , Y axis
Intervals = [0.2 0.4 0.6 1]; %Intervals Values , X axis
figure; %Show figure
hold on; %hold graph on screen
bar(Intervals,NumberTets); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('Tetrahedra Quality');  %label of axis x
ylabel('#Tetrahedra');  %label of axis y
axis([0 1 0 Tetracount ]); %Axis scaling 
title('Tetrahedral Mesh Quality'); %Title of graph



% --------------------------------------------------------------------
function MeshSimplices_Callback(hObject, eventdata, handles)
% View Mesh Simplices
global TetraDT; %global variable tetrahedral triangulation
figure; %show figure
cnames = {'V1','V2','V3','V4'}; % table column names
uitable('Data',sortrows(TetraDT(:,:)),'ColumnName',cnames,'Units','normalized'); %create sorted table simplices in figure



% --------------------------------------------------------------------
function VerticesInfo_Callback(hObject, eventdata, handles)
% View Mesh Vertices
global TetraCoordinates; % global variable triangulation points
figure; %show figure
cnames = {'X','Y','Z'}; % table column names
uitable('Data',sortrows(TetraCoordinates(:,:)),'ColumnName',cnames,'Units','normalized'); %create sorted table coordinates  in figure


% --------------------------------------------------------------------
function rtetra_Callback(hObject, eventdata, handles)
%Load regular Tetrahedron from mat file
tic; %start timer to measure performance
load t_regular;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display regular tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinenement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text



% --------------------------------------------------------------------
function ValueMenuItem_Callback(hObject, eventdata, handles)
% Local Longest Edge Refinement Algorithm by Input Value
%Input value by user
answer = inputdlg({'All Tetrahedra LE > Value will be refine:'},'Input Value');

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Value Input','Error Window','error');

elseif(value>=0) 

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
IterationIndex=[]; %init variable
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition if LE Distance > Input Value
 if(x<=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end    
  
    
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Value Input','Error Window','error');    
end    

end
 


% --------------------------------------------------------------------
function VertexMenuItem_Callback(hObject, eventdata, handles)
% Local Longest Edge Refinement Algorithm by Vertex ID
%Input Vertex ID by user
answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
[vertexID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
  
[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertexID>0  & vertexID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
 TV = vertexAttachments(trep,vertexID); %Return tetrahedra indices attached to specified vertex
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
    
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end    
  
    
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
 
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function EdgeMenuItem_Callback(hObject, eventdata, handles)
% Local Longest Edge Refinement Algorithm by Edge
%Input Edge by Vertex1 ID and Vertex2 ID
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
 %get Vertex1 ID   
[vertex1ID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    

 %get Vertex2 ID   
[vertex2ID status] =str2num(answer{2}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end  


[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertex1ID>0  & vertex1ID<=VertexCount & vertex2ID>0 & vertex2ID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
  %Test if Vertices are joined by Edge
  edge=isEdge(trep,vertex1ID,vertex2ID);
  
  if(edge==false)
     % Handle when vertices are not joined by edge
    msgbox('Vertices are not joined by Edge','Error Window','error');
    return;  
  end    
    
 TV = edgeAttachments(trep,vertex1ID,vertex2ID); %Return tetrahedra indices attached to specified edge defined by vertices
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
    
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end    
  
    
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
 
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function ThreeTMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ThreeTMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tetraSelect_Callback(hObject, eventdata, handles)
% hObject    handle to tetraSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function UniformRefinement_Callback(hObject, eventdata, handles)
% 3T-LE Uniform Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
  
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Longest Edge
  
  SortedDistance = sort(Distance,'descend'); %Sort distances in descending order
  SecondLE=SortedDistance(2); %get second longest edge distance
  
  %Check secondary distance value to corresponding edge given by vertexes id
  %Check the case when second longest edge distance is equal to primary
  %longest edge
  
  %When Secondary Longest Edge is equal to Edge1
  if(SecondLE==edge1)&&(j1~=1  ||  j2~=2)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,1),:)+TetraCoordinates(TetraDT(i,2),:))/2;
      j3=1; %save indices of secundary longest edge vertices
      j4=2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
  %When Secondary Longest Edge is equal to Edge2
  elseif(SecondLE==edge2) && (j1~=2  ||  j2~=3)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,2),:)+TetraCoordinates(TetraDT(i,3),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=3;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  %When Secondary Longest Edge is equal to Edge3
  elseif(SecondLE==edge3) && (j1~=3  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,3),:)+TetraCoordinates(TetraDT(i,1),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
      
  %When Secondary Longest Edge is equal to Edge4
  elseif(SecondLE==edge4) && (j1~=2  ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,2),:)+TetraCoordinates(TetraDT(i,4),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
 
   %When Secondary Longest Edge is equal to Edge5
  elseif(SecondLE==edge5) &&( j1~=3 ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,3),:)+TetraCoordinates(TetraDT(i,4),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  
  %When Secondary Longest Edge is equal to Edge6
  elseif(SecondLE==edge6) && (j1~=4  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,4),:)+TetraCoordinates(TetraDT(i,1),:))/2;
      j3=4; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  end %End of distance comparison with edges and computation of secondary midpoint 
      
  %Composing the new three tetrahedra
  % 3T-LE partition is composed of 2 subdivision patterns
  % Pattern 1: Longest Edge share a vertex with secondary longest edge
  % Pattern 2: Longest Edge is opposed to secondary longest edge
  %Applyinng proper pattern according to secondary longest edge position
    
  if(j1==j3 || j1==j4 || j2==j3 || j2==j4) %Pattern 1 is apply 
     
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
   Vertex3ID=TetraDT(i,j3);
   
  %Finding Vertex 4 ID
  VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)& VertexIDs ~= TetraDT(i,j3)& VertexIDs ~= TetraDT(i,j4));
  Vertex4ID = VertexIDs(b(1));
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
  Vertex3ID=TetraDT(i,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  %Finding Vertex 2ID
  VertexIDs =TetraDT(i,:);
   [a b] = find(VertexIDs ~=Vertex3 & VertexIDs ~= Vertex4 & VertexIDs ~= TetraDT(i,j3)& VertexIDs ~= TetraDT(i,j4));
  Vertex2ID = VertexIDs(b(1)); 
  
  
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  else %Pattern 2 is apply
  
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(i,j3);   
  
  %Finding Vertex 4 ID
  %Vertex 4 ID is one longest edge vertex id
  Vertex4ID=TetraDT(i,j2);
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(i,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
   %Finding Vertex 2ID , vertex 2 ID is one longest edge vertex
    Vertex2ID =TetraDT(i,j1);
    
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
    
  end    
  
        
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Longest Edge Vertexes ID for checking neighbor tetrahedra
  SecVertex1ID =TetraDT(i,j3);
  SecVertex2ID =TetraDT(i,j4);
  
  SecondaryEdge =[SecVertex1ID SecVertex2ID]; %Concatenate Secondary Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
   
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  TetraDT(1:Tetracount,:)=[];

    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function LeTrisectionMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to LeTrisectionMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tetMenu_Callback(hObject, eventdata, handles)
% hObject    handle to tetMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function UniformTriSection_Callback(hObject, eventdata, handles)
% Uniform Longest Edge Trisection Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate First Equidistant Point for longest edge trisection
    midP1=TetraCoordinates(TetraDT(i,j1),:)*2/3+TetraCoordinates(TetraDT(i,j2),:)/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge trisection
    midP2=TetraCoordinates(TetraDT(i,j1),:)/3+TetraCoordinates(TetraDT(i,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
 %Composing The New Three Tetrahedra by Trisection of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j2);
  
  Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  
  %Composing New Tetrahedron 3
   Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
     
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  TetraDT(1:Tetracount,:)=[];

    
  %Algorithm Assure-Conformity of the tet mesh
      while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
           SelectedTetraIndex=[]; %init variable
           flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
           
           %Calculate LEPP
           %Sequential Search for finding neighbors tetrahedra set
           [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
              for k=1:Tetracount
                  VertexIDs =TetraDT(k,:);
                  indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
                  indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
                  
                  if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                     SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                     flagHasNeighbor =true;
                  end    
                  
              end   
           
         if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
            SurroundingEdgeSet(1,:)=[]; 
         end 
             
        %Perform Longest Edge Trisection to selected Tetrahedra
        [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
        
    for i=1:Tetcount %iterate over each selected tetrahedron
            x=SelectedTetraIndex(i,1); %get selected tetrahedra index
      %Calculate Edge Length
     edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
     edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
     edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
     edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
     edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
     edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
    
     Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
     
     [z,d] = max(Distance,[],2); %Obtain Max Distance
     
     %Saving Original Edge Order
     % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
     % Edge number:      1    2    3    4    5    6
     
     V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
     
     [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
             
       % Calculate First Equidistant Point for longest edge trisection
        midP1=TetraCoordinates(TetraDT(x,j1),:)*2/3+TetraCoordinates(TetraDT(x,j2),:)/3;
        
      
      % Add new point into TetraCoordinates Vertex Matrix , checking if point
      % is duplicated in data structure
      
      [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
      
      if(isempty(r)==false)
      %First MidPoint ID
      MidPID=r;      
      else
      TetraCoordinates=[TetraCoordinates;midP1];
      %First MidPoint ID
      [row column] =size(TetraCoordinates);
      MidPID=row;
      end    
      
      % Calculate Second Equidistant Point for longest edge trisection
        midP2=TetraCoordinates(TetraDT(x,j1),:)/3+TetraCoordinates(TetraDT(x,j2),:)*2/3;
        
      
      % Add new point into TetraCoordinates Vertex Matrix , checking if point
      % is duplicated in data structure
      
      [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
      
      if(isempty(r)==false)
      %Second MidPoint ID
      SecMidPID=r;      
      else
      TetraCoordinates=[TetraCoordinates;midP2];
      %Second Midpoint ID
      [row column] =size(TetraCoordinates);
      SecMidPID=row;
      end    
      
     %Composing The New Three Tetrahedra by Trisection of Longest Edge
      
      % Composing New Tetrahedron 1
      %Finding Vertex 1 ID
      Vertex1=TetraDT(x,j1);
      
      %Finding Vertex 3 and 4 ID
        VertexIDs =TetraDT(x,:);
        [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
      
        Vertex3 = VertexIDs(b(1));
        Vertex4 = VertexIDs(b(2));
      
      Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
      
      %Updating Data Structure with tet1
      TetraDT =[TetraDT;Tetrahedron1];
      
      %Composing New Tetrahedron 2
      %Finding Vertex 1 ID
      Vertex1=TetraDT(x,j2);
      
      Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
        %Updating Data Structure with tet2
      TetraDT =[TetraDT;Tetrahedron2];
      
      
      %Composing New Tetrahedron 3
       Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
      %Updating Data Structure with tet3
      TetraDT =[TetraDT;Tetrahedron3];
     
        
      %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
      Vertex1ID =TetraDT(x,j1);
      Vertex2ID =TetraDT(x,j2);
      
      Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
      SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
      
      
    end
    
      %iterate over each old selected tetrahedron updating Data Structure deleting old ones
       if (isempty(SelectedTetraIndex)==0)
        TetraDT(SelectedTetraIndex(:),:)=[];
       end
        
              
     end   
  
 

tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths). Etha
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX.Whiteh
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.Ratio
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.Solid Angle
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function centroid_Callback(hObject, eventdata, handles)
% hObject    handle to centroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tetSelection_Callback(hObject, eventdata, handles)
% hObject    handle to tetSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function UniformDelaunayCentroid_Callback(hObject, eventdata, handles)
%Delaunay Uniform Centroid Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
   tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
   centroid=tetrahedron_centroid_3d(tetra); 
   TetraCoordinates=[TetraCoordinates;centroid];
   
end

%Removed Duplicated points if they exist
TetraCoordinates=unique(TetraCoordinates,'rows');

%Apply Delaunay Refinement Algorithm 
dt = DelaunayTri(TetraCoordinates);

%Update Global Data Structure
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;

%Check Delaunay Bug , Colinearity of points produces degeneracy
% [Tetracount vertexNumber]= size(TetraDT); %number of element

% SelectedTetraIndex=[]; %init variable
%for i=1:Tetracount %iterate over each tetrahedron
 %  tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
  % volume=tetrahedron_volume_3d(tetra); %Calculate Volume 
  % if(volume==0) %if tet is degenerate
   %  SelectedTetraIndex =[SelectedTetraIndex;i]; %store tetrahedra index in data structure next to deletion
  % end    
   
%end

%if (isempty(SelectedTetraIndex)==0) %Delete Degenerate Tets by delaunay colinearity
    %msgbox('Warning:Point Colinearity produces degeneracy , deleting degenerated tets','Warning Window','warn');%Show Warning Message 
   % TetraDT(SelectedTetraIndex(:),:)=[];
%end
    

tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function EvolutionG_Callback(hObject, eventdata, handles)
% Evolution Quality Graph
global refineIteration;%Number of refinement steps
%Quality by the four function of geometry package
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;

[row iterationLevel]=size(refineIteration);
[row qualityValues]=size(meanQualityValues);


%Remove zero element from refineIteration
 if(refineIteration(1)==0 && iterationLevel>1)
refineIteration(1)=[];
 end
 
%Remove zero element from meanQuality Values
 if(meanQualityValues(1)==0 && qualityValues>1)
meanQualityValues(1)=[];
 end
 
[row qualityValues]=size(meanQualityValues2); 
 
%Remove zero element from meanQuality2 Values
 if(meanQualityValues2(1)==0 && qualityValues>1)
meanQualityValues2(1)=[];
 end
  
 
[row qualityValues]=size(meanQualityValues3); 
 
%Remove zero element from meanQuality3 Values
 if(meanQualityValues3(1)==0 && qualityValues>1)
meanQualityValues3(1)=[];
 end 

  
[row qualityValues]=size(meanQualityValues4); 
 
%Remove zero element from meanQuality4 Values
 if(meanQualityValues4(1)==0 && qualityValues>1)
meanQualityValues4(1)=[];
 end
 
figure; %Show figure
hold on; %hold graph on screen
plot(refineIteration,meanQualityValues,refineIteration,meanQualityValues2,refineIteration,meanQualityValues3,refineIteration,meanQualityValues4,'LineWidth',2); %plot bar graph
grid on; % Turn on grid lines for this plot
xlabel('# of Refinements','FontSize',12);  %label of axis x
axis([0 iterationLevel 0 1 ]); %Axis scaling
legend('Etha','Whiteh','Ratio','Solid Angle'); %Legend of graph , four function for calculating quality
title('Mesh Quality Evolution','FontSize',12); %Title of graph


% --------------------------------------------------------------------
function cubeTag_Callback(hObject, eventdata, handles)
%Load Cube from mat file
tic; %start timer to measure performance
load cubo;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Cube with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinenement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;

% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function barycentric_Callback(hObject, eventdata, handles)
% hObject    handle to barycentric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tetselect_Callback(hObject, eventdata, handles)
% hObject    handle to tetselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function GlobalRefine_Callback(hObject, eventdata, handles)
% Uniform 4T-Barycentric Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element

tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
  
   tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
   centroid=tetrahedron_centroid_3d(tetra); %Calculate Centroid Point
   TetraCoordinates=[TetraCoordinates;centroid]; %Add point to data structure
 
   %Finding Centroid Vertex ID
  [row column] =size(TetraCoordinates);
  CentroidID=row;
   
  %Composing the 4 new Tetrahedra joining the centroid point with the faces of the initial Tet
  % Face1 =Vertices 1,3,4
  % Face2= Vertices 1,2,3
  % Face3= Vertices 3,2,4
  % Face4= Vertices 1,2,4
      
  % Composing New Tetrahedron 1
  %Face compose with vertex 1 ,3 and 4
  Vertex1=TetraDT(i,1);
  Vertex2 =TetraDT(i,3);
  Vertex3 =TetraDT(i,4);
    
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 CentroidID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Face compose with vertex 1 ,2 and 3
  Vertex1=TetraDT(i,1);
  Vertex2 =TetraDT(i,2);
  Vertex3 =TetraDT(i,3);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
  %Composing New Tetrahedron 3
  %Face compose with vertex 3 ,2 and 4
  Vertex1=TetraDT(i,3);
  Vertex2 =TetraDT(i,2);
  Vertex3 =TetraDT(i,4);
  
  Tetrahedron3 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
 
  %Composing New Tetrahedron 4
  %Face compose with vertex 1 ,2 and 4
  Vertex1=TetraDT(i,1);
  Vertex2 =TetraDT(i,2);
  Vertex3 =TetraDT(i,4);
  
  Tetrahedron4 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];
   
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  TetraDT(1:Tetracount,:)=[];


  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function capTet_Callback(hObject, eventdata, handles)
%Load Cap Tetrahedron from mat file
tic; %start timer to measure performance
load capTet;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Cap Tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function liujoeTet_Callback(hObject, eventdata, handles)
%Load LiuJoe Tetrahedron from mat file
tic; %start timer to measure performance
load liujoeTet;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display LiuJoe Tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;

% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function needleTet_Callback(hObject, eventdata, handles)
%Load Needle Tetrahedron from mat file
tic; %start timer to measure performance
load needleTet;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Needle Tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function sliverTet_Callback(hObject, eventdata, handles)
%Load Sliver Tetrahedron from mat file
tic; %start timer to measure performance
load sliverTet;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Sliver Tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function beyTet_Callback(hObject, eventdata, handles)
%Load Bey Tetrahedron from mat file
tic; %start timer to measure performance
load tbey;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Bey Tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function wedgeTet_Callback(hObject, eventdata, handles)
%Load Wedge Tetrahedron from mat file
tic; %start timer to measure performance
load wedgeTet;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Wedge Tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;

% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function tetmesh_Callback(hObject, eventdata, handles)
%Load Sample Tetrahedral Mesh from mat file
tic; %start timer to measure performance
load tetmesh;
% Update global data structure from file loaded
global TetraDT;
global TetraCoordinates;
TetraDT=tet;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Sample Tetrahedral Mesh with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text
%print('-djpeg','-r300','InitialMesh');

% --------------------------------------------------------------------
function midpoint_Callback(hObject, eventdata, handles)
% hObject    handle to midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tets_Callback(hObject, eventdata, handles)
% hObject    handle to tets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function globalRefine_Callback(hObject, eventdata, handles)
%Delaunay Uniform Midpoint Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element

tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
%Inserting Midpoints on each tetrahedron edge

% Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
% Edge number:      1    2    3    4    5    6
 
  % Calculate Mid Point for each edge of tet
  % Edge 1 Midpoint
  midP1=(TetraCoordinates(TetraDT(i,1),:)+TetraCoordinates(TetraDT(i,2),:))/2;
  TetraCoordinates=[TetraCoordinates;midP1];   
  
  % Edge 2 Midpoint
  midP2=(TetraCoordinates(TetraDT(i,2),:)+TetraCoordinates(TetraDT(i,3),:))/2;
  TetraCoordinates=[TetraCoordinates;midP2];   
  
  % Edge 3 Midpoint
  midP3=(TetraCoordinates(TetraDT(i,3),:)+TetraCoordinates(TetraDT(i,1),:))/2;
  TetraCoordinates=[TetraCoordinates;midP3];
  
  % Edge 4 Midpoint
  midP4=(TetraCoordinates(TetraDT(i,2),:)+TetraCoordinates(TetraDT(i,4),:))/2;
  TetraCoordinates=[TetraCoordinates;midP4]; 
  
   % Edge 5 Midpoint
  midP5=(TetraCoordinates(TetraDT(i,3),:)+TetraCoordinates(TetraDT(i,4),:))/2;
  TetraCoordinates=[TetraCoordinates;midP5]; 
  
  % Edge 6 Midpoint
  midP6=(TetraCoordinates(TetraDT(i,4),:)+TetraCoordinates(TetraDT(i,1),:))/2;
  TetraCoordinates=[TetraCoordinates;midP6]; 
   
end

%Removed Duplicated points if they exist
TetraCoordinates=unique(TetraCoordinates,'rows');

%Apply Delaunay Refinement Algorithm 
dt = DelaunayTri(TetraCoordinates);

%Update Global Data Structure
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;

%Check Delaunay Bug , Colinearity of points produces degeneracy
%[Tetracount vertexNumber]= size(TetraDT); %number of element

 %SelectedTetraIndex=[]; %init variable
%for i=1:Tetracount %iterate over each tetrahedron
  % tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
  % volume=tetrahedron_volume_3d(tetra); %Calculate Volume 
   %if(volume==0) %if tet is degenerate
     %SelectedTetraIndex =[SelectedTetraIndex;i]; %store tetrahedra index in data structure next to deletion
   %end    
   
%end

%if (isempty(SelectedTetraIndex)==0) %Delete Degenerate Tets by delaunay colinearity
   % msgbox('Warning:Point Colinearity produces degeneracy , deleting degenerated tets','Warning Window','warn');%Show Warning Message 
   % TetraDT(SelectedTetraIndex(:),:)=[];
%end
    

tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function circumsphere_Callback(hObject, eventdata, handles)
% hObject    handle to circumsphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uniform_Callback(hObject, eventdata, handles)
%Delaunay Uniform Circumsphere Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
   tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
  [radious circum]=tetrahedron_circumsphere_3d(tetra); 
   TetraCoordinates=[TetraCoordinates;circum];
   
end

%Removed Duplicated points if they exist
TetraCoordinates=unique(TetraCoordinates,'rows');

%Apply Delaunay Refinement Algorithm 
dt = DelaunayTri(TetraCoordinates);

%Update Global Data Structure
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;

%Check Delaunay Bug , Colinearity of points produces degeneracy
%[Tetracount vertexNumber]= size(TetraDT); %number of element

 %SelectedTetraIndex=[]; %init variable
%for i=1:Tetracount %iterate over each tetrahedron
  % tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
  % volume=tetrahedron_volume_3d(tetra); %Calculate Volume 
   %if(volume==0) %if tet is degenerate
    % SelectedTetraIndex =[SelectedTetraIndex;i]; %store tetrahedra index in data structure next to deletion
   %end    
   
%end

%if (isempty(SelectedTetraIndex)==0) %Delete Degenerate Tets by delaunay colinearity
    %msgbox('Warning:Point Colinearity produces degeneracy , deleting degenerated tets','Warning Window','warn');%Show Warning Message 
   % TetraDT(SelectedTetraIndex(:),:)=[];
%end
    

tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function longest_Callback(hObject, eventdata, handles)
% hObject    handle to longest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
%Delaunay Uniform Longest Edge Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element

tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
   %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    TetraCoordinates=[TetraCoordinates;midP];
     
end

%Removed Duplicated points if they exist
TetraCoordinates=unique(TetraCoordinates,'rows');

%Apply Delaunay Refinement Algorithm 
dt = DelaunayTri(TetraCoordinates);

%Update Global Data Structure
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;

%Check Delaunay Bug , Colinearity of points produces degeneracy
%[Tetracount vertexNumber]= size(TetraDT); %number of element

 %SelectedTetraIndex=[]; %init variable
%for i=1:Tetracount %iterate over each tetrahedron
   %tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
   %volume=tetrahedron_volume_3d(tetra); %Calculate Volume 
   %if(volume==0) %if tet is degenerate
     %SelectedTetraIndex =[SelectedTetraIndex;i]; %store tetrahedra index in data structure next to deletion
   %end    
   
%end

%if (isempty(SelectedTetraIndex)==0) %Delete Degenerate Tets by delaunay colinearity
    %msgbox('Warning:Point Colinearity produces degeneracy , deleting degenerated tets','Warning Window','warn');%Show Warning Message 
   % TetraDT(SelectedTetraIndex(:),:)=[];
%end
    

tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
%Delaunay Uniform Insphere Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
   tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
  [radious insphere]=tetrahedron_insphere_3d(tetra);
  insphere=[insphere(1) insphere(2) insphere(3)]; %Converting to proper format , from row vector to column vector
  TetraCoordinates=[TetraCoordinates;insphere];
   
end

%Removed Duplicated points if they exist
TetraCoordinates=unique(TetraCoordinates,'rows');

%Apply Delaunay Refinement Algorithm 
dt = DelaunayTri(TetraCoordinates);

%Update Global Data Structure
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;

%Check Delaunay Bug , Colinearity of points produces degeneracy
%[Tetracount vertexNumber]= size(TetraDT); %number of element

 %SelectedTetraIndex=[]; %init variable
%for i=1:Tetracount %iterate over each tetrahedron
   %tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
   %volume=tetrahedron_volume_3d(tetra); %Calculate Volume 
   %if(volume==0) %if tet is degenerate
     %SelectedTetraIndex =[SelectedTetraIndex;i]; %store tetrahedra index in data structure next to deletion
   %end    
   
%end

%if (isempty(SelectedTetraIndex)==0) %Delete Degenerate Tets by delaunay colinearity
    %msgbox('Warning:Point Colinearity produces degeneracy , deleting degenerated tets','Warning Window','warn');%Show Warning Message 
   % TetraDT(SelectedTetraIndex(:),:)=[];
%end
    

tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
%Delaunay Uniform Longest Edge Trisection Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate First Equidistant Point for longest edge trisection
    midP1=TetraCoordinates(TetraDT(i,j1),:)*2/3+TetraCoordinates(TetraDT(i,j2),:)/3;
     
    TetraCoordinates=[TetraCoordinates;midP1]; %Add point to Data Structure
  
  % Calculate Second Equidistant Point for longest edge trisection
    midP2=TetraCoordinates(TetraDT(i,j1),:)/3+TetraCoordinates(TetraDT(i,j2),:)*2/3;
    
   TetraCoordinates=[TetraCoordinates;midP2];%Add point to Data Structure
   
end

%Removed Duplicated points if they exist
TetraCoordinates=unique(TetraCoordinates,'rows');

%Apply Delaunay Refinement Algorithm 
dt = DelaunayTri(TetraCoordinates);

%Update Global Data Structure
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;

%Check Delaunay Bug , Colinearity of points produces degeneracy
% [Tetracount vertexNumber]= size(TetraDT); %number of element

%  SelectedTetraIndex=[]; %init variable
% for i=1:Tetracount %iterate over each tetrahedron
%    tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
%    volume=tetrahedron_volume_3d(tetra); %Calculate Volume 
%    if(volume==0) %if tet is degenerate
%      SelectedTetraIndex =[SelectedTetraIndex;i]; %store tetrahedra index in data structure next to deletion
%    end    
   
% end

% if (isempty(SelectedTetraIndex)==0) %Delete Degenerate Tets by delaunay colinearity
    %msgbox('Warning:Point Colinearity produces degeneracy , deleting degenerated tets','Warning Window','warn');%Show Warning Message 
%     TetraDT(SelectedTetraIndex(:),:)=[];
% end
    

tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function FourTMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to FourTMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Uniform_Callback(hObject, eventdata, handles)
% 4T-LE Uniform Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
  
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
        
  %Composing the new four tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(i,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID Vertex3];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
   
  %Vertex1ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex1ID Vertex4];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
     
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(i,j2);
     
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID Vertex3];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
   %Composing Tet#4
   %Vertex1ID is the same as tet#3 since tet3 and tet4 are neighbors 
   Tetrahedron4 =[LEVertexID Vertex1ID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
         
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
   
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  TetraDT(1:Tetracount,:)=[];

    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info 


% --------------------------------------------------------------------
function eightTetrahedra_Callback(hObject, eventdata, handles)
% hObject    handle to eightTetrahedra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tetselection_Callback(hObject, eventdata, handles)
% hObject    handle to tetselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
% 8T-LE Uniform Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
  
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
   
   %Calculate Third Mid point of secondary edge  
  thirdMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(i,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==thirdMidP(1) & TetraCoordinates(:,2)==thirdMidP(2) & TetraCoordinates(:,3)==thirdMidP(3));
  
  if(isempty(r)==false)
  %Finding Third Vertex ID
  ThirdVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;thirdMidP];
  %Finding Third Vertex ID
  [row column] =size(TetraCoordinates);
  ThirdVertexID=row;
  end     
  
     
   %Calculate Fourth Mid point of secondary edge  
  fourMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(i,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fourMidP(1) & TetraCoordinates(:,2)==fourMidP(2) & TetraCoordinates(:,3)==fourMidP(3));
  
  if(isempty(r)==false)
  %Finding Fourth Vertex ID
  fourVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fourMidP];
  %Finding Fourht Vertex ID
  [row column] =size(TetraCoordinates);
  fourVertexID=row;
  end  
  
       
   %Calculate Fifth Mid point of secondary edge  
  fiveMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(i,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fiveMidP(1) & TetraCoordinates(:,2)==fiveMidP(2) & TetraCoordinates(:,3)==fiveMidP(3));
  
  if(isempty(r)==false)
  %Finding Fifth Vertex ID
  fiveVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fiveMidP];
  %Finding Fifht Vertex ID
  [row column] =size(TetraCoordinates);
  fiveVertexID=row;
  end  
  
  
  %Calculate Sixth Mid point of secondary edge  
  sixMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(i,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==sixMidP(1) & TetraCoordinates(:,2)==sixMidP(2) & TetraCoordinates(:,3)==sixMidP(3));
  
  if(isempty(r)==false)
  %Finding Sixth Vertex ID
  sixVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;sixMidP];
  %Finding Sixht Vertex ID
  [row column] =size(TetraCoordinates);
  sixVertexID=row;
  end  
  
  %Composing the new Eight Tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(i,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID ThirdVertexID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  Tetrahedron2 =[LEVertexID SecLEVertexID ThirdVertexID Vertex3];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
        
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID fourVertexID];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  %Composing Tet#4
   Tetrahedron4 =[LEVertexID fourVertexID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
  %Composing Tet#5
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
  Vertex1ID=TetraDT(i,j2);
  Tetrahedron5 =[LEVertexID Vertex1ID SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet5
  TetraDT =[TetraDT;Tetrahedron5];  
  
  %Composing Tet#6
 
  Tetrahedron6 =[LEVertexID Vertex3 SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet6
  TetraDT =[TetraDT;Tetrahedron6];  
  
  %Composing Tet#7
 
  Tetrahedron7 =[LEVertexID Vertex1ID SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet7
  TetraDT =[TetraDT;Tetrahedron7];  
  
    
  %Composing Tet#8
 
  Tetrahedron8 =[LEVertexID Vertex4 SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet8
  TetraDT =[TetraDT;Tetrahedron8];  
  
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
  
  %Saving Third Edge Vertexes ID for checking neighbor tetrahedra
  ThirdEdge =[Vertex3 Vertex1ID]; %Concatenate Third Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;ThirdEdge];
  
  
  %Saving Four Edge Vertexes ID for checking neighbor tetrahedra
  FourEdge =[Vertex4 Vertex1ID]; %Concatenate Four Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FourEdge];
  
  %Saving Fifth Edge Vertexes ID for checking neighbor tetrahedra
  FifthEdge =[Vertex3 Vertex2ID]; %Concatenate Fifth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FifthEdge];
  
  %Saving Sixth Edge Vertexes ID for checking neighbor tetrahedra
  SixthEdge =[Vertex4 Vertex2ID]; %Concatenate Sixth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SixthEdge];
  
    
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  TetraDT(1:Tetracount,:)=[];

    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info 


% --------------------------------------------------------------------
function delaunayMesh_Callback(hObject, eventdata, handles)
%Generate Delaunay Mesh from random points
tic; %start timer to measure performance
points = rand(50,3); %generating 50 random points
dt = DelaunayTri(points); %Triangulating by Delaunay Algorithm
% Update global data structure after Delaunay Triangulation of random points
global TetraDT;
global TetraCoordinates;
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Delaunay Mesh with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;

% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text 


% --------------------------------------------------------------------
function LE4section_Callback(hObject, eventdata, handles)
% hObject    handle to LE4section (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
% Uniform Longest Edge 4-section Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate First Equidistant Point for longest edge 4-section
    midP1=1/4*(TetraCoordinates(TetraDT(i,j1),:)*3+TetraCoordinates(TetraDT(i,j2),:));
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge 4-section
    midP2=2/3*midP1+TetraCoordinates(TetraDT(i,j2),:)/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
    % Calculate Third Equidistant Point for longest edge 4-section
    midP3=1/3*midP1+TetraCoordinates(TetraDT(i,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
  
  if(isempty(r)==false)
  %Third MidPoint ID
  thirdMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP3];
  %Third Midpoint ID
  [row column] =size(TetraCoordinates);
  thirdMidPID=row;
  end 
  
  
  
 %Composing The New Four Tetrahedra by 4-section of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
   %Composing New Tetrahedron 3
   Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j2);
  
  %Composing New Tetrahedron 4
   Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];
    
  
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  TetraDT(1:Tetracount,:)=[];

    
  %Algorithm Assure-Conformity of the tet mesh
   while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
         SelectedTetraIndex=[]; %init variable
         flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
         
         %Calculate LEPP
         %Sequential Search for finding neighbors tetrahedra set
         [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
            for k=1:Tetracount
                VertexIDs =TetraDT(k,:);
                indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
                indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
                
                if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                   SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                   flagHasNeighbor =true;
                end    
                
            end   
         
       if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
          SurroundingEdgeSet(1,:)=[]; 
       end 
           
      %Perform Longest Edge 4-section to selected Tetrahedra
      [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
      
  for i=1:Tetcount %iterate over each selected tetrahedron
          x=SelectedTetraIndex(i,1); %get selected tetrahedra index
    %Calculate Edge Length
   edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
   edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
   edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
   edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  
   Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
   
   [z,d] = max(Distance,[],2); %Obtain Max Distance
   
   %Saving Original Edge Order
   % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
   % Edge number:      1    2    3    4    5    6
   
   V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
   
   [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
   
   % Calculate First Equidistant Point for longest edge 4-section
     midP1=1/4*(TetraCoordinates(TetraDT(x,j1),:)*3+TetraCoordinates(TetraDT(x,j2),:));
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge 4-section
     midP2=2/3*midP1+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
     % Calculate Third Equidistant Point for longest edge 4-section
     midP3=1/3*midP1+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
   
   if(isempty(r)==false)
   %Third MidPoint ID
   thirdMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP3];
   %Third Midpoint ID
   [row column] =size(TetraCoordinates);
   thirdMidPID=row;
   end 
   
   
   
  %Composing The New Four Tetrahedra by 4-section of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
    %Composing New Tetrahedron 3
    Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
   
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   %Composing New Tetrahedron 4
    Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
   %Updating Data Structure with tet4
   TetraDT =[TetraDT;Tetrahedron4];
     
   
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];        
    
    
  end
  
    %iterate over each old selected tetrahedron updating Data Structure deleting old ones
     if (isempty(SelectedTetraIndex)==0)
      TetraDT(SelectedTetraIndex(:),:)=[];
     end
      
            
   end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths). Etha
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX.Whiteh
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.Ratio
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.Solid Angle
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
 


% --------------------------------------------------------------------
function Untitled_18_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_19_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_20_Callback(hObject, eventdata, handles)
%Delaunay Uniform Longest Edge 4-Section Refinement Algorithm
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance

for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
    %Calculate First Equidistant Point for longest edge 4-section
    midP1=1/4*(TetraCoordinates(TetraDT(i,j1),:)*3+TetraCoordinates(TetraDT(i,j2),:));
         
    TetraCoordinates=[TetraCoordinates;midP1]; %Add point to Data Structure
  
   % Calculate Second Equidistant Point for longest edge 4-section
    midP2=2/3*midP1+TetraCoordinates(TetraDT(i,j2),:)/3;
        
   TetraCoordinates=[TetraCoordinates;midP2];%Add point to Data Structure
   
    % Calculate Third Equidistant Point for longest edge 4-section
    midP3=1/3*midP1+TetraCoordinates(TetraDT(i,j2),:)*2/3;
    
   TetraCoordinates=[TetraCoordinates;midP3];%Add point to Data Structure
   
   
end

%Removed Duplicated points if they exist
TetraCoordinates=unique(TetraCoordinates,'rows');

%Apply Delaunay Refinement Algorithm 
dt = DelaunayTri(TetraCoordinates);

%Update Global Data Structure
TetraDT=dt.Triangulation;
TetraCoordinates=dt.X;

%Check Delaunay Bug , Colinearity of points produces degeneracy
% [Tetracount vertexNumber]= size(TetraDT); %number of element

%  SelectedTetraIndex=[]; %init variable
% for i=1:Tetracount %iterate over each tetrahedron
%    tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
%    volume=tetrahedron_volume_3d(tetra); %Calculate Volume 
%    if(volume==0) %if tet is degenerate
%      SelectedTetraIndex =[SelectedTetraIndex;i]; %store tetrahedra index in data structure next to deletion
%    end    
   
% end

% if (isempty(SelectedTetraIndex)==0) %Delete Degenerate Tets by delaunay colinearity
    %msgbox('Warning:Point Colinearity produces degeneracy , deleting degenerated tets','Warning Window','warn');%Show Warning Message 
%     TetraDT(SelectedTetraIndex(:),:)=[];
% end
    

tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


% --------------------------------------------------------------------
function Untitled_21_Callback(hObject, eventdata, handles)
% Local 3T-LE Refinement Algorithm by Input Value
%Input value by user
answer = inputdlg({'All Tetrahedra LE > Value will be refine:'},'Input Value');

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Value Input','Error Window','error');

elseif(value>=0) 

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
IterationIndex=[]; %init variable
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition if LE Distance > Input Value
 if(x<=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
       
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Longest Edge
  
  SortedDistance = sort(Distance,'descend'); %Sort distances in descending order
  SecondLE=SortedDistance(2); %get second longest edge distance
  
  %Check secondary distance value to corresponding edge given by vertexes id
  %Check the case when second longest edge distance is equal to primary
  %longest edge
  
  %When Secondary Longest Edge is equal to Edge1
  if(SecondLE==edge1)&&(j1~=1  ||  j2~=2)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,1),:)+TetraCoordinates(TetraDT(i,2),:))/2;
      j3=1; %save indices of secundary longest edge vertices
      j4=2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
  %When Secondary Longest Edge is equal to Edge2
  elseif(SecondLE==edge2) && (j1~=2  ||  j2~=3)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,2),:)+TetraCoordinates(TetraDT(i,3),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=3;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  %When Secondary Longest Edge is equal to Edge3
  elseif(SecondLE==edge3) && (j1~=3  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,3),:)+TetraCoordinates(TetraDT(i,1),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
      
  %When Secondary Longest Edge is equal to Edge4
  elseif(SecondLE==edge4) && (j1~=2  ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,2),:)+TetraCoordinates(TetraDT(i,4),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
 
   %When Secondary Longest Edge is equal to Edge5
  elseif(SecondLE==edge5) &&( j1~=3 ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,3),:)+TetraCoordinates(TetraDT(i,4),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  
  %When Secondary Longest Edge is equal to Edge6
  elseif(SecondLE==edge6) && (j1~=4  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(i,4),:)+TetraCoordinates(TetraDT(i,1),:))/2;
      j3=4; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  end %End of distance comparison with edges and computation of secondary midpoint 
      
  %Composing the new three tetrahedra
  % 3T-LE partition is composed of 2 subdivision patterns
  % Pattern 1: Longest Edge share a vertex with secondary longest edge
  % Pattern 2: Longest Edge is opposed to secondary longest edge
  %Applyinng proper pattern according to secondary longest edge position
    
  if(j1==j3 || j1==j4 || j2==j3 || j2==j4) %Pattern 1 is apply 
     
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
   Vertex3ID=TetraDT(i,j3);
   
  %Finding Vertex 4 ID
  VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)& VertexIDs ~= TetraDT(i,j3)& VertexIDs ~= TetraDT(i,j4));
  Vertex4ID = VertexIDs(b(1));
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
  Vertex3ID=TetraDT(i,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  %Finding Vertex 2ID
  VertexIDs =TetraDT(i,:);
   [a b] = find(VertexIDs ~=Vertex3 & VertexIDs ~= Vertex4 & VertexIDs ~= TetraDT(i,j3)& VertexIDs ~= TetraDT(i,j4));
  Vertex2ID = VertexIDs(b(1)); 
  
  
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  else %Pattern 2 is apply
  
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(i,j3);   
  
  %Finding Vertex 4 ID
  %Vertex 4 ID is one longest edge vertex id
  Vertex4ID=TetraDT(i,j2);
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(i,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
   %Finding Vertex 2ID , vertex 2 ID is one longest edge vertex
    Vertex2ID =TetraDT(i,j1);
    
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
    
  end    
  
        
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Longest Edge Vertexes ID for checking neighbor tetrahedra
  SecVertex1ID =TetraDT(i,j3);
  SecVertex2ID =TetraDT(i,j4);
  
  SecondaryEdge =[SecVertex1ID SecVertex2ID]; %Concatenate Secondary Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
     
  %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Value Input','Error Window','error');    
end    

end
 


% --------------------------------------------------------------------
function Untitled_22_Callback(hObject, eventdata, handles)
% Local 3T-LE Refinement Algorithm by Vertex ID
answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
[vertexID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
  
[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertexID>0  & vertexID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
 TV = vertexAttachments(trep,vertexID); %Return tetrahedra indices attached to specified vertex
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
    
 % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Longest Edge
  
  SortedDistance = sort(Distance,'descend'); %Sort distances in descending order
  SecondLE=SortedDistance(2); %get second longest edge distance
  
  %Check secondary distance value to corresponding edge given by vertexes id
  %Check the case when second longest edge distance is equal to primary
  %longest edge
  
  %When Secondary Longest Edge is equal to Edge1
  if(SecondLE==edge1)&&(j1~=1  ||  j2~=2)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,1),:)+TetraCoordinates(TetraDT(y,2),:))/2;
      j3=1; %save indices of secundary longest edge vertices
      j4=2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
  %When Secondary Longest Edge is equal to Edge2
  elseif(SecondLE==edge2) && (j1~=2  ||  j2~=3)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,2),:)+TetraCoordinates(TetraDT(y,3),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=3;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  %When Secondary Longest Edge is equal to Edge3
  elseif(SecondLE==edge3) && (j1~=3  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,3),:)+TetraCoordinates(TetraDT(y,1),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
      
  %When Secondary Longest Edge is equal to Edge4
  elseif(SecondLE==edge4) && (j1~=2  ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,2),:)+TetraCoordinates(TetraDT(y,4),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
 
   %When Secondary Longest Edge is equal to Edge5
  elseif(SecondLE==edge5) &&( j1~=3 ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,3),:)+TetraCoordinates(TetraDT(y,4),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  
  %When Secondary Longest Edge is equal to Edge6
  elseif(SecondLE==edge6) && (j1~=4  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,4),:)+TetraCoordinates(TetraDT(y,1),:))/2;
      j3=4; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  end %End of distance comparison with edges and computation of secondary midpoint 
      
  %Composing the new three tetrahedra
  % 3T-LE partition is composed of 2 subdivision patterns
  % Pattern 1: Longest Edge share a vertex with secondary longest edge
  % Pattern 2: Longest Edge is opposed to secondary longest edge
  %Applyinng proper pattern according to secondary longest edge position
    
  if(j1==j3 || j1==j4 || j2==j3 || j2==j4) %Pattern 1 is apply 
     
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
   Vertex3ID=TetraDT(y,j3);
   
  %Finding Vertex 4 ID
  VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)& VertexIDs ~= TetraDT(y,j3)& VertexIDs ~= TetraDT(y,j4));
  Vertex4ID = VertexIDs(b(1));
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
  Vertex3ID=TetraDT(y,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  %Finding Vertex 2ID
  VertexIDs =TetraDT(y,:);
   [a b] = find(VertexIDs ~=Vertex3 & VertexIDs ~= Vertex4 & VertexIDs ~= TetraDT(y,j3)& VertexIDs ~= TetraDT(y,j4));
  Vertex2ID = VertexIDs(b(1)); 
  
  
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  else %Pattern 2 is apply
  
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(y,j3);   
  
  %Finding Vertex 4 ID
  %Vertex 4 ID is one longest edge vertex id
  Vertex4ID=TetraDT(y,j2);
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(y,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
   %Finding Vertex 2ID , vertex 2 ID is one longest edge vertex
    Vertex2ID =TetraDT(y,j1);
    
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
    
  end    
  
        
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];  
  
  %Saving Secondary Longest Edge Vertexes ID for checking neighbor tetrahedra
  SecVertex1ID =TetraDT(y,j3);
  SecVertex2ID =TetraDT(y,j4);
  
  SecondaryEdge =[SecVertex1ID SecVertex2ID]; %Concatenate Secondary Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];  
   
end   

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)
% Local 3T-LE Refinement Algorithm by Edge
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
 %get Vertex1 ID   
[vertex1ID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    

 %get Vertex2 ID   
[vertex2ID status] =str2num(answer{2}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end  


[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertex1ID>0  & vertex1ID<=VertexCount & vertex2ID>0 & vertex2ID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
  %Test if Vertices are joined by Edge
  edge=isEdge(trep,vertex1ID,vertex2ID);
  
  if(edge==false)
     % Handle when vertices are not joined by edge
    msgbox('Vertices are not joined by Edge','Error Window','error');
    return;  
  end    
    
 TV = edgeAttachments(trep,vertex1ID,vertex2ID); %Return tetrahedra indices attached to specified edge defined by vertices
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Longest Edge
  
  SortedDistance = sort(Distance,'descend'); %Sort distances in descending order
  SecondLE=SortedDistance(2); %get second longest edge distance
  
  %Check secondary distance value to corresponding edge given by vertexes id
  %Check the case when second longest edge distance is equal to primary
  %longest edge
  
  %When Secondary Longest Edge is equal to Edge1
  if(SecondLE==edge1)&&(j1~=1  ||  j2~=2)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,1),:)+TetraCoordinates(TetraDT(y,2),:))/2;
      j3=1; %save indices of secundary longest edge vertices
      j4=2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
  %When Secondary Longest Edge is equal to Edge2
  elseif(SecondLE==edge2) && (j1~=2  ||  j2~=3)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,2),:)+TetraCoordinates(TetraDT(y,3),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=3;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  %When Secondary Longest Edge is equal to Edge3
  elseif(SecondLE==edge3) && (j1~=3  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,3),:)+TetraCoordinates(TetraDT(y,1),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
      
  %When Secondary Longest Edge is equal to Edge4
  elseif(SecondLE==edge4) && (j1~=2  ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,2),:)+TetraCoordinates(TetraDT(y,4),:))/2;
      j3=2; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
 
   %When Secondary Longest Edge is equal to Edge5
  elseif(SecondLE==edge5) &&( j1~=3 ||  j2~=4)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,3),:)+TetraCoordinates(TetraDT(y,4),:))/2;
      j3=3; %save indices of secundary longest edge vertices
      j4=4;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  
  %When Secondary Longest Edge is equal to Edge6
  elseif(SecondLE==edge6) && (j1~=4  ||  j2~=1)
     %Calculate Mid point of secondary longest edge
      secMidP=(TetraCoordinates(TetraDT(y,4),:)+TetraCoordinates(TetraDT(y,1),:))/2;
      j3=4; %save indices of secundary longest edge vertices
      j4=1;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end
  
  end %End of distance comparison with edges and computation of secondary midpoint 
      
  %Composing the new three tetrahedra
  % 3T-LE partition is composed of 2 subdivision patterns
  % Pattern 1: Longest Edge share a vertex with secondary longest edge
  % Pattern 2: Longest Edge is opposed to secondary longest edge
  %Applyinng proper pattern according to secondary longest edge position
    
  if(j1==j3 || j1==j4 || j2==j3 || j2==j4) %Pattern 1 is apply 
     
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
   Vertex3ID=TetraDT(y,j3);
   
  %Finding Vertex 4 ID
  VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)& VertexIDs ~= TetraDT(y,j3)& VertexIDs ~= TetraDT(y,j4));
  Vertex4ID = VertexIDs(b(1));
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID
  Vertex3ID=TetraDT(y,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  %Finding Vertex 2ID
  VertexIDs =TetraDT(y,:);
   [a b] = find(VertexIDs ~=Vertex3 & VertexIDs ~= Vertex4 & VertexIDs ~= TetraDT(y,j3)& VertexIDs ~= TetraDT(y,j4));
  Vertex2ID = VertexIDs(b(1)); 
  
  
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  else %Pattern 2 is apply
  
  %Composing Tet#1
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(y,j3);   
  
  %Finding Vertex 4 ID
  %Vertex 4 ID is one longest edge vertex id
  Vertex4ID=TetraDT(y,j2);
  
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  %Vertex1 ID and Vertex2 ID is Primary and Secondary Longest Edges Midpoint Vertex ID
  %Finding Vertex 3 ID , Vertex 3 ID is one secondary longest edge vertex
  %id
  Vertex3ID=TetraDT(y,j4);
  
  %Vertex4ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex3ID Vertex4ID];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  %Composing Tet#3
  
  %Vertex1 ID is Primary Longest Edge Midpoint Vertex ID
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
   %Finding Vertex 2ID , vertex 2 ID is one longest edge vertex
    Vertex2ID =TetraDT(y,j1);
    
  Tetrahedron3 =[LEVertexID Vertex2ID Vertex3 Vertex4];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
    
  end    
  
        
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Longest Edge Vertexes ID for checking neighbor tetrahedra
  SecVertex1ID =TetraDT(y,j3);
  SecVertex2ID =TetraDT(y,j4);
  
  SecondaryEdge =[SecVertex1ID SecVertex2ID]; %Concatenate Secondary Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
   
end             

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_24_Callback(hObject, eventdata, handles)
%  Local LE Trisection Refinement Algorithm By Value
answer = inputdlg({'All Tetrahedra LE > Value will be refine:'},'Input Value');

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Value Input','Error Window','error');

elseif(value>=0) 

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
IterationIndex=[]; %init variable
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition if LE Distance > Input Value
 if(x<=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
  
   % Calculate First Equidistant Point for longest edge trisection
    midP1=TetraCoordinates(TetraDT(i,j1),:)*2/3+TetraCoordinates(TetraDT(i,j2),:)/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge trisection
    midP2=TetraCoordinates(TetraDT(i,j1),:)/3+TetraCoordinates(TetraDT(i,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
 %Composing The New Three Tetrahedra by Trisection of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j2);
  
  Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  
  %Composing New Tetrahedron 3
   Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
     
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
 %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end

    
  %Algorithm Assure-Conformity of the tet mesh
  while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
        SelectedTetraIndex=[]; %init variable
        flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
        
        %Calculate LEPP
        %Sequential Search for finding neighbors tetrahedra set
        [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
           for k=1:Tetracount
               VertexIDs =TetraDT(k,:);
               indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
               indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
               
               if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                  SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                  flagHasNeighbor =true;
               end    
               
           end   
        
      if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
         SurroundingEdgeSet(1,:)=[]; 
      end 
          
     %Perform Longest Edge Trisection to selected Tetrahedra
     [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
     
 for i=1:Tetcount %iterate over each selected tetrahedron
         x=SelectedTetraIndex(i,1); %get selected tetrahedra index
   %Calculate Edge Length
  edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
  edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
  edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 
  Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
  
  [z,d] = max(Distance,[],2); %Obtain Max Distance
  
  %Saving Original Edge Order
  % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
  % Edge number:      1    2    3    4    5    6
  
  V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
  
  [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
          
    % Calculate First Equidistant Point for longest edge trisection
     midP1=TetraCoordinates(TetraDT(x,j1),:)*2/3+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge trisection
     midP2=TetraCoordinates(TetraDT(x,j1),:)/3+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
  %Composing The New Three Tetrahedra by Trisection of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
     %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
   
   %Composing New Tetrahedron 3
    Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
  
     
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
   
   
 end
 
   %iterate over each old selected tetrahedron updating Data Structure deleting old ones
    if (isempty(SelectedTetraIndex)==0)
     TetraDT(SelectedTetraIndex(:),:)=[];
    end
     
           
  end       
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Value Input','Error Window','error');    
end    

end

 


% --------------------------------------------------------------------
function Untitled_25_Callback(hObject, eventdata, handles)
%Local LE Trisection Refinement Algorithm By Vertex ID
%Input Vertex ID by user
answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  
    

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
[vertexID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
   msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
 
   

[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices


if(vertexID>0  && vertexID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
 TV = vertexAttachments(trep,vertexID); %Return tetrahedra indices attached to specified vertex
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP1=TetraCoordinates(TetraDT(y,j1),:)*2/3+TetraCoordinates(TetraDT(y,j2),:)/3;
       
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge trisection
    midP2=TetraCoordinates(TetraDT(y,j1),:)/3+TetraCoordinates(TetraDT(y,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
 %Composing The New Three Tetrahedra by Trisection of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j2);
  
  Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  
  %Composing New Tetrahedron 3
   Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
     
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end    

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
  while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
        SelectedTetraIndex=[]; %init variable
        flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
        
        %Calculate LEPP
        %Sequential Search for finding neighbors tetrahedra set
        [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
           for k=1:Tetracount
               VertexIDs =TetraDT(k,:);
               indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
               indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
               
               if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                  SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                  flagHasNeighbor =true;
               end    
               
           end   
        
      if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
         SurroundingEdgeSet(1,:)=[]; 
      end 
          
     %Perform Longest Edge Trisection to selected Tetrahedra
     [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
     
 for i=1:Tetcount %iterate over each selected tetrahedron
         x=SelectedTetraIndex(i,1); %get selected tetrahedra index
   %Calculate Edge Length
  edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
  edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
  edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 
  Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
  
  [z,d] = max(Distance,[],2); %Obtain Max Distance
  
  %Saving Original Edge Order
  % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
  % Edge number:      1    2    3    4    5    6
  
  V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
  
  [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
          
    % Calculate First Equidistant Point for longest edge trisection
     midP1=TetraCoordinates(TetraDT(x,j1),:)*2/3+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge trisection
     midP2=TetraCoordinates(TetraDT(x,j1),:)/3+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
  %Composing The New Three Tetrahedra by Trisection of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
     %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
   
   %Composing New Tetrahedron 3
    Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
  
     
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
   
   
 end
 
   %iterate over each old selected tetrahedron updating Data Structure deleting old ones
    if (isempty(SelectedTetraIndex)==0)
     TetraDT(SelectedTetraIndex(:),:)=[];
    end
     
           
  end       
 
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

  for i=1:Tetracount
          tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
         
          quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
          quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
          quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
          quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
          
  end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info


else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end

    



% --------------------------------------------------------------------
function Untitled_26_Callback(hObject, eventdata, handles)
%Local LE Trisection Refinement Algorithm By Edge
%Input Edge by Vertex1 ID and Vertex2 ID
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
 %get Vertex1 ID   
[vertex1ID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end   

%get Vertex2 ID   
[vertex2ID status] =str2num(answer{2}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end  

[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertex1ID>0  & vertex1ID<=VertexCount & vertex2ID>0 & vertex2ID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
  %Test if Vertices are joined by Edge
  edge=isEdge(trep,vertex1ID,vertex2ID);
  
  if(edge==false)
     % Handle when vertices are not joined by edge
    msgbox('Vertices are not joined by Edge','Error Window','error');
    return;  
  end    
    
 TV = edgeAttachments(trep,vertex1ID,vertex2ID); %Return tetrahedra indices attached to specified edge defined by vertices
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
       
  % Calculate First Point of longest edge
    midP1=TetraCoordinates(TetraDT(y,j1),:)*2/3+TetraCoordinates(TetraDT(y,j2),:)/3;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge trisection
    midP2=TetraCoordinates(TetraDT(y,j1),:)/3+TetraCoordinates(TetraDT(y,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
 %Composing The New Three Tetrahedra by Trisection of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j2);
  
  Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  
  %Composing New Tetrahedron 3
   Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
     
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end  

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
 %Algorithm Assure-Conformity of the tet mesh
  while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
        SelectedTetraIndex=[]; %init variable
        flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
        
        %Calculate LEPP
        %Sequential Search for finding neighbors tetrahedra set
        [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
           for k=1:Tetracount
               VertexIDs =TetraDT(k,:);
               indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
               indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
               
               if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                  SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                  flagHasNeighbor =true;
               end    
               
           end   
        
      if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
         SurroundingEdgeSet(1,:)=[]; 
      end 
          
     %Perform Longest Edge Trisection to selected Tetrahedra
     [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
     
 for i=1:Tetcount %iterate over each selected tetrahedron
         x=SelectedTetraIndex(i,1); %get selected tetrahedra index
   %Calculate Edge Length
  edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
  edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
  edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 
  Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
  
  [z,d] = max(Distance,[],2); %Obtain Max Distance
  
  %Saving Original Edge Order
  % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
  % Edge number:      1    2    3    4    5    6
  
  V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
  
  [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
          
    % Calculate First Equidistant Point for longest edge trisection
     midP1=TetraCoordinates(TetraDT(x,j1),:)*2/3+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge trisection
     midP2=TetraCoordinates(TetraDT(x,j1),:)/3+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
  %Composing The New Three Tetrahedra by Trisection of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
     %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
   
   %Composing New Tetrahedron 3
    Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
  
     
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
   
   
 end
 
   %iterate over each old selected tetrahedron updating Data Structure deleting old ones
    if (isempty(SelectedTetraIndex)==0)
     TetraDT(SelectedTetraIndex(:),:)=[];
    end
     
           
  end     
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_27_Callback(hObject, eventdata, handles)
% By Value Longest Edge 4-section Refinement Algorithm
%Input value by user
answer = inputdlg({'All Tetrahedra LE > Value will be refine:'},'Input Value');

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Value Input','Error Window','error');

elseif(value>=0) 

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
IterationIndex=[]; %init variable
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition if LE Distance > Input Value
 if(x<=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
 % Calculate First Equidistant Point for longest edge 4-section
  midP1=1/4*(TetraCoordinates(TetraDT(i,j1),:)*3+TetraCoordinates(TetraDT(i,j2),:));
    
 
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge 4-section
    midP2=2/3*midP1+TetraCoordinates(TetraDT(i,j2),:)/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
    % Calculate Third Equidistant Point for longest edge 4-section
    midP3=1/3*midP1+TetraCoordinates(TetraDT(i,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
  
  if(isempty(r)==false)
  %Third MidPoint ID
  thirdMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP3];
  %Third Midpoint ID
  [row column] =size(TetraCoordinates);
  thirdMidPID=row;
  end 
  
  
  
 %Composing The New Four Tetrahedra by 4-section of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
   %Composing New Tetrahedron 3
   Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j2);
  
  %Composing New Tetrahedron 4
   Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];
    
  
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end
    
    
  %Algorithm Assure-Conformity of the tet mesh
   while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
         SelectedTetraIndex=[]; %init variable
         flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
         
         %Calculate LEPP
         %Sequential Search for finding neighbors tetrahedra set
         [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
            for k=1:Tetracount
                VertexIDs =TetraDT(k,:);
                indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
                indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
                
                if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                   SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                   flagHasNeighbor =true;
                end    
                
            end   
         
       if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
          SurroundingEdgeSet(1,:)=[]; 
       end 
           
      %Perform Longest Edge 4-section to selected Tetrahedra
      [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
      
  for i=1:Tetcount %iterate over each selected tetrahedron
          x=SelectedTetraIndex(i,1); %get selected tetrahedra index
    %Calculate Edge Length
   edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
   edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
   edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
   edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  
   Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
   
   [z,d] = max(Distance,[],2); %Obtain Max Distance
   
   %Saving Original Edge Order
   % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
   % Edge number:      1    2    3    4    5    6
   
   V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
   
   [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
   
   % Calculate First Equidistant Point for longest edge 4-section
     midP1=1/4*(TetraCoordinates(TetraDT(x,j1),:)*3+TetraCoordinates(TetraDT(x,j2),:));
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge 4-section
     midP2=2/3*midP1+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
     % Calculate Third Equidistant Point for longest edge 4-section
     midP3=1/3*midP1+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
   
   if(isempty(r)==false)
   %Third MidPoint ID
   thirdMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP3];
   %Third Midpoint ID
   [row column] =size(TetraCoordinates);
   thirdMidPID=row;
   end 
   
   
   
  %Composing The New Four Tetrahedra by 4-section of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
    %Composing New Tetrahedron 3
    Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
   
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   %Composing New Tetrahedron 4
    Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
   %Updating Data Structure with tet4
   TetraDT =[TetraDT;Tetrahedron4];
     
   
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];        
    
    
  end
  
    %iterate over each old selected tetrahedron updating Data Structure deleting old ones
     if (isempty(SelectedTetraIndex)==0)
      TetraDT(SelectedTetraIndex(:),:)=[];
     end
      
            
   end   
   
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Value Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_28_Callback(hObject, eventdata, handles)
% Vertex ID LE 4-section Refinement Algorithm
%Input Vertex ID by user

answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
[vertexID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
  
[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertexID>0  & vertexID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
 TV = vertexAttachments(trep,vertexID); %Return tetrahedra indices attached to specified vertex
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
  % Calculate First Equidistant Point for longest edge 4-section
  midP1=1/4*(TetraCoordinates(TetraDT(y,j1),:)*3+TetraCoordinates(TetraDT(y,j2),:));
    
   
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge 4-section
    midP2=2/3*midP1+TetraCoordinates(TetraDT(y,j2),:)/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
    % Calculate Third Equidistant Point for longest edge 4-section
    midP3=1/3*midP1+TetraCoordinates(TetraDT(y,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
  
  if(isempty(r)==false)
  %Third MidPoint ID
  thirdMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP3];
  %Third Midpoint ID
  [row column] =size(TetraCoordinates);
  thirdMidPID=row;
  end 
  
  
  
 %Composing The New Four Tetrahedra by 4-section of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
   %Composing New Tetrahedron 3
   Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j2);
  
  %Composing New Tetrahedron 4
   Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];
    
  
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

 %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
  
 %Algorithm Assure-Conformity of the tet mesh
   while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
         SelectedTetraIndex=[]; %init variable
         flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
         
         %Calculate LEPP
         %Sequential Search for finding neighbors tetrahedra set
         [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
            for k=1:Tetracount
                VertexIDs =TetraDT(k,:);
                indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
                indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
                
                if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                   SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                   flagHasNeighbor =true;
                end    
                
            end   
         
       if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
          SurroundingEdgeSet(1,:)=[]; 
       end 
           
      %Perform Longest Edge 4-section to selected Tetrahedra
      [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
      
  for i=1:Tetcount %iterate over each selected tetrahedron
          x=SelectedTetraIndex(i,1); %get selected tetrahedra index
    %Calculate Edge Length
   edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
   edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
   edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
   edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  
   Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
   
   [z,d] = max(Distance,[],2); %Obtain Max Distance
   
   %Saving Original Edge Order
   % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
   % Edge number:      1    2    3    4    5    6
   
   V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
   
   [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
   
   % Calculate First Equidistant Point for longest edge 4-section
     midP1=1/4*(TetraCoordinates(TetraDT(x,j1),:)*3+TetraCoordinates(TetraDT(x,j2),:));
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge 4-section
     midP2=2/3*midP1+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
     % Calculate Third Equidistant Point for longest edge 4-section
     midP3=1/3*midP1+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
   
   if(isempty(r)==false)
   %Third MidPoint ID
   thirdMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP3];
   %Third Midpoint ID
   [row column] =size(TetraCoordinates);
   thirdMidPID=row;
   end 
   
   
   
  %Composing The New Four Tetrahedra by 4-section of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
    %Composing New Tetrahedron 3
    Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
   
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   %Composing New Tetrahedron 4
    Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
   %Updating Data Structure with tet4
   TetraDT =[TetraDT;Tetrahedron4];
     
   
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];        
    
    
  end
  
    %iterate over each old selected tetrahedron updating Data Structure deleting old ones
     if (isempty(SelectedTetraIndex)==0)
      TetraDT(SelectedTetraIndex(:),:)=[];
     end
      
            
   end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_29_Callback(hObject, eventdata, handles)
% By Edge Longest Edge LE 4-Section Refinement Algorithm
%Input Edge by Vertex1 ID and Vertex2 ID
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
 %get Vertex1 ID   
[vertex1ID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    

 %get Vertex2 ID   
[vertex2ID status] =str2num(answer{2}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end  


[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertex1ID>0  & vertex1ID<=VertexCount & vertex2ID>0 & vertex2ID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
  %Test if Vertices are joined by Edge
  edge=isEdge(trep,vertex1ID,vertex2ID);
  
  if(edge==false)
     % Handle when vertices are not joined by edge
    msgbox('Vertices are not joined by Edge','Error Window','error');
    return;  
  end    
    
 TV = edgeAttachments(trep,vertex1ID,vertex2ID); %Return tetrahedra indices attached to specified edge defined by vertices
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate First Equidistant Point for longest edge 4-section
    midP1=1/4*(TetraCoordinates(TetraDT(y,j1),:)*3+TetraCoordinates(TetraDT(y,j2),:));
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge 4-section
    midP2=2/3*midP1+TetraCoordinates(TetraDT(y,j2),:)/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
    % Calculate Third Equidistant Point for longest edge 4-section
    midP3=1/3*midP1+TetraCoordinates(TetraDT(y,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
  
  if(isempty(r)==false)
  %Third MidPoint ID
  thirdMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP3];
  %Third Midpoint ID
  [row column] =size(TetraCoordinates);
  thirdMidPID=row;
  end 
  
  
  
 %Composing The New Four Tetrahedra by 4-section of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
   %Composing New Tetrahedron 3
   Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
  %Finding Vertex 1 ID
  Vertex1=TetraDT(y,j2);
  
  %Composing New Tetrahedron 4
   Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];
    
  
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
  
 %Algorithm Assure-Conformity of the tet mesh
   while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
         SelectedTetraIndex=[]; %init variable
         flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
         
         %Calculate LEPP
         %Sequential Search for finding neighbors tetrahedra set
         [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
            for k=1:Tetracount
                VertexIDs =TetraDT(k,:);
                indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
                indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
                
                if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                   SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                   flagHasNeighbor =true;
                end    
                
            end   
         
       if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
          SurroundingEdgeSet(1,:)=[]; 
       end 
           
      %Perform Longest Edge 4-section to selected Tetrahedra
      [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
      
  for i=1:Tetcount %iterate over each selected tetrahedron
          x=SelectedTetraIndex(i,1); %get selected tetrahedra index
    %Calculate Edge Length
   edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
   edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
   edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
   edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
   edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  
   Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
   
   [z,d] = max(Distance,[],2); %Obtain Max Distance
   
   %Saving Original Edge Order
   % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
   % Edge number:      1    2    3    4    5    6
   
   V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
   
   [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
   
   % Calculate First Equidistant Point for longest edge 4-section
     midP1=1/4*(TetraCoordinates(TetraDT(x,j1),:)*3+TetraCoordinates(TetraDT(x,j2),:));
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge 4-section
     midP2=2/3*midP1+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
     % Calculate Third Equidistant Point for longest edge 4-section
     midP3=1/3*midP1+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP3(1) & TetraCoordinates(:,2)==midP3(2) & TetraCoordinates(:,3)==midP3(3));
   
   if(isempty(r)==false)
   %Third MidPoint ID
   thirdMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP3];
   %Third Midpoint ID
   [row column] =size(TetraCoordinates);
   thirdMidPID=row;
   end 
   
   
   
  %Composing The New Four Tetrahedra by 4-section of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   Tetrahedron2 =[MidPID SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
    %Composing New Tetrahedron 3
    Tetrahedron3 =[SecMidPID thirdMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
   
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   %Composing New Tetrahedron 4
    Tetrahedron4 =[thirdMidPID Vertex1 Vertex3 Vertex4];
   %Updating Data Structure with tet4
   TetraDT =[TetraDT;Tetrahedron4];
     
   
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];        
    
    
  end
  
    %iterate over each old selected tetrahedron updating Data Structure deleting old ones
     if (isempty(SelectedTetraIndex)==0)
      TetraDT(SelectedTetraIndex(:),:)=[];
     end
      
            
   end   
  
 
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_30_Callback(hObject, eventdata, handles)
% 4T-LE By Value Refinement Algorithm
%Input value by user
answer = inputdlg({'All Tetrahedra LE > Value will be refine:'},'Input Value');

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Value Input','Error Window','error');

elseif(value>=0) 

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
IterationIndex=[]; %init variable
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition if LE Distance > Input Value
 if(x<=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
       
  % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
        
  %Composing the new four tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(i,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID Vertex3];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
   
  %Vertex1ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex1ID Vertex4];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
     
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(i,j2);
     
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID Vertex3];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
   %Composing Tet#4
   %Vertex1ID is the same as tet#3 since tet3 and tet4 are neighbors 
   Tetrahedron4 =[LEVertexID Vertex1ID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
         
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
  %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Value Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_31_Callback(hObject, eventdata, handles)
% 4T-LE Vertex ID
%Input Vertex ID by user
answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
[vertexID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
  
[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertexID>0  & vertexID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
 TV = vertexAttachments(trep,vertexID); %Return tetrahedra indices attached to specified vertex
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
      
  % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
        
  %Composing the new four tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(y,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID Vertex3];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
   
  %Vertex1ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex1ID Vertex4];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
     
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(y,j2);
     
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID Vertex3];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
   %Composing Tet#4
   %Vertex1ID is the same as tet#3 since tet3 and tet4 are neighbors 
   Tetrahedron4 =[LEVertexID Vertex1ID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
         
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
   
end

   %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_32_Callback(hObject, eventdata, handles)
% Local 4T-LE Refinement Algorithm by Edge
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
 %get Vertex1 ID   
[vertex1ID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    

 %get Vertex2 ID   
[vertex2ID status] =str2num(answer{2}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end  


[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertex1ID>0  & vertex1ID<=VertexCount & vertex2ID>0 & vertex2ID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
  %Test if Vertices are joined by Edge
  edge=isEdge(trep,vertex1ID,vertex2ID);
  
  if(edge==false)
     % Handle when vertices are not joined by edge
    msgbox('Vertices are not joined by Edge','Error Window','error');
    return;  
  end    
    
 TV = edgeAttachments(trep,vertex1ID,vertex2ID); %Return tetrahedra indices attached to specified edge defined by vertices
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
    
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
      
        
  %Composing the new four tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(y,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID Vertex3];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
   
  %Vertex1ID is the same as tet#1 since tet1 and tet2 are neighbors
  
  Tetrahedron2 =[LEVertexID SecLEVertexID Vertex1ID Vertex4];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
     
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(y,j2);
     
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID Vertex3];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
   %Composing Tet#4
   %Vertex1ID is the same as tet#3 since tet3 and tet4 are neighbors 
   Tetrahedron4 =[LEVertexID Vertex1ID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
         
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
   
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end

  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_33_Callback(hObject, eventdata, handles)
% 8T-LE By Value Refinement Algorithm
%Input value by user
answer = inputdlg({'All Tetrahedra LE > Value will be refine:'},'Input Value');

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Value Input','Error Window','error');

elseif(value>=0) 

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
IterationIndex=[]; %init variable
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition if LE Distance > Input Value
 if(x<=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
    % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(i,j1),:)+TetraCoordinates(TetraDT(i,j2),:))/2;
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
   
   %Calculate Third Mid point of secondary edge  
  thirdMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(i,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==thirdMidP(1) & TetraCoordinates(:,2)==thirdMidP(2) & TetraCoordinates(:,3)==thirdMidP(3));
  
  if(isempty(r)==false)
  %Finding Third Vertex ID
  ThirdVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;thirdMidP];
  %Finding Third Vertex ID
  [row column] =size(TetraCoordinates);
  ThirdVertexID=row;
  end     
  
     
   %Calculate Fourth Mid point of secondary edge  
  fourMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(i,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fourMidP(1) & TetraCoordinates(:,2)==fourMidP(2) & TetraCoordinates(:,3)==fourMidP(3));
  
  if(isempty(r)==false)
  %Finding Fourth Vertex ID
  fourVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fourMidP];
  %Finding Fourht Vertex ID
  [row column] =size(TetraCoordinates);
  fourVertexID=row;
  end  
  
       
   %Calculate Fifth Mid point of secondary edge  
  fiveMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(i,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fiveMidP(1) & TetraCoordinates(:,2)==fiveMidP(2) & TetraCoordinates(:,3)==fiveMidP(3));
  
  if(isempty(r)==false)
  %Finding Fifth Vertex ID
  fiveVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fiveMidP];
  %Finding Fifht Vertex ID
  [row column] =size(TetraCoordinates);
  fiveVertexID=row;
  end  
  
  
  %Calculate Sixth Mid point of secondary edge  
  sixMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(i,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==sixMidP(1) & TetraCoordinates(:,2)==sixMidP(2) & TetraCoordinates(:,3)==sixMidP(3));
  
  if(isempty(r)==false)
  %Finding Sixth Vertex ID
  sixVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;sixMidP];
  %Finding Sixht Vertex ID
  [row column] =size(TetraCoordinates);
  sixVertexID=row;
  end  
  
  %Composing the new Eight Tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(i,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID ThirdVertexID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  Tetrahedron2 =[LEVertexID SecLEVertexID ThirdVertexID Vertex3];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
        
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID fourVertexID];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  %Composing Tet#4
   Tetrahedron4 =[LEVertexID fourVertexID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
  %Composing Tet#5
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
  Vertex1ID=TetraDT(i,j2);
  Tetrahedron5 =[LEVertexID Vertex1ID SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet5
  TetraDT =[TetraDT;Tetrahedron5];  
  
  %Composing Tet#6
 
  Tetrahedron6 =[LEVertexID Vertex3 SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet6
  TetraDT =[TetraDT;Tetrahedron6];  
  
  %Composing Tet#7
 
  Tetrahedron7 =[LEVertexID Vertex1ID SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet7
  TetraDT =[TetraDT;Tetrahedron7];  
  
    
  %Composing Tet#8
 
  Tetrahedron8 =[LEVertexID Vertex4 SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet8
  TetraDT =[TetraDT;Tetrahedron8];  
  
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
  
  %Saving Third Edge Vertexes ID for checking neighbor tetrahedra
  ThirdEdge =[Vertex3 Vertex1ID]; %Concatenate Third Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;ThirdEdge];
  
  
  %Saving Four Edge Vertexes ID for checking neighbor tetrahedra
  FourEdge =[Vertex4 Vertex1ID]; %Concatenate Four Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FourEdge];
  
  %Saving Fifth Edge Vertexes ID for checking neighbor tetrahedra
  FifthEdge =[Vertex3 Vertex2ID]; %Concatenate Fifth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FifthEdge];
  
  %Saving Sixth Edge Vertexes ID for checking neighbor tetrahedra
  SixthEdge =[Vertex4 Vertex2ID]; %Concatenate Sixth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SixthEdge];
  
 %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end
    
   
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Value Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_34_Callback(hObject, eventdata, handles)
%Vertex ID 8T-LE Uniform Refinement Algorithm
%Input Vertex ID by user
answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
[vertexID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
  
[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertexID>0  & vertexID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
 TV = vertexAttachments(trep,vertexID); %Return tetrahedra indices attached to specified vertex
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
   
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
   
   %Calculate Third Mid point of secondary edge  
  thirdMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(y,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==thirdMidP(1) & TetraCoordinates(:,2)==thirdMidP(2) & TetraCoordinates(:,3)==thirdMidP(3));
  
  if(isempty(r)==false)
  %Finding Third Vertex ID
  ThirdVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;thirdMidP];
  %Finding Third Vertex ID
  [row column] =size(TetraCoordinates);
  ThirdVertexID=row;
  end     
  
     
   %Calculate Fourth Mid point of secondary edge  
  fourMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(y,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fourMidP(1) & TetraCoordinates(:,2)==fourMidP(2) & TetraCoordinates(:,3)==fourMidP(3));
  
  if(isempty(r)==false)
  %Finding Fourth Vertex ID
  fourVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fourMidP];
  %Finding Fourht Vertex ID
  [row column] =size(TetraCoordinates);
  fourVertexID=row;
  end  
  
       
   %Calculate Fifth Mid point of secondary edge  
  fiveMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(y,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fiveMidP(1) & TetraCoordinates(:,2)==fiveMidP(2) & TetraCoordinates(:,3)==fiveMidP(3));
  
  if(isempty(r)==false)
  %Finding Fifth Vertex ID
  fiveVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fiveMidP];
  %Finding Fifht Vertex ID
  [row column] =size(TetraCoordinates);
  fiveVertexID=row;
  end  
  
  
  %Calculate Sixth Mid point of secondary edge  
  sixMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(y,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==sixMidP(1) & TetraCoordinates(:,2)==sixMidP(2) & TetraCoordinates(:,3)==sixMidP(3));
  
  if(isempty(r)==false)
  %Finding Sixth Vertex ID
  sixVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;sixMidP];
  %Finding Sixht Vertex ID
  [row column] =size(TetraCoordinates);
  sixVertexID=row;
  end  
  
  %Composing the new Eight Tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(y,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID ThirdVertexID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  Tetrahedron2 =[LEVertexID SecLEVertexID ThirdVertexID Vertex3];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
        
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID fourVertexID];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  %Composing Tet#4
   Tetrahedron4 =[LEVertexID fourVertexID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
  %Composing Tet#5
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
  Vertex1ID=TetraDT(y,j2);
  Tetrahedron5 =[LEVertexID Vertex1ID SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet5
  TetraDT =[TetraDT;Tetrahedron5];  
  
  %Composing Tet#6
 
  Tetrahedron6 =[LEVertexID Vertex3 SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet6
  TetraDT =[TetraDT;Tetrahedron6];  
  
  %Composing Tet#7
 
  Tetrahedron7 =[LEVertexID Vertex1ID SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet7
  TetraDT =[TetraDT;Tetrahedron7];  
  
    
  %Composing Tet#8
 
  Tetrahedron8 =[LEVertexID Vertex4 SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet8
  TetraDT =[TetraDT;Tetrahedron8];  
  
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
  
  %Saving Third Edge Vertexes ID for checking neighbor tetrahedra
  ThirdEdge =[Vertex3 Vertex1ID]; %Concatenate Third Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;ThirdEdge];
  
  
  %Saving Four Edge Vertexes ID for checking neighbor tetrahedra
  FourEdge =[Vertex4 Vertex1ID]; %Concatenate Four Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FourEdge];
  
  %Saving Fifth Edge Vertexes ID for checking neighbor tetrahedra
  FifthEdge =[Vertex3 Vertex2ID]; %Concatenate Fifth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FifthEdge];
  
  %Saving Sixth Edge Vertexes ID for checking neighbor tetrahedra
  SixthEdge =[Vertex4 Vertex2ID]; %Concatenate Sixth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SixthEdge];
  
    
end   
    
  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
              
 end   
  
tElapsed=toc; %stop timer

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_35_Callback(hObject, eventdata, handles)
% 8T-LE By Edge
%Input Edge by Vertex1 ID and Vertex2 ID
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  


%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
 %get Vertex1 ID   
[vertex1ID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    

 %get Vertex2 ID   
[vertex2ID status] =str2num(answer{2}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end  


[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element

if(vertex1ID>0  & vertex1ID<=VertexCount & vertex2ID>0 & vertex2ID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
  %Test if Vertices are joined by Edge
  edge=isEdge(trep,vertex1ID,vertex2ID);
  
  if(edge==false)
     % Handle when vertices are not joined by edge
    msgbox('Vertices are not joined by Edge','Error Window','error');
    return;  
  end    
    
 TV = edgeAttachments(trep,vertex1ID,vertex2ID); %Return tetrahedra indices attached to specified edge defined by vertices
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance

for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(y,1),1)-TetraCoordinates(TetraDT(y,2),1)).^2+(TetraCoordinates(TetraDT(y,1),2)-TetraCoordinates(TetraDT(y,2),2)).^2+(TetraCoordinates(TetraDT(y,1),3)-TetraCoordinates(TetraDT(y,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,3),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,3),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(y,2),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,2),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,2),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(y,3),1)-TetraCoordinates(TetraDT(y,4),1)).^2+(TetraCoordinates(TetraDT(y,3),2)-TetraCoordinates(TetraDT(y,4),2)).^2+(TetraCoordinates(TetraDT(y,3),3)-TetraCoordinates(TetraDT(y,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(y,4),1)-TetraCoordinates(TetraDT(y,1),1)).^2+(TetraCoordinates(TetraDT(y,4),2)-TetraCoordinates(TetraDT(y,1),2)).^2+(TetraCoordinates(TetraDT(y,4),3)-TetraCoordinates(TetraDT(y,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
           
  % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(y,j1),:)+TetraCoordinates(TetraDT(y,j2),:))/2;
    
    %Performing Longest Edge Bisection
    
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex ID
  LEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex ID
  [row column] =size(TetraCoordinates);
  LEVertexID=row;
  end    
  
  %Finding Secondary Edge Vertexes opposite to Longest Edge
   VertexIDs =TetraDT(y,:);
    [a b] = find(VertexIDs ~=TetraDT(y,j1)& VertexIDs ~= TetraDT(y,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
    
   %Calculate Mid point of secondary edge opposite to Longest Edge 
  secMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(Vertex4,:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==secMidP(1) & TetraCoordinates(:,2)==secMidP(2) & TetraCoordinates(:,3)==secMidP(3));
  
  if(isempty(r)==false)
  %Finding Secondary Vertex ID
  SecLEVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;secMidP];
  %Finding Secondary Vertex ID
  [row column] =size(TetraCoordinates);
  SecLEVertexID=row;
  end     
   
   %Calculate Third Mid point of secondary edge  
  thirdMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(y,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==thirdMidP(1) & TetraCoordinates(:,2)==thirdMidP(2) & TetraCoordinates(:,3)==thirdMidP(3));
  
  if(isempty(r)==false)
  %Finding Third Vertex ID
  ThirdVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;thirdMidP];
  %Finding Third Vertex ID
  [row column] =size(TetraCoordinates);
  ThirdVertexID=row;
  end     
  
     
   %Calculate Fourth Mid point of secondary edge  
  fourMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(y,j1),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fourMidP(1) & TetraCoordinates(:,2)==fourMidP(2) & TetraCoordinates(:,3)==fourMidP(3));
  
  if(isempty(r)==false)
  %Finding Fourth Vertex ID
  fourVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fourMidP];
  %Finding Fourht Vertex ID
  [row column] =size(TetraCoordinates);
  fourVertexID=row;
  end  
  
       
   %Calculate Fifth Mid point of secondary edge  
  fiveMidP=(TetraCoordinates(Vertex3,:)+TetraCoordinates(TetraDT(y,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==fiveMidP(1) & TetraCoordinates(:,2)==fiveMidP(2) & TetraCoordinates(:,3)==fiveMidP(3));
  
  if(isempty(r)==false)
  %Finding Fifth Vertex ID
  fiveVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;fiveMidP];
  %Finding Fifht Vertex ID
  [row column] =size(TetraCoordinates);
  fiveVertexID=row;
  end  
  
  
  %Calculate Sixth Mid point of secondary edge  
  sixMidP=(TetraCoordinates(Vertex4,:)+TetraCoordinates(TetraDT(y,j2),:))/2;
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==sixMidP(1) & TetraCoordinates(:,2)==sixMidP(2) & TetraCoordinates(:,3)==sixMidP(3));
  
  if(isempty(r)==false)
  %Finding Sixth Vertex ID
  sixVertexID=r;      
  else
  TetraCoordinates=[TetraCoordinates;sixMidP];
  %Finding Sixht Vertex ID
  [row column] =size(TetraCoordinates);
  sixVertexID=row;
  end  
  
  %Composing the new Eight Tetrahedra
  %Composing Tet#1
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
   Vertex1ID=TetraDT(y,j1);
    
  Tetrahedron1 =[LEVertexID SecLEVertexID Vertex1ID ThirdVertexID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing Tet#2
  Tetrahedron2 =[LEVertexID SecLEVertexID ThirdVertexID Vertex3];
  
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
    
  %Composing Tet#3
        
  Tetrahedron3 =[LEVertexID Vertex1ID SecLEVertexID fourVertexID];
  
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];  
      
  %Composing Tet#4
   Tetrahedron4 =[LEVertexID fourVertexID SecLEVertexID Vertex4];
  
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];  
  
  %Composing Tet#5
  %Vertex1 ID is one longest edge vertex 
  %Finding Vertex 1 ID
  Vertex1ID=TetraDT(y,j2);
  Tetrahedron5 =[LEVertexID Vertex1ID SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet5
  TetraDT =[TetraDT;Tetrahedron5];  
  
  %Composing Tet#6
 
  Tetrahedron6 =[LEVertexID Vertex3 SecLEVertexID fiveVertexID];
  
  %Updating Data Structure with tet6
  TetraDT =[TetraDT;Tetrahedron6];  
  
  %Composing Tet#7
 
  Tetrahedron7 =[LEVertexID Vertex1ID SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet7
  TetraDT =[TetraDT;Tetrahedron7];  
  
    
  %Composing Tet#8
 
  Tetrahedron8 =[LEVertexID Vertex4 SecLEVertexID sixVertexID];
  
  %Updating Data Structure with tet8
  TetraDT =[TetraDT;Tetrahedron8];  
  
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(y,j1);
  Vertex2ID =TetraDT(y,j2);
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
  %Saving Secondary Opposite Edge Vertexes ID for checking neighbor tetrahedra
  SecondaryEdge =[Vertex3 Vertex4]; %Concatenate Secondary Opposite Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SecondaryEdge];
  
  
  %Saving Third Edge Vertexes ID for checking neighbor tetrahedra
  ThirdEdge =[Vertex3 Vertex1ID]; %Concatenate Third Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;ThirdEdge];
  
  
  %Saving Four Edge Vertexes ID for checking neighbor tetrahedra
  FourEdge =[Vertex4 Vertex1ID]; %Concatenate Four Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FourEdge];
  
  %Saving Fifth Edge Vertexes ID for checking neighbor tetrahedra
  FifthEdge =[Vertex3 Vertex2ID]; %Concatenate Fifth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;FifthEdge];
  
  %Saving Sixth Edge Vertexes ID for checking neighbor tetrahedra
  SixthEdge =[Vertex4 Vertex2ID]; %Concatenate Sixth Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;SixthEdge];
  
    
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end 
    
  %Algorithm Assure-Conformity of the tet mesh
 while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
       SelectedTetraIndex=[]; %init variable
       flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
       
       %Calculate LEPP
       %Sequential Search for finding neighbors tetrahedra set
       [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
          for k=1:Tetracount
              VertexIDs =TetraDT(k,:);
              indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
              indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
              
              if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                 SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                 flagHasNeighbor =true;
              end    
              
          end   
       
     if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
        SurroundingEdgeSet(1,:)=[]; 
     end 
         
    %Perform Longest Edge Bisection to selected Tetrahedra
    [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
    
for i=1:Tetcount %iterate over each selected tetrahedron
        x=SelectedTetraIndex(i,1); %get selected tetrahedra index
  %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [z,d] = max(Distance,[],2); %Obtain Max Distance
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
   % Calculate Mid Point of longest edge
    midP=(TetraCoordinates(TetraDT(x,j1),:)+TetraCoordinates(TetraDT(x,j2),:))/2;
 
  %Performing Longest Edge Bisection
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP(1) & TetraCoordinates(:,2)==midP(2) & TetraCoordinates(:,3)==midP(3));
  
  if(isempty(r)==false)
  %Finding Vertex 2 ID
  Vertex2=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP];
  %Finding Vertex 2 ID
  [row column] =size(TetraCoordinates);
  Vertex2=row;
  end   
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(x,:);
    [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2));
  Vertex3 = VertexIDs(b(1));
  Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(x,j2);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
    
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(x,j1);
  Vertex2ID =TetraDT(x,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
  
end

  %iterate over each old selected tetrahedron updating Data Structure deleting old ones
   if (isempty(SelectedTetraIndex)==0)
    TetraDT(SelectedTetraIndex(:),:)=[];
   end
    
          
 end   
  
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];


%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_36_Callback(hObject, eventdata, handles)
%Local 4T-Barycentric Refinement Algorithm by value
%Input value by user
answer = inputdlg({'All Tetrahedra LE > Value will be refine:'},'Input Value');
 
%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Value Input','Error Window','error');
 
elseif(value>=0) 
 
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element
 
tic; %start timer for measuring performance
IterationIndex=[]; %init variable
 
for i=1:Tetracount %iterate over each tetrahedron
      %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 
 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition if LE Distance > Input Value
 if(x<=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    

   tetra= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
   centroid=tetrahedron_centroid_3d(tetra); %Calculate Centroid Point
   TetraCoordinates=[TetraCoordinates;centroid]; %Add point to data structure
 
   %Finding Centroid Vertex ID
  [row column] =size(TetraCoordinates);
  CentroidID=row;
   
  %Composing the 4 new Tetrahedra joining the centroid point with the faces of the initial Tet
  % Face1 =Vertices 1,3,4
  % Face2= Vertices 1,2,3
  % Face3= Vertices 3,2,4
  % Face4= Vertices 1,2,4
      
  % Composing New Tetrahedron 1
  %Face compose with vertex 1 ,3 and 4
  Vertex1=TetraDT(i,1);
  Vertex2 =TetraDT(i,3);
  Vertex3 =TetraDT(i,4);
    
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 CentroidID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Face compose with vertex 1 ,2 and 3
  Vertex1=TetraDT(i,1);
  Vertex2 =TetraDT(i,2);
  Vertex3 =TetraDT(i,3);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
  %Composing New Tetrahedron 3
  %Face compose with vertex 3 ,2 and 4
  Vertex1=TetraDT(i,3);
  Vertex2 =TetraDT(i,2);
  Vertex3 =TetraDT(i,4);
  
  Tetrahedron3 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
 
  %Composing New Tetrahedron 4
  %Face compose with vertex 1 ,2 and 4
  Vertex1=TetraDT(i,1);
  Vertex2 =TetraDT(i,2);
  Vertex3 =TetraDT(i,4);
  
  Tetrahedron4 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];

  %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];

   
end
 
  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end
    
 
 
  
tElapsed=toc; %stop timer
 
 
tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color
 
%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;
 
%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];
 
%Calculating Mean Quality Value
 
[Tetracount vertexNumber]= size(TetraDT); %number of element
 
for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end
 
meanValue=mean(quality); %Mean Value of Quality
 
%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];
 
 
meanValue=mean(quality2); %Mean Value of Quality 2
 
%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];
 
meanValue=mean(quality3); %Mean Value of Quality 3
 
%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];
 
meanValue=mean(quality4); %Mean Value of Quality 4
 
%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];
 
 
%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
 
% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
 
%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
         
else
 msgbox('Wrong Value Input','Error Window','error');    
end    
 
end
 


% --------------------------------------------------------------------
function Untitled_37_Callback(hObject, eventdata, handles)
% Local 4T- Baricentric Refinement Algorithm by Vertex ID
%Input Vertex ID by user
answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');
 
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  
 
 
%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
[vertexID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
  
[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element
 
if(vertexID>0  & vertexID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
 TV = vertexAttachments(trep,vertexID); %Return tetrahedra indices attached to specified vertex
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance

 
for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
     tetra= [TetraCoordinates(TetraDT(y,1),:)', TetraCoordinates(TetraDT(y,2),:)', TetraCoordinates(TetraDT(y,3),:)', TetraCoordinates(TetraDT(y,4),:)'];
   centroid=tetrahedron_centroid_3d(tetra); %Calculate Centroid Point
   TetraCoordinates=[TetraCoordinates;centroid]; %Add point to data structure
 
   %Finding Centroid Vertex ID
  [row column] =size(TetraCoordinates);
  CentroidID=row;
   
  %Composing the 4 new Tetrahedra joining the centroid point with the faces of the initial Tet
  % Face1 =Vertices 1,3,4
  % Face2= Vertices 1,2,3
  % Face3= Vertices 3,2,4
  % Face4= Vertices 1,2,4
      
  % Composing New Tetrahedron 1
  %Face compose with vertex 1 ,3 and 4
  Vertex1=TetraDT(y,1);
  Vertex2 =TetraDT(y,3);
  Vertex3 =TetraDT(y,4);
    
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 CentroidID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Face compose with vertex 1 ,2 and 3
  Vertex1=TetraDT(y,1);
  Vertex2 =TetraDT(y,2);
  Vertex3 =TetraDT(y,3);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
  %Composing New Tetrahedron 3
  %Face compose with vertex 3 ,2 and 4
  Vertex1=TetraDT(y,3);
  Vertex2 =TetraDT(y,2);
  Vertex3 =TetraDT(y,4);
  
  Tetrahedron3 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
 
  %Composing New Tetrahedron 4
  %Face compose with vertex 1 ,2 and 4
  Vertex1=TetraDT(y,1);
  Vertex2 =TetraDT(y,2);
  Vertex3 =TetraDT(y,4);
  
  Tetrahedron4 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];

  
end
 
  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
 
  
tElapsed=toc; %stop timer
 
 
tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color
 
%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;
 
%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];
 
%Calculating Mean Quality Value
 
[Tetracount vertexNumber]= size(TetraDT); %number of element
 
for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end
 
meanValue=mean(quality); %Mean Value of Quality
 
%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];
 
 
meanValue=mean(quality2); %Mean Value of Quality 2
 
%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];
 
meanValue=mean(quality3); %Mean Value of Quality 3
 
%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];
 
meanValue=mean(quality4); %Mean Value of Quality 4
 
%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];
 
%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
 
% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
 
%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    
 
end


% --------------------------------------------------------------------
function Untitled_38_Callback(hObject, eventdata, handles)
% Local 4T-Baricentric Refinement Algorithm by Edge
%Input Edge by Vertex1 ID and Vertex2 ID
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');
 
global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points  
 
 
%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  return;  
else
 %get Vertex1 ID   
[vertex1ID status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end    
 
 %get Vertex2 ID   
[vertex2ID status] =str2num(answer{2}); %Convert String to number
if ~status
    % Handle empty vertex id returned for unsuccessful conversion
    msgbox('Wrong Vertex ID Input','Error Window','error');
    return;
end  
 
 
[VertexCount vertexColumn] =size(TetraCoordinates); %number of vertices
[Tetracount vertexNumber]= size(TetraDT); %number of element
 
if(vertex1ID>0  & vertex1ID<=VertexCount & vertex2ID>0 & vertex2ID<=VertexCount) %if vertex ID is in range
  %Load global data structure into TriRep object
  trep = TriRep(TetraDT,TetraCoordinates);
  %Test if Vertices are joined by Edge
  edge=isEdge(trep,vertex1ID,vertex2ID);
  
  if(edge==false)
     % Handle when vertices are not joined by edge
    msgbox('Vertices are not joined by Edge','Error Window','error');
    return;  
  end    
    
 TV = edgeAttachments(trep,vertex1ID,vertex2ID); %Return tetrahedra indices attached to specified edge defined by vertices
 tetSelection =TV{:}; %Convert Cell Array to matrix
 [row tetColumn]=size(tetSelection);
 
tic; %start timer for measuring performance

 
for i=1:tetColumn %iterate over each selected tetrahedron
    y=tetSelection(1,i); %get selected tetra index
    tetra= [TetraCoordinates(TetraDT(y,1),:)', TetraCoordinates(TetraDT(y,2),:)', TetraCoordinates(TetraDT(y,3),:)', TetraCoordinates(TetraDT(y,4),:)'];
   centroid=tetrahedron_centroid_3d(tetra); %Calculate Centroid Point
   TetraCoordinates=[TetraCoordinates;centroid]; %Add point to data structure
 
   %Finding Centroid Vertex ID
  [row column] =size(TetraCoordinates);
  CentroidID=row;
   
  %Composing the 4 new Tetrahedra joining the centroid point with the faces of the initial Tet
  % Face1 =Vertices 1,3,4
  % Face2= Vertices 1,2,3
  % Face3= Vertices 3,2,4
  % Face4= Vertices 1,2,4
      
  % Composing New Tetrahedron 1
  %Face compose with vertex 1 ,3 and 4
  Vertex1=TetraDT(y,1);
  Vertex2 =TetraDT(y,3);
  Vertex3 =TetraDT(y,4);
    
  Tetrahedron1 =[Vertex1 Vertex2 Vertex3 CentroidID];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Face compose with vertex 1 ,2 and 3
  Vertex1=TetraDT(y,1);
  Vertex2 =TetraDT(y,2);
  Vertex3 =TetraDT(y,3);
  
  Tetrahedron2 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
 
  %Composing New Tetrahedron 3
  %Face compose with vertex 3 ,2 and 4
  Vertex1=TetraDT(y,3);
  Vertex2 =TetraDT(y,2);
  Vertex3 =TetraDT(y,4);
  
  Tetrahedron3 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
 
  %Composing New Tetrahedron 4
  %Face compose with vertex 1 ,2 and 4
  Vertex1=TetraDT(y,1);
  Vertex2 =TetraDT(y,2);
  Vertex3 =TetraDT(y,4);
  
  Tetrahedron4 =[Vertex1 Vertex2 Vertex3 CentroidID];
  %Updating Data Structure with tet4
  TetraDT =[TetraDT;Tetrahedron4];
   
 
end
 
  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(tetSelection)==0)
  TetraDT(tetSelection(:),:)=[];
  end
    
  
tElapsed=toc; %stop timer
 
 
tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color
 
%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;
 
%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];
 
%Calculating Mean Quality Value
 
[Tetracount vertexNumber]= size(TetraDT); %number of element
 
for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end
 
meanValue=mean(quality); %Mean Value of Quality
 
%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];
 
 
meanValue=mean(quality2); %Mean Value of Quality 2
 
%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];
 
meanValue=mean(quality3); %Mean Value of Quality 3
 
%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];
 
meanValue=mean(quality4); %Mean Value of Quality 4
 
%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];
 
 
%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
 
% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
 
%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    
 
end


% --------------------------------------------------------------------
function Untitled_39_Callback(hObject, eventdata, handles)
% Quality Comparison Evolution of selected refinement algorithms
% Show Dialog
qualityComparisonDialog();


% --------------------------------------------------------------------
function Untitled_40_Callback(hObject, eventdata, handles)
%Load Well Shaped Tetrahedron from Coordinates
tic; %start timer to measure performance
Tes=[1 2 3 4];
X=[0.0 0.0 0.0
     4.0 2.0 2.0
     1.0 5.0 0.0
     0.5 0.5 5.0];
% Update global data structure from tet coordinates
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display well shaped tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinenement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function Untitled_41_Callback(hObject, eventdata, handles)
%Load Rectangular Tetrahedron from Coordinates
tic; %start timer to measure performance
Tes=[1 2 3 4];
X=[0.0 0.0 0.0
     4.0 0.0 0.0
     0.0 4.0 0.0
     0.0 0.0 4.0];
% Update global data structure from tet coordinates
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Rectangular tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinenement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function Untitled_42_Callback(hObject, eventdata, handles)
%Load Distorted Tetrahedron from Coordinates
tic; %start timer to measure performance
Tes=[1 2 3 4];
X=[0.0 0.0 0.0
     0.5 0.0 0.0
     1.0 5.0 2.0
     0.5 0.5 5.0];
% Update global data structure from tet coordinates
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Distorted tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinenement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function Untitled_43_Callback(hObject, eventdata, handles)
%Load Equilateral Tetrahedron from Coordinates
tic; %start timer to measure performance
Tes=[1 2 3 4];
X=[0.0 0.0 0.0
     3.4641016151377544 0.0 0.0
     1.7320508075688772 3.0 0.0
     1.7320508075688772 1.0 2.8284271247461900];
% Update global data structure from tet coordinates
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Equilateral tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinenement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function Untitled_44_Callback(hObject, eventdata, handles)
%Load Quasi-Equilateral Tetrahedron from Coordinates
tic; %start timer to measure performance
Tes=[1 2 3 4];
X=[0.0 0.0 0.0
     3.46 0.0 0.0
     1.73 3.0 0.0
     1.73 1.0 2.83];
% Update global data structure from tet coordinates
global TetraDT;
global TetraCoordinates;
TetraDT=Tes;
TetraCoordinates=X;
tElapsed=toc; %stop timer 

tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Display Quasi-Equilateral tetrahedron with removed face color
set(handles.RefinementMenuItem,'Enable','on'); % Enable Refinenement Menu Item
set(handles.QualityMenuItem,'Enable','on'); % Enable Quality Menu Item
set(handles.ViewMenuItem,'Enable','on'); %Enable View Menu Item

%Initialization of global variable for refinement level and mean quality values 
global refineIteration; %Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2; 
global meanQualityValues3; 
global meanQualityValues4;
refineIteration=0; %init to zero
meanQualityValues=0; %init to zero
meanQualityValues2=0;
meanQualityValues3=0;
meanQualityValues4=0;


% Updating GUI
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info
set(handles.tetLabel,'Visible','on'); % Enable Visible static text

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info
set(handles.vertLabel,'Visible','on'); % Enable Visible static text

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String concatenation
set(handles.timeLabel,'String',text); % Update Time info
set(handles.timeLabel,'Visible','on'); % Enable Visible static text


% --------------------------------------------------------------------
function Untitled_45_Callback(hObject, eventdata, handles)
% Local LE-Trisection Refinement Algorithm By Quality Threshold
answer = inputdlg({'All Tetrahedra Quality Value (Etha) < Value will be refine:'},'Input Value');

%Check if answer is empty , user click Cancel Button
if (isempty(answer)==true)
  %Do Nothing  
else
[value status] =str2num(answer{1}); %Convert String to number
if ~status
    % Handle empty value returned for unsuccessful conversion
    msgbox('Wrong Quality Value Input','Error Window','error');

elseif(value>=0 && value<=1) 

global TetraDT; %global variable tetrahedral triangulation
global TetraCoordinates; % global variable triangulation points
[Tetracount vertexNumber]= size(TetraDT); %number of element


tic; %start timer for measuring performance
SurroundingEdgeSet =[]; %preallocating for improving performance
IterationIndex=[]; %init variable
for i=1:Tetracount %iterate over each tetrahedron
    %Calculate Edge Length
 edge1 =sqrt((TetraCoordinates(TetraDT(i,1),1)-TetraCoordinates(TetraDT(i,2),1)).^2+(TetraCoordinates(TetraDT(i,1),2)-TetraCoordinates(TetraDT(i,2),2)).^2+(TetraCoordinates(TetraDT(i,1),3)-TetraCoordinates(TetraDT(i,2),3)).^2);   
 edge2 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,3),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,3),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,3),3)).^2);
 edge3 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 
 edge4 =sqrt((TetraCoordinates(TetraDT(i,2),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,2),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,2),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge5 =sqrt((TetraCoordinates(TetraDT(i,3),1)-TetraCoordinates(TetraDT(i,4),1)).^2+(TetraCoordinates(TetraDT(i,3),2)-TetraCoordinates(TetraDT(i,4),2)).^2+(TetraCoordinates(TetraDT(i,3),3)-TetraCoordinates(TetraDT(i,4),3)).^2); 
 edge6 =sqrt((TetraCoordinates(TetraDT(i,4),1)-TetraCoordinates(TetraDT(i,1),1)).^2+(TetraCoordinates(TetraDT(i,4),2)-TetraCoordinates(TetraDT(i,1),2)).^2+(TetraCoordinates(TetraDT(i,4),3)-TetraCoordinates(TetraDT(i,1),3)).^2); 

 Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
 
 [x,d] = max(Distance,[],2); %Obtain Max Distance
 
 % Check Condition IF Quality Value by Etha < Input Value
    tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
    quality=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).Etha Quality Measure.
 
 if(quality>=value) %Skip Tetrahedra , Jump to next iteration if true
   continue   
 end    
 
 %Saving Original Edge Order
 % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
 % Edge number:      1    2    3    4    5    6
 
 V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
 
 [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
         
  
   % Calculate First Equidistant Point for longest edge trisection
    midP1=TetraCoordinates(TetraDT(i,j1),:)*2/3+TetraCoordinates(TetraDT(i,j2),:)/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
  
  if(isempty(r)==false)
  %First MidPoint ID
  MidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP1];
  %First MidPoint ID
  [row column] =size(TetraCoordinates);
  MidPID=row;
  end    
  
  % Calculate Second Equidistant Point for longest edge trisection
    midP2=TetraCoordinates(TetraDT(i,j1),:)/3+TetraCoordinates(TetraDT(i,j2),:)*2/3;
    
  
  % Add new point into TetraCoordinates Vertex Matrix , checking if point
  % is duplicated in data structure
  
  [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
  
  if(isempty(r)==false)
  %Second MidPoint ID
  SecMidPID=r;      
  else
  TetraCoordinates=[TetraCoordinates;midP2];
  %Second Midpoint ID
  [row column] =size(TetraCoordinates);
  SecMidPID=row;
  end    
  
 %Composing The New Three Tetrahedra by Trisection of Longest Edge
  
  % Composing New Tetrahedron 1
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j1);
  
  %Finding Vertex 3 and 4 ID
    VertexIDs =TetraDT(i,:);
    [a b] = find(VertexIDs ~=TetraDT(i,j1)& VertexIDs ~= TetraDT(i,j2)); 
  
    Vertex3 = VertexIDs(b(1));
    Vertex4 = VertexIDs(b(2));
  
  Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
  
  %Updating Data Structure with tet1
  TetraDT =[TetraDT;Tetrahedron1];
  
  %Composing New Tetrahedron 2
  %Finding Vertex 1 ID
  Vertex1=TetraDT(i,j2);
  
  Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
    %Updating Data Structure with tet2
  TetraDT =[TetraDT;Tetrahedron2];
  
  
  %Composing New Tetrahedron 3
   Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
  %Updating Data Structure with tet3
  TetraDT =[TetraDT;Tetrahedron3];
  
     
  %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
  Vertex1ID =TetraDT(i,j1);
  Vertex2ID =TetraDT(i,j2);
  
  Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
  SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
  
 %Store for loop index to update data structure
  IterationIndex=[IterationIndex;i];
  
end

  %iterate over each old tetrahedron updating Data Structure deleting old ones
  if(isempty(IterationIndex)==0)
  TetraDT(IterationIndex(:),:)=[];
  end

    
  %Algorithm Assure-Conformity of the tet mesh
  while(isempty(SurroundingEdgeSet)==0) %while there exits at least one surrounding edge
        SelectedTetraIndex=[]; %init variable
        flagHasNeighbor =false; %flag variable to test if tetrahedra has neighbor     
        
        %Calculate LEPP
        %Sequential Search for finding neighbors tetrahedra set
        [Tetracount vertexNumber]= size(TetraDT); %number of element in data structure
           for k=1:Tetracount
               VertexIDs =TetraDT(k,:);
               indic = find(VertexIDs ==SurroundingEdgeSet(1,1)); 
               indic2 =find(VertexIDs ==SurroundingEdgeSet(1,2));
               
               if(isempty(indic)==0 & isempty(indic2)==0) %if it is neighbor tetrahedra       
                  SelectedTetraIndex =[SelectedTetraIndex;k]; %store tetrahedra index in data structure , next to refinement
                  flagHasNeighbor =true;
               end    
               
           end   
        
      if (flagHasNeighbor ==false) %if no neighbor tetrahedra exist
         SurroundingEdgeSet(1,:)=[]; 
      end 
          
     %Perform Longest Edge Trisection to selected Tetrahedra
     [Tetcount column]= size(SelectedTetraIndex); %number of selected tetrahedra
     
 for i=1:Tetcount %iterate over each selected tetrahedron
         x=SelectedTetraIndex(i,1); %get selected tetrahedra index
   %Calculate Edge Length
  edge1 =sqrt((TetraCoordinates(TetraDT(x,1),1)-TetraCoordinates(TetraDT(x,2),1)).^2+(TetraCoordinates(TetraDT(x,1),2)-TetraCoordinates(TetraDT(x,2),2)).^2+(TetraCoordinates(TetraDT(x,1),3)-TetraCoordinates(TetraDT(x,2),3)).^2);   
  edge2 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,3),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,3),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,3),3)).^2);
  edge3 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
  edge4 =sqrt((TetraCoordinates(TetraDT(x,2),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,2),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,2),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge5 =sqrt((TetraCoordinates(TetraDT(x,3),1)-TetraCoordinates(TetraDT(x,4),1)).^2+(TetraCoordinates(TetraDT(x,3),2)-TetraCoordinates(TetraDT(x,4),2)).^2+(TetraCoordinates(TetraDT(x,3),3)-TetraCoordinates(TetraDT(x,4),3)).^2); 
  edge6 =sqrt((TetraCoordinates(TetraDT(x,4),1)-TetraCoordinates(TetraDT(x,1),1)).^2+(TetraCoordinates(TetraDT(x,4),2)-TetraCoordinates(TetraDT(x,1),2)).^2+(TetraCoordinates(TetraDT(x,4),3)-TetraCoordinates(TetraDT(x,1),3)).^2); 
 
  Distance = [edge1 edge2 edge3 edge4 edge5 edge6];
  
  [z,d] = max(Distance,[],2); %Obtain Max Distance
  
  %Saving Original Edge Order
  % Edge per vertex: 1-2, 2-3, 3-1, 2-4, 3-4, 4-1
  % Edge number:      1    2    3    4    5    6
  
  V=sparse([1 2 3 2 3 4],[2 3 1 4 4 1],[1 2 3 4 5 6],6,6);
  
  [j1,j2]=find(V==d); %find the corresponding vertexes indexes to longest edge
          
    % Calculate First Equidistant Point for longest edge trisection
     midP1=TetraCoordinates(TetraDT(x,j1),:)*2/3+TetraCoordinates(TetraDT(x,j2),:)/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP1(1) & TetraCoordinates(:,2)==midP1(2) & TetraCoordinates(:,3)==midP1(3));
   
   if(isempty(r)==false)
   %First MidPoint ID
   MidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP1];
   %First MidPoint ID
   [row column] =size(TetraCoordinates);
   MidPID=row;
   end    
   
   % Calculate Second Equidistant Point for longest edge trisection
     midP2=TetraCoordinates(TetraDT(x,j1),:)/3+TetraCoordinates(TetraDT(x,j2),:)*2/3;
     
   
   % Add new point into TetraCoordinates Vertex Matrix , checking if point
   % is duplicated in data structure
   
   [r,c]=find(TetraCoordinates(:,1)==midP2(1) & TetraCoordinates(:,2)==midP2(2) & TetraCoordinates(:,3)==midP2(3));
   
   if(isempty(r)==false)
   %Second MidPoint ID
   SecMidPID=r;      
   else
   TetraCoordinates=[TetraCoordinates;midP2];
   %Second Midpoint ID
   [row column] =size(TetraCoordinates);
   SecMidPID=row;
   end    
   
  %Composing The New Three Tetrahedra by Trisection of Longest Edge
   
   % Composing New Tetrahedron 1
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j1);
   
   %Finding Vertex 3 and 4 ID
     VertexIDs =TetraDT(x,:);
     [a b] = find(VertexIDs ~=TetraDT(x,j1)& VertexIDs ~= TetraDT(x,j2)); 
   
     Vertex3 = VertexIDs(b(1));
     Vertex4 = VertexIDs(b(2));
   
   Tetrahedron1 =[Vertex1 MidPID Vertex3 Vertex4];
   
   %Updating Data Structure with tet1
   TetraDT =[TetraDT;Tetrahedron1];
   
   %Composing New Tetrahedron 2
   %Finding Vertex 1 ID
   Vertex1=TetraDT(x,j2);
   
   Tetrahedron2 =[Vertex1 SecMidPID Vertex3 Vertex4];
     %Updating Data Structure with tet2
   TetraDT =[TetraDT;Tetrahedron2];
   
   
   %Composing New Tetrahedron 3
    Tetrahedron3 =[MidPID SecMidPID Vertex3 Vertex4];
   %Updating Data Structure with tet3
   TetraDT =[TetraDT;Tetrahedron3];
  
     
   %Saving Longest Edge Vertexes ID for checking neighbor tetrahedra
   Vertex1ID =TetraDT(x,j1);
   Vertex2ID =TetraDT(x,j2);
   
   Edge =[Vertex1ID Vertex2ID]; %Concatenate Longest Edge Vertices ID
   SurroundingEdgeSet =[SurroundingEdgeSet;Edge];
   
   
 end
 
   %iterate over each old selected tetrahedron updating Data Structure deleting old ones
    if (isempty(SelectedTetraIndex)==0)
     TetraDT(SelectedTetraIndex(:),:)=[];
    end
     
           
  end       
 
  
tElapsed=toc; %stop timer


tetramesh(TetraDT,TetraCoordinates,'FaceColor','none'); % Displays each tetrahedron defined as a mesh with removed face color

%Updating refine level and mean quality for quality evolution graph
global refineIteration;%Number of refinement steps
global meanQualityValues; %Quality mean values
global meanQualityValues2;
global meanQualityValues3;
global meanQualityValues4;

%Increasing Refine Level variable
[row iterationColumn]=size(refineIteration);
iterationNumber =refineIteration(1,iterationColumn)+1;
refineIteration=[refineIteration iterationNumber];

%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];


meanValue=mean(quality2); %Mean Value of Quality 2

%Concatenating Quality 2 Mean Value
meanQualityValues2=[meanQualityValues2 meanValue];

meanValue=mean(quality3); %Mean Value of Quality 3

%Concatenating Quality 3 Mean Value
meanQualityValues3=[meanQualityValues3 meanValue];

meanValue=mean(quality4); %Mean Value of Quality 4

%Concatenating Quality 4 Mean Value
meanQualityValues4=[meanQualityValues4 meanValue];

%Updating GUI after refinement algorithm
% Visualization of # Tetrahedra
[Tetracount vertexNumber]= size(TetraDT); %number of element
text = ['Tetrahedra: ' int2str(Tetracount)]; %String concatenation
set(handles.tetLabel,'String',text); % Update #Tetrahedra info

% Visualization of # Vertices
[Vertexcount vertexNumber]= size(TetraCoordinates); %number of vertices
text = ['Vertices: ' int2str(Vertexcount)]; %String Concatenation
set(handles.vertLabel,'String',text); % Update #Vertices info

%Visualization of Time in seconds
text = ['Time Elapsed: ' num2str(tElapsed) ' Secs']; %String Concatenation
set(handles.timeLabel,'String',text); % Update Time info
    
else
 msgbox('Wrong Quality Value Input','Error Window','error');    
end    

end


% --------------------------------------------------------------------
function Untitled_46_Callback(hObject, eventdata, handles)
%Show info about the software
msgbox('Mesh-i. A 3D meshing software for Quality Assessment in tetrahedral meshes','Software Info','help');    
