function [ meanQuality ] = EightTLE(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure)
%Longest-Edge Refinement Algorithm (8T-LE)
  meanQualityValues=[];


switch refineType

    case 1
     
    for z=1:iterationNumber    
% 8T-LE Uniform Refinement Algorithm
[Tetracount vertexNumber]= size(TetraDT); %number of element
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
  

switch qMeasure

    case 1
       
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
       
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue]; 
    case 2
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
      
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
       
end

meanValue=mean(quality2); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

    case 3
    
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
             
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
               
end

meanValue=mean(quality3); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];
         
    case 4
          
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality4); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

end

    end
    case 2
        
    for z=1:iterationNumber    
%8T-LE By Value Refinement Algorithm
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
[Tetracount vertexNumber]= size(TetraDT); %number of element
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
 
 
switch qMeasure

    case 1
       
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
       
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue]; 
    case 2
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
      
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
       
end

meanValue=mean(quality2); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

    case 3
    
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
             
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
               
end

meanValue=mean(quality3); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];
         
    case 4
          
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality4); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

end

 

else
 msgbox('Wrong Value Input','Error Window','error');    
end    

end
 
    end


    case 3

for z=1:iterationNumber        
%Vertex ID 8T-LE Uniform Refinement Algorithm
%Input Vertex ID by user
answer = inputdlg({'All Tetrahedra attach to Vertex ID will be refine:'},'Input Vertex ID');
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
 

switch qMeasure

    case 1
       
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
       
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue]; 
    case 2
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
      
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
       
end

meanValue=mean(quality2); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

    case 3
    
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
             
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
               
end

meanValue=mean(quality3); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];
         
    case 4
          
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality4); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

end
    
 
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end

end
        
    case 4
      
   for z=1:iterationNumber  
% 8T-LE By Edge
%Input Edge by Vertex1 ID and Vertex2 ID
answer = inputdlg({'Enter Vertex1 ID:','Enter Vertex2 ID:'},'Input Edge by Vertex ID');
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
  
switch qMeasure

    case 1
       
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
        quality(i)=tetrahedron_quality3_3d(tet); %QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
       
end

meanValue=mean(quality); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue]; 
    case 2
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
       
      
        quality2(i)=tetrahedron_quality2_3d(tet); %QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
       
end

meanValue=mean(quality2); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

    case 3
    
    %Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
             
        quality3(i)=tetrahedron_quality1_3d(tet); %3.0 times the ratio of the radius of the inscribed sphere divided by that of the circumscribed sphere.
               
end

meanValue=mean(quality3); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];
         
    case 4
          
%Calculating Mean Quality Value

[Tetracount vertexNumber]= size(TetraDT); %number of element

for i=1:Tetracount
        tet= [TetraCoordinates(TetraDT(i,1),:)', TetraCoordinates(TetraDT(i,2),:)', TetraCoordinates(TetraDT(i,3),:)', TetraCoordinates(TetraDT(i,4),:)'];
        quality4(i)=tetrahedron_quality4_3d(tet); %sine of half the minimum of the four solid angles.
        
end

meanValue=mean(quality4); %Mean Value of Quality

%Concatenating Quality Mean Value
meanQualityValues=[meanQualityValues meanValue];

end
   
else
 msgbox('Wrong Vertex ID Input','Error Window','error');    
end    

end
      
   
    end           
end    
  
 meanQuality=meanQualityValues; 



end

