function [ meanQuality ] =FourTBarycentric(TetraDT,TetraCoordinates,iterationNumber,refineType,qMeasure)
%(4T-Barycentric)
  meanQualityValues=[];


switch refineType

    case 1
     
    for z=1:iterationNumber    
% Uniform 4T-Barycentric Refinement Algorithm
[Tetracount vertexNumber]= size(TetraDT); %number of element

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
 [Tetracount vertexNumber]= size(TetraDT); %number of element
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
% Local 4T- Baricentric Refinement Algorithm by Vertex ID
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
% Local 4T-Baricentric Refinement Algorithm by Edge
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

