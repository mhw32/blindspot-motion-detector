%---uh---getting rid of redundant code like tobes and Just, and finding
%more time efficient coding in general
obj=mmreader('Trial3.avi');
a=read(obj);
frames=get(obj,'numberOfFrames');
d1=figure('PaperSize',[600 600],'visible','off');
d2=figure('visible','off');
scrsz = get(0,'ScreenSize');% retrieve screensize of computer
d3=figure('Position',[1 1 scrsz(3)/.82 scrsz(4)/1],'visible','off');%create the position and size of figure d3
winsize = get(d3,'Position');%initiate winsize, which will determine the position and size of the movie,
winsize(1:2) = [0 0]; %adjust size of window to include whole figure window
A=moviein(frames,d3,winsize);%create movie matrix A
values=0; %intiate the values variable

xx=zeros(1,frames-1);
countit=0;
for tony=1:frames-1
    xx(1,tony)=tony;
end
yy=zeros(1,frames-1);%velocity values stored here
zz=zeros(1,frames-1);%acceleration values stored here
zz(1,1)=0;

im1 = a(:,:,:,1);
VerticalLength=size(im1(:,1,1),1);
HorizontalLength=size(im1(1,:,1),2);
Area=VerticalLength*HorizontalLength;
x=zeros(1,Area);
y=zeros(1,Area);
for i = 1:HorizontalLength
 for j= 1:VerticalLength
     x(1,(VerticalLength.*(i-1)+j ))=i;
y(1,(VerticalLength.*(i-1)+j ))= j;
end
end
NoOutliers=zeros(1,Area);
anglebisector=-0.9435*pi;
InitialAngleRange=pi/4;
accelerationcounter=0;

if anglebisector<-pi+InitialAngleRange
    anglebisector=anglebisector+2*pi;
end
for k=94 %loop
    if rem(k+3,5)==0
        countit=countit+1
        continue;
    end
%     
im1 = a(:,:,:,k);
im2=a(:,:,:,k+1);
[u,v]=HS(im1,im2);
displayImg=im1;
%--------------create first image (video with vectors)--------------------------------------------


d1=figure('visible','off','PaperSize',[600 600]); 

set(d1,'Position',[scrsz(3)/3.5 scrsz(4)/3 scrsz(3)/1.75 scrsz(4)/1.75])
set(d1,'Menubar','none');

subplot('position',[0.58,0.05,0.39,.68]);
if sum(sum(displayImg))~=0
        imshow(displayImg,[0 255]);
        hold on;% holds image onto figure d1, so that the next image(in this case, the vector arrows) does
        %not replace the image
    end


    rSize=5;


  


% Enhance the quiver plot visually by showing one vector per region
for i=1:size(u,1)
    for j=1:size(u,2)
        if floor(i/rSize)~=i/rSize || floor(j/rSize)~=j/rSize
            u(i,j)=0;
            v(i,j)=0;
        end
    end
end
Length=floor(size(u(1,:),2)); %horizontal
Height=floor(size(u(:,1),1));% vertical
u=reshape(u,[1,Area]);
v=reshape(v,[1,Area]);
timesaver=find(u~=0 & v~=0); %a matrice, size based on #of non-zero arrows, indices contains the indice number of non-zero indices 
arrowlengthsquareit=sqrt(u(1,timesaver).^2+v(1,timesaver).^2); %returns the vector's magnitude of non-zero indices
findthreshold=sort(arrowlengthsquareit,'descend');%sort from BIGGEST, vector magnitudes
threshold=findthreshold(1,floor(size(timesaver(1,:),2)/2));%take the value at the 50th percentile for all possible vectors, and sets it as THRESHOLD

arrowlength=sqrt(u.^2+v.^2);
validarrows=sqrt(u(find(arrowlength>threshold)).^2+v(find(arrowlength>threshold)).^2);
AverageMagnitude=mean(validarrows);

s=std(validarrows);

johnny=atan2(-v,u);%contains the angle of the vectors from the interval [-pi, pi]
%johnny contains the angle of the vectors, but the interval can be changed.
  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured

if  anglebisector<-pi+InitialAngleRange | anglebisector>pi-InitialAngleRange % if anglebisector>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
    %the first quartile(NE). to make these values evaluated, add 2pi to those values, as seen in line 200. 
    %transporse interval if necessary
    
        if anglebisector>pi-InitialAngleRange
            anglebisector=anglebisector-2*pi;
        end
        joy=find(johnny<anglebisector+InitialAngleRange);
    for yikes = 1: size(joy(1,:),2)
            johnny(1,joy(1,yikes))=johnny(1,joy(1,yikes))+2*pi; %transposes interval
        
    end
end
if anglebisector<-pi+InitialAngleRange
    anglebisector=anglebisector+2*pi;%convert angle back to higher angle
end
UpperInitialRange=anglebisector+InitialAngleRange;
LowerInitialRange=anglebisector-InitialAngleRange;
joseph=find(arrowlength>threshold & johnny>LowerInitialRange & johnny<UpperInitialRange);% find indice values containing all the car vectors
alltheangles=johnny(1,joseph);%gets the actual angles of the car vectors
RevisedAngleRange=std(alltheangles);%get the standard deviation
CarVectors=find(arrowlength>threshold & johnny>LowerInitialRange & johnny<UpperInitialRange); %stores indices values'
BackGroundVectors=find( johnny<LowerInitialRange &arrowlength>threshold | johnny>UpperInitialRange & arrowlength>threshold);%stores indices values'

thetally=2;
   while size(CarVectors(1,:),2)+size(BackGroundVectors(1,:),2)>5000%if there are more than 5000 arrows 
    
threshold=findthreshold(1,floor(size(timesaver(1,:),2)/(4*(thetally-1))));%take the 75%biggest and set it as the THRESHOLD so that less arrows appear
CarVectors=find(arrowlength>threshold & johnny>LowerInitialRange & johnny<UpperInitialRange); %stores indices values'
BackGroundVectors=find(johnny<LowerInitialRange &arrowlength>threshold | johnny>UpperInitialRange & arrowlength>threshold);%stores indices values'
 thetally=thetally+2  
   end
   UpperRevisedRange=anglebisector+RevisedAngleRange;
   LowerRevisedRange=anglebisector-RevisedAngleRange;
carvectorstwo=find(arrowlength>threshold & johnny> LowerRevisedRange & johnny< UpperRevisedRange);
    
   %-----start of determining initial box location code------------------------------------------
presence=0;

horizontalcenter=0;
verticalcenter=0;
itsthecounter=0;

if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.1
    presence=1;
res1=reshape(u,[VerticalLength,HorizontalLength]);
res2=reshape(v,[VerticalLength,HorizontalLength]);
thewhy=reshape(johnny,[VerticalLength,HorizontalLength]);
indices=find(sqrt(res1.^2+res2.^2)>threshold & thewhy>3*pi/4);
horizontal=ceil(indices/VerticalLength);
horizontalcenter=floor(sum(horizontal)'/size(horizontal)');
vertical=ceil(rem(indices+VerticalLength,VerticalLength));
verticalcenter=floor(sum(vertical)'/size(vertical)');

       
%-----end of determining initial box location
%code--------------------------------------------
   
%--------------------start-determine SIZE of box-----------------
HorizontalParameter=floor(Length*size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)/2);
if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.4
    HorizontalParameter=floor(Length*0.4);
end

VerticalParameter=floor(Height*size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)/2);
 
if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.4
    VerticalParameter=floor(Height*0.4);
end
%-----------------------end-determine SIZE of box-----------------------



%---start--make box location more accurate ---
for standarddeviationunit=1:5
NoOutliers=zeros(1,Area);
LeftBorder=floor(horizontalcenter(:,1)-HorizontalParameter*(6-standarddeviationunit));%4 standard deviations from box
RightBorder=floor(horizontalcenter(:,1)+HorizontalParameter*(6-standarddeviationunit));
TopBorder=floor(verticalcenter(:,1)-VerticalParameter*(6-standarddeviationunit));
BottomBorder=floor(verticalcenter(:,1)+VerticalParameter*(6-standarddeviationunit));
if LeftBorder<1
    LeftBorder=1;
end
if RightBorder>Length
    RightBorder=Length;
end
if TopBorder<1
    TopBorder=1;
end
if BottomBorder>Height
    BottomBorder=Height;
end
%--start--eliminate outliers
    for kimchi=LeftBorder:RightBorder%horizontal 
   for kimmy=TopBorder:BottomBorder%vertical
       if johnny(1,VerticalLength*(kimchi-1)+kimmy)>LowerRevisedRange & johnny(1,VerticalLength*(kimchi-1)+kimmy)< UpperRevisedRange
       NoOutliers(1,VerticalLength*(kimchi-1)+ kimmy)=1;
       end
   end
    end
     %--end--eliminate outliers-----------
     
%--remake box center-------------------------------------------
horizontalcenter=0;
verticalcenter=0;
itsthecounter=0;

NoOutliers=reshape(NoOutliers,[VerticalLength,HorizontalLength]);
%calculate box center

indices=find(sqrt(res1.^2+res2.^2)>threshold & thewhy>3*pi/4 & NoOutliers==1);
if size(indices)>0
horizontal=ceil(indices/VerticalLength);
horizontalcenter=floor(sum(horizontal)'/size(horizontal)');
vertical=ceil(rem(indices+VerticalLength,VerticalLength));
verticalcenter=floor(sum(vertical)'/size(vertical)');

end
end
%--end--remake box center first time---
%---end---make box location more accurate 1---


if(presence==1 )
Box=zeros(1,Area);



    %---start-check if box fits in picture---
   if horizontalcenter-HorizontalParameter<1 
       LeftBorder=1; end%what if parameters of box are out of range? default to a box hugging
    %the side of the frame
    if horizontalcenter+HorizontalParameter>HorizontalLength
        RightBorder=HorizontalLength; end
        if verticalcenter-VerticalParameter<1
            TopBorder=1; end
        if verticalcenter+VerticalParameter>VerticalLength
            BottomBorder=VerticalLength;end
    %---end-check if box fits in picture
    
    %---start--check if box is valid...and if car is actually present
    presenceTestTwo=0;
    for kimchi=LeftBorder:RightBorder%horizontal 
   for kimmy=TopBorder:BottomBorder%vertical
if (johnny(1,VerticalLength*(kimchi-1)+kimmy)>LowerRevisedRange & ...
        johnny(1,VerticalLength*(kimchi-1)+kimmy)< UpperRevisedRange ...
        & arrowlength(1,VerticalLength*(kimchi-1)+ kimmy)>threshold)
presenceTestTwo=presenceTestTwo+1;

end
   end
    end
     %---end--check if box is valid...and if car is actually present
    Isthereabox=0;
    if presenceTestTwo>size(carvectorstwo(1,:),2)/3%if # car vectors in box is not even 1/5 of total car vectors, we can assume that a car isn't in the picture
        Isthereabox=1;
    end
    %make vertical sides of box
for kimchi=LeftBorder:RightBorder-LeftBorder:RightBorder%horizontal 
   for kimmy=TopBorder:BottomBorder%vertical
Box(1,VerticalLength*(kimchi-1)+ kimmy)=2;
u(1,VerticalLength*(kimchi-1)+ kimmy)=2;
v(1,VerticalLength*(kimchi-1)+ kimmy)=2;

   end
end
%make horizontal sides of box
for kimchi=LeftBorder:RightBorder%horizontal 
   for kimmy=TopBorder:BottomBorder-TopBorder:BottomBorder%vertical
Box(1,VerticalLength*(kimchi-1)+ kimmy)=2;
u(1,VerticalLength*(kimchi-1)+ kimmy)=2;
v(1,VerticalLength*(kimchi-1)+ kimmy)=2;

   end
end

%----------------end of box code--------------
testvectors=find(   Box~=0  );
end%end of statement if a box radius was found
end
quiver(x(carvectorstwo), y(carvectorstwo),u(carvectorstwo),v(carvectorstwo),0, 'color', 'g', 'linewidth', 2);

%quiver(x(CarVectors), y(CarVectors),u(CarVectors),v(CarVectors),0, 'color', 'g', 'linewidth', 2);
quiver(x(BackGroundVectors),y(BackGroundVectors),u(BackGroundVectors),v(BackGroundVectors),0, 'color', 'r', 'linewidth', 2);
if(presence==1 && Isthereabox==1)
%quiver(x(testvectors),y(testvectors),u(testvectors),v(testvectors),0,'color','c','linewidth',2);
end
% 0 = won't autoscale to biggest arrow




%--------------------end of video(withvectors)----------------------------------------%
%-------------start construction of rose plot-----------------------------------------%





backgroundfrequency=  (johnny(find(johnny<LowerInitialRange &arrowlength>threshold &  (arrowlength-AverageMagnitude)/s<2 ...
    | johnny>UpperInitialRange & arrowlength>threshold &  (arrowlength-AverageMagnitude)/s<2)));%background
carfrequency= (johnny(find(arrowlength>threshold &  johnny>LowerInitialRange & johnny<UpperInitialRange  &  (arrowlength-AverageMagnitude)/s<2)));
backgroundfrequencysize=size(backgroundfrequency(1,:),2);
carfrequencysize=size(carfrequency(1,:),2);




hax = axes('Position', [-0.073, 0, .75, .75]);



if   backgroundfrequencysize*3>carfrequencysize || backgroundfrequencysize*3==carfrequencysize
h=rose(backgroundfrequency, 30); 
xd = get(h, 'XData') ;
yd = get(h, 'YData') ;
p = patch(xd, yd, 'r','EdgeColor',[1 1 0]) ;
hold on;
hg=rose(carfrequency,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'g','EdgeColor','b') ;


hold off;
elseif backgroundfrequencysize*3<carfrequencysize
hg=rose(hax,carfrequency,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'g','EdgeColor','b') ;
hold on;
h=rose(hax,backgroundfrequency, 30);
xd = get(h, 'XData') ;
yd = get(h, 'YData') ;
p = patch(xd, yd, 'r','EdgeColor',[1 1 0]) ;
hold off;
end
%---end construction of rose plot===
   A(k-countit).cdata=imcapture(d1);
  tally=0;
total=0; % 

NoOutliers=reshape(NoOutliers,[1,Area]);

joey=find(johnny>3*pi/4 & NoOutliers==1);%matrice containing index of car vectors in johnny
tally=size(joey(1,:),2);%# of car vectors

grrr=arrowlength(1,joey);
total=sum(grrr);

vectoraverage=total/tally; %average vector in certain direction
if size(joey(1,:),2)==0
    vectoraverage=0;
end
pixels=vectoraverage*9/10; %convert vectorlength to pixels
%window length (bottom portion of closest left car window) has 75 pixels, it is 81.3 cm.
actuallength=pixels*81.3/75; %converts from pixels to centimeters
actuallength=actuallength/2.54;%convert from centimeters to inches
actuallength=actuallength/12;%convert from inches to feet
actuallength=actuallength/5280;%convert from feet to mile
if rem(k+2-countit,2)==0
velocitytwo=actuallength*get(obj,'FrameRate')*3600;%convert from framerate units to hours;
   yy(1,k-countit)=velocitytwo;
acceleration=(velocitytwo-velocityone)*get(obj,'FrameRate')*3600;

zz(1,k-1-countit)=acceleration;

elseif rem(k+2-countit,2)==1
     velocityone=actuallength*get(obj,'FrameRate')*3600;
        yy(1,k-countit)=velocityone;
     if k==1
continue;
     end
acceleration=(velocityone-velocitytwo)*get(obj,'FrameRate')*3600;
zz(1,k-1-countit)=acceleration;


end
capturescreen
end%stop for loop
%plotvelocity=figure();
%plot(xx,yy)

%plotacceleration=figure();
%plot(xx,zz)

%moviefigure=figure();
%set(moviefigure,'Position',[scrsz(3)/16 scrsz(4)/18 scrsz(3)/0.82 scrsz(4)/1.1],'Menubar','none')
%movie(moviefigure,A,30,3,winsize)

if k>9000

triforce=figure();
scrsz = get(0,'ScreenSize');% retrieve screensize of computer
set(triforce,'Position',[scrsz(3)/18 scrsz(4)/17 scrsz(3)/0.88 scrsz(4)/1.15],'Menubar','none')
d3=figure('Position',[scrsz(3)/7 scrsz(4)/7 scrsz(3)/2 scrsz(4)/2]);%create the position and size of figure d3
winsize = get(d3,'Position');%initiate winsize, which will determine the position and size of the movie,
winsize(1:2) = [0 0]; %adjust size of window to include whole figure window
movie(triforce,A,30,3,winsize)

end
