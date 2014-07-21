obj=mmreader('Trial1.avi');
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
%close all
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
anglebisector=-pi;
accelerationcounter=0;
for k=44%loop
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
u=reshape(u,[1,Area]);
v=reshape(v,[1,Area]);
timesaver=find(u~=0 & v~=0); %a matrice, size based on #of non-zero arrows, indices contains the indice number of non-zero indices 
arrowlengthsquareit=sqrt(u(1,timesaver).^2+v(1,timesaver).^2); %returns the vector's magnitude of non-zero indices
findthreshold=sort(arrowlengthsquareit,'descend');%sort from BIGGEST, vector magnitudes
threshold=findthreshold(1,floor(size(timesaver(1,:),2)/2));%take the value at the 50th percentile for all possible vectors, and sets it as THRESHOLD

%calculate means

tobes=zeros(1,Area);%initiates tobes with matrix size of 1 x Area.
%tobes is used as a threshold factor to determine which vectors are greater
%than threshold
for trikk=1:Area
    if sqrt(u(1,trikk).^2+v(1,trikk).^2)>threshold 
        tobes(1,trikk)= 1; % find greater than threshold
    end
end
joes=atan2(-v,u);%contains the angle of the vectors from the interval [-pi, pi]
johnny=joes;%johnny contains the angle of the vectors, but the interval can be changed.
  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured

if anglebisector>7/4*pi || anglebisector<-3/4*pi  % if anglebisector>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
    %the first quartile(NE). to make these values evaluated, add 2pi to those values, as seen in line 200. 
    for yikes = 1: size(johnny(1,:),2)
        if johnny(1,yikes)<anglebisector+pi/4
            johnny(1,yikes)=johnny(1,yikes)+2*pi; %transposes interval
        end
    end
end
carvectors=find(tobes~=0 & johnny>pi & johnny<9*pi/8);
CarVectors=find(tobes~=0 & johnny>3*pi/4 & johnny<7*pi/8); %stores indices values'
carvectorss=find(tobes~=0 & johnny>7*pi/8 & johnny<pi);
carvectorsss=find(tobes~=0 & johnny>9*pi/8 & johnny<5*pi/4);
BackGroundVectors=find(tobes~=0 & johnny<3*pi/4);%stores indices values'

quiver(x(carvectorss), y(carvectorss),u(carvectorss),v(carvectorss),0, 'color', 'k', 'linewidth', 2);

quiver(x(carvectorsss), y(carvectorsss),u(carvectorsss),v(carvectorsss),0, 'color', 'w', 'linewidth', 2);
quiver(x(carvectors), y(carvectors),u(carvectors),v(carvectors),0, 'color', 'b', 'linewidth', 2);

quiver(x(CarVectors), y(CarVectors),u(CarVectors),v(CarVectors),0, 'color', 'g', 'linewidth', 2);
quiver(x(BackGroundVectors),y(BackGroundVectors),u(BackGroundVectors),v(BackGroundVectors),0, 'color', 'r', 'linewidth', 2);
% 0 = won't autoscale to biggest arrow




%--------------------end of video(withvectors)----------------------------------------%
%-------------start construction of rose plot-----------------------------------------%




%background
backgroundfrequency=  (johnny(find(tobes~=0 &  johnny<3*pi/4  )));%background
backgroundfrequencysize=size(backgroundfrequency(1,:),2);
carfrequency= (johnny(find(tobes~=0 & johnny>3*pi/4  & johnny<7*pi/8)));
carfrequencysize=size(carfrequency(1,:),2);
carfrequencytwo= (johnny(find(tobes~=0 & johnny>pi  & johnny<5*pi/4)));
carfrequencysizetwo=size(carfrequencytwo(1,:),2);
carfrequencythree=(johnny(find(tobes~=0 & johnny>7*pi/8 & johnny<pi)));
carfrequencysizethree=size(carfrequencythree(1,:),2);
carfrequencyfour=(johnny(find(tobes~=0 & johnny>9*pi/8 & johnny<5*pi/4)));
carfrequencysizefour=size(carfrequencythree(1,:),2);

hax = axes('Position', [-0.073, 0, .75, .75]);



if   backgroundfrequencysize*3>carfrequencysize || backgroundfrequencysize*3==carfrequencysize
h=rose(backgroundfrequency, 30); 
xd = get(h, 'XData') ;
yd = get(h, 'YData') ;
p = patch(xd, yd, 'r') ;
hold on;
hg=rose(carfrequency,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'g') ;
hold on;
hg=rose(carfrequencytwo,30);
xyz = get(hg, 'XData') ;
yzz = get(hg, 'YData') ;
ptt = patch(xyz, yzz, 'b','EdgeColor','b') ;
hold on;
hg=rose(carfrequencythree,30);
xyz = get(hg, 'XData') ;
yzz = get(hg, 'YData') ;
ptt = patch(xyz, yzz, 'k','EdgeColor','b') ;
hold on;
hg=rose(carfrequencyfour,30);
xyz = get(hg, 'XData') ;
yzz = get(hg, 'YData') ;
ptt = patch(xyz, yzz, 'w') ;

hold off;
elseif backgroundfrequencysize*3<carfrequencysize
hg=rose(hax,carfrequency,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'g') ;
hold on;
h=rose(hax,backgroundfrequency, 30);
xd = get(h, 'XData') ;
yd = get(h, 'YData') ;
p = patch(xd, yd, 'r') ;
hold on;

hg=rose(carfrequencytwo,30);
xyz = get(hg, 'XData') ;
yzz = get(hg, 'YData') ;
ptt = patch(xyz, yzz, 'b','EdgeColor','b') ;
hold on;
hg=rose(carfrequencytwo,30);
xyz = get(hg, 'XData') ;
yzz = get(hg, 'YData') ;
ptt = patch(xyz, yzz, 'k','EdgeColor','b') ;
hold on;
hg=rose(carfrequencyfour,30);
xyz = get(hg, 'XData') ;
yzz = get(hg, 'YData') ;
ptt = patch(xyz, yzz, 'w') ;
hold off;
end
%---end construction of rose plot===


capturescreen
end

%stop for loop
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



