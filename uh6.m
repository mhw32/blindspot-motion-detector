obj=mmreader('trial1.avi');
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
close all
xx=zeros(1,frames-1);
countit=0;
%threshold=30;
for tony=1:frames-1
    xx(1,tony)=tony;
end
yy=zeros(1,frames-1);
zz=zeros(1,frames-1);

im1 = a(:,:,:,1);
korean=size(im1(:,1,1),1);
thekoran=size(im1(1,:,1),2);
uv=korean*thekoran;
x=zeros(1,uv);
y=zeros(1,uv);
for i = 1:thekoran%352
 for j= 1:korean%288
     x(1,(korean.*(i-1)+j ))=i;
y(1,(korean.*(i-1)+j ))= j;
%x(1,(288.*(i-1)+j ))=i;
%y(1,(288.*(i-1)+j ))= j;
end
end
anglebisector=-pi;
for k= 33%loop
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
%figure(d1);

d1=figure('visible','off','PaperSize',[600 600]); 

set(d1,'Position',[scrsz(3)/3.5 scrsz(4)/3 scrsz(3)/1.75 scrsz(4)/1.75])
set(d1,'Menubar','none');

subplot('position',[0.58,0.05,0.39,.68]);
if sum(sum(displayImg))~=0
        imshow(displayImg,[0 255]);
        hold on;% holds image onto figure d1, so that the next image(in this case, the vector arrows) does
        %not replace the image, causing it to disappear.
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
u=reshape(u,[1,uv]);
v=reshape(v,[1,uv]);
timesaver=find(u~=0 & v~=0); %a matrice, size based on #of non-zero arrows, indices contains the indice number of non-zero indices 
arrowlengthsquareit=sqrt(u(1,timesaver).^2+v(1,timesaver).^2); %returns the vector's magnitude of non-zero indices
findthreshold=sort(arrowlengthsquareit,'descend');%sort from BIGGEST, vector magnitudes
threshold=findthreshold(1,floor(size(timesaver(1,:),2)/2));%take the 500th biggest magnitude, and set it as the THRESHOLD


%calculate means
arrowlength=sqrt(u(find(sqrt(u.^2+v.^2)>threshold)).^2+v(find(sqrt(u.^2+v.^2)>threshold)).^2);
dfdf=mean(arrowlength);

s=std(arrowlength);

Just=zeros(1,uv);
for trikk =1:uv
    if (sqrt(u(1,trikk).^2+v(1,trikk).^2)-dfdf)/s<2%if z-score is less than two(if its not greater than two standard deviation units)
        Just(1,trikk)=2;
    end
end
tobes=zeros(1,uv);%initiates tobes with matrix size of 1 x uv.
%tobes is used as a threshold factor to determine which vectors are greater
%than threshold
for trikk=1:uv
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

CarVectors=find(tobes~=0 & johnny>3*pi/4); %stores indices values'
BackGroundVectors=find(tobes~=0 & johnny<3*pi/4);%stores indices values'

thetally=2;
   while size(CarVectors(1,:),2)+size(BackGroundVectors(1,:),2)>5000%if there are lots of arrows
    
threshold=findthreshold(1,floor(size(timesaver(1,:),2)/(4*(thetally-1))));%take the 75%biggest and set it as the THRESHOLD

tobes=zeros(1,uv);
for trikk=1:uv
    if sqrt(u(1,trikk).^2+v(1,trikk).^2)>threshold 
        tobes(1,trikk)= 1; % find greater than threshold
    end
end


CarVectors=find(tobes~=0 & johnny>3*pi/4); %stores indices values'
BackGroundVectors=find(tobes~=0 & johnny<3*pi/4);%stores indices values'
 thetally=thetally+2  

   end

%--------------------start-determine SIZE of box-----------------
lickit=floor(Length*size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)/2);
if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.4
    lickit=floor(Length*0.4);
end

suckit=floor(Height*size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)/2);
 
if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.4
    suckit=floor(Height*0.4);
end
%-----------------------end-determine SIZE of box

%-----start of determining box location code------------------------------------------
presence=0;

horizontalcenter=0;
verticalcenter=0;
itsthecounter=0;

if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.1
    presence=1;
jerk=reshape(u,[korean,thekoran]);
berk=reshape(v,[korean,thekoran]);
thewhy=reshape(johnny,[korean,thekoran]);
for jimmy=1:korean
    for hendrix=1:thekoran
        if sqrt(jerk(jimmy,hendrix)^2+berk(jimmy,hendrix)^2)>threshold && thewhy(jimmy,hendrix)>3*pi/4
            horizontalcenter=horizontalcenter+hendrix;
            verticalcenter=verticalcenter+jimmy;
            itsthecounter=itsthecounter+1;
        end
    end
end
HorizontalCenter=floor(horizontalcenter/itsthecounter);
VerticalCenter=floor(verticalcenter/itsthecounter);

       
%-----end of determining box location
%code--------------------------------------------
%---start--make box location more accurate---
NoOutliers=zeros(1,uv);
LeftBorder=HorizontalCenter-lickit*2;
RightBorder=HorizontalCenter+lickit*2;
TopBorder=VerticalCenter-suckit*2;
BottomBorder=VerticalCenter+suckit*2;
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
%--start--find points two radius units away from box
    for kimchi=LeftBorder:RightBorder%horizontal 
   for kimmy=TopBorder:BottomBorder%vertical
       NoOutliers(1,korean*(kimchi-1)+ kimmy)=1;
   end
    end
     %--end--find points two radius units away from box
%--remake box center
horizontalcenter=0;
verticalcenter=0;
itsthecounter=0;

NoOutliers=reshape(NoOutliers,[korean,thekoran]);
for jimmy=1:korean
    for hendrix=1:thekoran
        if sqrt(jerk(jimmy,hendrix)^2+berk(jimmy,hendrix)^2)>threshold && thewhy(jimmy,hendrix)>3*pi/4 &&NoOutliers(jimmy,hendrix)==1
            horizontalcenter=horizontalcenter+hendrix;
            verticalcenter=verticalcenter+jimmy;
            itsthecounter=itsthecounter+1;
        end
    end
end

HorizontalCenter=floor(horizontalcenter/itsthecounter);
VerticalCenter=floor(verticalcenter/itsthecounter);
%--remake box center

%---end---make box location more accurate---
if(presence==1 )
%--------------beginning of box code---------
Box=zeros(1,uv);

    %---start-check if box fits in picture---
   if HorizontalCenter-lickit<1 
       HorizontalCenter=lickit+1; end%what if parameters of box are out of range? default to a box hugging
    %the side of the frame
    if HorizontalCenter+lickit>thekoran
        HorizontalCenter=thekoran-lickit; end
        if VerticalCenter-suckit<1
            VerticalCenter=suckit+1;end
        if VerticalCenter+suckit>korean
            VerticalCenter=korean-suckit;end
    %---end-check if box fits in picture
    
    %---start--check if box is valid...and if car is actually present
    detectortally=0;
    for kimchi=HorizontalCenter-lickit:HorizontalCenter+lickit%horizontal 
   for kimmy=VerticalCenter-suckit:VerticalCenter+suckit%vertical
if (johnny(1,korean*(kimchi-1)+ kimmy)>3*pi/4 & tobes(1,korean*(kimchi-1)+ kimmy)~=0)
detectortally=detectortally+1;

   end
   end
    end
     %---end--check if box is valid...and if car is actually present
    Isthereabox=0;
    if detectortally>size(CarVectors(1,:),2)/10%if # car vectors in box is not even 1/10 of total car vectors, we can assume that a car isn't in the picture
        Isthereabox=1;
    end
    %make vertical sides of box
for kimchi=HorizontalCenter-lickit:lickit*2:HorizontalCenter+lickit%horizontal 
   for kimmy=VerticalCenter-suckit:VerticalCenter+suckit%vertical
Box(1,korean*(kimchi-1)+ kimmy)=2;
u(1,korean*(kimchi-1)+ kimmy)=2;
v(1,korean*(kimchi-1)+ kimmy)=2;

   end
end
for kimchi=HorizontalCenter-lickit:HorizontalCenter+lickit%horizontal 
   for kimmy=VerticalCenter-suckit:suckit*2:VerticalCenter+suckit%vertical
Box(1,korean*(kimchi-1)+ kimmy)=2;
u(1,korean*(kimchi-1)+ kimmy)=2;
v(1,korean*(kimchi-1)+ kimmy)=2;

   end
end

%----------------end of box code--------------
testvectors=find(   Box~=0  );
end%end of statement if a box radius was found
end

quiver(x(CarVectors), y(CarVectors),u(CarVectors),v(CarVectors),0, 'color', 'g', 'linewidth', 2);
quiver(x(BackGroundVectors),y(BackGroundVectors),u(BackGroundVectors),v(BackGroundVectors),0, 'color', 'r', 'linewidth', 2);
if(presence==1 && Isthereabox==1)
quiver(x(testvectors),y(testvectors),u(testvectors),v(testvectors),0,'color','c','linewidth',2);
end
% 0 = won't autoscale to biggest arrow




%--------------------end of video(withvectors)----------------------------------------%
%-------------start construction of rose plot-----------------------------------------%





zombies=  (johnny(find(tobes~=0 &  johnny<3*pi/4  & Just==2)));%background
zombiesmyass= (johnny(find(tobes~=0 & johnny>3*pi/4  & Just==2)));
zombiessize=size(zombies(1,:),2);
zombiesmyasssize=size(zombiesmyass(1,:),2);


%end
%-------------end construction of rose plot--------------------------------------------%


hax = axes('Position', [-0.073, 0, .75, .75]);



if   zombiessize*3>zombiesmyasssize || zombiessize*3==zombiesmyasssize
h=rose(zombies, 30); 
xd = get(h, 'XData') ;
yd = get(h, 'YData') ;
p = patch(xd, yd, 'r','EdgeColor',[1 1 0]) ;
hold on;
hg=rose(zombiesmyass,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'g','EdgeColor','b') ;


hold off;
elseif zombiessize*3<zombiesmyasssize
hg=rose(hax,zombiesmyass,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'g','EdgeColor','b') ;
hold on;
h=rose(hax,zombies, 30);
xd = get(h, 'XData') ;
yd = get(h, 'YData') ;
p = patch(xd, yd, 'r','EdgeColor',[1 1 0]) ;
hold off;
end
   A(k-countit).cdata=imcapture(d1);
   %-----------------------------start of real acceleration
   %start of acceleration code
  
  tally=0;
total=0; % 
  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured

NoOutliers=reshape(NoOutliers,[1,uv]);

joey=find(johnny>3*pi/4 & NoOutliers==1);%matrice containing index of car vectors in johnny
tally=size(joey(1,:),2);%# of car vectors

grrr=sqrt(u(1,joey).^2+v(1,joey).^2);
total=sum(grrr);

vectoraverage=total/tally; %average vector in certain direction
pixels=vectoraverage*9/10; %convert vectorlength to pixels
% 29 frames= 1 second
%window length (bottom portion of closest left car window) has 75 pixels, it is 81.3 cm.
actuallength=pixels*81.3/75; %converts from pixels to centimeters
speedconversion=actuallength*get(obj,'FrameRate') %converts from per frame to per second 
   yy(1,k-countit)=speedconversion;
   
   
   if k~=frames-1
       
       
    if rem(k+4,5)~=0
im1=a(:,:,:,k+1);
im2=a(:,:,:,k+2);
    else
im1=a(:,:,:,k+2);
im2=a(:,:,:,k+3);
    end
[u,v]=HS(im1,im2);

for i=1:size(u,1)
    for j=1:size(u,2)
        if floor(i/rSize)~=i/rSize || floor(j/rSize)~=j/rSize
            u(i,j)=0;
            v(i,j)=0;
        end
    end
end

u=reshape(u,[1,uv]);
v=reshape(v,[1,uv]);


timesaver=find(u~=0 & v~=0); %a matrice, size based on #of non-zero arrows, indices contains the indice number of non-zero indices 
arrowlengthsquareit=sqrt(u(1,timesaver).^2+v(1,timesaver).^2); %returns the vector's magnitude of non-zero indices
findthreshold=sort(arrowlengthsquareit,'descend');%sort from BIGGEST, vector magnitudes
threshold=findthreshold(1,500);%take the 500th biggest magnitude, and set it as the THRESHOLD


for trikk=1:uv
    if sqrt(u(1,trikk).^2.*v(1,trikk).^2)>threshold 
        tobes(1,trikk)= sqrt(u(1,trikk).^2+v(1,trikk).^2);
    end
end
%calculate means
arrowlength=sqrt(u(find(sqrt(u.^2+v.^2)>threshold)).^2+v(find(sqrt(u.^2+v.^2)>threshold)).^2);
%means=mean(sqrt(u(arrowlength).^2+v(arrowlength).^2));
%totalvectos=size(arrowlength(1,:),2);
dfdf=mean(arrowlength);

s=std(arrowlength)

Just=zeros(1,uv);
for trikk =1:uv
    if (sqrt(u(1,trikk).^2+v(1,trikk).^2)-dfdf)/s<2%if z-score is less than two(if its not greater than two standard deviation units)
        Just(1,trikk)=2;
    end
end
tobes=zeros(1,uv);%initiates tobes with matrix size of 1 x uv.
%tobes is used as a threshold factor to determine which vectors are greater
%than 15

johnny=atan2(-v,u);
if anglebisector>7/4*pi || anglebisector<-3/4*pi  % if anglebisector>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
      for yikes = 1: size(johnny(1,:),2)
        if johnny(1,yikes)<anglebisector+pi/4
            johnny(1,yikes)=johnny(1,yikes)+2*pi; %transposes interval
        end
    end
end

total=0;
tally=0;


joey=find(johnny>3*pi/4);%matrice containing index of car vectors in johnny
tally=size(joey(1,:),2);%# of car vectors

grrr=sqrt(u(1,joey).^2+v(1,joey).^2);
total=sum(grrr);

vectoraverage=total/tally; %average vector in certain direction
pixels=vectoraverage*9/10; %convert vectorlength to pixels
% now i need to go downstairs and grab that picture of car
% dimension...great
%actually, lets just make life easy for now,50 pixels = 1 inch
% 29 frames= 1 second
%window length (bottom portion of closest left car window) has 75 pixels, it is 81.3 cm.
actuallength=pixels*81.3/75; %converts from pixels to centimeters
secondpicturespeedconversion=actuallength*get(obj,'FrameRate') %converts from per frame to per second 
acceleration=(secondpicturespeedconversion-speedconversion)*get(obj,'FrameRate')
zz(1,k-countit)=acceleration;
   end%stop acceleration loop
capturescreen
end%stop for loop

%close all
%values=values/(frames-31)

%thetrifroce=figure();
%plot(xx,yy)
%zelda=figure();
%plot(xx,zz)
%triforce=figure();
%set(triforce,'Position',[scrsz(3)/16 scrsz(4)/18 scrsz(3)/0.82 scrsz(4)/1.1],'Menubar','none')
%movie(triforce,A,30,3,winsize)

%implay(D);
%implay(G);
%implay(DF);


if k>9000
triforce=figure();
scrsz = get(0,'ScreenSize');% retrieve screensize of computer
set(triforce,'Position',[scrsz(3)/16 scrsz(4)/18 scrsz(3)/0.82 scrsz(4)/1.1],'Menubar','none')
d3=figure('Position',[scrsz(3)/7 scrsz(4)/7 scrsz(3)/2 scrsz(4)/2]);%create the position and size of figure d3
winsize = get(d3,'Position');%initiate winsize, which will determine the position and size of the movie,
winsize(1:2) = [0 0]; %adjust size of window to include whole figure window
movie(triforce,A,30,3,winsize)

end