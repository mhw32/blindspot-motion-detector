obj=mmreader('trial1.avi');
a=read(obj);
frames=get(obj,'numberOfFrames');
d1=figure('PaperSize',[600 600],'visible','off');
d2=figure('visible','off');
scrsz = get(0,'ScreenSize');% retrieve screensize of computer
d3=figure('Position',[scrsz(3)/7 scrsz(4)/7 scrsz(3)/2 scrsz(4)/2],'visible','off');%create the position and size of figure d3
winsize = get(d3,'Position');%initiate winsize, which will determine the position and size of the movie,
winsize(1:2) = [0 0]; %adjust size of window to include whole figure window
A=moviein(frames,d3,winsize);%create movie matrix A
values=0; %intiate the values variable
close all
xx=zeros(1,74);
countit=0;
threshold=10;
for tony=1:74
    xx(1,tony)=tony;
end
yy=zeros(1,74);
zz=zeros(1,74);

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

for k= 1:frames-2%loop
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
lickit=floor(size(u(1,:),2)*1/5);%1/5 horizontal
suckit=floor(size(u(:,1),1)*1/5);%1/5 vertical
u=reshape(u,[1,uv]);
v=reshape(v,[1,uv]);




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
%than threshold
for trikk=1:uv
    if sqrt(u(1,trikk).^2+v(1,trikk).^2)>threshold 
        tobes(1,trikk)= sqrt(u(1,trikk).^2+v(1,trikk).^2); % find greater than threshold
    end
end
joes=atan2(-v,u);
arrowlengthsquared = u.^2 + v.^2;
NorthWestVectors = find(tobes~=0 & joes>pi/2 & joes<pi   & Just==2);% looks for NorthWest vectors that are greater than threshold
NorthEastVectors = find(tobes~=0 & joes>0 & joes<pi/2   & Just==2);% looks for NorthEast vectors that are greater than threshold
SouthWestVectors = find(tobes~=0 & joes>-1*pi & joes<-1*pi/2 & Just==2);% looks for SouthWest vectors that are greater than threshold
SouthEastVectors = find(tobes~=0 & joes>-1*pi/2 & joes<0 & Just==2);% looks for SouthEast vectors that are greater than threshold


%-----start of determining box location code------------------------------------------
dork=-1;
nork=-1;
horizontalcenter=0;
verticalcenter=0;
itsthecounter=0;
if size(SouthWestVectors(1,:),2)>50%counter to make sure there are 105 arrows
    nork=0;
    bork=0;
jerk=reshape(u,[korean,thekoran]);
berk=reshape(v,[korean,thekoran]);
thewhy=reshape(joes,[korean,thekoran]);
for jimmy=1:korean
    for hendrix=1:thekoran
        if sqrt(jerk(jimmy,hendrix)^2+berk(jimmy,hendrix)^2)>threshold && thewhy(jimmy,hendrix)>-pi && thewhy(jimmy,hendrix)<-pi/2
            horizontalcenter=horizontalcenter+hendrix;
            verticalcenter=verticalcenter+jimmy;
            itsthecounter=itsthecounter+1;
        end
    end
end
anumber=floor(horizontalcenter/itsthecounter);
somenumber=floor(verticalcenter/itsthecounter);

       
%-----end of determining box location
%code--------------------------------------------
if(nork~=-1 && bork~=-1)
%--------------beginning of box code---------
Box=zeros(1,uv);

    %make vertical sides of box
    
   if anumber-lickit<1 
       anumber=lickit+1; end%what if parameters of box are out of range? default to a box hugging
    %the side of the frame
    if anumber+lickit>thekoran
        anumber=thekoran-lickit; end
        if somenumber-suckit<1
            somenumber=suckit+1;end
        if somenumber+suckit>korean
            somenumber=korean-suckit;end
            
    
for kimchi=anumber-lickit:lickit*2:anumber+lickit%horizontal 
   for kimmy=somenumber-suckit:somenumber+suckit%vertical
Box(1,korean*(kimchi-1)+ kimmy)=2;
u(1,korean*(kimchi-1)+ kimmy)=2;
v(1,korean*(kimchi-1)+ kimmy)=2;

   end
end
for kimchi=anumber-lickit:anumber+lickit%horizontal 
   for kimmy=somenumber-suckit:suckit*2:somenumber+suckit%vertical
Box(1,korean*(kimchi-1)+ kimmy)=2;
u(1,korean*(kimchi-1)+ kimmy)=2;
v(1,korean*(kimchi-1)+ kimmy)=2;

   end
end

%----------------end of box code--------------
testvectors=find(   Box~=0  );
end%end of statement if a box radius was found
end


%lines 58-69: used to determine which vectors are NW, NE,SW,SE
quiver(x(NorthWestVectors), y(NorthWestVectors),u(NorthWestVectors),v(NorthWestVectors),0, 'color', 'r', 'linewidth', 2);%
%line 88 shows blue northwest vectors
%the 0 standardizes the arrow length.

quiver(x(NorthEastVectors), y(NorthEastVectors),u(NorthEastVectors),v(NorthEastVectors),0, 'color', 'r', 'linewidth', 2);
%line 92 shows pink northeast vectors
quiver(x(SouthWestVectors), y(SouthWestVectors),u(SouthWestVectors),v(SouthWestVectors),0, 'color', 'g', 'linewidth', 2);
%line 94 shows green southwest vectors

quiver(x(SouthEastVectors), y(SouthEastVectors),u(SouthEastVectors),v(SouthEastVectors),0, 'color', 'r', 'linewidth', 2);
%line 97 shows red southeast vectors
if(nork~=-1 && bork~=-1)
quiver(x(testvectors),y(testvectors),u(testvectors),v(testvectors),0,'color','c','linewidth',2);
end
% 0 = won't autoscale to biggest arrow

im3=imcapture(d1);





%--------------------end of video(withvectors)----------------------------------------%
%-------------start construction of rose plot-----------------------------------------%




%need to create certain angle range...copy code from uh3.m

%%%--- start copied code from uh3.m


%%%----end of copied code from uh3.m




zombies=  (joes(find(tobes~=0 &  joes>-1*pi/2  & Just==2)));

zombiesmyass= (joes(find(tobes~=0 & joes>-1*pi & joes<-1*pi/2  & Just==2)));
zombiessize=size(zombies(1,:),2);
zombiesmyasssize=size(zombiesmyass(1,:),2);


%end
%-------------end construction of rose plot--------------------------------------------%

%----start construction of subplot ---------------------------------------------------%
d3=figure('visible','off')
set(d3,'Position',[scrsz(3)/3.5 scrsz(4)/3 scrsz(3)/1.75 scrsz(4)/1.75])
set(d3,'Menubar','none');

%subplot('position',[0,-0.05,0.5,0.95]);
%imshow(im5.cdata(33:367,33:367,:));
%hold on;
subplot('position',[0.5,-0.05,0.59,.88]);
imshow(im3,[]);


hax = axes('Position', [-0.08, 0, .75, .75]);



if   zombiessize*3>zombiesmyasssize | zombiessize*3==zombiesmyasssize
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
%set(p,'FaceColor','g');
%set(d2,'position',[0 0 400 400]);
%set(hax,'XTick',[])
%DF(k)=imcapture(d3);

   A(k-countit).cdata=imcapture(d3);
   %-----------------------------start of real acceleration
   %start of acceleration code
  johnny=joes;
  tally=0;
total=0; % 
desiredangle=-pi*.95;  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured

if desiredangle>3/4*pi | desiredangle<-3/4*pi  % if desiredangle>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
    %the first quartile(NE). to make these values evaluated, add 2pi to those values, as seen in line 200. 
    for yikes = 1: size(johnny(1,:),2)
        if johnny(1,yikes)<desiredangle+pi/4
            johnny(1,yikes)=johnny(1,yikes)+2*pi; 
        end
    end
end
        
for yantis=1:size(johnny(1,:),2)
    if johnny(1,yantis)<desiredangle+pi/4+2*pi % lines 206-207 make the range pi/4 units around specified direction.
        if johnny(1,yantis)>desiredangle-pi/4+2*pi 
            tally=tally+1; % tallies number of vectors in this direction
            total=total+sqrt((u(1,yantis)).^2+(v(1,yantis)).^2); % totals vector values
        end
    end
end
vectoraverage=total/tally; %average vector in certain direction
pixels=vectoraverage*9/10; %convert vectorlength to pixels
% 29 frames= 1 second
%window length (bottom portion of closest left car window) has 75 pixels, it is 81.3 cm.
actuallength=pixels*81.3/75; %converts from pixels to centimeters
speedconversion=actuallength*get(obj,'FrameRate') %converts from per frame to per second 
   yy(1,k-countit)=speedconversion;
   if k~=frames-1
       im1 = a(:,:,:,k+1);
im2=a(:,:,:,k+2);
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
for trikk=1:uv
    if sqrt(u(1,trikk).^2.*v(1,trikk).^2)>threshold 
        tobes(1,trikk)= sqrt(u(1,trikk).^2+v(1,trikk).^2);
    end
end
       %first end of acceleration
   %start of acceleration code
  
        

  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured
total=0;
tally=0;
joes=johnny;
if desiredangle>3/4*pi | desiredangle<-3/4*pi  % if desiredangle>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
    %the first quartile(NE). to make these values evaluated, add 2pi to those values, as seen in line 200. 
    for yikes = 1: size(joes(1,:),2)
        if joes(1,yikes)<desiredangle+pi/4
            joes(1,yikes)=joes(1,yikes)+2*pi; 
        end
    end
end
        
for yantis=1:size(joes(1,:),2)
    if johnny(1,yantis)<desiredangle+pi/4+2*pi % lines 206-207 make the range pi/4 units around specified direction.
        if johnny(1,yantis)>desiredangle-pi/4+2*pi 
            tally=tally+1; % tallies number of vectors in this direction
            total=total+sqrt((u(1,yantis)).^2+(v(1,yantis)).^2); % totals vector values
        end
    end
end
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


end%stop for loop

%close all
%values=values/(frames-31)

%thetrifroce=figure();
%plot(xx,yy)
%zelda=figure();
plot(xx,zz)
triforce=figure();

set(triforce,'Position',[scrsz(3)/18 scrsz(4)/17 scrsz(3)/0.88 scrsz(4)/1.15],'Menubar','none')

movie(triforce,A,30,3,winsize)

%implay(D);
%implay(G);
%implay(DF);