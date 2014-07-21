obj=mmreader('Pcar 7_good (2).avi');
a=read(obj);
frames=get(obj,'numberOfFrames');
d1=figure('PaperSize',[600 600]);
d2=figure();
scrsz = get(0,'ScreenSize');% retrieve screensize of computer
d3=figure('Position',[1 scrsz(4)/2 scrsz(3)/1.75 scrsz(4)/1.75]);%create the position and size of figure d3

winsize = get(d3,'Position');%initiate winsize, which will determine the position and size of the movie,
winsize(1:2) = [0 0]; %adjust size of window to include whole figure window
A=moviein(frames,d3,winsize);%create movie matrix A
values=0; %intiate the values variable
firstframe=30;%just an integer...not used in this code
close all
fudge=-5;
d3=figure('Resize','off');
xx=zeros(1,frames-1);
for tony=1:frames-1
    xx(1,tony)=tony;
end
yy=zeros(1,frames-1);
zz=zeros(1,frames-1);
for k= 78%loop
    
%     
im1 = a(:,:,:,k);
im2=a(:,:,:,k+1);
[u,v]=HS(im1,im2);

displayImg=im1;
%--------------create first image (video with vectors)--------------------------------------------
figure(d1);
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
 uv = size(u(:,1),1).*size(u(1,:),2);
korean=size(u(:,1),1);%vertical
thekoran=size(u(1,:),2);%horizontal
koreans=floor(size(u(:,1),1)*3/5); %vertical side
thekorans=floor(size(u(1,:),2)*3/5);%thekorans is the longer horizontal side
kolea=floor(size(u(1,:),2)*4/5);
lickit=floor(size(u(1,:),2)*1/5);%1/5 horizontal
elbow=floor(size(u(1,:),2)*4/5); %4/5 horizontal
suckit=floor(size(u(:,1),1)*1/5);%1/5 vertical
backpack=floor(size(u(:,1),1)*4/5);%4/5 vertical
u=reshape(u,[1,uv]);
v=reshape(v,[1,uv]);


Why=zeros(1,uv);        
for anthony=1:uv	
	if u(1,anthony)>0
Why(1,anthony)=300;
	end
end

for kim=1:uv
	if v(1,kim)>0
Why(1,kim)=Why(1,kim)+600;
	end
end


%calculate mean
mean=0;
totalvectos=0;
for trikk=1:uv
    if sqrt(u(1,trikk).^2.*v(1,trikk).^2)>15%vectorlength has to be greater than 15
        mean=mean+sqrt(u(1,trikk).^2.*v(1,trikk).^2);% add vectorlength
        totalvectos=totalvectos+1;%count vectorlength
        %IMPORTANT NOTE: [sqrt(u(1,trikk).^2.*v(1,trikk).^2)>15] is the vectorlength
    end
end

actualmean=mean/totalvectos;

sun=0;
for trikk=1:uv
    if sqrt(u(1,trikk).^2.*v(1,trikk).^2)>15
        sun=sun+(sqrt(u(1,trikk).^2.*v(1,trikk).^2)-actualmean)^2;% vectorlength subract average vectorlength
    end
end
standarddeviation=0;

standarddeviation=sqrt(sun/(totalvectos-1));

Just=zeros(1,uv);
for trikk =1:uv
    if (sqrt(u(1,trikk).^2.*v(1,trikk).^2)-actualmean)/standarddeviation<1%if z-score is less than two(if its not greater than two standard deviation units)
        Just(1,trikk)=2;
    end
end
%for pith=1:101376
 %   if Why(1,pith)==600
  %      if sqrt(u(1,pith).^2.*v(1,pith).^2)>15
   %     squat=sqrt(u(1,pith).^2.*v(1,pith).^2);
    %    end
    %end
%end
tobes=zeros(1,uv);%initiates tobes with matrix size of 1 x uv.
%tobes is used as a threshold factor to determine which vectors are greater
%than 15
for trikk=1:uv
    if sqrt(u(1,trikk).^2.*v(1,trikk).^2)>15 
        tobes(1,trikk)= sqrt(u(1,trikk).^2.*v(1,trikk).^2); % find greater than 15
    end
end

arrowlengthsquared = u.^2 + v.^2;
NorthWestVectors = find(tobes~=0 & Why==0   & Just==2);% looks for NorthWest vectors that are greater than 15
NorthEastVectors = find(tobes~=0 & Why==300 & Just==2);% looks for NorthEast vectors that are greater than 15
SouthWestVectors = find(tobes~=0 & Why==600 & Just==2);% looks for SouthWest vectors that are greater than 15
SouthEastVectors = find(tobes~=0 & Why==900 & Just==2);% looks for SouthEast vectors that are greater than 15


%-----start of determining box location code------------------------------------------
dork=-1;
nork=-1;
if size(SouthWestVectors(1,:),2)>105%counter to make sure there are 105 arrows
    
for tork=lickit:elbow%1/5 horizontal to 4/5 horizontal center
     bork=suckit;
     tttttttt=tork
     while bork<backpack%1/5 vertical to 4/5 vertical center%%goes down to right
        counter=0;
        sdsd=3
        captainplanet=0;
         for Bike=bork-suckit+1:bork+suckit%vertical
            for Trike=tork-lickit+1:tork+lickit% horizontal 
           
                if (tobes(1,(Trike-1)*288+Bike)~=0 & Why(1,(Trike-1)*288+Bike)==600 )%vectorlength has to be greater than 15
             
             counter=counter+1;
                end
            end
         end
         %starts with bork=suckit= 57, tork=lickit =70
         %lets say lickit is 70
         %lets say suckit is 57
        %%-------------------start code for efficient searching-----------
        for computer=bork-suckit+2:bork+suckit+1                   %vertical
          for laptop= tork-lickit:tork+lickit                    %horizontal
        if( Why(1,laptop*288+computer)==600 & tobes(1,laptop*288+computer)~=0)%find first green vector
             bork=computer+suckit+2;%determine new center
              captainplanet=1
        break
        end
         end
          if captainplanet==1
             break
          end
        end
        %%-------------------end code for efficient searching-------------
   if counter>105%%if a valid box is found, stop searching
        dork=bork+1;%y coordinate(vertical)
        break
   end
    
    bork=bork+1
     end
   if counter>105%if a valid box is found, stop searching
        nork=tork;%x coordinate(horizontal)
     break
  end
end

       
%-----end of determining box location
%code--------------------------------------------
if(nork~=-1 && bork~=-1)
%--------------beginning of box code---------
Box=zeros(1,uv);

    %make vertical sides of box
for kimchi=nork-lickit:lickit*2:nork+lickit%horizontal 
   for kimmy=dork-suckit:dork+suckit%vertical
Box(1,288*(kimchi)+ kimmy)=2;
u(1,288*(kimchi)+ kimmy)=2;
v(1,288*(kimchi)+ kimmy)=2;

   end
end
for kimchi=nork-lickit:nork+lickit%horizontal 
   for kimmy=dork-suckit:suckit*2:dork+suckit%vertical
Box(1,288*(kimchi)+ kimmy)=2;
u(1,288*(kimchi)+ kimmy)=2;
v(1,288*(kimchi)+ kimmy)=2;

   end
end

%----------------end of box code--------------
testvectors=find(   Box~=0  );
end%end of statement if a box radius was found
end
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

im3=getframe(d1);





%--------------------end of video(withvectors)----------------------------------------%
%-------------start construction of rose plot-----------------------------------------%

nonzeroarrows = find(sqrt(arrowlengthsquared)>15);

Z=zeros(1,uv);        
for anthony=1:uv	
	if u(1,anthony)>0
Z(1,anthony)=300;
	end
end

for kim=1:uv
	if v(1,kim)>0
Z(1,kim)=Z(1,kim)+600;
	end
end

for gee = 1:uv
	if u(1,gee)==0
u(1,gee)=1/(9*10^99);
end
end
anthonykim=sqrt((atan(v/u))*(atan(v/u)));

pants=zeros(1,uv);
for barbeque=1:uv % (NW)
    if u(1,barbeque)<0; % if direction is west
        if v(1,barbeque)<0; % and if direction is going north
            d=v(1,barbeque)/u(1,barbeque);
pants(1,barbeque)=2.*(0.5*pi-sqrt(atan(d)*atan(d)));
        end
    end
end

for barbequesteak=1:uv % (SW)
    if u(1,barbequesteak)<0; % if direction is west
        if v(1,barbequesteak)>0;
% and if direction is going south
pants(1,barbequesteak)=pi;
        end
    end
end
 
for barbequesteak=1:uv % (W)
    if u(1,barbequesteak)<0; % if direction is west
        if v(1,barbequesteak)==0;
 
pants(1,barbequesteak)=pi;
        end
    end
end
    

for bbq=1:uv   % (SE)
    if u(1,bbq)>0; % if direction is east
        if v(1,bbq)>0; % and if direction is going south
            
pants(1,bbq)=2*pi-2.*sqrt(atan(v(1,bbq)/u(1,bbq))*atan(v(1,bbq)/u(1,bbq)));

        end
    end
end



Zee=Z;
%point = Zee(find(tobes~=0 & Just==2))+ sqrt(arrowlengthsquared(find(tobes~=0&Just==2));
%[ncounts, binc]=hist(point);
joes=atan(v./u).*atan(v./u);
zombies= pants(find(tobes~=0 & Why==600 & Just==2))+ sqrt(joes(find(tobes~=0 & Why==600 & Just==2)));

zombiesmyass= pants(find(tobes~=0 & Why~=600 & Just==2))+ sqrt(joes(find(tobes~=0 & Why~=600 & Just==2)));
zombiessize=size(zombies(1,:),2);
zombiesmyasssize=size(zombiesmyass(1,:),2);


%end
%-------------end construction of rose plot--------------------------------------------%

%----start construction of subplot ---------------------------------------------------%
figure(d3);


%subplot('position',[0,-0.05,0.5,0.95]);
%imshow(im5.cdata(33:367,33:367,:));
%hold on;
subplot('position',[0.5,0.04,0.59,.93]);
imshow(im3.cdata,[]);


hax = axes('Position', [-0.08, 0, .75, .75]);
if   zombiessize*3>zombiesmyasssize | zombiessize*3==zombiesmyasssize
h=rose(zombies, 30); 
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'g') ;
hold on;
hg=rose(zombiesmyass,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'r') ;
hold off;
elseif zombiessize*3<zombiesmyasssize

 hg=rose(hax,zombiesmyass,30);
xy = get(hg, 'XData') ;
yz = get(hg, 'YData') ;
pt = patch(xy, yz, 'r') ;
hold on;
h=rose(hax,zombies, 30);
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'g') ;
hold off;
end
%set(p,'FaceColor','g');
%set(d2,'position',[0 0 400 400]);
%set(hax,'XTick',[])
DF(k)=getframe(d3);

   A(:,k)=getframe(d3,winsize);
   %-----------------------------start of real acceleration
   %start of acceleration code
   johnny=pants+joes;

tally=0;
total=0; % 
desiredangle=pi*1.05;  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured

if desiredangle>7*pi/4  % if desiredangle>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
    %the first quartile(NE). to make these values evaluated, add 2pi to those values, as seen in line 200. 
    for yikes = 1: size(johnny(1,:),2)
        if johnny(1,yikes)<1-desiredangle
            johnny(1,yikes)=johnny(1,yikes)+2*pi; 
        end
    end
end
        
for yantis=1:size(johnny(1,:),2)
    if johnny(1,yantis)<desiredangle+pi/4 % lines 206-207 make the range pi/4 units around specified direction.
        if johnny(1,yantis)>desiredangle-pi/4 
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
speedconversion=actuallength*29 %converts from per frame to per second 
   yy(1,k)=speedconversion;
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

Why=zeros(1,uv);        
for anthony=1:uv	
	if u(1,anthony)>0
Why(1,anthony)=300;
	end
end

for kim=1:uv
	if v(1,kim)>0
Why(1,kim)=Why(1,kim)+600;
	end
end

Just=zeros(1,uv);
for trikk =1:uv
    if (sqrt(u(1,trikk).^2.*v(1,trikk).^2)-actualmean)/standarddeviation<1%if z-score is less than two(if its not greater than two standard deviation units)
        Just(1,trikk)=2;
    end
end

tobes=zeros(1,uv);%initiates tobes with matrix size of 1 x uv.
%tobes is used as a threshold factor to determine which vectors are greater
%than 15
for trikk=1:uv
    if sqrt(u(1,trikk).^2.*v(1,trikk).^2)>15 
        tobes(1,trikk)= sqrt(u(1,trikk).^2.*v(1,trikk).^2);
    end
end
       %first end of acceleration
   %start of acceleration code
   johnny=pants+joes;


%Why=zeros(1,uv);        

tally=0;
total=0; % 
desiredangle=pi*1.05;  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured

if desiredangle>7*pi/4  % if desiredangle>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
    %the first quartile(NE). to make these values evaluated, add 2pi to those values, as seen in line 200. 
    for yikes = 1: size(johnny(1,:),2)
        if johnny(1,yikes)<1-desiredangle
            johnny(1,yikes)=johnny(1,yikes)+2*pi; 
        end
    end
end
        
for yantis=1:size(johnny(1,:),2)
    if johnny(1,yantis)<desiredangle+pi/4 % lines 206-207 make the range pi/4 units around specified direction.
        if johnny(1,yantis)>desiredangle-pi/4 
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
secondpicturespeedconversion=actuallength*29 %converts from per frame to per second 
acceleration=(secondpicturespeedconversion-speedconversion)*29
zz(1,k)=acceleration;
   end%stop acceleration loop



end%stop for loop

%close all
%values=values/(frames-31)
%triforce=figure('Position',[1 scrsz(4)/2 scrsz(3)/1.75 scrsz(4)/1.75])
%thetrifroce=figure();
%plot(xx,yy)
%zelda=figure();
%plot(xx,zz)
%triforce=figure();

%movie(triforce,A,30,3,winsize)

%implay(D);
%implay(G);
%implay(DF);