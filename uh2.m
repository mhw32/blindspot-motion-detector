obj=mmreader('Pcar 7_good (2).avi');
a=read(obj);
frames=get(obj,'numberOfFrames');
%d1=figure('PaperSize',[600 600]);
%d2=figure();
scrsz = get(0,'ScreenSize');
%d3=figure('Position',[1 scrsz(4)/2 scrsz(3)/1.75 scrsz(4)/1.75])
d4=figure('visible','off');
winsize = get(d4,'Position');
winsize(1:2) = [0 0];
A=moviein(frames,d4,winsize);
values=0;
firstframe=30;
A(1).cdata=[5];

for k= 1
    
    
im1 = a(:,:,:,k);
im2=a(:,:,:,k+1);
[u,v]=HS(im1,im2);

displayImg=im1;
%--------------create first image (video with vectors)--------------------------------------------
%figure(d1);
d1=figure('visible','off','PaperSize',[600 600]);
    if sum(sum(displayImg))~=0
        imshow(displayImg,[0 255]);
        hold on;
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
korean=size(u(:,1),1);
thekoran=size(u(1,:),2);


%quiver(u, v, scale, 'color', 'b', 'linewidth', 2);
u=reshape(u,[1,uv]);
%u(1,1:10000)=15;
%u(1,1:101372)=0;
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



arrowlengthsquared = u.^2 + v.^2;
NorthWestVectors = find(sqrt(arrowlengthsquared)>15 & Why==0);% looks for NorthWest vectors that are greater than 15
NorthEastVectors = find(sqrt(arrowlengthsquared)>15 & Why==300);% looks for NorthEast vectors that are greater than 15
SouthWestVectors = find(sqrt(arrowlengthsquared)>15 & Why==600);% looks for SouthWest vectors that are greater than 15
SouthEastVectors = find(sqrt(arrowlengthsquared)>15 & Why==900);% looks for SouthEast vectors that are greater than 15

%zeroarrows=find(sqrt(arrowlengthsquared)==0);
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
bryced=quiver(x(NorthWestVectors), y(NorthWestVectors),u(NorthWestVectors),v(NorthWestVectors), 'color', 'b', 'linewidth', 2);%
%line 88 shows blue northwest vectors


gyce=quiver(x(NorthEastVectors), y(NorthEastVectors),u(NorthEastVectors),v(NorthEastVectors), 'color', 'm', 'linewidth', 2);
%line 92 shows pink northeast vectors
mice=quiver(x(SouthWestVectors), y(SouthWestVectors),u(SouthWestVectors),v(SouthWestVectors), 'color', 'g', 'linewidth', 2);
%line 94 shows green southwest vectors

rice=quiver(x(SouthEastVectors), y(SouthEastVectors),u(SouthEastVectors),v(SouthEastVectors), 'color', 'r', 'linewidth', 2);
%line 97 shows red southeast vectors



im3=imcapture(d1);
%im34=imshow(im3.cdata(33:367,33:367,:));




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
point = Zee(find(sqrt(arrowlengthsquared)>15))+ sqrt(arrowlengthsquared(find(sqrt(arrowlengthsquared)>15)));
[ncounts, binc]=hist(point);
joes=atan(v./u).*atan(v./u);
zombies= pants(find(sqrt(arrowlengthsquared)>15))+ sqrt(joes(find(sqrt(arrowlengthsquared)>15)));
d2=figure('visible','off');
rose(zombies, 50); 

set(d2,'position',[0 0 400 400]);
im5=imcapture(d2);

%-------------end construction of rose plot--------------------------------------------%

%----start construction of subplot ---------------------------------------------------%

d3=figure('visible','off','Position',[1 scrsz(4)/2 scrsz(3)/1.75 scrsz(4)/1.75]);

subplot('position',[0,0,0.35,0.69]);
imshow(im3);
%imshow(im5.cdata(33:367,33:367,:));

%subplot('position',[0.2,0,0.34,.23]);
%imshow(im5,[]);
%imshow(im3.cdata,[]);
%DF(k)=imcapture(d3,'img',[300 400])
%IMG = IMCAPTURE(H,'img',[MROWS NCOLS])
   A(k).cdata=imcapture(d3)
end
close all
%values=values/(frames-31)
%triforce=figure();
%plot(x,y)
%movie(triforce,A,30,3,winsize)
%implay(F)


%implay(D);
%implay(G);
%implay(DF);