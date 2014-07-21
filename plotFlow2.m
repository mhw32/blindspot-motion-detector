function [tryce]=plotFlow2(u, v, imgOriginal, rSize, scale)
% Creates a quiver plot that displays the optical flow vectors on the
% original first frame (if provided). See the MATLAB Function Reference for
% "quiver" for more info.
%
% Usage:
% plotFlow(u, v, imgOriginal, rSize, scale)
%
% u and v are the horizontal and vertical optical flow vectors,
% respectively. imgOriginal, if supplied, is the first frame on which the
% flow vectors would be plotted. use an empty matrix '[]' for no image.
% rSize is the size of the region in which one vector is visible. scale
% over-rules the auto scaling.
%
% Author: Mohd Kharbat at Cranfield Defence and Security
% mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
% Published under a Creative Commons Attribution-Non-Commercial-Share Alike
% 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% October 2008
% Rev: Jan 2009

figure('visible','off');

if nargin>2
    if sum(sum(imgOriginal))~=0
        imshow(imgOriginal,[0 255]);
        hold on;
    end
end
if nargin<4
    rSize=5;
end
if nargin<5
    scale=3;
end

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
%v(1,1:101372)=0;
%v(1,10454)=226;
%v(1,34233)=300;
%v(1,74463)=342;

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
nonzeroarrows = find(sqrt(arrowlengthsquared)>15 & Why==0);
nonzeroarrowsz = find(sqrt(arrowlengthsquared)>15 & Why==300);
nonzeroarrowsd = find(sqrt(arrowlengthsquared)>15 & Why==600);
nonzeroarrowst = find(sqrt(arrowlengthsquared)>15 & Why==900);
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

quiver(x(nonzeroarrows), y(nonzeroarrows),u(nonzeroarrows),v(nonzeroarrows), 'color', 'b', 'linewidth', 2);

quiver(x(nonzeroarrowsz), y(nonzeroarrowsz),u(nonzeroarrowsz),v(nonzeroarrowsz), 'color', 'm', 'linewidth', 2);

quiver(x(nonzeroarrowsd), y(nonzeroarrowsd),u(nonzeroarrowsd),v(nonzeroarrowsd), 'color', 'g', 'linewidth', 2);

quiver(x(nonzeroarrowst), y(nonzeroarrowst),u(nonzeroarrowst),v(nonzeroarrowst), 'color', 'r', 'linewidth', 2);

%set(gca,'YDir','reverse');
 uv = size(u(:,1),1).*size(u(1,:),2);
korean=size(u(:,1),1);
thekoran=size(u(1,:),2);

x=zeros(1,uv);
y=zeros(1,uv);
for i = 1:thekoran
 for j= 1:korean
x(1,(korean.*(i-1)+j ))=i;
y(1,(korean.*(i-1)+j ))= j;
end
end

u=reshape(u,[1,uv]);
v=reshape(v,[1,uv]);
arrowlengthsquared = u.^2 + v.^2;
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
            k=v(1,bbq)/u(1,bbq);
pants(1,bbq)=2*pi-2.*sqrt(atan(k)*atan(k));

        end
    end
end


%plotting
%directions = atan2(v(nonzeroarrows),u(nonzeroarrows));
%hist(directions)
%plot(u(nonzeroarrows),v(nonzeroarrows),'.');

%point = sqrt(arrowlengthsquared(nonzeroarrows));

Zee=Z;
point = Zee(find(sqrt(arrowlengthsquared)>15))+ sqrt(arrowlengthsquared(find(sqrt(arrowlengthsquared)>15)));
[ncounts, binc]=hist(point);
%figure;
%bar(binc, ncounts);
%title('Data Analysis of Vector Magnitude and Direction');
%datacursormode on;
%axis on;
%xlabel('Direction/Magnitude of Vector: [NW(+0), NE(+300), SW(+600), SE(+900)]/[Left(small), Right(big)]');
%ylabel('Frequency of Vector');
joes=atan(v./u).*atan(v./u);
zombies= pants(find(sqrt(arrowlengthsquared)>15))+ sqrt(joes(find(sqrt(arrowlengthsquared)>15)));


steal=figure('visible','off');

rose(zombies, 50); 


tryce=getframe(steal);

title('Circular Data Analysis of Vector Magnitude and Direction');