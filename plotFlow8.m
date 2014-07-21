function [trololol,speedconversion] = plotFlow6(u, v, imgOriginal, rSize, scale )
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
 uv = size(u(:,1),1).*size(u(1,:),2); % initialize uv, which is the product of the horizontal  and vertical length of the image pixels
HorizontalLength=size(u(:,1),1); %initialize HorizontalLength, which is the value of the horizontal length pixels
VerticalLength=size(u(1,:),2);%intialize VerticalLength, which is the value of the vertical length pixels


u=reshape(u,[1,uv]);% u is now a 2-D matrix that has one row, and uv columns
v=reshape(v,[1,uv]);% v is now a 2-D matrix that has one row, and uv columns

Why=zeros(1,uv);        %initialize Why, which is used to distinguish vectors in the picture with color
for anthony=1:uv	
	if u(1,anthony)>0 % if the u vector points toward east
Why(1,anthony)=300;% add 300 to 'Why' for that east vector, if it points toward west, nothing happens
	end
end

for kim=1:uv
	if v(1,kim)>0 %if the v vector points toward south
Why(1,kim)=Why(1,kim)+600;% add 600 to 'Why' for that south vector, if it points north, nothing happens
	end
end
%all SouthEast vectors should have 600+300=900 value.
%all SouthWest vectors should have 600+0=600 value.
%all NorthEast vectors should have 0+300=300 value.
%all NorthWest vectors should have 0+0=0 value.


arrowlengthsquared = u.^2 + v.^2;
NorthWestVectors = find(sqrt(arrowlengthsquared)>15 & Why==0);% looks for NorthWest vectors that are greater than 15
NorthEastVectors = find(sqrt(arrowlengthsquared)>15 & Why==300);% looks for NorthEast vectors that are greater than 15
SouthWestVectors = find(sqrt(arrowlengthsquared)>15 & Why==600);% looks for SouthWest vectors that are greater than 15
SouthEastVectors = find(sqrt(arrowlengthsquared)>15 & Why==900);% looks for SouthEast vectors that are greater than 15
x=zeros(1,uv);%initialize x, creates the matrix 
y=zeros(1,uv);%intialize y, creates the matrix
for i = 1:VerticalLength
 for j= 1:HorizontalLength
     x(1,(HorizontalLength.*(i-1)+j ))=i; %82-83 produces the coordinates of x and y for all pixels.
y(1,(HorizontalLength.*(i-1)+j ))= j;
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

[trololol]=getframe; 


uv = size(u(:,1),1).*size(u(1,:),2);
HorizontalLength=size(u(:,1),1);
VerticalLength=size(u(1,:),2);

x=zeros(1,uv);
y=zeros(1,uv);
for i = 1:VerticalLength
 for j= 1:HorizontalLength
x(1,(HorizontalLength.*(i-1)+j ))=i;
y(1,(HorizontalLength.*(i-1)+j ))= j;
end
end

u=reshape(u,[1,uv]);
v=reshape(v,[1,uv]);
arrowlengthsquared = u.^2 + v.^2;


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

joes=atan(v./u).*atan(v./u);
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






vectoraverage=total/tally %average vector in certain direction
pixels=vectoraverage*9/10 %convert vectorlength to pixels
% now i need to go downstairs and grab that picture of car
% dimension...great
%actually, lets just make life easy for now,50 pixels = 1 inch
% 29 frames= 1 second
%window length (bottom portion of closest left car window) has 75 pixels, it is 81.3 cm.
actuallength=pixels*81.3/75 %converts from pixels to centimeters
speedconversion=actuallength*29 %converts from per frame to per second 


