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

arrowlengthsquared = u.^2 + v.^2;

hockey=0;
counterstick=0;
Why=zeros(1,uv);        
for anthony=1:uv	
	if u(1,anthony)>0
        if v(1,anthony)>0 %Southeast
            hockey=hockey+1;
            counterstick=counterstick+sqrt((u(1,anthony)).^2+(v(1,anthony)).^2);
        end
	end
end
monkey=0;
countertrap=0;
for zelda=1:uv
    if u(1,zelda)>0
        if v(1,zelda)<0 %northeast
            monkey=monkey+1;
            Why(1,zelda)=Why(1,zelda)+300;
            countertrap=countertrap+sqrt((u(1,monkey)).^2+(v(1,monkey)).^2);
        end
    end
end

silversurfer=0;
reptile=0;
for tektite=1:uv
    if u(1,tektite)<0
        if v(1,tektite)>0 % southwest
            silversurfer=silversurfer+1;
            Why(1,tektite)=Why(1,tektite)+600;
            reptile=reptile+sqrt((u(1,tektite)).^2+(v(1,tektite)).^2);
        end
    end
end

linktothepast=0;
town=0; % vector sum
for yantis=1:uv
    if u(1,yantis)<0
        if v(1,yantis)<0 %northwest
            linktothepast=linktothepast+1;
            Why(1,yantis)=Why(1,yantis)+900;
            town=town+sqrt((u(1,yantis)).^2+(v(1,yantis)).^2);
        end
    end
end
town
linktothepast



%Try NW avg
vectoraverage=town/linktothepast %average vector in certain direction
pixels=vectoraverage*9/10
% now i need to go downstairs and grab that picture of car
% dimension...great
%actually, lets just make life easy for now,50 pixels = 1 inch
% 29 frames= 1 second
%window length (bottom portion of closest left car window) has 75 pixels, it is 81.3 cm.
actuallength=pixels*81.3/75
speedconversion=actuallength*29 % 





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
bryced=quiver(x(nonzeroarrows), y(nonzeroarrows),u(nonzeroarrows),v(nonzeroarrows), 'color', 'b', 'linewidth', 2);



gyce=quiver(x(nonzeroarrowsz), y(nonzeroarrowsz),u(nonzeroarrowsz),v(nonzeroarrowsz), 'color', 'm', 'linewidth', 2);

mice=quiver(x(nonzeroarrowsd), y(nonzeroarrowsd),u(nonzeroarrowsd),v(nonzeroarrowsd), 'color', 'g', 'linewidth', 2);

rice=quiver(x(nonzeroarrowst), y(nonzeroarrowst),u(nonzeroarrowst),v(nonzeroarrowst), 'color', 'r', 'linewidth', 2);

[trololol]=getframe;

%----------plotFlow2 copy stuff below------------------------------%





