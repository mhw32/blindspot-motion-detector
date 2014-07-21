function plotFlow2(u, v, imgOriginal, rSize, scale)
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

figure();

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
 

%quiver(u, v, scale, 'color', 'b', 'linewidth', 2);
u=reshape(u,[1,101376]);
%u(1,1:10000)=15;
%u(1,1:101372)=0;
v=reshape(v,[1,101376]);
%v(1,1:101372)=0;
%v(1,10454)=226;
%v(1,34233)=300;
%v(1,74463)=342;


arrowlengthsquared = u.^2 + v.^2;
nonzeroarrows = find(sqrt(arrowlengthsquared)>15);
%zeroarrows=find(sqrt(arrowlengthsquared)==0);
x=zeros(1,101376);
y=zeros(1,101376);
for i = 1:352
 for j= 1:288
x(1,(288.*(i-1)+j ))=i;
y(1,(288.*(i-1)+j ))= j;
end
end

quiver(x(nonzeroarrows), y(nonzeroarrows),u(nonzeroarrows),v(nonzeroarrows), 'color', 'b', 'linewidth', 2);
%quiver(x(zeroarrows), y(zeroarrows),u(zeroarrows),v(zeroarrows), 'color', 'k', 'linewidth', 2);
set(gca,'YDir','reverse');