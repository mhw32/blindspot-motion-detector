obj=mmreader('Pcar 7_good (2).avi');
a=read(obj);
frames=get(obj,'numberOfFrames');
%for k = 1:frames-1
%I(k).cdata = a(:,:,:,k);
%I(k).colormap = [];
%j=k+1;
%I(j).cdata=a(:,:,:,k);
%I(j).colormap=[];
%end

values=0;
firstframe=30;
x=zeros(1,frames-1);
for tony=1:frames-1
    x(1,tony)=tony;
end
y=zeros(1,frames-1);
z=zeros(1,frames-1);
for k= firstframe:frames-1
    
    
%im1=imread('F13.tif');
 %  im2=imread('F14.tif');
im1 = a(:,:,:,k);
im2=a(:,:,:,k+1);
[u,v]=HS(im1,im2);

displayImg=im1;
%[tryce]=plotFlow2(u, v, displayImg, 5, 5);

[trololol,speedconversion]=plotFlow8(u,v,displayImg,5,5);
y(1,k)=speedconversion;
values=values+speedconversion
%z(1,k)=values
%-------------------
%
%

%[trololol]=plotFlow5(u,v,displayImg,5,5);


 
 


%tree=figure();



%set(clock,'Position', [0.15 0.15 0.77 0.77])

%subplot('position',[-0.1,0.2,0.7,0.7])
%tike=imshow(tryce.cdata,[])
%ok=subplot('position', [0.6,0.6,0.4,0.4])
%sike=imshow(trololol.cdata,[])

%set(ok,'Position',[0.5 0.5 0.5 0.5])
%A(k)=getframe(tree)

 %subplot(2,1,1), plot(geez)
%subplot(2,1,2), plot(tryce)
end
close all
values=values/(frames-31)

plot(x,y)

%implay(F)

%implay(D);
%implay(G);
%implay(A);

% h1 =figure 1
%loop
%plot rose in h1
% im1=getframe(h1)
%im1 trim function=im1(height, width, sw: ew)


