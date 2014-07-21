im1=imread('fr 3_15.tif');
im2=imread('fr 4_15.tif');
subplot(121)
imshow(im1);
subplot(122);
imshow(im2);


%optical flow parameters
alpha =1;
ite = 100;
uInitial=0;
vInitial=0;
displayFlow=1;
displayImg=im1;

[u, v] = HS(im1(1:181,1:321,:), im2, alpha, ite, uInitial, vInitial, displayFlow, displayImg);