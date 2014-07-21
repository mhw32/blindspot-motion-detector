aviobj = avifile('c:\video_out.avi', 'compression', 'divx', 'fps', 30);

for i = 1:numframes
 %get a video frame
 mov = aviread('c:\video.avi',i);
 %get image data
 im = mov.cdata;

 %process image
 [...add code here...]

 %output results
 imshow(im_processed,[]);
 %save movie
 aviobj = addframe(aviobj, im_processed);

end