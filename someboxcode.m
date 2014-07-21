
Length=floor(size(u(1,:),2)); %horizontal
Height=floor(size(u(:,1),1));% vertical
u=reshape(u,[1,uv]);
v=reshape(v,[1,uv]);
timesaver=find(u~=0 & v~=0); %a matrice, size based on #of non-zero arrows, indices contains the indice number of non-zero indices 
arrowlengthsquareit=sqrt(u(1,timesaver).^2+v(1,timesaver).^2); %returns the vector's magnitude of non-zero indices
findthreshold=sort(arrowlengthsquareit,'descend');%sort from BIGGEST, vector magnitudes
threshold=findthreshold(1,floor(size(timesaver(1,:),2)/2));%take the 500th biggest magnitude, and set it as the THRESHOLD


%calculate means
arrowlength=sqrt(u(find(sqrt(u.^2+v.^2)>threshold)).^2+v(find(sqrt(u.^2+v.^2)>threshold)).^2);
dfdf=mean(arrowlength);

s=std(arrowlength);

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
        tobes(1,trikk)= 1; % find greater than threshold
    end
end
joes=atan2(-v,u);%contains the angle of the vectors from the interval [-pi, pi]
johnny=joes;%johnny contains the angle of the vectors, but the interval can be changed.
  %initialize variable desiredangle with value of 1.05pi. Used to determine angle direction that will be measured

if anglebisector>7/4*pi || anglebisector<-3/4*pi  % if anglebisector>7pi/4, that means values greater than 2pi will be evaluated. These values are basically
    %the first quartile(NE). to make these values evaluated, add 2pi to those values, as seen in line 200. 
    for yikes = 1: size(johnny(1,:),2)
        if johnny(1,yikes)<anglebisector+pi/4
            johnny(1,yikes)=johnny(1,yikes)+2*pi; %transposes interval
        end
    end
end

CarVectors=find(tobes~=0 & johnny>3*pi/4); %stores indices values'
BackGroundVectors=find(tobes~=0 & johnny<3*pi/4);%stores indices values'

thetally=2;
   while size(CarVectors(1,:),2)+size(BackGroundVectors(1,:),2)>5000%if there are lots of arrows
    
threshold=findthreshold(1,floor(size(timesaver(1,:),2)/(4*(thetally-1))));%take the 75%biggest and set it as the THRESHOLD

tobes=zeros(1,uv);
for trikk=1:uv
    if sqrt(u(1,trikk).^2+v(1,trikk).^2)>threshold 
        tobes(1,trikk)= 1; % find greater than threshold
    end
end


CarVectors=find(tobes~=0 & johnny>3*pi/4); %stores indices values'
BackGroundVectors=find(tobes~=0 & johnny<3*pi/4);%stores indices values'
 thetally=thetally+2  

   end

%--------------------start-determine SIZE of box-----------------
lickit=floor(Length*size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)/2);
if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.4
    lickit=floor(Length*0.4);
end

suckit=floor(Height*size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)/2);
 
if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.4
    suckit=floor(Height*0.4);
end
%-----------------------end-determine SIZE of box

%-----start of determining box location code------------------------------------------
presence=0;

horizontalcenter=0;
verticalcenter=0;
itsthecounter=0;

if size(CarVectors(1,:),2)/size(BackGroundVectors(1,:),2)>0.1
    presence=1;
res1=reshape(u,[korean,thekoran]);
res2=reshape(v,[korean,thekoran]);
thewhy=reshape(johnny,[korean,thekoran]);
for sigma=1:korean
    for gamma=1:thekoran
        if sqrt(res1(sigma,gamma)^2+res2(sigma,gamma)^2)>threshold && thewhy(sigma,gamma)>3*pi/4
            horizontalcenter=horizontalcenter+gamma;
            verticalcenter=verticalcenter+sigma;
            itsthecounter=itsthecounter+1;
        end
    end
end
HorizontalCenter=floor(horizontalcenter/itsthecounter);
VerticalCenter=floor(verticalcenter/itsthecounter);

       
%-----end of determining box location
%code--------------------------------------------
%---start--make box location more accurate---
NoOutliers=zeros(1,uv);
LeftBorder=HorizontalCenter-lickit*2;
RightBorder=HorizontalCenter+lickit*2;
TopBorder=VerticalCenter-suckit*2;
BottomBorder=VerticalCenter+suckit*2;
if LeftBorder<1
    LeftBorder=1;
end
if RightBorder>Length
    RightBorder=Length;
end
if TopBorder<1
    TopBorder=1;
end
if BottomBorder>Height
    BottomBorder=Height;
end
%--start--find points two radius units away from box
    for kimchi=LeftBorder:RightBorder%horizontal 
   for kimmy=TopBorder:BottomBorder%vertical
       NoOutliers(1,korean*(kimchi-1)+ kimmy)=1;
   end
    end
     %--end--find points two radius units away from box
%--remake box center
horizontalcenter=0;
verticalcenter=0;
itsthecounter=0;

NoOutliers=reshape(NoOutliers,[korean,thekoran]);
for sigma=1:korean
    for gamma=1:thekoran
        if sqrt(res1(sigma,gamma)^2+res2(sigma,gamma)^2)>threshold && thewhy(sigma,gamma)>3*pi/4 &&NoOutliers(sigma,gamma)==1
            horizontalcenter=horizontalcenter+gamma;
            verticalcenter=verticalcenter+sigma;
            itsthecounter=itsthecounter+1;
        end
    end
end
if itsthecounter~=0
HorizontalCenter=floor(horizontalcenter/itsthecounter);
VerticalCenter=floor(verticalcenter/itsthecounter);
else%we dont change center
    
end
%--remake box center

%---end---make box location more accurate---
if(presence==1 )
%--------------beginning of box code---------
Box=zeros(1,uv);

    %---start-check if box fits in picture---
   if HorizontalCenter-lickit<1 
       HorizontalCenter=lickit+1; end%what if parameters of box are out of range? default to a box hugging
    %the side of the frame
    if HorizontalCenter+lickit>thekoran
        HorizontalCenter=thekoran-lickit; end
        if VerticalCenter-suckit<1
            VerticalCenter=suckit+1;end
        if VerticalCenter+suckit>korean
            VerticalCenter=korean-suckit;end
    %---end-check if box fits in picture
    
    %---start--check if box is valid...and if car is actually present
    detectortally=0;
    for kimchi=HorizontalCenter-lickit:HorizontalCenter+lickit%horizontal 
   for kimmy=VerticalCenter-suckit:VerticalCenter+suckit%vertical
if (johnny(1,korean*(kimchi-1)+ kimmy)>3*pi/4 & tobes(1,korean*(kimchi-1)+ kimmy)~=0)
detectortally=detectortally+1;

   end
   end
    end
     %---end--check if box is valid...and if car is actually present
    Isthereabox=0;
    if detectortally>size(CarVectors(1,:),2)/10%if # car vectors in box is not even 1/10 of total car vectors, we can assume that a car isn't in the picture
        Isthereabox=1;
    end
    %make vertical sides of box
for kimchi=HorizontalCenter-lickit:lickit*2:HorizontalCenter+lickit%horizontal 
   for kimmy=VerticalCenter-suckit:VerticalCenter+suckit%vertical
Box(1,korean*(kimchi-1)+ kimmy)=2;
u(1,korean*(kimchi-1)+ kimmy)=2;
v(1,korean*(kimchi-1)+ kimmy)=2;

   end
end
for kimchi=HorizontalCenter-lickit:HorizontalCenter+lickit%horizontal 
   for kimmy=VerticalCenter-suckit:suckit*2:VerticalCenter+suckit%vertical
Box(1,korean*(kimchi-1)+ kimmy)=2;
u(1,korean*(kimchi-1)+ kimmy)=2;
v(1,korean*(kimchi-1)+ kimmy)=2;

   end
end

%----------------end of box code--------------
testvectors=find(   Box~=0  );
end%end of statement if a box radius was found
end

quiver(x(CarVectors), y(CarVectors),u(CarVectors),v(CarVectors),0, 'color', 'g', 'linewidth', 2);
quiver(x(BackGroundVectors),y(BackGroundVectors),u(BackGroundVectors),v(BackGroundVectors),0, 'color', 'r', 'linewidth', 2);
if(presence==1 && Isthereabox==1)
quiver(x(testvectors),y(testvectors),u(testvectors),v(testvectors),0,'color','c','linewidth',2);
end