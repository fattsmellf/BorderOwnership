function BOmodel_final

%Model of the interactions between BO tuned cells (e.g. in V4) and V1.
%The 'ramps' variable controls whether a a prioori sawtooth linear ramp is used for the connections cheme (ramps = 1) or a post-hoc connections cheme based on teh nosie correlation data (ramps = 0)
ramps = 0;

%Master paramters
nPix = 200; %pixels in complete field
wdeg = 25;  %Size of complete field
fgsz = 4;  %Size of figure
pixperdeg = nPix./wdeg;
halfhole = 2; %size of n/u center
halfspan = 6;

%% Co-ordinate system and figure logical matrix
XV = linspace(-wdeg/2,wdeg/2,nPix);
[X,Y] = meshgrid(linspace(-wdeg/2,wdeg/2,nPix+1));
X = X(1:end-1,1:end-1);
Y = Y(1:end-1,1:end-1);
nu = 1;
if nu
    %u shape
    figcx = (X>(-fgsz/(fgsz/halfhole)) & X<(fgsz/(fgsz/halfhole))) & (Y>(-fgsz/(fgsz/halfspan)) & Y<(fgsz/(fgsz/halfhole)));
else
    %simple square
    figcx = (X>(-fgsz/(fgsz/halfhole)) & X<(fgsz/(fgsz/halfhole))) & (Y>(-fgsz/(fgsz/halfhole)) & Y<(fgsz/(fgsz/halfhole)));
end
figcx2 = (X>(-fgsz/(fgsz/halfspan)) & X<(fgsz/(fgsz/halfspan))) & (Y>(-fgsz/(fgsz/halfspan)) & Y<(fgsz./(fgsz/halfspan)));
figcx = figcx2-figcx;
figure,imagesc(XV,XV,figcx)

%% Make a checkerboard
check = ones(size(figcx));
xo = [-16:8:16];
y = [-10:8:14];
for j = 1:length(y)
    if j==1|j==3
        x = xo+4;
    else
        x = xo;
    end
    for i = 1:length(x)
        check((X>x(i)-2)&(X<=x(i)+2)&(Y>y(j)-4)&(Y<=y(j)+4)) = 0;
    end
end
figure,imagesc(XV,XV,check)


%% V4 cells
%4 degrees square, spaced throughout field
%4 maps, one for each BO tuning
%RF is either 1 or -1, when multiplied by the figure logical matrix they give
%stronger responses when the figure is on a certain side of the RF
clear V4
V4(:,:,1) = [ones(2*pixperdeg,4*pixperdeg).*-1;ones(2*pixperdeg,4*pixperdeg)]; %0 DOWN
V4(:,:,2)  = V4(:,:,1)'; %270
V4(:,:,3)  = flipud(V4(:,:,1)); %180
V4(:,:,4)  = fliplr(V4(:,:,2)); %90

%Convolve with figure
clear M C
for q = 1:4
    M(:,:,q) = conv2(figcx,rot90(V4(:,:,q),2),'same');
    M(:,:,q) = 2.*abs(M(:,:,q))+M(:,:,q);
    buf = conv2(check,rot90(V4(:,:,q),2),'same');
    C(:,:,q) = 2.*abs(buf);
end

%Image the V4 BO responses
titl{1} = 'Pref Figure DOWN';
titl{2} = 'Pref Figure RIGHT';
titl{3} = 'Pref Figure UP';
titl{4} = 'Pref Figure LEFT';
%U shape
figure
for q = 1:4
    subplot(2,2,q),imagesc(XV,XV,M(:,:,q)),axis off,title(titl{q}),drawnu([0,0,0]),axis equal
end
%Checkerboard
figure
for q = 1:4
    subplot(2,2,q),imagesc(XV,XV,C(:,:,q)),axis off,title(['Check ',titl{q}]),drawcheck([0,0,0]),axis equal
end

%%
%Now make connection scheme between V4 and V1
%V4 cells are assumed to connect to V1 cells according to their BO
%preference. For example a cell that prefers 'FIGURE DOWN' will
%preferentially connect to V1 cells located in teh downwards direction from
%the center of the V4 RF. This is modelled as a connection gradient, taken
%from the noise correlation data or using linear ramps (controlled by teh
%ramps flag)

%Co-ordiantes for V4
v4sz = 4;
xv = linspace(-8,8,v4sz*2*pixperdeg+1);
[x,y] = meshgrid(xv);
x = x(1:end-1,1:end-1);
y = y(1:end-1,1:end-1);
xv = xv(1:end-1);


if ramps

    %Make linear ramps in different directions
    orientation = 90;
    ramp_0 = cos(orientation*pi/180)*x - sin(orientation*pi/180)*y;
    orientation = 180;
    ramp_90 = cos(orientation*pi/180)*x - sin(orientation*pi/180)*y;
    orientation = 270;
    ramp_180 = cos(orientation*pi/180)*x - sin(orientation*pi/180)*y;
    orientation = 0;
    ramp_270 = cos(orientation*pi/180)*x - sin(orientation*pi/180)*y;

    % %Normalize from -1 to 1
    ramp_0 = 2.*(ramp_0-min(min(ramp_0)))./(max(max(ramp_0))-min(min(ramp_0)))-1;
    ramp_90 = 2.*(ramp_90-min(min(ramp_90)))./(max(max(ramp_90))-min(min(ramp_90)))-1;
    ramp_180 = 2.*(ramp_180-min(min(ramp_180)))./(max(max(ramp_180))-min(min(ramp_180)))-1;
    ramp_270 = 2.*(ramp_270-min(min(ramp_270)))./(max(max(ramp_270))-min(min(ramp_270)))-1;

    hsz = size(x,1)./2;
    ramp_0(1:hsz,:) = 1-ramp_0(1:hsz,:);
    ramp_180(1:hsz,:) = -1-ramp_180(1:hsz,:);
    ramp_90(:,1:hsz) = 1-ramp_90(:,1:hsz);
    ramp_270(:,1:hsz) = -1-ramp_270(:,1:hsz);
    ramp_0(hsz+1:end,:) = -1-ramp_0(hsz+1:end,:);
    ramp_180(hsz+1:end,:) = 1-ramp_180(hsz+1:end,:);
    ramp_90(:,hsz+1:end) = -1-ramp_90(:,hsz+1:end);
    ramp_270(:,hsz+1:end) = 1-ramp_270(:,hsz+1:end);

    w(:,:,1) = ramp_180;
    w(:,:,2) = ramp_270;
    w(:,:,3) = ramp_0;
    w(:,:,4) = ramp_90;
else
    %USe regression weights
    z = hypot(x,y);
    a = abs(atan2(y,x));
    B = [ -0.0508   -0.0404    0.0119    0.2387];  %Angle, ecc, interaction, offset from noise correlation results
    buf = B(1).*a+B(2).*z+B(4); %Not including intercept
    buf(z>4) = 0; %Set conenctions greater than 4 degrees distance to be zero.
    w(:,:,1) = rot90(buf,3); %down
    w(:,:,2) = buf; %right
    w(:,:,3) = rot90(buf,1); %up
    w(:,:,4) = rot90(buf,2); %left
end

%Visualize connection scheme
if 1
    figure,
    for q = 1:4
        subplot(2,2,q),imagesc(xv,xv,w(:,:,q)),axis off,title(titl{q}),colorbar,axis equal,axis off
    end
    figure,surf(xv,xv,w(:,:,q)),shading interp
    title('Connection Scheme')
end

%Calculate the connection strength for each class of BO cells, this is the
%response of the cell to the figure, multiplied by the connection strength
%ramp for that class.
for q = 1:4
    Mcs(:,:,q) = conv2(M(:,:,q),w(:,:,q),'same');
    Ccs(:,:,q) = conv2(C(:,:,q),w(:,:,q),'same');
end

if 1
    %Visualize the connection strength map
    figure,
    for q = 1:4
        subplot(2,2,q),imagesc(XV,XV,Mcs(:,:,q)),axis off,title(['FB ',titl{q}]),colorbar,drawnu([0,0,0]),axis equal,axis off
    end
    figure,
    for q = 1:4
        subplot(2,2,q),imagesc(XV,XV,Ccs(:,:,q)),axis off,title(['FB Chk',titl{q}]),colorbar,drawcheck([0,0,0]),axis equal,axis off
    end
end

%% Average across BO directions
output = mean(Mcs,3);
output_check = mean(Ccs,3);
figure,imagesc(XV,XV,output),axis off,colorbar,drawnu([0,0,0]),axis equal,axis off,title('U')
figure,imagesc(XV,XV,output_check),axis off,colorbar,drawcheck([0,0,0]),axis equal,axis off,title('Check')


%% Checkboard comparison through cross-section
%Blue = figure, green = check, red = modulation
ix = XV>-8&XV<8;
figure,plot(XV(ix),output(80,ix),'b'),hold on,plot(XV(ix),output_check(80,ix),'g'),plot(XV(ix),output(80,ix)-output_check(80,ix),'r')
yl = get(gca,'YLim');
xlim([-12 12])
hold on,line([-6 -6],yl)
line([-2 -2],yl)
line([2 2],yl)
line([6 6],yl)

%%%%%%%%%%
function drawnu(col)

%Adds a N/U to a graph

halfhole = 2; %size of n/u center
halfspan = 6;

left = -halfspan;
ileft = -halfhole;
iright = halfhole;
right = halfspan;
top = -halfspan;
bot = halfspan;
bail = halfhole;

cx = [left left ileft ileft iright iright right right left];
cy = [bot top top bail bail top top bot bot];

hold on
for j = 1:length(cx)-1
    h = line([cx(j) cx(j+1)],[cy(j) cy(j+1)]);
    set(h,'Color',col)
end

return

%%%%%%%%%%%%%%%%%%%
function drawcheck(col)

%Adds a chekcerboard to a graph
x = -18:4:18;
y = -6:8:10;
hold on
h =xline(x);
set(h,'Color',col)
h =yline(y);
set(h,'Color',col)


return


