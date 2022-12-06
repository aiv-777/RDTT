function varargout = surface1(varargin)
% SURFACE2 Application M-file for surface.fig
%    FIG = SURFACE2 launch surface GUI.
%    SURFACE2('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 10-Oct-2022 20:56:25

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
        
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.
% --------------------------------------------------------------------
function varargout = mnuParameters_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.set(HpbCalc,'Enable','Off');
global Htb1 Htb2 bContCalcMade;
Htb1=findobj(gcbf,'Tag','tb1');
Htb2=findobj(gcbf,'Tag','tb2');
bContCalcMade=0;
dial1;
% --------------------------------------------------------------------
function varargout = mnuParView_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.set(HpbCalc,'Enable','Off');
dial2;

% --------------------------------------------------------------------
function varargout = mnuBegin_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.

global R L Xmax amXnod amRnod mnAmNodes Node AxisXmax AxisYmax AxisXmin AxisYmin;
global maxL Field Coordcontour GridX GridR GStepX GStepR; 
global bContCalcMade bContPaintMade DraftAxisLevel;
global amIsoline bUniformGrid Yzero;

if bContPaintMade==0,
    errordlg('Error: No contour!','Error');
    return;
end

if bContCalcMade==1,
  hold on;
  cla reset;
  %newplot;
  %colormap('copper');
  set(gca,'DataAspectRatio',[1 1 1]);
  Vcont=linspace(maxL/amIsoline,maxL-0.05,amIsoline);
  [C,h]=contour(transpose(Field),Vcont,'-r');
  Coordcontour=C;
  set(gca,'XTickLabel',' ');
  set(gca,'YTickLabel',' ');
  set(gca,'DataAspectRatio',[GStepR GStepX 1]);
  return
end

%building rectangular grid with regular step

if bUniformGrid==1,
   amXnod=round(amRnod*(AxisXmax-AxisXmin)/(AxisYmax-AxisYmin));
end
GridX=linspace(AxisXmin,AxisXmax,amXnod);  % axial coordinates of grid nodes
GridR=linspace(AxisYmin,AxisYmax,amRnod);     % radial coordinates of grid nodes
GStepX=GridX(2)-GridX(1);   %step in axial direction
GStepR=GridR(2)-GridR(1);   %step in radial direction
Field=ones(amXnod,amRnod).*(-10000000);
bContCalcMade=0;

j1=0; %BV-arrays counter
% length(Node)
% Node(:).Type
% pause;
for i=1:mnAmNodes
   % making the matrix of contour' burning vertexes
   if Node(i).Type==1 & i>=2,
      j1=j1+1;
      BVX1(j1)=Node(i).X;
      BVY1(j1)=Node(i).Y;
      BVX0(j1)=Node(i-1).X;
      BVY0(j1)=Node(i-1).Y;
   end
end
BVX=[BVX0,BVX1];
BVY=[BVY0,BVY1];

L=AxisXmax;
R=AxisYmax-AxisYmin;
maxL=0;

Vrtcl=(abs(BVX1-BVX0)<(ones(1,length(BVX1)).*(0.001*R)));
nonVrtcl=not(Vrtcl);
BVXtemp=(BVX1-BVX0).*not(BVX1==BVX0)+(BVX1==BVX0);
K=((BVY1-BVY0)./BVXtemp)+(Vrtcl.*10e+10);
K2=K.^2;

%making martrices for node location procedure


Hwait=waitbar(0,'Nodes selection: please, wait...');

GX=repmat(GridX,amRnod,1);
GY=repmat(GridR',1,amXnod);
PImatrix=repmat(pi,amRnod,amXnod);
size(GX)
size(GY)
%checking: does the node belong to grain area
AlfaSum=zeros(amRnod,amXnod);

for i1=2:mnAmNodes
    NDX0=repmat(Node(i1-1).X,amRnod,amXnod);
    NDY0=repmat(Node(i1-1).Y,amRnod,amXnod);
    NDX1=repmat(Node(i1).X,amRnod,amXnod);
    NDY1=repmat(Node(i1).Y,amRnod,amXnod);
    [Alfa1,RTemp]=cart2pol((NDX0-GX),(NDY0-GY));
    [Alfa2,RTemp]=cart2pol((NDX1-GX),(NDY1-GY));
    CND1=(abs(Alfa2-Alfa1)>PImatrix);
    CND2=Alfa1<zeros(amRnod,amXnod);
    Alfa1=Alfa1+CND1.*CND2.*PImatrix.*2;
    Alfa2=Alfa2+CND1.*not(CND2).*PImatrix.*2;
    AlfaSum=AlfaSum+Alfa2-Alfa1;
    
    %          [Alfa1,RTemp]=cart2pol((Node(i1-1).X-GridX(i)),(Node(i1-1).Y-GridR(j)));
    %          [Alfa2,RTemp]=cart2pol((Node(i1).X-GridX(i)),(Node(i1).Y-GridR(j)));
    %          if abs(Alfa2-Alfa1)>pi,
    %             if Alfa1<0,
    %                Alfa1=Alfa1+2*pi;
    %             else
    %                Alfa2=Alfa2+2*pi;
    %             end
    %          end
    %          dAlfa=Alfa2-Alfa1;
    %          %dAlfa=dAlfa
    %          AlfaSum=AlfaSum+dAlfa;
    waitbar(i1/(mnAmNodes-1),Hwait);
end
delete(Hwait);
Hwait=waitbar(0,'Calculations: please,wait...');

%main cycle
for i=1:amXnod
   for j=1:amRnod
      if abs(abs(AlfaSum(j,i))-2*pi)<=0.05,%node is belonging
         %plot(GridX(i),GridR(j),'-+g','Markersize',2);
         %minL=1000*L;
         %calculating of minimal distance to vertexes
         x0=GridX(i);
         y0=GridR(j);
         DSTV=sqrt((BVX-x0).^2+(BVY-y0).^2);
         minDSTV=min(min(DSTV));
         %calculating of minimal distance to lines
         XH=(x0-K.*(BVY0-y0)+K2.*BVX0)./(K2+1);
         YH=BVY0+K.*(XH-BVX0);
         XH=XH.*nonVrtcl+BVX0.*Vrtcl;
         YH=YH.*nonVrtcl+Vrtcl.*y0;
         COND=(((XH-BVX0).*(XH-BVX1)+(YH-BVY0).*(YH-BVY1))<0);
         DSTL=sqrt((XH-x0).^2+(YH-y0).^2);
         DSTL=DSTL.*COND+not(COND).*10e+10;
         minDSTL=min(min(DSTL));
 
         Field(i,j)=min(minDSTL,minDSTV);
      end
      %pause(0.00001);
   end
   waitbar(i/amXnod,Hwait);
end
delete(Hwait);
maxL=max(max(Field));set(handles.text2,'string',maxL);
hold on;
cla reset;
%newplot;
%colormap('copper');
set(gca,'DataAspectRatio',[1 1 1]);
Vcont=linspace(maxL/amIsoline,maxL*0.99,amIsoline);
[C,h]=contour(transpose(Field),Vcont,'-r');
C=C*(AxisYmax-AxisYmin)/amRnod;
set(gca,'XTickLabel',' ');
set(gca,'YTickLabel',' ');
set(gca,'DataAspectRatio',[GStepR GStepX 1]);
Coordcontour=C;
bContCalcMade=1;
% pause;

% --------------------------------------------------------------------

function varargout = mnuDiagram_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.

global Field AxisXmax AxisXmin AxisYmax AxisYmin Node mnAmNodes;
global GStepX GStepR GridX GridR maxL CurL;
global bContCalcMade bContPaintMade DraftAxisLevel;
global Regime amDiagNodes;

%calculating of burning diagram

if bContPaintMade==0,
    errordlg('Error: No contour!','Error');
    return;
end
if bContCalcMade==0,
    mnuBegin_Callback;
end

j=0;
for i=2:mnAmNodes
   if Node(i).Type==0,
       j=j+1;
       ArmorSect(j).Beg.X=Node(i-1).X;
       ArmorSect(j).Beg.Y=Node(i-1).Y;
       ArmorSect(j).End.X=Node(i).X;
       ArmorSect(j).End.Y=Node(i).Y;
     
       ASX1(j)=Node(i).X;
       ASY1(j)=Node(i).Y;
       ASX0(j)=Node(i-1).X;
       ASY0(j)=Node(i-1).Y;
   end
end
mnAmArmorSections=j;

Vrtcl=(abs(ASX1-ASX0)<(ones(1,length(ASX1)).*(0.001*(AxisXmax-AxisXmin))));
nonVrtcl=not(Vrtcl);
ASXtemp=(ASX1-ASX0).*not(ASX1==ASX0)+(ASX1==ASX0);
K=((ASY1-ASY0)./ASXtemp)+(Vrtcl.*10e+10);
K2=K.^2;





mnAmContours=amDiagNodes;
CrosSectLen=max(GStepX,GStepR); %compare criterion
FValStep=maxL/mnAmContours;
CurFVal=FValStep*0.5;
j=0;

Hwait=waitbar(0,'Diagram building: please, wait...');

% cla reset;
%     
% Vax=[AxisXmin AxisXmax AxisYmin AxisYmax];
% axis(Vax);
% hold on;

CurCont(1).X=0;
CurCont(1).Y=0;
CurCont(1).Type=0;

while CurFVal<=maxL
   %building contour
   V=[CurFVal maxL*0.9999];
   C=contourc(transpose(Field),V);
   
   k1=1;
   SubContNum=1;
   mnAmContPoints=C(2,k1);
   SubContVal=C(1,k1);
   while abs(SubContVal-CurFVal)<0.0001*CurFVal
      k1=k1+1;
      SubContInfo(SubContNum).Beg=k1;
      SubContInfo(SubContNum).End=(k1+mnAmContPoints-1);
      for i=k1:(k1+mnAmContPoints-1)
         CurCont(i).X=AxisXmin+(C(1,i)-1)*GStepX;
         CurCont(i).Y=AxisYmin+(C(2,i)-1)*GStepR;
         CurCont(i).Type=0;
      end
      k1=k1+mnAmContPoints;%next point after values
      SubContNum=SubContNum+1;
      mnAmContPoints=C(2,k1);
      SubContVal=C(1,k1);
   end%while
   AmSubCont=SubContNum-1;
   
   %plotting initial contour
%    pause(0.0001);
%    cla;
%    for i=2:mnAmNodes
%     tmpX(1)=Node(i-1).X;
%     tmpY(1)=Node(i-1).Y;
%     tmpX(2)=Node(i).X;
%     tmpY(2)=Node(i).Y;
%     if Node(i).Type==1,
%        plot(tmpX,tmpY,'-r*','Markersize',2);
%     else
%        plot(tmpX,tmpY,'-k*','Markersize',2);
%     end
%    end 
%    set(gca,'DataAspectRatio',[1 1 1]);
   
   
   Sq=0;
   bCalcIt=1;
   %contour square calculating
   for SubContNum=1:AmSubCont
   for i=SubContInfo(SubContNum).Beg+1:SubContInfo(SubContNum).End
       %Does this contour-section belong to armoring surface?
       x0=CurCont(i).X;         %current contour point
       y0=CurCont(i).Y;
       
       XH=(x0-K.*(ASY0-y0)+K2.*ASX0)./(K2+1);
       YH=ASY0+K.*(XH-ASX0);
       XH=XH.*nonVrtcl+ASX0.*Vrtcl;
       YH=YH.*nonVrtcl+Vrtcl.*y0;
       COND=(((XH-ASX0).*(XH-ASX1)+(YH-ASY0).*(YH-ASY1))<0);
       DSTL=sqrt((XH-x0).^2+(YH-y0).^2);
       DSTL=DSTL.*COND+not(COND).*AxisXmax;
       
       REZ=(DSTL<CrosSectLen);
       if sum(REZ)==0
          CurCont(i).Type=1;
       else
          CurCont(i).Type=0;
       end
       
%        for i1=1:mnAmArmorSections
%           
%           xa=ArmorSect(i1).Beg.X;  %current arm-surf section
%           ya=ArmorSect(i1).Beg.Y;
%           xc=ArmorSect(i1).End.X;
%           yc=ArmorSect(i1).End.Y;
%           
%           if abs(xa-xc)>0.0001*(AxisXmax-AxisXmin),
%              k=(yc-ya)/(xc-xa);
%              xh=(x0-k*(ya-y0)+(k^2)*xa)/(1+k^2);
%              yh=ya+k*(xh-xa);
%           else
%              xh=xa;
%              yh=y0;
%           end
%           cond=(xh-xa)*(xh-xc)+(yh-ya)*(yh-yc);
%           if cond<0,
%              CurL=sqrt((xh-x0)^2+(yh-y0)^2);
%           else
%              CurL=AxisXmax;
%           end
%           %Is the point near arm-surface? ----main condition
%           if CurL<=CrosSectLen,
%               CurCont(i).Type=0;
%               break;
%           else
%               CurCont(i).Type=1;
%           end
% 
%        end %armor-nodes cycle
       if CurCont(i).Type==0 & CurCont(i-1).Type==0,
           bCalcIt=0;
       else
           bCalcIt=1;
           switch Regime
           case 1 
               r1=CurCont(i-1).Y;
               r2=CurCont(i).Y;
               tempDist=sqrt((r2-r1)^2+(CurCont(i).X-CurCont(i-1).X)^2);
               SinAlfa=(r2-r1)/tempDist;
               if SinAlfa==0,
                   Sq=Sq+2*pi*r1*abs(CurCont(i).X-CurCont(i-1).X);
               else
                   Sq=Sq+pi*(r2^2-r1^2)/SinAlfa;
               end
           case 2
               x1=CurCont(i-1).X;
               x2=CurCont(i).X;
               y1=CurCont(i-1).Y;
               y2=CurCont(i).Y;
               tempDist=sqrt((x1-x2)^2+(y1-y2)^2);
               Sq=Sq+tempDist;
           end %switch
       end
       
       %plotting contour section
      
%        tmpX(1)=CurCont(i-1).X;
%        tmpY(1)=CurCont(i-1).Y;
%        tmpX(2)=CurCont(i).X;
%        tmpY(2)=CurCont(i).Y;
%        
%        if bCalcIt==1,
%           plot(tmpX,tmpY,'-r*','Markersize',2);
%        else
%           plot(tmpX,tmpY,'-k*','Markersize',2);
%        end
%        %pause;
       
   end %subcontour points cycle
   pause(0.0001);
   end %subcontours cycle
   j=j+1;
   Diag(j)=Sq;
   CurFVal=CurFVal+FValStep;
   %pause;
   waitbar(CurFVal/maxL);
end
delete(Hwait);

diagram;
cla reset;
%Vax=[0 AxisXmax-AxisXmin 0 AxisYmax-AxisYmin];
%axis(Vax);
tempX=linspace(0,maxL,length(Diag));
plot(tempX,Diag,'-r','LineWidth',1.5);[DXFname,Pname]=uiputfile('*.xls');filename=strcat(Pname,DXFname);P={'d, mm'; 'S, mm^2'};xlswrite(filename,transpose(P));
sheet=1; xlRange1='A2'; xlswrite(filename,transpose(tempX),sheet,xlRange1); xlRange='B2'; xlswrite(filename,transpose(Diag),sheet,xlRange);
% --------------------------------------------------------------------
function varargout = mnuPaintCont_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbContour.

global R L Xmax amXnod amRnod Node mnAmNodes;
global AxisXmax AxisXmin AxisYmax AxisYmin;
global bContPaintMade bContCalcMade AxPosDef;
global Xr Yr Xk Yk
Node=[];
%set(HpbCalc,'Enable','Off');

%R=0.6; %chamber internal radius
%L=2*R; %chamber length
%Xmax=L/R;

cla;
Vax=[AxisXmin AxisXmax AxisYmin AxisYmax];
axis(Vax);
hold on;
bContCalcMade=0;
bContPaintMade=0;

%first point
NodeNum=1;
[X1,Y1]=ginput(1);
Node(1).X=X1(1);
Node(1).Y=Y1(1);
plot(X1(1),Y1(1),'-r*','Markersize',2);
%second point
[X1,Y1,Butn]=ginput(1);
while not(Butn==27) 
   NodeNum=NodeNum+1;
   Node(NodeNum).X=X1(1);
   Node(NodeNum).Y=Y1(1);
   tmpX=[Node(NodeNum-1).X Node(NodeNum).X];
   tmpY=[Node(NodeNum-1).Y Node(NodeNum).Y];
   
   if Butn==1,
      Node(NodeNum).Type=1;
      plot(tmpX,tmpY,'-r*','Markersize',2);
   elseif Butn==2 | Butn==3,
      Node(NodeNum).Type=0;
      plot(tmpX,tmpY,'-k*','Markersize',2);
   else
      Node(NodeNum).Type=0;
      plot(tmpX,tmpY,'-k*','Markersize',2);
   end
   [X1,Y1,Butn]=ginput(1);
end
NodeNum=NodeNum+1;
Node(NodeNum).X=Node(1).X;
Node(NodeNum).Y=Node(1).Y;
tmpX=[Node(NodeNum-1).X Node(NodeNum).X];
tmpY=[Node(NodeNum-1).Y Node(NodeNum).Y];
Node(NodeNum).Type=Node(NodeNum-1).Type;
if Node(NodeNum-1).Type==1,
    plot(tmpX,tmpY,'-r*','Markersize',2);
else
    plot(tmpX,tmpY,'-k*','Markersize',2);
end;
 %for i=1:mnAmNodes
%    if Node(i).X>AxisXmax,
%        AxisXmax=Node(i).X;
%    end
%    if Node(i).X<AxisXmin,
%        AxisXmin=Node(i).X;
%    end
%    if Node(i).Y>AxisYmax,
%        AxisYmax=Node(i).Y;
%    end
%    if Node(i).Y<AxisYmin,
%        AxisYmin=Node(i).Y;
%    end
%end

bContPaintMade=1;
mnAmNodes=NodeNum;


% --------------------------------------------------------------------
function varargout = pbPaintCont_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbPaintCont.
% Stub for Callback of the uicontrol handles.pbContour.

global R L Xmax amXnod amRnod Node mnAmNodes;
global AxisXmax AxisXmin AxisYmax AxisYmin;
global bContPaintMade bContCalcMade AxPosDef;
Node=[];
%set(HpbCalc,'Enable','Off');

%R=0.6; %chamber internal radius
%L=2*R; %chamber length
%Xmax=L/R;

cla;
Vax=[AxisXmin AxisXmax AxisYmin AxisYmax];
axis(Vax);
hold on;
bContCalcMade=0;
bContPaintMade=0;

%first point
NodeNum=1;
[X1,Y1]=ginput(1);
Node(1).X=X1(1);
Node(1).Y=Y1(1);
plot(X1(1),Y1(1),'-r*','Markersize',2);
%second point
[X1,Y1,Butn]=ginput(1);
while not(Butn==27) 
   NodeNum=NodeNum+1;
   Node(NodeNum).X=X1(1);
   Node(NodeNum).Y=Y1(1);
   tmpX=[Node(NodeNum-1).X Node(NodeNum).X];
   tmpY=[Node(NodeNum-1).Y Node(NodeNum).Y];
   if Butn==1,
      Node(NodeNum).Type=1;
      plot(tmpX,tmpY,'-r*','Markersize',2);
   elseif Butn==2 | Butn==3,
      Node(NodeNum).Type=0;
      plot(tmpX,tmpY,'-k*','Markersize',2);
   else
      Node(NodeNum).Type=0;
      plot(tmpX,tmpY,'-k*','Markersize',2);
   end
   [X1,Y1,Butn]=ginput(1);
end
NodeNum=NodeNum+1;
Node(NodeNum).X=Node(1).X;
Node(NodeNum).Y=Node(1).Y;
tmpX=[Node(NodeNum-1).X Node(NodeNum).X];
tmpY=[Node(NodeNum-1).Y Node(NodeNum).Y];
Node(NodeNum).Type=Node(NodeNum-1).Type;
if Node(NodeNum-1).Type==1,
    plot(tmpX,tmpY,'-r*','Markersize',2);
else
    plot(tmpX,tmpY,'-k*','Markersize',2);
end;
 %for i=1:mnAmNodes
%    if Node(i).X>AxisXmax,
%        AxisXmax=Node(i).X;
%    end
%    if Node(i).X<AxisXmin,
%        AxisXmin=Node(i).X;
%    end
%    if Node(i).Y>AxisYmax,
%        AxisYmax=Node(i).Y;
%    end
%    if Node(i).Y<AxisYmin,
%        AxisYmin=Node(i).Y;
%    end
%end

bContPaintMade=1;
mnAmNodes=NodeNum;



% --------------------------------------------------------------------
function varargout = figMain_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the figure handles.fi
global amXnod amRnod amDiagNodes bUniformGrid;
global bContCalcMade bContPaintMade AxPosDef AmFrames;
global Regime amIsoline DraftAxisLevel ArcDiskret;
global AxisXmax AxisXmin AxisYmax AxisYmin;

amXnod=60;
amRnod=40;
amDiagNodes=25;

bContCalcMade=0;
bContPaintMade=0;

AxPosDef=[41 30 722 452];
AxisXmin=0;
AxisYmin=0;
AxisYmax=1;
AxisXmax=AxPosDef(3)/AxPosDef(4);

Vax=[AxisXmin AxisXmax AxisYmin AxisYmax];
%axis(Vax);

Regime=1;

bUniformGrid=1;
AmFrames=30;
amIsoline=20;
DraftAxisLevel=0.5;
ArcDiskret=0.1;

%[X,map]=imread('but5.bmp');
%RGB=ind2rgb(X,map);
%image(RGB);
% --------------------------------------------------------------------
function varargout = mnuReadDXF_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbContour.

global mnAmNodes Node AxisXmax AxisYmax AxisXmin AxisYmin;
global bContPaintMade bContCalcMade DraftAxisLevel AxPosDef AxisCoef; 
global ArcDiskret AxisDistance Regime;
global LinType bInvert; %data for reading dialog window

HaxMain=findobj(gcbf,'Tag','axMain');
%set(HaxMain,'Position',AxPosDef);
Htb1=findobj(gcbf,'Tag','tb1');
Htb2=findobj(gcbf,'Tag','tb2');


[DXFname,Pname]=uigetfile('*.dxf');
filename=strcat(Pname,DXFname);
dxf=fopen(filename,'rt');
sTemp=fgetl(dxf);
sTemp=fgetl(dxf);
%[sTemp]=textread(filename,'%s')
pause(0.01);
bContPaintMade=0;
bContCalcMade=0;
%from vb
Cont(1).Kind='a';
Cont(1).LinType=0;
Cont(1).X(1)=0;
Cont(1).Y(1)=0;
Cont(1).Z(1)=0;

Hwait=waitbar(0.01,'Please, wait...');
%-----------------------reading text strings-------------------------

while strcmp(sTemp,'ENTITIES')==0
  sTemp=fgetl(dxf);
end

i = 0;
bStop = 0;
bRead = 0;
LinType(1).name='';%burning surface
LinType(1).amount=0;
LinType(2).name='';%armoring
LinType(2).amount=0;
LinType(3).name='';%axis type
LinType(3).amount=0;

bPolyline=0; %for reading: is polyline belongs to burning or armor. surf
nPolylinType=0;
while (feof(dxf)==0 & bStop==0)
  sCodeLine=fgetl(dxf);
  sValLine=fgetl(dxf);
  nCodeLine = str2num(sCodeLine);
  switch nCodeLine
     case 0
        switch sValLine
           case 'LINE'
               i = i + 1;
               bRead = 1;
               Cont(i).Kind = sValLine;
           case 'ARC'
               i = i + 1;
               bRead = 1;
               Cont(i).Kind = sValLine;
           case 'ENDSEC'
               bStop = 1;
               bRead = 0;
           case 'POLYLINE'
               bPolyline=1;
               nVert = 0;
               bRead = 0;
           case 'VERTEX'
               i = i + 1;
               nVert = nVert + 1;
               bRead = 1;
               Cont(i).Kind = sValLine;
               Cont(i).X(2) = nVert;
               Cont(i).LinType=nPolylinType;
           otherwise
               bRead = 0;
           end
     case 6
        if (strcmp(LinType(1).name,'')==1 | strcmp(LinType(1).name,sValLine)==1)
            LinType(1).name=sValLine;
            LinType(1).amount=LinType(1).amount+1;
        elseif (strcmp(LinType(2).name,'')==1 | strcmp(LinType(2).name,sValLine)==1)
            LinType(2).name=sValLine;
            LinType(2).amount=LinType(2).amount+1;
        elseif (strcmp(LinType(3).name,'')==1 | strcmp(LinType(3).name,sValLine)==1)
            LinType(3).name=sValLine;
            LinType(3).amount=LinType(3).amount+1;
        end
        
        if bPolyline==1,
            nPolylinType = sValLine;
            bPolyline=0;
        end
        if bRead==1
            Cont(i).LinType = sValLine;
        end
     case 10
        if bRead==1
            Cont(i).X(1) = str2num(sValLine);
        end
     case 20
        if bRead==1
            Cont(i).Y(1) = str2num(sValLine);
        end
     case 30
        if bRead==1
            Cont(i).Z(1) = str2num(sValLine);
        end
     case 11
        if bRead==1
            Cont(i).X(2) = str2num(sValLine);
        end
     case 21
        if bRead==1
            Cont(i).Y(2) = str2num(sValLine);
        end
     case 31
        if bRead==1
            Cont(i).Z(2) = str2num(sValLine);
        end
     case 40 %arc, radius
        if bRead==1
            Cont(i).X(2) = str2num(sValLine);
        end
     case 50 %arc, start angle
        if bRead==1
            Cont(i).Y(2) = str2num(sValLine)*pi/180;
        end
     case 51 %arc, end angle
        if bRead==1
            Cont(i).Y(3) = str2num(sValLine)*pi/180;
        end
     otherwise
        %aaa
    waitbar(i/50);    
    end
end
mnAmFigures = i;
fclose(dxf);

%sorting linetypes found in draft
maxAL=0;
for i=1:3
    if LinType(i).amount>maxAL
        maxAL=LinType(i).amount;
        maxIndex=i;
    end
end
%the majority of lines has burning linetype
LinType(4)=LinType(1);
LinType(1)=LinType(maxIndex);
LinType(maxIndex)=LinType(4);

%-----------------making Line record-structure------------------------

Line(1).Burn=0;  %1 means it's burning part; 0-armoring, -1 - axis
Line(1).Type='';
Line(1).Beg.X=0; %part begining
Line(1).Beg.Y=0;
Line(1).End.X=0; %part end
Line(1).End.Y=0;
Line(1).Tag=0;

bVertex=0; %for correct reading of vertex chain
ArcSection=ArcDiskret*pi;
j=0;
for i=1:mnAmFigures
    switch Cont(i).Kind
    case 'LINE'
      bVertex=0;
      j=j+1;
      if strcmp(Cont(i).LinType,LinType(1).name)==1
          Line(j).Burn=1;
      else
          Line(j).Burn=0;
      end
      Line(j).Type=Cont(i).LinType;
      
      Line(j).Beg.X=Cont(i).X(1);
      Line(j).Beg.Y=Cont(i).Y(1);
      Line(j).End.X=Cont(i).X(2);
      Line(j).End.Y=Cont(i).Y(2);
      Line(j).Tag=0;
    case 'ARC'
      bVertex=0;
      if Cont(i).Y(2)>Cont(i).Y(3), %correction
          Cont(i).Y(2)=Cont(i).Y(2)-2*pi;
      end
      Ang=Cont(i).Y(2);
      while Ang+ArcSection<Cont(i).Y(3)
          j=j+1;
          if strcmp(Cont(i).LinType,LinType(1).name)==1
             Line(j).Burn=1;
          else
             Line(j).Burn=0;
          end
          Ang2=Ang+ArcSection;
          Line(j).Beg.X=Cont(i).X(1)+Cont(i).X(2)*cos(Ang);
          Line(j).Beg.Y=Cont(i).Y(1)+Cont(i).X(2)*sin(Ang);
          Line(j).End.X=Cont(i).X(1)+Cont(i).X(2)*cos(Ang2);
          Line(j).End.Y=Cont(i).Y(1)+Cont(i).X(2)*sin(Ang2);
          Line(j).Tag=0;
          Ang=Ang2;
      end%while
      %last line building
      j=j+1;
      if strcmp(Cont(i).LinType,LinType(1).name)==1
         Line(j).Burn=1;
      else
         Line(j).Burn=0;
      end
      Ang2=Cont(i).Y(3);
      Line(j).Beg.X=Cont(i).X(1)+Cont(i).X(2)*cos(Ang);
      Line(j).Beg.Y=Cont(i).Y(1)+Cont(i).X(2)*sin(Ang);
      Line(j).End.X=Cont(i).X(1)+Cont(i).X(2)*cos(Ang2);
      Line(j).End.Y=Cont(i).Y(1)+Cont(i).X(2)*sin(Ang2);
      Line(j).Tag=0;
    case 'VERTEX'
      if bVertex==1, %this is not a begining of polyline
          j=j+1;
          if strcmp(Cont(i).LinType,LinType(1).name)==1
             Line(j).Burn=1;
          else
             Line(j).Burn=0;
          end
          Line(j).Beg.X=Cont(i-1).X(1);
          Line(j).Beg.Y=Cont(i-1).Y(1);
          Line(j).End.X=Cont(i).X(1);
          Line(j).End.Y=Cont(i).Y(1);
          Line(j).Tag=0;
      end
      bVertex=1;
    end%switch
    waitbar((j+mnAmFigures)/50);
end%for
mnAmLines=j;

%for i=1:mnAmLines
%     plot(Line(i).Beg.X,Line(i).Beg.Y,'-mo','Markersize',5);
%     plot(Line(i).End.X,Line(i).End.Y,'-k+','Markersize',5);
%end

%-------------representation contour into Node array-------------------

Condit=sqrt((Line(1).End.X-Line(1).Beg.X)^2+(Line(1).End.Y-Line(1).Beg.Y)^2)*0.00001; %if difference between two 
%nodes less than Condit, - it means they are the same node

%moving burning line to first place (for correct reading of axis)
for i=1:mnAmLines,
    if Line(i).Burn==1 & not(i==1),
        tmpRec=Line(1);
        Line(1)=Line(i);
        Line(i)=tmpRec;
        break;
    end
end
waitbar(1);

minX=10000000;
minY=10000000;
maxX=0;
maxY=0;

nLineNum=1;
Line(1).Tag=1;
k=1; %current node counter
bContourIsClosed=0;
Node(k).X=Line(1).End.X;
Node(k).Y=Line(1).End.Y;
Node(k).Type=Line(1).Burn;
for i=1:mnAmLines %current line counter  
  for j=1:mnAmLines
    if Line(j).Tag==0,
      if abs(Line(j).Beg.X-Node(k).X)<=Condit & abs(Line(j).Beg.Y-Node(k).Y)<=Condit,
          nLineNum=j;
          Line(nLineNum).Tag=1;
          k=k+1;
          Node(k).X=Line(nLineNum).End.X;
          Node(k).Y=Line(nLineNum).End.Y; 
          Node(k).Type=Line(nLineNum).Burn;
          break;
      elseif abs(Line(j).End.X-Node(k).X)<=Condit & abs(Line(j).End.Y-Node(k).Y)<=Condit,
          nLineNum=j;
          Line(j).Tag=1;
          k=k+1;
          Node(k).X=Line(j).Beg.X;
          Node(k).Y=Line(j).Beg.Y;
          Node(k).Type=Line(nLineNum).Burn;
          break;
      end %if
    end %if
  end%for
  
  if Node(k).X>maxX,
      maxX=Node(k).X;
  end
  if Node(k).X<minX,
      minX=Node(k).X;
  end
  if Node(k).Y>maxY,
      maxY=Node(k).Y;
  end
  if Node(k).Y<minY,
      minY=Node(k).Y;
  end
  %checking: does this node close the contour
  if abs(Line(1).Beg.X-Node(k).X)<=Condit & abs(Line(1).Beg.Y-Node(k).Y)<=Condit,
      k=k+1;
      Node(k).X=Node(1).X;
      Node(k).Y=Node(1).Y;
      Node(k).Type=Node(1).Type;
      bContourIsClosed=1;
      break; %the contour is closed
  end
end%for
mnAmNodes=k;
delete(Hwait);
%-------------------plotting contour----------------------------------

if bContourIsClosed>=0,%==1,
    %search for axis
    Regime=2; %default value - the contour represens flat grain
    AxisIndex=-1;
    for i=1:mnAmLines
        if not(Line(i).Tag==1) & abs(Line(i).Beg.Y-Line(i).End.Y)<0.01 & Line(i).Burn==0,
           DraftAxisLevel=Line(i).Beg.Y;
           %search for axis lintype
           if strcmp(Line(i).Type,LinType(2).name)==1
              LinType(4)=LinType(2);
              LinType(2)=LinType(3);
              LinType(3)=LinType(4);
           end
           Line(i).Burn=-1; %it's an axis
           Regime=1; %the contour represents cylindric-shaped grain
           AxisIndex=i;
        end
    end
    if Regime==1,
       set(Htb1,'Value',1.0);
       set(Htb2,'Value',0.0);
       for i=1:mnAmNodes
           Node(i).X=Node(i).X-minX;
           Node(i).Y=Node(i).Y-DraftAxisLevel;
       end
       maxY=maxY-DraftAxisLevel;
       minY=minY-DraftAxisLevel;
       maxX=maxX-minX;
       minX=0;
    else
       set(Htb2,'Value',1.0);
       set(Htb1,'Value',0.0);
       for i=1:mnAmNodes
           Node(i).X=Node(i).X-minX;
           Node(i).Y=Node(i).Y-minY;
       end
       maxY=maxY-minY;
       minY=0;
       maxX=maxX-minX;
       minX=0;
    end

    XSize=maxX-minX;
    YSize=maxY-minY;
    
    AxisXmax=maxX+0.07*XSize;
    AxisXmin=-0.07*XSize;
    AxisYmax=maxY+0.07*YSize;
    AxisYmin=minY-0.07*YSize;
    
    %plotting contour
    cla reset;
    
    Vax=[AxisXmin AxisXmax AxisYmin AxisYmax];
    axis(Vax);
    hold on;
    
    tmpX=zeros(2);
    tmpY=zeros(2);
    
    for i=2:mnAmNodes
        tmpX(1)=Node(i-1).X;
        tmpY(1)=Node(i-1).Y;
        tmpX(2)=Node(i).X;
        tmpY(2)=Node(i).Y;
        if Node(i).Type==1,
            plot(tmpX,tmpY,'-r','LineWidth',1);
        else
            plot(tmpX,tmpY,'-k','LineWidth',1);
        end
        %pause(0.2);
    end
%     length(Node)
%     Node(:).X    
    set(gca,'DataAspectRatio',[1 1 1]);
    %set(HaxMain,'Xlim',[AxisXmin AxisXmax]);
    %set(HaxMain,'Ylim',[AxisYmin AxisYmax]);
    %set(HaxMain,'Position',AxPos);
    
    bContPaintMade=1;

else %it means that draft is wrong 
   warndlg('Error: Contour is not closed!','Error');
   bContPaintMade=0;
end

dialread; %show dialog
    
% --------------------------------------------------------------------
function varargout = mnuSingleSurface_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.set(HpbCalc,'Enable','Off');
global Field DistForPlot maxL GStepX GStepR;
global bContPaintMade bContCalcMade;

if bContPaintMade==0,
    errordlg('Error: No contour!','Error');
    return;
end
if bContCalcMade==0,
    mnuBegin_Callback;
end

[X1,Y1]=ginput(1);
X1=round(X1);
Y1=round(Y1);
DistForPlot=Field(X1,Y1);
if DistForPlot<0,
    DistForPlot=0;
end
hold on;
cla reset;
%newplot;
colormap('copper');
Vcont=[DistForPlot maxL];
contourf(transpose(Field),Vcont);
set(gca,'XTickLabel',' ');
set(gca,'YTickLabel',' ');
set(gca,'DataAspectRatio',[GStepR GStepX 1]);

% --------------------------------------------------------------------
function varargout = mnuSavePicture_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbReadDXF.
[DXFname,Pname]=uiputfile('*.jpg');
filename=strcat(Pname,DXFname)
HaxMain=findobj(gcbf,'Tag','axMain');
saveas(gcbf,filename,'jpg');


% --------------------------------------------------------------------
function varargout = mnuIsoline_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbReadDXF.
global Field amIsoline maxL GStepX GStepR;
global bContPaintMade bContCalcMade;

if bContPaintMade==0,
    errordlg('Error: No contour!','Error');
    return;
end
if bContCalcMade==0,
    mnuBegin_Callback;
else
  hold on;
  cla reset;
  %newplot;
  %colormap('copper');
  set(gca,'DataAspectRatio',[1 1 1]);
  Vcont=linspace(maxL/amIsoline,maxL-0.05,amIsoline);
  [C,h]=contour(transpose(Field),Vcont,'-r');
  set(gca,'XTickLabel',' ');
  set(gca,'YTickLabel',' ');
  set(gca,'DataAspectRatio',[GStepR GStepX 1]);
end
% --------------------------------------------------------------------
function varargout = pbReadDXF_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbReadDXF.
%mnuReadDXF_Callback
surface1('mnuReadDXF_Callback',gcbo,[],guidata(gcbo));

% --------------------------------------------------------------------
function varargout = pbPaintCont_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbPaintCont.
[X,map]=imread('but1.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);



% --------------------------------------------------------------------
function varargout = pbReadDXF_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbReadDXF.
[X,map]=imread('but2.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);


% --------------------------------------------------------------------
function varargout = pbCalcPar_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbCalcPar.
[X,map]=imread('but3.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);


% --------------------------------------------------------------------
function varargout = pbBurnCalc_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbBurnCalc.
[X,map]=imread('but4.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);


% --------------------------------------------------------------------
function varargout = pbDiagCalc_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbDiagCalc.
[X,map]=imread('but5.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);


% --------------------------------------------------------------------
function varargout = axes1_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the axes handles.axes1.

% --------------------------------------------------------------------
function varargout = tb1_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.tb1.
[X,map]=imread('tbut1.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB,'Value',1.0);
% --------------------------------------------------------------------
function varargout = tb1_Callback(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.tb1.
global Regime;
Htb2=findobj(gcbf,'Tag','tb2');
set(Htb2,'Value',0.0);
Regime=1;
% --------------------------------------------------------------------
function varargout = tb2_Callback(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.tb1.
global Regime;
Htb1=findobj(gcbf,'Tag','tb1');
set(Htb1,'Value',0.0);
Regime=2;
% --------------------------------------------------------------------
function varargout = tb2_CreateFcn(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.tb1.
[X,map]=imread('tbut2.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB,'Value',0.0);

% --------------------------------------------------------------------
function varargout = pbClean_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbClean.
global AxPosDef bContPaintMade bContCalcMade;
cla reset;
bContPaintMade=0;
bContCalcMade=0;
HaxMain=findobj(gcbf,'Tag','axMain');
set(HaxMain,'Position',AxPosDef);set(handles.text2,'string',0);
% --------------------------------------------------------------------
function varargout = pbClean_CreateFcn(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbClean.
 [X,map]=imread('but6.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);
% --------------------------------------------------------------------
function varargout = pbViewPar_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbViewPar.
 [X,map]=imread('but7.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);



% --------------------------------------------------------------------
function varargout = pbSingleSurface_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbSingleSurface.
surface1('mnuSingleSurface_Callback',gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
function varargout = pbAnimation_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pbAnimation.
global Field DistForPlot maxL GStepX GStepR;
global bContPaintMade bContCalcMade AmFrames;

if bContPaintMade==0,
    errordlg('Error: No contour!','Error');
    return;
end
if bContCalcMade==0,
    mnuBegin_Callback;
end

AnimStep=maxL/AmFrames;
DistForPlot=0;

colormap('copper');
hold on;    
cla reset;

while DistForPlot<=maxL   
%cla;
Vcont=[DistForPlot maxL];
contourf(transpose(Field),Vcont);
DistForPlot=DistForPlot+AnimStep;
set(gca,'XTickLabel',' ');
set(gca,'YTickLabel',' ');
set(gca,'DataAspectRatio',[GStepR GStepX 1]);
pause(0.001);
end


% --------------------------------------------------------------------
function varargout = pbSingleSurface_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbSingleSurface.
[X,map]=imread('but9.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);



% --------------------------------------------------------------------
function varargout = pbAnimation_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbAnimation.
[X,map]=imread('but8.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);
% --------------------------------------------------------------------
function varargout = pbHelp_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbAnimation.
[X,map]=imread('but10.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);
% --------------------------------------------------------------------
function varargout = tbZoom_CreateFcn(h, eventdata, handles, varargin)
% Stub for CreateFcn of the uicontrol handles.pbAnimation.
[X,map]=imread('but11.bmp');
RGB=ind2rgb(X,map);
set(gcbo,'CData',RGB);

% --------------------------------------------------------------------
function varargout = tbZoom_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.tbZoom.
v=get(h,'Value');
if v==1,
    zoom on;
else
    zoom off;
end


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bContCalcMade bContPaintMade DraftAxisLevel;
global Node Coordcontour amIsoline;
if bContCalcMade==1
    [DXFname,Pname]=uiputfile('*.dxf');
    [v,p]=size(Coordcontour); %получаем общее количество точек в переменной p
    for g=1:amIsoline %присваиваем нули всем координатам будущих точек
     for l=1:p
     Xcoor(g,l)=0;
     Ycoor(g,l)=0;
     end
    end
    j=2; %начинаем со второй точки, первая вне изолиний, 
         %j-номер точки из Coordcontour
    for gg=1:amIsoline
        A(1,gg)=0;
    end
    for i=1:amIsoline;
        %рисуем первую замкнутую линию
              Xcoor(i,1)=Coordcontour(1,j);
              Ycoor(i,1)=Coordcontour(2,j);
              a=j; %фиксируем координаты точки, с которыми сравниваем последнюю точку
              j=j+1;
              Xcoor(i,2)=Coordcontour(1,j);  %записываем вторую точку
              Ycoor(i,2)=Coordcontour(2,j);
              k=2;
              XY=11*Xcoor(i,k)+Ycoor(i,k); %фиксируем первую точку, с которой сравниваем
        %Сначала делаем матрицу для 1-ой изолинии с помощью цикла while, 
        %пока координаты (обе) последней точки изолинии не совпадут с
        %координатами первой точки в этой матрице
        %Затем сохраняем эту матрицу
        %Переходим к следующей изолинии, взяв за первую точку координаты
        %предыдущей первой точки плюс координаты последней точки предыдущей
        %изолинии
       
                 while XY~=11*Coordcontour(1,a)+Coordcontour(2,a)
                    j=j+1; %номер точки из Coordcontour
                    k=k+1; %номер записываемой точки
                    Xcoor(i,k)=Coordcontour(1,j+1);
                    Ycoor(i,k)=Coordcontour(2,j+1);
                    XY=11*Xcoor(i,k)+Ycoor(i,k);
                 end %нарисовали 1 замкнутую линию, дальше рисуем 2, проверяя каждую точку
                    ki=1;
                    A(1,ki)=k;
                    q=0;
                    while q==0 %количество совпадающих точек среди замкнутых линий  
                        j=j+3; %скачок через 2 точки, затем заново тот же цикл для проверки,
                               %что очередная замкнутая линия не имеет общих точек
                        if j>=p
                            break
                        end
                        t=j;   %переменная для нахождения промежутка очередной замкнутой линии
                        tt=t;  %фиксируем 1 точку для нахождения промежутка замкнутой линии
                        t=t+1;
                        CC=11*Coordcontour(1,t)+Coordcontour(2,t);
                        while CC~=11*Coordcontour(1,tt)+Coordcontour(2,tt)
                           t=t+1;
                           CC=11*Coordcontour(1,t)+Coordcontour(2,t);
                        end
                        for ii=tt:t
                        
                            for iii=1:k
                               if Coordcontour(1,ii)>=(Xcoor(i,iii)-0.1)&&Coordcontour(1,ii)<=(Xcoor(i,iii)+0.1)&&Coordcontour(2,ii)>=(Ycoor(i,iii)-0.1)&&Coordcontour(2,ii)<=(Ycoor(i,iii)+0.1)
                                                                %сравниваем 
                                                                %взятую 
                                                                %точку со 
                                                                %всеми точками
                                                                %предыдущих
                                                                %замкнутых
                                                                %линий
                               q=q+1;
                               end
                            end
                        end
                          if q>0
                             break
                          end
                          j=t;
                          
                          
                          
                          for ii=tt:j
                          k=k+1;
                          Xcoor(i,k)=Coordcontour(1,ii);
                          Ycoor(i,k)=Coordcontour(2,ii);
                          end
                          ki=ki+1;%счётчик замкнутых линий
                          A(1,ki)=k;
                          j=j-1;
                          if t>=p
                              break
                          end
                    end
                 
                    
                    
               
                  

        
  
       %сохраняем полученную изолинию (возможно, из нескольких замкнутых
       %линий)
       isymbol=num2str(i);
       filename=strcat(Pname,isymbol,'-',DXFname);
       FID = dxf_open(filename);
       FID = dxf_set(FID,'Color',[1 0 0]);
       aa=1;
       for ii=1:ki
           xx=1;
           for w=aa:A(1,ii)
               Xcoor1(1,xx)=Xcoor(i,w); %вытаскиваем из матрицы вектор, чтобы 
               Ycoor1(1,xx)=Ycoor(i,w); %не было лишних нулей
               xx=xx+1;
           end
           [m,n]=size(Xcoor1);
           if n>=xx
              for hh=xx:n
                  Xcoor1(:,xx)=[];
                  Ycoor1(:,xx)=[];
              end
           end
           Z=zeros(1, xx-1);
           FID = dxf_polyline(FID, transpose(Xcoor1(1,:)), transpose(Ycoor1(1,:)), transpose(Z));
           
               Xcoor1(:,:)=[];
               Ycoor1(:,:)=[];
           
           aa=A(1,ii)+1;
       end
       dxf_close(FID);
       
    end
elseif bContPaintMade==1 && bContCalcMade==0
[DXFname,Pname]=uiputfile('*.dxf');
filename=strcat(Pname,DXFname);
FID = dxf_open(filename);
FID = dxf_set(FID,'Color',[1 0 0]);
[m,n]=size(Node);
for i=1:n;
Xnode(i)=Node(i).X;
Ynode(i)=Node(i).Y;
end
Z=zeros(n, 1);
FID = dxf_polyline(FID, transpose(Xnode), transpose(Ynode), Z);
dxf_close(FID);
else
    disp('Нет контура');
end


% --- Executes during object creation, after setting all properties.
