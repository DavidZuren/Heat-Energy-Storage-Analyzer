function varargout = HeatStorageAnalyzer(varargin)
% HEATSTORAGEANALYZER MATLAB code for HeatStorageAnalyzer.fig
%      HEATSTORAGEANALYZER, by itself, creates a new HEATSTORAGEANALYZER or raises the existing
%      singleton*.
%
%      H = HEATSTORAGEANALYZER returns the handle to a new HEATSTORAGEANALYZER or the handle to
%      the existing singleton*.
%
%      HEATSTORAGEANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HEATSTORAGEANALYZER.M with the given input arguments.
%
%      HEATSTORAGEANALYZER('Property','Value',...) creates a new HEATSTORAGEANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HeatStorageAnalyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HeatStorageAnalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HeatStorageAnalyzer

% Last Modified by GUIDE v2.5 19-Aug-2021 14:08:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HeatStorageAnalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @HeatStorageAnalyzer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before HeatStorageAnalyzer is made visible.
function HeatStorageAnalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HeatStorageAnalyzer (see VARARGIN)

MatTypeMenu_Callback(hObject, eventdata, handles)
set(handles.tspanEdit,'string','20E6')
handles.Ac = 1;
handles.Line = 1;

handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes HeatStorageAnalyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = HeatStorageAnalyzer_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function PhaseTransitionModel(hObject, eventdata, handles)
clc;
hold off
syms x xi(t) L(t) a1(t) b1(t) a2(t) b2(t) TH TC
syms a_1 a_2
syms T_f L_f rho_l rho_s kl ks Cl Cs

%Initial conditions
f = waitbar(0,'Setting up');
pause(1)
rho_l = str2double(get(handles.rho_lEdit,'string'));
rho_s = str2double(get(handles.rho_sEdit,'string'));
Cl = str2double(get(handles.ClEdit,'string'));
Cs = str2double(get(handles.CsEdit,'string'));
kl = str2double(get(handles.klEdit,'string'));
ks = str2double(get(handles.ksEdit,'string'));
L_f = str2double(get(handles.L_fEdit,'string'));
T_f = str2double(get(handles.T_fEdit,'string'));
TH = str2double(get(handles.THEdit,'string'));
TC = str2double(get(handles.TCEdit,'string'));
L_inic = str2double(get(handles.L_inicEdit,'string'));
xi_inic = str2double(get(handles.xi_inicEdit,'string')); 

%Goodman
Tl(x,t)=a1(t)*(x-xi(t))+b1(t)*(x-xi(t))^2+T_f;
Ts(x,t)=a2(t)*(x-xi(t))+b2(t)*(x-xi(t))^2+T_f;

%Boundary conditions
eq1=Tl(0,t)==TH;
eq2=Ts(L(t),t)==TC;

waitbar(0.2,f,'Building matrices');
pause(0.5)
%Motion equation of the interface
deqInt=L_f*rho_l*diff(xi(t),t)==-kl*subs(diff(Tl(x,t),x),x,xi(t))+ks*subs(diff(Ts(x,t),x),x,xi(t));

%Heat equations
deqQ1=rho_l*Cl*simplify(diff(int(Tl(x,t),x,0,xi(t)),t)-Tl(xi(t),t)*diff(xi(t),t))-kl*int(diff(diff(Tl(x,t),x),x),x,0,xi(t))==0;
deqQ2=rho_s*Cs*simplify(diff(int(Ts(x,t),x,xi(t),L(t)),t)+Ts(xi(t),t)*diff(xi(t),t))-ks*int(diff(diff(Ts(x,t),x),x),x,xi(t),L(t))==0;

%Conservation of mass equation
deqM = rho_l*(diff(xi(t),t))+rho_s*(diff(L(t),t)-diff(xi(t),t))==0;

%Solving for a1 and a2 with boundary conditions
a_1=solve(subs(eq1,a1(t),a_1),a_1);
a_2=solve(subs(eq2,a2(t),a_2),a_2);

%Substituting in differential equations
deqInt1=subs(deqInt,{a1(t),a2(t)},{a_1,a_2});
deqInt11=collect(subs(deqInt1,{L_f rho_l kl ks TC TH T_f},{L_f rho_l kl ks TC TH T_f}),diff(xi(t),t));
deqInt111 = collect(lhs(deqInt11)-rhs(deqInt11),diff(xi(t),t))==0;

%Matrix
NOPER_deqInt111 = nop(op(1,deqInt111));
A(1,:) = [op(1,op(1,deqInt111))/diff(xi(t),t) 0 0 0];
B(1,1) = -op(NOPER_deqInt111,op(1,deqInt111));

%Substituting in differential equations
deqQ1_1=subs(deqQ1,{a1(t)},{a_1});
deqQ1_11=collect(subs(deqQ1_1,{rho_l Cl kl TH T_f},{rho_l Cl kl TH T_f}),[diff(xi(t),t) diff(b1(t),t)]);
deqQ1_111=collect(lhs(deqQ1_11)-rhs(deqQ1_11),[diff(xi(t),t) diff(b1(t),t)] )==0;

%Matrix
NOPER_deqQ1_111 = nop(op(1,deqQ1_111));
A(2,:) = [op(1,op(1,deqQ1_111))/diff(xi(t),t) op(2,op(1,deqQ1_111))/diff(b1(t),t) 0 0];
B(2,1) = -op(NOPER_deqQ1_111,op(1,deqQ1_111));

%Substituting in differential equations
deqQ2_1=subs(deqQ2,{a2(t)},{a_2});
deqQ2_11=collect(subs(deqQ2_1,{rho_s Cs ks TC T_f},{rho_s Cs ks TC T_f}),[diff(xi(t),t) diff(b1(t),t)]);
deqQ2_111=collect(lhs(deqQ2_11)-rhs(deqQ2_11),[diff(xi(t),t) diff(b2(t),t)] )==0;

%Matrix
NOPER_deqQ2_111 = nop(op(1,deqQ2_111));
A(3,:) = [op(1,op(1,deqQ2_111))/diff(xi(t),t) 0 op(2,op(1,deqQ2_111))/diff(b2(t),t) 0];
B(3,1) = -op(NOPER_deqQ2_111,op(1,deqQ2_111));

%Substituting in differential equations
deqM_1=collect(subs(deqM,{rho_l rho_s},{rho_l rho_s}),[diff(xi(t),t) diff(L(t),t)]);
deqM_11=collect(lhs(deqM_1)-rhs(deqM_1),[diff(xi(t),t) diff(L(t),t)])==0;

%Matrix
NOPER_deqM_11 = nop(op(1,deqM_11));
if rho_l~=rho_s
    A(4,:) = [op(1,op(1,deqM_11))/diff(xi(t),t) 0 0 op(2,op(1,deqM_11))/diff(L(t),t)];
else
    A(4,:) = [0 0 0 op(1,op(1,deqM_11))];
end
B(4,1) = 0;

syms y1 y2 y3 y4;
AA = simplify(subs(A,{xi(t), b1(t), b2(t), L(t)},{y1, y2, y3, y4}));
BB = simplify(subs(B,{xi(t), b1(t), b2(t), L(t)},{y1, y2, y3, y4}));

%Defining a1, a2, b1, b2
syms aa1 aa2 bb1 bb2;

%Liquid temperature distribution
Tliq = T_f+aa1*(x-xi_inic)+bb1*(x-xi_inic)^2;
%Solid temperature distribution
Tsol = T_f+aa2*(x-xi_inic)+bb2*(x-xi_inic)^2;

waitbar(0.3,f,'Applying BC');
pause(0.5)
%%____________________BOUNDATY CONDITIONS ____________________%%
%Isotermal BC for the liquid
eq3 = subs(Tliq,x,0)==TH;
%Adiabatic BC on Tliq for t=0 to find b1
eq4 = subs(diff(Tliq,x),x,0)==0;
%Isotermal BC for the solid
eq5 = subs(Tsol,x,L_inic)==TC;
%Adiabatic BC on Tsol for t=0 to find b2
eq6 = subs(diff(Tsol,x),x,L_inic)==0;

%Solving for a1 and b1
[aa1,bb1] = solve([eq3,eq4],[aa1,bb1]);
aa1=double(aa1);
bb1=double(bb1);
%Solving for a2 and b2
[aa2,bb2] = solve([eq5,eq6],[aa2,bb2]);
aa2=double(aa2);
bb2=double(bb2);

%Solving the model
waitbar(.4,f,'Solving model (RK)');
y0=[xi_inic;bb1;bb2;L_inic];
Fun = @(t,y) double(subs(AA\BB,{y1, y2, y3, y4},{y(1), y(2), y(3), y(4)}));
tspan = [0,str2double(get(handles.tspanEdit,'string'))];
[t,y]=ode45(Fun,tspan,y0);

waitbar(0.5,f,'Evaluating Limits');
pause(0.4)

xi_lim = (rho_s*ks*(T_f-TC)*0+kl*(TH-T_f)*(rho_s*L_inic-(rho_s-rho_l)*xi_inic))/(rho_l*kl*(TH-T_f)+rho_s*ks*(T_f-TC));
L_lim = ((kl*(TH-T_f)+ks*(T_f-TC)))*(rho_s*L_inic-xi_inic*(rho_s-rho_l))/(rho_l*kl*(TH-T_f)+rho_s*ks*(T_f-TC));

%Mass of the system
Ac = handles.Ac;
M = rho_l*xi_inic*Ac+rho_s*(L_inic-xi_inic)*Ac;
q = -(1/M)*(kl*rho_l*(TH-T_f)+ks*rho_s*(T_f-TC));

waitbar(0.6,f,'Evaluating Masses');
pause(0.4)
%Simpli RK & t
its=224;
MRK=zeros(its,4);
tsim=zeros(its,1);
cc=round(linspace(1,length(t),its));
for i=1:its
     MRK(i,:)=y(cc(i),:);
     tsim(i)=t(cc(i));
end
%Mass of Liquid t
Mlt=zeros(its,1);
for i=1:its
    limsup=MRK(i,1);
    Mlt(i)=Ac*int(rho_l,x,0,limsup);
end
%Mass of Solid t
Mst=zeros(its,1);
for i=1:its
    Mst(i)=M-Mlt(i);
end

waitbar(0.65,f,'Evaluating energy limits');
pause(0.4)
%_______________________Energy functions______________________________
%Since beta = 0
deltaU1lim = double(Cl*Ac*(rho_l*(TH*xi_inic+(1/2)*(q/kl)*xi_inic^2)-int((T_f+aa1*(x-xi_inic)+bb1*(x-xi_inic)^2)*rho_l,x,0,xi_inic)));
deltaMslim = rho_s*Ac*(L_inic-xi_inic-L_lim+xi_lim);
Xa = (deltaMslim/(Ac*rho_s))+xi_inic;
deltaU2lim = double(Cs*deltaMslim*T_f-Cs*Ac*int((T_f+aa2*(x-xi_inic)+bb2*(x-xi_inic)^2)*rho_s,x,xi_inic,Xa));
deltaU3lim = double(Cl*Ac*(rho_l*(TH*(xi_lim-xi_inic)+(1/2)*(q/kl)*(xi_lim^2-xi_inic^2)))-Cl*deltaMslim*T_f);
deltaU4lim = double(Cs*Ac*(rho_s*(((kl/ks)*(TH-T_f)+T_f)*(L_lim-xi_lim)+(1/2)*(q/ks)*(L_lim^2-xi_lim^2)))-Cs*Ac*int(rho_s*(T_f+aa2*(x-xi_inic)+bb2*(x-xi_inic)^2),x,Xa,L_inic));
Qs_lim = deltaU1lim + deltaU2lim + deltaU3lim + deltaU4lim;
Qf_lim = L_f*(deltaMslim);
Qt_lim = Qs_lim + Qf_lim;

waitbar(0.7,f,'Evaluating U1');
Tl=subs(Tl,a1,a_1);
%Energy U1
DeltaU1t=zeros(its,1);
for i=1:its-1
    DeltaU1t(i)=double(Ac*Cl*rho_l*int((subs(Tl,{b1 xi},{MRK(i+1,2) MRK(1+i,1)})),x,0,MRK(1,1)))-double(Ac*Cl*rho_l*int((subs(Tl,{b1 xi},{MRK(1,2) MRK(1,1)})),x,0,MRK(1,1)));
    Percent = (i/its)*100;
    waitbar(0.7,f,sprintf('Evaluating U1: %.1f %',Percent))
end
DeltaU1t(its)=DeltaU1t(its-1);

waitbar(0.75,f,'Evaluating U2');
%Energy U2
xa=zeros(its,1);
for i=1:its-1
    xa(i)=(Mst(1)-Mst(i+1))/(Ac*rho_s)+MRK(1,1);
end
xa(its)=xa(its-1);
Ts=subs(Ts,a2,a_2);
DeltaU2t=zeros(its,1);
for i=1:its-1
    DeltaU2t(i)=(Mst(1)-Mst(i+1))*Cs*T_f-(Ac*Cs*rho_s*int((subs(Ts,{b2 xi L},{MRK(1,3) MRK(1,1) MRK(1,4)})),x,MRK(1,1),xa(i+1)));
    Percent = (i/its)*100;
    waitbar(0.75,f,sprintf('Evaluating U2: %.1f %',Percent))
end
DeltaU2t(its)=DeltaU2t(its-1);

waitbar(0.8,f,'Evaluating U3');
%Energy U3
DeltaU3t=zeros(its,1);
for i=1:its-1
    DeltaU3t(i)=double(Ac*Cl*rho_l*int((subs(Tl,{b1 xi},{MRK(i,2) MRK(i,1)})),x,MRK(1,1),MRK(i+1,1)))-(Mst(1)-Mst(i+1))*Cl*T_f;
    Percent = (i/its)*100;
    waitbar(0.8,f,sprintf('Evaluating U3: %.1f %',Percent))
end
DeltaU3t(its)=DeltaU3t(its-1);

waitbar(0.85,f,'Evaluating U4');
%Energy U4
DeltaU4t=zeros(its,1);
for i=1:its-1
    DeltaU4t(i)=double(Ac*Cs*rho_s*int((subs(Ts,{b2 xi L},{MRK(i+1,3) MRK(i+1,1) MRK(i+1,4)})),x,MRK(i+1,1),MRK(i+1,4))-Ac*Cs*rho_s*int((subs(Ts,{b2 xi L},{MRK(1,3) MRK(1,1) MRK(1,4)})),x,xa(i+1),MRK(1,4)));
    Percent = (i/its)*100;
    waitbar(0.85,f,sprintf('Evaluating U4: %.1f %',Percent))
end
DeltaU4t(its)=DeltaU4t(its-1);

waitbar(0.9,f,'Solving energies');
pause(0.4)
%Total Sensible Energy 
Qs=DeltaU1t+DeltaU2t+DeltaU3t+DeltaU4t;

%Total Latent Energy
Qf=zeros(its,1);
for i=2:its-1
    Qf(i)=(Mst(1)-Mst(i))*L_f;
end
Qf(its)=Qf(its-1);

%Total Energy
Qt=Qs+Qf;

waitbar(0.95,f,'Evaluating temperature profiles');
%________________________Temperature profiles_____________________________
%Crear matriz de z (tiempos a graficar)
dimen=6;
tiempos=its;
Mz1=zeros(1,tiempos);
cm=round(linspace(1,its,tiempos));
for i=1:tiempos
    Mz1(i)=tsim(cm(i));
end
Mz=ones(1,dimen,tiempos);
for i=1:tiempos
    Mz(:,:,i)=Mz1(i)*ones(1,dimen);
end
%crear matriz de x (puntos en x a graficar)
Mxl=zeros(1,dimen,tiempos);
for i=1:tiempos
    Mxl(:,:,i)=linspace(0,MRK(cm(i),1),dimen);
end
Mxs=zeros(1,dimen,tiempos);
for i=1:tiempos
    Mxs(:,:,i)=linspace(MRK(cm(i),1),MRK(cm(i),4),dimen);
end
%crear maztriz de y (temperaturas evaluadas en cada x y en 1 tiempo)
Myl=zeros(1,dimen,tiempos);
for j=1:tiempos
    Y=(subs(Tl,{b1 xi},{MRK(cm(j),2) MRK(cm(j),1)}));
    for i=1:dimen
        x=Mxl(:,i,j);
        Myl(1,i,j)=double(subs(Y));
    end
end
Mys=zeros(1,dimen,tiempos);
for j=1:tiempos
    Y2=(subs(Ts,{b2 xi L},{MRK(cm(j),3) MRK(cm(j),1) MRK(cm(j),4)}));
    for i=1:dimen
        x=Mxs(:,i,j);
        Mys(1,i,j)=double(subs(Y2));
    end
end

PCGN=38E6; %J/m^3
PCGLP=46E6; %J/m^3
VGN=(Qt_lim/PCGN)*1000; %L
VGLP=(Qt_lim/PCGLP)*1000; %L
PGN=4.16; %$/L
PGLP=10.99; %$/L

GGN=(VGN*PGN);
set(handles.NatGasText,'string',GGN)
%fprintf('El costo de utilizar gas natural es:%8.3f\n',GGN)
GGLP=VGLP*PGLP;
set(handles.LPgasText,'string',GGLP)
%fprintf('El costo de utilizar gas LP es:%8.3f\n',GGLP)
nc=double(0.98*Qt_lim);
uts = find((abs(nc-Qt)/nc)*100<0.05,1,'First');
ts=tsim(uts)/3600;
%fprintf('El tiempo de recarga es:%8.3f horas\n',ts)
set(handles.ChT,'string',ts)

waitbar(1,f,'Plotting');
pause(2)

set(handles.CurveTypeMenu,'Enable','on')
set(handles.LineWidth,'Enable','on')

handles.tsim = tsim;
handles.MRK = MRK; 
handles.xi_lim = xi_lim;
handles.L_lim = L_lim;
handles.tspan = tspan;
handles.Mlt = Mlt;
handles.Mst = Mst;
handles.its = its;
handles.rho_l = rho_l;
handles.rho_s = rho_s;
handles.Qs = Qs;
handles.Qs_lim = Qs_lim;
handles.Qf = Qf;
handles.Qf_lim = Qf_lim;
handles.Qt = Qt;
handles.Qt_lim = Qt_lim;
handles.Mz = Mz;
handles.Mxl = Mxl;
handles.Myl = Myl;
handles.Mxs = Mxs;
handles.Mys = Mys;
guidata(hObject, handles);
Plot(hObject, eventdata, handles)
close(f)

function Plot(hObject, eventdata, handles)
tsim = handles.tsim;
MRK = handles.MRK;
tspan = handles.tspan;
L_lim = handles.L_lim;
axes(handles.axes1)
plot(tsim,MRK(:,4));grid on;xlabel('t(s)');ylabel('L(m)','rotation',90);
hold on
plot(tspan,ones(size(tspan)) * L_lim);
legend('L(t)','L lim')
grid minor

% --- Executes on selection change in CurveTypeMenu.
function CurveTypeMenu_Callback(hObject, eventdata, handles)
CurveType = get(handles.CurveTypeMenu,'value');
tsim = handles.tsim;
MRK = handles.MRK;
tspan = handles.tspan;
xi_lim = handles.xi_lim;
L_lim = handles.L_lim;
Mlt = handles.Mlt;
Mst = handles.Mst;
its = handles.its;
Ac = handles.Ac;
rho_l = handles.rho_l;
rho_s = handles.rho_s;
Line = handles.Line;
Qs = handles.Qs;
Qs_lim = handles.Qs_lim;
Qf = handles.Qf;
Qf_lim = handles.Qf_lim;
Qt = handles.Qt;
Qt_lim = handles.Qt_lim;
Mz = handles.Mz;
Mxl = handles.Mxl;
Myl= handles.Myl;
Mxs = handles.Mxs;
Mys = handles.Mys;
hold off
axes(handles.axes1)
if CurveType == 1
    plot(tsim,MRK(:,4),'LineWidth',Line);grid on;xlabel('t(s)');ylabel('L(m)','rotation',90);
    hold on
    plot(tspan,ones(size(tspan)) * L_lim,'LineWidth',Line);
    legend('L(t)','L lim')
elseif CurveType == 2
    plot(tsim,MRK(:,1),'LineWidth',Line);grid on;xlabel('t(s)');ylabel('\xi(m)','rotation',90);
    hold on
    plot(tspan,ones(size(tspan)) * xi_lim,'LineWidth',Line);
    legend('\xi(t)','\xi lim')
elseif CurveType == 3
    plot(tsim,Mlt,'LineWidth',Line);grid on;xlabel('t(s)');ylabel('M(kg)','rotation',90);
    hold on
    plot(tsim,rho_l*(xi_lim)*Ac*ones(its,1),'LineWidth',Line);
    plot(tsim,Mst,'LineWidth',Line);
    plot(tsim,rho_s*(L_lim-xi_lim)*Ac*ones(its,1),'LineWidth',Line);
    legend('Ml(t)','Ml lim','Ms(t)','Ms lim')
elseif CurveType == 4
    plot(tsim,Qs,'LineWidth',Line);grid on;xlabel('t(s)');ylabel('Qs(J)','rotation',90);
    hold on
    plot(tsim,Qs_lim*ones(its,1),'LineWidth',Line);
    legend('Qs(t)','Qs lim')
elseif CurveType == 5
    plot(tsim,Qf,'LineWidth',Line);grid on;xlabel('t(s)');ylabel('Qf(J)','rotation',90);
    hold on
    plot(tsim,Qf_lim*ones(its,1),'LineWidth',Line);
    legend('Qf(t)','Qf lim')
elseif CurveType == 6
    plot(tsim,Qt,'LineWidth',Line);grid on;xlabel('t(s)');ylabel('Qt(J)','rotation',90);
    hold on
    plot(tsim,Qt_lim*ones(its,1),'LineWidth',Line);
    legend('Qt(t)','Qt lim')
elseif CurveType == 7
    for i=1:its
        plot3(Mz(:,:,i),Mxl(:,:,i),Myl(:,:,i),'c','LineWidth',Line)
        xlabel('t [s]'); ylabel('x [m]'); zlabel('T [°K]');
        set ( gca, 'ydir', 'reverse' )
        hold on
        plot3(Mz(:,:,i),Mxs(:,:,i),Mys(:,:,i),'b','LineWidth',Line)
    view(3)
    end
end
grid minor
set(handles.CalculateBtn,'Enable','on')

function exp2 = op(n,exp1)
s=children(exp1);
s = [s{:}];
exp2 = s(n);
function n = nop(exp1)
s=children(exp1);
s = [s{:}];
n = size(s,2);

function THEdit_Callback(hObject, eventdata, handles)
function THEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TCEdit_Callback(hObject, eventdata, handles)
function TCEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function L_inicEdit_Callback(hObject, eventdata, handles)
function L_inicEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xi_inicEdit_Callback(hObject, eventdata, handles)
function xi_inicEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rho_lEdit_Callback(hObject, eventdata, handles)
function rho_lEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rho_sEdit_Callback(hObject, eventdata, handles)
function rho_sEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ClEdit_Callback(hObject, eventdata, handles)
function ClEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CsEdit_Callback(hObject, eventdata, handles)
function CsEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function klEdit_Callback(hObject, eventdata, handles)
function klEdit_CreateFcn(hObject, eventdata, handles)
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function L_fEdit_Callback(hObject, eventdata, handles)
function L_fEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function T_fEdit_Callback(hObject, eventdata, handles)
function T_fEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in EditBtn.
function EditBtn_Callback(hObject, eventdata, handles)
set(handles.MatTypeMenu,'Enable','on')
set(handles.rho_lEdit,'Enable','on')
set(handles.rho_sEdit,'Enable','on')
set(handles.ClEdit,'Enable','on')
set(handles.CsEdit,'Enable','on')
set(handles.klEdit,'Enable','on')
set(handles.ksEdit,'Enable','on')
set(handles.L_fEdit,'Enable','on')
set(handles.T_fEdit,'Enable','on')
set(handles.THEdit,'Enable','on')
set(handles.TCEdit,'Enable','on')
set(handles.L_inicEdit,'Enable','on')
set(handles.xi_inicEdit,'Enable','on')
set(handles.CalculateBtn,'Enable','on')
set(handles.EditBtn,'Enable','off')
set(handles.tspanEdit,'Enable','on')
    
function CalculateBtn_Callback(hObject, eventdata, handles)
set(handles.rho_lEdit,'Enable','off')
set(handles.rho_sEdit,'Enable','off')
set(handles.ClEdit,'Enable','off')
set(handles.CsEdit,'Enable','off')
set(handles.klEdit,'Enable','off')
set(handles.ksEdit,'Enable','off')
set(handles.L_fEdit,'Enable','off')
set(handles.T_fEdit,'Enable','off')
set(handles.THEdit,'Enable','off')
set(handles.TCEdit,'Enable','off')
set(handles.L_inicEdit,'Enable','off')
set(handles.xi_inicEdit,'Enable','off')
set(handles.CalculateBtn,'Enable','off')
set(handles.EditBtn,'Enable','on')
set(handles.tspanEdit,'Enable','off')
%set(hanldles.LineWidth,'Value',1)
pause(0.2)
PhaseTransitionModel(hObject, eventdata, handles)

% --- Executes on selection change in MatTypeMenu.
function MatTypeMenu_Callback(hObject, eventdata, handles)
MatType = get(handles.MatTypeMenu,'Value');
if MatType == 1
    set(handles.rho_lEdit,'string','1344.6')
    set(handles.rho_sEdit,'string','1505')
    set(handles.ClEdit,'string','2730')
    set(handles.CsEdit,'string','1270')
    set(handles.klEdit,'string','0.36')
    set(handles.ksEdit,'string','0.52')
    set(handles.L_fEdit,'string','237600')
    set(handles.T_fEdit,'string','363.15')
    set(handles.THEdit,'string','473.15')
    set(handles.TCEdit,'string','353.15')
    set(handles.L_inicEdit,'string','1')
    set(handles.xi_inicEdit,'string','0.005593')
    set(handles.tspanEdit,'string','20E6')
elseif MatType == 2
    set(handles.rho_lEdit,'string','2096')
    set(handles.rho_sEdit,'string','2192')
    set(handles.ClEdit,'string','1500')
    set(handles.CsEdit,'string','1430')
    set(handles.klEdit,'string','0.8')
    set(handles.ksEdit,'string','1')
    set(handles.L_fEdit,'string','105000')
    set(handles.T_fEdit,'string','496.15')
    set(handles.THEdit,'string','780')
    set(handles.TCEdit,'string','480')
    set(handles.L_inicEdit,'string','1')
    set(handles.xi_inicEdit,'string','0.005228')
    set(handles.tspanEdit,'string','5E6')
elseif MatType == 3
    set(handles.rho_lEdit,'string','2545')
    set(handles.rho_sEdit,'string','2545')
    set(handles.ClEdit,'string','1130')
    set(handles.CsEdit,'string','1016')
    set(handles.klEdit,'string','225')
    set(handles.ksEdit,'string','215')
    set(handles.L_fEdit,'string','396000')
    set(handles.T_fEdit,'string','933.52')
    set(handles.THEdit,'string','1173')
    set(handles.TCEdit,'string','873.25')
    set(handles.L_inicEdit,'string','1')
    set(handles.xi_inicEdit,'string','0.2')
    set(handles.tspanEdit,'string','2E4')
end

function MatTypeMenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ksEdit_Callback(hObject, eventdata, handles)
function ksEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CurveTypeMenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tspanEdit_Callback(hObject, eventdata, handles)
function tspanEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function LineWidth_Callback(hObject, eventdata, handles)
Line = get(handles.LineWidth,'Value');
handles.Line = Line;
guidata(hObject, handles);
CurveTypeMenu_Callback(hObject, eventdata, handles)

function LineWidth_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end