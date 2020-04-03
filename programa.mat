% Continuous time system with multiplicative noise
% dx  = (A(\alpha(t)) + D(\beta(t)))x + B1*W(t) + B2(\alpha(t))u(t)
% z(t)= C1*(\alpha(t))x + D11(\alpha(t))w(t)
% synthesis condition
%
% Example:
% clear all,clc
% authors: elmerl@utfpr.edu.br 
%          ricfow@dt.fee.unicamp.br 
%          avargas@utfpr.edu.br 
function output = sintesis_robusto_LPV_sys_continuo(A,D,B1,B2,C1,D11,D21,dotbounds,xi,deg)
%
order = size(A{1},2);
p     = size(B1{1},2);
l     = size(B2{1},2);
r     = size(C1{1},1);
N1    = size(A,2);
N2    = size(D,2);

% system matrices
Ai    = rolmipvar(A,'A',N1,1);
Di    = rolmipvar(D,'D',N2,1);
Di    = fork(Di,'Dib');
B1i   = rolmipvar(B1,'B1',N1,1);
B2i   = rolmipvar(B2,'B2',N1,1);
C1i   = rolmipvar(C1,'C1',N1,1);
D21i  = rolmipvar(D21,'D21',N1,1);
D11i  = rolmipvar(D11,'D11',N1,1);
mu    = sdpvar();
%
% system variables
if deg(1)>0 & deg(2)>0
     Wa   = rolmipvar(order,order,'W','sym',[N1 N2],deg);
     dotW = diff(Wa,'dpoly',dotbounds);
elseif deg(1)>0
     Wa   = rolmipvar(order,order,'W','sym',N1,deg(1));
     dotW = diff(Wa,'dpoly',dotbounds{1});
     dotW = fork(dotW,'dpoly',[1 2],[1 3]);
elseif deg(2)>0
     Wa   = rolmipvar(order,order,'W','sym',N2,deg(2));
     dotW = diff(Wa,'dpoly',dotbounds{2});
     dotW = fork(dotW,'dpoly',[1 2],[1 4]);
else
     Wa = rolmipvar(order,order,'W','sym',[N1 N2],[0 0]);
     dotW = zeros(order);
end

X     = rolmipvar(order,order,'X','full',[0 0],[0 0]);
Z     = rolmipvar(l,order,'Z','full',[0 0],[0 0]);

Y1    = rolmipvar(order,order,'Y1','full',[N1 N2],deg);
Y2    = rolmipvar(order,order,'Y2','full',[N1 N2],deg);
Y3    = zeros(order,r);
Y4    = rolmipvar(order,order,'Y3','full',[N1 N2],deg);
Y5    = zeros(order,p);
Y6    = rolmipvar(order,order,'Y4','full',[N1 N2],deg);

T11  = Ai*X + X'*Ai' + B2i*Z + Z'*B2i' - dotW;
T12  = Wa - X' + xi*(Ai*X + B2i*Z);
T13  = X'*C1i' + Z'*D21i';
T14  = Wa - Y1';
T15  = B1i;
T16  = Y1'*Di';
T22  = -xi*(X + X');
T23  = xi*(X'*C1i' + Z'*D21i');
T24  = -Y2';
T25  = zeros(order,p);
T26  = Y2'*Di';
T33 = -mu*eye(r);
T34 = -Y3';
T35 = D11i;
T36 = Y3'*Di';
T44 = -(Y4 + Y4');
T45 = -Y5;
T46 = -Y6' + Y4'*Di';
T55 = -eye(p);
T56 = Y5'*Di';
T66 = -Wa + Di*Y6 + Y6'*Di';
T   = [T11  T12  T13 T14 T15 T16;
       T12' T22  T23 T24 T25 T26;
       T13' T23' T33 T34 T35 T36;
       T14' T24' T34' T44 T45 T46;
       T15' T25' T35' T45' T55 T56;
       T16' T26' T36' T46' T56' T66];

LMIs  = [T<=0,Wa>=0];

opt   = sdpsettings('verbose',0,'solver','sedumi');%mosek
sol   = solvesdp(LMIs,mu,opt);

output.tempo  = sol.solvertime;
output.p      = min(checkset(LMIs));
output.V      = size(getvariables(LMIs),2);
output.L      = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(LMIs(i),1);
end

delta=min(checkset(LMIs));
if delta > -1e-6
    output.feas = 1;
    output.Wab  = double(Wa);
    output.hinf = sqrt(double(mu));
    output.K    = double(Z)*inv(double(X));
    %fprintf('estavel!\n');
else
    output.feas = 0;
    %fprintf('instavel!\n');
end
%compile
