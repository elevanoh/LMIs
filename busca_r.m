clc;
k1        = [1.6 2.4]';
k2        = 9.6; % mudei aqui
beta      = 1;
d         = [0.1 0.45]; %mudei aqui
for i=1:2
    D{i}  = blkdiag(d(i),-d(i),0,0);
end
ind = 0;
for i=1:2
    ind     = ind + 1;
    A{ind}  = beta*[ 0.0      0.0              1.0   0.000;
        0.0      0.0              0.0   1.000;
        -k1(i)   k1(i)            -0.2  0.200;
        k1(i)/2  -(k1(i)+k2)/2 0.1   -0.15];
    B1{ind} = [0 0 0 0.5]';    
    B2{ind} = [0 0 1 0]';    
    C1{ind} = [1 0 0 0;0 0 0 0];
    D11{ind}= [0 0]';
    D21{ind}= [0 1]';
end

bestGamma = 1e10;
xis         = logspace(-2,0.89,50);
r = [0.1 1 5 10];
Hinf2      = [];
for g=1:1
    for i=1:size(r,2)
        y         = [];
        for xi=xis,
            dotbounds = {r(i)*[-1 1;-1 1],r(i)*[-1 1;-1 1]};
            deg       = [g 0];
            output = sintesis_robusto_LPV_sys_continuo(A,D,B1,B2,C1,D11,D21,dotbounds,xi,deg);
            if output.feas == 0
                output.hinf = 0;
            end
            bestGamma=output.hinf;
            bestXi   = r(i);
            fprintf('ri=%.2f xi=%.2f -> gamma=%.4f\n',r(i),xi,output.hinf);
            y        = [y;output.hinf];
        end
        Hinf2 = [Hinf2 y];
    end
end

semilogx(xis,Hinf2);
grid;
xlabel xi
ylabel \gamma
legend('r=0.1','r=1','r=5','r=10');
hold off
%print('-depsc2','Fig_IFAC2020_Exemplo_01')
