function out = Lugar_Raizes(n,N,Fronteira,aleatorios,d)


% Matrizes do sistema
for i=1:N
   A{i} = rand(n); 
end

eig_Fronteira= [];
eig_Rand     = [];

% Fronteira do politopo
for i=1:N
    for j=i+1:N
        for k=0:1/Fronteira:1
            alphas      = zeros(1,N); 
            alphas(1,i) = k;
            alphas(1,j) = 1-k;
            Ai          = alphas(1,i)*A{i} + alphas(1,j)*A{j};
            eig_max     = max(real(eig(Ai)));
            
            
            Ai            = Ai -(eig_max+d)*eye(n);
            eig_Fronteira = [eig_Fronteira;eig(Ai)];
        end
    end
end


% Interior Politopo
for i=1:aleatorios
    alphas = rand(1,N);
    alphas = alphas./sum(alphas);
    Ai     = zeros(n);
    for j=1:N
        Ai       = Ai + A{j}*alphas(1,j);
        eig_max  = max(real(eig(Ai)));
        Ai       = Ai -(eig_max+d)*eye(n);
    end
    m_i      = max(real(eig(Ai)));
    eig_Rand = [eig_Rand;eig(Ai)];
end

plot(real(eig_Rand),imag(eig_Rand),'o');
grid
