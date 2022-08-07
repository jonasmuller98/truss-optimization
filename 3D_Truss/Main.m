%% Particle Swarm Optimization - Truss
% Jonas Müller Gonçalves
clear all;clc;
% rng(2);
tic
% nrodadas = 10;
% for rodar = 1:nrodadas

%% Inicialização das Variáveis
population = 10; npart=population;   % População total de Partículas
variables = 1; nvar=variables;     % Variáveis de Projeto
lambda1 = 2;    % Lambda 1
lambda2 = 2;    % Lambda 2
omega = 0.5;    % Propriedade de Inércia

dominio = [6.426e-5 0.045238];%.*ones(nvar,1); % Domínio das nvar variaveis
    
tol = 1e-8;     % Diferença entre o xgbest anterior e o novo xgbest - erro
erro_contador = 0;  % contador de iteracoes com o verificador de erro ativado
nerros = 600;   % Número máximo de iterações com o erro ativado
contador = 1;
%% Inicialização dos Vetores de Projeto
% Distribuindo as Partículas
for icol = 1:nvar
    for ivar = 1:npart
        position(ivar,icol) = abs(dominio(icol,1)+rand*(dominio(icol,1)-dominio(icol,2)));
    end
end
ct1=1;ct2=1;    % contadores auxiliares

% Checando se o chute inicial está dentro do domínio
for ipart=1:npart
    for ivariab=1:nvar
        if position(ipart,ivariab) < dominio(ivariab,1)
            fora_dominio_menor(ct1,ivariab) = position(ipart,ivariab); ct1=ct1+1;
            position(ipart,ivariab) = dominio(ivariab,1);  % Induz o retorno da partícula que passou
        elseif position(ipart,ivariab) >  dominio(ivariab,2)
            fora_dominio_maior(ct2,ivariab) = position(ipart,ivariab); ct2=ct2+1;
            position(ipart,ivariab) = dominio(ivariab,2);   % Induz o retorno da partícula que passou
        end
    end
end

% Definindo as Velocidades Iniciais
velocity = zeros(npart,nvar);   % inicializa as velocidades

%% Cálculo Inicial da Melhor Partícula
for ires=1:size(position,1)
    [node_disp,Tensoes,Massa,SecForce,Pcrit] = FEM(position(ires,:));    % Define os valores de f(x,y)
    [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit); % Verifica as condições de contorno
    contador = contador +1;    % Contador de Chamada da Função
    resultados(ires,1) = Massa * fp1 * fp2 * fp3;   % Imposição da Restrição
end

[resg,ord]=sort(resultados);    % Ordena os valores do menor para o maior

for ixgbest=1:nvar
    xgbest(ixgbest) = [position(ord(1,1),ixgbest)];   % Define o melhor global como o primeiro do vetor ordenado
end
melhorglobal=xgbest;    % Cria um valor armazenado do melhor global
ordglobal = ord(1,1);   % Replica qual a linha do melhor global

xlbest = 0.04523.*ones(npart,nvar); % Inicia Vetor de melhor local

% Variáveis Auxiliares
parada=0; iter=0; pos_old=position; v_old=velocity; cont1=1; cont2=1;


%% Otimização
while parada == 0
    iter = iter+1; % Atualiza o Contador de Iterações
    r1=rand; r2=rand; % R's aleatórios
    
    % Cálculo das Novas Velocidades e Posições de cada partícula
    for i=1:npart
        for j=1:nvar
            v(i,j) = omega*v_old(i)+lambda1*r1*(xlbest(i,j) - pos_old(i,j)) + lambda2*r2*(xgbest(j) - pos_old(i,j));
            pos(i,j) = abs(pos_old(i,j) + v(i,j));
        end
    end
    cont1=1;cont2=1;    % Contadores Auxiliares
    
    % Verificação do Domínio das Novas Partículas
    for ipart=1:npart
        for ivariab=1:nvar
            if pos(ipart,ivariab) < dominio(ivariab,1)
                fora_dominio_menor(ct1,ivariab) = pos(ipart,ivariab); 
                ct1=ct1+1;    % armazena os valores que foram menores que o dominio (somente para conferencia)
                pos(ipart,ivariab) = dominio(ivariab,1);  % Induz o retorno da partícula que passou
            elseif pos(ipart,ivariab) >  dominio(ivariab,2)
                fora_dominio_maior(ct2,ivariab) = pos(ipart,ivariab); 
                ct2=ct2+1;    % armazena os valores que foram maiores que o dominio (somente para conferencia)
                pos(ipart,ivariab) = dominio(ivariab,2);  % Induz o retorno da partícula que passou
            end
        end
    end
    
% Atualiza o valor de xlbest
for j=1:nvar
    for i=1:npart
        [node_disp,Tensoes,fxlbest2,SecForce,Pcrit] = FEM(xlbest(i,:));
        [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit);
        fxlbest = fxlbest2 * fp1 * fp2 * fp3;
        contador = contador +1;

        [node_disp,Tensoes,fposi2,SecForce,Pcrit] = FEM(pos(i,:));
        [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit);
        fposi = fposi2 * fp1 * fp2 * fp3;
        contador = contador +1;
            
        if fxlbest > fposi   % Define o melhor global como o primeiro do vetor ordenado
            xlbest(i,:) = pos(i,:);
        else
            
        end
    end
end

% Armazena os resultados de f(xlbest)
    for ilc=1:npart
        [node_disp,Tensoes,fxlbest3,SecForce,Pcrit] = FEM(xlbest(ilc,:));
        [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit);
        fxlbest_aux = fxlbest3 * fp1 * fp2 * fp3;
        resultados_loc(ilc,1) = fxlbest_aux;
        contador = contador +1;
    end
    
    % Ordena os valores de f(xlbest) do menor para o maior
    [rloc,iloc] = sort(resultados_loc);
    
    % Armazena o melhor x local
    for imelhorl=1:nvar
        melhorlocal(imelhorl) = [xlbest(iloc(1),imelhorl)];
    end
    
    [node_disp,Tensoes,fmelhorglobal2,SecForce,Pcrit] = FEM(melhorglobal(1,:));
    [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit);
    fmelhorglobal = fmelhorglobal2 * fp1 * fp2 * fp3;
    contador = contador +1;

    [node_disp,Tensoes,fmelhorlocal2,SecForce,Pcrit] = FEM(melhorlocal(1,:));
    [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit);
    fmelhorlocal = fmelhorlocal2 * fp1 * fp2 * fp3;
    contador = contador +1;

    % Confere se o melhor x local é melhor que o melhor x global
    if fmelhorglobal > fmelhorlocal
        melhorglobal(1,:) = melhorlocal(1,:);
    end
    
    [node_disp,Tensoes,fxgbest2,SecForce,Pcrit] = FEM(xgbest(1,:));
    contador = contador +1;
    [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit);
    fxgbest = fxgbest2 * fp1 * fp2 * fp3;
   
    melhor_disp = node_disp;
    melhor_tensao = Tensoes;    
    melhor_massa = fxgbest2;
    melhor_area = xgbest;
    massas_conv(iter,1) = fxgbest2;
    
    erro = fxgbest - fmelhorglobal;   % Verifica o erro
    xgbest = melhorglobal;    % Atualiza o xgbest 
    v_old = v;  % Atribui a velocidade antiga como a nova velocidade calculada
    pos_old = pos;  % Atribui a posição antiga como a nova posição calculada (já com o domínio verificado)
    melhores(iter,:) = melhorglobal; % Armazena todos melhores globais
    
    % Condição de Parada
    if erro<=tol % se o erro calculado é menor que a tolerancia, atualiza o contador
        erro_contador = erro_contador+1;
        if erro_contador>nerros % se o contador atinge o numero maximo de iteracoes consecutivas com erro<=tol
            parada=1;
            
        end
    else 
        erro_contador=0;
    end
    
end
%     MASSAS(rodar,1) = melhor_massa;
% end
tempo = toc;

%% Pós Processamento

% Para fazer essa figura, ativar o laço mais externo - linha 6,7,180,181
% figure(1)
% plot(1:nrodadas,MASSAS,'-o','linewidth',1.4); grid minor;title('Convergência da Simulação','FontSize',17);xlabel('Nº de Simulações','FontSize',16);ylabel('Massa [kg]','FontSize',16)
% ytickformat('%.2f');

figure(2);
plot(1:iter,massas_conv,'linewidth',1.4); grid minor;title('Convergência de Resposta','FontSize',17);xlabel('Nº de Iterações','FontSize',16);ylabel('Massa [kg]','FontSize',16)
ytickformat('%.2f');
 %% Verificação dos Limites:

disp(' --- Verifição dos Limites --- ');

disp('Área: 6.426e-5 < A < 0.045238');
disp(['Melhor Área Obtida: [',num2str(melhor_area),']']);
disp(' ');

disp('Desloc Nodal: -0.058 < u < 0.058');
disp(['Max Obtido: [',num2str(max(melhor_disp(:,1))),' m]; Min Obtido: [',num2str(min(melhor_disp(:,1))),' m]']);
disp(' ');

disp('Tensão: -173.36e6 < sigma < 173.36e6');
disp(['Max Obtida: [',num2str(max(melhor_tensao)),' Pa]; Min Obtida: [',num2str(min(melhor_tensao)),' Pa]']);
disp(' ');

disp(['Melhor Massa Obtida: [',num2str(melhor_massa),' kg]']);
disp(' ');

disp(' --- Usando a melhor área no FEM --- ')
[node_disp,Tensoes,Massa,Pcrit,SecForce] = FEM(melhor_area);
disp(['Massa utilizando a melhor area: [',num2str(Massa),' kg]']);
disp(['Max tensao utilizando a melhor area: [',num2str(max(Tensoes)),' Pa]; Min: [',num2str(min(Tensoes)),' Pa]']);
disp(['Max disp utilizando a melhor area: [',num2str(max(node_disp(:,1))),' m]; Min: [',num2str(min(node_disp(:,1))),' m]']);
disp(' ');

disp(' --- Dados --- ');
disp(['Número de chamadas da função FEM: [',num2str(contador),']']);
disp(['Número de Iterações até a convergência: [',num2str(iter),']']);
disp(['Tempo de Solução: [',num2str(tempo),' s]']);
