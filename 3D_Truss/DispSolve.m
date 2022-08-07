function [deslocamentos] = DispSolve(K,F,nndof,fixdegrees)

activeDof=setdiff([1:nndof]',[fixdegrees]);    % Verifica os GDL ativos
U=K(activeDof,activeDof)\F(activeDof);  % Deslocamentos
deslocamentos=zeros(nndof,1); % Cria a matriz de deslocamentos para todos GDL
deslocamentos(activeDof)=U; % Insere o valor dos deslocamentos nos GDL ativos

end