function [node_disp,Tensoes,Massa,Pcrit,SecForce] = FEM(area_anterior)
Dados.coordenadas = [0 0 0;0 9.144 0;9.144/2 9.144/2 0;9.144 0 0;9.144 9.144 0;9.144+9.144/2 9.144/2 0;9.144*2 0 0;9.144*2 9.144 0];
Dados.elementos = [1 4;1 3;2 3;2 5;3 5;4 5;3 4;4 7;4 6;5 6;5 8;6 8;7 8;6 7];
Dados.geometria = [area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior;area_anterior];
Dados.material = [68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6;68947*10^6];
Dados.gdlfixos = [1 2 3 4 5 6 9 12 15 18 21 24]';
Dados.cargapontual = [0 0 0 0 0 0 0 0 0 0 -444.82*10^3 0 0 0 0 0 0 0 0 -444.82*10^3 0 0 0 0]';
Dados.cargadistrib=[];

[node_disp,deslocamentos] = FEMSolver(Dados);
% drawingMesh(Dados,10,deslocamentos);
d = sqrt(4*area_anterior/pi);
I = (pi*d^4)/32; % Momento de Inércia [m^4] 
[Tensoes,Pcrit] = stresses3Dtruss(Dados,deslocamentos,I);
k = 1; % k para pino-pino (flambagem)
E = 68947*10^6;
rho = 2700; % Densidade [kg/m^3]
Lbar = 6*9.144 + 4*9.144*sqrt(2); % comprimento total de barras [m]
SecForce = Tensoes.*area_anterior;
Massa = rho*area_anterior*Lbar;  % massa total das barras [m^2]

end

