function [node_disp,deslocamentos,F] = FEMSolver(Dados)
coordinates=Dados.coordenadas; elements=Dados.elementos; geometry=Dados.geometria;
materials=Dados.material; fixdegrees=Dados.gdlfixos; 
pointload=Dados.cargapontual;
distload=Dados.cargadistrib;

%% Verificação das Dimensões
nelem = size(elements,1);
nnode = size(coordinates,1);
nndof = 3*nnode;

%% Rigidez
[K] = StifMatrix(nndof,nelem,elements,coordinates,geometry,materials);

%% Forças
[F] = ForceVector(pointload,distload,coordinates);
[deslocamentos] = DispSolve(K,F,nndof,fixdegrees);

Fr=K*deslocamentos;
reactions = Fr(fixdegrees);

%% Segregando os Deslocamentos nos Respectivos Nós
d_nd = length(deslocamentos)/3;
i=1;j=2;k=3;

for ie = 1:d_nd
    Disp_X(ie,1) = deslocamentos(i);
    Disp_Y(ie,1) = deslocamentos(j);
    Disp_Z(ie,1) = deslocamentos(k);
    i=i+3; j=j+3; k=k+3;
end

node_disp = [Disp_X Disp_Y Disp_Z];
end