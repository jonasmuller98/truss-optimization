function [Tensoes,Pcrit] = stresses3Dtruss(Dados,deslocamentos,I)

coordinates=Dados.coordenadas; elements=Dados.elementos; geometry=Dados.geometria;
materials=Dados.material; fixdegrees=Dados.gdlfixos; 
nelem = size(elements,1);
pointload=Dados.cargapontual; distload=Dados.cargadistrib;

% fprintf('Tensão nos elementos\n')
ff=zeros(nelem,6); 
% format long
k = 1;

for e=1:nelem
    indice=elements(e,:) ;
    elementDof=[3*indice(1)-2 3*indice(1)-1 3*indice(1)...
        3*indice(2)-2 3*indice(2)-1 3*indice(2)] ;
    x1=coordinates(indice(1),1);
    y1=coordinates(indice(1),2);
    z1=coordinates(indice(1),3);
    x2=coordinates(indice(2),1);
    y2=coordinates(indice(2),2);
    z2=coordinates(indice(2),3);
    L = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
    CXx = (x2-x1)/L; CYx = (y2-y1)/L; CZx = (z2-z1)/L;
    u=deslocamentos(elementDof);
    E = materials(e,1);
    member_stress(e)=E/L*[-CXx -CYx -CZx CXx CYx CZx]*u;
%     fprintf('%3d %12.8f\n',e, member_stress(e));
    Pcrit(e) = (pi^2*E*I)/((k*L)^2);
end

Tensoes = member_stress;
