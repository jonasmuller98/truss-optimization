function [K] = StifMatrix(nndof,nelem,elements,coordinates,geometry,materials)
%% Função para construção da matriz de rigidez
K = zeros(nndof,nndof);

for i=1:nelem
    ielem = elements(i,1:2);  % Identifica os nós do elemento em questão
    
    igdl=[3*ielem(1)-2 3*ielem(1)-1 3*ielem(1)...
        3*ielem(2)-2 3*ielem(2)-1 3*ielem(2)]; % Vetor localização
    
    x1=coordinates(ielem(1),1);y1=coordinates(ielem(1),2);z1=coordinates(ielem(1),3);
    x2=coordinates(ielem(2),1);y2=coordinates(ielem(2),2);z2=coordinates(ielem(2),3);
    
    L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);   % Tamanho do elemento
    A = geometry(i);
    E = materials(i);
    
    Cx = (x2 - x1)/L;
    Cy = (y2 - y1)/L;
    Cz = (z2 - z1)/L;

    T = [Cx*Cx Cx*Cy Cx*Cz ;
        Cy*Cx Cy*Cy Cy*Cz ;
        Cz*Cx Cz*Cy Cz*Cz];
    
    K(igdl,igdl)=K(igdl,igdl)+E*A/L*[T -T ; -T T];    
end


end