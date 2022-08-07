function [F] = ForceVector(pointload,sideload,coordinates)
F=pointload;

for i = 1 : size(sideload,1)
    indice=[sideload(i,1) sideload(i,2)];
    LM=[3*indice(1)-2 3*indice(1)-1 3*indice(1)...
        3*indice(2)-2 3*indice(2)-1 3*indice(2)]; % Vetor localização
    
    x1=coordinates(sideload(i,1),1);
    y1=coordinates(sideload(i,1),2);
    z1=coordinates(sideload(i,1),3);
    x2=coordinates(sideload(i,2),1);
    y2=coordinates(sideload(i,2),2);
    z2=coordinates(sideload(i,2),3);
    l = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1)); % Compr. do elemento
    
    CXx = (x2-x1)/l;
    CYx = (y2-y1)/l;
    CZx = (z2-z1)/l;
    
    R = [CXx CYx CZx zeros(1,3); zeros(1,3) CXx CYx CZx]; % Matriz de transformação de coordenadas
    
    force=l/6*[2 1; 1 2]*[sideload(i,3); sideload(i,3)];
    
    F(LM)=F(LM) + R'*[force];
end

end