function drawingMesh(Dados,scaleFact,displacements)
coordinates=Dados.coordenadas; elements=Dados.elementos;
nelem=size(elements,1);

for e=1:nelem
          indice=elements(e,1:2);
          local=[indice(1)*3-2 indice(1)*3-1 indice(1)*3 indice(2)*3-2 indice(2)*3-1 indice(2)*3]; % local: � o vetor de localiza��o de cada GDL na Matriz de Rigidez
          x=coordinates(indice,1); % Coordenadas em x dos n�s do elemento e.
          y=coordinates(indice,2); % Coordenadas em y dos n�s do elemento e.
          z=coordinates(indice,3); % Coordenadas em z dos n�s do elemento e.
          DeslElem=displacements(local); % Este vetor cont�m apenas os deslocamentos do elemento e no sistema global.
          DeslocX=DeslElem([1 4]);
          DeslocY=DeslElem([2 5]);
          DeslocZ=DeslElem([3 6]);
          x1=x+DeslocX*scaleFact; % Coordenadas em x dos n�s do elemento e ap�s deforma��o.
          y1=y+DeslocY*scaleFact; % Coordenadas em y dos n�s do elemento e ap�s deforma��o.
          z1=z+DeslocZ*scaleFact; % Coordenadas em y dos n�s do elemento e ap�s deforma��o.
          plot3(x,y,z,'--','linewidth',1.5);hold on;
          plot3(x1,y1,z1,'b','linewidth',1.5);hold on;
      end
      grid minor
      legend('Estrutura Original','Estrutura Deformada')
      hold off
      xlabel('X','FontSize',14);ylabel('Y','FontSize',14);zlabel('Z','FontSize',14);
      title('Representa��o da Estrutura','FontSize',16);
      drawnow
end

