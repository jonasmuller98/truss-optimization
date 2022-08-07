function [fp1,fp2,fp3] = Verifica(node_disp,Tensoes,SecForce,Pcrit)

for inode = 1:size(node_disp,1)
    if node_disp(inode,1) > 0.058 || node_disp(inode,1) < -0.058
        fp1 = exp(abs(node_disp(inode,1)));
    elseif node_disp(inode,2) > 0.058 || node_disp(inode,2) < -0.058
        fp1 = exp(abs(node_disp(inode,2)));
    else
        fp1=1;
    end
end

for itensao = 1:size(Tensoes,1)
    if Tensoes(itensao) > 173.36e6 || Tensoes(itensao) < -173.36e6
        fp2 = exp(abs(Tensoes(itensao)));
    else
        fp2=1;
    end
end

for ipcrit = 1:3:size(SecForce,1)
    if SecForce(ipcrit) < 0    % condição de compressao
        if abs(SecForce(ipcrit)) > Pcrit
            fp3 = exp(abs(SecForce(ipcrit)));
        else
            fp3=1;
        end
    else
        fp3=1;
    end
end


end
