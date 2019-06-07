
% This function reads two images, finds their SIFT features, and
%   displays lines connecting the matched keypoints.  A match is accepted
%   only if its distance is less than distRatio times the distance to the
%   second closest match.
%


function [puntosMatch,im1,im2] = match(umbral,representaImagen, image1, image2)

%umbral: Umbral para la correspondencia (distRatio)
%representaImagen:
%   1: Representa la primera imagen con los puntos encontrados
%   2: Representa la segunda imagen con los puntos encontrados
%   3: Representa ambas imagenes con una linea que une los puntos casados
%image1: Nombre de la primera imagen en formato .pgm
%image1: Nombre de la segunda imagen en formato .pgm
%puntosMatch: Puntos casados. En cada fila x1,y1,x2,y2
%im1: Matriz imagen1
%im2: Matriz imagen2


% Find SIFT keypoints for each image
[im1, des1, loc1] = sift(image1);
[im2, des2, loc2] = sift(image2);

% For efficiency in Matlab, it is cheaper to compute dot products between
%  unit vectors rather than Euclidean distances.  Note that the ratio of 
%  angles (acos of dot products of unit vectors) is a close approximation
%  to the ratio of Euclidean distances for small angles.
%
% distRatio: Only keep matches in which the ratio of vector angles from the
%   nearest to second nearest neighbor is less than distRatio.
distRatio = umbral;   

% For each descriptor in the first image, select its match to second image.
des2t = des2';                          % Precompute matrix transpose
for i = 1 : size(des1,1)
   dotprods = des1(i,:) * des2t;        % Computes vector of dot products
   [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

   % Check if nearest neighbor has angle less than distRatio times 2nd.
   if (vals(1) < distRatio * vals(2))
      match(i) = indx(1);
   else
      match(i) = 0;
   end
end


%GENERA LA MATRIZ DE RESULTADOS
%Por filas los puntos casados
%Por columnas (x,y) primera imagen, (x,y) segunda imagen

%Numero de puntos casados
num = sum(match > 0);
puntosMatch=zeros(num,4);
indice=1;

%Bucle para cada caracteristica encontrada
for i = 1: size(des1,1)
    %Caracteristica casada
    if (match(i) > 0)
        puntosMatch(indice,1)=loc1(i,1);
        puntosMatch(indice,2)=loc1(i,2);
        puntosMatch(indice,3)=loc2(match(i),1);
        puntosMatch(indice,4)=loc2(match(i),2);
        indice=indice+1;
    end
end

    
%REPRESENTA LA IMAGEN


% Ambas imagenes con correspondencias
if (representaImagen==3)
    
    % Create a new image showing the two images side by side.
    % Select the image with the fewest rows and fill in enough empty rows
    %   to make it the same height as the other image.
    rows1 = size(im1,1);
    rows2 = size(im2,1);

    if (rows1 < rows2)
        im1(rows2,1) = 0;
    else
        im2(rows1,1) = 0;
    end

    % Now append both images side-by-side.
    im3 = [im1 im2];
    
    % Show a figure with lines joining the accepted matches.
    figure('Position', [100 100 size(im3,2) size(im3,1)]);
    colormap('gray');
    imagesc(im3);
    hold on;
    cols1 = size(im1,2);
    for i = 1: size(des1,1)
        if (match(i) > 0)
            line([loc1(i,1) loc2(match(i),1)+cols1], ...
                [loc1(i,2) loc2(match(i),2)], 'Color', 'c');
        end
    end
    hold off;
end

%Una sola imagen con puntos detectados
if ((representaImagen==1)||(representaImagen==2))
    
    if (representaImagen==1)
        im3=im1;
    else
        im3=im2;
    end
     
    % Show a figure with lines joining the accepted matches.
    figure('Position', [100 100 size(im3,2) size(im3,1)]);
    colormap('gray');
    imagesc(im3);
    hold on;
    
    for i = 1: size(des1,1)
        if (match(i) > 0)
            if (representaImagen==1)
                posX=loc1(i,1);
                posY=loc1(i,2);
            else
                posX=loc2(match(i),1);
                posY=loc2(match(i),2);
            end           
            
            for (j=-2:2)
                line([posX-2 posX+2],[posY+j posY+j], 'Color', 'c');
            end
        end
    end
    hold off;
end

end







