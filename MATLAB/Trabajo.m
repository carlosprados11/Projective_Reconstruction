clc
clear

%% Correspondences to obtain the fundamental matrix

%Dimensiones imagenes
NUM_FIL=1500;
NUM_COL=2000;

%Parametros
umbral=0.45;
representaImagen=1;
cmap = colormap('jet');

[puntosMatch,im1,im2]=match(umbral,representaImagen,'Escena1_Imagen1.pgm','Escena1_Imagen3.pgm');
num=size(puntosMatch,1);
fprintf('Puntos casados %d \n', num);


%% Filter trougthout RANSAC

% Calculo de la matriz fundamental (status 0 if not error)
[F,inliernsIndex,status] = estimateFundamentalMatrix(...
    puntosMatch(:,[1 2]),puntosMatch(:,[3 4]),'Method',...
    'RANSAC','NumTrials',50,'DistanceThreshold',2,...
    'DistanceType','Sampson','Confidence',95);

% Show the disparity
k=1;
for j=1:length(inliernsIndex)
    if inliernsIndex(j)==1
        newPoints1(k,:) = puntosMatch(j,[1 2]);
        newPoints2(k,:) = puntosMatch(j,[3 4]);
        k = k + 1;
    end
end
Muestra_Disparidad(newPoints1,newPoints2,umbral,im1);

% Percentage of inliners
inl = 0;
for j=1:length(inliernsIndex)
    if inliernsIndex(j)==1
        inl=inl+1;
    end
    error(j) = [puntosMatch(j,[3 4]) 1]*F*[puntosMatch(j,[1 2]) 1]';
end

inl = 100*inl/length(inliernsIndex);
fprintf('Porcentaje de inliners: %2.2f \n', inl);

% Median and variance error
M = mean(abs(error'));
V = var(abs(error'));

%% Filter troughtout measure error

rech_ant = 100;
rech = 2;
rango = 0;

% All of the matched points are used inicially
for i=1:length(puntosMatch)
    inlierns_sin(i) = 1;
end

for j=1:length(inlierns_sin)
    newPoints1(k,:) = puntosMatch(j,[1 2]);
    newPoints2(k,:) = puntosMatch(j,[3 4]);
end

% While there are eliminate samples
while rech_ant>rango
    
    % Calculate the fundamental matrix without noise elimination
    F_sin = estimateFundamentalMatrix(newPoints1,newPoints2);

    % Filter of error value
    rech = 0;
    
    % Calculate the error
    clearvars error
    ite = 1;
    for j=1:length(inlierns_sin)
        if inlierns_sin(j)==1
            error(ite) = [puntosMatch(j,[1 2]) 1]*F_sin*[puntosMatch(j,[3 4]) 1]';
            ite = ite + 1;
        end
    end
    
    M = mean(abs(error'));
    
    ite = 1;
    for j=1:length(inlierns_sin)
        if inlierns_sin(j)==1

            err = error(ite);
            ite = ite + 1;
            if abs(err)>2.5*M
                % Suprime value
                rech = rech + 1;
                inlierns_sin(j) = 0;
            end

        end
    end
    
    % Choose the correct values
    k=1;
    clearvars newPoints1 newPoints2
    for j=1:length(inlierns_sin)
        if inlierns_sin(j)==1
            newPoints1(k,:) = puntosMatch(j,[1 2]);
            newPoints2(k,:) = puntosMatch(j,[3 4]);
            k = k + 1;
        end
    end
    
    rech_ant = rech;
    
end

% Show the new disparity
Muestra_Disparidad(newPoints1,newPoints2,umbral,im1);

% New Fundamental Matrix
Fn = estimateFundamentalMatrix(newPoints1,newPoints2);


%% Correspondences to obtain the proyective reconstruction

% Calculate the epipole
ey = (F(1,2)*F(3,1)-F(3,2)*F(1,1))/(F(2,2)*F(1,1)-F(1,2)*F(2,1));
ex = -(F(3,1)+F(2,1)*ey)/F(1,1);

%Parametros
umbral=0.65;
representaImagen=1;

% Points linked
[puntosMatch,im1,im2]=match(umbral,representaImagen,'Escena1_Imagen1.pgm','Escena1_Imagen3.pgm');
num=size(puntosMatch,1);
fprintf('Puntos casados %d \n', num);

% RANSAC for eliminate wrong correspondences
[F,inliernsIndex,status] = estimateFundamentalMatrix(...
    puntosMatch(:,[1 2]),puntosMatch(:,[3 4]),'Method',...
    'RANSAC','NumTrials',50,'DistanceThreshold',2,...
    'DistanceType','Sampson','Confidence',95);

% Calculate the depth and position
for i=1:length(puntosMatch)
    x(i) = -1;
    y(i) = -1;
    z(i) = -1;
    if inliernsIndex(i)==1
        dist_x = ex + (puntosMatch(i,3)-puntosMatch(i,1));
        dist_y = ey + (puntosMatch(i,4)-puntosMatch(i,2));
        z(i) = sqrt((dist_x)^2+(dist_y)^2);
        x(i) = puntosMatch(i,1);
        y(i) = puntosMatch(i,2);
    end
end

%% Representate in 2D the image

% Select the min and max values of depth
min = 1000000;
max = -1000;
for i = 1:length(puntosMatch)
    if inliernsIndex(i)==1
        if z(i)<min
            min = z(i);
        elseif z(i)>max
            max = z(i);
        end
    end
end

% Assignment of a color for each depth
for i = 1:length(puntosMatch)
    if inliernsIndex(i)==1
        depth(i) = (z(i)-min)/(max-min);
        depth(i) = round(depth(i)*length(cmap));
        if depth(i)==0
            depth(i) = 1;
        end
            
    end
end
    
% Plotting the image
figure('Position', [100 100 size(im1,2) size(im1,1)]);
colormap('gray');
imagesc(im1);

hold on;
for i = 1:length(puntosMatch)
    if inliernsIndex(i)==1
        
        posX = puntosMatch(i,1);
        posY = puntosMatch(i,2);         

        for j = -3:3
            if ((posX-3)>0)&&((posX+3)<2001)&&((posY+j)>0)&&((posX+j)<1501)
                line([posX-3 posX+3],[posY+j posY+j],'Color',cmap(depth(i),:));
            end
        end
    end
end
hold off;

%% Representate in 3D the image

% Plotting
figure
hold on;
for i = 1:length(puntosMatch)
    if inliernsIndex(i)==1
        scatter3(x(i),y(i),z(i),'MarkerEdgeColor','k','MarkerFaceColor',cmap(depth(i),:))
    end
end
hold off;

