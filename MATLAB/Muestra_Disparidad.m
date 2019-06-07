function Muestra_Disparidad(newPoints1,newPoints2,umbral,img)

    % Muestra la disparidad
    for i=1:length(newPoints1)
        C(i) = newPoints1(i,1);
        D(i) = newPoints2(i,1);
        A(i) = newPoints1(i,2);
        B(i) = newPoints2(i,2);
    end

    % Disparidad
    figure('Position', [100 100 size(img,2) size(img,1)]);
    colormap('gray');
    imagesc(img);
    title(['Umbral: ',num2str(umbral)])
    for i=1:2:length(newPoints1)
        hold on
        plot([C(i) D(i)],[A(i) B(i)],'-b')
    end
    axis([0 2000 0 1500])
end

