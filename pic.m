function pic(output)
    close all

    fig1 = figure(1)
    h1 = plot3(output(:,2),output(:,3),output(:,4));
    hold on
    h2 = plot3(output(:,5),output(:,6),output(:,7));
    legend([h1, h2], {'p', 'pr'});

    fig2 = figure(2)
    subplot(3,1,1)
    p1 = plot(output(:,1),output(:,2));
    hold on
    p2 = plot(output(:,1),output(:,5));
    legend([p1,p2],{'x','xr'});
    
    subplot(3,1,2)
    p1 = plot(output(:,1),output(:,3));
    hold on
    p2 = plot(output(:,1),output(:,6));
    legend([p1,p2],{'y','yr'});

    subplot(3,1,3)
    p1 = plot(output(:,1),output(:,4));
    hold on
    p2 = plot(output(:,1),output(:,7));
    legend([p1,p2],{'z','zr'});

    frame1 = getframe(fig1); 
    frame2 = getframe(fig2); 
    img1 = frame2im(frame1); 
    img2 = frame2im(frame2);
    imwrite(img1,'./imgs/Figure of 3D trajectory tracking.png');
    imwrite(img2,'./imgs/Figure of position tracking.png');
end