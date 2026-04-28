function [plot_handles] = initParameterFigures(elements, recordVideo, xn)
figure,
plot_handles.plot_b_1 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_1.b0, 'facecolor', 'interp'); 
title('Reconstructed b');
view(0,90)  
colorbar
shading interp;
figure,
plot_handles.plot_s_1 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_1.s0, 'facecolor', 'interp'); 
title('Reconstructed s');
view(0,90)  
colorbar
shading interp;
figure,
plot_handles.plot_eta_1 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_1.eta0, 'facecolor', 'interp'); 
title('Reconstructed \eta');
view(0,90)  
colorbar
shading interp;

figure,
plot_handles.plot_b_2 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_2.b0, 'facecolor', 'interp'); 
title('Reconstructed b');
view(0,90)  
colorbar
shading interp;
figure,
plot_handles.plot_s_2 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_2.s0, 'facecolor', 'interp'); 
title('Reconstructed s');
view(0,90)  
colorbar
shading interp;
figure,
plot_handles.plot_eta_2 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_2.eta0, 'facecolor', 'interp'); 
title('Reconstructed \eta');
view(0,90)  
colorbar
shading interp;

figure,
plot_handles.plot_b_3 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_3.b0, 'facecolor', 'interp'); 
title('Reconstructed b');
view(0,90)  
colorbar
shading interp;
figure,
plot_handles.plot_s_3 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_3.s0, 'facecolor', 'interp'); 
title('Reconstructed s');
view(0,90)  
colorbar
shading interp;
figure,
plot_handles.plot_eta_3 = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_3.eta0, 'facecolor', 'interp'); 
title('Reconstructed \eta');
view(0,90)  
colorbar
shading interp;

if recordVideo == true
    % write each parameter into a separate video 

end


end

