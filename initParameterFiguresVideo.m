function [plot_handles]  = initParameterFiguresVideo(elements, xn, NewtonIterations, CGIterations, plot_handles)

set(0, 'CurrentFigure', plot_handles.video.fig);

plot_handles.video.plot_b = subplot(2,3,4,'Parent',plot_handles.video.fig);
plot_handles.video.b_axes = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_2.b0, 'facecolor', 'interp','Parent', plot_handles.video.plot_b); 
str = sprintf('Recon. b (Newton: %d, CG: %d)', NewtonIterations, CGIterations);
title(str);
view(0,90)  
colorbar
shading interp;

plot_handles.video.plot_s = subplot(2,3,5,'Parent',plot_handles.video.fig);
plot_handles.video.s_axes = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),xn.refState_2.s0, 'facecolor', 'interp', 'Parent', plot_handles.video.plot_s); 
title('Reconstructed s');
view(0,90)  
colorbar
shading interp;

plot_handles.video.plot_eta = subplot(2,3,6,'Parent',plot_handles.video.fig);
plot_handles.video.eta_axes = trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2),projectParameters(xn).refState_2.eta0, 'facecolor', 'interp', 'Parent', plot_handles.video.plot_eta); 
title('Reconstructed \eta');
view(0,90)  
colorbar
shading interp;

end

