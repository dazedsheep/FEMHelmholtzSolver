function updateVideoFigure(s, b, eta, plot_handles, NewtonIterations, CGIterations)

set(plot_handles.video.b_axes, 'CData', b);
set(plot_handles.video.s_axes, 'CData', s);
set(plot_handles.video.eta_axes, 'CData', eta);

str = sprintf('Recon. b (Newton: %d, CG: %d)', NewtonIterations, CGIterations);
title(plot_handles.video.plot_b, str);

drawnow;

writeVideo(plot_handles.video.video_handle, getframe(plot_handles.video.fig));

end

