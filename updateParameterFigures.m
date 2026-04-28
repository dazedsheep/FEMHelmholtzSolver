function updateParameterFigures(plot_handles,xn)

set(plot_handles.plot_b_1, 'CData', xn.refState_1.b0);
set(plot_handles.plot_s_1, 'CData', xn.refState_1.s0);
set(plot_handles.plot_eta_1, 'CData', xn.refState_1.eta0);

set(plot_handles.plot_b_2, 'CData', xn.refState_2.b0);
set(plot_handles.plot_s_2, 'CData', xn.refState_2.s0);
set(plot_handles.plot_eta_2, 'CData', xn.refState_2.eta0);

set(plot_handles.plot_b_3, 'CData', xn.refState_3.b0);
set(plot_handles.plot_s_3, 'CData', xn.refState_3.s0);
set(plot_handles.plot_eta_3, 'CData', xn.refState_3.eta0);

drawnow;
pause(0.1);


end

