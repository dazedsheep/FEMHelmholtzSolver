function [utt] = dtt(u, timeMeshh)
utt = zeros(size(u));
utt = [ (u(3,:) - 2*u(2,:) + u(1,:)) / timeMeshh^2;
          (u(3:end,:) - 2*u(2:end-1,:) + u(1:end-2,:)) / timeMeshh^2;
          (u(end,:) - 2*u(end-1,:) + u(end-2,:)) / timeMeshh^2 ];
end

