function [modDomain, modBoundary, obs] = forwardOperatorWestervelt( ...
    elements, observation, timeMeshh, L, M, u, eta, b, s, gamma)

Nt = size(u,1);
Np = size(u,2);

%% Observation
obs = u(:,observation);

%% Time derivatives (vectorized, no circshift)
u_t  = diff(u,1,1) ./ timeMeshh;
u_t  = [u_t; u_t(end,:)];          % Neumann extrapolation

uh   = b.'.*u - eta.'.*u.^2;
uh_t = diff(uh,1,1) ./ timeMeshh;
uh_t = [uh_t; uh_t(end,:)];

uh_tt = diff(uh_t,1,1) ./ timeMeshh;
uh_tt = [uh_tt; uh_tt(end,:)];

%% Laplacian (factorize once if reused often)
A = M\L;
Laplaceu   = (A * u.').';
Laplaceu_t = (A * u_t.').';

%% Domain residual
modDomain = zeros(Nt,Np);
modfull = uh_tt - s.' .* Laplaceu - Laplaceu_t;
modDomain(:,elements.interiorIdx) = modfull(:,elements.interiorIdx);

%% --------- FAST BOUNDARY GRADIENT ---------
persistent Gx Gy

if isempty(Gx)

    p = elements.points;
    t = elements.tri(:,1:3);

    x = p(:,1); 
    y = p(:,2);

    nt = size(t,1);
    np = size(p,1);

    % triangle vertex indices
    i1 = t(:,1); i2 = t(:,2); i3 = t(:,3);

    % coordinates
    x1 = x(i1); y1 = y(i1);
    x2 = x(i2); y2 = y(i2);
    x3 = x(i3); y3 = y(i3);

    % triangle areas
    area = 0.5*((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));

    % gradients of basis functions
    b1 = (y2-y3)./(2*area);
    b2 = (y3-y1)./(2*area);
    b3 = (y1-y2)./(2*area);

    c1 = (x3-x2)./(2*area);
    c2 = (x1-x3)./(2*area);
    c3 = (x2-x1)./(2*area);

    % assemble sparse gradient matrices
    rows = repmat((1:nt)',1,3);

    Gx_elem = sparse(rows(:), t(:), [b1;b2;b3], nt, np);
    Gy_elem = sparse(rows(:), t(:), [c1;c2;c3], nt, np);

    % average to nodes (mass lumping style)
    Mmap = sparse(repmat((1:nt)',3,1), t(:), 1, nt, np);
    weights = Mmap' * abs(area);

    Gx = (Mmap' * (abs(area).*Gx_elem)) ./ weights;
    Gy = (Mmap' * (abs(area).*Gy_elem)) ./ weights;

end

% compute gradient for ALL time steps at once
Ux = (Gx * u.').';
Uy = (Gy * u.').';

nb = elements.boundaryIdx;

normal_x = elements.boundaryNormals(:,1).';
normal_y = elements.boundaryNormals(:,2).';

normalGrad = Ux(:,nb).*normal_x + Uy(:,nb).*normal_y;

modBoundary = zeros(Nt,Np);
modBoundary(:,nb) = gamma.*u(:,nb) + normalGrad;

end