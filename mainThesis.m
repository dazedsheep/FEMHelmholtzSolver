
%% point sources
clearvars
% speef of sounde
speed_of_sound = 1540; % m/s


% signal period or cener frequency
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [1/2,1/2];
brad = 1/2;
domain = [bcenter, brad];


% point scatterers and their domain
values = [20000, 22000];
radii = [0, 0];
centers = [3/4, 3/4; 1/4, 3/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound, we set b=0, neglecting the diffusitivity of sound
kappa = w/speed_of_sound;

nHarmonics = 7;
[boundaryIndices, elements, U, F, coupling] = solveForward(speed_of_sound, w, kappa, centers, radii, values, domain, nHarmonics);

%% disc sources
clearvars
% speef of sounde
speed_of_sound = 1540; % m/s

% signal period or cener frequency
T = 10^-4;
w = 2*pi*1/T;

% our domain
bcenter = [1/2,1/2];
brad = 1/2;
domain = [bcenter, brad];


% point scatterers and their domain
values = [300, 300];
radii = [0.03, 0.03];
centers = [3/4, 3/4; 1/4, 3/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound, we set b=0, neglecting the diffusitivity of sound
kappa = w/speed_of_sound;

nHarmonics = 7;
[boundaryIndices, elements, U, F, coupling] = solveForward(speed_of_sound, w, kappa, centers, radii, values, domain, nHarmonics);


%% inspect what happens to the source by the (incomplete forward)coupling

% first harmonic right hand side
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(F(2,:)), 'facecolor', 'interp'); shading interp;
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), abs(F(2,:)), 'facecolor', 'interp'); shading interp;

% second hermonic right hand side
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(F(3,:)), 'facecolor', 'interp'); shading interp;
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), abs(F(3,:)), 'facecolor', 'interp'); shading interp;

% coupling
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(coupling(1,:)), 'facecolor', 'interp'); shading interp;
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(coupling(2,:)), 'facecolor', 'interp'); shading interp;
%%

t_step = 0.005;
tind = 0:t_step:1;
n_t = size(tind,2);

% construct p(t,x)
z = repmat(exp(-1i.*w.*T.*tind)', 1, size(U,2));
pC = zeros(size(z,1), size(z,2));

for k=1:nHarmonics
    pC = pC + repmat(U(k,:),size(z,1),1).*repmat(exp(-1i.*k.*w.*tind.*T)', 1, size(U,2));
end

P = real(pC);

%% let's take a look at p(t,x); playback the periodic solution

figure;
for i=1:n_t
    trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), P(i,:), 'facecolor', 'interp'); shading interp;
    view([0 90]);
    drawnow
    pause(0.1)
end

%% extract harmonics from (complex) solution
% this allows to check whether we have constructued the solution correctly

% DFT
for k=1:nHarmonics
    H(k,:) = 1/(n_t).*sum(repmat(exp(1i.*k.*w.*tind.*T)', 1, size(pC,2)).*pC,1);
    figure, subplot(2,1,1);
    trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(k,:)), 'facecolor', 'interp'); shading interp;
    subplot(2,1,2);
    trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(H(k,:)), 'facecolor', 'interp'); shading interp;
end



%% extract harmonics from (real) solution

% DFT
for k=1:nHarmonics
    Hr(k,:) = 1/(n_t).*sum(repmat(exp(1i.*k.*w.*tind.*T)', 1, size(P,2)).*P,1);
end

% there is a noticeable difference between real(H) and Hr, obviously!

%% inspect specific harmonics
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(1,:)), 'facecolor', 'interp'); shading interp;
title("Real part of p_1(x).")
xlabel('x');
ylabel('y');
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(2,:)), 'facecolor', 'interp'); shading interp;
title("Real part of p_2(x).")
xlabel('x');
ylabel('y');
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(3,:)), 'facecolor', 'interp'); shading interp;
title("Real part of p_3(x).")
xlabel('x');
ylabel('y');

%% try to reconstruct with known centers, only intensities

% find points nearest to center on our "grid"
for i=1:size(centers,2)

    [v,pcentersIdx(i)] = min(sum((elements.points - centers(:,i)').^2,2)); 

end

pcenters = elements.points(pcentersIdx,:)';

% normals of boundary nodes
normals = elements.points(elements.bedges(:,1),:) - bcenter;
normals = normals./sqrt(sum(normals.^2,2));

bFunc = @(x,d,K) exp(1i.*K.*(x*d'));

% random directions for our basis wavelets
nD = 500;

nKappa = nHarmonics-1;
vKappa = (1:nKappa).*kappa; % harmonic frequencies

dirs = rand(nD,2)*2 - 1;
dirs = dirs./sqrt(sum(dirs.^2,2));
elements.bedgesLength = zeros(size(elements.points,1),1);

% compute the length of the edges for numeric integration
for l = 1:size(elements.tri,1)
    c = intersect(elements.tri(l,:), elements.bedges(:,1), 'stable');
    if size(c,1)>1
        elements.bedgesLength(c(1),1) = sqrt(sum((elements.points(c(2), :) - elements.points(c(1),:)).^2));
    end
end

ip = 1 + dirs*normals';

for i = 1:nKappa
    for j=1:size(dirs,1)
        for k = 1:size(pcenters,2)
            A(i, j, k) = exp(1i.*vKappa(i).*pcenters(:,k)'*dirs(j,:)');
        end
        % int on \partial \Omega = \sum
        % weight_i*u(x_i)*<\nabla(w_xi),normal_i> dx
        B(i,j) = 1i*vKappa(i)*sum(elements.bedgesLength(elements.bedges(:,1),1).*U(i,(elements.bedges(:,1)))'.*exp(1i.*vKappa(i).*elements.points(elements.bedges(:,1), :)*dirs(j,:)').*ip(j,:)');
    end
end

Amn = reshape(permute(A,[2 1 3]), size(dirs,1)*nKappa, size(pcenters,2));
Bmn = reshape(permute(B,[2 1]), size(dirs,1)*nKappa, 1);

intensities = inv(Amn'*Amn)*Amn'*Bmn;
abs(intensities)

%% 
m = 2;
f = zeros(size(elements.points,1),1);

g = f;
h = f;

f(pcentersIdx) = abs(intensities')./(coupling(m,pcentersIdx));

tt = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, 1/speed_of_sound, f.*m^2.*kappa^2.*coupling(m,:), h, g, size(elements.points,1));
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(tt), 'facecolor', 'interp'); shading interp;
err = sqrt(sum(abs(tt - U(2,:)').^2,1));
%[boundaryIndices, elements, Uapprox, Fapprox, couplingapprox] = solveForward(speed_of_sound, w, kappa, centers, radii,  abs(intensities'./coupling(1,pcentersIdx)), domain, nHarmonics);
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(m,:)), 'facecolor', 'interp'); shading interp;


%% playground :)

%% plain gradient descent for identifying sources and there values
m = 2;
sourcePoints = [];
sourceValues = [];

f = zeros(size(elements.points,1),1);
y = U(2,elements.bedges(:,1));
g = f;
h = f;
f(elements.bedges(:,1)) = -1i.*m.*kappa.*y;
% K*(-y) 
%te = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));

%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(te), 'facecolor', 'interp'); shading interp;
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(2,:)), 'facecolor', 'interp'); shading interp;


for j = 1:2
    
    te = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));

    [v, pc] = max(abs(te));
    center = elements.points(pc,:)';
    sourcePoints = [sourcePoints center];
    % try to maximize the centers source value ||Kx - y||^2, where x are our
    % source values
    v = 1;
    F = constructF(elements, sourcePoints, zeros(size(sourcePoints,2),1), [sourceValues v]);
    v = [sourceValues v];
    fun = @(x) solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, 1/speed_of_sound, constructF(elements, sourcePoints, zeros(size(sourcePoints,2),1), [x]).*m^2.*kappa^2.*coupling(m,:), g, g, size(elements.points,1));
    sigma = 10^-4;
    for k = 1:200
        lRate = 1000./min(abs(coupling(m,:)));
        % gradient
        h = 0.001;
        x = v;
        for i=1:size(x,2)
            e_i = zeros(1,size(x,2));
            e_i(i) = 1;
            y0 = fun(x-h*e_i);
            z = (y0(elements.bedges(:,1)).' - y);
            val0 = z*z';
    
    %         y0 = fun(x);
    %         z = (y0(elements.bedges(:,1)).' - y);
    %         val0 = z*z';
            y1 = fun(x+h*e_i);
            z = (y1(elements.bedges(:,1)).' - y);
            val1 = z*z';
            df(i,k) = 1/(2*h) *(val1 - val0);  
        end
        
        direction = -df(:,k);
        y_current = fun(x);
        loss = (y_current(elements.bedges(:,1)).' - y)*(y_current(elements.bedges(:,1)).' - y)';
        
        % line search / amijo conditions
        while 1
            y_new = fun(x + lRate*direction');
            loss_new = (y_new(elements.bedges(:,1)).' - y)*(y_new(elements.bedges(:,1)).' - y)';
            if (loss_new < loss + sigma*lRate*(direction'*direction))
                break;
            end
            lRate = lRate/10;
        end
    
        v = x - lRate*df(:,k)';
    end
    
    sourceValues = v;
    f(elements.bedges(:,1)) = 1i.*m.*kappa.*(y_new(elements.bedges(:,1)).' - y);


end
% uu = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, 1/speed_of_sound, F, g, g, size(elements.points,1));
% f(elements.bedges(:,1)) = -1i.*m.*kappa.*(uu(elements.bedges(:,1)).' - y);
% te2 = solveHelmholtzVectorizedTmp(elements, m*w, m*kappa, -1/speed_of_sound, zeros(size(elements.points,1),1), f, g, size(elements.points,1));

%% try out the predicted sources
[boundaryIndices, elements, Uapprox, Fapprox, couplingapprox] = solveForward(speed_of_sound, w, kappa, sourcePoints, [0,0], sourceValues, domain, nHarmonics);

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(Uapprox(2,:)), 'facecolor', 'interp'); shading interp;
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(Uapprox(3,:)), 'facecolor', 'interp'); shading interp;