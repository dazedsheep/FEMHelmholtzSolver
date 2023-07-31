
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
gamma = 1;
nHarmonics = 7;
[boundaryIndices, elements, U, F, coupling] = solveForward(speed_of_sound, w, gamma, kappa, centers, radii, values, domain, nHarmonics);

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
values = [200, 400];
radii = [0.08, 0.03];
centers = [3/4, 3/4; 1/4, 3/4];

% actually kappa = w/(sqrt(c^2 + i*w*b), b accounts for the diffusitivity
% of sound, we set b=0, neglecting the diffusitivity of sound
kappa = w/speed_of_sound;

nHarmonics = 7;
gamma = 1;

[boundaryIndices, elements, U, F, coupling] = solveForward(speed_of_sound, w, gamma, kappa, centers, radii, values, domain, nHarmonics);


%% inspect what happens to the source by the (incomplete forward)coupling

% first harmonic right hand side
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(U(2,:)), 'facecolor', 'interp'); shading interp;
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


%% try out the predicted sources
[boundaryIndices, elements, Uapprox, Fapprox, couplingapprox] = solveForward(speed_of_sound, w, kappa, sourcePoints, [0,0], sourceValues, domain, nHarmonics);

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(Uapprox(2,:)), 'facecolor', 'interp'); shading interp;
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(Uapprox(3,:)), 'facecolor', 'interp'); shading interp;

%% this validates the linearity of the solution operator for m=2 :)

% F(\eta^\sim) - F(\eta) = F(\eta^\sim - \eta)
eta_sim = constructF(elements, centers(1), radii(1), values(1)*0.4);
eta = constructF(elements, centers(1), radii(1), values(1)*0.8);

m = 2;
f = zeros(size(elements.points,1),1);
f(elements.bedges(:,1)) = -(1i.*m.*w.*1/speed_of_sound + gamma).*U(m,elements.bedges(:,1));
g = zeros(size(elements.points,1),1);
f1 = solveHelmholtzVectorizedTmp(elements, 0, 0, m*kappa, 0, -(eta + eta_sim).*1/4.*m^2.*kappa^2.*coupling(m-1,:).', f, g, size(elements.points,1));
f2 = solveHelmholtzVectorizedTmp(elements, 0, 0, m*kappa, 0, -eta.*1/4.*m^2.*kappa^2.*coupling(m-1,:).', f, g, size(elements.points,1));
v = f1 - f2;

fdiff = solveHelmholtzVectorizedTmp(elements, 0, 0, m*kappa, 0, -(eta_sim).*1/4.*m^2.*kappa^2.*coupling(m-1,:).', zeros(size(elements.points,1),1), g, size(elements.points,1));
%figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(te), 'facecolor', 'interp'); shading interp;

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(fdiff), 'facecolor', 'interp'); shading interp;

sum(fdiff(elements.bedges(:,1)) .* U(2,elements.bedges(:,1))')



%%
f = zeros(size(elements.points,1),1);
y = -U(2,:).';
f(elements.bedges(:,1)) = -(1i.*m.*w.*1/speed_of_sound + gamma).*y(elements.bedges(:,1));
u = solveHelmholtzVectorizedTmp(elements, 0, 0, m*kappa, 0, zeros(size(elements.points,1),1), f, g, size(elements.points,1));
adj = 1/4.*kappa.^2.*m.^2.*coupling(m-1,:).'.*conj(u);

sum((eta_sim - eta).* adj)
%% adjoint of the solution operator, this works only for m = 2, because after m = 2 the solution operator is not linear anymore
m = 2;
y = U(2,:);
f = zeros(size(elements.points,1),1);
h = f;
g = f;

adjPDEsol = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -conj(y), h, g, size(elements.points,1));
adj = coupling(m-1,:)'.*conj(adjPDEsol).*1/4.*m.^2.*kappa.^2;

eta = constructF(elements, centers, radii, values);
testp = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -eta.*coupling(m-1,:).'.*1/4.*m.^2.*kappa.^2, h, g, size(elements.points,1));
interior = testp.*U(m,:).';
%error = sum(real(testp.*U(m,:)') - real((adjPDEsol.*coupling(m-1,:).'.*1/4.*m.^2.*kappa.^2).*eta),1)

%% formulate Landweber iteration for m=2, the linear inverse problem
% we reconstruct in L^2(\Omega;C) so we will use the forward oprator
% F(\eta.*kappa^2.*m^2.*coupling)=u solving the helmholtz equation for m=2
% this gives us an adjoint operator in the Hilbert space setting and we
% obtian by F^*(y) = a where a is the solution to the adjoint pde
m = 2;
f = zeros(size(elements.points,1),1);
g = f;
h = f;
f_x_k = f;
N = 2000;
eta = zeros(size(elements.points,1), N);
error = zeros(N,1);
%eta(:,1) = constructF(elements, centers, radii, [1 1]);
f_x_k = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -eta.*coupling(m-1,:).'.*1/4.*m.^2.*kappa.^2, h, g, size(elements.points,1));

stepsize = 1/mean(abs(coupling(m-1,:)));   % (initial) step size should be dependent on the damping of the medium for ultrasonic waves...

for i = 1:N
    h = zeros(size(elements.points,1),1);
    error(i) = real(f_x_k - U(m,:).')'*real(f_x_k - U(m,:).');
    f = (f_x_k - U(m,:).');

    adjointStatePDE = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, conj(f), h, g, size(elements.points,1));
    adjointState = conj(adjointStatePDE);
    
    eta(:, i+1) = eta(:,i) + stepsize*(adjointState);
    %f = 1/4.*eta(:,i+1).*m^2.*kappa^2.*coupling(m-1,:).';
    f = eta(:,i+1);
    f_x_k =solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -f, h, g, size(elements.points,1));
    %figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(f_x_k), 'facecolor', 'interp'); shading interp;
    
    while( real(f_x_k - U(m,:).')'*real(f_x_k - U(m,:).') > error(i))
        stepsize = stepsize/2;
        eta(:, i+1) = eta(:,i) + stepsize*(adjointState);
        f = eta(:,i+1);
        f_x_k =solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -f, h, g, size(elements.points,1));
    end
   
end

%     while( real(f_x_k - U(2,:).')'*real(f_x_k - U(2,:).') > error(i))
%         stepsize = stepsize/2;
%         eta(:, i+1) = eta(:,i) + stepsize*real(adjointState);
%         f = 1/4.*eta(:,i+1).*m^2.*kappa^2.*coupling(m-1,:).';
%         f_x_k =solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -f, h, g, size(elements.points,1));
%     end

my_eta = eta(:,i+1);
for j=1:size(my_eta,1)
    if abs(coupling(m-1,j)) > 0
        reconEta(j) = my_eta(j)/(coupling(m-1,j)*m^2*kappa^2*1/4);
    end
end
reconstructedEta = abs(reconEta);

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), reconstructedEta, 'facecolor', 'interp'); shading interp;
%% inner to outer value given
m = 2;
y = U(m,:);
h = zeros(size(elements.points,1),1);
g = h;
eta = constructF(elements, centers, radii, values);
innerPDE = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -conj(y), h, g, size(elements.points,1));

uu = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -eta.*coupling(m-1,:).'.*1/4.*m^2.*kappa^2, h, g, size(elements.points,1));
testh = zeros(size(elements.points,1),1);
testh(elements.bedges(:,1)) = -conj(y(elements.bedges(:,1)));
outerPDE = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, zeros(size(elements.points,1),1), testh, g, size(elements.points,1));

sqrt(sum(abs(innerPDE - outerPDE).^2))


y(elements.bedges(:,1))*y(elements.bedges(:,1))'
sum(-eta.*coupling(m-1,:).'.*1/4.*m^2.*kappa^2.*(outerPDE))

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(outerPDE), 'facecolor', 'interp'); shading interp;

%% data only given on the boundary, still the linear case of m=2 - TODO
m = 2;
f = zeros(size(elements.points,1),1);
g = f;
h = f;
f_x_k = f;
N = 145;
eta = zeros(size(elements.points,1), N);
boundaryElements = elements.bedges(:,1);
error = zeros(N,1);
%eta(:,1) = constructF(elements, centers, radii, [1 1]);
f_x_k = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -eta.*coupling(m-1,:).'.*1/4.*m.^2.*kappa.^2, h, g, size(elements.points,1));
f_x_k = f_x_k(boundaryElements);
y = U(m,boundaryElements).';
stepsize = 1/mean(abs(coupling(m-1,:)));   % (initial) step size should be dependent on the damping of the medium for ultrasonic waves...
stepsize = 1;
for i = 1:N
    h = zeros(size(elements.points,1),1);
    error(i) = real(f_x_k(boundaryElements) - y)'*real(f_x_k(boundaryElements) - y);
    f = (y - f_x_k(boundaryElements));
    h(boundaryElements) = -conj(f);
    adjointStatePDE = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, zeros(size(elements.points,1),1), h, g, size(elements.points,1));
    adjointState = conj(adjointStatePDE);
    % normalize
    adjnormalized = adjointState./(norm(adjointState,2));
    
    eta(:, i+1) = eta(:,i) - stepsize*(adjnormalized);
    f = eta(:,i+1);
    %eta(:, i+1) = eta(:,i) - stepsize*(conj(1/4 * kappa^2*m^2*coupling(m-1,:).' .* adjointStatePDE));
    %f = 1/4 * kappa^2*m^2*coupling(m-1,:).'.*eta(:,i+1);
    h = zeros(size(elements.points,1),1);

    f_x_k = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -f, h, g, size(elements.points,1));
    %figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(f_x_k), 'facecolor', 'interp'); shading interp;
    
    while( real(f_x_k(boundaryElements) - y)'*real(f_x_k(boundaryElements) - y) > error(i))
        stepsize = stepsize/2;
        eta(:, i+1) = eta(:,i) - stepsize*(adjnormalized);
        f = eta(:,i+1);
        f_x_k = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -f, h, g, size(elements.points,1));
    end
   
end

%     while( real(f_x_k - U(2,:).')'*real(f_x_k - U(2,:).') > error(i))
%         stepsize = stepsize/2;
%         eta(:, i+1) = eta(:,i) + stepsize*real(adjointState);
%         f = 1/4.*eta(:,i+1).*m^2.*kappa^2.*coupling(m-1,:).';
%         f_x_k =solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -f, h, g, size(elements.points,1));
%     end

my_eta = eta(:,i+1);
for j=1:size(my_eta,1)
    if abs(coupling(m-1,j)) > 0
        reconEta(j) = my_eta(j)/(coupling(m-1,j)*m^2*kappa^2*1/4);
    end
end
reconstructedEta = abs(reconEta);

figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), reconstructedEta, 'facecolor', 'interp'); shading interp;