
%% Solver for the periodic Westervelt equation
clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1480; % m/s

% signal period or center frequency
T = 10^-4;
omega = 2*pi*1/T;

% our domain
bcenter = [0,0];
brad = 1;
domain = [bcenter, brad];
% non linearity parameter of our domain (water = 5)
sourceValueDomain = 0;

% point scatterers and their domain
values = [9, 0];
refractionIndex = [1, 1];
% linear case
%values = [0, 0];
radii = [0.1, 0.1];
centers = [0, 0; 0.5, 0.2];

diffusivity = 10^(-9);

minHarmonics = 5; % minimum number of harmonics
nHarmonics = 5; % maximum number of harmonics

% impdeance boundary conditions --> massDensity cancels
% higher frequencies are taken into account later
beta = 1/(speed_of_sound);
gamma = 10^(-9);

meshSize = 0.005;

excitationPoints = [0;0.8];
% typical ultrasound pressure is 1MPa at a frequency of 1 MHz, lower
% frequency -> lower pressure!
pressure = 10*10^6;
excitationPointsSize = [0.01];
%excitationPower(1,1) = referencePressure;
%excitationPower(1,2:nHarmonics) = 0;
nExcitationHarmonics = 1;

[elements] = initializeMultiLeveLSolver(meshSize, domain);

% plot the positions of excitation(s) and source(s)
objects = getGridPointsLE(elements, [centers excitationPoints], [radii excitationPointsSize]);
figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), objects, 'facecolor', 'interp'); shading interp;
xlabel("x [m]");
ylabel("y [m]");

% construct nonlinearity
f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
% construct all space dependent wave numbers for all harmonics
kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);

% build a gaussian source
%source = 1./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*referencePressure.*gaussianSource(elements, excitationPoints, 0.6);
% build a point source (regularized dirac)

% the source needs to be scaled by omega^2, this matches the SI units
source = omega^2./(speed_of_sound.^2 + 1i .* omega .* diffusivity).*createPointSource(elements, excitationPoints, meshSize/2).*pressure;  
excitation = zeros(size(elements.points,1),nHarmonics);
excitation(:,1) = source;

%excitation = 1i./(speed_of_sound.^2 + 1i .* (1:nHarmonics) .* omega .* diffusivity).*excitation;
[cN, U, F] = solveWesterveltMultiLevel(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 10^(-12));
H = U;
U = squeeze(U(cN,:,:));

% difference between harmonics figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), abs(real(H(6,3,:))-real(H(14,3,:))), 'facecolor', 'interp'); shading interp;

%%
% calc the solution(s) for the multi level harmonic approximation
t_step = 1/(4/T);
%t_step = 0.0005;
% simulate 1000 steps
%tind = 0:t_step:10000*t_step;
tind = 0:0.01:0.02;
n_t = size(tind,2);

% compute only for n iterations
iter = 1;

% construct p(t,x)
z = repmat(exp(-1i.*omega.*T.*tind)', 1, size(U,2));
P = zeros(iter, size(z,1), size(z,2));
o = 0;
for j=(min(size(H,1),cN) - iter):min(size(H,1),cN)
    o = o + 1;
    rpC = zeros(size(z,1), size(z,2));
    for k=1:j
        rpC = rpC + squeeze(repmat(H(j,k,:),size(z,1),1)).*repmat(exp(1i.*k.*omega.*tind.*T).', 1, size(U,2));
    end
    p_L_2(o) = sqrt(sum(sum(real(rpC).^2, 2).^2,1));
    p_L_2_t(o) = sqrt(sum(real(rpC(1,:)).^2, 2));
    P(o,:,:) = real(rpC);
   
end

%figure, plot(abs(p_L_2(2:(size(p_L_2,2))) - p_L_2(1:(size(p_L_2,2)-1))))

%% points along a line segment (cut through)
% lsegStart = [-1/2;-1/2];
% lsegEnd = [2/4;2/4];
% due to the mesh, the peak  of the source may be a little off
[v,peakidx] = max(P(1,1,:));

lsegStart = elements.points(peakidx,:)';
lsegEnd = [0;-0.7];

npoints = 4000; 

unitVec = (lsegEnd - lsegStart)/norm((lsegEnd - lsegStart));
points = lsegStart + unitVec.*norm((lsegEnd - lsegStart),2)/npoints .* (0:(npoints-1));

idxList = zeros(npoints,1);

for j=1:size(points,2)
        [v,pcenterIdx] = min(sum((elements.points - points(:,j)').^2,2)); 
        idxList(j) = pcenterIdx;    
end

coords = [elements.points(idxList,1)';elements.points(idxList,2)'];

t = 1;

pPoint = squeeze(P(iter,t,idxList));

% find the points in space which are directly affected by nonlinearity
idxNonlinearity = find((sqrt(dot(points - centers(:,1) ,points - centers(:,1))) < radii(1))==1);

pointdist = sqrt(dot(points - excitationPoints,points - excitationPoints));
% smooth the data for plotting
smdata = smoothdata(pPoint, 'gaussian', 50);
figure, plot(pointdist, smdata)
hold on
plot(pointdist(idxNonlinearity), smdata(idxNonlinearity),'r')
xlabel('x');
ylabel('Pa');


%%
figure, plot(pointdist, smdata)
hold on
figure, plot(pointdist, 10*log10(abs(smdata)/referencePressure))
xlabel('x');
ylabel('Pa');

%% pick a point and look at the frequency domain
% therefore we sample at 2*f_max which is in our case 2*nHarmonics*1/T
point = [0;0.1];

[v,idx] = min(sum((elements.points - point(:)').^2,2)); 

node = [elements.points(idx,1);elements.points(idx,2)];

% sampling frequency in time
Fs = 1/T * 2 * (nHarmonics+1);

N = 2000;
pC = zeros(1,N);
for m=1:nHarmonics
    pC = pC + U(m,idx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
end

% time = (0:(N-1))*1/Fs;
% timescale = 10^3;
% figure, plot(time*timescale,real(pC))
% ylabel("Acoustinc Pressure [Pa]");
% xlabel("Time [ms]");

window = hanning(N);
freq =  Fs/N*(0:(N/2));
y = abs(fft(window'.*real(pC)))/N;
y1 = y(1:N/2+1);
y1(2:end-1) = 2*y1(2:end-1);
shiftedTF = fftshift(fft(real(pC)))/N;
TFdB = 10*log10(y1/pressure);
fscaling = 10^3;
figure, plot(freq./fscaling, TFdB(1:N/2+1))
xlabel("Frequency [kHz]")
ylabel("P/P0 [dB]")

%%
Fs = 1/T * 2 * (nHarmonics+1)*20;

N = 3000;
pC = zeros(1,N);
for m=1:nHarmonics
    pC = pC + U(m,idx) .* exp(1i.*m.*omega.*(0:(N-1))*1/Fs);
end

time = (0:(N-1))*1/Fs;
timescale = 10^3;
figure, plot(time*timescale,real(pC))
ylabel("Acoustinc Pressure [Pa]");
xlabel("Time [ms]");

%% solution with Green's function for the linear case
H = 1i/4 .* besselh(0,omega.*pointdist./speed_of_sound);
figure,plot(real(H))
figure,plot(pointdist,real(H))
%%
% calc the solution(s) for the multi level harmonic scheme starting at 1
t_step = 1/(4/T);
%t_step = 0.0005;
% simulate 1000 steps
%tind = 0:t_step:10000*t_step;
tind = 0:0.01:20;
n_t = size(tind,2);

% compute only for n iterations
maxHarmonics = 1;

% construct p(t,x)
z = repmat(exp(-1i.*w.*T.*tind)', 1, size(U,2));
P = zeros(maxHarmonics, size(z,1), size(z,2));
o = 0;
for j=1:maxHarmonics
    o = o + 1;
    rpC = zeros(size(z,1), size(z,2));
    for k=1:j
        rpC = rpC + squeeze(repmat(H(j,k,:),size(z,1),1)).*repmat(exp(1i.*k.*w.*tind.*T).', 1, size(U,2));
    end
    p_L_2(o) = sqrt(sum(sum(real(rpC).^2, 2).^2,1));
    p_L_2_t(o) = sqrt(sum(real(rpC(1,:)).^2, 2));
    P(o,:,:) = real(rpC);
end

%figure, plot(abs(p_L_2(2:(size(p_L_2,2))) - p_L_2(1:(size(p_L_2,2)-1))))

%%
L2_error = zeros(nHarmonics-1,1);

for i = 1:(nHarmonics-1)
    L2_error(i) = sum(sum(((P(i,:,:)) - P(nHarmonics,:,:)).^2));
end

figure, plot(log10(L2_error));
xlabel('Harmonics/Iterations');
ylabel('L2 Error to last iteration (log10)')



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

%% build pressure from multiharmonic expansion
t_step = 1/(nHarmonics/T);
%t_step = 0.0005;
% simulate 1000 steps
%tind = 0:t_step:10000*t_step;
tind = 0:t_step:1;
n_t = size(tind,2);

% construct p(t,x)
z = repmat(exp(-1i.*w.*T.*tind)', 1, size(U,2));
pC = zeros(size(z,1), size(z,2));

for k=1:nHarmonics
    pC = pC + repmat(U(k,:),size(z,1),1).*repmat(exp(1i.*k.*w.*tind.*T)', 1, size(U,2));
end

P = real(pC);


%% pick a point and look at the wave
point = [1/2;1/2];
[v,pcenterIdx] = min(sum((elements.points - point').^2,2));

coords = [elements.points(pcenterIdx,1);elements.points(pcenterIdx,2)];

pPoint = P(iter,:,pcenterIdx);

%% let's take a look at p(t,x); playback the periodic solution
figure;
for i=1:1:n_t
    trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), P(iter,i,:), 'facecolor', 'interp'); shading interp;
    view([0 90]);
    drawnow
    pause(1)
end
%% for a single point in space
t_step = 1/(2.5*nHarmonics/T);
tind = 0:t_step:10000*t_step;
n_t = size(tind,2);
x = 20000;
p = zeros(1,n_t);
for k=1:nHarmonics
    p = p + U(k,x).*exp(1i.*k.*w.*tind);
end
p = real(p);

% fft in order to show the n harmonics
fs = 1/t_step;
L = size(p,2);
Y = fft(p);
f = fs*(0:(L/2))/L;
figure, plot(f, abs(Y(1:(L/2)+1)/L))

%%
tvec = (0:size(P,1)-1)*1/(nHarmonics*1/T);
figure,  plot(1/(nHarmonics*1/T) * (0:(size(P,1)/2))/size(P,1), abs(Y(1:size(P,1)/2+1)/size(P,1)))

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
f_x_k = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -eta.*coupling(m,:).'.*1/4.*m.^2.*kappa.^2, h, g, size(elements.points,1));

stepsize = 1/mean(abs(coupling(m,:)));   % (initial) step size should be dependent on the damping of the medium for ultrasonic waves...

for i = 1:N
    h = zeros(size(elements.points,1),1);
    error(i) = (f_x_k - U(m,:).')'*(f_x_k - U(m,:).');
    f = (f_x_k - U(m,:).');

    adjointStatePDE = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, conj(f), h, g, size(elements.points,1));
    adjointState = conj(adjointStatePDE);
    
    eta(:, i+1) = eta(:,i) + stepsize*(adjointState);
    %f = 1/4.*eta(:,i+1).*m^2.*kappa^2.*coupling(m-1,:).';
    f = eta(:,i+1);
    f_x_k =solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -f, h, g, size(elements.points,1));
    %figure, trisurf(elements.tri(:,1:3), elements.points(:,1), elements.points(:,2), real(f_x_k), 'facecolor', 'interp'); shading interp;
    
    while( (f_x_k - U(m,:).')'*(f_x_k - U(m,:).') > error(i))
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
        reconEta(j) = abs(my_eta(j))*2;
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
f_x_k = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -eta.*coupling(m,:).'.*1/4.*m.^2.*kappa.^2, h, g, size(elements.points,1));
f_x_k = f_x_k(boundaryElements);
y = U(m,boundaryElements).';
stepsize = 1/mean(abs(coupling(m,:)));   % (initial) step size should be dependent on the damping of the medium for ultrasonic waves...
stepsize = 1;
for i = 1:N
    h = zeros(size(elements.points,1),1);
    error(i) = (f_x_k(boundaryElements) - y)'*(f_x_k(boundaryElements) - y);
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
    
    while( (f_x_k(boundaryElements) - y)'*(f_x_k(boundaryElements) - y) > error(i))
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

%% data only given on the boundary, still the linear case of m=2, partial data
percentageBoundaryData = 0.2;
m = 2;
f = zeros(size(elements.points,1),1);
g = f;
h = f;
f_x_k = f;
N = 145;
eta = zeros(size(elements.points,1), N);
numBoundaryElements = floor(size(elements.bedges,1)*percentageBoundaryData);
boundaryElements = elements.bedges(1:numBoundaryElements,1);
error = zeros(N,1);
%eta(:,1) = constructF(elements, centers, radii, [1 1]);
f_x_k = solveHelmholtzVectorizedTmp(elements, m*w, gamma, m*kappa, 1/speed_of_sound, -eta.*coupling(m,:).'.*1/4.*m.^2.*kappa.^2, h, g, size(elements.points,1));
%f_x_k = f_x_k(boundaryElements);
y = U(m,boundaryElements).';
stepsize = 1/mean(abs(coupling(m,:)));   % (initial) step size should be dependent on the damping of the medium for ultrasonic waves...
for i = 1:N
    h = zeros(size(elements.points,1),1);
    error(i) = (f_x_k(boundaryElements) - y)'*(f_x_k(boundaryElements) - y);
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
    
    while( (f_x_k(boundaryElements) - y)'*(f_x_k(boundaryElements) - y) > error(i))
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