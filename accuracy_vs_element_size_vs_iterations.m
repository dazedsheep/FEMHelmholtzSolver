%% runtime evals

clearvars

massDensity = 1000; %kg/m^3
speed_of_sound = 1450; % m/s

% signal period or center frequency 5 kHz

omega = 2*pi*5*10^3;

% our domain
bcenter = [0,0];
brad = 0.05;
domain = [bcenter, brad];
% nonlinearity parameter of our domain (water = 5)
sourceValueDomain = 5;

% point scatterers and their domain
values = [0];
refractionIndex = [1];
% linear case
%values = [0, 0];
radii = [0.05];
centers = [0; 0];

diffusivity = 10^(-9);

gamma = 1;

beta = 1/speed_of_sound;

meshSteps = 6;

initialMeshsize = 0.0025;

meshdivfactor = 2;


excitationPoints = [0;0];

excitationSize = 0.01;

pressure = 12*10^5;
meshSizes = flip( 0.0002:0.00025:0.0025 );
iters = 5:5:30;
% figure, trisurf(source_elements.tri(:,1:3), source_elements.points(:,1), source_elements.points(:,2), source, 'facecolor', 'interp'); shading interp;
%%
for i = 1:length(meshSizes)
    for j=1:length(iters)
        minHarmonics = iters(j); % minimum number of harmonics
        nHarmonics = iters(j); % maximum number of harmonics
        meshSize = meshSizes(i);
        [elements] = initializeMultiLeveLSolver(meshSize, domain);
        f = constructF(elements, massDensity, speed_of_sound, refractionIndex, centers, radii, values, sourceValueDomain, true);
        % construct all space dependent wave numbers for all harmonics

        kappa = constructKappa(elements, diffusivity, speed_of_sound, omega, refractionIndex, centers, radii, values, nHarmonics);
        pointSource = createPointSource(elements, excitationPoints, excitationSize);
        [area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');

        source = -exp(1i*pi/2).*pressure.*pointSource./s;
        excitation = zeros(size(elements.points,1),nHarmonics);
        excitation(:,1) = source;

        tic
        [cN, U, F] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, 0);
        runTime(i,j) = toc;
        numIterations(i,j) = cN;
        H = U;

        % calculate the difference of consecutive solution over the iterations
        for k=2:minHarmonics
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(H(k,:,:)))) - squeeze(real(sum(H(k-1,:,:))))).^2 )' );
            L2errors(i,j,k-1) = sqrt(int);
        end

        for k=1:minHarmonics
            [~, int] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(k,:,:))))).^2)');
            L2Norms(i,j,k) = sqrt(int);
            pollutionVal(i,j,k) = meshSize*abs(kappa(1,1))*k;
        end

        U = squeeze(U(cN,:,:));
        numDOF(i,j) = size(elements.points,1);
        degen(i,j)= max(max(abs(f.* squeeze(abs(squeeze(real(sum(U(cN,:,:)))))))));
    end
end
%% plots for the runtime
styles = {'--o','--+', '--*', '--x', '--square', '--diamond' };
figure,
t = tiledlayout(1,1);
ax1 = axes(t);

for j = 1:length(iters)
    plot(ax1,meshSizes,flip(runTime(:,j)), styles{j}, 'LineWidth',1.3,'MarkerSize',7);
    leg{j} = sprintf('N=%d iterations', iters(j));
    hold on
end

legend(leg, 'Interpreter','latex','Location','northwest');

yscale("log")
xlabel('Element size, $h$', 'Interpreter','latex');
ylabel('Execution time (s)', 'Interpreter','latex');
ax2 = axes(t);
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.YScale = 'log';
xlim (ax2,[numDOF(1,1),numDOF(10,1)]);
ylim (ax2,[0,10^3])
xlabel(ax2, 'DOF', 'Interpreter','latex')
grid on; grid minor;
%% tabular data for N=30
meshSizes(2:2:end)'
numDOF(2:2:end,6)
numElements(2:2:end,6)-1
%L2errors(1:2:end,6,29)
runTime(2:2:end,6)
%L2Norms(1:2:end,6,29)

%%
for i = 1:length(meshSizes)
        meshSize = meshSizes(i);
        [elements] = initializeMultiLeveLSolver(meshSize, domain);
        numNodes(i) = size(elements.points,1);
        numElements(i) = size(elements.tri,1);
end