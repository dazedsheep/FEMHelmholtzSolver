%% This is a greedy search algorithm that tries to find avlues for \delta
% We assume a fixed B/A value (here it is set to 5) and fixed source and
% that the nonlinearity spans the whole domain

% do not forget to initialise the parpool!

clearvars

fixedError = 10^-3; % L2 error of two consecutive iterations, we stop if we reach this

massDensity = 1000; %kg/m^3

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

minHarmonics = 3; % minimum number of harmonics
nHarmonics = 30; % maximum number of harmonics

gamma = 1;

excitationPoints = [0;0];
excitationSize = 0.01;

initialMeshSize = 0.0005;

speed_of_sound = 50;

% signal period or center frequency
T = 10^-2;

% create the parpool
% parpool(12);
% for this setup, the solutions start to blow up between 2000 and 3000
pressure = 3000;

beta = 1/speed_of_sound;
meshSize = initialMeshSize;

%
meshSize = 0.0005;
c = 100:100:1800;
freq = 500:100:(1/(10^-4));
pressure = 1000:1000:5*10^5;
[C,F] = meshgrid(c,freq);

% now we check the condition N*abs(\kappa)*meshSize < 1
[elements] = initializeMultiLeveLSolver(meshSize, domain);
numSims = 0;
skipped = 0;
for i=1:size(C,1)
    for j=1:size(F,2)
        % compute kappa (we neglect the diffuse part for screening the
        % permutations)
        omega = 2*pi*F(i,j,1);
        kappa = omega/sqrt((C(i,j,1))^2);
        if (meshSize * abs(kappa) * nHarmonics > 1)
            C(i,j,k) = 0; % 0 speeds will be skipped
            skipped = skipped + 1;
        else
            numSims = numSims + 1;
        end
    end
end

% prepare the source (our mesh is fixed over all iterations)
% normalise the regularised dirac
pointSource = createPointSource(elements, excitationPoints, excitationSize);
[area, s] = integrate_fun_trimesh(elements.opoints, elements.otri, pointSource.');

startPressureIncrease = 100000;
startPressure = 1000;
idM = zeros(size(C));
for i=1:size(C,2)
    for j=1:size(F,1)

        % skip prescreened value ranges that would lead to the
        % pollution effect
        if (C(i,j) < 1)
            break;
        end
        pressure = startPressure;
        while true
            f = constructF(elements, massDensity, C(i,j), refractionIndex, centers, radii, values, sourceValueDomain, true);
            omega = 2*pi*F(i,j,k);
            % construct all space dependent wave numbers for all harmonics
            kappa = constructKappa(elements, diffusivity, C(i,j), omega, refractionIndex, centers, radii, values, nHarmonics);

    

            source = -exp(1i*pi/2).*pressure.*pointSource./s;
            excitation = zeros(size(elements.points,1),nHarmonics);
            excitation(:,1) = source;

            [cN, U, G] =  solveWesterveltMultiLevelMT(elements, omega, beta, gamma, kappa, excitation, f, nHarmonics, minHarmonics, false, fixedError);

            % check if we encountered the case where the solution begins to
            % blow up due to the selected parameters (actually due to the
            % pressure/energy)

            [~, e1] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN-1,:,:)))) - squeeze(real(sum(U(cN-2,:,:))))).^2 )' );
            [~, e2] = integrate_fun_trimesh(elements.points', elements.otri,  squeeze(abs(squeeze(real(sum(U(cN,:,:)))) - squeeze(real(sum(U(cN-1,:,:))))).^2 )' );

            if e2 > e1 || isnan(e2)
                idM(i,j) = pressure;
            else
                lastWorkingPressure = pressure;
            end
        end
    end
end
