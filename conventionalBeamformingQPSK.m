%% Wireless Communication Project 2022/2023 - Politecnico di Milano
%% Giacomo Sguotti 10667547
% The following script performs Conventional beamforming using a linear 
% antenna array in uplink with 3D geometry, a variable number of fixed 
% vehicles, a single LoS ray (Free Space Loss) and a QPSK modulated
% narrowband signal modulating a variable.
% All the parameters and the geoemtry of the scenario can be easily
% modified in the first section of the code (Scenario Description).
close all
clear all

%% Scenario description
% Define parameters
Pars.fc = 1e9; % Carrier frequency
Pars.c = physconst('Lightspeed'); % speed of light

Pars.lambda = Pars.c/Pars.fc; % wavelength (derived from above parameters)

Pars.BSsize = 6; % Size of the linear antenna array
Pars.BSspacing = 1/2; % Spacing between antennas in terms of wavelengths

Pars.SNR = 1000; % Signal to Noise ratio of all the signals expressed in [dB]

Pars.nbBits = 4; % Number of bits transmitted by each vehicle (must be multiple of 2)

% Define geometry of the problem (xyz coordinates)
Geometry.BSPos=[0,0,100];    % Position of macrocell BS
Geometry.V1Pos=[25,-25,0];   % Position for Vehicle 1
Geometry.V2Pos=[-70,80,0];   % Position for Vehicle 2
Geometry.V3Pos=[0,50,-10];    % Position for Vehicle 3

% Define the antenna array
Geometry.BSArray = phased.ULA('NumElements',Pars.BSsize,'ElementSpacing',Pars.lambda*Pars.BSspacing, 'ArrayAxis','x');

Geometry.BSAntennaPos = getElementPosition(Geometry.BSArray);

% Visualize the scenario
nbVehicles = VisualizeScenario(Geometry);

%% Generate a QPSK signal for the vehicles

qpskmod = comm.QPSKModulator('BitInput', true, 'SymbolMapping', 'Gray');
nbAntennas = Geometry.BSArray.NumElements; % nb. of antennas at BS 

% Generate Pars.nbBits bits to be trasmitted for each vehicle
input_seqs = randi([0 1], Pars.nbBits, nbVehicles);

% Variable preallocation
sent_signals = zeros(ceil(Pars.nbBits/2), nbVehicles);
los_attenuations = zeros(nbVehicles);
attenuated_signals = zeros(ceil(Pars.nbBits/2), nbVehicles);
attenuated_signals_with_awgn = zeros(ceil(Pars.nbBits/2), nbVehicles);
steering_vectors = zeros(nbAntennas, nbVehicles);
weights = zeros(nbAntennas, nbVehicles);

for vehicleIndex = 1:1:nbVehicles
    % Generate the QPSK signal
    
    sent_signals(:,vehicleIndex) = qpskmod(input_seqs(:, vehicleIndex));
    
    % LoS attenuation
    startPositionField = strcat('V', num2str(vehicleIndex), 'Pos');
    distanceToBS = DistanceBetweenTwoPoints(Geometry.BSPos, Geometry.(startPositionField));
    
    los_attenuations(vehicleIndex) = sqrt(2)*Pars.lambda/(4*pi*distanceToBS);
    attenuated_signals(:,vehicleIndex) = los_attenuations(vehicleIndex).*sent_signals(:,vehicleIndex);
    
    % Add AWGN noise according to the provided SNR
    attenuated_signals_with_awgn(:,vehicleIndex) = awgn(attenuated_signals(:,vehicleIndex), Pars.SNR, 'measured');
    
    % Compute the steering vector
    steering_vectors(:,vehicleIndex) = getSteeringVector(Pars, Geometry, vehicleIndex);
    
    % Compute conventional beamforming weights (w = 1/N * s)
    weights(:, vehicleIndex) = (1/nbAntennas).*steering_vectors(:,vehicleIndex);
    
end


%% Beamformer output

% the signal received by each antenna is the sum of signals received by all sources
signals_received = zeros(height(sent_signals), nbAntennas);
for vehicleIndex = 1:1:nbVehicles
    signals_received(:,:) = signals_received(:,:) + attenuated_signals_with_awgn(:,vehicleIndex) * steering_vectors(:,vehicleIndex).';
end

ys_bf = zeros(nbVehicles, height(sent_signals));
for vehicleIndex = 1:1:nbVehicles
    ys_bf(vehicleIndex, :) = weights(:,vehicleIndex)'*signals_received(:,:).';
end

%% Results
qpskdemod = comm.QPSKDemodulator('BitOutput', true, 'SymbolMapping', 'Gray');
for vehicleIndex = 1:1:nbVehicles
    
    figure
    sgtitle(strcat("Vehicle ", int2str(vehicleIndex)));

    % tx signal on subplot 1 (useful if Pars.nbBits is low for debug
    % reasons)
    subplot(3,1,1);
    scatter(real(sent_signals(:,vehicleIndex)), imag(sent_signals(:,vehicleIndex)));
    hold on;
    scatter([1/sqrt(2), 1/sqrt(2), -1/sqrt(2), -1/sqrt(2)], [1/sqrt(2), -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)], "x");
    title('Tx signal'); 
    xlim([-1 1]);
    ylim([-1 1]);
    xline(0);
    yline(0);
    
    
    % rcv by antenna 1 on subplot 2
    y1 = signals_received(:,1)/steering_vectors(1,vehicleIndex);
    out_seq_ant1 = qpskdemod(y1);
    [nb_err_ant1, ber_ant1] = biterr(out_seq_ant1, input_seqs(:,vehicleIndex));
    
    subplot(3,1,2);
    scatter(real(y1), imag(y1) ,'b');
    hold on;
    scatter(los_attenuations(vehicleIndex).*[1/sqrt(2), 1/sqrt(2), -1/sqrt(2), -1/sqrt(2)], los_attenuations(vehicleIndex).*[1/sqrt(2), -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)], "x");
    title('Received by single antenna(1)');
    xlim(los_attenuations(vehicleIndex).*[-5 5]);
    ylim(los_attenuations(vehicleIndex).*[-5 5]);
    xline(0);
    yline(0);
    xlabel(['NumErr: ' num2str(nb_err_ant1) ' BER ' num2str(ber_ant1)])
    
    
    % Beaforming result on subplot 3
    out_seq_bf = qpskdemod(ys_bf(vehicleIndex, :).');
    [nb_err_bf, ber_bf] = biterr(out_seq_bf, input_seqs(:,vehicleIndex));
    
    subplot(3,1,3);
    scatter(real(ys_bf(vehicleIndex,:)), imag(ys_bf(vehicleIndex,:)));
    hold on;
    scatter(los_attenuations(vehicleIndex).*[1/sqrt(2), 1/sqrt(2), -1/sqrt(2), -1/sqrt(2)], los_attenuations(vehicleIndex).*[1/sqrt(2), -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)], "x");
    title('Beamforming result');
    xlim(los_attenuations(vehicleIndex).*[-5 5]);
    ylim(los_attenuations(vehicleIndex).*[-5 5]);
    xline(0);
    yline(0);
    xlabel(['NumErr: ' num2str(nb_err_bf) ' BER ' num2str(ber_bf)])
    
    % Print a table containing the value of Array Pattern Function in
    % correspondance of each vehicle for each weight
    fprintf('\n');
    disp(strcat('using weight vector of V', int2str(vehicleIndex)));
    disp('Vehicle | Array Pattern');
    disp('--------|--------------')
    for interferer=1:1:nbVehicles
            disp(strcat(int2str(interferer), '       | ', num2str(norm(weights(:,vehicleIndex)'*steering_vectors(:,interferer)))));
    end
    

    
end

%% Functions 
% DistanceBetweenTwoPoints funct
% Given two points in space, it returns their distance
function distance = DistanceBetweenTwoPoints(point1, point2)
    distance = norm(point1-point2);
end

% steering_vector funct
% Compute and returns the steering vector, given Pars and Geometry structs
% and the index of the vehicle.
% The phases are calculated starting from the transmitter.
function steering_vector = getSteeringVector(Pars, Geometry, vehicleIndex)
    startPositionField = strcat('V', num2str(vehicleIndex), 'Pos');
    nbAntennas = Geometry.BSArray.NumElements;
    steering_vector = zeros(nbAntennas,1);
    for antenna_index = 1:1:nbAntennas
        prop_delay = DistanceBetweenTwoPoints(Geometry.(startPositionField), Geometry.BSPos'+Geometry.BSAntennaPos(:,antenna_index)) / Pars.c;
        steering_vector(antenna_index,1) = exp(-1i*2*pi*Pars.fc*prop_delay);
    end
end

% VisualizeScenario funct
% Plot the Geometry scenario and, in a separate subplot the antenna array.
% It returns the number of vehicles from the positions, used later in code.
function nbVehicles = VisualizeScenario(Geometry)
    
    figure;
    subplot(2,1,1);

    % Base Station
    plot3(Geometry.BSPos(1),Geometry.BSPos(2), Geometry.BSPos(3), '^b', 'DisplayName','BS')
    text(Geometry.BSPos(1),Geometry.BSPos(2), Geometry.BSPos(3), 'BS', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    
    grid on
    hold on

    % Vehicles
    vehicleIndex = 1;
    while (isfield(Geometry, strcat('V', num2str(vehicleIndex), 'Pos')))
        % Start Position
        startPositionField = strcat('V', num2str(vehicleIndex), 'Pos');
        plot3(Geometry.(startPositionField)(1),Geometry.(startPositionField)(2), Geometry.(startPositionField)(3), 'v', 'DisplayName',strcat('Vehicle', num2str(vehicleIndex)))
        text(Geometry.(startPositionField)(1),Geometry.(startPositionField)(2), Geometry.(startPositionField)(3), strcat('V', num2str(vehicleIndex), 'start'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')

        % Update vehicle index
        vehicleIndex = vehicleIndex + 1;
    end
    
    nbVehicles = vehicleIndex - 1;
    disp(strcat('nb. vehicles: ', num2str(nbVehicles)))
    
    title('Scenario')
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    legend
    
    subplot(2,1,2);
    grid on
    hold on
    
    % Base Station
    for antenna_index = 1 : size(Geometry.BSAntennaPos,2)
        x_coord = Geometry.BSAntennaPos(1,antenna_index) + Geometry.BSPos(1);
        y_coord = Geometry.BSAntennaPos(2,antenna_index) + Geometry.BSPos(2);
        z_coord = Geometry.BSAntennaPos(3,antenna_index) + Geometry.BSPos(3);
        name = strcat('A', int2str(antenna_index));
        plot3(x_coord, y_coord, z_coord, '^b', 'DisplayName',name)
        text(x_coord, y_coord, z_coord, name, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    end

    title("Antenna Array in detail")
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    legend

end