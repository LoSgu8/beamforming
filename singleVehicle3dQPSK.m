% Conventional beamforming in uplink with 2D geometry, one fixed vehicle
% a single LoS ray (Free Space Loss) and a QPSK modulated narrowband signal
close all
clear all

%% Scenario description
% Define parameters
Pars.fc = 1e9; % Carrier frequency
Pars.c = physconst('Lightspeed');
Pars.lambda = Pars.c/Pars.fc;

Pars.BSsize = 4;
Pars.BSspacing = 1/2;

Pars.SNR = 0; %[dB]

% Define geometry of the problem (xyz coordinates)
Geometry.BSPos=[0,0,25];            % Position of macrocell BS
Geometry.V1PosStart=[70,-100,0];   % Start position for Vehicle 1

% Define the antenna array
Geometry.BSArray = phased.ULA('NumElements',Pars.BSsize,'ElementSpacing',Pars.lambda*Pars.BSspacing, 'ArrayAxis','x');

Geometry.BSAntennaPos = getElementPosition(Geometry.BSArray);

% Visualize the scenario
VisualizeScenario(Geometry);

%% Generate a QPSK signal for the vehicle

% Generate the QPSK signal for the vehicle
qpskmod = comm.QPSKModulator('BitInput', true, 'SymbolMapping', 'Gray');

nb_bits = 6;
input_seq = randi([0 1], nb_bits, 1);
sent_signal = qpskmod(input_seq);

% LoS attenuation
Geometry.DistV1Start = DistanceBetweenTwoPoints(Geometry.BSPos,Geometry.V1PosStart);

los_attenuation = sqrt(2)*Pars.lambda/(4*pi*Geometry.DistV1Start);

attenuated_signal = los_attenuation.*sent_signal;

% Add AWGN noise according to the provided SNR
attenuated_signal_with_awgn = awgn(attenuated_signal, Pars.SNR, 'measured');

% Initialize and compute the Steering vector
N = Geometry.BSArray.NumElements; % nb. of antennas at BS 
steering_vector = getSteeringVector(Pars, Geometry);

% Compute conventional beamforming weights (w = 1/N * s)
w = (1/N).*steering_vector;

% Compute the array pattern
arrayPattern = w'*steering_vector;

% Beamformer output
signals_received = attenuated_signal_with_awgn * steering_vector.';

y_bf = w'*signals_received.';

%% Results
qpskdemod = comm.QPSKDemodulator('BitOutput', true, 'SymbolMapping', 'Gray');

figure
% tx signal
subplot(3,1,1);
scatter(real(sent_signal), imag(sent_signal));
hold on;
scatter([1/sqrt(2), 1/sqrt(2), -1/sqrt(2), -1/sqrt(2)], [1/sqrt(2), -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)], "x");
title('Tx signal');
xlim([-1 1]);
ylim([-1 1]);
xline(0);
yline(0);


% rcv by antenna 1
y1 = signals_received(:,1)/steering_vector(1);
out_seq_ant1 = qpskdemod(y1);
[nb_err_ant1, ber_ant1] = biterr(out_seq_ant1, input_seq);

subplot(3,1,2);
scatter(real(y1), imag(y1) ,'b');
hold on;
scatter(los_attenuation.*[1/sqrt(2), 1/sqrt(2), -1/sqrt(2), -1/sqrt(2)], los_attenuation.*[1/sqrt(2), -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)], "x");
title('Received by single antenna(1)');
xlim(los_attenuation.*[-5 5]);
ylim(los_attenuation.*[-5 5]);
xline(0);
yline(0);
xlabel(['NumErr: ' num2str(nb_err_ant1) ' BER ' num2str(ber_ant1)])


% Beaforming result
out_seq_bf = qpskdemod(y_bf.');
[nb_err_bf, ber_bf] = biterr(out_seq_bf, input_seq);

subplot(3,1,3);
scatter(real(y_bf), imag(y_bf));
hold on;
scatter(los_attenuation.*[1/sqrt(2), 1/sqrt(2), -1/sqrt(2), -1/sqrt(2)], los_attenuation.*[1/sqrt(2), -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)], "x");
title('Beamforming result');
xlim(los_attenuation.*[-5 5]);
ylim(los_attenuation.*[-5 5]);
xline(0);
yline(0);
xlabel(['NumErr: ' num2str(nb_err_bf) ' BER ' num2str(ber_bf)])

%% Functions 
% DistanceBetweenTwoPoints funct
function distance = DistanceBetweenTwoPoints(point1, point2)
    distance = norm(point1-point2);
end

% steering_vector funct: Compute and returns the steering vector, given the
% Pars and Geometry struct, the phases are calculated from the transmitter
function steering_vector = getSteeringVector(Pars, Geometry)
    nbAntennas = Geometry.BSArray.NumElements;
    steering_vector = zeros(nbAntennas,1);
    for antenna_index = 1:1:nbAntennas
        prop_delay = DistanceBetweenTwoPoints(Geometry.V1PosStart, Geometry.BSPos'+Geometry.BSAntennaPos(:,antenna_index)) / Pars.c;
        steering_vector(antenna_index,1) = exp(-1i*2*pi*Pars.fc*prop_delay);
    end
end

% VisualizeScenario funct
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
    while ( isfield(Geometry, strcat('V', num2str(vehicleIndex), 'PosStart')) )
        % Start Position
        startPositionField = strcat('V', num2str(vehicleIndex), 'PosStart');
        plot3(Geometry.(startPositionField)(1),Geometry.(startPositionField)(2), Geometry.(startPositionField)(3), '^r', 'DisplayName','Vehicle')
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