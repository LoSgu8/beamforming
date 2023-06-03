# Conventional Beamforming for Uplink Communication in 3D Wireless Networks
This MATLAB script performs conventional beamforming using a linear antenna array in an uplink scenario with 3D geometry. It simulates a wireless communication system with multiple vehicles transmitting QPSK modulated narrowband signals to a base station. The script allows for easy modification of parameters and scenario geometry.
![scenario plot](https://github.com/LoSgu8/beamforming/assets/17884896/25424a17-a309-42dc-b3c8-bbc9df216499)

## Scenario Description

The script starts with a section where you can define the parameters and geometry of the scenario. Here are the parameters you can modify:

- **Pars.fc**: Carrier frequency in Hz.
- **Pars.c**: Speed of light.
- **Pars.lambda**: Wavelength derived from carrier frequency.
- **Pars.BSsize**: Size of the linear antenna array.
- **Pars.BSspacing**: Spacing between antennas in terms of wavelengths.
- **Pars.SNR**: Signal-to-noise ratio of all signals in dB.
- **Pars.nbBits**: Number of bits transmitted by each vehicle (must be a multiple of 2).

The geometry of the problem is defined by the positions of the base station (BS) and the vehicles. You can modify the positions of the following entities:

- **Geometry.BSPos**: Position of the macrocell BS.
- **Geometry.VxPos**: Position for Vehicle x.

The script uses a phased linear antenna array (**Geometry.BSArray**) with the specified parameters. The array geometry is visualized in a separate subplot.
## Signal Generation and Processing

The script generates a QPSK-modulated signal for each vehicle. The number of antennas at the base station is determined from the array size. The generated signals are attenuated based on the distance between the vehicles and the base station, and additive white Gaussian noise is added according to the specified SNR.

The script computes the direction of arrival (DoA) for each vehicle and introduces the channel phase. It then computes the steering vector and conventional beamforming weights. The received signals at each antenna are combined using beamforming weights to obtain the beamformer output.
## Results

The script plots the results for each vehicle in separate figures. Each figure contains three subplots:
 
1. The transmitted QPSK signal.
2. The received signal at a single antenna (after attenuation).
3. The beamforming result.

For each subplot, the script computes the bit error rate (BER) by comparing the received signal with the transmitted signal. The number of errors and the BER are displayed in the subplots.

The script also prints a table showing the Array Pattern Function values for each vehicle using the corresponding weight vector.
## Functions

The script includes several helper functions:

- **DistanceBetweenTwoPoints**: Computes the distance between two points in space.
- **computeDoA**: Computes the direction of arrival (AoA and ZoA) for a vehicle.
- **computeChannelPhase**: Computes the phase introduced by the channel for a vehicle.
- **VisualizeScenario**: Visualizes the scenario by plotting the positions of the base station and vehicles, along with the antenna array.

These functions are used in the main script to calculate distances, angles, and visualize the scenario.

Please note that this script assumes the availability of the phased array toolbox in MATLAB.
