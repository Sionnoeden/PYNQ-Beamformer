# An Ultrasound Beamformer Runing on PYNQ (_Still under debugging_)

- This project develops a portable ultrasound beamformer runing on PYNQ platform.
- This project aims to provide a back-end data processing system for wearable ultrasound devices.

## What you need for implementing this project

- PYNQ-Z2 board with PYNQ image (version 2.6.2)
- Vitis development tool kit (version 2020.1)
- Matlab (version R2018a)

## How to implement this project

1. Download the latest version of Matlab verification files and run the code to see what can the algorithm do.
2. Include the code in [HLS] folder in Vitis HLS, synthesize it and export the design as Vivado IP.
3. Download the .tcl script in [Vivado] folder and use the command below in Vivado platform to create Vivado project. (You may need to manually add a custom IP directory)
```sh
source ./Beamforming.tcl
```
4. Generate bitstream(.bit) and .hwh file using Vivado.
5. Login to the PYNQ system via explore.
6. Copy the .bit file, .hwh file and the jupyter notebook project in [PYNQ] folder to the same directory under PYNQ system.
7. Run the jupyter notebook project step by step.

## Acknowledgment
Many thanks to _Xilinx_ for providing the hardware platform and related technological support. Thank Yuxiang Ma([@Yuxiang-Ma]) for his help in developing the original algorithm. Thank Jian Li([@simlaugh]) and Yinfeng Tang([@onealtom]) for their help in developing peripheral jupyter notebook code.
   
   [PYNQ]: <https://github.com/Sionnoeden/Ultrasound/tree/main/PYNQ>
   [Vivado]: <https://github.com/Sionnoeden/Ultrasound/tree/main/Vivado>
   [HLS]: <https://github.com/Sionnoeden/Ultrasound/tree/main/HLS>
   [@Yuxiang-Ma]: <https://github.com/Yuxiang-Ma>
   [@simlaugh]: <https://github.com/simlaugh>
   [@onealtom]: <https://github.com/onealtom>
