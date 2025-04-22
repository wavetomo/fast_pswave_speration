# fast_pswave_speration
* 1. The main program is wavefront_phase_space_wavenumber_wavefield_decomposition 
* note：The program gives three examples, such as selecting example 1 for wavefield decomposition, 
* you need to comment out the parameter reading of other examples, and then run the program to get the result.

* 2. Homogeneous model wavefield data
* vxfile='homogeneous_model_snapshot_600_600_0700ms.vx';
* vzfile='homogeneous_model_snapshot_600_600_0700ms.vz';

Layered model wavefield data
* vxfile='layered_model_snapshot_460_860_1100ms.vx';
* vzfile='layered_model_snapshot_460_860_1100ms.vz';

Hess model wavefield data
* vxfile='Hess_model_snapshot_860_1560_2600ms.vx';
* vzfile='Hess_model_snapshot_860_1560_2600ms.vz';

* 3. Function description in the main program

* read_matrix：read a matrix from a file
* modpad2d：pading 2D model parameter for absorbing boundary condition
* generate_wavenumber: generate horizontal and vertical wavenumbers
* derivate1_fd8: the 8th order accuracy finite difference to calculate the first derivative
* derivate2_fd8: the 8th order accuracy finite difference to calculate the second derivative
* perclip: control the colorbar range
