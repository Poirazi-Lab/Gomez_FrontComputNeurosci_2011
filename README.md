This is the readme file for:

Gomez Gonzalez JF, Mel BW, Poirazi P. Distinguishing Linear vs. Non-Linear Integration in CA1 Radial Oblique Dendrites: It's about Time. Front Comput Neurosci. 2011;5:44. doi: 10.3389/fncom.2011.00044

### Instructions:

1. Copy all the folders in your computer.
2. Compile the mechanisms (in folder `../mechanism`) with command `nrnivmodl` (model is currently linux/unix only)
3. This step is necesary only if some mechanism is changed by the user. In this version, the neurons are tuned with the current mechanisms.(Gómez et al, Front Comput Neurosci 2011)

   For tunning the neuron:
   a) Select the cell that you want to tune in file `../experiment/tune-synapses/Simulation_gradient.hoc` in the variable `ID_cell=1, 2, ...`
   - 1 = n123
   - 2 = n125
   - 3 = n128
   - 4 = n129
   - 5 = n130
   b) Execute the tunning (`sh ../experiment/tune-synapses/run_simulation_gradient`) and to copy output file in `../junio/tuned-nXXX-trunk/`
4. Run simulation:
   a) Select the cell that you want to study in file `proximalCA1_1.hoc` in the variable `ID_cell=1, 2...`
   - 1 = n123
   - 2 = n125
   - 3 = n128
   - 4 = n129
   - 5 = n130
    b) Select the Parameters inside `../experiment/junio/ParameterFile.hoc`
    c) Design the protocols inside `../experiment/junio/ProtocolsFile.hoc`
    d) Execute `../experiment/junio/run_ProtocolsFile.hoc`
