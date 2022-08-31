a = load('test_phantom_patient1_bd10.mat');
sca_x = a.sca_x;
sca_z = a.sca_z;
ampl = a.ampl;
figure
scatter(sca_x, sca_z, 2, ampl)
