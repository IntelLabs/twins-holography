This code is used for hologram computer using wirtinger optimization. There are 2 important functions:

1) Wirtinger_test(mode = 'a', fname = 'test_images/fl1_gr'): used  for the actual computation. the mode paremeter is used to specify the propagation model. Currently only anisotripc propgation (mode='a') and the free space propagation (mode ='e') are supported. the Fname paremeter is the filename for  the target image.

2)run_analysis(path_iso, path_aniso, path_map='test_phases/Ply_r_114', fname = 'test_images/fl1_gr') : used for results analysis and PSNR computation.


To run it do : python pytorch_wirtinger.py
