from __future__ import print_function
import torch
import PIL
import PIL.Image
import math
import numpy as np
import scipy.io as sio
from scipy.io import savemat
from scipy.io import loadmat
from skimage.metrics import peak_signal_noise_ratio as psnr1
from datetime import datetime
import os.path
from os import path
import gc
from schedules import ExponentialSchedule, ConstantSchedule

#torch.cuda.device(1)
torch.set_printoptions(precision=10)
if torch.cuda.is_available():
	device = torch.device("cuda:0")
	print('device GPU name: ',torch.cuda.get_device_name(0))
	torch.set_default_tensor_type(torch.cuda.FloatTensor)
else:
	device = torch.device("cpu")
	torch.set_default_tensor_type(torch.FloatTensor)
dtype = torch.float32

def srgb2rgb(srgb) :
    rgb = np.zeros_like(srgb, dtype = srgb.dtype)
    rgb[srgb<0.04045] = srgb[srgb<0.04045]/12.92
    rgb[srgb>=0.04045]=np.power(((srgb[srgb>=0.04045]+0.055)/1.055),2.4);
    return rgb
    
def rgb2srgb(rgb) :
    srgb = np.zeros_like(rgb, dtype = rgb.dtype)
    srgb[rgb<0.0031308]=rgb[rgb<0.0031308]*12.92;
    srgb[rgb>=0.0031308]=1.055*np.power(rgb[rgb>=0.0031308], 1.0/2.4)-0.055;
    return srgb
    
def fftshift(x, axes=None):  
    if axes is None:
        axes = tuple(range(x.dim()))
        shift = [dim // 2 for dim in x.shape]
    else:
        shift = [x.shape[ax] // 2 for ax in axes]
    return torch.roll(x, shift, axes) 

def ifftshift(x, axes=None):  
    if axes is None:
        axes = tuple(range(x.dim()))
        shift = [-(dim // 2) for dim in x.shape]
    else:
        shift = [-(x.shape[ax] // 2) for ax in axes]
    return torch.roll(x, shift, axes) 
  
# load phase maps
def load_map(name):
    D = loadmat(name)
    Ply = D['Ply'].astype(np.float32)
    return Ply

Ply = torch.tensor(load_map('test_phases/Ply_r_114'), dtype=dtype ).to(device)

def Kpl(P, y2):
    return P*Ply[y2,0]+Ply[y2,1]

def free_prop_aniso( Uin, f, wvl, d1, direction=1):
    P=Uin[0]
    L =Uin[1]
    Kr =Uin[2]
       
    pix_loc = 0;  
    M= P.size(1) # assume rect grid
    N = P.size(0)
    k = 2.0*np.pi/wvl
    fY = torch.arange(-N/2,N/2,1, dtype=dtype, device= device)/(N*d1)
    fX = torch.arange(-M/2,M/2,1, dtype=dtype, device= device)/(M*d1)

    [x2, y2] = torch.meshgrid(wvl * f * fX, wvl * f * fY, indexing='xy');
    [x0, y0] = torch.meshgrid((torch.arange(-M/2,M/2,1, dtype=dtype, device=device)+pix_loc)*d1, 
                              (torch.arange(-N/2,N/2,1, dtype=dtype, device=device)+pix_loc)*d1, indexing='xy');           

    Uout1=torch.zeros([P.size(0),P.size(1)],dtype=torch.complex64, device= device)

    Uin1all=torch.zeros([P.size(0),P.size(1)],dtype=torch.complex64, device= device)
    Uin20=L*(torch.cos(np.pi/wvl/f*(torch.square(x0)+torch.square(y0)))+1j*torch.sin(np.pi/wvl/f*(torch.square(x0)+torch.square(y0))))

    for j in range(0,N): 
        Py = Kpl(P,j)
        Uy = torch.cos(Py)+1j*torch.sin(Py)
        y = wvl*f*fY[j]
        Uin1all=(torch.cos(-1*np.pi/wvl/f*2*y*y0)+1j*torch.sin(-1*np.pi/wvl/f*2*y*y0))   		
        Uinall1=Uin20*Uin1all*Uy
        Uout1[j,:] =d1/np.sqrt(N)*torch.sum(Uinall1,0)
        val = 1/(1j*wvl*f)*d1*(fftshift(torch.fft.fft(fftshift(Uout1[j,:] ), norm='ortho'))) 
        Uout1[j,:]=val	
         
    Uin5=torch.cos(np.pi/wvl/f*(torch.square(x2) + torch.square(y2)))+1j*torch.sin(np.pi/wvl/f*(torch.square(x2) + torch.square(y2)))
    Uout = Uout1*Uin5
    Sout = torch.abs(torch.tensor(1/(1j*wvl*f)))*d1*d1

    return Uout, Sout

def ft2(g, delta1=None, delta2=None):
    if(delta1==None):
        delta1 = 1;
    if(delta2==None):
        delta2 = delta1;
    G=fftshift(g)
    G = torch.fft.fft2(G, norm="ortho")
    G = fftshift(G) * delta1*delta2;
    return G
    
def ift2(G, delta_f1=None, delta_f2=None):
    N = G.size(0);
    M = G.size(1);
    if(delta_f1==None):
        delta_f1 = 1/(N*1);
    if(delta_f2==None):
        delta_f2 = delta_f1*N/M;
    g=ifftshift(G)   
    g = torch.fft.ifft2(g, norm="ortho")
    g = ifftshift(g)* M*N * delta_f1*delta_f2;
    return g

def free_prop_s( Uin, f, wvl, d1, dir1=None ):
	pix_loc = 0;
	M = Uin.size(1); # assume rect grid
	N = Uin.size(0);
	k = 2*np.pi/wvl; # optical wavevector
	fY = torch.arange(-N/2,N/2,1, device= device) / (N * d1);
	fX = torch.arange(-M/2,M/2,1, device= device) / (M * d1);
	# observation plane coordinates
	[x2, y2] = torch.meshgrid(wvl * f * fX, wvl * f * fY, indexing='xy')
	[xp, yp] = torch.meshgrid((torch.arange(-M/2,M/2,1)+pix_loc)*d1, (torch.arange(-N/2,N/2,1)+pix_loc)*d1,indexing = 'xy')

	if (dir1==None)or (dir1>=0):
		Uout = 1/(1j*wvl*f)*ft2(Uin * (torch.cos(k/(2*f)*(torch.square(xp) + torch.square(yp)))+1j*torch.sin(k/(2*f)*(torch.square(xp) + torch.square(yp)))),d1);
		Sout = torch.abs(torch.from_numpy(np.array((1/(1j*wvl*f))*d1*d1)));
	else:
		Uout = torch.div((1j*wvl*f)*ift2(Uin, 1/(N * d1), 1/(M * d1)), (torch.cos(k/(2*f)*(torch.square(xp) + torch.square(yp)))+1j*torch.sin(k/(2*f)*(torch.square(xp) + torch.square(yp)))));
		Sout = torch.square(torch.abs(torch.from_numpy(np.array(1j*wvl*f*(1/d1)))));
                 
	return [ Uout, Sout ]

# SLM and other important parametrs of the holographic setup
def mk_hs():
    hs={} #dictionary
    hs['f'] = 1.0 #1m
    hs['wvl'] = 638e-9
    hs['slm_pp'] = 4.25e-6 # compound photonics SLM spec
    return hs

########################################################################################################training
def Wirtinger_Algorithm(Tamp,Lamp, iteration_num,P0, prop_fn, Pmap = None, lr0 = None):    
    hs=mk_hs()	
    A = P0
    H = P0.size(0)
    lr = .50 

    s = 1.0
    A=torch.autograd.Variable(A, requires_grad=True)    
    optvars = [{'params': A}]
    optimizer = torch.optim.Adam(optvars, lr=lr)
    C=torch.zeros([Tamp.size(0),Tamp.size(1)])  
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time) 	
    for i in range(0, iteration_num):

        if lr0 is not None:
            for g in optimizer.param_groups:
                g['lr'] = lr0.get(i)

        optimizer.zero_grad()         
        Ac=torch.remainder(A, 2*np.pi);
        if(prop_fn==free_prop_aniso):
            B=[Ac,Lamp,torch.ones([1])] 
            [C,L] = prop_fn(B,hs['f'], hs['wvl'], hs['slm_pp'],1) 			
        else:
            if Pmap is not None:
                Ac=Pmap(Ac, int(H/2))
            Bv = Lamp*(torch.cos(Ac)+1j*torch.sin(Ac))
            [C,L] = prop_fn(Bv,hs['f'], hs['wvl'], hs['slm_pp'],1)                
        Cpred=torch.abs(C)/L  
        loss=0.7*torch.nn.MSELoss()((s*Cpred)** 2, Tamp** 2)+ 0.3*torch.nn.MSELoss()((s*Cpred), Tamp)         
        loss.backward()
        optimizer.step()
        #if (i%10==0):
        if (i%1==0 and Pmap==None):
            print('Err(',i,'): loss=', loss.item(), ' scale=', s); 

    Ci = (torch.abs(C)/L)** 2
    Ci=torch.clamp(Ci, 0.0, 1.0)
    err = loss;
	
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    return [A, Ci, err]
 
 
def Wirtinger_test(mode = 'a', fname = 'holoanders_half'):
    # cleanup GPU memory
    gc.collect()    
    torch.cuda.empty_cache()
    target_image = PIL.Image.open(fname+'.png').convert('L')
	
    target_arr=np.array(target_image) 
    im=target_arr/255.0 
    im = srgb2rgb(im) # ATTN
    Tin = torch.tensor(im, dtype=dtype, device=device)    
    P0=torch.tensor(np.random.rand(target_arr.shape[0], target_arr.shape[1])*2*math.pi, dtype=dtype, device=device)

    Lin  =  torch.ones(Tin.size(), dtype= dtype, device = device);
    Tamp = torch.sqrt(Tin)  
    beam_nrm = torch.sum(torch.sum(torch.abs(Lin))); # normalize beam profile
    image_nrm = torch.sum(torch.sum(torch.abs(Tin)));    
    Ntarget=torch.sqrt(image_nrm/beam_nrm)
    Lamp = torch.sqrt(Lin)*Ntarget
	
    Pmap = None
    lr = None
    iterations = 200
    if mode[0]=='a':
        print("propagation model: free_prop_aniso")
        filename = fname+"_Wirtinger_free_prop_aniso.png"; 
        prop_fn=free_prop_aniso;
        iterations = 3000
        lr = ExponentialSchedule(0.5/2, 0.5/64, 3000)
    else:
        if mode[0]=='e':
            Pmap = Kpl
            iterations = 2000
            lr = ConstantSchedule(0.5/16)
        print("propagation model: free_prop_s")
        filename = fname+"_Wirtinger_free_props.png"; 
        prop_fn=free_prop_s;              

    print("number of iterations: ",iterations)
    
    #cuda
    [A, Ci, err00] = Wirtinger_Algorithm(Tamp,Lamp,iterations, P0,prop_fn, Pmap, lr); 
    #end cuda   

    Ires_np=Ci.detach().cpu().numpy();
    Ires_np = rgb2srgb(Ires_np) # ATTN
    P = A.detach().cpu().numpy();
    mdic = {"Ci": Ires_np}
    mdic["P"] = P;
    mdic["Lamp"] = Lamp.detach().cpu().numpy();    
    sio.savemat('Ci_'+mode+'_'+fname+'.mat', mdic);
    
    Ires_np=(Ires_np * 255).round().astype(np.uint8)
    print('     err :' , err00.cpu().detach().numpy())
    print('     psnr :' , psnr1(target_arr, Ires_np))
    img = PIL.Image.fromarray(Ires_np)	
    img.save(filename)  
    
    return 0

def run_analysis(path_iso, path_aniso, path_map='test_phases/Ply_r_114', fname = 'test_images/fl1_gr') :
    # load image
    target_image = PIL.Image.open(fname+'.png').convert('L')
    target_image = np.array(target_image).astype(np.float32)
    target_image /= 255.0
    # load aniso phase
    global Ply
    Ply = torch.tensor(load_map(path_map), dtype=dtype ).to(device)

    # load isotropic solution
    Di=loadmat(path_iso)
    # load anisotropic solution
    Da = loadmat(path_aniso)
    
    # define SLM properties
    hs = mk_hs()
    
    # run prop models
    # iso solution, aniso prop
    Pc=torch.remainder(torch.tensor(Di["P"], dtype = dtype), 2*np.pi);
    B=[Pc, torch.tensor(Di["Lamp"], dtype = dtype), torch.ones([1])] 
    [Cia,Lia] = free_prop_aniso(B,hs['f'], hs['wvl'], hs['slm_pp'],1)
    Cia = torch.abs(Cia)/Lia
    # iso solution, iso prop
    Bv = torch.tensor(Di["Lamp"], dtype = dtype)*(torch.cos(Pc)+1j*torch.sin(Pc))
    [Cii,Lii] = free_prop_s(Bv,hs['f'], hs['wvl'], hs['slm_pp'],1)                
    Cii = torch.abs(Cii)/Lii
    # aniso solution, aniso prop
    Pc=torch.remainder(torch.tensor(Da["P"], dtype = dtype), 2*np.pi);
    B=[Pc, torch.tensor(Da["Lamp"], dtype = dtype), torch.ones([1])] 
    [Caa,Laa] = free_prop_aniso(B,hs['f'], hs['wvl'], hs['slm_pp'],1)
    Caa = torch.abs(Caa)/Laa
    # aniso solution, iso model
    Bv = torch.tensor(Da["Lamp"], dtype = dtype)*(torch.cos(Pc)+1j*torch.sin(Pc))
    [Cai,Lai] = free_prop_s(Bv,hs['f'], hs['wvl'], hs['slm_pp'],1)                
    Cai = torch.abs(Cai)/Lai
    
    Call = [Cii, Cia, Cai, Caa]
    Cdef = ['iso-iso', 'iso-aniso', 'aniso-iso', 'aniso-aniso']
    for i in range(len(Call)):
        C = Call[i].cpu().numpy()
        I=C**2 
        I=rgb2srgb(I) # ATTN
        Imax = np.amax(I)
        print('Solution-model: ', Cdef[i], '     psnr :' , psnr1(image_true=target_image*255, image_test=I*255, data_range=255), 'max = ', Imax)
        img = PIL.Image.fromarray(np.clip(I*255/Imax, 0, 255).round().astype(np.uint8))
        filename = 'cmp-'+fname+'-'+Cdef[i]+'.png'
        img.save(filename)
        
###########################################################################################################################	
if __name__== "__main__":
    fnames = ['test_images/fl1_gr', 'test_images/f20_gr', 'test_images/f28_gr']
    # Ply = torch.tensor(load_map('test_phases/Ply_r_114'), dtype=dtype ).to(device)
    # Ply = torch.tensor(load_map('test_phases/Ply_r_100'), dtype=dtype ).to(device)
    Ply = torch.tensor(load_map('test_phases/Ply_rx_100'), dtype=dtype ).to(device)
    for i in range(3):
        for fn in fnames:
            Wirtinger_test('a20'+str(i), fn)
            Wirtinger_test('i20'+str(i), fn)
            Wirtinger_test('e20'+str(i), fn)
    # then one can run analysis to calculate images and get PSNR
    # example: run_analysis('Ci_i200_fl1_gr', 'Ci_a200_fl1_gr', 'test_phases/Ply_rx_100', 'fl1_gr')
        
        