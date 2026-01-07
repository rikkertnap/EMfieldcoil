import numpy as np
from scipy.special import ellipk, ellipe
import argparse as argp
import matplotlib.pyplot as plt

# Physical constants

mu0 = 4*np.pi*1e-7

def Aphi_loop(rho, z, R, I=1.0):
    """
    Exact azimuthal vector potential of a single circular loop
    """
    k2 = 4*R*rho / ((R + rho)**2 + z**2)
    k2 = np.clip(k2, 1e-12, 1-1e-12) #  clipped away from 0 and 1 avoids elliptic-integral divergences
    k = np.sqrt(k2)

    pref = mu0 * I / (np.pi * k) * np.sqrt(R / rho)
    return pref * ((1 - k2/2)*ellipk(k2) - ellipe(k2))


def B_loop_zeroaxis(z, R, I=1.0):
    """
    Exact magnetic field components (B_rho, B_z) of a single loop on=axis limit
    """
    denom = ((R**2 + z**2)**(3/2.0))
    Bz =  mu0*I *( R**2)/ denom
    Brho= 0

    return Brho, Bz 

def B_loop(rho, z, R, I=1.0):
    """
    Exact magnetic field components (B_rho, B_z) of a single loop
    """
    k2 = 4*R*rho / ((R + rho)**2 + z**2)
    k2 = np.clip(k2, 1e-12, 1-1e-12)  # k2 clipped away from 0 and 1 avoids elliptic-integral divergences
    k = np.sqrt(k2)

    denom = np.sqrt((R + rho)**2 + z**2)

    Brho = (
        mu0*I*z / (2*np.pi*rho*denom)
        * (-ellipk(k2)
           + ((R**2 + rho**2 + z**2)/((R - rho)**2 + z**2))
           * ellipe(k2))
    )

    Bz = (
        mu0*I / (2*np.pi*denom)
        * (ellipk(k2)
           + ((R**2 - rho**2 - z**2)/((R - rho)**2 + z**2))
           * ellipe(k2))
    )
   

    return Brho, Bz


def Ephi_field(RHO,Z,R,N,L,I0,omega):
    """
        Computes Induced electric field (amplitude)
    """

    z_turns = np.linspace(-L/2, L/2, N)
    Aphi = np.zeros_like(RHO)

    for zt in z_turns:
        Aphi += Aphi_loop(RHO, Z - zt, R, I=1.0)

    # Electric field amplitude
    Ephi = omega * I0 * Aphi

    return(Ephi)


def B_field(RHO,Z,R,N,L,I0,omega):
    """
        Computes Induced magnetic field 
    """

    z_turns = np.linspace(-L/2, L/2, N)

    Brho = np.zeros_like(RHO)
    Bz = np.zeros_like(RHO)
    for zt in z_turns:
        Brhoval, Bzval = B_loop(RHO, Z - zt, R, I=1.0)
        Brho += Brhoval
        Bz += Bzval 
        
    # Magnetic field amplitude
    Brho = I0 * Brho
    Bz= I0 * Bz

    return(Brho, Bz)

def make_plot_Efield(RHO,Z,Ephi):
    """ 
        Plot induced E-field map 
    """
    plt.figure(figsize=(6, 5))
    plt.contourf(RHO, Z, Ephi, levels=40)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced azimuthal electric field $|E_\phi|$ (V/m)")
    plt.colorbar(label=r"$\mathrm{V/m}$")
    plt.tight_layout()
    plt.show()



def make_plot_Bfield(RHO,Z,Brho,Bz):
    """ 
        Plot induced B-field map 
    """
    plt.figure(figsize=(6, 5))

    Babs=np.sqrt(Brho*Brho+Bz*Bz)
    plt.contourf(RHO, Z, Babs, levels=40)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced magentic field $|B|$ (T)")
    plt.colorbar(label=r"$\mathrm{T}$")
    plt.tight_layout()
    plt.show()

def make_plot_Brhofield(RHO,Z,Brho):
    """ 
        Plot induced B-field map 
    """
    plt.figure(figsize=(6, 5))

    plt.contourf(RHO, Z, Brho, levels=40)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced radial magentic field $B_\rho$ (T)")
    plt.colorbar(label=r"$\mathrm{T}$")
    plt.tight_layout()
    plt.show()

  
def make_plot_Bzfield(RHO,Z,Bz):
    """ 
        Plot induced B-field map 
    """
    plt.figure(figsize=(6, 5))

    plt.contourf(RHO, Z, Bz, levels=40)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced  azimuthal magentic field $B_z$ (T)")
    plt.colorbar(label=r"$\mathrm{T}$")
    plt.tight_layout()
    plt.show()

  
def make_plot_Bzcenterfield(RHO,Z,Brho,Bz):
    """ 
        Plot induced B-field map 
    """
    plt.figure(figsize=(6, 5))

    plt.plot(Z, Bz)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced  azimuthal magentic field $B_z$ (T)")
    plt.tight_layout()
    plt.show()


def save_plot_Efield(RHO,Z,Ephi,dirname):
    """ 
        save Plot induced E-field map 
    """
    fig = plt.figure(figsize=(6, 5))
    plt.contourf(RHO, Z, Ephi, levels=40)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced azimuthal electric field $|E_\phi|$ (V/m)")
    plt.colorbar(label=r"$\mathrm{V/m}$")
    plt.tight_layout()
    # plot.show()
    result_data_dir=dirname+"/"
    output_file=result_data_dir+"Ephi_field.pdf"
    fig.savefig(output_file)

def save_plot_Brhofield(RHO,Z,Brho,dirname):
    """ 
        save plot induced B-field map 
    """
    fig = plt.figure(figsize=(6, 5))

    plt.contourf(RHO, Z, Brho, levels=40)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced radial magentic field $B_\rho$ (T)")
    plt.colorbar(label=r"$\mathrm{T}$")
    plt.tight_layout()
    #plt.show()
    result_data_dir=dirname+"/"
    output_file=result_data_dir+"Bpho_field.pdf"
    fig.savefig(output_file)

  
def save_plot_Bzfield(RHO,Z,Bz,dirname):
    """ 
        Plot induced B-field map 
    """
    fig = plt.figure(figsize=(6, 5))

    plt.contourf(RHO, Z, Bz, levels=40)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced  axial magentic field $B_z$ (T)")
    plt.colorbar(label=r"$\mathrm{T}$")
    plt.tight_layout()
    #plt.show()
    
    result_data_dir=dirname+"/"
    output_file=result_data_dir+"Bz_field.pdf"
    fig.savefig(output_file)


def save_plot_Bzcenterfield(RHO,Z,Brho,Bz,dirname):
    """ 
        Plot induced B-field map 
    """
    fig = plt.figure(figsize=(6, 5))

    plt.plot(Z, Bz)
    plt.xlabel(r"$\rho\ \mathrm{(m)}$")
    plt.ylabel(r"$z\ \mathrm{(m)}$")
    plt.title(r"Induced  axail magentic field $B_z$ (T)")
    plt.tight_layout()
     #plt.show()
    
    result_data_dir=dirname+"/"
    output_file=result_data_dir+"Bzcdnter_field.pdf"
    fig.savefig(output_file)



def parse_cmdline():  
    """ 
        Get coil parameters from cmdline
    """ 
    
    ap = argp.ArgumentParser(
        description="EM field of a coil in quasi-static approx."
    )

    ap.add_argument("--R", type=float, default="0.134", help="Coil radius (m) default %(default)s m")
    ap.add_argument("--N", type=int,   default="500",  help="Number of turns default %(default)s")
    ap.add_argument("--L", type=float, default="0.02", help="Axial width (m) default %(default)s m")
    ap.add_argument("--I0",type=float, default="0.05", help="Current amplitude (A) default %(default)s A")
    ap.add_argument("--f", type=float, default="200.0", help="Frequency (Hz) default %(default)s Hz")
    ap.add_argument("--save", action="store_true", help="Save figures")
    args = ap.parse_args()

    return(args)

def print_args(args):
    args_dict = vars(args) # makes a dictonary 
    for key, value in args_dict.items():
        print(f"{key}: {value}")


def main():

    # get coil parameters
    args = parse_cmdline()
    
    # print arguments
    print_args(args)

    N = args.N                # Number of turns  
    R = args.R                # Coil radius
    L = args.L                # axial width (m)
    I0 = args.I0              # current amplitude (A)
    f = args.f                # frequency (Hz)
    
    # angular frequency 
    omega = 2.0 * np.pi * f 

    # Cylindrical grid (rho, z)
    rho = np.linspace(1e-3, 0.35, 140) # rho > 0 avoids the cylindrical singularity at rho=0
    z   = np.linspace(-0.35, 0.35, 180)
    RHO, Z = np.meshgrid(rho, z)
    
    # EM fields
    Ephi = Ephi_field(RHO,Z,R,N,L,I0,omega)
    Brho, Bz = B_field(RHO,Z,R,N,L,I0,omega)
   
    # Em fields center 
    RHOc, Zc = np.meshgrid(0.001, z) # rho > 0 avoids 
    Brhoc, Bzc = B_field(RHOc,Zc,R,N,L,I0,omega)


    # output 
    if not args.save:
        make_plot_Efield(RHO,Z,Ephi)
        #make_plot_Bfield(RHO,Z,Brho,Bz)
        make_plot_Brhofield(RHO,Z,Brho)
        make_plot_Bzfield(RHO,Z,Bz)
        make_plot_Bzcenterfield(RHOc,Zc,Brhoc,Bzc) 

    if args.save: 
        # location output
        dirname="/Users/rna878/ipythonnotebook/EMfieldcoil"
        save_plot_Efield(RHO,Z,Ephi,dirname)
        save_plot_Brhofield(RHO,Z,Brho,dirname)
        save_plot_Bzfield(RHO,Z,Bz,dirname)
        save_plot_Bzcenterfield(RHOc,Zc,Brhoc,Bzc,dirname) 





if __name__ == "__main__":
    main()


