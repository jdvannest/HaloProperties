import pynbody,pickle,sys,pymp,warnings,argparse,os
import numpy as np

def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser()
parser.add_argument("-s","--simulation",choices=['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329'],required=True)
parser.add_argument("-o","--overwrite",action='store_true')
parser.add_argument("-p","--path",required=True)
parser.add_argument("-n","--numproc",required=True,type=int)
args = parser.parse_args()


#########################################################################
#Change the name of the property you're collecting (used for file naming)
property = 'ObservedShape'
image_dir = args.path+f'Images/ShapeImages/'
#########################################################################


#Load in Simulation Info (paths and halo numbers)
Sims = pickle.load(open('SimulationInfo.pickle','rb'))
#Load in the existing data file, or create new/overwrite if necessary
file_name = args.path+f'Marvel_DCJL.{property}.pickle'
if args.overwrite:
    os.system('rm '+file_name)
    Data = {}
    for sim in Sims:
        Data[sim] = {}
    print('Writing new Data File...')
else:
    try:
        Data = pickle.load(open(file_name,'rb'))
        print('Data File loaded.')
    except:
        print('No Data File found. Writing new one...')
        Data = {}
        for sim in Sims:
            Data[sim] = {}


sim  = args.simulation
#Skip if this simulation if it has already been done (unless overwriting)
if len(Data[sim])>0:
    print(f'{sim} already completed')
else:
    print(f'Loading {sim}...')
    #Create shared dictionary for multiprocessing data and shared 1-element
    # array to act as scalar for traking progress
    SimData = pymp.shared.dict()
    prog=pymp.shared.array((1,),dtype=int)
    #Load in simulation
    s = pynbody.load(Sims[sim]['path'])
    s.physical_units()
    h = s.halos()
    myprint(f'{sim} Loaded.',clear=True)
    print('\tWriting: 0%')
    #Begin parallelized analysis
    with pymp.Parallel(args.numproc) as pl:
        for i in pl.xrange(len(Sims[sim]['halos'])):
            #Load in specific halo
            halonum = Sims[sim]['halos'][i]
            halo = h[halonum]

            #The shared dictionary object does not allow writing to sub-dictionaries
            #So use 'current' as a 'middle-man' step that will be inserted as a whole
            #into the shared dictionary at the end
            current={}

            ########################################################
            #Replace this block with the desired property to collect
            from pynbody.derived import rxy,az
            import matplotlib.pylab as plt 
            from math import pi,degrees
            from scipy.optimize import curve_fit
            from pynbody.plot.sph import image
            from numpy.linalg import eig, inv
            from matplotlib.patches import Ellipse
            from numpy import sin,cos
            plt.rcParams.update({'text.usetex':False})
            #Sersic Function which is fit to the Surface-Brightness Pofiles of galaxy
            def sersic(r, mueff, reff, n):
                return mueff + 2.5*(0.868*n-0.142)*((r/reff)**(1./n) - 1)
            #Creates a g-band SB profile from the b- and v-band data
            def gbandprofile(b,v):
                B = np.array(b)
                V = np.array(v)
                return(V + 0.6*(B-V) - 0.12)
            #Conversion functions for image pixels to sim values in kpc
            def pix2kpc(pix,width):
                return(pix/1000.*width-(width/2.))
            def kpc2pix(kpc,width):
                return(int((kpc+(width/2.))/width*1000))
            
            #Ellipse functions from https://stackoverflow.com/questions/13635528/fit-a-ellipse-in-python-given-a-set-of-points-xi-xi-yi
            def EllipseFit(x,y):
                x = x[:,np.newaxis]
                y = y[:,np.newaxis]
                D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
                S = np.dot(D.T,D)
                C = np.zeros([6,6])
                C[0,2] = C[2,0] = 2; C[1,1] = -1
                E, V =  eig(np.dot(inv(S), C))
                n = np.argmax(np.abs(E))
                return( V[:,n] )
            def EllipseCenter(ellipse):
                b,c,d,f,g,a = ellipse[1]/2,ellipse[2],ellipse[3]/2,ellipse[4]/2,ellipse[5],ellipse[0]
                num = b*b-a*c
                x0=(c*d-b*f)/num
                y0=(a*f-b*d)/num
                return np.array([x0,y0])
            def EllipseAngle(ellipse):
                b,c,d,f,g,a = ellipse[1]/2,ellipse[2],ellipse[3]/2,ellipse[4]/2,ellipse[5],ellipse[0]
                return 0.5*np.arctan(2*b/(a-c))
            def EllipseAxes(ellipse):
                b,c,d,f,g,a = ellipse[1]/2,ellipse[2],ellipse[3]/2,ellipse[4]/2,ellipse[5],ellipse[0]
                up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
                down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
                down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
                res1=np.sqrt(up/down1)
                res2=np.sqrt(up/down2)
                return np.array([res1,res2])
            
            #Imaging and Fitting Function
            def FitImage(centered_halo, plot_axis, legend=True):
                #Create a smoothed SB profile to determine Reff and SB @ 2Reff
                try:
                    Rvir = pynbody.analysis.halo.virial_radius(centered_halo)
                except:
                    Rvir = 10
                prof = pynbody.analysis.profile.Profile(centered_halo.s,type='lin',min=.25,max=Rvir,ndim=2,nbins=int((Rvir-0.25)/0.1))
                #gband = gbandprofile(p['sb,b'],p['sb,v'])
                sb = prof['sb,v']
                #Smooth the sb profile
                smooth = np.nanmean(np.pad(sb.astype(float),(0,3-sb.size%3),mode='constant',constant_values=np.nan).reshape(-1,3),axis=1)
                try:
                    y = smooth[:np.where(smooth>32)[0][0]+1]
                except:
                    y = smooth
                #Create the x-data for the sb profile
                x = np.arange(len(y))*0.3 + 0.15
                x[0] = 0.05
                #Remove any NaN's from profile
                if True in np.isnan(y):
                    x = np.delete(x,np.where(np.isnan(y)==True))
                    y = np.delete(y,np.where(np.isnan(y)==True))
                #Perform Sersic fit on smoothed SB profile (if not all NaN's)
                if len(x)>0:
                    #Initial guess values for curve_fit
                    r0 = x[int(len(x)/2)]
                    m0 = np.mean(y[:3])
                    par,ign = curve_fit(sersic,x,y,p0=(m0,r0,1),bounds=([10,0,0.5],[40,100,16.5]))
                    #find the array index closest to 2*Reff (Reff is par[1] from the Sersic function)
                    ind = np.where(abs(prof['rbins']-2*par[1])==min(abs(prof['rbins']-2*par[1])))[0][0]
                    #Get the v-band luminosity density at 2*Reff
                    v = prof['v_lum_den'][ind]
                    #Create a stellar particle filter around the halo and generate image
                    width = 12*par[1]
                    sphere = pynbody.filt.Sphere(width*np.sqrt(2)*1.01)
                    im = image(s[sphere].s,qty='v_lum_den',width=width,subplot=plot_axis,units='kpc^-2',resolution=1000,show_cbar=False)
                    #Set image axes to -6 to 6 Reff and plot + marker at center
                    plot_axis.set_xlim([-width/2,width/2])
                    plot_axis.set_ylim([-width/2,width/2])
                    plot_axis.scatter(0,0,marker='+',c='k')
                    #Find the image indices if isophote (defined as having v-lum within tolerance of SB@2Reff)
                    inds,tolerance = [[[],[]],.01]
                    while len(inds[0])==0 and tolerance<.1:
                        inds = np.where((im>v*(1-tolerance)) & (im<v*(1+tolerance)))
                        tolerance+=.01
                    if len(inds[0])==0:
                        return(np.nan)
                    else:
                        #If isophote exists, redifine to be within 1% of sb @ 2Reff
                        inds = np.where((im>v*.99) & (im<v*1.01))
                        #Convert isophote image coordinates into kpc cooordinates and plot it over the image
                        iso_y = pix2kpc(inds[0],width)
                        iso_x = pix2kpc(inds[1],width)
                        plot_axis.scatter(iso_x,iso_y,c='r',s=.5**2)
                        #Fit Ellipse to isophote, and store axis ratio
                        E = EllipseFit(iso_x,iso_y)
                        cen = EllipseCenter(E)
                        phi = EllipseAngle(E)
                        a,b = EllipseAxes(E)
                        #Plot the ellipse fit on the image and set image title to axis ratio
                        plot_axis.add_patch(Ellipse(cen,2*a,2*b,angle=degrees(phi),facecolor='None',edgecolor='orange'))
                        plot_axis.plot([-a*cos(phi)+cen[0],a*cos(phi)+cen[0]],[-a*sin(phi)+cen[1],a*sin(phi)+cen[1]]
                                ,linewidth=.5,color='orange')
                        plot_axis.plot([-b*cos(phi+pi/2)+cen[0],b*cos(phi+pi/2)+cen[0]],[-b*sin(phi+pi/2)+cen[1]
                                ,b*sin(phi+pi/2)+cen[1]],linewidth=.5,color='orange')
                        plot_axis.set_title(f'{round(min([b,a])/max([b,a]),3)}')
                        if legend: 
                            plot_axis.plot([-100,-100],[-100,-100],c='r',label=r'Isophote (2 R$_{eff}$)')
                            plot_axis.plot([-100,-100],[-100,-100],c='orange',label='Ellipse Fit')
                            plot_axis.legend(loc='lower left',ncol=2,prop={'size':10})
                        #Return axis ratio
                        return(min([b,a])/max([b,a]))
                else:
                    return(np.nan)
            
            #Load in the halo, and use imaging/fitting function
            if len(halo.s)>0:
                #Create 2 panel figure, and plot faceon image in left pane and sideon image in right 
                f,ax = plt.subplots(1,2,figsize=(15,6))
                pynbody.analysis.angmom.faceon(halo)
                faceon = FitImage(halo,ax[0])
                pynbody.analysis.angmom.sideon(halo)
                sideon = FitImage(halo,ax[1],legend=False)
                current['b/a'] = faceon
                current['c/a'] = sideon
                f.savefig(image_dir+f'{sim}.{halonum}.png',bbox_inches='tight',pad_inches=.1)
            else:
                current['b/a'] = np.nan
                current['c/a'] = np.nan
            ########################################################
            

            #Print progress to terminal
            myprint(f'\tWriting: {round(float(prog[0]+1)/float(len(Sims[sim]["halos"]))*100,2)}%',clear=True)
            
            #with pl.lock:  <----- the 'lock' block probably isn't necessary since the only thing that might
            #                    encounter a race condition is the progress counter, which is only for the UI
            #Insert the current halo's data into the shared dictionary
            SimData[str(halonum)] = current
            prog[0]+=1
            del current

    #Transfer the data from the shared dictionary to the main Data File (saving the shared
    #dictionary just results in garbage data for some reason)
    for halo in Sims[sim]['halos']:
        Data[sim][str(halo)] = SimData[str(halo)]

    #Update the Data File
    out = open(file_name,'wb')
    pickle.dump(Data,out)
    out.close()
    print(f'File updated with {sim}.')