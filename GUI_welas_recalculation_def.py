#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 18:14:01 2022

@author: rasmus
"""
from tkinter import filedialog
import miepython
import numpy as np
import codecs
import os
from glob import glob
#import time
from matplotlib import pyplot as plt
#from matplotlib.ticker import FormatStrFormatter
from tkinter import messagebox, Toplevel, Label, SOLID, LEFT,  Canvas, NW, ttk
from PIL import Image, ImageTk

#%%


def browseFiles(self):
      filenames = filedialog.askopenfilenames(initialdir = "/",title = "Select a File",filetypes = (("Text files","*.txt"),("all files","*.*")))
      with open('filenames_file.txt', 'w') as f:
            f.write(str(filenames))
      
      filenames = list(filenames)
      filenames.sort()
      
      file_name_sep =[]  
      for  number,filename in enumerate(filenames):
          name = filename.split("/")
          file_name_sep.append(str(number+1)+ str("\t")+ name[-1]+str("\n"))
      #print(file_name_sep)
          
      #sort = list(file_name_sep)
      #sort.sort()
      #file_paths = lambda file_name_sep: "".join(file_name_sep)
      #file_paths = lambda filenames: "".join( sort(list("".join(i.split("/")[-1] ,"/n"))) for i in filenames)
      #file_paths = lambda sort: "".join(tuple(i +"\t"+str(number+1)+"\n") for number,i in enumerate(sort))
      
      for number,filename in enumerate(filenames):
          if number==0:
              file_name = filename.split("/")
              path = file_name[:-1]
              path = "".join(i +"/" for i in path)
              break
      
  
      #self.Canvas_Label = ttk.Label(self.frame.scrollable_frame, text="Selected files from Folder: \n"+str(path)+" :\n" +str("".join(file_name_sep))).grid(row=0, column=0)
      self.Canvas_Label['text'] = "Selected files from Folder: \n"+str(path)+" :\n" +str("".join(file_name_sep))
      #return(file_paths)
  
    
#%%  
    
def file_error():
    try:
        with open('filenames_file.txt', 'r') as f:
            f = f.read()
            f = f[ : -1]
            filenames = f.split(",")
            if len(filenames)==2 and "" in filenames[-1]:
                del filenames[-1]
         
        
                
        for number,filename in enumerate(filenames):
            filename = filename[ 2: -1]
            f=codecs.open(filename,'r',encoding="ISO-8859-1")
            
            dm = []     # mid of diameter bin
            dN = []     # number of particles in bin
        
            # indicators if lines for du, dm, do, dN have been found
            foundX   = False
            founddN = False
            
            for number,line in enumerate(f):
        
                if ('X [µm]' in line) & (foundX == False):
                    dm = line.replace('X [µm]','').strip().split('\t')
                    dm = 1e-6*np.asarray([float(x) for x in dm])
                    foundX = True
        
                if ('dN [P]' in line) & (founddN == False):
                    dN = line.replace('dN [P]','').strip().split('\t')
                    dN = np.asarray([float(x) for x in dN])
                    founddN = True
            
                if  foundX & founddN:
                    break
                if number == 100:
                    messagebox.showinfo(message="File error, file propably does not have the right format.")
                
            
            f.close()
            if  not foundX & founddN:
                messagebox.showinfo(message="File error, file propably does not have the right format.")
  
    except ValueError as e:
        messagebox.showerror(message='File error, file propably does not have the right format:\n"{}"'.format(e))
       
        
#%%       
        
            
def error_variable_entry_normal(self):
        try:
            
            old_ref = float(self.t_present_ref_index.get())
            new_res = int(self.t_new_resolution.get())
            new_ref = float(self.t_new_ref_index.get())
            
        except ValueError as e:
            messagebox.showerror(message='Variable entry error: The format of the given variable is propably faulty. \
                                 The format should be float for old refraction index, e.g. 1.59, \
                                     integer for new resolution, e.g. 59, \
                                         float for new refraction index, e.g. 1.38\n \
                                             Error:"{}"'.format(e))

 #%% 

  
def error_variable_entry_custom(self):
        try:
            
            old_ref = float(self.t_present_ref_index.get())
            new_res = int(self.t_new_resolution.get())
            ref_index_dry = float(self.t_ref_index_dry.get())
            ref_index_solvent = float(self.t_ref_index_solvent.get())
            diameter_dry = float(self.t_diameter_dry.get())
            density_dry = int(self.t_density_dry.get())
            density_solvent = int(self.t_density_solvent.get())
            
        except ValueError as e:
            messagebox.showerror(message='Variable entry error: The format of the given variable is propably faulty. \
                                 The format should be float for old refraction index, e.g. 1.59, \
                                     integer for new resolution, e.g. 59, \
                                         float for new refraction index dry, e.g. 1.54, \
                                             float for new refraction index solvent, e.g. 1.33 \
                                                 float for dry diameter, e.g. 400e-9 \
                                                     integer for dry density, e.g. 2160 \
                                                         integer for solvent density, e.g. 1000\n \
                                                             Error:"{}"'.format(e))



#%%
    
def chooser_silvio(self):
    
    mWelas = float(self.t_present_ref_index.get())
    nDs_get = int(self.t_new_resolution.get())
    mNew = float(self.t_new_ref_index.get())
    
    
    with open('filenames_file.txt', 'r') as f:
        f = f.read()
        f = f[ : -1]
        filenames = f.split(",")
        if len(filenames)==2 and "" in filenames[-1]:
            del filenames[-1]

            
    
    
    for number,filename in enumerate(filenames): 
        #creates a new directory in the directory of the original files
        if number==0:
            path_tuple = filename.split("/")
            path_tuple = path_tuple[1:]
            path_tuple = path_tuple[:-1]
            path_tuple = [ '/' + word  for word in path_tuple]
            path_tuple.append("/Files_with_new_refraction_index_and_resolution")
            path = "".join(path_tuple)
            face_id = 1
            while os.path.exists(path+"%s.1" % face_id):
                face_id += 1
            os.mkdir(path+"%s.1" % face_id)
        else: 
            break
            
       
    for number,filename in enumerate(filenames):
        
        filename = filename[ 2: -1]
        f=codecs.open(filename,'r',encoding="ISO-8859-1")
        
        dm = []     # mid of diameter bin
        dN = []     # number of particles in bin
        dNneu = []
    
        # indicators if lines for du, dm, do, dN have been found
        foundX   = False
        founddN = False
        
        for number,line in enumerate(f):
    
            if ('X [µm]' in line) & (foundX == False):
                dm = line.replace('X [µm]','').strip().split('\t')
                dm = 1e-6*np.asarray([float(x) for x in dm])
                foundX = True
    
            if ('dN [P]' in line):
                dNneu = line.replace('dN [P]','').strip().split('\t')
                dNneu = np.asarray([float(x) for x in dNneu])
                if not founddN:
                    dN = dNneu
                else:
                    dN = dN + dNneu
                founddN = True
        
            #if  foundX & founddN:
            #    break
        
        f.close()
        if not foundX & founddN:
            messagebox.showinfo(message="File error, file propably does not have the right format.")
     
        #mWelas = 1.59           # refractive index set in Welas evaluation
        #mNew   = 1.38           # refractive index for new evaluation
        
        
        nLambdas  = len(dN)#59#234#256         # number of wavelengths
        lambdaMin = 200.0e-9    # minimum wavelength [m]
        lambdaMax = 2500e-9     # maximum wavelength [m]
        
        nDs  = len(dN)#59#234#256              # number of diameters
        dMin = 100.0e-9         # minimum diameter [m]
        dMax = 10.0e-6          # maximum diameter [m]
        
        muRes = 1               # angle resolution [deg]
        muMin = 78.0            # minimum angle [deg]
        muMax = 102.0           # maximum angle [deg]
        
        # do not change the stuff below
        d = np.logspace(np.log10(dMin),np.log10(dMax),nDs)  # diameter range
        lmb = np.linspace(lambdaMin,lambdaMax,nLambdas)     # wavelength range
        theta = np.arange(muMin,muMax+1,muRes)              # angle range
        mu = np.cos(theta*np.pi/180.)                       # angle range in rad
        
        
        
        # specific spectral emission of wavelengths of xenon arc lamp
        c = 2.997e8     # speed of light [m/s]
        h = 6.626e-34   # Planck quantum [J/s]
        k = 1.381e-23   # Boltzmann constant [J/K]
        T = 5500.0      # black body temperature [K]
        
        c1 = 2.0*np.pi*h*c**2.0
        c2 = h*c/k
        
        specEm = c1/(pow(lmb,5.0)*(np.exp(c2/(lmb*T))-1))
        
        # initialise original and new transfer functions; original transfer function
        # uses the ref index specified as mWelas; one new trans. func. uses the refr.
        # index specified as mNew; second trans. func. calculates refr. index based on
        # mass fraction of water in the particle from the refractive indices of NaCl
        # and H2O
        
        transFuncOrig = np.zeros(nDs)
        transFuncNew = np.zeros(nDs)
        
        # loop over all diameters of the diameter range specified above and calculate
        # the intensities for each diameter for all three transfer functions
        for i in range(0,nDs):
        
           
        
            # loop over all wavelengths in range specified above and sum up intensities
            # for the current diameter
            for j in range(0,nLambdas):
                # optical parameter
                x = np.pi * d[i] / lmb[j]
        
                # calculate intensity for current diameter-wavelength-combination and
                # refr. index of Welas
                S1, S2 = miepython.mie_S1_S2(mWelas, x, mu)
                S11 = (pow(abs(S1),2.0) + pow(abs(S2),2.0))
                sigma_sca = miepython.mie(mWelas,x)[1]
                transFuncOrig[i] += np.sum(S11)*specEm[j]*(d[i]*sigma_sca)**2.0
        
                # calculate intensity for current diameter-wavelength-combination and
                # new refr. index
                S1, S2 = miepython.mie_S1_S2(mNew, x, mu)
                S11 = (pow(abs(S1),2.0) + pow(abs(S2),2.0))
                sigma_sca = miepython.mie(mNew,x)[1]
                transFuncNew[i] += np.sum(S11)*specEm[j]*(d[i]*sigma_sca)**2.0
        
                   
        
  
        
        # initialise arrays for new diameters
        dmNew = np.zeros(len(dm))
        
        dmNew_klein = np.zeros(nDs_get)
        d_klein = np.logspace(np.log10(dMin),np.log10(dMax),nDs_get)  # diameter range
        # loop over all diameter bin mids of the distribtion from Welas read in above
        for i in range(len(dm)):
        
            # calculate the intensity for this diameter with original refr. index
            intensityOrig = 0
            for j in range(0,nLambdas):
                x = np.pi * dm[i] / lmb[j]
        
                S1, S2 = miepython.mie_S1_S2(mWelas, x, mu)
                S11 = (pow(abs(S1),2.0) + pow(abs(S2),2.0))
                sigma_sca = miepython.mie(mWelas,x)[1]
                intensityOrig += np.sum(S11)*specEm[j]*(dm[i]*sigma_sca)**2.0
        
            # interpolate the diameters of new and new+adjusted transfer functions from
            # the original intensity
            dmNew[i] = np.interp(intensityOrig,transFuncNew,d)
            
        dmNew_klein = np.interp(d_klein,dmNew,dN)
        while len(dmNew_klein) <len(dm):
            dmNew_klein= np.append(dmNew_klein,np.nan)
            d_klein= np.append(d_klein,np.nan)
            
        safe_file = np.stack((dm, dN, dmNew, dN, d_klein, dmNew_klein), axis=-1)
        path_tuple = filename.split("/")
        old_filename = path_tuple[-1]
        with open(path+"%s.1" % face_id+"/New_n" +str(mNew)+"_"+str(nDs)+"bins_"+ old_filename, 'w') as f:
            f.write(str("Old_Particle_Diameter_[nm]\tOld_Number_Concentration\tNew_Particle_Diameter_[nm]_(old resolution)\tNew_Number_Concentration_(old resolution)\tNew_Particle_Diameter_[nm]_(new resolution)\tNew_Number_Concentration_(new resolution)\n"))
            for item in safe_file:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}\t{item[5]}\n")
        
    
    
#%%        
     
        
     
#old_ref, new_res, ref_index_dry, ref_index_solvent, diameter_dry, density_dry, density_solvent        
def chooser_silvio_custom(self):
    mWelas = float(self.t_present_ref_index.get())
    nDs_get = int(self.t_new_resolution.get())
    mNaCl = float(self.t_ref_index_dry.get())
    mH2O = float(self.t_ref_index_solvent.get())
    dpdry = float(self.t_diameter_dry.get())
    rhoSol = int(self.t_density_dry.get())
    rhoLiq = int(self.t_density_solvent.get())
    
    with open('filenames_file.txt', 'r') as f:
        f = f.read()
        f = f[ : -1]
        filenames = f.split(",")
        if len(filenames)==2 and "" in filenames[-1]:
            del filenames[-1]

            
    
    
    for number,filename in enumerate(filenames): 
        #creates a new directory in the directory of the original files
        if number==0:
            path_tuple = filename.split("/")
            path_tuple = path_tuple[1:]
            path_tuple = path_tuple[:-1]
            path_tuple = [ '/' + word  for word in path_tuple]
            path_tuple.append("/Files_with_new_refraction_index_and_resolution")
            path = "".join(path_tuple)
            face_id = 1
            while os.path.exists(path+"%s.1" % face_id):
                face_id += 1
            os.mkdir(path+"%s.1" % face_id)
        else: 
            break
            
    for number,filename in enumerate(filenames):
        filename = filename[ 2: -1]
        f=codecs.open(filename,'r',encoding="ISO-8859-1")
        
        dm = []     # mid of diameter bin
        dN = []     # number of particles in bin
        dNneu = []
    
        # indicators if lines for du, dm, do, dN have been found
        foundX   = False
        founddN = False
        
        for number,line in enumerate(f):
    
            if ('X [µm]' in line) & (foundX == False):
                dm = line.replace('X [µm]','').strip().split('\t')
                dm = 1e-6*np.asarray([float(x) for x in dm])
                foundX = True
    
            if ('dN [P]' in line):
                dNneu = line.replace('dN [P]','').strip().split('\t')
                dNneu = np.asarray([float(x) for x in dNneu])
                if not founddN:
                    dN = dNneu
                else:
                    dN = dN + dNneu
                founddN = True
        
            #if  foundX & founddN:
            #    break
        
        f.close()
        if not foundX & founddN:
            messagebox.showinfo(message="File error, file propably does not have the right format.")
      
        #mWelas = 1.59           # refractive index set in Welas evaluation
        #mNaCl  = 1.54           # refractive index of NaCl
        #mH2O   = 1.33           # refractive index water
        #dpdry  = 400e-9         # dry diameter
        #rhoSol = 2160           # density of NaCl
        #rhoLiq = 1000           # density of H2O
        
        nLambdas  = len(dN)#59#234#256         # number of wavelengths
        lambdaMin = 200.0e-9    # minimum wavelength [m]
        lambdaMax = 2500e-9     # maximum wavelength [m]
        
        nDs  = len(dN)#59#234#256              # number of diameters
        dMin = 100.0e-9         # minimum diameter [m]
        dMax = 10.0e-6          # maximum diameter [m]
        
        muRes = 1               # angle resolution [deg]
        muMin = 78.0            # minimum angle [deg]
        muMax = 102.0           # maximum angle [deg]
        
        # do not change the stuff below
        d = np.logspace(np.log10(dMin),np.log10(dMax),nDs)  # diameter range
        lmb = np.linspace(lambdaMin,lambdaMax,nLambdas)     # wavelength range
        theta = np.arange(muMin,muMax+1,muRes)              # angle range
        mu = np.cos(theta*np.pi/180.)                       # angle range in rad
        
        mSol = rhoSol*np.pi*(dpdry**3.0)/6.0    # dry particle mass
        
        # specific spectral emission of wavelengths of xenon arc lamp
        c = 2.997e8     # speed of light [m/s]
        h = 6.626e-34   # Planck quantum [J/s]
        k = 1.381e-23   # Boltzmann constant [J/K]
        T = 5500.0      # black body temperature [K]
        
        c1 = 2.0*np.pi*h*c**2.0
        c2 = h*c/k
        
        specEm = c1/(pow(lmb,5.0)*(np.exp(c2/(lmb*T))-1))
        
        # initialise original and new transfer functions; original transfer function
        # uses the ref index specified as mWelas; one new trans. func. uses the refr.
        # index specified as mNew; second trans. func. calculates refr. index based on
        # mass fraction of water in the particle from the refractive indices of NaCl
        # and H2O
        
        transFuncOrig = np.zeros(nDs)
        transFuncNewAdj = np.zeros(nDs)
        
        # loop over all diameters of the diameter range specified above and calculate
        # the intensities for each diameter for all three transfer functions
        t_mAdj = np.zeros( (nDs,) )
        for i in range(0,nDs):
        
            # calculate the liquid mass and solid mass fraction in the particle
            mLiq = max(0.0, rhoLiq*np.pi*((d[i])**3.0 - dpdry**3.0)/6.0)
            alpha = mSol/(mSol+mLiq)
        
            # use rule of mixture to calculate adjusted refr. index
            mAdj = alpha*mNaCl + (1.0-alpha)*mH2O
            t_mAdj[i] = mAdj
        
            # loop over all wavelengths in range specified above and sum up intensities
            # for the current diameter
            for j in range(0,nLambdas):
                # optical parameter
                x = np.pi * d[i] / lmb[j]
        
                # calculate intensity for current diameter-wavelength-combination and
                # refr. index of Welas
                S1, S2 = miepython.mie_S1_S2(mWelas, x, mu)
                S11 = (pow(abs(S1),2.0) + pow(abs(S2),2.0))
                sigma_sca = miepython.mie(mWelas,x)[1]
                transFuncOrig[i] += np.sum(S11)*specEm[j]*(d[i]*sigma_sca)**2.0
        
                
        
                # calculate intensity for current diameter-wavelength-combination and
                # new mass ratio adjusted refr. index
                S1, S2 = miepython.mie_S1_S2(mAdj, x, mu)
                S11 = (pow(abs(S1),2.0) + pow(abs(S2),2.0))
                sigma_sca = miepython.mie(mAdj,x)[1]
                transFuncNewAdj[i] += np.sum(S11)*specEm[j]*(d[i]*sigma_sca)**2.0     
        
 
        # initialise arrays for new diameters
        dmNewAdj = np.zeros(len(dm))
        
        dmNewAdj_klein = np.zeros(nDs_get)
        d_klein = np.logspace(np.log10(dMin),np.log10(dMax),nDs_get)  # diameter range
        # loop over all diameter bin mids of the distribtion from Welas read in above
        for i in range(len(dm)):
        
            # calculate the intensity for this diameter with original refr. index
            intensityOrig = 0
            for j in range(0,nLambdas):
                x = np.pi * dm[i] / lmb[j]
        
                S1, S2 = miepython.mie_S1_S2(mWelas, x, mu)
                S11 = (pow(abs(S1),2.0) + pow(abs(S2),2.0))
                sigma_sca = miepython.mie(mWelas,x)[1]
                intensityOrig += np.sum(S11)*specEm[j]*(dm[i]*sigma_sca)**2.0
        
            # interpolate the diameters of new and new+adjusted transfer functions from
            # the original intensity
            #dmNewAdj[i] = np.interp(intensityOrig,transFuncNewAdj,d) #from silvio, should be wrong according to dokumentation
            dmNewAdj[i] = np.interp(intensityOrig,transFuncNewAdj,d)
        
        dmNewAdj_klein = np.interp(d_klein,dmNewAdj,dN)
        while len(dmNewAdj_klein) <len(dm):
            dmNewAdj_klein = np.append(dmNewAdj_klein,np.nan)
            d_klein = np.append(d_klein, np.nan)
            
        safe_file = np.stack((dm, dN, dmNewAdj, dN, d_klein, dmNewAdj_klein), axis=-1)
        path_tuple = filename.split("/")
        old_filename = path_tuple[-1]
        
        with open(path+"%s.1" % face_id+"/New_custom_ref_index_"+str(nDs_get)+"bins_"+ old_filename, 'w') as f:
            f.write(str("Old_Particle_Diameter_[nm]\tOld_Number_Concentration\tNew_Particle_Diameter_[nm]_(old resolution)\tNew_Number_Concentration_(old resolution)New_Particle_Diameter_[nm]_(new resolution)\tNew_Number_Concentration_(new resolution)\n"))
            for item in safe_file:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}\t{item[5]}\n")
 
#%%        
 
def Show_Plot(self):   
    mWelas = float(self.t_present_ref_index.get())
    new_res = int(self.t_new_resolution.get())
    selected = int(self.t_select_show_data.get())
    nDs_get = int(self.t_new_resolution.get())


    
    with open('filenames_file.txt', 'r') as f:
        f = f.read()
        f = f[ : -1]
        filenames = f.split(",")
        if len(filenames)==2 and "" in filenames[-1]:
            del filenames[-1]

  
    
    for number,filename in enumerate(filenames): 
        #creates a new directory in the directory of the original files
        if number==0:
            path_tuple = filename.split("/")
            path_tuple = path_tuple[1:]
            path_tuple = path_tuple[:-1]
            path_tuple = [ '/' + word  for word in path_tuple]
            path_tuple.append("/Files_with_new_refraction_index_and_resolution")
            path = "".join(path_tuple)
            paths = glob(path+'*/')
            raw_list =[]
            #print(paths)
            for i in paths:
                if i.find('Files_with_new_refraction_index_and_') != -1:
                    raw_list.append(i)
                else:
                    continue
                
            raw_list.sort()
            #print(raw_list)
            directory = raw_list[-1]
           
            file_names=[]
            for filename in filenames:
                
                path_tuple = filename.split("/")
                path_tuple = str(path_tuple[-1])
                file_names.append(path_tuple[:-1])
                
                
            file_names.sort()
            file = file_names[selected-1]
            txt_files = [f for f in os.listdir(directory) if f.endswith(file)]
            if len(txt_files) != 1:
                raise ValueError('should be only one txt file in the current directory')

            filename = txt_files[0]

            #f=codecs.open(str(directory)+'/' + "New_custom_ref_index_"+str(nDs_get)+"bins_"+str(file),'r',encoding="cp1250", errors="ignore")#"ISO-8859-1")
            
            x_ges =[]
            with open(directory + str(filename), 'r') as f:
                f = f.readlines()
                for lines in f:
                    lines[:-2]
                    x = lines.split("\t")
                    x_ges.append(x)
                    

        else: 
            break
        x_ges =x_ges[1:]
        x_ges = np.array(x_ges,float)
        
        #print(x_ges)
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8),sharex=False,sharey=False)
        
        ax.plot(x_ges[:,0]*10e8, x_ges[:,1], "-y",label = 'orig. (ref.index=' + str(mWelas) + ')')
        ax.plot(x_ges[:,0]*10e8, x_ges[:,1], "yo")
        ax.plot(x_ges[:,2]*10e8, x_ges[:,3],"-b", label = 'new (old res.' + str(len(x_ges[:,2])) + ')')
        ax.plot(x_ges[:,2]*10e8, x_ges[:,3],"bo")
        ax.plot(x_ges[:,4]*10e8, x_ges[:,5],"-k", label = "new (res.=" + str(new_res)+")")
        ax.plot(x_ges[:,4]*10e8, x_ges[:,5],"ko")
        ax.grid(True,which='both',axis='both')
        ax.set_title("particle size-number distributions", fontsize = 22)
        ax.set_xlabel('particle diameter [nm]', fontsize=18)
        ax.set_xscale('log')
        ax.set_ylabel('Particle number concentration', fontsize =18)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.legend(fontsize=18,loc='upper right')
        plt.tight_layout()
        
        face_id = 1
        while os.path.exists(str(directory)+"refrac_plot_{0}.1".format(face_id)+str(file)+".png"):
            face_id += 1
        Image_name = str(directory)+"refrac_plot_{0}.1".format(face_id)+str(file)+".png" 
        fig.savefig(fname=str(directory)+"refrac_plot_{0}.1".format(face_id)+str(file)+".png" )
    
        
        #global second
        
        # Creating a second Level
        second = Toplevel()
        second.title(file) # Rename this
        second.geometry("910x610")
        
        #Create a canvas
        canvas= Canvas(second, width= 900, height= 600)
        canvas.pack()
        # Creation of image object
        
        img= (Image.open(Image_name))
        #Resize the Image using resize method
        resized_image= img.resize((900,600), Image.ANTIALIAS)
        new_image= ImageTk.PhotoImage(resized_image)
        
        #Add image to the Canvas Items
        canvas.create_image(10,10, anchor=NW, image=new_image)

        

#%%
class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "10", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)
    
    
#%%

def Show_Old_Plot(self):   
    selected = int(self.t_select_show_old_data.get())

    with open('filenames_file.txt', 'r') as f:
        f = f.read()
        f = f[ : -1]
        filenames = f.split(",")
        if len(filenames)==2 and "" in filenames[-1]:
            del filenames[-1]

    for number,filename in enumerate(filenames): 
        #creates a new directory in the directory of the original files
        if number==0:
            path_tuple = filename.split("/")
            path_tuple = path_tuple[1:]
            path_tuple = path_tuple[:-1]
            path_tuple = [ '/' + word  for word in path_tuple]
            #path_tuple.append("/Files_with_new_refraction_index_and_resolution")
            path = "".join(path_tuple)
            
            file_names=[]
            for filename in filenames:
                
                path_tuple = filename.split("/")
                path_tuple = str(path_tuple[-1])
                file_names.append(path_tuple[:-1])
                
                
            file_names.sort()
            file = file_names[selected-1]
            

            f=codecs.open(str(path)+'/' + str(file),'r',encoding="cp1250", errors="ignore")#"ISO-8859-1")
            
            dm = []     # mid of diameter bin
            dN = []     # number of particles in bin
            dNneu = []
            alle_dN = []
        
            # indicators if lines for du, dm, do, dN have been found
            foundX   = False
            founddN = False
            
            for number,line in enumerate(f):
                #print(line)
        
                if ('X [µm]' in line) & (foundX == False):
                    dm = line.replace('X [µm]','').strip().split('\t')
                    dm = 1e-6*np.asarray([float(x) for x in dm])
                    foundX = True
        
                if ('dN [P]' in line):
                    dNneu = line.replace('dN [P]','').strip().split('\t')
                    dNneu = np.asarray([float(x) for x in dNneu])
                    
                    array_dN = np.array(dNneu)
                    alle_dN.append(array_dN)
                    if not founddN:
                        dN = dNneu
                    else:
                        dN = dN + dNneu
                    founddN = True
            
                
            
            f.close()
            if not foundX & founddN:
                messagebox.showinfo(message="File error, file propably does not have the right format.")
         
        else:
            break
        
        dN = np.asarray(dN)
        
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8),sharex=False,sharey=False)
        
        #ax.step(dm*1000, norm_cumulative_array, linewidth=2.5,label = 'cumulative (ref.index=' + str(mWelas) + "resolution=" +str(len(dm))+ ')')
        for number, i in enumerate(alle_dN):
            norm=i/sum(i)
            ax.step(dm*10e8, norm, linewidth=1, label= str(number)+". time interval")
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel("diameter X in [nm]", fontsize=18)  # Add an x-label to the axes.
        ax.set_ylabel("Normalised particle number concentration dN", fontsize=18)  # Add a y-label to the axes.
        ax.set_title("Size number concentration of  measured particles particel",fontsize = 22)  # Add a title to the axes.
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.legend(fontsize=18,loc='upper right')
        plt.tight_layout()
        
        
        
        face_id = 1
 
    
        while os.path.exists(str(path)+"/Old_data_plot_{0}.1".format(face_id) + str(file)+".png"):
            face_id += 1
        Image_name = str(path)+"/Old_data_plot_{0}.1".format(face_id) + str(file)+".png"
        fig.savefig(fname=str(path)+"/Old_data_plot_{0}.1".format(face_id) +str(file)+".png")
        #global second
        
        # Creating a second Level
        second = Toplevel()
        second.title(file) # Rename this
        second.geometry("910x610")
        
        #Create a canvas
        canvas= Canvas(second, width= 900, height= 600)
        canvas.pack()
        # Creation of image object
        
        img= (Image.open(Image_name))
        #Resize the Image using resize method
        resized_image= img.resize((900,600), Image.ANTIALIAS)
        new_image= ImageTk.PhotoImage(resized_image)
        
        #Add image to the Canvas Items
        canvas.create_image(10,10, anchor=NW, image=new_image)

    
    
