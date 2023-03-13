#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 17:23:42 2022

@author: rasmus
"""

from tkinter import Tk, Label, Button, Entry, LEFT, Canvas, ttk
from GUI_welas_recalculation_def import browseFiles, chooser_silvio,chooser_silvio_custom,\
    file_error, error_variable_entry_normal, error_variable_entry_custom, \
        CreateToolTip, Show_Plot, Show_Old_Plot


file_paths = ""

class MyFirstGUI:
    def __init__(self, master):
        self.master = master
        master.title("Recalculation of refraction index and resolution for welas2300 measurements")

        
        self.lbl_Entry =            Label(master, text='Please insert relevant variables:')
        self.present_ref_index =    Label(master, text='Present refraction index of the data:')
        self.new_resolution =       Label(master, text='New resolution of the data:') 
        self.new_ref_index =        Label(master, text='New refraction index of the data:') 
        
        class ScrollableFrame(ttk.Frame):
            def __init__(self, container, *args, **kwargs):
                super().__init__(container, *args, **kwargs)
                canvas = Canvas(self, width=1350, height=300)
                scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
                self.scrollable_frame = ttk.Frame(canvas)

                self.scrollable_frame.bind(
                    "<Configure>",
                    lambda e: canvas.configure(
                        scrollregion=canvas.bbox("all")
                    )
                )

                canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

                canvas.configure(yscrollcommand=scrollbar.set)

                canvas.pack(side="left", fill="both", expand=True)
                scrollbar.pack(side="right", fill="y")

        self.frame = ScrollableFrame(master)
        self.Canvas_Label = ttk.Label(self.frame.scrollable_frame, text="")
        self.Canvas_Label.grid(row=0, column=0)
        
        self.progression =          Label(master, text='Select recalculated data to show:',justify=LEFT)
        self.show_old_data_label =  Label(master, text="Select file to examine:",justify=LEFT)
        
    
            

        self.t_present_ref_index =  Entry(bd=3)
        self.t_present_ref_index.insert(-1, '1.59')
        self.t_new_resolution =     Entry(bd=3)
        self.t_new_resolution.insert(-1, '59')
        self.t_new_ref_index =      Entry(bd=3)
        self.t_new_ref_index.insert(-1, '1.38')
        self.t_select_show_data =      Entry(bd=3, width=4)
        self.t_select_show_data.insert(-1, '1')
        self.t_select_show_old_data =      Entry(bd=3, width=4)
        self.t_select_show_old_data.insert(-1, '1')
        
        
        
        
        
        self.close_button =                 Button(master,\
                                                   text="Close",\
                                                   command=master.quit)
            
        self.file_button =                  Button(master,\
                                                   text="file explorer",\
                                                   command=lambda: [browseFiles(self),\
                                                                    self.show_old_data_label.place(x=50,y=375),\
                                                                        self.Show_Old_Plot_button.place(x=500,y=375),\
                                                                            self.t_select_show_old_data.place(x=280,y=375)])
            
        self.datein_umwandeln_button =      Button(master,\
                                                   text="Dateien umwandeln",\
                                                   command=lambda: [chooser_silvio(self),\
                                                                    file_error,error_variable_entry_normal(self),\
                                                                        self.progression.place(x=700,y=375),\
                                                                                  self.t_select_show_data.place(x=930, y=375),\
                                                                                      self.Show_Plot_button.place(x=1050, y=375)])
            
        
        self.Show_Old_Plot_button =             Button(master,\
                                                    text="Show Old Data",\
                                                    command=lambda: Show_Old_Plot(self))   
            
        self.Show_Plot_button =             Button(master,\
                                                   text="Show Plot",\
                                                   command=lambda: Show_Plot(self))
                
            
            
        self.L_ref_index_dry =        Label(master, text='Refraction index of the dry particles:')      
        self.L_ref_index_solvent =        Label(master, text='Refraction index of the solvent:')
        self.L_diameter_dry =        Label(master, text='Diameter of the dry particles:')      
        self.L_density_dry =        Label(master, text='Density of the dry particles:') 
        self.L_density_solvent =        Label(master, text='Density of the solvent:')     
        
        
        self.t_ref_index_dry =   Entry(bd=3)
        self.t_ref_index_dry.insert(-1,"1.54")
        self.t_ref_index_solvent =   Entry(bd=3)
        self.t_ref_index_solvent.insert(-1,"1.33")
        self.t_diameter_dry =   Entry(bd=3)
        self.t_diameter_dry.insert(-1,"400e-9")
        self.t_density_dry =   Entry(bd=3)
        self.t_density_dry.insert(-1,"2160")
        self.t_density_solvent =   Entry(bd=3)
        self.t_density_solvent.insert(-1,"1000")
        
        
            
        self.custom_ref_index_button =      Button(master,\
                                                   text="Custom refraction index",\
                                                   command=lambda: [self.L_ref_index_dry.place(x=700, y=75),\
                                                                            self.t_ref_index_dry.place(x=930, y=75),\
                                                                                self.L_ref_index_solvent.place(x=700, y=125),\
                                                                                    self.t_ref_index_solvent.place(x=930, y=125),\
                                                                                        self.L_diameter_dry.place(x=700, y=175),\
                                                                                            self.t_diameter_dry.place(x=930, y=175),\
                                                                                                self.L_density_dry.place(x=700,y=225),\
                                                                                                    self.t_density_dry.place(x=930,y=225),\
                                                                                                        self.L_density_solvent.place(x=700,y=275),\
                                                                                                            self.t_density_solvent.place(x=930,y=275),\
                                                                                                                self.datein_umwandeln_custom_button.place(x=930,y=325)])
                                                                                                                                                                       
            
        
        
        
        #print(ref_index_dry)
        
        self.datein_umwandeln_custom_button =     Button(master,\
                                                         text="Dateien umwandeln with custom refraction index",\
                                                         command=lambda: [chooser_silvio_custom(self),\
                                                                          file_error,error_variable_entry_custom(self),\
                                                                              self.progression.place(x=700,y=375),\
                                                                                  self.t_select_show_data.place(x=930, y=375),\
                                                                                      self.Show_Plot_button.place(x=1050, y=375)])
        
        
        
        self.lbl_Entry.place(x=50, y=25)
        self.close_button.place(x=1300, y=25)
        self.present_ref_index.place(x=50, y=75)
        self.t_present_ref_index.place(x=280, y=75)
        self.new_resolution.place(x=50, y=125)
        self.t_new_resolution.place(x=280, y=125)
        self.new_ref_index.place(x=50, y=175)
        self.t_new_ref_index.place(x=280, y=175)
        
        self.custom_ref_index_button.place(x=500,y=175)
        self.file_button.place(x=50, y=325)
        self.datein_umwandeln_button.place(x=280,y=325)
        self.frame.place(x=50, y=425)
        
        
        
        
        CreateToolTip(self.present_ref_index,\
                      text = 'Present refraction index of the data:\n\
                          Normaly refraction index of latex calibration particles 1.59,\n\
                              other possibilities are destilled water with 1.33,\n\
                                  sodium chloride with 1.54.')
                                  
        CreateToolTip(self.new_resolution,\
                      text = 'New resolution of the data:\n\
                          New resolution schould be as good or worse than input resolution.\n\
                              Normal Resolutions are 59bins (32bit), 117bins (64bit) and 234 bins (128bit).')
                              
        CreateToolTip(self.new_ref_index,\
                      text = 'New refraction index of the data:\n\
                          Should be choosen if the refraction index is not dependent on the size of the particle.\n\
                              Possibilities are destilled water with 1.33,\n\
                                  sodium chloride with 1.54, etc..')
                                  
        CreateToolTip(self.datein_umwandeln_button,\
                      text = 'Creates a new directory in the directory of the choosen files,\n\
                          recalculates the choosen files with the given spezifications and \n\
                              saves the choosen files with adaped names in the new directory.\n\
                                  The format of the new files is changed. Measuremnts got cumulated bin-wise, \n\
                                      with the changes from 100sec to 100sec step beeing lost.\n\
                                      This could be changed in the sourcecode without much trouble if needed.')
                                  
        CreateToolTip(self.custom_ref_index_button,\
                      text = 'Opens Menu to define a refraction index that is dependent on the size of the particle.\n\
                          The assumption is, that a dry salt particle gets inserted into the measurement chamber \n\
                              and than grows by absorbing a solvent, normaly water')
                              
        CreateToolTip(self.L_ref_index_dry,\
                      text = "New refraction index of the the dry aerosol particle, normaly a salt:\n\
                              Possibile is for example sodium chloride with 1.54, etc..")
                              
        CreateToolTip(self.L_ref_index_solvent,\
                      text = "New refraction index of the the solvent, normaly a destillated water:\n\
                              Possibile is for example destillated water with 1.33, etc..")
                              
        CreateToolTip(self.L_diameter_dry,\
                      text = "Diameter of the dry particle, normaly a salt crystall, in meter:\n\
                              Possibile is for example a dyameter of 400nm , wittten as 400e-9.")
                              
        CreateToolTip(self.L_density_dry,\
                      text = "Density of the dry particle, normaly a salt crystall, in kg per cubic-meter:\n\
                              Possibile is for example a density of 2160 kg per cubic-meter for sodium chloride.")
                              
        CreateToolTip(self.L_density_solvent,\
                      text = "Density of the solvent, in kg per cubic-meter:\n\
                              Possibile is for example a density of 1000 kg per cubic-meter for destillated water.")
                              
        CreateToolTip(self.datein_umwandeln_custom_button,\
                      text = 'Creates a new directory in the directory of the choosen files,\n\
                          recalculates the choosen files with the given spezifications for the custom refraction index and \n\
                              saves the choosen files with adaped names in the new directory.\n\
                                  The format of the new files is changed. Measuremnts got cumulated bin-wise, \n\
                                      with the changes from 100sec to 100sec step beeing lost.\n\
                                      This could be changed in the sourcecode without much trouble if needed.')
            
        CreateToolTip(self.progression,\
                      text = "Select recalculated data to show:\n\
                              Type in the number of the file displyed below, counted from the top, as an integer.\n\
                                  For example 1, 2, or 3.")
        
        CreateToolTip(self.Show_Plot_button,\
                      text = "Shows the selected Plot in a pop-up window.")
            
        CreateToolTip(self.show_old_data_label,\
                      text = "Select Old data file to show:\n\
                              Type in the number of the file displyed below, counted from the top, as an integer.\n\
                                  For example 1, 2, or 3.")
            
        CreateToolTip(self.Show_Old_Plot_button,\
                      text = "Shows the selected Plot in a pop-up window.")
        
        pad=3
        self._geom='20bd=30x200+0+0'
        master.geometry("{0}x{1}+0+0".format(
            master.winfo_screenwidth()-pad, master.winfo_screenheight()-pad))
        master.bind('<Escape>',self.toggle_geom)  
        

    
        
    def toggle_geom(self,event):
        geom=self.master.winfo_geometry()
        print(geom,self._geom)
        self.master.geometry(self._geom)
        self._geom=geom

root = Tk()
my_gui = MyFirstGUI(root)
root.mainloop()