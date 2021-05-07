################################ change params ########################################################
# Material's name
file_name                                      = "WSe"

# pixel size of the image (#A)
pixel_size                                     = 0.2074

#image size of the image (#pixel)
image_size                                     = 1024

#metal site atom number
metal_atom                                     = 74

#chalcogen site atom number
chalcogen_atom                                 = 34

#lattice constant a of your sample(#A)
lattice_constant_a                             = 3.3

#If there is subsitution on metal site, input the doped atom number or input zero
doped_metal_atom                               = 34

#If there is subsitution on metal site, input the concentration(0-1) or input zero
metal_atom_concentration                       = 0.02

#If there is vacancy on metal site, input the concentration(0-1) or input zero
metal_atom_vacancy_concentration               = 0.005

#If there is subsitution on chalcogen site, input the doped atom number or input zero
doped_chalcogen_atom                           = 74

#If there is two subsitution type on chalcogen site, input the concentration(0-1) or input zero
chalcogen_atom_concentration_two_subsititution = 0.001

#If there is one subsitution type on chalcogen site, input the concentration(0-1) or input zero
chalcogen_atom_concentration_one_subsititution = 0.001

#If there is one vacancy type on chalcogen site, input the concentration(0-1) or input zero
chalcogen_atom_concentration_one_vacancy       = 0.05

#If there is two vacancy type on chalcogen site, input the concentration(0-1) or input zero
chalcogen_atom_concentration_two_vacancy       = 0.02


voltage                   = 80     #Voltage(#kV)
Cs3_param_mean            = 0      #mean of Cs3(#um)
Cs3_param_std             = 0.0005 #the std of Cs3(#mm)
Cs5_param_mean            = 0      #mean of Cs5(#mm)
Cs5_param_std             = 0      #std of Cs5(#mm)
df                        = 0      #df
aperture                  = 25.2   #aperture
ADF_angle_min             = 63     #ADF angle min
ADF_angle_max             = 200    #ADF angle max
Source_size_param_mean    = 1.05   #mean of Source size(#Ang)
Source_size_param_std     = 0.1    #std of Source size(#Ang)
defocus_spread_param_mean = 0      #mean of defocus spread(#A)
defocus_spread_param_std  = 0.2    #std of defocus spread(#A)
probe_current_param_mean  = 30     #mean of probe current(#A)
probe_current_param_std   = 2      #std of probe current(#A)
dwell_time                = 20     #dwell time(#us)


################################ DO NOT MODIFY ########################################################
import numpy as np
import time
import os

sample_param_dic = {}
sample_param_dic['file_name'] = file_name
sample_param_dic['pixel_size'] = pixel_size
sample_param_dic['image_size'] = image_size
sample_param_dic['metal_atom'] = metal_atom
sample_param_dic['chalcogen_atom'] = chalcogen_atom
sample_param_dic['lattice_constant_a'] = lattice_constant_a
sample_param_dic['doped_metal_atom'] = doped_metal_atom
sample_param_dic['metal_atom_concentration'] = metal_atom_concentration
sample_param_dic['metal_atom_vacancy_concentration'] = metal_atom_vacancy_concentration
sample_param_dic['doped_chalcogen_atom'] = doped_chalcogen_atom
sample_param_dic['chalcogen_atom_concentration_two_subsititution'] = chalcogen_atom_concentration_two_subsititution
sample_param_dic['chalcogen_atom_concentration_one_subsititution'] = chalcogen_atom_concentration_one_subsititution
sample_param_dic['chalcogen_atom_concentration_one_vacancy'] = chalcogen_atom_concentration_one_vacancy
sample_param_dic['chalcogen_atom_concentration_two_vacancy'] = chalcogen_atom_concentration_two_vacancy

sample_param_dic['lattice_constant_b'] = sample_param_dic['lattice_constant_a']*np.sqrt(3)
sample_param_dic['lattice_constant_c'] = 1
sample_param_dic['rep_x'] = int(sample_param_dic['image_size']*sample_param_dic['pixel_size']/sample_param_dic['lattice_constant_a'])
sample_param_dic['rep_y'] = int(sample_param_dic['image_size']*sample_param_dic['pixel_size']/sample_param_dic['lattice_constant_b'])
sample_param_dic['rep_z'] = 1
sample_param_dic['metal_dopant_different'] = sample_param_dic['doped_metal_atom'] - sample_param_dic['metal_atom']
sample_param_dic['chalcogen_dopant_different'] = sample_param_dic['doped_chalcogen_atom'] - sample_param_dic['chalcogen_atom']
sample_param_dic['dopant_conc_two_subsitutions_high'] = sample_param_dic['chalcogen_atom_concentration_two_subsititution']
sample_param_dic['dopant_conc_one_subsitution_high'] = sample_param_dic['chalcogen_atom_concentration_one_subsititution']+sample_param_dic['dopant_conc_two_subsitutions_high']
sample_param_dic['dopant_conc_one_vacancy_high'] = sample_param_dic['dopant_conc_one_subsitution_high']+sample_param_dic['chalcogen_atom_concentration_one_vacancy']
sample_param_dic['dopant_conc_two_vacancies_high'] = sample_param_dic['dopant_conc_one_vacancy_high'] + sample_param_dic['chalcogen_atom_concentration_two_vacancy']
sample_param_dic['supercell_a'] = sample_param_dic['lattice_constant_a'] * sample_param_dic['rep_x']
sample_param_dic['supercell_b'] = sample_param_dic['lattice_constant_b'] * sample_param_dic['rep_y']
sample_param_dic['supercell_c'] = sample_param_dic['lattice_constant_c'] * sample_param_dic['rep_z']


EM_param_dic = {}
EM_param_dic['voltage'] = voltage
EM_param_dic['Cs3_param'] = (Cs3_param_mean,Cs3_param_std)
EM_param_dic['Cs5_param'] = (Cs5_param_mean,Cs5_param_std)
EM_param_dic['df'] = df
EM_param_dic['aperture'] = aperture
EM_param_dic['ADF_angle_min'] = ADF_angle_min
EM_param_dic['ADF_angle_max'] = ADF_angle_max
EM_param_dic["Higher_order"] = 'END'
EM_param_dic['Source_size_param'] = (Source_size_param_mean, Source_size_param_std)
EM_param_dic['defocus_spread_param'] = (defocus_spread_param_mean,defocus_spread_param_std)
counting_noise = 'y'
EM_param_dic['counting_noise'] = counting_noise
EM_param_dic['probe_current_param'] = (probe_current_param_mean,probe_current_param_std)
EM_param_dic['dwell_time'] = dwell_time


def generate_files(sample_param_dic,EM_param_dic,file_num):
    '''
    generate .xyz and .param files and one batch file(.bat)
    Input: sample_param_dic, EM_param_dic, file_num(int)
    Ruturn: None
    '''
    filesuffix = '_incostem_'
    xtalnm = sample_param_dic['file_name']
    image_size = sample_param_dic['image_size']
    rep_x, rep_y, rep_z = sample_param_dic['rep_x'], sample_param_dic['rep_y'], sample_param_dic['rep_z']
    a, b, c = sample_param_dic['lattice_constant_a'], sample_param_dic['lattice_constant_b'], sample_param_dic['lattice_constant_c']
    atomZ1, atomZ2 = sample_param_dic['metal_atom'], sample_param_dic['chalcogen_atom']
    atom1=np.array([atomZ1, 0.000000, 0.000000, 1.797500, 1, 0.08])
    atom2=np.array([atomZ2, a/2, b/6, 0.000000, 1, 0.08])
    atom3=np.array([atomZ2, a/2, b/6, 3.595000, 1, 0.08])
    atom4=np.array([atomZ1, a/2, b/2, 1.797500, 1, 0.08])
    atom5=np.array([atomZ2, a, b*2/3, 0.000000, 1, 0.08])
    atom6=np.array([atomZ2, a, b*2/3, 3.595000, 1, 0.08])
    supercell_a, supercell_b, supercell_c = sample_param_dic['supercell_a'], sample_param_dic['supercell_b'], sample_param_dic['supercell_c']
    metal_doped_prob = sample_param_dic['metal_atom_concentration']
    metal_vacancy_prob = sample_param_dic['metal_atom_vacancy_concentration']
    metal_dopant_different = sample_param_dic['metal_dopant_different']
    chalcogen_dopant_different = sample_param_dic['chalcogen_dopant_different']
    dopant_conc_two_subsitutions_high = sample_param_dic['dopant_conc_two_subsitutions_high']
    dopant_conc_one_subsitution_high = sample_param_dic['dopant_conc_one_subsitution_high']
    dopant_conc_one_vacancy_high = sample_param_dic['dopant_conc_one_vacancy_high']
    dopant_conc_two_vacancies_high = sample_param_dic['dopant_conc_two_vacancies_high']
    voltage = EM_param_dic['voltage']
    Cs3_param = EM_param_dic['Cs3_param']
    Cs5_param = EM_param_dic['Cs5_param']
    df = EM_param_dic['df']
    aperture = EM_param_dic['aperture']
    ADF_angle_min = EM_param_dic['ADF_angle_min']
    ADF_angle_max = EM_param_dic['ADF_angle_max']
    Higher_order = 'END'
    Source_size_param = EM_param_dic['Source_size_param']
    defocus_spread_param = EM_param_dic['defocus_spread_param']
    counting_noise = 'y'
    probe_current_param = EM_param_dic['probe_current_param']
    dwell_time = EM_param_dic['dwell_time']

    # open is an default function to deal with files in python
    Batch_File_name = 'Batch_'+str(file_num)+'files_'+xtalnm+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'.bat'
    file_batch = open(Batch_File_name,'w+')
    for N in range(file_num):

        #image files
        filename = xtalnm+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_'+str(N)
        fid = open(filename+'.xyz', 'w+')

        #Label files, depends on how many defects you want
        filename_metal_Doped = xtalnm+'_metal_Doped'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_'+str(N)
        fid_metal_Doped = open(filename_metal_Doped+'.xyz', 'w+')
        filename_metal_vacancy = xtalnm+'_cmetal_vacancy'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_'+str(N)
        fid_metal_vacancy = open(filename_metal_vacancy+'.xyz', 'w+')
        filename_2Doped = xtalnm+'_2Doped'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_'+str(N)
        fid_2Doped = open(filename_2Doped+'.xyz', 'w+')
        filename_1Doped = xtalnm+'_1Doped'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_'+str(N)
        fid_1Doped = open(filename_1Doped+'.xyz', 'w+')
        filename_1vacancy = xtalnm+'_1vacancy'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_'+str(N)
        fid_1vacancy = open(filename_1vacancy+'.xyz', 'w+')
        filename_2vacancy = xtalnm+'_2vacancy'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_'+str(N)
        fid_2vacancy = open(filename_2vacancy+'.xyz', 'w+')

        #Write commands into the batch file
        batch_str = './incostem<'+filename+'.param \n'
        batch_metal_Doped_str = './incostem<'+filename_metal_Doped+'.param \n'
        batch_metal_vacancy_str = './incostem<'+filename_metal_vacancy+'.param \n'
        batch_2Doped_str = './incostem<'+filename_2Doped+'.param \n'
        batch_1Doped_str = './incostem<'+filename_1Doped+'.param \n'
        batch_1vacancy_str = './incostem<'+filename_1vacancy+'.param \n'
        batch_2vacancy_str = './incostem<'+filename_2vacancy+'.param \n'
        file_batch.write(batch_str)
        if metal_doped_prob>0:
            file_batch.write(batch_metal_Doped_str)
        if metal_doped_prob>0:
            file_batch.write(batch_metal_vacancy_str)
        if sample_param_dic['chalcogen_atom_concentration_two_subsititution']>0:
            file_batch.write(batch_2Doped_str)
        if sample_param_dic['chalcogen_atom_concentration_one_subsititution']>0:
            file_batch.write(batch_1Doped_str)
        if sample_param_dic['chalcogen_atom_concentration_one_vacancy']>0:
            file_batch.write(batch_1vacancy_str)
        if sample_param_dic['chalcogen_atom_concentration_two_vacancy']>0:
            file_batch.write(batch_2vacancy_str)

        #Write the header into each xyz file
        DateString = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
        xyzheader=xtalnm+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_Created_at_'+DateString+'\n'
        xyzheader_metal_Doped=xtalnm+'_metal_Doped'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_Created_at_'+DateString+'\n'
        xyzheader_metal_vacancy=xtalnm+'_metal_vacancy'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_Created_at_'+DateString+'\n'
        xyzheader_2Doped=xtalnm+'_2Doped'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_Created_at_'+DateString+'\n'
        xyzheader_1Doped=xtalnm+'_1Doped'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_Created_at_'+DateString+'\n'
        xyzheader_1vacancy=xtalnm+'_1vacancy'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_Created_at_'+DateString+'\n'
        xyzheader_2vacancy=xtalnm+'_2vacancy'+filesuffix+str(rep_x)+'_'+str(rep_y)+'_'+str(rep_z)+'_Created_at_'+DateString+'\n'
        fid.write(xyzheader)
        fid_metal_Doped.write(xyzheader_metal_Doped)
        fid_metal_vacancy.write(xyzheader_metal_vacancy)
        fid_2Doped.write(xyzheader_2Doped)
        fid_1Doped.write(xyzheader_1Doped)
        fid_1vacancy.write(xyzheader_1vacancy)
        fid_2vacancy.write(xyzheader_2vacancy)

        # write xyz coordinates
        fid.write('{:.4f} {:.4f} {:.4f} \n'.format(supercell_a,supercell_b,supercell_c))
        fid_metal_Doped.write('{:.4f} {:.4f} {:.4f} \n'.format(supercell_a,supercell_b,supercell_c))
        fid_metal_vacancy.write('{:.4f} {:.4f} {:.4f} \n'.format(supercell_a,supercell_b,supercell_c))
        fid_2Doped.write('{:.4f} {:.4f} {:.4f} \n'.format(supercell_a,supercell_b,supercell_c))
        fid_1Doped.write('{:.4f} {:.4f} {:.4f} \n'.format(supercell_a,supercell_b,supercell_c))
        fid_1vacancy.write('{:.4f} {:.4f} {:.4f} \n'.format(supercell_a,supercell_b,supercell_c))
        fid_2vacancy.write('{:.4f} {:.4f} {:.4f} \n'.format(supercell_a,supercell_b,supercell_c))


        for i in range(rep_x):
            for j in range(rep_y):
                for k in range(rep_z):
                    prob_metal1 = np.random.rand()
                    metal_vacancy1 = False
                    if prob_metal1<metal_doped_prob:
                        final_atom1 = atom1 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0]) + np.array([metal_dopant_different,0,0,0,0,0])
                        fid_metal_Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom1[0],final_atom1[1],final_atom1[2],final_atom1[3],final_atom1[4],final_atom1[5]))
                    elif metal_doped_prob<prob_metal1<metal_doped_prob+metal_vacancy_prob:
                        final_atom1 = atom1 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        fid_metal_vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom1[0],final_atom1[1],final_atom1[2],final_atom1[3],final_atom1[4],final_atom1[5]))
                        metal_vacancy1 = True
                    else:
                        final_atom1 = atom1 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])

                    prob_metal2 = np.random.rand()
                    metal_vacancy2 = False
                    if prob_metal2<metal_doped_prob:
                        final_atom4 = atom4 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0]) + np.array([metal_dopant_different,0,0,0,0,0])
                        fid_metal_Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom4[0],final_atom4[1],final_atom4[2],final_atom4[3],final_atom4[4],final_atom4[5]))
                    elif metal_doped_prob<prob_metal2<metal_doped_prob+metal_vacancy_prob:
                        fid_metal_vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom4[0],final_atom4[1],final_atom4[2],final_atom4[3],final_atom4[4],final_atom4[5]))
                        metal_vacancy2 = True
                    else:
                        final_atom4 = atom4 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])

                    prob_chalcegon1 = np.random.rand()
                    one_vacancy1, two_vacancy1 = False, False
                    if prob_chalcegon1<dopant_conc_two_subsitutions_high:
                        final_atom2 = atom2+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])+np.array([chalcogen_dopant_different,0,0,0,0,0])
                        fid_2Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom2[0],final_atom2[1],final_atom2[2],final_atom2[3],
                                                                                      final_atom2[4],final_atom2[5]))
                        final_atom3 = atom3+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])+np.array([chalcogen_dopant_different,0,0,0,0,0])
                        fid_2Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom3[0],final_atom3[1],final_atom3[2],final_atom3[3],
                                                                                      final_atom3[4],final_atom3[5]))
                    elif dopant_conc_two_subsitutions_high<prob_chalcegon1<dopant_conc_one_subsitution_high:
                        final_atom2 = atom2+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])+np.array([chalcogen_dopant_different,0,0,0,0,0])
                        fid_1Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom2[0],final_atom2[1],final_atom2[2],final_atom2[3],
                                                                                      final_atom2[4],final_atom2[5]))
                        final_atom3 = atom3+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])

                    elif dopant_conc_one_subsitution_high<prob_chalcegon1<dopant_conc_one_vacancy_high:
                        final_atom2 = atom2+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        fid_1vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom2[0],final_atom2[1],final_atom2[2],final_atom2[3],
                                                                                      final_atom2[4],final_atom2[5]))
                        final_atom3 = atom3+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        one_vacancy1 = True
                    elif dopant_conc_one_vacancy_high<prob_chalcegon1<dopant_conc_two_vacancies_high:
                        final_atom2 = atom2+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        fid_2vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom2[0],final_atom2[1],final_atom2[2],final_atom2[3],
                                                                                      final_atom2[4],final_atom2[5]))
                        final_atom3 = atom3+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        fid_2vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom3[0],final_atom3[1],final_atom3[2],final_atom3[3],
                                                                                      final_atom3[4],final_atom3[5]))
                        two_vacancy1 = True
                    else:
                        final_atom2 = atom2 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        final_atom3 = atom3 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])

                    prob_chalcegon2 = np.random.rand()
                    one_vacancy2, two_vacancy2 = False, False
                    if prob_chalcegon2<dopant_conc_two_subsitutions_high:
                        final_atom5 = atom5+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])+np.array([chalcogen_dopant_different,0,0,0,0,0])
                        fid_2Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom5[0],final_atom5[1],final_atom5[2],final_atom5[3],
                                                                                      final_atom5[4],final_atom5[5]))
                        final_atom6 = atom6+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])+np.array([chalcogen_dopant_different,0,0,0,0,0])
                        fid_2Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom6[0],final_atom6[1],final_atom6[2],final_atom6[3],
                                                                                      final_atom6[4],final_atom6[5]))
                    elif dopant_conc_two_subsitutions_high<prob_chalcegon2<dopant_conc_one_subsitution_high:
                        final_atom5 = atom5+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])+np.array([chalcogen_dopant_different,0,0,0,0,0])
                        fid_1Doped.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom5[0],final_atom5[1],final_atom5[2],final_atom5[3],
                                                                                      final_atom5[4],final_atom5[5]))
                        final_atom6 = atom6+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])

                    elif dopant_conc_one_subsitution_high<prob_chalcegon2<dopant_conc_one_vacancy_high:
                        final_atom5 = atom5+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        fid_1vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom5[0],final_atom5[1],final_atom5[2],final_atom5[3],
                                                                                      final_atom5[4],final_atom5[5]))
                        final_atom6 = atom6+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        one_vacancy2 = True
                    elif dopant_conc_one_vacancy_high<prob_chalcegon2<dopant_conc_two_vacancies_high:
                        final_atom5 = atom5+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        fid_2vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom5[0],final_atom5[1],final_atom5[2],final_atom5[3],
                                                                                      final_atom5[4],final_atom5[5]))
                        final_atom6 = atom6+ i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        fid_2vacancy.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom6[0],final_atom6[1],final_atom6[2],final_atom6[3],
                                                                                      final_atom6[4],final_atom6[5]))
                        two_vacancy2 = True
                    else:
                        final_atom5 = atom5 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])
                        final_atom6 = atom6 + i*np.array([0,a,0,0,0,0])+j*np.array([0,0,b,0,0,0])+k*np.array([0,0,0,c,0,0])

                    if metal_vacancy1 == False:
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom1[0],final_atom1[1],final_atom1[2],final_atom1[3],
                                                                                      final_atom1[4],final_atom1[5]))



                    if one_vacancy1 == False and two_vacancy1 == False:
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom2[0],final_atom2[1],final_atom2[2],final_atom2[3],
                                                                                          final_atom2[4],final_atom2[5]))
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom3[0],final_atom3[1],final_atom3[2],final_atom3[3],
                                                                                          final_atom3[4],final_atom3[5]))
                    elif one_vacancy1 == True:
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom2[0],final_atom2[1],final_atom2[2],final_atom2[3],
                                                                                          final_atom2[4],final_atom2[5]))

                    elif two_vacancy1 == True:
                        pass

                    if metal_vacancy2 == False:
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom4[0],final_atom4[1],final_atom4[2],final_atom4[3],
                                                                                      final_atom4[4],final_atom4[5]))

                    if one_vacancy2 == False and two_vacancy2 == False:
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom5[0],final_atom5[1],final_atom5[2],final_atom5[3],
                                                                                          final_atom5[4],final_atom5[5]))
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom6[0],final_atom6[1],final_atom6[2],final_atom6[3],
                                                                                          final_atom6[4],final_atom6[5]))
                    elif one_vacancy1 == True:
                        fid.write('{:.0f} {:.6f} {:.6f} {:.6f} {:.0f} {:.2f} \n'.format(final_atom5[0],final_atom5[1],final_atom5[2],final_atom5[3],
                                                                                          final_atom5[4],final_atom6[5]))
                    else:
                        pass

        #-1 as End of the file and close all xyz files
        fid.write('-1')
        fid_metal_Doped.write('-1')
        fid_metal_vacancy.write('-1')
        fid_2Doped.write('-1')
        fid_1Doped.write('-1')
        fid_1vacancy.write('-1')
        fid_2vacancy.write('-1')

        # close the files
        fid.close()
        fid_metal_Doped.close()
        fid_metal_vacancy.close()
        fid_2Doped.close()
        fid_1Doped.close()
        fid_1vacancy.close()
        fid_2vacancy.close()

        #Generate parameter files
        fid_param = open(filename+'.param','w+')
        fid_metal_Doped_param = open(filename_metal_Doped+'.param','w+')
        fid_metal_vacancy_param = open(filename_metal_vacancy+'.param','w+')
        fid_2Doped_param = open(filename_2Doped+'.param','w+')
        fid_1Doped_param = open(filename_1Doped+'.param','w+')
        fid_1vacancy_param = open(filename_1vacancy+'.param','w+')
        fid_2vacancy_param = open(filename_2vacancy+'.param','w+')
        #Wirte the xyz filename as the header
        fid_param.write(filename+'.xyz \n1 1 1\n')
        fid_metal_Doped_param.write(filename_metal_Doped+'.xyz \n1 1 1\n')
        fid_metal_vacancy_param.write(filename_metal_vacancy+'.xyz \n1 1 1\n')
        fid_2Doped_param.write(filename_2Doped+'.xyz \n1 1 1\n')
        fid_1Doped_param.write(filename_1Doped+'.xyz \n1 1 1\n')
        fid_1vacancy_param.write(filename_1vacancy+'.xyz \n1 1 1\n')
        fid_2vacancy_param.write(filename_2vacancy+'.xyz \n1 1 1\n')
        #Save the image as .tif files
        fid_param.write('Image'+filename+'.tif')
        fid_metal_Doped_param.write('metal_Doped_'+filename+'.tif')
        fid_metal_vacancy_param.write('metal_vacancy_'+filename+'.tif')
        fid_2Doped_param.write('2Doped_'+filename+'.tif')
        fid_1Doped_param.write('1Doped_'+filename+'.tif')
        fid_1vacancy_param.write('1vacancy_'+filename+'.tif')
        fid_2vacancy_param.write('2vacancy_'+filename+'.tif')
        #Set parameters

        #get different str lines in .param files
        #Image and defect map have the same paramter file
        # except conting noise
        image_size_str = '\n'+str(image_size)+' '+str(image_size)
        Cs3 = np.random.normal(Cs3_param[0],Cs3_param[1])
        Cs5 = np.random.normal(Cs5_param[0],Cs5_param[1])
        STEM_Param_str = '\n'+str(voltage)+' '+str(Cs3)+' '+str(Cs5)+' '+str(df)+' '+str(aperture)
        ADF_str = '\n' + str(ADF_angle_min) + ' '+ str(ADF_angle_max)
        High_order_str = '\n'+Higher_order
        Source_size = np.random.normal(Source_size_param[0],Source_size_param[1])
        source_size_str = '\n'+str(Source_size)
        defocus_spread = np.random.normal(defocus_spread_param[0],defocus_spread_param[1])
        Defocus_str = '\n'+str(defocus_spread)
        probe_current = np.random.normal(probe_current_param[0],probe_current_param[1])
        noise_str = '\n'+str(probe_current)+' '+str(dwell_time)

        GeneralParam = image_size_str+STEM_Param_str+ADF_str+High_order_str+source_size_str+Defocus_str+'\ny'+noise_str+'\n-1'
        GeneralParam_defect = image_size_str+STEM_Param_str+ADF_str+High_order_str+source_size_str+Defocus_str+'\nn'+'\n-1'

        fid_param.write(GeneralParam)
        fid_metal_Doped_param.write(GeneralParam_defect)
        fid_metal_vacancy_param.write(GeneralParam_defect)
        fid_2Doped_param.write(GeneralParam_defect)
        fid_1Doped_param.write(GeneralParam_defect)
        fid_1vacancy_param.write(GeneralParam_defect)
        fid_2vacancy_param.write(GeneralParam_defect)
        #close param files
        fid_param.close()
        fid_metal_Doped_param.close()
        fid_metal_vacancy.close()
        fid_2Doped_param.close()
        fid_1Doped_param.close()
        fid_1vacancy_param.close()
        fid_2vacancy_param.close()
    #close batch files
    file_batch.close()

