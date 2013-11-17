# Dataset: '2011-06-06_DV_8x1_Dicr_xFRET10b_100pM_1'

ph2 = d.ph_times_m[1]
ph2_donor = ph2[~d.red[1]]
ph2_d_u32 = ph2_donor.astype('uint32')
fi = open(name+'_ch2_u32_TS.dat', 'w')
fi.write(ph2_d_u32.tostring()) # write in little-endian uint32 format
fi.close()
