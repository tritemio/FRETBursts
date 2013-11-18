"""
Process a batch of 4-spot FRET measurements. Obsolete.
"""


#from subprocess import Popen
from glob import glob

fname_list = [
#'Dualview_8x1_BS_200nm_100X_36uW_0.dat',
#'Dualview_8x1_BS_SSDNA_100pM_30mW_2.dat',
#'Dualview_8x1_Cy3Cy5Dicr_SSDNA_500pM_80mW_3.dat',
#'Dualview_8x1_Cy3Cy5Dicr_SSDNA_100pM_80mW_1.dat',
#'DV_8x1_Dicr_DSDNA_tmr_atto_100pM_80mW_2.dat',
#'2011-05-04/DV_8x1_Dicr_dsdna_100pM_160mW_12.dat',
#'2011-05-04/DV_8x1_Dicr_dsdna_100pM_240mW_14.dat',
#'2011-05-06/DV_8x1_Dicr_dsdna_10pM_75mW_2.dat',
#'2011-05-16/DV_8x1_Dicr_dsdna_10pM_1.dat',
#'2011-05-17/DV_8x1_Dicr_ssdna_q100pM_3.dat',
#'2011-05-18/DV_8x1_Dicr_ssdna_q100pM_sample2_4.dat',
#'2011-05-18/DV_8x1_Dicr_ssdna_tmr_alexa_100pM_2.dat',
#'2011-05-18/DV_8x1_Dicr_high_fret_no1_100pM_1.dat'
#'2011-05-18/DV_8x1_Dicr_high_fret_no1_10pM_2.dat'
'2011-05-27/DV_8x1_Dicr_xFRET3b_100pM_0.dat',
'2011-05-27/DV_8x1_Dicr_xFRET7b_100pM_1.dat',
'2011-06-06/DV_8x1_Dicr_xFRET5b_100pM_0.dat',
'2011-06-06/DV_8x1_Dicr_xFRET10b_100pM_1.dat',
'2011-06-07/DV_8x1_Dicr_xFRET10b_33pM_0.dat',
'2011-06-07/DV_8x1_Dicr_xFRET3b_33pM_3.dat',
'2011-06-07/DV_8x1_Dicr_xFRET5b_33pM_2.dat',
'2011-06-07/DV_8x1_Dicr_xFRET7b_33pM_1.dat'
]
BT=0.11

fname_list = glob('data/2011-06-06/*.dat'); BT=0.11

for fname in fname_list:
    print("\n >>>>>>>>>>>>>>> DATA: "+fname+"\n\n")
    #p = Popen("ipython" + " -pylab burst.py", shell=True)
    #sts = os.waitpid(p.pid, 0)[1]
    ip = get_ipython()
    ip.magic("run -i burst.py "+fname)
    close('all')
    del d
    del ph_times
    del ph_times_m
    del mburst

