#/usr/bin/env python
#*****************************************************************************************************************************
#Program to construct PFT map from 2010 to 850 for ORCHIDEE forcing (CMIP6 runs) using Hurtt's data set
#*****************************************************************************************************************************
#Note, we consider rangeland as natural vegetation
#import modules
#from __future__ import division
import sys
import cdms2, MV2, numpy as np
import ctypes
import cdtime, cdutil
#import  vcs
import numpy.ma as ma
#***************************************************************************************
#***************************************************************************************
#inialize vcs canvas to plot the variables for debugging purpose
#bg = not False
#x=vcs.init()
#x.setcolormap("rainbow")
#gm = vcs.createboxfill()

cdms2.setNetcdfShuffleFlag(1) ## where value is either 0 or 1, 0 for compressed NetCDf file writing
cdms2.setNetcdfDeflateFlag(1) ## where value is either 0 or 1
cdms2.setNetcdfDeflateLevelFlag(4) ## where value is a integer between 0 and 9 included

#***************************************************************************************
#***************************************************************************************
#  C3/C4 separation we use Koppen Geiger (KG) climate classification and follow Table 3 in Poulter et al. 2011 paper that you already to split the ~30 KG zones into 5 separate zones.
kg_file='/home/basc/dnarayan/LandUsemap_ORCHIDEE/Natasha/koppen_files/koppen_fill_025deg_mode-from-addo.dat'
kg_biome = np.loadtxt(kg_file, dtype=int)
pft_biome = np.zeros_like(kg_biome)
pft_biome = np.where(kg_biome == -9999, kg_biome, pft_biome)
print 'pft_biome',pft_biome
# Tropical
pft_biome = np.where( ( ((kg_biome >= 0) & (kg_biome <= 4)) | (kg_biome == 6) ), 1, pft_biome )
# Temperate (warm)
pft_biome = np.where( ((kg_biome == 5) | (kg_biome == 7) | (kg_biome == 8) | (kg_biome == 11) | (kg_biome == 14)), 2, pft_biome )
# Temperate (cold)
pft_biome = np.where( ((kg_biome == 9) | (kg_biome == 10) | (kg_biome == 12) | (kg_biome == 13) | (kg_biome == 15) | (kg_biome == 16)), 3, pft_biome )
# Boreal (warm)
pft_biome = np.where( ((kg_biome == 17) | (kg_biome == 21) | (kg_biome == 25)), 4, pft_biome )
# Boreal (cold)
pft_biome = np.where( ( ((kg_biome >= 18) & (kg_biome <= 20)) | ((kg_biome >= 22) & (kg_biome <= 24)) | (kg_biome >= 26) ), 5, pft_biome ) 
#C3 climatic zones are pft_biome =3,4,5 and c4 climatic zones are pft_biome=1,2
pft_biome = MV2.array((pft_biome))
pft_biome = MV2.masked_equal(pft_biome,-9999)
#**************************************************************************************************
#**************************************************************************************************

# Parameters
# Interpolation years
year_start = 2010
year_end = 1850

nyears = (year_start - year_end) + 1 # 1 is added to include end  year
# Latitude and longitude range
latrange = (90.,-90.)
lonrange = (-180, 180, 'co')
#select time period
time	= ('1850-01-01', '2010-01-01')
time_trans = ('1849-01-01', '2009-01-01') #time for transitions
print '>>>>>>>>>>>>> Number of years pft map to be created now <<<<<<<<<<<<<<<', nyears

#*************************************************************************************************
#*************************************************************************************************
# Input / Output files directory

esa_pft_dir = '/home/basc/dnarayan/LandUsemap_ORCHIDEE/Natasha/ctl2_13July2016/'

hurt_states_dir = '/home/basc/dnarayan/LandUsemap_ORCHIDEE/LUH2v2h/'

output_dir = '/home/basc/dnarayan/LandUsemap_ORCHIDEE/Program_LUmap_Dev/Program_Final_ESA_LUH2v2h/'

esa_pft_file = esa_pft_dir+'ESACCI-LC-PFT_biome_map_aggreg-0.25deg_2010_v1.6.1_upgraded_CW_legend_TUWienKG.nc' #map from natasha 

hurt_states_file = hurt_states_dir+'states1850_2015deflate.nc'
# Open specified Files
#ESA map file
esa_file 	= cdms2.open(esa_pft_file) 

#Hurtts states file
hurt_file 	= cdms2.open(hurt_states_file)

#Open transition files from Hurtt Data LUHv0.4_historical
transitionf 	= cdms2.open(hurt_states_dir+'transitions1850_2015_deflate.nc')
staticf 	= cdms2.open(hurt_states_dir+'staticData_quarterdeg.nc')

# output ORCHIDEE PFTmap (ORCmap)

pft_file_hist19 = output_dir +'ORCHIDEE_19PFTmap_1850to2010_cmpi6_LUH2v2h.nc'
pft_file_hist13 = output_dir +'ORCHIDEE_13PFTmap_1850to2010_cmpi6_LUH2v2h.nc'

#**************************************************************************************************
#**************************************************************************************************
# Specify Required Variables
#------------------------------
#ESA PFT map variable 'maxvegetfrac'
esa_var = 'maxvegetfrac'

# states variable names From Hurtt Data set
primf		= 'primf'
secdf		= 'secdf'
primn		= 'primn'
secdn		= 'secdn'
c3ann		= 'c3ann'
c4ann		= 'c4ann'
c3per		= 'c3per'
c4per		= 'c4per'
c3nfx		= 'c3nfx'
pastr		= 'pastr'          
rangel		= 'range'	
urban		= 'urban'
icwtr		= 'icwtr'

#store variables 
primf	= hurt_file(primf,time=time,latitude=latrange, longitude=lonrange)
secdf	= hurt_file(secdf,time=time,latitude=latrange, longitude=lonrange)
primn	= hurt_file(primn,time=time,latitude=latrange, longitude=lonrange)
secdn	= hurt_file(secdn,time=time,latitude=latrange, longitude=lonrange)
c3ann	= hurt_file(c3ann,time=time,latitude=latrange, longitude=lonrange)
c4ann	= hurt_file(c4ann,time=time,latitude=latrange, longitude=lonrange)
c3per	= hurt_file(c3per,time=time,latitude=latrange, longitude=lonrange)
c4per	= hurt_file(c4per,time=time,latitude=latrange, longitude=lonrange)
print 'c4per',c4per.shape
c3nfx	= hurt_file(c3nfx, time=time,latitude=latrange, longitude=lonrange)
pastr	= hurt_file(pastr,time=time,latitude=latrange, longitude=lonrange)
rangel	= hurt_file(rangel,time=time,latitude=latrange, longitude=lonrange) 
print 'range land',rangel.shape
urban	= hurt_file(urban,time=time,latitude=latrange, longitude=lonrange)  
print 'urban',urban.shape
icwtr	= staticf(icwtr,latitude=latrange, longitude=lonrange)  #constant over time

#Store Transition variables

primf_to_c3ann = transitionf('primf_to_c3ann',time=time_trans,latitude=latrange, longitude=lonrange)
primf_to_c4ann = transitionf('primf_to_c4ann',time=time_trans,latitude=latrange, longitude=lonrange)
primf_to_c3per = transitionf('primf_to_c3per',time=time_trans,latitude=latrange, longitude=lonrange)
primf_to_c4per = transitionf('primf_to_c4per',time=time_trans,latitude=latrange, longitude=lonrange)
primf_to_c3nfx = transitionf('primf_to_c3nfx',time=time_trans,latitude=latrange, longitude=lonrange)
primf_to_pastr = transitionf('primf_to_pastr',time=time_trans,latitude=latrange, longitude=lonrange)
primf_to_urban = transitionf('primf_to_urban',time=time_trans,latitude=latrange, longitude=lonrange)
primf_to_range = transitionf('primf_to_range',time=time_trans,latitude=latrange, longitude=lonrange)

secdf_to_c3ann = transitionf('secdf_to_c3ann',time=time_trans,latitude=latrange, longitude=lonrange)
secdf_to_c4ann = transitionf('secdf_to_c4ann',time=time_trans,latitude=latrange, longitude=lonrange)
secdf_to_c3per = transitionf('secdf_to_c3per',time=time_trans,latitude=latrange, longitude=lonrange)
secdf_to_c4per = transitionf('secdf_to_c4per',time=time_trans,latitude=latrange, longitude=lonrange)
secdf_to_c3nfx = transitionf('secdf_to_c3nfx',time=time_trans,latitude=latrange, longitude=lonrange)
secdf_to_pastr = transitionf('secdf_to_pastr',time=time_trans,latitude=latrange, longitude=lonrange)
secdf_to_urban = transitionf('secdf_to_urban',time=time_trans,latitude=latrange, longitude=lonrange)
secdf_to_range = transitionf('secdf_to_range',time=time_trans,latitude=latrange, longitude=lonrange)
print 'secdf to all reading is done',secdf_to_range.shape
primn_to_c3ann = transitionf('primn_to_c3ann',time=time_trans,latitude=latrange, longitude=lonrange)
primn_to_c4ann = transitionf('primn_to_c4ann',time=time_trans,latitude=latrange, longitude=lonrange)
primn_to_c3per = transitionf('primn_to_c3per',time=time_trans,latitude=latrange, longitude=lonrange)
primn_to_c4per = transitionf('primn_to_c4per',time=time_trans,latitude=latrange, longitude=lonrange)
primn_to_c3nfx = transitionf('primn_to_c3nfx',time=time_trans,latitude=latrange, longitude=lonrange)
primn_to_pastr = transitionf('primn_to_pastr',time=time_trans,latitude=latrange, longitude=lonrange)
primn_to_urban = transitionf('primn_to_urban',time=time_trans,latitude=latrange, longitude=lonrange)
primn_to_range = transitionf('primn_to_range',time=time_trans,latitude=latrange, longitude=lonrange)
print 'prim non forest to all reading is done', primn_to_pastr.shape

secdn_to_c3ann = transitionf('secdn_to_c3ann',time=time_trans,latitude=latrange, longitude=lonrange)
secdn_to_c4ann = transitionf('secdn_to_c4ann',time=time_trans,latitude=latrange, longitude=lonrange)
secdn_to_c3per = transitionf('secdn_to_c3per',time=time_trans,latitude=latrange, longitude=lonrange)
secdn_to_c4per = transitionf('secdn_to_c4per',time=time_trans,latitude=latrange, longitude=lonrange)
secdn_to_c3nfx = transitionf('secdn_to_c3nfx',time=time_trans,latitude=latrange, longitude=lonrange)
secdn_to_pastr = transitionf('secdn_to_pastr',time=time_trans,latitude=latrange, longitude=lonrange)
secdn_to_urban = transitionf('secdn_to_urban',time=time_trans,latitude=latrange, longitude=lonrange)
secdn_to_range = transitionf('secdn_to_range',time=time_trans,latitude=latrange, longitude=lonrange)
print 'secd non forest to all reading is done', secdn_to_pastr.shape


c3ann_to_secdf = transitionf('c3ann_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
c3ann_to_secdn = transitionf('c3ann_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)

c4ann_to_secdf = transitionf('c4ann_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
c4ann_to_secdn = transitionf('c4ann_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)

c3per_to_secdf = transitionf('c3per_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
c3per_to_secdn = transitionf('c3per_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)

c4per_to_secdf = transitionf('c4per_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
c4per_to_secdn = transitionf('c4per_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)

c3nfx_to_secdf = transitionf('c3nfx_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
c3nfx_to_secdn = transitionf('c3nfx_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)

pastr_to_secdf = transitionf('pastr_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
pastr_to_secdn = transitionf('pastr_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)

urban_to_secdf = transitionf('urban_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
urban_to_secdn = transitionf('urban_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)

range_to_secdf = transitionf('range_to_secdf',time=time_trans,latitude=latrange, longitude=lonrange)
range_to_secdn = transitionf('range_to_secdn',time=time_trans,latitude=latrange, longitude=lonrange)
print 'All transition variables are read', urban_to_secdf.shape

#store PFT variable from present day map
pft		= esa_file(esa_var,latitude=latrange, longitude=lonrange)
print 'Present day() pft map read',pft.shape
pft		= pft(squeeze=1)
print 'pft squeeze',pft.shape
# Variables reading is done here

#**************************************************************************************************
#**************************************************************************************************
#==================================================================================================
#Get the axis names for further stamping to new map
set_axes = pft.getAxisList(('lat','lon'))
lat	= pft.getLatitude()
lon	= pft.getLongitude()
nlat	= len(lat)
nlon    = len(lon)
print 'length of latitude', nlat
nlon	= len(lon)
print 'length of longitude', nlon
pft		= np.nan_to_num(pft)
pft[pft==1.e+20]=0
#=====================================================================================================
#Create a variable 'pft_back' to store the final PFT values after calculations (nyears,levels,lat,lon)
#for nyears, 19 veg types,0.25 deg lat lon resolution
pft_back		= np.zeros((nyears,19,nlat,nlon)) 
print 'pft_back', pft_back.shape
#create Time axis
time = cdms2.createAxis(np.arange(year_end,year_start+1,1))
time.designateTime(calendar = cdtime.MixedCalendar)
time.id = 'time_counter'
time.units = 'years since 0-1-1'
print 'time',time.shape
#create Level axis

lev = cdms2.createAxis(range(1,20))    #PFT numbers 1 to 19 
lev.designateLevel()
lev.id = 'veget'
lev.units ='-'

#lat lon axis values are taken directly from ESA map
pft_back = cdms2.createVariable(pft_back,axes=[time]+[lev]+set_axes,id='maxvegetfrac')
print '--------------------------------------------------->>>>>>   4d variable',pft_back.shape
pft_back.long_name ='Vegetation types'
pft_back.units = '-'
pft19_axislist_pft = pft_back.getAxisList(('veget', 'lat','lon'))
#==========================================================================================
#for nyears, 13 veg types,0.25 deg lat lon resolution
pft_13		= np.zeros((nyears,13,nlat,nlon)) 
print 'pft_13', pft_13.shape
#create Time axis
time = cdms2.createAxis(np.arange(year_end,year_start+1,1))
time.designateTime(calendar = cdtime.MixedCalendar)
time.id = 'time_counter'
time.units = 'years since 0-1-1'
print 'time',time.shape
#create Level axis

lev_13 = cdms2.createAxis(range(1,14))    #PFT numbers 1 to 13
lev_13.designateLevel()
lev_13.id = 'veget'
lev_13.units ='-'

#lat lon axis values are taken directly from ESA map
pft_13 = cdms2.createVariable(pft_13,axes=[time]+[lev_13]+set_axes,id='maxvegetfrac')
print '--------------------------------------------------->>>>>>   4d variable',pft_13.shape
pft_13.long_name ='Vegetation types'
pft_13.units = '-'

pft13_axislist_pft = pft_13.getAxisList(('veget', 'lat','lon'))
timax_slide = cdms2.createAxis(np.asarray([2010], np.float32))
timax_slide.designateTime(calendar = cdtime.MixedCalendar)
timax_slide.id = 'time_counter'
timax_slide.units = 'years since 0-1-1'
#********************************************************************************************************
#********************************************************************************************************
#Get the LUH masking
#mask_file	= cdms2.open('esa_mask.nc')
#esa_mask 	= mask_file('Mask_only',latitude=latrange, longitude=lonrange)
luh_mask 	= MV2.getmask(primf[0,:,:])
print 'luh_mask', luh_mask.shape
icwtr		= MV2.masked_array(icwtr,mask=luh_mask)
#Total land frac in ORCmap  (Sum of all land units in ORCmap is 1)
orc_pft = MV2.masked_array(pft)
for n1 in range(0,15):
	orc_pft[n1,:,:] = MV2.masked_array(orc_pft[n1,:,:],mask=luh_mask)
	#print n1
lndf_esa = MV2.minimum(MV2.sum(orc_pft,axis=0),1)
lndf_esa.id = 'lndfrac_esa'
lndf_esa.long_name = 'Sum of all ORCmap PFTs'


#Total land frac in Hurtt data set (sum of all land units in LUH is one after including icwtr)
lndf_hurt = primf[nyears-1,:,:] + secdf[nyears-1,:,:] + primn[nyears-1,:,:] + secdn[nyears-1,:,:] + c3ann[nyears-1,:,:] + c4ann[nyears-1,:,:] + c3per[nyears-1,:,:] + c4per[nyears-1,:,:] + c3nfx[nyears-1,:,:] + pastr[nyears-1,:,:]+rangel[nyears-1,:,:]+urban[nyears-1,:,:]+icwtr
lndf_hurt.id = 'landfrac_hurtt'
lndf_hurt.long_name = 'Sum of all land units in Hurtt Data 2010'
lndf_hurt = np.around(lndf_hurt,decimals=5)
lndf_hurt = MV2.masked_array(lndf_hurt,mask=luh_mask)
#x.plot(lndf_hurt[::-1,:],gm,bg=bg)
#x.png('lndf_hurt.png')
#x.clear()
#==================================================================================================
#Present day croplands from LUH data
c3ann2010 	= c3ann[nyears-1,:,:]
c4ann2010 	= c4ann[nyears-1,:,:]
c3per2010 	= c3per[nyears-1,:,:]
c4per2010 	= c4per[nyears-1,:,:]
c3nfx2010 	= c3nfx[nyears-1,:,:]
pastr2010	= pastr[nyears-1,:,:]
rangel2010  	= rangel[nyears-1,:,:]
urban2010	= urban[nyears-1,:,:]


#**************************************************************************************************
#**************************************************************************************************
#Impose Hurtt's present day anthropogenic area onto ORCmap and reduce/expand natural vegetation
#=============================================================================================
print 'Present day map creation starts here',year_start

#Baresoil in present ORCmap
bsoil2010	= MV2.minimum(orc_pft[0,:,:],1.)+np.where(lndf_esa==0,1.0,0)
#x.plot(bsoil2010[::-1,:],gm,bg=bg)
#x.png('bsoil2010.png')
#x.clear()
#Natural vegetation types in present day map
nveg2010 	= orc_pft[1:11,:,:]
#Total Natural vegetation in present day ORCmap; 
totnveg2010	= MV2.minimum(MV2.sum(nveg2010,axis=0),1.0)
# Save natural vegetation to check
totnveg2010.id = 'natveg'
totnveg2010.units = '-'

#total anthropogenic area from Hurtt's data in 2010
totanth_hurt2010= c3ann2010 + c4ann2010 + c3nfx2010 + c3per2010 + c4per2010 + pastr2010 + urban2010 + rangel2010

#Convert managed pasture as anthropogenic c3/c4 grass using koppen climate conditions
anthgr2010	= pastr2010
c3_anthgr2010 = orc_pft[11,:,:]*0.
c4_anthgr2010 = orc_pft[12,:,:]*0.
totcrop2010 = orc_pft[11,:,:] + orc_pft[12,:,:]
#c3angrf2010 = MV2.array(MV2.where(totcrop2010==0.,0.,orc_pft[11,:,:]/totcrop2010),dtype='f')
#c4angrf2010 = MV2.array(MV2.where(totcrop2010==0.,0.,orc_pft[12,:,:]/totcrop2010),dtype='f')
#print 'anthgr',anthgr2010.shape
c3_anthgr2010 = MV2.where(((anthgr2010!=0) & (pft_biome==3) | (pft_biome==4)|(pft_biome==5)),anthgr2010,0)
c4_anthgr2010 = MV2.where(((anthgr2010!=0) & (pft_biome==1) | (pft_biome==2)),anthgr2010,0)
#print 'c3_anthgrass',c3_anthgr2010.shape
#print 'c4_anthgrass',c4_anthgr2010.shape
print 'Anthropogenic grass partitioning is done here'
totanthgr2010 = c3_anthgr2010 + c4_anthgr2010
#c3/c4 natural grassland in ORCmap
c3gr2010	= orc_pft[9,:,:]
c4gr2010	= orc_pft[10,:,:]
totgr2010   = c3gr2010 + c4gr2010
#present day ratio of c3/c4 grass
c3grf2010	= MV2.where(totgr2010==0.,0.,c3gr2010/totgr2010)
c4grf2010	= MV2.where(totgr2010==0.,0.,c4gr2010/totgr2010)

#Derive natural vegetation after imposing anthropgenic area
totnveg2010_new = 1-(bsoil2010+icwtr)-totanth_hurt2010

#if there is not enough natural vegetation for anthropogenic area then reduce baresoil  
excess_anth		= MV2.where(totnveg2010_new>0.,0.,totnveg2010_new)
re_bsoil		= bsoil2010 + excess_anth
bsoil_now		= MV2.where(re_bsoil<=0.,0.,re_bsoil)
#x.plot(bsoil_now[::-1,:],gm,bg=bg)
#x.png('bsoil_now.png')
#x.clear()

#Update the new natural vegetation
totnveg2010_new 	= 1-(bsoil_now+icwtr)-totanth_hurt2010
totnveg2010_new		= MV2.where(totnveg2010_new<=0.,0,totnveg2010_new)

#Store all PFTs into one vaariable 'pft2010'
pft2010			= pft_back[nyears-1,:,:,:]
#baresoil
pft2010[0,:,:]  	= bsoil_now
#Crop PFTs impose
pft2010[11,:,:] 	= c3ann2010  
pft2010[12,:,:] 	= c4ann2010    
pft2010[13,:,:] 	= c3_anthgr2010
pft2010[14,:,:] 	= c4_anthgr2010	
pft2010[15,:,:] 	= c3per2010
pft2010[16,:,:] 	= c4per2010
pft2010[17,:,:] 	= c3nfx2010 
pft2010[18,:,:] 	= urban2010

#Adjust Natural vegetation i.e.expand/reduce
#before adjusting make sure not to overlap vegetation on baresoil
totveg2010 	 	= totnveg2010_new  + totanth_hurt2010
totveg2010_notOK	= MV2.masked_values(totveg2010-lndf_hurt,0.).mask
no_nvegmask		= MV2.masked_values(totveg2010,0.).mask
for npft in range(10):
	print 'npft',npft
	pft_n		= nveg2010[npft,:,:]
	pft_tmp= MV2.choose(totveg2010_notOK,(MV2.choose(no_nvegmask,(pft_n*totnveg2010_new/MV2.where(totnveg2010==0,1,totnveg2010),0.)),pft_n*(totnveg2010_new/MV2.where(totnveg2010==0,1,totnveg2010))))     
	print 'Reduction/ Expansion of Natural vegetation is done'
	pft2010[npft+1,:,:] = pft_tmp
	print 'pft_back assignment is done'

#Since we consider rangeland as Natural grassland partition; add here to natural grassland after partitioning c3/c4
c3_range2010 = c3gr2010*0.0
c4_range2010 = c4gr2010*0.0
c3_range2010 = MV2.where(((rangel2010!=0) & (pft_biome==3) | (pft_biome==4)|(pft_biome==5)),rangel2010,0)
c4_range2010 = MV2.where(((rangel2010!=0) & (pft_biome==1) | (pft_biome==2)),rangel2010,0)
pft2010[9,:,:] = pft2010[9,:,:] + c3_range2010
pft2010[10,:,:] = pft2010[10,:,:] + c4_range2010

pft_sum2010 = MV2.sum(pft2010, axis=0)+icwtr
print 'Sum of all PFT is',pft_sum2010.shape
#x.plot(pft_sum2010[::-1,:],gm,bg=bg)
#x.png('pft_sum2010.png')
#x.clear()

#In case if sum of all PFTs is greater than 1, reduce all PFTs
re_pft2010 = MV2.masked_greater(pft_sum2010 - lndf_hurt, 0.).mask
for npft in range(19):
	#print 'npft',npft
	pft_tmp1 = pft2010[npft,:,:]
	pft_tmp1 = MV2.choose(re_pft2010, (pft_tmp1, pft_tmp1*lndf_hurt/pft_sum2010))
	pft2010[npft,:,:] = pft_tmp1

#In case if sum of all PFTs is less than 1, increase all PFTs
inc_pft2010 = MV2.masked_less(pft_sum2010 - lndf_hurt, 0.).mask
for npft in range(19):
	#print 'npft',npft
	pft_tmp2 = pft2010[npft,:,:]
	pft_tmp2 = MV2.choose(inc_pft2010, (pft_tmp2, pft_tmp2*lndf_hurt/pft_sum2010))
	pft2010[npft,:,:] = pft_tmp2

newpftsum = MV2.sum(pft2010, axis=0)+icwtr
sys.stdout.flush()
pft_back[[nyears-1],:,:,:] = pft2010
#x.plot(newpftsum[::-1,:],gm,bg=bg)
#x.png('newpftsum.png')
#x.clear()
#if (pft_sum2010.all()==1):
#   print 'Sum of all PFTs is equal to 1'
#else:
#   print 'Sum of all PFTs is not equal to 1'
#Store the present day map into file
timax_slide[0] = 2010
g01 = cdms2.open('19PFTmap_2010_cmpi6_LUH2v2h.nc','w')
g01.write(pft_back[nyears-1,:,:,:],axes = [timax_slide] + pft19_axislist_pft)
g01.close()

#==========================================================================================================================
#Merge 19 PFTs into 13 PFTs for ORCHIDEE standard version model
#19 PFT names: 0-Baresoil, 1-TropBEf, 2-TropBRf,3-TempNEf, 4-TempBEf,5-TempBSf
#6-BorNEf, 7-BorBSf, 8-BorNSf, 9-C3Natural grass, 10-C4Natural grass
#11-C3 annual crops, 12-C4 annual crops, 13-C3 anthropogenic grass
#14-C4 anthropogenic grass,15-c3 perennial crops, 16-C4 perennial crops
#17-c3 nitrogen fixing crop, 18-urban
#Pasture and urban land is treated as Natural grassland
#c3 and c4 perennial crops are treated as c3/c4 annual crops

pft_13[[nyears-1],0,:,:] = pft2010[0,:,:]
pft_13[[nyears-1],1:9,:,:] = pft2010[1:9,:,:]
#we add urban to natural grassland
c3u2010 = totveg2010*0.0
c4u2010 = totveg2010*0.0
c3u2010 = MV2.where(((pft2010[18,:,:]!=0) & (pft_biome==3) | (pft_biome==4)|(pft_biome==5)),pft2010[18,:,:],0)
c4u2010 = MV2.where(((pft2010[18,:,:]!=0) & (pft_biome==1) | (pft_biome==2)),pft2010[18,:,:],0)
pft_13[[nyears-1],9,:,:] = pft2010[9,:,:]+c3u2010+pft2010[13,:,:]
pft_13[[nyears-1],10,:,:] = pft2010[10,:,:]+c4u2010+pft2010[14,:,:]
pft_13[[nyears-1],11,:,:] = pft2010[11,:,:]+pft2010[15,:,:]+pft2010[17,:,:]
pft_13[[nyears-1],12,:,:] = pft2010[12,:,:]+pft2010[16,:,:]
g11 = cdms2.open('13PFTmap_2010_cmpi6_LUH2v2h.nc','w')
g11.write(pft_13[nyears-1,:,:,:],axes = [timax_slide] + pft13_axislist_pft)
g11.close()
print '>>>>>>>>>>>>>>>>>Present day PFT map creation is done here'

#**************************************************************************************************
#**************************************************************************************************
print '>>>>>>>>>>>>>>>>>Start Backward reconstruction<<<<<<<<<<<<<<<<<<<<<<<<<'

#ratio of tree PFTs in present day map
tree2010    = pft2010[1:9,:,:]
tottree2010 = MV2.sum(tree2010,axis=0)
treef_2010  = tree2010/tottree2010
treef_2010  = np.nan_to_num(treef_2010)
print 'treef_2010',treef_2010.shape
mask_notreef2010 = MV2.masked_equal(tree2010/tottree2010,0.)
#x.plot(mask_notreef2010[6,:,:],gm,bg=bg)
#x.png('mask_notreef2010_6.png')
#x.clear()

#Apply year to year transitions
#=============================================================================
for i in range(nyears-1):
	print 'nyears', nyears
	print 'Calculating for the years',2010-i, 2010-(i+1)
	nyr = nyears-1
	sys.stdout.flush()
	bsoil_p   = pft_back[nyr-i,0,:,:]
	tree_p    = pft_back[nyr-i,1:9,:,:]
	#print 'Tree_p',tree_p.shape
	tottree_p = MV2.sum(tree_p,axis=0)
	c3gr_p 	  = pft_back[nyr-i,9,:,:]
	c4gr_p    = pft_back[nyr-i,10,:,:]
	totgr_p   = c3gr_p + c4gr_p
	
	#Read LUH2 anthropogenic area in the year before present
	j = i+1
	c3ann_b = c3ann[nyr-j,:,:] 
	c4ann_b = c4ann[nyr-j,:,:] 
	c3per_b = c3per[nyr-j,:,:] 
	c4per_b = c4per[nyr-j,:,:] 
	c3nfx_b = c3nfx[nyr-j,:,:]
	pastr_b = pastr[nyr-j,:,:]
	urban_b = urban[nyr-j,:,:]
	rangel_b = rangel[nyr-j,:,:]
	#c3/c4 anthropogenic grass partitioning from pasture
	c3angr_b = c3_anthgr2010*0.0
	c4angr_b = c4_anthgr2010*0.0
	c3angr_b = MV2.where(((pastr_b!=0) & (pft_biome==3) | (pft_biome==4)|(pft_biome==5)),pastr_b,0)
	c4angr_b = MV2.where(((pastr_b!=0) & (pft_biome==1) | (pft_biome==2)),pastr_b,0)
	print 'Anthropogenic grass partitioning is done here for the year before present'
	#total anthropogenic land in the year before present
	totanth_b = c3ann_b + c4ann_b + c3per_b + c4per_b + c3nfx_b + c3angr_b + c4angr_b + urban_b + rangel_b
	#total transition from anthropogenic to forest
	anth_to_forest = urban_to_secdf[nyr-i,:,:] + pastr_to_secdf[nyr-i,:,:] + c3nfx_to_secdf[nyr-i,:,:] + c4per_to_secdf[nyr-i,:,:] + c3per_to_secdf[nyr-i,:,:] + c4ann_to_secdf[nyr-i,:,:] + c3ann_to_secdf[nyr-i,:,:]+range_to_secdf[nyr-i,:,:]
	anth_to_forest = anth_to_forest

	#total transition from anthropogenic area to non forest
	anth_to_nforest = urban_to_secdn[nyr-i,:,:] + pastr_to_secdn[nyr-i,:,:] + c3nfx_to_secdn[nyr-i,:,:] + c4per_to_secdn[nyr-i,:,:] + c3per_to_secdn[nyr-i,:,:] + c4ann_to_secdn[nyr-i,:,:] + c3ann_to_secdn[nyr-i,:,:]+range_to_secdn[nyr-i,:,:]
	anth_to_nforest = anth_to_nforest

	#total transition from forest to anthropogenic area
	forest_to_anth =  primf_to_c3ann[nyr-i,:,:] + primf_to_c4ann[nyr-i,:,:] + primf_to_c3per[nyr-i,:,:] + primf_to_c4per[nyr-i,:,:] + primf_to_c3nfx[nyr-i,:,:] + primf_to_pastr[nyr-i,:,:] + primf_to_urban[nyr-i,:,:]+primf_to_range[nyr-i,:,:]+secdf_to_c3ann[nyr-i,:,:] + secdf_to_c4ann[nyr-i,:,:] + secdf_to_c3per[nyr-i,:,:] + secdf_to_c4per[nyr-i,:,:] + secdf_to_c3nfx[nyr-i,:,:] + secdf_to_pastr[nyr-i,:,:] + secdf_to_urban[nyr-i,:,:] + secdf_to_range[nyr-i,:,:]
	forest_to_anth = forest_to_anth

	#total transition from non forest  to anthropogenic area
	nforest_to_anth =  secdn_to_c3ann[nyr-i,:,:] + secdn_to_c4ann[nyr-i,:,:] + secdn_to_c3per[nyr-i,:,:] + secdn_to_c4per[nyr-i,:,:] + secdn_to_c3nfx[nyr-i,:,:] + secdn_to_pastr[nyr-i,:,:] + secdn_to_urban[nyr-i,:,:]+secdn_to_range[nyr-i,:,:]+primn_to_c3ann[nyr-i,:,:] + primn_to_c4ann[nyr-i,:,:] + primn_to_c3per[nyr-i,:,:] + primn_to_c4per[nyr-i,:,:] + primn_to_c3nfx[nyr-i,:,:] + primn_to_pastr[nyr-i,:,:] + primn_to_urban[nyr-i,:,:]+primn_to_range[nyr-i,:,:]
	nforest_to_anth = nforest_to_anth

	#apply the net transition to forests (add forest loss/remove forest gain in the previous year)
	tottree_b = tottree_p +(forest_to_anth-anth_to_forest)
	#apply the net transition to natural grassland(add loss/remove natural grass gain in the previous year)
	totngr_b   = c3gr_p + c4gr_p + (nforest_to_anth-anth_to_nforest)
	# if not enough area to decrease natural vegetation then encroah baresoil	
	excess_anthb = MV2.where(tottree_b>0.,0.,tottree_b) + MV2.where(totngr_b>0.,0.,totngr_b)
	re_bsoil_b = bsoil_p + excess_anthb
	bsoil_bnow  = MV2.where(re_bsoil_b<0,0.,re_bsoil_b)
	print 'bsoil_bnow',bsoil_bnow.shape
	#x.plot(re_bsoil_b[::-1,:],gm,bg=bg)
	#x.png('re_bsoil_b.png')
	#x.clear()
	tottree_b = MV2.where(tottree_b<0.,0.,tottree_b)
	totngr_b   = MV2.where(totngr_b<0.,0.,totngr_b)
	totnveg_b = tottree_b + totngr_b 
#******************************************************************************
	#Final pft_back after transitions
	pft_back[nyr-j,0,:,:]  = MV2.masked_array(bsoil_bnow,mask=luh_mask)
	pft_back[nyr-j,11,:,:] = MV2.masked_array(c3ann_b,mask=luh_mask)
	pft_back[nyr-j,12,:,:] = MV2.masked_array(c4ann_b,mask=luh_mask)
	pft_back[nyr-j,13,:,:] = MV2.masked_array(c3angr_b,mask=luh_mask)
	pft_back[nyr-j,14,:,:] = MV2.masked_array(c4angr_b,mask=luh_mask)
	pft_back[nyr-j,15,:,:] = MV2.masked_array(c3per_b,mask=luh_mask)
	pft_back[nyr-j,16,:,:] = MV2.masked_array(c4per_b,mask=luh_mask)
	pft_back[nyr-j,17,:,:] = MV2.masked_array(c3nfx_b,mask=luh_mask)
	pft_back[nyr-j,18,:,:] = MV2.masked_array(urban_b,mask=luh_mask)

	#c3/c4 natural grass partitioning
	c3gr_b = totngr_b*0.0
	c4gr_b = totngr_b*0.0
	c3gr_b = MV2.where(((totngr_b!=0) & (pft_biome==3) | (pft_biome==4)|(pft_biome==5)),totngr_b,0)
	c4gr_b = MV2.where(((totngr_b!=0) & (pft_biome==1) | (pft_biome==2)),totngr_b,0)
	#################################
	pft_back[nyr-j,9,:,:]  = c3gr_b
	pft_back[nyr-j,10,:,:] = c4gr_b
	#Adjust Natural vegetation i.e.expand/reduce
	#before adjusting make sure not to overlap vegetation on baresoil
	totveg_back 	 	= totnveg_b  + totanth_b
	totveg_back_notOK	= MV2.masked_values(totveg_back-lndf_hurt,0.).mask
	no_vegmask_back		= MV2.masked_values(totveg_back,0.).mask
	for npft in range(0,8):
		#print 'npft pft_back',npft
		pft_n		= pft_back[nyr-j,npft+1,:,:]
		#locate the points where new forest areas exists
		newTREE	= MV2.masked_where(tottree_p!=0,tottree_b)
		newTREE = MV2.masked_equal(newTREE,0.)
		indx_newTREE = zip(*ma.where(newTREE.mask==False))
		len_indx = len(indx_newTREE)
		print 'No of New TREE grids at year %i' %(nyr-j), len_indx
		for s in range(len(indx_newTREE)):
			#print 's',s
			lat_idx=indx_newTREE[s][0]
 			lon_idx=indx_newTREE[s][1]
			#print 'lat_idx',lat_idx
			#print 'lon_idx',lon_idx
			#Find the average fraction of nearest 5*5 pixels			
			il = max(0,lat_idx-2)
			ir = min(len_indx,lat_idx+2)
			jl = max(0,lon_idx-2)
			jr = min(len_indx,lon_idx+2)
			#print 'il,ir',il,ir
			#print 'jl,jr',jl,jr
			data_tmp = ma.mean(mask_notreef2010[npft,il:ir,jl:jr])
			#print 'data_tmp',data_tmp
			nomask_data  = np.nan_to_num(data_tmp)
			#print 'nomask_data',nomask_data
			treef_2010[npft,lat_idx,lon_idx]  = nomask_data
			#print 'treef_2010[npft,indx_newTREE[s]]',treef_2010[npft,lat_idx,lon_idx]
			if ((nomask_data==0) and (lat_idx>=266 and lat_idx<=454)):
				treef_2010[npft,lat_idx,lon_idx]= MV2.where(npft<=1,0.5,0)
			elif ((nomask_data==0) and (lat_idx<=160)):
				treef_2010[npft,lat_idx,lon_idx]= MV2.where(npft>=5,1/3,0)
			else:
				treef_2010[npft,lat_idx,lon_idx]= MV2.where(1<npft<=5,1/3,0)
				
		pft_tmp3 = MV2.choose(totveg_back_notOK,(MV2.choose(no_vegmask_back,(tottree_b*treef_2010[npft,:,:],0.)),tottree_b*treef_2010[npft,:,:]))		
		print 'Tree fraction classification is done'
		pft_back[[nyr-j],npft+1,:,:]	= MV2.masked_array(pft_tmp3,mask=luh_mask)
		print 'pft_back assignment is done'

	#Sum of all 19 new PFTs
	pft_backSUM = MV2.sum(pft_back[nyr-j,:,:,:], axis=0)+icwtr
	#x.plot(pft_backSUM[::-1,:],gm,bg=bg)
	#x.png('pft_backSUM.png')
	#x.clear()
	#check if sum of all PFTs is one or not
	#if (pft_backSUM==1):
	#	print 'Sum of all PFTs is equal to 1'
	#else: print 'Sum of all PFTs is not equal to 1'
	#In case if sum of all PFTs is greater than 1, reduce all PFTs
	re_pftback = MV2.masked_greater(pft_backSUM - lndf_hurt, 0.).mask
	for npft in range(19):
		#print 'npft',npft
		pft_b = pft_back[nyr-j,npft,:,:]
		pft_b = MV2.choose(re_pftback, (pft_b,pft_b*lndf_hurt/pft_backSUM))
        	pft_back[[nyr-j],npft,:,:] = pft_b
	
	newpftSUM = MV2.sum(pft_back[nyr-j,:,:,:], axis=0)+icwtr
	#In case if sum of all PFTs is less than 1, increase all PFTs
	inc_pftback = MV2.masked_greater(lndf_hurt-newpftSUM, 0.).mask
	for npft in range(19):
		#print 'npft',npft
		pft_b1 = pft_back[nyr-j,npft,:,:]
		pft_b1 = MV2.choose(inc_pftback, (pft_b1,pft_b1*lndf_hurt/newpftSUM))
        	pft_back[[nyr-j],npft,:,:] = pft_b1
	
	inc_pftSUM = MV2.sum(pft_back[nyr-j,:,:,:], axis=0)+icwtr
	#x.plot(inc_pftSUM[::-1,:],gm,bg=bg)
	#x.png('inc_pftSUM.png')
	#x.clear()
	#Store into File from 2009 to 1850
	fileyear = 2010-j
	timax_slide[0] = 2010-j
	g30 = cdms2.open('19PFTmap_%i_cmip6_LUH2v2h.nc' %fileyear,'w')#
	g30.write(pft_back[nyr-j,:,:,:],axes = [timax_slide] + pft19_axislist_pft)
	g30.close()
	#Merge 19 PFTs into 13 PFTs for ORCHIDEE standard version model	
	#*******************************************************************************
	#*******************************************************************************
	#19 PFT names: 0-Baresoil, 1-TropBEf, 2-TropBRf,3-TempNEf, 4-TempBEf,5-TempBSf
	#6-BorNEf, 7-BorBSf, 8-BorNSf, 9-C3Natural grass, 10-C4Natural grass
	#11-C3 annual crops, 12-C4 annual crops, 13-C3 anthropogenic grass
	#14-C4 anthropogenic grass,15-c3 perennial crops, 16-C4 perennial crops
	#17-c3 nitrogen fixing crop, 18-urban
	pft_13[[nyr-j],0,:,:] = pft_back[nyr-j,0,:,:]
	pft_13[[nyr-j],1:9,:,:] = pft_back[nyr-j,1:9,:,:]
	#we add urban to natural grassland
	c3u_b = totngr_b*0.0
	c4u_b = totngr_b*0.0
	c3u_b = MV2.where(((pft_back[nyr-j,18,:,:]!=0) & (pft_biome==3) | (pft_biome==4)|(pft_biome==5)),pft_back[nyr-j,18,:,:],0)
	c4u_b = MV2.where(((pft_back[nyr-j,18,:,:]!=0) & (pft_biome==1) | (pft_biome==2)),pft_back[nyr-j,18,:,:],0)
	pft_13[[nyr-j],9,:,:] = pft_back[nyr-j,9,:,:]+c3u_b+pft_back[nyr-j,13,:,:]
	pft_13[[nyr-j],10,:,:] = pft_back[nyr-j,10,:,:]+c4u_b+pft_back[nyr-j,14,:,:]
	pft_13[[nyr-j],11,:,:] = pft_back[nyr-j,11,:,:]+pft_back[nyr-j,15,:,:]+pft_back[nyr-j,17,:,:]
	pft_13[[nyr-j],12,:,:] = pft_back[nyr-j,12,:,:]+pft_back[nyr-j,16,:,:]
	pft13SUM = MV2.sum(pft_13[nyr-j,:,:,:], axis=0)
	
	g40 = cdms2.open('13PFTmap_%i_cmip6_LUH2v2h.nc' %fileyear,'w')#
	g40.write(pft_13[nyr-j,:,:,:],axes = [timax_slide] + pft13_axislist_pft)
	g40.close()
	#x.plot(pft13SUM[::-1,:],gm,bg=bg)
	#x.png('pft13SUM.png')
	#x.clear()
#Store into File from present day 2010 to 1850 
#g1 = cdms2.open(pft_file_hist19,'w')#
#g1.write(pft_back)
#g1.close()
#g2 = cdms2.open(pft_file_hist13,'w')#
#g2.write(pft_13)
#g2.close()
print '>>>>>>>>>>>>>>>>> Historical PFT map reconstruction is done here<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'	
