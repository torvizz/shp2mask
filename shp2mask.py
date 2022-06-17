#!/opt/anaconda3/envs/tvz/bin/python3.7


######### /home/tovrizz/anaconda3/envs/tvz/bin/python
##############################################
### shapefile to mask based on given array ###
##############################################

import os
import sys
import xarray as xr
import geopandas as gpd
from rasterio import features
from affine import Affine

def check_lat_lon_dim_names(xr_obj):

	if 'latitude' in xr_obj.dims:
		lat_key = 'latitude'
	elif 'lat' in xr_obj.dims:
		lat_key = 'lat'
	else:
		raise ValueError('"lat" or "latitude" dimension name not finded')

	if 'longitude' in xr_obj.dims:
		lon_key = 'longitude'
	elif 'lon' in xr_obj.dims:
		lon_key = 'lon'
	else:
		raise ValueError('"lon" or "longitude" dimension name not finded')

	return lat_key, lon_key



def shp2mask(
	shapely_feature,
	xr_obj,
	maskvalue = 1,
	fillvalue = 0,
	all_touched = False,
	**kwargs
):

	lat_key, lon_key = check_lat_lon_dim_names(xr_obj)

	shape = (xr_obj.dims[lat_key], xr_obj.dims[lon_key])

	lonres = (xr_obj[lon_key][1] - xr_obj[lon_key][0]).data.tolist()
	ullon = xr_obj[lon_key][0].data.tolist() - (lonres/2)

	latres = (xr_obj[lat_key][1] - xr_obj[lat_key][0]).data.tolist()
	ullat = xr_obj[lat_key][0].data.tolist() - (latres/2)
	
	gt = Affine.from_gdal(
		ullon,
		lonres,
		0,
		ullat,
		0,
		latres
	)

	### mantido caso de algum erro no novo jeito implementado acima
	# gt = Affine.from_gdal(
	# 	xr_obj[lon_key][0].data.tolist(),
	# 	(xr_obj[lon_key][1] - xr_obj[lon_key][0]).data.tolist(),
	# 	0,
	# 	xr_obj[lat_key][0].data.tolist(),
	# 	0,
	# 	(xr_obj[lat_key][1] - xr_obj[lat_key][0]).data.tolist()
	# )

	mask_array = [
		features.rasterize(
			[(shapely_feature, maskvalue)],
			out_shape = shape,
			transform = gt,
			fill = fillvalue,
			all_touched = all_touched,
			# dtype = np.uint8
		)
	]

	mask = xr.Dataset(
		{'mask':((lat_key, lon_key),mask_array[0])}, 
		{lat_key:xr_obj.coords[lat_key], lon_key:xr_obj.coords[lon_key]}
	)

	return mask

if __name__ == "__main__":

	help_str = '\nshp2mask.py --shp SHP_PATH --in_nc NETCDF_INPUT_PATH\
	--out_nc NETCDF_OUTPUT_PATH (--bbox N,S,E,W (FLOAT_NUMBERS) --maskvalue INTEGER_NUMBER\
	--fillvalue INTEGER_NUMBER --all_touched INTEGER_NUMBER)\
		\n\n\
		Parameters\n\
		----------\n\n\
		--shp : path\n\
		shapefile path to generate mask\n\n\
		--in_nc : path\n\
		netcdf file path of which 2d array shape will be used\n\
		to generate mask\n\n\
		--out_nc : path\n\
		netcdf file path in which the file will be written\n\n\
		--bbox : opicional, float numbers separated by commas\n\
		coordinates N,S,E,W. if provided, the in_nc is cuted\n\
		before of the mask generation\n\n\
		--maskvalue : opicional, int, default = 1\n\
		pixel value contained in the shapefile\n\n\
		--fillvalue : opicional, int, default = 0\n\
		pixel value out of shapefile\n\n\
		--all_touched : opicional, int, default = 0\n\
		if 0 (false), default, determines that only the\n\
		pixels whose center is contained in the geometry\n\
		are added to the mask. if 1 (true), determines\n\
		that all pixels touched by the geometry boundary\n\
		are added to the mask.'

	if '--help' in sys.argv:
		print(help_str)
		sys.exit()

	arg_keys = [a[2:] for a in sys.argv[1:] if a.startswith('--')]
	arg_vals = [a for a in sys.argv[1:] if not a.startswith('--')]

	args = dict(zip(arg_keys, arg_vals))


	if 'maskvalue' in args.keys():
		if args['maskvalue'][-1].isdigit():
			args['maskvalue'] = int(args['maskvalue'])
		elif args['maskvalue'].lower() == 'nan':
			args['maskvalue'] = np.nan
		else:
			raise ValueError('maskvalue must be numeric or "nan"')

	if 'fillvalue' in args.keys():
		if args['fillvalue'][-1].isdigit():
			args['fillvalue'] = int(args['fillvalue'])
		elif args['fillvalue'].lower() == 'nan':
			args['fillvalue'] = np.nan
		else:
			raise ValueError('fillvalue must be numeric or "nan"')

	if 'all_touched' in args.keys():
		if args['all_touched'] != '0' and args['all_touched'] != '1':
			raise ValueError('all_touched must be 0 to False or 1 to True')
		else:
			args['all_touched'] = bool(int(args['all_touched']))


	# if 'mode' not in args.keys():
	#     args['mode'] = 'gen_mask'

	# if args['mode'] == 'gen_mask':

	in_nc = xr.open_dataset(args.pop('in_nc'))
	# if (in_nc.longitude >180).any().data.tolist():
	#     in_nc['longitude'] = in_nc.longitude-360

	if 'bbox' in args.keys():
		N, S, E, W = [float(c) for c in args.pop('bbox').split(',')]

		bbox_dict = {}
		for k in in_nc.dims.keys():
			if k.startswith('lon'):
				bbox_dict[k] = slice(W, E)
			elif k.startswith('lat'):
				if (in_nc[k][0]-in_nc[k][1]) < 0:
					bbox_dict[k] = slice(S, N)
				else:
					bbox_dict[k] = slice(N, S)

		in_nc = in_nc.sel(bbox_dict)
	# in_nc = in_nc.sel(latitude=slice(S, N), longitude=slice(W, E))

	shp = gpd.read_file(args.pop('shp'))
	shp = shp.geometry.values[0]

	mask = shp2mask(shp, in_nc, **args)
		
	mask.to_netcdf(args['out_nc'])

# elif:

#     in_nc = xr.open_dataset(args.pop('in_nc'))
#     in_nc['longitude'] = in_nc.longitude-360

#     bbox_dict = {}
#     for k in in_nc.dims.keys():
#         if k.startswith('lon'):
#             bbox_dict[k] = slice(W, E)
#         elif k.startswith('lat'):
#             bbox_dict[k] = slice(S, N)

#     in_nc = in_nc.sel(bbox_dict)
#     # in_nc = in_nc.sel(latitude=slice(S, N), longitude=slice(W, E))

#     shp = gpd.read_file(args.pop('shp'))
#     shp = shp.geometry.values[0]

#     mask = shp2mask(shp, in_nc, **args)

#     out_nc = in_nc.where(mask.mask==1,drop=True)

# run scripts2/shp2mask.py --shp /media/Dados/Shapes/Bacias_Costa_do_Brasil/Talude_Santos.shp --in_nc /media/Dados/Trans_Dados/BS/CCMP/Y1988/M01/CCMP_Wind_Analysis_19880101_V02.0_L3.0_RSS.nc --out_nc testefinal.nc --bbox -20,-30,-39,-49

# args = {}
# args['shp'] = '/media/Dados/Shapes/Bacias_Costa_do_Brasil/Plataforma_Santos.shp'
# args['in_nc'] = '/media/Dados/Trans_Dados/BS/CCMP/Y1988/M01/CCMP_Wind_Analysis_19880101_V02.0_L3.0_RSS.nc'
# args['bbox'] = '-20,-30,-39,-49'
# args['maskvalue'] = '2'
# args['fillvalue'] = 'Na'
# args['all_touched'] = '1'

# for k, v in args.items():

#     if k == 'in_value':
#         args[k] = int(v)

#     else:
#         args[] = 1

# def shp2mask(
#     shapely_feature = args['shp'],
#     xr_obj = args['in_nc'],
#     in_value = args['inside_value'],
#     out_value = args['outside_value'],
#     all_touched = args['all_touched']
# ):


#############################################
# lembrete:
# verificar se funções com kwargs aceitam dicionários vazios para poder 
# retiras as entradas principais da função e passar um dict vazio caso
# não tenha sido passado nenhuma entrada não essencial

# f, ax = plt.subplots(1,3, figsize=(10,6))

# v18s1.sel().plot(ax=ax[0])
# v18s2.sel().plot(ax=ax[1])
# tfabio.wspd.sel().plot(ax=ax[2])
# names=['centro contido', 'todos que tocam', 'fabio']
# for a, n in zip(ax, names):
#     a.set_ylim([-30,-20])
#     a.set_xlim([-50,-40])
#     a.set_title(n)
#     s1.plot(ax=a, facecolor='none',ec='k')
# f.tight_layout()






