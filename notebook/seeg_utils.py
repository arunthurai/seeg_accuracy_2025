import os
import re
import pandas as pd
import numpy as np
import math
from statistics import NormalDist, mean
from collections import ChainMap

##HANDLING SLICER FCSV
def determineFCSVCoordSystem(input_fcsv,overwrite_fcsv=False):
	# need to determine if file is in RAS or LPS
	# loop through header to find coordinate system
	coordFlag = re.compile('# CoordinateSystem')
	verFlag = re.compile('# Markups fiducial file version')
	headFlag = re.compile('# columns')
	coord_sys=None
	headFin=None
	ver_fin=None

	with open(input_fcsv, 'r') as myfile:
		firstNlines=myfile.readlines()[0:3]

	for row in firstNlines:
		row=re.sub("[\s\,]+[\,]","",row).replace("\n","")
		cleaned_dict={row.split('=')[0].strip():row.split('=')[1].strip()}
		if None in list(cleaned_dict):
			cleaned_dict['# columns'] = cleaned_dict.pop(None)
		if any(coordFlag.match(x) for x in list(cleaned_dict)):
			coord_sys = list(cleaned_dict.values())[0]
		if any(verFlag.match(x) for x in list(cleaned_dict)):
			verString = list(filter(verFlag.match,  list(cleaned_dict)))
			assert len(verString)==1
			ver_fin = verString[0].split('=')[-1].strip()
		if any(headFlag.match(x) for x in list(cleaned_dict)):
			headFin=list(cleaned_dict.values())[0].split(',')


	if any(x in coord_sys for x in {'LPS','1'}):
		df = pd.read_csv(input_fcsv, skiprows=3, header=None)

		if df.shape[1] != 13:
			df=df.iloc[:,:14]

		df[1] = -1 * df[1] # flip orientation in x
		df[2] = -1 * df[2] # flip orientation in y

		if overwrite_fcsv:
			with open(input_fcsv, 'w') as fid:
				fid.write("# Markups fiducial file version = 4.11\n")
				fid.write("# CoordinateSystem = 0\n")
				fid.write("# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n")

			df.rename(columns={0:'node_id', 1:'x', 2:'y', 3:'z', 4:'ow', 5:'ox',
								6:'oy', 7:'oz', 8:'vis', 9:'sel', 10:'lock',
								11:'label', 12:'description', 13:'associatedNodeID'}, inplace=True)

			df['associatedNodeID']= pd.Series(np.repeat('',df.shape[0]))
			df.round(6).to_csv(input_fcsv, sep=',', index=False, lineterminator="", mode='a', header=False, float_format='%.6f')

			print(f"Converted LPS to RAS: {os.path.dirname(input_fcsv)}/{os.path.basename(input_fcsv)}")
	return coord_sys,headFin


def df_to_fcsv(input_df, output_fcsv):
	with open(output_fcsv, 'w') as fid:
		fid.write("# Markups fiducial file version = 4.11\n")
		fid.write("# CoordinateSystem = 0\n")
		fid.write("# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n")

	out_df={'node_id':[],'x':[],'y':[],'z':[],'ow':[],'ox':[],'oy':[],'oz':[],
		'vis':[],'sel':[],'lock':[],'label':[],'description':[],'associatedNodeID':[]
	}

	for idx,ifid in input_df.iterrows():
		out_df['node_id'].append(idx+1)
		out_df['x'].append(ifid.iloc[0])
		out_df['y'].append(ifid.iloc[1])
		out_df['z'].append(ifid.iloc[2])
		out_df['ow'].append(0)
		out_df['ox'].append(0)
		out_df['oy'].append(0)
		out_df['oz'].append(0)
		out_df['vis'].append(1)
		out_df['sel'].append(1)
		out_df['lock'].append(1)
		out_df['label'].append(str(ifid.iloc[3]))
		out_df['description'].append(str(ifid.iloc[4]))
		out_df['associatedNodeID'].append('1')

	out_df=pd.DataFrame(out_df)
	out_df.to_csv(output_fcsv, sep=',', index=False, lineterminator="", mode='a', header=False, float_format = '%.3f')

def sorted_nicely(lst):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	sorted_lst = sorted(lst, key = alphanum_key)

	return sorted_lst

def determine_groups(iterable, numbered_labels=False):
	values = []
	for item in iterable:
		temp=None
		if re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item):
			temp = "".join(list(re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item)[0]))
		elif '-' in item:
			temp=item.split('-')[0]
		else:
			if numbered_labels:
				temp=''.join([x for x in item if not x.isdigit()])
				for sub in ("T1","T2"):
					if sub in item:
						temp=item.split(sub)[0] + sub
			else:
				temp=item
		if temp is None:
			temp=item

		values.append(temp)

	vals,indexes,count = np.unique(values, return_index=True, return_counts=True)
	values_unique = [values[index] for index in sorted(indexes)]

	return values_unique,vals

##METRICS

def euclidianDistanceCalc(xyz_planned, xyz_actual):
	##greydon/jon
  #print(xyz_planned.ndim)
	if xyz_planned.ndim>1:
		euc_dist=[]
		for ipoint in range(xyz_planned.shape[0]):
			plan_act_diff = xyz_planned[ipoint] - xyz_actual[ipoint]
			euc_dist.append(math.sqrt(sum(plan_act_diff**2)))
	else:
		plan_act_diff = xyz_planned - xyz_actual
		euc_dist = math.sqrt(sum(plan_act_diff**2))
	return euc_dist


def euclidean_distance(point1, point2):
    return np.linalg.norm(point1 - point2)

def radialDistanceCalc(pt, xyz_entry, xyz_target):
    if xyz_entry.ndim > 1:
        dist3d = []
        for ipoint in range(xyz_entry.shape[0]):
            x1_minus_pt = pt[ipoint] - xyz_entry[ipoint]
            x2_minus_x1 = xyz_target[ipoint] - xyz_entry[ipoint]

            sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)
            sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

            mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

            dist3d.append(np.sqrt(
                (sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1))
    else:
        x1_minus_pt = pt - xyz_entry
        x2_minus_x1 = xyz_target - xyz_entry

        sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)
        sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

        mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

        dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 -
                         (mydotprod * mydotprod))/sumsq_x2_minus_x1)
    return dist3d


def ptLineAngleCalc(pt, x_entry, x_target):
    if x_entry.ndim > 1:
        deg_angle = []
        for ipoint in range(x_entry.shape[0]):
            try:
                x1_minus_pt = pt[ipoint] - x_entry[ipoint]
                x2_minus_x1 = x_target[ipoint] - x_entry[ipoint]

                sumsq_x1_minus_pt = sum(x1_minus_pt**2)
                sumsq_x2_minus_x1 = sum(x2_minus_x1**2)

                # sum of products of elements
                mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

                rad_angle = math.acos(
                    mydotprod/(np.sqrt(sumsq_x1_minus_pt)*np.sqrt(sumsq_x2_minus_x1)))
                deg_angle.append(math.degrees(rad_angle))
            except:
                deg_angle.append(np.nan)
                print(f"Check point {ipoint}")
    else:
        x1_minus_pt = pt - x_entry
        x2_minus_x1 = x_target - x_entry

        sumsq_x1_minus_pt = sum(x1_minus_pt**2)
        sumsq_x2_minus_x1 = sum(x2_minus_x1**2)

        # sum of products of elements
        mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

        rad_angle = math.acos(
            mydotprod/(np.sqrt(sumsq_x1_minus_pt)*np.sqrt(sumsq_x2_minus_x1)))
        deg_angle = math.degrees(rad_angle)
    return deg_angle

# Line angle calculation


def lineLineAngleCalc(a_entry, a_target, b_entry, b_target):
    if a_entry.ndim > 1:
        deg_angle = []
        for ipoint in range(a_entry.shape[0]):
            try:
                vectorA = a_target[ipoint] - a_entry[ipoint]
                vectorB = b_target[ipoint] - b_entry[ipoint]

                sumsq_vectorA = sum(vectorA**2)
                sumsq_vectorB = sum(vectorB**2)

                mydotprod = sum(vectorA*vectorB)

                rad_angle = math.acos(
                    mydotprod/(np.sqrt(sumsq_vectorA)*np.sqrt(sumsq_vectorB)))
                deg_angle.append(math.degrees(rad_angle))
            except:
                deg_angle.append(np.nan)
                print(f"Check point {ipoint}")
    else:
        vectorA = a_target - a_entry
        vectorB = b_target - b_entry

        sumsq_vectorA = sum(vectorA**2)
        sumsq_vectorB = sum(vectorB**2)

        mydotprod = sum(vectorA*vectorB)

        rad_angle = math.acos(
            mydotprod/(np.sqrt(sumsq_vectorA)*np.sqrt(sumsq_vectorB)))
        deg_angle = math.degrees(rad_angle)
    return deg_angle

def normDir(direction):
  return np.array(direction) / np.linalg.norm(np.array(direction))

def mag_vec(P1, P2):
  if isinstance(P1, list):
    P1 = np.array(P1)
  if isinstance(P1, list):
    P2 = np.array(P2)
  DirVec = P2-P1
  MagVec = np.sqrt([np.square(DirVec[0]) + np.square(DirVec[1]) + np.square(DirVec[2])])
  return MagVec

def norm_vec(P1, P2):
  if isinstance(P1, list):
    P1 = np.array(P1)
  if isinstance(P2, list):
    P2 = np.array(P2)
  DirVec = P2-P1
  MagVec = np.sqrt([np.square(DirVec[0]) + np.square(DirVec[1]) + np.square(DirVec[2])])
  NormVec = np.array([float(DirVec[0] / MagVec), float(DirVec[1] / MagVec), float(DirVec[2] / MagVec)])
  return NormVec

def mean_confidence_interval(data, confidence=0.95):
	dist = NormalDist.from_samples(data[~np.isnan(data)])
	z = NormalDist().inv_cdf((1 + confidence) / 2.)
	h = dist.stdev * z / ((len(data) - 1) ** .5)
	return dist.mean - h, dist.mean + h