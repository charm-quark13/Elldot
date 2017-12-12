import numpy
import gc
from sys import maxsize
import sys
import os
import matplotlib
#matplotlib.use('pdf')
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d


# utitlity functions
def PlotPHM(filename, phm_matrix = None, inverse = False, fmt = "pdf"):
	"""
	Given a particle-hole matrix or file, plot the particle-hole map.
	
	Parameters
	----------
	filename : file which contains the matrix of the particle-hole map or the file name to write the map.
	phm_matrix: the particle-hole matrix. It should be a numpy array
	inverse : turn the sign of the particle-hole map. (it is OK in the linear response theory).
	fmt: output format
	
	Returns
	-------
	no return. But there will be a new file named filename + '.pdf' which contains the plot in pdf format.
	"""
	try:
		phm = numpy.loadtxt(filename)
	except:
		if(type(phm_matrix) != None):
			phm = phm_matrix
		else:
			print("Error: please provide a file or a matrix as input")
	if(inverse): phm = -phm
	Lx = phm.shape[0]
	Ub = (numpy.ndarray.max(phm) + numpy.std(phm))/2
	#fig = pyplot.figure()
	cdict1 = {'red':  ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
			'green': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
			'blue':  ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0))}
	phm_standard = LinearSegmentedColormap('PH_map_standard', cdict1)
	pyplot.imshow(  phm   , interpolation = None, extent = (0.5,Lx+0.5,0.5,Lx+0.5), cmap = phm_standard, vmin = -Ub, vmax = Ub, origin='lower')
	pyplot.colorbar()
	filenameext = filename.find(".")
	if(filenameext == -1):
		file = filename
	else:
		file = filename[:filenameext]
	file = file + "."
	try:
		pyplot.savefig(file + fmt)
		pyplot.close()
		print(file + fmt, " is saved.")
	except:
		print("Error: the format is not support by matplotlib")



def Savebov(r0, nstep, data, filename):
	"""
	Save a 3D object (e.g. wave function) into binary file with the bov file (can be opened in VisIt to visualize).

    Parameters
    ----------
    r0: orient of the 3D object
	nstep: grid step
	data: the matrix
	filename: output filename. There will be filename.bin and filename.bov in the directory.

    Returns
    -------
    No return
    """	
	nx, ny, nz = data.shape
	lbox = nstep*numpy.array((nx,ny,nz))
	d2 = data.reshape(nx*ny*nz, order='F')
	d2.tofile(filename+".bin")
	bovf = open(filename+".bov", 'w')
	bovf.write("TIME:   0.000000000000E+00" + "\n")
	bovf.write("DATA_FILE: " + filename + ".bin" + "\n")
	bovf.write("DATA_SIZE: " + str(nx) + "   " + str(ny) + "   " + str(nz) + "\n")
	bovf.write("DATA_FORMAT:  DOUBLE" + "\n")
	bovf.write("VARIABLE:  density response" + "\n")
	bovf.write("DATA_ENDIAN:  LITTLE" + "\n")
	bovf.write("CENTERING:  nodal" + "\n")
	bovf.write("BRICK_ORIGIN:   " + str(r0[0]) + "   " + str(r0[1]) + "   " + str(r0[2]) + "\n")
	bovf.write("BRICK_SIZE:   " + str(lbox[0]) + " " + str(lbox[1]) + " " + str(lbox[2]) + "\n")
	bovf.write("BYTE_OFFSET: 0" + "\n")



def ShrinkMap(filename, reduce_every):
    """
    Given a file with matrix elements stored, plot the particle-hole map.

    Parameters
    ----------
    filename : file which contains the matrix of the particle-hole map.

    Returns
    -------
    no return. But there will be a new file named filename + 'small' which contains the smaller map.
    """
    phm = numpy.loadtxt(filename)
    small = numpy.empty((int(phm.shape[0]/reduce_every), int(phm.shape[0]/reduce_every)))
    print("new size is ", small.shape)
    for i in range(0, int(phm.shape[0]/reduce_every)):
        for j in range(0, int(phm.shape[0]/reduce_every)):
            small[i,j] = numpy.sum(phm[i*reduce_every:(i+1)*reduce_every, j*reduce_every:(j+1)*reduce_every])
    numpy.savetxt(filename + '_small', small)



def MakeGauge(filegeo, filebov, filegeoimg):
    """
    Given the geometry, lattice information and molecular image. Add the gauge to the image for PHM reading.

    Parameters
    ----------
    filegeo (must be xyz format): Geometry file in xyz format.
    filebov (must be bov format): VisIt file (.bov format) of a density, wave function or anything to extract the lattice information.
    filegeoimg (must be png format): A png image with the geometry contained. The background must be white for the program to understand where the molecule is.

    Returns
    -------
    no return. But there will be a new file named gauge.pdf which contains the plot in pdf format.
    """
    geo = numpy.genfromtxt(filegeo, skiprows = 2)
    geomax = numpy.max(geo[:,1])
    geomin = numpy.min(geo[:,1])
    fp = open(filebov)
    l1 = fp.readline()
    l1 = fp.readline()
    l1 = fp.readline().split()
    nx = int(l1[1])
    ny = int(l1[2])
    nz = int(l1[3])
    l1 = fp.readline()
    l1 = fp.readline()
    l1 = fp.readline()
    l1 = fp.readline()
    l1 = fp.readline().split()
    lx = float(l1[1])
    ly = float(l1[2])
    lz = float(l1[3])
    gaugex = numpy.arange(lx, -lx, -2*lx/(nx-1))
    numpy.append(gaugex, -lx)
    idmin = abs(gaugex - geomin).argmin()
    if(gaugex[idmin] - geomin > 0): --idmin
    idmax = abs(gaugex - geomax).argmin()
    if(gaugex[idmax] - geomax > 0): --idmax
    image = pyplot.imread(filegeoimg)
    avgimage = numpy.mean(image, axis = (0,2))
    iml = numpy.where(avgimage != 1)[0][0]
    imr = numpy.where(avgimage != 1)[0][-1]
    id0 = numpy.ceil((idmax - idmin)/(imr - iml)*(-iml) + idmin)
    idn = numpy.ceil((idmax - idmin)/(imr - iml)*(image.shape[1] - imr) + idmax)
    idy = int((idn - id0 + 1)/image.shape[1]*image.shape[0])
    fig = pyplot.imshow(image, extent = (id0, idn, 1, idy))
    fig.axes.get_yaxis().set_visible(False)
    #pyplot.tick_params(axis='y', which='both', bottom='off', top='off', labelbottom='off')
    pyplot.savefig('gauge.pdf')



#Gaussian post process
def GaussianCubeWFLoader(file, skip_first = 0, skip_last = 0, first_orb_in_file = 5):
	"""
	Load wave functions from Gaussian cube file.
	Storing the wave functions are expensive in memory. Choose the skip parameters wisely.

    Parameters
    ----------
    file: the wave function file
	skip_first: first n orbitals to ignore (default 0)
	skip_last: last n orbitals to ignore (default 0)
	first_orb_in_file: first orbital number in the cube file (default 1). Sometimes cube file doesn't hold all orbitals.

    Returns
    -------
    wave functions in a list. if the ignore parameters are set then the corresponding item has 0 length
    """
	f = open(file)
	f.readline()
	f.readline()
	line = f.readline()
	natoms = -int(line.split()[0])
	r0 = []
	r0.append(float(line.split()[1]))
	r0.append(float(line.split()[2]))
	r0.append(float(line.split()[3]))
	ngrid = []
	dr = []
	line = f.readline()
	ngrid.append(int(line.split()[0]))
	dr.append(float(line.split()[1]))
	line = f.readline()
	ngrid.append(int(line.split()[0]))
	dr.append(float(line.split()[2]))
	line = f.readline()
	ngrid.append(int(line.split()[0]))
	dr.append(float(line.split()[3]))
	for iatom in range(0, natoms):
		f.readline()
	line = f.readline()
	nwf = int(line.split()[0])
	print("grid start at: ", r0)
	print("number of grid: ", ngrid)
	print("grid space: ", dr)
	print("number of wf: ", nwf)
	print("will skip the first ", skip_first, " orbitals")
	print("will skip the last ", skip_last, " orbitals")
	ng = ngrid[0]*ngrid[1]*ngrid[2]
	data = []
	for iwf in range(0, nwf):
		if(iwf in range(0, skip_first) or iwf in range(nwf - skip_last, nwf)):
			data.append([])
		else:
			data.append(numpy.empty([ng], dtype = numpy.float32))
	for line in f:
		if(float(line.split()[0]) >= 0.999):
			pass
		else:
			break
	i = 0
	iwf = 0
	print("the first line of the cube is")
	print(line)
	for word in line.split():
		if(iwf < skip_first):
			pass
		else:
			data[iwf][i] = float(word)
		iwf += 1
	try:
		for line in f:
			for word in line.split():
				if(iwf in range(0, skip_first) or iwf in range(nwf - skip_last, nwf)):
					pass
				else:
					data[iwf][i] = float(word)
				if(iwf + 1 == nwf):
					iwf = 0
					i += 1
				else:
					iwf += 1
	except:
		print("error when reading line")
		print(line)
		print("to the orbital storage")
		print(iwf)
		print("at index")
		print(i)
	wf = [[]]*nwf
	for iwf in range(skip_first,nwf - skip_last):
		wf[iwf] = numpy.reshape(data[iwf], ngrid)
	#Make the first wf has index 0
	for i in range(0, first_orb_in_file - 1):
		wf.insert(0,[])
	return r0, dr, ngrid, wf



def GaussianCoeffLoader(file, cutoff = 0.0001):
	"""
	load TD (casida) coefficients
	and return a set with element of "from" "to" "coeff"
	"""
	minorbital = maxsize
	maxorbital = 1
	f = open(file)
	excitations = []
	for line in f:
		key = line.split()
		if1 = "Excited" in key
		if2 = "State" in key
		theexcitation = []
		coeffsum = 0.0
		if(if1 and if2):
			while True:
				tdline = f.readline()
				if(tdline.find("->") == -1 and tdline.find("<-") == -1):
					break
				flag = 0
				sep = tdline.find("->")
				if(sep == -1):
					sep = tdline.find("<-")
					flag = 1
				data = tdline[sep+2:].split()
				if (flag == 0):
					print(data[:])
				if(abs(float(data[1])) < cutoff):
					continue
				orb1 = int(tdline[0:sep])
				if(orb1 < minorbital):
					minorbital = orb1
				orb2 = int(data[0])
				if(orb2 > maxorbital):
					maxorbital = orb2
				print(flag)
				theexcitation.append([orb1, flag, orb2, float(data[1])])
				coeffsum += theexcitation[-1][-1]**2
			print("this excitation with sum: ", coeffsum)
			excitations.append(theexcitation)
	print("lowest excitation orbital : ", minorbital)
	print("highest excitation orbtial: ", maxorbital)
	print(excitations[0])
	return excitations

# def Slicer(func, direction = 'x'):
	# dir = 0
	# if(direction == 'x'):
		# dir = 0
	# elif(direction == 'y'):
		# dir = 1
	# elif(direction == 'z'):
		# dir = 2
	# else:
		# print("error")
	# dim = func.shape
	# wf_part = numpy.zeros([dim[dir]])
	# for ip in range(0,dim[dir]):
		# if(direction == 'x'):
			# wf_part[ip] = numpy.sum(func[ip,:,:])
		# elif(direction == 'y'):
			# wf_part[ip] = numpy.sum(func[:,ip,:])
		# elif(direction == 'z'):
			# wf_part[ip] = numpy.sum(func[:,:,ip])
	# return wf_part

def LoadGeoXYZ(file, take_hydrogen = True):
	f = open(file)
	line = f.readline()
	natoms = int(line)
	print(natoms)
	line = f.readline()
	print("comment of the file:", line)
	geo = []
	for line in f:
		if(line.strip() == ""):
			continue
		word = line.split()
		if(take_hydrogen):
			pass
		else:
			if(word[0].lower() == "h"):
				natoms -= 1
				continue
		geo.append([word[0], float(word[1]), float(word[2]), float(word[3])])
	print("load ", natoms, "atoms")
	return geo

def xyzswitch(file, ind1, ind2):
	f = open(file, 'r')
	lines = []
	for line in f:
		lines.append(line)
	f.close()
	assert(ind1 >= 0)
	assert(ind2 >= 0)
	temp = lines[ind1 + 2]
	lines[ind1 + 2] = lines[ind2 + 2]
	lines[ind2 + 2] = temp
	f = open(file, 'w')
	for line in lines:
		f.write(line)
	f.close()

def Partition(option, ngrid):
	"""
	Generate partition template.

    Parameters
    ----------
    option: list, the details of partition.
		option[0] is the partition style ("box"/"atomic")
		option[1:], the details of the partition.
			For "box" partition, one can use, [nboxx, nboxy, nboxz] and -1 for slice along the direction.
			For "atomic" partition: option[1] = geometry of molecule, option[2] = r0, option[3] = step

    Returns
    -------
    number of partitions and 3D array with partition index
    """
	Weight = 0.0
	PartInd = 0.0
	npart = 0
	if(option[0].lower() == "box"):
		option = option[1:]
		if(option[0] <= 0):
			option[0] = ngrid[0]
		if(option[1] <= 0):
			option[1] = ngrid[1]
		if(option[2] <= 0):
			option[2] = ngrid[2]
		npart = numpy.prod(option)
		Weight = numpy.ones(ngrid)
		PartInd = numpy.zeros(ngrid)
		PartInd.fill(-1)
		lbins = numpy.array((ngrid[0]//option[0],ngrid[1]//option[1],ngrid[2]//option[2]))
		print(lbins)
		for z in range(0, ngrid[2] - ngrid[2]%option[2]):
			for y in range(0, ngrid[1] - ngrid[1]%option[1]):
				for x in range(0, ngrid[0] - ngrid[0]%option[0]):
					PartInd[x,y,z] = x//lbins[0] + y//lbins[1]*option[0] + z//lbins[2]*option[0]*option[1]
	elif(option[0].lower() == "atomic"):
		atoms = option[1]
		r0 = option[2]
		lstep = option[3]
		if(atoms == None):
			raise ValueError('Cannot make this partition without geometry info.')
		if(r0 == None or lstep == None):
			raise ValueError('Cannot make this partition without knowing the box info.')
		natoms = len(atoms)
		atoms_grid = numpy.zeros((natoms,3))
		npart = natoms
		Weight = numpy.ones(ngrid)
		PartInd = numpy.zeros(ngrid)
		for i in range(0, natoms):
			atoms_grid[i,0] = (atoms[i][1] - r0[0])/lstep[0]
			atoms_grid[i,1] = (atoms[i][2] - r0[1])/lstep[1]
			atoms_grid[i,2] = (atoms[i][3] - r0[2])/lstep[2]
		dist = numpy.zeros((natoms))
		for z in range(0, ngrid[2]):
			for y in range(0, ngrid[1]):
				for x in range(0, ngrid[0]):
					for i in range(0, natoms):
						dist[i] = numpy.linalg.norm([x,y,z] - atoms_grid[i,:])
					PartInd[x,y,z] = numpy.argmin(dist)
	else:
		raise ValueError('There is no such partition mode by now.')
	return npart, PartInd

def PartitionApply(npart, PartInd, func):
	func_part = numpy.zeros((npart))
	for i in range(0, npart):
		func_part[i] = numpy.sum(func*(PartInd == i))
	return func_part

def CalDensResp(nks, wf, coeff):
	for i in range(0, len(wf)):
		if(len(wf[i]) == 0):
			pass
		else:
			break
	# split each orbital response by the order of coefficients
	respbuf = numpy.zeros([15, nks, wf[i].shape[0], wf[i].shape[1], wf[i].shape[2]])
	for excit in coeff:
		respbuf[-int(numpy.log10(abs(excit[3]))),excit[0]-1,:,:,:] += excit[3]*wf[excit[0]-1]*wf[excit[2]-1]
	resp = numpy.zeros([nks, wf[i].shape[0], wf[i].shape[1], wf[i].shape[2]])
	for oom in range(14,-1,-1):
		resp[:,:,:,:] += respbuf[oom,:,:,:,:]
	del respbuf
	return resp

def BuildPHM(orbdens0, orbdens1, iter_from = 0):
	PHM = numpy.zeros([orbdens0.shape[1], orbdens0.shape[1]])
	for iwf in range(iter_from,orbdens0.shape[0]):
		PHM += numpy.outer(orbdens0[iwf,:], orbdens1[iwf,:])
	return PHM


# class MolIndexReorder(object):
	# """
	# """
	# def __init__(self, file_geo, take_hydrogen = True):
		# self.geo = LoadGeoXYZ(file_geo, take_hydrogen)
		# self.coord = numpy.zeros((len(self.geo),3))
		# self.labels = numpy.zeros(len(self.geo), dtype = 'U5')
		# for i in range(0, len(self.geo)):
			# self.coord[i,0] = self.geo[i][1]
			# self.coord[i,1] = self.geo[i][2]
			# self.coord[i,2] = self.geo[i][3]
			# self.labels[i] = str(i)
	
	# def ViewMol(self):
		# self.fig = pyplot.figure()
		# self.ax = self.fig.add_subplot(111, projection='3d')
		# for i in range(0, len(self.geo)):
			# self.ax.scatter(self.coord[i,0], self.coord[i,1], self.coord[i,2])
			# self.ax.text(self.coord[i,0], self.coord[i,1], self.coord[i,2], self.labels[i])
		# pyplot.show()
		# pyplot.close()


class GaussianPostProcess(object):
	"""
	A class constructed from Gaussian output and compute the particle-hole map with your choice of partition mathod.
	"""
	def __init__(self, file_cube, file_output, nks, file_geo, skip_first = 0, skip_last = 0, first_orb_in_file = 1, take_hydrogen = True, cutoff = 0.005):
		self.file_cube = file_cube
		self.file_output = file_output
		self.file_geo = file_geo
		self.nks = nks
		print("the cube file is:      ", file_cube)
		print("the output file is:    ", file_output)
		print("the geometry file is:  ", file_geo)
		print("There are", nks, "Kohn-Sham orbitals")
		print("Loading the wave functions from cube file...")
		self.r0, self.dr, self.ngrid, self.wf = GaussianCubeWFLoader(file_cube, skip_first = skip_first, skip_last = skip_last, first_orb_in_file = first_orb_in_file)
		self.vol = numpy.prod(self.dr)
		print("done")
		if(file_geo == None):
			pass
		else:
			self.geo = LoadGeoXYZ(file_geo, take_hydrogen)
			self.natoms = len(self.geo)
		print("Loading the excitation info and coefficients...")
		self.coeffs = GaussianCoeffLoader(file_output, cutoff)
		self.first_wf = first_orb_in_file - 1 + skip_first
		self.flag_phmforapoint = 0
	
	def CalculateDensityResponse(self, whichexcitation = 'All'):
		excitation_list = []
		if(whichexcitation == 'All'):
			excitation_list = numpy.arange(len(self.coeffs), dtype = int)
		else:
			excitation_list = whichexcitation
		print("Will study the excitations in the following list:\n")
		print(excitation_list)
		print("loop for the excitations\n")
		for i in excitation_list:
			print("excitation", i)
			orbdens1 = CalDensResp(self.nks, self.wf, self.coeffs[i])
			dn1 = numpy.sum(orbdens1, axis = 0)
			Savebov(self.r0, self.dr, dn1, "dn" + str(i+1))
			norm = dn1[dn1 > 0].sum()
			print("norm = " + str(norm) + "\n")
			print("\n")
		print("done")
	
	def UpdatePartition(self, partition_choice, par_info = None):
		print("making partition template...")
		if(partition_choice.lower().strip() == "box"):
			self.npart, self.par = Partition([partition_choice.lower().strip(), par_info[0], par_info[1], par_info[2]], self.ngrid)
		elif(partition_choice.lower().strip() == "atomic"):
			self.npart, self.par = Partition([partition_choice.lower().strip(), self.geo, self.r0, self.dr], self.ngrid)
		else:
			raise ValueError('There is no such partition mode by now.')
		print("done")
		print("partition mode is:       ", partition_choice)
		print("number of partitions is: ", self.npart)
		print("creating Part[|psi0|^2]")
		self.orbdens0_part = numpy.zeros((self.nks, self.npart))
		for i in range(0, self.nks):
			print(i, end="")
			if(len(self.wf[i]) == 0):
				continue
			else:
				orbdens0 = self.wf[i]**2
				self.orbdens0_part[i,:] = PartitionApply(self.npart, self.par, orbdens0)*self.vol
		print("done")
	
	def CalculatePHM(self, whichexcitation = [0,1,2,3,4,5,6,7,8,9,10], file_type = "pdf", ifnormalize = True):
		excitation_list = []
		if(whichexcitation == 'All'):
			excitation_list = numpy.arange(len(self.coeffs), dtype = int)
		else:
			excitation_list = whichexcitation
		print("Will study the excitions in the following list:\n")
		print(excitation_list)
		print("loop for the excitations\n")
		self.PHM = []
		orbdens1_part = numpy.zeros([self.nks,self.npart])
		for i in excitation_list:
			print("excitation", i)
			orbdens1 = CalDensResp(self.nks, self.wf, self.coeffs[i])
			dn1 = numpy.sum(orbdens1, axis = 0)
			Savebov(self.r0, self.dr, dn1, "paper-ex" + str(i+1))
			norm = 1
			if(ifnormalize):
				norm = dn1[dn1 > 0].sum()
			print("norm = " + str(norm) + "\n")
			for iks in range(self.first_wf, self.nks):
				print(iks, end="")
				orbdens1_part[iks,:] = PartitionApply(self.npart, self.par, orbdens1[iks,:,:,:])*self.vol
			self.PHM.append(1.0/norm*BuildPHM(orbdens1_part, self.orbdens0_part, self.first_wf))
			PlotPHM( 'paper-PHM' + str(i+1), phm_matrix = self.PHM[-1], fmt = file_type)
			print("\n")
		print("done")
	
	def PHMForApoint(self, excitation, OriginLabels = [-1], DestiLabels = [-1], SumOriginOrDesti = "origin"):
		"""
		Produce a PHM(r,r') to PHM(r,r0') or PHM(r0, r')
		
		Parameters
		----------
		excitation: int
		            index of the excitation. 1 is the first excitation.
		
		OriginLabels: list. 
		            study the CT from the locations in the list. [-1] for all locations
		
		DestiLabels: list. 
		            study the CT to the locations in the list. [-1] for all locations
		
		SumOriginOrDesti: string "origin" or "destination"
					sum (integrate) over origin or destination.
		
		Returns
		-------
		No return. Results are stored as bov file.
		"""
		excitation -= 1 #coefficients start from 0
		if(OriginLabels == [-1]):
			#OriginLabels = numpy.arange(len(self.npart)).tolist()
			OriginLabels = numpy.arange(self.npart).tolist()
		mask1 = ((self.par) == OriginLabels[0])
		for i in OriginLabels:
			mask1 += ((self.par) == i)
		if(DestiLabels == [-1]):
			#DestiLabels = numpy.arange(len(self.npart)).tolist()
			DestiLabels = numpy.arange(self.npart).tolist()
		mask2 = ((self.par) == DestiLabels[0])
		for i in DestiLabels:
			mask2 += ((self.par) == i)
		orbdens1 = CalDensResp(self.nks, self.wf, self.coeffs[excitation])
		res = numpy.zeros((self.ngrid[0], self.ngrid[1], self.ngrid[2]))
		for i in range(0, self.nks):
			print(i, end="")
			if(len(self.wf[i]) == 0):
				continue
			else:
				orbdens0 = self.wf[i]**2
				if(SumOriginOrDesti == "origin"):
					res += numpy.sum(mask1*orbdens0)*(mask2*orbdens1[i,:,:,:])*self.vol
				else:
					res += (mask1*orbdens0)*numpy.sum(mask2*orbdens1[i,:,:,:])*self.vol
		Savebov(self.r0, self.dr, res, "PHMForApoint" + str(self.flag_phmforapoint))
		self.flag_phmforapoint += 1
		return

gs = GaussianPostProcess("paper.cub", "dimer.log", 128, "geo-pap.xyz", first_orb_in_file=1, cutoff = 0.001)
gs.UpdatePartition('atomic')
gs.CalculatePHM()
gs.PHMForApoint(2)
