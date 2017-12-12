def GaussianCubeWFLoader(file='dft-canon.cube', skip_first = 0, skip_last = 0, first_orb_in_file = 1):
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
        natoms = int(line.split()[0])
        print(natoms)
        r0 = []
        r0.append(float(line.split()[1]))
        r0.append(float(line.split()[2]))
        r0.append(float(line.split()[3]))
        print(r0)
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
        print(ngrid)
        print(dr)
        for iatom in range(0, natoms):
                f.readline()
        line = f.readline()
        print(line)
        nwf = int(line.split()[0])
        print(nwf)

GaussianCubeWFLoader()
