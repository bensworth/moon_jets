##############################################################
### File to post-process individual HDF5 files saved from C++
### simulations and combine into a single HDF5 file, including
### all simulation data as individual datasets in a group. The
### dataset name includes indices representing the particle
### radius and initial speed, and correspong to the datasets
### "radii" and "speeds".
##############################################################
import pdb
import os
import h5py
import numpy as np

# NOTE: HDF5 groups work like dictionaries, and datasets work like NumPy arrays

# Parameters specifying which data to parse
angular_dist = "uniform"
max_angle_load = 15

# Loop over all hdf5 files in the data directory. In this loop,
# we are simply collecting unique particle radii and speeds to
# enumerate, and making sure attributes are stored
speed_set = set()
radii_set = set()
flux_data_dims = np.array([-1,-1])
speed_data_dims = np.array([-1,-1])
size_data_dims = np.array([-1,-1])
altitudes = None
inclinations = None
if angular_dist == "uniform":
    prefix = "./raw_data_uni/"
else:
    prefix = "./raw_data_cos2/"
for filename in os.listdir(prefix):
    if filename.endswith(str(max_angle_load) + ".hdf5"):

        # Load HDF5 file, get speed and radius attributes
        f = h5py.File(prefix + filename, 'r')
        attrs = [attr for attr in f.attrs.keys()]   # Attribute dictionary keys
        if "maximum angle" in attrs:
            max_ang = int(f.attrs["maximum angle"][0])
            if max_angle_load != max_ang:
                print("Loaded max angle does not match filename!")
                pdb.set_trace()
                raise ValueError("Bad data.")
        else:
            print("Warning! Dataset without radius attribute!")
            pdb.set_trace()
            raise ValueError("Incomplete data.")
        if "particle radius" in attrs:
            radius = f.attrs["particle radius"][0]
        else:
            print("Warning! Dataset without radius attribute!")
            pdb.set_trace()
            raise ValueError("Incomplete data.")
        if "initial speed" in attrs:
            speed = f.attrs["initial speed"][0]
        else:
            print("Warning! Dataset without speed attribute!")
            pdb.set_trace()
            raise ValueError("Incomplete data.")
        if "uniform angular distribution" in attrs:
            if angular_dist is None:
                angular_dist = "uniform"
            elif angular_dist != "uniform":
                raise ValueError("Inconsistent angular distributions.")
        else:
            angular_dist = "cos^2"

        # Store speed and radius attributes
        speed_set.add(speed)
        radii_set.add(radius)
        
        # Check altitudes are same across data sets
        if altitudes is None:
            try: 
                altitudes = np.array(f['altitudes'])
            except:
                print("Warning! Dataset without altitudes key!")
                pdb.set_trace()
                raise ValueError("Incomplete data.")
        else:
            try: 
                test = np.max(np.abs(altitudes - np.array(f['altitudes'])))
                if test > 1e-5:
                    print("Warning! Inconsistent altitudes across files!")
                    pdb.set_trace()
                    raise ValueError("Inconsistent data.")
            except:
                print("Warning! Dataset without altitudes key!")
                pdb.set_trace()
                raise ValueError("Incomplete data.")

        # Check inclinations are same across data sets
        if inclinations is None:
            try: 
                inclinations = np.array(f['inclinations'])
            except:
                print("Warning! Dataset without inclinations key!")
                pdb.set_trace()
                raise ValueError("Incomplete data.")
        else:
            try: 
                test = np.max(np.abs(inclinations - np.array(f['inclinations'])))
                if test > 1e-5:
                    print("Warning! Inconsistent inclinations across files!")
                    pdb.set_trace()
                    raise ValueError("Inconsistent data.")
            except:
                print("Warning! Dataset without inclinations key!")
                pdb.set_trace()
                raise ValueError("Incomplete data.")

        f.close()

# Check to make sure files were loaded, if not quit here
if altitudes is None:
    raise ValueError("No files loaded!")

# Convert altitudes to m
altitudes *= 1000.0

# Convert speeds to m/s
speed_arr = np.array(np.sort(list(speed_set)))
speed_arr *= 1000.0

# Convert radii to m
radii_arr = np.array(np.sort(list(radii_set)))
radii_arr *= 1e-6

# New HDF5 files that will aggregate all others as groups
# - Use one for flux data and one for distribution data
if angular_dist == "uniform":
    hh = h5py.File('EncFluxData_max-angle'+str(max_angle_load)+'_uni.hdf5', 'w')
    ha = h5py.File('EncDistData_max-angle'+str(max_angle_load)+'_uni.hdf5', 'w')
else:
    hh = h5py.File('EncFluxData_max-angle'+str(max_angle_load)+'.hdf5', 'w')
    ha = h5py.File('EncDistData_max-angle'+str(max_angle_load)+'.hdf5', 'w')

# Store set of all speeds and radii as sorted numpy arrays;
# add as datasets to complete HDF5 file
hh.create_dataset("grid altitudes", data=altitudes)
hh.create_dataset("grid inclinations", data=inclinations)
hh.create_dataset("speeds m per s", data=speed_arr)
hh.create_dataset("radii m", data=radii_arr)
hh.attrs["num altitudes"] = altitudes.shape[0]
hh.attrs["num inclinations"] = inclinations.shape[0]
hh.attrs["num speeds"] = speed_arr.shape[0]
hh.attrs["num radii"] = radii_arr.shape[0]
if angular_dist == "uniform":
    hh.attrs["angular dist"] = "uniform"
else:
    hh.attrs["angular dist"] = "cos^2"
ha.create_dataset("grid altitudes", data=altitudes)
ha.create_dataset("grid inclinations", data=inclinations)
ha.create_dataset("speeds m per s", data=speed_arr)
ha.create_dataset("radii m", data=radii_arr)
ha.attrs["num altitudes"] = altitudes.shape[0]
ha.attrs["num inclinations"] = inclinations.shape[0]
ha.attrs["num speeds"] = speed_arr.shape[0]
ha.attrs["num radii"] = radii_arr.shape[0]
if angular_dist == "uniform":
    ha.attrs["angular dist"] = "uniform"
else:
    ha.attrs["angular dist"] = "cos^2"

alt_dists = np.zeros((altitudes.shape[0], speed_arr.shape[0], radii_arr.shape[0]))

# Loop over all hdf5 files in the data directory and copy flux 
# profiles to a single HDF5 file 
hh.create_group("flux profiles")
for filename in os.listdir(prefix):
    if filename.endswith(str(max_angle_load) + ".hdf5"):

        # Load HDF5 file, get speed and radius attributes
        f = h5py.File(prefix + filename, 'r')
        radius = f.attrs["particle radius"][0]*1e-6
        speed = f.attrs["initial speed"][0]*1000.0

        # Find index of this files speed and radius
        tempr = np.where(radii_arr == radius)[0][0]
        temps = np.where(speed_arr == speed)[0][0]

        # Copy datasets and attributes from f to a group in hh
        f.copy(f["flux"], hh["flux profiles"], "flux_r" + str(tempr) + "_s" + str(temps))

        # Set a;titude distribution for given speed and size 
        alt_dists[:,tempr,temps] = f["altitude distribution"]

        f.close()
    else:
        continue

# Create group of altitude distributions, where the ith
# distribution is a 2d array of distribution values in
# speed and size at grid altitudes[i]. Each one is labeled
# as, e.g., a5 for the 5th altitude in grid altitudes
# attribute. 
ha.create_group("altitude distributions")
for i in range(0,altitudes.shape[0]):
    f.copy(alt_dists[i,:,:], ha["altitude distributions"], "a" + str(i))

hh.close()
ha.close()