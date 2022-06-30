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
angular_dist = "cos2"
max_angle_load = 25

# Loop over all hdf5 files in the data directory. In this loop,
# we are simply collecting unique particle radii and speeds to
# enumerate, and making sure attributes are stored
speed_set = set()
flux_data_dims = np.array([-1,-1])
speed_data_dims = np.array([-1,-1])
size_data_dims = np.array([-1,-1])
altitudes = None
inclinations = None
prefix = "./gridtype1/"  # For large data in gridtype1 folder


for filename in os.listdir(prefix):
    # if filename.endswith(str(max_angle_load) + ".hdf5"):
    if filename.endswith(".hdf5"):     # For files where I hadnt yet added angle to filename

        # Load HDF5 file, get speed and radius attributes
        f = h5py.File(prefix + filename, 'r')
        attrs = [attr for attr in f.attrs.keys()]   # Attribute dictionary keys
        # if "maximum angle" in attrs:
        #     max_ang = int(f.attrs["maximum angle"][0])
        #     if max_angle_load != max_ang:
        #         print("Loaded max angle does not match filename!")
        #         pdb.set_trace()
        #         raise ValueError("Bad data.")
        # else:
        #     print("Warning! Dataset without max angle attribute!")
        #     pdb.set_trace()
            # raise ValueError("Incomplete data.")
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

# New HDF5 files that will aggregate all others as groups
# - Use one for flux data and one for distribution data
if angular_dist == "uniform":
    ha = h5py.File('EncAltitudeData_max-angle'+str(max_angle_load)+'_uni.hdf5', 'w')
else:
    ha = h5py.File('EncAltitudeData_max-angle'+str(max_angle_load)+'.hdf5', 'w')

# Store set of all speeds and radii as sorted numpy arrays;
# add as datasets to complete HDF5 file
ha.create_dataset("grid altitudes", data=altitudes)
ha.create_dataset("grid inclinations", data=inclinations)
ha.create_dataset("speeds m per s", data=speed_arr)
ha.attrs["num altitudes"] = altitudes.shape[0]
ha.attrs["num inclinations"] = inclinations.shape[0]
ha.attrs["num speeds"] = speed_arr.shape[0]
if angular_dist == "uniform":
    ha.attrs["angular dist"] = "uniform"
else:
    ha.attrs["angular dist"] = "cos^2"

alt_dists = np.zeros((altitudes.shape[0], speed_arr.shape[0]))

# Loop over all hdf5 files in the data directory and copy flux 
# profiles to a single HDF5 file 
for filename in os.listdir(prefix):
    # if filename.endswith(str(max_angle_load) + ".hdf5"):
    if filename.endswith(".hdf5"):

        # Load HDF5 file, get speed and radius attributes
        f = h5py.File(prefix + filename, 'r')
        speed = f.attrs["initial speed"][0]*1000.0

        # Find index of this files speed and radius
        temps = np.where(speed_arr == speed)[0][0]
        print(temps)

        # Set altitude distribution for given speed and size 
        alt_dists[:,temps] = np.array(f["residence_time"])

        f.close()
    else:
        continue

ha.create_dataset("altitude densities", data=alt_dists)
ha.close()


