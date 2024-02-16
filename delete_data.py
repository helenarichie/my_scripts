import os
import re
import glob

delete = False
rename = True

datadir = "/ix/eschneider/helena/data/cloud_wind/2024-02-07/hdf5/"

convert = lambda text: int(text) if text.isdigit() else text.lower()
alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

if delete:
    for i, filename in enumerate(sorted(glob.glob(os.path.join(datadir, "*.h5")), key=alphanum_key)):
        # print(i)
        # print(filename)
        if i % 2 != 0:
            os.remove(os.path.join(datadir, filename))
            continue

if rename: 
    for i, filename in enumerate(sorted(glob.glob(os.path.join(datadir, "*.h5")), key=alphanum_key)):
            os.rename(os.path.join(datadir, filename), os.path.join(datadir, str(i) + ".h5"))