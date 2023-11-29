from hconfig import *
from csv import writer

date = "2023-10-08"
new = True

basedir = f"/ix/eschneider/helena/data/cloud_wind/{date}/" # crc
logdir = os.path.join(basedir)

reference_frame_shift = []
with open(os.path.join(logdir, "output.log")) as f:
    for i, line in enumerate(f):
        if line.startswith("Average cloud velocity"):
            line = line.split(" ")
            reference_frame_shift.append(float(line[-2]))

if new:
    f = open(os.path.join(logdir, "frame_shift.csv"), "w")
    f.close()

with open(os.path.join(logdir, "frame_shift.csv"), "a") as f:
    writer_obj = writer(f)
    for i, v in enumerate(reference_frame_shift):
        writer_obj.writerow([v])
    f.close()