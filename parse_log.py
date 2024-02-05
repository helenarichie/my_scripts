from hconfig import *
from csv import writer

################################################################
date = "2024-02-04"
################################################################

basedir = f"/ix/eschneider/helena/data/testing/{date}/"
csvdir = os.path.join(basedir, "csv/")

f = open(os.path.join(csvdir, "sputter_hot.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "sputter.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "mass_cloud.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "mass_dust.csv"), "w")
f.close()
f = open(os.path.join(csvdir, "time.csv"), "w")
f.close()

sputter_hot, sputter = [], []
total_sputter = [sputter_hot, sputter]
strings_sputter = ["Mixed sputtered mass: ", "Hot sputtered mass: "]
mass_cloud, mass_dust = [], []
total_mass = [mass_cloud, mass_dust]
strings_mass = ["Cloud mass: ", "Dust mass: "]
times = []

i = 0
l = 0
m = 0

with open(os.path.join(basedir, "output.log")) as f:
    for line in f:
        if line.startswith("** "):
            line = line.lstrip("** ").rstrip(" \n")
            line = line.split("  ")
            for j, entry in enumerate(line):
                entry = entry.lstrip(strings_sputter[j])
                total_sputter[j].append(float(entry))
                if j == 0:
                    with open(os.path.join(csvdir, "sputter_hot.csv"), "a") as f_txt:
                        writer_obj = writer(f_txt)
                        writer_obj.writerow([total_sputter[j][i]])
                        f_txt.close()
                if j == 1:
                    with open(os.path.join(csvdir, "sputter.csv"), "a") as f_txt:
                        writer_obj = writer(f_txt)
                        writer_obj.writerow([total_sputter[j][i]])
                        f_txt.close()
            i += 1

        elif line.startswith("@@ "):
            line = line.lstrip("@@ ").rstrip(" \n")
            line = line.split("  ")
            for j, entry in enumerate(line):
                entry = entry.lstrip(strings_mass[j])
                total_mass[j].append(float(entry))
                if j == 0:
                    with open(os.path.join(csvdir, "mass_cloud.csv"), "a") as f_txt:
                        writer_obj = writer(f_txt)
                        writer_obj.writerow([total_mass[j][l]])
                        f_txt.close()
                if j == 1:
                    with open(os.path.join(csvdir, "mass_dust.csv"), "a") as f_txt:
                        writer_obj = writer(f_txt)
                        writer_obj.writerow([total_mass[j][l]])
                        f_txt.close()
            l += 1

        elif line.startswith("n_step"):
            line = line.split(" ")
            times.append(float(line[6]))
            with open(os.path.join(csvdir, "time.csv"), "a") as f_txt:
                writer_obj = writer(f_txt)
                writer_obj.writerow([times[m]])
                f_txt.close()

            m += 1