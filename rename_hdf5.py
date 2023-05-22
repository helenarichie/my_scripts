import os
import re

collection = "/Users/helenarichie/Desktop/flatiron/2023-04-26/"

convert = lambda text: int(text) if text.isdigit() else text.lower()
alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

for i, filename in enumerate(sorted(os.listdir(collection), key=alphanum_key)):
    if filename.endswith(".png"):
        #print(filename)
        #print(str(i) + "_proj.png")
        os.rename(os.path.join(collection, filename), os.path.join(collection, str(i) + "_proj.png"))