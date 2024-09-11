import h5py

h5_name = "39.h5"

f_h5 = h5py.File(h5_name, "r")
head = f_h5.attrs
nx, ny, nz = head["dims"][0], head["dims"][1], head["dims"][2]
dx, dy, dz = head["dx"][0], head["dx"][0], head["dx"][0]
f_h5.close()

with open("0.xmf", "w") as f_xmf:
    # Open the file and write the XML description of the mesh.
    f_xmf.write("<?xml version=\"1.0\" ?>\n")
    f_xmf.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
    f_xmf.write("<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n")
    f_xmf.write(" <Domain>\n")
    f_xmf.write("   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n")
    f_xmf.write(f"     <Topology TopologyType=\"3DCORECTMesh\" NumberOfElements=\"{nx+1} {ny+1} {nz+1}\"/>\n")
    f_xmf.write("     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n")
    f_xmf.write(f"       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n")
    f_xmf.write(f"        0 0 0\n")
    f_xmf.write("       </DataItem>\n")
    f_xmf.write(f"       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n")
    f_xmf.write(f"        {dx} {dy} {dz}\n")
    f_xmf.write("       </DataItem>\n")
    f_xmf.write("     </Geometry>\n")
    f_xmf.write("     <Attribute Name=\"density\" AttributeType=\"Scalar\" Center=\"Cell\">\n")
    f_xmf.write(f"       <DataItem Dimensions=\"{nx} {ny} {nz}\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n")
    f_xmf.write(f"        {h5_name}:/density\n")
    f_xmf.write("       </DataItem>\n")
    f_xmf.write("     </Attribute>\n")
    f_xmf.write("   </Grid>\n")
    f_xmf.write(" </Domain>\n")
    f_xmf.write("</Xdmf>\n")
    f_xmf.close()