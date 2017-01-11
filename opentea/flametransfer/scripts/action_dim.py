from XDR2 import Dataset

ds = Dataset("dataset.xml")

dim_list=["x","y"]
if ds.getValue("dim_choice") == "threed":
    dim_list.append("z")
ds.setValue(dim_list,"dim_list")

ds.save2file("out_dataset.xml")
