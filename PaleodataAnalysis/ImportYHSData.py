import os
import arcpy

# Run XYTableToPoint with safety checks in place

# Explicitly define a file geodatabase in your script's directory so it always has a valid location
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
gdb_name = "PaleomagData.gdb"
gdb_path = os.path.join(script_dir, gdb_name)

# Create the file geodatabase if it doesn't already exist
if not arcpy.Exists(gdb_path):
    arcpy.management.CreateFileGDB(script_dir, gdb_name)
    print(f"Created new file geodatabase at: {gdb_path}")

output_feature_class = os.path.join(gdb_path, "YhsPathFeatureClass")
    
try:
    yhs_data_name = "YHS_path.csv"
    yhs_data_path = os.path.join(script_dir, yhs_data_name)
    arcpy.management.XYTableToPoint(
        in_table= yhs_data_path,
        out_feature_class=output_feature_class,
        x_field="long",
        y_field="lat",
        coordinate_system=arcpy.SpatialReference(4326) # WGS 1984
    )

    print(f"\nSuccess! YHS Feature class created inside geodatabase:")
    print(output_feature_class)

except Exception as e:
    print(f"\nArcPy Error encountered: {e}")