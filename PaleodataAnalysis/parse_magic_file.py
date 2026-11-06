import os
import pandas as pd
import numpy as np
import arcpy
# Python.exe needs to be in ArcGIS:  C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe

def parse_magic_file(filepath):
    """Parses a multi-section MagIC v3 tab-delimited file into a dictionary of pandas DataFrames."""
    tables = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    sections = content.split('>>>>>>>>>>')
    for section in sections:
        lines = [line.strip() for line in section.strip().split('\n') if line.strip()]
        if not lines:
            continue
        
        table_name = None
        data_start_idx = 0
        for i, line in enumerate(lines):
            if line.startswith('tab delimited'):
                table_name = line.split('tab delimited')[1].strip()
                data_start_idx = i + 1
                break
        
        if table_name and data_start_idx < len(lines):
            header = lines[data_start_idx].split('\t')
            rows = [line.split('\t') for line in lines[data_start_idx+1:]]
            valid_rows = [r for r in rows if len(r) == len(header)]
            if valid_rows:
                tables[table_name] = pd.DataFrame(valid_rows, columns=header)
    return tables

# ==========================================
# 1. PARSE & CLEAN DATA
# ==========================================
# Ensure this matches the exact name of your downloaded MagIC text file
file_path = "magic_contribution_18469.txt" 

magic_tables = parse_magic_file(file_path)

if 'locations' in magic_tables:
    df = magic_tables['locations']
    
    # Cast numerical fields
    df['lon_w'] = pd.to_numeric(df['lon_w'], errors='coerce')
    df['lat_s'] = pd.to_numeric(df['lat_s'], errors='coerce')
    df['dir_dec'] = pd.to_numeric(df['dir_dec'], errors='coerce')
    df['dir_inc'] = pd.to_numeric(df['dir_inc'], errors='coerce')

    # Convert MagIC 0-360 longitude to standard -180 to 180 for ArcGIS
    df['x_lon'] = df['lon_w'].apply(lambda lon: lon - 360 if lon > 180 else lon)
    df['y_lat'] = df['lat_s']

    # Calculate Tectonic Rotation Angle relative to stable North America (ref_dec = 0.0)
    ref_dec = 0.0 
    df['rotation_angle'] = df['dir_dec'].apply(
        lambda d: (d - ref_dec + 180) % 360 - 180 if pd.notnull(d) else np.nan
    )

    # Write cleaned CSV using an absolute path
    cleaned_csv = os.path.abspath("cleaned_wells_coe_data.csv")
    df.to_csv(cleaned_csv, index=False)
    print(f"Successfully processed {len(df)} records. Saved to: {cleaned_csv}")

    # ==========================================
    #  INGEST INTO ARCGIS PRO
    # ==========================================

    # 1. Ensure ArcGIS can overwrite an existing feature class if you run this multiple times
    arcpy.env.overwriteOutput = True

    # 2. Drop any rows where latitude or longitude ended up as NaN/null (prevents XYTableToPoint crashes)
    df = df.dropna(subset=['x_lon', 'y_lat'])

    # Re-save the cleaned CSV so it only contains valid spatial rows
    cleaned_csv = os.path.abspath("cleaned_wells_coe_data.csv")
    df.to_csv(cleaned_csv, index=False)

    # 3. Handle workspace paths safely
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        default_gdb = aprx.defaultGeodatabase
        output_feature_class = os.path.join(default_gdb, "Wells_Coe_Tectonic_Rotations")
        print("Targeting active ArcGIS Pro project geodatabase.")
    except Exception:
        output_feature_class = "Wells_Coe_Tectonic_Rotations"
        print("Targeting default workspace environment.")

    # 4. Run XYTableToPoint with safety checks in place

    arcpy.env.overwriteOutput = True

    # Drop any rows with missing coordinates
    df = df.dropna(subset=['x_lon', 'y_lat'])
    cleaned_csv = os.path.abspath("cleaned_wells_coe_data.csv")
    df.to_csv(cleaned_csv, index=False)

    # Explicitly define a file geodatabase in your script's directory so it always has a valid location
    script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
    gdb_name = "PaleomagData.gdb"
    gdb_path = os.path.join(script_dir, gdb_name)

    # Create the file geodatabase if it doesn't already exist
    if not arcpy.Exists(gdb_path):
        arcpy.management.CreateFileGDB(script_dir, gdb_name)
        print(f"Created new file geodatabase at: {gdb_path}")

    output_feature_class = os.path.join(gdb_path, "Wells_Coe_Tectonic_Rotations")

    # Run XYTableToPoint with the fully qualified path
    try:
        arcpy.management.XYTableToPoint(
            in_table=cleaned_csv,
            out_feature_class=output_feature_class,
            x_field="x_lon",
            y_field="y_lat",
            coordinate_system=arcpy.SpatialReference(4326) # WGS 1984
        )
        print(f"\nSuccess! Feature class created inside geodatabase:")
        print(output_feature_class)
    except Exception as e:
        print(f"\nArcPy Error encountered: {e}")

