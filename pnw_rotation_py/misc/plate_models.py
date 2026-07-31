import json
import gplately
import pygplates

# All credit for this amazing work goes to Google Gemini

def track_and_export_models():
    # 1. Fetching the GPlately data models
    print("Fetching Müller et al. (2019) model...")
    server_2019 = gplately.DataServer("Muller2019")
    rotation_model_2019, _, _ = server_2019.get_plate_reconstruction_files()
    
    print("Fetching Matthews et al. (2016) model...")
    server_2016 = gplately.DataServer("Matthews2016")
    rotation_model_2016, _, _ = server_2016.get_plate_reconstruction_files()
    
    print("Fetching Seton et al. (2012) model...\n")
    server_2012 = gplately.DataServer("Seton2012")
    rotation_model_2012, _, _ = server_2012.get_plate_reconstruction_files()
    
    # 2. Setup your specific Yellowstone Hotspot (YHS) parameters
    yhs_lat = 44.43
    yhs_long = -110.67
    plate_id = 101       # North American Plate
    anchor_plate_id = 0  # Global spin axis / deep mantle frame
    
    yhs_point = pygplates.PointOnSphere(yhs_lat, yhs_long)
    timeline_ma = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
    
    models = {
        "Muller2019": {"obj": rotation_model_2019, "name": "Müller et al. (2019)", "features": []},
        "Matthews2016": {"obj": rotation_model_2016, "name": "Matthews et al. (2016)", "features": []},
        "Seton2012": {"obj": rotation_model_2012, "name": "Seton et al. (2012)", "features": []}
    }
    
    print(f"=========================================================================================")
    print(f" KINEMATIC & YHS PATH ANALYSIS FOR NORTH AMERICAN PLATE (ID {plate_id})")
    print(f"=========================================================================================")
    
    for time_ma in timeline_ma:
        print(f"\nReconstruction Time: {time_ma} Ma")
        print("-" * 115)
        
        for key, model_data in models.items():
            try:
                # Query finite rotation for Euler Pole
                finite_rot = model_data["obj"].get_rotation(time_ma, plate_id, anchor_plate_id)
                p_lat, p_lon, p_ang = finite_rot.get_lat_lon_euler_pole_and_angle_degrees()
                
                # Reconstruct the YHS crustal point backward
                reconstructed_point = finite_rot.get_inverse() * yhs_point
                y_lat, y_lon = reconstructed_point.to_lat_lon()
                
                # Print merged terminal view
                print(f"  {model_data['name']:<23} -> POLE [Lat: {p_lat:>8.4f}° | Lon: {p_lon:>9.4f}° | Ang: {p_ang:>8.4f}°] | YHS [Paleolat: {y_lat:>8.4f}° | Paleolon: {y_lon:>9.4f}°]")
                
                # Store structural GeoJSON Feature data attributes
                geojson_feature = {
                    "type": "Feature",
                    "geometry": {
                        "type": "Point",
                        "coordinates": [float(y_lon), float(y_lat)] # GeoJSON requires [Longitude, Latitude]
                    },
                    "properties": {
                        "Time_Ma": float(time_ma),
                        "Model": model_data["name"],
                        "Euler_Lat": float(p_lat),
                        "Euler_Lon": float(p_lon),
                        "Euler_Ang": float(p_ang)
                    }
                }
                model_data["features"].append(geojson_feature)
                
            except Exception as e:
                print(f"  {model_data['name']:<23} -> Error: {e}")
                
    # 3. Export collected tracking items to individual GeoJSON files
    print(f"\n=========================================================================================")
    print(" EXPORTING GEOSPATIAL VECTOR FILES FOR QGIS")
    print(f"=========================================================================================")
    for key, model_data in models.items():
        filename = f"yhs_track_{key}.geojson"
        geojson_output = {
            "type": "FeatureCollection",
            "features": model_data["features"]
        }
        with open(filename, "w") as f:
            json.dump(geojson_output, f, indent=4)
        print(f" Successfully exported: {filename}")

if __name__ == "__main__":
    track_and_export_models()
