Attributes {
    Depth {
        String long_name "Depths of Standard Layer";
        String standard_name "depth";
        String units "m";
        String positive "down";
        String axis "Z";
    }
    Latitude {
        String long_name "Latitude in common grid";
        String standard_name "latitude";
        String units "degrees_north";
    }
    Longitude {
        String long_name "Longitude in common grid";
        String standard_name "longitude";
        String units "degrees_east";
    }
    h {
        String long_name "bathymetry";
        String standard_name "sea_floor_depth_below_mean_sea_level";
        String units "meter";
        String positive "down";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude";
        String field "bath, scaler";
        Float64 _FillValue -99999.0;
    }
    mask {
        String long_name "mask";
        String standard_name "sea_binary_mask";
        Float64 flag_values 0.0, 1.0;
        String flag_meanings "land, water";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude";
    }
    time {
        String long_name "time since initialization";
        String units "seconds since 2018-01-01 00:00:00";
        String calendar "gregorian";
        String field "time, scalar, series";
        UInt32 _ChunkSizes 512;
    }
    Times {
        String time_zone "UTC";
        UInt32 _ChunkSizes 1, 26;
        DODS {
            Int32 strlen 26;
            String dimName "DateStrLen";
        }
    }
    zeta {
        String long_name "free-surface referenced to mean sea level";
        String standard_name "sea_surface_elevation";
        String units "meters";
        String time "time";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude time";
        String field "free-surface, scalar, series";
        Float32 _FillValue -99999.0;
        UInt32 _ChunkSizes 1, 777, 760;
    }
    zetatomllw {
        String long_name "free-surface referenced to mean lower low water";
        String units "meters";
        String time "time";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude time";
        String field "free-surface, scalar, series";
        Float32 _FillValue -99999.0;
        UInt32 _ChunkSizes 1, 777, 760;
    }
    u_eastward {
        String long_name "eastward momentum component";
        String standard_name "eastward_sea_water_velocity";
        String units "meters second-1";
        String time "time";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude Depth time";
        String field "u_eastward, scalar, series";
        Float32 _FillValue -99999.0;
        UInt32 _ChunkSizes 1, 7, 311, 304;
    }
    v_northward {
        String long_name "northward momentum component";
        String standard_name "northward_sea_water_velocity";
        String units "meters second-1";
        String time "time";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude Depth time";
        String field "v_northward, scalar, series";
        Float32 _FillValue -99999.0;
        UInt32 _ChunkSizes 1, 7, 311, 304;
    }
    temp {
        String long_name "potential temperature";
        String standard_name "sea_water_temperature";
        String units "Celsius";
        String time "time";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude Depth time";
        String field "temperature, scalar, series";
        Float32 _FillValue -99999.0;
        UInt32 _ChunkSizes 1, 7, 311, 304;
    }
    salt {
        String long_name "salinity";
        String standard_name "sea_water_salinity";
        String units "PSU";
        String time "time";
        String grid "grid";
        String location "face";
        String coordinates "Longitude Latitude Depth time";
        String field "salinity, scalar, series";
        Float32 _FillValue -99999.0;
        UInt32 _ChunkSizes 1, 7, 311, 304;
    }
    NC_GLOBAL {
        String title "SSCOFS";
        String institution "School for Marine Science and Technology";
        String source "FVCOM_4.4.7";
        String history "model started at: 01/10/2026   16:04";
        String references "http://fvcom.smast.umassd.edu, https://github.com/FVCOM-GitHub/FVCOM";
        String Conventions "CF-1.0";
        String CoordinateSystem "GeoReferenced";
        String CoordinateProjection "proj=utm +ellps=WGS84 +zone=10";
        String Tidal_Forcing "Tidal Forcing Time Series Title: sscofs lateral open boundary netCDF file";
        String River_Forcing "THERE ARE 38 RIVERS IN THIS MODEL.
RIVER INFLOW IS ON THE nodes WHERE TEMPERATURE AND SALINITY ARE specified IN THE MODEL.
THE FOLLOWING RIVER NAMES ARE USED:
14144700
14144700
08MF005
08MF005
12200500
12200500
12167000
12167000
12150800
12150800
14211720
14211720
12213100
12213100
12201500
12201500
12113000
12113000
12089500
12089500
12101500
12101500
12061500
12061500
12080010
12080010
12031000
12031000
12045500
12045500
12048000
12048000
12054000
12054000
14220500
14220500
14243000
14243000";
        String GroundWater_Forcing "GROUND WATER FORCING IS OFF!";
        String Surface_Heat_Forcing "FVCOM variable surface heat forcing file:
FILE NAME:sscofs.t15z.20261001.hflux.nowcast.nc
SOURCE:FVCOM grid (unstructured) surface forcing
Unknown start date meta data format";
        String Surface_Wind_Forcing "FVCOM variable surface Wind forcing:
FILE NAME:sscofs.t15z.20261001.met.nowcast.nc
SOURCE:FVCOM grid (unstructured) surface forcing
Unknown start date meta data format";
        String Surface_PrecipEvap_Forcing "SURFACE PRECIPITATION FORCING IS OFF";
    }
    DODS_EXTRA {
        String Unlimited_Dimension "time";
    }
}
