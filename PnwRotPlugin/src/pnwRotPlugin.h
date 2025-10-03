#ifndef _QGIS_pnwRotationPlugin_WORLD_H_
#define _QGIS_pnwRotationPlugin_WORLD_H_

/* MSVC workarounds */
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923132169163975144   // pi/2
#endif

#include "qgisplugin.h"
#include "qgsmapcanvas.h"
#include "qgisinterface.h"
#include "qgsvectorlayer.h"
#include "qgsmessagelog.h"
#include <iostream>
#include <QAction>
#include <QApplication>
#include "qgsVectorDataProvider.h"
#include "qgssinglesymbolrenderer.h"
#include "qgsmaplayeractionregistry.h"
#include "qgssymbol.h."
#include <QVariant>
#include <qgslogger.h> // For logging potential errors


class pnwRotationPlugin : public QObject, public QgisPlugin
{
   Q_OBJECT

public:
   /// @brief Constructor.
   /// @param qgis_if The Qgis interface.
   explicit pnwRotationPlugin(QgisInterface* qgis_if);

   /// @brief Destructor
   virtual ~pnwRotationPlugin() = default;

   /// @brief Called when the plugin is loaded.
   virtual void initGui() override;

   /// @brief Called when the plugin is unloaded.
   virtual void unload() override;

public slots:
   void clear_menu_button_action();
   void rot_menu_button_action();
   void yhs_menu_button_action();

private:
   QgisInterface* m_qgis_if;

   /// The actions in the QGIS menu bar.
   QAction *m_clear_menu_action;
   QAction *m_display_rot_menu_action;
   QAction *m_yhs_menu_action;

   QgsVectorLayer *m_rotSrcLayer = NULL;
   QgsVectorLayer *m_rotDestLayer = NULL;
   QgsVectorLayer *m_yhsDestLayer = NULL;

   QList<QgsField> m_fieldList;   
   QgsFeatureList m_rotFeatureList;
   QgsFeatureList m_yhsFeatureList;
   std::vector<QString> m_fieldNames;
   std::vector<std::vector<double>> m_rot_data;
   QgsPolylineXY m_line;

   
   bool m_verbose = true;
   bool m_layers_setup = false;
   bool m_rotDataLoaded = false;
   bool m_yhsDataLoaded = false;
   double m_lat, m_lon;

   // All valocity units match the Zeng data: mm/yr 
   // YHS current center
   const double EARTH_RADIUS = 6371000; // meters
   const double YHS_lat = 44.43;
   const double YHS_lon = -110.67;
   // YHS NA Plate Velocity: WSW @ 46 mm/yr is best estimate (WSW ~247.5 degrees)
   const double NA_Vel_N = -0.3826 * 46; // cos(247.5) * 46 mm/Y
   const double NA_Vel_E = -0.9239 * 46; // sin(247.5) * 46 mm/Y
   const double detlaT = 1E6; // 1 million year intervals

   bool setupLayers();
   bool loadRotData();
   bool setupRotLayer();
   bool setupYhsLayer();
   void displayRotData(QgsFeatureList& featureList);
   void displayYhsData(QgsPolylineXY& list);
   double getFeatureAttrubute(QgsFeature &feature, int index);
   bool setFeatureAttribute(QgsFeature &feature, int index, double value);
   QgsPointXY getClosestRotEntry(QgsPointXY endPoint);
   void clear_display_data();

   // Meters N,E to lat, Lon
   double latitudeFromDisatnce(double d);
   double longitudeFromDistance(double latitude, double d);
};

#endif
