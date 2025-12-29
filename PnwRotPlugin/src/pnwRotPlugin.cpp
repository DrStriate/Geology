#pragma warning(disable : 4996)
#define sqr(x) ((x) * (x))
#include "pnwRotPlugin.h"
#include <cmath>
#include <qgslinesymbol.h>
#include <QtWidgets>
#include <qaction.h>
#include <cstdio>

namespace
{
   const QString s_name = QStringLiteral("Pnw Rotation Analysis");
   const QString s_description = QStringLiteral("Pnw rotation analysis plugin");
   const QString s_category = QStringLiteral("Plugins");
   const QString s_version = QStringLiteral("Version 1.2.3");
   const QString s_icon = QStringLiteral(":/plugin.svg");
   const QString s_rotDatadestLayerName = QStringLiteral("nshm2023_GPS_velocity");
   const QgisPlugin::PluginType s_type = QgisPlugin::UI;
   const QString s_yhsDestLayerName = "YHS movement";
   const QString s_rotDestLayerName = "PNW rotation";
}

QGISEXTERN QgisPlugin *classFactory(QgisInterface *qgis_if)
{
   std::cout << "::classFactory" << std::endl;
   return new pnwRotationPlugin(qgis_if);
}

QGISEXTERN const QString *name()
{
   return &s_name;
}

QGISEXTERN const QString *description()
{
   return &s_description;
}

QGISEXTERN const QString *category()
{
   return &s_category;
}

QGISEXTERN const QString *version()
{
   return &s_version;
}

QGISEXTERN const QString *icon()
{
   return &s_icon;
}

QGISEXTERN int type()
{
   return s_type;
}

QGISEXTERN void unload(QgisPlugin *plugin)
{
   std::cout << "::unload" << std::endl;
   delete plugin;
}

pnwRotationPlugin::pnwRotationPlugin(QgisInterface *iface) : QgisPlugin(s_name, s_description, s_category, s_version, s_type), m_qgis_if(iface)
{
   m_qgis_if = iface;

   // Set initial state point
   m_NA_Vel_N = cos(NA_Bearing / 180.0 * M_PI) * NA_Speed;  // cos(247.5) * 46 mm/Y
   m_NA_Vel_E = sin(NA_Bearing / 180.0 * M_PI) * NA_Speed; // sin(247.5) * 46 mm/Y
   const pState yhsState0 = {YHS_lon, YHS_lat, m_NA_Vel_E, m_NA_Vel_N};
   m_pYhsState.push_back(yhsState0);

   m_passes = 0;
}

void pnwRotationPlugin::unload()
{
   // TODO - need to remove the actions from the menu again.
   // Get the QgsProject instance
   QgsProject *project = QgsProject::instance();

   // Find the layer by name
   // mapLayersByName() returns a list, even if there is only one match
   QList<QgsMapLayer *> foundLayers = project->mapLayersByName(s_yhsDestLayerName);

   // Check if a layer was found
   if (!foundLayers.isEmpty())
   {
      QgsMapLayer *layerToDelete = foundLayers.first();
      QString layerId = layerToDelete->id();

      // Remove the layer from the project.
      // The QgsProject takes care of memory management after removal.
      project->removeMapLayer(layerId);

      // After this, the layer pointer (`layerToDelete`) is no longer valid.
      // It is a dangling pointer, so do not use it again.
   }
   else
   {
      qDebug() << "Layer not found:" << s_yhsDestLayerName;
   }
}

void pnwRotationPlugin::initGui()
{
   QgsMessageLog::logMessage(QString("initGui"), name(), Qgis::MessageLevel::Info);

   if (m_qgis_if == nullptr)
   {
      QgsMessageLog::logMessage("failed to get the handle to QGIS", name(), Qgis::MessageLevel::Info);
      return;
   }

   // add display rot data action to the menu
   m_display_rot_menu_action = new QAction(QIcon(""), QString("Display rot data"), this);
   connect(m_display_rot_menu_action, SIGNAL(triggered()), this, SLOT(rot_menu_button_action()));
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_display_rot_menu_action);

   // add clear action to the menu
   m_clear_menu_action = new QAction(QIcon(""), QString("Clear data"), this);
   connect(m_clear_menu_action, SIGNAL(triggered()), this, SLOT(clear_menu_button_action()));
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_clear_menu_action);

   // add YHS action to the menu
   m_yhs_menu_action = new QAction(QIcon(""), QString("Move YHS"), this);
   connect(m_yhs_menu_action, SIGNAL(triggered()), this, SLOT(yhs_menu_button_action()));
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_yhs_menu_action);
}

bool pnwRotationPlugin::setupLayers()
{
   if (m_layers_setup)
      return true;

   if (!loadRotData())
   {
      QgsMessageLog::logMessage("failed to load rot data from layer", name(), Qgis::MessageLevel::Info);
      return false;
   }
   if (!setupRotLayer())
   {
      QgsMessageLog::logMessage("failed to setup layer", name(), Qgis::MessageLevel::Info);
      return false;
   }
   if (!setupYhsLayer())
   {
      QgsMessageLog::logMessage("failed to setup layer", name(), Qgis::MessageLevel::Info);
      return false;
   }
   m_layers_setup = true;
   QgsMessageLog::logMessage("Layers set up ", name(), Qgis::MessageLevel::Info);

   return true;
}

bool pnwRotationPlugin::setupRotLayer() // Rot layer must be loaded in qgis first (so not at qgis launch)
{
   QgsMessageLog::logMessage(QString("Setup dest layer "), name(), Qgis::MessageLevel::Info);

   QString providerName = "memory";     // Or "ogr" for file-based data
   QString uri = "Point?crs=epsg:4326"; // Example URI for memory provider
   if (!m_rotDestLayer)
      m_rotDestLayer = new QgsVectorLayer(uri, s_rotDestLayerName, providerName);
   if (!m_rotDestLayer->isValid())
   {
      qDebug() << "Could not instantiate plugin target layer";
      delete m_rotDestLayer; // Clean up if creation failed
      return false;
   }

   // copy over symbology from rot layer
   QgsFeatureRenderer *sourceRenderer = m_rotSrcLayer->renderer();
   if (sourceRenderer)
   {
      QgsFeatureRenderer *clonedRenderer = sourceRenderer->clone();
      if (clonedRenderer)
      {
         m_rotDestLayer->setRenderer(clonedRenderer);
      }
   }

   // copy rot data fields over
   if (m_rotDestLayer->dataProvider()->addAttributes(m_fieldList))
      m_rotDestLayer->updateFields();

   return true;
}

void pnwRotationPlugin::displayRotData(QgsFeatureList &featureList)
{
   // copy feature data over
   m_rotDestLayer->dataProvider()->addFeatures(featureList);

   QgsProject::instance()->addMapLayer(m_rotDestLayer);

   m_rotDestLayer->triggerRepaint();
}

bool pnwRotationPlugin::setupYhsLayer() // Rot layer must be loaded in qgis first (so not at qgis launch)
{
   QgsMessageLog::logMessage(QString("Setup YHS layer "), name(), Qgis::MessageLevel::Info);

   QString providerName = "memory";     // Or "ogr" for file-based data
   QString uri = "LineString?crs=epsg:4326"; // Example URI for memory provider
   if (!m_yhsDestLayer)
      m_yhsDestLayer = new QgsVectorLayer("LineString?crs=EPSG:4326", s_yhsDestLayerName, "memory");
   if (!m_yhsDestLayer->isValid())
   {
      qDebug() << "Could not instantiate plugin target layer";
      delete m_yhsDestLayer; // Clean up if creation failed
      return false;
   }

   QgsSingleSymbolRenderer *renderer = dynamic_cast<QgsSingleSymbolRenderer *>(m_yhsDestLayer->renderer());
   QgsLineSymbol *lineSymbol = dynamic_cast<QgsLineSymbol *>(renderer->symbol());
   lineSymbol->setWidth(0.5);
   lineSymbol->setColor(QColor (255, 0, 0)); // Set color to blue

   m_line << QgsPointXY(YHS_lon, YHS_lat); // Set initial point to current YHS

   return true;
}

void pnwRotationPlugin::yhs_menu_button_action()
{
   if (!setupLayers())
      return;

   pState p_last = m_pYhsState.at(m_line.length() - 1);
   while (p_last.lon > longitudeLimit)
   {
      // Get closest rotation vector velocity from rotation field
      QgsFeature rotFeature = getClosestRotEntry(p_last.lon, p_last.lat);
      double rotVe = getFeatureAttrubute(rotFeature, 2);
      double rotVn = getFeatureAttrubute(rotFeature, 3);

      // Get state change
       double deltaLon = longitudeFromDistance(p_last.lat, p_last.ve * detlaT);
      double deltaLat = latitudeFromDisatnce(p_last.vn * detlaT);
      double deltaVe = rotVe; // mm/yr
      double deltaVn = rotVn; // mm/yr

      pState p_next = {
          p_last.lon + deltaLon,
          p_last.lat + deltaLat,
         p_last.ve + deltaVe,
          p_last.vn + deltaVn};

      m_pYhsState.push_back(p_next);

      ////////////////// Update polyline layer line
      QgsPointXY nextPoint(p_next.lon, p_next.lat);
      m_line << nextPoint;

      if (m_verbose)
      {
         if (m_line.length() == 2)
         {
            QString entryDataString("NA Move \t");
            entryDataString += "E: " + QString::number(YHS_lon) + " deg,\t";
            entryDataString += "N: " + QString::number(YHS_lat) + " deg\t";
            entryDataString += "(E: " + QString::number(m_NA_Vel_E * detlaT / 1E6) + " km\t";
            entryDataString += "N: " + QString::number(m_NA_Vel_N * detlaT / 1E6) + " km)";
            QgsMessageLog::logMessage(entryDataString, name(), Qgis::MessageLevel::Info);
         }
         QString entryDataString("Rot Move \t");
         entryDataString += "E: " + QString::number(deltaLon) + " deg,\t";
         entryDataString += "N: " + QString::number(deltaLat) + " deg\t";
         entryDataString += "(E: " + QString::number(deltaVe * detlaT / 1E6) + " km\t";
         entryDataString += "N: " + QString::number(deltaVn * detlaT / 1E6) + " km)";
         QgsMessageLog::logMessage(entryDataString, name(), Qgis::MessageLevel::Info);
      }

      displayYhsData(m_line);

      ///////////////////// Update Rot layer vector
      if (m_verbose)
         printFeature(rotFeature, QString("YHS Rot: "));

      m_rotFeatureList2.push_back(rotFeature);
      displayRotData(m_rotFeatureList2);

      p_last = p_next;
      m_passes++;
   }
   QgsMessageLog::logMessage(QString("Completed. Passes = ") + QString::number(m_passes), name(), Qgis::MessageLevel::Info);
}

void pnwRotationPlugin::displayYhsData(QgsPolylineXY &line)
{
   QgsFeature feature(m_yhsDestLayer->fields());
   feature.setGeometry(QgsGeometry::fromPolylineXY(line));

   // Add feature to the layer
   m_yhsDestLayer->startEditing();
   m_yhsDestLayer->addFeature(feature);
   m_yhsDestLayer->commitChanges();
   QgsProject::instance()->addMapLayer(m_yhsDestLayer);

   m_yhsDestLayer->triggerRepaint();
}

void pnwRotationPlugin::rot_menu_button_action()
{
   if (!setupLayers())
      return;

   clear_display_data();

   displayRotData(m_rotFeatureList);
}

void pnwRotationPlugin::clear_menu_button_action()
{
   QgsMessageLog::logMessage(QString("Clearing yhs data"), name(), Qgis::MessageLevel::Info);
   clear_display_data();
}

void pnwRotationPlugin::clear_display_data()
{
   // Clear rotDestLayer
   if (!m_rotDestLayer)
   {
      QgsMessageLog::logMessage("Error: Provided m_rotDestLayer is null.", name(), Qgis::MessageLevel::Info);
      return;
   }

   m_rotDestLayer->startEditing();
   QgsVectorDataProvider *dataProvider = m_rotDestLayer->dataProvider();

   if (!dataProvider)
   {
      QgsMessageLog::logMessage("Error: Data provider not found for the m_rotDestLayer.", name(), Qgis::MessageLevel::Info);
      return;
   }

   if (!dataProvider->truncate())
   {
      QgsMessageLog::logMessage(QString("Failed to truncate m_rotDestLayer '%1'.").arg(m_rotDestLayer->name()), name(), Qgis::MessageLevel::Info);
      return;
   }

   m_rotDestLayer->commitChanges();
   m_rotDestLayer->triggerRepaint();

   // Clear yhsDestLayer
   if (!m_yhsDestLayer)
   {
      QgsMessageLog::logMessage("Error: Provided m_yhsDestLayer is null.", name(), Qgis::MessageLevel::Info);
      return;
   }

   m_line.clear();
   m_line << QgsPointXY(YHS_lon, YHS_lat); // Current YHS hotspot

   m_yhsDestLayer->startEditing();
   QgsVectorDataProvider *yhsDataProvider = m_yhsDestLayer->dataProvider();

   if (!yhsDataProvider)
   {
      QgsMessageLog::logMessage("Error: Data provider not found for the m_yhsDestLayer.", name(), Qgis::MessageLevel::Info);
      return;
   }

   if (!yhsDataProvider->truncate())
   {
      QgsMessageLog::logMessage(QString("Failed to truncate m_yhsDestLayer '%1'.").arg(m_yhsDestLayer->name()), name(), Qgis::MessageLevel::Info);
      return;
   }

   m_yhsDestLayer->commitChanges();
   m_yhsDestLayer->triggerRepaint();
}

bool pnwRotationPlugin::loadRotData()
{
   if (m_rotDataLoaded)
      return true;

   QgsMessageLog::logMessage(QString("loadRotData"), name(), Qgis::MessageLevel::Info);

   // get rotation data layer
   QgsMapLayer *mapLayer = QgsProject::instance()->mapLayersByName(s_rotDatadestLayerName).value(0);
   m_rotSrcLayer = dynamic_cast<QgsVectorLayer *>(mapLayer);
   if (!m_rotSrcLayer)
   {
      QString loadErrorMsg = QString("Could not load source layer ") + s_rotDatadestLayerName;
      QgsMessageLog::logMessage(loadErrorMsg, name(), Qgis::MessageLevel::Info);
      return false;
   }

   // get rot data fields
   QgsFields fields = m_rotSrcLayer->fields();
   for (const QgsField &field : fields)
      m_fieldList.append(field);

   // gety rot features
   QgsFeatureIterator featureIt = m_rotSrcLayer->getFeatures();
   QgsFeature feature;
   while (featureIt.nextFeature(feature))
   {
      m_rotFeatureList << feature;

      if (m_verbose)
         printFeature(feature, QString("Loaded "), fields.size());
   }

   m_rotDataLoaded = true;
   return true;
}

bool pnwRotationPlugin::setFeatureAttribute(QgsFeature &feature, int index, double value)
{
   QVariant variant;

   if (index < 4)
      variant = QVariant::fromValue(value);
   else if (index < 7)
      variant = QVariant::fromValue(0.0);
   else
      variant = QVariant::fromValue(QString("SJW"));
   if (!feature.setAttribute(index, variant))
   {
      qDebug() << "Could not set feature attribute";
      return false;
   }
   return true;
}

double pnwRotationPlugin::getFeatureAttrubute(QgsFeature &feature, int index)
{
   QVariant attributeValueByIndex = feature.attribute(index);
   return attributeValueByIndex.toDouble();
}

double pnwRotationPlugin::latitudeFromDisatnce(double distanceN)
{
   double latitude = atan(distanceN / (EARTH_RADIUS * 1000)) * 180.0 / M_PI;
   return latitude;
}

// Function to calculate new longitude after moving eastward
// distance is in mm
double pnwRotationPlugin::longitudeFromDistance(double latitude, double distance)
{

   double latitudeRadians = latitude * M_PI / 180.0;

   // Calculate the radius of the parallel of latitude at the starting latitude
   double radiusOfParallel = EARTH_RADIUS * std::cos(latitudeRadians);

   // Calculate the change in longitude in radians
   double deltaLongitudeRadians = distance / (radiusOfParallel * 1000);

   return deltaLongitudeRadians * 180.0 / M_PI;
}

QgsFeature pnwRotationPlugin::getClosestRotEntry(double lon, double lat)
{
   // get rot features
   QgsFeatureIterator featureIt = m_rotSrcLayer->getFeatures();
   QgsFeature feature;
   QgsFeature closestFeature;

   double Ve = 0, Vn = 0;
   double p_lon, p_lat;
   double minDist = 1e10;
   while (featureIt.nextFeature(feature))
   {
      double f_lon = getFeatureAttrubute(feature, 0);
      double f_lat = getFeatureAttrubute(feature, 1);
      double dist = (sqr(f_lon - lon) + sqr(f_lat - lat));
      if (dist < minDist)
      {
         minDist = dist;
         closestFeature = feature;
      }
   }
   return closestFeature;
}

void pnwRotationPlugin::printFeature(QgsFeature feature, QString label, int fields)
{
   QString entryDataString = label;
   for (int i = 0; i < fields; ++i)
   {
      QVariant attributeValueByIndex = feature.attribute(i);
      entryDataString += attributeValueByIndex.toString() + "\t";
   }
   QgsMessageLog::logMessage(entryDataString, name(), Qgis::MessageLevel::Info);
}
