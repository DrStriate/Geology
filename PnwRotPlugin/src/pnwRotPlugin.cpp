#include "pnwRotPlugin.h"
#include <cmath>

namespace
{
   const QString s_name = QStringLiteral("Pnw Rotation Analysis");
   const QString s_description = QStringLiteral("Pnw rotation analysis plugin");
   const QString s_category = QStringLiteral("Plugins");
   const QString s_version = QStringLiteral("Version 1.2.3");
   const QString s_icon = QStringLiteral(":/plugin.svg");
   const QString s_rotDatadestLayerName = QStringLiteral("nshm2023_GPS_velocity");
   const QgisPlugin::PluginType s_type = QgisPlugin::UI;
   const QString s_destLayerName = "PnwPluginData";
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
}

void pnwRotationPlugin::unload()
{
   // TODO - need to remove the actions from the menu again.
   // Get the QgsProject instance
   QgsProject *project = QgsProject::instance();

   // Find the layer by name
   // mapLayersByName() returns a list, even if there is only one match
   QList<QgsMapLayer *> foundLayers = project->mapLayersByName(s_destLayerName);

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
      qDebug() << "Layer not found:" << s_destLayerName;
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
   if (!setupDestLayer())
   {
      QgsMessageLog::logMessage("failed to setup layer", name(), Qgis::MessageLevel::Info);      
      return false;
   }
   m_layers_setup = true;
   QgsMessageLog::logMessage("Layers set up ", name(), Qgis::MessageLevel::Info);

   return true;
}

void pnwRotationPlugin::yhs_menu_button_action()
{
   if (!setupLayers())
      return;

   QgsFeature feature(m_fieldList);

   if (m_yhsFeatureList.size() == 0) // load first yhs loc/motion
   {
      // Create a point geometry
      QgsPoint pointGeometry(YHS_lon, YHS_lat);
      setFeatureAttribute(feature, 0, YHS_lon);
      setFeatureAttribute(feature, 1, YHS_lat);
      QgsGeometry geometry = QgsGeometry::fromPointXY(pointGeometry);
      feature.setGeometry(geometry);

      double deltaLat = latitudeFromDisatnce(NA_Vel_N * detlaT);
      double deltaLon = longitudeFromDistance(YHS_lat, NA_Vel_E * detlaT);
      setFeatureAttribute(feature, 2, deltaLon);
      setFeatureAttribute(feature, 3, deltaLat);
   }
   else // add incremental move
   {

   }
   if (m_verbose)
   {
      QString entryDataString;
      for (int i = 0; i < m_fieldList.length(); ++i)
      {
         QVariant attributeValueByIndex = feature.attribute(i);
         entryDataString += attributeValueByIndex.toString() + "\t";
      }
      QgsMessageLog::logMessage(entryDataString, name(), Qgis::MessageLevel::Info);
   }
   
   m_yhsFeatureList.push_back(feature);
   displayData(m_yhsFeatureList);
}

void pnwRotationPlugin::rot_menu_button_action()
{
   if (!setupLayers())
      return;

   clear_display_data();

   displayData(m_rotFeatureList);
}

void pnwRotationPlugin::clear_menu_button_action()
{
   QgsMessageLog::logMessage(QString("Clearing yhs data"), name(), Qgis::MessageLevel::Info);
   m_yhsFeatureList.clear();
   clear_display_data();
}

void pnwRotationPlugin::clear_display_data()
{
   if (!m_DestLayer)
   {
      QgsMessageLog::logMessage("Error: Provided m_DestLayer is null.", name(), Qgis::MessageLevel::Info);
      return;
   }

   m_DestLayer->startEditing();
   QgsVectorDataProvider *dataProvider = m_DestLayer->dataProvider();

   if (!dataProvider)
   {
      QgsMessageLog::logMessage("Error: Data provider not found for the m_DestLayer.", name(), Qgis::MessageLevel::Info);
      return;
   }

   if (!dataProvider->truncate())
   {
      QgsMessageLog::logMessage(QString("Failed to truncate m_DestLayer '%1'.").arg(m_DestLayer->name()), name(), Qgis::MessageLevel::Info);
      return;
   }

   m_DestLayer->commitChanges();
   m_DestLayer->triggerRepaint();
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
      {
         QString entryDataString;
         for (int i = 0; i < fields.count(); ++i)
         {
            QVariant attributeValueByIndex = feature.attribute(i);
            entryDataString += attributeValueByIndex.toString() + "\t";
         }
         QgsMessageLog::logMessage(entryDataString, name(), Qgis::MessageLevel::Info);
      }
   }

   m_rotDataLoaded = true;
   return true;
}

bool pnwRotationPlugin::setupDestLayer() // Rot layer must be loaded in qgis first (so not at qgis launch)
{
   QgsMessageLog::logMessage(QString("Setup dest layer "), name(), Qgis::MessageLevel::Info);

   QString providerName = "memory";     // Or "ogr" for file-based data
   QString uri = "Point?crs=epsg:4326"; // Example URI for memory provider
   if (!m_DestLayer)
      m_DestLayer = new QgsVectorLayer(uri, s_destLayerName, providerName);
   if (!m_DestLayer->isValid())
   {
      qDebug() << "Could not instantiate plugin target layer";
      delete m_DestLayer; // Clean up if creation failed
      return false;
   }

   // copy over symbology from rot layer
   QgsFeatureRenderer *sourceRenderer = m_rotSrcLayer->renderer();
   if (sourceRenderer)
   {
      QgsFeatureRenderer *clonedRenderer = sourceRenderer->clone();
      if (clonedRenderer)
      {
         m_DestLayer->setRenderer(clonedRenderer);
      }
   }

   // copy rot data fields over
   if (m_DestLayer->dataProvider()->addAttributes(m_fieldList))
      m_DestLayer->updateFields();
   
   return true;
}


void pnwRotationPlugin::displayData(QgsFeatureList& featureList)
{
   // copy feature data over
   m_DestLayer->dataProvider()->addFeatures(featureList);

   QgsProject::instance()->addMapLayer(m_DestLayer);

   m_DestLayer->triggerRepaint();
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

double pnwRotationPlugin::latitudeFromDisatnce (double distanceN)
{
   double latitude = atan(distanceN / EARTH_RADIUS) * 180.0 / M_PI;
   return latitude;
}

// Function to calculate new longitude after moving eastward
double pnwRotationPlugin::longitudeFromDistance(double latitude, double distance) {

    double latitudeRadians = latitude * M_PI / 180.0;

    // Calculate the radius of the parallel of latitude at the starting latitude
    double radiusOfParallel = EARTH_RADIUS * std::cos(latitudeRadians);

    // Calculate the change in longitude in radians
    double deltaLongitudeRadians = distance / radiusOfParallel;

    return deltaLongitudeRadians * 180.0 / M_PI;
}
