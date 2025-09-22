#include "pnwRotPlugin.h"

namespace
{
   const QString s_name = QStringLiteral("Pnw Rotation Analysis");
   const QString s_description = QStringLiteral("Pnw rotation analysis plugin");
   const QString s_category = QStringLiteral("Plugins");
   const QString s_version = QStringLiteral("Version 1.2.3");
   const QString s_icon = QStringLiteral(":/plugin.svg");
   const QString s_rotDatadestLayerName = QStringLiteral("nshm2023_GPS_velocity");
   const QgisPlugin::PluginType s_type = QgisPlugin::UI;
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
   m_rotDataLoaded = false;
}

void pnwRotationPlugin::unload()
{
   // TODO - need to remove the actions from the menu again.
}

void pnwRotationPlugin::initGui()
{
   QgsMessageLog::logMessage(QString("initGui"), name(), Qgis::MessageLevel::Info);

   if (m_qgis_if == nullptr)
   {
      QgsMessageLog::logMessage("failed to get the handle to QGIS", name(), Qgis::MessageLevel::Info);
      return;
   }

   // add an action to the menu
   m_display_rot_menu_action = new QAction(QIcon(""), QString("Display rot data"), this);
   connect(m_display_rot_menu_action, SIGNAL(triggered()), this, SLOT(display_menu_button_action()));
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_display_rot_menu_action);

   // add an action to the menu
   m_clear_menu_action = new QAction(QIcon(""), QString("Clear data"), this);
   connect(m_clear_menu_action, SIGNAL(triggered()), this, SLOT(clear_menu_button_action()));
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_clear_menu_action);
}

void pnwRotationPlugin::display_menu_button_action()
{
   if (!loadRotData(s_rotDatadestLayerName, m_verbose))
   {
      QgsMessageLog::logMessage("failed to load rot data from layer", s_rotDatadestLayerName, Qgis::MessageLevel::Info);
      return;
   }
   if (!displayData())
   {
      QgsMessageLog::logMessage("failed to write rot data to layer", name(), Qgis::MessageLevel::Info);
      return;
   }
   return;
}

void pnwRotationPlugin::clear_menu_button_action()
{
}

bool pnwRotationPlugin::loadRotData(QString rotDataLayer, bool verbose)
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
      m_featureList << feature;

   m_rotDataLoaded = true;
   return true;
}

bool pnwRotationPlugin::displayData()
{
   QgsMessageLog::logMessage(QString("displayData"), name(), Qgis::MessageLevel::Info);

   QString destLayerName = "PnwPluginData";
   QString providerName = "memory";     // Or "ogr" for file-based data
   QString uri = "Point?crs=epsg:4326"; // Example URI for memory provider
   QgsVectorLayer *m_DestLayer = new QgsVectorLayer(uri, destLayerName, providerName);
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

   // copy feature data over
   m_DestLayer->dataProvider()->addFeatures(m_featureList);

   QgsProject::instance()->addMapLayer(m_DestLayer);
   return true;
}
