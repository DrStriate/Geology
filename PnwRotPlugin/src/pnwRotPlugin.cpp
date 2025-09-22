#include "pnwRotPlugin.h"

namespace
{
   const QString s_name = QStringLiteral("Pnw Rotation Analysis");
   const QString s_description = QStringLiteral("Pnw rotation analysis plugin");
   const QString s_category = QStringLiteral("Plugins");
   const QString s_version = QStringLiteral("Version 1.2.3");
   const QString s_icon = QStringLiteral(":/plugin.svg");
   const QString s_rotDataLayerName = QStringLiteral("nshm2023_GPS_velocity");
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
   m_load_menu_action = new QAction(QIcon(""), QString("Load rot data"), this);
   connect(m_load_menu_action, SIGNAL(triggered()), this, SLOT(load_menu_button_action()));
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_load_menu_action);

   // add an action to the menu
   m_clear_menu_action = new QAction(QIcon(""), QString("Clear data"), this);
   connect(m_clear_menu_action, SIGNAL(triggered()), this, SLOT(clear_menu_button_action()));
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_clear_menu_action);
}

void pnwRotationPlugin::load_menu_button_action()
{
   if (!loadRotData(s_rotDataLayerName, m_verbose))
   {
      QgsMessageLog::logMessage("failed to load rot data from layer", name(), Qgis::MessageLevel::Info);
   }
}

void pnwRotationPlugin::clear_menu_button_action()
{
}

bool pnwRotationPlugin::loadRotData(QString rotDataLayer, bool verbose)
{
   // going to suck features from rot layer and write to new layer
   QgsMapLayer *mapLayer = QgsProject::instance()->mapLayersByName(s_rotDataLayerName).value(0);
   QgsVectorLayer *vlayer = dynamic_cast<QgsVectorLayer *>(mapLayer);
   if (!vlayer)
   {
      QString loadErrorMsg = QString("Could not load source layer ") + s_rotDataLayerName;
      QgsMessageLog::logMessage(loadErrorMsg, name(), Qgis::MessageLevel::Info);
      return false;
   }

   QString layerName = "PnwPluginData";
   QString providerName = "memory"; // Or "ogr" for file-based data
   QString uri = "Point?crs=epsg:4326"; // Example URI for memory provider
   QgsVectorLayer *mypLayer = new QgsVectorLayer(uri, layerName, providerName);
   if (!mypLayer->isValid())
   {
      qDebug() << "Could not instantiate plugin target layer";
      delete mypLayer; // Clean up if creation failed
      return false;
   }

   // copy fields over
   m_fields = vlayer->fields();
   QList<QgsField> fieldList;
   for (const QgsField &field : m_fields)
      fieldList.append(field);
   if (mypLayer->dataProvider()->addAttributes(fieldList))
      mypLayer->updateFields();

   // copy features over
   QgsFeatureIterator featureIt = vlayer->getFeatures();
   QgsFeature feature;
   QgsFeatureList featureList;
   while (featureIt.nextFeature(feature))
      featureList << feature;
   mypLayer->dataProvider()->addFeatures(featureList);

   if (mypLayer->isValid())
   {
      QgsProject::instance()->addMapLayer(mypLayer);
      return true;
   }
   else
   {
      QgsMessageLog::logMessage("failed to create ", name(), Qgis::MessageLevel::Info);
      return false;
   }
}
