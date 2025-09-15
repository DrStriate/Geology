#include "pnwRotPlugin.h"

namespace {
   const QString s_name = QStringLiteral("Pnw Rotation Analysis");
   const QString s_description = QStringLiteral("Pnw rotation analysis plugin");
   const QString s_category = QStringLiteral("Plugins");
   const QString s_version = QStringLiteral("Version 1.2.3");
   const QString s_icon =  QStringLiteral( ":/plugin.svg" );
   const QgisPlugin::PluginType s_type = QgisPlugin::UI;
}

QGISEXTERN QgisPlugin* classFactory(QgisInterface* qgis_if)
{
   std::cout << "::classFactory" << std::endl;
   return new pnwRotationPlugin(qgis_if);
}

QGISEXTERN const QString* name() {
   return &s_name;
}

QGISEXTERN const QString* description() {
   return &s_description;
}

QGISEXTERN const QString* category() {
   return &s_category;
}

QGISEXTERN const QString* version() {
   return &s_version;
}

QGISEXTERN const QString* icon() {
   return &s_icon;
}

QGISEXTERN int type() {
   return s_type;
}

QGISEXTERN void unload(QgisPlugin* plugin) {
   std::cout << "::unload" << std::endl;
   delete plugin;
}

pnwRotationPlugin::pnwRotationPlugin(QgisInterface* iface) : QgisPlugin(s_name, s_description, s_category, s_version, s_type), m_qgis_if(iface) {
}

void pnwRotationPlugin::unload() {
   // TODO - need to remove the actions from the menu again.
}

void pnwRotationPlugin::initGui() {
   std::cout << "pnwRotationPlugin::initGui" << std::endl;

   // add an action to the menu
   m_menu_action = new QAction(QIcon(""), QString("Pnw rotation analysis"), this);
   connect(m_menu_action, SIGNAL(triggered()), this, SLOT(menu_button_action()));
   if (m_qgis_if == nullptr) {
      std::cout << "failed to get the handle to QGIS" << std::endl;
      return;
   }
   m_qgis_if->addPluginToMenu(QString("&PnwRotationPlugin"), m_menu_action);
}

void pnwRotationPlugin::menu_button_action() {
   QgsMessageLog::logMessage(QString("Menu clicked!"), QString("Pnw Rotation Plugin"), Qgis::MessageLevel::Info);
}

