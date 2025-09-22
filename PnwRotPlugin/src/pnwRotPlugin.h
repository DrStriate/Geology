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
   void display_menu_button_action();

private:
   QgisInterface* m_qgis_if;

   /// The actions in the QGIS menu bar.
   QAction *m_clear_menu_action;
   QAction *m_display_rot_menu_action;

   bool m_rotDataLoaded;
   bool loadRotData(QString rotDataLayer, bool verbose);
   bool displayData();

   QgsVectorLayer *m_rotSrcLayer;
   QgsVectorLayer *m_DestLayer;
   // Rot layer data
   QList<QgsField> m_fieldList;   
   QgsFeatureList m_featureList;

   std::vector<QString> m_fieldNames;
   std::vector<std::vector<double>> m_rot_data;
   
   const bool m_verbose = true;
};

#endif
