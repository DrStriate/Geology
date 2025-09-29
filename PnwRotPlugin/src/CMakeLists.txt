#C:\Users\stria\AppData\Local\Programs\OSGeo4W\apps\Qt5\lib\cmake\Qt5
# set(Qt5_DIR ${OSGEO4W_ROOT}/apps/Qt5/lib/cmake/Qt5)

add_library(pnwRotationPlugin MODULE
  pnwRotPlugin.cpp
)

target_link_libraries(pnwRotationPlugin
  ${LIBS}
)

target_include_directories(pnwRotationPlugin PUBLIC
  ${INCLUDES}
)
