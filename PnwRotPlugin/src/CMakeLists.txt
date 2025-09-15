add_library(pnwRotationPlugin MODULE
  pnwRotPlugin.cpp
)

target_link_libraries(pnwRotationPlugin
  ${LIBS}
)

target_include_directories(pnwRotationPlugin PUBLIC
  ${INCLUDES}
)
