cmake_minimum_required(VERSION 3.16)

SET(MAIN_PROJ "PNWRotation")
project(${MAIN_PROJ} CXX)

SET(ROOT_DIR ${CMAKE_SOURCE_DIR})
SET(BIN_DIR ${ROOT_DIR}/bin)
SET(LIB_DIR ${ROOT_DIR}/lib)
set(CMAKE_CXX_STANDARD 20)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# Debug by default
# Turn on Release with "-DCMAKE_BUILD_TYPE=Release"
IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Debug )
ENDIF()
message(">>>>>>>>>>>>>>>>> Building ${CMAKE_BUILD_TYPE}")

if (MSVC)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
endif()

############### Main  #################
add_executable(${MAIN_PROJ})

target_link_libraries(${MAIN_PROJ} PRIVATE  
    ${OS_LIBRARIES}
#    ${CUDA_LIBRARIES}
)

if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux" ) 
target_link_libraries(${MAIN_PROJ} PRIVATE 
#    ${LIB_DIR}/lib${CORE_LIB}.so 
    $ENV{OpenCV_LIB}/libopencv_core.so
    $ENV{OpenCV_LIB}/libopencv_highgui.so
    $ENV{OpenCV_LIB}/libopencv_imgproc.so
    $ENV{OpenCV_LIB}/libopencv_imgcodecs.so)
else ()
target_link_libraries(${MAIN_PROJ} PRIVATE 
#   ${LIB_DIR}/${CORE_LIB}.lib 
    $ENV{OpenCV_LIB}/opencv_world420d.lib)
endif ()

target_sources(${MAIN_PROJ} PRIVATE
    main.cpp
)

target_include_directories(${MAIN_PROJ} PUBLIC
 #   ${ROOT_DIR}/Common 
 #   ${ROOT_DIR}/Core
    ${ROOT_DIR}/../../eigen-3.4.0
    $ENV{OpenCV_INC}
)

#SET_property(TARGET ${MAIN_PROJ} PROPERTY CUDA_ARCHITECTURES 61 75 87) 


