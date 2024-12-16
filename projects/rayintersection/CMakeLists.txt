#--------------------------------------------------------------------------
# spacegame project
#--------------------------------------------------------------------------

PROJECT(rayintersection)
FILE(GLOB project_headers code/*.h)
FILE(GLOB project_sources code/*.cc)

SET(files_project ${project_headers} ${project_sources})
SOURCE_GROUP("rayintersection" FILES ${files_project})

ADD_EXECUTABLE(rayintersection ${files_project})
TARGET_LINK_LIBRARIES(rayintersection core render)
ADD_DEPENDENCIES(rayintersection core render)

IF(MSVC)
    set_property(TARGET rayintersection PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
ENDIF()