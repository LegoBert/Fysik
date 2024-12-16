#--------------------------------------------------------------------------
# spacegame project
#--------------------------------------------------------------------------

PROJECT(collision)
FILE(GLOB project_headers code/*.h)
FILE(GLOB project_sources code/*.cc)

SET(files_project ${project_headers} ${project_sources})
SOURCE_GROUP("collision" FILES ${files_project})

ADD_EXECUTABLE(collision ${files_project})
TARGET_LINK_LIBRARIES(collision core render)
ADD_DEPENDENCIES(collision core render)

IF(MSVC)
    set_property(TARGET collision PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
ENDIF()