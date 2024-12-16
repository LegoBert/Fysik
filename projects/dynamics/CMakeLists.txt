#--------------------------------------------------------------------------
# spacegame project
#--------------------------------------------------------------------------

PROJECT(dynamics)
FILE(GLOB project_headers code/*.h)
FILE(GLOB project_sources code/*.cc)

SET(files_project ${project_headers} ${project_sources})
SOURCE_GROUP("dynamics" FILES ${files_project})

ADD_EXECUTABLE(dynamics ${files_project})
TARGET_LINK_LIBRARIES(dynamics core render)
ADD_DEPENDENCIES(dynamics core render)

IF(MSVC)
    set_property(TARGET dynamics PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
ENDIF()