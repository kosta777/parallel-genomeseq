include(ExternalProject)

add_library(Eigen INTERFACE)

ExternalProject_Add(
        EigenProject
        URL ${PROJECT_SOURCE_DIR}/cmake/eigen-3.3.7.zip
        URL_MD5 0d9c8496922d5c07609b9f3585f00e49
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

add_dependencies(Eigen EigenProject)
target_include_directories(Eigen
        SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/Eigen
        )