add_library(yocto_texture yocto_texture.h yocto_texture.cpp)

set_target_properties(yocto_texture PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yocto_texture PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yocto_texture yocto)
