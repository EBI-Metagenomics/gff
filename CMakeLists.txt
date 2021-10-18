cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(
  gff
  VERSION 0.0.3
  LANGUAGES C)

# Generate compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(external/cmake-boilerplate/cbp.cmake)
cbp_set_rpath()

add_library(gff src/error.c src/gff.c src/fsm.c src/tok.c src/strlcpy.c)
add_library(GFF::gff ALIAS gff)

cbp_hide_functions(gff)
cbp_set_warning_flags()
cbp_ci_build_option()
cbp_generate_export_header(gff EXPORT_FILE)

target_include_directories(
  gff
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(gff PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(gff PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
target_compile_features(gff PRIVATE c_std_11)
target_compile_definitions(
  gff PRIVATE $<$<BOOL:${WIN32}>:_CRT_NONSTDC_NO_DEPRECATE>)
target_compile_definitions(gff PRIVATE $<$<BOOL:${WIN32}>:__extension__>)
target_compile_options(gff PRIVATE ${WARNING_FLAGS})
target_compile_options(gff PRIVATE "-fverbose-asm")
target_link_libraries(gff INTERFACE $<$<BOOL:${UNIX}>:m>)

cbp_install(gff MERGED)
install(FILES ${EXPORT_FILE} DESTINATION include/gff)

enable_testing()
add_subdirectory(test)
