cmake_minimum_required(VERSION 3.20)
project(CheMPAS_TMT VERSION 1.0.0 LANGUAGES C Fortran)
include(GNUInstallDirs)
include(CTest)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set(TMT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../TMT" CACHE PATH "TMT source checkout")
if(NOT EXISTS "${TMT_SOURCE_DIR}/CMakeLists.txt")
  message(FATAL_ERROR "Set TMT_SOURCE_DIR to the TMT checkout")
endif()
if(DEFINED TMT_BUILD_FORTRAN AND NOT TMT_BUILD_FORTRAN)
  message(FATAL_ERROR "CheMPAS-A requires TMT_BUILD_FORTRAN=ON")
endif()

# Keep the ordinary TMT installation intact. Build an additional, isolated
# Fortran interface for applications that also link MUSICA's MIEM component.
set(chempas_testing ${BUILD_TESTING})
set(BUILD_TESTING OFF)
add_subdirectory("${TMT_SOURCE_DIR}" tmt)
set(BUILD_TESTING ${chempas_testing})

file(GLOB_RECURSE namespace_inputs CONFIGURE_DEPENDS
  "${TMT_SOURCE_DIR}/fortran/*.F90" "${TMT_SOURCE_DIR}/include/*.h")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  ${namespace_inputs} "${CMAKE_CURRENT_SOURCE_DIR}/namespace.py")
find_package(Python3 REQUIRED COMPONENTS Interpreter)
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/namespace.py"
  "${TMT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/namespaced"
  RESULT_VARIABLE namespace_result)
if(NOT namespace_result EQUAL 0)
  message(FATAL_ERROR "Could not generate the isolated TMT interface")
endif()

get_target_property(tmt_sources tmt SOURCES)
list(TRANSFORM tmt_sources PREPEND "${TMT_SOURCE_DIR}/")
add_library(chempas_tmt SHARED ${tmt_sources})
set_target_properties(chempas_tmt PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED YES C_EXTENSIONS NO)
target_include_directories(chempas_tmt PRIVATE "$<TARGET_PROPERTY:tmt,INCLUDE_DIRECTORIES>")
target_compile_definitions(chempas_tmt PRIVATE "$<TARGET_PROPERTY:tmt,COMPILE_DEFINITIONS>")
target_compile_options(chempas_tmt PRIVATE "$<TARGET_PROPERTY:tmt,COMPILE_OPTIONS>"
  -include "${CMAKE_CURRENT_BINARY_DIR}/namespaced/namespace.h")
target_link_options(chempas_tmt PRIVATE "$<TARGET_PROPERTY:tmt,LINK_OPTIONS>")
target_link_libraries(chempas_tmt PRIVATE "$<TARGET_PROPERTY:tmt,LINK_LIBRARIES>")

get_target_property(fortran_sources tmt_fortran SOURCES)
set(namespaced_fortran)
foreach(source IN LISTS fortran_sources)
  if(IS_ABSOLUTE "${source}")
    file(RELATIVE_PATH source "${TMT_SOURCE_DIR}" "${source}")
  endif()
  list(APPEND namespaced_fortran "${CMAKE_CURRENT_BINARY_DIR}/namespaced/${source}")
endforeach()
add_library(chempas_tmt_fortran STATIC ${namespaced_fortran})
set_target_properties(chempas_tmt_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/chempas_modules"
  POSITION_INDEPENDENT_CODE ON)
target_include_directories(chempas_tmt_fortran PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/chempas_modules>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/chempas-tmt>")
target_include_directories(chempas_tmt_fortran PRIVATE "${TMT_SOURCE_DIR}/include")
target_compile_options(chempas_tmt_fortran PRIVATE "$<TARGET_PROPERTY:tmt_fortran,COMPILE_OPTIONS>")
target_link_libraries(chempas_tmt_fortran PUBLIC chempas_tmt)

if(APPLE)
  set_target_properties(chempas_tmt PROPERTIES INSTALL_RPATH "@loader_path")
else()
  set_target_properties(chempas_tmt PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()
install(TARGETS chempas_tmt chempas_tmt_fortran
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/chempas_modules/"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/chempas-tmt)
configure_file(chempas-tmt-fortran.pc.in chempas-tmt-fortran.pc @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/chempas-tmt-fortran.pc"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if(BUILD_TESTING)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(MUSICA REQUIRED IMPORTED_TARGET musica-fortran=0.16.5)
  add_executable(test_coexistence test_coexistence.F90)
  target_link_libraries(test_coexistence PRIVATE chempas_tmt_fortran PkgConfig::MUSICA)
  add_test(NAME tmt_miem_coexistence COMMAND test_coexistence
    "${TMT_SOURCE_DIR}/tests/configs/decay.json"
    "${TMT_DEPENDENCY_DIR}/turbotuv/examples/tuv_5_4.json")
  set_tests_properties(tmt_miem_coexistence PROPERTIES
    WORKING_DIRECTORY "${TMT_DEPENDENCY_DIR}/turbotuv"
    ENVIRONMENT "TMT_CACHE=${CMAKE_CURRENT_BINARY_DIR}/test-cache")
endif()
