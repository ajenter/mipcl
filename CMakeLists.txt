cmake_minimum_required(VERSION 3.12)

get_filename_component(BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}" ABSOLUTE)

file(COPY "lib" DESTINATION "${BUILD_DIR}")
file(COPY "include" DESTINATION "${BUILD_DIR}")

set (MIPCL_PATH "${BUILD_DIR}" CACHE PATH "Path to MIPCL binary dir" FORCE)
file(APPEND "${CMAKE_BINARY_DIR}/options" " -DMIPCL_PATH=\"${MIPCL_PATH}\"")
message("MIPCL_PATH=\"${MIPCL_PATH}\"")
