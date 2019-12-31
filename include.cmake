include("3rdparty/find_libs.cmake")

message("MIPCL_PATH: ${MIPCL_PATH}")

set (MIPCL_INCLUDE_DIR ${MIPCL_PATH}/include CACHE INTERNAL "" FORCE)
message("MIPCL_INCLUDE_DIR: ${MIPCL_INCLUDE_DIR}")

set (MIPCL_LIBRARY_DIR ${MIPCL_PATH}/lib CACHE INTERNAL "" FORCE)
message("MIPCL_LIBRARY_DIR: ${MIPCL_LIBRARY_DIR}")

find_libs(MIPCL_DLL SHARED LIB_NAMES "mipcl" PREFIXES "lib" PATHS ${MIPCL_LIBRARY_DIR})
message("MIPCL_DLL: ${MIPCL_DLL}")
add_library(Mipcl SHARED IMPORTED)

if (WIN32)
    find_libs(MIPCL_LIB STATIC LIB_NAMES "mipcl" PREFIXES "lib" PATHS ${MIPCL_LIBRARY_DIR})
    message("MIPCL_LIB: ${MIPCL_LIB}")
    set_target_properties(Mipcl PROPERTIES
        IMPORTED_LOCATION ${MIPCL_DLL}
        IMPORTED_IMPLIB ${MIPCL_LIB}
        INTERFACE_INCLUDE_DIRECTORIES ${MIPCL_INCLUDE_DIR})
else()
    set_target_properties(Mipcl PROPERTIES
        IMPORTED_LOCATION ${MIPCL_DLL}
        INTERFACE_INCLUDE_DIRECTORIES ${MIPCL_INCLUDE_DIR})
endif()

function(installMipclDll dependingTarget)
    if(WIN32)
        add_custom_command(TARGET ${dependingTarget} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${MIPCL_DLL}
            $<TARGET_FILE_DIR:${dependingTarget}>)
    endif()
endfunction()
