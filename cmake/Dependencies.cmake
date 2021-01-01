# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(DownloadExternal)



if(NOT TARGET igl)
    download_libigl()
    add_subdirectory(${LearnFEM_EXTERNAL}/libigl EXCLUDE_FROM_ALL)
endif()

#if(NOT TARGET glad)
#    download_glad()
#    add_subdirectory(${LearnFEM_EXTERNAL}/glad EXCLUDE_FROM_ALL)
#endif()

#if(NOT TARGET glfw)
#    download_glfw()
#    add_subdirectory(${LearnFEM_EXTERNAL}/glfw EXCLUDE_FROM_ALL)
#endif()

