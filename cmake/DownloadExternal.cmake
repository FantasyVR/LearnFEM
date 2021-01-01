include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(EXTRA_OPTIONS "")
endif()


function(custom_download_project name)
    download_project(
            PROJ         ${name}
            SOURCE_DIR   ${LearnFEM_EXTERNAL}/${name}
            DOWNLOAD_DIR ${LearnFEM_EXTERNAL}/.cache/${name}
            QUIET
            ${EXTRA_OPTIONS}
            ${ARGN}
    )
endfunction()
################################################################################

## glad
function(download_glad)
    custom_download_project(glad
		GIT_REPOSITORY https://github.com/libigl/libigl-glad.git
		GIT_TAG        09b4969c56779f7ddf8e6176ec1873184aec890f
	)
endfunction()

# glfw
function(download_glfw)
    custom_download_project(glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.3
    )
endfunction()


# libigl
function(download_libigl)
    custom_download_project(libigl
            GIT_REPOSITORY https://github.com/libigl/libigl.git
            GIT_TAG        3cb4894eaf8ea4610467189ca292be349425d44b
            )
endfunction()

# eigen3
function(download_eigen3)
    custom_download_project(eigen3
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.3.9
            )
endfunction()
