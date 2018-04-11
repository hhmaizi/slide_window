
include_directories(
	./
	./factor/
	./utility/
)

# add add_executable
set(src_list vins_simple.cpp
# main.cpp
	estimator.cpp
	parameters.cpp
# 	vins_simple.cpp
	factor/pose_local_parameterization.cpp
	factor/stereo_projection_factor.cpp
	utility/utility.cpp
	utility/Converter.cc
)	

set(src_prjlist estimator.cpp
	parameters.cpp
	factor/pose_local_parameterization.cpp
	factor/stereo_projection_factor.cpp
	utility/utility.cpp
	utility/Converter.cc
)	

add_executable(vins_simple ${src_list}
)

target_link_libraries(vins_simple ${OpenCV_LIBS} ${CERES_LIBRARIES}) 

# endif (GFLAGS)
## example cmakefile in pose_graph_3d
# if (GFLAGS)
#   add_executable(pose_graph_3d pose_graph_3d.cc)
#   target_link_libraries(pose_graph_3d ceres ${GFLAGS_LIBRARIES})
# endif (GFLAGS)
add_subdirectory(tests)