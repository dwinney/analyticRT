# BUILD LIBRARY 
file(GLOB SRC "core/*.cpp" 
              "plotting/*.cpp") 

include_directories(core)
include_directories(plotting)

ROOT_GENERATE_DICTIONARY(G__ANALYTICRT
    data_set.hpp
    plot.hpp
    plotter.hpp
    trajectory.hpp
    LINKDEF ANALYTICRTLinkDef.h
)
add_library(ANALYTICRT SHARED ${SRC} G__ANALYTICRT.cxx )
target_link_libraries( ANALYTICRT ROOT::Core ROOT::MathMore )
target_link_libraries( ANALYTICRT Boost::boost)

install( TARGETS ANALYTICRT
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libANALYTICRT_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

# BUILD EXECUTABLE
add_executable( analyticRT "cling/analyticRT.cpp" )
target_link_libraries( analyticRT ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( analyticRT ANALYTICRT )