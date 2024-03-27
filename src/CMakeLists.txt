# BUILD LIBRARY 
file(GLOB SRC "*.cpp" ) 
include_directories($ENV{ANALYTICRT}/src)

ROOT_GENERATE_DICTIONARY(G__ANALYTICRT
    plot.hpp
    plotter.hpp
    trajectory.hpp
    isobar.hpp
    amplitude.hpp
    iterable.hpp
    data_set.hpp
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