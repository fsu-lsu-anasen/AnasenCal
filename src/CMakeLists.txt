add_subdirectory(CalDict)

add_executable(AnasenCal)

target_include_directories(AnasenCal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/CalDict SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS})

target_sources(AnasenCal PRIVATE
    ChannelMap.h
    ChannelMap.cpp
    DataCalibrator.h
    DataCalibrator.cpp
    DeadChannelMap.h
    DeadChannelMap.cpp
    EnergyCalibrator.h
    EnergyCalibrator.cpp
    GainMatcher.h
    GainMatcher.cpp
    main.cpp
    MapChecker.h
    MapChecker.cpp
    ParameterMap.h
    ParameterMap.cpp
)

target_link_libraries(AnasenCal PUBLIC CalDict ${ROOT_LIBRARIES})

set_target_properties(AnasenCal PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ANASENCAL_BINARY_DIR})