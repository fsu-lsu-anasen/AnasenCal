add_library(CalDict SHARED)

target_include_directories(CalDict SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

ROOT_GENERATE_DICTIONARY(cal_dict DataStructs.h LINKDEF LinkDef_AnasenEvent.h MODULE CalDict)

target_sources(CalDict PRIVATE
    DataStructs.h
    DataStructs.cpp
)

target_link_libraries(CalDict ${ROOT_LIBRARIES})
set_target_properties(CalDict PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${ANASENCAL_LIBRARY_DIR})

add_custom_command(TARGET CalDict POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_BINARY_DIR}/libCalDict_rdict.pcm
    ${ANASENCAL_LIBRARY_DIR}/libCalDict_rdict.pcm
)