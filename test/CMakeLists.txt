find_package (Boost COMPONENTS unit_test_framework)
if(${Boost_FOUND})
  add_executable(qtest-bin "")
  set_target_properties(qtest-bin PROPERTIES OUTPUT_NAME qtest)

  target_compile_features(qtest-bin PUBLIC ${CPP_FEATURES})
#  target_sources(qtest-bin PRIVATE "${CMAKE_CURRENT_LIST_DIR}/main.cpp")

 # include("${CMAKE_CURRENT_LIST_DIR}/AIM/CMakeLists.txt")

  target_sources(qtest-bin
    PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/static_buffer_test.cpp"
    #"${CMAKE_CURRENT_LIST_DIR}/derivtest.cpp"
    # "${CMAKE_CURRENT_LIST_DIR}/static_test_rrfld.cpp"
 
    )

  add_subdirectory("extras/rapidcheck")
  add_subdirectory("extras/rapidcheck/extras/boost_test")

  target_link_libraries(qtest-bin PUBLIC quest ${Boost_LIBRARIES} rapidcheck rapidcheck_boost_test)
endif(${Boost_FOUND})
