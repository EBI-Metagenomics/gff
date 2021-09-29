find_package(hope REQUIRED)

file(COPY assets DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(ASSETS ${CMAKE_CURRENT_BINARY_DIR}/assets)

set(TMPDIR ${CMAKE_CURRENT_BINARY_DIR}/tmp)
file(MAKE_DIRECTORY ${TMPDIR})

function(gff_add_test name srcs)
  add_executable(${name} ${srcs})
  target_link_libraries(${name} PRIVATE GFF::gff)
  target_link_libraries(${name} PRIVATE HOPE::hope)
  target_compile_options(${name} PRIVATE ${WARNING_FLAGS})
  target_compile_features(${name} PRIVATE c_std_11)
  add_test(NAME ${name} COMMAND ${name})
  target_compile_definitions(${name} PUBLIC "ASSETS=\"${ASSETS}\"")
  target_compile_definitions(${name} PUBLIC "TMPDIR=\"${TMPDIR}\"")
endfunction()

gff_add_test(gff_read gff_read.c)
gff_add_test(gff_write gff_write.c)
