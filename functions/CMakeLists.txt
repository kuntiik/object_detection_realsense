set(functions_source_files
  functions.cpp
  functions.hpp
  )

set(make_groups_souce_files
   make_groups.cpp
   make_groups.hpp
   )
set(bounding_box_source_files
    bounding_box.cpp
    bounding_box.hpp
    )
add_library(functions ${functions_source_files})
add_library(make_groups ${make_groups_souce_files})
add_library(bounding_box ${bounding_box_source_files})
