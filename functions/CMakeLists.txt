set(functions_source_files
  functions.cpp
  functions.hpp
  )

set(make_groups_souce_files
   make_groups.cpp
   make_groups.hpp
   )
add_library(functions ${functions_source_files})
add_library(make_groups ${make_groups_souce_files})
