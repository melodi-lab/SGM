FILE(REMOVE_RECURSE
  "CMakeFiles/create-docs"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/create-docs.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
