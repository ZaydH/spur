add_subdirectory (component_types)

SET(SOURCES
        cached_assignment.h
        component_cache.inc
        containers.h
        primitive_types.h
        sampler_tools.h
        solver_config.h
        structures.h
        alt_component_analyzer.cpp
        component_cache.cpp
        component_management.cpp
        instance.cpp
        model_sampler.cpp
        rand_distributions.cpp
        solver.cpp
        stack.cpp
        statistics.cpp
    )


ADD_LIBRARY(libSpur ${SOURCES} component_types)
target_link_libraries (libSpur component_types)
