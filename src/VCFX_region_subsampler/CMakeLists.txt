add_executable(VCFX_region_subsampler VCFX_region_subsampler.cpp VCFX_region_subsampler.h)
target_link_libraries(VCFX_region_subsampler PRIVATE vcfx_core)
