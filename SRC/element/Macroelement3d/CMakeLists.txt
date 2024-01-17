#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
  PRIVATE
    CohesiveSurface.cpp
    DamageShearInterface.cpp
	GambarottaLagomarsinoModel.cpp
	GenericDamagePlasticityShear.cpp
	Macroelement3d.cpp
	NoTensionSection3d.cpp
	WrappedMaterial.cpp
  PUBLIC
    CohesiveSurface.h
    DamageShearInterface.h
	GambarottaLagomarsinoModel.h
	GenericDamagePlasticityShear.h
	Macroelement3d.h
	NoTensionSection3d.h
	WrappedMaterial.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

