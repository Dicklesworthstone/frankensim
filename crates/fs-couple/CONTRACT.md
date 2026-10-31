

## Loaded thermal Voce member with physical rate dependence

ThermoplasticSolid composes selected affine-Q Voce hardening with linear or power-law Perzyna kinetics in either selection order. Each physical interval supplies one actual dt to the shared constitutive return; thermal substeps and fixed-point retries retain the incoming history. The six-tetrahedron held-member regression in tests/solid_plastic_enthalpy/voce.rs independently solves axial flow and caloric balance at n=1 and n=3, checks the full storage/physical/numerical energy split, and obtains the same accepted state with changed thermal substeps and enthalpy relaxation. Heating changes tension and lateral shape through the original spatial Newton and enthalpy owners. This is a bounded small-strain synthetic integration check, not empirical calibration, live-runtime cancellation qualification or finite-deformation plasticity.
