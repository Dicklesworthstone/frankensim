
## Positive finite-interval tensor-Prony interference

The shared tensor-Prony/plane-stress spectrum window now retains bounded real
temporal field modes for eligible Uniform intervals. It combines original
binary-scaled constitutive channels before Kelvin/Maxwell heat squares and
independent stress:strain-rate work products. Centered coefficients are exposed
read-only by both shared owned/borrowed temporal factor APIs so constitutive
scales and original relaxation times need not first fit ordinary DD amplitudes.
Bounds and profile/parity requirements remain explicit.

For admitted 3..32 actual-real frequencies with maximum absolute half-width
at most0.5, at most17 modes replace the quadratic pair Gram cache. Preparation
and retained storage are O(F*M); material point work is O(F*M*(branches+1)*6),
with fixed point scratch and no query heap allocation. Original wider/count
fallbacks remain. The original spring increment is assembled before
DeltaU=Deltaq*(q_start+Deltaq/2), including equilibrium bias, instead of
subtracting bright stored endpoints. Endpoint reporting remains physical.
These are unweighted physical heat/work intervals; shaped detector means do
not redefine constitutive energy balance.

Validation: the original Kelvin three-tone rate-null regression returned
5.254869241218064e-46 for a centered1e-10s interval whose physical leading heat
is2.53125e-51. The new retained-field result and independently assembled
original-viscosity work agree with81*T^5/320. Eighteen focused actual-source
tests pass (eleven existing and seven new), including original extreme
viscosity/finite work, Maxwell memory below ordinary f64 range with its
original tau, anisotropic real ODE/quadrature closure with equilibrium bias,
precise source clocks, and shared owned/borrowed coefficient parity/bounds.
The shared temporal library and material projection/window test graph compile
with warnings denied. The full native workspace and public TensorProny
integration suite were not run. Finite-mode tail bounds exclude roundoff and
do not certify relative error at arbitrary cancellation nulls.
