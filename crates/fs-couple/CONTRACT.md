


## Free two-component translation on the original nonlinear ALE owner

With `experimental-stokes-motion`,
`stokes_motion::translating::planar::PlanarTranslationMotion` frees both world
translation components on the SAME original nonlinear translating flow owner.
Each actual native rigid impulse and its pressure/viscous reaction close one
coupled two-component endpoint. The entire original domain follows the body's
finite translation; reference operators and material transport are retained.
No fitted drag, apparent rigid mass or independent x/y fluid update is used.

A bounded local damped Newton search obtains the full cross-coordinate response
from two original flow probes. Authored velocity boxes and dimensional probe
sizes are explicit. The reusable native two-by-two LU owner solves the scaled
response; backtracking uses fixed-base physical residual scaling. Actual probe
velocity changes, including represented cross-coordinate rounding, determine
the correction. Every base, probe and rejected trial counts against the same
native-flow-solve budget. The local method does not guarantee root existence,
uniqueness or global convergence; original native nonlinear gates remain active.

Acceptance closes both ACTUALLY stored body momentum increments, and separately
checks represented body/mesh drift and mechanical work. Springs and dashpots
act in fixed world coordinates. Their stored energy, physical damping, signed
interface/external work and backward-Euler numerical loss remain separate;
numerical loss is not heat. Only pressure/viscous stress acts on the impermeable
body wall; advected momentum and kinetic energy retain their native accounts.
Original component-local pressure datums and material density remain binding.
Orientation is fixed, with the complete zero-work restraint couple about the
moving COM reported explicitly. Initial spin/out-of-plane momentum refuses.

The accepted native ALE candidate publishes once, together with body momentum
and pose, after final cancellation admission. Budgets, bounds, source/datum
refusal, work failure and cancellation preserve both continuing physical owners
and their shared clock. Exterior callbacks describe the translating outer
boundary at endpoint world coordinates; independently stationary obstacles,
rotation, remeshing and a whole-fluid energy/passivity theorem are not implied.

Four focused regression groups specify changing-direction momentum with finite
fluid inertia and common boosts; an independent off-diagonal world-force/spring
oracle and complete uniform-container energy balance; simultaneous open pressure
and tangential viscous recoil; and budget/datum/bounds/final-cancellation rollback.
New source and tests passed pinned Rustfmt parsing. Native compilation and test
execution remain unverified in the delivery environment.
