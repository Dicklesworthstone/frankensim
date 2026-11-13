# fs-solid CONTRACT

## Native directional material tangents and local diagonals

Prepared finite-strain tetrahedral tangent actions call each original law's
`consistent_tangent_action`. Native hyperelastic laws evaluate a directional
nested-dual derivative without forming a full Hessian. Adapter defaults retain
the same owner's complete consistent tangent and held history. Reference
changes transform increments linearly before the original-volume Piola
pullback, retaining directions below the identity's ULP.

Retained local Hessians expose a signed diagonal without global or element
stiffness assembly or further constitutive calls. Supports remain present;
caller-owned scratch is reusable and output publishes only after complete
finite-value and cancellation checks. The diagonal includes the original
finite-strain geometric contribution already present in the Piola Hessian.
This changes evaluation/storage, not the material law or integration scheme.
Focused tangent-action/reference and native diagonal regressions cover these
contracts; no whole-solver speedup or cross-ISA timing claim follows.

## Progressive ductility in native fiber-frame equilibrium

`FiberLaw::DuctileJ2` and its owned history dispatch the original progressive
axial material response. Authored and original source-resolved Johnson--Cook
or tabulated initiation retain one physical duration and temperature.
`FiberMaterial::with_ductile_evolution` requires the exact original native J2
backbone and matching source card, fixed context and actual temperature;
reference density and original mechanical receipts survive. Supplied work
capacity, attenuation profile and receiving length remain authored inputs.

`frame2d::ductile_damage::FrameDuctileDamagePath` owns complete native mechanics,
progressive histories and a physical clock. It supports original load control,
bordered displacement control and Newmark motion. Prepared durations equal
each accepted physical interval; Newton and element iterations remain trials
from the same committed history. Other timed material families refuse.
The original untimed public entry points continue to reject timed fibers.

Each progressive band's length must equal its original Lobatto tributary
length `w*L`, with volume `area*w*L`. This establishes the local work-capacity
normalization, not a mesh-objective crack trajectory or localization theorem.
The explicit one-dimensional law degrades both tension and compression.
Original stress and full possibly negative tangent enter section equilibrium.
Progressive sections use a row-scaled actual inverse without the legacy
determinant floor. A singular failed section refuses; no residual stiffness
or automatic member removal is fabricated. A failed fiber can remain with
zero response while surviving fibers carry the section.

Accepted local increments retain their original material energy, so observation
requires no second return map. Fixed-order compensated sums separate stored
energy, plastic dissipation, damage release and numerical loss. These accounts
cover progressive fibers only; they are not extra terms to add to the native
dynamic work balance and do not deposit numerical loss as heat.
Late domain, mechanism, budget or cancellation failure rolls back the entire
proposed block, including earlier increments, inertia and physical time.

The four `fiber_ductile_evolution` and five `frame_ductile_damage` regression
cases exercise assembled stress/tangents, exact local completion, surviving
load paths, source handoff, physical rate, work, native motion, restart and
transactional refusals. Native execution remains subject to the real workspace
build; these cases alone assert no experimental material or component validation.
Finite-strain necking, crack closure/contact, separated frame topology and
automatic transfer of a singular member remain outside this implementation.

## Ductile initiation on native frame fibers and sourced tetrahedra

`frame2d::ductile::FrameDuctilePath` owns the original `Frame2d`, explicit
member/section/fiber calibrations, accepted initiation histories and physical
clock. `static_load`, `displacement_control` and `newmark_with_load` reuse the
native observed-step owners; every accepted physical increment is evaluated.

Each point consumes its original accumulated equivalent plastic strain and
signed physical axial stress, including reversed plastic flow. Section
resultants, total strain and backstress-shifted effective stress are never
substituted. Static increments require explicit physical duration; Newmark
inertia and the initiation strain rate share exactly the same `dt`.

Supported fibers are native rate-independent J2, Voce, tabulated hardening and
Chaboche with one through eight independent branches. Menegotto--Pinto,
Mander and timed/thermal laws cannot acquire invented plastic accumulation.
Virgin initiation cannot be attached to an already plastically loaded point.

Source calibrations require the caller's original `FiberMaterial` handoff.
Its coefficients must match the native fiber, and its card, state, fixed
context and typed absolute temperature must match the initiation calibration.
Raw frame fibers retain coefficients rather than card identity; equal
coefficients cannot reconstruct source provenance that the caller discarded.

A late point, step, source, domain, budget or cancellation refusal preserves
the complete input frame, velocity, acceleration, initiation states and clock.
The reported initiation candidate minimizes the endpoint-held plastic-increment
fraction deterministically; fractions from different fibers do not locate the
first physical event. The original interval is only an event bracket.

`DuctileHistory` also admits source-resolved Johnson--Cook initiation through
its existing J2 tetrahedral/source and thermoplastic candidate interfaces.
The actual criterion temperature must equal the native endpoint's resolved
source temperature. Temperature-dependent calibration alone cannot qualify a
mechanical law frozen at another temperature; original endpoint bindings remain.

The regression cases in `frame_ductile` and `ductile_source` exercise independent
Prager reversals, local bending despite cancelling resultants, rate/clock
semantics, native motion/work preservation, replay and transactional refusals.

Frame gauges qualify only their selected points at fixed temperatures under
native small-axial-strain assumptions; the declared strain bound is at most
0.05. Endpoint quadrature supplies initiation only, without necking, progressive
softening, finite-strain forming, fracture energy or separated members.

## Two-sided normal contact on native ductile fragments

`FiniteDuctileFractureMotion::advance_with_face_contact` resolves frictionless dwell or unloaded separation on every non-tied vertex pair of the retained released faces. The existing bounded native normal-complementarity solver probes the original finite backward-Euler owner with equal/opposite nodal loads. All probes borrow the same accepted material histories, consistent mass, velocity, temperature and physical duration; only the original final response matching the solver's force/path receipt can be published.

All participating rows must fit the native solver's maximum of 21 contacts and the caller's stricter budget. Numerical force seeds and natural-map scales neither evolve history nor add compliance. Final forces must be nonnegative with exactly zero tangential force. Every paired initial/final gap must be nonnegative. Loaded pairs must stay matched, active and normally at rest; incoming closing impacts and loaded finite sliding return typed refusals. Both oriented triangle normals must remain in the declared current release-normal cone.

The unchanged native tetrahedral orientation-sweep owner checks the straight nodal chord against element collapse under an explicit work budget. Nonnegative endpoint gaps bound the matching fixed-normal P1 trace over that chord. This does not certify unrelated surface collisions, arbitrary fragment rotation, finite sliding, friction, impact restitution, or exact continuous trajectories.

Actual pair-force work, its signed negative endpoint defect and summed absolute defect are reported separately from physical plastic/damage loss, native material defects and backward-Euler numerical kinetic loss. Original external dead-load work is independently integrated and checked against the native contact-inclusive load receipt. The complete external/support/material/inertia/contact identity is gated. Summed absolute endpoint force couples bound angular-momentum error without cancellation between pairs. Contact defects are never deposited as heat or fresh fracture energy.

Three focused G0/G3 tests in `tests/finite_ductile_fracture/contact.rs` exercise rotated, genuinely failed native specimens under exact equal/opposite contact loading; mass-driven free separation with independent backward-Euler kinetic/work formulas; and transactional loaded-sliding refusal. The transferred history, original mass and separate work accounts remain observable. Pinned rustfmt parsing and bounded static API/algebra review completed; Cargo typechecking and test execution remain unverified because of the workspace storage and dependency-build memory limits.

## Finite ductile failure to native deformable fragments

`hyper3::plastic::ductile::fracture` releases explicitly selected conforming P1 seams only after both incident original material points have exact finite-support completion, zero retained compressive storage, compatible original reference-band direction, and traction-free current facets. Current Cauchy traction and area-scaled first-Piola traction are admitted independently. A rounded exponential damage value cannot authorize a cut.

The zero-duration transfer duplicates original nodal stars. Element/local-corner order, deformation, multiplicative plastic history, original material owner, characteristic length, physical stress, and cumulative plastic, damage, and numerical-work accounts are copied unchanged. The cut creates no additional fracture-energy credit, impulse, or material return. Boundary conditions and loads must use the receiving node numbering.

Both cellwise density and the complete P1 nodal-density field retain their original consistent mass. Exact copied motion is checked by element-order mass, current first moment, linear momentum, angular momentum, and kinetic inventories; native CSR kinetic roundoff has a separate energy gate. Graded density is not replaced by cell means or a lumped/rigid-fragment mass. Released faces retain their current release normals and element/corner identities through later explicit cuts.

`FiniteDuctileFractureMotion` continues the original finite backward-Euler solve and checks paired openings and oriented current normals at each accepted endpoint. Contact-free continuation refuses closing endpoints and departure from the admitted normal cone. These are endpoint guards, not a swept-path or general collision proof. There is no automatic crack discovery, arbitrary sliding/rotating-fragment contact, erosion, general remeshing, or calibrated alloy failure claim.

Damage-increment continuation remains necessary to resolve the softening work: a permitted virgin-to-terminal jump can exhaust the attenuation capacity while its endpoint work and defect are both zero. Endpoint-defect control alone cannot detect that unresolved path. The transfer neither manufactures missing physical work nor reclassifies the attenuation diagnostic as heat or fracture loss.

Three focused G0/G3 tests in `tests/finite_ductile_fracture.rs` construct a stress-free rigid rotation followed by a resolved progressive loading path, then check both density modes, exact history/state transfer, native fragment motion, and invalid-release/closure/cancellation refusals. Source and test files were parsed by pinned rustfmt; Cargo typechecking and execution remain unverified because the workspace is full and the real dependency build hit the memory limit.

## Progressive ductile rigid-plane tooling

Static and moving `FiniteDuctileTetSolid` / `FiniteDuctileTetMotion` endpoints
now expose `advance_with_plane_contact`. The shared
`hyper3::plastic::plane_contact` owner reuses native sustained normal
complementarity, original force probes, all-vertex plane gaps and the original
tetrahedron orientation sweep. Thermal plastic tooling delegates to the same
owner while retaining its existing public names and complete energy arithmetic.

Every force probe starts from the same accepted deformation, velocity,
plastic/damage histories and physical duration. The actual native nonsymmetric
equilibrium or momentum solve supplies its response. Only the final response
matching the original contact receipt can be published. Original reference mass,
material identities and quadrature remain unchanged. Explicit candidates are
bounded by the native contact solver; every original vertex must remain outside
every tool, so omitted candidates cannot silently permit penetration.

The native material accounts retain storage, plastic loss, damage release and
signed endpoint-work defects. Moving specimens additionally retain physical
kinetic change and backward-Euler numerical loss. Tool work, body contact work,
original dead-load work and the signed normal endpoint defect are separate;
the sum of absolute contact defects has its own gate and never becomes heat.
A rejected material/contact/geometry/work/cancellation proposal publishes no
partial specimen or reaction state.

Three focused source tests cover independent compressive finite stress and
damage, rotation, holds and tool release; analytical rate-dependent loading and
held relaxation with original consistent mass; distinct contact-arrival,
material and kinetic accounts; omitted penetrating vertices; and transactional
refusals. Parsing/formatting and bounded compatibility review passed. Runtime
and type-check execution remain unverified in the exhausted shared workspace.
Virgin native histories have no previously stored damage energy to release;
subsequent damage growth releases the separately retained active storage.

This is frictionless fixed-normal translating-plane tooling, not self-contact,
arbitrary curved tools, an impact/impulse law, crack creation, erosion or ductile
fragment transfer. Straight-path gap/orientation admission is not an exact-flow
CCD or temporal-error estimate. Unstable force-controlled branches can still
require physical refinement, inertia or another continuation formulation.

## Progressive ductility on original finite H8 cells

`FiniteDuctileHexSolid` binds eight independent original `FiniteDuctileJ2`
histories and explicit immutable characteristic lengths per cell. Standard and
MeanDilatation formulations reuse the original PreparedHex reference geometry,
compatible/Fbar kinematics, and complete material/geometric/cross-quadrature
force pullback. The thermal wrapper retains its own scalar potential. The
ductile owner exposes force and the complete nonsymmetric derivative, and uses
the native stress-only residual path. Physical corner, centre and quadrature
Jacobian checks remain independent of mean dilatation.

One actual interval drives native overstress and initiation at every point.
Uniform, per-point and per-cell held temperatures preserve coverage and material
domain checks. Acceptance verifies original committed histories, stored energy,
plastic dissipation, damage release, support/dead-load work, and separate native
material and projection endpoint defects. The defect gate bounds the sum of
absolute volume-weighted point material defects plus absolute cell projection
defects; signed cancellation cannot hide an unresolved increment. A refusal
preserves the original specimen.

Characteristic lengths are authored model inputs bound to the original
mesh/material assignment. Arbitrary-mesh fracture objectivity, topology/erosion,
contact handoff, creep coupling, heat transport, arclength continuation and
experimental qualification are outside this component. Physical Jacobian
sampling does not prove global injectivity of a warped cell.

Five focused source tests cover affine analytic finite response and viscous
loading/holds in both formulations, the complete nonaffine 24-by-24 tangent,
independent point calibration, physical/projection work and noncancelling
refusal, coverage, point limits, corner folding and cancellation. Pinned
rustfmt parsing/formatting, scoped diff checks and bounded static review passed.
Cargo type checking and runtime tests remain unverified: earlier real builds
were killed in asupersync, and the shared workspace now has exhausted storage
and insufficient memory for another compile. No executed numerical claim is
made for this extension.

## Progressive finite ductile sheet motion

`hyper2d::plastic::ductile::dynamics::FiniteDuctileSheetMotion` advances the
original P1/Q1 progressive ductile sheet through the native backward-Euler
momentum owner. One positive physical interval drives inertia, rate-dependent
plastic flow and initiation; trial/Newton retries retain the same old history.
This is two-component in-plane membrane inertia. Transverse shear and thickness
are eliminated through the full quasistatic plane-stress solve, without
independent thickness or rotary inertia.

Each element's original thickness and volumetric density produce one retained
consistent reference mass. Neither damage nor the solved current thickness
removes mass. Acceleration loads use that same mass action; authored body forces
integrate the original reference shape functions and volumes. The existing
uniform-thickness path retains its accumulation expressions. The shared
increment owner adds inertial residual forces and support reactions once; the
sheet operator adds only the matching inertial tangent.

Complete material point histories and receipts, physical plastic dissipation,
released damaged storage, signed material endpoint-work defects, kinetic change,
and backward-Euler numerical kinetic loss remain separate. Independent loss and
work gates precede publication of the complete accepted motion. Material,
equilibrium, work-budget or cancellation refusals leave the prior sheet,
velocity, reference inventory and clock unchanged.

Focused source tests cover heterogeneous thickness/density mass and moment,
free acceleration, independent affine kinetic work during rate-dependent
softening and held relaxation, a loaded free deformation degree of freedom,
transactional refusals, and uniform-versus-elementwise mass parity. Formatting
and parser checks passed. Execution remains unverified: both offline Cargo
checks were terminated by SIGKILL during the real asupersync dependency, including
a one-job, one-compiler-thread retry, before compiling these changes. No passing
Rust test result is claimed.

This membrane path does not supply contact, impulses, mass transfer, independent
out-of-plane bending motion, thermal evolution, topological separation, or a
temporal error estimator. The directional crack-band and local regular-branch
limits of the progressive sheet/material remain in force.

## Progressive ductile native shell endpoints

The existing CST/DKT owner provides `solve_ductile_with_energy` and
`advance_ductile_backward_euler_with_energy` for original `DuctileJ2Law`
sections. Native surface/thickness quadrature, transverse equilibrium,
nonsymmetric damaged tangents, support reactions, reference drilling stiffness
and consistent mass are reused. One complete endpoint publishes only after
separate material, drilling, external-work and kinetic balances pass.

Prepared layers use the same physical duration; the dynamic duration matches
inertia exactly. Rejected damage increments, owner/time/gauge mismatches, work
defects and cancellation leave input checkpoints unchanged. Caller physical
cutback refreshes prepared intervals while retaining material/history ownership
and band length. Plastic loss, damage release, numerical material loss and
drilling/inertial numerical losses remain distinct; no heat partition is inferred.

Five focused cases in `ductile_shell_evolution` cover layer bending/softening,
unloading/reversal/restart, backward-Euler work partitions, successful caller
cutback and transactional refusal. This remains small-strain/small-rotation
shell mechanics. Finite forming, delamination, fracture topology, erosion,
generated contact, time-accuracy estimation and measured-alloy qualification
are outside its claim.

## Progressive finite ductile sheets

`hyper2d::plastic::ductile::FiniteDuctileSheet` advances the original P1/Q1
reference continuum beyond initiation using full three-dimensional progressive
ductile history, solved transverse tractions and the nonsymmetric condensed
tangent. The native sheet reference geometry, quadrature, force/tangent assembly
and shared Newton/FGMRES increment owner are reused. No scalar damage potential
or symmetric stiffness is substituted.

Each immutable `SheetCrackBandGeometry` binds the selected in-plane reference
normal, polygon and width (reference area divided by projected transverse span).
The original element thickness supplies volume, never the in-plane band width.
Every material length must match its retained geometry. This is an explicit
projected-footprint regularization, without arbitrary-mesh crack-path or
localization-objectivity claims.

Original point histories, actual point temperatures and one physical duration
remain frozen across solver iterations and retries. Accepted physical time
advances once. Complete material storage, physical plastic loss, damage release,
signed endpoint-work error and transverse residual work remain separate. Work
gates bound absolute point defects before publishing the complete endpoint.
Displacement control can traverse falling load branches; force-controlled
instability can refuse without consuming history.

Six G1/G3 cases in `finite_ductile_sheet` cover independent P1/Q1 rotated
biaxial softening/thickness/work, damaged free-DOF equilibrium, full assembled
tangent differences, unloading/replay, rate-dependent holds and geometric/
work/time/cancellation refusals. Bending, wrinkling, evolving thermal
coefficients, separation and contact are not supplied by this in-plane owner.

## Rate-dependent progressive finite solid damage

The finite ductile tetrahedral static, backward-Euler and adaptive path owners
preserve the explicitly selected native Perzyna kinetics. Each trial receives
its actual physical subinterval; retained overstressed states are validated
through the complete ductile owner. Adaptive initial-equilibrium inspection
uses the retained accepted Piola stress and consumes no relaxation interval.
Held deformation can therefore continue plastic flow and ductile degradation.

The existing finite ductile caloric coupling receives the same native physical
plastic loss and damage release through its solid owner. Reference inertia,
numerical kinetic loss and signed endpoint-work defects keep their separate
accounts. `finite_ductile::rate` adds force-controlled hold/replay and dynamic
continuation regressions, including exact replay of accepted physical
subintervals. Source-calibrated viscosity, evolving thermal kinetics and
post-failure component/contact transfer remain outside this selection.

## Pointwise ductile initiation in layered plastic shells

`ductile::shell` gates the native CST/DKT quasistatic and backward-Euler
solvers at every facet, surface and thickness material point. Source-bound
J2, Voce and Chaboche histories use read-only recovery of full physical stress
and retain equivalent plastic accumulation through reversal. Admission requires
the original material-card/state/context binding and the actual typed fixed
source temperature. The physical step duration drives both native inertia
and initiation-rate evaluation.

All mechanical, velocity and initiation histories publish atomically.
Membrane and bending resultants never replace local stress. Uniform fields
admit each original physical layer once; heterogeneous fields check the
calibration at each original point. The native element/material-point budgets
and cancellation checks remain active. These adapters supply initiation gates;
they do not add damage-softened stiffness, changing-temperature constitutive
response, crack topology or finite-strain shell forming.

G0/G1/G3 regressions in `ductile_shell` cover opposing plastic bending layers
with cancelling membrane forces, all 18 material points of a layered facet,
unload/reversal/restart, atomic refusal and refinement, heterogeneous source
bindings, typed temperature mismatch and native backward-Euler identity/rate.

## Ductile initiation in finite plastic sheet forming

`DuctileHistory::advance_finite_sheet` and `advance_thermoplastic_sheet` gate
the original P1/Q1 sheet endpoint at every material point. The complete native
Piola/deformation pair supplies spatial stress, including solved transverse
stretch and shear. Neither in-plane Piola components nor element-average
resultants substitute for triaxiality or the Lode invariant. The original
plastic metric, accumulated plastic strain, thickness and energy accounts
remain with the sheet solver; observation performs no second return map.

Explicit Perzyna and simultaneous Norton/J2 sheets retain their actual physical
interval. A mismatched ductility interval refuses, and independent creep
accumulation never consumes plastic ductility. Pure-creep thermal families
refuse before solving. Thermal sheets use actual endpoint quadrature
temperatures, not element means. Initiation, material/domain, clock, history or
cancellation refusals publish neither a mechanical nor an initiation endpoint.

The synthetic G1/G3 tests in `ductile_sheet` check independent equibiaxial
Hencky stress/thickness and initiation, rotated P1/Q1 specimens, rate-sensitive
flow and held-deformation relaxation, creep/plastic separation, nonuniform
point temperatures and transactional retries. This is the original local,
endpoint-held initiation criterion, not a forming-limit/necking solution,
progressive sheet fracture or experimental qualification of an alloy.

## Source-refreshed thermal strain with prescribed support dynamics

`TetTransientRun3::advance_thermal_with_supports_candidate` and
`advance_thermal_with_support_accelerations_candidate` join changing thermal
eigenstrain to the existing full consistent-mass support/Newmark equations.
They retain the same fixed components and reference inertia; refreshed source
identities, thermal fields and complete stiffness remain bound to the returned
run. Changing stiffness rebuilds the native reduced operator and retains full
free/support coupling, including consistent mass and both endpoint stiffnesses.
Stiffness-proportional damping uses the positive midpoint operator
`eta*M + zeta*(K_old+K_next)/2`. Its support work uses the reaction conjugate to
that original midpoint equation; physical endpoint acceleration and reactions
are separately recovered for observation and subsequent advancement.

The full endpoint thermal load enters reactions and the free equation.
Numerical work pairs the mean of accepted and candidate thermal forces with the
same full displacement as external and support work. Physical shifted elastic
storage separately balances signed thermal material work, external work,
support work, signed tensor-parameter work and declared Rayleigh loss. Full
coefficient work differences the original material tensors before strain
products, including fully prescribed bodies, without assembling a full matrix.
No reciprocal heat debit is inferred.
Ordinary later support steps consume the newly accepted field and reactions.

Both endpoint conventions share native support admission and predictors.
Cancellation, source, strain, displacement, energy and resource refusals leave
the input checkpoint intact. Focused regressions in `transient_supports3`
cover an independent damped free-component equation, heating and cooling,
deterministic replay, fully prescribed motion, retained thermal reaction/work
and refused candidates. Changing-tensor regressions add an independent damped
two-coordinate equation, softening/stiffening, ordinary-step continuation and
fully constrained parameter work at zero boundary displacement.

## Thermal strain rates in native viscous tensile fracture

`ThermalBrittleTetMotion` advances total deformation and source-integrated free
strain through the original timed crack-band adapters. Both native backward
Euler and midpoint retain separate parameter work, total-strain work,
eigenstrain work, intrinsic fracture and physical viscous loss. The total mean
stress enters both strain work ports; viscous loss is subtracted once from the
complete physical balance and never added to cumulative fracture history.

Every interval explicitly declares affine free strain. Actual endpoint stress
uses the native P1 projection of accepted velocity minus that imposed rate,
including the distinct centered endpoint velocity. Temperature endpoints do
not supply an inferred instantaneous derivative. The explicit initial-rate
constructor admits thermal kinematics without advancing work, history or time;
the existing constructor declares held initial temperature.

Accepted scalar free-strain rates survive zero-time normal/friction impulses
and completed-fracture topology transfer. Post-impulse stress is readmitted at
the resulting actual velocity; original reference mass stays frozen. Sustained
and contact-time balances carry physical viscous loss. Untimed static evolution
and elastic-only prehistory refuse surviving viscous bonds. Source, rate,
energy, geometry and cancellation gates publish no partial successor.

Focused native regressions cover matched free expansion, restrained heating and
cooling, independent inertial solutions for both schemes, endpoint versus chord
rates, rupture/reclosure, source/refusal atomicity and retained thermal rates
through fracture/contact. The coefficient-first thermal split remains explicit;
there is no inferred heat conversion, caloric feedback or unresolved temperature
trajectory guarantee.

## Physical fracture viscosity in native tetrahedral dynamics

`BrittleTetMotion` endpoint and midpoint advances use the original duration-bound
crack-band/cohesive law when explicit fracture-plane viscosity is selected.
Total material work includes the exact surviving-bond viscous integral, while
`viscous_dissipation_j`, intrinsic fracture history, contact loss and numerical
kinetic/endpoint defects remain separate. Adaptive, sustained-contact and
coefficient-temperature work accounts carry the same physical viscous loss.

The native P1 operator projects accepted nodal velocity to instantaneous stress
after initialization, centered steps and contact impulses. A secant-rate path
traction is not published as a midpoint endpoint-velocity observation. Source,
rate, energy, budget and cancellation refusals publish no partial motion.
Complete failure removes bond viscosity and retains native compression and
topology/momentum transfer. Elastic-only tensile predecessors and the untimed
thermal-eigenstrain owner cannot silently omit selected viscosity.

Coefficients remain explicit Pa s/m inputs without inferred bulk calibration
or arbitrary-width objectivity. Focused native regressions use independent
normal/shear oscillator and surviving-bond work oracles, checking actual
endpoint stresses, distinct losses, rupture, reclosure and refusal atomicity.

## Regional source temperatures in retained finite-strain membrane motion

`OwnedFiniteViscoMembraneMotion::with_temperatures` evaluates one original
engineering-length source per triangle at its prescribed absolute temperature.
Every endpoint keeps the original datum and fixed source frame, including
cooling and revisits. Exact per-triangle source/card coverage is required;
unsupported temperatures and missing or foreign paths refuse.

Preparation copies the original patch once and replaces only its sheet
endpoints and cached constitutive laws. Native geometry, gradients, topology,
grain, original corner density, mass and strain envelope are not rebuilt.
The native viscous laws retain their equilibrium coefficients, branch energies,
kinetics and solver controls; only the admitted natural map changes.

The existing owned reference transaction retains displacement, velocity, time
and viscous metrics while restoring transverse balance. Subsequent native
advancement resolves the changed global forces and physical relaxation. Signed
reference input, pressure-volume actuation, numerical defects and relaxation
heat retain separate accounts. Prescribed temperatures do not constitute a
reciprocal caloric evolution.

`tests/visco_membrane_temperature.rs` covers constrained regional heating and
subsequent relaxation, both volume treatments, graded inventory, original-datum
revisits and source/work/cancellation atomicity.

## Spatial indentation hardness from accepted J2 history

`SmallStrainTetProblem<ResolvedJ2Plasticity>::observe_hardness` evaluates each
accepted element through its independently sourced `ResolvedJ2Hardness` path.
Uniform and heterogeneous material fields use the original Newton solution,
accumulated plastic strain and exact mechanical sources. Hardness is neither
averaged across regions nor inferred from yield stress. Native constitutive
admission and the prepared hardness interval are checked before a complete
field is returned; late refusal or cancellation leaves all inputs unchanged.
The caller accepts its native mechanical candidate and hardness field together.

The observer allocates one scalar per element, borrows source receipts and
performs no second return map. Typed endpoints remain available from the original
material observer. It does not attest an arbitrary solution's origin, certify
unobserved intermediate mechanical states, evolve contact geometry, or add heat.
`tests/plastic3_hardness.rs` follows loading, unloading and reversals in a bonded
two-material native solid, exports original hardness tests and passes the local
field to the existing Archard owner without averaging hardness. Prescribed
surface pressure/slip quadrature is a coupling fixture, not a contact solve or
measured material qualification. Late source/range refusal remains atomic.

## Spatial indentation hardness from complete cyclic history

`SmallStrainTetProblem::observe_cyclic_hardness` observes the original prepared
hardness source on each native source-bound Voce or Chaboche element. Uniform
and heterogeneous fields preserve the accepted Newton solution, all independent
backstress branches and their exact named source order. The material observer
checks complete native history, actual nonlinear yield and original support;
no additional return map or per-element source cloning occurs. Only a complete
scalar field is returned, with indexed refusal or cancellation leaving inputs
unchanged. The existing J2 observer remains available with its original API.

`tests/plastic3_hardness_cyclic.rs` continues two separately sourced solid
coupons through nonproportional loading and reversal under physical boundary
tractions, checks the native nonsymmetric solve and complete branch histories,
and observes each local hardness and original indentation test. A native Voce
solve checks the uniform field. Late source, support and branch-state failures
publish no partial field. The fixtures do not certify intermediate trajectories,
coupled thermal history, measured cyclic-hardness accuracy or fatigue life.

The additional `export_hardness_field` and `export_cyclic_hardness_field` methods
return the material owner's compact `AcceptedHardnessField` in native element
order. They pass each ACTUAL element law alongside its accepted strain/history,
retain one fixed original source and poll the native Cx during construction.
Different original element laws refuse instead of becoming a homogeneous
half-space source. Every point keeps its own alpha/H with no averaging, native
history cloning or per-element source copy. Existing J2/cyclic test targets
exercise real multi-element exports, complete branch retention, late history
refusal and heterogeneous-law mismatch. This is frozen external/preconditioning
hardness transfer: the consumer supplies an explicit original cell mapping and
its original contact geometry, without transferring residual stress, plastic
displacements or duplicating the bulk indentation deformation/work.

## Original indentation hardness at native Hex8 points

`SmallStrainHexProblem::observe_hardness` and `observe_cyclic_hardness`
extend the original J2 and source-bound Voce/Chaboche observers to all eight
material points per brick. Uniform, per-element and per-quadrature hardness
sources remain explicit. Each actual law, accepted assumed strain and complete
history is admitted in native element/Gauss order, including the B-bar strain
used by the mechanical solver. No element-average hardness, extra return map or
per-point source-receipt clone is introduced; scalar output costs one f64 per
point. Allocation budgets, source/support mismatch and cancellation refuse the
complete field atomically.

`export_hardness_field` and `export_cyclic_hardness_field` retain every point's
own alpha/H in the original compact `AcceptedHardnessField`. A shared source is
stored once; different actual original laws refuse homogeneous export. External
contact still supplies its explicit point map and original geometry. Residual
stress, displacement, bulk work and thermal energy are not transferred.

`tests/plastic_hex_hardness.rs` covers Standard/B-bar quadrature, heterogeneous
laws, an independent first-load radial-return result, elastic unloading,
reversal, complete cyclic branches and late atomic refusals. Its Archard fixture
uses each local hardness and shows why averaging H changes removal; it does not
qualify measured material behavior or solve reciprocal bulk/contact loading.

## Original thermal hardness at native Hex8 points

`ThermalHexProblem::bind_hardness_field` pins the original mesh, formulation,
source-family order and one exact indentation test per cell.
`ThermalHexHardnessBinding::at_accepted` observes every actual caloric point
temperature, accepted assumed strain and complete native J2/Voce/Chaboche
history. The existing thermal family admits the complete state, and the
original hardness source evaluates its supported H(T, alpha) without another
constitutive solve. Source resolution retains the original mechanical,
caloric and hardness selections. The typed field preserves receipts and native
element/Gauss order; all points must succeed before publication.

Observation changes no native history, stress, heat or energy account. It
neither averages point temperatures nor substitutes a yield/hardness
correlation. The full typed endpoints allocate source receipts. Unsupported
temperature/history ranges, inconsistent stored temperatures, changed owners,
budgets and cancellation refuse atomically. The binding covers accepted
endpoints, not intermediate trajectories, annealing or unimplemented material
families.

`tests/plastic_hex_thermal_hardness.rs` exercises independent point temperatures
and plastic histories under loading, heating/cooling and reversal, Standard/B-bar
formulations, complete Chaboche branches, exact test selection, unchanged native
energy and late-point refusals. These synthetic fixtures verify implementation
and source use; they do not establish measured thermal-hardness accuracy.

## Owned finite-strain membrane reference continuation

`OwnedFiniteViscoMembraneMotion::change_reference` delegates to the original
native reference transaction, then retains its complete accepted target patch,
laws and histories without internal borrows. Original coordinates, topology,
source/model, material axes, thickness, corner inventory and volume treatment
remain subject to the native admission. Caller-supplied patch storage is shared;
uniform material fields remain uniform and accepted displacement/velocity
buffers survive.

The native signed reference input, transverse loading, exact-volume constraint
work and absolute endpoint-defect gate remain distinct. No elapsed time, heat,
impulse or kinetic-energy change is created. Retained viscous metrics continue
through the ordinary native motion solver after temporary target law storage
has gone away. A late cancellation publishes no candidate.

`tests/visco_membrane_owned_reference.rs` covers regional and uniform reference
events after nonproportional history, both volume treatments and mass images,
graded original inventory, native receipt parity, subsequent motion and refusal
atomicity. This ownership adapter supplies no reciprocal thermal closure.

## Retained thermal porous compaction

`hyper3::plastic::porous::thermal::FiniteThermalPorousTetSolid` carries the original pressure-sensitive plastic metric, accumulated flow, total deformation and temperature across native total-Lagrangian P1 tetrahedral equilibria. Each element retains its complete thermal family or original source binding. The old history is admitted at its accepted temperature before the endpoint-temperature return, so temperature-induced softening can cause further compaction under retained loads. The native total-Piola tangent includes thermal stretch.

Original reference geometry and mass persist through thermal expansion and plastic volume change. Current density is explicitly `rho0 / det(F)`; plastic volume is separately available as `Jp`. Optional backward-Euler inertia reuses the original consistent P1 mass operator, including an explicitly supplied nodal density field. Numerical kinetic loss is admitted separately from physical caloric heat. Source densities are not inferred by the mechanical owner.

The shared finite thermoplastic energy carrier separates yield heat, signed reversible caloric exchange, free/internal energy storage, coefficient storage, support/dead-load work, and signed endpoint work defects. The complete mechanical proposal is transactional; temperature/history, energy, late inertia and cancellation refusals leave the accepted owner unchanged. The `hyper_porous_thermal` regressions specify free thermal expansion, independently solved force-driven hot compaction and unloading, exact graded-inertia accounts, late rejection and retry.

This is a dry associated porous law on a compressible displacement-only P1 solid. Its original volumetric-locking, constitutive calibration and time-discretization limits remain. Supplied temperatures belong to an external caloric owner; numerical work defects are never relabelled as physical heat.

## Native thermal tensile motion and conservative fracture replay

`expansion::initiation::dynamics::ThermalTensileElasticTetMotion` continues
accepted elastic thermal prehistory with explicit velocity, physical time and
step budget. Density from the same original selected source bundle is bound
once to reference geometry and cached as consistent P1 mass. Later current
density queries cannot replace that inventory or supply an impulse.

`probe_midpoint` uses the existing centered inertia and nonlinear P1 owner.
A borrowed original elastic response evaluates C_next times the mean of OLD
mechanical strain and NEW mechanical strain, with analytic endpoint tangent
C_next/2. Integrated free expansion enters the actual total-to-mechanical
conversion. True physical endpoint stresses are restored before shared native
strain-chord inspection; work-conjugate mean stresses remain separate.

The distinct `Onset` and `Activation` results retain the old motion and the
exact proposed source, endpoint supports, ordered interval-MEAN forces and
solver settings. Activation precedes an elastic solve and has no predictor or
fraction API. Both wrappers keep the static endpoint-force replay owner private.
Their common handoff retains the SAME reference density/mass Arcs, actual
velocity, kinetic storage, clock and consumed steps after original source,
expansion, stress and storage continuity checks. Native thermal brittle
midpoint replay performs the complete interval and consumes one step.

Static and inertial elastic consumers share coefficient/eigenstrain/material
energy recovery. Independent gates admit mean response work, mechanical work,
eigenstrain work, coefficient work, elastic storage, equilibrium and kinetic
balance. Native numerical kinetic loss must be zero. A refused prediction,
handoff or replay publishes no physical successor. The coefficient-first
thermal ordering is not a global second-order or continuous-event theorem;
strain-chord fractions are not physical times. Source changes, thermal paths,
reversals and force means still require physical resolution.

Two native integration regressions specify a free-coordinate solve with
nonzero old strain/velocity and changed stiffness/density/expansion, then
cooling-induced plane discovery and softening replay. A second journey checks
coefficient-triggered fracture in a translating, expanding body. Independent
bilinear-work and full consistent-P1 kinetic formulas test retained density
when current values differ, consumed budgets, separate energy ports and
atomic refusal/replay. Focused formatting, module-path, API/physics review and
scalar calculations completed; native Rust execution remains unavailable
without the required runners and pinned constellation checkouts.

## Thermal tensile onset with original-source fracture handoff

`expansion::initiation::ThermalTensileElasticTetSolid` retains virgin elastic
prehistory on the original P1 mesh and selected thermal sources. The native
equilibrium solver evaluates original elasticity at TOTAL strain minus
integrated free expansion. Tensile detection consumes MECHANICAL strain:
free expansion is stress-free, while restraint or incompatible expansion can
select a crack plane. Reference temperature, source claims, expansion curve
and nonthermal query coordinates remain authoritative.

Each endpoint first observes the held old mechanical strain under the next
coefficients. An already-activated state now returns a distinct
`ThermalTensileTetProbe::Activation`, containing actual source-bound
`ResolvedTensileActivation` observations. It supplies no elastic displacement
predictor, strain fraction, onset temperature or event time. Finite original
supports/forces and native Newton settings are admitted before this early
return. Ordinary strain-chord `Onset` retains its existing diagnostics and
strict native path admission; repeated principal tension still requires an
explicit physical direction hint. Subcritical steps allocate no per-element
activation or event observation arrays.

Both triggers share one private source/state/replay request. Selected planes
and actual V/A widths are prepared at BOTH exact source endpoints. Handoff
seeds only the OLD accepted state and checks original source/reference,
geometry, stress, storage and undamaged mixed-mode history. Native
`ThermalBrittleTetSolid::advance` then owns coefficient-induced release,
joint displacement/eigenstrain evolution and irreversible fracture. No
prediction or preload work is accepted during handoff. Cooling or stronger
coefficients cannot erase damage or accumulated fracture energy.

Elastic acceptance shares the same quadratic material-work recovery with
inertial consumers. Signed coefficient work, total-strain mechanical work,
eigenstrain work, storage and numerical endpoint defects remain separate.
This is the original coefficient-first endpoint model, not an integration of
an unseen temperature trajectory: source extrema/reversals and physical
increments need resolution. The tensile-only criterion, fixed receiving
planes and small-strain limits remain; no heat equation or caloric conversion
is supplied.

Six native regression definitions exercise free expansion, restrained-cooling
onset/softening, separate coefficient work and activation, repeated principal
tension, analytically reconstructed coefficient-induced release, nonhealing
and atomic source/reference/request refusals. Source/API/physics review,
focused formatting and independent scalar work calculations completed.
Native compilation and execution remain unperformed because required
distributed runners and pinned constellation checkouts are unavailable.

## Owned membrane histories for evolving native kinetics

`membrane::visco::owned::OwnedFiniteViscoMembraneMotion` captures an accepted
native membrane once, retaining an owned original patch/law binding and shared
accepted displacement, Maxwell metrics, stress, pressure, resultants and velocity.
Borrowed native views and candidate clones share immutable buffers. All physical
advances delegate to the original native membrane momentum/material solve and
retain its physical heat, pressure work and numerical kinetic-loss accounts.

`with_kinetics` changes only finite positive branch relaxation times. Complete
native law equality pins equilibrium energy, stretch domain, ordered branch
energies, natural reference, solver controls and future law fields. It retains
every accepted physical buffer and ledger exactly, uses the same original patch,
and performs no transverse Newton solve. It supplies no heat, reference work,
impulse or elapsed time. Temperature adapters must shift their immutable original
spectrum; no public snapshot restoration can bind arbitrary histories to laws.

`tests/visco_membrane_owned.rs` covers ownership beyond original source
lifetimes, exact native pressure-interval agreement for both volume treatments
and mass images, regional kinetics with retained memory, nonkinetic refusals,
and late cancellation/work refusal.

## Conservative midpoint tensile onset and native fracture replay

`initiation::dynamics::TensileElasticTetMotion::probe_midpoint` advances the
original elastic prehistory through the native centered inertia and P1 Newton
owners. Its borrowed original elastic backbone supplies exact path-mean stress
and the analytic endpoint derivative C/2. Physical endpoint stresses are
restored before the shared tensile detector and elastic-storage observation;
the accepted step exposes mean stresses separately.

Actual interval-mean dead forces and complete support reactions do work on the
original displacement increment. Elastic work/storage, mean response work,
equilibrium work, retained kinetic storage and the complete energy balance have
independent admission gates. Native numerical kinetic loss must be exactly
zero; numerical damping is never charged as physical fracture or heat. The
same cached consistent reference mass and source density survive each step.

`TensileTetMidpointOnset` privately retains the original motion, physical
interval, ordered mean forces, endpoint supports and solver settings. It
exposes no backward-Euler or static-load replay. Discovered planes and geometric
widths pass through the original source-bound seed; the OLD accepted motion
transfers without impulse, new fracture work or consumed time. Replay uses
`BrittleTetMotion::advance_midpoint` with that exact mean-force request. A
refused solve, handoff or work account leaves the original motion usable.

The native detector still observes an affine strain chord, not a certified
continuous-time event. Exact supplied-path work does not remove oscillator
frequency warping or resolve unseen reversals. Physical timestep refinement,
the original tensile-only criterion and fixed receiving planes remain explicit
model limits. Two native regression definitions cover nonzero-preload and
nonzero-velocity oscillation with independent phase/energy formulas, followed
by inertial onset, conservative handoff and native softening replay with moving
supports and original duplicate force entries.

## Native inertial tensile initiation and conservative brittle replay

`initiation::dynamics::TensileElasticTetMotion` continues an accepted elastic
solid with its original source density, cached consistent reference mass,
actual velocity, clock and step allowance. Initial strains must be strictly
elastic under the retained source-bound tensile criterion. The existing native
backward-Euler owner solves each physical candidate endpoint; shared original
tensile detection then inspects every accepted-to-candidate strain chord.

An entirely elastic chord admits the motion only after independent material
storage/work, endpoint force/support work and native kinetic accounts pass.
Reversible elastic work, endpoint quadrature and backward-Euler numerical
kinetic loss remain separate; none becomes fracture dissipation or heat.
A crossing returns a disposable `TensileTetMotionOnset`, never an accepted
over-strength elastic body.

The onset retains the original physical interval, exact ordered endpoint
supports and applied forces, and solver settings. `resolve_bands` reuses
native discovered planes, original source bundles and geometry-derived widths.
`handoff` seeds the OLD accepted strain through the existing source/stress/
storage continuity gates and retains the SAME mass operator, density inventory,
velocity, kinetic storage, clock and consumed steps. `replay` re-solves the
original interval through native nonlinear brittle motion. Source mismatch,
work refusal or cancellation leaves the original motion and prediction usable.

This is small-strain, fixed-plane, tensile-only, first-order endpoint dynamics.
A detector fraction parameterizes the affine STRAIN CHORD, not physical event
time. Continuous reversals and oscillations require resolved physical intervals;
a smaller interval is solved again using its actual loads and support values.
No continuous-event certificate, new constitutive law or adaptive crack-path
claim is inferred.

Two native regressions specify consistent-mass free acceleration and an
inertial crossing that selects a y-plane from a dormant x-plane source. They
include an already consumed elastic step, a shorter physical re-solve that
differs from interpolation, analytic softening replay, retained momentum/
mass/time/budget, changed-density source refusal and late-refusal retry.


## Prescribed natural-strain events on installed viscoelastic membranes

`membrane::visco::FiniteViscoMembrane::change_reference` and the matching
motion method apply prescribed natural-strain or kinetic-time changes to the
original accepted Maxwell states. Original coordinates/connectivity, material
axes, per-corner reference inventory, source bundle/model, thickness and
full/exact-volume choice survive. Native reference transactions retain every
viscous metric, accumulated loss, displacement and clock. Motion additionally
retains velocity and the complete selected consistent/lumped mass image, with
zero impulse, kinetic change and heat.

Fixed-full-deformation reference work, exact-volume pressure actuation,
transverse loading and endpoint splitting defect have distinct accounts. The
absolute sum of point defects is gated before publication. Transverse plane
stress is restored locally; a subsequent native advance resolves changed global
forces. This is a prescribed zero-duration split event, not a continuous growth
law or automatic thermal coefficient evolution. The native step and motion-step
receipts also retain physical interval loss per original triangle before
cumulative-ledger rounding, permitting conservative thermal handoff without
turning integration defects into heat.

`tests/visco_membrane_reference.rs` covers native agreement after relaxation,
natural-shape release, regional changes, retained graded-inventory momentum and
atomic refusals. These spatial integration tests require the full workspace;
the underlying material reference transactions have focused native numerical
coverage in `fs-material/tests/visco_membrane_reference.rs`.

## Shared-face constraints in brittle crack-contact modes

Adjacent crack faces can repeat the same paired-node displacement constraint.
The native modal projection now reduces exact signed copies of each supported
covector before QR, including reversed pair/direction and aliases exposed by
held coordinates. No numerical rank tolerance identifies an alias; unresolved
linear combinations and near dependencies still refuse in the original QR.

`original_constraints()` retains every admitted physical row, including rows
identically satisfied by supports. `constraint_aliases()` records its exact
signed representative or support-only status. `constraints()` remains the
deterministic independent basis, so existing force-shape indexing stays valid.
Every lifted mode and nullspace vector is checked against all original rows.

Reaction increments are aggregate multipliers determined by the free equations:
each equals the signed sum over its original aliases. They do not determine
individual face pressures or the allocation between contact and support
reactions. `individual_reactions_determined()` exposes that distinction.
Normal, sticking, sourced and thermal entry points share this reduction;
original mass and frozen material stiffness are unchanged.

The incoming scalar-row cap applies before reduction, and dense admission
includes retained witnesses and alias mappings. Focused native-projection
regressions specify exact spring/mass spectra, aggregate force balance,
reversed signs, common sticking components, support-induced aliases and empty
physical spectra. Non-alias dependencies are explicitly rejected.


## Original interface sources in brittle sticking modes

`frozen_sourced_sticking_contact_modes` on mechanical and thermal brittle
motion uses the existing original material-interface friction sampler and
the same exact contact-mode reduction as explicit-law sticking. Bindings
retain ordered incident material identities and an explicit Static coefficient
role. Every loaded original face/vertex samples its own actual tangent speed
in the declared frame, even when subsequent constraint reduction identifies
duplicate covectors. Open and tied vertices do not query coefficients or
construct inactive tangent frames.

The original source owner refreshes temperature/pressure-dependent claims,
accepted distance, texture and slip consistently. Source temperature and
pressure remain explicit external conditions; nodal preload does not invent
pressure. Returned samples are observations: finite-step temperature, slip,
distance, physical time, body state and fracture history do not advance.
All samples return with the fully admitted spectrum, so a late source, mode,
budget or cancellation refusal publishes neither output.

Source coefficients determine whether the held preload fits strictly inside
the native static ellipse; they add no tangent spring or damping. The thermal
route retains accepted mechanical strain and its original reference mass.
Equilibrium and finite-amplitude sticking remain outside this observation.
Focused regressions specify source-temperature loss of sticking, exact
explicit-law recovery, measured speed rather than query placeholders,
wrong-role/material refusal, late-query atomicity and inactive-source omission.


## Static-friction modes of closed brittle cracks

`frozen_sticking_contact_modes` on mechanical and thermal brittle motion
reuses the same exact contact projection with the native directional static
friction ellipse. Explicit face/vertex preloads use the original normal and
two-tangent frame. Native `EllipticForceCoordinates` evaluates the physical
force in its work-conjugate disk coordinates; a caller-specified fractional
interior margin excludes the slip boundary. Loaded, tangentially resting
points with two positive principal coefficients contribute two exact tangent
constraints in addition to the normal constraint.

Two zero coefficients preserve the identical frictionless normal pencil and
require zero tangential preload. One zero coefficient retains the native
sustained law's explicit refusal; no regularization invents resistance on a
free axis. Open and tied pairs require zero force in all components and do
not evaluate inactive friction inputs. Bounds count every scalar constraint,
including both tangents, before allocation and reduction.

The ellipse selects admissible contact membership, not a spring stiffness.
Neither preload magnitude nor ellipse coordinates alter physical constrained
frequencies while the same sticking subspace remains admitted. Original
material stiffness, accepted damage, mechanical thermal strain, reference
mass and time remain with their native owners. Coefficients are explicit
model inputs; no new source authority, equilibrium, finite-amplitude sticking,
frictional heating, sliding trajectory or modal damping is inferred.

G3 regressions specify complete paired-trace equality, spectrum invariance
under a change of tangent frame, native strong/weak ellipse-axis admission,
and exact recovery of the frictionless pencil at zero coefficients. G0 cases
exclude the slip boundary, singular sustained ellipses, moving tangent traces
and insufficient scalar-row bounds. Thermal coverage checks original mass
normalization and unchanged source/fracture histories in both compressive and
fully failed normal-opening material branches.

## Preloaded brittle crack-contact vibration

`BrittleTetMotion::frozen_contact_modes` and the corresponding thermal method
admit an explicit normal preload at every retained face vertex. Strictly open
pairs require zero preload. Touching, tangentially matched, normally resting
independent pairs require a strictly positive preload; unloaded touching,
penetration and unresolved modes refuse. Tied crack-front nodes require zero
preload and introduce no independent constraint. Preloads select a conditional
contact branch; neither supplied values nor interval-mean sustained forces
prove equilibrium.

Native paired normal rows impose `J phi = 0` on the same free coordinates as
the original frozen bulk stiffness and retained consistent reference mass.
Native `fs-la` Householder QR builds the complete nullspace `Z`;
`fs-modal::eigh_gen_dense` solves `Z^T K Z, Z^T M Z`. Lifted shapes retain
zero supported coordinates and original mass normalization. Contact-force
increments are recovered in the original physical normal directions and
independently checked against `K phi - lambda M phi = J^T delta_f`.
No penalty stiffness, constraint inertia, damping, heat or eigenvalue clipping
is introduced. Frictionless tangential motion remains unconstrained.

This complete dense reduction is capped at 128 free displacement coordinates
before contact reduction. Explicit row and dense-entry budgets include lifted
full-node shapes; sparse assembly retains its own bounds. Exact signed aliases
enter one independent supported row; all original witnesses and their mappings
remain available. Support-only rows retain witnesses but no independent free
reaction. Force shapes are aggregate free-equation multipliers when aliases
exist. Other unresolved rank refuses without arbitrary constraint removal.
A fully constrained displacement space has an empty physical spectrum.

Mechanical and thermal owners share the same contact reduction. The thermal
owner retains its accepted mechanical-strain branch, current coefficients and
original reference mass. Damage, temperature, geometry and contact membership
remain fixed. Equilibrium, a finite amplitude that keeps every preload positive,
finite sliding, geometric prestress stiffness and new collision discovery are
outside this conditional observation.

G1 regressions use genuinely fractured and subsequently compressed mirror
tetrahedra with original consistent mass: three mechanisms and eigenvalues
`5E/3`, `20E/3` after normal pairing. They check lifted momentum/reaction
balance, mass normalization, preload-scale independence, thermal density
retention, and the fully failed open-material zero-stiffness branch despite
loaded touching geometric faces. G0 cases cover mode switches, malformed
preloads, moving pairs, cancellation and allocation bounds.

## Inertial motion of finite-strain viscoelastic membranes

`membrane::visco::dynamics::FiniteViscoMembraneMotion` advances the original
membrane configuration and every Maxwell branch with one backward-Euler
physical interval. The same matrix-free Newton equation combines the native
frozen-history material potential with reference inertia. Explicit lumped or
consistent P1 mass retains the original per-triangle corner areal density;
neither deformation nor transverse thickness changes inventory. Consistent
mass keeps moving-support cross terms without a global dense mass or tangent.

Velocity follows the actual accepted displacement increment. The complete
candidate reports physical relaxation heat separately from signed material
endpoint-work error and nonnegative numerical kinetic loss. Their independent
budgets can refuse a converged step without changing accepted motion or history.
Current-area pressure retains the existing prescribed-perimeter requirement.
This is first-order isothermal motion; time refinement remains necessary.
`tests/visco_membrane_dynamics.rs` covers graded-mass translation/acceleration,
held-strain relaxation, moving-support inertia, rejection atomicity and a
three-rung trajectory/kinetic-loss time refinement.

## Retained finite-strain membrane loading, creep and relaxation

`membrane::visco::FiniteViscoMembrane` binds original per-triangle Maxwell laws
to the existing source-bound hyperelastic patch. Equilibrium energy, natural
reference, thickness, grain, mass, volume treatment and local plane-stress
settings remain those of the original sheets. Explicit branch parameters do
not acquire source qualification from this binding. Each triangle retains its
own full or exact-volume remainder history through physical time intervals.

The existing membrane kinematics, reference-area virtual work, geometric
stiffness, matrix-free Newton solve and current-area pressure are reused.
An exact repeated-strain cache within each immutable history/interval avoids
repeating constitutive Hessians during Krylov products; it cannot reuse a
response across a changed history or timestep. Prescribed supports are
eliminated and the existing pressure-perimeter requirement remains in force.

Complete endpoint proposals separate native storage/loading/relaxation loss,
incremental potential, endpoint force work, pressure path work and signed
integration defect. Failed solves, work gates and cancellation publish no
partial displacement, history or elapsed time. `tests/visco_membrane.rs`
checks native material agreement, rotated histories/reactions, pressure creep,
regional spectra, natural references and refusal/cancellation atomicity.

## Bounded torsional shear and layer resultants on native rods

`MaterialRod` and `SpatialMaterialRod` extend their original
`segment_torsional_stress` with `segment_torsional_stress_bounded` and
`segment_torsional_layer_response` on their accepted native Cosserat strains.
They reuse the canonical circular-section point law and cached per-layer G,
including actual material-temperature refresh. Coordinates use the segment's
geodesic-midpoint material frame. Spatial rods subtract their original
intrinsic curvature, so stress-free reference twist creates no shear stress.
Layer torque and energy sum to the original native torsional resultant and
storage on each actual reference interval. There is no re-solve, source refit,
geometry update or independent copied constitutive history.

The material field regressions compare recovered torque/stress with a loaded
shaft's analytical solution, native energy on non-unit lengths, heterogeneous
series compliance, explicit surface limits and stress-free intrinsic twist.
These observations are only the Saint-Venant
torsional contribution; transverse shear, axial/bending stress, end/interface
effects, warping, plasticity and a failure criterion are outside their scope.
The caller supplies an explicit small-engineering-shear applicability bound.

## Frozen-history brittle vibration operators

`BrittleTetMotion::frozen_vibration_pencil` and
`ThermalBrittleTetMotion::frozen_vibration_pencil` expose sparse `K, M`
operators from the original accepted histories. Native cohesive integrity
supplies conservative unloading storage: complete tensile/shear failure
retains physical zero stiffness, compression retains its original penalty,
and integrity is never reconstructed as `1 - damage`. The existing P1 tangent
assembler supplies work-conjugate stiffness, while the exact retained
consistent reference mass supplies inertia. Temperature-dependent current
density does not replace reference mass. The thermal owner uses accepted
mechanical strain to choose opening/compression, holding thermal eigenstrain
and current coefficients fixed.

Explicit incremental supports select ascending free coordinates of both
operators. Count bounds, owner/shape checks and cancellation propagate before
returning a complete immutable observation. Extraction does not advance
damage, temperature or time, consume motion steps, or add artificial residual
stiffness, modal loss or fracture heat. Strictly open released-face pairs may
be observed; touching, closed or tolerance-ambiguous independent traces need
a contact tangent and refuse. Tied crack-front nodes remain tied. A damaged
material exactly on its opening/compression switch also refuses.

This is a conditional frozen-history linearization. The caller establishes
equilibrium when equilibrium frequencies are required. A state on a damage
envelope can evolve under arbitrarily small outward loading; this observation
does not establish a finite vibration amplitude, evolving-fracture stability,
geometric prestress stiffness, thermoelastic damping, contact-free trajectories,
new collision pairs or acoustic radiation. Existing `fs-modal` pencil
consumers remain the eigensolver boundary.

The focused `brittle_tet_modes` and thermal-motion `modes` regressions specify
analytic frequency loss, exact failed null modes, retained compression,
tensor-work symmetry, consistent mass/support projection, current thermal
stiffness with original inertia, immutable history and refusal behavior.
Native mechanical and thermal fragment-motion fixtures also specify separate
open-fragment operators and rejection at touching arrival.

## Native rod torsional stress observations

`MaterialRod::segment_torsional_stress` and the same method on
`SpatialMaterialRod` recover a layer's local torsional shear field from the
accepted native geodesic-midpoint elastic twist and original source section.
The spatial owner subtracts intrinsic reference twist before recovery, so an
unloaded pretwisted rod has zero torsional stress. The point is specified in
the segment's local [y,z] coordinates; reported tensors share that material
frame. These allocation-free observations do not alter equilibrium or history.
Existing uniform/world-couple, heterogeneous-segment and intrinsic-shape tests
exercise the adapter against analytic torsion. The underlying model remains
small-strain circular torsion and does not resolve end effects or add the
separate axial, bending and transverse-force stress contributions.

## Source-held and thermal frictional trajectories

The same complete-duration friction driver now accepts the original material
interface sources and the native thermal brittle body. Impacts query the
caller's explicitly selected static or kinetic role at each active row's actual
incoming slip speed. Sustained sticking queries the same original source's
Static role at the actual resting speed. The last accepted source candidate is
retained by complete paired-row identity; numerical localization probes never
publish their source samples.

Returned source records identify the accepted physical event or interval and
carry the original selected claims, coefficients and source state. Fixed external
conditions and texture frames stay explicit. Physical travel, slip and the last
finite-step contact temperature remain unchanged for these admitted modes.
Neither force-derived pressure nor contact heat is inferred. A refusal after an
accepted internal prefix discards the entire body/source transaction and leaves
the caller's original histories available for retry.

Thermal trajectories hold the current admitted bulk source snapshots over each
call and use the original thermal midpoint for every candidate. Reference matter,
total/mechanical strains, expansion reference and complete thermal work receipts
remain with the actual body. Explicit thermal source transitions belong between
calls; no intermediate material data is interpolated for localization. Both
fixed-coefficient and sourced thermal entries reuse the same friction policy,
constrained locator and complete physical work admission as the mechanical body.

The source regression uses an inelastic oblique capture under the explicit
kinetic role followed by static shear exceeding kinetic capacity, with exact
native mass/energy expectations and late whole-duration rollback. The thermal
regression continues a genuinely warmed fractured body under static shear,
retaining its original reference density and signed moving-support work across
multiple accepted intervals. These are manufactured implementation checks,
not experimental material validation or a general finite-sliding claim.

## Frictional time evolution after brittle fracture

The native mechanical released-face trajectory can now advance a complete
requested duration through simultaneous directional-friction impacts, sustained
static load transfer, unloading and free separation. The existing normal and
thermal impact entries retain the same time loop through a typed contact-policy
seam; the original midpoint, material, consistent-mass, solver and contact owners
remain responsible for every physical response.

Constrained localization keeps only the initially matched, static friction rows
in the force solve while scanning every retained released pair for arrival.
Each trial uses its actual prefix of the original affine support/load path.
Discarded trials retain no material successor and spend only the bounded native
response allowance. Only the accepted prefix consumes a physical step; impulses
consume an event allowance without advancing the clock. A failed full duration
returns no partially advanced body.

Directional force warm starts are physical forces matched by the complete
paired-row identity. Each current ellipse maps those forces into its own
work-conjugate numerical coordinates; transformed solver forces are never
reused as physical tractions. The complete work receipt separately retains
normal and tangential constraint work, the sum of each component's absolute
work, fracture loss, impact loss, support work and native kinetic increments.
Static constraint work is a bounded numerical defect, not heating or damping.

This path retains the existing fixed paired-face, small-strain and source-domain
limits. Loaded finite sliding, new-pair discovery and thermal allocation of
impact loss remain outside it. The native sustained ellipse still refuses a
single zero principal coefficient. Numerical safe/unsafe endpoint brackets do
not establish global continuous collision detection. Caller-selected interval
resolution and explicit load/source knots remain required.

Focused G0/G1/G3 journeys exercise native oblique capture followed by sticking,
affine load reversal followed by release, and arrival while an earlier face
remains loaded. Their force, impulse, internal-mode and work expectations come
from the original P1 consistent-mass formulas; they also cover late complete-time
budget refusal and deterministic retry.

## Sourced static friction and thermal fracture load transfer

The sustained-friction source child binds each retained face vertex to an
explicit original ordered interface and accepted history. It requires the
Static coefficient role and actual incident material identities, sharing the
impact adapter's identity and coefficient-sampling helpers. Only initially
active, matched, static rows sample their own actual initial relative speed.
The original source owner retains selected claims, source authority, complete
query/domain checks, accepted distance, slip and texture.

This is a source-held interval law: sampled coefficients and caller-supplied
external contact conditions are fixed for the accepted sticking interval.
Returned samples are the prepared coefficient/source candidates, not a new
finite-step thermal/friction-history integration. Accepted travel/slip and the
last finite-step contact temperature remain unchanged. Parameter changes,
reversals, finite sliding and thermal destinations require explicitly resolved
additional physical paths. Source preparation and every native force response
remain private until the whole native body interval passes; refusals publish
neither source nor body state.

The thermal `contact::friction::sustained` adapter provides current-source and
explicit `at_materials` transitions, each with explicit or sourced friction.
Both use the same generic contact result and original thermal midpoint.
Every trial replays the same accepted mechanical history and the same resolved
thermal successor. Native total/mechanical strain, thermal expansion reference,
reference density/mass, storage, fracture and kinetic accounts retain their
original owners. Signed parameter and eigenstrain work remain independent of
external/support work and bounded static contact work. No heat or additional
fracture debit is inferred from the constraint.

Focused native journeys exercise static capacity above the source's kinetic
coefficient, original source/path/history retention, unsupported-role and
insufficient-capacity refusal, deterministic retry, and a bulk thermal
transition under normal/shear load. The thermal reference independently checks
force, signed eigenstrain work, support/dead-load work, frozen original mass,
native history and transactional source transfer.

## Sustained static friction between released brittle fragments

The `contact::friction::sustained` owner advances the original nonlinear
brittle midpoint with simultaneous normal and tangential mean contact forces.
It reuses the existing force/configuration solver and work-conjugate elliptic
coordinates. All same-body and normal/shear coupling comes from complete
native responses; no rigid proxy, impulse-over-time conversion, penalty or
additional compliance supplies the force. Shared force assembly, response
replay and external/contact energy admission also serve the existing normal
contact path. A public generic result retains each original motion/energy type.

Only initially active, matched pairs at rest in all three relative velocity
components may carry load. Loaded pairs must remain matched, active and static
at both native endpoints with bounded physical tangential displacement.
All other non-tied pairs remain in the solve with zero friction so unexpected
loading or arrival is detected and refused. Unloaded pairs may separate freely.
The native sustained ellipse admits two positive coefficients or equal
coefficients including zero; its unsupported single-zero directional case is
refused without artificial resistance.

Every trial replays one accepted material/mass/history state. Final original
force replay, strict nonnegative gaps, both-endpoint complementarity, physical
ellipse capacity, actual path/force residuals and endpoint speeds gate
publication. Separate normal/tangential and absolute contact work, physical/
solver coordinate work, native/external work and absolute individual pair
couples must pass. Static constraint work is a bounded numerical residual,
not dissipation or heat. One accepted interval consumes one physical step;
a refusal publishes no partial material or motion state.

Focused native tests use exact P1 consistent mass and independent work/momentum
inventories. They cover timestep-independent force with deforming apex motion,
moving-support shear transfer, actual unloaded release, rotated partial support
with normal/tangent coupling, insufficient static capacity, bounded refusal and
replay. The model remains small-strain and fixed-pair: sustained sliding, new
spatial contact discovery and unlocalized impact/arrival require further modes.

## Source-bound friction on mechanical and thermal brittle fragments

The `contact::friction::source` adapter binds an explicit ordered interface,
accepted source history, complete external conditions and static/kinetic role
to each retained face vertex. Both principal interface sources must name the
actual material states of the incident elements. Source face IDs must match
the retained face order; texture coordinates are explicitly declared in the
same contact frame. The original source sampler resolves each active vertex
at its own actual incoming tangential speed and retains selected claims,
authority and validity-domain checks. Declared accumulated-distance coordinates,
native slip, transported texture and the last finite-step temperature history
stay at their accepted values. Inactive and tied vertices do not query sources.

Source samples remain private until the simultaneous contact, native momentum,
actual endpoint and work audits all pass. The returned event contains both
the complete mechanical successor and the active source candidates indexed by
original face/vertex identity. Input histories are immutable, including after a
late source or mechanical refusal. Pressure, time, travel, temperature change
and impact heat are never inferred from an impulse. The coefficient role
selects one event ellipse; it does not infer an intrastep static/kinetic switch.

The thermal `expansion::dynamics::fracture::contact::friction` child provides
both explicit and source-bound entry points through the same physical kernel.
Thermal material identity is checked through the actual element's retained band
source. Accepted bulk source bindings, thermal expansion reference, mechanical
strain and fracture history, displacement, stored energy, fracture loss,
reference density/mass, physical time and step budgets survive exactly; only
velocity and kinetic storage change. Contact-source coefficient refresh is not
a bulk-temperature update or thermal allocation.

Focused native tests cover actual and distinct per-vertex slip queries,
accepted history/path retention, wrong material identity, inactive queries,
late refusal and deterministic retry. The thermal journey uses cooled source
density with retained original mass, directional free slip/sticking, independent
P1 momentum/work inventories and subsequent native motion. All fixed-pair
geometry and contact-law limits of the frictional-impact owner remain in force.

## Native frictional impact of released brittle fragments

The `brittle::solid::dynamics::fracture::contact::friction` owner accepts
explicit per-vertex directional friction ellipses and one declared tangent
direction per retained face. Every active non-tied pair enters one simultaneous
normal/tangential solve, including initially separating pairs. Sparse physical
paired rows assemble the complete `J M^-1 J^T` from original consistent-mass
responses, with all normal/shear and intercontact blocks retained. One combined
impulse then updates the original body. Contact/scalar-row caps bound quadratic
scratch; only one full-body response is retained at a time.

The original elliptic Coulomb owner retains its complementarity, maximum-
dissipation, passivity, rank and determinacy limits. Equal principal coefficients
recover circular Coulomb friction; a zero principal coefficient is exact free
slip. The shared original physical-ellipse audit independently checks the actual
body endpoint, not merely a predicted contact velocity. Native impulse and work,
moving-support reaction work, total relative work, and the sum of absolute
individual contact couples must agree before publication. A velocity-only
determination does not establish unique local impulses, pressure or heat.

Accepted events retain physical configuration, source/material/fracture history,
reference mass, physical clock and step allowance. Any refusal leaves the
accepted motion available. Penetration and active tangentially separated paired
witnesses refuse. This is instantaneous contact of retained matching faces; it
does not supply finite-sliding geometry, new-pair discovery, sustained frictional
forces, or thermal allocation of the combined impact loss. Explicit coefficient
inputs carry no inferred material-source authority.

Focused manufactured tests use exact P1 simplex mass and independent momentum,
angular-momentum and kinetic inventories. They cover frictionless parity,
sliding, sticking, an exact free axis, physical texture under tangent-frame
rotation, moving-support work and a rotated partially constrained specimen whose
normal impulse depends on friction through a nonzero normal/tangent mobility.

## Brittle arrivals on an already loaded contact trajectory

`advance_toward_released_face_contact_with_sustained` locates a new released-face
arrival while initially resting, active, matched pairs retain the original
nonlinear sustained-force response. Initially open or separating pairs receive
no contact force during the search. Every complete candidate is scanned against
all released openings; only an overlap on an unforced pair narrows the sampled
endpoint bracket. Participating-row, material, solver, mode and work refusals
remain terminal. No free-motion substitute, force across an open gap, state
interpolation or position/velocity correction supplies an accepted prefix.

The existing affine-load sampler gives each trial its actual prefix force mean
and prescribed endpoint. All trials replay one accepted material/mass/history
state. Warm starts map forces by complete paired-row identity, never subset
position. Endpoint attempts and every native mechanical response have separate
whole-search budgets, including discarded candidates and final force replays.
Only the accepted prefix spends a physical timestep. An initially active
approach returns a zero-time witness without a force solve.

Accepted paths retain strict nonnegative paired openings, the original released
face guard, native support work, and all normal complementarity, signed/absolute
work and force-couple gates. The reported bracket concerns sampled native
endpoints and the accepted affine trace, not globally earliest continuous
collision detection. Fixed paired faces, small strain and frictionless normal
response remain the declared model.

G1/G3 tests in `tests/brittle_tet_fragment_contact/sustained/approach.rs` use a
shared-mass, two-face fractured specimen. Independent consistent-mass equations
check constrained acceleration, aggregate reactions at duplicated rows, contact
time, external work and retained fracture history, with bounded atomic refusal
and retry.

## Brittle time continuation through impact, sustained contact and release

`advance_released_face_time_with_sustained` extends the shared native time
driver with constrained approach, simultaneous impact, loaded dwell and unloaded
release. Active approaching pairs receive the original impact response first.
Initially resting matched pairs retain their native nonlinear sustained forces
while every released pair is checked for arrival. A new arrival shortens the
constrained interval; the subsequent simultaneous impact includes all active
pairs, and the remainder recomputes the physical contact subset. An old contact
may release when the shared consistent-mass response requires it.

Every trial samples forces and prescribed supports from the original global
affine path. Previous accepted force receipts seed later solves by full paired-row
identity, and seeds are cleared after impact. They are numerical iterates, never
material history or prescribed forces. Every mechanical response, including
discarded constrained endpoints, Newton probes and final replay, spends the
complete native-solve budget. Only accepted intervals consume physical steps.
All original source, material, mass, geometry, mode and work gates remain active.

Each force interval retains its original total-load receipt beside external-only
physical work. The shared driver preserves native energy types and separate
parameter/eigenstrain accounts for its thermal caller. Compensated kinetic,
elastic, fracture, applied/support and impact-loss accounts retain their physical
roles. Signed and absolute sustained normal work have their own whole-duration
gate; this numerical account is never added as external work, dissipation or
heat. A late refusal returns no partial successor and leaves the original motion
available for retry.

Force/approach activation and normal-speed tolerances must agree. The chosen
maximum interval remains temporal resolution, not an adaptive error certificate.
Fixed paired faces, small strain and frictionless normal response bound this
entry. Sampled endpoint brackets do not establish globally earliest continuous
collision detection; finite sliding and new-pair discovery remain outside it.

G1/G3 cases in `tests/brittle_tet_fragment_contact/sustained/time.rs` retain
loaded intervals, accelerated impact/dwell and affine load reversal. Additional
`sustained/mixed_time.rs` cases independently reconstruct a shared-mass arrival,
the unilateral impact that releases a different corner, the changed sustained
force subset, affine global-time loads, physical work and final motion. Late
interval/response refusal and deterministic full retry retain the predecessor.
Existing mechanical and thermal contact tests exercise the same shared kernels.

## Sustained contact between native brittle fragments

`BrittleTetMotion::advance_released_face_sustained` solves simultaneous
frictionless normal forces through the original nonlinear brittle midpoint,
using `fs_contact::sustained::solve`. Every force response starts from the same
accepted material history, frozen consistent mass and physical clock. The
shared contact owner retains its Newton–Krylov/FGMRES, native force-response
differentiation, final replay and force-space feasibility admission. All
non-tied released pairs participate within its explicit 21-contact hard cap.
The numerical normal-map scale supplies conditioning only, with no physical
spring, compliance, stored energy or impulse-to-force conversion.

A private midpoint candidate seam defers only released-face opening admission
during unpublished signed-force probes. Every original material/source,
constitutive history, inertia, Newton and work gate still runs. Public midpoint
retains its original opening guard. Only the final freshly replayed contact
candidate, independently checked for strict nonnegative actual gaps and normal
forces, can publish one physical motion step. The zero-friction embedding uses
only normal loads and requires exactly zero final algebraic tangent forces.

Loaded pairs must remain active, normally at rest and tangentially matched at
both ends of the native affine displacement trace. Approaching impacts,
open-to-loaded arrivals and unresolved modes refuse instead of receiving
fictitious contact damping. Unloaded pairs may separate. Independent gates
check normal force-gap complementarity at both endpoints, signed and absolute
per-pair normal work, nodal/paired work agreement, the native applied-load work
split, internal resultant force and the sum of absolute pair force couples.
Normal contact work remains bounded numerical work, not fracture loss or heat.

The result retains both the original total-load midpoint receipt and an
external-load-only receipt including native support reaction work. Original
compressive material response, source ownership, cumulative fracture loss and
reference matter remain present. An optional private normal-force seed supports
subsequent numerical solves without becoming physical history. Trial response,
contact-space and native solid budgets stay explicit; refusal/cancellation
leaves the predecessor unchanged. This is fixed-pair, small-strain normal dwell
and release, not finite sliding, friction, new-pair discovery, global continuous
collision detection or a continuous-trajectory accuracy certificate.

G1/G3 regressions in `tests/brittle_tet_fragment_contact/sustained.rs` use
independent consistent-P1 mass and compressive stiffness calculations for closed
loaded motion, retained elastic storage after fracture, outward release, moving
support work and bounded atomic refusal/retry. Their compression branch prevents
a mass-only or rigid-fragment replacement from satisfying the mechanics checks.

## Complete brittle-fragment time intervals through impact

`BrittleTetMotion::advance_released_face_time` composes the original native
midpoint, bounded released-face approach and simultaneous consistent-mass
impact owners over a complete requested physical duration. Each accepted
motion prefix retains its actual material history, physical endpoint state and
cached reference mass. Every instantaneous event gathers all active paired
rows, retains the native coupled contact determinacy, and consumes no physical
timestep. Only a complete duration returns its final motion; a late refusal
cannot publish internal partial progress.

Affine nodal forces are resampled from the original global interval at every
actual prefix and remainder. Affine support paths retain the caller's absolute
endpoint; their initial velocities must match the prescribed slope within the
explicit velocity allowance. The same support DOFs retain their current
velocities through impact, with the native signed reaction work included.
Motion, impact, native-trial and contact-row budgets remain explicit. Repeated
events separated by no more than the localization-time allowance refuse as
unresolved, rather than manufacturing a persistent-contact impulse sequence.

Compact records retain each accepted midpoint work receipt and each original
contact/body-impulse receipt. Compensated whole-duration sums retain material
work, fracture loss, applied/support mean work, signed impulse-support work,
native polarized kinetic increments and separate impact loss. Complete signed
work balance and the sum of absolute local work defects must both fit their
whole-duration gate; mean-stress and cumulative-fracture accounts are checked
again. Contact loss is mechanical dissipation, not inferred heat. Native
midpoint numerical kinetic loss remains zero.

The maximum free-flight interval is an explicit caller-selected temporal
resolution, not an adaptive error estimate. The approach owner's guarantee is
the nonpenetrating native affine-displacement trace and sampled endpoint
bracket, not global continuous collision detection or a convergence-order
claim through impacts. Fixed matching material faces and small-strain source
admission remain in force. Persistent contact forces, finite sliding,
friction, new fragment pairs and automatic heat allocation are not supplied.

G1/G3 regressions in `tests/brittle_tet_fragment_contact/time.rs` exercise a
complete ballistic collision and remainder against independent consistent-mass
impulse, momentum and energy oracles; original-global-time affine-force means
across multiple native intervals; moving-support reaction work; and late
budget/persistent-contact/cancellation refusal with immutable deterministic
retry. The original four simultaneous-impact tests and three approach tests
are retained as the same focused native integration lane.

## Native brittle-fragment contact localization

`brittle::solid::dynamics::fracture::contact::approach` advances an original
`BrittleTetMotion` toward approaching released faces without requiring a caller
to choose the impact time exactly. Every trial replays the same accepted native
material history through `advance_midpoint`, its cached consistent mass, Newton
solve and physical-work admission. Affine applied forces are resampled at each
actual prefix mean; prescribed coordinates follow the accepted-to-target path.
No state interpolation, gap correction, extra material law or impulse is added.

Only an explicit native released-face closure or a strictly negative paired gap
narrows the unsafe endpoint. A successful full endpoint returns directly;
otherwise bounded bisection retains a nonnegative native prefix until an actual
approaching witness lies within both the positive activation-distance allowance
and the safe/unsafe endpoint-time allowance. Every non-tied pair is checked
before returning a witness, so another pair's penetration cannot be hidden.
Even penetration inside the existing opening roundoff tolerance is rejected.
Initially active approach returns a zero-time witness without spending a motion
step; all accepted nonzero prefixes spend exactly one. Native attempts have a
hard count, cancellation is polled, and any refusal leaves the original clock,
velocity, history and step allowance available for deterministic retry.

The geometric claim is nonpenetration on the accepted affine-displacement
trace. The reported time bracket concerns sampled native endpoint admission;
it is not a globally first collision certificate for a continuous or
nonmonotone trajectory. Callers choose temporal resolution and resolve forcing
knots. Persistent contact forces, finite sliding and new fragment-pair discovery
remain outside this entry. A witness is only an approach observation; the
original simultaneous impact owner still collects every active contact row.

G1/G3 regressions in `tests/brittle_tet_fragment_contact/approach.rs` compare
ballistic closure and constant-acceleration contact against independent roots,
check actual prefix means under affine transverse forces, and cover finite
precontact motion, initial approach/separation, strict tiny-overlap rejection,
bounded refusal, cancellation and retained deterministic retry. These are
manufactured native mechanics checks, not experimental validation or global CCD.

## Simultaneous impacts on released brittle faces

`BrittleTetMotion::impact_released_faces` binds the retained matching-face
kinematics to the original simultaneous normal-impact owner in fs-contact.
Every active non-tied vertex pair participates in one complete J M^-1 J^T
solve using the original supported consistent mass. At most 64 active rows are
admitted; one native sparse response is retained at a time. Redundant rows keep
the original PSD velocity/work witness without inventing unique local pressure.

The resulting signed nodal impulses are applied together through the original
mass impulse owner. Actual endpoint restitution/complementarity, linear/angular
impulse, individual contact couples, retained kinetic storage, signed support
work and total contact loss must pass before publication. All position, source,
damage, elastic storage, fracture-loss history, density, time and step budgets
remain unchanged. Impact loss is explicit and is not inferred thermal heat.

Contact requires nonnegative gaps and matching tangential witness positions;
penetration and finite-sliding geometry refuse rather than being repaired.
Restitution is an explicit model input. The caller locates the instantaneous
event; this supplies no global collision search, persistent contact force or
friction. Native motion retains its released-face opening guard afterward.

## Completed brittle fracture and native fragment continuation

`BrittleTetSolid::split_fracture` releases explicitly selected conforming faces
only after BOTH incident source-bound bands have completely failed, align with
the physical face and carry no compressive contact or unresolved face traction.
Original element/local-corner order, coordinates, V/A widths, strain, stress,
source ownership, stored energy and cumulative fracture loss are preserved.
Every previously dissipated band inventory remains charged once; a selected face
does not replace it with a fresh Gc*A. This is connectivity release of the native
damaged bulk, not conversion to intact elastic fragments with new cohesive laws.

`BrittleTetMotion::split_fracture` additionally duplicates actual velocity,
reconstructs the original consistent mass on split nodal stars and checks native
kinetic storage against independent P1 mass/momentum/energy inventories. Original
density, physical time and used/remaining step budgets survive. A cut is neither
an impulse nor a physical timestep. Subsequent cuts remap retained face guards
by original element-local corners, preserving earlier nonzero openings.

All native static, backward-Euler, midpoint and temperature continuation keeps
the released matching-face guard: affine fixed-normal opening is checked at all
three paired vertices before publication. Closing endpoints require contact and
refuse atomically. This is not arbitrary collision detection, finite sliding,
loaded partial-damage surface insertion or a mesh-objective crack-path theorem.
Six native structural regressions cover actual continued solves and motion,
conserved material/inertia, thermal history, sequential partial fronts, closure
and cancellation/refusal retry, supplementing five material-release checks.

## Temperature-dependent brittle midpoint continuation

`BrittleTetMotion::advance_midpoint_at_temperature` reuses the original private
source/history temperature transfer, then completes the native mean-stress
midpoint mechanical interval before publishing any successor. It retains the
original consistent mass and density snapshot through that parameter stage,
along with velocity, kinetic storage, clock and remaining physical-step budget.
A changed queried density cannot silently create matter or momentum.

The combined first-law gate accounts for signed parameter work, external
mean-force work, actual storage change, thermal plus mechanical fracture loss,
and the native stable kinetic increment. It adds no backward-Euler loss or
endpoint-work correction. Source, domain, nonlinear, cancellation and late work
refusals retain the pre-transfer motion for retry. Temperature is a coefficient
change at fixed accepted strain followed by motion at those new coefficients;
this ordered split does not infer thermal expansion, caloric heat or global
second-order accuracy through an unresolved temperature history.

## Work-conjugate brittle midpoint motion

`BrittleTetMotion::advance_midpoint` uses the original analytic crack-band
path-mean stress, retained consistent reference mass and native midpoint
inertia through the existing P1 Newton/FGMRES assembly. Applied loads and support
reactions are interval means; displacement, velocity, constitutive history and
published solid stresses remain physical endpoints. Independent replay and
whole-mesh gates check material work, cumulative fracture loss, mean-force work,
kinetic continuity, equilibrium work and the complete first law. No endpoint
quadrature correction or numerical kinetic loss is relabeled as physical loss.

Unavailable classical mean tangents require an explicit solver option. The
native finite-difference owner then evaluates only the affected element's
six-component actual mean response at fixed history, with at most seven material
trials using a cached base. It preserves domain/cancellation refusals and uses
the general sparse solver; all ordinary endpoint paths remain strict. This is
centered, fixed-plane small-strain motion. Exact supplied-path work does not
resolve unseen reversals, certify temporal error, remove frequency warping or
initialize newly inserted crack surfaces in loaded bodies.

## Continued brittle seam growth with retained physical histories

`TetFractureTopology::refine` adds explicitly selected original faces to an
already split mesh. Every new nodal star maps to exactly one PREVIOUS split
star, preserving independent crack-face displacement and velocity instead of
looking them up on the unsplit mesh. Original elements, quadrature layout,
material bindings and the retained seam prefix stay fixed. The existing native
splitter and consistent P1 inventory owner verify unchanged mass, linear/angular
momentum and kinetic energy; no lumping, averaging, deletion or healing occurs.
Budgets cover the complete result, and refused/cancelled candidates leave the
accepted topology and fields reusable.

The retained-interface admission checks both complete matching P1 maps against
the actual reference geometry and verifies unchanged physical points, axes,
areas, shape order, surface roles and borrowed law/source owners on old seams.
It uses the native P1 frame arithmetic with linear seam work and no mesh-sized
temporary allocation. This permits retained point histories to remain with
their physical locations; it neither rebinds opaque solver checkpoints nor
initializes newly inserted cohesive laws in a loaded solid. New-surface
traction/history/energy admission and boundary/load reassignment remain with
the receiving solver; automatic crack discovery is not claimed.

## Partial brittle crack fronts and retained reference trace order

Matching P1 cohesive traces admit one or two shared crack-front vertices. Their
identical paired shape columns cancel before sparse assembly, so tied-front
motion contributes no jump, cohesive work or spurious force. A completely tied
face has no opening mode and refuses; nonmatching trace admission is unchanged.
`TetFractureTopology` therefore admits partial explicit seams without cutting
unselected faces. Ordered matching constructors preserve the supplied physical
A-side vertex order, quadrature identities and local axes through global node
renumbering. Callers retaining interface histories must preserve that reference
order and their original material/source bindings; constructing geometry does
not initialize or transfer loaded cohesive histories.

## Nonuniform source-bound Cosserat reference intervals

`SpatialMaterialRod::try_nonuniform_straight` and
`RodReferenceGeometry::try_nonuniform` admit an explicit positive material
interval for each original section. `RodMaterialSource::try_new_nonuniform`
and `try_new_nonuniform_with_isotropic_parameters` retain that partition
through source resolution and thermal refresh. Cumulative reference stations
must represent every interval within the existing dimensionless geometry
tolerance (64 machine epsilons); an unrepresentable tiny interval is refused
instead of introducing artificial prestress or changing its authored inventory.
Existing uniform constructors retain their original division and station
arithmetic.

Each interval drives the original shear/extension and curvature measures,
analytic elastic force and tangent, immutable reference mass and rotary
inertia, gravity and hydrostatic volume, and prepared dynamic displacement
gradients. `RodReferenceMass::segment_lengths_m()` exposes the complete
partition. The coupled axial, transverse, torsional and bending loss adapters
use the corresponding section resistance divided by its own interval, and
finite-cell enthalpy admission checks the matching original mass and volume.
A source refresh cannot silently replace a partition with another having
the same total length.

The focused rod regressions cover independent axial/torsional series
compliance, graded bending convergence, reference mass and inertia,
objectivity, uniform parity, station representability, thermal continuation,
all four dissipative adapters and caloric inventory. The existing bounded
slender-rod model, frozen section geometry and dense static solve remain the
declared scope; this does not add adaptive remeshing or a history-transfer law.


## Exact matrix-free finite-strain tangent actions

`hyper3::prepared::PreparedHyperTetProblem::tangent_apply` streams the original
per-element constitutive Hessian into `delta F`, `delta P` and nodal force
increments without assembling element or global stiffness matrices. The caller
reuses full-DOF output and scratch buffers; cancellation or a domain/arithmetic
refusal leaves output unchanged. Relative increments retain the exact translation
null mode and do not subtract identity from a rounded deformation gradient.

`linearized_tangent` freezes the original 9-by-9 Hessian per reference tetrahedron
under an explicit payload byte allowance. Repeated fallible forward and transpose
actions do not reevaluate the material. Original heterogeneous, Gent and thermal
reference laws remain the constitutive owners. The actual numerical transpose is
applied without silently symmetrizing the material entries. The operator contains
all original DOFs; support elimination, follower loads, gravity and inertia remain
separate. No new mesh accuracy, locking, throughput or empirical claim is made.

`tests/hyper3_tangent_action.rs` checks assembled-operator parity, independent force
differences, exact rigid translation, frame covariance, very small increments,
thermal-reference response, cached evaluation counts, transpose pairing, Gent
chain stiffening and transactional refusal. Existing prepared-force tests remain
the original reference-geometry regression seam.


## Explicit brittle through-seams and conservative fragment transfer

`brittle::topology::TetFractureTopology` lowers explicitly selected internal
tetrahedral faces into separate exterior traces for the existing P1 cohesive
owner. Per-node tetrahedral-star connectivity duplicates only the nodes needed
by the requested cuts. Unselected faces keep their shared DOFs; original
element/local-corner order, reference coordinates, seam IDs and explicit A/B
roles remain unchanged. Fragment IDs follow surviving volume-face connectivity.
The existing matching-trace quadrature supplies actual triangle areas.

Construction requires original `TetLoadGeometry3` admission and explicit
node/face/seam budgets. Exterior, nonmanifold, repeated or wrong-owner faces
refuse; original disconnected node stars and still-tied crack-front vertices
refuse because the current cohesive trace owner requires distinct side DOFs.
Construction and transfer are cancellable and publish no partial result.

Motion transfer duplicates nodal displacement and velocity exactly. Independent
before/after reconstructions preserve consistent P1 mass, first mass moment,
linear momentum, angular momentum using current positions, and kinetic energy.
Fragment inventories and split nodal row-sum masses remain inspectable; those
row sums do not replace the consistent mass operator. Element material/source
bindings and history blocks retain their original borrowed slices and order,
without averaging, resetting damage or rebinding opaque solver checkpoints.

The six pure-core tests use independent degree-two tetrahedral quadrature,
consistent-versus-lumped kinetic energy, rotation/origin transformations, history
identity, surviving-face continuity and refusal/cancellation checks.
`tests/brittle_topology.rs` exercises admitted geometry, native P1 mass/stiffness,
force-driven cohesive failure, separated motion and compressive reclosure.

This is an explicit through-seam preparation/transfer path. It does not convert
a loaded smeared crack-band state into a fresh cohesive surface. Such insertion
still needs traction/history/fracture-work transfer, and evolving crack fronts
need a trace owner supporting tied front nodes. Existing cohesive mechanics owns
subsequent fracture dissipation, inertia and reclosure; no element mass is deleted.

## Directional brittle crack bands in native solid solves

`brittle::TetCrackBandGeometry::from_element` consumes the original admitted
`loads3::TetLoadGeometry3` and one explicitly supplied reference normal.
It retains the original element index/connectivity and volume, normalizes the
direction, and computes the convex tetrahedron's orthogonal footprint from
all four faces: `A_projected = sum |n dot area_vector| / 2`. The width
`h = V / A_projected` makes complete pure-mode fracture work
`V*(Gc/h) = Gc*A_projected`. This footprint is not the area of an arbitrary
interior cut plane. Finite nonzero directions at either floating-point extreme
are scaled before normalization; zero/nonfinite directions, invalid element
indices and unrepresentable geometry refuse. No allocations or matrices are
created by this bounded geometry query.

`material_settings` carries that exact width and normal into the original
`fs-material::brittle::BrittleCrackBand`. Its borrowed `SmallStrainLaw`
uses the existing tetrahedral and HEX8 constitutive assembly, free-coordinate
equilibrium, nonsymmetric tangent solve and accepted-history publication.
Per-element fields keep distinct widths/laws; no mean material or new solver
is introduced. The fully integrated HEX8 path retains its original eight
quadrature histories and requires its own explicit geometric band width.

The six regressions in `tests/brittle_crack_band.rs` specify geometric
rotation/translation/scaling and winding admission, a free-sided tetrahedral
tension/softening/unloading/failure/reclosure path, rejected-step isolation,
distinct per-element widths, and native HEX8 fracture-work integration.
The uniaxial expectations independently condense the transverse strains;
loaded-node reaction is `sigma*V/L`, not a guessed face force. These tests
are manufactured-model checks, not experimental material validation.

This delivery retains reference geometry and inventory. It creates no
separated fragments, new contact surfaces or transferred topology/history.
The band plane stays fixed and kinematics remain small. Local positive
softening does not remove structural snap-back or the existing solver's
convergence/refusal limits. The wider MR63 post-failure transfer scope remains
open.

## Source-bound transverse shear in native frames

`Frame2d::add_material_shear_element` explicitly enables an elastic Timoshenko
member using the original `ConcentricShearSection`, including its authored
shear correction and every original layer/source. The original force-based
state determination adds the complementary shear flexibility
`(GA L)^-1 [0,1,1] [0,1,1]^T`, with GA the corrected stiffness. Accepted
shear resultant, effective engineering strain and reversible energy are
observed through that same law. Existing Euler–Bernoulli members retain their
original behavior.

Consistent reference inertia integrates the static Timoshenko translation and
independent section-rotation fields with exact positive four-point quadrature.
The same native stiffness and inertia feed static equilibrium, modes, Newmark
and harmonic/pole preparation. This is a convergent beam approximation;
neither exact continuum dynamics nor resolved three-dimensional shear stress,
shear failure or thermal feedback is claimed. Explicit shear viscosity is
available through the source-bound frame strain-rate path below. Legacy
Euler–Bernoulli strain-trace consumers, including longitudinal Maxwell/thermal
adapters, refuse shear-enabled members rather than silently misinterpreting
their curvature. Frozen source-temperature and small-strain limits still apply.

The focused `frame2d_shear` target supplies analytic cantilever work/compliance,
load reversal, pure bending, rigid-motion inertia, a shared modal/harmonic/
Newmark oscillator, the Euler–Bernoulli limit, a three-rung simply-supported
Timoshenko frequency ladder, and source-refusal rollback checks.

## Graded reference inventory and consistent membrane modes

`MembranePatch::with_reference_areal_density` retains three positive P1
reference areal densities [kg/m²] in each original triangle's connectivity
order. This explicit inventory overrides the source sheet's scalar areal
mass for inertia and acceleration loads; it does not change its original
elastic law, grain, thickness or source authority. Incident triangles can
retain different values at a shared node. Original scalar fields preserve
the existing lumped arithmetic exactly.

Positive row sums integrate density times each shape exactly, preserving
mass and first moments. `MembraneMassImage::Consistent` integrates the full
cubic product of P1 density and both P1 velocity shapes, retaining inter-node
inertia. The same kernel supplies `reference_mass_apply`,
`solve_vibration_with_inertia` and `modes_from_tangent_with_inertia`. The
native modal owner still uses the original material/geometric/follower
stiffness, support elimination, complete spectrum and work-conjugate pressure
ports. Existing modal and nonlinear-motion entries keep their explicit lumped
approximation; no mass or frequency is fitted. Mass matrices depend only on
reference geometry and the retained inventory.

`current_corner_areal_density_kg_m2` observes rho_area_reference/J on the
original admitted P1 surface deformation, with J its current/reference area.
It never writes that current density back into reference inertia. This is
areal mass, not inferred bulk density or through-thickness rotary inertia.
Density transport, absorption, evaporation and evolving material inventory
are outside this fixed-reference model. Density admission is not independent
measurement authority or empirical qualification.

`tests/membrane/density.rs` adds G1/G3 exact cubic matrix, mass/centroid and
acceleration-load moments; reference/current area conservation; equal-total-
mass profiles with different native supported modes; off-diagonal modal
inertia and pressure work; partition/permutation invariance; original uniform
bits through modes and nonlinear motion; and typed input/cancellation refusal.

## Sparse consistent reference inertia for nonlinear membrane motion

`MembraneMotion::new_with_sparse_consistent_inertia` extends the existing
nonlinear momentum-Verlet lane to explicitly admitted sparse mass blocks,
beyond the dense entry's 512-free-node scalar leaf. The original exact P1
reference-density element kernel supplies every physical mass contribution;
deterministic COO assembly preserves element accumulation order before
reference row sums scale solver coordinates. No density averaging, diagonal
replacement, resampled inertia or substitute force law enters this path.
For positive P1 triangle density, the exact quadratic forms obey
`D/5 <= M <= D`, with D the positive reference row sums; assembly and stationary
support elimination preserve this bound. This controls the scaled physical
mass independently of density contrast in exact arithmetic. Nonfinite values
or a numerical factor without strictly positive inertia refuse.

One original `fs-sparse::direct::SymbolicLdlt` analysis and LDLT factorization
is prepared per distinct Cartesian support pattern, using the caller's
existing Natural, AMD or nested-dissection ordering. Nested dissection receives
the original reference node coordinates. Equal patterns and cloned checkpoints
share immutable numeric factors and their prepared solve permutations. The
existing `LdltFactor::solve_in_place` acts directly on the motion's retained
axis-major output slice; warm inverse actions allocate no mass matrix, gather
buffer or solve workspace. `sparse_inertia_statistics()` reports distinct
factor storage and per-inverse factor nonzeros, counting shared factors once
for storage and once per active axis for work. These are counts and original
numeric-owner payload bounds, not measured throughput.

`SparseConsistentMassBudget` explicitly admits free nodes and original staged
element contributions before sparse assembly. Retained symbolic bytes,
predicted factor nonzeros and the existing no-delay structural flop estimate
are checked after synchronous analysis and before numeric factorization.
The original direct owner's strict numeric-payload cap includes its fronts,
updates, factors, scratch and prepared solve permutations. Actual factor
nonzeros are checked again before accepting the prepared motion; two
triangular traversals plus linear permutation/scaling work bound each warm
inverse action. At most three different factors are retained. These per-factor
limits are not a whole-process RSS cap, a symbolic-construction peak-memory
cap, or a wall-clock guarantee. Cancellation polls cover assembly and the
boundaries of original synchronous COO/ordering/analysis/factor/solve leaves;
those leaves do not provide mid-operation interruption.

The same physical initial momentum, kinetic energy and polarized energy
increment remain in use. Stationary support impulses retain the original
`M_sf * delta_v` contribution and endpoint force kicks. The source force,
accepted-state publication, checkpoint replay and failure/cancellation
atomicity retain their original owners. The lumped and bounded dense public
entries retain their previous arithmetic and selection behavior.

Four focused regressions in `tests/membrane/sparse_inertia.rs` assert
all-ordering dense parity for graded nonlinear supported motion, full-space
support impulse balance, sparse factor sharing and uniform acceleration on a
625-node graded sheet, all explicit allowance refusals, late cancellation and
checkpoint replay, and finite kinetic energy across 560 orders of reference
mass ratio. These manufactured cases do not establish measured-material
fidelity, mesh/time convergence, automatic stability or real-time performance.

## Consistent reference inertia in nonlinear membrane motion

`MembraneMotion::new_with_consistent_inertia` uses the original nonlinear
momentum-Verlet map with the exact density-weighted P1 reference mass.
Each Cartesian component prepares its scalar free-node block through the
existing cancellable `fs-la` Cholesky owner. Positive reference row sums
scale solver coordinates before factorization; all off-diagonal physical
mass terms remain present. Equal support patterns share one immutable
factor. The explicit `ConsistentMassBudget` and hard 512-node per-axis leaf
cap bound dense preparation and triangular solves. The explicit sparse entry
above owns larger sparse inverse preparations; automatic stability bounds
are not claimed by either entry.

Initial velocities become physical nodal momentum through the same original
consistent matrix, then are recovered with the prepared inverse; factorization
roundoff can change their last bits. Kinetic energy and its polarized increment
use that inverse consistently. The original lumped constructor retains its
arithmetic and avoids the consistent image's extra work buffers.

Stationary supports retain their original values, but their reaction includes
`M_sf * delta_v` as well as the original force kicks. The summed supported-row
weights are prepared once; delta_v is recovered from accepted delta momentum
so a resolved weak impulse does not disappear under a large velocity datum.
The same integrator workspace is reused. Density moments, factors and source
laws are never resampled during motion; the new inverse path adds no warm
allocation. Failure and cancellation publish no candidate state or cached force.

Four further G1/G3 regressions in `tests/membrane/density.rs` cover the
independent nonlinear graded-mass period with second-order time refinement,
physical support impulse/kinetic balances, free linear/angular momentum and
time reversal, checkpoint/cancellation and work limits, and finite kinetic
energy for light/heavy sheet inventories across 560 orders of mass ratio.
No measured-material, real-time throughput or continuum-error qualification
is inferred from these manufactured cases.

## Exact pressure traction over finite membrane displacement

`pressure_segment_load_n` integrates the original current-area pressure force
over a straight nodal-displacement segment. It accepts the same scalar or
material P1 pressure data and retains the original pressure moments. The area
vector is quadratic along that path; the same analytic mean is shared with
`pressure_work_j`. Dotting the mean nodal load with the complete displacement
therefore reproduces the original path work up to floating-point accumulation.
For uniform pressure and fixed perimeter, unit pressure supplies a discrete
gradient of the exact swept volume.

`pressure_segment_tangent_apply` differentiates that mean force with respect
to the AFTER endpoint, keeping the BEFORE endpoint and pressure fixed. At
coincident endpoints it is one half of the instantaneous follower tangent.
A coupled midpoint coordinate with `q1=2*q_mid-q0` must apply its own factor
of two. Both endpoint areas, complete finite data and cancellation are admitted;
these operators do not certify intermediate constitutive states or make
nonuniform pressure conservative.

Two focused G1/G3/G4 pressure-kernel regressions compiled and passed using the
actual production `pressure.rs` with geometry-only fixture containers. They
cover cubic fixed-rim swept-volume work that midpoint geometry misses, exact
Gauss integration, reversed paths, signed P1 force derivatives, the coincident
endpoint factor, rigid translation and refusal. Complete material/motion
integration remains a separate dependency-limited validation obligation.

## Spatially resolved follower pressure on nonlinear membranes

`MembranePressure::TriangleLinear` borrows three signed pressure-jump samples
for each original triangle, in that triangle's original node order. Samples
remain fixed in material/barycentric coordinates during a call; a spatial fluid
or hydrostatic source must resample its Eulerian field between calls. Different
triangle traces are admitted. No pressure interpolation from a foreign mesh,
EOS, fluid inertia or constitutive pressure derivative is inferred.

The original pressure owner integrates the P1 pressure and displacement bases
exactly: each local nodal force is the current oriented doubled-area vector
times `(2*p_i+p_j+p_k)/24`. Its matrix-free tangent differentiates that same
area vector. Straight-drift work integrates its quadratic area evolution
analytically with the same nodal pressure weights. This retains pressure moments
and local virtual work which replacing a field by its mean would discard.
Nonuniform material pressure can perform nonzero work on a closed displacement
loop even with a stationary outer perimeter. Only a uniform fixed-perimeter
load has the existing conservative pressure-volume interpretation.

`pressure_load_n`, `pressure_tangent_apply`, `pressure_work_j`, and both
native `step_with_pressure` methods accept the new pressure type or a scalar
through `From<f64>`. Constant triangle triples retain the original uniform
arithmetic. Dynamics admits the immutable field before trial work, then reuses
its prepared pressure and drift storage through both current-area kicks.
All original inertia images, support coupling, pressure-work accounts, material
damping and transactional publication remain shared.

G1/G3 regressions define exact force and moment under rigid motion, independent
tangent differences and Gauss path-work quadrature, a nonconservative fixed-rim
loop, a nonlinear loaded pyramid with lumped and consistent inertia, signed
material loss, uniform parity and source/cancellation rollback. The two
pressure-kernel regressions compiled and passed against the actual production
`pressure.rs`, with geometry-only fixture containers replacing the surrounding
material setup. The full membrane and damped-motion integration regressions
are authored but have not been executed against the complete fs-solid dependency
closure in this environment. This evidence does not establish that broader
integration, continuum convergence or a stability certificate.

## Nonlinear source-bound membrane motion

`membrane::dynamics::MembraneMotion::new` advances the original total-Lagrangian
membrane energy with `fs-time::verlet_momentum_step_with_displacement`. Reference lumped inertia,
grain, source laws and finite-strain geometric stiffening retain their existing
owners. Only free coordinates are integrated; stationary support values are
eliminated exactly, and inconsistent initial support motion refuses. Cloning
retains a complete deterministic checkpoint.

The two momentum kicks use the original material gradient. The first reuses
the accepted gradient of the immutable patch; the endpoint evaluates the full
nonlinear law once. Geometry, force and integrator work buffers are allocated
at preparation and reused. Failed state/domain checks or cancellation cannot
publish displacement, momentum, time, energy or an invalid cached force.

Each step reports stored and kinetic energy, held-dead-load work over the
resolved kinematic drift, support impulse from the same kicks, and the signed
energy/work defect. Work uses the original `h M^-1 p_half` increment, without
subtracting rounded absolute positions. Kinetic changes use the polarized
quadratic form so a small resolved transfer survives a large energy background;
strain-storage/integration error remains in the signed numerical defect. Stationary
supports do no work. This second-order symplectic explicit lane introduces no
physical damping and does not claim exact energy conservation, automatic time
step stability, real-time throughput, spatial convergence or experimental
validation. Moving supports, contact and material-history
evolution are outside this entry. The original native constitutive costs and
admission domains still apply.

`MembraneMotion::step_with_pressure` adds a held signed pressure difference
through the original current-area triangle owner at both force kicks. It uses
the selected lumped, dense-consistent or sparse-consistent reference inertia;
pressure changes neither inventory nor the material law. Positive pressure
acts along the current reference-oriented normal, initially +z. Loads at
stationary supports are retained in their impulse account. Zero pressure
preserves the original entry's arithmetic.

Pressure work has its own returned account, separate from dead loads, stored
energy and the signed numerical defect. The native exact quadratic-area
integral uses the original resolved kinematic drift, retaining physical work
even below the absolute-position ULP. One prepared node buffer serves both
endpoint forces and full-space drift; warm pressure steps introduce no global
matrix or pressure-work allocation. With fixed perimeter this is conservative
pressure-volume work. A free perimeter is also admitted for the positional
force map, with explicit path-dependent work and no symplecticity claim. No
gas inventory, acoustic stiffness, fluid inertia, continuum-error certificate
or physical reservoir evolution is inferred.

G1/G3 regressions in `tests/membrane/dynamics.rs` compare a four-triangle sheet
against the independent nonlinear Duffing energy-integral period on three time
rungs, including pure geometric restoring force at zero prestretch. Additional
checks cover source stiffness/mass scaling, linear and angular momentum,
time reversal, held-load work/support impulse, checkpoint replay and refusal
rollback. Ballistic work is also checked below the position ULP at translated
origins and against a large orthogonal kinetic background. These are
manufactured mechanics checks, not measured drum sounds.
Pressure regressions additionally compare changing-area kicks and exact drift
work, a signed-pressure nonlinear pyramid under both reference mass images,
complete support impulse, time reversal, refusal rollback and a damped scalar
source law with separate physical heat.

## Original interface friction in quasistatic obstacle contact

`contact::ContactProblem::solve_timed_from_with_sources` and
`ContactTrajectory::advance_with_sources` bind every active node to an original
`SourceThermalCoulomb` interface. Normal load, tangent frame and source coefficient
refresh together in the existing outer lag. Inner Newton solves retain the
existing convex frozen Huber potential. Final acceptance checks the refreshed
physical residual and returns the actual per-node source queries and law IDs.

Inputs use SI and the native plane-strain unit-depth convention: nodal normal
reaction is N/m, and an explicitly supplied tributary contact width gives Pa.
Temperature coordinates remain prescribed; sliding speed is the actual relative
tangential displacement over the physical interval, including surface motion.
`solve_timed_from_with_source_history` admits explicit per-node accumulated
sliding distance. Every lag queries distance at the same accepted history plus
the candidate's absolute relative tangential increment. Final acceptance therefore
closes distance-dependent running-in coefficients on the actual equilibrium slip;
Newton and lag rounds never add physical distance. `ContactTrajectory` owns this
history, starts at zero or an explicit `new_with_sliding_history` restart, and
commits distance, displacement, time and work together. Accepted resolved reversals
add path length; normal motion and inactive endpoints do not. The solution returns
the per-node distance increment, and final source samples retain the queried
distance. Missing active bindings, source-domain failures, unsupported pressure
coordinates and history-dependent sources on the stateless entry refuse. No
fallback coefficient supplies unavailable material data.

This lane requires equal static and kinetic source values. Huber regularization
permits creep and supplies no exact static hysteresis, thermal feedback or heat
allocation. Sliding history includes regularized creep and uses endpoint tangent
projection: it is not an exact arc-length integral and cannot reconstruct
within-step contact activation, tangent rotation or slip reversals. Distance,
secant speed and endpoint work require temporal refinement. Source failures and
unconverged lags publish no trajectory state or work. Constant sources preserve the original
coefficient path; focused G1/G3/G4 contact tests exercise native equilibria,
local temperature/pressure/speed response, final-query closure and rollback.

The same final Huber traction supplies `friction_work_per_node` in original
nodal order; its ordered sum is the existing `friction_work`. Inactive nodes
and disabled friction contribute zero. SI plane-strain work is J/m of physical
depth, so conversion to heat requires an explicit depth and heat fraction.
`ContactTrajectory` can be cloned for a private coupled trial while retaining
the immutable original problem and all accepted distance, displacement and
work. This permits joint mechanics/thermal acceptance without aging the
accepted interface during thermal iterations.

## Graded reference density on native HEX8 solids

`linear3::hex8::Hex8ElasticProblem::assemble_with_reference_density` accepts
eight positive finite corner densities in kg/m³ per element. Density uses the
original nonnegative Q1 interpolation, with discontinuities between elements
permitted. The existing shared isoparametric owner integrates consistent mass
and total inventory; Standard and B-bar stiffness retain their original paths.
Constrained elimination does not remove mass from total-inventory reporting.

Three Gauss points in each direction integrate `rho*N_a*N_b*det(J)` on a
trilinear map: its coordinatewise degree is at most five. The existing mass
rule therefore needs no increased quadrature work. An identical uniform
material-density profile retains the scalar path bitwise. The resulting mass
pencil feeds native modes and dynamics without replacing the graded field by
its mean. Reference inventory is held fixed through subsequent motion; this
does not introduce density transport, thermal inventory refresh, finite-strain
kinematics or global map-injectivity certification.

The `hex8_elastic` G0/G1/G3/G4 regressions compare every matrix coefficient
against independent monomial integration on affine and warped bricks, retain
centroid and constrained-mass moments, verify native eigenvalue scaling, and
exercise invalid density/count and cancellation refusals.

`hyper3::plastic::hex::dynamics::FiniteThermoplasticHexMotion::try_new_with_reference_density`
carries the same positive corner-density representation into native finite
thermoplastic motion, point mobility and impulses. Its independent 27-point
reference mass rule is exact on warped trilinear cells; the original eight
constitutive points, Standard/mean-dilatation formulation and all histories
remain unchanged. The explicit profile and immutable consistent mass survive
heating, plastic flow, finite deformation and impulses. Reported cell densities
are volume-weighted inventory means, never substitutes for the graded inertia.
`reference_acceleration_load_n` projects an original-node acceleration through
that same mass, preserving density-weighted force and moment at supported nodes.
The scalar constructor shares this exact reference rule; a constant explicit
profile reproduces its mass bits on both affine and warped bricks.
`hyper_plastic_hex_dynamics/density.rs` checks independent
polynomial inventory, first/kinetic moments, equal-mass redistribution changing
supported mobility, hot-history continuation and atomic refusal. This is
authored reference data, not material calibration or mass transport.

## Prestressed Standard and B-bar HEX8 vibration

`linear3::hex8::Hex8ElasticProblem::solve_thermal_tangent` and
`solve_thermal_tangent_with_nodal_forces` expose the existing conservative
initial-stress pencil for eight-node bricks. `Hex8PrestressedAssembly` names
the shared immutable response/pencil carrier. Thermal strain and explicit
dead forces first enter the original equilibrium owner; an independently
supplied stress or stale response cannot enter this path.

The shared prestress kernel integrates the accepted compatible or B-bar
stress at every original stiffness point. B-bar uses the same volume-weighted
strain projection as equilibrium and recovery; geometric stiffness uses the
physical basis gradients. Neither raw compatible stress nor element-average
stress replaces the accepted quadrature field. Reference consistent mass,
elastic stiffness, geometry and homogeneous incremental supports retain their
existing owners. Fixed-DOF forces retain their reaction meaning.

Compression remains signed, and fully constrained bodies return `NoFreeDofs`.
This is a frozen-reference initial-stress approximation, not a finite-strain
tangent, follower-load model, stability certificate or live thermal update.
The existing `hex8_thermal` tests check independent anisotropic polynomial
geometric work, signed thermal/mechanical shifts, unchanged mass and supports,
reactions, zero-stress recovery, cancellation and fully fixed refusal.

## Native plane elastic stress, storage and support reactions

`LinearProblem::response` and `response_with_element_scales` recover original
P1/Q1 quadrature strains, full three-dimensional Cauchy stress, elastic storage
and pre-elimination nodal residuals without assembling a stiffness matrix.
Standard and B-bar response use the same strain columns, original loads and
supports as the solve. Plane stress retains relaxed thickness strain; plane
strain retains transverse stress. B-bar reports compatible in-plane strain
separately from its two-normal assumed strain.

Forces and integrated storage are N/m and J/m of reference thickness. Body and
traction loads on constrained rows enter support reactions before elimination;
free-row residuals remain visible. Support entries are physical reactions only
for an equilibrated displacement satisfying the original prescriptions.
Recovery subtracts a uniform nodal translation before differentiating, and
polls cancellation between elements and loaded edges without publishing partial
results. Ordinary/scaled solves and recovery share material, geometry and
support admission. This adds no finite-strain, calibrated strain-domain or
continuum-convergence claim to the existing linear-kinematic model.

`tests/linear_response.rs` covers analytical plane-stress/plane-strain extension,
shear, heterogeneous moduli, B-bar energy derivatives against the original
stiffness, support force/moment balance, translation and typed refusals.

## Static elasticity preserves the original force scale

The plane P1/Q1 and native P1 tetrahedral, TET10, HEX8, HEX20 and HEX27 static
thermal/mechanical solves use `fs_solver::checked_cg` with their original
stiffness and Jacobi preconditioner. The shared solver normalizes only its
iteration right-hand sides, restores physical displacement, and checks the
returned result against the original force equation. A small recursive residual
alone cannot publish a solution. Up to two defect corrections share the original
total iteration cap; true-residual tolerances and source/material/load identities
retain their existing meanings. Context-bearing solves poll cancellation through
the same owner, including before publication.

Plane prescribed displacements enter the free force through the existing exact
column elimination. Their auxiliary identity rows remain homogeneous during
iteration, so metres never mask a force residual; the prescribed values are
restored exactly afterward. Plane-stress Lame reduction is evaluated directly,
without squaring the modulus or passing through a singular three-dimensional
bulk modulus. No stiffness, strain, constitutive law, thermal load, reaction,
or energy is rescaled physically.

Regression cases cover signed traction and stress-unit rescaling over ordinary,
tiny and large finite magnitudes, exact zero load, nonzero supports, and original
thermal/mechanical stress, energy and support reactions across element families.
They keep large-force strains small by scaling the modulus with the force. These
are numerical/model identities, not measured extreme-material applicability or
an additional finite-deformation claim.

## Native interval work in deposited strips

`deposition::LayeredStrip` integrates each original material point's
`work_conjugate_stress_pa` into axial force and bending moment for work. Actual
endpoint resultants and reactions remain separate outputs. Prescribed strain
and curvature therefore support the original `ExactStrainRamp` history without
inventing backward-Euler numerical dissipation or replacing physical dashpot
heat. New layers retain their stress-free reference at the accepted shape;
existing branch strains and stored energy survive birth unchanged.

`MeanForceN` and `MeanMomentNm` prescribe interval means through the original
affine mean-stress tangent and require exact-ramp integration in every layer.
They do not prescribe a constant load throughout the interval. Existing
`ForceN`/`MomentNm` retain backward-Euler endpoint-work semantics and explicitly
refuse exact-ramp sources. Material coefficients remain frozen at the proposed
temperature; no exact continuously varying-temperature solution is claimed.
The existing bridge's endpoint dead-load convention is unchanged.

G1/G3 regressions check analytic elastic/Maxwell ramp work and heat, mean-load
response, later stress-free birth and relaxation, incompatible load-history
refusals, and unchanged backward-Euler force/work identities.

## Sustained normal and Coulomb forces on exact complete solids

`ExactLinearSolid::advance_sustained_contacts` advances the original held-force
trajectory under simultaneous normal and two-tangent Coulomb forces. The shared
`fs-contact::sustained::solve` owner closes actual endpoint gaps and accumulated
slip; every response starts from the same accepted native state. Original full
consistent inertia, stiffness and optional physical viscosity remain present.
Prepared matrices are borrowed, and native matrix contractions reuse bounded
linear scratch without another exponential, factorization or modal truncation.

Actual point traces and their transposes use one canonical contact-ID order.
Independent work checks separate external work, prescribed-surface work,
normal contact work, relative frictional loss and the original integrated
viscosity. Friction is not added again to the native viscosity account and no
thermal destination is inferred. The actual complete energy change and native
cumulative accounts must both pass before the original sole acceptance.

Initial impacts, negative actual gaps, inconsistent loaded dwell, unresolved
endpoint slip reversals and invalid sticking/sliding directions refuse.
A mandatory caller-owned whole-path admission checks interior nonpenetration,
fixed geometry, contact modes and any necessary determinacy evidence; endpoint
convergence alone is insufficient. This bounded local constant-force interval
does not discover events, claim global uniqueness, or resolve a changing force
history. Cancellation, budgets, source errors and late admission failures leave
physical state, preparation, clock and cumulative work unchanged.

Focused G1/G3/G4 regressions in `tests/solid_sustained_contact.rs` cover static
support and breakaway, circular two-axis sliding, translating surfaces and
source work, full cross mass, simultaneous contacts and ordering, analytical
elastic/viscous motion and continuation, impacts/stopping and transactional
refusal. Runtime verification is recorded with the implementation commit.

`CompleteSpectralSolid::advance_sustained_contacts` shares this same contact
solve and full-coordinate physical admission. Its original complete basis and
cached oscillator coefficients retain every elastic and rigid mode; no spectral
cutoff, contact compliance, matrix clone or new factorization is introduced.
Five reusable linear vectors supply each read-only force response. The ordinary
affine advance and these responses share the original basis/mode kernels; the
final native replay retains its original interval and cumulative work gates.
An approximate basis admitted during preparation cannot bypass those gates.
Native operators remain borrowed and each response's work is explicitly bounded.

The path witness retains the actual dense or complete-spectrum owner, and the
error preserves its native error type. Existing dense callers retain their
default types. Source work, moving-surface work, relative friction and physical
viscosity keep one shared audit. The conservative complete-spectrum owner adds
no constitutive heat. Six focused `spectral_sustained_contact` regressions cover
sticking/breakaway, circular sliding, moving surfaces and boosts, off-diagonal
consistent mass, all-mode continuation and basis order/sign changes, and atomic
refusals including the original physical work gate. Native execution is reported
with the implementation commit.

`SupportedSpectralSolid::advance_sustained_contacts` extends the same shared
owner to authored support-acceleration histories. Contact traces and external
loads use every represented coordinate, including prescribed supports. Their
actual velocities enter slip and contact work; loads at those rows change the
original reactions without inventing mobility or discarding physical work.
The native cubic free source retains both K_fc and consistent M_fc coupling.

Each response is an unpublished original supported interval with all free modes,
continuous displacement moments, original full-row reactions, and independent
support-actuator work. No prepared matrix is cloned or refactored. Native free
and full energy gates remain active; contact/path admission and cancellation
precede one infallible publication of both physical/modal states, support motion
and cumulative accounts. Explicit response-product and simultaneous scratch
limits include this original cubic staging. Ordinary supported advance uses
the same private trial/commit seam with unchanged physical formulas.

The shared receipt now exposes `support_work_j` separately from the opposing
contact surface's work, external work and relative friction. It comes from the
original reaction integral and enters complete energy once. Four focused
`supported_sustained_contact` regression groups cover friction on prescribed
coordinates under boosts, independently integrated accelerating-support work
with cross mass and load partition, native continuation, and transactional
refusal including the original whole-supported energy gate. Geometry/event
discovery and allocation of frictional loss to thermal owners remain separate.

### Distinct static and kinetic sustained friction

`advance_sustained_static_kinetic_contacts` on `ExactLinearSolid`,
`CompleteSpectralSolid` and `SupportedSpectralSolid` resolves separate static
capture and kinetic sliding strengths. Each contact's existing coefficient is
the static bound; a same-order kinetic slice must satisfy `0 <= mu_k <= mu_s`.
Any nonzero incoming relative tangential speed selects kinetic friction before
the first force response. A resting trace first attempts static capture. If
that coupled solution has resolved tangential travel or endpoint speed, the
trace switches to its kinetic bound and all coupled forces are solved again
from the same original native predecessor. Each contact switches at most once.

Every re-solve shares the original response-count and native-product budgets,
prepared matrices, source/support inputs and reusable response scratch. No
intermediate trial state, clock, energy or support account is published. The
receipt reports the selected coefficients and number of regime re-solves in
canonical contact-ID order; `solve.response_evaluations` and feasibility
refinements include every attempt. The final solve retains its Newton report.
Equal static/kinetic coefficients retain the original single-cone result.

The final physical mode/work gates and mandatory whole-path callback apply to
the selected friction law. The callback receives those actual coefficients.
Interior stopping, reversal, loading changes or recapture still require event
resolution; this held-force interval does not replace that geometric/temporal
owner. Relative frictional loss stays separate from physical viscosity,
surface work, support-actuator work and any downstream thermal allocation.
Focused native regressions cover static hold, breakaway acceleration/work,
pre-existing slip, translating surfaces, support reactions and actuator work,
equal-coefficient equivalence and atomic budget/invalid-input/path refusals.

### Located Coulomb stopping, recapture and reversal

`advance_collinear_coulomb_events` on `ExactLinearSolid` and
`CompleteSpectralSolid` advances through a located sliding stop and either
static recapture or reversed kinetic sliding. `CollinearCoulombContact` supplies
one initially closed, compressively loaded contact, a fixed tangent axis,
constant surface velocity and separate static/kinetic bounds. The complete
external load stays constant for the original native interval. Initial normal
motion must lie within the existing velocity allowance; any actual residual is
carried, and resolved inward motion requires the impact owner.

Each original contact mass-response column must lie exactly in the admitted
stiffness and viscosity nullspaces, with its full consistent-mass equation
independently checked. The spectral response must excite only original zero
eigenvalue modes. This admits affine contact velocities while retaining every
other original elastic and damped mode. Non-null contact responses, detachment,
resolved transverse motion and additional unresolved reversals refuse. This is
a bounded rigid-trace event path, not a general elastic-contact event locator.
The supported extension below retains constant support acceleration; general
support histories remain with the existing supported interval owner.

Incoming slip uses kinetic strength immediately. Normal dwell and its kinetic
boundary give a force warm start and analytic stopping time; every actual
segment still uses the original shared force-law solve and fresh replay. A
resting trace or located stop first tries the static cone, then re-solves with
kinetic strength if capture fails. At most two physical phases and three force
solves share one response/product budget. Discarded capture trials count. No
phase changes the accepted state or advances the native clock independently.

The original held external-load flow is evaluated once. Contact force impulses
and first time moments provide its rigid-coordinate correction without another
exponential, a velocity reset or a fictitious impact loss. Segment integration
retains the low component of endpoint time differences. The actual stop-speed
residual remains visible. Normal extrema check both penetration and opening;
tangent extrema check static excursions and unresolved reversal travel against
both geometric and work allowances. Small endpoint or net-displacement errors
cannot conceal a resolved interior friction loop.

The receipt records actual phases, stopping time, force solves, gross sliding
travel, full consistent momentum change, external work, prescribed-surface work,
normal work and relative friction work. A reversal's two signed work segments
remain distinct, even when its net travel cancels. The dense owner's original
viscosity integral is unchanged because contact corrections lie in C's
nullspace; friction loss is not added to viscosity or assigned to heat.
Complete-coordinate energy, native interval/cumulative gates, the mandatory
geometry/feature callback and final cancellation all precede one original atomic
acceptance. The spectral owner updates actual modal and physical states together.

Focused native tests cover stopping and repeat-step capture, moving surfaces,
static hysteresis, reversal loss, off-diagonal mass, unchanged elastic/viscous
spectators and modal continuation, unsupported traces, global budgets, geometry
refusal and cancellation. Interior normal opening and hidden static/kinetic
out-and-back excursions have direct regression checks.

### Coulomb events with accelerating native supports

`SupportedSpectralSolid::advance_collinear_coulomb_events` extends the same
located stop/recapture/reversal resolver to prescribed supports whose retained
acceleration stays constant during the original interval. Contact rows and
external loads use all represented coordinates. The three free contact
responses must remain independent rigid directions of the original full
stiffness; other free elastic modes retain their original polynomial forcing.
Changing support acceleration or a nonrigid contact trace refuses atomically.

For each contact force component, the original free mass response is
`B = [M_ff^-1 J_f^T; 0]`. Its additional support reaction is
`R_c = M_cf B_f - J_c^T`; the full original equation `M B = J^T + R` is checked.
Thus consistent cross mass and contact loads directly applied to prescribed
coordinates both remain active. Constant support acceleration contributes
`-R_c dot a_c` to the contact acceleration. The full momentum audit includes
the independently derived actuator reaction impulse.

Actuator work adds the actual phase integral `integral R_c F dot v_c dt` to the
original native support-reaction work. Each support's quadratic displacement
supplies that integral, including either side of a friction reversal. Free
source work separately retains `-M_fc a_c`; contact work on supported rows is
removed from the free subsystem's port without dropping it from the full
solid. Neither source is inferred from an energy difference. The event receipt
exposes `support_work_j` independently of prescribed-surface work and friction.

The corrected native polynomial trial retains actual modal/physical position,
velocity, endpoint acceleration and all three continuous displacement moments.
Original free energy, full support reactions, free-row force residual, interval
and cumulative work gates precede geometry admission and one atomic publication
of both owners. The original supported advance uses the same final gate. Native
tests independently integrate actuator reactions for capture and reversal,
check cross-mass momentum and common velocity boosts, preserve a driven elastic
spectator and subsequent modal continuation, and reject changing acceleration,
nonrigid traces, spent budgets, geometry changes and cancellation. A direct
piecewise-polynomial oracle checks the three force-induced displacement moments.

### Sustained directional Coulomb ellipses

`advance_sustained_elliptic_contacts` on `ExactLinearSolid`,
`CompleteSpectralSolid` and `SupportedSpectralSolid` resolves an oriented
two-principal-coefficient friction ellipse through the same original held-force
flow. The shared fs-contact coordinates use physical force `f = A q` and
conjugate point trace `A^T J`, preserving every original normal/tangent and
cross-contact mechanical coupling. The existing simultaneous force solve,
native inertia/stiffness/viscosity and complete support work remain authoritative.
The extra staging is bounded O(contact count * original DOFs); no exponential,
factorization, modal cutoff or artificial contact compliance is added.

The mandatory callback receives original physical contacts, forces, actual
gap/slip and original native motion. Before calling it, the adapter reconstructs
the original transpose load and checks physical relative and moving-surface
work against the transformed response. Sticking path/velocity tolerances are
conservatively tightened by the smallest tangent scale. Transformation,
physical-work, budget and callback refusals all precede native atomic acceptance.
The result explicitly separates physical forces/paths from the transformed
solver receipt; native motion and work accounts retain their original units.

Equal coefficients preserve the existing isotropic response exactly, including
the frictionless case. Unequal coefficients must both be positive; a single
free-slip axis is refused rather than assigned an artificial coefficient.
Principal directions and coefficients are held for the interval. This lane
does not discover geometry/events, infer temperature or wear evolution, or
allocate mechanical frictional loss to heat. The caller still admits the full
trajectory and fixed directional law. Analytical point regressions specify
directional breakaway, static holding, oblique maximum dissipation, translating
surface work, tangent-frame covariance, isotropic reduction and atomic refusal.
Runtime verification is recorded with the implementation commit.

### Static capture and kinetic breakaway with directional friction

`advance_sustained_static_kinetic_elliptic_contacts` on the same dense,
complete-spectral and supported owners accepts a static ellipse and a kinetic
principal-coefficient pair on each contact's unchanged material axes. Each
kinetic strength lies between zero and its static counterpart. The two axes
can weaken by different proportions; the full work-conjugate map changes with
the selected ellipse, preserving its physical shape and original trace.

Any actual incoming tangential motion selects kinetic strength immediately.
An initially resting contact first attempts static capture. Resolved slip or
endpoint velocity discards that numerical trial, switches the contact to its
kinetic ellipse and re-solves all coupled forces from the same accepted native
predecessor. Selection occurs before physical mode admission, so an intermediate
static trial cannot wrongly reject the eventual kinetic trajectory. Each
contact switches at most once, and all coordinate staging, selection, response
evaluations and native product work share bounded budgets. No prepared matrix,
accepted checkpoint or physical clock is cloned or advanced during retries.

The returned `friction` array exposes actual physical principal strengths and
axes in canonical contact-ID order. The transformed receipt counts every force
response and feasibility refinement, including discarded static trials, and
reports the number of coupled regime re-solves. All original endpoint, work,
support-reaction and mandatory whole-path checks precede one native acceptance.
Exact isotropic pairs recover the original circular static/kinetic trajectory
bit for bit on the admitted fixtures. Rank-one ellipses, interior stopping,
reversal, changing texture or coefficients, and thermal allocation remain
outside this held-law interval's domain; existing event/source owners retain
those responsibilities.

G1/G3/G4 native regressions cover nonproportional breakaway and subsequent
kinetic continuation, oblique maximum dissipation, mixed static/kinetic contacts
with full normal/shear and cross-contact inertia, moving support and opposing
surface work under velocity boosts, elastic spectator/modal continuation,
isotropic equivalence, shared response budgets and late refusal/cancel/retry.

## Native quadratic Maxwell force, contact and pressure drift

`MixedViscoTetMotion::force_kick_candidate` integrates the accepted original
mixed Piola stress and authored nodal/reference-volume loads through the original
curved P2 quadrature and supported consistent mass. Signed corner weights, all
midsides, heterogeneous quadrature material and actual load moments remain
present. Compensated accumulation retains small loads between large opposing
inputs. The kick retains separate source, internal and support impulse/work
accounts; it changes velocity while sharing immutable geometry, pressure,
stress and Maxwell-history buffers with subsequent contact impulses.

`contact_drift_candidate` retains the post-contact velocity and prescribes every
displacement to its represented kinetic drift. The original mixed owner then
solves continuous P1 pressure with the original accepted histories frozen during
the solve. Physical storage remains remainder plus `U(J(p))`, including when the
weak pressure projection differs from pointwise `U(det F)`. Reactions from these
prescribed material-evaluation coordinates are not physical support impulses or
work. The private force-phase identity rejects skipped, superseded or already
advanced phases; public diagnostic mutation cannot change continuation.

Native relaxation, source displacement work, source/force splitting, endpoint
work defect and momentum residuals remain separate. The common drift receipt's
`native_loading_work_j` includes native remainder loading and reversible bulk
storage change. A noncancelling numerical-work allowance gates the split; no
numerical loss becomes heat. This first-order finite-compliance entry requires
stationary held velocities, caller-owned complete path admission and refinement.
Exact incompressibility explicitly refuses this unconstrained kinetic subflow.

Fully prescribed `OriginalMaterial` endpoints assemble only the original pp
tangent, retaining global pressure indices, scaling, factorization and native
residual/domain gates. Unused displacement Hessians and up blocks are omitted;
public full tangent evaluation and final native support evaluation are unchanged.
Existing assembly admission remains conservative. New physical endpoints own new
immutable buffers; velocity-only phases no longer copy every material history.

Focused native regressions specify curved heterogeneous force/support torque,
signed body loads, complete history replay, independent P1 pressure projection
and physical bulk storage, phase binding and rollback. Additional pressure tests
compare the original full block and residual derivatives for quadratic and
log-squared bulk laws. Exact-rational moment checks and pinned source formatting
passed; native Cargo compilation and Rust test execution remain unverified under
the workspace resource blocker. No runtime speedup is measured here.

## Original quadratic Maxwell momentum during native motion

`MixedViscoTetMotion::advance_with_momentum` advances the same original P2/P1
Maxwell solve and returns its physical momentum account about an authored origin.
The native consistent mass supplies actual linear and angular increments. The
exact external nodal vector already assembled by the mixed owner supplies applied
impulse and endpoint moment, including signed corner weights, all six midsides,
curved reference quadrature, distributed body force and repeated nodal loads.
Original support reactions supply the corresponding external support exchange.
Accounting performs no load reassembly, mass reconstruction or additional history
update. The existing `advance` follows the identical accepted-state path.

Backward Euler's angular integration defect remains explicit and separate from
physical torque and momentum residuals. A consumer must bound that numerical
defect as well as the residual; it cannot declare angular conservation from the
work identity alone. The same accounting applies to the native finite-compliance
and exact-pressure regimes. Exact endpoint volume nevertheless does not make a
backward-Euler displacement chord an instantaneous volume-tangent velocity;
pressure-constrained impact continuation retains its own incoming-velocity gate.

Focused native regressions specify affine and curved reference volumes/centroids
under uniform acceleration, original load torque, nonzero numerical angular drift,
signed supported quadratic loads on a deformed body and one retained Maxwell
history update. Independent polynomial moment arithmetic and source formatting
were checked; native Rust execution remains unverified under the resource blocker.

## Graded reference inertia in loaded mixed-solid vibration

`MixedUpProblem::solve_vibration_with_reference_density` and
`solve_incompressible_vibration_with_reference_density`, including their
`_from_with_reference_density` entries, retain four positive P1 reference-density
values per original Tet10. The cached reference quadrature integrates that field
against the signed consistent P2 mass before finite-compliance condensation or
exact-pressure nullspace projection. Changing a density distribution can change
the modes while total mass and mean density remain fixed. Pressure still has
zero inertia; equilibrium, material tangent, supports and gauge admission keep
their original owners.

Both returned spectra expose the original corner field, reference-volume cell
means and total reference mass. Their `reference_acceleration_load_n` applies
the retained full mass to authored nodal acceleration, including prescribed
coordinates and signed corner weights. Modal forcing selects the free force
entries; exact-pressure modal forcing additionally applies the original
nullspace transpose. The reduced mass does not replace that full load projection.
The full operator is assembled once under the existing mass budget and retained.
Graded curved quadratic maps require six Duffy points per direction; affine
maps use four. Uniform profiles preserve the original scalar arithmetic.

Focused analytical regressions in `tests/mixed_up/reference_density.rs` check
loaded finite-stretch modal frequencies at equal total mass, exact-pressure
gauge and mode normalization, full gravity force/centroid moments, every affine
and curved mass coefficient, scalar parity and explicit admission refusals.
These entries cover the original dead-load problem; existing dedicated gravity
and follower-pressure adapters retain their separately declared density APIs.

## Graded inventory in quadratic Maxwell motion

`MixedViscoTetMotion::try_new_with_reference_density` binds four positive P1
corner densities per original tetrahedron to the native signed consistent P2
mass. The original cached quadrature supplies geometry, displacement shapes and
barycentric density samples. Density and mass remain immutable across native
finite-compliance and exact-incompressible motion, impulses, force kicks and
drift. Graded affine maps require the existing four-point Duffy rule; graded
curved maps require six points per direction. Constant profiles preserve the
original scalar mass path's arithmetic.

`reference_corner_density_kg_m3` returns the complete retained profile.
`reference_density_kg_m3` reports reference-volume-weighted cell means; it does
not replace the density field in physical inertia. `reference_acceleration_load_n`
applies the same full consistent mass to authored nodal acceleration, including
held coordinates and signed corner entries. A held gravity or acceleration load
therefore preserves the original body's force resultant and torque without
changing material inventory as the body deforms. Evolving force fields still
require their own force-law owner.

Independent barycentric polynomial integrals check every affine and curved
mass coefficient. Focused native fixtures exercise the graded centroid,
finite-compliance and exact-pressure acceleration, momentum balance, retained
Maxwell history through impulse, scalar equality and explicit admission refusals.

## Current-area pressure on finite-deformation membranes

Source-bound hyperelastic sheets can explicitly select exact incompressibility
with `ResolvedSheet::with_incompressible_hyperelastic_membrane`. The existing
3-D `MembranePatch` consumes the native constrained material potential and
tangent through the same material-axis Green strain as its finite-compliance
path. At unchanged natural reference, normal thickness contracts inversely
with area. Transverse shears relax for tilted fibres and reference areal mass
remains unchanged. No large
bulk penalty, extra stiffness or inferred bending model enters the solve.
Dead-load equilibrium, current-area pressure and loaded small-vibration modes
therefore share the same explicit material choice, supports and reference mass.

Membrane kinematics retains the relative displacement gradient H before forming
F=A+H in the cached material axes A. Engineering Green strain is evaluated from
sym(A^t H)+H^t H/2, preserving small axial strain and quadratic transverse strain
without subtracting nearly equal reference metrics. The undeformed reference
has exactly zero total strain; any natural-reference stress still belongs to
the original material. Elastic forces, matrix-free tangents, vibration admission,
viscosity and current-density recovery share this kinematics and strain-domain
admission. The `membrane` regressions check signed affine extension, quartic
transverse energy and its force/tangent, and the absence of artificial transverse
stiffness at rest. This is a numerical kinematics correction, not a new bending,
wrinkling, material calibration or finite-strain constitutive model.

The `hyperelastic_membrane` regressions differentiate the native assembled
energy/force/tangent for both choices, including tilted fibres. An independent
exact-incompressible Neo-Hookean pyramid potential supplies center displacement,
energy, support reactions and the loaded modal eigenvalue under fixed-rim
pressure. These are discrete manufactured controls; they do not establish
continuum convergence, material calibration, wrinkling or finite-amplitude
vibration. Full native test execution requires the complete workspace.

A source-bound hyperelastic sheet may also carry an explicit finite natural
reference, including three retained mean-thermal-expansion endpoints. The
same native membrane mechanics then consumes its pulled-back constitutive
energy, resultants and tangent. The installed mesh, reference thickness and
areal mass remain original; only the elastic part Fe=F G^-1 sets constitutive
stretch and fibre activation. Exact incompressibility acts on Fe and permits
total det(F)=det(G). Source temperatures/axis roles and the endpoint material
context stay in the shared immutable material carrier. Native consumer tests
cover tilted-fibre frame/tangent behavior and an independently derived heated
pressure-loaded pyramid equilibrium and modal stiffness with unchanged inertia.
This is fixed-endpoint mechanical preparation; it supplies no heat integrator,
evolving thermal work or added wrinkling/bending law.

`MembranePatch::solve_with_pressure` uses the original membrane material,
kinematics and Newton--Krylov owner with a uniform held pressure jump. Pressure
acts on current oriented triangle areas, including the exact matrix-free
follower-load tangent and loads at prescribed nodes. Positive pressure acts
along the reference-oriented normal (initially +z); it is a signed pressure
difference, not absolute thermodynamic pressure. Nonzero-pressure equilibrium
requires a consistently oriented surface and a fully prescribed perimeter.

`pressure_work_j` integrates current-area virtual work exactly along a straight
nodal-displacement segment. At fixed perimeter this is pressure times swept
volume. The equilibrium potential uses zero free displacement at the same
prescribed perimeter as its pressure-work datum. Moving-edge segment work is
path-dependent and is not advertised as a global potential.

`solve_vibration_with_pressure` solves that same equilibrium and uses
`K_material + K_geometric - dF_pressure/du` with the conserved reference mass.
Its additional pressure/volume port reuses the exact current-area load owner.
This holds the background pressure fixed; it does not infer gas compressibility,
reservoir volume, flow, radiation or a finite-amplitude dynamic trajectory.
The existing material-domain, cancellation, bounded solve and full-spectrum
admission rules still apply. G1/G3 tests in `tests/membrane/pressure.rs` check
analytic pyramid deflection/reactions/work, deformed-area frame covariance,
every pressure-force derivative and the loaded modal pencil with two interior
nodes (nonzero follower stiffness). No experimental membrane validation or
global injectivity/continuum stability certificate is claimed.

`MembranePatch::modes_from_tangent` is the shared lower modal owner for
composing physical models. It accepts a finite row-major full residual tangent
in an ordered unique free-DOF map, reuses original geometry/material admission
and reference mass, and retains the original symmetry, complete-positive-spectrum
and residual gates. Every eigenfrequency enclosure must also fit the common
support of the original areal-density and active elastic-law sources. This
admission applies before damping is attached; unsupported modes cannot escape
through the undamped entry. Per-area sheet laws do not borrow thickness support,
while bulk hyperelastic energy and volume viscosity retain its applicability.
Unrelated carried material channels do not narrow the spectrum. `MembraneModes`
owns the same pressure projection and
reconstruction as held-pressure `MembraneVibration`. This lower entry does not
certify force equilibrium or constitutive provenance of a supplied matrix;
the composing physical owner must supply the same accepted state, supports and
complete tangent. Existing held-pressure solve semantics are unchanged.

## Source-derived membrane viscosity and exact modal motion

`membrane::viscosity::MembraneViscosity` binds one reduced plane-stress
viscosity tensor per original sheet triangle, using the identical complete
material-state bundle, ordered warp/weft axes and reference thickness.
Engineering Green-strain rate gives objective resisting force and independently
integrated power over reference area times thickness. Rigid motions and
directional loss nullspaces are preserved. Flat prestressed transverse motion
has zero first-order material strain-rate loss despite positive geometric
stiffness. Admitted natural references use the original hyperelastic convention
Fe=F G^-1: rates transform by D=diag(g1^-2,g2^-2,(g1*g2)^-1), conjugate stress
by det(G) D, and heat by det(G). The original mesh volume and inertia remain
unchanged, and source Gram factors receive the same pullback. Identity references
preserve the original arithmetic. Thermal references require the explicit
`new_with_thermal_endpoints` binder: original full sheet bundles/card hashes and
ordered axes are checked, and the same six pinned viscosity claims are resolved
at the hyperelastic model's exact endpoint context. Endpoint passivity and band
admission run again; an expired source cannot switch to another hot claim.
`reference_laws()` retains reference receipts while `laws()` supplies the actual
endpoint coefficients. Reference mass/thickness remain unchanged. This binds a
fixed endpoint, not evolving temperature or a hot source/state replacement.

`MembraneModes` retains the exact accepted linearization displacement.
Projection keeps the complete modal resistance and source Gram rows, with
source-band admission covering every eigenfrequency enclosure. The opt-in
`viscosity::modal::PreparedMembraneDamping` reuses `fs-phs`'s exact held-force
quadratic flow and positive physical-heat owner. State is `[q, qdot]`; ports
pair generalized force with interval-mean modal velocity. Successful warm steps
reuse prepared buffers, and late cancellation publishes no partial outputs.
All retained modes and cross-mode terms remain; two states plus one force port
per mode limit this existing dense exact path to 42 complete modes. Larger
bases refuse rather than truncate. This is motion linearized at the original
geometry, not nonlinear damped integration or a guarantee that arbitrary
forcing stays within the declared source band.

`membrane::dynamics::viscosity::DampedMembraneMotion` also composes physical
loss with the original finite-amplitude `MembraneMotion`. It retains lumped,
dense consistent or sparse consistent reference inertia and all original
stationary supports. Each held-geometry source Gram channel has an exact
rank-one exponential and independently integrated nonnegative heat. Symmetric
forward/reverse channel sweeps around the original elastic Verlet step give
second-order split dynamics, with explicit time/energy-defect and source-row/
inverse-work budgets. There is no dense damping matrix or modal ceiling.
Lumped channel work is local and linear in triangle count; consistent channels
reuse the original mass inverse and retain M_cf delta_v in support impulses.
Prepared state/scratch buffers are reused; cancelled or refused trials publish
neither motion nor heat. Per-element and accumulated physical heat are separate
from the signed numerical defect. `step_with_pressure` composes those same
loss stages around the original current-area pressure step. Pressure work
occurs only during the elastic drift and is retained separately from physical
heat, including the pressure forces on supports. The scalar-source regression
checks both signs of pressure with lumped and consistent reference inertia.
This is not the full coupled exponential, a nonlinear stability/phase
certificate, fluid-reservoir evolution, material temperature evolution or an
inferred source-bandwidth qualification.

## Variable-amplitude critical-plane Miner damage on native J2 histories

The quasistatic and backward-Euler J2 fatigue owners now admit
`J2VariableAmplitudeGauge` alongside fixed-plane CMB, representative-cycle SWT,
and Fatemi-Socie gauges. Every accepted native endpoint contributes its actual
stress and total mechanical-strain tensor to the material-level
`CriticalPlaneMiner`; no external tensor history or reconstructed equivalent
stress is accepted.

Each gauge retains its original strain-life calibration/source binding and
declared finite plane set. Uncorrected, Morrow and Smith-Watson-Topper policies
reuse the existing bounded rainflow/Palmgren-Miner owner. Paired policies keep
their existing refusal for phase-lag, relaxation, or nested reductions that
would hide unmatched stress extrema. Sample and samples-times-planes work
budgets remain explicit.

Native mechanical advancement and every fatigue bank are one transaction.
Quasistatic block failure, inertial-step failure, cancellation, source mismatch,
or a late variable-amplitude gauge refusal restores the previously accepted J2
state and every fatigue history. Finalization is likewise atomic across all
banks. The variable-amplitude bank does not require or manufacture representative
cycle closure and is not reset when a separate SWT/Fatemi-Socie cycle window is
restarted.

This is a finite-plane multiaxial variable-amplitude life estimate under linear
Miner summation. It does not establish cyclic stabilization, optimize the
continuum plane sphere, infer unsampled extrema, model load-sequence interaction,
or evolve fatigue damage into stiffness, heat, fracture topology, or crack
growth. The existing fixed-temperature small-strain J2 and source-validity
limits remain authoritative.

## Fatigue on original frame fibers and material layers

`frame2d::fatigue::FrameFatiguePath` owns one native frame and fixed gauges at
existing member/Lobatto stations. Static load ramps, displacement-controlled
cycles and Newmark records observe every accepted step. The initial supplied
native state is observed once; damage before admission is not reconstructed.
Newton trials never enter a fatigue counter. A failure in mechanics, any gauge,
or cancellation publishes neither the requested block nor its fatigue histories.
The native frame and gauge banks are cloned once per requested block, and the
Newmark operator is prepared once per record. Bounded unresolved histories
survive blocks and path clones without retaining the entire loading record.

Metal fiber gauges use accepted `epsilon_0 - y*kappa` and the original committed
uniaxial stress. J2, Voce and Chaboche (one through eight backstress branches)
retain their actual plastic strain; Menegotto-Pinto retains its actual committed
hysteretic stress. No elastic reconstruction replaces plastic stress. The shared
bounded CMB counters support uncorrected, Morrow and paired SWT policies using
the actual accepted stress/strain endpoints. Raw native fiber laws do not carry
material-card provenance. An explicitly supplied `FiberMaterial` must match
the exact original law, and a source-backed CMB curve additionally checks its
original material card, state, typed operating coordinates and modulus receipt.
Authored curves remain empirical estimates. Concrete and time-dependent or
thermal fiber histories are outside this fixed-condition, rate-independent path.

Elastic gauges name a physical `[y,z]` point in one original concentric material
annulus. Each uses that layer's actual retained modulus, not an effective section
modulus, to observe signed normal stress `E_layer*(epsilon_0-y*kappa)`. The shared
S-N counter retains its declared mean-stress correction. Admission checks the
point's annulus, original source identity/state/query, material/process and
explicit typed absolute-temperature coordinate. This is an axial/bending normal-
stress estimate; it does not infer shear/torsion or multiaxial fatigue. Material
viscosity is refused because omitting its local stress would change the observed
stress history. Stress-free reference and source temperature remain fixed.

`FrameFatigueSummary` returns separate fiber strain-life and elastic stress-life
results together. `preview_finish` evaluates all endpoint half cycles without
changing the frame, counters or future continuation; `finish` publishes both
banks atomically. These observations do not introduce stiffness degradation,
heat generation, crack topology, a sequence-dependent damage law or a qualified
component service life. `examples/frame_fatigue.rs --demo` demonstrates the
native plastic loading path with deliberately synthetic CMB parameters.

## Critical-plane strain-life at original plastic shell material points

`plastic_shell::fatigue::ShellFatiguePath<L>` binds a native source-resolved
CST/DKT shell to finite-plane variable-amplitude strain-life gauges. The default
law is `ResolvedJ2Plasticity`; original Voce and one-through-eight-branch
Chaboche laws retain their complete native states through the same owner. A gauge
names an original element, one of its three surface quadrature points, and one
original layer/thickness quadrature point. Plane normals are fixed in that
facet's material frame. The complete accepted plane-stress strain, including
transverse contraction and signed plastic strain, is retained. Stress recovery
applies the original elastic map to accepted elastic strain; it performs no
return mapping or extra constitutive update, and does not subtract backstress
from physical fatigue stress. Native history admission checks every retained
branch at the observed material point. Membrane force and bending moment
resultants never replace the local tensor histories.

The existing critical-plane Miner owner counts each declared plane throughout
the complete observation and selects the greatest accumulated damage. It uses
the declared uncorrected, Morrow or paired SWT policy and retains their existing
phase/coincident-extrema refusals. Original material identity, typed absolute
temperature, modulus and source-backed CMB receipts must match the assigned
layer. Per-gauge storage/work limits and checked aggregate plane/residue/work
budgets apply before admission. These are bounded empirical predictions, not
new constitutive damage, stabilization or measured service-life claims.

Prepared geometry is shared across continuation and clones. The virgin native
zero-load state is observed once, followed by every absolute load-path target,
including holds. Static preload may continue under the original backward-Euler
inertia. Dynamic advancement retains that exact reference-mass owner, accepted
velocity and representable elapsed time; returning to statics after motion
begins refuses. Original physical kinetic and numerical-loss accounts pass
through unchanged. Native solution, fatigue banks and motion publish together;
source, solver, late-gauge and cancellation failures preserve the accepted step.
`preview_finish` evaluates a prefix without consuming it, while `finish`
requires completion of the declared path. Fixed source conditions, small
strain/rotation and native shell/reduction limits remain authoritative.

Three nonlinear regressions in `tests/plastic_shell_fatigue/cyclic.rs` exercise
borrowed Voce and one/two/eight-branch Chaboche sources against native exact
free-sided axial histories, including each tensor backstress and Morrow/SWT
damage. Original backward-Euler motion, physical energy, prefix continuation
and late-gauge rollback retain the full nonlinear history. Changing a branch
source refuses even when material labels, modulus and temperature match.
These synthetic comparisons establish implementation references, not measured
cyclic calibration or service-life qualification.

### Representative shear-cycle life on the native shell path

`ShellFatiguePath::new_with_fatemi_socie` admits original material-point
`ShellFatemiSocieGauge` windows alongside normal-strain Miner gauges, or alone.
The existing Fatemi-Socie owner evaluates engineering shear-strain paths with
its minimum-enclosing-circle amplitude, whole-cycle maximum normal stress and
original maximum-shear-plane selection. It consumes the complete accepted
facet-frame stress and mechanical strain, including the native plane-stress
contraction and every Voce/Chaboche history variable. No extra return map,
backstress subtraction or axial-to-torsional calibration substitution occurs.

Admission binds the original immutable material identity, typed absolute
temperature, shear modulus and initial yield reference; sourced calibrations
also match the exact card, operating query and selected E/nu/yield answers.
Combined gauge counts and checked aggregate window sample/projection/circle
budgets supplement each window's original limits. Ordinary native steps retain
shear histories in place, with accepted-prefix rollback on late refusal rather
than copying the entire retained trace at every endpoint.

Explicit restart seeds a new representative window at the current accepted
material point while preserving normal rainflow history, completed optional
shear Miner totals and native motion.
Full tensor stress/strain closure remains required. Static and backward-Euler
cycle endpoint methods accept mechanics and all observations together; failed
closure, calibration or cancellation preserves the previous solution, velocity,
time and counters. Finalization closes both kinds of observation atomically;
normal-history preview does not close or reinterpret a shear window.

Representative shear-life results remain separate from the optional fixed-plane
FS/Miner totals. Opting in counts each explicitly closed window once per plane;
repeated closure adds no damage and failed closure publishes no partial total.
Neither repetitions, arbitrary variable-amplitude shear cycles, stabilization
nor creep/thermal interactions are inferred. Focused synthetic regressions
in `tests/plastic_shell_fatigue/shear.rs` cover an independent engineering-shear
life formula, native conditioned Voce/Chaboche histories, original inertia and
clone continuation, late-gauge rollback, source mismatches and window budgets.

### Independently phased shell loading and native fatigue

`ShellFatigueLoading` supplies new reference prescribed values and dead nodal
forces/couples for the next original load-path entry. The original factor scales
both fields; unit factors admit arbitrary absolute endpoint histories. Values
follow the original support-DOF order, so changing the load path cannot replace
supports, geometry, source laws or the reference-inertia owner. New field counts,
DOFs and finite values are admitted before the prepared native solve; its
coordinate-scaling, overflow, rotation and material-domain checks still apply.

Static and backward-Euler loading methods reuse the existing prepared shell
geometry and accepted-state observation. Their cycle variants close shear
windows in the same transaction. Independent membrane, bending and shear phases
therefore produce actual native tensor paths without imported stress arrays or
a second time/constitutive solver. A failed endpoint or closure consumes no load
entry and preserves displacement, velocity, time, plastic state and all counters.

Three focused regressions in `tests/plastic_shell_fatigue/loading.rs` exercise
a nonproportional triangle in engineering-shear space with an independent
enclosing-circle/normal-stress/life oracle in both static and dynamic motion,
native independently forced support motion with nonunit factors and original
inertial energy, and invalid-field/late-closure rollback. Sample endpoints must
resolve physical extrema; these controls do not infer an unobserved load path
or establish that the complete mechanical state is periodic.

## Fatigue at fixed material points in spatial rods

`rod::material_field::fatigue::SpatialRodFatiguePath` observes chronological
accepted native rod equilibria at fixed segment/layer/`[y,z]` material points.
Each point retains its original concentric layer, material identity and typed
absolute-temperature coordinate. Local axial/bending stress uses the original
layer modulus and the native intrinsic-strain subtraction. Torsional shear uses
the original layer shear modulus and signed circular-section stress field. The
point remains attached to the material; it does not move to whichever surface
point currently has the largest stress. Fixed thermal expansion is removed at
each original layer, preserving residual stress from incompatible expansion.

Signed normal-stress gauges feed the shared bounded S-N/rainflow owner with the
declared mean-stress policy. An S-N scope names the layer material
"chemistry; phase", as fs-material's S-N source writes it and as the plate
and frame S-N gauges read it; a Findley scope names the chemistry. Findley
gauges retain the full local tensor and
both shear-vector directions over an explicitly closed representative cycle.
Findley reports a calibrated fatigue-limit utilization, not an invented cycle
life. Its finite plane search and enclosing-circle work remain bounded by the
original material owner. The original constructors retain an explicit
unassessed transverse-force limit: omitted forces remain inspectable, zero
requires vanishing transverse force, and a nonzero allowance is not a local
stress-error bound. `new_with_transverse_shear` and
`RodFatigueObservation::new_motion_with_transverse_shear` instead recover this
contribution through the same material-owned circular Saint-Venant field.
They require explicit original isotropic roles and homogeneous elastic
coefficients; no `Q/A` surrogate or heterogeneous averaging enters the tensor.
Both normal-stress life and Findley receive the combined physical shear at
their fixed material point. The sum of the separate whole-annulus transverse
and torsional strain peaks conservatively bounds the combined small-strain
domain. The native force, authored shear correction and storage remain with
the original rod; quasistatic section recovery introduces no stress-energy
replacement, cross-sectional inertia or end-effect solution.

The initial accepted state is sampled once. Explicit target-load blocks use the
original native equilibrium solver; intermediate Newton and continuation trials
are not fatigue samples. The caller must resolve physical extrema in the target
history. Mechanics and all fatigue banks publish together, so a late gauge
refusal or cancellation preserves the preceding accepted block. Unresolved
rainflow reversals survive continuation, cloning and Findley-window restarts.
Provisional endpoint evaluation leaves the live history available for further
loading. Material/source temperature and stress-free geometry remain fixed;
finite strain, plasticity, creep interaction, stiffness degradation and crack
growth are outside this elastic rod consumer.

`RodFatigueObservation::new_motion` instead starts at an explicitly supplied
native motion pose, with no phantom sample from the source's initialization
pose. It retains the immutable source owner, allowing constant-time admission
of subsequent observations. Only selected segment strains are recovered through
the original native kinematics. Source/temperature/geometry replacement refuses,
and a failed observation restores every counter and last published tensor. A
native time owner must join this observation to mechanical acceptance.

## Provisional observation of native fatigue histories

The quasistatic and inertial elastic/J2 fatigue owners expose `preview_fatigue`
on their accepted histories. It includes current endpoint half cycles and uses
the same fixed-plane life owners as finalization. J2 owners additionally expose
`preview_variable_amplitude_fatigue` for their critical-plane Miner banks.
Returned summaries describe a hypothetical finalized prefix;
the native trajectory and all live counters remain unchanged, so loading may
continue after a preview. No mechanical solve, extra stress/strain sample,
source requery or representative SWT/Fatemi–Socie cycle closure is performed.
Cancellation or unsupported provisional cycles leave the accepted state intact.
This is an estimate of the observed prefix, not a prediction of future loading
or a rule for changing stiffness, heat, fracture or remaining service life.

## Shear-strain Fatemi–Socie life from native J2 cycles

The shared quasistatic and inertial J2 fatigue observer now admits
`J2FatemiSocieGauge` alongside fixed-plane CMB and critical-plane SWT gauges,
or as the only fatigue bank. Both owners consume the original native stress
and total mechanical strain, including plastic history and unforced inertial
motion. Existing cycle closure, conditioning restart and finish operations
include both critical-plane banks; no alternative constitutive or motion
solver is introduced.

Admission binds the actual immutable material identity, temperature, native
shear modulus and initial yield stress. Source-backed calibrations additionally
match the original card, named state, typed operating coordinates and exact
selected E, nu and initial-yield answers. The evolving hardened yield surface
cannot replace the fitted initial-yield reference. Original material and
calibration evidence remain inspectable.

A native block and every fatigue bank publish together only after requested
closure succeeds. Failure restores mechanics, fixed/SWT snapshots, each
Fatemi–Socie accepted prefix, provisional result and current-window damage. Ordinary advances do
not copy the retained shear tensor histories. Restart samples the actual
accepted mechanical state and preserves fixed-plane counters and the global
sample budget. A decaying or ratcheting path that fails closure cannot silently
become a stabilized periodic-life result.

All four native owners, `J2FatiguePath`, `J2DynamicFatigueRun`,
`PlaneFatiguePath` and `ShellFatiguePath`, additionally expose
`enable_fatemi_socie_miner(maximum_cycles, cx)`,
`fatemi_socie_miner_summaries()` and `fatemi_socie_plane_damage(gauge)`.
Enable once before the first loading excursion of a fresh window, including
after explicitly discarded conditioning. Admission is atomic across gauges
and checks cancellation before publication; late enabling cannot reinterpret
an already observed trace.

The shared material owner adds one calibrated life fraction per explicitly
closed cycle to EACH fixed material plane, then selects maximum accumulated
damage. Changing loading direction never combines distinct per-cycle winning
planes. Every nonzero plane must fit the original finite torsional calibration;
zero shear contributes zero. The representative-cycle result retains its
original maximum-shear convention separately from the cumulative summary.
Native restart carries completed plane totals and seeds the next window from
actual accepted stress and strain. Closing the same window twice adds no
exposure. Failed mechanics, late-gauge calibration refusal, exhausted cycle
budget or cancellation retains all prior completed damage and native state.
Ordinary steps still do not copy the retained tensor trace.

Four additional native regressions exercise changing shear direction,
amplitude and opening stress, plane/shell restart, repeated closure, clone
continuation and cycle-budget rollback. Existing late-bank and inertial
rollback checks also compare per-plane totals. These are synthetic numerical
references; their source definitions do not constitute experimental evidence.

This extension predicts bounded empirical shear-controlled life and optionally
linear accumulated damage for explicitly resolved cycles at fixed source
conditions. It does not establish stabilization or experimental component
validity, automatically count multiaxial variable-amplitude cycles, model
load-order interactions, or evolve fatigue stiffness loss or crack geometry.

## Findley shear-path fatigue from native elastic and Newmark responses

`fatigue::TetFatiguePath::new_with_findley` and
`TetTransientRun3::with_fatigue_and_findley` attach Findley gauges alongside
fixed-plane S-N gauges, or use Findley gauges alone. Both feed the same shared
observer from the original native tensor stress recovery, including actual
inertial stress under zero external load and prescribed-support motion.
Physical tensor shear is retained in the declared mesh frame. The calculation
does not replace the elastic solve, Newmark step or material constitutive law.

A sourced Findley gauge requires the original immutable resolved material
bundle, with identity equal to the actual native elastic material identity.
Its original card, named material state and every calibration query coordinate
and quantity convention must match. The retained thermal state and fatigue
calibration must have the same temperature. A label or caller-supplied hash
alone cannot authorize a source-backed transfer; derived/oriented material
identities need their own explicit original receipt mapping before admission.

Explicit window closure evaluates the full stress path. A combined endpoint
advance and closure accepts native mechanics, S-N history and all Findley
gauges together. Any late failure or cancellation restores every observation
and leaves native motion unchanged. Ordinary steps append bounded tensor
samples without cloning the entire preceding Findley history. A restart begins
at the last actual accepted native stress, preserves S-N and mechanical history,
and does not reset the global sample budget. Completed Findley usage factors
are exposed independently from S-N/Miner summaries.

This is fixed-temperature, small-strain elastic fatigue-limit usage on resolved
closed representative blocks and declared candidate planes. It does not infer
periodicity for decaying vibration, count multiaxial variable-amplitude cycles,
predict finite life or couple fatigue to stiffness loss or crack propagation.
Synthetic native regressions cover rotating shear, unforced vibration,
source binding, restart and all-gauge rollback against the original solvers.

## Critical-plane SWT from accepted native J2 cycles

The quasistatic and inertial J2 fatigue owners can carry bounded critical-plane
cycle gauges alongside fixed-plane rainflow gauges, or use critical-plane
gauges alone. Both use the same private observer, original material-state
binding, CMB source checks and actual native stress/total-mechanical-strain
tensors. Plane projection and SWT inversion remain in `fs-material`; no
alternative constitutive or motion solver is introduced.

An explicit cycle-window restart samples the current accepted native solution
and discards only the preceding critical-plane window and its optional result. It permits
mechanical conditioning before the representative observation, preserving
displacement, velocity, plastic history and fixed-plane fatigue counters.
A closed stress/strain window is an admitted empirical input, not proof that
the evolving constitutive state has stabilized.

Ordinary advances retain every observation across processing chunks. Explicit
`advance_cycle` operations close the requested observation on a private native
candidate before publishing any motion or fatigue state. Nonclosure, nested
strain cycles, unsupported life, a late gauge failure and cancellation preserve
the previously accepted state. Finalization closes fixed and critical gauges
together. Critical-plane results are exposed separately from the existing
fixed-plane summaries.

This fixed-temperature, small-strain tensile-mode SWT path admits only resolved
representative cycles and declared candidate planes. Decaying free vibration,
virgin plastic transients or ratcheting paths that do not close cannot silently
become periodic-life estimates. Plane/time refinement and calibration
applicability remain the caller's obligations; there is no multiaxial
variable-amplitude counting, thermal-fatigue transfer, crack growth or
fatigue-to-stiffness/heat coupling in this extension.

## Conservative gravity on the current finite-strain configuration

With cellwise reference density, `hyper3::gravity` integrates a stationary
world-field potential, acceleration and analytic spatial gradient at four
current P1 tetrahedral quadrature points.
`PreparedHyperTetProblem::gravity_load` retains the original reference geometry;
the element-ordered densities are the same fixed reference densities consumed
by native `reference_mass_matrix`. Current volume never creates or removes mass.
Force, external-force tangent and potential use the same degree-two quadrature.
Uniform and affine acceleration reproduce the full consistent reference-mass
projection, including prescribed coordinates and material interfaces.

The problem and prepared `gravity_load_with_nodal_density` entries admit one
positive P1 reference density per original node. Their eight positive cubic
points (vertices of weight V/40 and face centroids of weight 9V/40) integrate
rho0*N_a*N_b exactly. Force, potential and analytic tangent share those fixed
reference mass weights, preserving graded center of mass and inertia moments.
Under admitted positive-J deformation this represents rho_current = rho0/det F: changed current
volume does not change inventory or replace the retained field by cell means.
Globally uniform nodal input retains the original four-point arithmetic.
`PreparedHyperTetProblem::gravity_potential_energy_change` integrates
field-owned endpoint differences with those same material points and fixed
weights. It sums differences in compensated precision before rounding, rather
than subtracting two independently rounded total potentials. The callback must
supply the original stationary potential difference, not midpoint force work.
This retains small gravitational work beside a large arbitrary datum when the
source supplies an analytic difference. It allocates neither nodal forces nor a
Newton tangent. Both cellwise and graded P1 density are supported; nonlinear
fields retain the same spatial-quadrature and source-domain boundaries.
G0/G3 unit regressions check exact P1 mass moments beside a 1e20 J/kg datum,
endpoint reversal, quartic potential differences and source refusal.

`solve_with_nodal_density_gravity` and its combined pressure entry use the same
native Newton, line search and reaction owner; the bound mass weights are
prepared once. Current vertices and face centroids must lie in the field domain.

`HyperTetProblem::solve_with_gravity` uses that owner at every Newton and trial
configuration, subtracts its analytic tangent from internal stiffness, includes
its potential in conservative globalization, and recovers reactions with all
supported weight. The existing load factor ramps gravity. The combined pressure
entry preserves current-area pressure and its conservative/nonsymmetric split.
If a large absolute gravity datum obscures a small potential increment, an
actual residual decrease is required; equilibrium tolerance is not relaxed.

`StationaryGravity3` can bind that same endpoint field to its own analytic
potential increment. `solve_with_gravity_source` accepts this source with either
elementwise or nodal reference density and optional native boundary pressure.
Conservative trials accumulate the source increment in compensated precision
in the existing force quadrature sweep, with no extra field sample, mass
binding or tangent assembly. That direct change replaces endpoint subtraction
in Armijo globalization; its tolerance scales with actual work, independently
of the potential datum. Sample-only bindings retain the original residual
safeguard. Domain refusals backtrack without switching laws, and cancellation
propagates. The source remains responsible for the consistency of its two
callbacks; this interface does not turn moving sources into static loads.

Callbacks must supply derivatives of one stationary scalar potential. Finite
values and gradient symmetry are checked; derivative consistency remains the
field owner's contract. Domain exits refuse or backtrack, never clamp. This is
an external field, not self-gravity, moving-source dynamics or a new material
law. Smooth nonlinear fields require spatial refinement; no quadrature error
certificate, singular-source resolution or mesh-independent result is inferred.
Focused G0/G1/G3 regressions in `tests/hyper3_gravity.rs` exercise heterogeneous
native mass, finite-difference derivatives, rotated geometry, analytic finite
stretch with pressure and support reactions, potential-datum invariance, and a
three-rung nonlinear quadrature ladder. These are synthetic numerical checks.
The nodal-density unit regressions additionally compare all cubic simplex
moments, actual shared reference-mass action, deformed gravity force/energy/
tangent and center-of-mass torque, nonlinear discrete-potential derivatives,
and scalar-path bits. `tests/hyper3_gravity/nodal_density.rs` adds the native
finite-stretch/pressure/reaction and public-admission regressions. The bounded
native source validation executes the three new assembly regressions and the
three shared-mass regressions; the full Newton integration regressions remain
unexecuted in the current dependency-constrained environment.

### Gravity-loaded finite-strain vibration

The complete-spectrum and sparse-window owners accept the same borrowed
stationary gravity field, either original elementwise or P1 nodal reference
density, and optional conservative closed pressure. Their common equilibrium
helper calls the native finite-strain load owner. At the accepted deformation,
the physical pencil is the original material tangent minus the pressure and
gravity external-force derivatives, each included once. Reference mass retains
the identical density basis used for weight; current volume and cell averaging
do not replace a graded material inventory. Both results expose the retained
gravity callback, including an explicitly supplied zero field.

The sparse path merges all three canonical rows directly and releases full
load matrices before entering the existing spectral owner. Signed negative
eigenvalues remain visible. Original field-domain, conservative-gradient,
pressure, material, cancellation and spectral-resource gates still apply;
these modes do not introduce evolving gravity, damping or a new spatial
discretization.

`tests/hyper3_gravity/vibration.rs` adds independent finite-height equilibrium
and modal-stiffness formulas for constant/graded density, uniform/affine gravity
and affine pressure, dense/sparse parity, unstable modes and refusals. The
unchanged production sparse merge passed bounded original-source tests against
independent dense arithmetic, including cancellation and nonfinite/symmetry
refusals. Full native equilibrium/API/eigenmode integration tests remain
unexecuted without the complete workspace dependency closure.

## Fatigue on source-bound inertial J2 motion

`plastic3::dynamics::fatigue::J2DynamicFatigueRun` owns accepted physical time,
displacement, velocity, native plastic history and bounded fatigue observation
while borrowing the original source-resolved cyclic plasticity problem. The
law parameter defaults to `ResolvedJ2Plasticity` and also admits native Voce
and Chaboche histories through `SourceBoundFatigueLaw`. It initializes through
the real zero-load solve, admits explicit velocity compatible with the original
supports, and freezes reference density from each original material bundle.
Every strictly increasing endpoint calls the existing `advance_backward_euler`
with actual endpoint loads. Native stress, reactions and inertial energy remain
inspectable; there is no alternate constitutive or dynamic solve.

Single steps and requested endpoint blocks publish only after every mechanical
increment and every fatigue gauge accepts. A late gauge failure, invalid time,
nonconvergence, cancellation or original strain-domain refusal leaves all prior
motion, plastic state, cycles and damage unchanged. Checkpoint clones retain
native motion and unresolved fatigue cycles together. Only the latest native
solution/step and bounded gauge residues are retained.

Quasistatic `J2FatiguePath` and inertial motion share the same private observer.
Uncorrected strain-life uses the scalar counter; named Morrow/SWT policies feed
actual simultaneous signed stress and total mechanical strain to the paired
counter, including its coincident-extrema and source-domain checks. A sourced
CMB curve must match the original material card, named state, external query
conventions and modulus claim. Caller-declared curves keep their explicit
identity/chemistry/process/fixed-temperature/modulus admission.

Backward-Euler numerical dissipation remains separately reported in the native
inertia receipt; it is not material heat or fatigue damage. Physical timestep
refinement and resolved turning points remain necessary. Miner accumulation
does not change stiffness, crack topology or thermodynamic state. The six native
regressions in `tests/plastic_dynamic_fatigue.rs` compare unforced vibration and
time refinement, cyclic J2/Prager response and inertial support reactions, named
mean-stress policies, continuation, source admission and atomic rollback.

## Strain-life from the original cyclic J2 constitutive path

`plastic3::fatigue::J2FatiguePath` owns the latest native solution and one
scalar or paired strain-life counter per fixed-plane gauge while borrowing the original
`SmallStrainTetProblem<L>`, with `ResolvedJ2Plasticity` as the default law and
source-bound Voce/Chaboche laws admitted without changing their native state.
The real virgin zero-load solve
is the first observation. Every original absolute load target goes through
`solve_from`; its plastic strain, accumulated isotropic history, Voce saturation,
Prager shift and every Chaboche recovery branch continue from the preceding
accepted state as applicable to the original law. The observer projects
the returned tensor-Voigt total mechanical strain and signed normal stress.
It does not substitute nondecreasing equivalent plastic strain, unsigned von
Mises stress or a second constitutive approximation.

Each calibration must match its original resolved material-state identity,
chemistry, process, fixed temperature and Young's modulus. Original material
receipts remain inspectable through the path. The admitted original cyclic
laws have no thermal eigenstrain; thermal-wrapped
laws and changing-temperature histories require a separately admitted transfer.
Every accepted element, including ungauged ones, must remain within the native
law's finite-state and componentwise small-strain envelope, including every
Chaboche branch's native history bounds. Physical fatigue stress is the original
solver's accepted stress tensor, with no backstress subtraction or second
constitutive evaluation. Complete source receipts, fixed modulus and initial
yield references bind the existing CMB, SWT and Fatemi-Socie calibrations.

`advance(count)` observes every intermediate target and commits the whole
requested block only after mechanics and all gauges accept. Cancellation,
nonconvergence, a later gauge refusal or an unsupported strain leaves the prior
mechanical and fatigue histories intact. Clones retain the complete current
state and unresolved cycles, without storing the full response trace. `finish`
requires completion of the declared observation and includes final half cycles
once. Processing windows do not create artificial cycles.

The CMB model and selected Uncorrected, Morrow or Smith-Watson-Topper policy
retain their empirical, fixed-plane scope. No stiffness, heat, fracture, failure topology, ratcheting
model or validated service life is inferred from Miner damage. The five focused
G0/G1/G3/G4 tests in `tests/plastic_fatigue.rs` compare the actual cyclic J2 coupon
with independent one-dimensional Prager equations and an analytical strain-life
damage sum, and check continuation, source admission, budgets and rollback.
Their data are synthetic numerical references, not held-out material evidence.

Four further regressions in `tests/plastic_cyclic_fatigue.rs` exercise Voce
against an independent scalar consistency equation and analytic strain-life
count, three-branch reversal against the native exact axial restriction,
variable-amplitude critical planes, original inertial stress/velocity/energy,
and late static/dynamic refusal preserving every branch and fatigue counter.
The existing fixed-plane, representative-cycle and variable-amplitude fatigue
algorithms retain their original calibration and closure limits; native cyclic
hardening alone does not validate ratcheting, stabilization or component life.


## Fatigue from accepted inertial stress histories

`loads3::transient::fatigue::TetTransientFatigueRun3` owns the original sparse
Newmark run and the same fixed-plane rainflow/Miner observer used by
`fatigue::TetFatiguePath`. `TetTransientRun3::with_fatigue` admits only step zero,
requires an explicit retained thermal field, checks every gauge's material-state
identity and calibration temperature, and samples the actual initial stress.
Every successful stationary- or moving-support advance recovers native endpoint
stress and samples every gauge once. Inertia, support motion, damping, thermal
eigenstrain, reactions and work remain in their original mechanical owners;
zero applied load does not imply zero vibration stress or zero fatigue.

Source-resolved S-N curves require the original immutable material bundle in
`TetFatigueGauge::material_source`. The shared observer verifies its actual
elastic identity, original material card, named state and complete fixed typed
query; only calibration life and R coordinates are excluded. Missing bundles,
foreign cards or states and different fixed physical context refuse before
observations begin. Both native paths expose the retained bundle through
`material_source(gauge)`. Caller-declared curves remain explicitly supported;
an optional bundle still has to be the actual elastic state. The synthetic
`tests/fatigue_source.rs` regressions cover static damage, native Newmark motion,
source refusals and complete checkpoint replay without changing mechanics.

`apply_load_jump` and `apply_support_acceleration_jump` expose the original
native finite-force and support-acceleration events at the accepted time.
They retain the native `LoadJump3` diagnostics and full mass/support coupling.
Displacement, velocity, temperature and stress stay fixed, so no extra stress
recovery, history copy, sample or cycle is needed. Open rainflow residues and
sample/step allowances remain unchanged. Native rejection or cancellation
publishes neither force nor acceleration; finalized observers refuse events.

The wrapper exposes no mutable run, external stress injection or material/thermal
reconfiguration. A late gauge refusal, stress-domain failure or cancellation
commits neither native mechanics nor any fatigue counter. Cloning retains both
accepted mechanics and unresolved cycles, with shared immutable matrices and
bounded state/residue copies. Read-only stress recovery adds no sample.
`finish` closes endpoint half cycles exactly once and prevents further steps;
processing chunks must not finalize their intermediate residues.

This is an elastic, fixed-temperature, fixed-plane empirical stress-life
estimate. Caller-supplied S-N/process/environment applicability is not inferred
from matching identity/temperature. Time and spatial refinement remain necessary
to resolve physical stress extrema. There is no multiaxial criterion, plastic
strain-life, changing-temperature law, remesh transfer, crack propagation,
stiffness degradation, damage-derived heat or service-life validation claim.

G0/G1/G3/G4 regressions in `tests/transient_fatigue.rs` cover unforced oscillator
stress and independently accumulated damage, a three-rung physical-period
reference, exact native mechanics/reactions, moving supports, checkpoint/chunk
continuation, finite load/support-acceleration events through unresolved
histories and exhausted sample allowances, complete rollback and source/history
refusals. All fatigue data
in those tests is explicitly synthetic; the regressions are implementation
checks, not experimental material calibration.

## Ordered Mode-I crack growth from native elastic cycles

`fatigue::crack_growth::TetCrackGrowthPath` recovers the original quasistatic
P1 stress, and `TetTransientRun3::with_crack_growth` attaches the same observer
to the actual initial and accepted sparse Newmark responses. The inertial owner
retains original mass, damping, thermal eigenstrain, prescribed-support motion,
reactions and work. Read-only response recovery adds no observation; there is no
mutable mechanical escape or external stress injection into either history.

Each gauge declares an uncracked nominal-stress location, fixed crack-plane
normal and either constant geometry or a bounded log-log geometry-factor table
with its admitted crack-length interval. Physical
tensor shear is converted by the existing tensor-basis owner. Normal traction
must be tensile; the norm of tangential plane traction must stay within the
explicit Mode-I numerical allowance at every sample. The current crack's Kmax
is checked on every sample, including static holds and unclosed excursions.
The fixed material identity and actual retained temperature must match the
Paris calibration scope. Caller-declared laws may carry an original material
bundle, which binds identity, chemistry, process and `T` without experimentally
validating those coefficients. A source-resolved Paris law requires the actual
immutable material bundle used by the native elastic owner. The shared static
and transient observer verifies the source card, named material state and
every fixed typed query coordinate. Only calibration R and squared Delta-K
coordinates are excluded; the runtime law enforces their source domains.
Custom source temperature axes remain intact. These checks retain original
calibration evidence without establishing validated component life.

An explicitly selected Walker fit derives its operating R from those same
native tensile stress extrema. Mean load and fixed thermal/residual prestress
therefore affect growth through the calibrated shared law while original
elastic stress, displacement, inertia, reactions and work remain unchanged.
The reference R and gamma remain inspectable in the retained law and, for a
source-resolved fit, in its original six material claims. This does not change
the fixed-temperature restriction or establish a changing-temperature transfer.

The inertial wrapper also exposes native `apply_load_jump` and
`apply_support_acceleration_jump` events. At the accepted time they update
force, acceleration and reactions through the original consistent-mass owner,
retaining its `LoadJump3` report. They preserve displacement, velocity,
temperature, stress, open excursions and closed but unconsumed cycles exactly.
There is no duplicate observation, observer clone or cycle closure, and no
step/sample allowance is consumed. Native failure or cancellation leaves both
owners unchanged; a finalized observer refuses further events.

The original immutable geometry and its context remain inspectable through
`geometry(gauge)`. The shared material owner evaluates the changing factor
during advancement and forecasting; Y is not frozen at the initial crack.
Constant geometry retains constant work, while tabulated evolution visits only
crossed intervals and checks all intermediate intensity extrema. Native
forecast cancellation propagates through these intervals. Table refinement
and geometric applicability remain caller obligations, and the original
mechanical stiffness, native stress and physical motion remain unchanged.

A cycle starts at a resolved physical extremum and contains one monotone
outward segment and one monotone return segment, with plateaus admitted.
Explicit closure advances exactly one full cycle through the shared analytical
`fs-material` Paris owner at the current accepted crack length. Peak-started
and valley-started cycles are both supported. Additional reversals, a return
past the starting extremum beyond the declared closure tolerance, or incomplete
closure refuse. The last actual endpoint anchors the next cycle; tolerances do
not snap stress, discard drift or declare an arbitrary waveform periodic.
No rainflow completion order is used as load order and no endpoint half cycle
or repeated spectrum is invented. Retained cycle storage is constant per gauge.

For an explicitly declared identical quasistatic loading block,
`TetCrackGrowthPath::close_repeated_cycle(full_cycles, &Cx)` closes one sampled
representative cycle and analytically integrates that many full cycles,
including the representative cycle. Work and retained storage per gauge are
independent of the repeat count; the original native solves are not repeated.
The entire block must remain inside the same material, stress-intensity,
geometry and exposure bounds. Counts must be whole cycles in `1..=2^53` for
exact conversion to the material owner's cycle coordinate. The caller declares
repetition under unchanged elastic, thermal and support fields; physical sample
time and sample count stay unchanged. No transient periodicity, thermal history,
plastic stabilization or crack-induced stiffness change is inferred. A refused
block preserves every crack and the representative window for a bounded retry.
The inertial wrapper intentionally has no repeated-cycle operation. The runnable
`fatigue_crack_growth --demo CYCLES` example uses three native solves and clearly
synthetic coefficients, and prints an independent analytical m=2 comparison.
Its optional `--varying-geometry` mode exercises a declared Y proportional to
the square root of crack length and compares the independent rational solution.
The `--finite-width` mode uses the bounded elastic centre-crack secant helper,
retaining full plate width, physical half-crack length and remote gross-stress
conventions. It requests a relative-Y interpolation tolerance and a point cap,
and reports the actual adaptive table resolution and direct secant factor.
Its material coefficients remain explicitly synthetic; the analytical
interpolation criterion is not a floating-point error certificate and does not
establish component or material validation.
The composable `--walker` flag selects explicitly synthetic reference-R/gamma
parameters and prints their independent rate multiplier. Constant-geometry,
varying-geometry and finite-width comparison references all include that
same multiplier; no hidden mean-stress coefficient is selected.

Both native owners expose a read-only `cycles_to_length` forecast for a gauge
whose representative cycle is fully sampled but not yet explicitly consumed.
It uses that original nominal stress range and the shared analytical inverse,
subject to every fit, geometry and remaining-exposure bound. The requested crack
length is a caller target, not an inferred failure length. Identical future
extrema are an explicit forecasting assumption; no future dynamics, component
life, cycle count or mechanical motion is silently accepted by this query.

An endpoint and requested closure publish together only after original
mechanics, every projection and every crack-growth update accept. A late gauge
failure or cancellation preserves native motion, physical time, all open
windows and crack lengths. Inertial checkpoint clones retain both histories.
A separate closure of an already accepted window changes growth only; a failed
closure leaves that accepted endpoint intact. Global gauge/sample bounds and
per-gauge cycle/block bounds remain in force. Finalization refuses incomplete
excursions and then prevents further samples. Static holds produce no crack
extension; positive ranges below the calibrated Paris regime refuse rather
than becoming an invented fatigue-arrest threshold.

This is fixed-temperature, small-strain, tensile Mode-I nominal-stress LEFM
postprocessing. Small-scale yielding, geometric applicability and adequate
spatial/time resolution of stress extrema remain explicit caller obligations.
The mesh has no advancing crack surface, and growth does not change stiffness,
thermal energy, constitutive state or topology. There is no mixed-mode growth,
plastic-zone correction, closure/retardation model, thermal history transfer,
crack initiation inference or validated component-life claim. Synthetic
regressions in `tests/crack_growth.rs` compare independent analytical growth,
actual inertial/support-driven stress and original mechanical responses, and
exercise chronological closure, finite drift, domain refusals, cancellation,
resource bounds and complete checkpoint rollback.

## Finite-strain traction-free plane stress

`HyperProblem::solve_plane_stress(PlaneStressSettings, &Cx)` runs the original
load-stepped Newton solver on P1 triangles and Q1 quadrilaterals while the native
`fs-material` energy relaxes all three transverse nominal tractions at each
quadrature point. Thickness contraction and anisotropic transverse shear are
solved, not prescribed. The in-plane force and Schur-complement tangent derive
from the same reduced energy; supports, dead tractions and reactions retain the
original reference-thickness convention. `plane_stress_states` exposes the
actual relaxed 3-D deformation and local equilibrium residual.

Reference gradients and quadrature are prepared once per solve. Energy-only
line searches use the material's three-direction local Newton path without
forming full constitutive Hessians. The global solver checks the true linear
residual, polls cancellation between quadrature points and bounded Krylov
chunks, and refuses nonconvergence, folded/inverted elements, invalid controls
and material-domain violations without returning a partial solution.

G0/G1/G3 tests in `tests/hyper_plane_stress.rs` cover an independently derived
finite Neo-Hookean uniaxial traction solve with lateral/thickness contraction,
nodal energy/force/tangent consistency on P1 and skew Q1 elements with tilted
fibres, rigid-motion invariance, and domain/cancellation refusals. This is a
thin in-plane continuum: no shell bending, wrinkling, global injectivity,
fully incompressible formulation, general stability or convergence certificate
is implied. Existing `solve()` and continuation retain plane-strain semantics.
`cargo run -p fs-solid --example hyperelastic_plane_stress` solves a thin coupon
under nominal traction and reports displacement, relaxed thickness and support
reaction. Its constitutive coefficients are explicitly authored model inputs.

## Directional small-strain forming through the native solid solvers

The existing generic `plastic2::SmallStrainPlaneProblem` and
`plastic3::SmallStrainTetProblem` consume `fs-material::plastic::hill::Hill48`
without replacing the element, quadrature, plane-stress reduction or Newton
owner. Each point retains its complete Hill plastic strain and work-conjugate
accumulation. Directional yielding therefore changes actual nodal deformation,
support reactions and solved transverse strain, including permanent width and
thickness changes after unloading. Fixed rotated material axes may be supplied
through `OrientedHill48`; this does not introduce finite-rotation kinematics.

G0/G1/G3 `tests/plastic_hill.rs` uses independent free-sided uniaxial flow
equations for P1/Q1 plane stress and native tetrahedra, both loading signs and
all three principal material directions. It also checks the condensed
algorithmic tangent against assembled-force differences on skew Q1 elements.
The original small-strain/rotation and locking limits remain binding. The
separate finite-plastic sheet owner retains its existing J2 history type;
these additions make no finite-strain anisotropic-sheet claim.

## Retained finite-plastic sheet forming

`hyper2d::plastic::FinitePlasticSheet` advances a P1/Q1 sheet through the
original three-dimensional finite J2, Perzyna overstress, Norton creep or
combined creep/J2 return map and shared plastic increment solver. The material
bound is `FinitePlaneStressMaterial`; kinetic owners are never reduced to
their rate-independent backbones. All three transverse nominal tractions relax at each point; the
Schur-complement material tangent enters the actual element tangent. The sheet
retains full accepted deformation, plastic history and original per-element
native/source-resolved law bindings. Geometry is prepared once, and every
weight is reference area times original element thickness. P1 has one exact
centroid evaluation; Q1 retains four independent Gauss-point histories.

Each endpoint starts from the accepted displacement and uses frozen history
through local/global Newton and backtracking. Free in-plane coordinates solve
equilibrium; recovered normal thickness is `h0 * F22`. Independent energy
checks distinguish elastic/hardening storage, native plastic/creep dissipation,
the actual reduced incremental potential, support/dead-load work and the signed
endpoint quadrature defect. Potential is checked against assembly independently;
it is not identified with storage plus physical loss for kinetic laws. The
inelastic-dissipation accessor sums the native irreversible accounts, while
signed work defects never become heat. An exhausted
solve, material-domain exit, folded element, energy refusal or cancellation
returns no partially updated sheet. `SheetAssemblyBudget` bounds topology and
free coordinates; cancellation is polled between material points. Assembly is
sequential and deterministic. The shared increment solver's component count
defaults to three; sheets select two solely for translation-free force scaling.

Every kinetic element must declare the same positive physical duration, held
fixed through local/global trials and endpoint audit. Static elements can share
that interval. `physical_time_s()` is `None` for a wholly static sheet and starts
at `Some(0.0)` for a kinetic sheet; it advances only with an accepted endpoint.
Repeated fixed-deformation or fixed-load endpoints therefore resolve physical
relaxation or creep. Accepted-state inspection retains legitimate viscous
overstress instead of imposing equilibrium yield or integrating a second time.
Construction refuses any inelastic initialization flow, and incompatible
element durations or an unrepresentable next clock value refuse transactionally.

G0/G1/G3/G4 regressions in `tests/finite_plastic_sheet.rs` exercise independent
Hencky/J2 uniaxial tension/compression formulas, free lateral/interior motion,
unloading with permanent shape and thickness, hardening versus heat, assembled
force/tangent/energy derivatives, heterogeneous original laws and thicknesses,
and transactional refusal/cancellation. The rate-sheet regressions in
`tests/finite_plastic_sheet_rate.rs` add kinetic holds, creep under free-coordinate
equilibrium, distinct physical/kinetic accounts, and retained clock/history.
Coefficients are synthetic. This is an isothermal in-plane continuum; it does not claim shell
bending, wrinkling, necking, contact, anisotropic yielding, fracture, calibrated
forming limits or a general mesh-convergence result.

## Thermally evolving finite-plastic sheets

`hyper2d::plastic::thermal::FiniteThermoplasticSheet` retains the original
borrowed `FiniteThermalJ2Material` families and uses the same prepared P1/Q1
geometry, condensed plane-stress assembly and increment solver as the
fixed-temperature sheet. `try_new` admits virgin identity deformation only at
each family's own reference temperature, using duration-free inspection.
Every `advance` receives an exact element-temperature field and positive
physical duration. Native thermal expansion, coefficient softening and selected
rate-independent, Perzyna, Norton or simultaneous Norton/J2 evolution remain
inside the original family; no mechanical backbone or averaged property field
replaces it. Old history is inspected at its accepted old temperature, so a
new hotter endpoint can legitimately cause flow while displacement is held.

Each element freezes its new temperature and actual duration through all
Newton/backtracking trials. The accepted temperature, time, full deformation,
traction-free thickness and inelastic history advance atomically. The endpoint
audit calls `checked_increment_for_step` from the original accepted deformation,
temperature and history, independently matching the solver's candidate history
and assembled potential. `evaluate_trial` accepts endpoint temperatures first,
then displacements, duration and tangent selection, without committing them.

`quadrature_energy` reuses the established original-volume
`FiniteThermoplasticElementEnergy` accounts in element/quadrature order;
`quadrature_volume_m3` supplies their original integration weights. The whole
sheet energy account retains free and mechanical internal-energy changes,
plastic/creep dissipation, reversible caloric transfer, coefficient and entropy
storage changes, actual kinetic potential, support/dead-load work and signed
endpoint defect separately. Total caloric transfer is native irreversible loss
plus reversible exchange; the endpoint quadrature defect never becomes heat.
Point balances, retained energy, assembled potential, external/internal work
and the complete thermal balance are gated before acceptance. Temperature is
caller-supplied; this is not a thermal-equation or caloric-material owner.

`tests/finite_thermoplastic_sheet.rs` covers free thermal expansion/cooling,
held-strain thermal softening, explicit-duration viscous relaxation, heterogeneous
original families and temperatures, native point/whole-sheet work accounts,
and transactional domain, duration and cancellation refusals. Coefficients
are synthetic and carry the same in-plane continuum limitations as the
fixed-temperature sheet.

### Temperatures at native sheet quadrature points

The element-temperature API remains available. The explicit
`try_new_with_quadrature_temperatures`,
`advance_with_quadrature_temperatures` and `evaluate_trial_for_step` paths
retain an independent temperature at every P1 centroid or Q1 Gauss point.
Every constitutive trial and accepted two-temperature audit uses that point's
own original family, temperature, deformation and plastic/creep history.
No element average enters the mechanical response. Equal temperatures within
an element reuse its frozen material value; the common prepared geometry,
assembly, condensed tangent and Newton solver are shared with static sheets.

`quadrature_temperatures_k` exposes the actual accepted point temperatures.
The established `temperatures_k` accessor retains one value per element;
after a pointwise endpoint it reports an original-volume weighted diagnostic
mean. The physical clock advances only after the complete endpoint is accepted.
The accepted response, physical plastic/creep heat, native incremental
potential, assembled work, free/internal storage and reversible heat are
audited separately before any state is published.

`tests/finite_thermoplastic_sheet_points.rs` adds independent scalar
Hencky/Voce/Perzyna and simultaneous Norton/J2 oracles for free-sided forming,
relaxation and heated dead-load creep. These tests solve free lateral and
interior displacements and compare full deformation, thickness, retained
history and thermal work accounts; their coefficients are synthetic.

### Conservative nodal temperature and caloric transfer

`quadrature_temperatures_from_nodes` interpolates one positive absolute
temperature per original node with the native P1/Q1 shape functions.
`advance_from_nodal_temperatures` feeds that point field into the same
physical-duration, frozen-history endpoint solve. Cached shapes are admitted
once with the immutable reference geometry; constant fields are preserved
exactly, and no node-to-element average replaces a point's temperature.

An accepted step's `nodal_caloric_heat_j` returns signed joules on the original
nodes through the transpose transfer, `Q_a = sum_q N_a(q) V0(q) h_q`.
The density `h_q` is native plastic/creep loss plus reversible caloric exchange.
Original quadrature volume includes the original thickness. Reversible cooling
remains negative; current volume, work-quadrature defects and kinetic potentials
never replace physical heat. The transpose conserves the integrated signed
energy and its pairing with every representable nodal temperature field.
The receiving thermal owner applies the returned energy once and retains its
own caloric storage and transport law.

`tests/finite_thermoplastic_sheet_transfer.rs` checks exact constant and affine
temperature reproduction on translated/sheared P1/Q1 meshes, independent signed
heat conservation and polynomial duality through cooling, forming and spatial
warming, and transactional input/cancellation refusals.

## Exactly incompressible finite-strain plane stress

`HyperProblem::solve_incompressible_plane_stress(PlaneStressSettings, &Cx)`
uses the same P1/Q1 geometry, integration, Newton solver and boundary loads as
the finite-bulk plane-stress path. At each point the native material reduction
enforces `det F = 1`, relaxes transverse shear and recovers the traction-free
pressure multiplier. Its reduced energy, forces and tangent do not use a bulk
penalty. The stored constitutive bulk coefficient can span finite scales
without changing this explicitly selected incompressible problem.

`incompressible_plane_stress_states` recovers physical first Piola stress,
volume-preserving 3-D deformation and the original material remainder's
pressure multiplier. This multiplier's gauge is not mean Cauchy pressure.
The corresponding potential and residual/tangent APIs retain reference-thickness
work conjugacy and reaction recovery. Original material stretch bounds apply to
the actual relaxed thickness; an allowed in-plane stretch does not override them.

G0/G1 tests in `tests/hyper_plane_stress.rs` independently recover incompressible
Neo-Hookean uniaxial contraction, stress, energy and support reaction on triangles
and quads; show bulk-coefficient independence through `1e300`; differentiate
assembled reduced energy/forces for tilted fibres on skew elements; and check
thickness-domain/cancellation refusals. The existing example accepts
`--incompressible`. This is the exact thin-continuum volume constraint, not a
general 3-D mixed formulation, shell bending model or arbitrary-mesh stability
certificate.

## Finite-strain Q1 mean dilatation

`HyperProblem::solve_mean_dilatation(&Cx)` adds an opt-in plane-strain Q1
formulation for quasi-incompressible compression. Its condensed element energy
integrates the original material remainder at all four Gauss points and the
original bulk energy at reference-volume-mean det F. Eight-DOF nested automatic
differentiation supplies the complete residual/Hessian, including geometric and
cross-quadrature terms. NH, MR, Ogden and full-invariant HGO retain their original
material split; no alternate constitutive coefficients are substituted.

Reference geometry is prepared once. The original load stepping/Newton owner
is reused, with bounded cancellation and an independently recomputed linear
residual. `mean_dilatation_states` distinguishes local det F, mean det F and
tension-positive condensed bulk Cauchy stress from the full constitutive stress.
Reference/current corners and every local material state must remain admissible;
an admissible mean cannot hide local inversion or a stretch-domain violation.

Focused G0/G1 regressions in `tests/mean_dilatation.rs` cover energy derivatives,
analytic finite compression/reactions on two meshes, an independent bulk-stiffness
identity, a nonuniform compressed-block reference, objectivity and refusal.
This mitigates the finite bulk penalty's excess local constraints. It supplies
no 3D/P1/TDNNS, fully incompressible, hourglass-stability, continuum-convergence
or contact guarantee. Existing `solve()` and continuation retain standard behavior.

## Static compression by spatial follower pressure

`HyperTetProblem::solve_with_face_pressures` consumes the existing canonical
face-pressure field in native finite-deformation equilibrium. Explicit zero
values select unloaded faces; positive values act inward on CURRENT area.
Loads, prescribed displacements and pressures share the original load ramp,
constitutive law, Newton assembly and support-reaction owner. The original
scalar uniform-pressure entry point retains its conservative solver behavior.

A partially loaded boundary generally has a nonsymmetric follower tangent and
path-dependent work. The new path retains the full analytic pressure Jacobian,
uses cancellable FGMRES with ILU(0) or absolute-diagonal preconditioning, and
globalizes against the force residual. It does not manufacture a pressure-volume
potential. The explicit linear iteration budget and independently recomputed
linear residual apply. All-zero face fields use the dead-load path exactly.

G1 tests in `tests/hyper3.rs` reproduce independent Neo-Hookean uniaxial
compression on two meshes, including lateral expansion, current-area reaction,
volume and stored energy. `tests/hyper3_pressure.rs` exercises a compressed,
sheared equilibrium with an actually nonsymmetric free tangent, reaction loads,
uniform-field equivalence, exact zero-field identity, invalid face coverage,
nonfinite values, cancellation and an exhausted linear budget. The existing
P1 volumetric-locking and stability limitations remain. Face loads are supplied
data, not a fluid equation of state or contact/platen law; no mesh-convergence
claim follows from an exactly representable affine test.

## Spatial follower pressure on finite-deformation exterior faces

`hyper3::pressure::TetPressure::PerFace` supplies finite signed pressures in the
canonical order of `ClosedTetPressureSurface::faces`. Explicit zeros select
unloaded faces. The original prepared geometry and force/tangent kernel serve
both a scalar uniform pressure and a face field. Positive pressure acts inward
on current face area; a patch keeps its own resultant, torque and generally
nonsymmetric external-force Jacobian. It is not assigned a volume potential.

Held face pressures have exact work along a straight nodal-displacement segment
to rounding, using Simpson integration of the quadratic P1 area vector. Common
translation work on an open patch remains present. Curved physical trajectories
require time refinement; this is not path independence or an injectivity proof.
Uniform pressure retains its original exact endpoint-volume work and static
solver path. Missing/extra/nonfinite face values and cancellation refuse before
returning a load. G0/G1/G3/G4 tests check localized outward orientation, independent
finite-difference tangents, path-dependent analytic work, uniform compatibility
and refusal. Dynamic consumption belongs to `fs-couple::finite_solid`.

### Hydrostatic and affine Eulerian pressure

`TetPressure::Affine` prescribes a stationary world-space field
`p(x) = p_origin + gradient dot (x - origin)` on the complete exterior. A uniform
fluid in uniform gravity uses `gradient = rho_fluid * gravity`; signed pressure
is differential traction data. The original load owner samples CURRENT node
positions and integrates `N_i p(x)` exactly on each P1 triangle. The resulting
nodal force uses weights `(2 p_i + p_j + p_k)/12` multiplying the inward area
vector, retaining both force and torque. The analytic Jacobian differentiates
the current area and spatial pressure, and is symmetric on a closed surface.
Its net force is `-current_volume * gradient`, so physical buoyancy follows
directly from the surface traction.

The potential is the volume integral of `p(x)`. Simpson integration of the
cubic force along the declared straight nodal segment gives its exact work to
rounding, including rigid translation through the field. A zero gradient uses
the existing uniform-pressure arithmetic. `solve_with_pressure_field` scales
both pressure and gradient with the native load factor and includes the same
conservative work in Newton globalization and support reactions. The existing
gravity-plus-pressure and ordinary finite-solid dynamic owners consume the
same field. The specialized energy/momentum dynamic pressure algorithm remains
limited to uniform pressure.

`tests/hyper3_pressure.rs` compares forces and work against independently
differentiated tetrahedral volume moments; checks buoyancy, torque, tangent
symmetry and finite differences, path independence, reference-anchor covariance,
uniform compatibility and refusal; and verifies a manufactured Neo-Hookean
equilibrium and support reactions using independent Piola stress. This is a
fully immersed affine field: free-surface clipping, varying fluid density,
cavitation, fluid evolution and global self-intersection remain outside this
load owner's contract.

### Partial immersion at a prescribed planar free surface

`TetPressure::FreeSurface` applies gauge pressure
`max(0, p_origin + gradient dot (current_x - origin))`; uniform-fluid
hydrostatics uses `gradient = rho_fluid * gravity`. Wet polygons are split into
at most two triangles in the original P1 barycentric chart. Exact integrals of
`N_i p` retain the nodal force and torque, while wet `N_i N_j` moments supply the
spatial-pressure tangent alongside the current-area derivative. The moving-cut
term vanishes at zero pressure. A whole face on the plane selects the dry-side
generalized tangent; no unique two-sided tangent is claimed at that state.

Conservative work uses exact endpoint moments of `integral_volume max(0,p) dV`.
Positive decompositions into at most three wet tetrahedra avoid subtracting an
almost equal dry volume or dividing by a small pressure gradient. Element work
increments retain changes in immersion without treating a crossing as a fixed
Simpson interval. The native Newton load ramp scales pressure and gradient
together, preserving the zero plane, and uses this work in globalization and
the same traction in reactions. Fully wet and zero-gradient limits retain the
existing affine or nonnegative uniform load; a currently dry body can enter
the fluid during a solve.

`tests/hyper3_pressure/free_surface.rs` checks independent truncated-power
volume moments, wet-volume buoyancy and torque, consistent tangents, thin caps,
transition work, extreme finite pressure ratios, refusal and a native equilibrium
that crosses a wet vertex. This is a prescribed exposure field on the prepared
boundary; fluid/free-surface evolution, capillarity and cavitation are not added.

## Finite-deformation static follower pressure

`hyper3::HyperTetProblem::solve_with_pressure` couples dead referential loads
and prescribed displacement with uniform follower pressure on the exact matching
complete closed boundary. The existing Newton owner evaluates current-area
forces, the analytic pressure tangent and conservative pressure-volume work at
every accepted or trial state. Force and volume share one geometry evaluation;
Armijo compares energy/work increments. Pressure ramps with load factor, and
support reactions include current pressure and applied dead loads. Positive
pressure compresses; negative pressure denotes outward traction. Invalid trial
geometry backtracks, cancellation aborts, and nonfinite reactions refuse.

This does not add an equation of state, open-patch loading, nonuniform pressure,
mixed incompressible elements or limit-point continuation. The focused static
tests in `tests/hyper3_pressure.rs` compare finite hydrostatic compression and
expansion, frame transformations, fixed-area traction differences, zero-pressure
equivalence and finite-deformation support reactions with independent formulas.
These synthetic regressions require native execution; they do not establish
experimental material calibration or mesh convergence.

## Conservative gravity on mixed P2/P1 solids

`MixedUpProblem::prepare_with_gravity` binds one positive reference density per
Tet10 element and the original stationary `GravityField3`, optionally composed
with uniform pressure on the complete closed material boundary. The native
mixed preparation and solver retain their pressure space, volumetric treatment,
pressure gauges, source domains and stable-KKT admission. Gravity acts only on
displacement; it adds no pressure mass or bulk penalty.

Current sampling uses the original cached quadratic shapes and positive
reference quadrature, with translation-stable relative coordinates. Force
and its spatial derivative enter the existing 34-by-34 element block and sparse
scatter directly. The same density, shapes and weights supply consistent modal
mass, including the signed quadratic shape contributions. Current det(F)
never multiplies the fixed material mass. General nonlinear fields and curved
geometry retain the declared reference quadrature approximation.

The wrapper exposes evaluation and native static solves, plus finite-compliance
and exact-incompressible vibration, each with optional accepted initial iterates.
Both modal paths consume the SAME loaded equilibrium and tangent; original
pressure condensation or exact-constraint projection and full residual gates
remain intact. Results retain the exact density/field through `HeldMixedGravity`.
Static results report gravitational and boundary-pressure potentials separately
from the original material saddle potential.

Field-domain and gradient errors retain element/point attribution. Callback
cancellation remains canonical and sticky within the native solve, participates
in its existing stop callback and prevents further field sampling. Refusals
publish no altered input state. This is external gravity for equilibrium and
small vibration; self-gravity, changing sources and gravity-loaded motion are
not introduced.

Native `mixed_gravity` regressions include independent rational Tet10 moments,
force/potential derivatives, harmonic-field frequency shifts, exact-mode bulk
independence, pressure gauges, combined-pressure reactions and source/cancellation
refusals. Bounded original-source kernel checks cover sampling and element
assembly; the complete static/modal integration tests remain unexecuted until
the workspace sibling dependencies are available.

## Native three-dimensional hyperelastic continuation

`HyperTetProblem::prepare_continuation` binds the original prepared P1
tetrahedra, hyperelastic material adapters, homogeneous supports and dead
nodal/reference-volume loads to the existing bordered pseudo-arclength owner.
`prepare_continuation_with_loads` additionally retains native current-area
pressure and stationary spatial gravity. The pressure surface must match the
exact reference geometry; original elementwise or P1 nodal density is bound
once into fixed reference gravity weights. No current-volume mass update or
cell-average replacement occurs.

Explicit positive displacement and force scales define q=u/L and
r=(f_int-lambda*f_ext(u))/F. Both the complete displacement tangent
(K_int-lambda*K_ext)*L/F and the current parameter derivative f_ext(u)/F
enter every bordered Newton linearization. Material, pressure and gravity
rows are merged directly in canonical sparse order without another global
COO buffer, sort or dense intermediate. The parameter scales the held force;
spatial origins and free-surface planes remain fixed even for negative lambda.
Nonuniform face patches retain their original nonsymmetric tangent.

Fixed coordinates are eliminated from the path metric, while all supported
dead, pressure and gravity loads remain in recovered reactions. `held_loads`
retains the exact original bindings, including explicit zero fields.
`start_from` admits an actual physical equilibrium against the same force
equation. `advance` returns a complete candidate, preserving its input even
when a later requested step fails or is cancelled. `response` recovers metres,
newtons, original material stored energy, det F and free-force residuals.

The public `PathState` retains dimensionless free coordinates and complete
trace, tangent and step control for split/resume. `PathResidual::load_vector_at`
admits state-dependent parameter derivatives while its default preserves
constant dead loads. Checkpoints run at corrector and controlled FGMRES
boundaries; dense factorizations remain bounded leaves. Original
element-attributed material errors survive the generic boundary. The existing
Cholesky branch diagnostic skips nonsymmetric or oversized tangents; it does
not interpret a skipped measurement as a change of definiteness.

Native `hyper3_continuation` tests include independent finite-strain force,
energy/reactions, scaling, pressure/gravity derivatives, hydrostatic load
turning points, field-domain refusals and atomic cancellation. These folds
verify the declared constitutive/loading equations; they do not validate
calibrated rubber behavior or geometric buckling. Generic original-source
tests additionally exercise nonlinear parameter amplitudes, follower folds,
bit-identical dead-load compatibility and sparse cancellation. Full native
assembly/API tests require the complete workspace. Moving supports, evolving
temperature/history, mixed incompressibility, branch uniqueness and continuum
stability or flutter certificates remain outside this path interface.

## Prepared finite-deformation reference geometry

`hyper3::prepared::PreparedHyperTetProblem` retains admitted reference shape
gradients and volumes while immutably borrowing the original mesh and material
field. Repeated force, energy and tangent evaluations use the same native
constitutive and assembly kernels; current displacement, positive Jacobian,
material-domain and finite-response checks remain active. New interval dead
loads are validated on every call and use the original reference volumes.
Preparation takes O(elements) storage and removes repeated reference geometry
construction; it does not change the P1 discretization or mitigate volumetric
locking. The finite-solid dynamic adapter shares this preparation across
candidate steps. G0/G1 checks cover changing finite stretches against closed-form
Neo-Hookean energy/force/tangent, exact native assembly, new load admission,
inversion, invalid dimensions, budgets and cancellation.

`stored_energy_interval_j` streams the original heterogeneous material field
and cached reference volumes once, admitting both total deformations and each
original endpoint stress/potential. It accumulates endpoint storage and the
source-owned potential increment separately with compensated products and
sums, without allocating a force vector or tangent. `stored_energy_change_j`
exposes that direct change. Custom and frozen inelastic adapters retain their
own potential semantics and default endpoint difference unless they supply an
increment capability; a generic stable-increment claim is not inferred.

### Objective energy/momentum force on the prepared mesh

`PreparedHyperTetProblem::energy_momentum_internal_force` and
`energy_momentum_tangent` evaluate an objective strain-space discrete gradient
between two configurations on the original P1 mesh. The mean right
Cauchy--Green metric is `Cbar=(C0+C1)/2`, represented by
`F*=chol(Cbar)^T`; it is not the metric of the mean deformation gradient.
The original material supplies its direct energy increment, Piola stress, tangent and endpoint-domain
admission. A symmetric second-Piola work correction, multiplied by midpoint F,
retains partition-of-unity force balance, midpoint internal torque balance and
the potential-work identity. No physical damping or constitutive replacement
is introduced.

The analytic generally nonsymmetric 9-by-9 endpoint tangent differentiates the
same Cholesky representative, material pullback and work quotient. A consumer
using midpoint displacement must multiply its action by two. Near
`||Delta C||^2 <= 64 epsilon ||Cbar||^2`, the stress takes the continuous limit;
measured nonlinear and energy tolerances remain authoritative. Original
reference geometry, material/reference adapters and consistent mass are retained.

Endpoint metric increments, the original potential-work correction and the
force/tangent pullback retain compensated arithmetic through final rounding.
The stiff-bulk regression resolves a 6e-11 Pa isochoric work increment and its
375 Pa shear tangent beside 1e12 Pa bulk stiffness, while preserving endpoint
reversal and midpoint torque. The original continuous-limit cutoff and
material/source domain gates remain unchanged.

Four exact native stress/tangent tests pass, covering work, internal force and
torque, all nine tangent directions, objectivity, endpoint reversal, the
continuous limit and source-domain refusal. Two 250-step native free-tumbling
P1 tetrahedron trajectories with the original generalized-alpha/Newton-FGMRES
owner retain energy and both momenta near roundoff for the admitted fixtures.
Full prepared/Cx/material-card integration regressions are authored but remain
unexecuted because the required sibling runtime dependencies are unavailable.

### Uniform closed pressure in the objective discrete gradient

`energy_momentum_force_with_uniform_pressure` and its endpoint tangent combine
the original frozen material with the pressure potential `p det(F)` before
the same objective strain-space algorithm. Positive pressure is inward. The
native prepared closed-pressure surface owns reference matching and finite-load
admission; its original endpoint-volume work remains independently measurable.
The analytic cofactor stress and Hessian preserve the original material/source
domain. Zero pressure delegates exactly to the unpressurized algorithm.

The work target combines the original material increment with `p Delta J`
before rounding, rather than subtracting two absolute `W+pJ` values. This
retains acoustic-scale elastic work at equal volume even under much larger
ambient pressure. The shared cubic determinant increment evaluates
`cof(M):D + det(D)/4` with compensated endpoint means and differences, retaining
the tails lost by a rounded midpoint. Its endpoint derivative remains the
original material stress plus `p cof(F_after)`; the pressure-coefficient force
uses the same volume increment. The material's separate increment capability
retains its source and reference adapter; the pressure decomposition does not
override the material's potential or its cancellation guarantees.

Uniform closed-surface pressure work accumulates oriented element volume
increments before applying pressure. Affine, free-surface and stratified work
use the same increment in the symmetric volume/pressure product difference.
Held per-face loading retains its original path-dependent force quadrature.
The geometry path normalizes physical edge matrices before the cubic operation
and retains original endpoint volume, orientation and source-domain admission.
It preserves small work beside a large unrelated mesh volume without changing
surface traction or pressure fields. Eight original stress tests and three
volume-kernel tests pass, including the previously failing 100 kPa small-shear
case, rounded-midpoint tails, scales, reversal and consistent tangents. Three
pressure-owner regressions also pass through the original geometry, clipping
and stratified algorithms with a non-cancelled context shell; full executor/Cx
and Cargo qualification remain outside that focused verification.

The symmetric algorithmic second-Piola stress closes pressure work and midpoint
torque together. Exact straight-segment averaging of face forces closes work
alone and can produce a midpoint torque defect; a nontrivial regression exposes
that difference. This is a finite-step algorithmic volume load. It approaches
native boundary traction as endpoints coalesce, but its finite-step nodal field
need not equal midpoint traction or cancel on interior nodes. Nonuniform face
pressure remains outside this energy/momentum entry.

All six production stress/tangent groups pass, including signed pressure,
endpoint tangent directions, objectivity/reversal, native face-traction limits
and the work-versus-torque counterexample. Native time integration also closes
energy and both momenta and shows resolved second-order refinement. Public
prepared/Cx integration is additionally qualified by the native spatial-pressure
tests below; full Cargo qualification remains blocked by missing siblings.

### Spatial hydrostatics in the objective discrete gradient

`energy_momentum_force_with_pressure` and its endpoint tangent admit the native
`Uniform`, `Affine` and `FreeSurface` fields. Each tetrahedron retains the exact
pressure potential `V * m(p_nodes)`, where `m` is the original affine or positive
clipped-volume mean. The product discrete gradient is
`m_bar * DG(V) + V_bar * DG(m)`. The volume term reuses the original objective
strain-space kernel; the mean term averages coordinate increments over all node
orders. Fully wet and dry elements take their affine and zero limits. Existing
uniform-pressure and zero-pressure arithmetic remain unchanged.

Field forces lie along the pressure-gradient axis. A free body therefore keeps
linear momentum transverse to that axis and angular momentum about it, while
buoyancy changes the momentum along the axis. No conservation of all momenta
is claimed for an external spatial field. The analytic endpoint tangent includes
both changing volume and spatial pressure, using native scalar AD of the same
clipped moment. Midpoint unknowns still require the original factor of two.

Rational divided differences retain coincident-pressure limits and thin wet
caps without subtracting endpoint potentials. An edge landing on another
coplanar vertex uses its finite complete-moment quotient; only a numerator
unresolved at `64 epsilon` relative to its endpoint moments takes the analytic
midpoint limit. Entirely coplanar mixed states select the dry-side generalized
derivative. Native work and nonlinear gates remain authoritative. These are
finite-step conservative algorithmic loads; native surface traction is recovered
as endpoints coalesce and the trajectory still needs time refinement.

Three native scalar tests and four complete prepared/Cx tests pass, including
independent endpoint work, every assembled tangent entry, admissible force and
torque, reversal, thin caps, coplanar arrival, uniform bits and cancellation.
The existing 22 pressure/surface/static-Newton regressions also pass. These
focused executions use original source and numerical dependencies; the full
workspace Cargo suite remains blocked by missing sibling repositories. This
entry supplies no fluid evolution, evolving free surface or conservative
potential for `PerFace` loads.

## Source-pack mechanical coupons and nonthermal dead loads

`linear3::eigenstrain` admits explicit sorted unique dead nodal forces in
assembly, equilibrium solve and input-bound stress recovery. It shares the
existing force validator and physical solver, while retaining generic initial
strain and requiring no fictitious temperature. Fixed-DOF loads remain in the
full force vector for support reactions. An empty force list preserves the
original eigenstrain-only identity and response exactly; changed loads refuse
recovery against an older solution.

`examples/material_coupon.rs` reads a normalized material pack and an explicit
bounded TSV scenario, then runs native tetrahedral face-load integration,
equilibrium and stress recovery. It resolves independent isotropic parameter
pairs or complete source-bound FPL orthotropy with explicit L/R/T orientation.
Source property roles, quantity-typed coordinates, source/derivation/orientation
identities and strain ceilings remain explicit. Templates exercise aluminum
5005, C52100 bronze and loblolly pine without inline material presets.

The aligned uniform traction fixture allows Poisson contraction and has an
independent axial compliance reference. Output reports source-bound mass,
displacement, applied/support forces, strain energy and residuals. Its
load-displacement product is F dot u, not a loading-history work integral.
Rotated general anisotropy is a coarse FE result. Neither a source data lookup
nor this linear model supplies strength, creep, fracture or mesh-error bounds.

G0/G1 regressions cover combined initial strain and dead force, exact affine
displacement/stress/energy, fixed-support reactions, empty-load equivalence,
malformed forces, stale-load refusal, signed coupon loads and strain admission.

## Exact partially submerged reference pressure loads

`loads3::fields::SurfaceField3::FreeSurfacePressure` integrates the positive
part of an authored P1 gauge-pressure field over each selected exterior face.
The caller supplies the original signed coefficients; negative values locate
the dry portion. The existing nonlinear pressure owner's exact clipped-triangle
kernel supplies consistent nodal forces without mesh refinement, face sampling,
or additional per-face allocations. Clipping nodal samples before interpolation
would load dry area and is not equivalent.

Uniform-liquid hydrostatics uses `p(x) = rho * gravity dot (x - surface_point)`.
Atmospheric pressure is already subtracted. The resulting case feeds the native
linear static solver and the existing nodal-load history path, preserving total
force, first moments and work against P1 displacement. The load remains fixed
on reference geometry; moving immersion, follower stiffness, fluid evolution,
capillarity and cavitation require their separate owners.

G0/G1/G3 regressions in `tests/load_fields3.rs` check exact wet-triangle nodal
work, a partially submerged wall's force and center of pressure on three meshes,
Archimedes buoyancy and torque, native support reactions, wet/dry limits, and
invalid coefficients and face assignments. These are discrete manufactured
hydrostatic controls, not free-surface-flow validation.

`SurfaceField3::StratifiedPressure` consumes the original prepared
`StratifiedHydrostaticPressure` directly on selected reference faces. The same
clipped-simplex integration used by follower pressure retains density interfaces,
density decreases and a continuous density slope before integrating each P1
nodal force. Three sampled pressures are not substituted for the reservoir:
that substitution changes hydrostatic force and center of pressure on a face
crossing a layer or free surface, and linearizes a quadratic pressure law.
No face allocations or alternate hydrostatic law are added. Nonnegative load
scaling acts after clipping; the original layer budget and cancellation polls
remain active. The assembled force stays a reference dead load and feeds the
existing static, reaction and nodal-history consumers without claiming evolving
immersion or fluid inertia.

Three additional G0/G1/G3 regressions in `tests/load_fields3.rs` cover exact
layered-wall force and first moment on three meshes, displaced-density buoyancy,
native support reactions, agreement with the follower owner at zero deformation,
the single-density limit, load scaling and invalid face assignments. The wall
case explicitly detects the error from replacing the reservoir by nodal samples.

## Material-bound spatial gravity loads

`loads3::history::VolumeHistory3::FixedAcceleration` accepts uniform, nodal P1,
or element-local P1 accelerations through `fields::AccelerationField3`. Values
are world-frame m/s² on the retained reference mesh. The borrowed solve problem
supplies the authoritative element density. Full consistent P1 mass acts before
support elimination, preserving interface density jumps, restrained-node loads,
resultants and first moments. The same scalar mass kernel supplies the original
elastic inertia. Typed uniform acceleration histories use this same integration
owner without allocating a second per-element force-density field.

The accepted sample enters the existing static, stress/reaction, support-history,
and sparse transient owners. G0/G1/G3 tests in `load_history3` compare against
full M*g, heterogeneous material mass, an analytic gravitational-potential
gradient, rotated reactions, initial dynamic acceleration and refusal behavior.
This is prescribed reference-sampled loading, not self-gravity, follower gravity,
or an accelerating-coordinate model. An inverse-square field sampled at vertices
retains P1 spatial approximation error and requires explicit resampling if its
physical owner changes position or time.

`TetLoadGeometry3::assemble_acceleration_with_nodal_reference_density` and
`TetLoadHistory3::at_with_nodal_reference_density` explicitly bind a positive
P1 reference-density field instead. They use the same exact cubic mass kernel
as graded elastic inertia, preserving density-weighted gravity moments and
the product of spatially varying density and acceleration. Full supported-node
forces remain available to the original static stress/reaction owner; surface
and point forces are superposed once. Uniform fields retain the original bits.
This entry projects a held force for static use; it does not replace the mass
in the separate `loads3` transient constructors. Three focused G1/G3 regressions
in `load_history3` cover independent barycentric integration, graded `M*a`,
weight/torque/reactions, scalar parity and refusals.

## Moving force/couple ports on the native shell trajectory

`moving_load::MovingShellModel` binds the original `shell::ShellPlate` geometry,
six-DOF spatial operators and fixed support membership once. A point or explicit
patch quadrature transfers Cartesian forces and free couples with the transpose
of its P1 material-point motion, including `offset cross force`. Full restrained
rows survive for reactions. The same rows sample displacement, director rotation
and material velocity; moving a load's attachment across triangles does not add
the convective traversal velocity to that material velocity or its power.

`MovingShellMotion` initializes acceleration from the original mass equilibrium
and reuses `shell_time::PreparedNewmark` with beta=1/4 and gamma=1/2. Retained old
and supplied new endpoint loads provide trapezoidal force and reference-moment
impulses. Support and damping exchanges, supplied work, and measured momentum
and energy residuals stay separate. Mass-proportional Rayleigh damping is a
momentum sink; it is not silently attributed to internal material heat. Force
projection is linear in point count, factors are reused, reaction evaluation
visits only restrained rows, and impulse reductions need no extra dense vectors.
The caller's point/step/horizon bounds and cancellation checks precede atomic
publication; duplicate, stale or future endpoint timestamps refuse.

This is the existing flat, small-rotation, estimate-only shell surrogate. Its
angular inventory and wrench moments use original reference positions and the
linearized director field. It is not a finite-rotation continuum shell,
follower-force/contact solver, automatic patch tracker, changing support topology,
or modal reduction. Time-varying load paths require temporal refinement. The
physical subset of bead `frankensim-b8bxd.9.3` is implemented; that broader
adapter's remaining features are not claimed complete. G0/G1/G3/G4 tests in
`tests/moving_load.rs` exercise independent wrench resultants, rigid motion,
reciprocal work, analytic accelerated translation, native damped/undamped moving
loads, support reactions, common-edge traversal, and atomic refusal/retry.

### Driven pin motion and actuator momentum

`start_driven` accepts the Cartesian displacement, velocity and acceleration of
all three original pinned nodes, in their original support-map order.
`advance_driven` accepts only new endpoint accelerations and advances the pin
velocities and displacements with the same average-acceleration kinematics as
the free shell. Prescribed motion enters the existing free solve as
`f_f - M_fc a_c - C_fc v_c - K_fc q_c`; the effective matrix and its prepared
factorization remain unchanged. Condensation visits the nine prescribed
translation columns, and all director rotations remain original free DOFs.

Reactions use the original full operators and retain prescribed inertial,
viscous and elastic terms. Material observations include the prescribed node
motion. Mechanical storage and linear/reference-angular momentum include both
free and prescribed DOFs; the vectors in `state()` remain free-DOF vectors,
while its energy fields describe the full shell. Signed `support_work_j` is
actuator work supplied to the shell, separate from applied-load work and
Rayleigh damping work. The reported residual is
`delta_energy + damping_work - applied_work - support_work`. The stationary API
retains zero support work, and identically zero driven inputs reproduce it.

A failed endpoint, nonfinite support input, or cancellation publishes neither
free nor prescribed candidate state. Support identities stay fixed; no finite
support-body dynamics, sliding constraints, or moving-contact law is implied.
Nonuniform or time-varying prescribed motion remains subject to the existing
small-rotation surrogate and requires temporal refinement. Focused native
tests in `tests/moving_load/driven.rs` cover analytic rigid acceleration and
actuator work, manufactured full-operator deformation with and without damping,
shuffled support order, changing accelerations, stationary parity, and late
cancellation with exact retry.

### Finite rigid payloads at material ports

`MovingShellModel::new_with_payloads` binds finite masses and complete symmetric
COM inertia tensors to immutable `ShellLoadPoint` attachments. The point offset
locates the payload center of mass; the original P1 translation plus director
rotation supplies its velocity. Each payload contributes the consistent form
`B^T diag(m I, I_COM) B` to the original full mass before support restriction.
This retains inter-node inertial coupling, the reference parallel-axis terms,
and reciprocal translation/rotation coupling. It is not nodal mass lumping.

Payload mass must be finite and positive. COM inertia must be finite, exactly
symmetric, and physically realizable: the normalized second moment
`trace(I_COM)/2 * identity - I_COM` must admit a positive-semidefinite
factorization, enforcing the principal-inertia triangle inequalities. Zero
COM inertia is admitted for a point mass. Negative or numerically unresolved
pivots refuse without clipping. Derived mass must remain finite and symmetric;
the original bounded algebraic diagnostics are refreshed for the augmented mass.

The original shell stiffness and damping are preserved exactly; payload mass
does not acquire fictitious Rayleigh drag. Full and reduced augmented inertia
feed the existing initial-equilibrium solve, prepared Newmark factors, driven
pin condensation, reactions, kinetic storage and linear/reference-angular
momentum inventory. Payload assembly visits bounded 18-DOF point stencils
once and adds no DOFs or per-step payload loops. Payload count and each endpoint
load count are separately capped by the model's admitted point bound.

Payload attachment and reference orientation remain fixed. Uniform gravity can
be supplied explicitly through the body-acceleration entries below; other
external payload forces use the existing COM force port. Applying body gravity
and the same payload weight through its point port would count that weight twice.
This remains a small-rotation attached-inertia model, without sliding,
collision, finite-rotation gyroscopic dynamics or finite support-carrier
motion. Empty payload input preserves the original path exactly. Focused
G0/G1/G3/G4 regressions in `tests/moving_load/payload.rs` assert independent COM
energy and angular-momentum identities, native forced acceleration, driven
support impulse/work, physical tensor admission, parity, and atomic retry.

### Density-consistent prescribed body acceleration

`project_with_body_acceleration`, `start_with_body_acceleration` and
`advance_with_body_acceleration` superpose bounded point wrenches with an
explicit uniform Cartesian acceleration in m/s². The original full augmented
mass acts on a repeated translation field with zero director acceleration:
`f_body = M_full * [g,0,g,0,...]`. Graded shell density, payload COM offset
moments and restrained rows remain in the same physical load. No separate
payload-weight loop, scalar density substitution or reduced-mass projection
is used. Zero acceleration preserves the original point-load path exactly.

Three rigid-translation mass columns are cached during model preparation and
refreshed once after payload attachment. Endpoint evaluation is linear in
full DOFs and consumes no point-load budget. The native initial-equilibrium,
prepared Newmark, impulse, support-reaction and work owners consume the
combined load without new factors or separate accounting. The driven variants
accept body and support accelerations independently; full-mass condensation
retains their physical difference and the original separate support work.

Every new endpoint supplies its acceleration explicitly. A held endpoint load
survives instantaneous impulses; the original `advance` methods supply no body
acceleration at their newly supplied endpoint. Nonfinite or unrepresentable
loads and cancellation leave accepted state, held load and reactions unchanged.
These are prescribed uniform loads on the retained reference inventory in the
existing small-rotation shell model, without self-gravity, spatially varying
fields, accelerating-coordinate terms or configuration-dependent gravity.

Six focused regressions in `tests/moving_load/body_acceleration.rs` cover
independent density-weighted mass/centroid/potential moments, offset payloads,
free fall, stationary reactions, falling driven pins, changing acceleration,
held forces through impulses, point-port parity, zero-field bits and atomic
refusal/retry for both native shell mass formulations.

### Localized instantaneous impulse ports

`MovingShellMotion::apply_impulses` accepts bounded physical linear and free
angular `ShellPointImpulse` inputs at the current accepted time. Forces and
impulses share the exact original material-point transpose, including offset
moments. The initialization mass LU is retained and reused to solve
`M_ff delta_v_f = J_f`; attached payload and inter-node inertia are preserved.
The event uses quadratic-cost actions and triangular solves, with no new
factorization, contact solver, dependency edge or generalized coordinates.

Original pin velocities remain continuous during an event. Their impulse on
the shell is `(M delta_v - J)_c`, and its work uses that held prescribed
velocity. Every receipt uses the velocity difference actually accepted after
floating addition, full original inertia, and mean material velocity.
Kinetic change is `v_mean^T M delta_v`, avoiding subtraction of large absolute
energies. A free-row backward-error check scaled by the impulse and absolute
mass-row products refuses unrepresentable kicks; it does not normalize or
report the ideal solve increment as delivered momentum.

Time, step, displacement, support motion and held endpoint forces stay fixed.
The original force equilibrium refreshes acceleration and ordinary support
reactions before one atomic publication. Invalid input, a failed physical
solve, an unrepresentable kick or late cancellation retains the accepted
motion. Viscosity supplies no impulse over a zero-duration event. The caller
owns event timing, collision/restitution laws and duplicate prevention; the
shell retains its existing reference-angular, small-rotation scope.

G0/G1/G3/G4 regressions in `tests/moving_load/impulse.rs` independently inventory
bare-shell and off-center anisotropic payload momentum/energy, driven-pin
reaction/work, zero-free-impulse support events, unchanged force continuation,
actual rounded increments, and atomic refusal/retry.

### Coupled point inverse inertia

`MovingShellMotion::point_inverse_inertia` prepares the complete physical
`B M_ff^-1 B^T` response for an ordered material-point set using the retained
native mass LU. Each point contributes six rows and columns: Cartesian linear
then angular components. Linear/free-angular impulses in N s / N m s map to
material/angular velocity increments in m/s / rad/s; the blocks therefore
retain their mixed physical units. `ShellPointInverseInertia::apply` applies
simultaneous impulses with every cross-point and offset coupling preserved.

The original force/impulse transpose and a shared observation sampler own
both maps. Finite payload inertia and prescribed-pin constraints stay in the
original mass solve. Constrained or redundant directions remain singular;
no response inversion, diagonal replacement, symmetrization or regularizing
stiffness is introduced. Prescribed support velocities contribute no affine
event increment. This is a linear response before accepted-velocity rounding;
actual event admission and support impulse/work remain with `apply_impulses`.

Checked dimensions and byte counts, the original point budget, and the
original `max_matrix_entries` bound dense allocation. Preparation requires
six native mass solves per point, with no new factorization. Callers can
reuse the immutable operator while the model and material points remain
fixed; it need not be rebuilt each step. Preparation, refusal and cancellation
leave accepted motion, loads and support reactions untouched. The existing
reference-frame, small-rotation scope and caller-owned contact law remain.

G0/G1/G3/G4 tests in `tests/moving_load/inverse_inertia.rs` compare complete
cross-point blocks against analytic nodal inertia, independently inventory
payload work, check simultaneous physical events and reciprocity, and retain
constrained/redundant null directions, dense budgets and cancellation purity.

### Native simultaneous normal impacts

`prepare_normal_impacts` binds fixed material points, unit separating normals
and native law controls to the original shell model. One directed normal
mass solve per contact constructs the complete `J_f M_ff^-1 J_f^T`, including
all cross-contact and finite-payload coupling. The retained LU, original
point transpose and shared material sampler own the response. The prepared
matrix is reusable; the native 64-contact cap, original point budget and
checked original matrix-entry/byte limits bound its allocation.

`apply_normal_impacts` receives explicit instantaneous plane positions and
surface velocities. It derives signed gaps and relative normal speeds from
actual shell observations; penetration refuses. Every geometry-active row,
including initially separating rows, enters one original
`fs-contact::multicontact::solve_normal_impacts` call. That owner supplies
Newton restitution, complementarity, global relative-work admission and
SPD/PSD determinacy. No pairwise impact sequence or restitution change is
introduced; zero restitution gives its perfectly inelastic projection.

One aggregate native impulse trial refreshes force equilibrium but retains
accepted motion until measured endpoint speeds, complementarity, full-mass
linear/reference-angular impulse, and kinetic/source/support work pass.
Obstacle work and original driven-pin work remain separate. Redundant
impulses retain the native `UniqueVelocity` classification: individual loads
and separate source/support work shares describe a checked representative,
while combined physical work is audited. Incompatible prescribed motion
still refuses through the native range/work admission.

Displacement, held forces, support motion, time and step count stay fixed.
Every geometry, law, impulse, accounting or cancellation refusal is atomic.
Zero-impulse events leave even acceleration and regular reactions unchanged.
Relative mechanical loss has no automatic heat/damage allocation. Event
timing and contact-set completeness remain caller-owned; this retains the
original reference-frame, small-rotation shell scope without CCD, penetration
repair, sustained forces or finite-rotation claims.

The three native groups in `tests/moving_load/normal_impact.rs` cover analytic
simultaneous offset impacts and global restitution refusal, independently
condensed anisotropic payload inertia with moving-plane/driven-pin work,
and actual post-trial audit refusal, cancellation/retry, redundant contact
motion, incompatible source work, model identity and no-impulse parity.

### Native simultaneous Coulomb shell impacts

`prepare_frictional_impacts` binds stable contact IDs, fixed material points,
orthonormal normal-first Cartesian frames, authored friction coefficients
and the existing `NormalVelocityLaw` to the original shell model. Three
directed mass solves per contact reuse its retained LU and exact point
kinematics/transpose. The complete contact-major `3N x 3N` mobility preserves
normal/tangent, cross-contact, finite-payload and support coupling, including
null directions. The native 21-contact cap, original point budget and checked
matrix-entry/byte limits bound preparation; the immutable result is reusable.

`apply_frictional_impacts` derives gaps and full relative velocity from actual
material-point observations and explicit instantaneous plane motion. One
`fs-contact::frictional_multicontact::solve_frictional_contacts` call owns the
non-associated Coulomb disk, Newton/target normal law, simultaneous stick/slip,
passivity and SPD/PSD determinacy. Geometry-active separating contacts remain
in the coupled solve. Stable IDs canonicalize native solving, world-impulse
scatter and work accumulation without replacing the original law or witnesses.

One aggregate native impulse trial remains private until its actual rounded
endpoint passes speed, normal complementarity, friction-cone and final-slip
checks with the native velocity, impulse and work-unit allowances. Actual
relative work retains the native energy allowance and work-scale roundoff
term. Full-mass linear/reference-angular momentum and kinetic energy are
checked with separate prescribed-plane and driven-support work. Redundant
contact impulses retain representative load/source shares and native total-work
authority; neither individual pressure nor per-contact heat is inferred.

Acceptance holds displacement, time, support motion and force inputs fixed,
refreshing native force acceleration/reactions. Any geometry, law, trial,
accounting or cancellation refusal preserves accepted motion. The caller owns
event discovery, contact completeness and friction-coefficient evidence.
The original small-rotation/reference-frame scope remains; this event consumer
does not implement sustained friction, accumulated partial slip or thermal loss.

The three physical groups in `tests/moving_load/frictional_impact.rs` establish
analytic off-center sticking and slip reversal, simultaneous shared-node
contact with input-order and tangent-frame invariance, and independent finite
payload momentum/energy with moving planes, driven pins, late refusal,
cancellation/retry and unchanged force-continuation behavior.

### Reciprocal Coulomb impacts between finite shell trajectories

`prepare_reciprocal_impacts` combines two original `MovingShellMotion` point
responses into the complete relative mobility
`J_A M_A,ff^-1 J_A^T + J_B M_B,ff^-1 J_B^T`. Both finite inertias, payloads,
supports and cross-contact blocks remain active. Each original body matrix is
retained for independent endpoint checks; their sum enters the same native
Coulomb multicontact owner. Preparation reuses both original mass LUs and checks
the aggregate storage for all three matrices against each model's matrix budget,
alongside existing point and native 21-contact limits.

Prepared pairs bind both model identities, original material witnesses and
fixed reference-world frames. Reference attachments must agree within the
explicit componentwise length allowance. At an event, accepted times must match,
actual deformed witnesses must agree along both tangents, and their observed
normal separation supplies the original activation/penetration law. Neither
attachment is relocated and no compensating couple is added. Nonzero common
displacement and prior deformation are admitted; this is the original linearized
shell interface, not a finite-rotation spatial-angular closure.

`apply_reciprocal_impacts` solves all paired contacts once and stages equal/opposite
point impulses through both native full-mass owners. Each body's actual rounded
point increment, momentum residual and energy residual is checked independently,
so opposing endpoint errors cannot cancel. Combined momentum excludes only the
actual external support impulses. The net commanded reference contact couple is
reported and independently gated, even when the reference-position allowance
would admit a finite offset. Both actual support-work accounts remain external;
there is no prescribed infinite-inertia obstacle supply.

Actual relative midpoint work, native Coulomb residuals/passivity, total kinetic
change and the combined support-work balance gate a final infallible double
acceptance. Any first/second trial refusal, accounting failure or cancellation
leaves both complete accepted trajectories and held inputs unchanged. Native
determinacy and total-work witnesses are retained; redundant individual loads and
support-work shares remain representative. Contact discovery, sustained friction,
finite rigid motion and per-contact heat remain separate owners.

The three groups in `tests/moving_load/reciprocal_impact.rs` check analytic unequal
finite-mass restitution/recoil, simultaneous off-center spin and tangential
transfer with ID/frame invariance, and independent two-body payload inventories.
Deformed trajectories at nonzero accepted time retain different driven-support
work and native force continuation. Permissive geometry cannot hide a reference
couple; strict late audits, final cancellation/retry and an unrepresentable
second-body recoil all preserve both accepted states.

## Momentum-stable shell integration

`shell_time::step_newmark` solves the effective mass equation for acceleration,
then adds complete displacement and velocity increments. Absolute positions are
never differenced to reconstruct acceleration; large rigid drift and short time
steps therefore do not inject artificial momentum. The original `fs-la`
partial-pivot LU owns the solve, without a dimensional absolute pivot cutoff.
For average-acceleration parameters and initially equilibrated acceleration,
the full consistent mass carries the trapezoidal force impulse and the discrete
force-work balance includes internal viscous dissipation. Other Newmark
parameters do not inherit those conservation or unconditional stability claims.
The G0/G3 `shell_time` tests cover free drift, applied impulses, off-diagonal
mass, internal spring/damper reactions, work and physical coefficient scaling.
The stiffness action retains the Newmark predictor increment separately from
the accepted absolute displacement, subtracting the static preload before its
small restoring force. This prevents a rounded absolute predictor from changing
the damped rate update about an equilibrated preload. The regression checks
linear superposition with coupled mass and nonproportional viscosity, including
strong damping. It does not recover displacement increments already below the
precision of the stored absolute position or certify a long unresolved motion.

`shell_time::PreparedNewmark` borrows fixed operators and factors the effective
mass once through the same `fs-la::LuWorkspace` owner. Repeated steps reuse that
factorization with quadratic-cost matrix-vector products and triangular solves;
the one-shot API shares its complete-increment advancement and admission rules.
The prepared system must be rebuilt when operators, timestep or Newmark
parameters change. It can continue a restored checkpoint or another trajectory
without refactoring, with bit-identical results to fresh one-shot solves.
Additional `shell_time` tests exercise full nonproportional damping with
force-work balance, prepared/fresh equality over a forced trajectory, exact
checkpoint continuation after rebuilding factors, and second-order convergence
to a manufactured forced response.

`PreparedNewmark::initialize` solves the complete initial mass equation
`M a = f - C v - K u` through the same LU owner and evaluates the actual
quadratic energies. `step_with_work` accepts both physical endpoint loads and
returns the accepted motion with independently evaluated applied work and
full-matrix damping work. Ordinary prepared `step` calls expose the same
per-interval quantities through `last_step_work`, inferring the old load from
the supplied acceleration equilibrium. Work uses midpoint power quadrature;
it is not an exact continuous-time heat integral. Symmetric PSD damping gives
nonnegative discrete heat, but this generic integrator does not certify PSD or
clip an active operator's signed work. The signed algorithmic energy defect
remains separate from material heat, including for non-default Newmark
parameters and inconsistent startup acceleration. Quadratic energy changes use
endpoint vector differences rather than cancellation-prone subtraction of
large cached energies. Accepted work is replaced per step, preserved on refusal,
and cleared on successful initialization; checkpoint continuation requires no
accumulated work state. G0/G1/G3 regressions cover first-step conservation and
passivity, full-operator force work, independently integrated manufactured heat,
algorithmic versus physical loss, and resolved energy changes on large preload.

## Purpose and layer

Layer: L3 (FLUX). The elasticity core (plan §8.2, bead tfz.13):
small-strain elasticity through finite-strain hyperelasticity, with a
MEASURED locking-free formulation, on both the body-fitted and the
CutFEM-on-SDF frontends. Constitutive laws live in fs-material; this
crate owns body-fitted kinematics, weak forms, boundary operators, and
the globalized Newton loop. Vector CutFEM assembly and solution are
owned canonically by fs-cutfem; fs-solid retains a legacy constructor facade
and a structural-time adapter over the original three-dimensional operators.

`cutfront::dynamics3::CutSecondOrderProblem3` binds Cartesian or adaptive native
elasticity to its exact `MaterialMass3` owner and fixed physical occupancies.
One free-coordinate map removes the same homogeneous strong support rows from
mass, stiffness, forcing and transposes. Adaptive hanging constraints remain in
the native owners. The existing `fs-time` generalized-alpha stepper owns the
accepted state, clock, nonlinear solve and cancellation transaction; the adapter
adds neither artificial support mass nor stabilization inertia. Native material
strain admission runs at initial, trial and accepted states. Its diagonal
preconditioner does not replace the consistent physical mass action.
Cancellation polls run inside native mass and strain callbacks and before/after
the original bounded, indivisible stiffness action.

The same adapter supplies fixed-configuration state/force VJPs and consistent
initial acceleration through `fs-time`. Geometry, material, occupancy, stiffness
scales and supports remain fixed; their derivatives, moving supports, physical
damping, prestress, support reactions and a continuum error certificate are not
inferred. Positive mass diagonals do not establish mass nonsingularity. Focused
G1/G3 regression sources in `tests/cutfront_dynamics3.rs` cover forced rigid
translation and work, clamped Cartesian/adaptive dense-pencil agreement,
forcing pullbacks and transactional refusals. Native execution remains unverified.

The `linear3` module is the generic body-fitted three-dimensional
small-strain rung. It assembles a sparse stiffness/consistent-mass pencil from
P1 tetrahedra plus uniform or per-element `fs-material` state points. Its
constitutive seam is a finite, exactly symmetric, positive-definite 6x6 tensor
in an explicitly documented Mandel basis; isotropic and oriented-orthotropic
constructors lower into that shared representation. Object names and material
presets are deliberately absent: downstream vibration and vibroacoustic
consumers receive the same `(K,M)` artifact for any admitted
geometry/material composition. Anisotropy is therefore supported at the
operator level, but material-card resolution of engineering/crystal constants
and texture/orientation fields remains a separate evidence-bearing admission
step; absence of that data is not permission to assume isotropy.

`TetLinearElasticProblem::assemble_with_nodal_reference_density` admits an
explicit positive P1 density field on the original shared nodes. It produces
the same native sparse pencil for `fs-modal`, with exact cubic integration of
`rho_0 N_a N_b`; a cell mean preserves total mass but not these inertia moments.
The scalar element kernel is shared with finite-strain `P1ReferenceMass`.
Density controls reference mass while the original material field controls
stiffness. Uniform nodal input preserves the original mass bits; connectivity,
support reduction, quality, budgets and cancellation share the original
assembly path. Invalid or unrepresentable density/mass refuses before any
pencil is returned. Total inventory includes supported nodes. The density
field is explicit authored input, not material-source resolution or a coupled
thermal/phase refresh; a closed body's reference field must be retained.
Native modal regression moves density toward the free vertex at identical
total mass/stiffness and verifies the analytic eigenvalue ratio `2/3`.

`HyperTetProblem::reference_mass_matrix_with_nodal_density` carries that same
exact P1 density integration into finite-strain mechanics. It assembles in the
original node/element order, retains supported coordinates, and admits finite
positive nodal density and representable mass under the native mesh, budget
and cancellation checks. The matrix is reference inventory, held through
deformation. It does not refresh density from current volume or infer a
constitutive field from inertia. G0/G1/G3 regressions in
`tests/hyper3_reference_density.rs` cover cubic moments, centroid, kinetic
energy against independent quadrature, positive definiteness, scalar bit
parity, shared-vertex/permuted assembly and refusal boundaries.

The native `hyper3::vibration` dense and sparse-window owners also accept this
original P1 field through their `with_nodal_density` variants, including initial
iterates and held closed-surface pressure. They reuse the same consistent-mass
assembly before nonlinear work, preserve the material and pressure tangents,
and retain every original nodal density. The existing cell-density accessor
reports cell means for inventory only; these means never replace the modal
mass. Neither deformed volume nor a refreshed material density changes this
closed-body inertia. Tests in `hyper3_reference_density.rs` specify exact
finite-stretch modal mass, equal-total-mass redistribution, density scaling,
uniform bit parity, dense/sparse agreement, the independent pressure-loaded
three-mode spectrum, and invalid-density/cancellation refusals. The original
dense DOF cap, sparse spectral budgets and signed-mode limitations remain.

Those same pressure entry points accept the original `TetPressure::Affine`
and `TetPressure::FreeSurface` fields as well as scalar uniform pressure.
The native equilibrium samples the held Eulerian field on the current surface;
the complete pencil subtracts that same owner's analytic follower derivative
exactly once. Field gradients, buoyancy and moving wet-area contributions are
retained with no added fluid inertia, density refresh or separate prestress
matrix. Both cellwise and original P1 nodal reference inertia remain available.
`held_pressure_field()` returns the exact bound surface and field; the existing
`held_pressure()` accessor continues to identify uniform loads. Per-face patches
refuse this symmetric spectral route even if their current matrix happens to
be symmetric: they do not supply a conservative loading-potential contract.
At a whole face on the zero-pressure plane, the native dry-side generalized
tangent is reported, not a unique two-sided harmonic law. Pressure evolution,
finite reservoirs, coupled added mass and continuum stability remain separate.
G1/G3/G4 regression definitions in `tests/hyper3_pressure/vibration.rs` cover
independent affine and clipped-volume modal stiffness, Archimedes waterplane
restoring work, original graded mass, dense/sparse agreement, uniform limits,
invalid fields, patch refusal and cancellation.

`tet10::QuadraticTetProblem::assemble_with_reference_density` carries graded
reference inertia into the same native quadratic stiffness/mass pencil. Each
element supplies four positive corner densities, interpolated with P1 reference
barycentric weights, including on curved physical maps. Material interfaces may
be discontinuous. The ten quadratic displacement shapes retain the full
consistent mass; their negative corner lobes are not used to interpolate density
or to lump mass. Total inventory is integrated before constrained rows are
removed. `assemble_thermoelastic_with_reference_density` changes stiffness with
the admitted thermal material field while retaining that original density.
No temperature-dependent current density refreshes a closed body's mass.

Four mass Duffy points per direction integrate affine graded quadratic mass
exactly (degree five). A non-default profile on a curved map requires six
(degree eight including the cubic Jacobian determinant); this guarantees the
mass quadrature degree, not stiffness accuracy or global map injectivity.
Existing geometry, work budgets and cancellation still apply. A profile equal
to the original scalar material density preserves the existing mass bits.
Native regression checks every mass entry against independent barycentric
moments, translation/rotation inertia, positive definiteness, curved quadrature,
and frozen mass under changed thermal tensors. The actual generalized modal
solve verifies that relocating density at identical total mass and stiffness
increases a supported edge's mass by 25% and reduces its eigenvalues by 20%.

`hex::HexElasticProblem::assemble_with_reference_density` extends original
reference inertia to HEX20 and HEX27. Each element supplies eight positive
corner densities; nonnegative Q1 weights interpolate that reference field
while the original quadratic displacement shapes retain the full consistent
mass. Interfaces may be discontinuous. The resulting native sparse pencil
retains total inventory before support reduction, density-weighted first
moments and rotational inertia. Thermal assembly can change local stiffness
through `assemble_thermoelastic_with_reference_density` while the same graded
reference inventory remains fixed.

The explicit mass quadrature is independent of stiffness quadrature. Equal
rules reuse one prepared point array. Three Gauss points per direction
integrate affine graded quadratic mass exactly; other new-profile geometries
conservatively require six for coordinatewise degree ten, including the
quadratic map's Jacobian determinant. Exact nodal affine relations determine
that admission. These bounds certify mass polynomial integration, not global
map injectivity or stiffness accuracy. An identical scalar profile with the
original rule preserves prior mass bits; a different requested mass rule is
honored. Native regressions compare every mass entry with independent monomial
integrals on affine and curved maps, retain first/rotational moments, expose a
five-point underintegration case, and check exact supported eigenvalue shifts
at equal total mass and stiffness. Existing geometry, budget and cancellation
admission remains in force.

`linear3::TetElasticMaterial::acoustic_response` exposes the shared
`fs-material::elastic::acoustic` law using the exact native Mandel stiffness,
density and material identity. Source-resolved isotropy, oriented orthotropy and
full tensors therefore derive directional wave speeds/impedances and hydrostatic
compressibility/strain from the same coefficients as assembled mechanics. The
caller declares an explicit strain ceiling, as for stress recovery, without
extending source validity. This is a frozen small-strain model query; static
moduli do not become measured ultrasonic data. Native regressions connect
source-selected bulk/shear constants and oriented full-tensor states to these
derived responses while preserving original identities.

The shared `linear3::source` adapter resolves explicitly selected isotropic
Young/Poisson, Young/shear or bulk/shear source pairs, or oriented
orthotropic elastic constants, and independently selected isotropic or three-axis
integrated expansion. Both `fs-mbd-material` and spatial
`fs-couple::solid_enthalpy` use this one lowering into the original Mandel
tensor and eigenstrain owners. Current stiffness does not require current
density: the caller supplies conserved reference mass per reference volume.
Reference resolution can select density and the full elastic schema atomically.
Every public lowering checks a complete matching property schema before
accessing engineering constants. Source integration, compliance admission and
proper-frame transforms remain in their original owners.

All three isotropic source choices lower through
`fs-material::state_point::ResolvedIsotropicElasticModuli`; derived E/nu values are not
inserted into the card or relabeled as measurements. The selected pair, source
receipts and reference inventory remain bound to the native material identity.
The source owner's Mandel tensor carries the retained K/G or E/G moduli directly
into mechanical, modal and thermal assembly. It never reconstructs them through
a rounded Poisson ratio near incompressibility or the auxetic limit. Canonical
E/nu input preserves its original tensor arithmetic. The same direct K/G tensor
construction serves the acoustic owner; this numerical correction does not
remove the existing displacement elements' volumetric locking.
Changing temperature re-resolves the same selected source schema before
rebuilding stiffness, thermal stress and vibration inputs, while reference mass
and inertia remain unchanged. G1/G3 tests exercise equivalent source pairs,
temperature-dependent compliance, thermal stress, exact mass retention and
source/identity refusals through the existing tetrahedral solve. This does not
extend the separate reciprocal affine thermoelastic potential's E/nu derivative
law: that consumer explicitly refuses alternative pairs until its thermal
derivatives are supported.

TET10 and HEX20/HEX27 expose `assemble_thermal_load` using the same
isoparametric geometry, Mandel strain operator, elastic tensors, constraints
and fully integrated stiffness rule as their mechanical assembly. Existing
`TetThermalStrainField` values supply total eigenstrain, uniform or constant per
element, bound to the exact elastic material state. The sign is
`K u = f_external + integral B^T C epsilon_free dV`. The existing `TetThermalLoad`
carrier retains both full support rows and the canonical reduced free-DOF order,
so the existing solver can compute compatible thermal displacement. Load identity
binds geometry, connectivity, source states, actual elastic coefficients,
quadrature and supports. Both stiffness and mass quadrature locations retain the
original geometry-admission checks. Cancellation and numerical/material refusal
publish no partial load. This does not add a thermal transport solver, within-cell
temperature interpolation, geometric prestress, or
an accuracy certificate for curved-map quadrature. `tests/quadratic_thermal.rs`
checks uniform anisotropic eigenstrain including Mandel shear, resultant/moment
balance, thermal work, actual free-expansion solves on all three element families,
source/state mismatch and constrained loads against analytic references.

The same quadratic APIs also accept `QuadraticThermalStrainField::AtQuadrature`.
Each element supplies one source-integrated Mandel eigenstrain at each of its
original stiffness quadrature points, in explicit element-major/rule order.
The field checks both the points-per-element value and complete sample count;
every sample must bind the exact element elastic material identity. The
elementwise input path retains its original arithmetic, load identity and
one-constitutive-contraction-per-element behavior. No temperature interpolation,
constant-alpha approximation or averaging of nonlinear expansion is inserted.
Elastic tensors remain frozen per element. Changing quadrature requires new
samples, and source-data and quadrature error remain the caller's obligations.

Pointwise thermal loading is shared by equilibrium, stress/energy recovery and
the initial-stress vibration pencil. A zero mean strain does not erase stored
energy or the correlation between stress and shape gradients. G1 tests in
`tests/quadratic_thermal_gradient.rs` verify exact compatible thermal-gradient
bending on TET10/HEX20/HEX27, the analytic restrained energy of a zero-mean
through-thickness gradient, and its signed geometric-stiffness quadratic form.
Count/rule and late source mismatches refuse without publishing a partial load.

`resolve_nodal_thermal_strains` connects actual nodal reference/current absolute
temperatures to those same TET10/HEX20/HEX27 mechanics. It uses each family's
native geometry/shapes and stiffness quadrature, interpolates temperature first,
then calls the existing source-integrated isotropic or oriented orthotropic
expansion owner. Nonconstant alpha is integrated over each local temperature
path; nodal eigenstrains or endpoint alpha times temperature change are not
substitutes. The original resolved expansion reference pins the card, selected
claims, nonthermal coordinates and quantity conventions. A source selection is
prepared once per element, and uniform nodal fields resolve their integral once
per element. Returned original source receipts and accepted strains remain in
native element/quadrature order and borrow directly into load/response/tangent
methods. There is no second thermal law, shape kernel, reference-mass update or
implicit elastic-tensor refresh. Invalid counts, source mismatch, unsupported
temperature paths, nonpositive/nonfinite temperatures and higher-order
interpolation overshoot outside source support refuse the whole result; no
clamping hides those exits. Cancellation is checked per element and point.
Temperature transfer from a different thermal mesh/space remains explicit,
and this prescribed-temperature bridge does not claim caloric or reciprocal
energy closure. `tests/quadratic_temperature_source.rs` checks nonconstant-alpha
integration with a nonuniform stress-free reference, real gradient bending and
joint geometry/material-frame rotation, exact constant-field load equivalence,
source/count/domain refusal and positive-nodal temperature overshoot.

`resolve_nodal_thermoelasticity` explicitly couples source-resolved local
elasticity to the same nodal-temperature bridge. `QuadraticThermoelasticSource`
binds the original elastic and expansion references to one card and physical
context, with independent declared tensor symmetries and orientations. Pinned
current elastic properties are resolved before assembly; current density is
never queried. Uniform element temperatures resolve each source law once.
`ThermomechanicalAtQuadrature` carries matching local materials and integrated
eigenstrains at the native stiffness points. The exact local C(T) enters the
actual displacement stiffness, thermal load, recovered stress/energy and signed
initial-stress stiffness. `assemble_thermoelastic` exposes that same stiffness
and original consistent mass. Every local density must equal the original
element density; heating cannot redefine the reference inventory. Local sample
counts, material/strain identity and source support are checked, and the existing
load identity includes every local tensor's identity and coefficient bits.
Legacy elementwise and strain-only field paths retain their previous behavior.

This is the mechanical tangent at prescribed, frozen temperature, not a
coupled temperature derivative, transient thermoelastic feedback or a caloric
energy closure. Source and quadrature validity remain explicit: variable C(T)
and prestress can require a higher integration rule. The independent polynomial
oracles in `tests/quadratic_thermoelastic_material.rs` check stiffness, thermal
work, stored energy, signed prestress, compatible thermal bending and unchanged
mass on all three quadratic families. They also check a stiffness-domain exit
despite admitted expansion data and late local material/count/density refusals.

The same TET10/HEX20/HEX27 problems expose `solve_thermal_response` for source-bound
thermal equilibrium on a frozen reference mesh. The existing `fs-solver::CgState`
and Jacobi preconditioner solve the assembled stiffness with bounded
iteration/cancellation windows and an independently recomputed true Euclidean
residual. Every connected component must have all six rigid motions restrained.
Fully fixed bodies have an exact empty solve, while their nonzero thermal stress,
elastic energy and support force are retained.

Recovery integrates actual `sigma = C (B u - epsilon_free)` and
`epsilon_elastic:sigma/2` at the same stiffness quadrature points, checking both
total and elastic Mandel strain against the caller's finite positive ceiling.
`QuadraticThermalResponse` contains compatible nodal displacement, explicitly
named volume-mean element strains/stresses, integrated energy, and full nodal
internal forces. Those forces are reactions **on the solid** at fixed DOFs and
unmodified equilibrium residuals at free DOFs, with no external mechanical load.
Element means do not estimate pointwise or peak stress; energy integrates the
varying strain field before averaging. The immutable result binds the thermal
load, actual material coefficients, mesh frame, controls and observations;
failure publishes no partial solution or field. Geometry, mass and density are
unchanged. This does not add live thermoelastic transport, geometric stiffness,
finite-strain behavior, stress extrapolation or quadrature-error certification.
`tests/quadratic_thermal_response.rs` exercises free expansion, a piecewise
nonuniform expansion field, clamped stress/energy/support reactions, repeatable
results and refusals for rigid modes, controls and applicability limits.

Quadratic problems also expose `assemble_thermal_load_with_nodal_forces` and
`solve_thermal_response_with_nodal_forces`. Their dead nodal forces use the same
strictly ordered, unique, finite full-DOF validation and combined-load identity
owner as P1. Thermal and mechanical forces enter one stiffness solve. Actual
stress still integrates `C(Bu - epsilon_free)`; `nodal_internal_force_n` retains
that full internal force, while `nodal_force_balance_n` subtracts applied dead
nodal forces. The latter gives reactions on the solid at fixed DOFs and unchanged
residuals at free DOFs, including applied loads directly on supports. The result
retains both the original thermal-load identity and the complete load identity.
`tests/quadratic_thermal_loads.rs` checks analytic thermal expansion plus end
traction on all three families, support-load accounting, and malformed loads.

`linear3::visco` (IRP M1.3-W2a, in progress) prepares P1, TET10, HEX20 or
HEX27 elements with caller-supplied, mesh-frame `fs-material::TensorProny`
laws and density.
It reuses the existing elastic mesh admission, equilibrium stiffness and
consistent mass, and the shared Mandel B operator. Complex stiffness retains
the full storage and loss tensors, including off-diagonal couplings. Every
frequency query must lie in every element's declared law band. Empty branches
reproduce the canonical elastic stiffness bitwise. Positive quadrature weights
propagate pointwise relative quadratic-form bounds to the assembled forms;
this alone is not a frequency-dependent modal-eigenvalue certificate.


`ViscoElements::prepare_harmonic_response` connects the complete material
spectrum and the original consistent mass to the shared `fs-modal` harmonic
owner. At each positive requested frequency it solves
`(K'(omega) - omega^2 M - i K''(omega)) u = f` under `exp(-i omega t)`;
the exact single-frequency identity `C_eff = K''/omega` does not replace
hereditary histories with memoryless damping. The full coupled tensor entries
remain present, and every element law must admit the frequency. Preparation
reuses one dense factor for multiple force phasors: O(n^3) factorization,
O(n^2) retained storage and solves, plus original quadrature strain checks.
Before publication, the scaled 2x2 Gram problem of real/imaginary Mandel
strain admits the maximum norm over the complete cycle at each source point.
Outputs retain displacement/velocity phasors, numerical backward error,
independent mean input work and constitutive heat, and their signed defect.
No dispersive stored-energy formula is inferred from storage stiffness.
Budgets include the three dense input matrices and shared harmonic payload;
source sparse assembly retains its existing limits. Unrepresentable effective
loss, source-domain violations, cancellation and singular pencils refuse.
This bounded dense route does not claim spatial convergence, unsampled strain
extrema, transient response, moving supports or experimental validation.
`tests/visco_harmonic.rs` specifies analytic Prony/Kelvin response and heat,
nonproportional coupled tensor/consistent-mass phasors, complete-cycle strain
refusal, deterministic factor reuse, lossless response, budget/band and
cancellation checks (G0/G1/G3/G4).

`SolidMeshRef::Hex` retains the original hexahedral connectivity and positive
tensor-product quadrature (at least three points per axis). The original
`HexElasticProblem` admits topology and sampled Jacobians and assembles the
reference elastic pencil; the same hex shape gradients feed full-tensor
spectra, projected operators, constitutive histories and integrated heat.
This includes curved HEX20/HEX27 geometry without a tetrahedral remesh.
Jacobian and small-strain checks apply at declared integration points;
they do not certify unsampled geometry or strain extrema.

The required `ViscoStateBudget` counts six internal strain scalars per branch
and stiffness quadrature point, before history allocation. A history borrows
its immutable prepared mesh/laws and starts relaxed at zero displacement.
Prescribed full nodal displacement ramps advance the constitutive owner's
exact exponential update; B-transpose returns complete force rows, including
support reactions, and reduced forces use the canonical free-DOF ordering.
Stored energy, exact ramp work and dashpot heat are integrated in joules.
Work is the material owner's time-integrated work, not endpoint force times
displacement. Candidate histories commit only after all points and assembled
forces/energies succeed; cancellation and invalid/nonfinite inputs preserve
the previous state. Cx is checked per element/point, including state copying.
No damping clipping, diagonal loss approximation or secondary constitutive
law is introduced. Determinism is same-host at fixed element/point order.

`tests/visco_elements.rs` names G0 uniform-strain tensor/element energy,
elastic parity, 4096 counter-stream multiaxial steps, one-branch Zener response,
band/state-cap/invalid-history refusals; G1 harmonic endpoint force on three
time rungs (observed order 2 +/- 0.2); G3 rotation invariance; and a G4
cancellation-preserves-state check. These fixtures are synthetic software
checks and require native execution.
The HEX20/HEX27 fixtures additionally check rotated multiaxial spectra,
work-conjugate full forces, exact material work/heat integration, canonical
elastic mass/stiffness parity, quadrature/state-budget refusals, and curved
HEX27 storage/loss convergence against an analytic spatial integral.

`ViscoElements::quasistatic_operator(dt, cx)` prepares a reusable sparse
consistent tangent for one positive duration; `ViscoElementState::step_quasistatic`
solves endpoint nodal loads and prescribed support motion through the shared
Jacobi-preconditioned CG owner. The original full Mandel Kelvin/Prony tensors,
quadrature and accepted branch histories remain intact. Each step solves free
displacement increments, independently checks the original internal-force
residual, and returns complete displacements, support reactions and constitutive
storage/work/heat. Support motion contributes its actual integrated material
work. All connected components must pass the existing rigid-motion constraint
admission, and every original material band must include DC for sustained-load
creep. Failed solves, source-domain violations and cancellation preserve the
accepted displacement and memory. This is fixed-material quasistatic evolution
without inertia: each affine displacement ramp has exact constitutive
integration, while a continuously prescribed load trajectory requires time-step
convergence. DC admission does not certify that trajectory's transient spectrum.
`tests/visco_quasistatic.rs` checks analytic Maxwell/Kelvin creep convergence,
independent heat, nonproportional coupling, moving supports, tangent derivatives
and transactional refusals.

`ViscoElements::dynamic_operator(dt, config, cx)` prepares consistent-mass
midpoint dynamics on the same P1/TET10/HEX20/HEX27 tensor histories. Its
`start(history, full_velocity, cx)` retains an accepted constitutive preload;
`ViscoDynamicState::advance(mean_force, cx)` consumes the actual interval-mean
nodal load. The solve uses `M + dt² K_mean/2`, where `K_mean` differentiates
the original exact strain-ramp mean force. Geometry, tensor assembly and CG
are shared with the existing owners. Endpoint stresses remain observations;
their trapezoidal average does not replace the work-conjugate material force.
The accepted interval independently checks momentum and
`external_work = change(kinetic + material_storage) + physical_heat`.
No energy defect is relabeled as heat. Exact per-step work and loss survive
large cumulative counters. Explicit step/free-DOF/iteration budgets, source
strain limits, finite arithmetic, cancellation and energy allowances are
checked before publishing motion, memory, work or clock changes.

The step also retains `kinetic_energy_change_j` through the shared original-mass
quadratic increment and `material_stored_energy_change_j` through the original
per-element constitutive increments. Applied and support work retain their
displacement/product tails through the independent residual. Relative energy
admission scales against interval work, storage changes and physical heat;
unchanged kinetic or material inventories cannot hide a local imbalance.
The finite-cell thermal consumer uses these same accepted increments. Focused
regressions include supported acceleration and Kelvin strain heating beside
large unchanged translation; full coupled execution is separately required.

Ordinary `start`/`advance` retain stationary supports and zero support work.
`ViscoElements::new_with_full_reference_mass` additionally retains the original
full consistent mass, including support rows and free/support couplings, for
P1/TET10/HEX8/HEX20/HEX27. The ordinary constructor keeps its original reduced
storage cost. `start_with_support_motion` and `advance_with_support_velocity`
admit explicit support velocities; displacement follows the same midpoint
kinematics as the free coordinates. Complete mean reactions include inertia
and original work-mean material forces. Their signed work enters the existing
independent energy gate separately from applied loads, eigenstrain work and
physical heat. Thermal forcing composes through the same interval, and source
continuation preserves full reference mass and accumulated support work.
Missing retained mass, malformed/free-DOF support input, unrepresentable motion,
source refusal and cancellation leave accepted motion and memory unchanged.
Free bodies are admitted through their original consistent mass, without
artificial constraints. Every
material band must include DC; transient spectral admission and temporal
convergence remain caller obligations. The material ramp is exact, while the
coupled inertial trajectory is a time discretization. `tests/visco_dynamics.rs`
checks elastic energy conservation, analytical Kelvin-oscillator convergence,
independent heat, nonproportional memory, applied work, free-body motion and
transactional refusals. The existing element suite also carries a shared
TET10/HEX20/HEX27 kinetic/material-energy regression.

`ViscoElements::modes` solves the storage-dispersive ordered eigenproblem
`lambda(K'(omega), M) = omega^2` using the existing dense fs-modal eigensolver
and an explicit frequency bracket, iteration/tolerance budget, and a dense
reference dimension cap (at most 256 free DOFs). For positive equilibrium and
PSD Prony branches, `K'(omega)/omega^2` decreases in Loewner order, so ordered
positive-mode roots can be bisected without assuming fixed-point contraction.
Endpoints must bracket the selected mode inside the common material band;
missing roots, exhausted iterations, and an eigenvalue whose certified
quadratic (Kato–Temple) radius exceeds the budget refuse. The radius comes
from `fs_modal::recentred_dense_radii`, because a dense eigenvector's linear
residual stays near `ε·λ_max` on a stiff mesh and is no measure of the
eigenvalue; the linear residual is still reported.
The output retains the mass-normalized shape, full loss-energy quotient,
numerical frequency interval and nonlinear residual. Numerical bisection is
not an outward-rounded modal certificate, and the loss quotient does not
certify a diagonal damping approximation. A one-coordinate P1/TET10 Zener
fixture compares the self-consistent frequency with an independent analytic
quadratic (2e-10 relative), and deliberately shifts it by more than a factor
of two from the equilibrium result.

This module does not yet supply the fractional-target modal-loss ladder,
certified modal block reduction, valve-card integration, or material/specimen
validation. Arbitrary
time-domain input is not a certification that its spectrum lies in the law's
band. Intrinsic material heat excludes support, contact and radiation losses.

`linear3::tet10` adds a separate isoparametric quadratic displacement/geometry
operator (IRP M1.3-W1a). `QuadraticTetProblem` consumes ten-node connectivity,
the existing uniform/per-element material tensor, fixed DOFs and assembly
budget, and a required `TetQuadrature`. Canonical order is four vertices, then
edges `(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)`. Shared edges must share exactly one
midside node; a midside node cannot also be a corner or represent another edge.
The signed Jacobian must be finite, positive and above the scaled quality floor
at every stiffness and mass quadrature point. The P1 arithmetic, thermal loads,
thermal prestress and stress-recovery paths remain unchanged.

`linear3::hex` (IRP M1.3-W1b) adds the fully integrated isoparametric HEX20
serendipity and HEX27 tensor-Q2 families. They use the **same** const-generic
Mandel/consistent-mass/COO assembler extracted from TET10, the same material
field, fixed-DOF map, `TetElasticAssembly` carrier and fs-modal eigen path.
Canonical node coordinates and Gmsh-12/17 permutations live in L2 fs-mesh's
`structured_hex`; fs-solid now has a runtime dependency on that owner. The
historical `maximum_tetrahedra` budget counts hex elements too. Family and
connectivity cardinality must agree. Shared corners/edges/HEX27 face centres
have one global node per topological entity; duplicate cells and conflicting
node roles refuse. This checks conformity, not arbitrary cell intersections.

`HexQuadrature` requires 3..=12 Gauss-Legendre points per direction for both
operators. A 2x2x2 rule refuses `UnderIntegrated`: a test-only contraction
counts its HEX20 hourglass modes and cross-checks full-rule stiffness. Mass
is always consistent, SPD on admitted test meshes; HEX20 corner row sums are
negative and row-sum lumping is not available. Every stiffness/mass Gauss
point checks signed Jacobian and the declared scaled floor. General curved
maps still require a quadrature ladder; 3x3x3 is exact for affine mass/stiffness,
not a curved-map accuracy certificate. The common kernel retains TET10's
ordered local triangle arithmetic and polls Cx per element and quadrature point.

`tests/hex.rs` specifies G0 partition/Kronecker/gradient identities (1e-14),
six rigid modes, exact mass/first moments (1e-14), quadratic bending completeness
(1e-10), topology/Jacobian/underintegration/cap/cap+1/cancellation refusals, and
Gmsh permutation checks. HEX27 reproduces quadratic physical fields under the
tested trilinear distortion; HEX20's loss there is logged, not qualified away
(Arnold, Boffi & Falk 2002, https://arxiv.org/abs/math/0005036).
G1 roller-face P waves use an fs-mms OrderGate on three x-refinement rungs.
G2 square free-free bars (L/h=50, nu=.3) gate the first three bending pairs and
first axial mode to Timoshenko/Cowper and c/(2L), respectively (0.1 percent),
after a three-rung fs-mms-fitted Richardson frequency estimate <=0.02 percent.
The Timoshenko characteristic is solved by deterministic Newton and checked
against independent free-end shear/moment residuals; Euler-Bernoulli is a
diagnostic only. The scalar oracle follows the cosh-scaled free-free equation
in Khasawneh & Segalman 2019, Table 7, https://arxiv.org/abs/1811.03222.
G3 compares the first six HEX27/TET10 modes on the same bar within the sum of
their three-rung discretization estimates plus propagated eigen residuals.
These are estimated discretization bounds, not certified continuum enclosures.
All new numerical gates require central native execution before qualification.
No mixed/incompressible stabilization, reduced-integration stabilization,
finite-strain hexes, inside-element material variation or global injectivity
claim is introduced. Distorted bending-dominated parts require HEX27/TET10;
HEX20's affine benchmark does not qualify it on arbitrary distorted geometry.

Stiffness and mass use the existing `fs-feec` positive-weight collapsed Duffy
rule. Required counts are respectively 3..=12 and 4..=12 per direction: affine
stiffness has degree two and consistent mass degree four, while the Duffy map
adds degree two. Mass is never row-sum lumped. Curved-element stiffness is
rational and needs an independent quadrature ladder; specifying a count does
not certify its error. `tet10_shapes` and `tet10_point` expose the same mapping
for consistent loads and independent norm integration. Assembly is ordered by
element and quadrature point, stores symmetric local triangles identically,
and polls the caller's `Cx` before allocation and at each quadrature point.
Mesh/free-DOF caps refuse before allocating maps or matrix storage; no partial
matrix is published after an error or cancellation. Replay is same-host
deterministic; no new cross-ISA claim is made.

`promote_p1_to_p2` gives each sorted global endpoint pair one midpoint ID,
independent of element order. It canonicalizes element orientation and applies
an optional caller-supplied source-surface projection only to boundary edges.
Straight-sided approximation of curved geometry is a separate error axis.
Positive sampled Jacobians do not prove injectivity between samples; projection
does not supply source-surface or specimen authority. Finite-strain and mixed incompressible formulations remain outside this
quadratic operator's scope. Its explicit initial-stress extension is described below. Thermal stress observations are limited to
the explicit volume averages returned by `solve_thermal_response`.

Valve reduction requires `ValveCardRequest::element_order`: `Tet4` or explicit
`Tet10 { elements, quadrature }`. Quadratic corner sets must match the original
per-element topology, and every observation/pressure-load node must refer to
the supplied complete quadratic mesh. No boundary condition or load is
implicitly interpolated from P1. `ValveCard` v2 retains the discretization;
its identity binds order, quadrature and full quadratic connectivity. The v2
decoder refuses v1 bytes rather than reinterpreting them. This changes card
identity, not the legacy P1 numerical operator or the rl-001 authored band.
Reduction remains an explicitly truncated exploratory consumer, not an IRP
acceptance rendering engine.

### Thermal equilibrium in pressure-loaded valve reduction

`reduce::thermal::reduce_pressure_loaded_valve` carries the original P1 or
TET10 thermal/dead-load equilibrium into the existing valve/reed consumer.
The material field, integrated eigenstrain and explicit dead nodal forces
use the original solid owners. The common reducer receives that equilibrium's
`K_elastic + K_geometric` pencil; both its modal spectrum and exact opening/
pressure compliance therefore change with the same actual stress field.
TET10 retains pointwise quadrature stress, not element-mean stress.

The existing LDLT inertia owner must admit the WHOLE tangent as positive
definite. Selecting a positive modal window cannot hide an unstable direction
outside it. A singular/unstable pencil refuses instead of converting negative
eigenvalues into frequencies. This is numerical admission within the original
factorization tolerances, not a continuum buckling certificate.

Reference geometry gradients, mass and homogeneous incremental supports remain
unchanged. The aperture gap, width and thickness are measured from the accepted
displaced nodes against the original fixed observation axes and opposing plane.
The original P1 geometry-update owner checks its displacement bound and moved
cell quality. TET10 uses that declared displacement bound and the existing
isoparametric point evaluator at both original quadrature rules. Positive sampled
Jacobians do not establish injectivity between points. Closed apertures, excessive
strain/displacement and incompatible source states refuse; no contact model is
selected implicitly.

The pressure-force-per-pascal field remains an explicit FIXED reference load.
Its work-conjugate modal areas use the same new modes; surface-normal updates,
follower stiffness, pressure-dependent equilibrium and changed load geometry
are not inferred. The loss factor remains authored. This is the conservative
small-displacement initial-stress approximation, not finite-strain dynamics,
a source loss-spectrum model or live modal-state transfer.

`ThermalPressureLoadedValve` retains the original equilibrium/tangent,
displaced observation nodes, ordinary pressure-loaded valve and a combined
identity. The existing `BeatingReedSource` ingress compiles its actual pressure
port after the original first-mode compliance admission. The ValveCard codec
is unchanged; its ordinary material/mesh fields keep their original meaning.
The complete thermal receipt must accompany that card to retain thermal input
provenance. Every refusal or cancellation returns no partially minted reduction.

G1 definitions in `reduce/tests/thermal.rs` use a six-tet unit cube with one
free corner displacement. Independent barycentric integrals give closed-form
P1 and TET10 mass, stiffness, thermal/dead-load displacement, signed geometric
stiffness and pressure compliance; the same result reaches the playing-reed
ingress. G3 checks zero-preload parity and an unstable disconnected direction
outside the selected positive window. G4 covers displacement, closed-aperture,
source mismatch and cancellation refusals. These are synthetic implementation
checks, not measured valve/material validation.

The focused `tests/tet10.rs` target specifies G0 shape, quadrature, consistent
mass, rigid-mode, distorted bending patch, rotated-tensor, refusal and budget
checks; G1 MMS (L2 order 3 and energy order 2, ±0.2), projected-cylinder geometry
and free-bar eigenvalue ladders; and the Cook/MacNeal–Harder G2 bands. The
`reduce::reduce_tests::rl_005_tet10_cantilever_card` regression specifies a
fivefold static-compliance improvement over P1 at matched DOFs on the rl-001
slab geometry. `scripts/e2e/irp_s1_structures.sh` retains their structured
receipts. These new probes require central native execution before any DEF-08
closure claim. The Euler-disc structural consumer and Bar-S1 descriptions
require their own explicit migration; adding this operator alone does not
establish those consumers' convergence or remove their P1 limitations.

The free-bar per-commit test computes r=1,2,4 and explicitly compares against
the historical ts1 r=8 diagnostic identified in the .58 bead; it does not
qualify that historical rung. `g1_free_free_bar_eigen_ladder_r8_qualification`
computes all four rungs (the first diagnostic cost 2101 s on ts1). Its gates
require rising observed order toward the theoretical four, finest order at
least 3.5, and eigen residual uncertainty at most 2% of each finest difference.
Successive Richardson limits are compared with an explicitly estimated
pre-asymptotic truncation allowance plus propagated eigen residual uncertainty.
Mass scaled by `(1+h²)` must fail the same trend/limit gate. MacNeal–Harder
alone uses `8 eps kappa_1_est ||b||_1` for the ordinary-SpMV residual, with a
five-step inverse-norm estimate and a local 1e-6 displacement falsifier. This
is a conditioning estimate, not a verified inverse bound; all other static
fixtures retain their original relative Euclidean residual limit of 1e-9.

`TetElasticMaterial::try_new_oriented_tensor` admits complete anisotropic
stiffness data into that same operator. `ElasticTensorBasis` requires explicit
tensor-shear, engineering-shear or Mandel notation, one of the documented
`xx,yy,zz,xy,yz,zx` / `xx,yy,zz,yz,zx,xy` orders, and a nonzero source-frame
identity. The caller also supplies a target-frame identity, a proper rotation
whose columns are source axes in the target frame, density and an upstream
material-state identity. A frame mapped to itself requires identity rotation.
The constructor checks full major symmetry in work-conjugate coordinates and
positive definiteness of the entire source and target law. Tensor-shear
matrices with normal/shear coupling are generally not ordinarily symmetric;
engineering matrices transform as `C_M = S C_E S`, with
`S = diag(1,1,1,sqrt(2),sqrt(2),sqrt(2))`. Negative off-diagonal coefficients
are admissible when the whole law is stable. This does not apply a stiffness
positivity rule to arbitrary off-diagonal blocks of coupled constitutive laws.

`ElasticTensorTransformReceipt` retains the exact source matrix and descriptor,
both frame identities, rotation and upstream material identity alongside the
executable target-frame material. The material's identity binds those fields,
density and every output coefficient under
`org.frankensim.fs-solid.elastic-tensor-frame-transform.v1`. Rotation reuses the
existing Mandel basis transform (equivalent to the four-index tensor law) and
its deterministic triangle storage; source asymmetry is refused before that
step. No missing coefficient, noisy-source symmetry projection or frame
calibration is inferred. G0/G3 tests compare all three notations and both
orders against a separate four-index rotation with a fully coupled synthetic
tensor, check whole-law instability and frame/identity refusals, and retain
the existing orthotropic tests. A G1 affine single-tet test compares actual
nodal forces and strain energy against direct tensor contraction. This is
software evidence, not measured crystal/wood validation. The receipt has no
portable codec yet. `from_resolved_elastic_tensor` now connects the complete
`fs-material::state_point::ElasticTensorStatePoint` directly to this constructor.
The shared basis descriptors are owned by L1 `fs-matdb` and re-exported here.
Source compiler v6 and normalized-pack v5 preserve each coefficient context;
the state resolver requires all 36 entries plus density. The transform binds
the bundle identity, hence every selected receipt and the actual query point.
The G1 tet regression also decodes a portable material card, resolves the
complete tensor, and reproduces the independently checked assembled operator.
Coefficient uncertainties remain in the upstream receipts; no independent
coefficient assumption, propagated tensor bound or empirical authority is added.

The committed `silicon-cubic-25c-nasa-rp1057` source bundle now supplies a
separate sourced witness in `xtask/tests/matdb_pack_cli.rs`. It compiles all
36 explicit engineering coordinates, persists and reloads the material card,
and resolves its density/tensor into the same oriented tetrahedral operator.
The G1/G3 test checks mass, nodal forces and strain energy at zero and 45-degree
orientations against an independent three-constant cubic stress calculation.
The matrix entries derive from three handbook constants at 25 C; density is
from a separate pure-crystal reference. This cross-source engineering example
retains unknown pressure and specimen/dopant conditions. It neither qualifies
a particular wafer nor supplies a thermal/electrical silicon model.

`TetThermalStrainState` and `assemble_thermal_load` add the matching generic
small-strain thermomechanical load boundary. An upstream material law supplies
the already-integrated, possibly anisotropic Mandel thermal strain between an
explicit positive reference and current temperature. Each state is bound to
the exact elastic material-state identity and a nonzero thermal-law identity;
the isotropic binder consumes fs-material's exact signed
`integral(alpha(T) dT)` only when its current endpoint and tangent elasticity
share the exact card and complete state point. Otherwise it refuses before
assembly.
the operator assembles `integral B^T C epsilon_free dV` on the same tetrahedra
and publishes both full nodal forces and the constrained free-DOF projection.
`solve_thermal_displacement` advances a Jacobi-PCG solve in cancellation-bounded
chunks and accepts only an independently recomputed true Euclidean residual;
it reconstructs the full zero-Dirichlet nodal displacement field for an
explicit later geometry-update transaction. That transaction is
`update_geometry_from_displacement`: it identity-checks the reference mesh,
enforces a caller displacement ceiling, preserves connectivity, and refuses
every degenerate or orientation-flipped tetrahedron before publishing the
updated fixed-topology chart.
Fully fixed homogeneous-Dirichlet meshes now admit their exact zero-displacement
thermal response with zero iterations and zero reduced residual. The thermal
load still validates all geometry/material/state inputs and retains its full
force vector; recovered stress, support reactions and stored elastic energy
can be nonzero. Matrix/modal `assemble` continues to refuse `NoFreeDofs`.
A G1 fully clamped isotropic heating test checks hydrostatic compression,
reactions, energy and unchanged geometry against closed-form values; invalid
thermal counts and degenerate geometry still refuse before the empty solve.
No material name, constant-alpha assumption, or Euler-disc geometry appears in
this layer.

`TetThermalStrainState::try_from_oriented_strain` connects a source-frame
symmetric strain to that same thermal operator. `StrainTensorBasis` declares
physical tensor shear, engineering strain shear (`2 epsilon_ij`), or Mandel
shear and either existing six-component order. These are strain conventions,
not stress-vector or elastic-matrix scaling. The numerical transform is the
second-order law `epsilon' = Q epsilon Q^T`. The admitted
`ElasticTensorTransformReceipt` supplies the target frame and elastic-state
identity; a source frame mapped to itself requires identity rotation. Identity
rotation preserves normalized Mandel bits. Nonfinite input, improper rotation,
shear overflow, complete nonzero-to-zero rotation underflow, and target strains
outside the declared Mandel component ceiling refuse.

`ThermalStrainTransformReceipt` retains the source six-vector, basis, both
frames, rotation and executable thermal state. State identity domain
`org.frankensim.fs-solid.thermal-strain-frame-transform.v1` binds the ordinary
thermal-state identity (temperatures, law, elastic state, target strain and
validity ceiling) plus the source vector and complete transform. The existing
constructor's identity preimage is unchanged. G1/G3 compares all six source
conventions against independent second-order rotation, trace/Frobenius
invariance, and four-index elastic stress contraction; actual thermal loads,
constrained displacement and fixed-topology geometry reproduce affine
anisotropic expansion. This is nominal, synthetic numerical evidence. The
upstream law must supply the integrated strain and admit its temperature path;
no coefficient integration, frame calibration, portable transform codec or
uncertainty propagation is added here.

`try_from_resolved_strain` connects the complete six-component
`fs-material::StrainTensorStatePoint` (in its `state_point` module) to this same
operator ingress. The selected card/claim/query-point bundle identity becomes
the upstream thermal-strain identity; callers still declare the strain to be
stress-free and integrated over the supplied temperature path. The G1/G3 test
now loads every convention through a portable material card before checking
the same independent stress, nodal force, displacement and geometry oracles.

`try_from_joint_strain` additionally retains admitted source covariance or its
exact unknown-correlation reason. With a fixed proper rotation Q and fixed
admitted stiffness C, it computes `Sigma_e = A Sigma_source A^T` and
`Sigma_stress = C Sigma_e C^T`, where A includes source shear scaling, order
and the second-order frame transform. Both outputs use target Mandel order;
strain covariance is dimensionless squared and free-stress covariance is Pa
squared. No off-diagonal covariance is discarded. The fixed 6x6 calculations
publish symmetric results or refuse nonfinite values, negative variances and
complete loss of a nonzero covariance. No eigenvalue clipping manufactures
positive semidefiniteness.

`JointThermalStrainTransformReceipt` holds the executable nominal receipt,
conditional covariance outcome and source joint-state identity. Domain
`joint-thermal-strain-transform.v1` binds all of them and the output covariance
bits. G3 compares all six conventions and both source policies with a correlated
finite ensemble rotated by an independent full-tensor calculation, including
covariance of actual thermal nodal loads. G0 retains NoBlock, PartialMembership
and UnstatedMarginal outcomes and checks arithmetic loss/overflow. These are
conditional linear second moments, not confidence intervals or guarantees that
all realizations remain in the admitted small-strain/material domain. Uncertain
stiffness, orientation, temperatures, their cross-correlations and nonlinear
uncertainty propagation remain unsupported; unknown covariance is never zero.

`recover_thermal_stress` consumes the accepted thermal displacement and its
exact mesh, thermal field, material coefficients and constraints. It computes
the actual P1 elastic strain `B u - epsilon_free`, Cauchy stress
`C (B u - epsilon_free)`, stored elastic energy `V/2 epsilon_elastic:sigma`,
and full nodal internal force `integral B^T sigma dV`. These forces are support
reactions on fixed DOFs and equilibrium residuals on free DOFs; residual entries
are retained. Free expansion produces zero stress up to numerical error, while
constrained expansion produces the corresponding physical stress and energy.
The caller names the reference mesh frame and supplies a finite positive
ceiling applied to both total and elastic Mandel strains. This is a declared
small-strain applicability bound, not a yield or material-validation certificate.

Recovery rechecks the original thermal-load identity and the actual density
and stiffness bits retained under `thermal-response-input.v1`; changing a law
while reusing its caller-supplied material ID therefore refuses. Existing v1
load and displacement identity formats are unchanged. The immutable
`TetThermalStressReport` binds those inputs, solution, frame, strain ceiling,
element outputs, nodal forces and energy under `thermal-stress-recovery.v1`.
Cancellation is checked before allocation and at every element boundary; no
partial field is published. Existing mesh/element budgets bound report storage.

`stress_in_frame` expresses an element's recovered stress by `Q sigma Q^T`
with an explicit `StressTensorBasis`. Physical and engineering stress shear
are both `sigma_ij`; Mandel shear is `sqrt(2) sigma_ij`. Engineering-strain
factor two must never be applied to stress. Both existing component orders
are supported. Proper rotation and nonzero frames are required; a frame mapped
to itself requires identity rotation. Nonfinite and complete nonzero-to-zero
transforms refuse. `StressTensorTransformReceipt` retains source stress,
report/element identity, both frames, target convention/order, rotation and
output, binding them under `stress-frame-transform.v1`.

G1/G3 checks free expansion and two supported bodies with distinct thermal
strains/stiffness against closed-form displacement, stress, reaction and
energy formulas. All six output conventions are compared with independent
full-tensor rotation and stress/strain work conjugacy. G0/G4 checks reused
material IDs with changed stiffness, changed thermal fields/constraints,
strain bounds, invalid frames and pre-cancellation. This is nominal linear
elastic recovery in the reference geometry, not stress-source pack import,
finite-deformation/Piola stress, yield/buckling admission, evolved-state
recoupling, propagated solved-stress uncertainty or measured validation.

`JointStressTensorTransformReceipt::try_from_resolved` consumes the immutable
joint stress state from fs-material. It transforms source stress and full
covariance through stress-specific order/shear normalization and a fixed proper
rotation into canonical target Mandel coordinates. Tensor and engineering stress
both store physical shear; the transform never applies engineering strain's
factor two. The existing bounded 6x6 covariance kernel computes `A Sigma A^T`
and preserves output symmetry without clipping or assuming independence.
Nonfinite normalization/rotation/covariance, negative variance and complete
nonzero-to-zero numeric loss refuse. Nonzero target frames and identity rotation
for a frame mapped to itself are required; identity rotation retains normalized
nominal bits. The work is one bounded tile and allocates no growing buffers.

The receipt retains source joint identity, original coordinates and nominal
stress, target frame, rotation and both numerical outputs under
`joint-stress-frame-transform.v1`. `StressTensorUncertainty` holds the complete
target covariance [Pa^2] or the unchanged source unknown reason. G3 compares
all six source conventions and both joint policies against a full-rank
twelve-realization ensemble transformed as full tensors and passed through
actual clamped thermal displacement/stress recovery and checks receipt tampering.
G0 retains NoBlock, PartialMembership and UnstatedMarginal and checks frame
errors, shear overflow and complete subnormal loss. These are synthetic numerical tests.
This is source-data covariance conditional on fixed orientation, not solved-field
uncertainty, a confidence/domain bound, material/frame calibration, propagation
of uncertain temperature/stiffness/orientation, or experimental validation.

## Material replacement on an accepted P1 trajectory

`loads3::transient::events::rebind_material_candidate` permits source-data
changes to density, stiffness and stress-free strain on the same reference
mesh and fixed-DOF set. It retains the accepted physical displacement, clock
and step count, solves `M_after_ff v_after = M_before_ff v_before`, and uses
the original initializer to rebuild operators and re-equilibrate acceleration
under the explicitly supplied new load/thermal endpoint. Subsequent ordinary
Newmark steps and stress recovery consume the actual retained new problem.
The event does not consume a time-advancement slot.

The conserved quantity is FREE canonical momentum. Stationary constraints
may exchange full nodal or total physical momentum; this operation reports
no support impulse. Moving-support trajectories and changed meshes/restraints
refuse. This is externally supplied material reconfiguration, not mass
transport or a continuous material transition.

Required signed work is `Delta(1/2 v'Mv) + Delta(sum V/2 e'C e)`, where
`e = strain - strain_free`. Kinetic and physical shifted-elastic changes
remain separate. Positive work must be supplied by the reconfiguration owner,
and negative work returned to it; neither is automatically caloric heat or
dissipation. An explicit absolute work budget and independent SI momentum
residual gate admit the candidate. Per-element old/new material and thermal
identities and reference masses remain in its report. Failure of source,
solve, momentum, strain, work or cancellation checks preserves the old run.

`tests/transient_events3.rs` exercises immutable source-card substitution
through the native solver: doubled density, tripled stiffness, changed
eigenstrain, exact scalar momentum/acceleration/work, subsequent motion,
zero-work equivalent replacement, retained nonzero clock/step history,
fully clamped storage, and late work/strain or cancellation refusal.
This does not differentiate the replacement map or validate a real material
transition.

Material replacement retains its kinetic increment through the original sparse
mass action and the actual old/new mass coefficients before combining external
work or applying its work bound. An unchanged large kinetic inventory cannot
erase a smaller density-reconfiguration exchange. Native impulse admission
likewise retains actual endpoint velocity differences through the original full
consistent mass, including restrained cross rows; that same mass action supplies
momentum and kinetic storage. Independent applied and signed support work retain
endpoint means and support-impulse low words through the work gate. Reports
expose retained kinetic/source/support or reconfiguration work for composition.
These numerical work defects remain separate from physical loss and heat.

## Instantaneous impulses on an accepted P1 trajectory

`loads3::transient::TetTransientRun3::apply_nodal_impulse` applies
explicit nodal impulses through the same sparse consistent reference mass as
Newmark stepping. It solves `M_ff delta_v_f = J_f`, retains supported velocity,
and recovers full support impulses as `M delta_v - J` on restrained rows.
Momentum gates use the actual representable velocity change. Moving supports
retain their impulse work; `average(v) dot M delta_v` is checked against applied
plus support work without subtracting large endpoint kinetic energies.

The event preserves displacement, time, step count, material, thermal state
and finite loads. Native acceleration and finite-force reactions are refreshed
for the new velocity and damping force before transactional publication.
It does not advance a strain/history law or generate damping heat at zero time.
Collision timing, restitution and the opposite partner remain with the caller.
`tests/transient_impulse3.rs` supplies direct consistent-mass, support-work,
continued-motion and refused-candidate regressions. These are local numerical
checks, not calibration of an impact law.

## Thermal material evolution in sparse P1 dynamics

`loads3::transient::thermal::advance_thermal_candidate` refreshes the actual
source-bound thermal strain, material identity and global Mandel stiffness.
Changed tensors rebuild the native endpoint stiffness using the existing P1
assembly; unchanged tensors reuse the original sparse operators. Both paths
preserve the accepted Newmark displacement, velocity, acceleration and clock.
Reference density, sparse mass, reference inventory, mesh connectivity and
stationary zero restraints stay fixed. Refreshed source counts, identities,
geometry, quality and strain limits pass the original assembly and response
checks. A fully clamped trajectory uses the existing zero-free-DOF dynamic
assembly and still reports stress and stored energy.

The ordinary Newmark audit retains its unshifted quadratic-energy convention,
averages old/new thermal nodal forcing, and records signed stiffness parameter
work `W_K = q_old' * (K_next-K_old) * q_next / 2`.
`ThermalTransientStep3` separately audits physical kinetic plus shifted elastic
energy, with `e = epsilon-epsilon_free`. Its independent element work is
`W_T = -sum V * (sigma_old+sigma_next)/2 : delta(epsilon_free)` and
`W_C = sum V/2 * e_old : (C_next-C_old) : e_next`. Both terms are needed when a
fully restrained solid expands and softens with zero displacement work.
Material parameter work is signed prescribed input; it is not dissipated heat.

Changing stiffness requires zero stiffness-proportional Rayleigh damping;
fixed reference mass damping remains supported. Changing endpoint Rayleigh
stiffness damping is not unconditionally dissipative and is refused.
Unchanged-stiffness trajectories retain their existing damping behavior.
The native and physical audits use the declared run tolerances and refreshed
source identities. These are exact discrete endpoint work identities;
continuum accuracy still requires resolving the prescribed material path and
checking timestep refinement. Failed assembly, motion, stress, energy or
cancellation checks leave the original trajectory unchanged.

G1/G3/G4 tests in `tests/transient3.rs` cover an analytical thermally driven
oscillator with timestep refinement, source-identity refresh, scalar softening
and stiffening with retained motion and mass damping, fully clamped thermal
storage with distinct expansion/stiffness work, changed-density and variable
stiffness-damping refusal, late strain refusal and cancellation. This is
one-way prescribed thermal/material loading. Reciprocal caloric feedback,
changing reference mass, moving thermal supports, an adjoint of the new
material inputs and measured validation remain outside this path.

## Prescribed support motion in sparse P1 dynamics

`loads3::transient::motion` extends the existing sparse average-acceleration
Newmark trajectory with absolute world-component support displacement. The
fixed DOF set, reference mesh, materials, thermal eigenstrain and Rayleigh
damping remain frozen. `begin_transient_with_supports` accepts full initial
displacement/velocity fields and explicit constrained acceleration; free
acceleration entries in that input must be zero because they are solved.
Initialization and every step include the complete partitioned load
`f_f + fT_f - M_fc a_c - C_fc v_c - K_fc u_c`, with
`C = eta M + zeta K`. Consistent-mass coupling is retained.

`advance_with_supports` accepts sorted, unique fixed-DOF displacement values.
Omitted fixed components mean zero at that endpoint. Support acceleration
and velocity follow the same beta=1/4, gamma=1/2 recurrence as the free state;
they are not finite differences with an unrelated integration convention.
Recovering acceleration from rounded displacement samples amplifies their
roundoff by the inverse squared step size; a smaller step is not an unlimited
accuracy guarantee for prescribed derivatives.
The displacement-controlled recurrence subtracts the absolute endpoints before
the velocity/acceleration predictor drift, retaining resolved drift even when
adding it to a nonzero displacement would round it away.

`advance_with_support_accelerations` instead accepts sorted, unique fixed-DOF
endpoint accelerations in m/s². It integrates support velocity and displacement
with the same Newmark map while retaining the supplied acceleration directly.
An omitted component has zero endpoint acceleration and continues its existing
motion. It does not reset displacement or stop a moving support. This route
avoids inverse-squared-step amplification of sampled displacement roundoff and
retains consistent-mass/damping/stiffness coupling into free acceleration.
Support work uses the accepted Newmark increment `h * (v_old + v_new)/2`, so
resolved kinetic gain still has reaction work when absolute displacement
cannot represent a tiny increment. The original admission and energy gates
remain in force; inputs describe absolute world Cartesian acceleration,
not fictitious forces in a rotating frame or impulsive velocity changes.

Full element actions reuse the canonical P1 stiffness/mass kernel. No dense
global operator is assembled. Fully prescribed meshes use empty reduced
matrices and zero solve iterations after full geometry/material validation;
public modal assembly still refuses `NoFreeDofs`.

`nodal_motion` includes the accepted support u/v/a. State energy includes all
nodes, and `TransientStep3::support_work_j` is the trapezoidal work of support
reactions **on the solid**. Mechanical work includes external loads on both
free and prescribed components; the frozen full thermal equivalent load has
its own work term. The acceptance identity is
`delta(kinetic + quadratic_strain_energy) + damping_loss = mechanical_work +
thermal_equivalent_work + support_work`. Quadratic strain energy remains
`u'Ku/2`; `recover_response` separately reports the actual thermally shifted
elastic energy and stress, full inertia/damping and support reactions.

Input admission, full nodal displacement limits, element limits, sparse solve
limits, per-node/element cancellation checkpoints and the energy/residual
gates precede atomic publication. A clone retains the prescribed kinematics,
reaction endpoint and immutable full thermal load as well as the free state.
Force-only advancement and reduced-only state adjoints
explicitly refuse this lane; they cannot omit
support terms while claiming a valid trajectory or gradient. Displacement
impulses, changing restraint axes/topology, adaptive event discovery,
temporal/continuum error certification and measured validation are not supplied.
Prescribed motion requires a physically resolved time-step ladder.

`loads3::history::dynamics` connects typed `SupportMotion3` displacement
histories to this same inertial trajectory. `begin_transient_with_support_history`
requires the initial restrained displacement to match the retained history
sample, with explicit initial velocity and constrained acceleration.
`advance_history_with_supports` samples all loads and support values together
at the run's next physical time before delegating one coupled Newmark step.
Value-table clamping or wrapping does not replace physical time for moving
frames, and rotating a payload does not rotate the fixed world restraint axes.
Malformed bindings, sampling failures and cancellation leave the accepted
trajectory intact. These adapters do not infer analytic derivatives, impulses
or an event schedule from displacement tables.

`SupportAcceleration3` supplies the parallel m/s² input convention through
`at_with_support_accelerations`,
`begin_transient_with_support_acceleration_history`, and
`advance_history_with_support_accelerations`. The same bounded, unique fixed-DOF
binding and typed world-vector sampler serves both conventions. Initial q/v
remain independent physical states; initial constrained acceleration is sampled
from the record and free acceleration is solved from the original full coupled
equation. Omitted fixed components integrate zero acceleration, retaining their
existing motion. Source values and all external forces share one physical clock.

`EventLoadHistory3::new_with_support_accelerations` reuses the finite load/frame
knot index and exact left/right samplers. Its acceleration initializer and
`advance_to` integrate the left acceleration through each interval, then apply
the right acceleration and external force together through
`apply_support_acceleration_jump`. That zero-duration operation retains q/v,
stored energy, clock and step count bit-for-bit; the original mass solve resets
free acceleration with the consistent M_fc term and refreshes support reactions.
There is no impulse, backward-smearing of a pulse edge, or event work. A late
right-side refusal rolls back the entire private interval. Operators are shared,
not reassembled, and the original adaptive comparison includes integrated fixed
coordinates even for fully prescribed bodies. Smooth records still need temporal
refinement; periodic records need the existing explicit finite unrolling. Source,
frame and event-time parameter gradients are not supplied by this adapter.

The same event owner retains borrowed native chirp laws during bounded setup.
Before fixed, adaptive or recorded advancement, their shared scenario limiter
restricts duration times the maximum instantaneous source frequency to 1/16
cycle. Direct pressure histories and translating/tilting frames referenced by
actual surface, volume or support vectors participate, including their ancestors
and continuing tones outside a sweep. Unused frames retain their authored knots
without imposing dynamic frequency limits. Each used chirp tilt also scales its
cycle allowance by `1/max(1, abs(amplitude_rad))`, bounding that individual tilt's
angular travel to `2*pi/16` radians per step. Frame reachability and the scaled
allowances are retained during setup; source admission does not allocate or
rescan the full histories per step. Cancellation and unrepresentable clock advances
refuse before publishing motion. Adaptive coarse/fine comparisons therefore
cannot silently skip complete cycles of an authored tone, while their original
motion-error, work, solve and lifetime-step gates still apply. This controls
the authored sources, not every frequency generated by nonlinear composition;
other smooth signal families and solution accuracy still require refinement.
The native `transient_events3` regressions check the independent center-of-mass
response to whole-cycle pressure forcing and frame-driven clock limits.

`EventLoadHistory3::new_with_supports` explicitly indexes those original
support, load and frame knots together under one schedule budget. Its support
initializer requires the same full initial kinematics. Event stepping uses
the exact left force before a knot and the right force after it, retaining
the native support Newmark recurrence. Finite force jumps keep support u/v/a,
include the original inertial/damping/stiffness lift in free acceleration,
and refresh reactions for subsequent support-work accounting; they add no
impulse or energy. A support-displacement jump, including one caused by a
frame jump, refuses without advancing accepted motion. The existing adaptive
two-half-step comparison now includes constrained-node kinematics as well as
free motion, including fully prescribed bodies. Support-driven runs compare
raw full/half disagreement for all kinematics: prescribed derivative errors
also enter free acceleration through the consistent mass. No Richardson order
is assumed for that recurrence, which can retain alternating acceleration
error. This local consistency indicator is not a derivative-error bound;
neither unknown events, exact support derivatives nor a global temporal bound
are inferred. Accepted work, clock and fork/restart semantics
remain those of the original trajectory.

`tests/transient_supports3.rs` supplies G1 exact mass/stiffness/both-Rayleigh
cross-term initialization, a closed-form accelerating-base oscillator with
second-order refinement, fully prescribed thermal/external/support work and
reaction checks, and G3/G4 fork, malformed-input and cancellation checks.
It also checks direct acceleration forcing against an independent scalar
partitioned equation with both Rayleigh cross terms, unchanged constrained
acceleration and exact kinetic/support work below displacement resolution,
unchanged prescribed displacement with a resolved velocity reversal, and
transactional malformed-input/cancellation handling.
`tests/transient_support_history3.rs` covers rigid translation with damping,
agreement with explicit core endpoints, physical frame time beyond a clamped
table, a base-excited oscillator with both held-force pulse edges, full-support
adaptive rejection, and transactional sampling, binding and cancellation refusals.
`tests/transient_support_acceleration_history3.rs` adds G1/G3/G4 checks for
finite acceleration switches with both Rayleigh terms, rigid pulse motion/work,
independent coupled base/body-acceleration convergence, physical frame time and
frame-tilt events, all-prescribed adaptive rejection, and late-event rollback.

### Prescribed-support motion adjoints

`loads3::transient::adjoint::support` supplies a separate full-state transpose
companion through `prescribed_adjoint`. `NodalMotionCotangent3` includes q/v/a
on every node, so direct objective dependence on constrained motion is retained.
`pullback_step` returns previous full-state cotangents, endpoint free-force
cotangents, and absolute support-displacement control cotangents in sorted
`fixed_dofs()` order. An omitted primal fixed row is still a zero displacement
control; varying it requires explicitly adding that row at the endpoint.
`pullback_initialization` reverses the actual initial coupled mass solve and
returns initial full q/v and constrained-a derivatives. Free initial a was
solved and its unused input cotangent is exactly zero.

The reduced solve reuses the original sparse transpose owner. A bounded
canonical element sweep supplies actual `M_fc^T lambda` and `K_fc^T lambda`
without assuming bitwise matrix symmetry. For `b=h²/4`, `d=h/2`, the support
endpoint seeds are `q_bar-K_fc^T lambda`,
`v_bar-(eta M_fc^T+zeta K_fc^T)lambda`, and `a_bar-M_fc^T lambda`.
The pullback then reverses the original prescribed Newmark recurrence, including
the previous support acceleration. Fully prescribed meshes require no sparse
solve but retain the kinematic transpose and direct objective seeds.

Full-DOF and element caps precede scratch allocation; independent sparse
transpose/solve limits remain explicit. Failed or cancelled pullbacks publish
no partial gradient and leave both the companion and primal run reusable. These
are derivatives of fixed-step discrete motion equations, not of PCG iterations,
acceptance, floating-point rounding, moving constraint axes, impulses, or step
times. Direct stress/reaction/work/energy-observation and material derivatives
are not supplied by this motion pullback. Prescribed derivatives retain the
forward recurrence's inverse-squared-step conditioning.

`pullback_acceleration_step` is the separate full-state transpose for direct
support acceleration controls. It returns previous full q/v/a, free-force and
endpoint prescribed-acceleration cotangents while retaining the same reduced
transpose solve and canonical cross-mass/stiffness actions. Reversing the
integrated support kinematics uses only `h` and `h²`, never their inverses.
The initial coupled acceleration pullback chains unchanged. The support-adjoint
tests check analytical gradients down to a 1e-8 s step and all ten initial,
force and acceleration controls against actual two-step trajectory perturbations.

`tests/transient_support_adjoint3.rs` checks the analytical one-free-DOF
mass/stiffness/both-Rayleigh transpose, every initial/force/support control
group through three native steps against centered differences with full nodal
objectives, the fully prescribed kinematic limit, fixed-index ordering and
malformed-input, budget and cancellation refusals.

`material_scale_pullback` now accepts prescribed-support runs at their exact
accepted endpoints, using actual full nodal q/v/a. For the corresponding
step/initialization free-force cotangent lambda, it contracts
`-lambda' d(M*a + (eta*M + zeta*K)*v + K*q - fT)` with zero constrained
equation multipliers. Canonical P1 element matrices and Mandel B are reused;
the complete `M_fc*a_c`, `C_fc*v_c` and `K_fc*q_c` terms are retained, including
the explicitly supplied initial support acceleration. Returned independent
per-element stiffness/density scale derivatives and global Rayleigh eta/zeta
derivatives compose with the existing full-motion pullback. A common material
scale sums its element entries; fully prescribed motion has zero implicit
material sensitivity.

The free thermal strain is frozen, but the stiffness derivative includes its
equivalent-force dependence. Do not count that stiffness-dependent thermal
load a second time. Geometry, support axes/authored kinematics and external
nodal forces are fixed. Parameter-dependent external forces, initial/support
controls and direct objective terms must be chained separately. This is not a
tensor-entry, thermal-law, shape, reaction or stress-observation derivative.
Independent full-DOF and element bounds precede scratch allocation; cancellation
or refusal returns no partial result and never changes accepted motion. The
reduced-state-only internal helper continues to refuse driven supports because
reduced q/v/a cannot recover their physical cross terms.

`tests/transient_support_material3.rs` checks the analytical initial coupled
mass/stiffness/damping terms, a three-step anisotropic heterogeneous thermal
trajectory against six independent material/damping finite-difference
directions at two parameter points, and the fully prescribed zero-gradient,
budget, cancellation and retry limits.

`adjoint::support::reaction` supplies physical support-force observation
pullbacks through `reaction_pullback`. The observed force is the force ON the
solid, `R_c=(M*a+(eta*M+zeta*K)*v+K*q-fT-f_external)_c`, in fixed world
Cartesian components. A caller supplies finite objective derivatives with
respect to sorted unique fixed DOFs, typically from a recovered reaction
residual. The result includes full motion seeds `(K' w,C' w,M' w)`, direct
fixed-force rows `-w`, and direct material/Rayleigh terms `+w' dR` at fixed
motion. The full canonical transpose and material contraction owners are
reused; no global full matrix, new equilibrium solve or finite difference is
introduced. `reaction_motion_pullback` on the prepared prescribed companion
supplies the same motion seed without a primal-state tape.

For a complete reaction objective, chain the motion seed through the existing
initialization/step/trajectory reverse, then add the direct material and
fixed-force terms. Fully prescribed motion has no implicit material derivative,
but its physical reactions retain nonzero direct thermal, inertia and damping
derivatives. The same `-w` applies independently to a full frozen thermal
equivalent force. Stiffness-dependent thermal effects already present in the
material term must not be counted again. The endpoint API supports stationary
and driven restraints, keeps full consistent-mass and damping coupling, and
checks independent full-DOF/element allowances before scratch allocation.
Malformed rows, nonfinite arithmetic, refusal and cancellation publish no
partial result and never change accepted motion. No moving-axis, shape,
thermal-law, stress-observation, physical-work or acceptance-gate derivative
is inferred.

`tests/transient_reaction_adjoint3.rs` checks analytic direct and implicit
initial reaction terms, the same observation seed on stationary and driven
supports, three-step reaction-target finite differences in seven independent
material/damping/support/free-force/fixed-force directions, fully prescribed
nonzero direct material sensitivities, and malformed/budget/cancellation retry
boundaries.

`adjoint::support::stress` adds pointwise elastic stress objectives through
`stress_pullback`. Rows contain sorted unique tetrahedron IDs and finite
`dJ/dsigma` in mesh-frame Mandel `[xx,yy,zz,xy,yz,zx]`; the final three stress
entries are `sqrt(2)` times physical tensor shears. The canonical P1 strain
operator and admitted constitutive tensor supply `B' C' w`, including fixed
DOFs, and the direct per-element stiffness-scale term `w' sigma`, where
`sigma=C*(B*q-epsilon_free)`. Velocity, acceleration, direct density and direct
Rayleigh terms are zero. No element-volume or time-quadrature factor is added;
the objective owner includes any such weights. For complete design gradients,
chain the full displacement seed through the matching accepted-state reverse
and add its implicit material terms to the direct stiffness contribution.

Geometry, material axes and free thermal strain remain fixed. Stationary,
driven and fully prescribed endpoints use their actual full displacements.
Full-DOF and element caps apply before scratch allocation; both total and
elastic Mandel strains must obey the same finite positive ceiling as stress
recovery, even for unobserved elements. Invalid rows, nonfinite arithmetic and
cancellation yield no partial gradient and preserve accepted motion.
`von_mises_limit` supplies the smooth dimensionless inequality
`(sigma_vm/limit)^2-1 <= 0` and its Mandel stress seed, with a defined zero
derivative at hydrostatic stress and a caller-supplied positive finite limit.
It does not infer a yield law, failure probability or physical admission.
`tests/transient_stress_adjoint3.rs` checks analytic thermal/shear terms without
volume weighting, anisotropic nodal and material finite differences at two
mesh scales, von-Mises derivatives, and refusal/cancellation/retry boundaries.

`adjoint::support::trajectory` records complete prescribed-support simulations
through `record_prescribed_controls`. For controls alone, its owning recorder
retains only accepted transition times and control kinds: the frozen linear
motion Jacobian needs no full primal tape or finite-difference replay.
Each successful original coupled step is recorded
atomically; failed sampling, support admission, solves, cancellation or storage
admission preserves both the accepted state and schedule. Cloning forks both.
The typed-history adapter samples one retained force/support endpoint before
advancing, using the same physical time and frame binding as the original run.

`reverse_with_observations` adds sparse full-nodal objective derivatives at any
accepted endpoint, including initialization, and reverses the complete schedule
once. It returns every initial/endpoint force row, each prescribed displacement
or acceleration control row, independent initial free q/all v/fixed a, and the sum of force
cotangents for the one frozen thermal equivalent load. Initial fixed q appears
only in the initial displacement-control row, avoiding duplicate control
contributions. Sparse observations preserve caller summation order with
O(steps + observations) indexing. Explicit step, returned-gradient and observation
caps precede allocation. Refused/cancelled reversals return no partial gradient.

`enable_material_recording` additionally retains the original full nodal q/v/a
at initialization and each accepted endpoint. It must be enabled once, before
the first finite step. Independent element/full-DOF limits and a checked
`samples * nodes * 9` scalar cap precede recording; failed advancement preserves
both the accepted motion and material tape. Cloning forks the retained tape.
The typed-history adapter uses the same bounded atomic recording path.

`reverse_design_with_observations` returns all existing control derivatives and
per-element stiffness/density plus global Rayleigh damping derivatives in one
reverse sweep. It contracts each residual multiplier against the exact retained
full motion, including the solved initial acceleration. No additional transpose
solve, replay or finite difference is needed. Material derivatives hold geometry,
restraint axes, authored support kinematics, free thermal strain and independent
external nodal loads fixed. The stiffness-induced thermal equivalent force is
already included; its dependence must not be chained again through the separate
frozen-material thermal-force cotangent. Parameter-dependent external loads and
direct objective terms remain caller-owned chains. Refusal/cancellation returns
no partial combined gradient and leaves the tape reusable.

`reverse_design_with_stress` adds sparse element-stress observations at any
recorded endpoint to those motion objectives. `ElementStressObservation3`
contains objective/Pa seeds in mesh-frame Mandel `[xx,yy,zz,xy,yz,zx]`; duplicate
element/state rows add in caller order, and callers include quadrature weights.
The same reverse sweep accumulates full `B' C' w` displacement seeds, direct
`w' sigma` stiffness terms, and implicit trajectory material/load/support terms,
including initialization. No additional equilibrium solves or stress tape are
needed. The combined motion/stress row cap and existing element/full-motion
budgets apply. Every state with stress rows obeys the finite positive total and
elastic strain ceiling, even when its seeds cancel to zero. Thermal strain,
geometry, frames, restraint axes and accepted times remain fixed. The smooth
`von_mises_limit` seed composes with this path for strength objectives; it does
not add plasticity or infer material failure limits. The three regressions in
`tests/transient_stress_trajectory3.rs` cover six material/damping/support/load
directions for mixed stress/strength/motion objectives, fully prescribed
intermediate-state controls, and observation-budget/strain/cancellation retries.

Callers supply actual residual/quadrature weights for their motion objectives
and chain returned world-component controls through authored waveforms. This
does not infer derivatives of frame transforms, knot/step times, interpolation,
geometry, direct control costs, physical work reports or acceptance.
`tests/transient_support_trajectory3.rs` supplies a multi-step least-squares
objective checked against all twelve independent initial/force/support control
directions, an analytical typed fully-prescribed ramp-amplitude chain rule,
and bounded recording, cloning, cancellation and refusal/retry checks.
`tests/transient_support_design3.rs` adds a thermal multi-step least-squares
objective observed at every endpoint, compared with six material/damping/support/
force directions at two parameter points, plus opt-in recording, storage refusal,
atomic failure, cancellation/retry, clone identity and fully prescribed limits.

The same recorder accepts `advance_with_support_accelerations`,
`apply_support_acceleration_jump` and atomic
`advance_to_support_acceleration_jump`. Displacement and acceleration controls
may alternate; each transition reverses its original kinematic and coupled
equilibrium owner. A finite jump preserves q/v, replaces acceleration through
the consistent mass solve, and records distinct left/right states at one time.
`pullback_acceleration_jump` discards all previous-a cotangents while retaining
the actual cross-mass, damping and stiffness terms. Acceleration controls have
separate `accelerations` rows; jump controls and forces carry `RightLimit`,
while ordinary endpoints carry `IntervalEndpoint`. Initial fixed acceleration
remains in `initial.a`, and solved free initial acceleration has no input gradient.

The typed smooth-history adapter and
`EventLoadHistory3::advance_recorded_accelerations_to` sample the original
frame-bound load/support histories, including exact left/right event limits.
The event adapter forks only the primal candidate; it never copies the existing
trajectory tape. Both transitions and their optional material motion are admitted
before the interval, and a late right-limit failure publishes neither. Recording
limits count transitions, including zero-duration jumps. Sparse motion, reaction
and stress observations can address either side independently; material and
Rayleigh gradients use the corresponding full accepted motion in the same reverse
sweep. Authored amplitude chains remain caller-owned, with frames and event times
fixed. `tests/transient_support_acceleration_trajectory3.rs` checks all eighteen
independent control/material/damping directions with reaction and stress
objectives, an analytic typed acceleration pulse, smooth-history parity, and
cancellation, budget and late-jump rollback.

## Public types and semantics

- `Mesh2` / `Patch`: structured body-fitted meshes — P1 triangles, Q1
  quads, mapped Q1 panels (`cooks_membrane`), with the four sides as
  named patches. `shapes_at`/`quad_points` expose the isoparametric
  machinery (chain rule with the transposed inverse Jacobian — the
  conformance battery pinned this exact slip on triangles and mapped
  panels).
- `LinearProblem` (plane strain/stress via `PlaneKind`, `lame`):
  small-strain elasticity with strong Dirichlet and dead-load traction
  patches; `Formulation::Standard` vs `Formulation::BBar` (Hughes'
  dilatation projection — element-average volumetric gradient,
  equivalent to condensed constant-pressure mixed elements).
- `HyperProblem` + `NewtonSettings`/`NewtonReport`: plane-strain
  finite strain; the 2D displacement gradient embeds into the 3D
  deformation gradient (F₃₃ = 1); residual `∫ P : ∇δu` and consistent
  tangent `∫ ∇δu : A : ∇Δu` come from fs-material's exact AD Piola
  stress and 9×9 tangent. Newton globalizes by Armijo backtracking on
  the potential energy (material refusals — det F ≤ 0 — halve the
  step) under linear load stepping; systems solve with MINRES
  (symmetric, indefinite-safe near buckling). `residual_and_tangent` /
  `potential_energy` are the public consistency probes.
- `CutElasticity`/`CutSolution`: source-compatible legacy facade over
  `fs_cutfem::CutElasticity` and its canonical solution type. The facade
  constructs `IsotropicElastic::new(E, ν, 1.0)`; the fixed strain limit
  supplies compatibility metadata only and is not evidence that a solved
  field remains inside that domain. Historical stabilization knobs retain
  their physical coefficients by translating
  `β_canonical = β_legacy(λ+2μ)/μ` and
  `γ_canonical = γ_legacy(λ+2μ)/μ` before delegation. Canonical assembly,
  graded-tree constraints, Nitsche and ghost terms, deterministic CG,
  retained topology, and error integration are not reimplemented here.
  The facade fixes the legacy solver controls at tolerance `1e-12` and
  `60_000` iterations. The removed implementation's post-hoc
  `1e-8` acceptance after targeting `1e-12` is not reproduced: a miss
  of the canonical `1e-12` gate fails closed. `BoundaryTraction`,
  `DesignBoxEdge`, and opaque checked `EdgeBand` are direct fs-cutfem
  re-exports. `solve_with_boundary_traction` delegates the typed support
  descriptor without inspecting or reimplementing it; the legacy callback
  field must be `None`, and an SDF crossing through supported load still
  refuses canonically.
  Physical-coefficient equivalence is claimed only when each active term's
  dimensionless translated coefficient is finite. Canonical validation refuses
  an overflowing Nitsche translation when `traction_free_interface` is false
  and any overflowing ghost translation. When `traction_free_interface` is
  true, the unused translated Nitsche value is deliberately ignored and
  carries no physical-coefficient claim. A future physical-coefficient input
  surface is required to represent active terms without a non-representable
  intermediate.
- `select_formulation(RegimeIndicators)`: the documented
  element-selection policy — B-bar when ν ≥ 0.45 or slenderness ≥ 5
  (fs-regime's indicators feed this; the locking battery measures the
  boundary).
- `stability` (bead tfz.15): `reduced_pencil` (free-DOF K and
  geometric stiffness K_G(σ₀) from the linear prebuckling state),
  `buckling_loads` — the pencil `(K + λK_G)φ = 0` solved MATRIX-FREE
  by fs-la LOBPCG on the Cholesky-reduced operator L⁻¹(−K_G)L⁻ᵀ
  (dense reduction fixture-gated at n ≤ 4096), modes −K_G-normalized.
  Both matrices are independently scaled before reduction; positive inverse
  loads have no dimensioned absolute cutoff. Transformed and original-pencil
  relative residuals must pass 1e-8 before publication, with original residuals
  retained in the result. Exhausted iterations refuse. Small pencils (n ≤ 64)
  whose requested block cannot fit LOBPCG use the native dense symmetric solver;
  invalid shapes, maps, nonfinite coefficients and asymmetry refuse explicitly.
  The reduced DOF map is component-independent. `buckling_loads_cancellable`
  polls inside the native Cholesky factorization and between bounded, resumable
  LOBPCG advances; resolved blocks stop before exhausting their iteration cap.
  `stability::solid3` exposes `proportional_buckling` on accepted P1 and
  isoparametric prestress assemblies (TET10/HEX8/HEX20/HEX27), reusing their
  original elastic and geometric matrices, supports and stress quadrature.
  It solves `K_e + lambda K_g`, not the already prestressed `K_t + lambda K_g`;
  its load factor scales the complete reference stress without counting it
  twice. Results retain the reference assembly identity and original xyz mode
  mapping, including fully supported nodes. Geometry, mass and constitutive
  tangent remain frozen. The proportional entry point scales the complete
  reference and does not admit a scaled material state.
  Accepted P1 and isoparametric assemblies also expose
  `base.incremental_buckling(&reference, ...)`: lambda=0 denotes the fixed base
  equilibrium and lambda=1 denotes the accepted unit-increment reference. The
  pencil is `K_e + K_g,base + lambda*(K_g,reference - K_g,base)`. Pairing binds
  original coordinates/connectivity/frame, actual original/incremental material
  tensors and source identities, homogeneous incremental supports, and native
  isoparametric formulation/quadrature. Thermal strain, dead loads and prescribed
  values may differ. Canonical CSR subtraction is linear in nonzeros. The complete
  base tangent must pass native SPD admission even for a zero increment; an
  unstable or singular base refuses. Results retain both accepted state
  identities and original xyz modes. These frozen-reference linear screens do
  not predict changing material laws, follower loads or postbuckling.
  `tests/buckling3.rs` supplies analytical constrained-tet and Standard/B-bar
  HEX8 references; `tests/buckling3_preload.rs` independently checks fixed thermal
  preloads plus prescribed increments, source/configuration pairing and unstable
  bases. `tests/buckling_scaling.rs` exercises numerical cancellation and odd
  three-component maps.
  `tests/buckling_scaling.rs` covers load/unit invariance, geometric work, small
  and full blocks, known noncommuting pencils and convergence refusal;
  `lambda_indicator` (Richardson) + `evidence_row`;
  `eigenvalue_derivative` (direct pencil derivative at frozen
  prebuckling stress), `group_stiffness` (per-group Young's-scale
  lever), and the clustered-eigenvalue policy: `ks_aggregate` /
  `ks_aggregate_derivative` — the conservative smooth KS lower
  envelope with softmax-weighted derivatives, because individual
  eigenvalues are not differentiable where branches cross (measured
  in the battery).
  `stability::sensitivity::element_scale_sensitivities[_cancellable]` computes
  total P1/Q1 element-modulus derivatives through the accepted reference
  equilibrium. `LinearProblem::solve_with_element_scales` and
  `reduced_pencil_with_element_scales` use the same finite positive per-element
  multipliers in elastic stiffness, prescribed-column lifting and recovered
  prestress; geometry, Poisson ratio and dead loads remain fixed. Direct elastic
  and geometric terms plus one native checked-CG adjoint per mode supply all
  element derivatives. The actual Standard/B-bar strain operator and its
  transpose retain full prescribed displacements and homogeneous incremental
  supports. No perturbed forward solve or third-order derivative tensor is
  needed. Returned scales, modes, gradients and true adjoint residuals remain
  paired; each adjoint must pass a 1e-8 physical relative-residual gate.
  One extra eigenpair checks the retained-spectrum boundary; relative gaps
  <=1e-6 refuse individual gradients rather than publishing a branch choice
  inside an unresolved cluster. The existing factorization/workspace caps
  remain. This is a fixed-geometry, fixed-load-parameterization sizing
  derivative, not a changing-support, follower-load or nonlinear stability
  derivative. `tests/buckling_sensitivity.rs` checks independently re-solved
  finite differences on heterogeneous triangles/distorted quads, plane stress/
  strain and Standard/B-bar, plus force/displacement homogeneity, mixed nonzero
  support lifting, cluster refusal and cancellation.
- `stability::sensitivity::aggregate` (bead tfz.15): conservative smooth KS
  total gradients through complete lowest-positive-mode packets. Interior
  repeated/crossing modes are admitted; a relative boundary gap <=1e-6 and
  missing requested positive modes refuse. The physical eigenvalues remain
  unchanged. One normalized native checked-CG adjoint sums the full material
  and prestress-redistribution pullbacks across the packet, including the
  Standard/B-bar transpose and original prescribed columns. Finite negative
  KS values are valid smoothing bias, not negative physical eigenvalues.
  Active interior clusters refuse sharpness that amplifies the native
  residual/roundoff resolution indicator above 1e-6 in log-weight ratios;
  this admission policy is not an eigenvalue forward-error certificate.
  Original eigenpair and adjoint residual gates, cancellation and fixed-load
  sizing scope remain. G1/G3 `tests/buckling_aggregate_sensitivity.rs` checks
  analytic repeated-column traces, independently re-solved total gradients,
  mixed/prescribed supports, near-incompressible B-bar, force-unit scaling,
  packet boundary/shortfall and excessive-sharpness refusal, and cancellation.
- `continuation` (bead tfz.15): `PathResidual` (generic equilibrium
  residual R(u,λ) = R_int − λF_ext; `HyperProblem` implements it for
  homogeneous clamps + ramped dead tractions), `advance` — Keller
  bordered pseudo-arclength with pivoted dense LU for total bordered dimension
  at most 1024 and matrix-free FGMRES above that cap. Both predictor and
  corrector solve the complete, row-equilibrated nonsymmetric border, so a
  simple fold does not require an inverse of the singular structural tangent.
  Corrector acceptance requires both physical force balance and the tangent-plane
  arclength constraint; a lost displacement increment cannot pass by returning
  to the previous equilibrium. Step reduction/growth respects the declared
  minimum/maximum, and failure preserves the complete pre-step checkpoint,
  limit-point events (λ̇ sign flips) and branch-point candidates
  (Cholesky definiteness flips without a λ̇ flip, gated n ≤ 1024);
  `PathState` — the checkpointable trajectory (clone to checkpoint,
  hand back to resume; split runs are BITWISE identical to straight
  runs); `switch_branch` — null-direction predictor pinning
  (`pending_switch`): the arc constraint forces the first step along
  the buckling mode so the corrector lands on the bifurcated branch
  (a state-only perturbation relaxes back to the fundamental branch;
  a basin-scale jump inverts elements — both measured failure modes).
- `rod` (bead tfz.14): geometrically exact Cosserat rods — nodal
  positions plus canonical `fs-ga::So3` frames updated
  MULTIPLICATIVELY through body/right group operations; strains Γ = Rᵀr′ − e₁ and
  κ = Log(Rᵢ⁻¹Rᵢ₊₁)/L₀ from RELATIVE quantities (rigid motions give
  exactly zero strain — energy-invariance battery). Statics: analytical native
  residual/tangent Newton with load stepping and frozen-Jacobian
  correction backtracking. The current LU is reused for trial solves;
  scaled displacement/rotation corrections choose accepted steps, while
  the original physical force/couple residual decides convergence.
  Constitutive and discretization validity retain their original
  admission and refinement requirements.
  Arbitrary energy callbacks retain an element-local finite-difference fallback;
  SE(3) dynamics under fs-time integrators remain successor scope.
  `Rod::is_stable_under_dead_load` independently checks the actual free residual
  and original analytical tangent of a uniform rod with a dead tip force and
  zero tip moment. Its symmetric second variation is diagonally equilibrated
  and must remain Cholesky-positive after a `64 * ndof * EPSILON` dimensionless
  margin. Negative, neutral or numerically unresolved tangents return false;
  unconverged states, couples, inadmissible parameters and more than 64 segments
  refuse. This immutable bounded query shares the raw rod's lack of cancellation
  context and does not change Newton convergence or its accepted state.
  `tests/rod_stability.rs` checks the clamped-free Euler threshold under refinement,
  force/length rescaling, rigid rotation and explicit refusal cases. These are
  local discrete stability screens, not dynamic or continuum stability proofs.
- `fiber` (bead tfz.14): fiber sections over fs-material `Uniaxial`
  cards (Menegotto–Pinto steel; Mander concrete with the
  compression-positive convention flipped at the section boundary);
  `respond`/`commit` split (pure response vs state commitment);
  `try_respond` uses the native checked Mander response, including finite
  parameters, envelope-consistent history and the inclusive ultimate compression
  strain. `try_commit` admits the complete section before publishing any proposed
  fiber histories; a later fiber's refusal cannot partially update earlier ones.
  The infallible `commit` convenience requires an admitted strain and panics
  before mutation on refusal. A virgin concrete fiber and a closed residual
  gap use the compression-side tangent; tensile gaps remain traction-free.
  `tests/fiber_material.rs` checks a concrete-only axial frame against the
  independent `r=2` envelope inverse, then verifies an over-crushing displacement
  refuses with the accepted frame state intact. Section tests exercise the
  inclusive bound, corrupted parameters/history and atomic commitment.
  `update_sections_batched` — the §15.2 hot loop through fs-la's
  batched Cholesky with a per-outlier direct fallback, throughput
  measured and ledgered; `rc_section` fixture (confined core +
  unconfined cover + steel layers).
- `fiber::material`: source-selected or explicitly authored Mander/MP laws
  share the existing state-point resolver and fiber engine. Every density and
  constitutive parameter is required; selected source bundles and receipts
  survive section construction and area scaling. Static source admission
  refuses unowned strain, stress, rate, frequency or history dependence.
  Selecting coefficients does not calibrate the cyclic model. `RcSection`
  partitions four-sided cover, core and two symmetric reinforcement bands;
  steel displaces concrete in both response and mass. Steel response is lumped
  at the declared band centroids, so this is smeared reinforcement, not resolved
  round-bar geometry. Geometry and source bindings stay immutable; consumer
  clones own independent constitutive histories. Area scaling changes mass and
  every fiber area together. The legacy `rc_section` remains an authored fixture
  with its original gross-concrete-plus-rebar convention. Focused G0/G1/G3 tests
  are in `tests/fiber_material.rs`; no measured-building validation is implied.
- `beamcol` (bead tfz.14): force-based (Spacone) cantilever element —
  EXACT equilibrium interpolation, five-point Gauss–Lobatto section
  integration (a point AT the hinge), per-section Newton with
  residual backtracking, virtual-work tip deflection; pushover and
  cyclic histories with committed, cloneable state (G4).
- `frame2d` (plan §15.2 step 3): planar frames of force-based
  distributed-plasticity beam-columns (`Frame2d`, `FrameElement`,
  `SectionModel::{Elastic, Fiber}`). Basic system `q = [N, Mᵢ, Mⱼ]` with
  EXACT equilibrium interpolation `D(ξ) = [N, (ξ − 1)Mᵢ + ξMⱼ]`, 3–7
  Gauss–Lobatto sections, Neuenhofer–Filippou state determination
  (compatibility correction + section Newton to compatibility AND section
  equilibrium; sub-incremented `v_c → v` in up to 64 substeps when the
  direct iteration fails), tangent `K_b = F⁻¹`. `Geometry::Linear` or
  `Geometry::Corotational` (exact chord rotation, consistent geometric
  stiffness, chord angle unwrapped against the committed angle). Analyses:
  `static_load` (load control), `displacement_control` (bordered pushover
  through softening), `natural_frequencies` (massless DOFs statically
  condensed), `newmark` (average acceleration, uniform base motion,
  Rayleigh damping on the initial stiffness via `Rayleigh::from_modes`) with
  a `DynamicHistory` energy ledger (input, kinetic, damping, internal work)
  and base shear. `newmark_with_load` retains an explicit constant nodal load
  throughout the base-motion record: `M a + C v + f_int = load - M r a_ground`.
  It includes that load in the initial acceleration and external-work ledger,
  so gravity-preloaded columns keep their compression and P–Δ response.
  `newmark` continues to mean zero sustained load; a preload is never inferred
  from internal forces. The energy audit uses the change in kinetic energy
  from the supplied initial state, alongside incremental internal work.
  Newmark solves directly for endpoint acceleration with tangent
  `M + gamma*h*C + beta*h^2*K_t`; position and velocity use the corresponding
  increments. It never recovers acceleration from two absolute displacements.
  Backtracking and the residual gate act on that same physical acceleration
  equation. Force work uses the resolved kinematic displacement increment,
  retaining momentum/work balance when the stored position cannot resolve it.
  `velocities()` and `accelerations()` expose the accepted base-relative
  generalized state, including angular components. Absolute displacement
  remains `f64`; compensated position or sub-ULP constitutive strain is not
  claimed. `tests/newmark_acceleration.rs` covers lumped and consistent mass,
  sub-ULP vibration, changed base forcing, damping and kinetic/internal work.
  Activating a fixed coordinate freezes its accepted displacement and clears
  its velocity/acceleration. This is a kinematic boundary reset; no resolved
  support-engagement impact or associated impulse/energy receipt is claimed.
  The initial Rayleigh reference tangent is frozen at the first admitted
  transient initialization (after any gravity/static preload) and retained
  across record boundaries, clones and retries after a partially accepted
  history. Static loading does not silently replace it. Adding a node or
  member invalidates the reference; `reset_rayleigh_reference` explicitly
  requests a new one without resetting mechanical or material state. Current
  point masses still enter the mass-proportional term, and supports restrict
  the same full reference operator. Explicit coefficients may change between
  records; that is a caller-requested change of the damping model.
  Every frame, including one without material viscosity, refuses nonfinite
  forcing/directions and negative/nonfinite Rayleigh coefficients before state
  mutation. For nonzero stiffness-proportional damping, the complete free-DOF
  C (including mass and material loss) must have nonnegative work within a
  floating-point sign-resolution check. Positive diagonal congruence scaling
  followed by the existing Cholesky owner with a scratch-only
  `64*eps*n_active` diagonal allowance admits semidefinite loss directions;
  negative diagonals and nonzero rows with zero diagonal are refused.
  The original C is never shifted, clipped or replaced in forces or work.
  This numerical admission is not an interval passivity certificate and does
  not suppress physical buckling in the stiffness/internal-force model.
  The focused frame tests cover nonlinear whole/split-history agreement,
  explicit reference changes, structural edits, negative-damping refusal,
  mass-stabilized damping and retained zero-loss directions.
  The `static_load` and `newmark` Newton steps backtrack
  (step 1, ½, … down to 2⁻⁹) on the free residual norm, because fiber laws
  are non-smooth at tensile-gap closure and load reversals. Concrete's virgin
  tangent includes its compression-side stiffness; subsequent trials still
  evaluate the original active branch, including zero stiffness in an open gap.
- `koiter` [F], feature `koiter-asymptotics` (off by default): FD
  energy expansion along the buckling mode at the critical state →
  a/b coefficients and `Bifurcation` classification, with the
  SAMPLED-CONTINUATION fallback oracle (the imperfect-geometry path)
  cross-checking imperfection tolerance.
- `contact` [F], feature `contact` (off by default): SDF-native
  log-barrier contact around `HyperProblem`, with `Barrier`,
  optional lagged smoothed-Coulomb `Friction`, `ContactProblem`,
  `ContactSolution`, and the frictionless `translation_gradient`
  adjoint probe. The solver requires a finite, intersection-free
  initial state; bad parameters, non-finite SDF values/gradients,
  invalid pins, and unsupported contact-adjoint inputs return
  `SolidError::InvalidInput`. `solve_from` holds the previous physical
  displacement fixed as the slip reference; `solve` uses the zero state.
  `solve_from_with_surface_displacement` additionally admits one finite
  world-coordinate material-displacement vector per node. Friction uses the
  relative increment `u - reference_u - surface_displacement`; the prescribed
  surface increment stays fixed and is not scaled by the load factor. The
  supplied SDF independently defines current obstacle geometry, so a conveyor
  can transport the solid without changing its planar chart. Reference-state
  feasibility is still checked against that chart; no surface trajectory or
  material correspondence is inferred.
  Numerical lag rounds refresh normal loads and tangent frames only,
  and success requires the refreshed friction residual to meet tolerance.
  The frozen-frame Huber energy, gradient, and Hessian agree in both
  sticking and sliding branches. Exact Newton directions are retained when
  finite, descending, and admitted by a true max-norm relative linear
  residual at most `sqrt(f64::EPSILON)`. Otherwise sliding search directions
  first add the Huber quadratic majorizer curvature
  `mu * normal_load / abs(slip)`. This positive search metric changes neither
  the physical Hessian nor the force law. A stiffness-scaled diagonal shift
  remains the fallback if the majorized system also fails direction admission;
  it preserves the physical force law, tangent, pin equations, and final
  residual gate. Both line-search acceptance branches
  require admitted material stresses as well as positive sampled gaps;
  an energy decrease alone cannot bypass the material's stretch domain.
  At satisfied pin equations, the residual-decrease shortcut also requires
  energy nonincrease within a scaled floating-point allowance, preventing
  alternating energy/residual merit cycles. A tiny admitted Newton correction
  may terminate despite unresolved potential differences only if its actual
  force residual already meets the unchanged tolerance; its coordinate size
  is bounded by `sqrt(EPSILON)` times the shortest reference edge in max norm.
  This numerical-resolution rule is not an energy-error certificate. Pending
  prescribed displacement equations retain their boundary-work path.
  `eps_v` is a displacement smoothing
  width; this regularized law permits small creep and does not claim exact
  set-valued static friction. Solutions expose separate normal and
  tangential reactions plus nonnegative endpoint relative friction work.
  Signed `friction_support_work = sum(f_t dot surface_displacement)` records
  supplied support work, with `f_t` the force on the solid. Endpoint products
  close `body_friction_work + friction_work = friction_support_work`.
  These products are not trajectory-integrated heat or actuator work;
  thermal deposition and rigorous continuous-trajectory CCD remain separate
  capabilities.
  `solve_timed_from` additionally admits a positive physical interval and a
  declared slip-speed regularization. Its displacement smoothing width is
  `dt * regularization_speed`, so subthreshold creep follows elapsed time
  rather than the number of equilibrium solves. Existing entry points retain
  their displacement-width semantics. `ContactTrajectory` binds one unchanged
  problem, holds accepted displacement and physical time, and accumulates body,
  relative and support work with backward-Euler endpoint quadrature. A failed
  solve or nonrepresentable time/work update changes none of that state.
  The cumulative ports retain `body + relative = support` up to roundoff;
  body work is independently evaluated from actual displacement increments.
  This is quasistatic velocity-regularized friction: it adds no inertia,
  recoverable tangential state, ideal static-friction constraint, thermal
  partition or continuous-trajectory accuracy claim. Changing physical loads
  and slip directions still require temporal refinement.
- `accept_scenario_bc`: fs-scenario integration — `Dirichlet` and
  `Traction` under `Physics::Elasticity` accepted; dimensioned-value
  resolution stays with the scenario consumer.
- `shell` (bead frankensim-b8bxd.9.1): a safe-Rust, deterministic,
  flat, single-patch triangular **estimate-only surrogate** for the Euler-disc
  base. It combines lumped inertia, CST isotropic plane-stress membrane
  elasticity, and normal-slope/rotation penalties. The in-plane operator
  reproduces constant-strain energy and equilibrium with the material's
  actual Poisson coupling, including auxetic ratios. The bending/rotation
  terms remain a surrogate, not a Kirchhoff--Love continuum discretisation.
  `ShellPlate::assemble` returns type-separated dense
  `ShellMatrix<Mass>`, `ShellMatrix<Stiffness>`, and optional Rayleigh
  `ShellMatrix<Damping>` operators in node-major `(x,y,z,rx,ry,rz)` order.
  Each triangle contributes exact affine plane-stress membrane energy plus
  normal-slope and rotation-gradient surrogate energy; translational and rotary inertia are lumped from the
  material card. `ShellSupport` applies three distinct, normal-aligned point
  pinned supports by reducing all translations at those nodes while retaining
  rotations; its normal identifies the plate side and it does not claim a
  frictional or rotational support law. Inputs carry model,
  source, and state IDs and fixed SI/right-handed-Cartesian interpretation.
  Assembly bounds node, element, dense-memory, deterministic
  validation-and-assembly work, and diagnostic sizes before allocation; it
  returns no partial operator. Validation charges finite-node scans,
  scale-relative coincident-node pairs, canonical face insertion, oriented
  edge incidence/adjacency, connectivity traversal, used-node marking, and
  triangle assembly against the declared work budget. Small assemblies report
  raw algebraic mass-array definiteness, stiffness-array nullity, and symmetry
  residue via a deterministic Jacobi eigensweep. Its raw mixed-unit spectral
  spread is diagnostic only, not a physical condition number or continuum
  certificate; larger valid assemblies explicitly report that those bounded
  algebraic diagnostics were not computed.
  Scale-relative coincident nodes, duplicate index faces, geometrically
  duplicated faces through coincident nodes, non-manifold edges (more than two
  incident faces), same-orientation shared edges, and components disconnected
  by an edge are refused before assembly; callers requiring several components
  must model separate plates.

- `reduce` module (music bead 3ez8g.3.3): the offline lip/reed reduce
  lab. `reduce_valve(request, cx) -> ValveCard` runs
  `TetLinearElasticProblem` -> `fs_modal::slice_window` -> per-mode
  re-normalization from the solver's `phi^T M phi = 1` convention to
  UNIT MEAN TRANSLATION of the orifice face (mass normalization would
  destroy the units). The mean includes every declared face node, with
  constrained components contributing zero; duplicate face nodes refuse.
  This mints physical `m_eff = phi^T M phi` [kg],
  `k_eff = lambda m_eff` [N/m], and `c = eta sqrt(k m)` from the card's
  AUTHORED loss factor (the lab does not invent Q values). The
  retained-compliance law: the modal compliance
  `sum (phi^T r)^2 / lambda` of the retained modes — `phi^T r` is the
  FORCE-pattern projection; the first build used the base-excitation
  `phi^T M r` and retained 1e-11 of the compliance, an executed
  falsification — is compared against the EXACT static compliance
  `r^T K^{-1} r` (sparse LDLT), and a reduction below the authored
  floor REFUSES (`RetainedComplianceTooLow`). Orifice geometry is
  MEASURED from the mesh: the face is the surface strip facing the
  opposing plane, gap = closest approach, width/channel-thickness =
  face extents. Opening and width axes must be perpendicular (normalized
  dot-product tolerance 1e-12); nonfinite opposing-plane offsets refuse
  before assembly. An oblique frame is not silently projected or rescaled.
  The card carries the provenance chain (mesh digest,
  material state-point receipt identity, source id) under a
  domain-separated schema-versioned content identity with
  canonical-bytes round-trip and tamper refusal. Loading requires at least
  one mode, bounds mode allocation by available input rows, and refuses
  trailing records after the identity. `beating_reed()`
  requires an explicit effective pressure area `A` [m²], conjugate to the
  mean-opening coordinate, and mints `fs_scenario::BeatingReed` with
  `P_c = k H / A`. No fixed pressure-face length is assumed. Distributed
  pressure-load projection is available through `reduce_pressure_loaded_valve`;
  authored areas supplied directly to the plain card remain research inputs.
  Empty cards and invalid converted parameters refuse.
  It also passes the card loss as explicit reed damping
  `zeta = eta/2`, preserving `c = eta sqrt(k m)` in the playing dynamics.
  This is the existing reduced viscous model, not a broadband loss-law fit.
  `reduce_pressure_loaded_valve` accepts one finite mesh-axis force-per-pascal
  vector [m²] per node. It projects the load onto each retained mode using
  `A = phi^T load / mean_face(phi)`. This is the same work-conjugate coefficient
  for generalized force `p A` and swept-volume flow `A v`; constrained DOFs
  do no virtual work. The loaded surface need not equal the opening-observation
  face. The returned `PressureLoadedValve` retains the card, signed per-mode
  areas and an identity covering card/load/projections, and converts its first
  mode without a separately authored area. Nonpositive first-mode area refuses
  the positive-closing reed conversion instead of taking an absolute value.
  The caller owns consistent surface-load assembly, pressure orientation and
  units. This does not discover wetted surfaces, solve a pressure field, or
  certify pressure-load modal convergence from the opening-compliance fraction.
  The existing card serialization remains structural; it does not serialize
  this pressure port. The cane integration assembles uniform-pressure loads
  from its planar underside triangles and exercises the resulting bore response.

### Three-dimensional finite-strain hyperelasticity (`hyper3`)

`HyperTetProblem { nodes_m, tetrahedra, material, prescribed_m,
nodal_forces_n, body_force_n_m3, budget, settings }.solve(cx)` is the 3-D
Total-Lagrangian counterpart of `hyper2d`: P1 tetrahedra (constant `F` per
element, so one card evaluation is exact quadrature), the `fs-material`
energy cards (Neo-Hookean, Mooney–Rivlin, decoupled Ogden, HGO/GOH;
AD Piola stress and consistent 9x9 tangent, with coalescence-safe spectral
derivatives for Ogden), DEAD loads (nodal forces, body force per reference
volume lumped to vertices), and linear ramping over `load_steps`. Each
Newton step solves `K_ff d = -r_f - K_fp dp` so pending prescribed motion
enters through the consistent linearization; trial states that invert an
element are halved; once the boundary data are applied an Armijo condition
on the potential `sum V W - f_ext . u` is enforced (residual decrease on a
non-descent direction). The gate is `||r_f|| <= relative_tolerance *
max(||f_int||, ||f_ext||, max|K_dd| max|u|)`; the stiffness term keeps it
meaningful at stress-free states (a rigid rotation, where every force
vanishes). Linear solves: ILU(0)-PCG, falling back to |diag|-Jacobi PMINRES
for an indefinite tangent; acceptance on the recomputed residual.
`HyperTetSolution` carries displacements, reactions in `prescribed_m`
order, stored energy, the smallest `det F`, and per-step residual histories,
backtracks, linear iterations and fallbacks. `internal_force` and
`tangent_matrix` are exposed for verification. Element quality and sizes are
admitted against `TetAssemblyBudget`; refusals are `HyperTetError`
(`InvalidInput`, `DegenerateElement`, `MaterialRefused`, `NewtonStalled`,
`LinearSolveFailed`, `Cancelled`).

Conservative static globalization integrates each material's direct finite
potential change in the existing trial force sweep. Dead-load work and the
combined potential difference use compensated arithmetic; the Armijo tolerance
depends on the actual changes and work, so an additive material energy datum
does not relax the decision. Uniform, affine, free-surface and stratified
pressure all use their original pressure owner's direct work. A bound gravity
source can supply its direct potential change in the trial load sweep. Legacy
gravity callbacks still expose endpoint potentials here: when their subtraction
cannot resolve the decision, an actual residual decrease is required. Sealed
cavities use original gas caloric increments and polarized enclosed-volume work
between their already admitted trial states. Custom material
adapters that retain the default endpoint subtraction
retain its precision limits. Force equilibrium and material/load admission are
unchanged. Numerical solve regressions retain the same nonlinear path with a
`1e20 Pa` material energy offset and resolve a pressure-loaded small cell beside
an immobile cell of approximately `1.7e20 m3` whose rounded global volume does
not change.

An already force-balanced static trial can be accepted when its final potential
decrement is unresolved. This check uses only the actual free internal/external
force norms, a lower bound on the native static normalization; supported
reactions and a stale tangent cannot relax it. All source/domain admissions and
the original tolerance remain required. The next outer iteration admits the
current tangent and certifies the full equilibrium gate. This avoids a
near-root line-search stall in weak cavity loading without changing the gas
entropy controls or treating numerical energy noise as physical work.

Backward-Euler globalization also consumes a direct increment of its original
held-predictor quadratic potential. With `q=(u-u_old)/dt-v_old`, compensated
polarization evaluates `Delta Psi=mean(q)^T M Delta q`, retaining small work
beneath a large predictor baseline. Both endpoint quadratic forms retain their
finite/nonnegative admission and previous floating-point reduction order.
Those endpoint values and the direct increment share one traversal of the
original consistent reference mass, including all signed quadratic-element
entries and cross terms. No lumped mass, copied tangent or second sparse sweep
is introduced. Physical kinetic energy, numerical loss, material dissipation
and the Backward-Euler scheme retain their existing meanings. Original-owner
tests cover weak predictor-relative work, common-offset invariance, independent
rational P2 mass work, reversal/subdivision, endpoint compatibility and refusal.

`HyperelasticLawField::Uniform` and `PerElement` use this same owner. A
per-element field follows tetrahedron connectivity order and must cover the
mesh exactly (`LawCountMismatch` otherwise). Each local law supplies every
trial energy, Piola stress and consistent tangent, including line-search and
final evaluations. Regional materials are not averaged. Fields may borrow
native `Hyperelastic` laws or `ResolvedHyperelastic` values retaining the exact
admitted constitutive model cards, including parameters, validity, source hashes
and provenance. Its explicit fixed external query must satisfy every source
bound and exact quantity descriptor. Static temperature, pressure, composition,
moisture and concentration conditions are supported; evolving mechanical axes
require their own owner and refuse. The context-free constructor only accepts
unconstrained source domains. The source wrapper delegates constitutive parameter
and native stretch checks to the existing hyperelastic owner; it introduces no
additional energy model.
The focused heterogeneous regression solves two source-defined Ogden layers in
series against an independent axial stress/energy formula and traction-balance
root, changes only the second layer's material data, and checks actual interface
motion and reaction. Uniform native/source/element fields agree bit for bit;
wrong cardinality and a single element's unsupported stretch refuse. Material
conditions remain fixed over the load path; thermal evolution, dynamic inertia
and source uncertainty propagation remain outside this static reduction.

Conformance (`tests/hyper3.rs`, all executed): the affine finite-strain
patch test is reproduced to 1.2e-14 m with the exact stored energy; the
compressible Neo-Hookean uniaxial stretch to 1.6 (rollers) matches the
closed-form lateral stretch and nominal stress to 1e-9 in 4-5 Newton
iterations per step; the tangent matches a central difference of the
internal force to 4e-11 (Neo-Hookean and Mooney–Rivlin); at the reference
state it equals the independently written `linear3` stiffness to 4e-16; a
rigid rotation of the boundary stores 1e-17 J with reactions below 1e-15 N;
inverted elements, repeated prescribed DOFs and cancellation refuse.

No-claims (`HyperTetProblem::solve`): static, isothermal, rate-independent; constant-strain
tetrahedra (volumetric locking of nearly incompressible cards is not
mitigated; no mixed or F-bar formulation); dead loads only (no follower
pressure); no contact; no limit-point continuation (drive by displacement
or refuse); no mesh-convergence claim from one run.

#### Complete-boundary current-area pressure (`hyper3::pressure`)

`ClosedTetPressureSurface::new` derives the entire exterior of the original
tetrahedral mesh through `TetLoadGeometry3`'s incidence, connectivity and quality
admission. It orients each exterior triangle outward and requires every
boundary edge to have exactly two oppositely directed incidences. Interior
faces receive no load. Duplicate/same-side tetrahedra, nonmanifold faces or
edges, open/inconsistent boundaries, unused vertices and work-budget overflow
refuse. The prepared owner retains exact reference-coordinate bits and ordered
connectivity for binding to a consumer's original geometry. Either reference
tet winding is admitted; every current tet must preserve its original sign.

For each outward triangle `(a,b,c)`, `force_n` places
`-p/6 * ((x_b-x_a) cross (x_c-x_a))` on each of its three nodes. Positive signed
pressure acts inward on the current area. The analytic `force_tangent` is the
external force Jacobian: perturbing a face vertex by `h` varies the cross
product by `h cross (x_b-x_c)`, `h cross (x_c-x_a)` or
`h cross (x_a-x_b)`. A consumer subtracts this Jacobian from its internal
residual tangent. On the complete closed boundary these assembled loads equal
`-p grad(V)` and the tangent equals `-p Hessian(V)`.

`volume_m3` sums current oriented tetrahedron volumes with finite arithmetic
and inversion/collapse refusal. `work_j` returns `-p (V_after - V_before)` for
pressure held between the two endpoints. This exact discrete volume difference
is not replaced by midpoint traction work; P1 volume is cubic in position.
Finite signed pressure, geometry and all returned arithmetic are admitted,
work is bounded by the reference node/tet/face limits, and cancellation is
checked during geometry and load operations. The owner changes no material,
reference mass or accepted trajectory and introduces no time solver.

`tests/hyper3_pressure.rs` specifies independent scalar-triple-product force
and volume oracles, finite-difference tangent and symmetry, rotated and
translated force/resultant/moment checks, a two-tet union with its interior
face removed, signed work, invalid/budget and cancellation refusal. Native
parsing passed. Full compilation and execution remain pending because the
shared executor killed the real `asupersync` prerequisite even during a
metadata-only check. No selected open patch, spatially varying pressure,
fluid equation of state, global injectivity certificate or continuous-time
inversion guarantee is claimed. Static `HyperTetProblem::solve` remains the
dead-load owner; finite-solid dynamics consumes this separate pressure kernel.

### Incremental small-strain inelasticity (`plastic3`)

`SmallStrainTetProblem { nodes_m, tetrahedra, laws, prescribed_m,
nodal_forces_n, load_path, budget, settings }.solve(cx)` follows a load path
on P1 tetrahedra for ANY `fs-material::SmallStrainLaw` (J2 radial return
with the Simo–Hughes algorithmic tangent, isotropic/orthotropic
elasticity). Tensor-component Voigt strain `eps = B u`, internal force
`sum V B^T W sigma`, tangent `sum V B^T W C_alg B` with `W = diag(1,1,1,2,2,2)`.
Each load factor (which may decrease or reverse) is one Newton solve from
the last converged state, prescribed increments through the consistent
linearization, residual-increase halving, and the same linear-solve policy
as `hyper3`. The static force gate uses element stress-force bounds on free
rows and actual free external loads, retaining the largest accepted scale
through each unloading increment. Fully fixed regions and rigid travel do
not weaken another region's equilibrium tolerance. Constitutive nonfinite
stress, force or tangent refuses before publishing state. Element states are updated by the law and
committed only after convergence. The solution carries displacements,
element strain/stress/states, and per-increment reactions, residual
histories, and counts of evolving elements. `internal_force` and
`tangent_matrix` take explicit committed states for verification.
The convergence check precedes tangent assembly. The element tangent kernel
precomputes its twelve constitutive strain-column contractions once on the
stack and reuses them for all test columns; no changed convergence order,
speedup factor or mesh-accuracy claim is inferred from this operation reduction.

`SmallStrainLawField::Uniform(&law)` and
`SmallStrainLawField::PerElement(&laws)` use this same numerical owner.
Per-element entries follow tetrahedron connectivity order and may come from
different material cards or spatially varying resolved conditions. Each
tetrahedron initializes its own law's state and uses that law for every trial
stress, consistent tangent, converged history update, and final stress.
A wrong field cardinality returns `LawCountMismatch` before assembly.
No effective-property averaging or separate heterogeneous solver is introduced.
Constitutive coefficients remain fixed throughout the load path.
Reported endpoint stress is evaluated from the same previous committed history
as the converged force and tangent, before publishing the next internal state.
An explicit-duration overstress law therefore advances exactly once per load
increment; reading the result does not perform an extra relaxation interval.
The same publication rule serves the backward-Euler inertial adapter.
`ResolvedJ2Plasticity` implements the native `SmallStrainLaw` interface directly,
so fields can borrow complete resolved material inputs and retain all source
receipts without extracting bare coefficients. This field alone does not add
thermal feedback, state remapping, or finite-strain mechanics.

The homogeneous numerical measurements below predate this field extension.

The source-selected `MixedPlasticStrainModuli` J2 extension uses this same
P1 owner and retained per-element plastic strain. Its Prager backstress changes
reverse yielding while the supplied isotropic and kinematic moduli can retain
the same monotonic tangent. A focused authored regression in
`tests/plastic3_material.rs` drives tension, unloading and reversal on the
original six-tetrahedron coupon and compares reactions and accumulated/signed
plastic strains with independent uniaxial equations. This introduces no
alternate solver, fatigue-life estimate or finite-strain plasticity claim;
the test's presence is not a claim of native execution.

Conformance (`tests/plastic3.rs`, executed): homogeneous uniaxial tension
with linear hardening (E = 1000, nu = 0.3, sigma_y = 1, H = 100) follows the
closed form through yield to eps = 0.006 (force 1.4545454547 vs
1.4545454545) and elastic unloading to the plastic strain leaves zero force
(2.6e-11), alpha = eps_p and the isochoric lateral contraction eps_p / 2 to
1e-8 relative, in 2-5 Newton iterations per increment; the algorithmic
tangent matches central differences of the internal force to 4.6e-10 at a
plastic state; the elastic-law tangent equals `linear3` to 4e-13 relative; a
cantilever driven by tip displacement yields progressively (1, 20, 33 newly
yielding elements) with a falling secant stiffness.

No-claims (`plastic3`): small strain and rotation; no implicit clock
(rate-dependent laws require explicit physical intervals);
constant-strain tetrahedra lock under confined plastic incompressibility
(confined limit loads are over-predicted; no B-bar/mixed formulation); dead
loads; no limit-load continuation; no mesh-convergence claim.

### Sub-yield creep on the native solid (`plastic3::creep`)

`SmallStrainTetProblem::solve_creep_from` consumes the accepted nodal
displacement and per-element `CreepState` through the original sparse Newton
owner. Uniform `NortonCreep` and source-bound `ResolvedNortonCreep` fields
retain the same constitutive flow, tangent, source receipts and reference
geometry. Exactly one absolute load endpoint represents one positive common
physical duration; the adapter refuses multi-entry paths or unequal element
durations. Newton and backtracking always borrow the old history. The published
stress is the converged trial stress evaluated before committing that history,
so observing a rate-dependent endpoint cannot advance another interval.

The adapter admits finite old/new material histories and bounds the full
displacement gradient by each law's small-strain limit, including its skew
part. Each checked constitutive increment is compared with the actual native
stress and committed state. Per-element accounts in J/m³ and reference
volumes are returned for conservative physical-heat transfer. The whole-solid
audit independently compares applied-force and support-reaction endpoint work
with integrated stress work, then with elastic storage change, physical creep
dissipation and backward-Euler numerical loss. Numerical loss is not heat.
All local and global gates and cancellation precede publication; caller-owned
displacements and histories are never mutated.

`advance_creep_backward_euler` couples that same checked material interval to
the original `plastic3::dynamics` consistent reference mass and sparse Newton
owner. Its explicit duration must equal every element's material duration;
inertia requires the absolute load endpoint `[1.0]`. The returned
`CreepTetDynamicStep` contains the checked material account, accepted velocity,
and separate inertial account. Applied and moving-support endpoint work balance
material work plus inertial work; the complete energy balance includes kinetic
change and backward-Euler kinetic loss. Neither kinetic loss nor the material
integration loss is caloric heat. Reference density and geometry remain the
original closed inventory. The scheme is first order and requires timestep
refinement; this extension does not establish a transient accuracy envelope.
`tests/plastic3_creep_dynamics.rs` checks repeated accepted steps against an
independent one-coordinate mass/creep recurrence, inertia in moving-support
work, stress-free rigid translation, mismatched clocks and cancellation.

Elastic coefficients and geometry must remain fixed across an interval.
Prescribed temperature may change Norton kinetics between intervals while
retaining creep history; elastic parameter work, thermal expansion, reciprocal
thermal evolution require their separate owners. This is a
small-strain, secondary-creep model, with no primary/tertiary creep, rupture,
confined-flow accuracy or material-population qualification claim.

`tests/plastic3_creep.rs` checks a loaded six-tetrahedron bar against analytical
constant-stress Norton creep, permanent set after unloading, held-length
relaxation against the backward-Euler recurrence and three timestep rungs,
per-element temperature-dependent histories and heat, and transactional
refusals including cancellation. These are manufactured implementation checks,
not measured validation.
The source-bound regression warms an accepted loaded trajectory, preserves
every pinned claim, and substitutes a material with twice the reference rate;
the actual bar extension follows the independent accumulated creep law.

### Simultaneous creep and yielding on native tetrahedra

`SmallStrainTetProblem::solve_creep_from` consumes the material owner's
`CreepLaw` and preserves its complete associated history. The same native
Newton solve, small-strain/rotation checks, endpoint stress/state comparison,
support/applied work and physical/numerical energy audit serve Norton creep
and simultaneous isotropic Norton/J2. `CreepTetStep` defaults to the original
creep state; combined steps retain separate creep and plastic memories.
Per-element accounts expose both physical heat channels and their total;
energetic hardening storage and backward-Euler loss remain excluded from heat.
The single physical interval, bounded work and immutable proposal semantics
are unchanged. `tests/plastic3_creep_plastic.rs` exercises an actual loaded
tetrahedron yielding and continuing to creep, unloading, held-shape relaxation,
analytical deformation/energy and preservation of both histories on refusal.
Finite deformation, inertia and noncoaxial kinematic hardening remain outside
this checked creep path.

### Permanent creep in native fiber members and frames

`FiberLaw::Norton` and `FiberState::Norton` consume the original material
owner's free-sided axial reduction, with signed permanent strain, accumulated
creep and optional original source receipts. The existing section integration,
force-based member compatibility and global frame Newton solver remain shared.
No effective modulus, alternate beam solver or inferred creep coefficient is
introduced. Explicit authored reference density in `FiberMaterial::norton`
does not acquire source authority from the law's separate coefficient receipts.

`Frame2d::solve_creep_step` admits one absolute load endpoint and one common
positive physical duration across the original Norton fiber points. Section
Newton, compatibility retries and global backtracking see accepted history;
the returned frame commits once. Stored forces and `creep_fiber` observations
do not evaluate another physical interval. `with_creep_laws` refreshes time
and prescribed temperature while retaining native geometry, history, elastic
and kinetic parameters, original card/claim identities and nonthermal context.

The existing untimed static, pushover and modal routes, and the inertial
Newmark route, refuse Norton fibers. They do not own the explicit creep clock
or its checked work account; an algorithmic creep tangent cannot silently
become an instantaneous modal stiffness or advance on a mismatched time step.

The checked step integrates each constitutive energy account using its original
fiber area, Lobatto weight and member length. Applied/support work is checked
independently against this material work and against storage change plus
physical creep loss and numerical loss. Force and moment gates use separate SI
tolerances. Physical heat and point volumes remain available to thermal
consumers; numerical loss is not heat. Immutable proposals, explicit fiber-work
budgets and cancellation checks preserve accepted history on refusal.

This path admits stationary inertia, linear geometry and Norton or joint
Norton/J2 sections.
It bounds nodal rotations, chord displacement and absolute quadrature curvature;
unsampled spatial extrema require refinement. Corotational geometric work,
inertia, thermal expansion, reciprocal temperature evolution, unsupported
inelastic fiber laws and measured-life qualification are outside this checked account.
`tests/frame_creep.rs` exercises constant-force extension/permanent set,
held-length tension relaxation, actual cantilever curvature/tip shape under
sustained moment, source-temperature continuation and atomic refusal/cancellation.
It also verifies that the untimed and inertial entry points cannot age fibers.

### Rod couples held in laboratory coordinates

`rod::WorldTipLoad` declares a world-held tip force in newtons and a free couple
about that moving tip in newton metres. `Rod`, `MaterialRod` and
`SpatialMaterialRod::solve_static_world` use the original energy, Lie-group
Newton and load stepping with frozen-Jacobian correction backtracking.
Each trial evaluates
`M_body = R_tip^T M_world`; the exact external contribution to the body/right
residual tangent is `-hat(M_body)`. It is included alongside both native
analytical internal tangents and the arbitrary-energy finite-difference fallback.
The constant couple is not replaced by the gradient of `-M dot Log(R)`.

At each Newton iterate, rod statics factors the current free-coordinate
Jacobian once. Backtracking reuses that LU and compares each trial's correction
norm with the current iterate's correction norm. Translations use the existing
translation scale; body/right rotation increments are dimensionless. A trial is
accepted when it satisfies the original physical residual tolerance or strictly
reduces the correction norm under that same frozen Jacobian. The raw force/couple
residual may rise during accepted iterations; its original tolerance remains
the sole equilibrium convergence gate. The Jacobian is refreshed at the next
Newton iterate. The returned residual vector contains final physical residuals
for completed load steps; a `NewtonStalled` history need not be monotone.

The unchanged finite-bending world-couple regression and the complete rod
material/source integration targets pass with the original loads, tolerances,
Newton limits and backtrack limits. This step-acceptance change preserves the
original constitutive model, reference inventory and material publication rules.

`TipLoad` continues to express a world-held force plus a body-held follower
couple. These are distinct physical load laws when a tip turns. World loading
composes with the existing source-resolved shear/torsional stiffness, spatial
materials, thermal sections, intrinsic geometry, gravity, tensile links and
fixed supports. `support_reactions_world` uses the identical load transform at
the accepted pose. The material owners retain atomic publication, source
admission and existing cancellation boundaries; no separate rod solver,
material averaging or load-induced change of reference mass is introduced.

G1/G3 tests check the load tangent against independent rotation differences,
laboratory-frame covariance, actual whole-rod force/torque balance after finite
bending, distinction from follower loading, fixed-support frame conversion,
source-bound torsional compliance and temperature-dependent shear modulus.
These remain bounded static equilibrium results. A residual is not a stability
or continuum-error certificate; no dynamic trajectory, accumulated torque work,
general conservative torque potential or material calibration is claimed.

### Intrinsically curved material rods

`SpatialMaterialRod::try_from_reference` binds `RodReferenceGeometry` to the
existing per-segment material sections. Both present and unloaded strain measures
use the same midpoint directors/SO(3) logarithm; elastic energy and resultants use
their difference. A cambered or twisted unloaded part therefore carries no
invented bending/torsional preload. The original Newton, support reactions,
reference-mass gravity and local thermal-section owners remain authoritative.
`RodMaterialSource::with_reference_geometry` preserves this datum, source pins,
and mass through temperature refresh; gravity uses the actual reference nodes.

Inputs declare an unloaded shape with equal material-coordinate intervals and
resolved neighboring rotations below pi. A loaded scan cannot supply unloaded
geometry without an independently justified reconstruction. Refine nodes,
frames and section assignments together. This remains bounded static slender-rod
mechanics with circular concentric sections and finite-difference tangents;
it does not supply dynamic stick/ribbon coupling, scan accuracy, contact,
self-intersection certification, or experimentally validated instrument audio.
G0/G1/G3 tests in `rod_material_field/reference.rs` cover intrinsic zero stress,
a three-rung circular-arc end-moment limit, material changes, gravity datum and
objectivity; `rod_material_source.rs` covers source/temperature continuity.

### Source-bound tensile members on material rods

`SpatialMaterialRod::with_tension_links` attaches up to 64 immutable
`MaterialTensionLink` bundles to existing nodes, including body-frame offsets.
Each keeps its original `ConcentricElasticSection`, natural length, strand count
and explicit tensile-strain limit. At the actual deformed endpoint distance L,
the complete bundle stores EA/(2 L0) max(L-L0,0)^2; slack links exert no compressive
force. The original rod energy/Newton/support-reaction owner sees this energy
once. Rotating offsets therefore transmit reciprocal forces and moments through
the same equilibrium, including on intrinsically cambered rods. Public link
observations report span, tension, storage and endpoint force/moment pullbacks.

The source-bound reference mass uses natural length, not the stretched span.
Explicit gravity includes half that mass at each endpoint and the corresponding
rotating offset potential and moment. Replacing a link rebuilds these weights
from retained inventories. `reference_mass()` remains the rod-only receipt;
`total_reference_mass_kg()` includes links and each link retains its own source.
`gravity_nodal_forces_n()` and `gravity_nodal_moments_body_nm()` expose the
matching loads. The separate link is held at its explicitly resolved material
state; rod heating does not silently refresh its properties or natural length.

This is bounded quasistatic straight-link mechanics with tensile linear strain,
not a sagging cable, dynamic bow ribbon, thermal link, rosin/contact law, damage,
or validated instrument model. Nonstatic material claims, invalid attachments,
unrepresentable inventory and trial strains outside the declared domain refuse.
G1/G3 tests in `rod_material_field/tension.rs` cover analytic parallel compliance,
energy-derived offset reactions, camber/material preload sensitivity, conserved
mass, gravity work, slack and refusal preservation. The existing
`fs-couple` `heated_material_rod` example accepts an explicit `--tension-link`
natural-length/offset/strand-count record and reports its physical observations;
its separately cold 300 K member is not a horsehair material preset.

### Shared material fluid loading on submerged rods

`SpatialMaterialRod::with_submerged_fluid` consumes the existing
`fs-material::FluidTransport` interface and the specimen's explicit gravity.
The uniform hydrostatic pressure resultant is `-rho_fluid V g`; its reference
volume uses the original layered section geometry rather than an effective
material density. `ReferenceDisplacementVolume` explicitly distinguishes flooded
central voids (annular material volume) from sealed cavities (outer envelope).
Solid mass, density receipts and native motion inertia remain unchanged.

The original shape-function integration retains total displaced volume and its
first moment, including clamped nodes. Attached strands use their natural-length
volume and strand count. Their rotating offsets contribute a separate buoyant
potential, force, moment and analytical tangent through the same rod Newton,
support-reaction and motion owners. Replacing fluid, gravity or links rebuilds
the load without accumulation; zero gravity retains an explicitly bound fluid.
The original material sources retain measurement authority; the shared numerical
transport snapshot does not invent a source receipt.

All represented bodies must be fully immersed with nonoverlapping displaced
volumes. This is an explicit fixed-reference-volume approximation: current
solid density does not redefine mass or displacement volume during thermal
resolution. No partial immersion, free surface, radial expansion, fluid
inventory in sealed cavities, capillarity, drag, added mass or compressible
hydrostatic profile is inferred. Material-fluid replacement is static specimen
construction, not a conservative time-dependent environment update.
Focused `rod_gravity` regressions cover density-driven equilibrium and neutral
buoyancy, layered/hollow volume and first moments, gas-pressure changes, thermal
reference inventory, eccentric strand work/moments, motion inertia, and atomic
load replacement/cancellation. These synthetic checks make no measured-fluid
or experimentally validated material-response claim.

## Invariants

1. Patch tests exact: linear displacement fields are reproduced to
   solver tolerance on P1, Q1, B-bar, and the CutFEM frontend
   (all-embedded fixture) — sol-001.
2. MMS optimal orders per implemented family × frontend: L2 ≈ 2,
   H1 ≈ 1 for P1/Q1/CutFEM-Q1 (sol-002).
3. Objectivity through the fs-material interface: W(RF) = W(F) and
   P(RF) = R·P(F) for both cards (sol-003).
4. Consistency merge gate: the assembled tangent equals the FD
   directional derivative of the residual, and the residual is the
   exact gradient of the potential energy (rel ≤ 1e-6) — sol-004.
5. Locking is MEASURED, not assumed: on thin near-incompressible
   bending, the standard element's error degrades by ≥10× from
   ν = 0.3 to ν = 0.4999 while B-bar stays within 2× (sol-005).
6. Cook's membrane envelope: converged tip at the literature reference
   point (48, 52) inside [23.5, 24.5] (plane stress, ν = 1/3), ≤2%
   self-convergence deviation; near-incompressible B-bar
   self-converges while the standard element's shortfall is logged
   (sol-006).
7. Newton robustness: large-strain bending and a buckling-adjacent
   compressed strip converge under load stepping + line search with
   fast terminal contraction; histories and backtracks are evidence
   (sol-007).
8. Determinism: BTree/insertion-ordered assembly, deterministic
   Krylov solvers; bit-identical across runs on a platform —
   including continuation trajectories (checkpoint/resume proven
   bitwise-equal in stab-002).
9. Buckling: the Euler strut's EXTRAPOLATED critical load lands
   within 3% of the analytic value with the Richardson indicator
   covering the raw fine-mesh gap (Q1 parasitic-shear inflation,
   reported, stab-001); pencil symmetry to 1e-10 and mode
   K-orthogonality to 1e-6.
10. Continuation: validated against the CLOSED-FORM von Mises truss
    (limit points within 5–10% of analytic, manifold deviation
    ≤ 1e-7, load control jumps where arclength traces, stab-002);
    the continuum tent traces the full snap-through Z-curve on one
    path — snap and snap-back limit points, load recovery beyond the
    load-control failure point (stab-003).
11. Branch handling: bifurcation detected within 15% of the pencil
    prediction on the compressed column; null-direction switching
    lands on the bent branch (transverse growth ≥ 20× measured,
    stab-004).
12. Flat plate-surrogate energy consistency: each assembled stiffness contribution is a
    positive weighted outer product, so the free operator is symmetric
    positive-semidefinite while positive lumped translational/rotary inertia
    makes the raw mass array positive-definite. A single free triangular panel has the six
    infinitesimal rigid motions in its stiffness null space; the focused
    `shell_operator` battery verifies this and independently reconstructs the
    quadratic energy from `Kx`. Its G1 membrane patch checks compare continuum
    plane-stress energy on three spatial rungs under normal, biaxial, shear,
    and combined strains, with positive and negative Poisson ratios and a
    jointly rotated geometry/displacement frame. Interior forces balance,
    and prescribed uniaxial extension with Poisson contraction reproduces the
    independent traction resultant. These checks qualify only the in-plane
    CST response, not the bending surrogate.

## Error model

`SolidError` teaching errors: `SolveFailed` (linear gate missed),
`NewtonStalled` (carries the residual history; repair = more load
steps), `MaterialRefused` (det F ≤ 0 escaped globalization; repair =
smaller steps), `UnknownPatch`, `UnsupportedBc` (fs-scenario mapping
outside Dirichlet/Traction × Elasticity), and `InvalidInput`
(non-finite/inconsistent caller data such as invalid contact
parameters, infeasible starts, or malformed adjoint seeds).
`InternalInvariant` distinguishes a corrupted lower-level solver graph from
repairable caller input.

The CutFEM facade maps only `CutFemError::SolveNotConverged` to
`SolidError::SolveFailed`, preserving its iteration count and relative
residual. A canonical `ConstraintCycle` maps to `InternalInvariant`; every
material, parameter, geometry, callback, topology-admission, and
unsupported-regime refusal, including an invalid typed edge band or a cut
supported segment, maps to `SolidError::InvalidInput` with the canonical
diagnostic. No canonical refusal is converted into a solution.

## Determinism class

Bit-deterministic across runs on a fixed platform (ordered assembly,
fs-solver CG/MINRES, no threading, no ambient state). CutFEM determinism
and topology ordering are inherited directly from fs-cutfem rather than
duplicated in this crate. Cross-ISA golden hashes are not yet recorded
(follow-up).

## Cancellation behavior

Bounded synchronous loops (assembly, capped Krylov iterations, capped
Newton with capped backtracking). Chunked Cx polling belongs to the
fs-exec driver (the L3 discipline; fs-feec/fs-cutfem precedent).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

`koiter-asymptotics` (off by default): the plan-flagged [F] Koiter
post-buckling module and its battery, per the Ambition-Tag gating
rule.

`contact` (off by default): the plan-flagged [F] SDF-native barrier
contact prototype and its battery, per the Ambition-Tag gating rule.

### Structural invariants (tfz.14)

12. Rod objectivity: superposed rigid motion leaves strain energy
    invariant to 1e-12 (str-001).
13. Elastica: large-deflection tip matches an in-test SHOOTING oracle
    (RK4 + bisection on the elastica BVP) within 0.015·L at
    PL²/EI = 2 (str-002).
14. Helical family: pure bending gives the EI/M circle (2%), pure
    torsion stays straight (1e-6) at rate M/GJ (2%), combined moments
    give the constant helix strain state (3%) — str-003.
15. RC hysteresis: per-cycle dissipation positive and growing with
    amplitude, peak moment inside the hand capacity band, section
    states bitwise-deterministic across runs (str-004).
16. Batched section updates are scalar-consistent to 1e-8 with
    ledgered throughput (str-005).
17. Force-based pushover softens past yield (secant < 0.5× initial),
    dissipates under reversal, and resumes from a checkpoint
    BITWISE (G4) — str-006.
18. Contact barrier calculus and solve hygiene (feature `contact`):
    barrier derivatives match finite differences, invalid/penetrating
    contact inputs return structured `InvalidInput` instead of
    panicking, accepted iterates keep positive gap, reactions balance
    the applied normal load, tight-clearance and curved-obstacle
    fixtures retain positive gap, lagged friction separates stick/slip
    cases, and the frictionless translation adjoint matches FD on a
    plane — cnt-001..007.

11. Reduce-lab gates (rl-001..rl-004): the thick cantilever-slab
    oracle (f above the slender Euler-Bernoulli line by the plate-strip
    factor, measured +11.4% inside an authored 0-35% band;
    strip-normalized effective-mass ratio 0.355 in (0.30, 0.42);
    `k = lambda m` at machine precision; the
    `k_ratio = f_ratio^2 * mass_ratio/0.25` self-consistency identity;
    retained compliance 0.98 > 0.95), the minted card driving the
    existing fs-phs valve island to its quasi-static equilibrium at
    1e-3 with supply defect < 1e-8, refusals (empty/out-of-range face,
    the retained-compliance floor via a window that excludes the
    compliant fundamental, a face behind the opposing plane), and the
    canonical-bytes round-trip + tamper refusal.

## Conformance tests

`tests/battery.rs`: sol-001 patch tests (body-fitted max nodal error
≤ 1e-9; CutFEM at the CG-tolerance floor ≤ 1e-7); sol-002 MMS orders
(P1 2.00/1.00, Q1 2.00/1.00, CutFEM 2.19/1.04); sol-003 objectivity
(energy rel ≤ 1e-10, stress rotation ≤ 1e-8, 8 random states × 2
cards); sol-004 tangent/energy consistency (≤ 1e-6); sol-005 locking
battery (standard degrades ≥10×, B-bar ≤2×); sol-006 Cook's envelope;
sol-007 Newton robustness (5-step large bending, 8-step compression,
terminal contraction ≥ 100×); the selection-policy probe. Unit tests:
mapped-mesh areas, selection thresholds. The CutFEM ownership probes
compare facade and direct canonical coefficients, nodal values,
residuals, iteration counts, active cells, and error norms bitwise after
the documented stabilization translation. A typed right-edge band likewise
compares coefficients, nodal values, and assembled-load compliance bitwise and
proves that an occupied legacy callback slot refuses. The probes also require
invalid material data, the uncertified near-incompressible regime, and an
empty domain to fail closed as `InvalidInput`.

`tests/stability.rs` (bead tfz.15): stab-001 Euler pencil (analytic
G2 + Richardson indicator); stab-002 von Mises truss closed-form
continuation oracle (+ bitwise resume); stab-003 continuum tent
snap-through with the self-calibrating load-control-failure probe;
stab-004 branch detection + null-direction switching; stab-005
eigenvalue-derivative gradient gate (rel ≤ 1e-3 vs frozen-K_G FD) and
the clustered-eigenvalue trap (min() kink ≈ 2 measured, KS aggregate
smooth with derivative matching FD to 1e-6, conservative).
`tests/koiter.rs` (feature-gated): stab-006 symmetric-stable
classification + imperfect-geometry sampled oracle.

`tests/contact.rs` (feature-gated): cnt-001 barrier calculus and
force divergence; cnt-002/003 intersection-free hard press and
reaction balance; cnt-004 tight-clearance and curved-obstacle gap
audits; cnt-005 lagged friction cone separation; cnt-006
frictionless contact translation adjoint versus FD; cnt-007
structured invalid-input returns; cnt-008 fixed-reference friction balance;
cnt-009 unpinned conveyor transport, common-translation invariance, motion
reversal and signed endpoint work closure, zero-motion equivalence, and invalid
surface-displacement refusal (G1/G3, force/translation tolerance 1e-9 and
endpoint work balance 1e-14); cnt-010 fixed-duration 1/2/4/8-step friction
against an independent force-balance creep oracle (1e-10), timed conveyor
reversal and cumulative work closure (1e-14), held-state zero work and
prescribed boundary work;
cnt-011 timed/default equivalence, invalid time/speed/product admission,
representable-clock progression and atomic nonlinear-failure behavior.

`tests/structural.rs` (bead tfz.14): str-001 objectivity; str-002
elastica vs shooting oracle; str-003 circle/twist/helix; str-004 RC
moment-curvature hysteresis + determinism; str-005 batched
consistency + throughput ledger; str-006 force-based pushover,
reversal dissipation, G4 bitwise resume.

`tests/frame2d.rs`: a single elastic element reproduces the Euler–Bernoulli
cantilever tip deflection, rotation and axial shortening to 1e-12; a
fixed-base portal matches the closed-form lateral stiffness
`(24EI_c/h³)(6ρ+1)/(6ρ+4)` to 1e-6; the corotational cantilever under tip
moment `πEI/L` rolls into the exact half circle (tip within 0.5% L, rotation
π to 1e-9); tip-mass frequency `√(3EI/mL³)` to 1e-6 and its softening
`ω²/ω₀² ≈ 1 − P/P_cr` at half the Euler load; RC pushover with exact
`Mᵢ = H·L` equilibrium bounded by an independent section moment envelope and
base-concentrated curvature; cyclic pushover dissipates energy over closed
cycles; sustained-gravity equilibrium for linear and corotational members,
loaded/unloaded linear superposition with the exact sustained-load work,
continued-state kinetic-energy accounting, and compressed-column vibration
against the gravity-loaded tangent oscillator; elastic SDOF Newmark against
the exact Duhamel response (< 0.1% of
peak) with energy balance < 1e-6; a 3-story 2-bay RC fiber frame with gravity
preload, modal Rayleigh damping and an 8 s synthetic record (energy balance <
2%); bit-identical replay; malformed-input refusals.

## No-claim boundaries

- `frame2d` members are planar Euler–Bernoulli (no shear deformation, torsion,
  warping or out-of-plane buckling); reference mass is lumped by default or
  consistently integrated when explicitly selected; no
  shear failure, bond-slip, joint-panel flexibility, P-Δ of leaning columns,
  or soil–structure interaction; the corotational formulation is large
  rotation / small strain. The dense LU global solve is sized for building
  frames of a few hundred DOFs, not for sparse large models. Fiber-section
  response inherits fs-material's Mander / Menegotto–Pinto cards and their
  no-claims.

- `linear3` claims infinitesimal-strain elasticity and work-conjugate thermal
  eigenstrain loading on a fixed, conforming, non-degenerate tetrahedral mesh.
  The thermal operator does not infer expansion coefficients, solve a transient
  response, remesh, establish mesh/stress convergence, or claim finite-strain
  validity. Its coordinate update is fixed-topology and inversion-gated, not a
  finite-strain constitutive update. Its static solve establishes only its reported linear
  residual on the supplied discretization. It does not claim plasticity,
  physical damping, phase change, fracture, moving interfaces, or topology
  evolution. Callers must escalate or refuse when the resolved material state
  leaves this solid rung.
- Its free-body stiffness correctly contains rigid modes. Removing, fixing,
  or slicing around those modes is a boundary-condition/eigenproblem decision;
  the assembler never adds hidden grounding stiffness.

- The facade inherits fs-cutfem's no-claim for clipping a supported traction
  segment cut by the SDF and for quadrature error bounds on arbitrary traction
  callbacks. It only delegates checked support; it does not certify callback
  smoothness or grow a second boundary assembler.
- Extreme `(E, beta, gamma)` combinations whose active physical stabilization
  products are finite but whose fs-cutfem dimensionless coefficients are not
  representable are not claimed. Nitsche-active and ghost-active terms fail
  closed pending a canonical physical-coefficient API; an unused Nitsche value
  in traction-free mode is ignored rather than validated.
- TDNNS-proper / weakly-imposed-symmetry FEEC stress elements: awaits
  the simplicial H(curl)/H(div) family bead (dcng); B-bar is the
  shipped, measured locking-free path with the same acceptance metric.
- 3D and shells; higher-order families (Qk/Pk, k ≥ 2); IGA frontend
  (bead tfz.9).
- Hyperelasticity on the CutFEM frontend (linear ships there; the
  finite-strain cut path composes this crate's Newton loop with
  fs-cutfem quadrature in a successor bead).
- The legacy CutFEM facade makes no independent coercivity,
  near-incompressible, quadrature, graded-tree, convergence, or
  cancellation claim. In particular, material ratios beyond
  fs-cutfem's certified compressible regime are refused rather than
  routed to fs-solid's body-fitted B-bar formulation.
- fs-opdsl-generated residual/adjoint paths (the DSL coverage bead);
  the hand path here passes the same consistency gates the DSL output
  must pass.
- Contact no-claims: mesh-mesh point-triangle fallback, rigorous
  interval CCD over curved SDFs, adaptive barrier-stiffness policy,
  frictional adjoints/VJPs, dynamics/time stepping, and fs-scenario
  contact-law coefficient wiring. Plasticity flow remains successor
  scope (fs-material's plastic module wires in later).
- Ogden cards use the existing HyperTetProblem material interface unchanged;
  `tests/hyper3.rs::ogden_card_solves_uniaxial_bar` compares reactions and
  lateral motion against an independent compressible material-point root on
  three affine P1 mesh rungs (1e-8 relative reaction band). This establishes
  card consumption and the affine patch response, with no P1 locking claim
  or intrinsic soft-tissue calibration.
- Rod SE(3) DYNAMICS (fs-time symplectic wiring) and analytic rod
  tangents; 3D frames/lattices assembled from rods (topo-truss's
  consumer pass); J2 continuum element wiring (the fs-material card
  ships with algorithmic moduli; the element battery is successor
  scope); published El Centro fiber-model plot digitization (the
  shipped envelope is the capacity band + self-consistency).
- fs-scenario dimensioned-value resolution (units plumbing stays with
  the scenario consumer; kind/physics validation ships here).
- Arbitrary-scale sparse-fill guarantees, spectrum-completeness certification,
  LOBPCG preconditioning hooks, and higher-dimensional/nonlinear buckling
  sensitivities. Bounded sparse positive-factor pencils are supplied by
  `stability::sparse`. Fixed-load P1/Q1 element-modulus total gradients are
  supplied by `stability::sensitivity`, including its complete-packet
  `aggregate` owner, and consumed by the dedicated `fs-topopt` sizing study.
- Plate/shell buckling references (needs bending elements — the
  shells bead); the 2D canonical here is the strut.
- Kirchhoff--Love/bending-continuum claims,
  physical condition-number authority, curved shells, multi-patch/mortar joins,
  NURBS geometry, transverse-shear
  (thick/Reissner--Mindlin) behavior, nonlinear/finite-rotation response,
  loads, time integration, constraint multipliers for oblique support normals,
  shell error enclosures, and NAFEMS qualification. The `shell` operator is a
  flat, linear spatial base only; fs-iga's dedicated surface-shell frontend
  remains the required successor for those claims.
- Koiter coefficients versus LITERATURE tables (v1 validates
  classification + the sampled oracle; quantitative a/b benchmarking
  awaits the shell fixtures where the classic tables live).
- The reduce lab is Estimate-authority end to end until measured
  tissue/cane packs exist: the material card in the battery is
  synthetic-authored (no soft-tissue elastic pack exists in
  data/matdb — the recorded upgrade path), the loss factor is an
  authored card field, and the "lip-shaped" fixture is an authored
  tapered slab. The lab claims reduction correctness, never tissue
  truth; claim-side promotion happens in the per-track gates beads.

### Continuation of accepted 3D inelastic history

`SmallStrainTetProblem::solve_from` starts the original native increment/Newton
loop from supplied accepted nodal displacements and element histories.
`solve` supplies virgin history to that same implementation. Every nonlinear
trial uses the previous converged element history, and the method returns a
complete immutable candidate only after the requested load path converges.
Wrong array sizes and non-finite displacements refuse before history copying;
cancellation or a later failed increment leaves the supplied state reusable.

The new problem's load factors are absolute multipliers of its own prescribed
displacements and dead forces. A continuation may release constraints, change
loads or rebind per-element coefficients. The caller owns mesh/material lineage
and law-specific history and parameter-work admission: this generic numerical
entry point neither transfers constitutive families nor certifies a thermal
coefficient change. It retains the existing small-strain and P1 limitations.

Focused G1/G3 regression sources cover split-path equivalence through reversal
and reloading, closed-form permanent strain after releasing the driven end,
elementwise yield reduction at fixed displacement, and retry after a later
increment refusal without consuming the accepted history.

## Quadratic thermal prestress and vibration

`linear3::QuadraticPrestressedAssembly` extends the reference initial-stress
approximation to TET10, HEX20 and HEX27. Each problem's
`solve_thermal_tangent` or `solve_thermal_tangent_with_nodal_forces` first solves
its source-bound thermal/dead-load equilibrium and returns that immutable response
with the existing `TetElasticAssembly` vibration pencil. An independently supplied
stress array or stale solution cannot enter this API.

TET10, HEX8 (Standard/B-bar), HEX20 and HEX27 also expose
`solve_thermal_response_with_prescribed_dofs` and
`solve_thermal_tangent_with_prescribed_dofs`. Sorted, unique full xyz prescriptions
specify absolute support displacement in metres on already-fixed components;
omitted components stay zero. The shared owner solves
`K_ff u_f = f_f + f_thermal,f - K_fc g` with a matrix-free lift using the original
quadrature, local `C(T)` and strain projection. Recovery uses the complete motion,
retaining physical internal-minus-applied reactions, elastic storage and signed
initial-stress stiffness. Static support motion leaves homogeneous incremental
vibration constraints and reference mass unchanged. Existing residual, strain,
source and geometry checks remain active. `tests/quadratic_thermal_prescribed.rs`
checks support/thermal superposition across all element families, pointwise
material variation, reactions, prestress work, warped anisotropic B-bar response
and invalid support declarations.

The correction is
`K_g[(a,i),(b,j)] = delta_ij integral grad(N_a)^T sigma_0 grad(N_b) dV`.
The same isoparametric map and stiffness quadrature evaluate the actual
`sigma_0 = C(Bu - epsilon_free)` at every point, through the shared response
kernel. Element-average stress is not substituted into this product. The
original reference consistent mass, geometry, material state and free-DOF order
remain unchanged. The pencil `(K_e + K_g, M)` feeds the existing modal/spectral
owner directly; its original elastic and signed geometric matrices remain
separately inspectable.

Tension stiffens and compression softens the relevant incremental motions.
Neither negative stiffness nor nonpositive eigenvalues may be clamped,
absolutized or presented as real vibration frequencies. This symmetric pencil
does not itself certify stability. It omits finite-strain kinematics, follower
loads, geometry-dependent boundary loads, live modal state transfer and
post-buckling. Fully fixed bodies have a valid static thermal response but no
vibration pencil and return `NoFreeDofs`.

Existing assembly, source-state, rigid-support, residual and pointwise
small-strain gates precede publication. Element/quadrature and sparse-row
checkpoints preserve cancellation. Identity binds the admitted equilibrium and
the complete signed stiffness/mass pencil.

G1/G3/G4 definitions in `tests/quadratic_thermal_tangent.rs` use all three
families to check zero-stress recovery, signed uniform thermal prestress,
unchanged mass/support maps, spectral shifts, and a manufactured varying-stress
quadratic field whose geometric energy differs from an element-average shortcut.
They also check the fully constrained and cancellation boundaries. These are
numerical approximation checks, not experimental material validation.

## Generic initial strain on the native tetrahedral owner

`linear3::eigenstrain` admits source-bound mesh-frame Mandel initial strain
without manufacturing a temperature. Principal strains rotate through the
existing second-order Mandel transform with explicit nonzero source/mesh frames.
Generic and thermal inputs share the original equivalent-load, PCG displacement,
stress recovery, geometry admission and prestress kernels; thermal source types
and thermal identity formats remain intact. Generic output carriers now name
the shared eigenstrain response; existing thermal names denote those same results.
`SolidStiffnessModel` is the existing source schema and lowering shared across
thermal and other evolving material conditions. No duplicate constitutive or
linear solver is introduced. Deformed charts expose their admitted tetrahedral
volumes for mass-conserving current density.

Generic initial strain is checked against exact current elastic identities;
stress/tangent recovery refuses a stale material/source/load/solution. Fully
constrained equilibrium returns zero displacement and actual reactions; a modal
pencil still requires free DOFs. The tangent uses the existing signed
reference-configuration initial-stress approximation, not a finite-deformation
or post-buckling theory. Input source work is external to this mechanical energy
account; no combined chemo/thermo-mechanical passivity claim is added. G1 tests
cover free directional deformation, constrained stress and energy, thermal/generic
numerical equivalence, frame rotation, stale input and cancellation.


## Source viscoelastic pressure-valve dynamics

The material-valve path in fs-solid::reduce::visco connects the existing
TetViscoField / TensorProny owner to the ordinary valve reduction, pressure
ports and BeatingReedSource interface. Uniform or per-element anisotropic
Kelvin-Voigt/Prony data supplies density, storage and loss. The same TET4 or
explicit-quadrature TET10 mesh, constraints, modal shapes, face normalization
and pressure work projection are used throughout. No material-name branch or
independent damping constant is inserted.

This is an explicitly admitted frozen-frequency reduction:
K = K'(omega_ref), C = K''(omega_ref)/omega_ref. All modes share ONE storage
pencil and its mass-orthogonal basis. The existing fs-modal damping certifier
examines the full projected symmetric C; the legacy `MaterialValve` refuses
a required coupled block instead of exporting independent decay rates. Full retained storage
and loss matrices are also compared with the frozen model at EVERY retained
modal frequency. Controls explicitly bound the resulting storage-energy and
absolute loss-factor operator errors, the modal loss factor and allocation.
These are finite-subspace, discrete-frequency checks. They do not establish
broadband/hereditary accuracy, coupling to omitted modes, nonlinear mechanics,
or a thermally coupled constitutive history.

The static opening and pressure comparisons still use the source EQUILIBRIUM
stiffness, not K'(omega_ref). Thus even a complete retained frozen model must
disclose any lost equilibrium compliance. Conversion to the legacy one-mode
playing image must additionally pass its existing playing-mode admission. Every law must admit DC, the reference frequency and each modal
sample. The wrapper retains the original field and approximation report;
identities bind constitutive tensors, branch times, density, source identities,
field ordering and the selected reference frequency.

Authored loss must be zero on this source-material path. The existing card
codec preserves each mode's actual damping_n_s_m, and reed lowering reads that
coefficient as c/(2*sqrt(k*m)). The legacy loss_factor field denotes only the
authored contribution; it is not total material loss. Exact legacy eta/2
rounding is retained only when the authored formula reproduces that mode's c.
The same source loss consequently changes the existing mechanical runtime's
energy decay after serialization, rather than merely changing metadata.

The focused visco_valve integration target covers analytic single-coordinate
mass/stiffness/pressure normalization, source-viscosity changes to actual
port-Hamiltonian decay, Prony changes to frequency and damping, equilibrium
compliance disclosure, codec-to-reed continuity, uniform/per-element identity,
and refusal of incompatible bands, double-counted loss, excessive dispersion
or modal loss requiring a coupled runtime.

### Complete selected material valve modes

`reduce::visco::reduce_modal_pressure_loaded_valve` explicitly exports a
`MaterialModalValve` for a coupled runtime. It reuses the same source assembly,
modal window, frozen-frequency dispersion checks and equilibrium comparators,
retains the original `ViscoElements`, every selected mass-normalized mode and
the complete projected viscosity, and permits the certifier's full blocks.
A mode with zero mean opening remains available to carry pressure work,
external force or coupled material loss. The legacy card and one-mode path
retain their refusal; no independent-mode conversion is exposed by this owner.

Signed mean-opening and pressure projections stay separate. Nodal force
projection and nodal displacement/velocity expansion use the same original
shapes and support map. Static diagnostics include all selected modes while
preserving the original equilibrium reference. Identity binds the complete
basis, full viscosity and signed ports. This is still an admitted frozen K/C
image, with no claim to broadband hereditary, omitted-mode, contact or
large-strain accuracy.

The added `visco_valve` regressions exercise an exactly zero-opening,
pressure-carrying mode in an anisotropic two-coordinate tetrahedron, full-C
work against the original element assembly, virtual-work reciprocity, legacy
refusal, source-dispersion admission and cancellation.

### Hereditary complete material valve modes

`reduce::visco::reduce_hereditary_pressure_loaded_valve` exports a separate
`HereditaryModalValve` from the original equilibrium pencil. It retains all
selected mass-normalized modes, including zero-opening modes, their signed
pressure/force observations, original tetrahedral laws and complete Kelvin C.
The same `ViscoElements::project_basis` supplies every quadrature strain map,
physical volume and tensor-Prony factor. This path does not freeze storage or
loss at a reference frequency, use modal loss ratios, or loosen the legacy
frozen/independent-mode admission.

Physical coordinates satisfy q=sqrt(total_mass)*x. Explicit projection and
factor-work caps apply before workspaces are allocated. The projected mass
and equilibrium stiffness must match the original modal pencil within the
caller-declared basis error. DC and the conservative complete instantaneous
frequency bound must lie in the original material bands.

Bit-identical relaxation times share one existing Householder-QR row reduction.
Column equilibration and original Gram reconstruction preserve all coupled and
weak physical directions without a relative rank cutoff. At most n channels
remain per distinct source time, independent of the number of mesh points.
This is energy/work/heat equivalent for wholly relaxed or unrelaxed initial
histories driven by the same retained coordinates; arbitrary independently
authored quadrature memories and per-element heat reconstruction are not
provided. Original material strain limits and maps remain attached and bound
into image identity. A conservative Mandel-norm bound admits ordinary runtime endpoints in O(n);
near a limit the original quadrature maps decide. Endpoint strain admission
is not an interval enclosure.
Source/modal/spatial convergence and sampled coupling accuracy remain separate.

Regressions compare full frequency-dependent matrices with the original
sparse material assembly, and compact-channel stress, energy, work and heat
with the original tensor history under signed shear ramps and holds. They
also cover strong dispersion refused by a frozen fit, dark modes, very weak
loss directions, source bands and explicit capacity refusals.

## Source tetrahedra projected into full-memory structural dynamics

`linear3::visco::ViscoElements::project_basis` projects an explicit displacement
basis using the original Tet4/Tet10 and Standard/B-bar HEX8 strain matrices, physical quadrature volumes,
material laws and constrained free-DOF map. The resulting
`ProjectedViscoSolid` retains the complete consistent mass matrix
`Phi^T M Phi`, each original quadrature-point tensor constitutive contribution,
and its owning element. It converts Mandel shear to tensor shear at the existing
`TensorProny::project_strain` boundary. All Kelvin cross terms and all original
Maxwell factor rows survive; no mass lumping, material averaging, fitted damping
or silent basis truncation is introduced.

The same retained basis maps nodal loads by its transpose and expands
displacement/velocity with the original supports. These maps preserve discrete
virtual work; fixed rows contribute no work. `SolidProjectionBudget` bounds
coordinate count, stored basis entries, quadrature points and scalar projection
work before constructing the result. Overflow, nonfinite maps, source refusal
and cancellation return errors without a partial projection. The immutable
projection borrows its original source assembly.

B-bar projections also retain each original compatible quadrature strain map.
The shared material owner checks it independently of the assumed work strain,
including after thermal eigenstrain or source-law rebuilding. This permits the
same `TensorSolidBody` and prepared material dynamics consumers without hiding
local overstrain or substituting compatible stiffness for the B-bar operator.
The extra geometric projection work is charged to the original scalar budget.
Focused HEX8 tests compare projected stiffness and consistent mass to the native
assembly, reject a strain hidden by B-bar averaging, and retain that refusal on
material rebuilding. The coupled-body regression checks actual motion and
element heat against an independent damped-oscillator solution.

A dynamic consumer must independently admit positive-definite, representable
projected inertia. Geometry, supports, basis and material state remain fixed.
This API does not claim spatial/modal convergence, arbitrary frequency validity,
nonlinear strain or thermally evolving geometry. Focused G0/G3 integration
regressions in `fs-couple/tests/tensor_solid.rs` cover Tet4/Tet10 shear and
consistent mass, nodal virtual work, nonorthogonal constrained bases, original
Maxwell prehistory and regional material-loss swaps. They use synthetic material
fixtures and do not establish empirical calibration.

## Continuing thermal-valve motion across source updates

`reduce::thermal::motion::ThermalValveModalMotion` owns the actual thermal
reduction and the existing `fs-modal` time-domain state. The reducer retains
its original mass-normalized mode shapes and exact pressure-load projections.
Pressure steps use the original exact ZOH owner; applied work is conjugate to
swept volume, and nodal observations reconstruct equilibrium plus vibration.

For affine Tet4 cells,
`reduce::thermal::reduce_pressure_loaded_valve_with_nodal_reference_density`
uses the exact original P1 density field in the native consistent mass and
initial-stress pencil. It retains that complete field on the thermal valve;
neither a cell mean nor endpoint constitutive density substitutes for it.
Quadratic requests on this entry point refuse explicitly. The original
scalar-density reduction keeps its existing arithmetic.

Candidate source updates preserve the exact sample clock and use the shared
mass-inner-product transfer on the same reference mesh, frame, supports and
consistent mass. Explicit P1 profiles must also agree at every original vertex:
equal reduced matrices do not identify density redistributed at constrained
nodes. The G1/G3 graded-inventory regression in `tests/thermal_valve_motion.rs`
checks exact cubic mass moments, uniform-field parity, thermal momentum
continuity and refusal of a density change hidden by the reduced matrix. New equilibrium, stiffness, mode count/order and supplied
pressure observations are admitted with explicit state-loss and signed
incremental-storage budgets. Updates and pressure steps publish only complete
candidates. Discarded state is not damping or heat; bore feedback, radiation
reaction and reciprocal thermoelastic heat are not inferred. Focused
one-tetrahedron tests use independent consistent mass, stiffness, thermal
equilibrium and oscillator references to check motion/work, source-update
continuity, changed pressure response and refused/cancelled candidates.


### Native finite-strain storage and reference-mass access

The `hyper3` owner exposes full referential dead loads, actual stored energy and
a consistent reference mass for dynamic consumers. Energy uses the same
deformation gradient, determinant/stretch admission and constitutive owner as
the native Piola force and tangent. Reference mass delegates to
`linear3::p1_tet_consistent_mass`, with one positive finite density per element;
unrepresentable element/assembled mass refuses. Consumers must retain that
mass during deformation instead of multiplying updated volume by a fixed
reference density. These accessors introduce no constitutive or time solver.

`FinitePlasticTetMotion::try_new_with_nodal_density` also admits an authored
positive P1 reference density on the original shared vertices. Its consistent
matrix integrates rho*N_a*N_b analytically through degree three and retains
the complete nodal inventory across the existing finite-plastic momentum
steps. The cell-density accessor reports cell means; the separate nodal
accessor preserves the actual distribution. `reference_acceleration_load_n`
uses this SAME mass action to project held nodal accelerations, including
uniform gravity, into the existing nodal-load input. No current-volume mass
refresh or one-quarter-per-node lumping is introduced. Uniform nodal density
delegates to the original cellwise mass path. Independent simplex moments,
rigid kinetic energy and native backward-Euler work checks cover the mass
operator; the finite-plastic regression exercises the native free-body Newton
step. This is not material-source inference or a new nonlinear-gravity
quadrature, deformation-field solver, or thermal mass-transfer model.


### Impulses with complete finite constitutive memory

`FinitePlasticTetMotion::impulse_candidate` and
`FiniteViscoTetMotion::impulse_candidate` share the original `ReferenceMass`
owner for zero-duration momentum transfer. Native checked CG solves
`M_ff delta_v_f = J_f` through borrowed sparse rows with diagonal
preconditioning; it builds neither a dense inverse nor a second mass matrix.
Both cellwise and graded nodal reference density retain their exact existing
mass. Caller-supplied held DOFs keep their current velocity, including nonzero
support motion. Repeated applied DOFs accumulate in input order.

`SolidImpulseReport` resolves linear and angular momentum at the actual
deformed positions, applied and support impulses, kinetic change, and signed
applied/support work. Physical gates use the actual representable velocity
change and include each free-DOF residual, so cancellation of net impulses
cannot conceal lost local momentum. Kinetic change is evaluated through
`average(v) dot M delta_v`; the work residual uses the independently
recomputed impulse defect instead of subtracting large work totals.

Only velocity changes. Original deformation, plastic metrics and accumulated
flow, every Maxwell branch, stress, storage, material time, law bindings and
reference inventory survive the event. Candidates publish only after solve,
workspace, momentum, work and cancellation gates. Subsequent native motion
continues the accepted material history. Impulses do not themselves generate
relaxation/plastic heat or select a contact/restitution law. The caller owns
the opposite impulse and any contact energy loss; mixed incompressibility
constraints are outside these displacement-only P1 entry points.
The `hyper_plastic_dynamics/impulse.rs` and `hyper_visco_dynamics/impulse.rs`
regressions use independent cubic-simplex mass, deformed torque and support-work
oracles, retained inelastic prehistory and subsequent native motion.

`FiniteThermoplasticTetSolid::impulse_candidate` extends this same account to
an already inertia-bound thermal tetrahedral solid; a quasistatic owner without
original reference inertia refuses. `FiniteThermoplasticHexMotion` supplies
the same event for both standard and mean-dilatation Hex8 formulations, using
a separate 3x3x3 reference-mass rule, exact for warped trilinear bricks, and
physical current positions. The eight material quadrature histories are
unchanged by this mass rule.
The thermal family bindings and reference temperatures, current temperature,
all plastic/creep and quadrature history, physical/assumed deformation, stress,
free energy and internal energy remain unchanged. The event advances no
material interval and deposits no heat. `tests/thermoplastic_impulse.rs`
checks hot plastically set bodies, independent reference-mass/torque/work
oracles, moving supports, state retention and subsequent native thermal motion.

### Graded reference inertia in finite thermoplastic tetrahedra

`FiniteThermoplasticTetSolid::with_nodal_inertia` binds a positive P1 reference
density through the original exact consistent-mass owner. The complete nodal
field and mass operator remain fixed through thermal softening, plastic/creep
history, finite deformation and instantaneous impulses. Reported cell means
describe inventory; they do not replace the graded inertia. Existing cellwise
binding retains its original arithmetic, including the uniform-density reduction.

`reference_acceleration_load_n` applies that same mass operator to an
original-node P1 acceleration field, retaining supported-node forces and
density-weighted moments. The caller supplies these held loads to the native
step once and owns any evolving field/potential. Heating does not requery
density to change the mass of a closed specimen, and numerical kinetic losses
remain separate from caloric heat. The retained field is available through
`reference_nodal_density_kg_m3`.

The focused `hyper_plastic_thermal/nodal_density.rs` regressions specify exact
simplex mass/gravity moments, heated plastic continuation, equal-total-mass
redistribution of contact mobility, impulse/history retention, the uniform
reduction, and atomic refusal/cancellation. This authored reference inventory
adds no new material-source or pressure-dependent density claim.

### Material-point mobility with inelastic memory

The finite plastic and viscoelastic motion owners expose `TetImpulsePoint`,
`point_kinematics`, `point_directional_inverse_mass_per_kg` and
`point_impulse_candidate`. One admitted P1 barycentric trace B determines
current position and velocity; its exact transpose distributes a world
impulse. Nonnegative unit-sum weights admit interior, face, edge and vertex
points, with only roundoff-sized sum normalization. Extrapolation refuses.
The directional mobility is `n^T B M_ff^-1 B^T n` for a unit world direction,
using the original cellwise or graded reference density and held velocity DOFs.
It borrows the existing checked-CG mass response without cloning constitutive
history or forming a dense inverse. A fully held trace has zero mobility.
Only the accepted impulse candidate copies the complete unchanged material
state. The native visco impulse regression checks an interior point against
independent cubic-simplex mass inversion, including negative consistent-mass
nodal responses, supported mobility, deformed torque and conjugate point work.

The same owners expose `point_mobility_matrix_per_kg` and
`point_impulses_candidate` for simultaneous contact. Every row/column pair of
`B M_ff^-1 B^T` retains its full coupling through the original supported mass.
One checked response is retained at a time; no dense mass inverse or material
history copy is needed to build mobility. Caller row and workspace limits
bound the matrix, trace maps and solve reports before allocation. The returned
matrix is not silently symmetrized or clamped; the contact-law owner admits it
with its own numerical tolerances. All point loads reach the compensated nodal
accumulator together and produce one native velocity, support and work account.

Hot Tet4 and Hex8 motion expose the same five point operations. Tet4 uses the
original barycentric trace; `HexImpulsePoint` supplies reference coordinates
in `[-1,1]^3` for the original trilinear trace. Physical
deformed coordinates determine position, velocity and torque even when the
Hex8 constitutive owner uses mean dilatation. Original density, thermal source
families, temperatures, quadrature histories and free/internal storage survive
unchanged. Tet4 requires its reference inertia to be bound first. Focused native
regressions compare complete two-point matrices with independent simplex and
tensor-product inverses, including oblique directions, nonvertex points,
supported responses, summed work and unchanged thermal/material state.

### Exact-incompressible pressure impulses

`mixed::conservative::ConservativeTetSolid::impulse_candidate` solves the
original P2/P1 velocity and pressure-impulse KKT equations at the current
configuration. The endpoint cofactor supplies the same weak volume derivative
B as the native mixed geometry. Original consistent mass supplies M; native
sparse LDL handles `[M_ff B_f^T; B_f 0]`. Explicit velocity, pressure-impulse
and work scales preserve a symmetric dimensionless system, with assembly and
factor budgets and a recomputed true residual. No second mass solve occurs.

The multiplier is an impulse in Pa s, separate from the stored algorithmic
pressure. Held velocities remain unchanged; full support impulses and work
remain explicit. The report separates external, internal pressure and support
contributions. Actual represented velocities must satisfy every original weak
volume-rate row, including a removed gauge row. Pressure contributes neither
net linear/angular momentum nor admitted work beyond explicit tolerances.
A gauge is allowed only for a verified connected constant-pressure null mode;
other rank deficiencies refuse. Incoming incompatible velocity refuses rather
than silently receiving a dissipative projection. In particular, an earlier
volume-constrained time step does not automatically guarantee endpoint tangent
velocity admissibility. The impulse preserves geometry, density, material
points, stored pressure and time, and adds external/support work to the
existing conservative trajectory ledger. `tests/mixed_conservative.rs` checks
independent P2 mass moments, pressure loading, rotated geometry, supported
gauge reactions and continuation. The exact-incompressible branch uses this
pressure solve; finite-compliance conservative formulations retain their
original mass-only response and supply no pressure impulse.

The shared native impulse owner now retains a prepared mass/pressure response.
It assembles the endpoint gradient and its transpose and factors the scaled KKT
system once. Every subsequent mobility column solves `B delta-v = 0` directly;
it never estimates an increment by subtracting large background velocities.
The original free-row residual, weak rate and pressure zero-work checks apply
to each response. The final event retains the original incoming-rate correction,
actual rounded endpoint checks and complete pressure/support account. One
factor may be reused only while its borrowed accepted geometry and supports
remain unchanged. No constitutive trial or history copy occurs for mobility.

Independent rational P2/P1 elimination fixtures check cross-point and cross-axis
responses, annihilation of pure pressure loads, retained factor reuse and the
pressure multipliers recovered by the original event. For a unit-mass affine
tetrahedron, centroid xx mobility is `17735/6528`, distinct from the mass-only
`95/32`; this difference participates in the contact law rather than being
deferred to a post-impact projection.

### Reference-aware finite-deformation response

`hyper3` evaluates energy, Piola force and tangent through the original
`HyperelasticMaterial` interface, preserving source reference-deformation
adapters. Bare laws retain the default native response. The actual total
deformation must have a finite positive Jacobian independently of the
material's elastic-gradient admission. Reference geometry and the consistent
reference mass remain fixed; a thermal adapter does not alter the mass mesh.

### Source-defined cohesive load paths

`cohesive::CohesivePatch::from_interface` consumes
`fs-material::cohesive::source::ResolvedBilinearCohesive`, checks the actual
ordered surface material assignments, and retains the original source bundle
through every accepted displacement/damage state. `p1::MaterialFacePair` and
`P1Interfaces::from_interfaces` lower the same source through the existing
matching-face quadrature and work-conjugate maps; both full and reduced patches
retain their source. No constitutive, quadrature or nonlinear solver is copied.

Swapping an interface dataset changes the traction/separation envelope, failure
point and area-integrated fracture work in the native equilibrium solver.
Conditions remain fixed during each mechanical increment. The explicit
`cohesive::temperature` transfer below advances the same sourced interface to
new temperatures without resetting history. Chemical evolution, adaptive
fracture topology and complete structural/thermal energy accounting are not
introduced. The `cohesive_source` integration tests
check an independently solved compliant pull, material substitution, exact P1
area/fracture accounting, surface-order refusal and atomic cancellation.

### Nonmatching bonded P1 surfaces

`cohesive::p1::nonmatching` binds independently triangulated planar surfaces
through their triangle intersections. Each intersection is triangulated and
integrated with positive degree-two quadrature. Distinct barycentric weights
on the two sides build one shared jump map; the existing cohesive operator uses
its transpose for forces and its consistent tangent for Newton iterations.
Undamaged linear traction is integrated exactly to floating-point precision,
and nonlinear damage retains discrete virtual work and force balance. Moment
and rigid-motion consistency hold on the common plane; any admitted off-plane
projection contributes error bounded by the supplied geometric tolerance.
The source-bound constructor retains the original ordered interface source at
every full and reduced quadrature point.

Construction checks exterior-face incidence, disjoint side DOFs, opposite-side
bodies, face quality, planarity, within-side overlap and full area coverage on
both sides. Work and storage limits bound face count, pair checks and quadrature
points, with cancellation throughout. A deterministic bounding-box sweep skips
disjoint pairs before they consume the pair budget. Worst-case preprocessing
remains quadratic for heavily overlapping boxes. Runtime
solves reuse the fixed sparse maps and original constitutive/solver path.

This is an explicit initially bonded, fixed small-strain reference interface.
The caller admits the relative common-plane projection/coverage tolerance.
Curved interfaces require separate locally planar regions; large sliding,
dynamic rebonding and adaptive fracture topology are not inferred. Nonlinear
softening still requires independent surface and load-path refinement. Tests
cover crossed-diagonal area and analytic linear energy, virtual work and rigid
motion, matching-mesh stiffness recovery, force-controlled softening/fracture/compressive
reclosure, source retention, malformed coverage and construction budgets.


### Exact native transient continuation across material lifetimes

`loads3::transient::checkpoint::TetTransientCheckpoint3` detaches accepted native
motion from its borrowed material problem. It retains original geometry, shared
assembled operators, motion/acceleration, thermal field, applied loads, support
state, controls and accepted step count. Binding requires the same mesh, fixed
components and complete material values/source identities. It neither
reassembles nor restarts dynamics; changed coefficients still use the existing
thermal/rebind owners. This immutable in-process continuation is not portable
serialization. G3 continuation, step-budget preservation and changed-material
refusal are exercised by the reciprocal dynamic thermoelastic tests in
`fs-couple`.

## Consistent reference inertia for planar frame vibration

`frame2d::FrameMassFormulation::ConsistentReference`, selected by
`Frame2d::with_mass_formulation`, integrates linear axial and cubic Hermite
transverse centerline motion on each original Euler–Bernoulli member. The
result retains translation/rotation and endpoint cross terms; separately
added point masses remain diagonal. `Frame2d::new` retains the original
lumped default, and a model with only nodal mass retains the diagonal path
even when consistent element inertia is selected.

Cross-section rotary inertia uses the same owner: material members derive
`sum(rho I)` from each original layer's density and centroidal area moment,
including the layer's own material state. `add_element_with_rotary_inertia`
accepts an explicit section mass moment per length [kg m] for raw-coefficient
members; `add_element` retains its declared zero-section-spin approximation.
Consistent mass integrates `rho I (d(v_dot)/dx)^2` exactly with the Hermite
slope (Rayleigh beam kinematics), including its endpoint and translation-spin
cross terms. Endpoint lumping assigns half the section spin inventory to each
node. Neither path counts centerline orbital inertia a second time.

`FrameInertiaProfile` and `add_element_with_inertia_profile` retain original
line mass [kg/m] and section spin inertia per length [kg m] at both
member endpoints in i-to-j order. `with_first_moment` also retains the signed
transverse density moment `integral rho*y dA` [kg]; positive y is the original
member's left normal, and spin is measured about its declared centerline.
Each moment is interpolated linearly on the original reference member.
Endpoint positive-semidefinite moment matrices enforce `first² <= mass*spin`
throughout the member. The centerline centroid remains the default when no
first moment is supplied. Consistent inertia uses closed-form endpoint-weighted
Gram matrices
for the existing linear axial, cubic Hermite transverse and quadratic slope
velocities. These integrate the complete degree-seven density/velocity and
degree-five spin/slope products exactly, without runtime quadrature or replacing
the field by a member mean. The signed term `-2*first*u_dot*v_dot'` integrates
the actual Euler--Bernoulli velocity `u_dot-y*v_dot'` with exact degree-four
products and no additional runtime quadrature. The resulting physical rigid
mass, first moment
and angular inertia retain the original density placement.

Uniform centered profiles dispatch to the original constant-coefficient arithmetic.
Explicit lumping uses the positive partition `L*(2*value_i+value_j)/6` at
endpoint i, preserving total inventory and its first axial moment. A
nonzero transverse first moment requires consistent inertia; diagonal lumping
refuses because it cannot retain the translation/rotation coupling. The complete
profile is immutable member data, available through `reference_member_inertia`;
cloning, constitutive history and motion retain it with the original geometry.
Profile and accumulated-inventory refusal occurs before the frame is changed.
`tests/frame2d_graded_mass.rs` checks independent rigid spatial moments under
rotation and endpoint reversal, polynomial kinetic energy, equal-total-mass
native modes and supported impulses, loaded/base-driven Newmark response under
time-step refinement, positive lumping, uniform-path bits and refusal atomicity.
Its eccentric-section cases independently integrate rigid spatial moments and
polynomial kinetic energy, check signed rotation from an axial impulse, and
reject nonphysical moments and a diagonal representation of eccentric inertia.

One reference mass owner supplies the generalized vibration pencil, the
initial acceleration solve, every Newmark inertial residual/tangent, Rayleigh
mass damping, uniform-base excitation, kinetic energy and support reactions.
Massless free coordinates still use the original static stiffness condensation;
the full massive block is reduced with the existing fs-la Cholesky owner before
the existing symmetric eigenvalue solve. Mass is prepared once per analysis,
not rebuilt within Newton iterations. `inertial_forces` applies the same M to
an explicit full acceleration vector, including consistent gravity/body loads
and their generalized nodal moments. Reaction records include supported inertia,
damping and applied loads, including the initial endpoint.

This is fixed reference-configuration inertia for linear or moderate-rotation
vibration. It does not include configuration-dependent mass, centrifugal or
gyroscopic terms, or shear deformation. Distributed body loads are
projected through the reference interpolation; one force-based element does
not certify spatial convergence of dynamic stress or distributed-load section
recovery. The existing dense building-scale solver and constitutive/kinematic
limits remain. Selecting a mass model does not alter the material law.

`tests/frame2d_mass.rs` checks exact rigid translation/rotation inertia in
rotated frames, two cantilever frequencies under mesh refinement, free-body
base acceleration and work, full supported inertia reactions, a loaded
two-coordinate history against independently solved generalized eigenvalues,
and unchanged nodal-only behavior. The one-element first bending-frequency
error falls from about 30% for endpoint lumping to below 0.5% with consistent
mass; higher-mode accuracy still requires member refinement.
Additional G0/G1 tests compare section spin against independent slope
quadrature, exact rigid-body rotational inertia and the analytic simply
supported Rayleigh spectrum under mesh refinement. Source-layer tests check
the full density-weighted spin inventory and its effect on the original
Newmark trajectory.

## Instantaneous native frame impulses

`Frame2d::apply_impulse` applies a full node-major generalized impulse through
the same reference mass used by Newmark: `M_ff * delta_v_f = J_f`. The original
lumped path stays diagonal; the consistent path factors the free mass once and
retains all endpoint, translation and rotation cross terms. Fixed coordinates
keep zero velocity. The returned `FrameImpulse::support_impulse` is the impulse
exerted by the supports **on the frame**, `M * delta_v - J` on fixed coordinates.
Supported consistent-inertia coupling is retained, including when the applied
impulse is entirely on free coordinates.

The event reports the actual represented velocity and generalized momentum
changes, endpoint kinetic energies, signed impulse work and kinetic-energy
change. Work uses the mean of the two endpoint velocities; the energy change
uses that mean dotted with `M * delta_v`, avoiding subtraction of large nearly
equal kinetic energies. Each free momentum equation and the work balance must
pass a scaled floating-point residual check. A jump lost when adding to a large
existing velocity is refused rather than reported as applied.

The caller supplies the complete finite `post_event_load`, Rayleigh coefficients
and the real `continuation_dt`. Any uniform-base inertial force belongs in that
full load vector. The method reequilibrates finite right-limit acceleration
from `M * a_plus = post_event_load - C * v_plus - f_internal` using the original
source viscosity and Rayleigh operators. `continuation_dt` admits those sources
against the following sampled clock; it never scales the impulse or advances
the material clock. Positions, section deformations and material histories are
unchanged. Velocity, acceleration and any newly captured Rayleigh reference
publish together only after every check succeeds.

All free coordinates must have positive inertia. Free massless coordinates
require a DAE velocity/acceleration projection that this event does not provide,
even when their applied impulse entries vanish. Fixed massless coordinates and
all-fixed support-only events are admitted. Timed constitutive laws retain the
existing untimed-operation refusal. This is native generalized momentum under
the frame's reference-inertia assumptions; it does not certify finite-rotation
world angular momentum or reconstruct an impact contact law.

`tests/frame2d_impulse.rs` checks independently integrated Hermite mass coupling,
support impulse, signed work, finite-force and Rayleigh right limits, Newmark
continuation, rotation covariance, lumped section spin, both all-fixed mass
formulations, massless refusal and failure atomicity. The material-viscosity
battery checks the original source loss in the event acceleration, independence
from admitted continuation clocks and refusal of inadmissible clocks.

### Simultaneous normal impacts on the native frame

`Frame2d::apply_normal_impacts` resolves one declared, fixed-configuration
frictionless contact set through `fs-contact::multicontact::solve_normal_impacts`.
Each `FrameNormalContact` supplies the derivative of its physical signed gap in
the original node-major coordinates and a prescribed obstacle normal speed.
The same covector maps velocities to gap rates and contact impulses to nodal
forces and moments. Translational weights are dimensionless; rotational weights
have metres. Node contacts use their unit normal, linear-reference member
contacts can use the original Hermite `reference_transverse_trace`, and internal
contacts include both signed traces in one row. Geometry and event discovery
remain caller-owned; penetration is refused without position repair.

The original free mass supplies every entry of `G M_ff^-1 G^T`. One native mass
factor serves streamed mobility columns, the combined impulse and finite
right-limit acceleration. Additional adapter scratch is `O(ndof + contacts^2)`;
no retained inverse, frame clone or `ndof`-by-contacts response matrix is needed.
All geometry-active contacts participate together, including initially
separating contacts another impulse can drive closed. The shared contact count,
iteration, PSD/rank, complementarity and global restitution/passivity admissions
remain authoritative. A redundant result preserves `UniqueVelocity` and its
checked impulse representative; individual contact loads and local source-work
allocations do not acquire a uniqueness claim.

The contact solution enters the same private native impulse trial used by
`apply_impulse`. Actual stored endpoint gap rates, complementarity and relative
work must agree with the shared law before publication. The receipt separates
work supplied by prescribed obstacles, native kinetic change, contact loss and
any explicitly allowed numerical positive relative work. The original-mass
work balance is `kinetic change = obstacle work + relative work`; stationary
frame supports contribute their full cross-mass reaction impulse and zero
work. Relative work is measured directly with compensated extended products,
not by subtracting large world-frame kinetic or obstacle energies. Cancellation
or any late endpoint, work, force-balance or continuation refusal preserves all
accepted frame state, including a previously uncaptured Rayleigh reference.

The existing positive-free-inertia, untimed-material and fixed-reference-inertia
boundaries remain. Contact loss has no automatic thermal allocation. This is
an instantaneous event response; sustained contact, friction, moving-support
integration and finite-rotation world-inertia claims remain separate owners.
`tests/frame2d_normal_impact.rs` checks analytic coupled contacts, inactive
endpoints, Hermite moments, support cross mass, moving-obstacle work, internal
momentum and velocity-boost invariance, finite-force/Newmark continuation,
redundant determinacy and atomic physical/numerical refusals.

### Simultaneous planar Coulomb impacts on the native frame

`Frame2d::apply_frictional_impacts` resolves declared planar contacts through
the existing non-associated `fs-contact::frictional_multicontact` owner. Each
`FrameFrictionalContact` supplies normal and tangent covectors at the same
physical point, in an orthonormal contact frame, with prescribed obstacle
velocities, friction coefficient, stable ID and the native normal endpoint law.
The same rows map frame velocities to contact rates and impulses to nodal
forces and moments. Geometry, physical unit directions and event discovery
remain caller-owned; penetration is refused without repair.

The full original-mass contact mobility couples every normal and tangent row.
One free-mass factor serves two streamed physical columns per contact, the
combined native impulse and finite right-limit acceleration. The shared owner's
second tangent is an exact zero coordinate: zero mobility, velocity and accepted
impulse, without artificial inertia or compliance. This embedding conservatively
retains the native PSD velocity/work qualification; a chosen representative
never grants unique local loads or individual source-work allocations. The
shared count, nonlinear evaluation, iteration, rank, determinacy and passivity
limits remain authoritative. Adapter scratch is `O(ndof + contacts^2)`, and
contact-ID-ordered reductions preserve replay under input permutation.

Normal complementarity and maximum-dissipation Coulomb sticking/sliding share
one solve. Actual stored frame velocities must then reproduce both components,
close the original-mass momentum and obstacle-work balances, and satisfy the
same physical Coulomb residual gates before the existing private impulse trial
publishes. The receipt separates obstacle work, kinetic change, actual relative
work and total impact loss. Relative work uses extended products directly;
large world-frame energies are not subtracted to infer a small loss. Loss is
not assigned to frictional heat. Cancellation, budget exhaustion or any late
endpoint, work, continuation or force-balance refusal preserves all accepted
state, including a newly staged Rayleigh reference.

The original stationary-support, positive-free-inertia, untimed-material and
reference-inertia boundaries remain. No position, member force, constitutive
history or material clock changes during the event. Finite post-event loads,
source viscosity and Rayleigh damping determine the actual right-limit
acceleration and subsequent Newmark continuation. Sustained contact, collision
detection, finite-rotation world inertia and heat allocation remain separate
owners.

`tests/frame2d_frictional_impact.rs` checks independent sliding/sticking solutions,
non-diagonal normal/tangent coupling, consistent-mass support moments, canonical
ID replay, inactive endpoints, internal linear/reference angular momentum,
moving-obstacle work and velocity boosts, finite-force/Newmark continuation,
the zero-friction limit, redundant-contact qualifications and atomic refusals.

### Sustained planar contact through native Newmark motion

`Frame2d::advance_sustained_contacts` advances one held-force interval through
the original `newmark_with_load` response and the existing
`fs-contact::sustained::solve` force/displacement Coulomb owner. Fixed normal
and tangent covectors at each physical point apply forces and moments through
the same transpose that measures contact motion. Stable IDs order every
response and work reduction. The planar embedding retains an exactly zero
auxiliary tangent; no spring, compliance, impact kick or alternate integrator
is introduced.

Every nonlinear probe starts from the same accepted native frame. All original
responses, including rejected probes and the final fresh replay, spend the
shared response budget. The final replay supplies the actual frame/history;
it is not replaced by an uncounted reconstruction. The native 21-contact limit,
256-coordinate/member hard bounds, explicit adapter workspace and native
solver budgets bound work. Each original Newmark interval is indivisible for
cancellation; polling surrounds it and continues through the contact solve
and physical audits.

Acceptance measures actual stored displacement, velocity and finite endpoint
acceleration. Original full-mass momentum includes integrated support reaction
and consistent-inertia cross terms. Actual displacement must agree with the
Newmark midpoint-velocity increment using separate physical metre and radian
allowances, without an absolute-world-coordinate rounding floor. Lost resolved
motion cannot supply fictitious contact gap, slip or work. Compensated audits
separate elastic storage, original material/Rayleigh damping, held external
work, moving-obstacle work, signed normal work and relative friction loss.
The whole balance and original native energy ledger must both close; contact
loss is not added again to native damping or assigned a thermal destination.

The receiving path owner must explicitly admit the complete physical interval,
fixed traces, contact-feature coverage and contact mode before any state
publishes. Loaded contacts dwell; interior-cone forces stick; sliding cannot
cross an unresolved stop. Resolved endpoint slip reversal refuses even if net
slip retains its initial sign. Endpoint force convergence alone establishes
neither full-path admissibility nor global existence/uniqueness; any additional
branch or determinacy evidence belongs to that mandatory admission. Initial
closing impacts and contact-arrival transitions require separate event
resolution. Actual negative final gap is always refused.

Scope is frozen linear geometry, elastic sections, untimed source laws,
positive inertia on every free coordinate, stationary frame supports and
held obstacle velocities/loads. Original source-viscosity and retained
Rayleigh-reference rules remain authoritative. Automatic discovery, changing
contact frames, general finite-rotation inertia, massless DAE projection,
interior stick/slip event location and thermal allocation remain separate.
Any source, budget, geometry, mode, balance, callback or cancellation refusal
preserves the complete accepted frame. The callback's own state must likewise
remain staged until the complete call succeeds.

`tests/frame2d_sustained_contacts.rs` checks analytic gravity support,
sticking/threshold/breakaway/sliding, native unknown-force recovery,
detachment, moving-belt and lifting-obstacle work, stop-time refusal,
post-stop/native continuation, original consistent-mass support impulse,
elastic/damping work, actual small translated/rotated increments, large-coordinate
lost-motion refusal and atomic path/cancellation/resource/domain refusals.

## Signed, mass-normalized frame modes and participation

`Frame2d::modal_analysis` exposes the complete signed current-tangent spectrum,
full node-major mode shapes and horizontal/vertical uniform-base participation.
Shapes satisfy `phi^T M phi = 1`; fixed coordinates are zero and massless
coordinates use the original static condensation map
`phi_r = -K_rr^-1 K_rt phi_t`. The massless stiffness must be positive definite:
an unstable or unconstrained massless subsystem cannot silently disappear
through condensation. Each massless and mass factor is prepared once.
The original fs-la Cholesky and symmetric Jacobi owners supply the algebra;
the former frame-only eigenvalue kernel is retired.

With `r_x/r_y` the full unit translational influence, including supported
coordinates, participation is `Gamma = phi^T M r` and effective mass is
`Gamma^2`. The same reference M is used by the frame transient and body-force
owners. Supported inertia may also contribute directly to reactions, so the
modal sum need not equal the complete `r^T M r` inventory.

Every mode retains its unmodified signed eigenvalue, a freshly recomputed
free-pencil backward residual and a mass-whitened residual. The reported
sign-resolution diagnostic is the latter plus
`64*eps*n*||L^-1 K_cc L^-T||_inf`. It is a floating-point estimate, not a
certified interval or stability theorem. Modes are explicitly oscillatory,
unstable, or unresolved near zero; unresolved does not certify a neutral mode.
Resolved negative eigenvalues expose their linearized exponential growth rate.
`natural_frequencies` now requires every mode to be resolved and oscillatory;
rigid or unstable frames must use the signed analysis instead of receiving
clipped zero frequencies.

The analysis is frozen at the current admitted material/displacement state;
it does not locate a nonlinear buckling load or refresh during a history.
Repeated eigenvalues define a mass-orthonormal subspace, not unique vectors.
`tests/frame2d_modes.rs` covers exact massless beam-shape restoration and
participation, full-mass orthogonality and reconstruction, free rigid motion,
and an axially compressed column crossing its Euler instability.

### Material-defined frame members

`Frame2d::add_material_element` consumes the existing
`fs-material::section::ConcentricElasticSection`. It derives axial rigidity,
bending rigidity and reference line mass together from the original layer
integrals and retains the full section through `material_section`. There are
no independent stiffness/mass overrides on this entry point. Heterogeneous
member assignments and perfectly bonded concentric layers retain their own
material cards, query states and receipts; no average material replaces them.

The original force-based static, signed-modal, reference-mass, gravity and
Newmark owners consume those coefficients. A sampled frequency/rate law is
refused, as is a finite frequency-applicability band that this general time
history cannot enforce. Coefficients use the shared material policy for the
held typed temperature, pressure, composition, moisture and concentration;
unowned mechanical/rate response or applicability still refuses. The same
policy reaches static rods, shear frames, tensile links and rod motion while
retaining each owner's existing frequency restrictions. This constructs a new
stress-free reference specimen; it does not implement thermal expansion, live
mass exchange, plasticity,
interlayer slip or material-derived damping. Existing planar Euler–Bernoulli
and reference-inertia assumptions still apply. Missing data and source/geometry
admission failures publish no partial member or added mass.

`tests/frame2d_material.rs` exercises independent cantilever compliance,
`sqrt(E/rho)` frequency scaling, integrated weight, layered and spatial
assignments, and a guided-tip oscillator's actual Newmark history and energy
balance. `cargo run -p fs-solid --example material_frame -- 1000 2e9 4000 2e9`
demonstrates data-only substitution through the same geometry and solvers.
The example inputs and regression data are synthetic, not measured alloy
qualification.

### Physical bond rates in native cohesive equilibrium

`CohesiveState::solve_over_time[_controlled]` uses the same accepted history,
Newton/FGMRES owner, work-conjugate jump/scatter map, force gate and atomic
publication as the original equilibrium path. Explicit dt gives the physical
constant jump rate; source-selected bond viscosity contributes to total endpoint
traction and its complete consistent tangent. Untimed solves refuse viscous
bonds. Intrinsic cohesive work/storage/fracture and physical viscous loss are
reported separately; no structural stability or unique-root claim is added.

The existing restrained thermal joint passes its actual interval to both bulk
Maxwell and cohesive owners. Its total bond path work includes analytically
integrated viscous loss; endpoint-force quadrature corrections and the combined
energy audit use that same work. Bulk loss, intrinsic fracture and bond viscosity
stay distinct. Loss is not silently deposited as heat, and source coefficients
are not refreshed across a damaged bond history.

### Source-bound frame strain-rate dissipation

`FrameViscousSection::try_new` binds explicit axial and bending viscosity
properties to the SAME original layered elastic section.
`try_new_with_shear` additionally retains an original
`ConcentricTransverseViscosity`: its shear stiffness, authored correction and
independent shear viscosity use the same Timoshenko interpolation as native
stiffness and consistent inertia. All participating source bands intersect;
an explicit zero viscosity does not remove shear elasticity.
`Frame2d::add_material_viscous_element` retains those sources, stiffness and
reference inertia together. Every native Newmark call includes the complete
physical damping matrix. Heterogeneous member losses retain directional,
rotational and inter-node terms; no common Rayleigh ratio is inferred.

`FrameMaterialDamping` retains the original sections and three positive
work-conjugate resistance ports per member. The native full matrix is assembled
from those exact ports, which are also available to reduced receivers. Each
port maps physical nodal displacement to an extension in metres; rotational
weights therefore carry metres. In a shear-deformable member, symmetric
curvature and constant shear share one displacement coordinate: their exact
integrated loss coefficients add in that port. The other two ports retain
axial and antisymmetric-curvature work. Public Euler–Bernoulli strain traces
still refuse these members. `reference_damping_n_s_m` and
`reference_layer_heat_fractions` expose the same native decomposition without
allocation or source queries: original layers receive axial eta*A, curvature
eta*I, and corrected transverse eta_s*A work separately. Symmetric bending and
shear share one physical port; antisymmetric curvature has no shear loss.
`with_material_section` retains every viscosity selector, the original
named/derived G role and the authored correction across source refreshes.
The finite-cell thermal consumer admits temperature-dependent viscosity with
fixed sourced elasticity. Live thermoelastic consumers can retain a
`FrameReferenceKinematics` value independently of later E/G/eta updates. It
supplies layer spring coefficients, damping, heat fractions, RMS curvature for
quadratic energy integration and distinct peak curvature for fibre bounds.
The original beta/shear fractions, correction and consistent mass stay fixed:
updated moduli act on one Galerkin basis, with source potential, entropy and
coefficient work owned by the complete thermal consumer. A new instantaneous
static interpolation is not substituted into an unchanged mass matrix.
`material_damping_energy` and
`rayleigh_damping_energy` are separate parts of the existing total damping work,
not additional heat sources.

This path requires fixed linear geometry and elasticity throughout the frame,
positive inertia on every free coordinate, explicit viscosity support from DC
through the actual sampled Nyquist band, and the complete original spectrum
below Nyquist even when a receiver retains fewer modes. Temperature and geometry
remain frozen. It adds no hereditary loss, thermal deposition or acoustic
radiation. `tests/frame2d_viscosity.rs` supplies continuous oscillator and
heterogeneous two-coordinate references, physical-port work equivalence,
continuation, and source/admission regressions. Its shear checks independently
reconstruct equilibrated section stress rates, integrate viscous virtual work,
and compare native shear-damped decay and harmonic phase/power with the
physical one-coordinate oscillator, including explicit zero shear viscosity.

## Complete frame material inventory for live thermoelastic consumers

`FrameMaterialDamping::covers_all_members` records whether every original frame
member has its material-bound elastic and viscous section, including members
fixed at both ends. A complete material thermal consumer must account for those
members' stored energy and physical layer inventories even when they do not
appear in a free modal coordinate. This flag adds no mass, stiffness or damping
law. The original operator still rejects massless free-coordinate condensation
and retains its complete reference geometry, source sections and native ports.

## Forced resonance on the native linear frame

`Frame2d::prepare_harmonic_response` connects the original reference mass,
elastic stiffness, full section viscosity and explicitly authored Rayleigh
matrix to `fs_modal::harmonic::PreparedHarmonicResponse`. The prepared factor
is reusable for force cases at one frequency. Its `free_dofs()` map preserves
node-major translations/rotations and their conjugate forces/moments; no mode
or off-diagonal inertia/loss term is discarded. The source frame is borrowed
immutably, including its accepted displacement and retained Rayleigh reference.

Peak phasors use `exp(-i omega t)`. Responses contain physical amplitude,
phase, velocity, numerical backward residual and independent cycle-averaged
input/viscous powers. Source viscosity and Rayleigh damping appear once in the
total loss; this is not automatic thermal deposition. Original source bands
admit the actual queried frequency, including explicit DC when supported.
Singular lossless poles refuse without a damping floor or artificial mass.

The scope is fixed linear geometry, elastic sections, stationary homogeneous
supports and positive inertia on every free coordinate. Consistent member mass
or source-derived section spin supplies physical rotational inertia. Original
rigid modes remain intact at positive forcing frequency. Source viscosity
assembly is shared with the native transient owner, without repeating its
natural-mode eigensolve or inventing a Nyquist clock for a harmonic query. The
public transient operator retains its complete-spectrum/timestep admission.
Massless viscous coordinates, moving supports, support reactions,
nonlinear periodic orbits and startup transients are not supplied by this API.
Roundoff-sized native stiffness asymmetry is admitted in diagonal-scaled
coordinates and averaged in scratch; larger skew refuses. Source coefficients
and accepted frame state are not rewritten.

Controls bound the full coordinate count and conservative numerical payload
before dense assembly. Cancellation polls surround native bounded preparations,
per-member stiffness assembly and the shared factorization/solve; an interrupted
native mass/material call finishes before its next checkpoint. G1/G3 tests in
`tests/frame2d_harmonic.rs` and `tests/frame2d_viscosity.rs` cover analytic axial
resonance, full consistent bending inertia, heterogeneous nonproportional loss,
rotational covariance, reciprocity, source-band boundaries, physical power,
singular limits, budgets, cancellation and unchanged accepted state.

### Prescribed harmonic support motion and reactions

`Frame2d::prepare_supported_harmonic_response` retains the original full
reference M/C/K before support elimination. Its support map is the frame's
original fixed DOFs; those translations and rotations can now have independently
phased displacement phasors. `TetElasticAssembly::prepare_supported_harmonic_response`
provides the same physical solve for the native P1/P2 tetrahedral and hexahedral
assembly carrier, with an explicit complete viscosity matrix or explicit zero
viscosity. Its prescribed supports are sorted global labels still present in
the assembly: assemble without eliminating the desired support coordinates
first. Previously eliminated cross terms and reactions cannot be recovered.

Both adapters use `fs_modal::harmonic::PreparedSupportedHarmonicResponse`:
`A_ff q_f = F_f - A_fc q_c`, `A = K - omega^2 M - i omega C`.
The original consistent cross inertia, material/Rayleigh cross damping and
relative support phases remain. Returned displacement/velocity spans every
represented coordinate. Reactions `R_c = (Aq-F)_c` are forces or moments BY the
supports ON the structure, including subtraction of applied support loads.
Independent applied and support work plus full symmetric-viscosity loss close
the cycle-power account; skew gyroscopic terms are nondissipative. Empty and
all-prescribed maps are supported. The original source is immutable.

G1/G3 checks in `tests/frame2d_harmonic.rs` and
`tests/solid_harmonic_supports.rs` cover analytic base excitation with consistent
cross mass, source-coordinate mapping, damping and support force/power, native
tetrahedral motion, DC support translation and total rigid-body inertial
reaction. Shared bounds and cancellation remain active. This is steady linear
harmonic response: it does not evolve startup transients, infer actuator
dynamics, deposit thermal energy, or restore coordinates already removed from
an assembly. Arbitrary boundary impedances use the separate stationary-support
API below; their composition with moving supports is not implied.

### Physical boundary loads on the native frame

`Frame2d::prepare_harmonic_response_with_port_impedance` takes the full
native node-major force map, including supported rows, and the original
frequency-dependent complex load impedance. The same map and its transpose
pair nodal force/moment with physical port flow. Reduction retains every free
coordinate and every supplied port; no modal basis or bare-frame inverse is
introduced. The complete pencil is
`K - omega^2 M - i omega C - i omega B Z B^T`, so attached inertia,
compliance and mutual losses can shift and damp the actual coupled resonance.
A physical load may resolve an otherwise singular bare-frame pole.

The shared harmonic owner admits the load's pointwise Hermitian part and
returns original port flows, opposing efforts and net absorbed boundary work.
That work remains separate from structural viscous power. The load owner
supplies frequency-band validity, physical units, causality and source accuracy;
pointwise admission is not a broadband certificate. Finite impedance at DC
does not encode a spring's `1/omega` pole. Supported ports have zero flow,
but mutual loading may give them nonzero effort; full support reactions are
outside this response. Separate load payload limits reserve the native-map
reduction before delegating the remaining budget to the shared owner.
The same original section-band, geometry, inertia and cancellation conditions
apply. Added `frame2d_harmonic` cases specify analytic attached-mass resonance,
lever-arm/mutual bending loads, physical work, supported ports and refusals.

### Complete damped poles of the native frame

`Frame2d::resonance_poles` returns all `2 * free_dofs.len()` rates of
`(s^2 M + s C + K) q = 0` under `exp(s t)`. A shared native linear-pencil
preparation now supplies both this query and the existing harmonic response,
retaining consistent inertia, section spin, every material-viscosity cross
term and the authored Rayleigh reference. No undamped modal projection or
diagonal loss fit is introduced. The existing scaled quadratic owner supplies
the numerical spectrum; exact diagonal pencils retain its analytic scalar
path and slow overdamped poles.

Fixed linear geometry, elastic sections, stationary supports and positive
inertia on every free coordinate have the same admission as native harmonic
response. All rigid and real relaxation roots remain visible. The query
explicitly selects the original constant Kelvin-Voigt/elastic coefficients
as a continuous law and returns the source's intersected real-frequency band.
It does not filter roots to that band or assert experimental validity at a
complex rate. Pure real roots are relaxation or rigid motion, not oscillatory
resonances; numerical eigenvalues alone do not establish a forced peak,
observability, a mode shape or a certified stability verdict.

Full-coordinate and numerical-payload limits precede dense assembly. An
explicit conservative arithmetic-work envelope and per-window QR sweep bound
limit the full quadratic solve. Existing native assembly and eigensolver
checkpoints propagate cancellation without publishing a partial spectrum or
changing accepted state. `tests/frame2d_poles.rs` adds independent scalar,
coupled-material, full-inertia, rigid-motion and refusal regressions. These
regression definitions require native integration execution; syntax and source
review alone do not claim that gate has passed.

### Native frame resonance with complete physical port-load states

`Frame2d::prepare_coupled_resonance` takes the same original M/C/K pencil and
an existing exact-linear `PortHamiltonian` or `DescriptorPortHamiltonian` load.
The original native force map B supplies both opposing frame effort `B*y_load`
and source input `B^T*v`. All mutual load terms, relaxation states, zero-output
internal modes and algebraic effort multipliers remain in one continuous
descriptor. Neither an isolated load impedance nor its descriptor mass is
inverted. Supported native rows have zero velocity; support reactions are
outside the returned force-port response.

State order is `[frame q, frame p=M*v, original load state]`; frame q/p follow
the original free-DOF map. The load differential prefix and trailing multipliers
remain unchanged. Returned effort entries for p are physical frame velocities.
The full native mass inverse uses equilibrated existing Cholesky solves and an
independent original-row check after symmetric rounding. Every reported mode or
AC response additionally checks `q_dot=v`, `p=M*v` and the original full frame
force equation with the two original port maps. No mass or damping diagonalization
is introduced. Shared pHS storage/loss observations include each physical owner
once; internal exchange is not an additional source, heat or radiation debit.

Homogeneous eigenmode residuals use each original coefficient row's one-norm
times the full original state/physical-velocity/load-effort infinity norm.
This normwise gate allows roundoff in a physically zero frame block of a dark
load mode, without an absolute floor. The separate kinematic and momentum rows
retain their own coefficients. Forced AC keeps componentwise force-relative
checks. These are numerical backward errors, not invariant physical condition
numbers or state-error bounds.

The prepared model reuses the existing `fs-phs` harmonic factor, seeded pole
refinement and contour spectrum owners. `resonance_pole` finds a local isolated
mode; `resonance_spectrum` searches only its declared Laplace-plane ellipse.
Numerical count/rank/quadrature agreement is not a rigorous completeness proof.
Repeated/defective modes can refuse. Stiffness can be semidefinite; no fictitious
spring removes an energy-null rigid coordinate. Infinite descriptor roots are
not reported as finite resonances, and a finite mode need not be observable.

The native source-viscosity band gates each actual AC frequency. For complex
rates, the API explicitly selects the same constant elastic/Kelvin-Voigt law as
`resonance_poles` and retains its original band as provenance; it neither filters
real relaxation roots nor claims complex-modulus validity from a real-frequency
band. The physical load owner must supply its own source-domain and continuation
admission. A periodic particular response does not erase an existing transient.

Full native-coordinate, port, combined-state, numerical-byte and arithmetic-work
bounds precede assembly. The prepared model and one AC solve share that envelope;
pole/contour work has the existing separate query options. Checkpoints bracket
bounded native mass/material and shared PSD-admission calls and poll inverse
columns, assembly members, coupling rows and original physical observations.
No source state or partial result is published on failure or cancellation.
Six G0/G1/G3/G4 regressions in `tests/frame2d_poles/coupled.rs` specify analytic
memory-shifted and hidden poles, descriptor inertance/drag, rigid position,
full consistent mutual inertia and port covariance, original material loss,
physical work and refusals. Native integration execution remains outstanding.

## Original material strain domains on frozen modal solid trajectories

`ViscoElements::validate_displacement` evaluates the existing element Mandel
strain maps at every original material quadrature point and invokes the source
law's strain admission. `MaterialModalValve::validate_modal_displacement` uses
the existing modal-to-node expansion before that same check. These observations
advance no material history and construct no alternate elastic operator. They
poll cancellation at source points and retain all accepted state on refusal.
TET10 observes its admitted quadrature; unsampled spatial extrema and maxima
between time samples remain separate accuracy obligations.

## Original exterior geometry on material valve images

Frozen and hereditary `MaterialModalValve` images retain one immutable
`MaterialValveSurface` prepared from the original Tet4 nodes and exterior
`TetLoadGeometry3` faces. Canonical node triples, physical face areas and
outward normals use the original solid owner, including reversed tetrahedron
winding. Preparation respects its existing cancellation and load budgets.
No geometry is rebuilt during stepping and no material response is fitted.
Tet10 images deliberately carry no P1 exterior surface: a consumer requiring
that trace must refuse rather than invent midside geometry or interpolate an
unqualified contact surface. Mechanical pressure and swept-volume maps retain
their original ownership independently of a consumer's contact-face selection.

## Continuing material-fluid rod loads

`PreparedRodMotion::body_loads` exposes the external conservative wrench at an admitted dynamic pose using the same gravity, hydrostatic displacement, attached-strand forces, and reference potential as the native force assembly. Nodal forces are expressed in world axes; pure attachment couples are expressed in the current body axes. Rotating attachment arms retain their actual mass and displaced-volume distributions, including cases where their centroids differ. `rebind_submerged_fluid` replaces an existing uniform fluid snapshot through `FluidTransport` while retaining the original solid reference inventory, inertial data, displaced-volume policy, reference geometry, and stationary supports. This model-level operation neither advances time nor accounts for environmental work; the continuing time owner must account for the potential change at its actual pose. These operations retain the existing full-immersion, nonoverlapping-volume, fixed-reference-volume approximation and infer no drag, added mass, partial immersion, free surface, solid compression, or reservoir dynamics.

## Native finite-rotation axial-extension covectors

`PreparedRodMotion::axial_extension_gradients_into` writes each original
segment's derivative of `L0*gamma_x` into caller-owned bounded storage. Each
endpoint row contains the world translational covector followed by its body
right-trivialized rotational covector. It reuses the exact geodesic-midpoint
kinematics and analytic virtual-work pullback of the native rod elastic owner,
without cloning a rod or introducing an approximate velocity difference.
The accessor validates the same pose/support dimensions, polls cancellation
per segment, and refuses the same principal-log cut as force assembly. Scratch
contents after refusal are unspecified; no accepted material or motion state
changes. Common rigid translation/spin lies in its nullspace, and finite
world rotations transform the force rows while preserving body torque rows.
Consumers must transform these material-node covectors through their actual
inertia and centroid coordinates before applying a physical loss law.

## Continuing source-temperature rod sections

`PreparedRodMotion::rebind_material_temperatures` refreshes section response through the original `RodMaterialSource`, retaining reference inventory, intrinsic geometry, supports, attached links, gravity and hydrostatic displacement. It admits the exact reference material description and verifies replay of current section claims, expansion paths and strain bounds before resolving a new temperature field. Original unheated reference sections may acquire their explicit thermal binding. The source owner re-evaluates E/rho/G and integrates expansion from the original stress-free reference; current density does not redefine closed-body mass or rotary inertia. `strain_energy_j` evaluates the original conservative strain-storage owner at an admitted dynamic pose. Section geometry, reference displaced volume and separately bound tension-link material states remain fixed.

`LayerIsotropicParameters` explicitly selects each layer's density key and E/nu,
E/G or K/G source pair. Material frames, spatial rods, tension links and prepared
rod motion consume the shared section integrals and check the original source
parents, preserving their state and frequency restrictions without requiring
canonical density/E/G names. `RodMaterialSource::try_new_with_isotropic_parameters`
accepts these bindings in segment/layer order. Temperature refresh retains their
roles, complete pinned requirements, geometry and reference mass/rotary inertia,
and integrates the same expansion source from the original stress-free reference.
`thermoelastic_laws` preserves those bindings through the existing potential:
scalar reference density and expansion, with affine direct E/G, affine E and
scalar nu, or scalar K/G. K/G with non-scalar parents, or E/nu with non-scalar
nu, refuses this thermal potential; prescribed-temperature section refresh
still evaluates supported source curves through the shared projector.
`tests/isotropic_sections.rs` and `tests/rod_material_source.rs` cover source-role
continuity, native mechanical response, reference inventory and these thermal
admission boundaries using synthetic data.

## Native rod strains and source thermoelastic lowering

`RodMaterialSource::thermoelastic_laws` lowers each original concentric segment/layer through the same pinned mechanical and thermal-expansion data into the complete material thermoelastic law. `PreparedRodMotion::segment_strain_states` returns the actual-pose axial, bending, shear and torsional strains and section energy through the existing native kinematic and energy owners. Neither path introduces a replacement rod model, effective material, or extra trajectory.

The reciprocal finite-enthalpy consumer in `fs-couple::rod_motion::enthalpy` uses these views to evaluate thermal exchange and refresh the original source-derived stiffness and integrated expansion while retaining reference mass/inertia, pose, momentum, supports, attached links, loads and clock. Existing source and small-strain refusals remain authoritative. Independent native-energy, nonuniform conduction and torsional-reaction regressions were added under `fs-couple/tests/rod_motion/enthalpy.rs`; native execution remains pending the unchanged dependency/environment build blocker.

## Original bulk and interface sources on thermal joints

`cohesive::thermal_joint::source::SourceJoint` binds the existing restrained-member
and cohesive solver to original `ResolvedThermoMaxwell` bulk sources and one ordered
`ResolvedBilinearCohesive` interface. Receiving member identities must match both
interface sides. Source temperature, strain and rate domains remain enforced by
the native material law; member memory, cohesive damage, force/extension
compatibility and work/storage checks remain owned by the existing joint.

The ordinary advance explicitly admits the unchanged interface query, including quantity
conventions, so evolving bulk temperature cannot silently freeze an unsupported
cohesive-state change. Three source integration tests cover independent series and
thermal-work references, data-only substitutions with relaxation and fracture
continuation, and source/domain/cancellation refusal. This remains a quasistatic
fixed-geometry series model with prescribed bulk temperatures and fixed external
interface state; no thermal transport, fracture-heat partition or
temperature-dependent cohesive-history transfer is inferred from that path.

`SourceJoint::advance_with_interface_temperature[_controlled]` explicitly
accepts a temperature successor from the same card, envelope and pinned claims.
Its normal series geometry consumes common-mode or mixed-mode pure-opening
data; the material extension does not supply shear equilibrium in this owner.
It stages the material owner's fixed-jump history transfer before the original
native equilibrium solve, retaining every bulk Maxwell strain, reference
geometry and source assignment. `TemperatureJointReport` separates signed
interface parameter work, fixed-jump damage release and the native mechanical
accounts. An independent old-to-new energy balance includes both stages before
joint state, source receipts, damage and physical time publish together.

Three focused source/native regressions add an analytical series-force and
parameter-work reference, hot strength/toughness loss followed by cooling
without healing, and source-domain/late-cancellation refusal with identical
retry. This is first-order prescribed-temperature composition; subdivide
combined thermal/load paths for accuracy. It does not solve thermal transport,
assign fracture heat, validate an unobserved continuous temperature path, or
refresh frozen-onset mixed-mode history.

The 34 cohesive constitutive tests, including the four temperature-transfer
tests, passed in an exact-physics slice. Source/native integration targets remain
unexecuted in this worker: compiling the unchanged pinned `asupersync` dependency
was SIGKILLed even with single-threaded metadata-only compilation.

## Source-bound relaxing adhesive bondlines in thermal series joints

`cohesive::thermal_joint::relaxation::SourceRelaxingJoint` consumes the original
`SourceMember` geometries and `ResolvedThermoMaxwell` bulk laws together with an
ordered `ResolvedViscoelasticBond`. It retains the adhesive's independent normal
and shear Maxwell spectra and original thickness; no Kelvin dashpot or bulk
modulus substitutes for adhesive memory. The axial consumer advances the normal
branch bank, while the original zero shear histories remain intact. Actual
receiving member identities must match both ordered interface sides. Combined
member and branch-state budgets are checked before accepted state is built.

Prescribed endpoint extension and endpoint force use the same scalar series
compatibility. Each original member supplies its native affine endpoint stress,
and the adhesive supplies its exact prepared straight-jump-ramp response. Their
affine equation is solved directly, followed by independent original material
updates, common-force and total-extension checks. Bulk Maxwell histories and
the selected backward-Euler or exact-strain-ramp integration survive each step.
All memories, displacement and time publish together; source, temperature,
strain/rate/time, energy and late-cancellation refusals leave the old state
retryable. The complete adhesive query remains fixed, including temperature;
bulk heating does not silently change adhesive coefficients or source claims.

The report separates exact adhesive work, stored change and physical relaxation
heat from original bulk heat and backward-Euler numerical dissipation. Signed
bond and bulk endpoint-work quadrature corrections close the common actuator
account, including exact-ramp bulk laws; neither correction is heat. Thermal
eigenstrain work and modulus-change work retain their original native meaning.
The combined balance uses constitutive increments rather than subtracting large
whole-joint energy inventories. This is a fixed-geometry, massless, small-strain
axial model; endpoint control requires time refinement. It adds no damage,
healing, thermal transport, heat destination, inertial or 3D peel model.

`tests/relaxing_adhesive_joint.rs` checks the independent standard-linear-solid
series relaxation time and constant-force creep retardation time, second-order
time refinement and integrated adhesive heat, unequal heated members with both
native bulk history integrators, signed work corrections, ordered source and
domain refusal, and identical retry after late cancellation.

## Original thermal eigenstrain in tensor-Prony projection

The existing tetrahedral tensor-Prony projector carries each explicitly supplied
thermal free strain into its original material quadrature instead of discarding
it or inventing an external force. Original geometry, physical strain columns,
consistent inertia, material frame and supports remain unchanged. Eigenstrain
admission uses the existing small-strain material domain. Its constitutive
energy and continuing Maxwell memory belong to the original projected material
and native coupled PHS owners; projection adds no heat or independent solver.
Source-integrated isotropic/oriented-orthotropic expansion and reciprocal heat
are exercised by the tensor coupling regressions, currently authored but unrun
in this environment.

## Shared J2 plasticity in fibers, beam-columns and frames

`FiberLaw::J2` consumes the original material owner's exact uniaxial stress
reduction, including plastic lateral contraction and separate isotropic/Prager
history. The existing section tangent, force-based beam-column and planar frame
solvers therefore support metal yielding, elastic unloading and permanent
curvature/deflection without a second structural plasticity law. Raw section
kernels retain their caller-admission convention; `Section::try_respond` and the
beam/frame nonlinear paths check finite response and the full compatible tensor
strain against the original J2 small-strain domain.

`FiberMaterial::resolve_j2` uses the same source requirements, explicit hardening
basis and coefficient conversion as tetrahedral solids, retaining density and
all original receipts at one fixed external state. RC reinforcement may select
the existing steel law or native J2; its replaced concrete area and reference
mass remain coherent. A source post-yield tangent does not invent a kinematic
hardening split. Temperature feedback, rate effects and calibrated cyclic steel
behavior are not inferred for this fixed-state axial reduction.

A force-based cantilever commits every section only after all sections converge;
a refused step preserves the whole accepted element for retry or reversal.
Focused regressions compare independent bilinear section and beam formulas,
source/mass consistency, full-tensor domain refusal and transactional history.
These additions do not provide finite-strain forming, plastic collapse
continuation, fatigue/damage or a complete nonlinear plastic shell solver.

Complete unloading in the native frame iteration retains its independent
compatibility and relative force checks while admitting a local constitutive
roundoff floor derived from section tangent and strain. This avoids rejecting
a physically equilibrated permanent bend solely because endpoint load is zero.

## Saturating metal fibers

`FiberLaw::Voce` and its retained axial history use the exact native
`UniaxialVoce` restriction in existing sections, force-based beam-columns, and
frame elements. `FiberMaterial::resolve_voce` admits independent saturation
magnitude/rate, linear isotropic/Prager coefficients and reference density in
one source bundle. Reinforced sections preserve the same concrete replacement
and constituent mass calculation. Checked section response admits the full
tensor strain including lateral contraction before a structural commit.

`fiber_plastic` checks the independent nonlinear moment/curvature relation,
consistent section tangent, source retention, unloading curvature, and matching
beam/frame permanent deflection. These fixed-condition, small-strain material
checks do not establish cyclic calibration, thermal recoupling, damage, or
finite-strain forming. Existing geometric choices keep their own model limits.

## Incremental plane continua

`plastic2::SmallStrainPlaneProblem` uses the existing `Mesh2` P1 triangles and
Q1 quadrilaterals with native uniform, per-element or per-quadrature laws.
Plane stress solves all transverse stress constraints with the original
material reduction and its condensed algorithmic tangent. Plane strain
retains the constrained transverse stress. Every quadrature point carries its
complete three-dimensional strain and committed constitutive history.

`SmallStrainLawField::PerQuadrature` borrows each original constitutive law in
the exact native integration order. P1 plane triangles require three laws,
Q1 planes four, HEX8 bricks eight, and P1 tetrahedra one. Spatial source
temperature, expansion, stiffness and material directions can therefore vary
inside an element without replacing them by element averages. The same point
law initializes history and supplies stress, tangent and accepted history;
B-bar projects kinematics while retaining these separate constitutive owners.
Pointwise symmetry declarations select the existing general or symmetric sparse
solve. Native fatigue gauges bind the source at their actual quadrature point.
Coefficients remain fixed within an interval; source refresh and retained
continuation use the existing `solve_from` semantics. This adds no thermal
transport, quadrature-error estimate or material interpolation model.

Focused `thermoelastic_plane` regressions exercise actual differential expansion
and restrained stresses under a within-element thermal gradient, with an
independent integrated energy solution for both planes and bricks. They also
check pointwise stress, tangent and exact P1/Q1 coverage. `plastic_hex` checks
separate point histories through reversal and unchanged uniform-field behavior.

Absolute chronological load factors support loading, unloading, reversal and
repeated physical intervals for timed laws. All Newton and line-search trials
reuse the incoming history. Only a converged increment returns updated history
and its matching endpoint stresses. `solve_from` leaves its borrowed accepted
state untouched on refusal. Reactions and free-row force convergence retain
their original mechanical meaning. Geometry, quadrature, matrix-entry and
load-path budgets are explicit, as are full tensor strain and small-rotation
bounds. Prepared geometry and B columns, reused C B products, and one accumulated
element matrix avoid repeated geometry work and duplicate quadrature insertion.
The existing tetrahedral sparse linear correction is shared without changing
its numerical acceptance criteria.

G0/G1/G3 `plastic2` tests cover P1/Q1 cyclic analytical references, confined
plane strain, skew-element global tangents, heterogeneous yielding, timed
relaxation, atomic continuation/refusal, geometry and rotation admission, and
Voce saturation. Standard displacement elements can lock in confined nearly
incompressible plastic flow; the opt-in Q1 plane-strain B-bar path below reduces
volumetric locking. Shell bending, finite-strain forming and limit-load
continuation remain outside this owner.

## Nonlinear backstress recovery in metal fibers and members

`FiberLaw::Chaboche` and its matching state use the exact free-sided axial
restriction of the native tensor law. Signed plastic strain, accumulated
equivalent plastic strain and independent recovering backstress remain separate
through section curvature reversals and member unloading. Section forces and
their tangent come from the same checked candidate; rejected domain exits leave
accepted fiber histories unchanged. Existing force-based beams and frame fiber
sections consume this law through their original section interface.

`FiberMaterial::resolve_chaboche` resolves density, Voce/Prager coefficients and
one independently sourced C/gamma branch from the same fixed-condition bundle.
All original property receipts remain available and density retains its original
mass meaning. The adapter does not infer nonlinear recovery from a monotonic
total-strain tangent or silently freeze evolving material coordinates.

The native fiber battery compares beam deflection and permanent set against an
independently stress-parametrized implicit-recovery curve, checks a reversed
section tangent and retained history, and exercises complete sourced properties
and mass. These are numerical constitutive/member checks, not calibration of
real cyclic behavior, finite-strain forming or a fatigue-life prediction.

## B-bar plastic flow and general constitutive Newton tangents

`SmallStrainPlaneProblem::formulation` admits `Standard` for existing P1/Q1
behavior and `BBar` for Q1 plane strain. The latter replaces only the full
three-dimensional volumetric strain by its element-quadrature mean, adding
`(mean(div u) - div u)/3` to each normal component. Local deviatoric strain,
plastic flow and separate quadrature histories survive unchanged. The same
projected B drives constitutive strain, internal force and algorithmic tangent;
the zz stress contributes to projected virtual work. Returned B-bar strains are
assumed material strains, whose zz component can be nonzero even though the
compatible displacement is plane strain. Both raw and projected strains, and
the raw rotation, must satisfy their independently declared admission limits.
B-bar plane stress and P1 are explicit refusals. History cannot be transferred
between formulations without a separately admitted projection.

Both native plane and tetrahedral increments preserve nonsymmetric constitutive
tangents. A law declaring `tangent_symmetric=false` selects existing FGMRES with
ILU(0) or absolute-diagonal preconditioning, bounded restart/iteration work,
in-Arnoldi cancellation checks and a recomputed true residual. The original
symmetric PCG/PMINRES path remains available. Nonlinear backstress recovery is
therefore usable under nonproportional load without false symmetrization; timed
tetrahedral dynamics inherits the same correction owner.

The focused plane battery includes constant plastic patches, skew-element
deviator/tangent checks, raw-domain refusal and a nearly incompressible pressure
cylinder against the independent plane-strain Lame displacement solution.
`tests/plastic_chaboche.rs` exercises actually nonsymmetric free tangents,
finite-difference force derivatives, triaxial/plane-stress/B-bar equilibrium,
retained nonproportional backstress and exhausted linear-budget refusal.
These tests do not establish a finite-strain, shear-locking or plastic-collapse
mesh-convergence claim.

## True-residual correction for confined continuum solves

The shared small-strain symmetric correction checks the actual Euclidean
residual independently of MINRES's preconditioner-norm recurrence estimate.
If the estimate stops early but the true residual misses the existing gate,
the solver may perform one fresh MINRES correction on `b - A x`. Its iteration
allowance is only the unused part of the original MINRES-path cap; no extra
allowance is created. Acceptance recomputes `b - A (x + correction)` relative
to the original right-hand side, retains only a better finite candidate, and
preserves the original linear gate and physical Newton force tolerance.
Cancellation or an exhausted allowance returns no accepted material history.

The nearly incompressible B-bar cylinder regression now also resolves a
12-by-48 mesh: a single bounded correction removes the observed mismatch
between the MINRES estimate and the independently checked residual. Its
displacement is compared with the same analytical Lame reference and the two
coarser meshes. A separate finite-load refusal test exercises an insufficient
linear allowance. This is a bounded residual correction, not a replacement
preconditioner or a claim that every ill-conditioned mesh can converge.



## Eight-node small-strain plasticity with full 3D B-bar

`plastic3::hex::SmallStrainHexProblem<L>` uses the original hexahedral geometry/Jacobian owner, FEEC two-point Gauss rule in each direction, and the shared native continuum increment/Newton solver. Uniform and per-element laws retain eight independent histories per cell. Absolute load paths support force/displacement control, unloading, reversal and continuation. The tetrahedral and inertial paths retain the same force scaling, corrections and old-history endpoint stress semantics after extraction into `plastic3::increment`.

Standard integration uses compatible strain. B-bar replaces only volumetric strain by the reference-volume-weighted element mean, with the full three-dimensional one-third correction. Identical projected B columns enter strain, virtual work and tangent; no material entries are symmetrized. Existing FGMRES handles genuinely nonsymmetric recovery tangents. Deviatoric quadrature remains fully integrated; there is no one-point stiffness or artificial hourglass penalty. Local 24-by-24 matrices accumulate on the stack before sparse insertion; no per-point matrix allocation is introduced.

Explicit node/element/point/free-DOF/matrix-entry/load-step budgets and sampled Jacobian checks precede assembly. Raw compatible strain and rotation are checked independently of projected strain. Histories commit only at converged increments; a refused candidate leaves caller-owned state unchanged. A rate law contributes one physical interval per increment, including endpoint stresses.

Focused verification covers cyclic axial response and permanent set, complete distorted-element nonproportional Chaboche tangent and force equilibrium, distinct point histories, heterogeneous yielding, actual-time overstress, refused-step retry, and geometry/domain/budget admission. The nearly incompressible axially constrained Lamé-cylinder displacement error falls from 98.3257% (standard 2x8x1) to 1.86313% (B-bar 2x8x1) and 0.500276% (B-bar 4x16x1) under the original physical force gate. This is small-strain volumetric-locking mitigation, not a general stable mixed pressure element, global injectivity certificate, finite-strain model, general shear-locking fix or plastic-limit convergence claim.



### Consistent inertia for generic inelastic Hex8 motion

`SmallStrainHexProblem::advance_backward_euler` and
`advance_backward_euler_with_corner_density` continue the same Standard or
B-bar solid with accepted displacement, velocity and eight native constitutive
histories per cell. Both use `plastic3::increment` and the existing
`BackwardEulerInertia`: the momentum residual, physical force scale and
support reactions include consistent mass, while the original material tangent
receives the same mass divided by the squared physical interval. Every Newton
proposal borrows the same old histories, and endpoint stresses are evaluated
before their single converged commit.

A separate three-point Gauss rule per coordinate integrates the original
trilinear shapes, warped reference geometry and either constant cell density
or a trilinear eight-corner profile. B-bar projection does not change reference
inventory. Positive corner profiles remain local to their original cells.
Combined material/mass admission requires at most 768 tangent entries per cell
before shared-node assembly; the 27 mass points add no material histories.
Prepared constitutive geometry is reused for the motion solve.

The returned `BackwardEulerHexStep` includes the original solid solution,
endpoint velocity and `InertialEnergy`. Inertial work equals kinetic change
plus the nonnegative backward-Euler kinetic defect, including moving-support
work and free/support mass coupling. The defect is numerical loss, never
plastic heat. This generic material trait exposes neither a duration nor an
energy account: a caller selecting a timed law must bind it to the same
physical interval and retain its original reference inventory. Material energy
and reciprocal caloric closure stay with their existing owners. This remains
a first-order, small-strain motion method without a finite-deformation,
contact, mass-transfer or temporal-accuracy certification claim.

Authored regressions in `tests/plastic_hex_dynamics.rs` exercise independent
force-driven shear, unloading and permanent set for rate-independent J2 and
timed overstress laws; exact warped-brick graded-inventory translation and
fully prescribed support work; and refused-step retry, density, endpoint and
combined matrix-budget admission. Their presence is not an executed-test claim.

## Exact-branch Chaboche fibers and structural source receipts

`FiberLaw::chaboche_branches` and `FiberMaterial::chaboche_branches` admit one through eight native axial Chaboche branches. One branch retains the existing variant; a bounded typed dispatcher selects the actual two-through-eight branch law and matching history. Response/commit allocate no heap storage and do not execute padded inactive branches. Enum storage is bounded by its largest variant; this is not a branch-proportional memory claim. Native return, tangent, tensor-compatible lateral flow and source coefficients are unchanged.

`FiberMaterial::resolve_chaboche_branches` resolves the ordered independent C_i/gamma_i roles, Voce/Prager data and density together. It retains the complete original bundle and an explicit immutable branch-index-to-source-role mapping, including N=1. Density names must differ from every constitutive role; a conflicting binding returns a typed refusal before lookup. The material-card owner independently rejects same-name claims with incompatible dimensions. Existing Norton, single-branch, Voce, J2 and RC variants remain available.

Checked section responses refuse material-domain and branch/history-shape mismatches without changing accepted history. Raw multi-branch responses return NaNs on these failures, including newly introduced single/multi Chaboche mismatches. Existing force-based beam and static-frame owners use the same retained fiber histories. RC reinforcement retains the original geometry-consistent density and line mass.

Seven focused regressions cover every supported branch count through cycling, independent three-scale bending/tangent/backstress formulas and reversal, failed-response isolation, complete ordered source metadata and RC mass, custom-density alias refusal, sourced beam/frame loaded and permanent deflection with reversal, and an independent one-interval Norton relaxation oracle. Existing fiber plasticity regressions remain required. These are fixed-condition small-strain cyclic laws, not thermal Chaboche history transfer, new frame kinematics, damage or fatigue calibration.


## Enclosed gas-cavity geometry on deforming solids

`hyper3::pressure::cavity::ClosedTetCavity` selects one connected closed inner
boundary of an already admitted tetrahedral solid. Enclosed gas volume is
separate from total material volume. Gas-outward orientation supplies the volume
and its derivatives while canonical solid-face pressure signs are preserved.
Evaluated geometry stores compact boundary-node `dV/du` and analytic
O(boundary faces) Hessian-vector actions. Positive differential gas pressure
produces `+p grad(V)`; held-pressure endpoint work is exactly `p ΔV`.
The generated canonical face field feeds the existing static and dynamic
pressure owners without a second solid solver or equation of state.
The cavity also exposes `adiabatic_gas_load`, binding an original
`ClosedGasParcel` at the undeformed cavity volume to the exact current cavity
geometry. Every evaluation follows the source-owned isentrope to the trial
volume and pulls the same pressure and acoustic bulk modulus through
`GasMechanicalLoad`. Its conservative contribution is
`Delta U_gas + p_ext Delta V`. This bridge uses no fixed gamma, does not freeze
pressure or mutate the caller parcel, and adds neither heat nor time integration.
The force Jacobian is a `TetCavityGasTangent`: sparse current-surface geometry
plus one retained acoustic rank-one term, `J = J_geometry - w w^T`, where
`w = sqrt(rho*c^2/V) grad(V)`. It retains all nonlocal gas coupling with
O(surface nonzeros + solid coordinates) storage and application work. Entry
access and `apply` include both terms; `apply_into` uses caller-owned output
and scratch without heap allocation, and leaves output unchanged on refusal.
The sparse `geometric_part` alone is not the complete gas stiffness.

The gas load now assembles its geometric term directly from the admitted
cavity's volume Hessian on selected inner faces. Assembly and matrix-free
Hessian actions share the original oriented area-gradient derivative. This
removes an all-boundary pressure vector and a second complete-mesh admission
from every gas load evaluation; sparse dimensions still cover the full solid
chart and unrelated rows are structurally zero. Zero differential pressure
removes the geometric term while retaining the original gas compressibility.
The focused oracle compares positive, zero and negative differential-pressure
tangents against the original pressure owner and the volume Hessian action.

The static acoustic operator retains independent factors as
`K_acoustic = sum_c w_c w_c^T`; merging factors would introduce unphysical
cross-chamber terms. Empty and single-factor actions keep the original
arithmetic. Multiple factors store ordered nonzero free-coordinate entries,
so Krylov application uses work and storage proportional to the active boundary
entries and allocates no temporary vectors. The Jacobi fallback includes the
full signed geometric-plus-acoustic diagonal before taking its absolute value.
Multiple-source force normalization adds each cavity's translation-invariant
bound on each free row before taking the norm. Focused operator regressions
check independent and opposing factors, action and transpose, diagonal,
source/free-chart refusal and the empty/single reductions.

`ClosedTetCavity::volume_change_m3` and held-pressure work use compensated
polarized boundary-cone determinants before subtracting any rounded volume
totals. Both original solid configurations and positive gas volumes remain
admitted. The gas load retains its already allocated compact boundary positions;
static trials integrate `m Delta u + p_ext Delta V` between these cached
endpoints. The original species calorics cancel common formation references,
with complete source model, inventory and auxiliary-coordinate identity checked.
This avoids both subtraction of reference-relative potentials and another
entropy solve or complete mesh admission. Thermodynamic endpoint accuracy still
depends on the declared source entropy controls.

`AdiabaticCavity` is a borrowed binding for the original source and geometry.
`solve_with_adiabatic_cavity_loads` composes it with stationary gravity and
optional native boundary pressure in one Newton solve. All loads share the
current displacement, stiffness, continuation factor and support reactions.
The additional boundary pressure is explicit: applying ambient pressure to the
complete material boundary together with cavity pressure minus that same ambient
pressure leaves gas pressure acting once inside and ambient pressure once outside.
The gas remains on its original reversible isentrope with fixed mass and
composition; this static composition adds no gas transport or physical time step.
Focused native regressions exercise original nonlinear calorics, direct volume
and gas work, source refusal and a gas-filled hollow solid under ambient pressure
and spatial gravity against independent equilibrium and gas-law balances.

`solve_with_adiabatic_cavities` extends that same native solve to an ordered
slice of independently sourced sealed chambers. Each retains its original mass,
composition, calorics, entropy controls and exterior pressure. Current forces,
geometric stiffness and independent acoustic factors contribute once; the
solid carries mechanical interaction between chambers. Multiple-cavity force
and potential reductions use compensated arithmetic. Every binding must match
the exact solid reference chart, and overlapping canonical boundary faces
refuse before gas evaluation, including duplicates made by independent geometry
preparation. A refusal leaves every borrowed inventory unchanged.

`HyperTetCavitiesSolution` returns the solid solution and the actual accepted
gas loads in input order. Their volumes, thermodynamic states, potentials and
tangents need no final entropy re-solve; only the final load-step states are
retained for publication. Empty and single-cavity inputs preserve their original
loading semantics. Native tests cover two mechanically connected hollow solids
with distinct variable heat capacities, an independently specified finite
compression with gravity and ambient pressure, and the second chamber's volume
and temperature response when the first is heated. Source ordering, inventory,
force balance, empty/single reductions and boundary ownership are checked.

For multiple chambers, a private immutable batch validates common reference
charts and boundary ownership once. Every Newton/Armijo trial then performs
the original complete-solid geometry admission once, using a canonical anchor
independent of source order. Each chamber still forms its current compact
positions from its own cached reference offsets and original displacement
offsets. It never subtracts already rounded global positions. All positive gas
volumes and original thermodynamic sources admit independently; a failure
publishes no partial batch. Empty and single-source paths retain the original
evaluation. This removes repeated complete-mesh preparation and allocation;
full-chart force and acoustic source vectors remain per cavity. No measured
wall-clock speedup is asserted. A distant-chamber regression checks retained
tiny volume work, independent equal-reference surfaces, source-order invariance
and bitwise agreement with standalone force, potential and tangent evaluations.

The same geometry and static equilibrium path accept native
`ClosedRealGasParcel` inventories through the sealed `AdiabaticCavityGas` seam.
`AdiabaticCavity`, `TetCavityGasLoad` and `HyperTetCavitiesSolution` retain the
concrete parcel type, defaulting to the existing ideal source. Peng–Robinson
endpoints retain the original temperature/density-dependent calorics, source
cards and homogeneous-branch admission. Gas storage increments include both
NASA and PR departures, while force tangents use the same-state `rho*c^2/V`;
`gamma*p/V` is not substituted for dense-fluid compressibility. Both sources
use the original `GasMechanicalLoad`, selected-face geometric tangent and
independent acoustic factors. There is no dynamic dispatch in boundary or
Krylov loops and no second equation of state or equilibrium implementation.

Four focused real-gas regressions use an independent scalar PR/Helmholtz
reference for finite-volume energy, entropy, pressure and acoustic bulk modulus;
check the full force tangent and potential work; solve a manufactured 4%
compression of the native hollow solid under ambient pressure; and heat one of
two mechanically connected dense-gas chambers. They retain distinct calorics,
large formation datums and original inventory/source identity, and exercise
source-domain refusal and cancellation. The manufactured gas is synthetic;
these checks do not supply empirical calibration, phase equilibrium, heat
transfer during deformation or time integration.


Preparation rejects duplicate, out-of-range, open, disconnected and outer-shell
selections. Each evaluation admits every material tetrahedron and requires
positive current cavity volume. Relative reference and displacement coordinates
preserve a large common translation. These local checks do not certify
self-intersection, global injectivity or contact. Geometry supports the admitted
boundary size. The cavity supplies that admitted size to the shared material
pullback and uses its cancellable force path, avoiding the default small-chart
256-coordinate limit. Spatial admission and Jacobian actions poll cancellation
in bounded tiles. Gas evolution, heat exchange, mixed incompressibility and
time integration remain with their respective owners.

Additional G1/G4 regressions refine the hollow specimen to 6,144 tetrahedra
and 390 cavity coordinates, compare matrix-free forces/stiffness/work against
source-resolved finite differences and analytic dilation, retain off-graph
acoustic coupling, and check output atomicity under refusal and cancellation.

Five focused regression groups cover independent volume/cofactor derivatives,
work and pressure signs, matrix-free tangent action and symmetry, objectivity,
actual native static-solver compression of a hollow twelve-tetrahedron specimen,
and invalid geometry/selection/cancellation. Independent analytical and NumPy
calculations verify the fixture and its tangent; the authored native Rust tests
have not been executed in this disk-constrained session.


### Simultaneous creep and yielding in native frames

`FiberLaw::NortonJ2` and `FiberState::NortonJ2` retain the original free-sided
joint axial law, separate signed creep/plastic strains and both accumulated
paths. `Frame2d::solve_creep_step` shares its existing section/member/global
Newton solve, reference quadrature, physical-duration checks, small-geometry
admission and work gates across Norton and joint Norton/J2 fibers. Repeated
section and equilibrium trials borrow the same accepted history; acceptance
advances both memories once.

`FrameCreepEnergy` reports separate physical creep and plastic dissipation.
Energetic hardening remains storage and numerical loss never becomes heat.
`creep_plastic_fiber` reads accepted histories without a constitutive update.
`with_creep_laws` accepts either axial family and preserves the original family,
elastic/hardening storage, kinetic coefficients, complete card/claim/query
identities and fixed external coordinates. Yield may refresh with temperature
only through its original pinned source. Authored and sourced histories cannot
be silently exchanged.

Untimed static, displacement, modal and inertial routes refuse both timed
families. Corotational geometry, finite deformation, thermal strain and reciprocal
heating remain outside this frame account. The `frame_creep_plastic` target
contains independent axial load/hold/unload references, nonuniform-temperature
cantilever curvature and tip-displacement references, and immutable-history
checks for unsupported routes, source replacement, inconsistent clocks and
cancellation. These manufactured fixtures are implementation checks, not
measured service-life or structural-safety qualification.

Verification for this increment: the complete native material/cut-cell/solid
source chain and both frame integration targets passed pinned-compiler metadata
checks. The joint material behavior passed its 22 focused executable regressions;
the subsequent expanded material suite passed 56/56. Native frame tests were
typechecked, not numerically executed. This is not the full Cargo/DSR gate.


### Thermal expansion and reciprocal heat in native creep intervals

The same generic `SmallStrainTetProblem::solve_creep_from` accepts
`ThermalCreep` with its complete temperature-bearing history. It retains the
original sparse equilibrium solve, common physical duration, total displacement
gradient bounds, frozen trial histories and endpoint constitutive consistency
gates. `CreepTetEnergy` additionally integrates the material's change of T times
mechanical entropy and signed reversible caloric heat over the original reference
tetrahedral volumes. The complete applied/support work balance includes these
terms exactly once. Plain Norton and Norton/J2 laws retain exact zero values for
both new accounts. Thermal transfer uses the per-element total caloric heat,
including reciprocal exchange; numerical loss and hardening storage remain separate.

The complete native source chain and the original/joint frame and tetrahedral
integration targets typechecked successfully. Those four native integration targets
were not numerically executed in this validation; the material behavior suite was.


## Native tetrahedra with thermally transported Chaboche history

The existing generic SmallStrainTetProblem directly consumes ThermalChaboche with its complete ChabocheState and nonsymmetric tangent. Each resolved old-to-new temperature law represents one physical solve_from increment; all Newton trials retain old state. For a rate-independent source, a new-to-new law recovers accepted endpoint tractions without repeating thermal transport. Timed sources retain the accepted stress; reevaluation represents another physical interval. tests/plastic_thermal_chaboche.rs verifies a six-tetrahedron held member through heating, cooling and compressive reversal, retained separate branch histories, free-row equilibrium and native per-element thermal energy. A nonproportional traction case independently checks the full nonsymmetric assembled tangent. Prescribed temperature changes reaction and loaded lateral displacement through actual constitutive flow. These tests establish a bounded small-strain prescribed-temperature capability, not an automatically heat-coupled Chaboche solid, finite-deformation or experimental qualification.


## Timed cyclic plasticity in the native tetrahedral consumer

The existing generic tetrahedral owner also advances source-pinned timed thermal Chaboche histories through heating, tensile/compressive holds and reversal. The timed regression in tests/plastic_thermal_chaboche.rs uses affine Voce Q(T), independent affine branch moduli and n=3 Perzyna kinetics. Independent endpoint flow, branch-recovery and viscous heat, coefficient storage and the first law are checked alongside force equilibrium and support work. Repeated frozen-history solves are identical; a shorter physical interval relaxes less. The test uses prescribed temperature and does not claim automatic Chaboche-state coupling in ThermoplasticSolid's heat owner, finite deformation or experimental validation.

## Retained finite-strain viscoelastic tetrahedra

`hyper3::visco::FiniteViscoTetSolid` consumes uniform or per-element original
`ViscoHyperelastic` laws. `try_new` admits the native reference geometry and
relaxed history. `advance(dt_s, supports, loads, body_load, settings, controls,
cx)` solves one actual quasistatic endpoint through the existing generic
`HyperTetProblem`, with its original force-equilibrium gate, consistent tangent,
linear solver and globalization. Every trial borrows the same accepted
`ViscoState` and physical interval; no trial or load-search iteration advances
viscous memory. Original reference geometry is prepared once and shared between
accepted candidates. Regional spectra are never averaged.

Displacement, complete element histories, Piola stress, physical storage and
elapsed time change together only after the native solve and physical work
gates pass. Dropping or refusing a candidate leaves the previous solid intact.
Absolute endpoint supports and loads permit continued loading, unloading,
stress relaxation and creep; `dt_s = 0` is the native elastic jump. One native
load endpoint is required. The caller controls physical interval subdivision
and explicit finite work tolerances.

The accepted energy record distinguishes physical storage before/after,
nonnegative native relaxation loss, native frozen-state loading work, reduced
incremental potential, actual internal/dead-load/support endpoint work, native
material balance and signed endpoint integration defect. The potential is used
only by Newton. Native work minus storage change and relaxation loss is checked
independently; assembled potential and support/load equilibrium work are also
checked. The signed endpoint defect is actual internal work minus physical
storage change and native loss. It has its own explicit allowance and is never
deposited as heat or hidden by clipping dissipation.

`tests/hyper_visco.rs` covers an affine relaxation patch with a free interior
node against native material-point histories, maintained-load creep with
physical work/loss accounts, recovery of the original elastic solver when there
are no Maxwell branches, and heterogeneous compatible histories with failed-step
rollback and deterministic retry. `cargo run -p fs-solid --example
viscoelastic_creep` reports extension, maintained force, storage, interval loss
and endpoint defect for an explicitly synthetic constrained-lateral coupon.

The base `FiniteViscoTetSolid::advance` path is isothermal quasistatic finite-bulk
P1 tetrahedra. Volumetric locking is not mitigated; this is not the Taylor–Hood
mixed/incompressible formulation from IRP M1.3-W2b. The motion owner below adds
inertia. Contact, evolving temperature, source-parameter resolution, general
time/mesh convergence and experimental qualification require their respective
owners and evidence. No closure of that full mixed-element bead is claimed.

### Native motion with retained viscoelastic memory

`hyper3::visco::dynamics::FiniteViscoTetMotion::try_new` binds the complete
accepted solid to explicit original reference density and initial nodal velocity.
Consistent `P1ReferenceMass` is prepared once and shared between candidates;
current volume and viscous deformation do not change mass. Previously relaxed
or deformed material histories and elapsed time are retained.

`try_new_with_nodal_density` uses the same native reference-mass owner for a
positive P1 nodal density. It retains the complete original field and exposes
cell means only as inventory diagnostics. `reference_acceleration_load_n`
projects acceleration through this exact matrix, including supported rows.
Candidates share the mass and preserve all accepted Maxwell histories; neither
viscous deformation nor a new temperature query changes the closed inventory.
G1/G3 regressions in `tests/hyper_visco_dynamics/nodal_density.rs` cover graded
gravity/free flight, a relaxed-memory oscillator with analytic recurrence,
kinetic-loss rollback, uniform bit parity and invalid-input refusal.

`advance` uses the same positive physical interval for the native branch update
and `BackwardEulerInertia`. The original hyper3 Newton owner assembles material
and inertial forces/tangents together, uses its physical free-row momentum
stopping rule, and includes consistent-mass coupling in support reactions.
Free bodies require no artificial supports. There is no alternate motion solver
or lumped-mass approximation. The quasistatic endpoint remains unchanged when
inertia is absent; its zero-time elastic jump remains available.

Tiny-strain energy subtraction can obscure a final Newton decrement. An
otherwise rejected inertial trial may be accepted only when its supports are
fully applied and its actual free-force residual already satisfies the original
convergence gate using the trial's native material/inertia/load normalization.
The next Newton iteration still certifies that same gate. No residual tolerance,
material-domain condition or physical work allowance is relaxed.

The returned step distinguishes physical stored-energy change, native viscous
loss, signed material endpoint defect, kinetic-energy change, nonnegative
backward-Euler numerical kinetic loss, and momentum/kinetic balance residuals.
Their sum closes actual applied/support endpoint work. Numerical kinetic loss
has its own finite caller allowance and is never material heat. A refused solve,
work gate, numerical-loss budget or cancellation returns no next motion and
leaves displacement, velocity, time, original mass and all viscous histories
unchanged. Dead loads are held over each interval; no follower/contact load,
mass transfer or automatic thermal coupling is introduced.

`advance_with_momentum` runs that SAME interval and adds an original-mass
`SolidMotionMomentum` account about a supplied world origin. It transfers the
already assembled dead-load vector from the native solve; it does not repeat
geometry assembly, freeze/update material history again, or solve another
momentum equation. Two sparse actions evaluate `M delta_v` and `M v_before`.
Actual represented endpoint positions give
`Delta L = sum((x_after-origin) cross M delta_v + delta_x cross M v_before)`.
Compensated products/reductions preserve small signed contributions. External
and support torque impulses use the same endpoint positions as the original
force equation. The second angular term is separately reported as backward-
Euler numerical drift, including coordinate rounding; it is not physical
torque. Consumers must bound both the balance residual and that numerical
defect, including accumulated absolute defects over a trajectory. Physical
relaxation, signed endpoint-work defect and numerical kinetic loss retain
their existing distinct accounts. The ordinary `advance` result is unchanged.
An independent graded-simplex regression checks moving supports, body/nodal
loads, actual angular inventory and nonzero numerical drift at two origins.

`hyper3::visco::dynamics::split` provides native force/contact/drift primitives
for loaded sustained-contact consumers. `force_kick_candidate` uses the accepted
original first Piola stresses and explicit referential loads, then applies their
combined impulse through the original supported consistent mass. Configuration,
material history and time remain fixed during the kick and any following contact
impulses. Immutable coordinate, Piola and Maxwell-history buffers are shared by
velocity-only candidates; new buffers accompany actual material/configuration
evolution. `contact_drift_candidate` advances coordinates with the actual retained
post-contact velocity and updates every Maxwell branch through the same native
constitutive endpoint evaluation used by the implicit interval owner.

The kick privately binds the original immutable accepted material state, mass,
loads, duration and stationary support set. Velocity-only native impulses retain
its new force-phase identity; a skipped or superseded force phase, different
configuration, clock, material history or mass refuses. Repeated nodal and body
loads accumulate with compensation before each resultant is rounded and scaled
by the interval, retaining a small imbalance between large opposed loads.
Editing public diagnostic fields cannot change the continuation or its physical
account. The private identity check is constant-time and requires no duplicate
constitutive-history snapshot.

Native relaxation loss, source displacement work, force splitting and the signed
material endpoint defect remain distinct. The signed integration account closes
storage plus kinetic change and physical loss; the sum of absolute phase defects
has an independent caller allowance. Kinetic drift retains the same velocity and
mass, so its kinetic change is exactly zero. One sparse mass action computes both
momentum inventories and the stable actual angular increment `sum(delta_x cross
M v)` with compensated products. Coordinate-rounding residue is bounded as a
numerical momentum error, never an external torque or heat source.

These primitives are a first-order explicit force/kinetic split with native
constitutive relaxation. The enclosing consumer must refine its timestep and
admit contact geometry; the primitives do not supply collision detection or a
general contact-existence guarantee. Failure or cancellation preserves the input
owner. The ordinary implicit `advance` retains its original integration rule.
Focused native tests specify graded supported recoil, load cancellation,
intervening contact impulses, constitutive replay, actual momentum inventory and
phase-binding refusal. Their execution remains unverified while the required
constellation sibling checkouts are unavailable.

Native plastic, viscoelastic and thermal Tet4 owners expose their immutable
original reference coordinates/connectivity and assembly bounds for geometry
consumers. These views borrow accepted owners and never rebuild material state.

`tests/hyper_visco_dynamics.rs` checks unconstrained rigid acceleration against
the independent discrete momentum solution; a deforming oscillator against
the small-strain scalar Maxwell/backward-Euler recurrence with the original
consistent mass; three timestep rungs against the continuous elastic oscillator;
and post-solve kinetic-budget rollback with deterministic retry. The comparison
explicitly bounds the small finite-strain correction. These are local model/time
checks, not a general finite-strain trajectory or spatial-convergence certificate.
`cargo run -p fs-solid --example viscoelastic_creep -- --dynamic` releases a
synthetic initially stretched coupon and reports displacement, velocity,
physical storage/kinetic energy, viscous loss and both numerical work terms.

## Finite-strain Taylor–Hood displacement/pressure solids

`hyper3::mixed::MixedUpProblem` solves P2 Tet10 displacements with one continuous
P1 pressure per original corner vertex. Uniform, per-element and per-quadrature
mixed material fields retain their original laws; the material owner supplies
all four consistent saddle-point tangent blocks. Exact incompressibility is
weakly imposed in that pressure space. Original finite-bulk mode preserves the
material's admitted volume/pressure dual, including its domain restrictions.

The original Tet10 topology/map and explicit Duffy rules prepare shared immutable
reference geometry. Four through twelve points per collapsed direction are
admitted: the nonlinear pressure constraint needs the quartic rule, not the
linear-stiffness minimum. Curved geometry and nonlinear energies still require
quadrature refinement. Tangents are accumulated as 34×34 element blocks and
scattered once per element; repeated Newton actions reuse the assembled matrix,
its original sparse LDL factor, and one symbolic analysis for the fixed pattern.
The previous numeric factor is released before its replacement is allocated.

The native Newton–Krylov owner supplies residual Armijo globalization. Required
length and stress scales define positive work-conjugate congruence scaling;
force [N] and constraint-volume [m³] residuals are not directly mixed. The same
physical force and complete pressure residuals are checked at acceptance.
Native LDL uses 2×2 pivots and requires inertia (n_free_u positive, n_free_p negative),
including for an initially equilibrated or zero-iteration endpoint. An unstable
or singular candidate refuses. Numeric factor payload, sparse assembly payload,
mesh/free-DOF counts and quadrature-point count have explicit caller caps.
Cancellation is checked during assembly and around the bounded native factor.

Pressure normally follows traction equilibrium. An optional exact-mode pressure
pin is accepted only when the selected component's constant-pressure field is a
true null mode on free displacement rows; the omitted constraint is also checked
at the accepted endpoint. Finite-bulk pins and physically overconstraining pins
refuse. The public raw evaluation retains constrained rows, reactions, physical
residuals, the full tangent and mixed saddle potential. That potential is not
silently interpreted as viscous stored energy.

`tests/mixed_up.rs` verifies independent uniaxial traction/stretch/pressure,
original bulk/shear ratios 1 through 1e8, full tangent derivatives and objectivity,
quartic constraint integration, pressure gauges, budget/cancellation boundaries,
and rejection of an already equilibrated unstable compressed bar. A non-affine
Cook panel has a three-rung refinement diagnostic and compares bulk/shear 1e6
against exact incompressibility. This bounded evidence is not the complete
Treloar/MMS/inf-sup campaign or a universal locking-free certificate; bead
`frankensim-irp-v5-instrument-realism-pexfp.61` remains open. The displacement-only
P1 hyper3 path retains its original volumetric-locking limitation.

## Retained finite-strain mixed viscoelastic solids

`hyper3::mixed::visco::MixedViscoTetSolid` owns the accepted displacement,
continuous P1 pressure, time, Piola stresses and one distinct
`ViscoRemainderState` per TET10 quadrature point. `MixedViscoLawField`
accepts uniform, per-element or per-quadrature original spectra without
averaging. Prepared P2/P1 geometry is shared across accepted endpoints.

`advance` solves one physical endpoint with accepted metrics and interval frozen
throughout the native mixed Newton solve, warm-starting displacement and pressure.
A zero interval means an instantaneous elastic jump. Numerical load continuation
must have exactly one endpoint; physical substeps remain explicit caller actions.
A returned candidate owns the complete new state. Cancellation, Newton refusal
or a failed work gate cannot mutate the accepted solid.

Physical storage is the native remainder plus `U(J(p))`, reconstructed from
pressure. It is not the original bulk energy at the quadrature deformation:
the mixed pressure constraint holds weakly. Exact incompressibility never
evaluates the unused bulk penalty, even when publishing accepted branch updates.
The report separates remainder/bulk storage, native viscous loss, frozen-history
loading work, incremental remainder and mixed saddle potentials, pressure
constraint work, actual internal/dead-load/support endpoint work, equilibrium
work residual and signed endpoint-work defect. Explicit
`FiniteViscoStepControls` gate the work residual and defect; neither numerical
defect nor saddle potential is counted as physical dissipation or heat.

`tests/mixed_visco.rs` verifies held-stretch relaxation against the native
material point, maintained-traction creep with an independent scalar equilibrium,
original-bulk storage and the non-aging elastic limit, heterogeneous spectra and
transactional work-gate refusal. Exact-mode relaxation is bitwise independent
of unused bulk moduli `1e3` and `1e300`. The runnable
`mixed_viscoelastic_creep` example applies exact quadratic-face loads and reports
stretch, pressure, storage, interval viscous loss and signed endpoint defect.
The standalone solid advance is isothermal and quasistatic; the retained motion
owner below adds consistent inertia. Neither path supplies thermal closure, a
general mesh inf-sup certificate or the full M1.3-W2b acceptance ladder.

## Consistent TET10/P1 mixed viscoelastic motion

`hyper3::mixed::visco::dynamics::MixedViscoTetMotion` binds an accepted mixed
solid to explicit original per-element density and initial nodal velocity.
Its immutable consistent reference mass is integrated from the same cached
quadratic shapes and positive reference weights. Signed corner/edge entries
are retained; row-sum lumping and deformation-dependent mass changes are not
used. The matrix is exact for affine TET10 geometry and uses the caller's
explicit quadrature approximation for curved geometry. Mass assembly has a
separate byte cap and cancellation checkpoints.

The existing native backward-Euler inertia enters the same mixed Newton solve:
`M*a` augments displacement forces and `M/dt^2` augments only the displacement
tangent. The pressure equations are unchanged and pressure carries no mass.
Residuals, final admission, support reactions and the KKT inertia check use the
complete incremental system. The reported mixed potential remains the material
saddle potential. Static entry points continue to omit inertia.

The original reference mass owner also supports signed consistent entries in
its tangent and uses `abs(M)` for its positive momentum reference bound; existing
positive P1 and Hex8 actions retain their accumulation order. Density-weighted
mass, not current volume, defines kinetic energy and momentum.

A returned motion candidate advances displacement, pressure, velocity, material
history and time together. Physical branch loss, signed material endpoint-work
defect and nonnegative backward-Euler kinetic loss remain separate. Existing
`FiniteViscoMotionControls` independently bound material work and numerical
kinetic loss. Refusal leaves the complete accepted endpoint unchanged.

Exact pressure enforces endpoint volume weakly. Velocity is the backward-Euler
displacement chord; there is no separate tangent-velocity projection. This is
a first-order convergent motion path with explicitly budgeted numerical damping,
not the instrument-render doctrine's zero-algorithmic-damping or discrete-gradient
integrator. Closed work accounts do not establish timestep accuracy or general
mesh inf-sup stability. Thermal evolution, follower loads and contact retain
their separate owners.

`tests/mixed_visco_dynamics.rs` checks exact-incompressible free-body
acceleration, a genuine quadratic shear mode against its independent scalar
backward-Euler recurrence and continuous solution on 20/40/80 intervals,
released finite-bulk versus exact-pressure viscoelastic motion, signed moving-
support reactions and complete rollback after numerical-loss refusal. Native
inertia unit tests check rational P2 coefficients, linear-velocity kinetic
moments, tangent finite differences, the absolute-mass bound, original positive
P1/Hex8 arithmetic, budget and cancellation. The existing creep example also
offers `--dynamic` to release the deformed body and report physical and
numerical losses independently.

### Adaptive physical intervals with retained finite-strain history

`MixedViscoTetMotion::advance_adaptive_pair_toward` advances one accepted pair
toward an explicit target time. It compares one original full interval with two
original half intervals from the same accepted displacement, pressure, velocity
and per-quadrature Maxwell history. Only the two-half trajectory is returned:
there is no extrapolation of displacement, pressure, metrics or material work.
The original consistent P2 mass and immutable reference geometry remain shared.

A caller callback supplies explicit support, nodal-force, body-force and pressure-
gauge snapshots at each actual candidate midpoint and endpoint. The captured end
snapshot is reused by coarse and fine solves. No load/support interpolation is
inferred. The caller must resolve source bandwidth and stop at discontinuities;
local sampling cannot detect skipped source cycles.

All original nodes, including supports, participate in displacement and velocity
checks. Every original quadrature point contributes objective Green-strain and
full Cauchy-stress checks, and every original branch contributes its symmetric
inverse viscous-metric error with the true Frobenius shear weights. Cauchy stress
is reevaluated at zero duration from the accepted native history and actual P1
pressure; positive-duration incremental stress is not substituted for physical
storage stress. Explicit SI absolute tolerances combine with change-relative
motion/strain/history scales and rotation-invariant endpoint stress magnitudes.
The largest local ratio controls first-order step-doubling refinement.

One call accepts exactly two physical history intervals and performs at most
three native endpoint solves per attempt, with an explicit cap of 32 attempts.
The raw full/fine disagreement is a local consistency indicator, not a global
error bound or a certificate for prescribed-support derivatives. Native source,
material, equilibrium and work failures propagate; this is accuracy adaptation,
not general failure recovery. Both fine reports retain actual material loss,
signed endpoint defect and numerical kinetic loss separately. Pair-wide sums of
absolute endpoint/balance defects and numerical kinetic loss must fit the same
supplied work allowances. No rejected coarse history or work enters acceptance,
and cancellation or any refusal leaves the original motion unchanged.

The focused adaptive regressions contain independent quadratic-shear recurrence
and time-refinement checks, scalar nonlinear uniaxial Maxwell log-strain updates,
physical pressure-inclusive stress and frame-invariance checks, actual callback
times under retries, pair-wide kinetic-budget refusal, and immutable rollback.
Authored native tests have not been executed in this disk-constrained session;
stdin rustfmt and independent scalar/equation checks are separate evidence.

### Instantaneous impulses with retained mixed viscoelastic memory

`MixedViscoTetMotion::impulse_candidate` transfers nodal impulses through the
original signed consistent P2 mass for `VolumetricTreatment::OriginalMaterial`.
The same shared mass owner supplies bounded CG, repeated-load accumulation,
actual free-row/linear/angular momentum gates, kinetic change and signed
support work. Held DOFs retain their current velocity, including nonzero
support motion. Angular momentum uses the current deformed geometry.

The complete mixed solid survives unchanged: continuous pressure, displacement,
all original quadrature Maxwell metrics and accumulated relaxation loss,
stresses, remainder and pressure-reconstructed bulk energy, physical clock,
source laws, reference mesh and density. Only velocity changes. A later native
`advance` continues these accepted histories. An instantaneous event creates no
material relaxation loss or heat, and the caller retains ownership of its
contact/restitution law and opposite body's impulse.

This mass-only entry returns the typed
`MixedViscoTetError::ExactIncompressibleImpulse` for an exact pressure
constraint. It does not assign pressure a mass or silently project an
incompatible velocity. All solve, input, budget, balance and cancellation
refusals leave the accepted motion unchanged.

`MixedViscoTetMotion::constrained_impulse_candidate` supplies the separate
exact-incompressible path through the existing shared P2/P1 pressure-impulse
owner. Each original accepted quadrature deformation supplies its endpoint
cofactor; the shared native mass/pressure KKT solve retains the original
supports, weak volume rows, physical scales and sparse budgets. Internal
pressure impulse has independently gated zero resultant, torque and work.
Its Pa s multipliers are reported separately and never overwrite physical
pressure or any material history. The finite-compliance formulation returns
the typed `FiniteCompliancePressureImpulse` refusal from this entry.

The incoming velocity must already meet the caller's weak-volume-rate allowance.
A velocity accepted through backward-Euler position chords is not automatically
tangent-compatible; the event refuses a large pre-existing rate defect instead
of concealing it through a dissipative projection. A pressure-impulse gauge is
admitted only for a true connected constant null mode. Shared pressure and
mass accounting gate the actual rounded successor before the same atomic
complete-history transfer used by the finite-compliance entry.

`MixedViscoTetMotion::constrained_impulse_response` exposes that prepared
exact-volume response at material points. Its directional and complete matrix
queries use the same ten signed quadratic trace weights as its single and
simultaneous point-impulse candidates. All response columns and the final native
event share one KKT factor. Construction checks the accepted weak velocity rate;
finite-compliance motion retains its separate typed mass-only entries. The
prepared borrow prevents a stale factor from following a changed accepted state.

`Tet10ImpulsePoint` names a point by its original cell and four barycentric
coordinates. The native quadratic shape owner supplies all ten signed weights
for current isoparametric position, velocity and transpose impulse distribution;
negative corner weights are retained. `point_kinematics` samples either pressure
formulation. `point_directional_inverse_mass_per_kg` uses the original supported
consistent mass response without a material trial or history clone, while
`point_impulse_candidate` delegates the complete-history nodal event. These two
mass-only operations explicitly refuse exact incompressibility. The shared P1
and P2 trace helpers retain conjugate point work and torque on curved geometry;
they do not replace the original mass by a point mass. The finite-compliance
entry is connected to the existing rigid impact owner in `fs-couple`.

The same finite-compliance owner exposes `point_mobility_matrix_per_kg` and
`point_impulses_candidate` for simultaneous quadratic contacts. The original
signed P2 trace and supported consistent mass retain every cross-point and
cross-direction entry. Caller row/workspace caps are admitted before allocation;
one mesh-sized response is retained at a time, without a dense mass inverse or
material-history copy. Loads remain separate until the original compensated
nodal admission combines them, followed by one complete native candidate.
Exact incompressibility retains its typed refusal on these mass-only entries.

The focused `mixed_visco_dynamics/impulse.rs` regressions use independent
rational P2 integrals and signed corner loads, deformed angular impulse,
moving-support work, nonzero retained pressure and relaxed Maxwell histories,
subsequent native motion, typed formulation refusal and transactional retry.
The exact-mode extension checks an independently integrated pressure impulse
on an already stretched, relaxed body, retained translation and Maxwell
continuation, and incompatible velocity/gauge/factor-budget refusals. A genuine
pressure-gauge regression also checks cancelling applied and moving-support
work with unchanged kinetic energy and retained relaxation history.
The material-point extension checks signed quadratic interpolation and loads,
supported mobility against independent P2 integrals, torque, work and retained
history. The point owner's polynomial fixture also checks curved reference
geometry and virtual work independently of the finite-element solve.
The simultaneous-point regression independently inverts a rational P2 corner/
edge mass block, including negative cross mobility and rotated directions,
and checks summed point work, torque, moving supports and retained histories.
The exact prepared-point fixture additionally checks a relaxed, pressure-loaded
body against independently tabulated centroid/corner impulse distributions,
matrix-predicted velocities and the original nodal constrained event.
Native execution remains pending the existing dependency-build blocker; this
addition is not a completed Cargo or DSR gate.

## Deformable cohesive dynamics with retained inertia

`cohesive::dynamics::CohesiveMotion` adopts an existing `CohesiveState`, its
original reduced consistent mass operator, explicit velocity and clock. It
retains the same elastic operator, sparse jump maps, patch areas, ordered source
bundles and accepted fracture histories. There is no second interface assembly
or independent damage state. Matching and nonmatching P1 traces use their
existing work-conjugate maps. The mass must use exactly the same reduced DOFs,
be symmetric positive definite and remain unchanged; the elastic operator must
be symmetric positive semidefinite and unchanged. These are requirements on
the original operator owner, not conclusions from dimension or energy samples.
Rigid-body modes do not require artificial elastic restraints in this path.
The bulk potential is `u^T K u/2` in the admitted reduced coordinates; an
additional prestress potential or moving essential boundary data requires its
own mechanics, inertia and support-work treatment.

For a physical interval `h`, the existing Newton/FGMRES owner solves
`M (v1-v0)/h + K (u0+u1)/2 + sum B^T A t_mean = f_mean`, with
`v1 = 2 (u1-u0)/h - v0`. Each trial starts from the same accepted material
history. `t_mean` is the original analytic vector integral along the straight
three-component jump path and includes surviving-bond viscosity. The endpoint
traction tangent is not substituted for the derivative of this path mean.
The original material path derivative supplies the action
`(2 M/h^2 + K/2) v + sum A B^T D(t_mean) B v`. Fixed-size local material blocks
are cached at each exact Newton iterate and reused throughout FGMRES; all
original sparse jump columns and native mass/elastic couplings remain present.
No global dense tangent or positive-definite softening approximation is used.
Velocity-regularized damaged-area friction differentiates the original
path-integrated damaged compression together with the full vector slip law.
Its local block retains normal/shear coupling, changing damaged area and the
transverse derivative of the saturated slip direction; zero slip in the linear
regularization branch has its finite shear slope. These contributions share
the same cached material block and leave the physical force and work unchanged.
Ideal Coulomb retains its independent reaction coordinates `lambda = k z`.
The augmented action differentiates the original natural map
`lambda - projection(lambda + k delta_g_t, mu P)` through both disk-interior
sticking and disk-exterior sliding, including the original pressure gradient.
Mechanical rows retain `A B_t^T d_lambda`; constraint rows retain both the
jump and reaction-coordinate blocks. The cache keys the complete augmented
iterate, and a changed static/kinetic coefficient starts a new solve/cache.
Zero-coefficient bindings add no coordinates and use the exact original native
mechanical action. No friction spring, regularization, or extra storage is added.
Explicitly unsupported material switches, the Huber speed corner and Coulomb
disk boundary retain the existing fallible numerical-difference adapter.
Residuals and analytic actions/preparation attempts share the original
operator-evaluation allowance; numerical fallback counts its actual residual
calls. The complete physical endpoint is independently replayed before
acceptance. A source, stencil, cancellation or work-budget refusal preserves
the original accepted motion and material history.

The G1 friction regression reconstructs a three-coordinate inertial reclosure
from independent momentum equations, including changing compression and both
linear-regularized and saturated shear. It reproduces native displacement,
velocity and friction work while remaining bit-identical under widely different
inactive numerical-difference scales. Retained partial damage and physical
viscosity retain their separate original loss accounts. This establishes the
solver path and derivative behavior, not a wall-clock speedup or specimen
calibration.

The ideal-Coulomb regressions also require exact zero-capacity equivalence and
static recapture after a stopping impulse, including finite-difference scales
which would be inadmissible as physical constitutive probes. The reaction stays
within the original static capacity with unchanged accepted damage/history and
zero frictional heat during the held interval. Independent natural-map derivative
checks cover pressure coupling, tangent rotation, and projection boundaries.

The step reports original kinetic, bulk and cohesive endpoint energies,
intrinsic fracture loss, physical bond-viscosity loss, optional damaged-area
friction loss, external force work and
prescribed-jump work separately. Prescribed work is
`sum A t_mean dot (jump1-jump0-B(u1-u0))`, so adoption of a previously loaded
solid neither guesses its old actuator offset nor drops that actuator's work.
The full gate is external plus prescribed work minus kinetic/bulk/cohesive
storage changes, fracture, viscosity and friction. No numerical residual is relabeled as
physical dissipation or routed to a thermal field. Signed roundoff-scale bulk
energies are retained rather than clamped, with materially negative values
refused under the explicit energy allowance.

`DynamicsStep::kinetic_energy_change_j` retains the actual fixed-mass quadratic
increment `mean(v) dot M delta(v)`. The original mass acts on compensated
endpoint differences; a nonzero low word gets its own action. Mean-product
layers are reduced before final rounding, so opposing transfers and a large
unchanged mode cannot erase small accepted storage changes. The usual exact
difference needs one mass action and two vector buffers. The mass owner's own
action accuracy, symmetry and purity remain required. The native interval
independently replays its original-to-final kinetic increment instead of
subtracting total inventories or summing rounded substep changes. Absolute
endpoint diagnostics, original force/constitutive work and physical losses
remain separate. G0 tests cover a coupled mass counter-mode and a 1.5-joule
exchange beside a large transferred inventory; the native interval regression
retains a small oscillator exchange beside unchanged free motion.

Native cohesive storage changes use each original material trial's incremental
account, integrated over its original patch area, rather than subtracting large
absolute endpoint inventories. A passive compressive preload cannot hide the
recoverable energy of small shear motion or loosen an absolute-only work gate.
The focused consistent-mass P1 regression compares that motion, velocity and
bond-energy increment with an independent midpoint oscillator solution while
the rounded absolute cohesive endpoint energies remain equal. The complete
interval gate retains those accepted increments across its substeps; the same
regression compares the interval with repeated native steps and their original
constitutive history.

Thermally driven series joints likewise combine the original member and bond
storage increments. Their member/constitutive gates, force compatibility and
work-conjugacy gates remain in force. A separate native spring/dashpot regression
retains a small bond storage gain and viscous loss beside an unchanged bulk
preload, with independently specified endpoint work and its signed quadrature
correction.

The accepted endpoint replay also reports `DynamicsPatchDissipation` in original
patch order, with actual area-integrated fracture, viscosity and friction losses.
Its sums use the same arithmetic as the global loss accounts. A downstream
thermal owner can identify the original patch/source through the read-only
`CohesiveState::system` accessor; these reports alone assign no heat fraction.

Complete failure removes tensile and shear resistance through the original
material law while retaining mass, velocity, history and clock. Original normal
compression remains available on reclosure. The default is the existing
frictionless penalty law; the explicit sourced damaged-area extension below
adds residual shear resistance. Neither is exact nonpenetration or newly
discovered contact. Accepted failure also releases the material's former
opening/shear separation and dashpot-rate restrictions; its independent normal
penetration bound and the structural displacement/velocity bounds still apply.
Native P1 regressions continue both free bodies' mass-weighted centres beyond
the adhesive range with unchanged momenta, and source-temperature regression
retains the detached displacement, original claims and spent fracture work.
Immutable
candidate publication makes late cancellation, material-domain refusal, failed
work gates and exhausted solves leave the entire accepted trajectory intact.

This is fixed-reference small-strain dynamics with fixed material/source
conditions during each mechanical increment, fixed frames and fixed topology.
The explicit temperature transfer below may advance the cohesive source and
history between increments. It does not update bulk coefficients, impose
moving reduced constraints, discover contact or provide a unique/stable
nonlinear root.
Prescribed interface offsets and reduced external forces are explicit inputs;
their physical construction and resolution of reversals, crack events and
softening branches remain caller responsibilities. Smooth elastic motion has
the midpoint second-order limit, but a closed energy account alone does not
establish temporal accuracy through fracture.

`tests/cohesive_dynamics.rs` specifies original P1 consistent-mass frequency
and refinement, independent scalar actuator/viscosity work, release of an
already softened solid, and source-bound three-component debonding of two
unrestrained deformable tetrahedra. Its postfailure checks retain each body's
momentum, drive physical compressive reclosure and retain the same sources and
damage. Late cancellation, a material-domain refusal and work-budget exhaustion
also retain the complete accepted trajectory. These focused native tests do
not constitute a dynamic crack-propagation, finite-strain or service-validation
campaign.

### Instantaneous cohesive momentum transfer

`CohesiveMotion::impulse_candidate` applies an explicitly supplied generalized
impulse in the original free-coordinate ordering by solving `M delta_v = J`
with the native checked CG owner and supplied SPD preconditioner. It uses
operator actions and bounded Krylov storage, without a dense inverse or a
replacement mass. The original immutable motion supplies displacement, damage,
source-bound viscosity, friction/contact bindings, clock, step count and controls.
Only an accepted velocity and its incremental event account are published.

The event independently rechecks `M (v_after - v_before) - J` using the actual
representable returned velocity. Its signed work is
`((v_before + v_after)/2) dot J`; kinetic change uses the corresponding mass
increment. Compensated accumulation and a work tolerance scaled to the event
avoid hiding a lost impulse in a large unchanged background energy. Numerical
workspace, solver iterations, cancellation and the original velocity limit are
explicit gates. Refusals leave the accepted state intact.

No displacement path or physical time elapses, so the event does not invent
fracture, viscous or frictional loss. The next physical interval retains those
original owners. An affected Coulomb tangential trace reconsiders inherited
sliding occupancy using that interval's actual incoming velocity, including its
prescribed rate; prior reactions and material history survive the impulse.
An unaffected trace retains its existing occupancy. Accepted finite-time replay
clears this pending classification, and temperature transfer preserves it.

The caller owns impulse generation, other-body reaction, restitution and any
impact energy loss, and adds the signed event work to its trajectory account.
These reduced-coordinate balances do not establish Cartesian angular momentum
without the original geometric trace. Moving support impulses, new collision
topology and automatic heat routing are outside this operation.

### Damaged-area friction in the original deformable dynamics

`CohesiveMotion::with_damaged_area_friction` binds one explicit regularization
speed per selected original patch. The kinetic coefficient is read from the
canonical `cohesive_residual_kinetic_friction` property already retained in
that patch's complete `ResolvedBilinearCohesive` source bundle. The receiving
owner requires a nonnegative dimensionless coefficient at fixed external
conditions; it never infers friction from cohesive strength or fracture energy.
Unselected patches keep their original response, and zero mu retains the
original arithmetic. Immutable shared configuration avoids per-step copying of
the patch coefficient array.

The existing material owner supplies the analytic path measure
`P = integral_0^1 d(s) Kc max(-gn(s),0) ds`, including changing damage and
compressive reclosure. With the actual straight-path tangential velocity `v`,
the mean resisting shear is `mu P v/max(|v|,v_reg)`. It is assembled through the
same `B^T A` map as the original bond and solved with the original mass,
elasticity, cohesive history, viscosity and Newton/FGMRES owner. A signed net
normal adhesive force or endpoint `d*p` is not substituted for that integral.
The original material's damage variable remains the sole damage history.

Friction work is the mean added shear dotted with the actual tangential jump,
integrated over original area. It is nonnegative, reported separately from
fracture and surviving-bond viscosity, and debited exactly once in the full
mechanical energy audit. Prescribed-interface work uses the total bond plus
friction traction, so a moving support supplies its actual contribution. No
heat destination or entropy claim is assigned implicitly. Source-temperature
transfer refreshes the same selected coefficient in the original successor
bundle while preserving regularization speed, damage, mass, velocity, time and
step count. Changing mu at fixed jump creates no friction work or heat.

This is velocity-regularized kinetic Coulomb on a declared damaged-area mixture,
not ideal static friction, a tangential elastic spring, roughness calibration or
new contact discovery. The speed is a numerical closure in m/s, independent of
the timestep; timestep refinement still resolves reversals and stopping. The
fixed-reference and original material-domain limitations remain in force.

`tests/cohesive_dynamics_friction.rs` compares actual native inertial response
against independent scalar equations for intact, partially damaged and fully
failed compressed bonds, including surviving-bond viscosity. It checks an
independent analytic integral while shear initiates and grows damage under
compression, rotated prescribed-slip work, zero-mu equivalence, low-speed
physical-time regularization, source-temperature substitution, missing/invalid
data and cancellation before final publication. These manufactured-data tests
establish implementation behavior, not measured adhesive-friction calibration.

### Static sticking and kinetic breakaway of reclosed deformable bonds

`CohesiveMotion::with_damaged_area_coulomb` selects unregularized static/kinetic
friction on the same original damaged area. It requires both canonical residual
friction properties in the original source bundle, dimensionless at fixed
external conditions, with `0 <= mu_kinetic <= mu_static`. Initial closed-contact
occupancy is explicit; accepted endpoint motion and occupancy determine whether
the next interval may test static capture. This owner is mutually exclusive
with the velocity-regularized extension. Neither introduces a second cohesive
history or changes the original damage, pressure or surviving-bond viscosity.

Two scaled tangential traction coordinates per nonzero-capacity selected patch
join the existing midpoint Newton/FGMRES system. The natural map is
`t - project_disk(t + k*delta_g, mu*P) = 0`, using the common `fs-material`
projection already used by rigid and finite-patch damaged-area owners. `P` is
the original analytic integral of damage times compressive pressure over the
actual straight cohesive path, including changing damage and closure. The
positive numerical `k` conditions the map and scaled unknowns; it is not an
interface spring and stores no energy. The native body mass, elasticity,
work-conjugate trace map and full coupling between patches remain in the solve.
An initially static candidate with resolved slip or endpoint velocity is re-solved at its sourced
kinetic capacity. All passes share the original nonlinear iteration and
physical-residual budgets; Krylov admission counts the augmented dimension.

The accepted physical audit checks the static reaction disk or the direction
and magnitude of kinetic traction. Sticking has exactly zero physical friction
loss only when slip and both endpoint velocities are unresolved; finite
tolerance-scale traction times slip remains in the work/energy
defect. Sliding loss is `A * mu_kinetic * P * |delta_g|`, and the same added
traction contributes to prescribed-interface work. Original fracture and bond
viscosity losses remain separate. Temperature transfer refreshes both original
coefficients while retaining damage, source identity, geometry, inertia,
velocity, clock and occupancy. Zero coefficients preserve original arithmetic.

This is a bounded straight-path rigid-friction specialization on fixed damaged
area, without tangential elastic slip storage or newly discovered contact.
Resolved reversal of tangential endpoint velocities while the pressure path is
active refuses the whole increment, including nonzero-net-slip reversals;
time-step/event refinement must resolve stopping before physical heat is
accepted. Nonlinear convergence does not establish time accuracy or uniqueness
of redundant static reactions, and singular or exhausted native solves refuse.

The Coulomb cases in `tests/cohesive_dynamics_friction/coulomb.rs` check static
hold and kinetic breakaway against independent inertial equations for partial
and complete damage, a rotated two-component reaction disk, changing damage
and opening-path pressure, actual prescribed work, source-temperature threshold
changes, physical kinetic heat below the displacement tolerance, zero-coefficient
parity, source/budget/cancellation refusal and both
zero- and nonzero-net-slip stopping refusals. Source data are manufactured;
these tests do not establish specimen calibration.

### Native cohesive intervals across Coulomb stopping and resticking

`CohesiveMotion::advance_interval` completes a requested duration under constant
external force and affine prescribed jumps using the original mechanical owner.
Each substep retains the original maximum step, constitutive laws, histories,
momentum equation and physical acceptance tolerances. Only the typed
`UnresolvedFrictionReversal { patch }` refusal triggers bracketing: shorter
candidates are solved from the same accepted predecessor until an admitted
candidate stops inside the original velocity-resolution band or resolves the
declared time bracket, then the remaining interval is advanced. Other refusals
propagate unchanged. A rejected trial's damage, motion
and losses are never adopted. No velocity projection, stopping impulse, extra
friction spring or numerical heating is introduced.

The original endpoint velocity tolerance determines static capture. Residual
kinetic energy stays in the actual native state and first-law account. This is
bounded continuation of the discrete straight-path model, not an exact stopping
time, earliest-event certificate, monotone-guard proof or unique softening-path
claim. Source conditions stay fixed over the interval; affine prescribed motion
and constant force are explicit restrictions. Native time resolution, conditioning,
source domains and unresolved events can still refuse.

Independent caps bound attempted solves, accepted substeps, cumulative physical
residual evaluations and nonlinear iterations. Each attempt reserves its native
allowance; successful work refunds unused counts while a refused solve charges
the whole reserved allowance. Per-residual material work retains the original
patch bound. Cancellation or exhaustion returns no partially accepted interval.
Memory retains only the current motion and best bracket candidate, not an entire
state trajectory. Original fracture, bond viscosity, friction and prescribed
work are accumulated separately, including per-patch totals. The interval's
endpoint first law is checked against the sum of original accepted tolerances.

The interval also retains each original patch's maximum normal jump over its
incoming endpoint and all accepted straight local paths. Rejected trials never
enter these bounds; a final closed endpoint cannot erase earlier opening. The
final remainder uses the exact requested clock coordinate and the unchanged
positive native timestep, avoiding an extra rounded time addition after a stop.
Ordinary single-step clock arithmetic is unchanged.

For resolved incoming slip, the augmented Coulomb solve starts its reaction
guess on the sourced kinetic disk at the incoming damaged compression and in
the actual displacement predictor's slip direction. This avoids a zero-reaction
line-search barrier for tiny sliding increments without another material trial
or solve. The original integrated path pressure, nonlinear equations, work
limits and physical audits still determine the accepted reaction and losses;
stationary capture retains its previous reaction predictor.

`tests/cohesive_dynamics_friction/interval.rs` compares an unaligned scalar stop
with independent constant-friction stopping distance and kinetic-energy loss,
then holds above kinetic capacity and breaks away above static capacity. It
also checks accumulated original fracture/viscosity/friction on prescribed
paths, native replay, work bounds, late cancellation and immutable rollback.
Path maxima are checked against native replay and incoming opening followed by
closure; nonzero-origin stopping retains the exact requested final clock.
A small native kinetic step is checked against constant-force displacement,
endpoint velocity and frictional-work formulas before completing the stop.

### Temperature-dependent adhesion in deformable equilibrium and dynamics

`CohesiveSystem::with_temperature_sources` prepares temperature successors on
the original patches. `CohesiveState::transfer_temperature` accepts those
coefficients only after matching the original elastic operator, DOFs, patch
identities, areas and complete sparse trace maps. Every ordered source retains
the original card, selected claims, material surfaces and nonthermal query
coordinates through `ResolvedBilinearCohesive::require_temperature_successor`.
Unsourced interfaces and substituted geometry/operators refuse.

Each quadrature point uses the existing material owner's fixed-jump coefficient
change followed by irreversible relaxation. Common-mode and frozen-physical-mode
BK histories retain their original interpretation; cooling cannot heal bonds or
rescale past fracture work. Area-integrated signed parameter work, additional
damage loss and storage change are reported separately and gated by
`parameter work - storage change - damage loss`. No displacement or dashpot
work occurs during transfer, and neither parameter work nor damage loss is
automatically assigned to heat. Bulk elasticity and thermal expansion remain
outside this transfer.

`CohesiveMotion::transfer_temperature` retains the original mass, velocity,
physical time, accepted step count and numerical controls. The next native
mechanical solve responds to the resulting interface tractions. The transferred
state is not asserted to be equilibrated at its new temperature. All sources,
histories and integrated accounts publish atomically; per-site and bounded
vector polls preserve cancellation without partial accepted state.

`tests/cohesive_temperature.rs` checks independent heated-bond storage/damage
identities, area integration across two original sites, cooling irreversibility,
native continuation to failure, changed-frequency inertial motion with retained
mass, trajectory-budget retention, source/geometry/operator refusal and late
cancellation. These manufactured source curves verify the integration; they
are not physical validation of a temperature-dependent adhesive specimen or a
coupled thermal-transport/heat-partition model.

## Shared original consistent-mass kinetic increments

Backward-Euler and centered reference inertia and the original cohesive
dynamics owner use `fs-math::quadratic::increment` with their unchanged mass
actions. Actual endpoint difference/mean tails and compensated reductions retain
small kinetic transfers beside large common motion or cancelling coordinate
transfers. The native reference-mass row accumulation order and force/velocity
trajectories remain unchanged; no lumping or substitute mass is introduced.

Backward Euler independently retains `delta-v^T M delta-v / 2` as numerical
kinetic loss. Centered inertia reports zero numerical loss. Independent native
force work uses actual displacement endpoints and remains double-double through
the kinetic residual, so kinematic rounding or operator error is still visible.
Neither numerical account supplies material heat. Absolute kinetic diagnostics
and existing admission tolerances retain their original semantics. Focused
owner fixtures cover coupled counter-motion and the separate BE/centered work
identities; their full integration execution remains pending the native
dependency/resource blocker. The shared arithmetic leaf has executed directly.

The native tetrahedral `loads3::transient` Newmark owner also evaluates kinetic
and unshifted strain increments through its original sparse operators; moving
supports use the original full element actions, including free/fixed mass
coupling. Source and reaction work remain independent and retained through the
gate. Changing stiffness contributes its original coefficient difference before
contraction with displacement. Interval transfers, rather than unchanged stored
energy, scale the relative work tolerance. Public increments and retained native
kinetic/work/loss getters support subsequent physical heat composition.

Thermal material refresh independently evaluates the original shifted element
storage and retains it with kinetic energy through physical admission. Rayleigh
loss remains the explicitly supplied damping model; no numerical residual becomes
heat. Regressions cover a resolved local acceleration beside a disconnected
coasting inventory and transactional refusal when a prescribed velocity rounds
away but support work is nonzero. These integration definitions await the full
native build; they do not inherit the arithmetic leaf's executed result.

## Conservative finite-strain elastic motion

`hyper3::mixed::conservative::{ConservativeTetDefinition, ConservativeTetSolid}`
retain quadratic TET10 displacement and velocity, physical stored energy, elapsed
time, original element material identities, and a shared consistent reference
mass. `FullMaterial` evaluates the original finite-strain energy through the
objective C-based discrete gradient in `fs-material`. Initial displacement and
velocity are explicit; preparation does not silently relax a supplied state.

`ConservativeTetSolid::impulse_candidate` applies zero-duration nodal impulses
through the same full consistent reference mass for `FullMaterial` and
`OriginalMaterial`. Current deformed positions supply the angular-momentum
arms. Held velocity coordinates retain their actual velocity, and the shared
impulse owner reports their reaction impulses and signed work. Displacement,
physical time, material points, original density, retained bulk volumes and
algebraic pressure remain unchanged. Applied and support impulse work enter
the original compensated cumulative account; the initial energy is not reset.
Subsequent native intervals therefore retain the entire physical work history.

The same public method retains the native pressure-impulse constraint solve
for `ExactIncompressible`. Compressible formulations use only `controls.balance`
and require no pressure gauge; pressure/rate arrays are empty and pressure
momentum/work is zero. `MixedImpulseReport::kkt_inertia` is `Some` only when the
native constrained factor exists, and `None` for the checked mass solve.
The caller owns restitution/contact physics and any opposite body's impulse.
G1/G3 fixtures in `tests/mixed_conservative_bulk/impulse.rs` compare independent
quadratic-simplex momentum and kinetic integrals, signed support impulses,
moving-support work, retained state, native continuation and refusal behavior.

`ConservativeTetSolid::point_kinematics` and `point_impulse_response` connect
all three formulations to the same native signed quadratic point owner. The
prepared response returns full cross-point/direction mobility and one summed
point-impulse candidate. Finite-compliance modes use their original mass solves;
exact mode reuses one pressure factor for all columns and the event. Both nodal
and point events share a single complete-state transfer and the original
compensated cumulative-work gate. Focused fixtures compare the three modes
against exact simplex responses, curved point torque, moving-support work,
retained pressure/bulk storage and subsequent conservative time evolution.
Native execution of these additions remains pending the existing dependency/
resource blocker; formatting and independent algebra are not Cargo validation.

Each interval uses centered kinematics, `v_after = 2 Delta u / dt - v_before`,
and force `2 M (Delta u - dt v_before) / dt^2`. The native consistent P2 mass
retains its signed entries and quadratic body-force weights. The inertia work
equals the kinetic-energy increment with zero algorithmic loss. This is a
separate `MidpointInertia` path; existing backward-Euler motion keeps its
original dissipative semantics.

`ConservativeTetDefinition::prepare_with_reference_density` binds four positive
P1 corner densities per original tetrahedron to that same signed consistent P2
mass. The cached barycentric pressure basis supplies the density samples; no
geometry reconstruction, row lumping or new time integrator is introduced.
The immutable field and assembled mass remain shared through accepted motion
for full-material, exact-incompressible and original finite-bulk formulations.
Four Duffy points per direction admit graded affine inertia; a varying profile
on a curved quadratic map requires six. Original constant cell fields retain
their scalar quadrature semantics and accumulation bits.

`reference_corner_density_kg_m3` exposes the retained field. For a graded body,
`reference_density_kg_m3` reports reference-volume-weighted cell means; physical
inertia and loads use the full field and mass operator. The body's
`reference_acceleration_load_n` projects nodal acceleration through its actual
full `M`, including prescribed coordinates and signed quadratic corner loads.
Uniform gravity and affine acceleration therefore retain the body's original
mass, force resultant and torque. The returned loads are held over a native
interval; evolving position-dependent force laws retain their own owner.

Independent barycentric polynomial integrals check every affine and curved
mass entry. Native deforming, translating and spinning trajectories check
energy and both momenta; acceleration fixtures exercise all three formulations
and constrained support reactions. Admission, scalar parity and cancellation
checks retain the original transactional boundary.

Assembly contracts each point's true generally nonsymmetric endpoint tangent
into one local displacement block before one element scatter. The native
Newton-FGMRES owner uses that actual Jacobian; only its symmetric part supplies
the sparse LDL preconditioner. The preconditioner is neither a substituted
elastic tangent nor an inertia-based stability certificate. Geometry, mass and
accepted material points are retained between candidate evaluations.

Dead nodal/body loads, moving Dirichlet supports and returned signed reactions
have explicit work accounts. Each accepted interval checks constitutive work,
equilibrium work, local total-energy balance and cumulative total-energy balance
against the caller's allowance in joules. Kinetic work and numerical loss remain
separate fields. Invalid domains, exhausted budgets, failed solves, failed work
gates or cancellation return no changed accepted state. One physical interval
requires one load step; Newton backtracking does not advance material history.

The mixed conservative owner carries `DiscreteGradientState::energy_change_pa`
and the original bulk increment through its existing quadrature sweep.
`ConservativeTetEnergy::stored_change_j` and `bulk_stored_change_j` expose these
physical increments even when the absolute before/after totals round equal.
Material, bulk, pressure, dead-load and support work use compensated products
and sums through the final balance. The cumulative gate sums each independently
computed work-minus-energy defect in double-double precision; neither a large
initial energy nor large accumulated through-work masks a small defect. Nodal
and point impulses add their original direct kinetic increment to the same
account, preserving all local work and atomic refusal gates.

Prescribed reference events integrate the material and bulk point-to-point
potential increments while retaining physical mixed volumes. Their source input
is defined by that potential change, so the event's local balance is zero by
construction and does not claim an independent thermal or entropy closure.
Previous cumulative defects remain intact. Both direct changes are exposed in
`ConservativeReferenceEnergy`, and any following native interval separately
audits mechanical and exact-volume actuation work. Natural-volume actuation
moments use the reference owner's direct determinant increment. Focused solve
regressions exercise weak shear under high bulk preload, small impulses under
large strain storage, and cancelling reference/modulus work; original full,
exact-constraint and finite-bulk motion and reference-event fixtures remain
covered by the same numerical owners.

G0/G1 evidence in `tests/hyper_conservative.rs` includes independent consistent
mass and original-energy checks for a freely translating, rotating and deforming
body; exact constant-acceleration motion; moving-support work; and rejected or
cancelled candidate rollback. A genuinely nonlinear quadratic shear benchmark
uses an independent quartic Ogden oscillator, independent scalar discrete steps,
and a time-refinement ladder against a converged continuous reference. The
focused native-source lane passed four tests: the three-rung time order was
1.99 and 2.00, and the 48-step free-body energy/momentum discrepancies were at
floating-point roundoff. The native inertia unit lane passed nine tests.

`FullMaterial` is displacement-only and a very large bulk modulus may still
lock. The exact and original finite-bulk mixed formulations below supply
separate pressure treatments. Nonlinear modal reduction, contact, viscoelastic
dissipation and the frequency-exact finite-rotation nonlinear split remain
separate. Centered nonlinear time integration has phase error even when its
discrete work balance closes; the complete-coordinate exact linear owner below
has its own scope. The bounded nonlinear fixtures do not satisfy the separate
long-horizon, nonlinear frequency-matching and complete port/checkpoint
obligations of bead `frankensim-irp-v5-instrument-realism-pexfp.62`.

## Conservative exact-incompressible P2/P1 motion

`ConservativeFormulation::ExactIncompressible` extends the retained elastic
motion owner with continuous corner-node P1 interval pressure and endpoint
constraints `integral N_p (J_after - 1) dV = 0`. The original bulk-free
remainder supplies physical storage and elastic work. No bulk penalty is
evaluated, subtracted or approximated; the pressure carries neither mass nor
stored energy. Original full-invariant fibre dependence is retained.

Momentum uses the objective discrete volume derivative, while the endpoint
constraint uses the endpoint cofactor. The two Jacobian cross-blocks are therefore
distinct at finite motion. The pressure contribution to the displacement block
includes the full derivative of the discrete cofactor. Native Newton-FGMRES
retains all these nonsymmetric derivatives. Its symmetric saddle preconditioner
uses the endpoint constraint derivative in both cross-blocks, keeping redundant
endpoint constraints singular instead of hiding them by averaging the crosses.

Preparation admits both the initial weak volume and the initial weak volume
rate `integral N_p cof(F_initial):grad_X(v_initial) dV`. Their explicit
controls are `maximum_relative_volume_residual` and
`maximum_initial_relative_volume_rate_per_s`, normalized by total reference
volume. Inconsistent initial velocities are refused without a pressure impulse
or silent projection. This is an initial compatibility check; subsequent
centered velocities are not claimed to satisfy a separate exact endpoint
tangency constraint.

`maximum_point_volume_work_defect` bounds each dimensionless local
discrete-volume identity. The accepted step separately reports actual pressure
work, its discrepancy from pressure dotted with the weak-volume change, and
the endpoint pressure residuals in cubic metres. Pressure work, local and
cumulative physical energy, and weak volume are all admitted before publishing
the new state. `algorithmic_pressure_pa()` returns the last interval's
algebraic work multipliers, initialized at zero; these are not endpoint Cauchy
pressures. `pressure_nodes()` identifies their sorted corner-node layout.

Independent endpoint pressure rows on the free displacement space are required.
The final factorization checks this even for a zero-iteration Newton solution.
There is no pressure-pin/gauge API; closed fully prescribed displacement cases
with redundant volume equations refuse. This method enforces weak P1 volume
constraints, not pointwise J=1; a coarse element can retain local volume
variation even when every weak moment closes. It is not a conservative
finite-bulk mixed formulation.

G0/G1 regressions in `tests/mixed_conservative.rs` use independently integrated
Mooney-Rivlin remainder energy, rational consistent P2 mass, momenta and weak
volume. The free nonaffinely strained, spinning body changes shape substantially
over 32 intervals while preserving energy and both momenta to roundoff and
normalized weak volume within `2e-15` in this fixture. Separate checks cover
bitwise trajectory independence from unused bulk moduli `60` and `1e300`,
incompatible initial volume/velocity, invalid rate controls, and redundant
pressure-row refusal without state mutation. Full-material regression tests
remain green. `examples/conservative_elastic_motion.rs` runs both formulations
and prints deformation, energy, work, sampled J, weak volume and interval
pressure with explicit SI controls. These bounded tests do not establish a
mesh-independent inf-sup constant or the remaining long-horizon/render contract.

## Continuous harmonic loads on complete structural coordinates

`ExactHarmonicSolid::prepare` and `prepare_viscous`
retain the same validated complete consistent
mass, stiffness, optional symmetric PSD physical viscosity, initial displacement
and initial velocity as the held/affine owner. The existing
`fs-phs::linear_exact::PreparedExactTones` supplies continuous prescribed forcing,
complete coupled propagation and independently integrated work/viscous loss.
There is no second exponential solver, modal truncation, fitted damping,
mass lumping or frequency retuning. Exact undamped resonance keeps its growing
transient instead of substituting a singular steady-state phasor response.

`ExactHarmonicSolid::advance` takes cosine/sine coefficients in N in free-coordinate-major,
caller-tone order. Frequency entries are angular frequencies in rad/s and may
include explicit DC, noncommensurate tones and coherent repeated frequencies.
Each pair describes `Fc*cos(omega*tau)+Fs*sin(omega*tau)` from that interval's
local start. `advance_force_phase` uses cached rotations to continue a
source after acceptance; it does not mutate the physical state or source phase.
Changing a next-interval pair is an explicit finite force jump, not an impulse.
Checkpoints must keep the caller's source coefficients with the accepted state.

Work and dissipation retain every cross-frequency term and initial transient;
they are not sums of cycle-averaged isolated-tone powers. The original `M/K`
stored/kinetic energy, interval work balance and compensated cumulative accounts
remain the publication gate. Invalid data, cancellation,
native refusal or energy refusal publishes no state or accounting change.
The separate harmonic type cannot be passed to held-force contact replay. Source
validation and physical state acceptance are shared without retaining an unused
held-force exponential. The original `ExactLinearSolid` remains `Send + Sync`;
the generic tone service's erased storage is retained only by `ExactHarmonicSolid`.

For N free coordinates and K supplied frequency slots, admission bounds the
native motion dimension `2*N+2*N*K+1` and independent power dimension
`4*N+2*N*K` by the existing 256-coordinate dense cap. The shared tone owner's
numerical-payload bound is reserved twice to cover its caches and the original
physical matrices, inverse mass and adapter state. Existing controls bound bytes,
native multiplication work and accepted intervals. Repeated frequencies may
share native source columns after this conservative caller-count admission.
Warm steps reuse all exponentials, work/heat factors and phase rotations; only
O(N) output vectors are allocated. Cancellation is checked before and after the
bounded native step, not inside every dense contraction.

Supports remain homogeneous and stationary. There is no missing constrained-row
reaction reconstruction, moving-support model, hereditary constitutive-law
replacement, sparse-scale performance claim or continuum accuracy claim.
Focused regressions in `tests/solid_linear_harmonic_dynamics.rs` cover analytic
resonant growth, a full mutual-mass/nonproportional-damping periodic response,
independently integrated transient beating work/heat, interval subdivision,
coherent duplicate cancellation, DC/rigid motion and transactional refusals.

## Exact complete-coordinate linear solid motion

`dynamics::ExactLinearSolid` advances every free coordinate of an admitted
`TetElasticAssembly` using the existing
`fs-phs::linear_exact::prepared::PreparedQuadraticFlow`. Full consistent mass
and stiffness are retained; no modal truncation, mass lumping or frequency
retuning is introduced. State is `[u,v]`, with original physical storage
`(u^T K u + v^T M v)/2`. Native Cholesky prepares the inverse mass and its
residual is independently admitted against the original matrix.

Preparation accepts explicit initial free-coordinate displacement/velocity,
one fixed positive clock, and `ExactLinearSolidControls`: complete-coordinate
capacity, dense preparation bytes, exponential multiplication work,
operator-error allowance, work allowance in joules and accepted-step cap.
The existing complete augmented dimension is bounded by `3*N+1 <= 256`.
The prepared exponential is reused for warm intervals. Zero-stiffness directions
remain in the state and retain rigid motion and constant-acceleration response.

Force inputs, initial fields and returned mean velocities follow the source
assembly's full `free_dofs()` order. The mean velocity is `Delta u/dt`.
Held-force work is independently `f dot Delta u`; it is not inferred from
the energy difference. Interval and accumulated original physical work balances
gate publication. Invalid operators, exhausted budgets, invalid force, work
refusal and cancellation preserve every accepted state field.

The dense, harmonic, complete-spectrum and prescribed-support owners retain
their original reference displacement and velocity. Both interval and
whole-window storage changes use the original K/M actions and the shared
`fs-math::quadratic` polarized increment, with retained low words before the
work balance is rounded. Rounded endpoint energy totals remain diagnostics;
subtracting them cannot erase a small transverse force, impulse or support
exchange beside a large unchanged kinetic store. Stiffness replacement carries
the coefficient-storage change at the fixed reference coordinates, separately
from parameter work at the current displacement. No tolerance is enlarged and
no missing kinetic change is filled from a source-work residual. Additional
reference storage and scratch are O(n), within the existing conservative
envelopes; original dense operator actions remain O(n²). These are arithmetic
and storage bounds, not a measured speed claim or an arbitrary-conditioning
accuracy certificate.

Focused G3 regressions in `solid_linear_dynamics`, `solid_spectral_dynamics`
and `solid_supported_dynamics` exercise sub-ULP diagnostic energy exchanges,
signed impulse return, support acceleration, physical-work refusal and retry.
`solid_spectral_update` checks a one-joule stiffness change beside a 5e15 J
unchanged store with a 1e-10 J work gate and exact return of that parameter work.

The supported reduced assembly represents homogeneous stationary constraints.
This API neither accepts a full nodal-force vector that silently drops
constrained forces nor reconstructs missing support reactions or moving-support
coupling. It supplies no material damping, nonlinear finite-rotation claim or
sparse exponential scalability claim. Exactness refers to the admitted spatial
pencil and held input, not continuum accuracy or arbitrary time-varying forcing.
The nonlinear conservative owner retains its separate phase limitation.

G0/G1/G5 evidence in `tests/solid_linear_dynamics.rs` passed all five native
tests. A 5 kHz oscillator at a 48 kHz clock runs 1,000,000 accepted intervals
with pitch error below `1e-6` cent and maximum relative energy drift
`5.03e-10`; the independent midpoint formula gives `-59.1076` cents.
A coupled non-diagonal pencil retains both 5 kHz and 11 kHz analytic frequencies
and held-force work. Rational P1 and signed P2 consistent-mass pencils, zero-K
motion, capacity/domain refusal, cancellation and repeatability are covered.
These tests construct independent pencils in the native assembly type; they
do not execute the original geometry assembly frontends. The native lane
compiled the complete original fs-phs crate and was rerun against the current
shared propagator. No full-workspace Cargo result is claimed.

## Original finite-bulk conservative P2/P1 motion

`ConservativeFormulation::OriginalMaterial` extends the same quadratic solid
owner with the ORIGINAL finite-compliance bulk energy. It accepts
`DiscreteBulkSettings`, an explicit relative weak-volume allowance and a local
discrete-volume-work allowance. Quadratic and log-squared laws use their own
native secant inverses. Actual geometric and retained bulk volumes must both
remain inside the caller's positive-volume branch; no clipping, averaged modulus
or surrogate quadratic law is used.

At each reference quadrature point, accepted `Jbar_before` is frozen.
The interval P1 pressure determines `Jbar_after` through
`p_alg = Delta U / Delta Jbar`. The spatial equation is
`integral N_p (det(F_after)-Jbar_after) dV = 0`. Actual geometric pressure
work uses the objective discrete cofactor. Its two distinct pressure/displacement
cross-blocks remain unchanged, and the pressure block is the negative original
secant compliance. The symmetric preconditioner retains that pressure block,
so a finite-compliance body with fully prescribed displacement can still solve
its pressure equations. This does not add a gauge to the exact-constraint mode.

Physical storage is the original remainder at actual F plus original U at
retained Jbar. `bulk_volume_ratios()` exposes those values in reference
element/quadrature order; `bulk_stored_energy_j()` exposes their storage.
The energy record separately reports bulk storage before/after, local bulk work
and its storage defect. Its constraint-work defect subtracts both bulk work and
pressure dotted with the weak-volume change. Bulk constitutive work, pressure
storage work, weak volume, complete physical energy and cumulative work must
all pass before the new displacement, velocity, pressure and volume history
are published together.

`prepare_with_pressure` takes a compatible P1 initial ENDPOINT bulk pressure.
Original `pressure_potential` reconstructs Jbar, and the initial weak-volume
equations must close against the supplied displacement. `prepare` uses zero
initial pressure and refuses incompatible initial geometry. No initial
displacement, velocity or volume is silently projected. In particular, arbitrary
geometric quadrature volumes are not copied into a mixed state: unresolved
pressure-space components would alternate under the quadratic secant.
Compressible initial velocity is admitted. Subsequent accepted volumes follow
the solved original history and are never reset to geometric determinants or
reprojected to endpoint pressure. Interval pressure thereafter remains an
algorithmic work multiplier, not endpoint Cauchy pressure.

G0/G1 tests in `tests/mixed_conservative_bulk.rs` passed all six native
cases, alongside all four full-material and three exact-constraint regressions.
Independent original-energy, rational consistent-mass, momentum and weak-volume
calculations cover freely deforming/spinning bodies, heterogeneous quadrature
moduli, fully prescribed dilation/compression, compatible prestressed
initialization and immutable domain/budget/cancellation refusal. Over 32 free
steps, maximum energy drift was `1.14e-13 J` for quadratic bulk and
`2.09e-13 J` for logarithmic bulk, with nonzero bulk-energy exchange.
The logarithmic trajectory ladder compares displacement, velocity AND retained
volume at 8/16/32 steps against a 128-step reference; observed orders were
1.961 and 2.057. This is bounded time-convergence evidence for the retained
nonlinear history, not a pointwise incompressibility or mesh-independent
stability claim.

`examples/conservative_elastic_motion.rs` now executes FullMaterial,
OriginalMaterial and ExactIncompressible with explicit controls and reports
physical bulk storage separately. All three 32-step native example runs passed.
Contact, material relaxation, moving-support linear exponentials and a
frequency-exact nonlinear finite-rotation split remain separate capabilities.

## Exact affine-force histories on linear solids

`ExactLinearSolid::prepare_affine_forces` prepares the same complete physical
pencil for a fixed clock and an explicit positive characteristic force in N.
`advance_affine(force_start_n, force_end_n, cx)` evolves the force linearly
over each interval. Native `PreparedQuadraticFlow::new_extended_input_controlled`
retains source coordinates `[u,v,f]`; a held force increment enters through
`I/dt`, avoiding a per-step force-slope calculation. These added force
coordinates carry no structural storage or inertia, and no new exponential
implementation or modal approximation is introduced.

Supplied work is independently integrated as `integral f(t) dot v(t) dt`
through the native quadratic observation. It is not inferred from energy change
or replaced by endpoint-average force dotted with displacement. Returned
`mean_velocity_mps` remains the kinematic `Delta u/dt`; only held-force
exchange can use it alone as a work-conjugate flow. The affine source endpoint,
physical interval/cumulative work and cancellation gates all pass before the
physical `[u,v]` state is accepted. Force jumps between intervals are explicit
and change acceleration without introducing a velocity impulse.

The affine complete dimension is `4*N+1 <= 256` (at most 63 free coordinates),
with a distinct explicit cold-memory estimate and the same multiplication cap.
The held route retains its original capacity. Wrong prepared mode, malformed
forces, unrepresentable finite force increments, source error, exhausted budgets
and cancellation refuse without changing accepted state. Exactness applies to
the declared affine force history; approximating another continuous history
still requires temporal refinement.

G0/G1 tests in `tests/solid_linear_force_history.rs` passed all four cases.
Coupled non-diagonal pencils are checked against analytic oscillator responses,
with work evaluated independently by integration by parts. Zero-stiffness
constant-jerk motion uses exact cubic/quadratic trajectories and a polynomial
work integral. One interval and four subintervals agree through continuous
ramps and finite force jumps. Affine-specific dimension/overflow boundaries,
wrong mode, cancellation and deterministic replay are covered. All five
existing held-force tests, including the million-step phase case, remain green.

## Exact coupled viscosity on linear solids

`dynamics::ExactLinearSolid::prepare_viscous` and
`prepare_viscous_affine_forces` retain the original complete free-coordinate
stiffness K and consistent mass M while admitting an explicit symmetric
positive-semidefinite viscosity C [N s/m] in the SAME free-DOF order. They solve
M vdot + C v + K u = f through the existing fixed-clock quadratic exponential.
C can include non-proportional coupling; no modal diagonalization or Rayleigh
assumption replaces its entries. This is a memoryless constitutive resistance
or boundary impedance, not the hereditary tensor-Prony material owner.

The original-dimensional C is checked for symmetry and negative eigenvalues
against the caller's explicit relative operator allowance, normalized by C's
own magnitude with no unit-sized floor. No eigenvalues are clipped. Dynamics
use R_vv = M^-1 C M^-1, retaining Q = diag(K,M) and the original physical
energy. An independent quadratic observation integrates D = integral v^T C v dt.
Held-force work remains f dot Delta u; affine work remains integral f(t) dot v(t) dt.
Interval and compensated cumulative balances enforce Delta H + D = W before
state publication. Negative integrated loss beyond the stated energy allowance
refuses. Energy reports expose interval and cumulative physical dissipation.

Zero viscosity retains the conservative trajectory and existing fixed-clock
semantics. Observation storage/work are included in the same explicit preparation
budgets. Clock, forces, overflow, cancellation and work/loss refusals leave all
accepted physical state, time and cumulative history unchanged. Stationary
homogeneous supports and the existing dense augmented-dimension limits remain;
this adds neither support reactions nor moving-support forcing.

Native `tests/solid_linear_viscous.rs` passes five focused tests: underdamped,
critical and overdamped analytic motion with independently integrated loss;
rigid drag and terminal velocity; coupled consistent M with noncommuting K/C;
affine force work and heat against continuous integrals and partitioned motion;
and invalid/indefinite viscosity, resource limits, zero-C equivalence and
refused/cancelled-step rollback. The five held-force and four affine-force
regressions also pass, including the existing million-step phase check.

## Complete-spectrum linear solid motion

`dynamics::spectral::CompleteSpectralSolid::prepare_with_basis` admits all n
generalized eigenpairs of the original n-coordinate TetElasticAssembly. The
n*n basis is row-major with modes in columns. Every supplied eigenvalue must
be finite and nonnegative and is retained verbatim. Missing modes, duplicate
columns, invalid source matrices or incompatible dimensions refuse. This is a
complete coordinate transformation; no spectral window, modal truncation,
tail substitute, diagonal mass or fitted stiffness enters the motion.

Preparation recomputes K phi - lambda M phi against the original matrices,
using the caller's explicit absolute eigen-residual allowance plus a relative
scale from that mode's actual K phi and lambda M phi. There is no unit-sized
floor. The full row-sum norm of Phi^T M Phi - I must satisfy a declared
tolerance strictly below one; initial u/v reconstruction has a separate
relative allowance. Symmetry is admitted and rounded as in the dense exact
owner. These are numerical residual checks, not interval-certified spectral
enclosures or experimental qualification.

Each fixed-clock mode uses exact oscillator/affine-source coefficients with
stable entire-function series near zero. Zero-frequency modes retain the
exact cubic displacement and quadratic velocity under force ramps without
inverse stiffness. Physical force projection uses Phi^T; all displacement
and velocity coordinates are reconstructed. The independent original-coordinate
work is f_start dot Delta u + Delta f dot (u_after - mean(u)), obtained from
integration by parts along the continuous trajectory, not from the energy
difference or endpoint-average force. Original K/M energy, interval and
compensated cumulative work are checked before accepting modal state, physical
state, time or history. Force jumps between intervals apply no impulse.

The 3n+1/4n+1 augmented-exponential limit does not apply to this owner. Explicit
admission caps cover free coordinates, simultaneous preparation bytes,
preparation and interval scalar products, and accepted steps. Preparation is
O(n^3); retained storage and warm intervals are O(n^2). Cancellation is checked
at rows and vector-sized arithmetic tiles; capacity, overflow, energy and
cancellation refusals preserve the accepted state. This conservative path has
zero physical dissipation and does not silently diagonalize a coupled viscous
operator. General coupled viscosity retains the separate dense exact owner.

The caller supplies the complete basis. The existing
`fs_modal::eigh_gen_dense` can produce all n modes and is used explicitly by
`examples/spectral_elastic_motion.rs`; its small, bounded cold eigensolve has
no internal cancellation polls. No automatic unpolled eigensolve is hidden in
the spectral owner, and no sparse million-coordinate scalability claim is made.

Native verification passes four `tests/solid_spectral_dynamics.rs` tests:
128 coupled coordinates against continuous analytical affine motion and work;
rigid and tiny-positive modes against an independent ODE power series and
partitioned ramps; invalid/incomplete bases including a 1e-30 stiffness
falsifier; and exact capacity boundaries, cancellation, energy refusal and
retry. `tests/solid_spectral_tet10.rs` runs actual 12-cell TET10 assembly,
retains signed consistent mass, obtains all 108 native modes, exercises the
lowest and highest modes, and independently checks original sparse K/M energy
and force work. Measured free-motion energy error is 1.345e-20 J in that
fixture. The example compiles and runs a synthetic nodal pulse/release through
the same production geometry, assembly, spectrum and dynamics. This is
represented-pencil temporal evidence, not mesh convergence or measured-material
validation. The previous nine held/affine and five viscosity regressions pass
alongside these new tests.

### Complete-spectrum viscous motion

`dynamics::spectral::DampedSpectralSolid::prepare_with_basis` adds an original
free-coordinate viscosity matrix C to the same complete M/K/basis admission.
Every column remains, including a rigid coordinate with nonzero physical drag.
The computed symmetric `Phi^T C Phi` must have finite nonnegative diagonal
rates. Each absolute off-diagonal row sum must be no larger than
`maximum_relative_operator_error` times that row's diagonal. A zero diagonal
requires an exactly zero row; a large rate elsewhere cannot mask a coupling to
a weakly damped mode. Negative rates and meaningful nonclassical viscosity
refuse, without clipping or a unit-sized floor. General coupled viscosity
continues to use `ExactLinearSolid`.

Held and affine force intervals use one prepared
`fs-modal::time::LinearForceTransition` per mode. Underdamped, critical,
overdamped and damped zero-stiffness coordinates share the same continuous
entire-function transition, with no division by frequency or pole separation.
Scalar force moments independently integrate supplied work, and prepared
positive velocity Gramians independently integrate viscous heat. Heat is never
an energy residual. Original M/K energies and compensated cumulative work and
heat must satisfy the caller's joule allowance before any state is published.
Impulses use the existing original-mass momentum audit and retain earlier heat
without introducing instantaneous dissipation. A zero C delegates to the
conservative arithmetic. The separate public type exposes no conservative
cubic, moving-support or sustained-contact operation that could omit damping.

The reported heat is the integral of the admitted classical diagonal law.
`relative_damping_projection_error()` reports the maximum off-diagonal row
sum divided by its diagonal. Along the represented trajectory, this bounds
the omitted symmetric modal cross-term work by that value times the reported
heat, before projection roundoff. It does not assert exact equality to an
arbitrary original C, an interval enclosure, or experimental material accuracy.
Preparation requires `96*n^2 + 2048*n` bytes and admits
`6*n^3 + 64*n^2 + 8192*n + 512*sum(s_i)` scalar-product terms, where each `s_i`
halves `max(omega_i*dt, c_i*dt, 1)` down to at most 1/4. Warm intervals retain
the `16*n^2 + 512*n` bound. Cancellation and every admission or arithmetic
refusal preserve accepted motion, time, work and heat.

`tests/solid_spectral_damping.rs` checks all 128 coupled consistent-mass
coordinates against independent continuous underdamped, critical, overdamped
and rigid affine solutions, with independent quadrature for work and heat.
Additional regressions cover terminal velocity and rigid heat, tiny high-Q
loss, impulse history, nonclassical/negative viscosity refusal, cancellation
and arithmetic rollback, and zero-C conservative bit equivalence.

### Updating complete-spectrum stiffness and viscosity without restarting motion

`DampedSpectralSolid::replace_parameters` installs new original-coordinate
K/C and a newly admitted complete eigenbasis on the owner's EXISTING mass,
free-DOF map, coordinate origin and fixed clock. No replacement assembly,
mass, geometry or constraint map is accepted. The new basis must satisfy the
original generalized-eigenpair, mass-norm and reconstruction checks. Physical
displacement and velocity remain bitwise unchanged and are reprojected through
the retained M into every new mode, including reordered/signed modes, rotated
repeated eigenspaces and transitions between true rigid and elastic modes.

The zero-duration event independently contracts signed
`W_parameter = u^T (K_after-K_before) u / 2` from the operator difference,
retaining product roundoff. A large unchanged storage inventory cannot erase
the reported exchange. Positive work is supplied by the parameter owner;
negative work is returned to it. The explicit `maximum_parameter_work_j`
limits the absolute magnitude of this physical exchange. It is separate from
the original numerical residual allowance. Changing C alone adds no work or
instantaneous heat. Original-coordinate before/after storage and cumulative
energy must still pass the existing joule gate before publication.

The original initial-energy datum, kinetic-energy bits, accepted step count,
clock, applied force/impulse work, viscous heat and all their compensation
survive. Parameter work has its own compensated cumulative account and
`cumulative_parameter_work_j()` accessor. Subsequent intervals and impulses
audit `W_force + W_parameter - Q_viscous - (H-H_initial)`. In particular,
switching C to zero retains past heat while using the conservative interval
arithmetic. Neither old loss nor signed parameter work is recast as the other.

The update admits `128*n^2+4096*n` simultaneous bytes for the live owner and
candidate, and `6*n^3+80*n^2+8192*n+512*sum(s_i)` preparation product terms.
Only a completely prepared and audited candidate replaces the live state;
budget, work, operator, reconstruction or cancellation refusal preserves it.
The operation represents caller-timed piecewise constant reconfiguration on
a fixed physical reference space. It does not model moving geometry/supports,
mass transport, equilibrium shifts, continuous coefficient evolution or
constitutive memory. Harmonic and multi-tone source owners retain immutable
underlying-solid views, preventing their prepared source kernels from being
invalidated through this API.

`tests/solid_spectral_update.rs` supplies analytical regressions for a coupled
consistent mass with degenerate rotations and signed/permuted modes, rigid
release and restiffening, zero-amplitude and damping-only work, prior heat
through zero-C affine motion and impulses, independent small parameter work
beside large unchanged storage, and refused/cancelled candidate continuity.

### Continuous harmonic forcing across the complete damped spectrum

`dynamics::spectral::HarmonicDampedSpectralSolid` propagates a prescribed
continuous frequency through every original free coordinate, using physical
cosine/sine force arrays in `free_dofs()` order. The complete M/K/basis and
classical viscosity admission remain those of `DampedSpectralSolid`. Each mode
uses the existing `fs-modal::time::PreparedHarmonicForce` motion and power
owner. Initial transients, exact resonant growth, critical/overdamped decay and
viscously damped rigid coordinates remain in the trajectory. There is no large
augmented solid exponential, mode truncation, Nyquist retuning or substitution
of sampled loads for the continuous forcing law.

`advance` takes quadratures measured from that interval's start. The pure
`advance_force_phase` rotates them through the same prepared phase, while
`advance_continuous` updates caller-owned source arrays only after the original
M/K and cumulative work/heat gates accept the motion. Cancellation, phase
overflow and other refusals preserve both source arrays and all physical
accounts. Source jumps do not imply momentum impulses; native impulses retain
past heat and consume neither a harmonic interval nor source phase. At zero
forcing frequency the cosine load delegates to held-force propagation, with
the sine coefficient physically inactive and DC arithmetic preserved bitwise.

Work includes both source quadratures and the initial transient. Viscous heat
is independently integrated with the shared positive Gramian; it is not
inferred from energy subtraction or cycle averages. The classical-law
projection residual retains the same meaning and boundary as the damped
affine owner. This API supplies one simultaneous frequency with arbitrary
spatial phase, not independent-tone heat summation. Repeated floating-point
source rotations accumulate roundoff; callers may form quadratures from an
absolute clock when that error matters. Temporal/observation resolution and
spatial/material fidelity remain separate requirements.

Nonzero-frequency preparation admits `96*n^2 + 4096*n` bytes and
`6*n^3 + 64*n^2 + 32768*n + 8192*sum(s_i)` product terms, where each `s_i`
halves `max(omega_i*dt, c_i*dt, Omega*dt, 1)` to at most 1/8. These envelopes
include the underlying damped preparation; warm work retains
`16*n^2 + 512*n`. `tests/solid_spectral_harmonic.rs` adds G1/G3/G4 regressions
for 128 coupled consistent-mass coordinates with analytic phase/work/heat,
exact and near resonance, the damped initial transient, DC rigid drag,
source/state refusal atomicity, and impulses following harmonic dissipation.

### Simultaneous tones and affine modulation on the complete spectrum

`dynamics::spectral::MultiToneDampedSpectralSolid` applies a fixed caller-ordered
frequency list to every original free coordinate. Physical force arrays are
DOF-major: `i*K+j` is tone j on free row i, in N. Preparation retains the
original complete M/K basis and the same classical-viscosity admission as
`DampedSpectralSolid`; no mass lumping or mode truncation is introduced.
The shared model-independent `PreparedModalTones` bank receives physical
modal damping rates, including positive drag at zero natural frequency.

Each mode has one total displacement/velocity state. Supplied work and
positive viscous heat retain every cross-frequency term and the initial
transient, using the existing paired-QR power owner. They are never sums of
isolated-tone dissipations. Exactly repeated frequencies coalesce coherently;
near frequencies remain distinct and retain beats. Zero-frequency sources
are physical DC forces. No acoustic or sample-Nyquist gate changes the source.

An explicit preparation flag selects constant or affine quadrature envelopes.
Affine endpoints share the interval-start carrier basis, so attack/release
and changing spatial phases act on the actual continuous force trajectory.
`advance_affine_continuous` publishes the rotated END coefficients as the next
start only after motion and original M/K work/heat gates accept. Constant
continuation has the same transactional rule. Cancellation, invalid arrays,
phase overflow and energy refusal leave physical state, caller source arrays
and all accounts unchanged. Impulses retain earlier heat and consume no source
phase or interval. The immutable `solid()` view exposes the existing state
and energy accounts without allowing a source-law bypass.

The overall solid byte cap reserves `96*n^2+4096*n+32*n*K` for physical
matrices and projection scratch, then clamps the aggregate native bank's
storage budget to the remaining bytes. The caller's smaller tone-storage cap
still applies. Work envelopes are additional: the solid controls cover its
original preparation and `16*n^2+512*n+P*n^2*K` warm terms (P=2 or 4), while
`ModalToneBudget` bounds the complete bank's native preparation and step
work scores. No budget is reset per mode. Preparation polls between native
cost rows, modes and tones; warm projections and mode evaluations poll by row.
One explicitly bounded native power preparation remains atomic.

`tests/solid_spectral_tones.rs` adds G1/G3/G4 cases for 128 coupled coordinates
against analytic phase and cross-frequency work/heat, affine manufactured
motion with direct original-C dissipation quadrature, coherent duplicates,
source/state rollback, global storage admission and impulses after loss.
Exact integration concerns the declared fixed carriers and their envelopes;
it does not resolve continuous peaks or remove aliasing from sparse outputs.

### Retaining continuous force sources through physical parameter changes

`HarmonicDampedSpectralSolid::replace_parameters` and
`MultiToneDampedSpectralSolid::replace_parameters` use the same fixed-mass
physical candidate and separate signed parameter-work account as the damped
owner. They also rebuild every exact source motion/work/heat kernel against
the new spectrum and admitted viscosity before publication. A late source
preparation or capacity refusal preserves the complete previous physical
state, parameter/force/heat histories and source caches.

The retained source frequencies, ordering, interval length and constant/affine
envelope mode do not change. Caller-owned force quadratures keep their phase
and amplitudes; the zero-duration parameter event consumes no source rotation.
Subsequent intervals use the new physical response with the same continuous
carrier. Exact resonance and simultaneous-source interference use their
existing scalar owners. At DC, the harmonic wrapper keeps its original
held-force arithmetic and lower preparation cost.

Nonzero single-frequency updates admit `128*n^2+8192*n` simultaneous bytes
and `6*n^3+80*n^2+32768*n+8192*sum(s_i)` preparation terms. Multi-tone
updates reserve `128*n^2+8192*n+32*n*K` physical bytes PLUS the old native
bank's payload/peak-scratch bound, then clamp the new bank's supplied
`ModalToneBudget` to the remaining global bytes. Thus both old and new
banks are covered during transactional preparation. Native preparation/step
work remains a separate aggregate allowance, with no per-mode reset.

`tests/solid_spectral_source_update.rs` adds independent resonant convolution
and affine multi-tone particular-plus-transient references, with continuous
quadrature for work and combined viscous heat. Cases retain earlier heat,
force work, clock and source phases through retuning, exercise DC after
release, and reject late source-cache preparation, cancellation and
insufficient simultaneous storage without publishing a partial update.


## Exact cubic force histories on the complete elastic spectrum

`CompleteSpectralSolid::prepare_cubic_forces_with_basis` explicitly prepares
piecewise cubic physical forces while retaining every generalized mode and the
original consistent mass and stiffness. `advance_cubic([F0, F1, F2, F3], cx)`
uses `F(s) = F0 + F1*s + F2*s^2 + F3*s^3` for normalized interval time
`s = t/dt`; every coefficient has units of newtons, not newtons per second.
Initial state, mode admission, the fixed time grid and force coordinates retain
the complete-spectrum contract above. No mode truncation, mass lumping, physical
damping, or nonlinear stiffness is introduced by this owner.

Cold preparation caches the entire oscillator functions through E8. A compensated
series handles zero and small phase, including exact rigid polynomial motion;
trigonometric recurrence handles larger phase. The fixed-clock step projects the
four force coefficients, advances all modal positions and velocities, and
reconstructs physical endpoint acceleration and three normalized displacement
moments `Ik = integral_0^1 s^k*u(dt*s) ds`, k = 0, 1, 2, in metres. These moments
allow support work and other polynomial coupling observations to be integrated
without subdividing the physical interval or evaluating extra force samples.

The applied work is independently evaluated as
`F0 dot (u1-u0) + sum_(r=1..3) Fr dot (u1-r*I_(r-1))`. It is compared with
original-coordinate strain plus kinetic energy both per interval and
cumulatively. The accepted state, moments, work accumulation and clock publish
only after all finite-value, work-budget and cancellation checks pass. Internal
read-only trial evaluation permits an enclosing physical coupling to impose its
own full-system gates before the same infallible state publication.

The explicit cubic capacities are `64*n^2 + 2048*n` bytes,
`4*n^3 + 32*n^2 + 4096*n` cold multiplication terms and
`24*n^2 + 1024*n` warm multiplication terms. Existing held and affine preparation
capacities and their arithmetic are unchanged. Source coefficients, candidate
moments and acceleration are finite-checked; refusal preserves accepted history.

`solid_spectral_polynomial` checks the zero/small/series-switch/multiple-period
phase regimes against independent continuous solutions, all 96 coupled
coordinates with a non-diagonal mass, independent physical work, affine
agreement, interval partitioning of all three moments, exact capacity boundaries,
and invalid, overflowing, cancelled and energy-refused candidate rollback.
Focused native source-closure validation also reruns the held, affine, viscous and
complete-spectrum tests. This is floating-point exact linear propagation subject
to admitted source/basis and balance tolerances, not an interval certificate or a
finite-rotation nonlinear dynamics claim.


## Moving supports with original consistent inertia and reaction work

`dynamics::supported::SupportedSpectralSolid` advances the conservative linear
solid with explicitly prescribed support motion. Supply an original assembly
that still contains the prescribed rows and columns. Support labels must be a
nonempty, strictly increasing proper subset of the assembly's represented global
DOF labels; at least one free coordinate is required. Initial displacement,
velocity and applied-force endpoints use the original represented order; support
accelerations and reactions use the supplied support order. Coordinates already
omitted by assembly remain homogeneous stationary supports, and their omitted
reaction rows cannot be recovered by this owner.

`partition_free_pencil` extracts the entire free/free K and M blocks for a native
complete generalized eigensolve. `prepare_with_basis` requires all free modes,
checks the whole source mass for positive definiteness and the whole source
stiffness for numerical positive semidefiniteness, and admits the same source
symmetry allowance as the complete-spectrum owner. A valid free block alone
cannot hide indefinite support inertia or a zero-diagonal indefinite stiffness
coupling. The full K check uses complete diagonal pivoting and the infinity norm
of the entire unresolved Schur residual; no eigenvalues are clipped and retained
operators are not replaced by the admission factors. These floating-point
operator checks are not interval-certified positivity claims.

Each `advance_affine_supports` linearly interpolates retained support acceleration
to the supplied endpoint and integrates it to continuous support velocity and
cubic displacement. It linearly interpolates the full physical applied force.
The free equation uses
`M_ff*a_f + K_ff*u_f = f_f - K_fs*u_s - M_fs*a_s`, retaining the ORIGINAL
consistent cross mass. This effective cubic force is advanced through every free
mode by the exact polynomial spectral kernel. There is no mass lumping, omitted
mode, velocity impulse or repeated dense eigensolve in an interval.

Reported endpoint reactions are the forces exerted BY prescribed supports ON the
solid, `r_s = (M*a + K*u - f)_s`. The original free-row equation residual must fit
its explicit absolute newton allowance. Strain and kinetic energies use full
original-coordinate K/M, including free/support cross terms. Applied work is
integrated from the continuous trajectory. Reaction work is separately integrated
using the exact displacement moments and twice integrating its inertial term by
parts; it is never inferred from the energy difference. Applied loads on prescribed
rows correctly change the reaction and the two work accounts without changing the
prescribed trajectory. The physical closure is
`Delta H = W_applied + W_support` both per interval and cumulatively.

The inner free candidate, original-row residual, independent full work balance,
finite values and physical endpoint energy all pass before either owner's accepted
state changes. Materially negative original strain energy is refused even if a
numerically tolerated stiffness residual produces matching negative support work.
The final cancellation checkpoint precedes infallible joint publication of modal
state, full motion, support acceleration, compensated work histories and clock.
No prepared matrix is cloned during an interval.

Separate whole-system caps admit n represented coordinates, preparation storage
`192*n^2 + 8192*n + 8*element_count` bytes, cold multiplication work
`12*n^3 + 128*n^2 + 16384*n` and warm multiplication work
`64*n^2 + 8192*n`. The inner cubic spectral controls also apply. Preparation is
O(n^3), while retained matrices and intervals are O(n^2). Cancellation is polled
through preparation, positivity factors and interval rows. The basis remains
caller-supplied; the example's small native dense eigensolve is an explicitly
bounded cold operation without internal cancellation polls.

Five `solid_supported_dynamics` regressions check an independently solved
moving-base oscillator with non-diagonal consistent mass, separately integrated
applied and actuator powers, uniform accelerating translation, partition invariance,
constrained-load invariance, original full-pencil rejection, capacities, invalid
sources, cancellation, energy refusal and bitwise rollback. The independent
oscillator's cumulative balance defect is below 3.886e-16 J.
`examples/supported_elastic_motion.rs` assembles an actual 30-coordinate TET10
specimen, prescribes 18 base-face coordinates and retains all 12 free modes. Two
opposite triangular acceleration pulses move the base 62.5 micrometres and stop it;
480 intervals at 48 kHz retain independently observed actuator work and original
physical energy within 9.672e-17 J. Geometry and material are synthetic. This is a
native assembly/solve/motion consumer, not spatial convergence or experimental
calibration evidence. Viscosity, hereditary loss, finite rotations, changing
supports/topology and impulsive support motion are outside this conservative owner.

## Physical boundary loading of hereditary harmonic solids

`ViscoElements::prepare_harmonic_response_with_port_impedance` assembles the
original material spectrum and consistent mass together with a caller-supplied
physical load `Z(omega)` and real power-conjugate force map `B`. The complete
pencil is `K' - omega^2 M - i K'' - i omega B Z B^T`. Reactive and mutual
entries remain present; a loaded resonance is solved before any bare-system
inversion. One factor is reused at that frequency.

The shared harmonic owner supplies pointwise Hermitian-PSD admission and
separate port/projection/payload limits. The load owner retains responsibility
for its frequency domain, units, causality and source accuracy. Returned
`port_load` contains the original flow, opposing effort and net absorbed
boundary work. This work is separate from the original Prony/Kelvin material
heat and is not relabeled as radiation or heat. The signed balance subtracts
both contributions. Original source-band and full-cycle strain admission,
structural payload budgets and cancellation still apply.

`tests/visco_harmonic.rs` checks a load-regularized bare resonance, an independent
coupled spring/mass/resistance boundary on a hereditary solid, separate work and
heat, exact zero-load response parity, and material/load/capacity/cancellation
refusals (G0/G1/G3/G4). These bounded numerical checks do not establish broadband
load passivity, a spatial convergence certificate or measured material accuracy.

### Original hereditary free-decay poles

`ViscoElements::prepare_resonance_poles` prepares the complete homogeneous
fixed-support spatial pencil in `[free displacement; original Maxwell strain]`
coordinates. Every six-component Mandel strain at every original quadrature
point and branch remains present, including equal-time and zero-energy
directions. The existing consistent mass, equilibrium stiffness and complete
Kelvin tensor are retained. Under `exp(s*t)`, the original equations are
`(K_eq+s*C_kv+s^2*M)u + sum(V*B^T*C_j*q_j)=0` and
`(1+s*tau_j)q_j-s*tau_j*B*u=0`. Keeping the branch coordinates makes this
pencil polynomial through material relaxation rates; eliminating them would
introduce poles into the operator and invalidate holomorphic contour counts.

Three real coefficient matrices and fixed two-sided scales are cached once.
Subsequent operator/analytic-derivative evaluations cost O(N^2), with all dense
factorization and contour iterations delegated to the existing `fs-la`
Newton/Beyn owners. Source coordinate, numerical payload, preparation-work and
per-query evaluation limits are admitted explicitly. The source budgets are
separate from the shared solver's iteration/storage/work limits.

Returned modes retain original coordinates and separately check mechanical
force and material-strain residual families against their original coefficient
envelopes. Negative real Laplace rate means decay; a purely real pole is
relaxation, and an unobservable memory state is not an audible resonance.
Eigenvectors have arbitrary common amplitude and are not admitted finite-strain
histories. The original real-frequency material-band intersection is retained
as a conservative `|s|` gate, including the complete contour region. Analytic
continuation of the declared Prony model is not experimental validation at
complex frequency or a material-data extrapolation certificate.

Local refinement needs a supplied structural seed and refuses initialization
exactly at an uncoupled material relaxation rate. Contour discovery keeps the
shared numerical conditioning, nested-quadrature, count/rank and simple-root
admissions. Two-moment Beyn cannot resolve more independent roots than its
probe rank; a quadratic augmented pencil may therefore require contour
subdivision even when every retained coordinate fits the resource budget.
Unresolved multiplicities/defective roots refuse. Numerical residual and
quadrature agreement do not certify spatial convergence, spectral completeness
or rigorous damping enclosures.

`tests/visco_poles.rs` checks an independently factorized scalar Prony cubic,
its oscillatory pair and real relaxation pole, a live contour quadrature ladder,
coupled tensor/Kelvin/consistent-mass roots against an independent determinant,
source/resource/cancellation refusals, lossless roots and full-spectrum rank
refusal (G0/G1/G3).

### Physical boundary stiffness, damping and inertia in hereditary solid modes

`ViscoElements::prepare_resonance_poles_with_port_load` adds a complete real
power-conjugate port map B and symmetric numerical-PSD boundary coefficients
K_port, C_port and M_port to the original hereditary spatial pencil. The
opposing effort is `K_port*d+C_port*d'+M_port*d''`, where `d=B^T*u`; the
body receives `-B*effort`. The mechanical pencil therefore receives
`B*(K_port+s*C_port+s^2*M_port)*B^T` before balancing. Its causal impedance
is `K_port/s+C_port+s*M_port`. No single-frequency impedance sample is
continued, and no boundary mass, mutual entry or source Maxwell coordinate is
discarded. This represents kinematically attached mass, springs and dampers;
an attachment with independent moving coordinates requires its own state.

The caller supplies the complete original free-coordinate force map, all
three port matrices, a declared real-frequency band, and additional port,
payload and preparation-work limits. Each coefficient is checked for exact
symmetry, finite values and numerical positive semidefiniteness under a
diagonal congruence, reusing the bounded shared fs-la Jacobi eigensolver.
Negative diagonal entries and nonzero rows with zero diagonals refuse before
numerical tolerance could conceal an active direction. The original matrices
are retained unchanged; singular passive directions are allowed. These
floating-point checks are not interval-certified passivity statements.

The load's additional payload envelope is `24*n*p+80*p^2+256*p` bytes. Its
preparation-work envelope is `3*(720*p^3+256*p^2)+6*n*p*(n+p)` scalar-work
units; every individual PSD admission also obeys fs-la's unpolled-work cap.
These limits are separate from the original source-pencil and Newton/Beyn
budgets. Port projection is prepared once; repeated pole and contour queries
reuse the original quadratic coefficients and fixed scales. Cancellation is
polled before/after bounded PSD calls and at projection/observation rows.

Each returned pole includes optional original port displacement, flow and
opposing effort, load storage `(d^H*K_port*d+flow^H*M_port*flow)/2`, boundary
heat `Re(flow^H*C_port*flow)`, net absorbed power `Re(flow^H*effort)`, and the
signed defect `power-2*Re(s)*storage-heat`. These are complex-envelope forms
at the mode's arbitrary normalization, not finite-cycle averages or admitted
finite-amplitude trajectories. Boundary storage and loss remain separate
from intrinsic material loss. The original mechanical and memory residuals
include the added load. The intersection of original material and load bands
gates `|s|` and the whole contour; analytic continuation is not measured-data
validation at complex frequency. Existing numerical-root, multiplicity,
completeness and spatial-accuracy limits remain unchanged.

`tests/visco_loaded_poles.rs` checks an independently factorized loaded Prony
cubic, original material memory ratios, independent boundary power/storage/heat,
two mutually coupled ports against an independently assembled rational
determinant with consistent cross mass, exact zero-load pole/vector parity,
singular passive directions, source-band admission, active and malformed
loads, exact resource boundaries, and cancellation (G0/G1/G3/G4). The full
native fs-solid library builds, all four loaded-pole tests pass, and the five
existing hereditary-pole regressions pass against the updated library.

## Full hereditary solid poles with independent dynamic physical loads

`ViscoElements::prepare_resonance_poles_with_dynamic_load` assembles the complete
original body/material/load pencil. The supplied `ViscoPoleDynamicLoad` retains
every coordinate of an admitted exact quadratic `PortHamiltonian`, including
mechanically dark states, full mutual storage/resistance and all original port
columns. For body displacement `u`, the same signed physical map `B` gives load
flow `v=B^T*u'` and opposing nodal force `-B*effort`. With the original load law
`z'=(J-R)Qz+Gv`, `effort=G^T*Qz`, the added equations are
`(sI-(J-R)Q)z-s*G*B^T*u=0` and the body receives `B*G^T*Qz`. No standalone load
resolvent, sampled impedance, material denominator, projection or modal damping
ratio replaces those rows. The joint polynomial is regular to evaluate at bare
load and Maxwell relaxation poles. Its bounded dense coefficients and fixed
scales are prepared once, while all root iterations remain with `fs-la`.

The existing complete-state dimension/payload/evaluation/solver budgets include
the appended load states. `ViscoPoleDynamicBudget` separately admits original
load state/port counts, payload, projection work and the shared numerical PSD
check. Exact finite quadratic storage, symmetric numerically PSD Q/R and the
original admitted pHS structure are required; nonlinear storage is refused.
Source and load bands are intersected and applied conservatively to `|s|` and
the entire contour. This is a declared causal continuation, not experimental
evidence outside measured frequencies. Preparation and observation poll `Cx`,
and no accepted material/attachment history is read, advanced or reset.

`physical_vector` is original free displacement, all original six-component
Mandel branch states, then the independent load's original state order.
`dynamic_load_range()` exposes the suffix. Each pole must pass separate
original coefficient-scaled mechanical, Maxwell and load-state residuals.
`resonance_pole_from_state` accepts a full authored seed, including dark load
or memory modes, without a displacement-only resolvent initialization.
The convenience displacement seed initializes load coordinates to zero and
retains its existing exact Maxwell initialization/refusal behavior.

The dynamic-load observation independently reports `H=z^H Qz/2`,
`D=(Qz)^H R(Qz)`, incoming port power, and `P-2 Re(s)H-D`. Body equilibrium,
consistent-inertia and every Maxwell storage term are evaluated separately,
as are Kelvin and Maxwell dissipation. The joint signed defect is
`2 Re(s)*(H_body+H_load)+D_body+D_load`. These are complex-envelope quantities
at the returned arbitrary mode amplitude, not cycle averages or fitted heat.
Signed numerical defects remain visible. Finite pencil roots need not be
observable mechanical resonances. Local/contour results remain numerical,
without interval enclosures, global completeness, nonlinear-load support,
measured material validation or a new spatial-convergence certificate.

G1/G3/G4 specimens in `visco_dynamic_poles.rs` verify analytical attached-mass
mode splitting, a seed exactly at the bare attachment pole, 128/256-node contour
agreement and honest 64-node count refusal, hereditary loaded roots against an
independent scalar eliminant, separate original body/load energy and heat,
signed mutual ports, mechanically dark original load states, nonlinear-law
refusal, resource/band refusal and cancellation. Native fallback verification
compiles the complete `fs-solid` crate against aligned cached dependencies;
it does not claim a DSR or full-workspace run.

### Forced response of the complete hereditary solid and dynamic load

`PreparedViscoPoles::harmonic_factor` factors its original conditioned
polynomial at `s = -i omega`, retaining every original free displacement,
six-component Maxwell strain state and independent quadratic pHS load state.
It accepts positive finite frequencies in the intersection of the original
material and attached load bands, including regular coupled problems at a bare
host resonance or bare attachment resonance. It never evaluates an isolated
load inverse, condenses away a load/material state, fits a replacement damping
law, advances accepted histories or introduces a homogeneous state.
`ViscoJointHarmonicFactor::solve` reuses that immutable shared `fs-la` complex LU
for multiple forcing columns and phases. Original mechanical, memory and load
rows are independently checked with the prepared physical residual tolerance;
up to three corrective solves reuse the same factor without changing that
tolerance. Genuine joint singularities and unrepresentable states refuse.

Force and state phasors are peak amplitudes for `exp(-i omega t)`. Every original
material quadrature point undergoes the same full-cycle strain-envelope
admission as `ViscoHarmonicOperator`; separate real/imaginary amplitude checks
cannot admit an excessive intermediate phase. The returned body storage uses
the original equilibrium stiffness, consistent mass and Maxwell strain energy;
material heat uses the original Kelvin tensors and every Maxwell relaxation
branch once. Static K/C/M and independent dynamic loads retain separately
observed original storage and dissipation. Body energy is never obtained by
subtracting a much larger attachment energy from a combined account.

All reported energies, heat rates and powers are cycle means. Applied work is
`Re(force^H velocity)/2`; absorbed interface work is
`Re(flow^H effort)/2`, positive into the load. Signed body, total and individual
load power defects remain visible. Power admission uses the unchanged physical
tolerance and the sum of apparent applied/interface work, actual dissipation
magnitudes and `2 omega` times original mean storage magnitudes. This physical
oscillatory work scale admits exact antiresonances with zero interface flow and
finite internal load motion; it adds no constant floor, clips no signed power,
and does not classify reactive storage as heat. Passivity restrictions remain
those of the original source and load preparation.

`ViscoJointHarmonicBudget` bounds complete coordinate count, numerical payload
and conservative factor/solve/refinement/observation work, separately from the
already admitted immutable source cache. Polling occurs during assembly, shared
LU/solves, original-row checks, quadrature admission and observation. Refusals
leave both the factor and accepted source histories unchanged. Determinism is
fixed-order same-build/ISA behavior; numerical budgets exclude allocator
metadata and process RSS. This is a dense bounded response path, not a sparse
large-model performance claim, an interval enclosure, a nonlinear harmonic
balance method, a zero-frequency solver or support for descriptor/nonquadratic
independent loads.

Focused G1/G3/G4 evidence in `tests/visco_joint_harmonic.rs` covers analytic
coupled-oscillator states at bare attachment antiresonance and bare host
resonance; full Maxwell/Kelvin and mutual load dissipation across five forced
frequencies; resolving body energy beside an attachment more than ten orders
larger; phase/repeated-RHS reuse; original full-cycle strain refusal; singular,
invalid-band, resource and cancellation refusal; and subsequent factor reuse.
Native compilation of the complete `fs-solid` crate and the four new tests,
nine existing `visco_harmonic` tests and four existing `visco_dynamic_poles`
tests passed. These are focused native checks, not a workspace or DSR gate.


### Original constrained dynamic attachments

`prepare_resonance_poles_with_descriptor_load` accepts original quadratic `DescriptorPortHamiltonian` loads through `ViscoPoleDescriptorLoad`. It retains `E=diag(I_diff,0)`, effort `e=[Q_diff*z_diff,z_alg]`, every differential state, every zero-storage multiplier, and every original Maxwell strain. Velocity input `B^T*u'` and opposing force `-B*G^T*e` preserve power conjugacy. No descriptor matrix, isolated load impedance or constraint is inverted or regularized.

Finite poles and forced harmonic response reuse the existing prepared pencil and shared numerical owners. Original state ranges distinguish differential states and reactions. Multipliers contribute to port power and full mutual resistance heat but never acquire stored energy. Every original algebraic row independently passes the physical residual gate, including same-factor harmonic refinement. An all-differential descriptor preserves ODE arithmetic and state order. Full source bands, state/port/payload/work limits and cancellation apply. Exact quadratic storage is required; descriptor index and infinite-spectrum completeness are not certified; redundant or singular constraints remain explicit refusals.

Native full `fs-solid` compilation and all 12 focused regressions passed: four descriptor tests, four original joint harmonic tests and four original dynamic-pole tests, using the actual source and cached production dependency artifacts. This focused validation is not a full workspace Cargo gate or a descriptor-index/infinite-spectrum certificate.


### Original continuous state-space boundaries in hereditary solid resonance

`ViscoPoleStateSpaceLoad` and `prepare_resonance_poles_with_state_space_load` retain the complete original real law `z'=A*z+B*v`, `effort=C*z+D*v+E*v'` with `v=force_map^T*u'` and opposing body force `-force_map*effort`. Full mutual D and derivative E terms enter the existing first- and second-degree mechanical coefficients. Every original state, including conjugate-pair coordinates and mechanically dark states, remains in the full polynomial; zero-state static D/E boundaries are supported. No isolated load impedance, source resolvent or derivative/direct matrix is inverted, fitted or condensed. Shared Newton/Beyn and reusable complex LU remain the only numerical solver owners.

The caller owns the original realization, state/port units, real-frequency support and stability/passivity evidence. Source and boundary bands are intersected and applied conservatively to |s|. Complete state/port/payload/preparation budgets and cancellation apply before copying the load. Each original boundary-state equation is checked independently within the physical residual/refinement gate. Original body Maxwell/Kelvin storage and heat remain separate from the boundary.

A generic state-space realization supplies no physical storage/dissipation witness. Pole and forced observations therefore expose original states, flows, efforts and signed per-port/total absorbed work, including negative mutual exchanges. They do not label that work as heat or invent boundary energy. Forced `average_total_storage_j` and `total_power_balance_residual_w` are now `Option<f64>`: the existing complete physical pHS/KCM paths return `Some`; an unqualified state-space boundary returns `None`. The original body energy/heat and body-plus-boundary-work balance remain independently computed and enforced. Numerical roots do not certify source passivity, experimental complex-frequency support or spectral/spatial completeness.

G1/G3/G4 regressions in `visco_state_space_loads` cover bare oscillator antiresonance with original derivative inertance, exact loaded free modes, signed mutual D/E work beside all Maxwell/Kelvin strains, original dark-state finite poles, finite-amplitude admission, resources, support, cancellation and factor reuse. Existing fully physical response tests explicitly unwrap complete energy accounts.

### Instantaneous momentum in complete linear-solid trajectories

`ExactLinearSolid::apply_momentum_impulse` and
`CompleteSpectralSolid::apply_momentum_impulse` accept resolved impulses in
the original free-DOF order. The dense owner retains the inverse already
computed and admitted for its source map; the complete-spectrum owner uses
`Phi^T J` and all basis columns, including zero-stiffness rigid modes. Each
event uses O(n²) operations and O(n) scratch, without a new factorization or
propagator rebuild. The dense owner retains one additional n² matrix within
its conservative preparation envelope. These are operation/storage bounds,
not measured performance results.

Both owners audit `M * (v_after - v_before)` using the actual represented
endpoint and an explicit absolute per-component impulse allowance. An admitted
approximate basis cannot bypass this physical audit. Signed work is
`J dot (v_before + v_after)/2`; its independent residual uses the momentum
residual directly, with compensated product accounting. Work enters the
original cumulative account once. Displacement, clock, accepted-step count,
viscosity/source preparation and earlier physical losses remain intact.
Zero impulses preserve physical/modal bits; cancellation, unrepresentable
kicks and failed momentum/work gates publish nothing.

These reduced pencils cannot recover discarded support reactions or claim a
full-body angular inventory from an incomplete DOF map. The API applies a
resolved impulse; it does not select contact laws, restitution or event times.
G1/G3/G4/G5 regression definitions cover off-diagonal mass, viscous/affine
continuation, a free tetrahedral mass's independent linear/angular balance,
mixed rigid/elastic spectral response and analytic continuation, zero identity,
rounded-basis work refusal and cancellation. Native execution remains outstanding.

### Impulses and instantaneous support velocity changes

`SupportedSpectralSolid::apply_momentum_impulse` accepts applied impulses
in the original represented-coordinate order and the prescribed support
velocities immediately after the event. Positions, clock, accepted-step count,
prepared operators and retained support acceleration remain unchanged.
Subsequent continuous intervals start from the new support and free velocities.

The original consistent mass supplies
`J_eff = J_f - M_fc * Delta v_c` to the complete free spectral owner.
Support reactions are `R_c = M_cf * Delta v_f + M_cc * Delta v_c - J_c`.
Both cross blocks remain active; applied impulses on support rows change the
actuator reaction without changing the prescribed motion. The actual represented
full velocity is independently checked against the original mass in every row,
including whether large opposing applied/support impulses can represent their
small required net transfer. A failed outer audit also discards the unpublished
inner spectral candidate.

Applied and support work are separately computed from their impulses dotted
with average pre/post velocity. Their signed cumulative accounts continue through
later native intervals. Compensated product accounting and direct momentum
increments avoid subtracting large momentum totals. The existing rounded
total-energy reference convention remains unchanged, so an exact zero event
does not invent a debt when one store is below the ULP of another.
The event borrows all prepared matrices and uses O(n²) operations with O(n)
scratch. These are algorithmic bounds, not measured performance.

G1/G3/G4 regression definitions cover independent consistent-mass reactions,
actuator work, analytical post-event continuation, uniform boosts, support-load
partition invariance, zero identity across separated energy scales, and late
full-row refusal/cancellation followed by an unchanged native continuation.
This API prescribes resolved impulses/support velocities; it does not choose
contact, restitution, event timing or a compliant actuator law. Native test
execution remains outstanding.

### Simultaneous normal impacts on complete solid trajectories

`ExactLinearSolid::apply_normal_impacts` and
`CompleteSpectralSolid::apply_normal_impacts` resolve a declared simultaneous
frictionless contact set through the existing `fs-contact` owner. Dimensionless
normal covectors use the original free-coordinate order. The prepared inverse
or complete basis supplies the full Delassus matrix; every rigid and elastic
mode remains represented. Restitution, complementarity, PSD/rank admission and
relative-work limits remain owned by the shared contact law.

One native impulse candidate checks the actual original-mass momentum change.
Actual endpoint normal speeds, complementarity, prescribed-surface work,
relative work and applied-work closure must pass before physical or modal
state commits. Contact loss is reported separately and is not added again to
native dissipation or assigned to heat. Displacements, clocks, prepared
operators and material history stay fixed. Explicit adapter scratch and
response-assembly limits accompany the contact solver's bounded iterations;
there is no prepared-matrix clone or refactorization.

A redundant result retains the native impulse-representative qualification:
unique velocity/relative work does not grant authority to individual contact
loads or source-work allocations. Geometry/event discovery, friction,
sustained contact and thermal loss partition remain separate owners.
Focused regression definitions cover simultaneous cross-contact response,
moving surfaces, rigid/elastic continuation, redundant contacts and boosts,
independent linear/angular impulse of a barycentric P1 tetrahedral point load,
and late endpoint refusal with unchanged-state retry and cancellation/budgets.

### Simultaneous normal impacts with prescribed solid supports

`SupportedSpectralSolid::apply_normal_impacts` composes the same complete
normal-contact owner with the original supported solid. Contact covectors use
all `represented_dofs()`, including prescribed coordinates. The accepted
support velocities remain fixed during the instantaneous event. Their real
velocities enter contact approach and work while their incremental mobility is
zero; the free block uses every original spectral mode, not a lumped mass or
sequential contact kick.

The original consistent mass then supplies each support reaction as
`M_cf Delta-v_f - J_c`. Cross mass, applied contact work, prescribed-obstacle
work and support-actuator work remain separate. A passive relative impact may
increase the solid's kinetic energy through physical work by a moving support.
The shared law's rank/determinacy classification remains authoritative:
redundant contact impulses and their source-work allocations are checked
representatives, not unique local loads or heat partitions.

The original support impulse implementation now exposes a private trial/accept
seam. Contact endpoint, complementarity and relative-work checks and the full
original momentum/energy gates precede publication of either the supported
state or its free spectral trajectory. Preparation, material history,
positions, support accelerations and clocks are retained; continuous motion
resumes from the accepted velocities. Scratch and response products use the
explicit shared contact limits; original impulse and native-solver budgets
retain their separate scopes. No prepared operator is cloned or refactored.

Focused regression definitions cover analytic cross-mass reactions and moving
support work followed by exact continuation, a moving opposing surface under a
common velocity boost, redundant/locked contact sets, cancellation/resource
refusal, and a late full-system work refusal with unchanged accepted state.
Geometry/event discovery, friction, sustained-contact integration and conversion
of contact loss to physical heat remain outside this instantaneous normal law.


### Simultaneous Coulomb impacts on complete solid trajectories

`ExactLinearSolid::apply_frictional_impacts` and
`CompleteSpectralSolid::apply_frictional_impacts` admit up to 21 declared
normal-first orthonormal contact frames through the existing simultaneous
Signorini–Coulomb owner. Each frame supplies three original-coordinate
covectors, a stable identity, an opposing-surface velocity, friction
coefficient and explicit normal endpoint law. Geometry owns the frame,
shape functions, event time and contact-set completeness.

The prepared inverse or complete spectral basis supplies every entry of
the full 3N-by-3N Delassus matrix. Normal/shear and intercontact coupling
remain in one original non-associated Coulomb solve. No sequential
normal/friction clipping, physical compliance, modal truncation, new
factorization or prepared-operator clone is introduced.

An original native impulse trial precedes publication. Actual represented
velocity independently checks predicted kinematics, normal complementarity,
the Coulomb cone, final-slip traction and tangential duality using the
original contact tolerances. Original-mass momentum and native work gates
remain active. Opposing-surface work, relative midpoint work and applied
work close independently with compensated arithmetic. Total impact loss
is not added again to native viscosity or assigned a heat destination.

Stable contact IDs determine transpose and work accumulation order.
Explicit adapter workspace and response-product limits accompany the
contact solver's separately bounded iterations and capacity. Cancellation
or any failed law, momentum, velocity or work gate preserves physical and
modal state, positions, clock, preparation and cumulative accounts.

The native determinacy and total-work witnesses remain authoritative:
redundant impulses are checked representatives, not unique individual
loads or local source/heat allocations. Geometry/event discovery,
sustained contact and physical heat partition remain separate owners.

Focused G1/G3/G4 regression definitions cover coupled sticking/sliding and
analytic all-mode continuation, independent P1 consistent-mass world
linear/angular impulse, redundant contacts under boosts and tangent-frame
rotation, and late native-readback refusal with unchanged-state retry.
Native compilation and test execution remain outstanding.

### Directional Coulomb impacts on complete solid trajectories

`ExactLinearSolid::apply_elliptic_impacts` and
`CompleteSpectralSolid::apply_elliptic_impacts` extend that same native impulse
path with two principal friction coefficients and their material axis in the
original tangent frame. Shared assembly retains the complete original mass
response, normal/shear and intercontact blocks, canonical transpose order,
private native momentum trial and physical work accounting. The circular entry
points retain their original scalar dispatch and arithmetic. The existing
`fs-contact` elliptic owner alone maps the directional cone to its scalar disk;
equal coefficients select its exact identity map, and a zero principal
coefficient gives physical free slip without an artificial traction floor.

Native endpoint readback independently reconstructs the original material
ellipse, forbidden free-axis traction, normal complementarity, final-slip
support gradient and tangential duality. It also checks actual native velocity
and relative work against the mapped solve. Loosening the adapter readback
allowance cannot silently loosen these original physical law gates. Actual
opposing-surface work and native applied work retain their separate balance.
All refusals leave physical/modal state and cumulative accounts unchanged.

The law's determinacy vectors and total-work witness remain explicitly in its
mapped coordinates; its physical impulse/velocity outputs remain in the
original contact frame. These are caller-declared fixed impact coefficients,
not inferred material calibration, intrastep static/kinetic switching,
sustained contact, event discovery or a frictional heat partition.

Four G1/G3/G4 regression groups in `tests/solid_frictional_impacts/elliptic.rs`
cover the isotropic limit on both complete owners, a manufactured coupled
stick/slide event with independent all-mode continuation, physical texture
rotation and free-axis frame covariance, and late native refusal/retry.
These integration tests require the real constellation runtime and native
solid owners; formatting alone does not establish that these gates pass.

The same owners' `apply_sourced_elliptic_impacts` methods bind original
`SourceThermalEllipticCoulomb` laws by contact ID and input order. Each active
contact samples both pinned principal sources from its actual native incoming
relative tangent velocity and explicit temperature/pressure/query/impact role.
The source sampler owns domain, units, authority, claim identity and original
history admission. Every accepted path coordinate, native slip, transported
texture and finite-step temperature survives the instantaneous source refresh.
Inactive rows remain mechanically present but do not query or alter history;
malformed gap/tolerance values cannot masquerade as inactive contacts.

The reusable callback entry performs queries in canonical ID order only after
native kinematic and mobility assembly. Returned histories remain private
through the complete native impulse and physical-law/work audit; all body and
history publication follows the final successful gate without fallible work.
Accepted histories need not implement `Clone`. Callback-owned memory/evaluation
budgets, interior mutation and external effects remain outside this candidate
transaction. No source pressure, time increment, travel or impact-heat partition
is inferred. Three G1/G3/G4 groups in the nested source tests cover actual
relative-speed sampling and explicit roles on both complete owners, inactive
history preservation, ID/order covariance, material refusal, invalid activation,
and late native refusal with unchanged-state retry. These share the same native
integration lane and original runtime requirements above.


### Simultaneous Coulomb impacts with moving solid supports

`SupportedSpectralSolid::apply_frictional_impacts` uses the same original
non-associated normal/two-tangent `fs-contact` law as the complete free solid
owners. Contact covectors include every represented coordinate. Prescribed
supports retain their accepted instantaneous velocities: they contribute to
actual approach/slip and work, while their incremental compliance is zero.
Normal and frictional supported events now share one private full-coordinate
mobility action using every free spectral mode.

The complete contact matrix is `J_f M_ff^-1 J_f^T`, retaining normal/tangent
and inter-contact coupling. Original consistent mass supplies
`R_c = M_cf Delta-v_f - J_c`, including direct contact impulses at prescribed
coordinates. Applied-contact, opposing-surface and support-actuator work are
separate. Actual final slip, normal/complementarity and Coulomb residuals,
relative passivity, original momentum and whole-solid energy all gate one
unpublished native impulse candidate before both owners accept.

The original trajectory keeps displacement, time, support acceleration,
prepared matrices, material state and prior cumulative work. Contact loss is
not assigned to heat or counted again as viscous dissipation. Native contact
determinacy and redundant-impulse representative qualifications remain intact;
a uniquely determined velocity does not establish unique local contact loads
or source-work partitions. Adapter workspace/product limits, native solver
limits and original supported-state limits retain their separate scopes.
No prepared operator is cloned or refactored.

Four focused G1/G3/G4 regression definitions cover analytical sliding and
sticking against moving supports, original cross-mass reaction and work,
independent elastic continuation, common-boost work transformation, locked
contact refusal, cancellation/budgets and a late original-energy refusal
followed by an unchanged-state retry. These four regressions pass in the
focused supported-dynamics native run described below.
Geometry/event discovery, sustained friction and physical loss destinations
remain separate owners.

### Directional material impacts with prescribed supports

`SupportedSpectralSolid::apply_elliptic_impacts` carries the original physical
friction ellipse through the same complete supported momentum trial as the
circular law. Both paths share native reaction, work and acceptance code.
Contact covectors include all represented coordinates; actual support velocity
enters incoming and outgoing slip, while incremental mobility is exactly zero
on prescribed coordinates. The complete free-mass response retains every
normal/shear and cross-contact block. Physical impulses return through the
original transpose, and support reactions retain both `M_cf Delta-v_f` and
direct contact loads `-J_c`. Contact work, opposing-surface work and support
actuator work remain separate. The original physical endpoint ellipse,
zero-axis forbidden traction, normal law, relative passivity and full native
momentum/energy gates all precede publication. Equal principal coefficients
retain the circular owner's exact arithmetic; a single zero coefficient gives
true free slip without an artificial traction floor.

`apply_elliptic_impacts_with` and `apply_sourced_elliptic_impacts` reuse the
original complete-solid source transaction. Active sources are sampled in
canonical contact-ID order at actual native incoming relative velocity;
history slices retain original input order, and inactive sources are neither
queried nor changed. The typed entry preserves both pinned principal claim
sets, original source applicability/authority, explicit temperature/pressure
and coefficient role. Accepted slip, accumulated distance, transported texture
and last finite-step temperature survive the instantaneous coefficient refresh.
Source histories need not implement `Clone`. They and the supported/free
trajectories publish only after every original law and full-support work gate;
source, cancellation, budget or late energy refusal leaves all accepted state
available for retry. Callback-owned allocation/evaluation and external side
effects retain the existing source entry's boundary.

The new nested tests in `tests/solid_supported_dynamics/frictional_impacts/elliptic.rs`
and its `source.rs` exercise the isotropic limit, oblique maximum-dissipation
traction, exact free-axis slip and frame rotation, two simultaneous coupled
stick/slide contacts, cross-mass and direct support reactions, independently
computed source work, subsequent original elastic motion, physical source
queries with transported material axes, and source/native refusal with retry.
All 24 tests in the original supported-dynamics integration target passed in a
focused native build, including these eight new cases and the existing normal,
circular-friction, impulse and continuous-support regressions. The run compiles
the original dynamics/source/constitutive modules and real `Cx`/`CancelGate`;
it does not establish a full Cargo workspace or DSR gate. These are numerical
model checks. The instantaneous entry does not infer
contact discovery, pressure, duration, sliding travel, intrastep static/kinetic
selection, a frictional heat partition or measured material calibration.


## Simultaneous reciprocal contact across complete native solid trajectories

`dynamics::contact_network` supplies one bounded, simultaneous Coulomb event across exclusively borrowed complete exact, spectral and supported solid trajectories. Geometry-owned body-local point maps use one common orthonormal contact frame; a finite pair automatically receives reciprocal transpose impulses. Independent original mass blocks feed the existing shared frictional-impact/`fs-contact` law owner, preserving complete cross-contact and normal/shear mobility without cloning or refactorizing prepared operators. All bodies must share the identical accepted event time. Explicit body, total-DOF and additional workspace limits supplement the existing adapter and solver capacities. Zero input blocks skip their native mass action.

Every original native trial, actual endpoint contact-speed and Coulomb-law audit, contact-work identity, whole-network applied-plus-support-work balance, and final cancellation check precedes publication of any body. Original per-body receipts preserve supported cross-mass reactions and prescribed support work separately from opposing-surface work and relative contact loss. No contact loss is silently allocated to heat or native viscosity. Geometry owns complete simultaneous contact discovery, common frames and physically consistent coincident-point interpolation; arbitrary algebraic point maps alone do not certify world linear or angular momentum. Determinacy and total-work witnesses retain the shared law owner's meaning.

Five focused regression groups cover finite-pair transfer and original trajectory continuation, a three-body coupled closure, consistent-mass tetrahedral world momentum at two origins, supported source work under a common velocity shift, and multi-owner refusal/retry. New files passed Rust formatting checks; edited parents passed syntax parsing. Small independent arithmetic checks covered analytic fixture impulses, work and tetrahedral mobility. Native compilation, Rust test execution, convergence and performance benchmarking remain outstanding in this resource-constrained environment.

## Original layer Maxwell memory in native frame resonance

`Frame2d::prepare_relaxing_resonance` binds authored longitudinal layer spectra
to the same complete native `PreparedFrameResonance` used for physical boundary
AC and poles. The original source resolver retains equilibrium modulus, every
delta_E/tau pair, selectors, receipts and held-environment applicability. The
explicit longitudinal declaration applies one scalar E(t) to axial strain and
both bending strains; it does not reinterpret a bending-only or torsional law.
The same reference member coordinates supply the geometry and conjugate force maps.
Their exact coefficients are delta_E*A/L, 3*delta_E*I/L^3 and delta_E*I/L^3.

Each positive original arm/channel has an energy-normalized elastic extension
z, with z_dot=-z/tau+sqrt(k)*extension_dot, opposing force sqrt(k)*z,
H=z^2/2 and heat=z^2/tau. No modulus is frozen at a resonance frequency and no
modal loss ratio is fitted. Original equilibrium stiffness enters once.
Different layer times remain distinct; exact-zero increments retain source
identity without a state, while positive arms with zero free-motion projection
retain their physical energy/decay states. Selected members with nonzero
Kelvin–Voigt material loss refuse an unspecified second material calibration.
Separately authored Rayleigh loss and other members' native material loss remain.

The optional physical boundary retains all original states, mutual structure,
ports and descriptor constraints. State order is
`[q,p,material_memory,original_boundary]`; no algebraic multiplier acquires
storage. Full consistent translation/rotation inertia and frame rigid modes
remain with the existing native coupled owner. The adapter adds no integrator,
factorization or eigensolver. Preparation charges constitutive composition and
native work/payload together, with conservative pre-admission of three states
per declared arm and cancellation between bounded operations.

Original material frequency bands are intersected and installed in the native
AC factor gate. Complex poles describe the declared continuous Maxwell law;
finite real-frequency source support is provenance, not experimental authority
for its complex continuation and not an Im(s) root filter. Local and requested
contour queries retain the shared solver's completeness qualifications. Boundary
source applicability remains with its owner. Fixed temperature/reference geometry,
perfect bonding, small strain and Euler–Bernoulli kinematics apply. This prepares
a linear constitutive model without importing transient history or depositing heat.

`harmonic_memory` reconstructs each original extension, elastic/viscous split,
opposing force, mean storage and physical material heat from the same peak native
phasors. The material account is already part of native total storage/dissipation;
it is not added again or confused with boundary work or Rayleigh damping.

Three native tests of the unmodified production coefficient/observation helper
passed against the existing production `fs_math` artifact: original channel/rate
integration, dynamic storage/loss, and all three analytic coupled Maxwell cubic
rates. Six full integration regression groups define source-bound pole/heat
closure, distinct layer times and zero identity, both bending traces under full
translation/rotary inertia and orientation changes, descriptor inertance/drag,
zero-only/empty spectra, and source/budget/cancellation refusals. New modules
passed pinned Rust formatting and syntax parsing. Full dependency compilation
and those integration executions remain outstanding in this constrained session;
the leaf checks are not a full native frame validation or performance benchmark.

## Independent shear Maxwell memory on the original Timoshenko basis

`Frame2d::prepare_relaxing_resonance_with_transverse_shear` composes original
longitudinal and independent transverse spectra in the same native resonance
and exact transient owner. Each selected Timoshenko member requires a
`FrameTransverseRelaxation` binding for every layer; an explicit empty spectrum
declares elastic G. Missing G, foreign/duplicate selections and nonzero selected
axial, bending or shear Kelvin loss refuse. Euler--Bernoulli members keep their
existing channels and refuse a transverse binding.

Original equilibrium E/G fixes the translation and section-spin interpolation,
inertia, and three geometric work ports. With original bending/shear fractions
beta and s, an E arm has coefficients `delta_E*A/L`,
`3*delta_E*I*beta^2/L^3`, `delta_E*I/L^3`; each G arm contributes
`kappa*delta_G*A*s^2/(4*L)` through the symmetric port. All original times and
layer identities survive. Neither complex frequency nor evolving memory refits
the spatial shape, and equilibrium stiffness is not counted twice. This fixed
Galerkin beam requires spatial refinement; it is not exact continuum dynamics.

`FrameStrainChannel::TransverseShear` distinguishes its physical memory storage
and heat in the existing harmonic/transient reports. Source-band intersection,
boundary states, cancellation and native balance retain their existing owners.
Preparation charges three states per declared E arm and one per G arm before
source evaluation, including declared zero increments. Four exact-source
coefficient/observation tests passed. Source and native integration regressions
cover corrected annular arms, explicit elastic channels, independent AC rigidity,
physical power and constrained shear-memory decay; those full targets await
the pinned dependencies and prescribed build runner. No live thermal update
or heat deposition is added to this fixed-source model.

## Native frame hereditary time histories

`PreparedFrameResonance::prepare_transient` and
`PreparedFrameRelaxation::prepare_transient` retain the original assembled
descriptor and prepare the existing complete `fs-phs` quadratic flow and its
integral observers for one fixed clock. Every native displacement, momentum,
original material memory and ordinary boundary state advances together under
held nodal forces. Off-diagonal consistent inertia, nonproportional viscosity,
Maxwell reaction and mutual boundary structure remain intact. No modal loss
ratio, fitted impedance or split held-velocity constitutive update replaces the
source. Algebraic boundary constraints explicitly refuse this ordinary time
path; AC/pole support is unchanged. Material time evolution requires an explicit
continuous-law support band `[0, f64::MAX]` rad/s, without extrapolating a finite
measured frequency band.

Native reports retain kinetic, equilibrium and attached storage, independently
integrated original frame and attached resistance heat, held-force work, zero
stationary-support work, a polarized storage change and original integrated
q/p residuals. Each Maxwell coordinate also retains its original element,
layer, source-arm and strain identity, elastic/viscous extension, opposing
force, storage and interval heat. Material and boundary heat are subaccounts
of attached heat, never additional energy charges. Explicit relaxed,
unrelaxed or per-coordinate initial history and complete-state restart retain
dark constrained memory. Time, motion and all memory observations publish
together after all gates; refusal leaves the accepted history unchanged.
Preparation has explicit complete dimension, byte and work bounds. Successful
warm steps allocate no memory and prepare no exponential. Reference geometry,
fixed temperature and stationary supports apply; no thermal destination is
invented by this mechanical runtime.

Five focused regression definitions cover all three analytic coupled Maxwell
decay rates with and without held forcing; dark memory and explicit startup;
full two-member consistent mass, nonproportional viscosity and mutual boundary
heat against independent physical-equation integration; refusal/retry; and
finite-frequency/algebraic boundary refusal. Three small native checks of the
unmodified production heat-form helpers, linked to the existing production
`fs_math` artifact, passed for nonproportional mass pullback, mutual Q/R cross
terms and exact source partition. These do not execute the complete prepared
flow. Full native integration execution remains outstanding in this
resource-constrained environment.

## Original element partitions of native viscoelastic heat

`ViscoForceStep.element_steps` exposes storage, mean mechanical/eigenstrain
work and true Kelvin–Voigt/Maxwell heat in original element connectivity order.
`ViscoElements` retains each element's mass and volume from the same
consistent-mass quadrature. Existing aggregate accounts are unchanged; element
receipts partition those accounts and are not added again. The native element
regression suite is authored; full Cargo execution remains outstanding.

## Translating-frame event histories

The original scheduled solid-load/support owner indexes dimensioned translating
frame histories alongside tilt histories and freezes their exact left values
when evaluating a scheduled jump. A newly introduced displacement schedule
cannot silently disappear from event timing. This does not reinterpret a vector
support displacement as a moving frame origin or change the native transient
integration law. The focused translating-frame case in `tests/transient_events3.rs`
checks retained finite schedule knots.


## Ductile-initiation gates on native inelastic continua

`ductile::DuctileHistory` retains the original immutable calibration field and
one initiation history per native material point. Its small-strain tetrahedral,
plane-stress/plane-strain and fully integrated/B-bar hexahedral entry points
propose exactly one native mechanical increment and gate every returned point.
A tetrahedral backward-Euler entry uses the same physical duration for inertia
and plastic strain rate. Full three-dimensional stress drives plane reductions;
neither element-averaged stress nor total strain replaces native point history.

Finite static/moving tetrahedra and thermoplastic tetrahedra/hexahedra use the
same owner through `ductile::finite`. Spatial Kirchhoff stress is reconstructed
as `P F^T`, with symmetry/finite checks and no repeated constitutive solve.
For B-bar hexes, the native assumed deformation is paired with its conjugate
stored Piola stress. Per-cell temperatures broadcast to their eight native
points without allocating an expanded temperature field.

The adapters read original accumulated plastic alpha for J2, Chaboche,
plane reductions and combined creep/plastic histories. Creep accumulation is
excluded. Native pure-creep thermal families explicitly refuse because their
alpha has a different meaning. Isothermal finite trait bounds and the
small-strain history trait exclude native pure-creep owners. Timed small-strain
laws must be prepared for the same physical duration supplied to the gate.

Native and ductile alpha must match exactly before every proposal; resetting
history after plastic loading or skipping increments refuses. Mechanics and
ductility publish together only after all point, source/domain and cancellation
checks pass. Initiation, native failure or invalid input leaves the borrowed
accepted state unchanged. Callers refine from that accepted pair to bracket
initiation. Reported crossing fraction refers to the held-condition PLASTIC
increment and is not a resolved event time.

`from_source_bound_j2_tet` binds every calibration to its actual source-resolved
J2 material, and `advance_source_bound_j2_tet` rechecks original card, named
state and fixed context at every subsequent increment. Changing material/context
does not silently preserve a source qualification. Generic authored numeric
entry points retain their explicitly supplied calibration assumptions.

`DuctileMaterialBinding` also admits already converged coupled tetrahedral
intervals through `validate_material_history` and
`check_material_tet_candidate`. Source-resolved laws recheck the original
card, named state and complete fixed context against both endpoint material
bundles. Explicit numerical Johnson--Cook and tabulated laws retain their
authored assumptions; there is no blanket unchecked source implementation.
The gate requires one native increment, actual endpoint temperatures and its
physical duration, reads the complete native plastic history, and performs no
second constitutive update or equilibrium solve. A thermal owner additionally
admits the immutable source reference before calling this candidate check.

Regression definitions exercise actual native uniaxial plastic response,
initiation refusal and refined retry, unloading/reversal/replay, all plane and
hex point histories, rotated finite deformation, creep separation, source/card
changes, cancellation and old-history mismatch. Static interface review checked
the original native calls and stress/deformation pairing. New Rust compilation
and test execution remain outstanding after the execution environment failed.
The gate stops at damage initiation: it supplies no stress erosion, fracture
energy, characteristic length, crack insertion, fragments or post-failure
load transfer, and does not close the broader regularized-fracture bead.

## Tensile strength on accepted material-frame states

`Frame2d::assess_tensile_strength` compares every original material member's
accepted basic forces with source-bound `AxialTensileCapacity` layer limits.
Signed axial force is retained; the actual section end moments are `-M_i` and
`M_j`. Their two assessments bound the native prismatic member's affine
axial/bending normal-stress field: compression can coexist with tensile outer
fibers, and opposite end moments cannot average away double-curvature demand.
The same observer consumes accepted static and Newmark states and retains each
layer's own modulus, geometry, strength and source conditions. It performs no
extra mechanical solve or tangent assembly.

Native Kelvin--Voigt members additionally use their actual accepted nodal
velocities and the same original reference columns as the physical damping
operator. Each layer keeps its independent axial and bending viscosity:
`sigma_mean_i = E_i epsilon + eta_ai epsilon_dot`, with outer-fiber amplitude
`r_i abs(E_i kappa + eta_bi kappa_dot)`. Signed elastic and viscous curvature
contributions combine before extrema, so unloading cancellation and a viscous
critical layer survive. Timoshenko members use their original bending fraction;
transverse shear rate is not substituted for curvature. Total material axial
force and bending moment are reported, but these aggregate resultants are never
redistributed through `E_i/EA` to recover layer stress. Their affine spatial
variation still admits the exact two-end envelope. Actual transient source-band
admission remains with native Newmark. No thermal, hereditary or Rayleigh
material stress is inferred by this observer.

Exact section/source matching, a preallocation layer budget and cancellation
checks make refusals transactional. Raw-coefficient/fiber-plastic members refuse,
since their actual stresses cannot be recovered from this source section law.
G1/G3/G4 regression definitions in
`tests/frame2d_tensile_strength.rs` cover signed compression/bending, layered and
spatially varying strength, rigid rotation, inertial double curvature, exact
budgets and foreign/unsourced refusal. Additional accepted-Newmark cases cover
distinct layer viscosities, reference orientation, signed unloading cancellation
and original Timoshenko curvature. This is a normal-stress observation;
it does not include transverse shear failure, plastic redistribution, buckling,
fracture evolution or a validated component allowable.

`PreparedFrameRelaxationTransient::assess_tensile_strength` observes the original
selected hereditary members from their accepted nodal motion and every native
Maxwell history. Equilibrium stress and each original layer/arm's signed stress
add before extrema. Its reconstruction uses the individual native channel
forces and original layer area/inertia, avoiding both aggregate-force
redistribution and loss from an underflowed diagnostic arm extension. A nonzero
canonical memory whose native force underflows refuses explicitly. Constrained
internal stresses remain observable even at zero nodal displacement, and decay
on their original relaxation times through the native transient owner.

The original Timoshenko bending fraction and two-end affine normal-stress
envelope are retained. Independent transverse-shear histories remain in the
mechanics but are not relabeled as longitudinal normal stress. Capacities must
match the original sections in selected-member order. A checked layer budget,
fallible result allocation and cancellation protect immutable observation;
work is linear in selected layers plus original memory coordinates, with no
secondary constitutive update or source lookup. Regression definitions in
`tests/frame2d_poles/relaxation/tensile.rs` cover distinct layer histories at
identical motion, signed stress relief, constrained bending stress and decay,
original shear-flexible interpolation, foreign sources and numerical refusal.
This observer covers selected Maxwell members and longitudinal tensile demand;
it does not claim shear interaction, fatigue or fracture evolution.

## Material-bound tensile strength on native solid stresses

`strength::TensileMaterialBinding` binds an explicitly selected isotropic
maximum-principal-tension criterion to the native source-derived elastic state.
Elasticity and strength must use the same immutable material card, physical
condition coordinates and quantity conventions; selecting additional strength
properties does not require an identical complete property-bundle identity.
Every assessed stress must retain the binding's actual elastic-state identity.
Uniform, per-element and per-point bindings support spatial heterogeneity and
temperature-dependent stiffness/strength without substituting a material-name
constant or a silently refreshed state.

`TetThermalStressReport::assess_tensile_strength` consumes actual constant P1
stresses, including the same reports returned by inertial trajectories.
`QuadraticThermalResponse` now retains native stiffness-quadrature stresses,
their local material identities and integration weights. Its tensile assessment
checks these samples individually: opposite bending stresses cannot cancel
through the existing element mean. Standard and B-bar samples use the same
strain as native force and energy recovery. Mandel shears are converted to
physical tensor components before the shared `fs-material` criterion.

The largest utilization and its first deterministic sample govern. Sample
budgets are checked before allocation, and each sample polls cancellation;
source mismatch, invalid counts or numerical refusal publishes no partial
assessment and changes no accepted mechanical state. The G1/G3/G4 tests in
`tests/tensile_strength.rs` cover signed loading, pure shear, zero-mean bending
hotspots, independently integrated energy, local thermal-strength changes and
transactional refusals. Results compare finite sampled stresses with a named
limit; they do not bound unsampled surface/interior stresses, supply anisotropic
or ductile failure laws, or predict crack growth, fracture energy, plastic
failure, fragmentation or service life. Finite-strain Piola stress is not
accepted as Cauchy stress by these small-strain adapters.

## Direction-dependent tensile assessment on native solid stresses

`strength::orthotropic::OrthotropicTensileMaterialBinding` binds the three
explicitly named `fs-material::strength::OrthotropicTensileLimit` sources to
the original source-derived orthotropic elastic state. Card, material
condition, complete query and typed axes must agree. The same validated
proper rotation and orientation identity construct the native
`TetElasticMaterial` and the cached stress pullback. A stress-principal frame
is never substituted for the material frame.

`TetThermalStressReport::assess_orthotropic_tensile_strength` and
`QuadraticThermalResponse::assess_orthotropic_tensile_strength` consume the
actual P1 or native stiffness-quadrature stress in its original oriented
constitutive state. Mandel stress is converted to tensor shear once and
pulled back as `Q^T sigma Q`, then the shared material owner compares each
positive normal component with that material axis's own strength. The native
response identity, oriented constitutive identity, source bundle, property
axis order and full local stress are retained in the observation identity.
Uniform, per-element and per-point bindings admit heterogeneous temperature
and orientation without averaging away tensile bending peaks. No tangent,
eigensolve or second mechanical solve is created during assessment.

The first sample at maximum directional utilization governs. Each sample
retains its governing material axis. Size budgets precede allocation and
every sample checks cancellation. Wrong source state, wrong orientation,
invalid field count, failed allocation or invalid arithmetic refuses the
whole observation while preserving the mechanical response. G1/G3/G4 tests
in `tests/orthotropic_tensile_strength.rs` exercise unequal directional
strengths, rotation, material-axis shear, compressive loading, zero-mean
Q1/Q2 bending, local thermal/orientation changes and atomic refusals.

These are the tensile branches of the explicitly selected maximum-stress
criterion, not a general composite failure law or a component safety
certificate. Compression/shear strength, interaction, delamination, damage,
fracture energy, fragmentation and stresses between retained samples remain
unassessed. Engineering ultimate strength is not converted to finite-strain
Cauchy failure strength.

## Finite-strain tensile assessment on prepared P1 solids

`HyperTetProblem<ResolvedTensileHyperelastic>::assess_tensile_strength` and the
prepared owner's same method use the original total F, Piola stress and Cauchy
push-forward from `hyper3::response`. The actual mechanics field also owns each
selected Cauchy-strength limit. Cached reference geometry is reused for repeated
observations; one result vector is allocated and no tangent or second solve is
formed. Every sample retains the unmodified native response and material/limit
identities. Only roundoff-scale Cauchy skew is admitted before comparing its
symmetric part; original stress components remain inspectable. Sample-count,
allocation, source/law-domain and cancellation refusals publish no partial field.
The manufactured cases in `tests/hyper3_tensile_strength.rs` distinguish Piola
from current-area stress, test hydrostatic extension/compression, rigid rotation,
heterogeneous materials and deterministic governing samples. These are local
stress/limit observations, without equilibrium, continuum-peak or fracture-
evolution authority. Caller-selected engineering strength is never silently
converted into a finite-strain rupture law.

## Retained source-bound brittle P1 solids

`brittle::solid::BrittleTetSolid` retains the original admitted
`TetLoadGeometry3`, each resolved brittle source/law and its native history.
Each fixed plane selects the projected footprint, and the law's band width must
match V/A within 64 relative machine epsilons. `advance` executes one native
`SmallStrainTetProblem::solve_from` endpoint with frozen accepted histories,
then checks the original constitutive trial, stress, irreversibility and work.
It integrates reference-volume storage and fracture loss while separately
reporting/gating material balance, cumulative loss, support/dead-load equilibrium
work and signed endpoint integration error. Only a complete accepted candidate
is returned; failure and cancellation preserve the original solid. Native
small-strain P1, fixed-plane and fixed-coefficient limits remain. This does not
create cracks or separated topology, infer contact, or convert numerical error
into heat. The manufactured cases in `tests/brittle_tet_solid.rs` cover softening,
unloading, complete Gc-times-area loss, compressive reclosure, independent source
strengths/widths, force/support work and rejected-candidate rollback.

## History-preserving brittle arclength continuation

`brittle::solid::continuation::BrittleTetContinuation` binds an accepted
source-defined P1 brittle solid to fixed support coordinates, a proportional
dead-force pattern, explicit SI displacement/force scales, and bounded step and
assembly work. The existing Keller bordered solver advances free normalized
coordinates and the load multiplier together, so the force multiplier may
decrease while tensile damage grows. Supports may retain a nonzero preload but
cannot move along this path. Original forces on supports remain in reactions;
duplicate nodal forces add in their input order.

Every residual and Jacobian evaluates the original material trial from the
same accepted damage history, including the native evolving-damage tangent.
A crate-private P1 evaluator retains the native reference geometry through all
corrector evaluations; no second constitutive or assembly kernel is introduced.
The native material selects unloading stiffness at its accepted envelope, so
subsequent predictors preserve the actual accepted secant direction. An explicit
physical direction can enter a damaged branch or request unloading. This
numerical direction never changes material state. Sharp constitutive corners
may need a smaller arc step or a supplied direction; no unique branch selection
or automatic event localization is claimed.

The final free-force residual is checked independently. Candidate strains,
stresses and histories then pass the same original material-work, cumulative
fracture, endpoint-work and released-face gates as ordinary `BrittleTetSolid`
equilibrium. A separate damage-increment allowance limits accepted jumps.
Refusal, exhausted assembly/step work and cancellation return no partial
history. Lowering the explicit maximum step permits a retry from the identical
checkpoint. The load-factor/deflection trace retains at most the admitted
step count plus its initial observation.

This is quasistatic, small-strain, fixed-plane equilibrium with immutable
material conditions and topology. Arc length supplies no physical time or
dynamic stability proof. Endpoint quadrature error is neither heat nor fracture
loss. Four regression groups in `tests/brittle_tet_solid/continuation.rs`
check the independent bilinear tensile response through the descending branch
and complete `Gc*A` loss, original-history unloading and replay, rejected
energy/damage/budget/cancelled proposals, and initial-force/support/scale
admission. Measured material validation and performance benchmarking remain
outside these manufactured implementation checks.

## Quasistatic normal recontact after tensile fracture

`brittle::solid::fracture::contact` continues supported, explicitly separated
P1 specimens through native quasistatic normal recontact. Existing
Newton–FGMRES solves native displacement equilibrium and `fs-contact`'s
unregularized Signorini natural map simultaneously, with one normal multiplier
per distinct full physical contact covector. It does not require an inverse
of an isolated softened fragment's stiffness. SI scales condition equations
without introducing physical contact compliance, stored energy or fracture loss.
Prepared P1 geometry, cached original tangents and reusable Krylov scratch
preserve the existing material/history owner; tangent products visit free rows.

Exact positively equivalent full rows at shared crack edges share an aggregate
force. Every original face/vertex/node/normal witness and the strictest original
opening allowance are retained, and incoming rows are budgeted before reduction.
Seeds are supplied per original non-tied pair and combined for each aggregate.
Opposite unilateral half-spaces, near duplicates and aliases that exist only
after eliminating supports are not merged. Aggregate force does not determine
individual face pressure, and no pressure allocation is inferred.

Publication replays force equilibrium directly in newtons and the original
constitutive, irreversible-history and work gates. Nonnegative normal force,
original opening allowances, physical complementarity, absolute pair work,
force/jump work consistency and contact couples must pass independently.
Separate external-work accounting removes internal contact work. Original
state and histories remain unchanged after refusal or cancellation. Explicit
budgets cover input rows, native assemblies, mixed workspaces and FGMRES history
across all configured restart cycles.

Loaded pairs must touch and remain tangentially matched at both endpoints;
open unloaded pairs carry zero force. Arbitrary closure from finite opening
requires a resolved preceding endpoint. Impact, unconstrained rigid modes,
finite sliding, collision discovery and equilibrium branch stability are
outside this fixed-trace small-strain endpoint. No position, force or history
projection is used.

`brittle_tet_fragment_contact/quasistatic.rs` specifies native analytical cases
for pressure transmission through fully softened tensile volumes, restored
compressive response, SI conditioning, common support translation and external
work, exact shared-edge aggregation, open unloaded interfaces and transactional
budget/cancellation/refusal. Native test execution is a separate obligation.

## Support-driven approach to released tensile-face contact

`BrittleTetSolid::advance_toward_released_face_contact` advances initially
open, support-driven released faces through bounded native equilibrium prefixes.
It returns the full open interval or an accepted nonnegative gap within explicit
activation and load-fraction tolerances. Every candidate replays the unchanged
original `advance_path` history, preserving native constitutive evolution,
source admission, damage subdivision and physical work. The endpoint-work
allowance scales with the actual accepted prefix fraction.

A private read-only endpoint admission checks nonnegative gaps at every
represented intermediate endpoint, including openings inside the original
negative roundoff allowance. Only finite negative released-face openings narrow
the sampled admissibility bracket. Other native material, source, strain,
equilibrium, work and cancellation refusals propagate without publishing a
partial state. No positions or constitutive histories are interpolated or
projected into contact.

The returned solid continues through the existing loaded released-face
equilibrium owner with its geometry, activation and contact-work gates intact.
This resolves the represented support-driven path; it does not certify global
continuous first contact or a monotone nonlinear branch. Initially active pairs
and unsupported free-fragment motion require their existing contact or dynamic
owners. Explicit prefix and native-path trial budgets bound repeated solves.

`brittle_tet_fragment_contact/quasistatic/approach.rs` covers analytic approach
followed by native pressure transmission, full-open-path receipt equivalence,
strict exclusion of the original negative roundoff allowance, and rollback on
work, initial-equilibrium, strain-domain, active-pair and cancellation refusals.

## Native pressure-sensitive porous compression

The existing generic `SmallStrainLaw` assembly and plastic-increment solver
accept `fs_material::plastic::porous::PorousPlasticity` directly. Its complete
volumetric plastic strain and accumulated work-conjugate strain remain in the
native per-quadrature history; no solid-specific return map or constitutive
parameter copy is introduced.

The focused `porous_compression` regression drives a P1 tetrahedron by a
hydrostatic force beyond yield, unloads to zero force, and reloads below the
previous peak. It checks equilibrium, retained permanent compaction and
elastic reloading against the material's closed-form hydrostatic solution.
This establishes the intended native interface and small-strain load path;
it does not extend the material to finite-strain foam crushing, contact or
densification.

## Temperature-dependent brittle solid equilibrium

`BrittleTetSolid::advance_at_temperature` transfers each original accepted
history through the material owner's checked temperature-successor API, then
executes the native P1 equilibrium endpoint with those new coefficients. Source
card/claim identity, nonthermal coordinates, Poisson ratio, crack plane, band
width and numerical bounds remain subject to the original material checks.
Transfer at held strain may initiate further damage; existing damage and
cumulative fracture loss cannot heal. Parameter work, transfer-induced fracture
loss and storage change are integrated separately from the subsequent mechanical
work, support work and signed endpoint integration defect. A combined exact-work
balance gates the complete result. Source mismatch, cancellation, failed Newton
iteration, loss of supported static capacity or a late energy refusal returns no
candidate and leaves the old material/history usable for retry.

The temperature path changes constitutive coefficients only. It introduces no
thermal expansion, caloric equation, heat conversion or entropy-production
claim; the caller resolves the ordered transfer/equilibrium split with physical
temperature and load increments. The added manufactured regression definitions
cover retained-damage compliance under load, damage onset at held strain,
source/overload/cancellation refusals and deterministic retries. Rust parsing and
formatting were checked; compilation and execution of these new tests remain
outstanding because the execution environment had no free disk space.


## Adaptive native continuation through progressive ductile softening

`ductile::evolution::DuctileTetEvolution` continues the ORIGINAL
`SmallStrainTetProblem::solve_from` through post-initiation softening with
`fs-material::ductile::evolution::DuctileJ2`. It reuses the native return map,
element assembly, Newton iteration and nonsymmetric FGMRES path. Each material
point retains its original immutable model, calibration and directional
characteristic length. The existing `TetCrackBandGeometry` supplies
h = V / projected area. Admission checks original connectivity and volume
exactly, derived area/width within normalization roundoff, and exact agreement
between the material length and the supplied band record. Agreement of those
geometric measures does not authorize history transfer to a changed mesh.

One call advances a finite positive physical interval with affine prescribed
displacements and dead forces, using one actual uniform or per-element held
temperature field. Rejected Newton/constitutive trials and excessive endpoint
numerical loss halve the proposed interval. Material damage increments retain
their own explicit refusal threshold. Accepted steps may grow up to twice the
previous fraction using measured damage and numerical-loss margins. Explicit
minimum fraction, maximum accepted steps and maximum total trials bound the
work; each retry starts from the same last accepted pair and consumes no
physical time until accepted. Load reversals must be separate intervals.

The adapter independently checks whole-solid external work from actual native
support reactions, applied forces and displacement increments. Original element
volumes integrate stored-energy change, plastic dissipation, damage release and
backward-Euler numerical loss. Original projected areas integrate the calibrated
equivalent-plastic-work diagnostic; that diagnostic is not charged again as
physical loss. Compensated totals and finite/sign gates preserve separate
accounts. The declared numerical-loss and external-work tolerances are checked
per accepted fraction and again over the whole interval. Numerical loss does
not become heat or fracture energy, and its budget is not a trajectory-error
certificate.

No constitutive solve is repeated to audit an equilibrated endpoint. Only the
final native full-field state and bounded scalar records are retained. Original
material/source cards are borrowed, not cloned per point or trial. An error,
budget exhaustion, minimum-interval refusal or cancellation returns no partially
advanced interval and leaves the caller's accepted displacement and histories
unchanged. Exact held-strain material replay prevents the native solver's
initial assembly from inventing plastic flow from yield roundoff.
The accepted displacement must reproduce each retained material strain using
the exact native assembly operator. This is checked once before subdivision,
without a constitutive trial. Inconsistent displacement/history pairs refuse
even when zero stress at the requested endpoint would hide their missing
initial energy from the external-work check. The dynamic adapter shares this
admission, and both reuse the admitted material-law buffer across retries.

Four regression definitions in `tests/ductile_evolution.rs` exercise actual
native tetrahedral assembly: post-initiation reaction loss with bounded damage,
an independently calculated shear/work oracle, whole-interval rollback and
deterministic restart, ordinary non-axis-aligned normal admission, foreign
geometry/width refusal, and heterogeneous widths with distinct histories and
projected-area work. Native API, lifetime and oracle review passed. These four
structural tests and the five motion/restart regressions below now pass in a
bounded actual-source integration build. It executes the original material,
P1 geometry and assembly, consistent reference mass, native Newton/Krylov
solvers, allocation arena, execution context and cancellation gates. Unrelated
module exports and upstream timestamp serialization are omitted from that
compilation root; no numerical or cancellation substitutes are used. All nine
regressions execute successfully. This is integrated structural execution,
not the full-workspace Cargo or DSR gate.

This is quasistatic reference-geometry small-strain continuation with fixed
mechanical coefficients. It adds no finite-deformation degradation, automatic
crack direction, crack insertion, geometric contact, erosion or fragment state.
Directional length scaling and adaptive damage steps do not establish arbitrary-
mesh fracture objectivity or experimentally calibrated component failure.
MR63's broader topology/contact/post-failure-transfer requirements remain open.

### Momentum continuation through ductile softening

`ductile::evolution::dynamics::DuctileTetMotion` binds the original admitted
ductile tetrahedral specimen to one immutable consistent reference-mass
inventory. `advance_from` continues its actual displacement, velocity and
plastic/damage histories through initiation and progressive softening using
the native backward-Euler momentum solve. Reference density is admitted once;
damage and temperature do not silently alter mass. Affine endpoint loads and
prescribed motion use each actual accepted physical subinterval, and that same
duration reaches the material law. Rejected trials retain the last accepted
motion and history; whole-interval refusal leaves all caller inputs unchanged.
Initial prescribed displacements must exactly match the declared starting
load factor, so a different support origin cannot introduce an accidental jump.

Complete external work includes native inertial support reactions. Independent
receipts retain material storage, plastic dissipation, damage release and
material numerical loss, alongside kinetic change and the consistent-mass
`1/2 delta_v M delta_v` numerical loss. Inertial work is removed from the material
work account exactly once. Material and kinetic numerical losses have separate
whole-interval budgets allocated across accepted fractions; excessive loss
refines the attempted interval. Compensated sums and per-fraction and final
balance checks retain these distinctions. Neither numerical loss is physical
heat or fracture work. Reused load/law buffers and the retained mass avoid
rebuilding that inventory for each retry.

The two source-bound tetrahedral regressions in `tests/ductile_source.rs` also
pass through original J2/history/source-binding code, including the prepared
Lode surface and transactional source-swap refusal. With the 20 core and 19
material-source tests, these 11 structural checks give 50 focused passes.

Regression definitions in `tests/ductile_evolution/dynamics.rs` cover an exactly
accelerated free tetrahedron and momentum-preserving kinetic-loss refinement,
an independently forced free shear DOF that solves consistent inertia and active
ductile softening together, imposed shear through damage initiation with
independently calculated stress and work under applied support forces, and
bounded/cancelled continuation with
exact deterministic replay. A shared restart regression refuses virgin history
at a deformed displacement before unloading and replays the correct unchanged
native state in both static and dynamic continuation. All five regressions
pass in the actual-source build described above. This remains small-strain/rotation,
held-temperature continuation with fixed reference geometry and directional
length; it adds no finite-deformation damage, caloric feedback, crack topology,
contact or fragment transfer. Loss budgets do not certify trajectory accuracy.

## Finite-deformation ductile tetrahedra

`hyper3::plastic::ductile::FiniteDuctileTetSolid` carries the original
multiplicative finite-J2 metric, accumulated plastic strain and shared ductile
history through progressive softening. The native reference tetrahedron,
deformation-gradient operator and dead/body-force assembly remain authoritative.
The original directional `TetCrackBandGeometry` fixes the reference length;
admission checks connectivity, volume, projected area and the material width.
Prepared material owners are borrowed, and reference elements are shared across
accepted states. Original source-bound finite laws retain their explicit
Kirchhoff-plastic-work-per-reference-area meaning and source receipts.

The damaged first-Piola response and its complete analytic 9-by-9 derivative
enter the existing nonpotential increment operator, Newton iteration and
nonsymmetric FGMRES solve. Geometry preparation uses the original finite
reference owner; no scalar potential or symmetric tangent is invented for the
damage response. Every trial starts from the same accepted finite plastic and
damage histories. Stress-only backtracking avoids full derivative assembly.
One physical duration and actual uniform or per-element held temperature reach
both native flow and initiation. Material, work, geometry or cancellation refusal
publishes none of the proposed displacement or histories.

Native forces and support reactions independently establish external work;
reference volumes and endpoint Piola stress establish internal material work.
Compensated accounts retain storage, nonnegative plastic dissipation and damage
release separately from signed finite-step work defects. The explicit defect
budget bounds the SUM of absolute element defects, so opposite rotation-induced
errors cannot cancel. Constitutive work consistency, assembled work, equilibrium
work and physical balance receive independent gates. Numerical work defects
never become physical heat or fracture energy. Exact held-deformation replay
retains the accepted history without repeating plastic flow or physical time.

All four focused regressions in `tests/finite_ductile.rs` pass in the bounded
actual-source build: independent affine reaction decay and work, a nonlinear
free-coordinate equilibrium through active damage, restart, proper finite
rotations with opposing signed defects, and geometric/damage/cancellation
refusals. The build executes the original material, reference geometry,
Newton/FGMRES, allocation and cancellation owners; unrelated modules and unused
hyperelastic generic defaults are omitted. This is not the full workspace
Cargo or constellation gate. The quasistatic endpoint advances one physical
interval; the bounded path API below automates subdivision. Evolving caloric
feedback, crack
topology, erosion and fragment/contact transfer remain outside this owner.
MR63's post-failure topology requirements remain open.

### Consistent momentum through finite ductile softening

`hyper3::plastic::ductile::dynamics::FiniteDuctileTetMotion` retains the complete
finite ductile specimen, physical nodal velocity and one original consistent
`P1ReferenceMass`. The supplied positive reference densities declare the mass
inventory without claiming material-card calibration. Mass and reference
density are admitted once and shared across accepted endpoints; neither damage
nor current volume changes the inventory. Free bodies require no artificial
supports. Nonzero initial velocity and moving prescribed supports are retained.

`try_new_with_nodal_density` admits the original positive P1 density field and
reuses the exact cubic reference-mass integration. Complete nodal values and
the prepared mass are shared across accepted endpoints; cell-average densities
are reporting values only. `reference_acceleration_load_n` uses that same full
mass to project a prescribed P1 acceleration into held nodal forces, including
supported nodes and the graded mass centroid. Neither the deformation Jacobian
nor the damage/integrity multiplier changes this original inventory. Current
density obeys rho_0/J in this closed body; `body_force_n_m3` remains force per
original volume. Position-dependent field evolution and potential work belong
to their external owner; no mass transfer or fragment removal is introduced.

The original `BackwardEulerInertia` supplies momentum forces, the mass tangent,
endpoint velocity and kinetic work. Material and momentum share the SAME
nonpotential Newton/FGMRES solve. The inertia residual and derivative enter
exactly once. One `DuctileStepConditions.duration_s` drives both momentum and
rate-dependent initiation; no independent constitutive clock is introduced.
All material, geometry, directional-band and signed-defect gates from the
quasistatic owner remain active.

Native support reactions include inertia. The original endpoint-work owner
removes inertial work from the material account exactly once; independent
checks retain storage, plastic dissipation, damage release, signed material
work defect, kinetic change and `1/2 delta_v M delta_v` numerical loss. Material
absolute-defect and kinetic-loss budgets are separate. Both component work
identities and their combined balance must pass before displacement, velocity
and histories publish together. Numerical losses never become heating or
fracture energy. Refusal leaves the original complete motion available for
deterministic retry or physical-interval refinement.

The four regressions in `tests/finite_ductile/dynamics.rs` cover nonzero-velocity
free acceleration and retained mass, an independently forced free coordinate
with simultaneous damage and inertia at two physical durations, moving-support
work, independent material/kinetic loss refusals, replay, refined retry and
cancellation rollback. The added density cases compare constant and graded
inventories with equal total mass, exact original nodal gravity shares and
cubic inertia/reactions through finite deformation and damage. The original
regressions passed their earlier actual-source build; the added graded-density
cases have not been executed in the disk-exhausted implementation environment.
This is one backward-Euler interval with held temperature and original reference inventory. It adds no
impulses, evolving caloric feedback or fragment
transfer, and the loss budgets do not certify trajectory accuracy.

### Adaptive finite ductile physical intervals

`hyper3::plastic::ductile::path` adds `advance_path` to the original finite
ductile solid and motion owners. One shared controller advances their existing
endpoint solves; it adds no constitutive, equilibrium or inertia solver.
Supports interpolate from the actual accepted coordinates to explicit final
coordinates. Dead nodal and referential body loads interpolate between explicit
initial and final values. Quasistatic initial free equations must satisfy the
original force/load assembly; moving specimens retain incoming velocity and
may begin with physical acceleration. Held uniform or per-element temperature
and the original geometry, material, directional bands and reference mass
remain fixed through the interval.

Every attempted fraction supplies its actual physical duration to both native
momentum and rate-sensitive initiation. Only typed Newton/linear failures,
material/controller damage-increment limits, summed absolute material-work
defects and numerical kinetic-loss excess trigger subdivision. Material
source/domain/history refusals and cancellation stay terminal. Failed material
evaluations retain their typed cause across the native Newton boundary;
successful evaluations cannot inherit a rejected backtracking trial's cause.
Rejected endpoints consume neither physical time nor accepted history.

Material absolute-defect, kinetic-loss and equilibrium-work allowances are
whole-interval budgets apportioned by actual subinterval size. Bounded growth
uses accepted damage and independent loss margins. Explicit minimum/maximum
fractions, total trials and accepted-record counts bound work and storage.
Compensated whole-interval accounts separately check material storage/loss,
external/support work, signed defects, absolute defects and kinetic work.
Neither signed cancellation nor another loss account can hide an exceeded
material or kinetic budget. A failed whole request publishes no private prefix.

Results retain the last complete native endpoint plus compact accepted records
and whole-interval accounts. Original endpoint results also retain scalar
`element_energy` receipts from their existing final material audit, allowing
coupled consumers to read physical losses without repeating a return map or
differencing cumulative histories. A path endpoint's per-element receipts
describe only its LAST substep. Load reversals require separate intervals;
this controller provides no arc-length method, trajectory-error certificate,
thermal evolution, crack topology or fragment transfer. MR63 remains open.

All four new `tests/finite_ductile/path.rs` regressions and the eight original
finite static/motion regressions pass in focused actual-source verification.
Independent references cover axisymmetric stress/work through softening,
accepted-partition replay, nonzero-velocity free acceleration and kinetic-loss
refinement, and a genuinely free rate-sensitive J2/Johnson--Cook coordinate
with consistent inertia at two physical durations. Initial-force/domain
refusal, failure after a privately accepted prefix, distinct loss-budget
exhaustion and cancellation preserve original state. This executes the original
material, structural, inertia, Newton/FGMRES and cancellation owners; it does
not establish the full-workspace Cargo or constellation gate.

## Adaptive brittle load intervals

`BrittleTetSolid::advance_path` reuses the original native endpoint solve and
retained histories across bounded affine physical load intervals. Initial dead
forces are explicit and must satisfy the original assembled free-force residual
before continuation. Subdivision limits the largest element damage increment
and allocates the whole-interval endpoint-work-error budget in proportion to
attempted interval size. The sum of absolute local defects is checked, so signed
cancellation cannot conceal numerical work error. Native Newton/linear failures
and these refinable endpoint gates may subdivide; source/domain/admission errors,
other work gates and cancellation are terminal. Trials and compact accepted
records are bounded, aggregate work uses compensated sums, and any failed whole
request leaves the original solid unchanged even after private partial progress.

The native solution describes the last substep; accepted records retain each
substep's work while the result reports complete-interval storage, material work,
fracture loss and external work. This is fixed-condition quasistatic load
continuation: reversals require separate intervals, and force-controlled limit
points still require another physical loading/control method. It is not an
arc-length solver, a physical time integrator or a crack-path convergence claim.
The focused regression definitions in `tests/brittle_tet_solid/path.rs` cover
complete Gc-times-area failure, retained unloading/reclosure, force histories and
transactional refusal. Rust parsing/formatting and independent bilinear arithmetic
checks passed; native compiled execution remains outstanding.

## Finite porous compression through the native tetrahedral solver

`fs_material::plastic::finite::porous::FrozenPorousTrial` is consumed by
the existing generic `HyperTetProblem` and its per-element
`HyperelasticLawField`. Its accepted material metric remains borrowed and
immutable throughout assembly, Newton and line searches. The caller checks
all returned material candidates before committing any history. The solver's
reported potential is the reduced incremental potential; physical stored
energy and plastic heat remain separate in the native material response.

The focused `hyper_porous` regression drives all six deformational degrees
of freedom of a P1 tetrahedron with dead referential force resultants.
It compresses to stretch 0.8 in all directions, unloads to retain permanent
volume loss, and reloads elastically. Independent hydrostatic formulas check
deformation, Piola/Cauchy stress, current/plastic volume, storage and physical
heat. This uses the original geometry, constitutive adapter and nonlinear
solver. Retained body stepping is supplied by `FinitePorousTetSolid` below;
the adapter adds no contact, remeshing, finite-porosity densification or
locking-free mixed element.

## Elastic prehistory and source-bound brittle initialization

`brittle::solid::initiation::TensileElasticTetSolid` retains original per-element
source bundles and advances their borrowed elastic backbones through the native
P1 solver. Its accepted elastic history keeps a fixed prescribed DOF set. Every
proposed affine endpoint is checked with the source-bound native tensile
detector; a crossing returns a disposable elastic predictor instead of accepting
an over-strength elastic state. Repeated principal tension requires an explicit
physical hint. Per-element events, numerical first-crossing brackets and the
original predictor remain inspectable, with bounded native onset work and Cx
polling.

`resolve_bands` derives each actual V/A width on the retained geometry. Predicted
crossings select their principal planes; elements without crossings retain the
explicit original source plane. `handoff` requires original sources, coefficients,
strain/work bounds and matching geometry, then seeds every receiving native law
at the OLD accepted elastic strain. Local stress/storage continuity and an
explicit whole-solid energy gate prevent a constitutive change from creating
work, loss or healing. Only the complete result is published; it can replay the
requested load with the existing single-endpoint or adaptive brittle solver.

This is an elastic-predictor, fixed-plane approximation. Later predicted planes
and crossings do not incorporate stress redistribution caused by earlier
fracture; those brackets are not accepted event times or a mesh-objective crack
path. Resolved mixed-mode damage at the seed refuses rather than being erased.
The new regression definitions cover physical Y-plane selection after elastic
preload, native softening after handoff, repeated-stress hints, heterogeneous
sources and transactional refusal/retry. New code parsed/formatted and passed
bounded API/analytical review; compilation and runtime execution remain pending.


## Retained source-bound brittle physical-time motion

`brittle::solid::dynamics::BrittleTetMotion` starts from the complete accepted `BrittleTetSolid`, including damaged histories and stored/fracture energy, plus actual nodal velocity and a bounded physical clock. It requires positive density from each original resolved material bundle, binds that snapshot to the supplied reference geometry, and assembles the native consistent P1 reference mass once. Every interval retains that exact mass inventory, source identity, fixed crack planes and irreversible history; damage cannot erase matter or momentum. No conversion from current/hot density to a different reference volume is inferred.

Each physical endpoint uses the existing P1 backward-Euler Newton owner and its original support/inertial coupling through a crate-private immutable-mass entry point. Static and dynamic brittle paths share one native constitutive-work recovery. Loads and support displacements are actual endpoint values; free-body motion requires no artificial supports. Above-peak tensile loads may evolve dynamically when the inertial tangent admits the endpoint; this is not a claim that arbitrary softening steps converge.

The returned receipt separates exact constitutive path work, recoverable storage, irreversible fracture dissipation, applied/support endpoint work, consistent-mass kinetic change, constitutive endpoint quadrature error, and backward-Euler numerical kinetic loss. Separate absolute gates bound equilibrium, physical/cumulative fracture work, retained kinetic storage, the native kinetic identity and both numerical accounts. Numerical damping is neither fracture work nor generated heat. A rejected solve, source, time, budget, cancellation or late energy gate preserves the entire old owner and physical clock.

Four native regression definitions cover original density/mass scaling and free momentum, an independently calculated above-static-peak tensile-softening endpoint with retained damage, moving-support work with concurrent applied forces, and transactional refusal/restart. New Rust was parsed/formatted with the pinned nightly formatter; compile/runtime execution is not claimed in the full-filesystem editing environment. This remains a fixed-plane small-strain, first-order physical-time slice; it does not create crack topology, fragment contact, evolving reference inventory or a time-accuracy certificate. Physical timestep refinement remains required to control the separately visible integration losses.


## Native physical-load progress before first tensile onset

`brittle::solid::initiation::first_event` adds `TensileTetOnset::advance_before_first_event`. The first lower onset bracket is a proposed load fraction, never an accepted finite-element state. The method interpolates original accepted applied forces and support coordinates toward the original requested target, then re-solves and re-probes through the same native elastic and source-bound tensile owners. Only a candidate with no tensile crossing is accepted. A bounded replay count, minimum physical load fraction, nonzero displacement progress, cancellation and native domain admission prevent unbounded or stationary continuation; numerical solve/onset-resolution refusals may back off, while material/input failures remain terminal.

The accepted elastic owner now retains original applied-force entries, including loads on prescribed coordinates. The result keeps the exact original target entries and ordering, duplicates, original solver/onset controls and physical direction hints. It also returns canonical remaining per-DOF force intervals from the actual accepted load to the original target. Support reactions are never substituted for applied loads. `probe_remainder` uses the retained original inputs, and reversible elastic work is reported as the original storage change without inventing endpoint work or heat. Pure displacement loading avoids dense force scratch allocation.

Four native regression definitions cover elastic preload/work continuity followed by brittle replay, duplicate force entries and independent applied support loads, bounded/cancelled/zero-progress refusal, and a retained physical hint for a degenerate principal spectrum. The complete new and modified Rust files parsed/formatted with the pinned formatter; native compile/runtime execution remains unperformed in the full-filesystem environment. This is quasistatic progress toward the first onset only. Later elastic-predictor crossings do not establish event order after damage redistributes stress; later plane selection, rotating cracks and fragment topology/contact remain with their own physical models.


## Temperature-dependent tensile failure in retained physical motion

`brittle::solid::dynamics::temperature::BrittleTetMotion::advance_at_temperature` performs the original checked temperature/history transfer at held strain, then advances the same specimen through the native physical-time inertial solver. The shared internal preparation preserves the static temperature path and exposes no public partially transferred state. Original card/claim identity, held nonthermal coordinates, plane/width, Poisson projection, constitutive bounds and irreversible history remain governed by the existing source owner.

The parameter stage retains the exact cached reference mass and original reference-density snapshot, velocity, kinetic storage, clock and step budget. A different density evaluated at a successor temperature does not change this closed reference inventory or supply an impulse. The result separates signed coefficient work, temperature-induced fracture, subsequent mechanical fracture, elastic storage, kinetic change and both numerical accounts. Combined full-work and cumulative-fracture gates use the native stable kinetic increment. Every source, solver, time/budget, cancellation or late energy refusal preserves the old motion and material state.

Four native regression definitions cover heating-induced tensile failure during rigid translation, frozen mass/force acceleration under temperature-varying density, cooling without healing, and late refusal/retry. The new Rust parsed/formatted with pinned Rustfmt, and independent coupon arithmetic supports the work signs; no compile or runtime execution is claimed in the full-filesystem environment. This is a first-order coefficient/history split followed by motion. It does not infer thermal eigenstrain, caloric heat, fracture heating, entropy production or an unresolved temperature trajectory.


## Adaptive physical time for retained tensile failure

`brittle::solid::dynamics::time::BrittleTetMotion::advance_time` completes an affine physical load duration using one native coarse step versus two native half steps from the same accepted material/history owner. It accepts the unmodified two-half-step trajectory; no Richardson-extrapolated displacement, velocity or irreversible state enters the material solver. Consistent reference mass stays cached. Compact records retain physical times, local estimates, damage increments and the two actual energy receipts without retaining complete intermediate fields.

Local indicators compare Euclidean displacement/velocity differences with increment-based relative scales, per-element damage, and the sum of absolute element differences in original incremental fracture dissipation. Separate gates limit either half-step's damage increase, whole-duration absolute constitutive quadrature defects, and whole-duration backward-Euler numerical kinetic loss. Applied forces and prescribed support coordinates follow the original affine physical interval, including actual endpoint values and native duplicate-force semantics. Native Newton/linear failures and the named finite numerical-work gates can refine; source/domain, geometry, physical-energy and cancellation failures remain terminal.

Explicit comparison, native-solve, retained-step and minimum-interval budgets bound work and memory. Coarse and rejected trials consume solve work but never accepted motion steps. Every accepted comparison advances two native physical steps; refusal of an incomplete duration exposes none of its internal partial progress. Compensated whole-duration work sums reuse the existing brittle path helper, and kinetic change sums the native stable increments instead of subtracting large kinetic stores. Full material, cumulative-fracture, equilibrium, kinetic and combined energy balances remain independently gated.

`advance_time_at_temperature` composes this trajectory with one original source/history transfer. Fixed-step and adaptive thermal motion share preparation and combined work accounting; every adaptive replay retains the same successor temperature and original reference mass. The mechanical substep damage cap excludes the explicit initial parameter-induced jump. A refused duration publishes neither intermediate motion nor the temperature change.

Four adaptive-motion regression definitions cover independent constant-force temporal convergence, exact native half-step history replay through tensile softening, moving-support/affine-force work, and bounded transactional refusal. A fifth integrated regression covers heated adaptive motion with temperature-varying queried density and late duration refusal. Pinned Rustfmt parsing/formatting and independent free-body/softening arithmetic were performed; native compile/runtime execution remains unperformed in the full-filesystem environment. Step doubling remains a local error indicator for this first-order fixed-plane small-strain model, not a global error certificate, unseen-reversal guarantee, crack-path convergence proof or adaptive thermal-path integration.


## Adaptive physical time with work-conjugate tensile motion

The midpoint owner now exposes `advance_midpoint_time` with explicit local motion, damage and fracture-work tolerances. A private scheme policy reuses the original adaptive backward-Euler subdivision, bounded solve/history budgets, original-state replays and compensated work summation. Existing backward-Euler controls and records retain their default public types and original endpoint-force/work semantics. Midpoint records retain genuine mean-work receipts and zero numerical kinetic loss.

One coarse native midpoint solve is compared with two native half steps from the same accepted history. Only the unmodified half-step trajectory can be accepted; irreversible damage is never extrapolated. Affine applied forces are sampled at the actual mean of each coarse or fine interval, while prescribed supports reach their actual absolute endpoints. Local indicators use increment-based displacement/velocity scales, per-element damage and noncancelling element fracture-work differences. Either accepted half step must satisfy the damage-increment limit. Only numerical Newton/linear refusals and excessive local indicators refine; native material/domain, physical-work and cancellation refusals remain terminal.

Whole-duration material, cumulative-fracture, mean-stress work, equilibrium, kinetic and total energy gates supplement the native per-step gates. Absolute accumulated mean-work defects cannot cancel over time. Tangent-fallback telemetry counts every completed coarse/fine solve, including discarded comparisons; partial failed-solve telemetry is unavailable. A refused duration publishes no intermediate clock, velocity or history.

Four focused regression definitions cover independent elastic-oscillator phase convergence without kinetic damping, affine-force mean sampling with its analytic position error and support work, exact accepted-half-step replay through tensile softening, and transactional work/minimum-interval/cancellation/native-refusal behavior. Changed Rust parsed/formatted with the pinned formatter and received bounded source/API review; compilation and native runtime tests were not run in the full-filesystem environment. Raw coarse/fine differences remain conservative local indicators, not a smooth-order guarantee across fracture switches or a global time-error certificate. Non-affine load/reversal knots and thermal-path resolution remain explicit caller responsibilities.


## Sourced thermal expansion drives native solid tensile failure

`brittle::solid::expansion::ThermalBrittleTetSolid` retains the original source-integrated expansion reference, native brittle histories, directional geometric band widths and P1 reference mesh. Virgin construction requires zero free thermal strain. Every successor must retain original source claims, reference, held nonthermal coordinates, projection and bounds. The returned native FE strain is total displacement strain; each original material state's strain is mechanical strain.

An endpoint first performs the original coefficient/history transfer at held mechanical strain, then uses the existing per-element small-strain Newton owner to solve total displacement and thermal eigenstrain together. All disposable proposals start from the same prepared mechanical history. This allows a free body to expand without accepting a fictitious fully restrained damage state; restraint instead creates actual mechanical tensile loading and irreversible fracture. No alternate finite-element assembler, nonlinear solver or fracture law is introduced.

Recovery replays the converged mechanical strain through the original exact path-mean stress and endpoint state. Whole-mesh accounts separate coefficient work, total-strain mechanical work, eigenstrain work, storage, coefficient-induced fracture, strain-induced fracture and cumulative physical loss. Mechanical plus eigenstrain work must equal original constitutive strain work. Support/dead-load endpoint work is checked against endpoint internal work, and its signed quadrature defect is reported separately from physical fracture. Transfer, work-split, material, cumulative-loss, equilibrium and full-balance gates pass before publication. Source, native-solve, cancellation and late energy refusals leave the original owner retryable.

Four native integration regression definitions cover freely expanding/cooling geometry across an alpha(T) knot, restrained cooling with independently integrated mode-I softening plus transverse elastic storage, reheating/unloading without healing, and source/reference/initial-state/native-solve/cancellation/late-work rollback. Changed Rust parsed/formatted with pinned Rustfmt and received bounded source/API/physics review; compile and native runtime execution remain unperformed in the full-filesystem environment.

This owner is quasistatic, small-strain and fixed-plane. Coefficient changes precede the coupled strain solve; thermal and mechanical reversals require resolved physical endpoints. A large free-cooling jump can put the initial Newton proposal on a softening branch even when a stress-free equilibrium exists, so the caller must resolve thermal increments. Thermal source increments remain explicit; physical motion is provided by the thermal dynamic owner below. Complete failed-face topology release is admitted separately. Automatic thermal substeps, partially loaded cuts, finite rotations and caloric quantities remain outside this owner.


## Retained physical motion with sourced thermal tensile fracture

`brittle::solid::expansion::dynamics::ThermalBrittleTetMotion` starts from an accepted original thermal solid and explicit nodal velocity. It binds positive density from those same selected source bundles to the original reference mesh once, assembles the existing consistent P1 mass, and retains that inventory across every temperature, deformation and fracture endpoint. A changed current-density query supplies neither new mass nor an impulse. Binding a chosen initial snapshot to reference geometry is explicit; the expansion reference alone does not infer a current/reference density conversion.

Static, backward-Euler and midpoint thermal endpoints now share one private coefficient/history preparation and native material-work recovery. Source changes first preserve the held mechanical history under the original temperature-transfer law. Displacement, thermal strain and inertia then evolve together through the existing P1 Newton solver; disposable proposals cannot accept an intermediate artificially restrained damage state. Recovery verifies the original endpoint history, validates the actual endpoint or mean solver response, and restores physical endpoint stress before exposing the accepted solution. Native material and source owners, fixed planes, geometric band widths, original budgets and cancellation remain authoritative.

`advance` uses the existing backward-Euler inertial owner. Its complete balance includes external endpoint work, signed coefficient work and eigenstrain work against elastic storage, irreversible fracture, native kinetic change, explicit numerical kinetic loss and the total-displacement endpoint quadrature defect. The last two are numerical accounts, never fracture heat or physical damping. `advance_midpoint` uses the existing centered inertia and original analytic path-mean material stress. Loads are interval-mean forces; support coordinates are endpoint displacements. The complete mean-work balance has no endpoint compensation or numerical kinetic loss. The original optional finite-difference tangent fallback differentiates this same mean response only when its classical derivative is unavailable, with bounded native telemetry.

Both schemes retain original source references, mechanical histories, mass, velocity, kinetic storage, physical clock and accepted-step budget until all material/work-split/cumulative-fracture, force equilibrium, kinetic-continuity and complete energy gates pass. Any source/domain, solver, time/budget, cancellation or late work refusal preserves the complete old motion for retry. Shared recovery integrates actual solver stress before endpoint restoration, so response work and physical endpoint stress are never silently interchanged.

Eight native regression definitions cover cooling-induced tensile softening with a genuinely free inertial degree of freedom and independent endpoint/mean-force analytical solutions; complete transverse elastic storage and irreversible loss; reference-mass retention under temperature-dependent density; matched thermal expansion with explicit support inertia; independent discrete oscillator phase and energy; and damaged-history/native-tangent/source/time/budget/cancellation/refusal replay. Four material regressions support the same-source density and mean-force adapter. Changed Rust parsed/formatted with pinned Rustfmt, bounded source/API/physics reviews completed, and independent scalar arithmetic checked the physical oracles. Native compilation and runtime tests remain unperformed in the full-filesystem environment.

The thermal path remains prescribed and explicitly resolved. Coefficients change first, then total/eigenstrain follow one straight increment under new coefficients; centered momentum does not upgrade that thermal split to a second-order coupled temperature integrator. Small-strain fixed-plane fracture, immutable reference inventory and caller-resolved thermal/mechanical reversals remain the admitted domain. Completed tensile cuts and fixed-configuration contact events are supplied by the thermal fracture continuation below. This owner does not infer heat capacity or entropy, convert fracture/numerical losses to heat, or supply automatic thermal-path adaptivity.


## Completed thermal tensile fractures retain their original bulk history

`ThermalBrittleTetSolid::split_fracture` releases explicitly selected conforming faces through the same native topology and failed-plane admission as `BrittleTetSolid`. Both incident bands must have completely failed on the selected normal and carry no unresolved traction or compressive contact. Reference coordinates, local element order, volumes, projected areas and V/A regularization must remain exactly unchanged.

The shared handoff keeps native MECHANICAL strain/history distinct from TOTAL finite-element displacement strain. Thermal source owners, integrated free-strain reference, physical endpoint stress, bulk storage and cumulative irreversible fracture survive the cut exactly. Existing fracture loss remains attached to each original volume; the operation performs no new surface work, adds no cohesive law and resets no damaged material. Old reaction records are cleared because global node numbering changes. Previous released faces are remapped through their original element corners.

All thermal static, backward-Euler and midpoint endpoints retain the same paired-face opening guard through common material recovery. A closing endpoint refuses before publication, including during a source-temperature change; the accepted owner remains available for retry. This is a small-strain fixed-reference release of fully completed cracks, without arbitrary collision detection, partial-damage insertion or finite sliding.

Native regression definitions cover restrained cooling to complete tensile failure while total displacement remains zero, retention of transverse storage and irreversible history through release, further source changes with matching total expansion, incomplete/compressed-face refusal and transactional closure refusal. Rustfmt parsing and bounded source/API/physical review completed; native compilation and runtime tests remain unperformed in the full-filesystem environment.


## Thermal fracture continues as native fragments and contacting bodies

`ThermalBrittleTetMotion::split_fracture` composes the original thermal bulk-history release with the same physical motion transfer as the nonthermal owner. Actual velocity is copied through the admitted topology map, and the retained original density snapshot reconstructs consistent P1 inertia on the split nodal stars. Independent mass, first-moment, linear/angular-momentum and kinetic inventories are retained alongside whole-body and per-fragment results. Native kinetic storage must agree before and after the representation change. The physical clock, reference inventory and used/remaining step allowances do not change; a cut is neither an impulse nor a timestep.

`ThermalBrittleTetMotion::impact_released_faces` uses the identical bounded simultaneous normal-impact kernel as `BrittleTetMotion`. Every active non-tied paired vertex participates in the coupled consistent-mass mobility solve, including cross terms and initially separating active rows. One combined native impulse updates physical velocity. Actual relative speeds, restitution/complementarity, passivity, linear/angular impulse moments, absolute contact couples and whole-body work must pass the existing physical gates. Held velocity components retain their current speeds and expose their external support work.

An accepted event preserves the original thermal source bindings, integrated eigenstrain reference, native mechanical state, total displacement, stored energy, cumulative tensile-fracture loss, mass, time and step budget. The material owner never passes through a substitute nonthermal constitutive state. Penetration, tangentially separated active witnesses, unresolved contact rank/budgets or late work failures refuse transactionally. Impact dissipation remains a contact account, without inferred heat or a fracture-history reset.

These APIs provide explicit completed-cut and fixed-configuration contact events. The thermal approach and sustained owners below supply contact localization and normal-force continuation. Finite sliding and automatically resolved temperature paths remain outside this extension. Regression definitions exercise retained thermal matter/history and motion across cuts, source-changing continuation, opening guards and simultaneous contact. Changed Rust parsed/formatted and received bounded cross-owner review; native compilation and runtime execution remain unperformed in the full-filesystem environment.


## Thermal-fragment contact localization at a held source state

`ThermalBrittleTetMotion::advance_toward_released_face_contact` uses the same bounded released-face locator as the original nonthermal motion. Every trial replays the accepted thermal body's CURRENT source snapshots, integrated free strain, native mechanical history and frozen reference mass through its own midpoint owner. Support coordinates follow the accepted-to-target affine path, and each shortened trial uses the actual prefix mean of the supplied affine forces. No intermediate temperature, source data or damaged material state is interpolated.

The result retains the original thermal midpoint step and work receipt, an optional approaching-face witness, the native trial count and successful-trial tangent telemetry. Initial active approach consumes no physical step. Only an accepted prefix advances the physical clock; discarded localization trials leave the predecessor unchanged and spend the explicit trial allowance.

Only explicit native released-face closure or a strict negative paired gap narrows the unsafe endpoint. Every non-tied pair is inspected before returning a witness, so a later penetration cannot hide behind an earlier active face. A nonnegative full endpoint returns directly; a bracketed result must meet both the activation-distance and endpoint-time allowances. Source, solver, work, budget and cancellation errors remain terminal and transactional. No penetrating candidate, projected gap or interpolated history is published.

This is bounded localization of the native affine displacement trace and sampled endpoints, not a global earliest-collision certificate for a nonmonotone continuous trajectory. Temperature changes and load reversals require separately resolved native intervals. Closed contacts driven inward at rest can require sustained forces and refuse this impact-approach entry.

Focused thermal regression definitions retain actual cooling-induced tensile failure, transverse prestress and the original density inventory while checking contact timing, actual prefix force work, history preservation and immutable refusal/retry. The changed Rust parsed/formatted and received bounded cross-owner/API and independent mechanical review. Native compilation and runtime execution remain unperformed in the full-filesystem environment.


## Complete thermal-fragment intervals through native impacts

`ThermalBrittleTetMotion::advance_released_face_time` advances a complete requested physical duration through native midpoint motion, bounded released-face localization, simultaneous consistent-mass impacts and the remaining motion. The accepted thermal source state is held fixed during the interval. Actual mechanical histories, physical stress, integrated free strain, retained reference density and mass stay with the original thermal owner throughout every internal step.

Thermal and nonthermal time entries share one interval driver. The original nonthermal sustained-contact path remains an optional native callback to that same driver; the thermal impact entry does not substitute impulses for sustained forces. Each accepted record retains its owner's actual native midpoint energy type. No alternate material law, mass solve or contact solver is introduced.

Applied affine forces and prescribed support paths are resampled from the original global interval at every actual prefix and remainder. Initial support velocities must agree with that path, and the same support DOFs retain their velocities during each instantaneous impact. Accepted motion consumes native physical steps; impacts consume only the event allowance. All localization and native response trials spend the complete solve budget. Repeated unresolved events inside the localization-time allowance refuse.

Whole-duration compensated accounts keep mechanical work on TOTAL strain, signed parameter work, eigenstrain work, original fracture loss, contact loss, applied/support mean work, support-impulse work and native kinetic increments distinct. Material balance includes parameter and eigenstrain work; complete balance adds those ports to external work before subtracting storage changes and fracture/contact losses. Actual thermal parameter/eigenstrain accounts are retained even though the current-source interval has no imposed temperature change. Signed and absolute local work defects and cumulative fracture history must pass the original gates before the complete successor is returned. Contact loss is not automatically allocated to heat.

Any late contact, source, geometry, solver, energy, budget or cancellation refusal leaves the original complete motion available for retry. Maximum interval size is caller-selected temporal resolution, not an adaptive error estimate or convergence-order guarantee through impacts. Persistent thermal contact forces, finite sliding, friction, new-pair discovery, global continuous collision detection and automatic temperature-path subdivision remain outside this entry.

Focused regression definitions exercise automatic impact and remainder, original-time affine forcing, retained thermal sources and density, separate contact/fracture energy, and bounded transactional refusal. Rustfmt parsing and source/API/physical review completed; native compilation and runtime tests remain unperformed in the full-filesystem environment.

## Sustained thermal-fragment contact and explicit temperature transitions

`ThermalBrittleTetMotion::advance_released_face_sustained` holds the admitted thermal source state while solving simultaneous frictionless normal forces through the original nonlinear midpoint body response. `advance_released_face_sustained_at_materials` admits an explicitly supplied successor source set in that same coupled response. Each force trial replays the same accepted mechanical history, original reference mass and thermal coefficient/eigenstrain transition. No disposable force proposal can commit restrained damage or replace the thermal owner with a nonthermal material state.

Both standalone thermal and nonthermal entries share the existing sustained-contact force kernel and `fs_contact::sustained::solve`. All non-tied released pairs participate within the original contact/response budgets and hard cap. The normal-map scale supplies numerical conditioning only. Normal forces are interval means, with no penalty spring, artificial compliance, stored contact energy or impulse-to-force conversion.

An internal midpoint candidate seam defers only released-face opening admission for unpublished signed-force probes. Public static, backward-Euler and midpoint endpoints retain their opening checks. Every source, material, history, inertia, native force-equilibrium and thermal work gate still runs during a force response. The final freshly replayed candidate must satisfy strict nonnegative actual openings and normal forces, original force/path identity, loaded resting/matching conditions, complementarity, paired/nodal work and internal resultant/couple gates before publication. Unloaded pairs can separate. Approaching impacts and unresolved mode changes refuse.

The native total-load midpoint receipt remains available beside an external-load-only receipt. Signed parameter work and thermal eigenstrain work remain in the complete physical balance when sources change. Native support reaction work is retained; signed and absolute normal constraint work have separate numerical allowances and are not fracture dissipation or inferred heat. Tensile damage and cumulative fracture loss persist while the original compression branch can recover elastic stress and storage.

The thermal model keeps its existing coefficient-first, then jointly solved total/eigenstrain increment. A supplied source endpoint does not establish exact unresolved temperature evolution, second-order thermal coupling, caloric heat or automatic path subdivision. Fixed matching traces, small strain, frictionless normal contact and original source domains remain the admitted scope.

Focused regression definitions check retained-reference-mass compression compliance and storage after tensile failure; thermally generated contact force with changing source ownership and signed eigenstrain work; unloading separation and transactional public-opening/contact-budget/cancellation refusal. Rustfmt parsing and bounded native API/physical review completed. Native compilation and runtime tests remain unperformed in the full-filesystem environment.

## Stratified hydrostatic loading on deformable tetrahedral solids

`StratifiedHydrostaticPressure` prepares a bounded sequence of positive
piecewise-constant fluid densities below one prescribed free surface.
`HydrostaticDensityLayer::top_depth_m` increases along gravity; the first
layer starts at zero and the final layer extends indefinitely. Increasing and
decreasing density profiles are admitted as prescribed data; no stability or
fluid-evolution claim is attached to a density inversion.

The physical pressure is continuous across the interfaces. With
`q = gravity dot (x - surface)`, each density jump contributes
`(rho_j-rho_previous) * max(q - |gravity| depth_j, 0)`.
Signed density increments multiply the clipped moment AFTER clipping; a
density decrease never reverses the wet half-space. Equal neighboring
densities are compressed at preparation. The explicit nonnegative
`TetPressure::Stratified::load_factor` scales pressure without moving any
interface.

The existing current-area surface owner integrates those moments against the
original P1 basis, including pressure-induced compression as well as buoyancy.
Its follower tangent uses the same density-weighted basis moments and the
shared gravity vector. No displaced material density or separate buoyancy
body force is substituted. Force, torque and deformation therefore derive
from the same closed-surface potential `integral_body p(x) dV`.

Endpoint work uses exact clipped-volume moments, including crossings of
multiple interfaces. Geometry and nodal depth fields are admitted once per
evaluation; each tetrahedral volume is formed once before its layer sum.
Fully wet terms cancel the layer datum through displacement increments.
Compensated sums retain signed profile contributions. The ordinary static
continuation and conservative-work globalization include the new field.

The native prepared energy-momentum owner combines the original clipped
means, discrete-gradient weights and their endpoint derivatives across all
layers before evaluating material stress once per element. The resulting
force and algorithmic tangent share the same endpoint potential. Existing
conservative-pressure vibration, continuation and finite-solid dynamic
consumers retain their original solvers and admission boundaries.

Preparation rejects missing/over-budget layers, nonpositive densities,
unordered depths, invalid coordinates and unrepresentable arithmetic.
Layer traversals poll cancellation in addition to the original geometry
checkpoints. Zero gravity and zero loading retain the unloaded path.
This remains a prescribed hydrostatic reservoir: no evolving fluid inventory,
wave pressure, drag, added mass, wetting or surface tension is inferred.

Six G0/G3 regression definitions in
`tests/hyper3_pressure/stratified.rs` cover independent displaced-volume
force/torque and frame covariance, independent volume-potential derivatives,
interface-crossing work, native static equilibrium and reactions, analytic
and algorithmic tangents, discrete-gradient work, density-decrease and
single-layer limits, budgets and cancellation. Native rustfmt parsed the
changed implementation and new regression source. Independent numerical
reference arithmetic was checked, but Cargo/DSR execution was unavailable;
no passing native regression or benchmark result is claimed here.

## Continuous density gradients in deformable-solid hydrostatics

`StratifiedHydrostaticPressure::with_density_gradient(beta)` adds or replaces
a nonnegative SI density slope on the admitted base layers:
`rho(d) = rho_layer(d) + beta*d`. Its gauge pressure is the existing sum of
layer ramps plus `|g|*beta*max(d,0)^2/2`, with depth measured along unit gravity.
Density drops between base layers remain prescribed data. Negative slopes are
refused for this unbounded reservoir; no fluid stability or mixing is inferred.
`density_at_depth_kg_m3` exposes the same physical density, and all sampled
columns are admitted before force, work or discrete-gradient publication.

The original closed tetrahedral surface remains the geometry owner. Exact
positive clipped-simplex quadratic moments supply traction and its consistent
follower tangent, including partial immersion. No wet-fraction sampling,
auxiliary volume mesh or per-layer material query is introduced. Static
equilibrium, load continuation and native energy-momentum integration consume
the same pressure field. The native discrete gradient differentiates the same
quadratic volume potential, including fully wet and coplanar limits, with its
after-endpoint tangent. Fully wet endpoint work uses polarized depth increments
to retain small motions without subtracting large absolute potentials.

A zero slope keeps the prepared piecewise-constant path, and force-only calls
skip tangent work. The existing layer budget, mesh admission, load scaling and
cancellation boundaries remain in force. The reservoir is prescribed: this
extension supplies neither evolving fluid density nor fluid inertia.

Five regressions in `tests/hyper3_pressure/stratified/gradient.rs` define
analytic Archimedes force, density-weighted torque and heave stiffness;
independent potential/force/tangent/work checks through changing immersion;
native finite-step work and endpoint tangents; exact coplanar emergence; and
zero, overflow and cancellation limits. All seven changed Rust files parsed
successfully with rustfmt. Independent arithmetic checks compared clipped
moments with divided differences and checked discrete work/reversal. Native
Cargo compilation and regression execution remain outstanding: DSR/RCH are
unavailable and the shared filesystem is full.

## Consistent reference inertia for the flat shell and moving loads

`ShellPlate::assemble_with_mass_formulation` selects either the original
`ShellMassFormulation::Lumped` inertia or `ConsistentReference`. The original
`assemble`, `MovingShellModel::new` and `new_with_payloads` entries retain
lumped inertia. Stiffness, geometry admission and deterministic element order
are shared between the formulations.

The consistent image integrates the existing P1 velocity fields exactly:
`T = integral_A (rho*h*|v_h|^2/2 + rho*h^3*|omega_h|^2/24) dA`.
Each triangle uses `integral_A N_i*N_j dA = A*(1 + delta_ij)/12`, retaining
cross-node translational and rotary mass. This preserves the total mass and
uniform-velocity inventory of the lumped image while resolving affine
velocity and rigid-rotation kinetic energy without vertex-quadrature bias.
The fixed local formula introduces no quadrature allocation or new mass solver.

`MovingShellModel::new_with_mass_formulation` and
`new_with_payloads_and_mass_formulation` carry that same native mass through
support restriction, prescribed support coupling, initial acceleration,
cached Newmark/mass solves, material-point inverse inertia and event impulses.
Payload COM and full rotational inertia are added once to the selected full
mass. Rayleigh damping uses the selected bare-shell mass; attached payloads
retain the existing absence of inferred payload damping.

Both formulations retain the shell's original isotropic rotary surrogate,
including artificial drilling inertia about the plate normal. This changes
the integration of that declared model, not its plate/shell constitutive
scope: no Kirchhoff--Love continuum, finite-rotation shell or
configuration-dependent inertia is introduced. Existing dense assembly
budgets, complete-result publication and cancellation boundaries remain.

Regressions in `tests/shell_operator.rs` compare kinetic energy with independent
degree-two triangle quadrature, uniform inventories, rotated/refined rigid
motion, support restriction and selected-mass damping. The existing moving-load
payload and event tests exercise both formulations. An independent analytic
triangle-inverse test verifies the complete post-impulse velocity field and
its kinetic energy, then continues the same state through passive Newmark
motion with the selected Rayleigh damping.

## Retained finite porous solids

`hyper3::plastic::porous::FinitePorousTetSolid` retains original per-element
native or source-resolved `FinitePorousMaterial` owners, non-isochoric plastic
metrics, accumulated work-conjugate flow, configuration and Piola stress.
Complete source cards and their fixed physical contexts remain borrowed across
load, unload and reload; constitutive trials never repeat source resolution.
Virgin reference preparation and every
absolute load/support endpoint reuse the existing tetrahedral geometry,
Newton, line search and prepared-element storage. Both porous and J2 owners
use the same mechanical endpoint/work implementation; neither infers plastic
volume from total volume nor normalizes the porous metric to unit determinant.

Every trial borrows the old accepted history. A candidate is published only
after all material-domain/yield, independent storage/potential, equilibrium
work and signed endpoint-work-defect gates pass. Rejection or cancellation
preserves the complete old solid. Original reference volume integrates energy;
physical `yield_stress * delta_alpha` loss remains separate from hardening
storage and the signed numerical endpoint defect. Prestressed finite rotations
can have negative endpoint-work defect without generating plastic heat.

Accepted total volume `J` and permanent plastic volume `Jp` are retained from
the native material response without conflation. The density query maps an
explicit original reference inventory to `rho_ref/J`; it does not use Jp as the
current volume or invent a source density calibration. The source-driven
`hyper_porous` cases check the original material/context binding, conserved
density, independent hydrostatic storage/dissipation/work, late work refusal
and deterministic retry, and per-element source retention without averaging.

The `hyper_porous_retained` regressions specify force-driven compression,
permanent volume loss on unloading, elastic reload against independent
hydrostatic formulas, heterogeneous pressure sensitivity, rigid-motion
covariance, physical work closure, refused-step retry and cancellation.

`FinitePorousTetSolid::advance_with_pressure` now drives this same accepted
history with current-area follower pressure on its exact original closed mesh.
It reuses the native pressure force, analytic geometric tangent, globalization
and support reactions for uniform, per-face and admitted spatial fields.
Uniform pressure is Cauchy traction: its exact held-interval work is
`-p * (V_after - V_before)`, not a frozen referential nodal force. Each material
trial still borrows the old accepted plastic metric; no Newton iteration
commits compaction. Pressure changes are resolved by successive held intervals.
The current-area virtual-work and load-stiffness convention agrees with the
[Abaqus Theory Guide, Pressure load stiffness](https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-pressloadsolids.htm).

The returned optional pressure account separates exact held-field work,
endpoint-force work, their signed geometric correction, and the total external
work defect. Dead-force work retains its original meaning, while the common
equilibrium residual includes pressure endpoint work. Acceptance bounds the
SUM of absolute material and pressure endpoint defects by the original work
budget, preventing cancellation from admitting an under-resolved step. Neither
defect becomes plastic heat. Per-face work follows the original straight
displacement-segment convention; it is not path-independent work. Focused
native regressions cover pressure-driven compaction, zero-pressure release,
elastic reload, exact volume work, patch/support accounting and atomic refusal.

This is bounded fixed-temperature, rate-independent associated porous
plasticity on displacement P1 tetrahedra. It does not infer evolving porosity,
cell-wall densification, asymmetric yield, damage, contact, pore-fluid pressure,
thermal coupling, inertia, locking-free accuracy or experimental calibration.

## Retained membrane and bending plasticity on native CST/DKT shells

`plastic_shell::SmallStrainShellProblem` consumes the original
`fs_plate::shell::ShellMesh` facet frames and `dkt_curvature_quadrature`.
Every original midside point carries an independent
`fs_material::plastic::section::LayeredPlaneStressState`, retaining all
three-dimensional transverse strains and material histories through thickness.
Native plane-stress return maps supply both stress resultants and the complete
consistent membrane/bending tangent. Engineering shear, asymmetric stacks,
yield-induced coupling and genuinely nonsymmetric native laws remain intact.

The existing solid increment owner supplies sparse corrections, line search,
frozen-history trials, physical load entries and atomic publication. Six
physical translation/axial-rotation coordinates are normalized by an explicit
positive length: `(u,L theta)` pairs with `(F,M/L)`. This preserves work and
makes the equilibrium residual a force norm. Output coordinates and reactions
retain their original metre/radian and newton/newton-metre meanings. Restart
admission checks complete section histories against the same reference geometry.
Small physical rotations, displacement slopes and material strains have
separate explicit limits; a large rigid spin cannot pass as zero strain.

The original drilling-spin penalty uses the virgin section shear stiffness,
stays independent of yielding, and is numerical elastic storage. Generic
force/history solves do not infer plastic heat. The lane does not claim
finite-strain forming, evolving contact, transverse-shear response, material
calibration, plastic-limit-load convergence or through-thickness error bounds.
Refine the explicit thickness integration and physical load path as required.

Five `plastic_shell` regressions cover the original full elastic stiffness
in flat/oblique frames, eccentric physical layer signs and virtual work,
free-DOF cantilever yielding and permanent shape after unloading, numerical
coordinate-length invariance, complete restart and transactional refusal/
cancellation. They are test definitions; execution requires the native build.


### Native plastic-shell work and separate physical dissipation

`plastic_shell::SmallStrainShellProblem::solve_with_energy` solves exactly one
physical load entry from virgin state or a retained checkpoint and publishes only after the native material,
section and whole-shell work gates pass. It uses the same prepared facet
operators as the mechanical solve, original area/thickness weights, and the
exact accepted native histories. Applied forces and moments plus prescribed
coordinate reactions supply the external endpoint work.

The account reports recoverable elastic/hardening storage, native physical
plastic/overstress/recovery dissipation, material numerical loss, and signed
transverse residual work separately. Virgin drilling stiffness supplies only
numerical elastic storage and its explicit backward-Euler endpoint loss; both
remain outside the physical material heat. Neither an equilibrium-work deficit
nor a stabilization term becomes an inferred heat source.

The same reference geometry, section assignments, coefficients and energetic
interpretation must own both endpoints. This mechanical account performs no
thermal transport or temperature evolution. Two added energy regressions define
yielded cantilever load/unload work against independent yield-increment heat
and pure-drilling storage/loss with zero material heat. Source/API review and
rustfmt checks are complete; native compilation and regression execution remain
unverified because the runner exhausted its available disk and memory.


## Oriented Hill48 compression, release and confinement

`tests/plastic_hill_compression.rs` extends the original directional cases
through `OrientedHill48` and the native P1 tetrahedral and P1/Q1 plane solvers.
Three principal material orientations and an independent 45-degree oracle
exercise compressive yield, elastic unloading, reversal and reaction balance.
Lateral unknowns remain free apart from rigid-motion constraints, so oblique
plastic shear and unequal width/thickness flow are solved rather than pinned
away. Comparisons convert original material-coordinate history through
`plastic_strain_spatial`. A rotated Q1 Standard/BBar affine comparison retains
plane-strain confinement and transverse plastic history.

`examples/anisotropic_compression.rs` resolves nine original scalar properties
through `plastic::hill::source::planar`, including explicit plastic r0/r45/r90
and both independent transverse shear yield ratios. The 100 x 20 x 10 mm
coupon is compressed and released in rolling, cross, thickness and 45-degree
orientations. CSV output separates stress, reaction, width/thickness strain,
plastic shear, Hill accumulation, permanent dimensions, lateral end offset
and residual released stress. Width and thickness use same-section gauges,
so plastic shear is not mislabeled as a width change. Run with
`cargo run -p fs-solid --example anisotropic_compression -- 0.004`.

These are synthetic constitutive examples, not measurements or calibrated
forming limits. The fixed-axis small-strain symmetric-yield model does not
add finite forging, contact, evolving texture, buckling, necking or damage.
The additional FE cases are authored and statically reviewed; their full
native Cargo execution remains outstanding under shared executor resource
constraints. Focused constitutive execution is not an FE integration pass.


### Retained plastic-shell momentum and physical kinetic work

`plastic_shell::dynamics` prepares immutable original layer-density moments and
one full consistent sparse mass before essential-boundary elimination.
Reference geometry, physical layer bounds and numerical coordinate length are
checked on reuse. Signed first moments retain translation/director coupling
for eccentric or unequal-density stacks. Reference inventory is actual
material mass; the numerical relative-drilling gauge adds no matter.

`PlasticShellDensity::PerElementCorners` supplies one positive original density
triplet for each physical layer of each facet. Density interpolates in the
original triangle and may jump between facets. The actual layer bounds give
three distinct zeroth, signed first and second density moments; the original
native shell quadrature retains their full P1 fields in momentum, support
reactions and physical/gauge kinetic work. No facet or layer mean replaces the
original distribution. Constant inputs keep their compact scalar inventory;
corner inputs retain their exact triplets and both expose an allocation-free
corner iterator. Sourced thermoplastic/contact motion compares every original
corner against its source density, rejecting inventories with matching layer
means but different spatial inertia. This does not add in-plane caloric
reconstruction, finite rotations or matter transfer. Existing dynamics tests
now exercise graded two-layer fields against exact cubic simplex/thickness
integrals and native moving-support equilibrium, plus source/count refusals.

`advance_backward_euler` uses the original shell constitutive history and
solid Newton owner with the original dimension-independent reference inertia.
The same mass supplies force, M/dt² tangent, free-row force scaling and full
support reactions. Coordinates and velocities retain physical translations
and axial rotations; the explicit numerical length transforms mass by
congruence. The physical load path contains exactly one factor-one endpoint.
Material trials always start from the accepted old history. Failure or
cancellation leaves the complete supplied checkpoint and reference mass intact.
The dynamic tangent budget includes both material and retained inertia entries.

The native one-pass kinetic integral reports physical and numerical drilling
storage, changes, backward-Euler losses and endpoint work separately. Both
partitions must agree with the mass actually used by the solve. Only one global
mass is retained; compact original facet moments drive kinetic partitioning
without extra sparse mass matrices. Nearly rigid gauge motion uses common
before/after point amplitudes rather than independently subtracting large terms.

`advance_backward_euler_with_energy` shares prepared geometry between momentum
and the original section/material energy audit. It gates applied/support work,
physical kinetics, numerical gauge kinetics, recoverable material storage,
physical dissipation and numerical losses before returning a candidate.
Supported timed native laws must expose the same physical duration as the
inertia interval. Neither numerical kinetic damping nor stabilization storage
is deposited as material heat. Temperature evolution and heat transport are
not performed. This remains first-order small-strain, small-rotation CST/DKT
motion with the native P1 midsurface/director kinetic field; finite forming,
contact impulses and a temporal-error or performance certificate are outside
this extension.

Four dynamics and three dynamic-energy regressions define unrestrained oblique
ballistic motion, consistent-mass moving-support equilibrium and reactions,
yielded release with coordinate-scale invariance, retained-mass/refusal and
cancellation behavior, yielded bending/unloading with independent yield heat,
physical/gauge kinetic separation, and material/inertial duration mismatch.
Rustfmt, scoped whitespace and source/API/physical-account reviews completed.
Native compilation and execution remain unverified under the shared runner's
disk and memory exhaustion; these are regression definitions, not pass results.


## Retained numerical drilling across thermal plasticity

`plastic_shell::drilling::PlasticShellReferenceDrilling` binds the original numerical spin-gauge weights to the reference mesh, physical layer bounds, coordinate length and drilling coefficient. Native static and backward-Euler entry points retain these weights across thermal/material-coefficient changes, reusing the original CST/DKT, Newton and consistent-inertia owners. Preparation with a retained gauge skips evaluation of the endpoint virgin material. Gauge endpoint work, storage and integration loss use one shared native integration kernel and are independently reported; none becomes physical heat. This retains the original small-strain/rotation and thickness quadrature scope, not a finite-forming claim.

`PerElementBorrowed` admits separately owned thermal sections without cloning material laws. Independent-mesh-lifetime constructors allow the gauge and mass to outlive temporary endpoint laws. The mass owner exposes its exact original per-layer density inventory: coupled consumers must match each layer, since equal total mass alone does not establish equal eccentric inertia. Three focused native regressions cover response under material softening, dynamic parity and temporary law lifetimes, incompatible references, incomplete fields and cancellation. These are numerical model checks, not experimental accuracy or throughput claims.


### Retained finite porous inertia and compaction

`hyper3::plastic::porous::dynamics::FinitePorousTetMotion` advances the original finite porous tetrahedra with the existing P1 consistent reference mass and backward-Euler momentum owner. Native equilibrium trials freeze the accepted non-isochoric plastic metric; admitted endpoints retain full hardening/compaction history, total and plastic Jacobians, displacement, velocity and original mass. Current density is reference density divided by total J, never plastic Jp. Original mass and density are prepared once and shared by later steps.

Physical yield dissipation, signed material endpoint-work defect and numerical kinetic loss retain independent work accounts and admission bounds; numerical loss is not heat. Prescribed support reactions include the full consistent-mass action, and free bodies require no artificial restraints. Refusal and cancellation leave the accepted state retryable. The current porous source schema has no density capability, so density is an explicitly supplied reference inventory, not a sourced calibration claim. Fixed-condition porous compression/release does not add evolving porosity, densification, thermal evolution, contact, mass transfer or pressure stabilization. Four native regressions cover independent hydrostatic and consistent-mass formulas, compacted elastic release, free motion and atomic refusal; their execution remains pending the shared dependency build resource blocker.

Both cellwise inventories and positive original P1 nodal density fields are
supported. The shared mass owner integrates `rho_ref * N_a * N_b` exactly;
nodal values and the full consistent matrix survive compaction without being
replaced by cell means. Density queries return current cell means through
`rho_ref/J`. `reference_acceleration_load_n` projects `M*a` through that same
matrix, including graded mass and support coupling.

Three additional `hyper_porous_dynamics` regressions specify source/context
retention, graded free-body acceleration, force-driven crush/unload/reload,
independent material and kinetic work, consistent support reactions, density,
late kinetic-budget refusal/retry, admission boundaries and cancellation.


## Primary and secondary creep in checked shell motion

The existing layered CST/DKT `solve_with_energy` and `advance_backward_euler_with_energy` paths accept native standalone or joint Norton creep through the shared section energy adapter. Primary strain resistance, signed creep strain and independent plastic history remain with their original constitutive owners. Every timed law must match the inertia interval. Layered stress/resultants, elastic/plastic-hardening storage and physical creep/plastic losses are integrated with the original thickness and surface rules; material/inertial numerical loss and drilling storage are never physical heat.

Three new `plastic_shell/creep` regressions define an actual loaded primary-creep membrane followed by clamped relaxation and permanent unloading, a jointly creeping/yielding bending cantilever, and prescribed creep motion with an independent P1 kinetic integral and duration-mismatch refusal. Isothermal sections do not evolve temperature; sourced thermal feedback is owned by `fs-couple::shell_plastic_enthalpy`. Existing small-strain/rotation, quadrature and timestep limits apply. The regressions were not executed in the implementation environment; no performance or experimental prediction claim follows from formatting/source checks.

## Full orthotropic strength interaction on native solid stresses

`strength::orthotropic::OrthotropicTsaiWuMaterialBinding` binds the shared
three-dimensional material criterion to its exact source-resolved, oriented
stiffness. The native P1 and stiffness-quadrature response owners expose
`assess_orthotropic_strength`; the existing tensile-only entry point uses the
same shared sampler. The sealed criterion type retains its own assessment,
ordered property interpretation and identity domain. No dynamic dispatch,
second mechanical solve or intermediate full stress report is introduced.

Actual Mandel stress is converted once to tensor shear and pulled back with
the original material rotation. Each element or quadrature point retains its
own source condition, orientation, complete normal/shear stress and selected
strength surface. Neither averaged stress nor principal-stress axes replace
the original samples. Missing/stale bindings, unresolved interaction convexity,
sample budgets and cancellation refuse before a result is returned.

Three new cases in `tests/orthotropic_tensile_strength/tsai_wu.rs` exercise a
triaxial tensile crossing that independent axis checks miss, rigidly rotated
native stress, local temperature/orientation at Q1 points, ordered-source
identity and atomic refusals. The existing P1/Q1/Hex20/Hex27 tensile cases
continue through the same sampler. These are manufactured model checks;
the observations add no failure-mode selection, damage evolution, fracture
energy, post-failure redistribution or continuum maximum bound.


### Sustained native plastic-shell contact

`plastic_shell::contact` couples the original sustained unilateral/Coulomb solver to native layered-shell equilibrium and backward-Euler momentum. Immutable material-point/plane witnesses retain IDs, physical patch areas, reference mesh and layer bounds. The existing gauge-free P1 shell point rows observe translations and physical offset rotations; their exact transpose applies forces and moments without directly loading numerical drilling. Original element geometry is prepared once per contact interval. Every force proposal restarts from the same accepted complete material history, and the returned endpoint must match the solver's final mechanical force/path replay.

Each contact reports solved force and pressure (normal force divided by explicit area), old/new gap, relative slip, body and obstacle work, physical Coulomb loss, normal complementarity work and numerical approach loss. Tangential and normal work defects have independent gates. Aggregate normal approach loss fn*g0 has its own finite budget, zero by default; it is never frictional or plastic heat. Retained drilling stiffness and reference mass stay unchanged, and failed solves, budgets, work gates or cancellation return no accepted state. Fixed frames and an affine point/plane gap path do not establish contact search, global CCD, enriched DKT interior trace accuracy, impact time or restitution. Three native regressions cover plastic pressing/sliding/release, analytical consistent-mass contact and atomic refusal; runtime execution remains pending the shared dependency build resource blocker.


### Identical shell-contact gap admission across retained intervals

Sustained shell contact evaluates final and restarted initial gaps through the same absolute original point/plane arithmetic. It does not accumulate a gap by subtracting endpoint displacements and then re-add it to an earlier rounded gap: that could accept a zero gap which recomputed as a negative gap on the next interval. Tangential increments and relative work remain unchanged, and geometry is never clamped. An eight-interval native J2 elastic load/unload/reload regression checks bit-identical carried gaps, independent force/displacement relations, unchanged plastic history and the explicit numerical approach budget. The regression is authored; native execution remains pending the shared build resource blocker.

## Native frame motion with one physical creep interval

`Frame2d::advance_creep_backward_euler` advances displacement, velocity,
acceleration and original standalone/joint/thermal creep fiber histories as one
immutable proposal. The native reference mass, force-based fiber assembly,
consistent tangent and Newton backtracking remain the mechanical owners. The
endpoint acceleration solves `M a1 + f_int(u1) = load`, with
`v1 = v0 + dt a1` and `u1 = u0 + dt v1`. Every trial and the independent fiber
work audit start from the same accepted history. There is no initial timed
constitutive evaluation that consumes an additional creep interval.

All fibers and inertia must have the exact same positive physical duration.
Separate force, moment, impulse and energy gates check the accepted endpoint.
Kinetic change uses average velocity contracted with the actual momentum
increment. Backward-Euler inertia loss `delta(v)^T M delta(v)/2` is explicit,
nonnegative and separate from native material numerical loss; neither is heat.
Physical creep/plastic work and reciprocal thermal transfer retain their native
partition. Failed trials and cancellation publish no changed frame.

Stationary supports, linear small-strain/rotation geometry and native
Euler–Bernoulli fibers are admitted with lumped or consistent reference inertia.
Attached section viscosity, transverse shear and retained Rayleigh damping
require separate work owners and refuse in this route. The existing untimed
and Newmark guards remain in force. This first-order dissipative integrator
requires physical timestep refinement; it does not establish finite-rotation
creep dynamics or coupled thermal transport. Focused analytical and history
regressions are in `tests/frame_creep_dynamics.rs`; native execution remains
pending an available prescribed build runner.

## Recovering primary creep in native frames

`FiberLaw::RecoveringNorton` carries the original free-sided recovering law and
all three axial histories through the native section and frame paths. Loaded
holds, clamped relaxation, unloaded dwell and reloading use the same constitutive
owner. `solve_creep_step` and `advance_creep_backward_euler` retain the actual
physical duration, consistent tangent and complete frozen trial state. Their
work gates delegate to the original tensor recovery account; resistance decay
supplies neither physical heat nor stored energy.

`Frame2d::creep_recovery_fiber` inspects accepted histories without advancing
creep. Rebinding preserves original kinetic/recovery coefficients, primary
resistance calibration, source claim identities and every fixed nonthermal
condition, while allowing a supported temperature or duration to change.
Authored density remains authored even when the law itself retains all thirteen
source receipts. Untimed static and Newmark paths continue to refuse timed
fibers. Four manufactured frame regressions specify load/hold/dwell/reload,
variable-duration momentum with both mass formulations, source refresh and
atomic refusal. Scoped parsing and source checks passed; native execution is
pending the prescribed build runner.

## Energetic Hill plastic shells

The existing generic plastic-shell equilibrium, energy and retained-inertia
owners consume native fixed-axis Hill laws through the material section-energy
adapter. `tests/plastic_shell/hill.rs` exercises directional membrane flow with
free width/shear, genuine zero-load release and permanent extension, reversed
oriented layer stacks with plastic membrane/bending coupling, and backward-Euler
motion with original reference inertia. Physical yield heat remains distinct
from hardening storage, numerical material loss, kinetic loss and drilling-gauge
energy. Local Hill histories survive every accepted point update.

Scope remains small strain and small rotation with an authored fixed material
metric and frame; these fixtures do not calibrate sheet forming or texture.

## Directional strength on native tensor-viscoelastic histories

`ViscoElements::prepare_orthotropic_strength` binds the existing directional
tensile or full three-dimensional Tsai-Wu criterion to each original source
element and its material frame. Immutable material-card identity and the
complete typed physical query must match the current native coefficients.
The prepared observer borrows that exact geometry/coefficient owner, retains
the complete original source receipts, and caches the inverse material frame.
A foreign or refreshed system requires its own prepared binding.

`PreparedViscoStrength::assess` reads each accepted native integration-point
stress, including equilibrium elasticity, Maxwell history, Kelvin stress from
the actual last strain ramp, and imposed thermal eigenstrain. It converts
work-conjugate Mandel shear to physical tensor shear exactly once, then pulls
the complete stress into the original material axes before applying the shared
criterion. Signed normals and shear remain visible; samples are never averaged.
The original P1, quadratic, hexahedral and B-bar point ordering survives, with
deterministic governing-point ties and explicit coefficient/strength identities.

Assessment performs no material query, constitutive update, force solve,
operator construction or displacement-based stress reconstruction. Only its
bounded output vector allocates; cancellation is checked per point and before
return. If an imposed thermal field is present, every actual sample temperature
must equal the current coefficient/strength query. Frozen coefficients do not
authorize a strength observation at another temperature.

A law-changing Kelvin parameter refresh intentionally has no native endpoint
rate. A retained transactional stress-observation flag therefore refuses
assessment until a successful positive-duration native step supplies complete
stress. Unchanged-law rebinds preserve the flag; source refresh without Kelvin
continues to have physically determined spring/Maxwell stress. Relaxed initial
states, clones and rejected trials preserve their existing semantics. Reference
mass, geometry, histories, work and material evolution are unchanged.

`visco_source/strength.rs` supplies four analytic native regressions: a
Maxwell-plus-Kelvin strength crossing missed by equilibrium-only checking and
subsequent fixed-displacement relaxation; rotated signed tensor stress and
shear-sensitive Tsai-Wu interaction; unaveraged thermal/B-bar quadrature with
stale-temperature refusal; and source refresh, actual resumed viscous stress,
owner mismatch, budgets and cancellation. These observe represented endpoint
stresses, not unsampled space/time peaks, rate-dependent strength calibration,
progressive damage or fracture evolution.


### Joint and thermal recovery throughout native member paths

The shared timed-fiber adapter also admits `RecoveringNortonJ2`,
`ThermalRecoveringNorton` and `ThermalRecoveringNortonJ2`. Section response,
static creep steps and backward-Euler motion keep their existing assembly,
mass, Newton and work owners. The joint fibers retain both permanent strain
paths and hardening storage; the thermal fibers additionally retain their
actual accepted temperature, original expansion receipts and reciprocal heat.
Rebinding never advances any history or silently replaces the recovery family,
kinetic calibration, energetic hardening, fixed source context or expansion law.

`creep_recovery_plastic_fiber`, `thermal_creep_recovery_fiber` and
`thermal_creep_recovery_plastic_fiber` expose accepted laws and histories without
another timed trial. Independent member regressions cover load/clamp/dwell/reload,
variable-duration joint momentum, heated yielding and an unloaded temperature
gradient bending a previously crept cantilever while local resistance recovers
at different rates. Physical heat, signed reciprocal exchange and material/
inertial numerical losses remain distinct. Temperature is prescribed between
intervals; this path reports its heat account but does not add caloric transport.
All original small-geometry, source, clock and atomic refusal gates remain.
Native compilation and test execution are still pending the prescribed runner.

## Independent boundary histories for native cyclic solid fatigue

`J2FatigueLoading` supplies reference prescribed displacements in the original
support order and dead nodal forces in the original mesh frame. The quasistatic
`J2FatiguePath::advance_with_loading` and `advance_cycle_with_loading` consume
the next original load factors; unit factors admit absolute endpoint histories.
This enables independently phased normal and shear loading while preserving
the original support DOFs, source-bound J2/Voce/Chaboche laws and complete
constitutive history. Additive duplicate nodal loads retain native semantics.

Every endpoint still calls the native `solve_from` owner. Its actual accepted
stress and total mechanical strain feed the same fixed-plane, representative
SWT, Fatemi-Socie and variable-amplitude banks. Changing boundary values does
not inject stress, erase plastic memory or imply cyclic stabilization. Original
paired Morrow/SWT rainflow refusals remain in force for unsupported phase-lag
histories; representative tensor-cycle models retain their own closure rules.

One reusable prescribed-field buffer serves an entire independent-loading
block. All input fields are checked before native work, and one original
mechanics/fatigue checkpoint covers the block. A later native solve, gauge,
closure or cancellation refusal publishes no partial state or samples. Ordinary
`advance` continues to use the original reference fields. Source validity,
fixed-temperature calibration, small-strain geometry, finite plane/sample/work
budgets and the need to resolve actual physical turning points are unchanged.

The physical-time `J2DynamicFatigueRun` accepts the same fields through
`J2DynamicFatigueLoadingTarget`, `advance_with_loading`,
`advance_cycle_with_loading` and `advance_to_with_loading`. Each target declares
a strictly increasing finite time and an absolute reference-field multiplier.
Scaling occurs once before the original backward-Euler solve with its unit
load factor. The previous full displacement, velocity and constitutive state
remain the initial conditions; reference density and support topology remain
fixed. Moving-support mass coupling, reactions, support work and the native
inertial energy receipt therefore accompany the actual fatigue samples.

Original proportional targets and independent fields share one native
transaction loop. A cheap borrowed endpoint iterator avoids allocating a
second target trace, and one prescribed buffer plus a force buffer sized for
the block are reused. A late mechanical, fatigue or closure refusal restores
all banks and leaves accepted time, velocity, state and last energy unchanged.
Backward-Euler numerical damping remains separate from physical dissipation
and empirical fatigue damage; interval refinement and resolved extrema remain
the caller's responsibility.

The static and dynamic loading regressions include an analytic nonproportional
shear-circle reference, native plastic continuation with independent forces,
an independent consistent-mass recurrence with moving supports, original
energy equality and late refusal/retry across every fatigue bank. Coefficients
are synthetic implementation references, not measured service-life evidence.

## Tabulated malleability in native structural consumers

The original small-strain continuum and plane-stress owners accept
`fs_material::plastic::tabulated::TabulatedJ2` through `SmallStrainLaw`.
Layered membrane/bending sections retain the same native plane-stress histories
and energetic hardening, including partly yielded thickness and springback.
`FiberLaw::Tabulated` carries its exact free-sided axial restriction through
existing section and beam/frame mechanics; checked responses and commits
refuse unsupported accumulated plastic strain without changing accepted fibers.
No structural consumer refits the curve or inserts missing material data.

The original finite tetrahedron and in-plane sheet owners accept the explicit
Kirchhoff table selected on `FiniteJ2`; their retained plastic metric, physical
work gates, reference volumes and traction-free thickness remain unchanged.
Existing locking, small-shell-rotation, finite-sheet geometry and material
applicability limits continue to apply. Synthetic regressions in
`tests/plastic_tabulated.rs` exercise native plane-stress loading/release,
permanent lateral and thickness strain, partial section yielding and bending
springback, force-based beam/frame release and atomic fiber refusal.
`tests/hyper_plastic.rs` and `tests/finite_plastic_sheet.rs` add finite
compression/release checks with solved free sides and thickness. These numerical
checks do not establish metal-specific calibration or finite shell forming,
contact geometry, necking or fracture limits.

## Native plane-continuum fatigue

`plastic2::fatigue::PlaneFatiguePath` carries shared finite-plane strain-life
fatigue through the original P1/Q1 nonlinear load path. It observes each accepted
integration point's complete signed three-dimensional constitutive strain and
physical stress: solved thickness contraction under plane stress, constrained
transverse stress under plane strain, and the actual assumed strain under B-bar.
J2, Voce and Chaboche retain their native histories and original return map.
Geometry is prepared once and shared by restart clones; only the latest accepted
solution, material samples and bounded rainflow residues are retained.

Each gauge binds its original material-state identity, chemistry, process,
absolute temperature and elastic modulus to the original strain-life calibration.
Source-backed curves additionally retain card, state and fixed query-coordinate
identity. The existing finite-plane CMB, Morrow and Smith-Watson-Topper laws own
cycle counting and life inversion; unsigned equivalent or instantaneous principal
stress does not replace the material-point tensor history. Caller-declared
normals are a finite sampling of candidate planes, not a continuous orientation
maximum.

The initial zero-load equilibrium is a real sample. Chronological blocks retain
every accepted intermediate endpoint, and a late mechanics, gauge, budget or
cancellation refusal rolls back the complete requested block. Prefix previews
include endpoint half cycles without closing the live history; finalization
requires the full declared path and commits every gauge atomically. Aggregate
plane/residue/projection-work admission bounds memory and work before a solve.

`tests/plastic_plane_fatigue.rs` supplies native P1/Q1 Prager cycle oracles,
plane-strain transverse-stress/SWT and B-bar checks, restart/chunk equivalence,
and source, budget, cancellation and late-refusal regressions. These synthetic
implementation checks do not establish transferable component service life.
The observer changes neither stiffness nor heat and adds no thermal/rate law,
crack growth, remeshing or second constitutive owner.

## Pointwise HEX8 thermal fields for native directional strength

`Hex8ElasticProblem::solve_thermal_response` and its dead-nodal-force variant
accept the original `QuadraticThermalStrainField` as well as existing uniform
and elementwise thermal inputs. For Standard strain, the unchanged native
isoparametric owner consumes each original stiffness-point eigenstrain and
matching local material tensor directly. Thermal stresses, support reactions,
energy and material/strength identities retain their actual point values;
a zero mean cannot conceal a local directional tensile peak.

This is a source-compatible public admission of the existing pointwise
implementation. B-bar still requires its original elementwise thermal field;
no within-element pressure/eigenstrain projection is invented. Mismatched
quadrature counts and unsupported B-bar point fields retain explicit refusal.
The native directional-strength tests cover Q1/Q2 zero-mean local peaks,
per-point material/temperature binding, full Tsai-Wu interaction and these
quadrature/formulation boundaries.

### Native frame full-stress Findley observation

`FrameFatiguePath::new_with_findley` adds fixed material-point stress windows to
its existing fiber and elastic S-N banks. Accepted Lobatto axial/bending strains
supply the original layer's normal stress; the canonical Timoshenko transverse
resultant supplies both local shear components through the shared homogeneous
isotropic circle/annulus Saint-Venant owner. The authored shear correction is
not a point-stress multiplier and recovered local storage is not added to the
native frame energy. Original material identity, typed temperature and available
calibration source pins are retained. Nonisotropic or heterogeneous elastic
sections, viscous members and unsupported source conditions refuse.

The selected annulus must satisfy the declared normal and engineering-shear
bounds even when the observed point lies on a zero-stress contour. Findley
fatigue-limit usage remains separate from Miner damage. Closed representative
windows can be finished or restarted explicitly. Gauge count, aggregate retained
samples and each original Findley work budget are bounded. Static load,
displacement control and Newmark all observe the same accepted native steps;
a late source, strain, budget or cancellation failure rolls back mechanics,
scalar counters and appended tensor samples together. Tensor rollback retains
only prefix lengths rather than copying all previous windows. Synthetic tests
in `tests/frame_fatigue/shear.rs` cover independent off-axis stresses, correction
invariance, Newmark equivalence, source/budget admission and transactional
rollback. This observation supplies no variable-amplitude multiaxial life,
plastic stress recovery, thermal recoupling or cross-section inertia correction.

### Native rod fatigue with bonded layered transverse shear

Explicit `new_with_layered_transverse_shear` and
`new_motion_with_layered_transverse_shear` constructors connect original
static and accepted-motion transverse resultants to the shared bonded
common-Poisson circular stress owner. Each original layer retains its E/G,
point, material source and fatigue calibration. Recovered transverse stress
joins the existing axial/bending/torsion tensor for signed S-N projections and
Findley windows. Existing source pins, combined whole-layer strain bounds and
atomic mechanical/history acceptance remain authoritative. Earlier homogeneous
and explicitly unassessed modes retain their original behavior. The new
independent two-layer native regression is in `tests/rod_fatigue.rs`; this
observation does not alter native stiffness, work, inertia or stored energy.

## Retained porous motion under prescribed follower pressure

`FinitePorousTetMotion::advance_with_pressure` combines the original current-area
surface load and analytic follower tangent with retained porous compaction and
consistent reference inertia in one native backward-Euler solve. The same
physical load fields, exact original mesh binding and declared-path pressure
work apply. Pressure is held over the requested interval; changing reservoir
pressure requires resolved successive intervals. Original reference mass and
accepted plastic history remain unchanged until every solve and work gate passes.

The complete balance keeps storage, physical plastic loss, kinetic change,
backward-Euler numerical kinetic loss, material endpoint error and pressure
endpoint correction separate. The original material-work allowance bounds the
sum of absolute material and pressure corrections; the numerical kinetic-loss
allowance remains independent. Numerical losses do not become material heat.
The pressure account's external defect excludes any inertial endpoint work.

The `hyper_porous_retained::dynamics` regressions specify free regular-tet
pressure compaction against an independent consistent-mass hydrostatic formula,
a nonsymmetric one-face pressure patch, support reactions with off-diagonal
mass coupling, retained history, exact pressure work, late refusal and retry.
The scope is prescribed-pressure motion at fixed material conditions; sealed
cavity and spatial-gravity inertia remain outside this entry point.

## Sealed cavities with retained porous compaction

`hyper3::plastic::porous::FinitePorousTetSolid::advance_with_cavities`
couples independently sealed, explicitly resolved gas cavities to the original
retained finite porous material history. It reuses the native cavity batch,
current-area force and geometric tangent, independent acoustic rank-one factors,
and original Newton/globalization owner. Every constitutive trial starts from
the same accepted plastic metric.

Each cavity binds the exact original solid reference chart and an unchanged gas
parcel at its undeformed cavity volume. Original ideal and Peng--Robinson source
owners supply calorics, phase-domain admission, entropy closure and acoustic
compressibility. Accepted deformed parcels are output evidence, not replacement
undeformed reference states. Distinct chamber inventories do not mix.

Per increment, gas work on the solid is `-delta U_gas` from the original source's
direct stored-energy increment; fixed exterior work is
`-p_exterior * delta V_gas` from retained compact geometry. Their compensated sum
is separate from physical porous yield loss. Each source records its endpoint
force work and signed work correction. The existing work budget bounds the sum
of absolute material, held-pressure and individual cavity corrections, so
opposing chamber errors cannot conceal an under-resolved increment. The common
equilibrium account includes all cavity endpoint work.

Gas endpoints, material histories and the complete next solid are returned only
after every source, material, equilibrium, work and cancellation gate passes.
Any refusal leaves accepted state and borrowed inventories unchanged. Empty
cavity input preserves the original advance path. Optional held boundary
pressure remains an explicit additional load; callers must account for its
boundary separately from a cavity's declared exterior pressure.

This is a quasistatic reversible adiabatic gas closure coupled to irreversible
solid compaction. It adds no gas transport, heat exchange, phase split, physical
time integration, unresolved pore-fluid model, evolving porosity, or measured
foam calibration. The `hyper_porous_cavity` regressions specify native pressure
feedback, permanent compaction and elastic unloading, direct variable-caloric
work with large energy datums, independent connected chambers, and atomic refusal.

## Native plane-continuum shear fatigue

`plastic2::fatigue::PlaneFatiguePath::new_with_shear` admits normal-strain
Miner gauges, Fatemi-Socie gauges or both on the original P1/Q1 load path.
Both banks observe the same accepted nonlinear solution and complete stress/
mechanical-strain tensors, including plane-stress thickness contraction,
plane-strain transverse stress and B-bar assumed strain. The original
J2/Voce/Chaboche histories retain sole ownership of constitutive integration.

`PlaneFatemiSocieGauge` binds material-state identity, chemistry/process,
typed absolute temperature, shear modulus and the initial-yield reference.
Source-backed shear fits additionally bind the original card, fixed typed
query and E/nu/initial-yield receipts. Axial fatigue coefficients and evolved
hardening strength cannot substitute for the separate torsional calibration
and its original initial-yield normalization.

The shared `FatemiSocieCycle` owns the full two-dimensional shear path,
maximum-shear-plane convention, whole-cycle normal-stress maximum and bounded
shear-life inversion. `restart_fatemi_socie_cycle` explicitly begins a window
at the accepted state after conditioning, preserving mechanics, lifetime
normal history and completed optional shear Miner totals.
`finish_fatemi_socie_cycle` requires complete stress/strain closure and exposes
the representative cycle's life separately. With explicitly enabled FS/Miner,
it credits exactly one cycle to each fixed plane and remains idempotent.
Repetitions and unresolved multiaxial cycle decompositions are not inferred.

`PlaneShearFatigueBudget` bounds retained paired tensor samples, search planes,
projection work and circle work across all gauges. Ordinary advancement stages
only current mechanics and bounded normal residues; a late failure truncates
shear histories to their saved accepted prefixes. It does not clone accumulated
shear tensors on each step. Whole-path finalization seals both banks atomically.

Focused regressions in `tests/plastic_plane_shear_fatigue.rs` cover native
P1/Q1 reductions and opening-stress effects, Prager shear reversal, conditioning
and nonclosure under isotropic hardening, source/G/yield binding, aggregate caps,
restart equivalence and late append/finalization rollback. They remain synthetic
implementation checks on the existing scalar-factor load path, with no new
independent multiaxis excitation or validated component-life claim.

## Native finite-deformation motion and stress power

`hyper3::response::HyperTetElementMotionResponse` observes supplied displacement
and velocity on the same original P1 tetrahedra and combined material response.
It reports reference velocity gradient, spatial velocity gradient, spatial
symmetric strain rate, Green-Lagrange strain rate, Jacobian and current-volume
rates, and both reference and spatial stress power. The reference contraction
is `V0 P : Fdot`; the spatial contraction is `V sigma : L`, retaining the full
stress tensor even for an explicitly supplied nonobjective material.

The problem and prepared problem expose bulk observations in original element
order. The prepared single-element callback is evaluated once per original
node, with fixed-size scratch and no mesh-sized displacement or velocity copy.
Velocity gradients use relative nodal velocities directly, preserving tiny
rates and exact uniform-velocity cancellation. Observation reuses the original
deformation-gradient/material owner, requests no tangent, advances no history,
and rejects nonfinite inputs/results and cancellation without publishing a
partial response. These are total geometric rates; they do not infer an
elastic/plastic/thermal split or heat production.

Eight independent component oracle cases passed for affine/noncoaxial motion,
time-dependent rotation, rigid spin, tiny rates, translated motion,
nonsymmetric stress and uniform velocity, with additional boundary refusals.
Five public integration tests also exercise original assembled-force power,
prepared/callback equivalence and cancellation; those complete tests remain
unexecuted in the available dependency environment. Component execution does
not claim validation of the complete material/assembly/runtime closure.

### Graded reference density in native plane viscoelastic motion

`hyper2d::visco::dynamics::FiniteViscoPlaneMotion::try_new_with_nodal_density`
retains positive original P1/Q1 nodal densities and prepares the consistent
`integral(rho_0 N_a N_b dV_0)` mass once. P1 uses analytic cubic simplex
moments; mapped Q1 uses 3-by-3 Gauss integration, including its varying
Jacobian. Locally constant fields retain the existing mass arithmetic exactly.
Reported cell means use original-volume weighting and never replace the
graded operator. Original thickness and the full nodal field stay immutable
through deformation, viscoelastic relaxation and material-reference changes.

Native momentum, support reactions, acceleration loads and kinetic work all
consume this same retained mass. Positive-density/geometry admission, bounded
assembly storage and cancellation remain enforced. G0/G1/G3 regressions in
`inertia/plane/tests.rs` and `tests/visco_plane_stress/dynamics.rs` cover
independent simplex/tensor polynomial moments, shifted mass centroid and
forcing, actual plane motion, reference-change preservation and constant-field
bit parity. This is authored reference inventory, with no inferred material
calibration, mass transfer, transverse inertia or current-density resampling.

## Dispersive structural harmonic frequency derivatives

`ViscoHarmonicOperator::with_frequency_derivative` retains the existing physical
factor and prepares its material frequency direction once. Original tensor-Prony
laws supply `dC'/domega` and the directly evaluated `d(C''/omega)/domega`; the
shared const-generic tensor assembly retains every original quadrature point,
Mandel cross term, B-bar column, physical measure and free-DOF map. Symmetric
entries are mirrored from one contraction. The primal spectral path retains
exact reuse of the original equilibrium stiffness when all storage tensors
equal their equilibrium values. Reference geometry, supports and mass are fixed.

The resulting `ViscoHarmonicFrequencyOperator` reuses the same harmonic factor
for each source. It passes `M_omega=0`, assembled `K'_omega` and effective-viscosity
slope, and the explicit frequency direction one to the shared total-parameter
solver. Thus `A_omega=K'_omega-2*omega*M-i*K''_omega`, including the actual
constitutive dispersion. Authored complex force slopes and complete attachment
and impedance slopes are optional; an absent port direction explicitly fixes
B and Z, while retaining their kinematic frequency action. No sampled material
slope, extra factorization or material law reconstruction is introduced.

The response includes original motion plus displacement/velocity slopes and
separately differentiated source work, intrinsic material heat and boundary
work. Its original displacement must pass the existing full-cycle material
strain admission at every retained quadrature point, including compatible
strain for B-bar. A derivative is not itself an additional strain cycle or a
certificate that neighboring frequencies satisfy the source domain. Refusal or
cancellation returns no partial response; the prepared factor remains reusable.

Derivative preparation admits exactly `24*n*n` additional retained dense source
bytes before assembly: zero mass slope, effective-viscosity slope and storage
slope. Original sparse assembly and constitutive temporaries keep their own
ownership. Per-solve tangent scratch/work uses `HarmonicParameterControls`.
The primal structural preflight now matches the shared dense owner's
`144*n*n+128*n` bytes plus its existing `24*n*n` dense source matrices.

`tests/visco_harmonic_frequency.rs` defines scalar closed-form dispersion and
heat slopes with complex source variation, full anisotropic/coupled spatial
derivatives with consistent mass, independently perturbed boundary/source
responses, and budget/cancellation/full-cycle refusals. The existing primal
one-coordinate resource boundary is updated from 280 to 296 bytes. These native
regressions remain unexecuted pending the real dependency graph's build resource
blocker; no measured material, modal-truncation or throughput claim follows.

## Thermal measured-curve forming and retained shape

The existing native finite thermoplastic tetrahedron and P1/Q1 plane-stress
sheet accept affine-temperature tabulated Kirchhoff hardening through the
shared `FiniteThermalJ2Material` family. All original plastic metric and
accumulation history, temperature, coefficient work, entropy, mechanical internal
energy and signed caloric transfer remain with those owners. The sheet still
solves transverse equilibrium and thickness locally; heating does not prescribe
a global scale or reset an existing plastic state.

Two structural regressions apply an unchanged nominal dead load through loading,
heating, cooling and release. Independent free-sided Hencky equations distinguish
nominal from Kirchhoff stress and predict further yielding on heating, a
hardening-knot crossing, lateral motion, sheet thickness and stress-free
permanent shape. They exercise one tetrahedron and both triangular and
quadrilateral sheets. Unsupported temperature or measured strain support refuses
without changing accepted history, geometry or energy.

These are synthetic constitutive/structural fixtures within the existing
displacement-solid and membrane domains, not alloy formability calibration.
Scoped Rust syntax, whitespace and independent scalar arithmetic were checked;
native Cargo regression execution remains unverified in the available shared
build environment.

## Physical-time tabulated malleability in fibers and frames

`FiberMaterial::tabulated_overstress` and its source-bound constructor retain the native free-sided tabulated Perzyna law, independent Prager history and actual physical duration. The source-bound route admits density from the same original material bundle. Immutable measured knots and source evidence are shared when the duration changes; the original elastic, hardening, kinetic and source calibration cannot change underneath accepted history.

The existing physical-time frame owner now admits this fiber law for quasistatic load/hold intervals and backward-Euler dynamics. Every equilibrium trial reuses the same accepted axial history and interval; a successful interval publishes the complete frame atomically. Loaded members can continue flowing, held members can relax stress, and bending members retain permanent curvature on release. Unsupported measured-strain continuation, source substitution, inconsistent durations and cancellation refuse without modifying accepted history. Clockless static, modal, Newmark and standalone pushover routes reject the timed law.

`FrameInelasticEnergy` preserves the existing creep and yield-offset fields and adds a separate physical overstress channel. Native hardening remains recoverable storage. Physical dissipation is the sum of yield, creep and overstress work exactly once; material and inertial integration losses remain numerical and never become caloric heat. Existing creep conversions retain their previous values and semantics.

Four structural regressions cover original-source force loading across measured knots, variable-duration hold relaxation and release; moment-driven bending and retained curvature; power-law inertial motion with both lumped and consistent reference mass; and clock, source, domain and cancellation refusals with usable retained history. The fixtures use independent scalar constitutive and momentum/work equations. Rust syntax, formatting, whitespace and those independent arithmetic checks passed. Native Cargo regression execution remains unverified in the resource-constrained shared environment. These tests establish model behavior only; they do not supply empirical alloy or forming calibration.


## Current-area pressure in retained thermal porous compaction

FiniteThermalPorousTetSolid::advance_with_pressure composes the original
current-area surface force and analytic follower tangent with the SAME thermal
porous advance and acceptance path. Supplied endpoint temperatures, old-temperature
history admission, original reference mass and optional backward-Euler inertia
remain with their original owners. Uniform pressure supplies -p delta V; other
admitted fields retain the surface owner's declared-path work. The prepared
pressure surface must bind the exact original mesh.

The step retains a separate FinitePorousPressureEnergy account. Material endpoint
defect, yield loss and reciprocal thermal exchange remain unchanged. Complete
work includes pressure path work and its geometric correction, plus kinetic
change and separately admitted numerical kinetic loss when present. Only the
original physical caloric heat enters a thermal owner; pressure work and
numerical defects are never deposited again as heat. The material-work allowance
bounds the sum of absolute material and pressure corrections without cancellation.

Three focused manufactured regressions specify thermal pressure continuation
and elastic release against an independent hydrostatic return, free inertial
compression against the consistent-mass formula, and translated pressure-patch
work with graded mass and support reactions. They also specify caloric identity,
late work/kinetic refusal, exact source/mesh admission, zero-pressure compatibility
and retry/cancellation without changing the accepted owner. The existing dry
porous, P1 locking, source calibration and prescribed-field limitations remain.

Native regression execution remains unverified on the exhausted shared build host.


### Retained porous dynamics with independently sealed cavities

`hyper3::plastic::porous::dynamics::FinitePorousTetMotion::advance_with_cavities`
advances one physical backward-Euler interval through the original native
material, cavity and consistent-mass owners. Current gas pressure, sparse
follower geometry, independent acoustic rank-one factors, porous tangents and
the full original `M/dt²` enter the same Newton solve. Original reference mass
survives compaction; current density remains `rho_ref/J`. Free bodies require
no artificial restraints, and prescribed-motion reactions include the full
mass action and cavity forces.

Each cavity retains its original undeformed-volume gas reference across
increments. Ideal and Peng--Robinson source owners retain their own calorics,
homogeneous-phase admission and acoustic derivatives. Accepted gas parcels are
output evidence, never replacement reference sources. The homogeneous
isentropic gas closure resolves neither gas momentum nor pressure waves and
does not transfer solid plastic heat into gas energy.

`FinitePorousCavitiesMotionStep` retains the material/momentum step, each
independent gas endpoint, original-source gas/exterior work, and the combined
signed endpoint-integration defect. Direct caloric and compact geometric
increments preserve small interval work without subtracting absolute gas
energies or reference-relative cavity potentials. Its whole-system identity is

`dead + support + held pressure + sum(exterior work)
 = delta(solid storage) + sum(delta(gas storage)) + plastic loss
   + delta(kinetic) + numerical kinetic loss
   + external_work_defect + step.balance_error`.

Gas storage, physical plastic loss, signed material/pressure/cavity endpoint
defects and numerical kinetic loss remain separate. The material work budget
bounds the sum of absolute constituent endpoint defects; the kinetic-loss
budget remains independent. Failure of any source, history, momentum, work or
cancellation gate preserves the original motion and all borrowed gas references.
An empty cavity list preserves the original pressure/dead-load dynamic path.
Gravity coupling remains outside this finite-inertia API.

The `hyper_porous_cavity/dynamics.rs` regressions specify a fully free hollow
shell with independent consistent-mass loads, analytic gas/solid compaction
and kinetic work; and retained release of two mechanically connected cavities
with different source calorics, independent support reactions, a large energy
datum, unchanged inventories, source order, empty-list parity and late kinetic
refusal/retry. Sources were parsed/formatted and the independent shell force,
mass and kinetic formulas were checked algebraically. Native Cargo compilation
and regression execution remain unverified under the shared resource blocker.

## Complete hereditary harmonic frequency tangents and driven linewidths

`ViscoJointHarmonicFactor::solve_with_frequency_derivative` differentiates the
original retained polynomial `P(s)=P0+s*P1+s²*P2` at `s=-i*omega`. Its additional
right-hand side is `f_omega+i*P1*x+2*omega*P2*x`; the existing complex LU solves
it and up to three corrections. Every original displacement, Mandel Maxwell
strain, independent load state and algebraic reaction remains present. Material
relaxation and Kelvin viscosity are differentiated through their original
equations. No effective damping is frozen, no neighbouring solves approximate
the slope, and no extra matrix is factored or inverted.

The response includes original state, displacement and physical velocity slopes,
with `v_omega=-i*(u+omega*u_omega)`. Original K/C/M, pHS, descriptor and continuous
state-space load owners differentiate their own displacement, flow and effort,
including derivative inertance and zero-storage constraint reactions. Signed
per-port and total absorbed-work slopes retain both product-rule terms. They
remain exchange rates; no boundary heat or storage is inferred for an
unqualified state-space realization. Primal energy accounts remain independently
computed by their existing owners; this API does not publish energy derivatives.
The existing harmonic owner's binary-normalized DD contraction retains reactive
and mutual cancellation before rounding work. The total is contracted from the
original products independently of the rounded per-port outputs. Every requested
per-port value and the total must be representable; either can explicitly refuse.

The complete original tangent equations are checked in mechanical, memory and
load families with coefficient-scaled residuals. Every original algebraic
constraint and each unqualified state-space row also passes its own gate. The
primal and tangent share the preparation's operator-evaluation limit. The new
`ViscoJointHarmonicFrequencyBudget` admits an additional `384*n+256*p` numerical
bytes and `1024*(n*n+n*p+p*p)+512*(n+p)` work before solving. Here `n` includes
all retained states and `p` all load ports; factor and primal-response resources
remain separately admitted. The extra storage is linear in coordinates, and
each tangent solve/audit is quadratic; this does not claim sparse large-model
performance. Polling and refusal leave source histories and the factor unchanged.

All original material and load bands still apply at each sampled frequency.
The full-cycle source strain check admits the primal physical response; the
tangent is not interpreted as another strain cycle. Geometry, supports,
material coefficients, boundary realization and attachment map are fixed.
An optional authored complex force slope is retained. This provides a local
frequency derivative, not neighbouring source-domain admission or a forward
error certificate near an ill-conditioned joint resonance.

`displacement_transfer` and `velocity_transfer` return the existing
`HarmonicTransferSample`. The shared `fs-modal::harmonic::resonance` owner remains
the only peak/notch and half-power refinement algorithm. Complex or dispersive
receivers use its existing complex projection functions on the exposed original
values and slopes. A bare attachment pole can be a regular joint response and
an exact forced-response notch; this path needs no singular standalone impedance.
An actual joint singularity still refuses. Searches must isolate the desired
continuous extremum/crossings and preserve the same physical source and receiver;
the measured response linewidth is not an inferred modal damping ratio.

Focused regressions extend `visco_joint_harmonic`, `visco_harmonic_frequency`,
`visco_descriptor_loads` and `visco_state_space_loads`. Their independent oracles
cover anisotropic tensor cross terms and original memory slopes, complex force
variation, mutual attachment work, a bare-load antiresonance, original constraint
reactions, a dark state and direct boundary inertance. A rational hereditary
oscillator with a physical K/C/M load exercises the existing peak and half-power
owner against independent stationary/half-level polynomials and agrees with the
separate `ViscoHarmonicFrequencyOperator` at its peak and both edges. Resource,
shared-evaluation, strain-domain and cancellation refusals preserve factor reuse.
These native tests remain unexecuted under the real dependency build's resource
blocker; scoped syntax and whitespace checks do not replace their execution.

## Temperature-dependent tabulated malleability in timed frames

`FiberMaterial::thermal_tabulated` admits an owned original-source thermal interval and derives reference density from its caloric inventory. `FiberLaw::ThermalTabulated` retains the native axial plastic state and accepted temperature. The existing physical-time quasistatic and backward-Euler frame owners consume it through the same force-based equilibrium, history commit and work audit. Rebinding admits new temperatures and duration only when the original calibration and accepted-temperature continuity agree; Newton trials share immutable preparation and never query sources or advance extra material time.

The thermal adapter reads old Helmholtz storage directly at the accepted old temperature, retains native endpoint storage, and reports yield and viscous work separately. Exact mechanical-entropy storage and signed reversible heat remain explicit. Native material and inertial numerical losses are never deposited as heat. The scalar history must agree with the native tensor/thermal increment before the frame is returned. Unsupported source/domain changes, clockless routes, incompatible intervals and cancellation leave the accepted frame unchanged.

`Frame2d::thermal_tabulated_fibers` exposes a bounded, cancellable inventory in the original integration order: member/section/fiber identity, reference volume, member reference mass, accepted axial strain and borrowed material/history. It requires complete thermal-fiber coverage and performs no equilibrium solve or constitutive advance. This supplies a direct correspondence to existing caloric cells without reconstructing a second frame.

Three G0/G1 regression sources exercise heated force-controlled extension and permanent shape, unequal-fiber-temperature bending and residual curvature, and heated power-law dynamics with both reference-mass formulations followed by held-strain relaxation. Independent scalar momentum, storage, entropy, heat and work calculations passed, as did scoped Rust syntax/format, whitespace and API checks. Native Cargo execution remains unverified because DSR/RCH are unavailable and shared build storage is full. The existing linear geometry, small-strain source domains, stationary-support and reference-mass assumptions remain binding.

## Full-stick vector adhesion on native P1 faces

`AdhesiveP1Patch3::with_full_stick_tractions` explicitly projects the original
mixed-mode JKR normal field and both full-stick tangential fields onto the same
exterior P1 traces. Local components follow the declared x tangent, y tangent
and inward normal. Exact sector force and first moments supply every nodal
component; singular edge tractions are integrated rather than sampled.
Normal-only projection remains the default.

The same explicit vector path accepts fresh elastic torsion profiles from the
original `SourceJkrTorsionContact`. Its azimuthal traction projects through the
same affine barycentric functions; original signed torque acts about the
declared inward normal. Projection gates compare the actual represented nodal
moment with this original nonzero couple. No force-normalized application
point, plastic-history field or additional mechanical solver is introduced.

`AdhesiveTractionProfile::jkr_torsion_equilibrium` supplies this same native P1
path with the original source's changing JKR area under prescribed approach
and twist. Normal pressure and resisting torque use that single evaluated
radius and its remaining normal opening, so peeling changes normal load and
twisting couple together. The torsional-peeling regression compares actual
native displacements, represented nodal work and independently reconstructed
spatial work against separately assembled analytic nodal forces at both source
endpoints. It retains the original source admission and prescribed-field scope.

The native transient owner retains the selected vector endpoint atomically,
including receiver refusal and supported-node motion. Its independent spatial
work reconstructs all three affine displacement components from the original
mesh and compares them with the represented nodal work. Existing support
reaction work stays separate. Tangential radial moments are analytic and
reuse bounded ring scratch; adding shear introduces no pressure quadrature.
The reference geometry and original source/full-stick applicability remain
fixed. Partial slip, follower contact and a new finite-body contact solve are
outside this explicit prescribed-traction extension.

The focused vector regressions specify independent closed-form nodal forces,
rotated-frame wrenches, actual Newmark motion under changing shear direction,
atomic receiver retry and tangential adhesion work on prescribed supports.
Elastic torsion cases independently recover four edge forces of magnitude
`|M|/(4L)`, rotated source couples, native Newmark motion against independently
assembled loads, and signed `M*delta_phi` work on actual P1 support rotation.

## Full original-quadrature harmonic material fields

`ViscoElements::prepare_harmonic_material_observer` caches the original free/support
maps and each element's original tensor-Prony harmonic spectrum. Repeated forcing
cases reuse this immutable observer without another mechanical factorization.
It accepts the original free-coordinate displacement phasor from either harmonic
solver and returns every original constitutive quadrature point: six-component
strain and stress, equilibrium/Maxwell storage, separate Kelvin/Maxwell heat and
signed stress-work closure. The reference quadrature volumes independently
integrate total stored energy, physical material heat and stress work.

Original compatible strain is retained alongside B-bar constitutive strain, and
both independently satisfy the full-cycle source limit. Original fixed supports
are zero, matching the harmonic mechanical owners. No point is subsampled and no
strain is reconstructed from a modal loss ratio. Preparation admits point count,
heap numerical payload and conservative work before allocating the cache; each
observation reuses one local nodal buffer and returns one complete field. Budget,
cancellation, band and strain refusals publish no partial field or material history.
This observes a supplied phasor, without certifying an external force equation or
adding startup transients, thermal evolution or a spatial error estimate.

Three native regressions cover analytic tetrahedral stress and spring energy,
agreement with the independent mechanical solver's input/heat, deterministic
cache reuse, source/budget/cancellation refusals, and every original compatible
and B-bar strain in a nonuniform HEX8 field. These integration regressions remain
unexecuted because the original runtime dependency build is resource-blocked;
the material observer's production arithmetic was exercised separately.

`ViscoHarmonicMaterialObserver::observe_window` uses the same original spatial
maps and reference quadrature to observe a finite physical interval. It retains
the actual phase, every compatible/B-bar strain and every Kelvin/Maxwell source
channel. One prepared material window per element supplies allocation-free point
evaluation. The returned field separately integrates start/end recoverable
storage, direct stored-energy change, Kelvin heat, Maxwell heat and stress work;
the signed work-minus-storage-minus-heat residual is never forced to zero. Direct
pointwise storage increments are summed rather than subtracting rounded total
endpoint stores. Uniform time intervals are required: a triangular detector
sensitivity is not unweighted physical heat or endpoint energy conservation.

The original supplied payload/work caps also apply to this operation. Its actual
returned point type is sized explicitly, alongside the existing spectrum/maps
and one local nodal buffer. Before any result allocation the window path admits
the conservative work bound `262144*sum(points*(1+branches)) +
8192*sum(1+branches)`. A cap that admits cycle means but refuses interval work
leaves the original mean observer usable. No extra constitutive history or
thermal solve is created, and cancellation/domain refusal publishes no field.
`visco_harmonic_window` checks affine tetrahedral and nonuniform HEX8 B-bar fields
against independent real-field Maxwell/Simpson calculations, physical full-cycle
and split-interval limits, and budget/profile/strain/cancellation refusals.
Those four new solid integration tests are present but remain unexecuted under
the shared runtime build constraint. Their material owner passed seven native
selected-source interval tests and eighteen existing numerical regressions.

### Coherent multitone material heating on original solid fields

`ViscoElements::prepare_harmonic_material_spectrum_observer` connects the same
original free/support maps and constitutive quadrature to the shared
`TensorPronyHarmonicSpectrum` material owner. It accepts up to 64 positive
source-admitted angular frequencies; `observe_window` consumes one complete
free-coordinate displacement vector per authored line. This applies to every
element family already served by `ViscoElements`, with every original compatible
strain independently retained and admitted alongside B-bar constitutive strain.

Equal-frequency fields combine before source admission and energy contractions.
Distinct-frequency sum and beat terms remain in actual finite-interval Kelvin
heat, Maxwell heat, recoverable storage and stress work. No common period or
independent-phase assumption is imposed. Each point exposes the original tone
order, reference-volume weight and physical material account. Volume integration
uses compensated totals and direct pointwise storage increments; stress work and
positive physical heat are computed independently of the signed balance defect.

Preparation admits authored line count, complete point count, numerical payload
and conservative quadratic-in-lines work before allocating spectra or fields.
One material temporal cache per element and one gathered nodal buffer serve all
its points; original frequency spectra are reused across forcing cases/windows.
Cancellation polls between material channels/branches as well as original points.
Refusal returns no partial field and changes no source or accepted history.

The material owner uses a conservative sum of distinct-line peak strain norms to
admit the complete waveform. It can refuse physically valid waveforms whose
separate peaks cannot coincide. This is settled linear viscoelastic response at
fixed material conditions, without startup memory, a DC preload, temperature
feedback, spatial convergence or an external-force solution certificate.
The G1/G3 tests in `visco_harmonic_window` compare affine tetrahedral and nonuniform
HEX8 B-bar fields with independent real Maxwell fields/Simpson quadrature, exercise
clock translation and split intervals, reject incompatible strain/budget/profile
inputs, preserve coherent duplicate cancellation and one-tone limits, and feed
the existing dispersive mechanical solutions directly into the shared observer.

### Thermal retained porous solids with independently sealed cavities

`hyper3::plastic::porous::thermal::FiniteThermalPorousTetSolid::advance_with_cavities`
closes the prescribed-temperature thermal-material/sealed-gas composition. It
freezes only the supplied endpoint solid temperatures, admits accepted plastic
history at its original temperature, and sends the original thermal porous
tangent, each cavity's follower geometry/acoustic rank-one tangent, optional held
current-area pressure and optional original-reference consistent inertia through
one existing native Newton solve. Cavity volume comes from total solid
displacement, not plastic volume or the thermally reduced material deformation.
Heating can therefore compact the material while fixed cavity geometry leaves
the independent gas states unchanged.

The returned `FiniteThermalPorousCavitiesStep` owns the complete original thermal
step, accepted gases in input-source order, the existing
`FinitePorousCavityEnergy` receipts and a complete signed external-work defect.
The dry and thermal paths share one direct original-source interval-work owner.
Original gas references remain at undeformed cavity volume across retained
increments; accepted deformed gas parcels are evidence, never replacement source
references. Stable direct gas caloric increments and compact geometric work
remain valid with large caloric energy datums and small retained increments.

The complete accepted identity is:

`dead + support + held-pressure work + sum(exterior work) =
material internal-energy change + original material caloric heat +
sum(gas stored-energy change) + kinetic change + numerical kinetic loss +
external-work defect + balance error`.

Original material caloric heat already contains physical yield loss and signed
reciprocal thermal exchange. Neither gas storage/work, pressure quadrature error,
constitutive endpoint defect nor numerical kinetic loss enters caloric heat.
The material work allowance bounds the sum of absolute material, held-pressure
and individual cavity corrections; independent sources cannot cancel that
acceptance budget. Backward-Euler numerical kinetic loss retains its own budget.
Empty cavity lists follow the original thermal pressure/dead-load branches.
Failures, including final work gates and cancellation, publish no material,
temperature, gas or momentum history.

Focused regressions are registered under
`tests/hyper_porous_thermal/cavity.rs`. A completely free heated hollow shell
combines analytic hydrostatic thermal compaction, independent consistent P1
mass forcing, current-area exterior pressure and a sealed gas; it checks real
Newton motion, source gas storage, material caloric heat and the full kinetic
balance. A connected two-cavity case uses different original heat-capacity
laws and reversed source order, heats at fixed total volume, releases
elastically, and checks direct interval energy with a 1e24 J/kg datum,
original-reference refusal, empty-path equivalence and atomic work/cancellation
refusals. Source parsed/formatted with pinned rustfmt; native tests were not run
in the resource-exhausted session.

This composition prescribes the solid temperature and keeps each gas
isentropic. It does not add gas-wall heat exchange, gas momentum/acoustic
propagation, pore-fluid transport, contact or evolving conductance. Existing
P1 locking, thermal-law domain, associated compaction and homogeneous gas-phase
admission limits remain explicit.

## Current-area pressure in finite thermal plastic shaping

`hyper3::plastic::thermal::FiniteThermoplasticTetSolid::advance_with_pressure`
uses the original closed tetrahedral pressure surface in every nonlinear trial.
Its follower force and tangent act on the current deformed area while the native
thermal J2/Voce/tabulated, creep or overstress material retains its complete
finite plastic history. The same prepared reference geometry and optional
consistent reference inertia remain the mechanical owners; a foreign pressure
surface is refused before solving.

The returned pressure account retains original path work, endpoint work and
their signed geometric correction. Uniform held pressure uses exact
`-p delta(V)`; other pressure fields retain their original potential or
declared displacement-path convention. Combined balance includes actual pressure
work and its correction once. The sum of absolute material and pressure
integration defects must fit the caller's work budget, so opposing defects
cannot hide each other. Physical plastic/creep and reversible thermal heat are
unchanged; neither pressure quadrature correction nor numerical inertia heats
the material. The existing unpressured entry point preserves its solve and
arithmetic path.

Two G1/G3 regression definitions use independent finite coaxial equations for
pressure-driven loading, heating and release with retained shape, pressure
volume work and reciprocal heat, and for motion with an independent P1 kinetic
integral, original mass and transactional refusal/retry. Independent fixture
arithmetic was checked in JavaScript. The execution host became unavailable;
Rust formatting, compilation and native regression execution were not performed.
This remains the native bounded P1 finite-strain model with held interval
pressure, explicit reference transport geometry and its original locking and
constitutive limits; contact tooling, mesh adaptation and experimental forming
validation are separate capabilities.


### Gas isentropes anchored at accepted deformed cavity geometry

`hyper3::pressure::cavity::CavityGasEpoch<G>` admits a physical gas parcel at
the actual accepted cavity geometry. This closes the native seam needed after
source-owned heat exchange at fixed current volume: the heated parcel starts
the next mechanical isentrope directly, without a fictitious entropy solve back
to undeformed solid volume. The caller retains the physical heat receipt and
owns coupled time integration and energy acceptance.

`try_new(Arc<ClosedTetCavity>, accepted_displacement_m, gas, cx)` admits the
complete solid mesh through the existing geometry owner and requires bit-exact
agreement between gas volume and represented accepted cavity volume. The epoch
retains the exact prepared-cavity identity, compact accepted positions and
volume in one shared immutable anchor. Newton clones share the anchor and
prepared cavity; they do not clone the full mesh or retain a full displacement
copy. Gas inventory, source model, auxiliary query axes and homogeneous-phase
admission remain with the original ideal or Peng--Robinson parcel.

A sealed reference-geometry hook keeps ordinary parcels' strict undeformed-volume
requirement. Only explicitly constructed epochs use an accepted deformed anchor.
Their force, acoustic compressibility, geometry tangent and direct stored-energy
increments delegate to the existing native owners. Exterior pressure potential
uses the original compensated compact volume increment from this physical
anchor, so the complete conservative potential remains
`U(current)-U(anchor) + p_exterior*(V(current)-V(anchor))`. Physical heat already
supplied to the anchored parcel is not replayed as mechanical work.

An accepted load returns an output-only epoch sharing the original anchor.
It may expose or yield its physical parcel through `gas()`/`into_gas()`, but
cannot silently replace a reference, even if its volume happens to equal the
anchor volume. A different prepared-cavity identity also refuses. To start a
new reference after physical heating, the caller explicitly constructs a fresh
epoch at the accepted displacement. Continuing the same isentrope instead uses
the original admitted reference. This explicit alternative applies to the
existing generic static, retained porous, thermal and consistent-inertia
cavity consumers without changing their work or history acceptance owners.

Focused G1/G3/G4 regressions in `tests/hyper3_cavity/epoch.rs` specify a heated
variable-Cv ideal gas at a deformed anchor, an admissible next isentrope whose
fictitious undeformed rebase would leave its temperature bracket, independent
force/acoustic/work derivatives, source and inventory preservation, exact
reference identity, output-only/fresh-binding behavior, complete-mesh refusal
and cancellation. The existing `hyper3_cavity/real_gas.rs` adds a heated dense
gas case using its independent PR Helmholtz/entropy oracle, including both
density-dependent caloric storage and the true acoustic bulk modulus. Both
gas-law tests retain a 1e24 J/kg energy datum to exercise direct increments.

The execution host was unavailable during implementation. A bounded pinned
rustfmt stdin attempt could not start; these new sources and tests have source
review and independent scalar-oracle checks, but no new formatting, compilation
or native-runtime pass is claimed. This seam itself adds no heat-transfer law,
simultaneous heat/motion integrator, gas momentum, phase split or contact model.

## Thermal cyclic fiber malleability

The native physical-time frame owner accepts source-bound thermal Chaboche
intervals with one through eight independently retained backstress branches.
Runtime dispatch preserves the exact native law, branch ordering and complete
history. Clones share material preparation; mechanical trials perform no source
query, branch aggregation or history reset.

Quasistatic and backward-Euler steps retain the actual physical interval and
original reference mass. Temperature changes transport each branch through
the material owner's strain-like history before the existing axial return.
Rebinding preserves original source/caloric identities, branch count and
ordering, and accepted temperature. Untimed structural routes refuse these
interval-bearing laws.

Yield-offset work, physical viscous overstress and nonlinear branch recovery
are reported separately. Their sum enters physical material heat once;
reversible thermal exchange remains signed, and material/inertial numerical
loss stays outside heat. The audit compares every independent plastic and
backstress component with the native tensor account in declared stress units.

The shared thermal inventory covers heterogeneous tabulated and cyclic fibers
with original quadrature volumes and member mass. Missing or mismatched
families refuse without omitting physical inventory. Three focused regression
definitions cover thermal reversal, unequal heating with eight branches, and
timed motion with both native mass formulations. Independent fixture arithmetic
and bounded API/conservation review were completed; native compilation,
execution and formatting remain unverified because the execution environment
became unavailable.

The native linear-geometry, small-strain, stationary-support and fixed-reference
mass limits remain in force. This supplies no finite cyclic-forming or
experimental calibration claim.


The admitted reference epoch additionally exposes `with_gas_at_anchor(gas, cx)`
for bounded nonlinear heating probes at the same accepted geometry. It requires
reference status and exact accepted-volume agreement, retains the same shared
anchor/prepared-cavity identity, and performs no complete-mesh scan, compact
geometry copy or thermodynamic reconstruction. Output-only epochs still refuse.
The original gas owner supplies the physical replacement and its heat/source
account; this operation infers neither heat nor inventory continuity. A changed
accepted geometry still requires `try_new`. Focused existing epoch assertions
cover heated same-anchor reuse, original-reference preservation, wrong-volume
and output-only refusal, and cancellation. The host remains unavailable; this
increment has source review but no new formatting or native-runtime claim.


### Sustained rigid-tool forming on retained thermal plastic history

With `contact`, `FiniteThermoplasticTetSolid::advance_with_plane_contact`
consumes the original `fs-contact::sustained` force-space complementarity solve.
Each force response calls the same native finite thermoplastic equilibrium from
accepted deformation, plastic metric, source family, temperature and optional
consistent reference inertia. Fixed-normal rigid planes translate over the
actual interval. There is no penalty stiffness, surrogate compliance, impulse
to force conversion, friction inference or competing constitutive update.

At most 21 explicit original node/plane rows enter the native solve. Every
original vertex is admitted against every plane at both endpoints, preventing
an omitted candidate from permitting penetration. Affine plane and nodal paths
then retain nonnegative gaps to their compensated arithmetic precision. The
original interval tetrahedron query independently excludes collapse throughout
the declared chord under an explicit work budget. This is a discrete geometric
scope, not exact-flow CCD, self-contact, curved tooling or automatic discovery.

Original dead-load work, body-contact work, prescribed-tool work and tool
reactions remain separate. Signed normal endpoint-work defects are bounded by
the sum of absolute row defects; they are not heat. The original material heat,
reference inertia, numerical kinetic loss and material endpoint defect retain
their owners. Independent load decomposition and complete energy gates apply
before publication. Physical arrival, unloading and material kinetics require
interval refinement. All accepted histories remain unchanged on refusal.

Three synthetic regression definitions cover heated compression/dwell/release,
rotated geometry and work ports, original inertia and normal-work refusal,
omitted-row penetration and cancellation. Rust parsing/formatting passed;
independent scalar fixture calculations do not substitute for native compilation
or test execution, which remain unverified in this environment.

## Temperature-dependent relaxation of original adhesive joints

`cohesive::thermal_joint::relaxation::SourceRelaxingJoint` now exposes
`advance_with_interface_temperature[_controlled]`. Its explicit successor is
guarded by `ResolvedViscoelasticBond::require_temperature_successor`: original
ordered surfaces, card and selected coefficient claims, indexed Maxwell
branches, thickness, use bounds and every nonthermal query coordinate remain
fixed. The ordinary fixed-query advance retains its original semantics.

The existing material temperature transfer retains physical spring strains,
total jump, elapsed adhesive time and accumulated mechanical/heat history. At
held jump, changing spring moduli supplies signed parameter work; changing only
relaxation times supplies exactly zero. The same native series compatibility
solve then advances the original members and retained adhesive memory under
prescribed endpoint force or extension. It borrows the transferred bond rather
than cloning all member histories into a staging assembly.

`TemperatureReport` separates interface parameter work and independently
evaluated held-jump storage change from mechanical motion, physical relaxation
heat and numerical endpoint-work corrections. The complete old-to-new balance
uses constitutive increments and compensated sums, not subtraction of large
whole-joint stored-energy inventories. Source, memory, geometry, displacement
and clock publish together only after local/global energy and final
cancellation gates; failed source admission or mechanics leaves the prior joint
retryable.

The existing `tests/relaxing_adhesive_joint.rs` includes G1 analytical
heating/cooling references for force and extension control, retained branch
memory, a kinetics-only heated series relaxation/time-refinement reference,
and G3/G4 source/nonthermal-coordinate refusal and identical retry after late
cancellation. Coefficients are manufactured inputs, not calibrated adhesive
specimens. This is a first-order staggered prescribed-temperature update of an
intact, massless, small-strain axial assembly; combined temperature/load paths
require refinement. It introduces no adhesive thermal expansion, caloric or
entropy model, thermal transport, automatic heat deposition, damage or healing.

### Complete finite creep and recovery on the original tetrahedral owner

`FinitePlasticTetSolid` and `FinitePlasticTetMotion` consume the material owner's `FiniteInelasticMaterial` protocol. The original hyper3 geometry, sparse Newton/globalization and consistent reference mass remain the sole mechanical owners. Native/source J2 callers retain their previous state through the blanket adapter; Norton, independently recovering Norton, Perzyna and combined Norton/J2 retain their complete typed histories. Every Newton and line-search trial borrows accepted history. Kinetic elements must share one positive physical interval, exposed by `duration_s()`, and motion requires that interval to equal its explicit timestep. Rejection leaves displacement, velocity, irreversible metric/path and resistance untouched.

The spatial account independently retains stored energy, integrated flow potential, native plastic/viscous and combined-creep heat channels, and signed endpoint quadrature defect. Newton's potential is checked against storage plus flow potential, not physical dissipation. `inelastic_dissipation_j()` sums the original physical channels; numerical endpoint/kinetic losses are excluded. The native sole-Norton convention still uses the plastic-named metric/alpha/heat slots for creep, while combined laws keep plastic and creep accumulations separate. This adds no recovery energy, rupture qualification or P1 locking-accuracy claim.

`tests/hyper_plastic_creep_recovery.rs` authors loaded free-sided recovery/dwell/reload against independent finite Hencky and quadratic kinetic equations, nonvirgin stress-free translation with consistent inertia, simultaneous creep/yielding with separate paths and heat, and clock refusals. Coefficients are manufactured. Native execution of this addition is not claimed; dependency compilation was killed before fs-solid was reached.

## Complete finite thermal recovery histories in structural consumers

The original `FiniteThermoplasticTetSolid` and step, pressure, sustained-plane
contact and impulse paths accept `FiniteThermalInelasticMaterial` with its full
associated state. Every mechanical/thermal trial freezes the same accepted
history and actual interval. Permanent creep, plastic hardening, independent
recovery resistance, deformation, temperature and optional momentum advance
atomically. The original Newton/tangent, reference-volume, consistent-mass and
work gates remain authoritative; there is no second structural solver.

Caloric heat is physical native inelastic loss plus signed reversible transfer.
Flow potential, signed endpoint quadrature defects and numerical kinetic loss
remain separate. Sole Norton retains the original plastic-named loss field;
combined families retain distinct plastic and creep channels. Zero-duration
impulses preserve the complete material memory without recovery or flow.

Synthetic free-coordinate thermal recovery regressions compare Hencky/Arrhenius
loading, hot stress-free dwell, reload and full energy accounts to independent
scalar equations; graded-mass impulse continuation and simultaneous yielding
retain their original owners. No measured calibration, life prediction, locking-
free P1 behavior or automatic remeshing/state transfer is claimed.

### Native source-bound thermal Hex8 continuation

The small-strain Standard/B-bar Hex8 owner now exposes
`plastic3::hex::thermal::ThermalHexProblem`. It binds one original
`ThermalSectionSourceFamily` per reference cell and eight independent specific
enthalpies/native histories per cell. Each physical interval resolves all eight
endpoint laws once, validates accepted strain/history at its **old** temperature,
and freezes those histories throughout the existing Newton/backtracking solve.
Original families, geometry, connectivity and formulation must remain bound
across continuation; no material or history remapping is inferred.

The static and consistent-mass backward-Euler entries use the same assumed
strain, projected tangent, quadrature and history owner. Timed source laws
receive exactly the supplied positive physical duration. Motion density comes
from each original source family; it is never recomputed from endpoint
temperature. Accepted stresses retain the old-history constitutive endpoint.

Before returning a candidate, the original source owner independently audits
every point and must reproduce the native history/stress. Reference-volume
point and whole-solid receipts retain physical heat (with creep as a
subaccount), reversible caloric transfer, free-energy change (with coefficient
work as a subaccount), change in temperature times mechanical entropy, and
material numerical loss. Actual applied/support work must close both mechanical
equilibrium and the complete energy balance. Dynamic support reactions include
consistent inertia; kinetic change and backward-Euler numerical kinetic loss
remain separate from material heat.

The returned displacement, native histories, point enthalpy/temperature,
energies and optional velocity are one complete candidate. Source-domain,
history, convergence, cancellation or energy refusal changes no accepted input.
A caloric owner must jointly accept this candidate and its heat transaction;
this adapter does not predict temperature or add heat transport, finite strain,
phase change, moving mass, or a new constitutive return map. Focused regressions
cover free/restrained heating, thermal softening and residual set, nonuniform
warped B-bar quadrature/tangents, consistent inertial support work, and
source-bound Hill-Voce rate loading/holds at distinct physical durations.

The native `SmallStrainHexProblem::reference_quadrature_volumes` accessor
reuses the complete original geometry/budget/cancellation admission and returns
eight reference-volume weights per element in native state order. It evaluates
no constitutive stress, tangent or state update. Caloric owners can therefore
bind warped-cell point masses without taking a fictitious material interval.

## Sparse buckling above the dense fixture bound

`stability::sparse` solves `(K + lambda G) phi = 0` using native bounded
sparse LDLᵀ, positive-factor whitening and
`fs-spectral::service::EigenService`. Existing
`buckling_loads[_cancellable]` dispatches above 4096 free coordinates;
smaller pencils retain the legacy path. The explicit sparse entry also
admits smaller cases fitting its native search block.

`SparseBucklingLimits` bounds input dimensions/entries, symbolic fill/work,
numeric factor payload/work/front size, per-application operator work and
total operator applications. The caller separately supplies the outer
iteration allowance. No global dense K, G or inverse is assembled. Dense
fronts remain smaller than n and below the configured limit. Native
cooperative sparse products and positive-factor triangular actions poll
cancellation through preparation, iteration and final publication.

The complete factor must have positive finite scalar pivots. Singular,
indefinite or unresolved bases refuse without regularization; 2×2 pivots
are disabled. G may be indefinite, tensile or zero. Largest algebraic
inverse-load roots select compressive modes. Independent matrix scaling
preserves units/reference-load multipliers. Whitened and original-pencil
relative residuals must pass 1e-8, and modes have unit negative-geometric work.

Both dense and sparse routes restore physical loads through binary
mantissa/exponent arithmetic. Original-pencil residuals use the actual
rounded load published to the caller, including its subnormal precision.
An accurately solved transformed root cannot excuse an inaccurate returned
physical factor. Exactly representable subnormal factors remain admissible;
unrepresentable residual coefficients or failed physical residuals refuse.

Resource/convergence exhaustion returns an error, never partial modes.
These remain numerical Ritz estimates: fill can exceed the budget, and
missed-mode exclusion, continuum convergence, nonlinear collapse and
postcritical stability are not certified. G1/G3 verification covers 5003
diagonal DOFs, small noncommuting signed pencils, small P1/Q1/B-bar continuum
comparisons, scale restoration, resource/cancellation refusals and exact
versus coarse subnormal factors. This is not a 100k-continuum performance
claim.

## Retained finite inelasticity on mixed P2/P1 tetrahedra

`hyper3::mixed::plastic::MixedPlasticTetSolid` advances isothermal,
finite inelastic loading on the existing quadratic-displacement,
continuous-linear-pressure tetrahedra. Uniform, per-element and per-quadrature
`MixedFiniteInelasticMaterial` bindings borrow the original law. The blanket
`FiniteJ2Material` implementation retains existing native/source-resolved J2
callers; explicit Perzyna, Norton, combined and recovering owners keep their
original kinetics. Every quadrature point retains the complete native associated
state, including the metric, permanent flow counters and independently recovering
resistance when selected. Linear, Voce and tabulated energetic hardening use
the same native return map and analytic logarithm derivative as displacement
plasticity; no bulk modulus is replaced and no large bulk penalty is subtracted
from an already rounded full stress or tangent.

`duration_s()` reports the common original kinetic interval, or `None` for a
fully rate-independent field. Every kinetic point must declare the same finite,
positive interval; a heterogeneous field cannot silently mix physical clocks.
Newton trials and accepted-state inspections reuse the same old complete state
and never advance recovery or flow independently of the accepted endpoint.

`OriginalMaterial` preserves the original log-squared bulk law on its locally
convex branch. `ExactIncompressible` explicitly imposes the existing weak volume
constraint with no bulk penalty or storage. Physical storage is the native
deviatoric/hardening remainder plus `U(J(p))` at the pressure-reconstructed
volume. Pointwise `U(det F)` is not added to that mixed energy. Both actual
deformation and reconstructed volume retain the material's original bounds.

One `advance` is one physical endpoint: solver `load_steps` must equal one and
Newton/line-search trials always use the same accepted histories. Native mixed
equilibrium, native old/new state admission, independently re-evaluated storage,
reduced-potential consistency, pressure-constraint work, support/dead-load work,
and the caller's signed endpoint-work-defect allowance must all pass before
publishing a successor. Failure or cancellation leaves the original state
available for retry or forking. Native admission preserves legitimate Perzyna
overstress and creep stress; equilibrium yield is enforced only by owners that
select it. Plastic/viscous dissipation, combined creep dissipation, flow
potential, bulk storage and numerical endpoint defects remain separate result
fields. `inelastic_dissipation_j()` sums only the two physical heat channels.
The native sole-Norton convention reports its creep heat in the historically
plastic-named channel, consistently with the existing displacement owner.

`dynamics::MixedPlasticTetMotion` binds an accepted solid to explicit original
reference density and nodal velocity. Its consistent P2 mass is assembled once
from the same admitted quadrature and shared by successors. Signed quadratic
mass couplings are retained; deformation never changes the mass inventory.
The existing mixed endpoint solver receives the native backward-Euler inertia
only in its displacement block, preserving the pressure equation and frozen
complete material histories. Rate-dependent material intervals must equal the
momentum interval exactly. Source ownership, pressure, displacement, history,
velocity and original mass are returned together only after all gates pass.

Dynamic reactions include inertia. Independent kinetic change, inertial work,
nonnegative numerical kinetic loss, physical inelastic heat and signed material
endpoint error remain distinct. The combined work identity is gated, as is the
caller-specified numerical kinetic-loss allowance. Accepted-endpoint weak volume
admission is inherited from the mixed solve; it does not assert an instantaneous
volume-rate constraint or time accuracy without refinement. Original-mass
acceleration loads, kinetic energy and current-geometry momentum use the shared
mass owner without another constitutive update.

`tests/mixed_plastic.rs` checks coaxial finite loading, free unloading, reverse
yielding, permanent deformation and pressure against independent scalar
Hencky/J2 equations at moderate and very high bulk modulus and at exact
incompressibility. It also checks quadrature bindings and transactional
refusal/cancellation/replay. `tests/mixed_plastic/bending.rs` adds a bounded
non-affine two-mesh plastic cantilever diagnostic and compares high finite bulk
response to exact incompressibility. `tests/mixed_inelastic.rs` adds finite-time
relaxation against independent scalar kinetic equations, complete recovery
memory, shared-clock admission and rejected-step/retry cases.
`tests/mixed_inelastic_dynamics.rs` covers free translation, consistent mass
acceleration, finite inelastic work accounts, kinetic clocks and transactional
motion acceptance. These manufactured cases are numerical
evidence, not material calibration, a universal inf-sup theorem, continuum
convergence, or an injectivity certificate between quadrature points. This
adapter's evolving-temperature counterpart is described below. Porous/ductile
damage, contact and follower loads retain their separate owners.

### Evolving temperature in retained mixed inelastic solids

`hyper3::mixed::plastic::thermal::MixedThermoplasticTetSolid` owns one original
thermal material family and supplied temperature per element, with complete
typed history at every native P2/P1 quadrature point. A virgin configuration
starts at each family's original reference temperature. Subsequent heating,
cooling and loading occur through actual accepted steps. Native and source-bound
plastic, overstress, creep, combined and recovering families use the same driver.

Old history and pressure are admitted at OLD temperature. The new temperature
and physical duration freeze each element's native law once; all its quadrature
points and every Newton trial reuse that view. Original source families and
prepared reference geometry are shared across successors. The existing mixed
solver supplies compatible displacement and continuous pressure without a new
return map or a replacement bulk modulus. Exact elastic incompressibility
constrains the thermal volume weakly, preserving the declared expansion.

The native mixed thermal account supplies physical free energy, pressure-dependent
entropy and mechanical internal energy. Plastic/creep heat, reversible exchange,
variational flow potential, weak-pressure work and signed endpoint error remain
separate. The error allowance gates the sum of absolute quadrature defects so
opposing local errors cannot hide by cancellation. `element_caloric_heat_j`
returns the integrated physical transfer for each original element; it includes
reversible exchange and excludes numerical loss. It must be deposited once by
an external caloric owner. Prescribed temperatures alone do not solve heat
transport or establish the pressure-coupled caloric capacity condition.

Optional `with_inertia` binds original density and nodal velocity to a cached
signed consistent P2 mass. The same physical interval drives material kinetics
and backward-Euler motion. Heating never changes reference mass. Reactions
include inertia, while kinetic change and numerical damping remain independent
of material heat. Equilibrium, native point identities, reduced potentials,
pressure work, total work, material-defect and kinetic-loss gates all precede
publication of temperature, displacement, pressure, history and optional velocity.
Failure or cancellation leaves the accepted state unchanged for retry.

`tests/mixed_thermoplastic.rs` specifies thermal expansion, clamped heating,
loaded temperature-dependent yielding/creep, original-mass motion, spatial
temperature binding and transactional refusal/retry against bounded independent
equations. The original quadrature, material, local bulk branch, timestep and
mesh-convergence limits remain in force.

### Integrated expansion in native elastic frames

`Frame2d::with_thermal_material_sections` binds the shared material owner's
`ConcentricThermoelasticSection` at every original Lobatto station. The native
force-based compatibility solve, support reactions, corotational initial-stress
tangent, modal stiffness and frozen-temperature Newmark mechanics consume the
same source-integrated strain and current elastic response. Reference member
mass, rotary inertia and geometry remain unchanged even when source rho(T)
changes. Every layer retains its original material/card selection, reference
temperature, geometry and expansion path; changing an inelastic history or a
separate shear/viscous operator through this operation refuses.

Residual/tangent trials allocate no layer recovery output. Checked recovery
reports per-layer thermal mismatch, signed outer-fiber axial stress ranges and
axial-plus-bending elastic storage. Layer bounds apply to total and elastic
outer-fiber strains. The existing resultant-only tensile-strength adapter
refuses these fields because a zero section resultant can conceal opposing
thermal layer stresses. Refresh/cancellation is transactional. The operation
refreshes internal compatibility at the retained nodal pose and does not itself
solve global supported equilibrium or advance thermal/mechanical time.

`tests/frame_thermal_elastic.rs` provides analytical G1/G3 definitions for free
and restrained spatial heating, current-modulus compliance and frequencies,
unchanged reference mass action, layered mismatch stress/storage, cooling,
negative expansion, rotation, source/reference refusal and cancellation. This
is prescribed-temperature, small-strain, circumferentially uniform concentric
expansion. It does not supply a heat equation, reciprocal thermal exchange,
radial expansion, cross-section thermal-gradient bending or a material yield
calibration. Existing linear harmonic/viscosity routes retain their explicit
section admission limits.


## Material-admitted harmonic support excitation

`ViscoElements::prepare_supported_harmonic_response` prepares prescribed
vibration on original global coordinates that remain represented in the source.
It shares the force-driven owner's original frequency-dependent pencil:
`A = K'(omega) - omega^2 M - i K''(omega)`. The native supported solver
retains the full consistent cross inertia and every Kelvin/Maxwell cross term,
factors only the free principal block, and independently computes reactions,
applied-force power, support-actuator power and intrinsic material heat.
Applied loads on supported coordinates are subtracted once from reactions.

Every returned displacement includes the prescribed coordinates and passes
whole-cycle constitutive and compatible B-bar strain admission at every original
quadrature point. Those same vectors enter the existing single-tone and coherent
multitone material observers without another solve or projection. All native
TET4/TET10/HEX8/HEX20/HEX27 assemblies use this one path. The source is immutable;
cancelled or refused cases publish no partial response and do not poison reuse.

Moving supports must be retained during source construction. Coordinates already
eliminated remain homogeneous stationary supports, and a request to move them
refuses. Empty/all-represented support maps are admitted. Frequencies must be
positive and admitted by every original source. Dense payload admission before
assembly is `344*n*n + 544*n` bytes; sparse/material temporaries retain native
source limits. No penalty, artificial damping, modal truncation or additional
constitutive fit is introduced.

`tests/visco_harmonic_supports.rs` contains independent P1 complex-pencil and
real-time Maxwell/Kelvin integral references, applied/support work closure,
empty/all-support limits, close-tone beat heating, source/strain/budget refusal,
cancellation and reuse cases. This is a particular periodic response at fixed
material conditions; startup, actuator dynamics, evolving material state,
spatial convergence and measured physical validation are separate concerns.


### Source-bound laminate thermal frames

`add_laminate_thermal_element` and `with_laminate_temperatures` feed the
original laminate cards and integrated per-layer expansion into native
force-based `Frame2d` equilibrium. The shared `PreparedLaminateBeam` reduction
is an aligned, narrow, free-sided Euler–Bernoulli approximation: each layer
relaxes transverse normal/shear stress independently. Mechanical or thermal
shear/twist coupling refuses per layer; cancellation in aggregate coefficients
does not admit an unsupported layer.

Source reference temperatures, geometry, density-centroid centerline, mass and
centered rotary inertia remain fixed while current elastic moduli, thermal
force and thermal bending moment refresh. Axial/bending coupling survives
mass centering. Native thermal section flexibility uses a scaled checked
inverse without an axial-stiffness determinant floor. Recovery includes layer
face stresses and positive residual elastic storage at zero net resultants.
Temperature changes preserve nodal pose/velocity and solve station
compatibility transactionally; heat and changed-storage accounting remain
caller-owned.

This is not a finite-width constrained plate, interlaminar-stress,
transverse-mode-inertia, torsion, finite-strain or reciprocal thermal model.
Source temperatures are piecewise constant by layer; continuous gradients
require through-thickness refinement. `tests/frame_laminate_thermal.rs`
defines independent free/restrained bimetal, nonuniform heating, inertia and
modal-refresh, thin-flexibility, layer-refinement and cancellation regressions.


## Held-eigenstrain CutFEM dynamics

`cutfront::dynamics3::CutSecondOrderProblem3::new_with_eigenstrain` extends the
original Cartesian/adaptive structural-time adapter to `M a + K q - b_free = f`.
The native equivalent force is cached once in free coordinates at the borrowed
source's fixed stiffness scales. Physical reference-density mass retains its
independent cell occupancy. Caller forcing is external only. Tangents and existing
state/force adjoints apply pure K and hold free strain fixed.

Initial, actual generalized-alpha stage and endpoint states use the original
thermal admission owner: total small-deformation and elastic calibration limits
are separate, checked at retained bulk points without allocating stress outputs.
Construction validates binding and policy without assuming the unexpanded state
is admissible. The original consistent initialization, split resistance, nonlinear
iteration, cancellation and accepted-state publication remain the sole time owners.
`cutfem_eigenstrain_dynamics` adds hot equilibrium, consistent acceleration,
translated trajectory, affine tangent/adjoint and rollback regression sources.
The field, geometry, tensor and scales stay fixed during the trajectory; evolving
temperature or prestress and moving supports require their own explicit
integration and are not silently approximated here.

### Frozen CutFEM prestress dynamics

`CutSecondOrderProblem3::with_reference_prestress` prepares signed native G from
the adapter's held eigenstrain and an admitted full reference displacement q_0.
Its resistance is `K q - b_free + G(q-q_0)` and its tangent/state-force transpose
is K+G. The geometric displacement difference is formed before multiplication,
and the thermal offset is removed before adding the geometric force. Thus a
nonzero heated reference gains no spurious Gq_0 force. Actual stage and endpoint
displacements retain the original total/elastic strain admission. Strong support
rows must be exactly zero; adaptive constraints and selected Nitsche geometric
flux stay with the native preparation owner. Physical consistent mass and cell
occupancy are unchanged; signed compression is never clipped.

Preparation bounds bulk plus selected support points and scalar geometric
entries, with two additional native-sized vectors for the reference and action
scratch. The immutable source borrow preserves material/scales identity, and
read-only `reference_prestress` allows modal queries to reuse that exact G.
Reference equilibrium remains caller-owned. This is a frozen initial-stress
approximation, not a finite-deformation/follower-load tangent, evolving thermal
coupling or stability certificate. State/force adjoints hold G and q_0 fixed.
The focused prestress child of `cutfem_eigenstrain_dynamics` defines Cartesian
and hanging-node signed-work, reference-balance, mass and transpose regressions,
an independent dense midpoint trajectory, and physical refusal/cancellation
rollback checks.


## Material-aware support vibration with physical port loading

`ViscoElements::prepare_supported_harmonic_response_with_port_impedance` combines
the original consistent inertia, dispersive storage and loss with complete
`-i omega B Z B^T` loading before the native free/support partition. Force-map
rows use all represented original DOFs, including every moving support. Global
support labels, original material frequency bands and complete-cycle strain
admission retain the existing source authority.

Responses report full physical-port flow, effort and absorbed boundary work
separately from Kelvin/Maxwell material heat and signed actuator work. Reactions
include every cross-support inertia, material and physical-load term. Structural
and port budgets remain separately admitted; failed load cases do not mutate the
cached factor or hereditary material state. The original material observer can
independently observe the same complete displacement field. Focused Tet4 cases
compare direct loaded equations, reaction, original-quadrature material heat and
port power across frequency, and retain domain, strain and cancellation refusals.


## Signed CutFEM prestress modes

`cutfront::prestress3::CutPrestressModes3` targets low algebraic Ritz estimates of
`(K_e + G) phi = lambda M phi` using the original Cartesian/adaptive stiffness,
prepared signed geometric stiffness and consistent reference-density inertia.
Only the free-coordinate mass is sparsely assembled and positively factored;
`fs-sparse::PositiveWhitening` and `fs-spectral::EigenService` own the transforms
and LOBPCG iteration. K/G stay matrix-free, with no coordinate probing or dense
production pencil. The mass factor is reused across queries. Strong/hanging
constraints, physical occupancy and frozen material-scale identity are retained.

Returned values retain their sign, including negative and near-zero values.
Each physical mode is mass-normalized, checked for M orthogonality and verified
against fresh native K/G/M actions in the mass-whitened residual norm, with
explicit absolute and relative tolerances. The spectral residual is also admitted.
DOF/cell/entry limits, sparse factor work/storage, spectral workspace and total
operator applications are bounded; cancellation publishes no partial packet.

These are residual-admitted estimates, not a complete-spectrum/inertia theorem,
critical-load factor, follower-load correction or finite-deformation tangent.
The caller separately establishes the intended prestress equilibrium. The four
`cutfem_prestress_modes` regression sources compare an independent dense Cholesky
oracle, signed/rigid modes, adaptive nonuniform mass, and resource/owner refusal.
The `cutfem_prestress_modes` example demonstrates a synthetic restrained bar with
explicit free extension and an independent native equilibrium check; it reports
signed squared-frequency estimates without a measured-material or buckling claim.


## Relaxed static displacement in coherent material observations

`ViscoHarmonicMaterialSpectrumObserver::observe_window_with_bias` uses one real
bias value per original represented coordinate. Original constitutive and
compatible B-bar kinematics recover and separately admit the full DC/AC strain
waveforms; already eliminated supports remain fixed. Point reports retain both
static strain fields. Nonzero physical strain requires original DC applicability;
exact zero strain, including rigid translations, preserves AC-only admission.
The material owner retains direct equilibrium storage/work cross terms without
adding Maxwell/Kelvin heat. Original quadrature, physical interval, cancellation
and sizeof-based point payload bounds remain. No static equilibrium, geometric
prestress, thermal strain or causal Maxwell memory is inferred.

## Causal harmonic evolution of original quadrature histories

`ViscoElementState::bind_harmonic_clock` retains an accepted original history
and a normalized DD source instant. A prepared causal harmonic material
observer advances that history through a prescribed full nodal waveform,
including moving supports and a real static bias. Each original integration
point supplies its retained raw-Mandel Maxwell state to the shared material
owner. Both constitutive and compatible B-bar strain fields pass original
source admission. No reduced/modal strain or settled-memory substitution is
made.

The endpoint carries original per-point memory, stress, full nodal
displacement and cumulative ledgers together with an independent interval
work/heat/storage field. Exact source identity, start displacement/strain
continuity and an advancing DD clock are admitted. Empty AC banks permit
natural memory relaxation about a held field. Payload/work bounds include
returned histories, temporal scratch, nodal reconstruction and original
per-point constitutive/compatible strain recovery. Refusal or cancellation
leaves the borrowed accepted history unchanged.

This is constitutive evolution under prescribed displacement, not a mechanical
equilibrium, kinetic-energy or support-work solve. Material laws remain fixed;
accepted thermal eigenstrain remains fixed with its original source snapshot.
Four G1/G3/G4 regression groups specify preloaded Tet4 memory and moving
supports against independent real-time oracles, original time-domain
continuation, free relaxation on large DD clocks, B-bar admission, and
source/continuity/budget/cancellation refusal. They were authored and reviewed
but not formatted, compiled or run after the execution environment became
unavailable.

### Frozen graded density in Tet10 thermal valve motion

`QuadraticTetProblem::solve_thermal_tangent_with_reference_density` and
`reduce_pressure_loaded_valve_with_tet10_reference_density` retain the native
positive P1 reference-density field as four corner values per original Tet10
cell. The same consistent mass enters the thermal prestress pencil and modal
motion; curved graded cells retain the native six-point mass-quadrature floor.
The complete cell/corner profile binds the prestress identity and must remain
identical during modal rebinding, including density moments hidden by supports.
Heating may change constitutive metadata without replacing this frozen matter
inventory. This does not introduce ten-node density interpolation, density
averaging, a new constitutive law, or a new integration rule.

## Held thermal strain in causal harmonic material histories

Causal harmonic evolution retains an accepted thermal snapshot with its
original material/source identities, reference temperatures and raw samples.
Constitutive static bias subtracts the original projected free strain;
compatible B-bar bias subtracts the raw point free strain. The complete total
and mechanical waveforms independently satisfy the original material norm
and each thermal sample's component ceiling. Pure free expansion preserves
AC-only constitutive admission when mechanical DC is zero.

The returned point biases describe mechanical strain. Full nodal displacement
and the original thermal snapshot remain with the endpoint history. No thermal
sample is rebound or averaged away. Frozen eigenstrain contributes zero new
eigenstrain work; cumulative eigenstrain and parameter work stay unchanged,
while actual retained Maxwell energy may become heat. Existing cancellation,
source identity and resource bounds apply.

Two additional regression groups specify free expansion with AC-only material
support, a future-phase thermal-ceiling refusal, retained cancellation state,
and nonuniform Hex8 B-bar memory/heat against independent original spring
solutions followed by native history continuation. Rustfmt and bounded source
review completed; these Cx-backed tests remain unexecuted because the original
runtime dependency closure is unavailable. Evolving free strain and material
coefficients require their original separate owners.

## Coherent nodal continuity for causal material histories

Causal solid observation combines equal-frequency full nodal phasors before
phase evaluation and tolerance scaling, then adds the fixed bias to the
resolved coherent field. Cancelled authored amplitudes cannot erase a retained
bias or hide an actual nodal jump. This check includes prescribed supports and
rigid translations, whose discontinuities need not appear in constitutive
strains. The explicit nodal-work bound includes the bounded carrier grouping.
`cancelled_rigid_carriers_preserve_bias_without_hiding_nodal_jumps` exercises
both refusal of a new rigid offset and preservation of an already accepted one.


### Graded reference density through hexahedral thermal prestress

`Hex8ElasticProblem` and `HexElasticProblem` expose `solve_thermal_tangent_with_reference_density` and `solve_thermal_tangent_with_reference_density_and_prescribed_dofs`. The positive Q1 eight-corner reference inventory reaches consistent prestressed mass through the native graded assembler; changed elastic-material density never refreshes that inventory. The prescribed-support entry retains physical nonzero support displacement and homogeneous incremental vibration supports. Hex8 preserves Standard/BBar and `max(3, stiffness_points_per_dir)` mass integration. Hex20/27 retain independent explicit mass quadrature and their original affine/curved admission. Density and mass-rule admission precede equilibrium. The complete profile and mass order bind prestress and configuration identities, including inventory whose moments are hidden by support elimination. Native thermal equilibrium, stress recovery and signed geometric stiffness retain their existing owners and limits; these entries add no finite-strain, follower-load or modal-transfer capability.

## Typed tensor stress histories on native reference surfaces

`loads3::SurfaceLoad3::Stress` maps a finite symmetric WORLD tensor to each
selected exterior face's traction `sigma * n_ref`, then uses the existing P1
surface-force integral. Tensile stress is positive; `sigma=-p I` is the same
load as positive inward pressure. The original geometry owner supplies normals
away from the incident tetrahedron, independent of connectivity winding.

`loads3::history::SurfaceHistory3::Stress` consumes the existing scenario tensor
carrier and shared history/frame evaluator. It admits dimension-only Pa or
exact semantic Stress in static/instantaneous form, with even parity and a
symmetric sampled tensor. Pressure semantic aliases, amplitude forms, skew
tensors and unsupported source/frame policies refuse. Both tensor indices
rotate; no source table is cloned or rescanned at each ordinary endpoint.
Static equilibrium, reactions, support histories and sparse Newmark trajectories
reuse this load path with their existing transactional publication and
cancellation. Event-aligned trajectories use the shared component-history
index, account for all nine components, retain exact held left limits and
register the tensor's frame ancestry. Periodic tables require an explicit
unrolled schedule, just as existing vector event histories do.

Focused regressions in `load_history3` specify an independent affine elastic
solution, closed-boundary force/moment balance, hydrostatic equivalence and
rotating stress resultants. `transient_history3` specifies independent nodal
forcing, an analytic ramp-driven oscillator, held stress jumps, event-budget
accounting and late-refusal/cancellation rollback. These tests have not executed:
Cargo typechecking and runtime verification were blocked when the build
environment disconnected during dependency preparation. The pre-recovery
source was syntax-checked with pinned rustfmt; the recovered event-scheduler
extension has source review only. This is uniform prescribed reference dead
loading in the small-strain owner; it does not infer follower geometry,
spatial stress interpolation or a finite-deformation Cauchy-to-Piola conversion.


### Consistent acceleration loads on graded reference inventory

`ReferenceAccelerationField::{Uniform, PerNode}` and the native `assemble_acceleration_load_with_reference_density` methods integrate `rho_ref N_i a_h dV_ref` on Tet10, Hex8, Hex20 and Hex27. The returned node-major xyz load includes supported rows and acceleration supplied on supported nodes. Quadratic shape-function signs remain intact. Uniform acceleration and identical nodal values produce identical loads; callers supply the physical acceleration and its sign in mesh Cartesian axes.

A shared matrix-free mass action uses the same positive P1/Q1 density interpolation, original geometry and native mass quadrature as consistent mass assembly. It constructs no global stiffness or mass matrix. Hex8 retains Standard/BBar-independent mass integration; curved graded Tet10/Hex20/Hex27 retain their existing quadrature floors. Material density metadata cannot replace an explicit reference profile. Existing mesh/work admission and cancellation apply, with finite acceleration and force checks. These are prescribed reference dead loads for the existing equilibrium, reaction and prestress owners; no follower geometry or transient acceleration law is inferred.

Four focused regression groups cover analytical force/centroid torque/rigid work on affine and curved maps, full consistent mass action and supported columns, actual native self-weight equilibrium/reactions at equal total mass, and quadrature/field/overflow/cancellation refusal. All four and fourteen existing Hex mass/prestress regressions passed in the native source-subset harness; full workspace Cargo validation was unavailable.

### Native reduced-plane pseudo-arclength continuation

`HyperProblem::prepare_continuation` binds the original elastic problem to
`hyper2d::continuation::PlaneContinuation`. The admitted formulations are
finite-bulk P1/Q1 plane stress, exact-volume P1/Q1 plane stress, and native
Q1 plane-strain mean dilatation. Each trial reuses the original prepared
reference geometry, point-material/reference owners, transverse or volumetric
condensation, complete energy Hessian, and native deformation-domain checks.
Frozen incremental inelastic potentials refuse; this is a held elastic path.

With positive physical scales `L,F` and original thickness `h`, free
coordinates are `q=u_f/L`, prescribed motion is `u_c=lambda*g_c`, and
dead tractions scale with `lambda`. The native equations are
`R=h*(f_int-lambda*f_dead)/F`, `R_q=h*L*K_ff/F`, and
`-R_lambda=h*(f_dead-K_fc*g_c)/F`. The support column is formed from
the complete current tangent before restriction. Prescribed coordinates stay
outside the free arclength metric. One exact-state cache shares native
constitutive assembly across residual, tangent and parameter-column requests.

The existing bordered continuation owner supplies predictor/corrector solves,
step control and resumable path state, including negative tangents and simple
load-control folds. Preparation requires nonempty free coordinates under an
explicit caller cap, finite representable scaling, and a nonzero free dead
force or prescribed-motion ramp. `start_from` admits a complete physical
equilibrium with its exact prescribed values; resume rechecks the original
force equation. A failed/cancelled advance returns no candidate and preserves
the input checkpoint. Context cancellation is polled through native assembly
and the existing continuation callbacks. Callbacks must remain pure while held.

Response retains original support reactions, explicit free-force residual,
stored elastic energy, and stored energy minus current dead-load work in
physical units. Moving-support work is not classified as dead-load potential.
Path events remain numerical diagnostics; no certified bifurcation, continuum
stability, wrinkling, global injectivity, evolved history, or updated reference
state is inferred.

`tests/hyper_plane_continuation.rs` specifies analytical biaxial energy,
reaction and dead-work checks across all five mesh/formulation combinations,
exact split resume, displacement/parameter finite differences, heterogeneous
held-reference parity, and frozen-history/domain/cancellation refusal.
A bounded native source check exercised the new support restriction and the
existing bordered continuation, CSR and dense factorization on an independent
moving-support von Mises reference: 40 derivative comparisons, a 140-step
limit-point crossing, exact 70+70 resume, and cancellation passed. The complete
HyperProblem integration regressions remain unexecuted in the dependency-
blocked workspace; no large-system performance claim is made.

## Scenario traction rows reach the native stress load owner

`SurfaceHistory3::scenario_traction` now binds an Elasticity/Traction vector or
3-by-3 stress tensor to explicitly selected native exterior faces. The existing
Elasticity/Neumann vector path remains available. Both declared Stress vectors
and tensors use the same quantity-admission function as scenario validation;
volume force and acceleration contracts retain their original dimensions.
Canonical tensors therefore reach the existing two-index frame transform,
reference-normal traction integral, static/reaction solve, transient histories
and event schedule without a second load assembler. Legacy scalar traction
rows still require an explicit physical interpretation outside this adapter.

The independent affine-equilibrium regression now round-trips a complete
scenario tensor row before native assembly and solution. The manual-force
comparison now enters through a semantic Stress vector boundary. These changes
retain the small-strain reference-surface scope and original cancellation and
publication behavior. Source review is complete; compilation and executable
regressions were unavailable after the runtime-dependency build disconnected
from the execution environment.

## Mixed finite-strain equilibrium paths through limit points

`hyper3::mixed::continuation::MixedContinuation` connects the original prepared
TET10/P1 material system to the existing bordered pseudo-arclength solver.
`MixedUpProblem::prepare_continuation` admits native hyperelastic sources;
`PreparedMixedUpProblem::into_continuation` and the original uniform-pressure
and stationary-gravity wrappers retain their already admitted geometry and loads.
The source must implement both `HyperelasticMaterial` and
`MixedFiniteStrainPotential`; frozen incremental inelastic potentials are not
substituted for a history-independent equilibrium path.

The original length and stress scales form diagonal variable scale D. With
W = stress_scale_pa * reference_volume_m3, free coordinates are D^-1 x,
the normalized residual is D r/W, and its tangent is the congruence D K D/W.
Pressure residuals remain volumes, paired with their pressure scale. Original
dead forces, held current-geometry pressure/gravity, prescribed displacements,
and an authored pressure gauge datum all ramp with lambda, exactly as in the
native static endpoints. The negative parameter derivative includes -D K_fc g/W
for every prescribed displacement or pressure datum; constraint rows receive
the same full tangent lifting. Held loads use their original quadrature and
geometry derivatives without another constitutive assembly.

A shared exact-point cache retains the residual, tangent and load column.
Every evaluation preserves native material and geometry domains. Both Newton
and continuation use one original constant-pressure-null-mode test. Every
accepted continuation step additionally checks ALL pressure equations,
including the gauge equation omitted from the algebraic solve. Thus an
incompatible prescribed volume cannot be hidden behind a pressure pin.
The complete requested advance is transactional with respect to its input
PathState; cancellation or any refused intermediate endpoint returns no
successor. Continuation state includes independent pressures and retains the
native deterministic restart representation.

Responses expose physical displacement, pressure, reactions, free-force
residuals, all volume residuals, sampled det(F), and the original MATERIAL
MIXED SADDLE potential. This potential is not relabeled as physical stored
energy, and external loading potentials are excluded. Preparation and path
growth obey the original assembly/factor payload caps, including retained
CSR, bounded dense/Krylov border storage, and candidate path records.

G0/G1 regressions in `tests/mixed_up/continuation.rs` cover an analytic
incompressible Ogden tensile fold and split restart, finite-bulk/exact mixed
tangent and load-column finite differences with moving supports plus native
pressure/spatial gravity, pressure-gauge shear motion, incompatible volume,
and cancellation on a cached point. The mixed saddle tangent is generally
indefinite even on a stable branch: generic definiteness events are not a
mechanical bifurcation certificate, and unverified branch-switch directions
are refused. No mesh/inf-sup convergence, spatial injectivity, stability,
evolving temperature/reference, or constitutive-history claim is added.

## Full forced harmonic response of finite-strain mixed equilibria

Both finite-compliance and exact-incompressible vibration owners expose
`prepare_harmonic_response`, returning a reusable native harmonic factor at a
held physical frequency. `solve` accepts peak force phasors in ALL original
node-major displacement coordinates, including support entries. It retains the
original full consistent reference mass, mixed material/geometric tangent,
follower-pressure or conservative-gravity tangent already in that equilibrium,
pressure ordering, and admitted gauge. Pressure has no mass. `solve` keeps zero
incremental support displacement. `solve_with_support_motion` additionally
accepts finite peak displacement phasors on strictly increasing original
prescribed DOFs; omitted supports have zero increment. These are perturbations
about the original absolute static motion, not a new equilibrium or support map.

Finite compliance uses the original condensed stiffness and `dp=R du`. Moving
supports add `p_g=H^-1 B_support g`, using the ORIGINAL pressure-compliance LDLT
factor retained from vibration preparation under its original factor-byte cap.
There is no second pressure factorization at a frequency or a forcing solve.
Exact incompressibility projects force through the COMPLETE retained nullspace, then
recovers pressure from the original forced momentum equation using one reused
native QR factor. Constraint-blocked loading can therefore produce nonzero
pressure with zero displacement; summing modal pressure shapes alone would miss
this contribution. A zero-dimensional admissible space remains a pressure-only
solve for fixed supports. For moving supports, the SAME scaled `B_free^T/s=Q1 R`
factor gives the minimum-Euclidean-norm lift `u0=Q1 R^-T(-B_support g/s)`; native
fs-la LU owns the small triangular-factor solve. The existing homogeneous
harmonic factor then solves for the correction in `u_free=u0+Z q`. All original
weak volume equations, including a gauge row omitted from QR, must accept the
complete displacement. Incompatible prescribed volume is refused even when the
independent equations have a solution. A zero-dimensional homogeneous space can
still move through this unique compatible lift. No modes are discarded and no
bulk penalty, normal equations, pressure inertia or physical damping is added.
The exact-constraint Galerkin mass and stiffness each compute one triangle and
mirror it after qualifying the original source's conservative symmetry. The
finite-compliance mass and completed Schur stiffness likewise qualify both
triangles before canonical storage. This preserves exact stored symmetry
without accepting a materially skew source or weakening the native harmonic
owner's symmetry gate. Physical residuals still use the original full operators.
Native harmonic admission refuses singular lossless poles.

Returned displacement, velocity and pressure use original physical coordinates.
Support reactions are `(Kuu-omega² M)du + Kup dp - authored_force`; full support
mass coupling is retained and an authored support force is subtracted once.
The complete dynamic action on the prescribed displacement and constraint lift
enters the free right-hand side, including `-omega² M_free,support g`, pressure
coupling and all original conservative load tangents. A static support force or
unlifted nullspace projection cannot substitute for this base excitation.
Every original pressure equation, including the gauge row, and every free
momentum equation must meet their respective backward-error budgets before
publication. Each physical equation block uses `max(abs(residual))` divided by
the maximum row sum of absolute uncancelled terms, matching the original mixed
vibration owner. Momentum and weak volume are normalized separately without an
absolute unit floor. An individually zero-action row cannot reject a valid
hydrostatic pressure response merely because both its residual and denominator
are floating-point noise; every row still contributes to the block audit.
Reported cycle-mean input power is `Re(F^H v)/2`, including applied
support-coordinate loads. Support power is separately `Re(R_support^H
v_support)/2`, positive for work supplied BY the support ON the solid. Their sum
is the signed lossless power-balance residual. These are independently computed
external work terms, not material heat or a passivity certificate. The nested
native harmonic response describes the homogeneous correction only; its force
and power do not replace the complete physical account when supports move.

`HarmonicResponseControls` counts free displacement coordinates BEFORE exact
constraint reduction. Its payload cap covers the adapter, QR, native harmonic
factor/preparation, affine lifting, and one physical response; the borrowed
original matrices and finite-compliance factor retain their original budgets
and are excluded from this extra payload count. With n free displacements,
m independent exact pressures, r=n-m (r=n for finite bulk), nu original
displacements and np pressures, the bound is `144*r²+128*r` native bytes plus
`8*r²+24*n*m+48*m²+256*(nu+np+n)+4096` adapter bytes. Native QR and the triangular
factor LU are bounded by the original 128-free-displacement cap and polled
before and after; the original sparse pressure solves retain their native
synchronous semantics with polls around each call. Sparse audits and
physical-vector work poll cancellation. The source
equilibrium and prepared factors remain immutable after all refusals.

This is a frozen-equilibrium small-perturbation response, not finite-amplitude
material admission, nonlinear transient motion, evolving source/history,
support-actuator dynamics, thermoelastic recoupling or physical damping.
`tests/mixed_up/harmonic.rs` defines G1/G3 native regression groups for the
loaded finite-bulk analytic response, hydrostatic force with zero motion,
graded inertia with pressure gauge and support reactions, and zero-dimensional
constraints plus budgets/cancellation. Moving-support regressions compare the
complete original mixed KKT solve against finite-bulk and exact affine-lift
responses, including graded cross-support inertia, DC and nonzero frequency,
zero-dimensional motion, signed actuator work, gauge-preserving rigid
translation, incompatible volume, support-map admission and cancellation.
These are discrete implementation
checks, not continuum convergence or physical calibration.
