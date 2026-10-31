# fs-vkplate: a drum kit from one solver

Every vibrating surface of a standard kit is one object here: a thin plate on a square
finite-difference grid. A drum head is a clamped plate with in-plane tension and the film's own
small bending stiffness; a cymbal is a free-edged plate with bending, a rest shape (a shallow
shell) and full von Karman coupling. Kick, snare, toms, hi-hat, crash and ride are parameter sets of the same
equations, couplings and output stage (`src/kit.rs`).

Why it exists: FrankenSim's modal percussion (`crates/fs-couple/examples/percussion`) keeps at
most 63 head modes, admits audio only for head bands inside 40-1640 Hz, and caps a cymbal shell
at 32 modes (an 18-inch crash stops near 400 Hz). A handful of lightly damped modes is a bell,
and that is what those renders sounded like. A grid solver avoids a prescribed modal cutoff,
but its usable bandwidth and dispersion still require spatial and temporal convergence checks.
Node count alone does not establish audio-band fidelity.

## Physics

- **Membranes**: `rho t w_tt = T Lap(w) - D Lap^2(w) + ...`, rim clamped, with the Mylar film's
  bending stiffness `D = E t^3 / 12(1 - nu^2)`, which takes over from tension for waves with
  `D k^2 > T` (above roughly 0.7 kHz on the slack 2-ply kick and floor-tom batters, 5 kHz on the
  snare and rack tom). Berger tension modulation (the mean stretch raises the tension) gives the pitch drop of a
  hard-struck tom.
- **Plates and shells**: bending energy averages four quarter-cell tensor energies with the full
  `(1 - nu)` twist term. Squaring a node-centred average of mixed curvature instead loses
  diagonal-wave stiffness. The quarter-cell form gives exactly `D Lap_h^2` in a uniform
  interior, second-order continuum consistency and zero bending force for rigid tilts;
  free edges are the natural boundary of that energy. Sourced bending-viscosity roots use
  the same quadrature. Von Karman coupling
  through the Airy stress function with a rest shape `z0(r)`:
  `Lap^2 F = -E t [l(w + z0) - l(z0)]` for uniform thickness,
  `l(u) = u_xx u_yy - u_xy^2`. The outer rim fixes the affine Airy
  gauge with `F = dF/dn = 0`. A hole has its own three affine trace
  coordinates; imposing the same zero trace on both boundaries incorrectly
  removes physical stress fields. Annuli use the complementary weak load
  `cof(Hess F) : [grad(w)⊗grad(w)/2 + sym(grad(z0)⊗grad(w))]`,
  retaining the boundary contribution and its exact discrete energy gradient.
  The three extra coordinates use a Schur complement of the narrow-band factor.
  Native and CUDA steppers consume the same quadratic strain rows.
  With varying thickness, the stress operator comes from the local complementary
  membrane energy `(F_xx^2 + F_yy^2 - 2 nu F_xx F_yy + 2(1+nu) F_xy^2)/(2 E t)`.
  Normal derivatives live at nodes and mixed derivatives at cell centres. This
  preserves the uniform-thickness operator and retains thickness gradients in
  the variable-coefficient compatibility equation; an average thickness cannot.
  The nonlinear force is the exact gradient of the discrete potential.
- **Loss**: either `s0 + s1 beta^2` per mode, fixed by two (frequency, T60) points,
  or an explicitly selected source-bound plane-stress Kelvin–Voigt bending tensor.
  The former is a phenomenological law, not measured material viscosity. Its gradient
  term is `2 s1 div(rho t grad(w_t))`: on a tapered surface, each grid edge uses
  the mean adjacent mass, giving equal and opposite internal forces and nonnegative
  dissipation. Multiplying the ordinary Laplacian by each node's mass instead would
  violate those properties. Uniform-mass surfaces retain their original operator;
  clock admission includes the mass-normalized damping bound on tapered surfaces.
- **Contacts** (stick, felt beater, snare wires, hi-hat rims): one-sided power laws with
  Hunt-Crossley loss, `f = K d^alpha (1 + mu d_dot)`, each solved implicitly as a monotone scalar
  equation, so a hickory tip on bronze is stable at audio rates. A zero-radius tip uses bilinear
  grid weights for both displacement and force; an incomplete boundary cell falls back to its
  nearest node. A finite striker face retains its equal-weight footprint.
- **Air**: a lumped adiabatic cavity couples a drum's two heads; each head's volume motion
  carries the radiation mass of a baffled piston (`8 rho a^3 / 3`, implicit); a monopole
  radiation resistance acts on the drum's NET volume flow (two heads moving the same way in
  space radiate as a weak dipole, which is why a tom's fundamental rings). An optional port
  models a vented front head.
- **Sound**: the Rayleigh integral of surface acceleration at each microphone, per-node delays.
  A microphone hears only surfaces facing it. Over a rigid floor (any plane, `floor_z`, `floor_y`)
  it also hears, by the floor's image, the surfaces facing that image.
  A cymbal is not in a baffle: its microphones apply the unbaffled-disc limit (a dipole,
  `(4/3pi) ka` times the baffled field, for `ka << 1`), joined to the baffled field for `ka >> 1`
  by a first-order high-pass with `tau = 4a/(3 pi c)`.

## Shared atmosphere

`Instrument::with_gas` and the `kit::*_with_gas` builders take one coherent
`fs_material::gas::GasState`. Its density and sound speed drive cavity bulk
modulus, port inertia and radiation resistance, head radiation mass, microphone
amplitude and travel time, and the unbaffled receiver filter. The exported
`InstrumentImage` retains that same state for external steppers, including CUDA.
Clock admission includes cavity stiffness and explicit gas damping, so changing
the gas cannot silently bypass the existing numerical stability check.

`live::Live::with_gas_record` and browser `KitPiece.withGas` consume the canonical
material parser's SI records, for example:

```text
293.15 50000 0
mixture 293.15 100000 mole 3 N2 0.7 O2 0.2 CO2 0.1
```

The first record means temperature [K], pressure [Pa], and relative humidity
[0..1]. Named mixtures also accept `mass` fractions; species and validity limits
belong to `fs-material`. Existing constructors resolve dry air at 293.15 K and
101325 Pa through that same owner. At fixed temperature, doubling ideal-gas
pressure doubles density and bulk modulus while leaving sound speed unchanged.

This is one prescribed stationary medium. Changing a live atmosphere starts a
new trajectory; it is not an in-place thermodynamic update of a ringing body.
The existing lumped-cavity and small-ka radiation approximations still apply.
Viscosity and conductivity remain in the retained state, but this solver does
not yet model viscothermal cavity loss or hi-hat squeeze films.

## Shared solid and striker materials

`materials::KitMaterials` binds three canonical
`fs_material::state_point::IsotropicElasticStatePoint` values: the primary
surface, secondary surface, and spherical striker. The complete source bundles
and imported material cards remain attached to the instrument and its exported
image. Use `kit::drum_with_materials`, `kit::cymbal_with_materials`, or
`kit::hihat_grid_mu_with_materials` with the same gas and geometry declarations.
Existing builders retain the standard kit presets.

The material density supplies each surface's actual mass. Young's modulus and
Poisson ratio supply its bending, Berger and von Karman operators. Striker mass
is its source density times the declared effective moving volume; the original
Hertz law uses both contact materials and the declared spherical tip radius.
Changing cards does not retune installed tension, rescale geometry to preserve
mass, or substitute a material-name switch in the physics code.

`KitMaterials::from_card_packs` accepts canonical `NormalizedMaterialCardPack`
bytes and an explicit `QueryPoint`. The live browser accepts those same packs
through `KitPiece.withMaterials` and the **Solid and striker materials** controls.
The primary pack applies to batter heads and upper cymbals; an omitted secondary
pack explicitly reuses it for resonant heads and lower hi-hats. The striker pack
replaces every preset beater, including the kick, with the declared spherical
impactor. `Live::with_materials_record` is the native counterpart. All live grid
trials, fallback rebuilds, overlapping strikes and copies retain the binding.
GPU plans also require matching solid and gas source bindings for every member.

Browser imports declare one finite positive absolute solid temperature in kelvin
and the source's temperature-axis name. A source that constrains that coordinate
must use the canonical typed absolute-temperature/static quantity. Native callers
can provide a general query point, including additional declared state coordinates.
Missing, conflicting or out-of-domain properties refuse; imported files are
bounded to 16 MiB in aggregate. These are new trajectories from rest, not an
energy-conserving replacement of a ringing specimen.

Loss and contact scope remains explicit. Drum heads built from imported cards
never keep the preset's T60 decay points, which were fitted to recordings. They
need a selected material bending-loss law, the head cards' film loss curve
(`KitMaterials::with_film_loss_curve`, below), or the declaration
`KitMaterials::with_stated_unknown_film_loss`, which enters head film and edge
loss as 0. Without one, `drum_with_materials` and every live drum refuse with
`materials::FILM_LOSS_REFUSAL`. A declared drum carries
`materials::FILM_LOSS_STATED_UNKNOWN_RECEIPT` in its notes and in `info_json`
(`filmLossStatedUnknown`): it omits the measured high-mode excess loss, so T60s
above the first two modes are upper bounds. Declared muffle dashpots remain. In
the browser, these are the `filmLoss: 'stated-unknown'` material option and
`KitPiece.withMaterialsStatedUnknownFilmLoss`. Cymbals and hi-hats without a
selected bending law still use their installed T60 decay points.

`KitMaterials::with_film_loss_curve(temperature_axis, frequency_axis, T)` reads each
head card's `loss_tangent_tensile` (`materials::FILM_LOSS_TANGENT`), a curve over a
typed cyclic-frequency axis, through the canonical resolver. Each head's bending
operator becomes a generalized-Maxwell (Prony) law fitted to that curve alone, never
to a recording (`film::fit`, `Body::attach_film_maxwell`):
- the fit is non-negative least squares on a fixed relaxation ladder (3 per decade, a
  decade beyond each end of the band);
- its receipt gives the measured sup error of tan delta against the curve (under 2 %
  for the PET card);
- the storage modulus equals the card's E at 632.5 Hz, the log-centre of Matsuzawa
  1959's 100-4000 c/s, and causality makes it rise with frequency (the receipt states
  E' at both ends of the band);
- it is stepped with fs-material's exact linear-ramp branch recursion
  (`GeneralizedMaxwell::step`) on the field K q. Its heat closes the discrete energy,
  and the leapfrog admission uses its glassy modulus.

Every grid mode therefore carries the tan delta at its own frequency times its
bending share, including the modes above the radiation table's shapes and the
`(0,1)` volume shape. The tension (a pre-stress) and the Berger in-plane stretching
stay elastic. The head radiation then carries air mass, radiation and the Stokes
viscous layer on every shape it retains (`Body::attach_head_radiation_with_viscous`),
and no film loss; each shape's receipt gives its load T60 (radiation, compact share
and viscous layer).

The load reads each shape's air mass and resistance at its loaded frequency on the
card's E, not on the law's E'. It does not iterate on the law's E'. On the PET snare E'
runs from -1.7 % to +2.8 % over the batter's retained shapes (181-3190 Hz). On the
resonant head it runs from -0.9 % to +3.9 % (343-5685 Hz). A shape's frequency
therefore differs from the one its load is read at by about half that times its bending
share: at most 1.4 % on the batter and 2.0 % on the resonant head (pure bending). Each
head's radiation scope receipt states its own figures. E' is +4.7 % at 8 kHz, where no
shape is retained.

The PET law keeps 11 branches, each updated on every node at every step. That is
0.5 % of the 3.5 mm snare render's profile. On the smallest live grids of the card
film route (not in the browser), it adds up to +36 % per output sample.

The drum's notes carry the law's receipts and
`materials::EDGE_LOSS_STATED_UNKNOWN_RECEIPT`, since no edge, hoop, shell or mount
loss law exists, and `info_json` has `filmLossCurve`. Only the native card route
uses it, offline and in `live`; the browser page's drums keep the stated unknown.

`packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk` is a sourced head-film
card, compiled with `xtask matdb-pack` from
`data/matdb/seed-v1/head-film-pet-mylar-dupont-matsuzawa1959`. It holds density
1390 kg/m³, Young's modulus 4.90 GPa (the mean of the datasheet's machine and
transverse directions) and Poisson ratio 0.38 for biaxially oriented PET, admitted
from 293.15 to 296.15 K, with material state
`polyethylene-terephthalate solid biaxially-oriented-film 0`. Its audio-band loss is a DERIVED
curve, `loss_tangent_tensile` at 293.15 K from 20 Hz to 20 kHz. It comes from Hitachi High-Tech's
DMA application brief No. 9: a 10 µm biaxially oriented PET film annealed at 130 C, not head gauge,
swept in tension at 0.5-10 Hz. The 1 Hz tan delta trace is read off the brief's plot at 1/3-octave
centres; it reads 0.0443 against the annotated 0.045 beta peak. The trace is then shifted to 20 C by
the beta relaxation's Arrhenius law, `1/T_eq = 1/293.15 K + (R / 76 kJ/mol) ln(f / 1 Hz)`, with the
brief's activation energy and no vertical shift. That gives 0.0183 at 100 Hz, 0.0254 at 1 kHz,
0.0298 at 3.15 kHz and 0.0343 at 10 kHz. These lie 1 to 3.3 decades beyond the brief's 10 Hz, so
the estimate is about ±30 %. The card also carries one composition-proxy reading for reference:
Kunugi et al. (1992) measured uniaxially drawn PET at 110 Hz and 25 C, giving a loss tangent of
0.022 along the draw and 0.029 across it (read off a figure). The kit does not consume it.

`packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk` is a sourced stick card,
compiled from `data/matdb/seed-v1/stick-hickory-shagbark-fpl-gtr282` with
material state `carya-ovata solid clear-straight-grained-12pct-moisture 0`. It
uses Wood Handbook Table 5-3a for shagbark hickory at 12 % moisture content:
- density 806.4 kg/m³, from specific gravity 0.72;
- Young's modulus 16,390 MPa, the longitudinal modulus: bending MOE 14,900 MPa
  plus the handbook's 10 % shear correction;
- Poisson ratio 0.371, a declared composition proxy. Hickory is absent from
  Table 5-2, so this is the mean μLR of its 16 hardwoods;
- E_T/E_L 0.05693, G_LT/E_L 0.06657, ν_LT 0.4668 and a transverse-plane
  Poisson ratio 0.518, the hardwood means of Tables 5-1 and 5-2, also
  composition proxies;
- side hardness 8,400 N (the load to embed an 11.28 mm ball to half its
  diameter) and compression perpendicular to the grain 12.1 MPa (the stress at
  the proportional limit), from the same Table 5-3a row.

It is admitted from 288.15 to 298.15 K, by Table 5-15. The other true
hickories' values are kept as the stated spread.

### The crash struck with the shaft

`crash-shaft` (examples/kit) strikes the crash's edge with the card stick's
shaft, at right angles (`shaft`). The stick is the Z5A shape estimate, turned
about its grip; the edge lands at its 0.33 m station. The stick's angle and
that landing point are declared, not sourced. The shaft crossing the edge is
an elliptic Hertz contact of crossed cylinders, with the wood indented across
its grain (Barber's anisotropic half-space with the Delafargue-Ulm cos 2φ).

Elastic, that contact's mean pressure would reach 119-245 MPa at 1-8 m/s. That is
above the card's side hardness read as a mean pressure, 84.1 MPa (8,400 N over
π (5.64 mm)²), so the wood crushes. The contact is a `Crush` (lib.rs):
- Hertz until the mean pressure reaches 84.1 MPa;
- deeper, that mean pressure on the Hertz area, which Thornton (1997) keeps
  proportional to the approach;
- unloading and reloading along Hertz from a permanent set.

Its potential only falls as the wood crushes, so the discrete-gradient solve
stays passive. The Hertz area is half the edge's truncation of the shaft, the
area the hardness test divides by at depth, so this is the lower of the two
readings of the same datum. The card's compression perpendicular to the grain,
a seventh of that pressure, is stated and not used.

The shaft's Hunt-Crossley loss is the card route's declared 0.2 s/m. Each hit
crushes fresh wood. Not modelled: friction along the edge and the stick's
bending. The GPU renderer refuses a crushing contact.

A striker card that carries all four elastic ratios (`TIP_ANISOTROPY`) has a
tip transversely isotropic about the grain. Its Hertz modulus is then the
Delafargue-Ulm (2004) axial indentation modulus: 7.245 GPa for this card,
against 19.006 GPa for an isotropic E_L tip. On the PET film card that gives
E* 3.200 GPa against 4.403 GPa, and a contact time 1.136 times longer. A card
with none of the ratios keeps the isotropic E / (1 - ν²); a card with only
some of them refuses. The drum's notes print both ends of this bracket. They
also print the contact's own limit: the head side of the Hertz pair is the film
as an elastic half-space, which a tensioned head film is not, and that error may
exceed the tip bracket.
Supports, snare-wire contact and hi-hat rim/friction laws remain separately
declared presets; the snare wires' steel is sourced on the card route (see
Snare wires). The striker's Hunt-Crossley loss is an independent input, not
inferred from elastic modulus. The head radiation-resistance estimate retains
its tension-only membrane-frequency approximation, using the new density; full
bending elasticity reaches the computed structural motion and Rayleigh pressure.
No measured broadband loss or stiff-head radiation qualification is inferred
from bulk elastic data. Imported materials bypass the original Mylar/bronze
grid-gain correction, so material comparisons keep the declared pressure-to-PCM
gain instead of being equalized to the original kit.

### Source-bound bending loss

`KitMaterials::with_bending_viscosity` accepts the existing canonical
`fs_material::plane_stress_viscosity::PlaneStressViscosity` for each surface and
its orientation relative to the plate grid. Each retained law must come from
the same material card and exact typed state point as the corresponding elastic
properties. `KitMaterials::from_packs_with_bending_loss` resolves those laws from
the same imported cards. The six required Pa·s properties are:

```text
plane_stress_viscosity_11   plane_stress_viscosity_22
plane_stress_viscosity_66   plane_stress_viscosity_12
plane_stress_viscosity_16   plane_stress_viscosity_26
```

They describe the symmetric reduced tensor acting on engineering strain rate
`[epsilon_dot_11, epsilon_dot_22, gamma_dot_12]`. Every entry, including zero,
must be declared. The canonical law admits a positive-definite tensor or an
explicit all-zero tensor; unsupported degeneracy refuses. Elastic Poisson ratio
does not supply a missing viscous tensor entry. At each curvature point the
solver uses the admitted rotation, actual local thickness cubed over twelve,
and engineering curvature rate `[w_xxt, w_yyt, 2 w_xyt]`. The sparse damping
matrix is assembled as `C = G^T G` from the same curvature stencils and area
quadrature as bending elasticity. Changing only loss data therefore changes
the actual motion and pressure while preserving structural mass and stiffness.

The browser's **Surface material viscosity** selection calls the
distinct `KitPiece.withMaterialBendingLoss` constructor. It declares a typed Hz
source coordinate and the two material orientations. An older WASM binary
without this constructor refuses the selection. Grid trials, fallbacks, clones
and exported GPU plans retain the complete law; GPU batches must match both
source identities and orientations. The ready material metadata reports the
two tensors, orientations, source-state identities and admitted angular-frequency
bands. Selecting this law removes the preset surface `s0` and `s1` entirely;
explicit zero viscosity does not restore the old T60 decay.

This is a fixed-temperature, broadband Kelvin–Voigt model. Coefficients must be
constant or temperature-only source curves at that declared temperature, and
the intersection of elastic and viscous source bands must cover DC through the
**internal integration Nyquist frequency**, including any oversampling. Merely
placing the query frequency inside a narrow source band does not admit a
broadband time step. Frequency-dependent sampled coefficients, narrower bands,
and clocks exceeding the supported stability margin refuse. No broadband
measurement claim is inferred from the six inputs.

The force uses backward velocity, `-C (w_n - w_{n-1}) / dt`, so this damping
discretization is first order in time while preserving the existing contact and
SAV compliance. The positive-energy condition is
`dt^2 lambda_K + 2 dt lambda_C < 4`; the material term uses a mass-normalized
absolute-row-sum bound. Kit clock selection evaluates that term at each candidate
internal rate. The structural part remains the existing rest-state eigenvalue
estimate, not a contact or large-deflection convergence certificate.

`Body::discrete_energy` includes the staggered correction
`-dt/4 v_half^T C v_half`. For an isolated linear body without other losses,
its step change is exactly `-dt ||G v_center||^2`, where
`v_center = (w_{n+1} - w_{n-1}) / (2 dt)`. Native and CUDA steppers accumulate
that material-loss contribution as a sum of squares, with no negative-loss
clipping. It is exposed by `Body::bending_loss_dissipated_j`, browser
`KitPiece.bendingLossJ`, and CUDA `Run::bending_loss_j`. The existing shell SAV
drift correction may spend the previous step's loss budget to restore its
auxiliary energy, so this counter is a constitutive dissipation contribution,
not a complete coupled energy audit. Support, contact, friction, snare and gas
losses are independently declared and are excluded from this counter. The
negative staggered correction is numerical stored energy; it is never reported
as heat. This route does not yet update temperature, add nonlinear membrane
viscosity, or implement a hereditary viscoelastic model.

### Head radiation

A drum built from imported head cards, with a sourced bending loss, a film loss
curve or a declared stated-unknown film loss, takes its heads' exterior loss and air mass from
`data/head-radiation-planar-baffled.tsv`, not from a fitted decay. The table holds `zhat_rs(ka) = Z_rs / (rho c a^2)`: the exterior
loading of clamped-membrane shapes `J_m(j_mn r/a) cos(m theta)` (with
`j_mn <= 20`) on a disk in an infinite rigid baffle. `fs_bem::modal_loading`
bakes it on the power-consistent Rayleigh route
(`crates/fs-bem/examples/membrane_head_radiation.rs`):

- A head of radius `a` in the instrument's gas sees
  `Z = rho c a^2 zhat(omega a / c)`. Density and sound speed follow the current
  gas record, live and offline.
- The ka grid is log-spaced from 0.02 to 1, then has uniform steps of 0.2 up
  to 12. A shape whose ka falls outside it, or in a supersonic span whose
  reactance is not mass-like, is refused for that shape alone. It is never
  clamped or extrapolated. The shape gets no exterior load, and the receipt
  names it as an upper bound. The summary gives the lowest in-vacuo frequency
  among refused shapes. A refused or unresolved (1,1) shape refuses the drum.
- Without a modelled head loss, only the `m >= 1` families load a head. Volume
  motion keeps the compact monopole, the piston air mass and the cavity. The
  baked `m = 0` family only prints `R_compact / R_baffled` per volume shape in
  the receipt.
- A head whose own loss is modelled (a film law on its bending operator,
  attached before the load, or the viscous face layer) retains every shape its
  grid resolves, so that the air mass and the viscous layer reach every resolved
  mode (`ModalAirLoad::new_with_losses`, `extend`). The film's loss is never in
  the load; the law is the plate's own:
  - A volume shape `(0,n)`, `n >= 2`, takes the table's baffled load less what
    the compact net-volume model already gives it (`R_c V^2` and the piston air
    mass times `V^2 / A^2`), so nothing is counted twice. `(0,1)` stays on the
    compact model alone.
  - A shape beyond the table's `j_mn <= 20` set takes the air mass of its
    untruncated evanescent field, `rho (integral of phi^2 dA) / sqrt(k_s^2 - k^2)`,
    and no radiation resistance, so its radiation T60 is an upper bound. The
    receipt compares that air mass with the table's own for `j_mn >= 15`.
  - It costs more. On the offline 3.5 mm snare a head carries 623 columns
    instead of 86, and the wires-on render takes about 8x longer. Live grids
    add only a few volume shapes, about +6.5 % per sample on the card routes.
- Each grid-resolved shape, in both parities, gets a narrowband load at its own
  loaded frequency. That frequency is a fixed point of
  `omega^2 (M + Ma(omega)) = K` with the head's Rayleigh-quotient `K` and `M`.
  The fixed point is refused beyond 100 iterations.
  - Resistance `Re Z_rr` acts as an explicit backward-velocity force.
  - Air mass `-Im Z_rr / omega` acts implicitly, through a Woodbury correction
    shared with the Berger and piston-mass columns.
- The in-family coupling is dropped. Its first-order back-reaction is printed
  per shape, and its largest value per head.
- Each shape's MAC against its family's Rayleigh-Ritz FD mode is printed. A
  shape below 0.9 is flagged.
- Shapes the grid cannot resolve (fewer than 6 spacings per wavelength) get no
  exterior load. They are named and counted as upper bounds.
- The table states its own checks:
  - probe values of every shape (a mismatch refuses);
  - midpoint samples, against which fs-vkplate measures its interpolation
    error per interval and shape, and prints each shape's own error;
  - the change between the two finest rungs of a 24/33/45-ring ladder;
  - per-family resolution of the small-ka resistance. Below it no resistance
    is applied, and the receipt prints the T60 bound.

Scope, also printed in the receipt:
- each head sits in its own infinite rigid baffle;
- no shell diffraction;
- no exterior coupling between the batter and resonant heads' `m >= 1` motion;
- no cavity loading of `m >= 1` modes.

The GPU and CUDA planners refuse a head with this load. Render it natively.

### Head viscous layer

Native card-built drums also apply the two-face Stokes boundary-layer loss
from `air_loss::plate_viscous_layer` to every retained head shape (with the
layer the load retains every shape the grid resolves, the volume shapes
(0,n>=2) included), evaluated at its loaded frequency and its own areal
density. (0,1) and the dropped shapes carry none, so their decays remain
upper bounds. Every admitted
shape must be below coincidence. This is a gas-owned loss, not a decay-time
fit; live keeps its existing installed decays.

### Snare wires

On the card route (`examples/drum_cards`, the offline examples' default) each
snare strand is a coil of sourced steel (`src/coil.rs`, ruling 7373). The
binding is `KitMaterials::with_snare_wire`. A drum without it keeps the fitted
strings and carries `kit::FITTED_WIRES_RECEIPT`: "fitted to a recording
(96ccf8f1b), not physics". That covers the presets, live, the browser and the
GPU.

- **Geometry:** Canopus Vintage Snare Wire, the maker's replica of the 1960s
  Slingerland wire. Its Japanese page gives 線径0.5mm、線外径1.35mm、ピッチ3.4mm:
  wire diameter 0.5 mm, coil outer diameter 1.35 mm, pitch 3.4 mm. The mean
  diameter is therefore 0.85 mm, and the helix angle is atan(p / (pi D)) =
  51.9 degrees.
- **Steel:**
  - `packs/steel-music-wire-string-heptonstall-ligo-p1400229.fsmatpk` is
    Heptonstall, Waller and Robertson, LIGO P1400229. It gives density
    7800 kg/m³, E 212 GPa, the three thermal coefficients, and the 1.4e-4
    bending loss tangent left after the thermoelastic peak is fitted out
    (8-3000 Hz). It is the engine-piano fixture: material state
    `steel solid music-wire-untreated 0`, at T = 293.15 K.
  - `packs/music-wire-spring-moduli-fuchs-1968.fsmatpk` is Fuchs (1968),
    Table 2. It gives the shear modulus, 79.29 GPa, compiled with material
    state `steel solid music-wire-cold-wound 0`. It is read on its untyped
    `temperature` axis at 293 K, the only temperature it holds.
- **Stiffness and mass:** over a turn, a bending moment loads the wire in
  bending and in torsion, so
  `B = 2 sin a / (cos^2 a / GJ + (1 + sin^2 a) / EI)` = 4.81e-4 N m². Wire
  torsion holds 24 % of that compliance. The mass per length is the wire's
  over sin a, 1.95e-3 kg/m; the fitted string had 4e-3. `Body::coil` adds
  `B / h^3 lap^2` to the strand's tension chain, with pinned ends. The
  bending operator therefore shares the chain's sine modes.
- **Loss** is applied exactly per mode at the mode's own discrete frequency,
  as `C = M sum 2 sigma_j phi_j phi_j^T` on the backward velocity. Its
  dissipation is counted in `Body::modal_loss_dissipated_j`. sigma is the sum
  of two parts:
  - **Material:** fs-material's exact thermoelastic series for the round wire,
    plus the card's excess loss tangent inside its band only. Both act on the
    wire-bending share of the mode's bending energy. Wire torsion loss is
    unmeasured and omitted.
  - **Gas:** the exact linear drag of the gas on the wire, across it
    (`air_loss::oscillating_cylinder`) and along it
    (`air_loss::sliding_cylinder`), turn-averaged. Linear Stokes flow
    understates the drag once the wire's amplitude nears its 0.25 mm radius,
    which a rattling wire exceeds.

  Every term is therefore a lower bound.

Contact with the head (`coil::turn_contact`, ruling 7520): each turn presses
the resonant film as an elliptical Hertz contact (`hertz::EllipticHertz`) of
the wire's surface, whose radii are the wire's 0.25 mm across it and, along
it, the helix's radius of curvature plus the wire's: (R^2 + (p/2 pi)^2)/R +
d/2 = 1.364 mm at the 51.9 degree helix. Steel E from the wire card (nu 0.3
assumed) on the film card's E and nu gives K 1.91e8 N/m^1.5 per turn on the
4.903 GPa PET; the wire's radius alone, a point contact, gives 1.18e8. A
3.5 mm node carries 1.03 turns at one approach. The film is a half-space
only while the ellipse is small against its 76 um: the semi-axes are 2.2 x
6.8 um at the 0.85 mN preload per turn, in range, and 11 x 33 um at 0.1 N per
turn. Above about 0.1 N per turn the law is outside its range. The contact's
loss, 2.0 s/m, is declared.

Still declared, not sourced, and printed with every coil: strand count,
spacing, length, tension, preload and the contact loss (`kit::Wires`). The
presets' tension, preload and contact law were adopted in a kit retuning
graded against the recording (96ccf8f1b); the card route's are 88cccbb78's
declared estimates (ruling 7486, `kit::card_drum`), but for the computed
contact stiffness.
The ends
are pinned with no strand-end or cord damping, so the lowest coil modes keep
their long material decay. On the 14x6 snare this leaves 63-125 Hz louder than
the recording.


## Numerics

Leapfrog in time. The nonlinear potential (Berger or von Karman) enters through a scalar
auxiliary variable `psi = sqrt(2V)` (Bilbao, Ducceschi and Zama): the update stays explicit,
needs one rank-one correction (Sherman-Morrison; with the air mass a 2x2 Woodbury solve), and
conserves a discrete energy, so the nonlinearity cannot destabilise the step. The plain explicit
von Karman update blew up a crash struck at 8 m/s within 7000 steps; this one stays bounded.
`psi` is advanced along the gradient only, so it drifts from `sqrt(2V)`; a shallow shell carries its
curvature stiffness in `V`, so that drift scales the shell's stiffness (on a live grid it wandered
between 0.5 and 1.5 within 2 s and the 8-16 kHz bands grew back after a second). Each step pulls `psi`
back to `sqrt(2V)` at the half step: lowering it only removes energy, and raising it spends at most
the energy the loss terms removed in the step before, so the discrete energy still never grows.
Power iteration estimates the larger of the linear and full rest-shell tangent stiffnesses.
For annuli, the rest tangent is the exact energy Hessian `h² E t Jᵀ A⁻¹ J`, where `J` is the
derivative of the weak strain load at rest. Its sparse Jacobian is assembled once per estimate;
this avoids repeated nonlinear evaluation at a tiny test displacement during instrument setup.
Rendering still uses the full nonlinear energy. The internal rate is the smallest supported
multiple of 48 kHz using at most 60 % of that estimate
(1, 2, 3, 4, 6, 8, 12, 16, 24 or 32 times); an insufficient maximum rate refuses the build.
Finer meshes can therefore use smaller steps without relaxing the stability margin;
this is not a contact or large-deflection convergence bound. Output uses a windowed-sinc filter.

The Airy operator is factored once (banded Cholesky, row-major ordering); both triangular sweeps
stream the one stored factor, which dominates a cymbal's cost on a fine grid. Contacts warm-start
their Newton solve from the previous step's force. Contacts on disjoint nodes (the snare wires, the
hi-hat rims) are solved four at a time, interleaved, each with the same operations it would take
alone, so the result is the one-at-a-time result bit for bit. A snare wire's solve stops when a
Newton step lands exactly on its own iterate (an exact fixed point) instead of bisecting away from
it: about 3 iterations per solve instead of 10, and the applied force is within a few ulps of the
root rather than 1e-13.

## Hi-hat

Two shells on one rod. The bottom shell sits on the stand's felt; the top shell hangs from the
clutch (felt washers gripping the annulus round its hole) on the pull rod, which a return spring
holds up and the pedal pulls down. Earlier constant-gap measurements sent about 62 N of an
80 N pedal load to the rims; that load split has not been requalified for the opposing geometry.
An earlier version also tied the top shell's centre to the ground through felt,
which took ~75 N of the pedal force, so the "closed" hat was barely pressed.

The shells face each other with opposite bell and bow profiles. The requested opening is the
undeformed face-to-face rim gap; it is no longer silently capped at 0.14 mm. Contact clearance
includes the local profiles, thicknesses, and up to 0.15 mm of two-lobe waviness on every grid.
Clutch axial stiffness/damping totals, angular stiffness, and rim traction are preserved under
grid refinement. A weighted affine grip separates translation, two tilts, and non-affine
deformation. Its condensed spring/dashpot form is positive semidefinite as a whole; individual
pair coefficients can be signed. It retains the 8 MN/m axial grip while allowing the cymbal
to rock against finite washer stiffness. The former point-to-rod links silently imposed an
angular clamp proportional to the squared attachment radius, which grew on coarse grids.

The angular estimate uses two circular bearings at the annulus mid-radius (13.25 mm):
each contributes `pi K_line r^3`, giving 2.412 N m/rad in total. `K_line = 1.65e5 N/m/m`
is a measured crash-washer value from [Samejima (2021)](https://doi.org/10.1250/ast.42.314),
transferred here as an estimate, not a measurement of this hi-hat or an audio fit. Angular
damping is unmeasured and set to zero. The lower cup uses one effective bearing
(1.206 N m/rad); its axial foundation coefficients (9e7 N/m^3, 2700 N s/m^3) are
integrated over the physical annulus, independently of the grid node count.
This is a small-slope impedance model, not resolved felt contact geometry.
The bodies' existing inertia carries the rocking motion.
These correct geometry and load semantics; acoustic improvement requires a separate comparison.
The offline and live constructors share `kit::HIHAT_RIM_LAW`: reference-cell
`K=2e9`, Hertz exponent 1.5 and Hunt-Crossley `mu` 0.08 s/m. Grid assembly scales
`K` with cell area. This remains a declared contact estimate, not a measured law
for these cymbals or an overall acoustic qualification.

The rims also rub (`Friction`; dry friction coefficient 0.5 in the presets, `RIM_FRICTION`, picked by
grade; 0.6 on the card route, `BRONZE_RIM_FRICTION`, below). Bending moves a plate's faces sideways by
half its thickness times its slope, so where the rims touch, the bottom shell's upper face and the top
shell's lower face slip by (t_b/2) grad w_b + (t_t/2) grad w_t, slopes from the grid's finite
differences. Each step, after the contacts, every rim contact's friction is the force that would stop
that slip (stick) unless it exceeds 0.5 times the contact's normal force. Sliding solves
`(C + lambda I) f = -d`, with `|f| = mu N`, so the force opposes the final slip even when
the tangent compliance is anisotropic; radially clipping the sticking force does not do this.
It acts through the same per-node compliances as normal contact and reduces free-mass
kinetic energy. This backward-slip discretization does not establish passivity of the
complete coupled oscillator. Normal contact reads the actual previous and predicted next
node positions, including shared-node updates, rather than a cached compression.

The centered normal solve uses `U(d) = K max(d,0)^(alpha+1)/(alpha+1)` and force
`[U(d_next)-U(d_prev)]/(d_next-d_prev) * max(1 + mu*(d_next-d_prev)/(2*dt), 0)`.
Here `d_prev` belongs to time step `n-1` and `d_next` to `n+1`; equal endpoints
use the continuous secant limit `K max(d,0)^alpha`. Stable divided differences
avoid subtracting nearly equal powers. For isolated free masses with `mu=0`,
this balances the discrete contact-potential and kinetic-energy changes through
compression and release. It is not a passivity theorem for the coupled shell,
SAV, friction and multiple-contact system.
Neither correction alone qualifies the timbre. The following grades used the earlier
geometry, gap and clock, and do not qualify the
corrected implementation. Graded against the Salamander closed hat, friction brought its octave-band decay inside
the real-vs-real ruler (0.32 octaves from the real stroke, two real strokes 0.47, 0.95 before) at
1.4 dB more octave-balance error (10.4 -> 11.8 dB: it is still far too dark, centroid 1.7 kHz against 11 kHz);
the open hat, whose rims seldom touch, is unchanged (balance 4.8 dB, decay 2.6 octaves).

## Validation (`cargo test --release`)

- clamped membrane: fundamental within 1.5 % of `2.405 c / (2 pi R)`;
- free disc: lowest elastic mode within 4 % of `lambda^2 = 5.253` (nu = 0.33);
- crash struck at 8 m/s: bounded, finite output;
- overlapping live strikes retain independent strikers; separated, receding slots are reused;
- a shallow shell's auxiliary variable keeps tracking `sqrt(2V)` (live-grid ride, 2 s);
- the stencil rows of the sparse products sum bit for bit as their CSR rows (five pieces, signed zeros);
- the four-lane contact solves match solving one contact at a time, force, penetration and every
  update bit for bit, with node-sharing and other-law contacts in the list;
- an empty matrix adds nothing, and a NaN penetration exerts no force, alone or queued;
- rim friction stops the slip exactly inside its cone, slides opposite final slip at mu N outside
  it, and contracts backward slip; prescribed cone solutions cover rotated anisotropic compliance;
- centered contact uses actual previous and predicted next positions after shared-node updates,
  in queued and scalar solves; isolated zero-loss free-mass collisions conserve discrete energy;
- every piece builds on the live grid, sounds, and the kick and snare stop being stepped once
  silent; a clone, taken settled or mid-ring, renders bit for bit what the original would (`live`).

## Real time (`live`)

`live::Live` is one piece on a grid with spacing multiplied by `scale`, struck anywhere on its
playing surface. It retains the preset's internal physics substeps and returns 48 kHz audio.
Substepped output uses the offline sinc kernel with 16 output frames of causal filter delay;
filter history and fixed 128-frame sleep windows survive arbitrary caller block sizes and clones.
It stops stepping once its motion and microphone output have decayed. On grids too coarse to put four nodes in
the 6.5-20 mm annulus round a cymbal's hole, the felt washers and the hi-hat clutch spread their
total stiffness and damping over the first ring of nodes instead. The live hi-hats now use the
same rim contact law as offline (`kit::HIHAT_RIM_LAW`, Hunt-Crossley `mu` 0.08 s/m),
replacing the stale live-only 3.0 s/m value. Their 5 mm stick tip still differs from the
offline hat's 4 mm tip. They rub with the same friction as offline. The following measurements predate
the bilinear point-contact, opposing-geometry and internal-clock corrections and do not
qualify their sound: at 4x / 4.5x the earlier closed hat's
500 Hz-1 kHz T60 falls from 1.3-2.8 s to 0.5-0.7 s, and a stroke puts
87-91 % of its energy in its first 50 ms (65-75 % without; the real closed hat 90-93 %). On grids of
22-25 mm it still carries little above ~3.5 kHz. A clone of a `Live` steps bit for bit as the original would, so the browser kit
settles each grid once: it plays a clone kept from before the timing stroke.
`apps/instrument-lab/engine-kit` compiles it for the browser.

`Live::with_hihat_pedal(scale, force_n)` constructs one hi-hat with a 2 mm rest
opening and an explicitly declared initial rod load. A nonzero load settles
through the ordinary microphone/decimator output path before playback: at least
0.5 s, then until the existing quiet gate passes, with at most 8 s extra.
The fixed closed preset uses the same preparation from its initial clock.
Nonfinite output or failure to become quiet within that budget returns an error;
preparation does not clear state or substitute temporary damping to force a pass.
`set_hihat_pedal(force_n, transition_seconds)` ramps force on that same rod;
opening or closing retains the ringing shells, contact and filter history.
The fixed `hihat-open` and `hihat-closed` presets keep their distinct rest gaps
and refuse pedal commands. These force trajectories are performance inputs.

Closed-hat rim solves now stop when a Newton update rounds to the same value;
the stick and open-hat rim policies are unchanged. In Chrome 155 on an Apple
Silicon Mac, 20 alternating pairs per grid (600 ms, 48 kHz, 128-frame blocks,
2.5 m/s at 0.42 radius) measured median paired render speedups of 1.109x at
4.5x and 1.024x at 3x. The largest PCM change was 7.46e-9; the open-hat
control was bit-identical. These measurements include bilinear point contacts.
This saves solver work; it does not fix the timbre or qualify production audio.
At 4.5x the closed hat still puts only 0.15% of its first 250 ms energy above
4 kHz, versus 91.65% in the reference, and the 3x grid still takes 719 ms to
render 600 ms. The four reference spectrogram comparisons remain failures.

The 2026-10-20 native live-grid A/B of the combined geometry/clock corrections is mixed.
At 4.5x and 2.5 m/s, closed-hat octave-balance distance to the original Salamander overhead
falls from 17.97 to 10.86 dB, but open-hat distance rises from 11.82 to 21.31 dB (first 250 ms).
Corrected open grids at 3x and 2.5x still measure 19.07 and 19.42 dB: refinement alone does
not restore the missing wash. The 3.5x crash improves in balance but worsens in attack/decay.
These are native diagnostics, not browser qualification or an overall acoustic pass; the kit
WASM package was held back. The old `rim_geometry` example assumes the pre-correction constant
gap and 48 kHz clock and remains a historical probe, not a test of the corrected geometry.

The subsequent finite-washer A/B isolates the missing rocking motion: at the true 2 mm
opening the rigid clutch produces no rim contacts and no bottom-shell energy; disabling rim
contacts gives identical PCM. The affine mount restores collisions without changing the
opening, material, axial grip or strike (2.5 m/s at 0.85 radius). Native open-hat diagnostics
against the same original overhead give the following top-mount-only isolation
(the lower cup still has its earlier distributed springs in this table):

| Grid spacing | Balance gap, first 250 ms | Energy above 4 kHz | Render duration |
| --- | --- | --- | --- |
| 24.75 mm, rigid clutch | 21.31 dB | 0.00552% | 1 s |
| 24.75 mm, affine mount | 10.40 dB | 0.954% | 1 s |
| 16.5 mm, affine mount | 13.07 dB | 2.725% | 1 s |
| 11 mm, affine mount | 11.45 dB | 2.648% | 1 s |
| 8.25 mm, affine mount | 8.30 dB | 6.790% | 0.6 s |
| 5.5 mm, affine mount | 5.91 dB | 11.412% | 0.35 s |

The real stroke has 52.693% above 4 kHz; another real stroke differs by 2.12 dB in
balance. Refinement is not monotone or shown converged, and the short fine-grid renders
cannot qualify the tail. On the 24.75 mm grid, open attack energy in the first 50 ms is
36.14% (real 38.54%), but its decay descriptor still fails. Closed-hat balance worsens
10.86 to 11.30 dB at that grid, while attack/decay improve; at 16.5 mm its balance improves
14.02 to 12.67 dB. These adverse results are retained. All 27 native tests pass, including
mount energy/passivity and axial/torque refinement regressions. This is a mechanical fix
with partial acoustic gains, not an overall fidelity pass or browser timing qualification;
the production kit remains held back.

Using finite impedance for both mounts subsequently improves the 24.75 mm open hat to
9.31 dB balance distance, 2.656% above 4 kHz, and 37.83% first-50-ms energy. Closed-hat
balance/decay worsen relative to the top-only arm (11.81 dB, 1.73 octaves), while remaining
better than the pre-mount geometry result. At 16.5 mm, both mounts give open balance
11.95 dB and 3.291% above 4 kHz; closed balance is 12.69 dB. The lower cup's physical
axial totals and both mounts' tilt torques are checked across grid sizes.
A pre-discrete-gradient fixed-settled-state rim-loss control (live `mu=3` versus offline `mu=0.08`) moves
open balance by at most 0.19 dB at these two grids and leaves its decay descriptor
unchanged. That negative result did not explain the missing wash. It predates the
centered-contact work correction and does not qualify the current shared 0.08 s/m law.
A zero-friction control also barely changes the open-hat result (at most 0.14 dB balance
at these grids). It substantially lengthens the closed attack and decay; dry friction is
retained rather than removed to improve one spectral score.

A full one-second 5.5 mm render with both finite mounts measures 6.30 dB balance distance
and 10.478% above 4 kHz. Its first-50-ms energy is 62.83% versus the reference's 38.54%,
and the 8 kHz band falls 20 dB in 170 ms versus 640 ms. Thus finer spatial resolution
does not close the tail or attack gap.

The SIMD WASM built from this source (SHA-256
`d93446abd4cf6c1b090e2328ca36d814f063dbbdf856ff7505c5f45f0bc4afd1`)
completed eight 600 ms buffered Chrome Worker checks with contiguous finite pressure,
including 96 kHz internal physics returning 48 kHz output for the 2x open hat and 7x snare.
The old level table was correctly ignored for the new hash. Browser open-hat balance
at 4.5x improves 11.60 to 9.32 dB, and energy above 4 kHz changes 2.596% to 2.661%.
This is not real-time qualification: a forced 3x live-clock probe dropped frames; buffered
render rates on that host were 0.40-0.49x at 3x and 0.086x for the 2x open hat.
The initial failed probe is retained separately from the 10-second-lead transport check.
A fresh level table and real-time qualification remain outstanding; this package is not deployed.

## Graded against a real kit

The kit's sizes follow the Salamander Drumkit (Alexander Holm, CC BY-SA 3.0): 18x18 kick,
14x6 snare, 12x7 rack tom, 14x14 floor tom, 14-inch hi-hat, Paiste 18-inch crash and 20-inch ride,
recorded with overhead microphones. `apps/instrument-lab/tools/kitgrade.py` compares every piece
band by band.

The presets (and so live, the browser and GPU) still carry the values a 2026-10-14 retuning
adopted because they graded closer to those recordings (96ccf8f1b): head tensions, strike points,
the stick's Hunt-Crossley loss, the wires' tension, preload and contact, the kick's felt crown
and pillow, the ride's shell and the hi-hat's rim loss, waviness and stick tip. Their
notes say so (`kit::FITTED_DRUM_RECEIPT`, `FITTED_CYMBAL_INTERFACE`, the ride's provenance). The
offline examples' card route keeps none of them (ruling 7486):
- **Tensions** are a drummer's tuning: both heads keep 88cccbb78's ratio and are scaled so the card
  drum's own lowest mode, read from a 0.1 m/s stroke 0.3-1.0 s after onset, sits at the late pitch
  of the Salamander strokes. That pitch is the lowest spectral peak within 20 dB of the strongest,
  30-400 Hz, where two dynamics agree: kick 40.6 Hz, snare 134.4 Hz, rack tom 133.9 Hz, floor
  tom 60.4 Hz. The batter tensions that put the card drums there are 375.5 (kick), 1072 (snare),
  1264 (rack tom) and 478 N/m (floor tom). The early pitch of a hard stroke is up to 15 Hz sharper
  (the stretch); kitgrade's "fundamental", the strongest component in the first 250 ms, reads that.
- **The wires' contact stiffness** is computed from the steel and film cards (ruling 7520, above).
- **The hi-hat's rim friction**, 0.6, is measured: Cu-10 at.% Sn (17 wt% Sn; B20 is 20 wt%) sliding
  on itself, dry in air (Buckley 1972, NASA TN D-6912, fig. 2: a polished 0.48 cm hemispherical
  rider, 250 g, 300 cm/min, 25 C, averaged over an hour; read to +/- 0.02). The other primary
  value, Morin 1834 (Tableau LXXVII), is 0.20 for an unnamed bronze sled on a bronze bed at about
  1.2 m/s, so primary values span 0.2-0.6 across contact geometry. Neither is a lathed B20 rim at
  a hi-hat's load and slip (ruling 7520, ledger
  `nasa-buckley-1972-cusn-self-friction`).
- **Everything else** is 88cccbb78's declared estimate: strike points; stick loss 0.2 s/m; wire
  tension 2.0 N, preload 0.25 N/m and contact loss 2.0 s/m; the kick's 0.12 kg
  beater with a 3 cm felt face (k 1e8 N/m^2.5, alpha 2.5, mu 0.3 s/m) and muffle dashpot 25 N s/m^3
  (no pillow is modelled);
  the ride's 1.60 mm bow, 2.20 mm bell and 3.0 s 10 kHz T60; flat hi-hat rims and a 5 mm stick
  tip; the open hat's 8 mm rim opening (`drum_cards::OPEN_HAT_OPENING_M`; the recording's pedal
  position is unknown, and the presets' 2 mm was picked by grade, bb1e5abf3).
- **The kick's microphone** (presets and card route alike, ruling 8318) is a declared overhead
  1 m above the shell's top edge, 0.1 m toward the drummer so that the batter passes the facing
  test, over a floor tangent to the hoops at a declared 0.02 m beyond the shell
  (`kit::KICK_HOOP_RISE_M`). The drum lies on its side, so the overhead and its floor image are
  both on the batter's side and the resonant head faces away from both: the overhead hears the
  batter only. The resonant head's radiation around the shell, and
  the dipole cancellation of the heads moving together, are not modelled
  (`kit::OVERHEAD_HEARS_BATTER_ONLY_RECEIPT`).
- **The hi-hat rims carry no contact loss** (`DECLARED_HIHAT_RIM_LAW`, mu 0; ruling 7865).
  88cccbb78's 3.0 s/m stood in for rim friction, which `BRONZE_RIM_FRICTION` now carries from
  data, and no source gives a separate loss at closure. So rim dissipation is the friction
  alone. The presets keep their fitted `HIHAT_RIM_LAW` (0.08 s/m).

88cccbb78's own loss points came from the recordings' band decay times. The card drums no longer
use head T60s; the cymbals' and hi-hat's T60 pairs still do.

## Estimates (not measurements)

Film and bronze properties are handbook values; head tensions are tuned; the stick's effective
tip mass (6-7 g), its Hertz constant (hickory on the target), the felt beater law (a 32 g crown in
point contact), strike points, snare-wire count, spacing, length, tension and preload (and, for the fitted strings, mass and T60), felt washer stiffness and
loss, and the hi-hat rim contact law, rest gap and waviness are declared estimates in
`src/kit.rs`. The lumped cavity has no cavity modes; there is no shell, hoop or rim
vibration, no in-plane (membrane) slip between hi-hat plates (only the bending slip rubs), and each piece is simulated alone (a kit
mix adds independent renders: the snare does not buzz when the kick is struck).

## Usage

```sh
cargo run --release --example kit -- crash crash.wav          # PIECE OUT.wav [SPEED_M_S] [SECONDS]
cargo run --release --example perform -- kick hits.csv kick.wav 17   # time_s,speed_m_s per stroke
cargo run --release --example livebench -- 5                   # real-time factor per piece at grid x5
cargo test --release
```
Live pieces: `kick snare snare-off tom-12 tom-14 hihat hihat-closed hihat-open crash ride`.
The offline `kit` and `perform` examples use the fixed hi-hat names; `kit` also renders `splash`
and `crash-shaft`, which has no default speed (give SPEED_M_S).

Their drums (`kick snare snare-off tom-12 tom-14`) are built from the sourced cards by default
(`examples/drum_cards/mod.rs`):
- the heads come from the PET film card and the stick from the hickory card, at the standard gas's
  293.15 K;
- the film loss is the film card's derived loss-tangent curve, as each head's generalized-Maxwell
  law on every grid mode (`with_film_loss_curve`); the notes carry its receipts, the source note
  `drum_cards::FILM_LOSS_SOURCE_NOTE`, and the edge loss's stated unknown;
- the stick's moving volume is `kit::STICK_TIP_KG` over the 800 kg/m³ at which that Z5A
  contact-effective mass was integrated, so the card's density sets the mass (7.30 g);
- the tip radius (5 mm) and Hunt-Crossley loss (0.2 s/m) are declared values;
- tensions, strike points, wires and the kick's beater and muffle dashpot are `kit::card_drum`'s
  (above);
- the kick keeps its spec's felt beater (`kit::kick_with_materials`) and says so, and its notes say
  that no pillow is modelled (`kit::KICK_PILLOW_RECEIPT`): only the declared muffle dashpot, an
  estimate, stands in for it.

Their cymbals and hi-hat take `kit::DECLARED_CYMBAL_INTERFACE` (and the hi-hat
`DECLARED_HIHAT_RIM_LAW`), and their ride is `kit::card_ride()`.

`--preset-heads` (anywhere on the command line) renders the recording-fitted presets instead, and
adds "heads: PRESET, recording-fitted T60 (opt-in --preset-heads); not a sourced head loss" to the
drums' notes.
