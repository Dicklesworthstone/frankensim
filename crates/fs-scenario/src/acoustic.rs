//! Acoustic-assembly *descriptions*: typed data that compose generic
//! physics. A guitar or clarinet is not a type — it is a prestressed
//! waveguide, a Bernoulli aperture, a frictional contact, a modal
//! radiator, and a viscothermal duct filled in together.

/// First-principles dry air (USSA 1976 constants live in fs-material).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AmbientGas {
    /// Static temperature [K].
    pub temperature_k: f64,
    /// Static pressure [Pa].
    pub pressure_pa: f64,
    /// Relative humidity as a fraction in `[0, 1]`. Explicit: never a
    /// hidden scenario default. ISO 9613-1 uses this; `0` is dry air.
    pub relative_humidity: f64,
}

impl AmbientGas {
    /// ISA sea-level dry air.
    #[must_use]
    pub const fn sea_level() -> Self {
        Self {
            temperature_k: 288.15,
            pressure_pa: 101_325.0,
            relative_humidity: 0.0,
        }
    }
}

/// Prestressed taut string (Kirchhoff–Carrier / modal).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrestressedString {
    /// Speaking length [m].
    pub length_m: f64,
    /// Tension [N].
    pub tension_n: f64,
    /// Linear density [kg/m].
    pub lin_density_kg_m: f64,
    /// Axial stiffness `E A` [N] (0 = linear wave equation).
    pub axial_stiffness_n: f64,
    /// Radiating width used by the compact observer [m].
    pub width_m: f64,
    /// Retained sine modes.
    pub n_modes: usize,
    /// Viscous modal damping ratio used when [`Self::rayleigh`] is `None`.
    pub damping_ratio: f64,
    /// Optional Rayleigh `ζ(ω) = α/(2ω) + βω/2` (air + internal).
    pub rayleigh: Option<RayleighParams>,
    /// Bending stiffness `E I` [N m²]. Zero is the ideal flexible string.
    /// Nonzero gives Fletcher inharmonicity `ω_n = n ω_1 √(1 + B n²)`.
    pub bending_stiffness_n_m2: f64,
    /// Second-polarization fractional detune. Zero keeps one polarization.
    /// A few 10⁻³ is a typical bridge-rocking split and produces beating.
    pub polarization_detune: f64,
    /// If true the attachment end is free (`φ(0) ≠ 0`) and realization
    /// Dirac-joins the waveguide to any plate/cavity. False is
    /// fixed-fixed sines with a one-way bridge force.
    pub moving_end: bool,
}

/// Rayleigh damping coefficients for a modal family.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayleighParams {
    /// Mass-proportional coefficient α [1/s].
    pub alpha_per_s: f64,
    /// Stiffness-proportional coefficient β [s].
    pub beta_s: f64,
}

/// One waveguide run: a cylinder when the radii match, a
/// truncated cone when they differ.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderSegment {
    /// Inlet radius [m].
    pub radius_m: f64,
    /// Axial length [m].
    pub length_m: f64,
    /// Outlet radius [m]. Equal to [`Self::radius_m`] is a cylinder.
    pub outlet_radius_m: f64,
}

impl CylinderSegment {
    /// Uniform cylinder.
    #[must_use]
    pub const fn cylinder(radius_m: f64, length_m: f64) -> Self {
        Self {
            radius_m,
            length_m,
            outlet_radius_m: radius_m,
        }
    }

    /// Linear radius taper (truncated cone).
    #[must_use]
    pub const fn taper(inlet_radius_m: f64, outlet_radius_m: f64, length_m: f64) -> Self {
        Self {
            radius_m: inlet_radius_m,
            length_m,
            outlet_radius_m,
        }
    }

    /// True when the end radii differ.
    #[must_use]
    pub fn is_taper(self) -> bool {
        (self.outlet_radius_m - self.radius_m).abs() > 1.0e-15 * (1.0 + self.radius_m)
    }
}

/// Far-end termination of a 1D waveguide.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaveguideEnd {
    /// Rigid closure.
    Closed,
    /// Unflanged open radiation (low-`ka` fit).
    UnflangedOpen,
    /// Flanged / baffled open radiation (`0.8216 a`, Rayleigh
    /// piston above the compact-`ka` ceiling).
    FlangedOpen,
}

/// Ordered 1D viscothermal duct (bore, muffler, HVAC run — same object).
#[derive(Debug, Clone, PartialEq)]
pub struct ViscothermalDuct {
    /// Inlet-first segments.
    pub segments: Vec<CylinderSegment>,
    /// Compact tone holes inserted after the named segment index.
    pub tone_holes: Vec<ToneHole>,
    /// Outlet termination.
    pub termination: WaveguideEnd,
    /// Locally reacting wall (`None` is rigid). Same mass-spring
    /// specific impedance the ODE `WallPin` is.
    pub wall: Option<LocallyReactingWall>,
}

/// Locally reacting duct wall: specific impedance
/// `Z' = r − iωσ + i K/ω` under `e^{-iωt}`. `None` on the
/// duct is a rigid wall.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocallyReactingWall {
    /// Surface density `σ` [kg/m²].
    pub surface_density_kg_m2: f64,
    /// Stiffness per unit area `K` [Pa/m].
    pub stiffness_pa_per_m: f64,
    /// Specific resistance `r` [Pa s / m].
    pub resistance_pa_s_per_m: f64,
}

/// One compact-limit tone hole on a 1D bore.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToneHole {
    /// Insert after this inlet-first cylinder index.
    pub after_segment: usize,
    /// Chimney radius [m].
    pub radius_m: f64,
    /// Chimney height [m].
    pub chimney_m: f64,
    /// Opening in `[0, 1]`. `0` is a sealed pad, `1` is fully
    /// open, and a fraction is the admittance mix
    /// `Y = σ Y_open + (1−σ) Y_closed` already used by the
    /// TMM `HoleState::Vent` and the ODE `AcousticTap`.
    pub open_fraction: f64,
}

/// Instantaneous point displacement of a taut string (the static
/// triangular shape of a pluck, not a sampled envelope).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pluck {
    /// Station as a fraction of speaking length in (0, 1).
    pub station_frac: f64,
    /// Transverse height [m].
    pub height_m: f64,
}

/// Prescribed inlet volume-velocity pulse (a blow, a loudspeaker, a fan).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumeVelocityPulse {
    /// Peak volume velocity [m³/s].
    pub peak_m3_s: f64,
    /// Pulse duration [s].
    pub duration_s: f64,
}

/// Compact listener placement for the radiation observer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Listener {
    /// Distance from the source compact origin [m].
    pub distance_m: f64,
}

/// A pressure-driven Bernoulli aperture (reed, vocal fold, valve).
///
/// Closing pressure is `k H / S` in the massless limit: the rest
/// opening shuts when the pressure drop reaches
/// [`Self::closing_pressure_pa`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeatingReed {
    /// Rest opening `H` [m].
    pub rest_opening_m: f64,
    /// Slit width [m].
    pub width_m: f64,
    /// Pressure drop that just closes the reed [Pa].
    pub closing_pressure_pa: f64,
    /// Blowing (mouth) pressure [Pa].
    pub blowing_pressure_pa: f64,
    /// Raised-cosine attack to the blowing pressure [s].
    pub attack_s: f64,
    /// Effective reed mass [kg]. Zero selects the quasistatic valve.
    pub mass_kg: f64,
    /// Reed stiffness [N/m]. Ignored when `mass_kg == 0`.
    pub stiffness_n_m: f64,
}

/// Regularized Stribeck friction at a station (bow, brake, fault).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BowStroke {
    /// Contact station as a fraction of speaking length in (0, 1).
    pub station_frac: f64,
    /// Normal force [N].
    pub normal_force_n: f64,
    /// Bow velocity [m/s] (signed).
    pub velocity_m_s: f64,
    /// Static friction coefficient.
    pub mu_static: f64,
    /// Dynamic friction coefficient.
    pub mu_dynamic: f64,
    /// Stribeck velocity scale [m/s].
    pub stribeck_m_s: f64,
}

/// A declared one-dimensional contact-path height spectrum.
///
/// This is surface geometry, not a bow or a rosin curve. Realization
/// asks `fs-tribo::surface_excitation` to sample it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactTexture {
    /// RMS height of the periodic profile [m].
    pub rms_height_m: f64,
    /// One-dimensional self-affine Hurst exponent in (0, 1).
    pub hurst_exponent: f64,
    /// Inclusive lowest retained spatial cycle.
    pub min_cycles: u32,
    /// Inclusive highest retained spatial cycle.
    pub max_cycles: u32,
    /// Deterministic phase seed.
    pub phase_seed: u64,
    /// Material-frame track length [m].
    pub track_length_m: f64,
    /// Linearized contact stiffness `dF_n / dh` [N/m].
    pub tangent_stiffness_n_m: f64,
}

/// A thin orthotropic plate (a panel, a soundboard, a bulkhead).
///
/// Modes are not data. Realization asks `fs-plate` + `fs-modal` for
/// certified eigenpairs and radiates those. Music is not a special case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThinPlate {
    /// Side length along material axis 1 [m].
    pub length_m: f64,
    /// Side length along material axis 2 [m].
    pub width_m: f64,
    /// Thickness [m].
    pub thickness_m: f64,
    /// Density [kg/m³].
    pub density_kg_m3: f64,
    /// Young's modulus along axis 1 [Pa].
    pub e1_pa: f64,
    /// Young's modulus along axis 2 [Pa].
    pub e2_pa: f64,
    /// Poisson ratio ν12.
    pub nu12: f64,
    /// In-plane shear modulus G12 [Pa].
    pub g12_pa: f64,
    /// Viscous modal damping ratio.
    pub damping_ratio: f64,
    /// How many certified modes to keep.
    pub n_modes: usize,
    /// If true the plate is a von Karman modal pHS (`fs-nlmodal`),
    /// not a linear modal bank. Isotropic SS uses sine modes;
    /// clamped or orthotropic bending uses FE-sampled displacement.
    pub geometric_nonlinearity: bool,
    /// Isotropic in-plane pretension [N/m]. Zero is unloaded.
    pub pretension_n_m: f64,
    /// If true, every edge is clamped (`w = ∇w = 0`). False is
    /// simply supported. Clamped von Karman samples DKT `w` and
    /// keeps the sine Airy membrane channel.
    pub clamped: bool,
}

/// A distributed unilateral obstacle under a taut span (a fretboard,
/// a reed lay, a snare wire, a cable against a stay).
#[derive(Debug, Clone, PartialEq)]
pub struct UnilateralObstacle {
    /// Sample stations as a fraction of span in (0, 1).
    pub stations: Vec<f64>,
    /// Gaps from the rest span to the obstacle [m]. Positive is clearance.
    pub gaps_m: Vec<f64>,
    /// Contact stiffness `K`.
    pub stiffness: f64,
    /// Power-law exponent `α >= 1`.
    pub alpha: f64,
    /// Kinetic friction coefficient at the contact. Zero is frictionless
    /// (normal-only). Nonzero composes `fs-tribo` on the same stations.
    pub mu_kinetic: f64,
    /// Hunt–Crossley internal-loss coefficient `χ` [s/m]. Zero is the
    /// elastic potential (restitution 1). Nonzero is a dissipative
    /// port force, not a gradient of `H`.
    pub internal_loss: f64,
    /// Caller provenance; never invented.
    pub provenance: String,
}

/// A lumped Helmholtz volume facing a radiating panel.
///
/// Geometry is data. A guitar body, a vented box, and a bottle are
/// the same object; realization composes `fs-phs` flow-driven
/// Helmholtz with the plate monopoles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelmholtzCavity {
    /// Cavity volume [m³].
    pub volume_m3: f64,
    /// Neck radius [m].
    pub neck_radius_m: f64,
    /// Physical neck length [m] (end correction is applied in the pHS).
    pub neck_length_m: f64,
}

/// One compact modal monopole (a body mode, a panel, a loudspeaker).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiatingPlate {
    /// Radiating area [m²].
    pub area_m2: f64,
    /// Effective modal mass [kg].
    pub mass_kg: f64,
    /// Natural frequency [Hz].
    pub frequency_hz: f64,
    /// Viscous damping ratio.
    pub damping_ratio: f64,
}

/// A described acoustic assembly. Empty of both string and duct is
/// refused at realization, not here — this type is data.
#[derive(Debug, Clone, PartialEq)]
pub struct AcousticAssembly {
    /// Ambient gas (T, p). Composition is dry air in v1.
    pub ambient: AmbientGas,
    /// Optional prestressed string.
    pub string: Option<PrestressedString>,
    /// Optional viscothermal duct.
    pub duct: Option<ViscothermalDuct>,
    /// Optional string pluck.
    pub pluck: Option<Pluck>,
    /// Optional bow stroke (may replace or accompany a pluck IC).
    pub bow: Option<BowStroke>,
    /// Optional inlet volume-velocity pulse.
    pub blow: Option<VolumeVelocityPulse>,
    /// Optional beating reed (replaces a prescribed blow on a duct).
    pub reed: Option<BeatingReed>,
    /// Optional radiating body driven by the string bridge force.
    pub soundboard: Option<RadiatingPlate>,
    /// Extra body modes (top + Helmholtz, bracing, …) driven by the
    /// same bridge force as [`Self::soundboard`].
    pub body_modes: Vec<RadiatingPlate>,
    /// Optional thin plate whose certified modes replace a named-Hertz
    /// soundboard. Driven by bridge force and, when a duct is present,
    /// by mouth pressure (structure–bore).
    pub plate: Option<ThinPlate>,
    /// Optional lumped Helmholtz volume facing the plate monopoles.
    /// `None` is a plate in free half-space.
    pub cavity: Option<HelmholtzCavity>,
    /// Distributed obstacles under the string (rattle, frets, a stay).
    pub obstacles: Vec<UnilateralObstacle>,
    /// Optional declared contact-path texture. When a bow is present
    /// this height trace modulates the normal load that friction sees.
    pub contact_texture: Option<ContactTexture>,
    /// Observer.
    pub listener: Listener,
    /// Output sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Realized duration [s].
    pub duration_s: f64,
}
