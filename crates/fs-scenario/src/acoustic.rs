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
}

impl AmbientGas {
    /// ISA sea-level air.
    #[must_use]
    pub const fn sea_level() -> Self {
        Self {
            temperature_k: 288.15,
            pressure_pa: 101_325.0,
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
}

/// Rayleigh damping coefficients for a modal family.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayleighParams {
    /// Mass-proportional coefficient α [1/s].
    pub alpha_per_s: f64,
    /// Stiffness-proportional coefficient β [s].
    pub beta_s: f64,
}

/// One cylindrical waveguide segment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderSegment {
    /// Radius [m].
    pub radius_m: f64,
    /// Axial length [m].
    pub length_m: f64,
}

/// Far-end termination of a 1D waveguide.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaveguideEnd {
    /// Rigid closure.
    Closed,
    /// Unflanged open radiation (low-`ka` fit).
    UnflangedOpen,
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
    /// Open to the exterior (`true`) or pad-closed (`false`).
    pub open: bool,
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
    /// If true and the section is isotropic, the plate is a von
    /// Karman modal pHS (`fs-nlmodal`), not a linear modal bank.
    pub geometric_nonlinearity: bool,
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
    /// Caller provenance; never invented.
    pub provenance: String,
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
