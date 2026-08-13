//! Acoustic-assembly *descriptions*: typed data that compose generic
//! physics crates. There is no instrument module here — a guitar or
//! clarinet is a description of strings, ducts, gas, and events.

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
    /// Viscous modal damping ratio (same for every mode in v1).
    pub damping_ratio: f64,
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
    /// Outlet termination.
    pub termination: WaveguideEnd,
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
    /// Optional inlet volume-velocity pulse.
    pub blow: Option<VolumeVelocityPulse>,
    /// Observer.
    pub listener: Listener,
    /// Output sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Realized duration [s].
    pub duration_s: f64,
}
