//! Free modes of the already selected causal shell material image, exp(s t).
//!
//! The original mass-normalized equilibrium coordinates, complete background
//! and interface viscosity, and EVERY selected substrate/film memory row enter
//! the existing fs-phs storage owner. Its shared Newton/Beyn analysis supplies
//! the poles. No modal damping fit, frequency freeze, new material lowering,
//! branch merging or extra eigenvalue algorithm occurs here.
//!
//! The existing film Prony lowering and constant-loss substrate continuation
//! remain approximations to their source laws. Selecting them authors a global
//! passive model; its complex roots are not experimental complex-frequency
//! validation. Frozen films and unselected anchored substrate loss refuse.
//! The original finite source band is retained beside the results. The shell
//! basis is unchanged, with no new spatial/modal convergence claim.

use super::harmonic::{Image, projection, quadratic};
use super::strain::CoatingStrainSource;
use crate::acoustic_realize::AcousticRealizeError;
use crate::shear_shell::{
    MaterialShellImage, ResolvedShearShell, assembly::MaterialShellAssemblyImage, poll,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicModel, HarmonicPole, HarmonicPoleError,
    HarmonicSpectrumError, MAX_HARMONIC_DIMENSION,
};
use fs_phs::{DescriptorPortHamiltonian, PhsError, PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage};

pub use super::harmonic::loaded::ShellHarmonicPort;
pub use fs_phs::harmonic::{BeynOptions, PoleOptions, ResonanceRegion};

/// Full-state forced response of the same causal shell and independent loads.
pub mod response;

/// An original causal impedance owner, connected through physical shell flows.
/// The load receives the selected velocity/volume flow and returns its opposing
/// force/pressure. Its storage must be declared exactly quadratic; all of its
/// states and mutual terms are retained. No isolated load impedance is inverted.
pub struct ShellResonanceLoad<'a> {
    pub system: PortHamiltonian,
    /// One original shell projection for every load input, in authored order.
    /// Signed, repeated and nonorthogonal projections are retained unchanged.
    pub ports: &'a [ShellHarmonicPort],
    /// Real-frequency support of the authored load [rad/s]. Selecting this
    /// continuous owner also selects its global causal continuation. Results
    /// retain the common source support without claiming measured complex poles.
    pub band_rad_s: (f64, f64),
}

/// An original constrained impedance, connected to physical shell flows.
/// Every differential state and algebraic reaction remains in the joint model.
/// Its exact quadratic storage acts only on differential coordinates; reaction
/// coordinates enter port effort and resistance but store no energy.
pub struct ShellResonanceDescriptorLoad<'a> {
    pub system: DescriptorPortHamiltonian,
    /// One original signed shell force/volume-flow projection per load input.
    pub ports: &'a [ShellHarmonicPort],
    /// Original real-frequency support of the explicitly selected continuation.
    pub band_rad_s: (f64, f64),
}

enum ShellLoadSource<'a> {
    Ode(ShellResonanceLoad<'a>),
    Descriptor(ShellResonanceDescriptorLoad<'a>),
}
impl ShellLoadSource<'_> {
    fn state_dim(&self) -> usize {
        match self {
            Self::Ode(load) => load.system.state_dim(),
            Self::Descriptor(load) => load.system.state_dim(),
        }
    }
    fn differential_dim(&self) -> usize {
        match self {
            Self::Ode(load) => load.system.state_dim(),
            Self::Descriptor(load) => load.system.differential_dim(),
        }
    }
    fn port_dim(&self) -> usize {
        match self {
            Self::Ode(load) => load.system.port_dim(),
            Self::Descriptor(load) => load.system.port_dim(),
        }
    }
    fn ports(&self) -> &[ShellHarmonicPort] {
        match self {
            Self::Ode(load) => load.ports,
            Self::Descriptor(load) => load.ports,
        }
    }
    fn band(&self) -> (f64, f64) {
        match self {
            Self::Ode(load) => load.band_rad_s,
            Self::Descriptor(load) => load.band_rad_s,
        }
    }
    fn storage(&self) -> &dyn Storage {
        match self {
            Self::Ode(load) => load.system.storage(),
            Self::Descriptor(load) => load.system.storage(),
        }
    }
    fn structure(&self) -> (&[f64], &[f64], &[f64]) {
        match self {
            Self::Ode(load) => load.system.structure(),
            Self::Descriptor(load) => load.system.structure(),
        }
    }
}

/// Original dynamic load coordinates and separate envelope energy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellResonanceLoadResponse {
    /// Original differential states followed by any algebraic reactions.
    /// Prepared source ranges identify the two blocks without condensation.
    pub state: Vec<C64>,
    /// Work-conjugate shell flow and load opposing effort, in load port order.
    pub flow: Vec<C64>,
    pub effort: Vec<C64>,
    /// Re(conj(flow)*effort) at each original port; positive into the load.
    pub absorbed_port_power_w: Vec<f64>,
    /// Original differential quadratic energy; reactions store none.
    pub stored_energy_j: f64,
    /// Full original resistance evaluated on [Q_diff*z_diff, z_alg].
    pub dissipated_power_w: f64,
    /// 2 Re(s) H_load + D_load - sum(absorbed_port_power_w).
    pub decay_power_residual_w: f64,
    /// Independent componentwise residual of every original load state row.
    pub relative_equation_residual: f64,
}

/// Cold complete-state preparation, including the shared physical admissions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellResonanceBudget {
    /// q,p for EVERY retained mode, every selected material memory row and
    /// every original load coordinate, including algebraic reactions.
    pub maximum_coordinates: usize,
    /// Conservative simultaneous numerical payload; excludes borrowed source
    /// storage, allocator overhead and process RSS.
    pub maximum_payload_bytes: usize,
    /// Conservative scalar work including shared PSD admissions and generator
    /// preparation. The eigensolver's per-query budgets are separate.
    pub maximum_preparation_work: u64,
    /// Original componentwise equation and independent energy/loss gate.
    pub physical_relative_tolerance: f64,
}

/// Original complete mode and source-specific physical observations.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellResonancePole {
    /// Complete fs-phs mode. State order is q0,p0,q1,p1,..., then original
    /// substrate z rows, then original coating z rows, then all load states.
    /// Momentum p equals modal
    /// velocity because these unchanged coordinates have unit mass.
    /// The numerical vector normalization does not admit a finite-amplitude
    /// material history or establish that this mode is acoustically observable.
    pub mode: HarmonicPole,
    pub displacement_m_sqrt_kg: Vec<C64>,
    pub velocity_m_sqrt_kg_per_s: Vec<C64>,
    /// Original physical force/displacement rows, in source port order.
    pub force_port_displacement_m: Vec<C64>,
    /// Mechanical kinetic and equilibrium strain envelope energy [J].
    pub mechanical_storage_j: f64,
    pub background_dissipation_w: f64,
    pub interface_dissipation_w: f64,
    /// In original selected substrate-row order: |c.q-z|²/2 and |c.q-z|²/tau.
    pub substrate_branch_storage_j: Vec<f64>,
    pub substrate_branch_dissipation_w: Vec<f64>,
    /// In original selected coating-row order. Equal-time factors can mix
    /// source films; these are branch accounts, not invented per-film shares.
    pub coating_branch_storage_j: Vec<f64>,
    pub coating_branch_dissipation_w: Vec<f64>,
    /// Present only for a complete dynamic load; no load state is condensed.
    pub load: Option<ShellResonanceLoadResponse>,
    /// Shell-only account: 2 Re(s) H_shell + D_shell + power_into_load.
    pub shell_decay_power_residual_w: f64,
    /// Complete joint envelope account: 2 Re(s) H + D. This is not a
    /// finite-cycle time average for a decaying real trajectory.
    pub decay_power_residual_w: f64,
    /// Maximum relative discrepancy of independent physical accounting.
    pub physical_relative_residual: f64,
    /// |Im(s)| lies in the original real-frequency support. Even true does
    /// not validate continuation to this complex frequency experimentally.
    /// A pure real relaxation pole has no oscillation and reports false.
    pub oscillation_within_source_band: bool,
}

impl ShellResonancePole {
    /// Original work-conjugate force-port velocities [m/s].
    pub fn force_port_velocity_m_s(&self) -> &[C64] {
        &self.mode.port_output[..self.force_port_displacement_m.len()]
    }
    /// Signed exterior normal-volume flow [m³/s]. A zero value is not silence.
    pub fn exterior_volume_velocity_m3_s(&self) -> C64 {
        self.mode.port_output[self.force_port_displacement_m.len()]
    }
    /// Independently signed outward-from-cavity volume flow [m³/s].
    pub fn cavity_volume_velocity_m3_s(&self) -> C64 {
        self.mode.port_output[self.force_port_displacement_m.len() + 1]
    }
}

/// Shared numerical contour discovery with the same physical mode observers.
/// Numerical count/rank agreement is not a rigorous completeness enclosure.
#[derive(Clone, Debug)]
pub struct ShellResonanceSpectrum {
    pub poles: Vec<ShellResonancePole>,
    pub unrefined_rates_per_s: Vec<C64>,
    pub argument_count: usize,
    pub resolved_rank: usize,
    pub quadrature_nodes: usize,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ShellResonanceError {
    Source(AcousticRealizeError),
    PhysicalModel(PhsError),
    Preparation(HarmonicError),
    Pole(HarmonicPoleError),
    Spectrum(HarmonicSpectrumError),
    Invalid(&'static str),
    Budget {
        quantity: &'static str,
        requested: u128,
        cap: u128,
    },
    PhysicalAccounting {
        relative: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for ShellResonanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SHELL-RESONANCE: {self:?}")
    }
}
impl std::error::Error for ShellResonanceError {}
impl From<AcousticRealizeError> for ShellResonanceError {
    fn from(error: AcousticRealizeError) -> Self {
        Self::Source(error)
    }
}

fn cap(quantity: &'static str, requested: u128, limit: u128) -> Result<(), ShellResonanceError> {
    if requested > limit {
        Err(ShellResonanceError::Budget {
            quantity,
            requested,
            cap: limit,
        })
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Prepared once from an immutable shell or monolithic shell-assembly image.
/// All spectral work delegates to the complete fs-phs continuous model.
pub struct PreparedShellResonance<'a> {
    image: Image<'a>,
    // Same immutable original geometry, physical modes and authored film limits.
    strain_source: CoatingStrainSource<'a>,
    model: HarmonicModel,
    dimension: usize,
    substrate_rows: usize,
    coating_rows: usize,
    load: Option<PreparedLoad>,
    band: (f64, f64),
    tolerance: f64,
}

struct PreparedLoad {
    states: usize,
    differential_states: usize,
    ports: usize,
    /// Original port-major physical shell projections.
    rows: Vec<f64>,
    q: Vec<f64>,
    structure: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
}

fn source_band(
    source: &ResolvedShearShell,
    gate: &CancelGate,
) -> Result<(f64, f64), ShellResonanceError> {
    poll(gate)?;
    if !source.coatings().is_empty() && source.coating_relaxation_options().is_none() {
        return Err(ShellResonanceError::Invalid(
            "shell poles require the selected causal coating continuation, not frozen films",
        ));
    }
    if source.substrate_relaxation_options().is_none()
        && source
            .regions()
            .iter()
            .any(|r| r.constant_loss_reference_omega_rad_s().is_some())
    {
        return Err(ShellResonanceError::Invalid(
            "anchored substrate loss requires its selected causal relaxation before pole analysis",
        ));
    }
    let mut band = (0.0_f64, f64::MAX);
    for b in source
        .regions()
        .iter()
        .map(|r| r.omega_band_rad_s())
        .chain(source.interfaces().iter().map(|r| r.omega_band_rad_s()))
        .chain(source.coatings().iter().map(|r| r.omega_band_rad_s()))
    {
        poll(gate)?;
        band.0 = band.0.max(b.0);
        band.1 = band.1.min(b.1);
    }
    if band.0 > band.1 {
        return Err(ShellResonanceError::Invalid(
            "shell source bands do not intersect",
        ));
    }
    Ok(band)
}

impl MaterialShellImage {
    /// Prepare poles of this SAME equilibrium-basis image and its complete
    /// selected memory. Bare Kelvin-Voigt images are supported. Films and any
    /// anchored constant-loss substrate require their explicit existing causal
    /// selections. No accepted motion, strain history or material fit is changed.
    pub fn prepare_coating_resonance(
        &self,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        self.prepare_coating_resonance_inner(None, budget, gate)
    }

    /// Couple the same causal shell image to every original dynamic load state.
    /// The selected force/pressure reactions enter the joint generator before
    /// any pole refinement or contour solve, including at a bare load resonance.
    pub fn prepare_coating_resonance_with_load(
        &self,
        load: ShellResonanceLoad<'_>,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        self.prepare_coating_resonance_inner(Some(ShellLoadSource::Ode(load)), budget, gate)
    }

    /// Retain an original constrained attachment beside every selected material
    /// memory coordinate. The shared descriptor owner keeps all constraint
    /// reactions with zero storage and no artificial inertia or regularization.
    /// Finite poles and forced responses use the same complete source model;
    /// descriptor index and infinite-spectrum completeness are not certified.
    pub fn prepare_coating_resonance_with_descriptor_load(
        &self,
        load: ShellResonanceDescriptorLoad<'_>,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        self.prepare_coating_resonance_inner(
            Some(ShellLoadSource::Descriptor(load)), budget, gate,
        )
    }

    fn prepare_coating_resonance_inner(
        &self,
        load: Option<ShellLoadSource<'_>>,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        let band = source_band(self.source(), gate)?;
        if self.coating_relaxation().map_or(0, |m| m.films().len())
            != self.source().coatings().len()
        {
            return Err(ShellResonanceError::Invalid(
                "every original coating projection is required",
            ));
        }
        PreparedShellResonance::new(
            Image::shell(self),
            CoatingStrainSource::Shell(self),
            band,
            load,
            budget,
            gate,
        )
    }
}
impl MaterialShellAssemblyImage {
    /// The same complete material dynamics on one monolithic tied-body basis;
    /// component interface and assembly glue viscosity enter exactly once.
    pub fn prepare_coating_resonance(
        &self,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        self.prepare_coating_resonance_inner(None, budget, gate)
    }

    /// Attach one complete multiport dynamic load to the monolithic tied-body
    /// image. Component and glue losses remain distinct from load dissipation.
    pub fn prepare_coating_resonance_with_load(
        &self,
        load: ShellResonanceLoad<'_>,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        self.prepare_coating_resonance_inner(Some(ShellLoadSource::Ode(load)), budget, gate)
    }

    /// Retain an original constrained attachment beside every selected material
    /// memory coordinate. The shared descriptor owner keeps all constraint
    /// reactions with zero storage and no artificial inertia or regularization.
    /// Finite poles and forced responses use the same complete source model;
    /// descriptor index and infinite-spectrum completeness are not certified.
    pub fn prepare_coating_resonance_with_descriptor_load(
        &self,
        load: ShellResonanceDescriptorLoad<'_>,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        self.prepare_coating_resonance_inner(
            Some(ShellLoadSource::Descriptor(load)), budget, gate,
        )
    }

    fn prepare_coating_resonance_inner(
        &self,
        load: Option<ShellLoadSource<'_>>,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellResonance<'_>, ShellResonanceError> {
        let mut band = (0.0_f64, f64::MAX);
        let mut films = 0usize;
        for part in &self.source().components {
            let b = source_band(&part.specimen, gate)?;
            band.0 = band.0.max(b.0);
            band.1 = band.1.min(b.1);
            films += part.specimen.coatings().len();
        }
        if self.coating_relaxation().map_or(0, |m| m.films().len()) != films {
            return Err(ShellResonanceError::Invalid(
                "every component coating projection is required",
            ));
        }
        PreparedShellResonance::new(
            Image::assembly(self),
            CoatingStrainSource::Assembly(self),
            band,
            load,
            budget,
            gate,
        )
    }
}

impl<'a> PreparedShellResonance<'a> {
    fn new(
        image: Image<'a>,
        strain_source: CoatingStrainSource<'a>,
        mut band: (f64, f64),
        load: Option<ShellLoadSource<'_>>,
        budget: ShellResonanceBudget,
        gate: &CancelGate,
    ) -> Result<Self, ShellResonanceError> {
        poll(gate)?;
        let m = image.modes.len();
        let substrate_rows = image.substrate.map_or(0, |s| s.relaxation_times_s().len());
        let coating_rows = image.films.map_or(0, |s| s.relaxation_times_s().len());
        let shell_dimension = 2 * m + substrate_rows + coating_rows;
        let load_dimension = load.as_ref().map_or(0, |load| load.state_dim());
        let load_ports = load.as_ref().map_or(0, |load| load.port_dim());
        let dimension =
            shell_dimension
                .checked_add(load_dimension)
                .ok_or(ShellResonanceError::Invalid(
                    "joint shell/load dimension overflow",
                ))?;
        if let Some(load) = &load {
            let load_band = load.band();
            if load_ports == 0
                || load.ports().len() != load_ports
                || !load_band.0.is_finite()
                || load_band.0 < 0.0
                || !load_band.1.is_finite()
                || load_band.1 < load_band.0
            {
                return Err(ShellResonanceError::Invalid(
                    "one physical shell port per load input and an ordered finite load band required",
                ));
            }
            band.0 = band.0.max(load_band.0);
            band.1 = band.1.min(load_band.1);
        }
        let ports = image.force_ports.len() + 2;
        let n = dimension as u128;
        // Covers physical J/R/Q admissions (including the shared bounded
        // Jacobi PSD checks), complete memory attachment and model assembly.
        let all_ports = ports as u128 + load_ports as u128;
        let bytes = 256 * n * n + 64 * n * all_ports + 256 * (n + all_ports);
        let cubic = if load.is_some() { 8192 } else { 2048 };
        let work =
            cubic * n * n * n + 64 * n * n * (all_ports + image.interfaces.len() as u128 + 8);
        if m == 0
            || m > crate::shear_shell::MAX_SHELL_MODES
            || dimension > MAX_HARMONIC_DIMENSION
            || band.0 > band.1
            || !budget.physical_relative_tolerance.is_finite()
            || budget.physical_relative_tolerance <= 0.0
            || budget.physical_relative_tolerance >= 1.0
            || image.background.len() != m * m
            || image.exterior.len() != m
            || image.cavity.len() != m
            || image.force_ports.iter().any(|p| p.len() != m)
            || image
                .interfaces
                .iter()
                .any(|(p, c)| p.len() != m || !c.is_finite() || *c < 0.0)
        {
            return Err(ShellResonanceError::Invalid(
                "shell resonance needs a complete bounded image and positive finite accuracy",
            ));
        }
        cap(
            "complete coordinates",
            n,
            budget.maximum_coordinates as u128,
        )?;
        cap(
            "preparation payload bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "preparation work",
            work,
            budget.maximum_preparation_work as u128,
        )?;
        let base = 2 * m;
        let mut storage = vec![0.0; base * base];
        let mut j = vec![0.0; base * base];
        let mut r = vec![0.0; base * base];
        let mut g = vec![0.0; base * ports];
        for (i, mode) in image.modes.iter().enumerate() {
            poll(gate)?;
            let (qi, pi) = (2 * i, 2 * i + 1);
            storage[qi * base + qi] = mode.omega * mode.omega;
            storage[pi * base + pi] = 1.0;
            j[qi * base + pi] = 1.0;
            j[pi * base + qi] = -1.0;
            for k in 0..m {
                r[pi * base + 2 * k + 1] = image.background[i * m + k];
            }
            for (p, row) in image.force_ports.iter().enumerate() {
                g[pi * ports + p] = row[i];
            }
            g[pi * ports + ports - 2] = image.exterior[i];
            g[pi * ports + ports - 1] = image.cavity[i];
        }
        for &(row, damping) in &image.interfaces {
            for (i, &a) in row.iter().enumerate() {
                poll(gate)?;
                for (k, &b) in row.iter().enumerate() {
                    r[(2 * i + 1) * base + 2 * k + 1] += (a * b) * damping;
                }
            }
        }
        if storage.iter().chain(&r).chain(&g).any(|v| !v.is_finite()) {
            return Err(ShellResonanceError::Invalid(
                "shell physical coefficients are not representable",
            ));
        }
        let storage =
            QuadraticStorage::new(storage, base).map_err(ShellResonanceError::PhysicalModel)?;
        poll(gate)?;
        let system = PortHamiltonian::new(base, ports, j, r, g, Box::new(storage))
            .map_err(ShellResonanceError::PhysicalModel)?;
        poll(gate)?;
        let mut branches = Vec::with_capacity(substrate_rows + coating_rows);
        for (rows, times) in memories(&image) {
            if rows.len() != times.len() {
                return Err(ShellResonanceError::Invalid(
                    "complete material rows and times must match",
                ));
            }
            for (row, &time) in rows.iter().zip(times) {
                poll(gate)?;
                if row.len() != m {
                    return Err(ShellResonanceError::Invalid(
                        "material rows must retain the full shell basis",
                    ));
                }
                let mut projection = vec![0.0; base];
                for (i, &c) in row.iter().enumerate() {
                    projection[2 * i] = c;
                }
                branches.push(RelaxationBranch {
                    projection,
                    stiffness: 1.0,
                    relaxation_time_s: time,
                });
            }
        }
        // This is the SAME internal-strain storage owner used by playback.
        let system = system
            .with_relaxation_branches(branches)
            .map_err(ShellResonanceError::PhysicalModel)?;
        poll(gate)?;
        let harmonic_budget = HarmonicBudget {
            max_dimension: budget.maximum_coordinates,
            max_dense_bytes: budget.maximum_payload_bytes,
            max_dense_work: budget.maximum_preparation_work,
            relative_residual_tolerance: budget.physical_relative_tolerance,
        };
        let (model, load) = if let Some(load) = load {
            let q = load
                .storage()
                .quadratic_part()
                .filter(|q| q.exact)
                .ok_or(ShellResonanceError::Invalid(
                    "dynamic shell loads require declared exactly quadratic storage",
                ))?;
            if q.matrix.len() != load_dimension * load_dimension
                || q.matrix.iter().any(|value| !value.is_finite())
            {
                return Err(ShellResonanceError::Invalid(
                    "complete finite original load storage matrix required",
                ));
            }
            let (j, r, g) = load.structure();
            let mut rows = Vec::with_capacity(load_ports * m);
            let mut columns = vec![0.0; shell_dimension * load_ports];
            for (p, port) in load.ports().iter().enumerate() {
                poll(gate)?;
                let row = match port {
                    ShellHarmonicPort::Force(index) => image
                        .force_ports
                        .get(*index)
                        .map(|row| row.as_slice())
                        .ok_or(ShellResonanceError::Invalid(
                            "unknown original shell force port",
                        ))?,
                    ShellHarmonicPort::ExteriorVolume => image.exterior,
                    ShellHarmonicPort::CavityVolume => image.cavity,
                };
                rows.extend_from_slice(row);
                for (i, &coefficient) in row.iter().enumerate() {
                    columns[(2 * i + 1) * load_ports + p] = coefficient;
                }
            }
            let prepared = PreparedLoad {
                states: load_dimension,
                differential_states: load.differential_dim(),
                ports: load_ports,
                rows,
                q: q.matrix,
                structure: j.iter().zip(r).map(|(j, r)| j - r).collect(),
                r: r.to_vec(),
                g: g.to_vec(),
            };
            poll(gate)?;
            let model = match load {
                ShellLoadSource::Ode(load) => {
                    let system = system
                        .with_parallel_load_columns(load.system, &columns)
                        .map_err(ShellResonanceError::PhysicalModel)?;
                    poll(gate)?;
                    HarmonicModel::from_ode(&system, harmonic_budget, || gate.is_requested())
                }
                ShellLoadSource::Descriptor(load) => {
                    let system = system
                        .with_parallel_descriptor_load_columns(load.system, &columns)
                        .map_err(ShellResonanceError::PhysicalModel)?;
                    poll(gate)?;
                    HarmonicModel::from_descriptor(&system, harmonic_budget, || gate.is_requested())
                }
            }.map_err(ShellResonanceError::Preparation)?;
            (model, Some(prepared))
        } else {
            (HarmonicModel::from_ode(&system, harmonic_budget, || gate.is_requested())
                .map_err(ShellResonanceError::Preparation)?, None)
        };
        poll(gate)?;
        Ok(Self {
            image,
            strain_source,
            model,
            dimension,
            substrate_rows,
            coating_rows,
            load,
            band,
            tolerance: budget.physical_relative_tolerance,
        })
    }

    pub const fn dimension(&self) -> usize {
        self.dimension
    }
    pub fn mode_count(&self) -> usize {
        self.image.modes.len()
    }
    /// Original equilibrium-basis frequencies, useful as mechanical seeds;
    /// these are not the final hereditary complex poles.
    pub fn equilibrium_angular_frequencies_rad_s(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.image.modes.iter().map(|mode| mode.omega)
    }
    pub const fn substrate_memory_count(&self) -> usize {
        self.substrate_rows
    }
    pub fn coating_memory_count(&self) -> usize {
        self.coating_rows
    }
    /// Every original load coordinate, including zero-storage reactions.
    pub fn load_state_count(&self) -> usize {
        self.load.as_ref().map_or(0, |load| load.states)
    }
    /// Original differential load suffix before its algebraic reactions.
    pub fn load_differential_range(&self) -> Option<core::ops::Range<usize>> {
        self.load.as_ref().map(|load| {
            let start = self.dimension - load.states;
            start..start + load.differential_states
        })
    }
    /// Original zero-storage reactions; empty for an ordinary ODE load.
    pub fn load_algebraic_range(&self) -> Option<core::ops::Range<usize>> {
        self.load.as_ref().map(|load| {
            self.dimension - load.states + load.differential_states..self.dimension
        })
    }
    /// Real-axis support retained from the original sources, not an enclosure
    /// of the selected globally passive continuation's complex-frequency error.
    pub const fn source_band_rad_s(&self) -> (f64, f64) {
        self.band
    }

    /// Lift a mechanical seed by p=s*q and z=c.q/(1+s*tau). Near an exact
    /// material relaxation denominator supply a complete independent state seed
    /// to `resonance_pole` instead; no memory direction is removed from it.
    /// Independently dynamic load states start at zero; callers can replace the
    /// appended entries to select a load-dominated mode, or use contour discovery.
    pub fn displacement_seed(
        &self,
        rate_per_s: C64,
        displacement: &[C64],
        gate: &CancelGate,
    ) -> Result<Vec<C64>, ShellResonanceError> {
        poll(gate)?;
        if !finite(rate_per_s)
            || displacement.len() != self.mode_count()
            || displacement.iter().any(|&q| !finite(q))
        {
            return Err(ShellResonanceError::Invalid(
                "finite complete mechanical pole seed required",
            ));
        }
        let mut state = vec![C64::ZERO; self.dimension];
        for (i, &q) in displacement.iter().enumerate() {
            state[2 * i] = q;
            state[2 * i + 1] = rate_per_s * q;
        }
        let mut at = 2 * self.mode_count();
        for (rows, times) in memories(&self.image) {
            for (row, &tau) in rows.iter().zip(times) {
                poll(gate)?;
                let denominator = C64::ONE + rate_per_s.scale(tau);
                if denominator.abs() <= 64.0 * f64::EPSILON * (1.0 + rate_per_s.abs() * tau) {
                    return Err(ShellResonanceError::Invalid(
                        "use a complete state seed at a material relaxation rate",
                    ));
                }
                state[at] = projection(row, displacement) / denominator;
                at += 1;
            }
        }
        if state.iter().any(|&x| !finite(x)) {
            return Err(ShellResonanceError::Invalid(
                "complete material pole seed is not representable",
            ));
        }
        poll(gate)?;
        Ok(state)
    }

    fn observation_budget(&self, count: usize) -> (u128, u128) {
        let n = self.dimension as u128;
        let p = self.image.force_ports.len() as u128
            + 2
            + self.load.as_ref().map_or(0, |load| load.ports) as u128;
        let interfaces = self.image.interfaces.len() as u128;
        (
            count as u128 * 64 * (3 * n + p + interfaces + 8),
            count as u128 * 16 * (n * n + n * (p + interfaces + 8)),
        )
    }

    /// Refine one simple mode of the selected complete continuous model. The
    /// complete state seed can select a real relaxing or zero-port-output mode.
    /// Internal Newton iterates use the explicitly authored global continuation.
    /// No finite-amplitude strain or rigorous frequency enclosure is claimed.
    pub fn resonance_pole(
        &self,
        rate_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        gate: &CancelGate,
    ) -> Result<ShellResonancePole, ShellResonanceError> {
        poll(gate)?;
        let (bytes, work) = self.observation_budget(1);
        cap(
            "pole observation bytes",
            bytes,
            options.max_dense_bytes as u128,
        )?;
        cap("pole observation work", work, options.max_work as u128)?;
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let pole = self
            .model
            .resonance_pole(rate_seed_per_s, state_seed, &numerical, || {
                gate.is_requested()
            })
            .map_err(ShellResonanceError::Pole)?;
        self.observe(pole, gate)
    }

    /// Discover simple modes inside a Laplace-plane ellipse, reusing the same
    /// complete source model and observers. Rank, count, quadrature, conditioning
    /// or resource refusal publishes no partial spectrum. The existing source
    /// fits/basis remain independent approximation axes.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        gate: &CancelGate,
    ) -> Result<ShellResonanceSpectrum, ShellResonanceError> {
        poll(gate)?;
        let (bytes, work) = self.observation_budget(options.probe_columns);
        cap(
            "spectrum observation bytes",
            bytes,
            options.max_dense_bytes as u128,
        )?;
        cap("spectrum observation work", work, options.max_work as u128)?;
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let result = self
            .model
            .resonance_spectrum(region, &numerical, || gate.is_requested())
            .map_err(ShellResonanceError::Spectrum)?;
        let mut poles = Vec::with_capacity(result.poles.len());
        for pole in result.poles {
            poles.push(self.observe(pole, gate)?);
        }
        poll(gate)?;
        Ok(ShellResonanceSpectrum {
            poles,
            unrefined_rates_per_s: result.unrefined_laplace_poles_per_s,
            argument_count: result.argument_count,
            resolved_rank: result.resolved_rank,
            quadrature_nodes: result.quadrature_nodes,
            argument_integral: result.argument_integral,
            coarse_argument_integral: result.coarse_argument_integral,
            moment_relative_change: result.moment_relative_change,
            minimum_reciprocal_condition: result.minimum_reciprocal_condition,
        })
    }

    fn observe(
        &self,
        mode: HarmonicPole,
        gate: &CancelGate,
    ) -> Result<ShellResonancePole, ShellResonanceError> {
        poll(gate)?;
        let m = self.mode_count();
        let displacement: Vec<_> = (0..m).map(|i| mode.state[2 * i]).collect();
        let velocity: Vec<_> = (0..m).map(|i| mode.state[2 * i + 1]).collect();
        let mechanical_storage_j: f64 = self
            .image
            .modes
            .iter()
            .enumerate()
            .map(|(i, body)| {
                0.5 * (velocity[i].norm_sq() + displacement[i].scale(body.omega).norm_sq())
            })
            .sum();
        let background_dissipation_w = quadratic(self.image.background, &velocity, gate)?;
        let mut interface_dissipation_w = 0.0;
        for &(row, c) in &self.image.interfaces {
            poll(gate)?;
            interface_dissipation_w += c * projection(row, &velocity).norm_sq();
        }
        let mut stores = [Vec::new(), Vec::new()];
        let mut losses = [Vec::new(), Vec::new()];
        let mut at = 2 * m;
        for (family, (rows, times)) in memories(&self.image).into_iter().enumerate() {
            stores[family].reserve(rows.len());
            losses[family].reserve(rows.len());
            for (row, &tau) in rows.iter().zip(times) {
                poll(gate)?;
                let square = (projection(row, &displacement) - mode.state[at]).norm_sq();
                stores[family].push(0.5 * square);
                losses[family].push(square / tau);
                at += 1;
            }
        }
        let shell_energy = mechanical_storage_j + stores.iter().flatten().sum::<f64>();
        let shell_heat = background_dissipation_w
            + interface_dissipation_w
            + losses.iter().flatten().sum::<f64>();
        let load = self
            .load
            .as_ref()
            .map(|load| load.observe(&mode.state[at..], &velocity, mode.laplace_pole_per_s, gate))
            .transpose()?;
        let load_energy = load.as_ref().map_or(0.0, |load| load.stored_energy_j);
        let load_heat = load.as_ref().map_or(0.0, |load| load.dissipated_power_w);
        let absorbed_power = load
            .as_ref()
            .map_or(0.0, |load| load.absorbed_port_power_w.iter().sum::<f64>());
        let interface_scale = load.as_ref().map_or(0.0, |load| {
            load.flow
                .iter()
                .zip(&load.effort)
                .map(|(u, y)| u.abs() * y.abs())
                .sum::<f64>()
        });
        let energy = shell_energy + load_energy;
        let heat = shell_heat + load_heat;
        let shell_decay_power_residual_w =
            2.0 * mode.laplace_pole_per_s.re * shell_energy + shell_heat + absorbed_power;
        let decay_power_residual_w = 2.0 * mode.laplace_pole_per_s.re * energy + heat;
        let energy_error = (energy - mode.stored_energy_j).abs()
            / (energy.abs() + mode.stored_energy_j.abs()).max(f64::MIN_POSITIVE);
        let heat_error = (heat - mode.dissipated_power_w).abs()
            / (heat.abs() + mode.dissipated_power_w.abs()).max(f64::MIN_POSITIVE);
        let decay_error = decay_power_residual_w.abs()
            / (2.0 * mode.laplace_pole_per_s.abs() * energy + heat).max(f64::MIN_POSITIVE);
        let shell_error = shell_decay_power_residual_w.abs()
            / (2.0 * mode.laplace_pole_per_s.abs() * shell_energy.abs()
                + shell_heat.abs()
                + interface_scale)
                .max(f64::MIN_POSITIVE);
        let load_error = load.as_ref().map_or(0.0, |load| {
            (load.decay_power_residual_w.abs()
                / (2.0 * mode.laplace_pole_per_s.abs() * load_energy.abs()
                    + load_heat.abs()
                    + interface_scale)
                    .max(f64::MIN_POSITIVE))
            .max(load.relative_equation_residual)
        });
        let relative = energy_error
            .max(heat_error)
            .max(decay_error)
            .max(shell_error)
            .max(load_error);
        if [
            energy,
            heat,
            decay_power_residual_w,
            shell_decay_power_residual_w,
            relative,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(ShellResonanceError::Invalid(
                "physical shell pole accounting overflow",
            ));
        }
        if relative > self.tolerance {
            return Err(ShellResonanceError::PhysicalAccounting {
                relative,
                tolerance: self.tolerance,
            });
        }
        let mut force_port_displacement_m = Vec::with_capacity(self.image.force_ports.len());
        for row in self.image.force_ports {
            poll(gate)?;
            force_port_displacement_m.push(projection(row, &displacement));
        }
        let oscillation = mode.laplace_pole_per_s.im.abs();
        let [substrate_branch_storage_j, coating_branch_storage_j] = stores;
        let [substrate_branch_dissipation_w, coating_branch_dissipation_w] = losses;
        poll(gate)?;
        Ok(ShellResonancePole {
            mode,
            displacement_m_sqrt_kg: displacement,
            velocity_m_sqrt_kg_per_s: velocity,
            force_port_displacement_m,
            mechanical_storage_j,
            background_dissipation_w,
            interface_dissipation_w,
            substrate_branch_storage_j,
            substrate_branch_dissipation_w,
            coating_branch_storage_j,
            coating_branch_dissipation_w,
            load,
            shell_decay_power_residual_w,
            decay_power_residual_w,
            physical_relative_residual: relative,
            oscillation_within_source_band: oscillation > 0.0
                && oscillation >= self.band.0
                && oscillation <= self.band.1,
        })
    }
}

impl PreparedLoad {
    fn observe(
        &self,
        state: &[C64],
        shell_velocity: &[C64],
        rate: C64,
        gate: &CancelGate,
    ) -> Result<ShellResonanceLoadResponse, ShellResonanceError> {
        let n = self.states;
        let m = shell_velocity.len();
        let mut effort_state = vec![C64::ZERO; n];
        let mut flow = vec![C64::ZERO; self.ports];
        let mut effort = vec![C64::ZERO; self.ports];
        for (i, h) in effort_state.iter_mut().enumerate() {
            poll(gate)?;
            *h = if i < self.differential_states {
                projection(&self.q[i * n..(i + 1) * n], state)
            } else {
                state[i]
            };
        }
        for p in 0..self.ports {
            poll(gate)?;
            flow[p] = projection(&self.rows[p * m..(p + 1) * m], shell_velocity);
            for (i, &h) in effort_state.iter().enumerate() {
                effort[p] = effort[p] + h.scale(self.g[i * self.ports + p]);
            }
        }
        let stored_energy_j = 0.5
            * state
                .iter()
                .zip(&effort_state)
                .take(self.differential_states)
                .map(|(z, h)| (z.conj() * *h).re)
                .sum::<f64>();
        let mut dissipated_power_w = 0.0;
        let mut relative_equation_residual = 0.0_f64;
        for i in 0..n {
            poll(gate)?;
            let loss_row = projection(&self.r[i * n..(i + 1) * n], &effort_state);
            dissipated_power_w += (effort_state[i].conj() * loss_row).re;
            let mut residual = if i < self.differential_states {
                rate * state[i]
            } else {
                C64::ZERO
            };
            let mut scale = if i < self.differential_states {
                rate.abs() * state[i].abs()
            } else {
                0.0
            };
            for (j, &h) in effort_state.iter().enumerate() {
                let a = self.structure[i * n + j];
                residual = residual - h.scale(a);
                scale += a.abs() * h.abs();
            }
            for (p, &u) in flow.iter().enumerate() {
                let b = self.g[i * self.ports + p];
                residual = residual - u.scale(b);
                scale += b.abs() * u.abs();
            }
            let relative = residual.abs() / scale.max(f64::MIN_POSITIVE);
            if !finite(residual) || !scale.is_finite() || !relative.is_finite() {
                return Err(ShellResonanceError::Invalid("dynamic load row overflow"));
            }
            relative_equation_residual = relative_equation_residual.max(relative);
        }
        let absorbed_port_power_w: Vec<_> = flow
            .iter()
            .zip(&effort)
            .map(|(u, y)| (u.conj() * *y).re)
            .collect();
        let decay_power_residual_w = 2.0 * rate.re * stored_energy_j + dissipated_power_w
            - absorbed_port_power_w.iter().sum::<f64>();
        if state
            .iter()
            .chain(&effort_state)
            .chain(&flow)
            .chain(&effort)
            .any(|&z| !finite(z))
            || absorbed_port_power_w
                .iter()
                .chain(&[stored_energy_j, dissipated_power_w, decay_power_residual_w])
                .any(|x| !x.is_finite())
        {
            return Err(ShellResonanceError::Invalid(
                "dynamic load accounting overflow",
            ));
        }
        poll(gate)?;
        Ok(ShellResonanceLoadResponse {
            state: state.to_vec(),
            flow,
            effort,
            absorbed_port_power_w,
            stored_energy_j,
            dissipated_power_w,
            decay_power_residual_w,
            relative_equation_residual,
        })
    }
}

fn memories<'a>(image: &Image<'a>) -> [(&'a [Vec<f64>], &'a [f64]); 2] {
    [
        image.substrate.map_or((&[][..], &[][..]), |m| {
            (m.coupling_rows_per_s(), m.relaxation_times_s())
        }),
        image.films.map_or((&[][..], &[][..]), |m| {
            (m.coupling_rows_per_s(), m.relaxation_times_s())
        }),
    ]
}
