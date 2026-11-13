//! Causal normal-incidence reflection from original oriented wall mechanics.
//!
//! Inputs and outputs are REAL instantaneous modal velocities, not peak phasors.
//! For Cartesian incident characteristics `a`, wall velocity `w` and outgoing
//! characteristics `b=w-a`, instantaneous powers are `a^T D a` and `b^T D b`,
//! with the original positive characteristic impedance D. Native Christoffel
//! modes diagonalize D for anisotropic solids. There is no averaging factor 1/2.
//!
//! Six original principal displacement/velocity coordinates evolve through the
//! native exact held-input quadratic flow. Its cached trajectory integrals keep
//! original pin heat separate from outgoing-wave energy. Numerical energy
//! defects are never assigned to resistance. Warm steps use fixed-size scratch.
//!
//! This is a stationary, locally reacting, linear support on a lossless elastic
//! half-space at normal wavevector incidence. It retains causal support memory, not an
//! oblique convolution, resolved coating, moving surface or spatial wave history.
//! Original strain admission bounds this interval's incoming/emitted waves;
//! previously emitted fields elsewhere in the half-space remain caller-owned.

use super::OrientedOscillatorSupport;
use fs_material::MaterialError;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::normal::NormalElasticCharacteristic;
use fs_material::elastic::acoustic::interface::surface::{
    ElasticSurfaceCondition, IsotropicElasticSurface,
};
use fs_material::elastic::acoustic::interface::{ElasticMode, ElasticWaveTrace};
use fs_math::{dd::Dd, det};
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, prepared::PreparedQuadraticFlow, split_quadratic,
};
use fs_phs::{
    PhsError, PortHamiltonian, PreparedStepError, QuadraticStorage, wall_specific_impedance,
};

/// Native dense preparation budget and a relative physical work tolerance.
#[derive(Clone, Copy, Debug)]
pub struct TransientElasticReflectionLimits {
    pub preparation: DenseExponentialBudget,
    pub energy_balance_tolerance: f64,
}

impl Default for TransientElasticReflectionLimits {
    fn default() -> Self {
        Self {
            preparation: DenseExponentialBudget {
                // PreparedQuadraticFlow budgets n+m+1=10. The original split
                // separately admits its doubled 20-coordinate block <=256.
                max_dimension: 10,
                max_dense_bytes: 4 * 1024 * 1024,
                max_multiply_terms: 10_000_000,
            },
            energy_balance_tolerance: 5e-10,
        }
    }
}

/// Original material/dynamics refusals or an unaccepted interval.
#[derive(Debug)]
pub enum TransientElasticReflectionError {
    Material(MaterialError),
    Dynamics(PhsError),
    Invalid(&'static str),
    Cancelled,
    EnergyBalance {
        residual_j_m2: f64,
        tolerance_j_m2: f64,
    },
}

impl From<MaterialError> for TransientElasticReflectionError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<PhsError> for TransientElasticReflectionError {
    fn from(value: PhsError) -> Self {
        Self::Dynamics(value)
    }
}
impl core::fmt::Display for TransientElasticReflectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "transient elastic material: {error}"),
            Self::Dynamics(error) => write!(f, "transient elastic mechanics: {error}"),
            Self::Invalid(what) => write!(f, "transient elastic reflection: {what}"),
            Self::Cancelled => f.write_str("transient elastic reflection cancelled"),
            Self::EnergyBalance {
                residual_j_m2,
                tolerance_j_m2,
            } => write!(
                f,
                "transient reflection work defect {residual_j_m2} J/m² exceeds {tolerance_j_m2} J/m²"
            ),
        }
    }
}
impl core::error::Error for TransientElasticReflectionError {}

type Result<T> = core::result::Result<T, TransientElasticReflectionError>;

/// Complete accepted original support motion and one physical clock.
/// Compatible preparations with different durations can continue this state.
#[derive(Clone, Debug)]
pub struct TransientElasticSupportState<'a> {
    material: &'a ElasticAcousticState,
    support: OrientedOscillatorSupport,
    surface_to_material: [[f64; 3]; 3],
    time_s: f64,
    displacement_m: [f64; 3],
    velocity_m_s: [f64; 3],
}

impl TransientElasticSupportState<'_> {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn principal_displacement_m(&self) -> [f64; 3] {
        self.displacement_m
    }
    pub const fn principal_velocity_m_s(&self) -> [f64; 3] {
        self.velocity_m_s
    }
}

/// Instantaneous original pin storage; no cycle averaging is applied.
#[derive(Clone, Copy, Debug, Default)]
pub struct TransientElasticSupportEnergy {
    pub kinetic_j_m2: [f64; 3],
    pub elastic_j_m2: [f64; 3],
    pub total_j_m2: f64,
}

/// Independent accounts on the exact same held-input trajectory [J/m²].
/// Small signed contraction roundoff is retained within the reported tolerance.
#[derive(Clone, Copy, Debug, Default)]
pub struct TransientElasticReflectionWork {
    pub initial: TransientElasticSupportEnergy,
    pub final_energy: TransientElasticSupportEnergy,
    /// Evaluated bilinearly, avoiding subtraction of large total energies.
    pub stored_energy_change_j_m2: f64,
    pub incoming_wave_energy_j_m2: f64,
    /// Integral of `(w-a)^T D (w-a)`, independently expanded from native
    /// quadratic and linear trajectory observations, not an energy deficit.
    pub outgoing_wave_energy_j_m2: f64,
    /// Integral of the original principal `sum r_i v_i²` only.
    pub support_dissipation_j_m2: f64,
    pub balance_tolerance_j_m2: f64,
}

impl TransientElasticReflectionWork {
    pub fn balance_residual_j_m2(self) -> f64 {
        sum([
            self.stored_energy_change_j_m2,
            self.outgoing_wave_energy_j_m2,
            self.support_dissipation_j_m2,
            -self.incoming_wave_energy_j_m2,
        ])
    }
}

/// Endpoint and interval-mean outgoing characteristics from one accepted step.
#[derive(Clone, Copy, Debug, Default)]
pub struct TransientElasticReflectionFrame {
    pub time_s: f64,
    pub duration_s: f64,
    /// Held REAL instantaneous modal velocities [m/s], in the prepared owner's
    /// channel order: legacy P/SV/SH or native ascending-speed quasi-modes.
    pub incoming_velocity_m_s: [f64; 3],
    pub initial_reflected_velocity_m_s: [f64; 3],
    pub final_reflected_velocity_m_s: [f64; 3],
    pub mean_reflected_velocity_m_s: [f64; 3],
    /// Actual final Cartesian wall velocity, x/y tangential and z normal.
    pub support_velocity_m_s: [f64; 3],
    /// Final traction on the original solid's +z surface [Pa].
    pub support_traction_pa: [f64; 3],
    /// Conservative original Mandel-strain bound throughout this interval.
    pub interval_strain_norm_bound: f64,
    pub work: TransientElasticReflectionWork,
}

/// Fixed-size reusable native contraction buffers. Contents are scratch until
/// a successful step publishes its complete state and frame.
#[derive(Default)]
pub struct TransientElasticReflectionWorkspace {
    original: [f64; 6],
    augmented: [f64; 10],
    balanced: [f64; 6],
    endpoint: [f64; 6],
    integrals: [f64; 5],
    frame: TransientElasticReflectionFrame,
}

/// One original material/support pair and a cached held-input duration.
pub struct PreparedTransientElasticReflection<'a> {
    material: &'a ElasticAcousticState,
    support: OrientedOscillatorSupport,
    rotation: [f64; 9],
    surface_to_material: [[f64; 3]; 3],
    incoming_polarizations: [[f64; 3]; 3],
    reflected_polarizations: [[f64; 3]; 3],
    principal_to_modes: [[f64; 3]; 3],
    impedance: [f64; 3],
    root_mass: [f64; 3],
    root_stiffness: [f64; 3],
    incoming_strain: [f64; 3],
    reflected_strain: [f64; 3],
    duration_s: f64,
    tolerance: f64,
    flow: PreparedQuadraticFlow,
}

impl<'a> PreparedTransientElasticReflection<'a> {
    pub fn prepare(
        material: &'a ElasticAcousticState,
        support: OrientedOscillatorSupport,
        duration_s: f64,
    ) -> Result<Self> {
        Self::prepare_controlled(
            material,
            support,
            duration_s,
            TransientElasticReflectionLimits::default(),
            &mut || false,
        )
    }

    /// Native exact preparation retains all six original support coordinates.
    /// The finite dense cache is admitted before its native matrix allocations.
    pub fn prepare_controlled(
        material: &'a ElasticAcousticState,
        support: OrientedOscillatorSupport,
        duration_s: f64,
        limits: TransientElasticReflectionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let surface =
            IsotropicElasticSurface::new(material, ElasticSurfaceCondition::TractionFree)?
                .prepare_incidence_cancellable(ElasticMode::P, 1.0, cancelled)?;
        let modes = [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh];
        let incoming = modes.map(|mode| {
            surface
                .incoming_mode(mode)
                .expect("normal-incidence modes are propagating")
        });
        let reflected = modes.map(|mode| surface.reflected_mode(mode));
        Self::prepare_modes(
            material,
            support,
            duration_s,
            limits,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            incoming,
            reflected,
            cancelled,
        )
    }

    /// Retain the original anisotropic Christoffel modes and physical frame.
    /// Input/output channels follow the characteristic's native ascending-speed
    /// order and SAME-polarization two-way convention, not legacy P/S labels.
    /// Normal wavevectors may carry energy laterally in an anisotropic solid.
    pub fn prepare_normal_characteristic(
        characteristic: &NormalElasticCharacteristic<'a>,
        support: OrientedOscillatorSupport,
        duration_s: f64,
    ) -> Result<Self> {
        Self::prepare_normal_characteristic_controlled(
            characteristic,
            support,
            duration_s,
            TransientElasticReflectionLimits::default(),
            &mut || false,
        )
    }

    pub fn prepare_normal_characteristic_controlled(
        characteristic: &NormalElasticCharacteristic<'a>,
        support: OrientedOscillatorSupport,
        duration_s: f64,
        limits: TransientElasticReflectionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let incoming = core::array::from_fn(|mode| {
            characteristic
                .incoming_mode(mode)
                .expect("three native normal modes")
        });
        let reflected = core::array::from_fn(|mode| {
            characteristic
                .reflected_mode(mode)
                .expect("three native normal modes")
        });
        Self::prepare_modes(
            characteristic.state(),
            support,
            duration_s,
            limits,
            characteristic.surface_to_material(),
            incoming,
            reflected,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_modes(
        material: &'a ElasticAcousticState,
        support: OrientedOscillatorSupport,
        duration_s: f64,
        limits: TransientElasticReflectionLimits,
        surface_to_material: [[f64; 3]; 3],
        incoming: [ElasticWaveTrace; 3],
        reflected: [ElasticWaveTrace; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !limits.energy_balance_tolerance.is_finite()
            || limits.energy_balance_tolerance <= 0.0
            || limits.energy_balance_tolerance >= 1.0
        {
            return Err(invalid(
                "positive finite duration and work tolerance in (0,1) required",
            ));
        }
        let pins = *support.principal_supports();
        for pin in &pins {
            // The native wall owner admits the original positive mass/stiffness
            // and nonnegative resistance. No frequency response is fitted here.
            wall_specific_impedance(pin, 1.0)?;
        }
        let impedance = incoming.map(|wave| 2.0 * wave.normal_power_factor_pa_s_m);
        if impedance.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(invalid(
                "normal characteristic impedances are unrepresentable",
            ));
        }
        let rotation = support.principal_to_surface().matrix().m;
        let incoming_polarizations = incoming.map(|wave| wave.velocity.map(|v| v.re));
        let reflected_polarizations = reflected.map(|wave| wave.velocity.map(|v| v.re));
        let principal_to_modes: [[f64; 3]; 3] = core::array::from_fn(|mode| {
            core::array::from_fn(|axis| {
                sum((0..3)
                    .map(|cart| incoming_polarizations[mode][cart] * rotation[3 * cart + axis]))
            })
        });
        let root_mass = pins.map(|pin| det::sqrt(pin.surface_density));
        let root_stiffness = pins.map(|pin| det::sqrt(pin.stiffness_per_area));
        // Energy-scaled ORIGINAL coordinates [sqrt(K)x, sqrt(M)v]. Q=I;
        // the coordinate change adds no physical storage or artificial pole.
        let mut q = vec![0.0; 36];
        let mut j = vec![0.0; 36];
        let mut r = vec![0.0; 36];
        let mut g = vec![0.0; 18];
        let mut heat = vec![0.0; 36];
        let mut radiation = vec![0.0; 36];
        let mut linear = vec![vec![0.0; 6]; 3];
        for axis in 0..3 {
            poll(cancelled)?;
            q[6 * axis + axis] = 1.0;
            q[6 * (3 + axis) + 3 + axis] = 1.0;
            let canonical = root_stiffness[axis] / root_mass[axis];
            if !canonical.is_finite() || canonical == 0.0 {
                return Err(invalid("original support rate is unrepresentable"));
            }
            j[6 * axis + 3 + axis] = canonical;
            j[6 * (3 + axis) + axis] = -canonical;
            linear[axis][3 + axis] = 1.0 / root_mass[axis];
            let loss = pins[axis].resistance / root_mass[axis] / root_mass[axis];
            if pins[axis].resistance > 0.0 && loss == 0.0 {
                return Err(invalid(
                    "original support resistance is below representable range",
                ));
            }
            heat[6 * (3 + axis) + 3 + axis] = loss;
            for other in axis..3 {
                let value = sum((0..3).map(|mode| {
                    principal_to_modes[mode][axis]
                        * impedance[mode]
                        * principal_to_modes[mode][other]
                })) / root_mass[axis]
                    / root_mass[other];
                if axis == other && value <= 0.0 {
                    return Err(invalid("original radiation resistance is unrepresentable"));
                }
                radiation[6 * (3 + axis) + 3 + other] = value;
                radiation[6 * (3 + other) + 3 + axis] = value;
                let total = value + if axis == other { loss } else { 0.0 };
                r[6 * (3 + axis) + 3 + other] = total;
                r[6 * (3 + other) + 3 + axis] = total;
            }
            for mode in 0..3 {
                g[3 * (3 + axis) + mode] =
                    (impedance[mode] / root_mass[axis]) * (2.0 * principal_to_modes[mode][axis]);
            }
        }
        if q.iter()
            .chain(&j)
            .chain(&r)
            .chain(&g)
            .chain(&heat)
            .chain(&radiation)
            .any(|v| !v.is_finite())
        {
            return Err(invalid("original support operator is unrepresentable"));
        }
        let system = PortHamiltonian::new(6, 3, j, r, g, Box::new(QuadraticStorage::new(q, 6)?))?;
        poll(cancelled)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: 6,
                max_omega_dt: 1.0,
            },
        )?;
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            duration_s,
            &[heat, radiation],
            &linear,
            limits.preparation,
            &mut *cancelled,
        )
        .map_err(|error| match error {
            PreparedStepError::Solver(error) => TransientElasticReflectionError::Dynamics(error),
            PreparedStepError::Cancelled => TransientElasticReflectionError::Cancelled,
        })?;
        poll(cancelled)?;
        Ok(Self {
            material,
            support,
            rotation,
            surface_to_material,
            incoming_polarizations,
            reflected_polarizations,
            principal_to_modes,
            impedance,
            root_mass,
            root_stiffness,
            incoming_strain: incoming.map(strain_norm),
            reflected_strain: reflected.map(strain_norm),
            duration_s,
            tolerance: limits.energy_balance_tolerance,
            flow,
        })
    }

    pub const fn material(&self) -> &'a ElasticAcousticState {
        self.material
    }
    pub const fn support(&self) -> &OrientedOscillatorSupport {
        &self.support
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Columns define the physical surface frame in the original material axes.
    pub const fn surface_to_material(&self) -> [[f64; 3]; 3] {
        self.surface_to_material
    }
    /// Outer index is the input channel; inner index is surface Cartesian x/y/z.
    pub const fn incoming_polarizations_surface(&self) -> [[f64; 3]; 3] {
        self.incoming_polarizations
    }
    pub const fn reflected_polarizations_surface(&self) -> [[f64; 3]; 3] {
        self.reflected_polarizations
    }
    pub const fn mode_impedances_pa_s_m(&self) -> [f64; 3] {
        self.impedance
    }
    pub fn workspace(&self) -> TransientElasticReflectionWorkspace {
        Default::default()
    }

    /// Initialize the six physical principal coordinates without erasing memory.
    /// The original stored-energy envelope must admit freely emitted waves.
    pub fn state(
        &self,
        time_s: f64,
        principal_displacement_m: [f64; 3],
        principal_velocity_m_s: [f64; 3],
    ) -> Result<TransientElasticSupportState<'a>> {
        if !time_s.is_finite()
            || principal_displacement_m
                .iter()
                .chain(&principal_velocity_m_s)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite original support coordinates and clock required",
            ));
        }
        let state = TransientElasticSupportState {
            material: self.material,
            support: self.support,
            surface_to_material: self.surface_to_material,
            time_s,
            displacement_m: principal_displacement_m,
            velocity_m_s: principal_velocity_m_s,
        };
        let energy = self.energy(&state)?;
        self.admit_strain([0.0; 3], energy.total_j_m2)?;
        Ok(state)
    }

    pub fn energy(
        &self,
        state: &TransientElasticSupportState<'_>,
    ) -> Result<TransientElasticSupportEnergy> {
        self.validate(state)?;
        let mut result = TransientElasticSupportEnergy::default();
        for axis in 0..3 {
            result.kinetic_j_m2[axis] =
                energy_square(self.root_mass[axis], state.velocity_m_s[axis])?;
            result.elastic_j_m2[axis] =
                energy_square(self.root_stiffness[axis], state.displacement_m[axis])?;
        }
        result.total_j_m2 = sum(result.kinetic_j_m2.into_iter().chain(result.elastic_j_m2));
        if !result.total_j_m2.is_finite() {
            return Err(invalid("original support storage overflows"));
        }
        Ok(result)
    }

    /// Advance one held real incident triple transactionally. Preparation owns
    /// every exponential and Gramian; this path allocates and factorizes nothing.
    /// Failure leaves the accepted physical state and clock unchanged.
    pub fn step_into<'w>(
        &self,
        state: &mut TransientElasticSupportState<'a>,
        incoming_velocity_m_s: [f64; 3],
        workspace: &'w mut TransientElasticReflectionWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<&'w TransientElasticReflectionFrame> {
        poll(cancelled)?;
        self.validate(state)?;
        if incoming_velocity_m_s.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite real held incident modal velocities required",
            ));
        }
        let time_s = state.time_s + self.duration_s;
        if !time_s.is_finite() || time_s <= state.time_s {
            return Err(invalid(
                "held interval must advance the finite physical clock",
            ));
        }
        let a = incoming_velocity_m_s;
        let initial = self.energy(state)?;
        let incoming_power = sum((0..3).map(|axis| {
            let root = det::sqrt(self.impedance[axis]) * a[axis];
            root * root
        }));
        let incoming_work = incoming_power * self.duration_s;
        if !incoming_work.is_finite()
            || (incoming_power > 0.0 && incoming_work == 0.0)
            || (incoming_power == 0.0 && a.into_iter().any(|value| value != 0.0))
        {
            return Err(invalid("incident interval energy is unrepresentable"));
        }
        // Passivity gives H(t)<=H0+integral Pin for EVERY point of the interval.
        // The resulting emitted-mode envelope is checked before advancing.
        let energy_bound = sum([initial.total_j_m2, incoming_work]);
        let interval_strain_norm_bound = self.admit_strain(incoming_velocity_m_s, energy_bound)?;
        for axis in 0..3 {
            workspace.original[axis] = self.root_stiffness[axis] * state.displacement_m[axis];
            workspace.original[3 + axis] = self.root_mass[axis] * state.velocity_m_s[axis];
        }
        self.flow.evaluate_into(
            &workspace.original,
            &a,
            &mut workspace.augmented,
            &mut workspace.balanced,
            &mut workspace.endpoint,
            &mut workspace.integrals,
        )?;
        poll(cancelled)?;
        let candidate = TransientElasticSupportState {
            material: self.material,
            support: self.support,
            surface_to_material: self.surface_to_material,
            time_s,
            displacement_m: core::array::from_fn(|axis| {
                workspace.endpoint[axis] / self.root_stiffness[axis]
            }),
            velocity_m_s: core::array::from_fn(|axis| {
                workspace.endpoint[3 + axis] / self.root_mass[axis]
            }),
        };
        let final_energy = self.energy(&candidate)?;
        let change = sum((0..6).map(|i| {
            let first = Dd::from_f64(workspace.original[i]);
            let last = Dd::from_f64(workspace.endpoint[i]);
            (Dd::from_f64(0.5) * (last - first) * (last + first)).to_f64()
        }));
        let integrated_principal = core::array::from_fn(|i| workspace.integrals[2 + i]);
        let integrated_wall = self.rotate(integrated_principal);
        let integrated_modal = self.principal_modes(integrated_principal);
        let supplied = sum((0..3).map(|i| 2.0 * self.impedance[i] * a[i] * integrated_modal[i]));
        let outgoing = sum([workspace.integrals[1], -supplied, incoming_work]);
        let heat = workspace.integrals[0];
        let scale = initial
            .total_j_m2
            .max(final_energy.total_j_m2)
            .max(change.abs())
            .max(incoming_work)
            .max(supplied.abs())
            .max(workspace.integrals[1].abs())
            .max(heat.abs())
            .max(outgoing.abs());
        let tolerance = (self.tolerance + 512.0 * f64::EPSILON) * scale;
        let work = TransientElasticReflectionWork {
            initial,
            final_energy,
            stored_energy_change_j_m2: change,
            incoming_wave_energy_j_m2: incoming_work,
            outgoing_wave_energy_j_m2: outgoing,
            support_dissipation_j_m2: heat,
            balance_tolerance_j_m2: tolerance,
        };
        let residual = work.balance_residual_j_m2();
        if [change, incoming_work, outgoing, heat, supplied, tolerance]
            .into_iter()
            .any(|v| !v.is_finite())
            || residual.abs() > tolerance
            || outgoing < -tolerance
            || heat < -tolerance
        {
            return Err(TransientElasticReflectionError::EnergyBalance {
                residual_j_m2: residual,
                tolerance_j_m2: tolerance,
            });
        }
        let initial_wall = self.rotate(state.velocity_m_s);
        let final_wall = self.rotate(candidate.velocity_m_s);
        let mean_wall = integrated_wall.map(|value| value / self.duration_s);
        let final_modal = self.principal_modes(candidate.velocity_m_s);
        let traction = self.modal_surface(core::array::from_fn(|mode| {
            self.impedance[mode] * (final_modal[mode] - 2.0 * a[mode])
        }));
        let incoming_surface = self.modal_surface(a);
        if initial_wall
            .into_iter()
            .chain(final_wall)
            .chain(mean_wall)
            .chain(traction)
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "emitted characteristics or support traction overflow",
            ));
        }
        let frame = TransientElasticReflectionFrame {
            time_s,
            duration_s: self.duration_s,
            incoming_velocity_m_s,
            initial_reflected_velocity_m_s: self.reflected_modes(initial_wall, incoming_surface),
            final_reflected_velocity_m_s: self.reflected_modes(final_wall, incoming_surface),
            mean_reflected_velocity_m_s: self.reflected_modes(mean_wall, incoming_surface),
            support_velocity_m_s: final_wall,
            support_traction_pa: traction,
            interval_strain_norm_bound,
            work,
        };
        poll(cancelled)?;
        *state = candidate;
        workspace.frame = frame;
        Ok(&workspace.frame)
    }

    fn validate(&self, state: &TransientElasticSupportState<'_>) -> Result<()> {
        if !core::ptr::eq(self.material, state.material)
            || self.support != state.support
            || self.surface_to_material != state.surface_to_material
        {
            return Err(invalid(
                "support state belongs to another original material, support or surface frame",
            ));
        }
        Ok(())
    }

    fn rotate(&self, principal: [f64; 3]) -> [f64; 3] {
        core::array::from_fn(|row| {
            sum((0..3).map(|axis| self.rotation[3 * row + axis] * principal[axis]))
        })
    }

    fn principal_modes(&self, principal: [f64; 3]) -> [f64; 3] {
        core::array::from_fn(|mode| {
            sum((0..3).map(|axis| self.principal_to_modes[mode][axis] * principal[axis]))
        })
    }

    fn modal_surface(&self, modal: [f64; 3]) -> [f64; 3] {
        core::array::from_fn(|cart| {
            sum((0..3).map(|mode| self.incoming_polarizations[mode][cart] * modal[mode]))
        })
    }

    fn reflected_modes(&self, wall: [f64; 3], incoming: [f64; 3]) -> [f64; 3] {
        core::array::from_fn(|mode| {
            sum((0..3).map(|cart| {
                self.reflected_polarizations[mode][cart] * (wall[cart] - incoming[cart])
            }))
        })
    }

    fn admit_strain(&self, incoming: [f64; 3], energy_bound: f64) -> Result<f64> {
        if !energy_bound.is_finite() || energy_bound < 0.0 {
            return Err(invalid(
                "complete support energy envelope is unrepresentable",
            ));
        }
        let energy_root = core::f64::consts::SQRT_2 * det::sqrt(energy_bound);
        let wall_bounds: [f64; 3] = core::array::from_fn(|mode| {
            let inverse_mass_norm = (0..3).fold(0.0, |norm, axis| {
                let projected = sum((0..3).map(|cart| {
                    self.reflected_polarizations[mode][cart] * self.rotation[3 * cart + axis]
                }));
                det::hypot(norm, projected / self.root_mass[axis])
            });
            energy_root * inverse_mass_norm
        });
        let incoming_surface = self.modal_surface(incoming);
        let bound = sum((0..3).map(|mode| {
            let incoming_projection = sum((0..3)
                .map(|cart| self.reflected_polarizations[mode][cart] * incoming_surface[cart]));
            self.incoming_strain[mode] * incoming[mode].abs()
                + self.reflected_strain[mode] * (wall_bounds[mode] + incoming_projection.abs())
        }));
        if !bound.is_finite() || bound > self.material.maximum_absolute_mandel_strain() {
            return Err(invalid(
                "incoming/emitted interval exceeds the original material strain domain",
            ));
        }
        Ok(bound)
    }
}

fn strain_norm(wave: ElasticWaveTrace) -> f64 {
    wave.strain_mandel_s_m
        .iter()
        .fold(0.0, |norm, value| det::hypot(norm, value.abs()))
}
fn energy_square(root: f64, amplitude: f64) -> Result<f64> {
    let value = root * (core::f64::consts::FRAC_1_SQRT_2 * amplitude);
    let energy = value * value;
    if !energy.is_finite() || (amplitude != 0.0 && energy == 0.0) {
        Err(invalid("original support energy is unrepresentable"))
    } else {
        Ok(energy)
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |value, term| value + Dd::from_f64(term))
        .to_f64()
}
fn invalid(what: &'static str) -> TransientElasticReflectionError {
    TransientElasticReflectionError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(TransientElasticReflectionError::Cancelled)
    } else {
        Ok(())
    }
}
