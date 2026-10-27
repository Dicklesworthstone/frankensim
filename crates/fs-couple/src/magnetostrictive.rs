//! A biased magnetostrictive coil/core drives real axial motion on fs-phs.
//!
//! A fixed-root rod uses u(z)=x*z/L and the material's original incremental
//! energy. The same geometry supplies m=rho*A*L/3, flux linkage lambda=N*A*B,
//! current i=L*H/N, and force A*stress. State [x,p,lambda] advances through the
//! existing port-Hamiltonian descriptor; there is no separate actuator gain or
//! mechanics integrator. B, H, lambda, i, strain and work are perturbations of
//! the prepared DC-biased reference. Bias preparation/supply is outside this
//! incremental account; this is not the absolute energy of the biased device.
//!
//! The coil is ideal, uniformly linked, with a prescribed fixed magnetic path
//! equal to the rod length. End/fringing/demagnetizing fields, coil resistance,
//! eddy loss, hysteresis, temperature changes and acoustic radiation are absent.
//! Electrical and mechanical loading can use the original two-port component.

use fs_material::magnetostrictive::{
    LongitudinalMagnetostrictive, MagnetostrictiveError, MagnetostrictiveResponse,
};
use fs_phs::{
    DescriptorPortHamiltonian, PhsError, PortHamiltonian, StepRecord, Storage, step_descriptor,
};

/// Caller-admitted perturbation envelope; not an inferred material strength or
/// guarantee that a measured nonlinear response is linear throughout the range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MagnetostrictiveLimits {
    /// Positive absolute axial strain increment bound, strictly below one.
    pub maximum_abs_strain: f64,
    /// Must be positive and smaller than |H_bias|; cannot cross zero bias.
    pub maximum_abs_field_increment_a_m: f64,
}
/// Source, native solve, model-domain, balance or cancellation refusal.
#[derive(Debug)]
pub enum MagnetostrictiveCouplingError {
    /// Original material-law refusal.
    Material(MagnetostrictiveError),
    /// Original descriptor-system or stepper refusal.
    Native(PhsError),
    /// Invalid input or unrepresentable physical quantity.
    Invalid(&'static str),
    /// Incremental storage and actual port work disagree.
    Balance {
        /// Storage change minus mechanical/electrical supply [J].
        residual_j: f64,
        /// Caller allowance [J], not a discretization-error bound.
        tolerance_j: f64,
    },
    /// Cancellation before atomic candidate publication.
    Cancelled,
}
impl From<MagnetostrictiveError> for MagnetostrictiveCouplingError {
    fn from(value: MagnetostrictiveError) -> Self {
        Self::Material(value)
    }
}
impl From<PhsError> for MagnetostrictiveCouplingError {
    fn from(value: PhsError) -> Self {
        Self::Native(value)
    }
}
impl core::fmt::Display for MagnetostrictiveCouplingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "magnetostrictive coupling: {self:?}")
    }
}
impl std::error::Error for MagnetostrictiveCouplingError {}
type Result<T> = core::result::Result<T, MagnetostrictiveCouplingError>;
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn invalid(message: &'static str) -> MagnetostrictiveCouplingError {
    MagnetostrictiveCouplingError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MagnetostrictiveCouplingError::Cancelled)
    } else {
        Ok(())
    }
}

/// Axis-3 material in an ideal N-turn coil. The reference length/area are the
/// already biased equilibrium geometry. Positive coil current makes +3 H;
/// reversing winding orientation requires reversing both conjugate electrical
/// coordinates. Signed source d controls the actual direct/converse response.
#[derive(Clone, Debug, PartialEq)]
pub struct MagnetostrictiveRod {
    material: LongitudinalMagnetostrictive,
    length_m: f64,
    area_m2: f64,
    turns: u32,
    volume_m3: f64,
    modal_mass_kg: f64,
    limits: MagnetostrictiveLimits,
}
impl MagnetostrictiveRod {
    /// Derive consistent inertia and reciprocal coil/mechanical storage from
    /// one material and reference geometry, without a fitted transducer gain.
    pub fn new(
        material: LongitudinalMagnetostrictive,
        length_m: f64,
        area_m2: f64,
        turns: u32,
        limits: MagnetostrictiveLimits,
    ) -> Result<Self> {
        let volume_m3 = length_m * area_m2;
        let modal_mass_kg = material.density_kg_m3() * volume_m3 / 3.0;
        if [
            length_m,
            area_m2,
            volume_m3,
            modal_mass_kg,
            1.0 / modal_mass_kg,
            limits.maximum_abs_strain,
            limits.maximum_abs_field_increment_a_m,
        ]
        .iter()
        .any(|&x| !positive(x))
            || turns == 0
            || limits.maximum_abs_strain >= 1.0
            || limits.maximum_abs_field_increment_a_m >= material.bias_field_a_m().abs()
        {
            return Err(invalid(
                "positive representable geometry/turns/inertia and sub-bias small-signal limits required",
            ));
        }
        let rod = Self {
            material,
            length_m,
            area_m2,
            turns,
            volume_m3,
            modal_mass_kg,
            limits,
        };
        if [
            rod.fixed_current_stiffness_n_m(),
            rod.fixed_flux_stiffness_n_m(),
            rod.clamped_inductance_h(),
            1.0 / rod.clamped_inductance_h(),
        ]
        .iter()
        .any(|&x| !positive(x))
            || !rod.coupling_wb_m().is_finite()
            || rod.material.piezomagnetic_coefficient_m_a() != 0.0 && rod.coupling_wb_m() == 0.0
            || !rod.bias_current_a().is_finite()
        {
            return Err(invalid(
                "derived magnetostrictive storage, coupling or bias current is unrepresentable",
            ));
        }
        Ok(rod)
    }
    /// Complete prepared material, including the DC-bias source context.
    pub const fn material(&self) -> &LongitudinalMagnetostrictive {
        &self.material
    }
    /// Reference length and admitted ideal magnetic path length [m].
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Reference rod and uniformly linked flux cross section [m^2].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Positive ideal winding turn count.
    pub const fn turns(&self) -> u32 {
        self.turns
    }
    /// Consistent mass for the axial shape u(z)=x*z/L [kg].
    pub const fn modal_mass_kg(&self) -> f64 {
        self.modal_mass_kg
    }
    /// Explicit admitted perturbation envelope.
    pub const fn limits(&self) -> MagnetostrictiveLimits {
        self.limits
    }
    /// DC current for the explicitly ideal magnetic path; not a dynamic port.
    pub fn bias_current_a(&self) -> f64 {
        self.material.bias_field_a_m() * self.length_m / f64::from(self.turns)
    }
    /// Mechanical stiffness under an ideal held incremental current [N/m].
    pub fn fixed_current_stiffness_n_m(&self) -> f64 {
        self.material.modulus_at_fixed_field_pa() * self.area_m2 / self.length_m
    }
    /// Mechanical stiffness with incremental flux linkage fixed [N/m].
    pub fn fixed_flux_stiffness_n_m(&self) -> f64 {
        self.material.modulus_at_fixed_flux_pa() * self.area_m2 / self.length_m
    }
    /// Incremental inductance with axial extension held fixed [H].
    pub fn clamped_inductance_h(&self) -> f64 {
        self.material.axially_clamped_permeability_h_m()
            * f64::from(self.turns).powi(2)
            * self.area_m2
            / self.length_m
    }
    /// Reciprocal theta in lambda=L_c*i+theta*x and F=k_H*x-theta*i [Wb/m=N/A].
    pub fn coupling_wb_m(&self) -> f64 {
        f64::from(self.turns) * self.area_m2 * self.material.stress_coefficient_t() / self.length_m
    }
    fn response(&self, extension: f64, linkage: f64) -> Result<MagnetostrictiveResponse> {
        let strain = extension / self.length_m;
        let response = self
            .material
            .response(strain, linkage / (f64::from(self.turns) * self.area_m2))?;
        if strain.abs() > self.limits.maximum_abs_strain
            || response.field_increment_a_m.abs() > self.limits.maximum_abs_field_increment_a_m
        {
            return Err(invalid(
                "magnetostrictive perturbation left its strain or field envelope",
            ));
        }
        Ok(response)
    }
    fn energy(&self, state: &[f64]) -> Result<f64> {
        let response = self.response(state[0], state[2])?;
        let kinetic = 0.5 * state[1] * (state[1] / self.modal_mass_kg);
        let constitutive = self.volume_m3 * response.incremental_energy_j_m3;
        let total = kinetic + constitutive;
        if !total.is_finite()
            || state[1] != 0.0 && kinetic == 0.0
            || response.incremental_energy_j_m3 != 0.0 && constitutive == 0.0
        {
            return Err(invalid("incremental coil/core energy is unrepresentable"));
        }
        Ok(total)
    }

    /// Native two-port component: state [x,p,Delta lambda], inputs [force,
    /// Delta voltage], outputs [velocity,Delta current]. The external composer
    /// owns whole-step source-domain and work acceptance. DC bias is NOT a port.
    pub fn port_hamiltonian(&self) -> Result<PortHamiltonian> {
        Ok(PortHamiltonian::new(
            3,
            2,
            vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0; 9],
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            Box::new(RodStorage(self.clone())),
        )?)
    }
    fn descriptor(&self, current_drive: bool) -> Result<DescriptorPortHamiltonian> {
        let n = if current_drive { 4 } else { 3 };
        let mut j = vec![0.0; n * n];
        let mut g = vec![0.0; n * 2];
        j[1] = 1.0;
        j[n] = -1.0;
        g[2] = 1.0;
        if current_drive {
            // lambda_dot=V, 0=-dH/dlambda+i; V is algebraic, never storage.
            j[2 * n + 3] = 1.0;
            j[3 * n + 2] = -1.0;
            g[7] = 1.0;
        } else {
            g[5] = 1.0;
        }
        Ok(DescriptorPortHamiltonian::new(
            n,
            3,
            2,
            j,
            vec![0.0; n * n],
            g,
            Box::new(RodStorage(self.clone())),
        )?)
    }
}
#[derive(Clone)]
struct RodStorage(MagnetostrictiveRod);
impl Storage for RodStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.0.energy(x).unwrap_or(f64::NAN)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        match self.0.response(x[0], x[2]) {
            Ok(r) => {
                out[0] = self.0.area_m2 * r.stress_increment_pa;
                out[1] = x[1] / self.0.modal_mass_kg;
                out[2] = self.0.length_m * r.field_increment_a_m / f64::from(self.0.turns);
            }
            Err(_) => out.fill(f64::NAN),
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        Some(if [a, b, mid].iter().all(|x| self.0.energy(x).is_ok()) {
            0.0
        } else {
            f64::NAN
        })
    }
}

/// Held incremental electrical boundary. Current-zero holds H at its DC bias;
/// voltage-zero conserves incremental flux. These are different mechanical limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoilDrive {
    /// Voltage perturbation across the ideal coil [V]; zero holds its flux.
    VoltageIncrement(f64),
    /// Current perturbation [A]. Initial flux must satisfy this ideal clamp.
    CurrentIncrement(f64),
}
/// Independent physical-work and ideal-current-constraint allowances.
#[derive(Clone, Copy, Debug)]
pub struct MagnetostrictiveSolve {
    /// Absolute incremental storage/work residual allowance [J].
    pub energy_tolerance_j: f64,
    /// Absolute initial/final current-constraint allowance [A].
    pub current_tolerance_a: f64,
}
/// Accepted mechanical and incremental magnetic state with immutable sources.
#[derive(Clone, Debug, PartialEq)]
pub struct MagnetostrictiveBody {
    rod: MagnetostrictiveRod,
    state: [f64; 3],
}
/// Complete joint proposal; inspecting it does not advance the accepted body.
pub struct MagnetostrictiveCandidate {
    body: MagnetostrictiveBody,
    step: MagnetostrictiveStep,
}
impl MagnetostrictiveCandidate {
    /// Proposed mechanical and magnetic endpoint.
    pub const fn body(&self) -> &MagnetostrictiveBody {
        &self.body
    }
    /// Native solver result and independently reconstructed work.
    pub const fn step(&self) -> &MagnetostrictiveStep {
        &self.step
    }
    /// Publish endpoint and work account together.
    pub fn into_parts(self) -> (MagnetostrictiveBody, MagnetostrictiveStep) {
        (self.body, self.step)
    }
}
/// Incremental reversible work; no physical heat or absolute bias energy claim.
pub struct MagnetostrictiveStep {
    /// Original descriptor record; current drive appends algebraic voltage.
    pub native: StepRecord,
    /// Held axial force times actual extension change [J].
    pub mechanical_supply_j: f64,
    /// Incremental voltage/current port work [J], excluding DC-bias preparation.
    pub incremental_electrical_supply_j: f64,
    /// Incremental storage change minus both source works [J].
    pub incremental_energy_residual_j: f64,
}
impl MagnetostrictiveBody {
    /// Admit perturbations about the prepared reference without projecting flux
    /// or discarding pre-existing mechanical/magnetic incremental energy.
    pub fn new(
        rod: MagnetostrictiveRod,
        extension_m: f64,
        velocity_m_s: f64,
        flux_linkage_increment_wb: f64,
    ) -> Result<Self> {
        let momentum = rod.modal_mass_kg * velocity_m_s;
        let state = [extension_m, momentum, flux_linkage_increment_wb];
        if !velocity_m_s.is_finite() || velocity_m_s != 0.0 && momentum == 0.0 {
            return Err(invalid("initial velocity or momentum is unrepresentable"));
        }
        rod.energy(&state)?;
        Ok(Self { rod, state })
    }
    /// Immutable original geometry and material binding.
    pub const fn rod(&self) -> &MagnetostrictiveRod {
        &self.rod
    }
    /// [Extension m, momentum kg*m/s, flux-linkage increment Wb].
    pub const fn state(&self) -> [f64; 3] {
        self.state
    }
    /// Tip extension from the prepared biased shape [m].
    pub const fn extension_m(&self) -> f64 {
        self.state[0]
    }
    /// Physical axial tip velocity [m/s].
    pub fn velocity_m_s(&self) -> f64 {
        self.state[1] / self.rod.modal_mass_kg
    }
    /// N times the incremental linked flux [Wb].
    pub const fn flux_linkage_increment_wb(&self) -> f64 {
        self.state[2]
    }
    /// Accepted stress, H and material energy increments.
    pub fn response(&self) -> Result<MagnetostrictiveResponse> {
        self.rod.response(self.state[0], self.state[2])
    }
    /// Coil-current increment about the prescribed DC bias [A].
    pub fn current_increment_a(&self) -> Result<f64> {
        Ok(self.rod.length_m * self.response()?.field_increment_a_m / f64::from(self.rod.turns))
    }
    /// Reversible perturbation energy plus actual kinetic energy [J].
    pub fn incremental_energy_j(&self) -> Result<f64> {
        self.rod.energy(&self.state)
    }

    /// Advance the coupled state with the original descriptor solver. Sources,
    /// geometry and bias remain fixed. Refusal/cancellation never modifies this
    /// accepted state. Work admission is not a time-discretization accuracy claim.
    pub fn step_candidate(
        &self,
        force_n: f64,
        drive: CoilDrive,
        duration_s: f64,
        controls: MagnetostrictiveSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MagnetostrictiveCandidate> {
        poll(&mut cancelled)?;
        let (value, current_drive) = match drive {
            CoilDrive::VoltageIncrement(v) => (v, false),
            CoilDrive::CurrentIncrement(i) => (i, true),
        };
        if !force_n.is_finite()
            || !value.is_finite()
            || !positive(duration_s)
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.current_tolerance_a)
        {
            return Err(invalid(
                "finite held drives and positive duration/work/current tolerances required",
            ));
        }
        if current_drive
            && (self.current_increment_a()? - value).abs() > controls.current_tolerance_a
        {
            return Err(invalid(
                "ideal current clamp conflicts with initial flux; switching needs an explicit circuit",
            ));
        }
        let system = self.rod.descriptor(current_drive)?;
        let mut x0 = self.state.to_vec();
        if current_drive {
            x0.push(self.rod.coupling_wb_m() * self.velocity_m_s());
        }
        let inputs = [force_n, value];
        if current_drive {
            let initial = step_descriptor(&system, &x0, &inputs, 0.0)?;
            if initial.x[..3] != self.state {
                return Err(invalid(
                    "current-clamp consistency changed conserved coordinates",
                ));
            }
        }
        poll(&mut cancelled)?;
        let native = step_descriptor(&system, &x0, &inputs, duration_s)?;
        let state = [native.x[0], native.x[1], native.x[2]];
        let body = Self {
            rod: self.rod.clone(),
            state,
        };
        let before = self.incremental_energy_j()?;
        let after = body.incremental_energy_j()?;
        let mechanical_supply_j = force_n * (state[0] - self.state[0]);
        let incremental_electrical_supply_j = if current_drive {
            value * (state[2] - self.state[2])
        } else {
            duration_s * value * 0.5 * (self.current_increment_a()? + body.current_increment_a()?)
        };
        let incremental_energy_residual_j =
            after - before - mechanical_supply_j - incremental_electrical_supply_j;
        if current_drive
            && (body.current_increment_a()? - value).abs() > controls.current_tolerance_a
        {
            return Err(invalid("endpoint violated held current constraint"));
        }
        if [
            before,
            after,
            mechanical_supply_j,
            incremental_electrical_supply_j,
            incremental_energy_residual_j,
            native.supplied,
            native.delta_h,
            native.dissipated,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "incremental magnetostrictive work is unrepresentable",
            ));
        }
        let tolerance = controls.energy_tolerance_j;
        if incremental_energy_residual_j.abs() > tolerance
            || native.balance_residual().abs() > tolerance
            || (after - before - native.delta_h).abs() > tolerance
            || (mechanical_supply_j + incremental_electrical_supply_j - native.supplied).abs()
                > tolerance
            || native.dissipated != 0.0
        {
            return Err(MagnetostrictiveCouplingError::Balance {
                residual_j: incremental_energy_residual_j,
                tolerance_j: tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(MagnetostrictiveCandidate {
            body,
            step: MagnetostrictiveStep {
                native,
                mechanical_supply_j,
                incremental_electrical_supply_j,
                incremental_energy_residual_j,
            },
        })
    }
}
