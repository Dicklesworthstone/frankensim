//! Reciprocal Coulomb momentum and heat transfer between Maxwellian species.
//!
//! This is a classical, unmagnetized, weak-coupling, weak-drift Landau moment
//! closure for the actual frozen populations. It is not independent drag
//! against stationary ions, Spitzer-Harm electrical conductivity, or a kinetic
//! solution for a distorted electron distribution. All unordered charged pairs
//! are retained, including different ionic species; self-collisions transfer
//! neither a species' bulk momentum nor its thermal energy to another pool.
//!
//! For relative one-component thermal variance `S=kB(Ta/ma+Tb/mb)`, reduced
//! mass `mu=ma*mb/(ma+mb)` and `C=|qa*qb|/(4*pi*eps0)`, the pair coefficients are
//! `K=(4*sqrt(2*pi)/3)*na*nb*C²*lnLambda/(mu*S^(3/2))` [kg/(m³ s)] and
//! `G=3*kB*K/(ma+mb)` [W/(m³ K)]. Thus `Fa=-K(ua-ub)=-Fb` and the zero-drift
//! thermal powers are `Qa=-G(Ta-Tb)=-Qb`. These are the Maxwellian Landau
//! moments, with finite masses retained: Fitzpatrick, Eqs.3.179 and3.183,
//! <https://farside.ph.utexas.edu/teaching/plasma/lectures/node42.html>;
//! NRL Plasma Formulary (2023), thermal equilibration, pp.33–34,
//! <https://www.nrl.navy.mil/News-Media/Publications/NRL-PlasmaFormulary/>.
//!
//! Coulomb logarithms are approximate leading-logarithmic prescriptions, not
//! exact screening calculations. The caller selects all-carrier STATIC Debye
//! screening or electron-only Debye screening. Both use the explicitly fixed
//! thermal convention `b90=C/(mu*S)` and `lnLambda=ln(lambda_D/b90)`; no fitted
//! constant, logarithm floor or quantum interpolation is inserted. The quantum
//! length `hbar/(2*mu*sqrt(S))` must be smaller than a declared fraction of b90.
//! This thermal convention and the static screening choices must not be
//! confused with NRL's velocity-dependent screening or its fitted logarithms.
//!
//! Leading weak-drift internal heating includes a signed correction to thermal
//! exchange. With `beta=kB(Ta-Tb)/[(ma+mb)*S]`, the drift terms through order
//! `w²=|ua-ub|²` are `Qa=K*w²*[mb/(ma+mb)+3*beta/2]` and the complementary Qb.
//! This follows by conditioning the center-of-mass velocity on the relative
//! Gaussian velocity in the same Rutherford/Landau energy moment: its mean
//! is `U+beta*(g-w)`, while `<1/|g|>=sqrt(2/pi)/sqrt(S)*[1-w²/(6*S)+...]`.
//! The correction reduces the zero-drift thermal transfer by `w²/(2*S)`.
//! It can make one drift contribution negative; neither pool is clipped.
//! Their sum is exactly the retained linear-drag dissipation `K*w²`.
//! Zero-drift exchange is separate and must be counted once by its energy owner.

use crate::electron_collision::{AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K, ELECTRON_MASS_KG};
use crate::equilibrium::composition::FrozenIonizedGasStateV1;
use crate::{SpeciesId, ThermochemErrorV1};
use fs_math::det::{exp, ln, sqrt};
use fs_qty::Temperature;
use std::sync::Arc;

mod evaluate;

const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;
// CODATA 2022; these are model constants, not an uncertainty enclosure.
const VACUUM_PERMITTIVITY_F_M: f64 = 8.854_187_818_8e-12;
const PLANCK_J_S: f64 = 6.626_070_15e-34;
const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// Explicit leading-log screening prescription, fixed for a retained model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoulombScreeningV1 {
    /// lambda_D^-2=sum_present n_s q_s²/(eps0 kB T_s), including both pools.
    /// This is static screening; dynamic dielectric screening is not solved.
    StaticAllChargedDebye,
    /// lambda_De^-2=n_e e²/(eps0 kB Te). A populated structural electron is
    /// required whenever any charged species is present.
    ElectronDebye,
}

/// Caller-selected validity limits, not experimental or numerical error bars.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoulombLimitsV1 {
    /// Smallest admitted ln(lambda_D/b90), strictly greater than one.
    pub minimum_coulomb_logarithm: f64,
    /// Smallest admitted (4*pi/3)*n_charged*lambda_D³, greater than one.
    pub minimum_particles_per_debye_sphere: f64,
    /// Largest z_max²e²/(4*pi*eps0*a*kB*T_min), with total-carrier Wigner-Seitz a.
    pub maximum_coulomb_coupling: f64,
    /// Largest n_s*[h/sqrt(2*pi*m_s*kB*T_s)]³; no inferred spin degeneracy.
    pub maximum_thermal_phase_space_density: f64,
    /// Largest kB*T_s/(m_s*c²), screening nonrelativistic thermal motion.
    pub maximum_thermal_energy_rest_energy_ratio: f64,
    /// Largest [hbar/(2*mu*sqrt(S))]/b90, screening classical scattering.
    pub maximum_quantum_to_classical_ratio: f64,
    /// Largest |ua-ub|/min(sqrt(kB*Ta/ma),sqrt(kB*Tb/mb)). Actual drift is
    /// checked by the pair's heating evaluator and by its momentum consumer.
    pub maximum_relative_drift_thermal_ratio: f64,
}
impl CoulombLimitsV1 {
    fn validate(self) -> Result<(), CoulombErrorV1> {
        if !self.minimum_coulomb_logarithm.is_finite()
            || self.minimum_coulomb_logarithm <= 1.0
            || !self.minimum_particles_per_debye_sphere.is_finite()
            || self.minimum_particles_per_debye_sphere <= 1.0
            || [
                self.maximum_coulomb_coupling,
                self.maximum_thermal_phase_space_density,
                self.maximum_thermal_energy_rest_energy_ratio,
                self.maximum_quantum_to_classical_ratio,
                self.maximum_relative_drift_thermal_ratio,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0 || *x >= 1.0)
        {
            return Err(CoulombErrorV1::Invalid(
                "finite classical Coulomb limits: minima>1, maxima in (0,1)",
            ));
        }
        Ok(())
    }
}

/// An original source, physical domain or representability refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum CoulombErrorV1 {
    Invalid(&'static str),
    Domain(&'static str),
    Species {
        species: SpeciesId,
        source: ThermochemErrorV1,
    },
    Arithmetic(&'static str),
    AllocationRefused,
    Cancelled,
}
impl core::fmt::Display for CoulombErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Maxwellian Coulomb model: {self:?}")
    }
}
impl std::error::Error for CoulombErrorV1 {}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Carrier {
    index: usize,
    mass: f64,
    log_mass: f64,
    log_density: f64,
    log_charge: f64,
    electron: bool,
}

/// Immutable source and population binding; clones share carrier metadata.
/// Construction allocates bounded metadata. Successful evaluations and pair
/// visits allocate no trial storage, with at most128 source species and8128 pairs.
#[derive(Clone, Debug, PartialEq)]
pub struct CoulombModelV1 {
    initial_state: FrozenIonizedGasStateV1,
    carriers: Arc<[Carrier]>,
    screening: CoulombScreeningV1,
    limits: CoulombLimitsV1,
}
impl CoulombModelV1 {
    /// Bind actual populations, charge axis, particle masses and original source
    /// support. No LTE projection, seed, neutral collision or transport source
    /// is inserted. Initial temperatures must pass the same evaluation gates.
    pub fn new(
        state: &FrozenIonizedGasStateV1,
        screening: CoulombScreeningV1,
        limits: CoulombLimitsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CoulombErrorV1> {
        limits.validate()?;
        poll(&mut cancelled)?;
        let mut carriers = Vec::new();
        carriers
            .try_reserve_exact(state.amounts().len())
            .map_err(|_| CoulombErrorV1::AllocationRefused)?;
        for (index, (model, amount)) in state
            .thermodynamics()
            .models()
            .iter()
            .zip(state.amounts())
            .enumerate()
        {
            poll(&mut cancelled)?;
            let z = state
                .thermodynamics()
                .charge_vector()
                .get(index)
                .expect("retained source charge axis");
            if z == 0 || amount.value() == 0.0 {
                continue;
            }
            if z.unsigned_abs() > (1u128 << 53) {
                return Err(CoulombErrorV1::Invalid(
                    "charge number cannot be represented exactly",
                ));
            }
            let mass = model.molar_mass().value() / AVOGADRO_PER_MOL;
            let density = state
                .number_density_per_m3(index)
                .expect("retained source amount axis");
            let electron = state.electron_index() == Some(index);
            if !positive(mass)
                || !positive(density)
                || (electron
                    && ((mass - ELECTRON_MASS_KG) / ELECTRON_MASS_KG).abs() > 16.0 * f64::EPSILON)
            {
                return Err(CoulombErrorV1::Invalid(
                    "source particle mass, physical electron mass or density",
                ));
            }
            carriers.push(Carrier {
                index,
                mass,
                log_mass: ln(mass),
                log_density: ln(density),
                log_charge: ln((z as f64).abs()) + ln(ELEMENTARY_CHARGE_C),
                electron,
            });
        }
        if screening == CoulombScreeningV1::ElectronDebye
            && !carriers.is_empty()
            && !carriers.iter().any(|c| c.electron)
        {
            return Err(CoulombErrorV1::Invalid(
                "electron-only screening requires populated structural electrons",
            ));
        }
        let model = Self {
            initial_state: state.clone(),
            carriers: carriers.into(),
            screening,
            limits,
        };
        model.evaluate(
            state.heavy_temperature(),
            state.electron_temperature(),
            &mut cancelled,
        )?;
        Ok(model)
    }
    /// Original source-bound populations and initial temperatures.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial_state
    }
    /// Exact source bank, amounts and volume binding; temperatures may evolve.
    pub fn is_bound_to_populations(&self, state: &FrozenIonizedGasStateV1) -> bool {
        self.initial_state.thermodynamics() == state.thermodynamics()
            && self.initial_state.amounts() == state.amounts()
            && self.initial_state.volume() == state.volume()
    }
    pub const fn limits(&self) -> CoulombLimitsV1 {
        self.limits
    }
    pub const fn screening(&self) -> CoulombScreeningV1 {
        self.screening
    }
    /// Evaluate the complete state at actual Th and Te, without allocating.
    pub fn evaluate(
        &self,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoulombSummaryV1, CoulombErrorV1> {
        self.evaluate_with_pairs(heavy_temperature, electron_temperature, |_| {}, cancelled)
    }
    /// Visit every canonical unordered pair and return the complete summary.
    /// A visitor receives provisional local values: consumers must discard their
    /// accumulation if a later source/domain/cancellation check returns Err.
    pub fn evaluate_with_pairs(
        &self,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        visitor: impl FnMut(CoulombPairV1),
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoulombSummaryV1, CoulombErrorV1> {
        self.evaluate_inner(heavy_temperature, electron_temperature, visitor, cancelled)
    }
}

/// Fully evaluated canonical pair from one bound model and actual temperatures.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoulombPairV1 {
    first_index: usize,
    second_index: usize,
    friction: f64,
    thermal: f64,
    coulomb_logarithm: f64,
    relative_variance: f64,
    maximum_speed_squared: f64,
    first_mass_fraction: f64,
    second_mass_fraction: f64,
    thermal_drift_coefficient: f64,
}
impl CoulombPairV1 {
    pub const fn first_index(&self) -> usize {
        self.first_index
    }
    pub const fn second_index(&self) -> usize {
        self.second_index
    }
    /// Symmetric K_ab=K_ba [kg/(m³ s)], not an independent drag frequency.
    pub const fn friction_coefficient_kg_m3_s(&self) -> f64 {
        self.friction
    }
    /// Symmetric G_ab=G_ba [W/(m³ K)] for zero relative drift.
    pub const fn thermal_conductance_w_m3_k(&self) -> f64 {
        self.thermal
    }
    pub const fn coulomb_logarithm(&self) -> f64 {
        self.coulomb_logarithm
    }
    pub const fn relative_thermal_variance_m2_s2(&self) -> f64 {
        self.relative_variance
    }
    pub const fn maximum_relative_speed_squared_m2_s2(&self) -> f64 {
        self.maximum_speed_squared
    }
    /// Internal heating correction from the ACTUAL squared relative speed.
    /// This excludes zero-drift heat transfer ±G(Ta-Tb). Its signed pool values
    /// add to K*w²; the thermal correction redistributes energy with zero sum.
    pub fn drift_heating(
        &self,
        relative_speed_squared: f64,
    ) -> Result<CoulombDriftHeatingV1, CoulombErrorV1> {
        if !relative_speed_squared.is_finite() || relative_speed_squared < 0.0 {
            return Err(CoulombErrorV1::Invalid(
                "finite nonnegative pair relative speed squared",
            ));
        }
        if relative_speed_squared > self.maximum_speed_squared {
            return Err(CoulombErrorV1::Domain(
                "pair drift exceeds declared Maxwellian thermal-speed limit",
            ));
        }
        let total = self.friction * relative_speed_squared;
        let first_friction = self.second_mass_fraction * total;
        let second_friction = self.first_mass_fraction * total;
        let correction = self.thermal_drift_coefficient * total;
        let first = first_friction + correction;
        let second = second_friction - correction;
        if [
            total,
            first_friction,
            second_friction,
            correction,
            first,
            second,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || (relative_speed_squared > 0.0 && total == 0.0)
        {
            return Err(CoulombErrorV1::Arithmetic("Coulomb drift heating"));
        }
        Ok(CoulombDriftHeatingV1 {
            first,
            second,
            first_friction,
            second_friction,
            correction,
            total,
        })
    }
}

/// Signed per-pool weak-drift heating and its physically distinct components.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoulombDriftHeatingV1 {
    first: f64,
    second: f64,
    first_friction: f64,
    second_friction: f64,
    correction: f64,
    total: f64,
}
impl CoulombDriftHeatingV1 {
    pub const fn first_power_density_w_m3(&self) -> f64 {
        self.first
    }
    pub const fn second_power_density_w_m3(&self) -> f64 {
        self.second
    }
    pub const fn first_friction_power_density_w_m3(&self) -> f64 {
        self.first_friction
    }
    pub const fn second_friction_power_density_w_m3(&self) -> f64 {
        self.second_friction
    }
    /// Equal/opposite change from the zero-drift thermal-transfer moment.
    pub const fn thermal_exchange_correction_to_first_w_m3(&self) -> f64 {
        self.correction
    }
    pub const fn total_power_density_w_m3(&self) -> f64 {
        self.total
    }
}

/// Complete trial diagnostics. Screening/logarithm/occupation are absent for
/// exactly zero carriers; zero coupling and conductance then have exact meaning.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoulombSummaryV1 {
    pair_count: usize,
    electron_heavy_thermal: f64,
    debye_length: Option<f64>,
    minimum_logarithm: Option<f64>,
    debye_occupation: Option<f64>,
    maximum_coupling: f64,
    maximum_phase_space: f64,
    maximum_relativistic: f64,
    maximum_quantum: f64,
}
impl CoulombSummaryV1 {
    pub const fn pair_count(&self) -> usize {
        self.pair_count
    }
    /// Sum of G for every electron/heavy pair [W/(m³ K)], counted once.
    pub const fn electron_heavy_thermal_conductance_w_m3_k(&self) -> f64 {
        self.electron_heavy_thermal
    }
    pub const fn debye_length_m(&self) -> Option<f64> {
        self.debye_length
    }
    pub const fn minimum_coulomb_logarithm(&self) -> Option<f64> {
        self.minimum_logarithm
    }
    pub const fn particles_per_debye_sphere(&self) -> Option<f64> {
        self.debye_occupation
    }
    pub const fn maximum_coulomb_coupling(&self) -> f64 {
        self.maximum_coupling
    }
    pub const fn maximum_thermal_phase_space_density(&self) -> f64 {
        self.maximum_phase_space
    }
    pub const fn maximum_thermal_energy_rest_energy_ratio(&self) -> f64 {
        self.maximum_relativistic
    }
    pub const fn maximum_quantum_to_classical_ratio(&self) -> f64 {
        self.maximum_quantum
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CoulombErrorV1> {
    if cancelled() {
        Err(CoulombErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn from_log(log: f64, name: &'static str) -> Result<f64, CoulombErrorV1> {
    let value = exp(log);
    if positive(value) {
        Ok(value)
    } else {
        Err(CoulombErrorV1::Arithmetic(name))
    }
}
fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let maximum = a.max(b);
    maximum + ln(1.0 + exp(a.min(b) - maximum))
}

#[cfg(test)]
mod tests;
