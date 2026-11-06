//! Isotropic test-ion energy moments in a classical Maxwellian background.
//!
//! These are the finite-mass Landau/Fokker–Planck moments, including energy
//! diffusion, not the loss of squared mean velocity. For a background species
//! b, x = m_b E/(m_t k_B T_b), Psi = erf(sqrt(x)) - 2 sqrt(x) exp(-x)/sqrt(pi),
//! and nu0 = 4 pi n_b C_tb² ln(Lambda)/(m_t² v³), C_tb = q_t q_b/(4 pi eps0):
//!
//! A_b = -2 E nu0 [(m_t/m_b) Psi - Psi'];
//! D_b = 2 (m_t/m_b) E k_B T_b nu0 Psi, where D is HALF the energy variance rate.
//!
//! Each prepared target freezes its Coulomb logarithm at the background thermal
//! reference S = k_B(Th/m_t + T_b/m_b), using the parent's screening and cutoff.
//! This is a leading-log closure, not velocity-dependent screening or a quantum
//! collision integral. A frozen logarithm preserves exact Maxwellian detailed
//! balance: A = D' + D(1/(2E) - 1/(k_B T_b)). The zero-energy limit heats rather
//! than absorbs particles. A kinetic consumer must retain diffusion and account
//! separately for particle number, charge, energy and background backreaction.
//!
//! Reference: NRL Plasma Formulary (2023), “Collisions and Transport”, pp. 31–35,
//! https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf
//! The test population must be dilute; no neutral stopping, strong coupling,
//! relativistic motion, collective instability or anisotropic transport is added.

use super::*;
use fs_math::det::erf;
use fs_qty::Energy;

/// Original source-bound test ion and admitted Maxwellian background model.
#[derive(Clone, Debug, PartialEq)]
pub struct TestParticleCoulombV1 {
    background: CoulombModelV1,
    projectile_index: usize,
    projectile_mass: f64,
    projectile_log_charge: f64,
    maximum_energy_rest_energy_ratio: f64,
    maximum_energy: Energy,
}

impl TestParticleCoulombV1 {
    /// Select a positive heavy ion by its original source-bank index. An absent
    /// product is allowed; its original source support is checked when preparing.
    pub fn new(
        background: CoulombModelV1,
        projectile_index: usize,
        maximum_energy_rest_energy_ratio: f64,
    ) -> Result<Self, CoulombErrorV1> {
        let state = background.initial_state();
        let charge = state.thermodynamics().charge_vector().get(projectile_index);
        if state.electron_index() == Some(projectile_index)
            || !charge.is_some_and(|z| z > 0 && z.unsigned_abs() <= (1u128 << 53))
            || !positive(maximum_energy_rest_energy_ratio)
            || maximum_energy_rest_energy_ratio >= 1.0
        {
            return Err(CoulombErrorV1::Invalid(
                "original positive heavy projectile and energy/rest-energy limit in (0,1)",
            ));
        }
        let projectile_mass = state.thermodynamics().models()[projectile_index]
            .molar_mass()
            .value()
            / AVOGADRO_PER_MOL;
        if !positive(projectile_mass) {
            return Err(CoulombErrorV1::Invalid("original projectile particle mass"));
        }
        let maximum_energy = Energy::new(from_log(
            ln(maximum_energy_rest_energy_ratio)
                + ln(projectile_mass)
                + 2.0 * ln(SPEED_OF_LIGHT_M_S),
            "maximum nonrelativistic projectile energy",
        )?);
        let projectile_log_charge =
            ln(charge.expect("validated projectile charge") as f64) + ln(ELEMENTARY_CHARGE_C);
        Ok(Self {
            background,
            projectile_index,
            projectile_mass,
            projectile_log_charge,
            maximum_energy_rest_energy_ratio,
            maximum_energy,
        })
    }

    pub const fn background(&self) -> &CoulombModelV1 {
        &self.background
    }
    pub const fn projectile_index(&self) -> usize {
        self.projectile_index
    }
    pub const fn maximum_energy_rest_energy_ratio(&self) -> f64 {
        self.maximum_energy_rest_energy_ratio
    }
    pub const fn projectile_mass_kg(&self) -> f64 {
        self.projectile_mass
    }
    pub const fn maximum_energy(&self) -> Energy {
        self.maximum_energy
    }

    /// Prepare the retained original populations at current pool temperatures.
    pub fn prepare(
        &self,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedTestParticleCoulombV1, CoulombErrorV1> {
        self.prepare_point(
            self.background.initial_state().amounts(),
            self.background.initial_state().volume(),
            heavy_temperature,
            electron_temperature,
            cancelled,
        )
    }

    /// Refresh all screening, density and target coefficients once per physical
    /// stage without allocating. Source identity and species order are retained.
    /// As with the parent's raw-point API, the caller owns hybrid charge/energy
    /// closure: thermal-background amounts may exclude a separate fast-ion bank.
    pub fn prepare_point(
        &self,
        amounts: &[Amount],
        volume: Volume,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedTestParticleCoulombV1, CoulombErrorV1> {
        let summary = self.background.evaluate_point(
            amounts,
            volume,
            heavy_temperature,
            electron_temperature,
            &mut cancelled,
        )?;
        let state = self.background.initial_state();
        let source = &state.thermodynamics().models()[self.projectile_index];
        source
            .evaluate_properties(heavy_temperature)
            .map_err(|source_error| CoulombErrorV1::Species {
                species: source.species().clone(),
                source: source_error,
            })?;
        let log_mt = ln(self.projectile_mass);
        let log_kt = ln(BOLTZMANN_J_PER_K) + ln(heavy_temperature.value());
        if log_kt - log_mt - 2.0 * ln(SPEED_OF_LIGHT_M_S)
            > ln(self
                .background
                .limits()
                .maximum_thermal_energy_rest_energy_ratio)
        {
            return Err(CoulombErrorV1::Domain(
                "projectile thermal-reference relativity limit",
            ));
        }
        let mut prepared = PreparedTestParticleCoulombV1 {
            targets: [Target::EMPTY; MAX_EQUILIBRIUM_SPECIES_V1],
            target_count: 0,
            temperatures: [heavy_temperature, electron_temperature],
            maximum_energy: self.maximum_energy,
            minimum_logarithm: None,
        };
        let log_4pi = ln(4.0 * core::f64::consts::PI);
        for (index, amount) in amounts.iter().copied().enumerate() {
            poll(&mut cancelled)?;
            let Some(target) = Carrier::from_point(state, index, amount, volume)? else {
                continue;
            };
            let pool = usize::from(target.electron);
            let log_kb_tb = ln(BOLTZMANN_J_PER_K) + ln(prepared.temperatures[pool].value());
            let log_s = log_add(log_kt - log_mt, log_kb_tb - target.log_mass);
            let mass_ratio =
                self.projectile_mass.min(target.mass) / self.projectile_mass.max(target.mass);
            let log_mu = log_mt.min(target.log_mass) - ln(1.0 + mass_ratio);
            let log_c = self.projectile_log_charge + target.log_charge
                - log_4pi
                - ln(VACUUM_PERMITTIVITY_F_M);
            let log_b90 = log_c - log_mu - log_s;
            let log_quantum = ln(PLANCK_J_S) - log_4pi + 0.5 * log_s - log_c;
            let log_minimum = match self.background.cutoff() {
                CoulombShortRangeCutoffV1::ClassicalClosestApproach => log_b90,
                CoulombShortRangeCutoffV1::QuantumDiffraction => log_b90 + log_quantum.max(0.0),
            };
            let logarithm = ln(summary.debye_length_m().ok_or(CoulombErrorV1::Arithmetic(
                "populated target without background screening",
            ))?) - log_minimum;
            if !logarithm.is_finite()
                || logarithm < self.background.limits().minimum_coulomb_logarithm
            {
                return Err(CoulombErrorV1::Domain(
                    "test-ion Coulomb logarithm below limit",
                ));
            }
            if self.background.cutoff() == CoulombShortRangeCutoffV1::ClassicalClosestApproach
                && log_quantum > ln(self.background.limits().maximum_quantum_to_classical_ratio)
            {
                return Err(CoulombErrorV1::Domain(
                    "test-ion quantum diffraction exceeds classical limit",
                ));
            }
            // nu0 = rate_scale / x^(3/2); all fixed logarithms are paid once.
            let log_rate = log_4pi + target.log_density + 2.0 * log_c + ln(logarithm)
                - 2.0 * log_mt
                - 1.5 * (ln(2.0) + log_kb_tb - target.log_mass);
            let log_r = log_mt - target.log_mass;
            let log_drag = ln(2.0) + log_r + log_rate;
            prepared.targets[prepared.target_count] = Target {
                pool,
                energy_scale: from_log(log_r + log_kb_tb, "test-ion target energy scale")?,
                heating_scale: from_log(log_drag + log_kb_tb, "test-ion heating scale")?,
                drag_scale: from_log(log_drag, "test-ion energy drag scale")?,
                velocity_scale: from_log(
                    log_rate + log_add(0.0, log_r),
                    "test-ion velocity drag scale",
                )?,
            };
            prepared.target_count += 1;
            prepared.minimum_logarithm = Some(
                prepared
                    .minimum_logarithm
                    .map_or(logarithm, |previous| previous.min(logarithm)),
            );
        }
        poll(&mut cancelled)?;
        Ok(prepared)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Target {
    pool: usize,
    energy_scale: f64,
    heating_scale: f64,
    drag_scale: f64,
    velocity_scale: f64,
}
impl Target {
    const EMPTY: Self = Self {
        pool: 0,
        energy_scale: 0.0,
        heating_scale: 0.0,
        drag_scale: 0.0,
        velocity_scale: 0.0,
    };
}

/// Fixed-size prepared workspace. Every energy evaluation is bounded by the
/// number of populated charged targets (at most 128), with no allocation.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedTestParticleCoulombV1 {
    targets: [Target; MAX_EQUILIBRIUM_SPECIES_V1],
    target_count: usize,
    temperatures: [Temperature; 2],
    maximum_energy: Energy,
    minimum_logarithm: Option<f64>,
}

/// Per-projectile energy moments; pool order is [heavy, electron]. Positive A
/// heats the projectile and removes energy from that background pool. D enters
/// `d_t f = -d_E(A f) + d_E²(D f)`; it is not the full variance rate B = 2D.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TestParticleEnergyCoefficientsV1 {
    pub mean_energy_drift_j_s: [f64; 2],
    pub energy_diffusion_j2_s: [f64; 2],
    /// Decay rate of mean velocity, distinct from kinetic-energy loss.
    pub velocity_drag_frequency_s: [f64; 2],
}

impl PreparedTestParticleCoulombV1 {
    pub const fn temperatures(&self) -> [Temperature; 2] {
        self.temperatures
    }
    pub const fn target_count(&self) -> usize {
        self.target_count
    }
    pub const fn minimum_coulomb_logarithm(&self) -> Option<f64> {
        self.minimum_logarithm
    }
    pub const fn maximum_energy(&self) -> Energy {
        self.maximum_energy
    }

    /// Full finite-mass energy drift and diffusion, including E=0. The model is
    /// nonrelativistic and rejects energies above its explicit rest-energy bound.
    pub fn evaluate(
        &self,
        energy: Energy,
    ) -> Result<TestParticleEnergyCoefficientsV1, CoulombErrorV1> {
        let e = energy.value();
        if !e.is_finite() || e < 0.0 {
            return Err(CoulombErrorV1::Invalid(
                "finite nonnegative test-ion kinetic energy",
            ));
        }
        if e > self.maximum_energy.value() {
            return Err(CoulombErrorV1::Domain("test-ion kinetic/rest-energy limit"));
        }
        let mut sums = [[0.0; 2]; 3];
        let mut corrections = [[0.0; 2]; 3];
        let inv_sqrt_pi = 1.0 / sqrt(core::f64::consts::PI);
        for target in &self.targets[..self.target_count] {
            let x = e / target.energy_scale;
            if !x.is_finite() {
                return Err(CoulombErrorV1::Arithmetic("test-ion thermal speed ratio"));
            }
            let exponential = exp(-x);
            let (f, ef) = if x <= 0.25 {
                // Psi/x^(3/2) is regular at zero. Twelve series terms have
                // truncation below 2e-17 on this interval; no subtractive loss.
                let mut term = 1.0;
                let mut series = term;
                for k in 1..12 {
                    term *= -x / f64::from(k) * f64::from(2 * k + 1) / f64::from(2 * k + 3);
                    series += term;
                }
                let f = (4.0 / 3.0) * inv_sqrt_pi * series;
                (f, e * f)
            } else {
                let u = sqrt(x);
                let psi = erf(u) - 2.0 * inv_sqrt_pi * u * exponential;
                // E*F is evaluated independently to retain the high-energy
                // energy-loss tail even when F itself is extremely small.
                ((psi / x) / u, target.energy_scale * (psi / u))
            };
            let values = [
                target.heating_scale * (2.0 * inv_sqrt_pi * exponential) - target.drag_scale * ef,
                target.heating_scale * ef,
                target.velocity_scale * f,
            ];
            for (moment, value) in values.into_iter().enumerate() {
                if !value.is_finite() || (moment != 0 && value < 0.0) {
                    return Err(CoulombErrorV1::Arithmetic("test-ion energy moment"));
                }
                let sum = &mut sums[moment][target.pool];
                let next = *sum + value;
                corrections[moment][target.pool] += if sum.abs() >= value.abs() {
                    (*sum - next) + value
                } else {
                    (value - next) + *sum
                };
                *sum = next;
            }
        }
        for (sum, correction) in sums.iter_mut().flatten().zip(corrections.iter().flatten()) {
            *sum += correction;
            if !sum.is_finite() {
                return Err(CoulombErrorV1::Arithmetic("summed test-ion energy moment"));
            }
        }
        Ok(TestParticleEnergyCoefficientsV1 {
            mean_energy_drift_j_s: sums[0],
            energy_diffusion_j2_s: sums[1],
            velocity_drag_frequency_s: sums[2],
        })
    }
}

#[cfg(test)]
mod tests;
