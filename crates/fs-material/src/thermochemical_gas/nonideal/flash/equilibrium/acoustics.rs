//! Equilibrium acoustic response of a multicomponent liquid/vapor system.
//!
//! Every perturbation retains the original total component inventories while
//! re-equilibrating phase compositions and amounts. Pressure, temperature,
//! chemical potentials and velocity equilibrate during the wave. This explicit
//! homogeneous-equilibrium limit supplies no relaxation rate or attenuation.
//!
//! A fixed nine-PT stencil resolves the center and two centered step sizes in
//! each coordinate. Original entropy increments avoid subtracting arbitrary
//! species datums. Refinement, phase count, Maxwell reciprocity and stability
//! must all pass before a response is published. The refinement estimate is
//! numerical evidence, not an interval-certified error bound. Bubble/dew
//! crossings, critical singularities and PT-indeterminate coexistence refuse.

use super::super::{CaloricSample, E, ThermochemicalPengRobinsonModel, poll, positive};
use super::ResolvedThermochemicalPhaseState;
use crate::MaterialError;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::det;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::eos::peng_robinson::{
    flash::PengRobinsonFlashOptionsV1,
    phase_screen::{PengRobinsonPhaseScreenOptionsV1, PengRobinsonPhaseScreenV1},
    stability::PengRobinsonStabilityOptionsV1,
};

pub use crate::thermochemical_gas::nonideal::equilibrium::acoustics::EquilibriumThermodynamicDerivatives;

/// Center plus two signed perturbations at two resolutions in each coordinate.
pub const PHASE_ACOUSTIC_PT_EVALUATIONS: usize = 9;

/// Bounded phase work and a caller-declared differentiable neighborhood.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseEquilibriumAcousticControls {
    /// Original phase screening and fugacity policy for all nine PT samples.
    pub phase: PengRobinsonPhaseScreenOptionsV1,
    /// Coarse centered temperature step divided by center T; in (0, 0.05].
    pub relative_temperature_step: f64,
    /// Coarse centered pressure step divided by center p; in (0, 0.05].
    pub relative_pressure_step: f64,
    /// Maximum normalized step-refinement estimate and Maxwell discrepancy.
    /// Derivatives use their magnitude or the corresponding ideal-gas scale,
    /// whichever is larger, so a vanishing cross derivative remains resolvable.
    /// Derived heat capacities and pressure responses must also converge.
    pub relative_tolerance: f64,
}
impl Default for PhaseEquilibriumAcousticControls {
    fn default() -> Self {
        Self {
            phase: PengRobinsonPhaseScreenOptionsV1 {
                stability: PengRobinsonStabilityOptionsV1 {
                    stationarity_tolerance: 1e-12,
                    tangent_plane_tolerance: 1e-12,
                    ..Default::default()
                },
                flash: PengRobinsonFlashOptionsV1 {
                    log_k_tolerance: 1e-12,
                    ..Default::default()
                },
            },
            relative_temperature_step: 1e-3,
            relative_pressure_step: 1e-3,
            relative_tolerance: 1e-4,
        }
    }
}
impl PhaseEquilibriumAcousticControls {
    /// Check finite steps, accuracy and original phase-solver work limits.
    pub fn validate(self) -> Result<(), E> {
        self.phase.validate().map_err(E::Phase)?;
        if [
            self.relative_temperature_step,
            self.relative_pressure_step,
            self.relative_tolerance,
        ]
        .into_iter()
        .any(|x| !positive(x))
            || self.relative_temperature_step > 0.05
            || self.relative_pressure_step > 0.05
            || self.relative_tolerance > 0.01
        {
            return Err(E::Invalid(
                "invalid equilibrium acoustic derivative controls",
            ));
        }
        Ok(())
    }

    /// Conservative cubic-EOS ceiling for the complete response construction.
    /// Consumers of the returned acoustic state perform no phase solves.
    pub fn maximum_eos_evaluation_bound(self) -> Result<usize, E> {
        self.validate()?;
        Ok(PHASE_ACOUSTIC_PT_EVALUATIONS
            * self
                .phase
                .maximum_eos_evaluation_bound()
                .map_err(E::Phase)?)
    }
}

/// Original resolved phase state and its explicitly chosen equilibrium response.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseEquilibriumAcoustics {
    source: ResolvedThermochemicalPhaseState,
    acoustic: FluidAcousticState,
    derivatives: EquilibriumThermodynamicDerivatives,
    isobaric_heat_capacity_j_kg_k: f64,
    refinement_error: f64,
    maxwell_residual: f64,
    eos_evaluation_upper_bound: usize,
}
impl PhaseEquilibriumAcoustics {
    /// Center state, unchanged component inventory and complete source receipts.
    #[must_use]
    pub const fn source(&self) -> &ResolvedThermochemicalPhaseState {
        &self.source
    }
    /// Coherent density, sound speed, isentropic bulk modulus and impedance.
    #[must_use]
    pub const fn state(&self) -> FluidAcousticState {
        self.acoustic
    }
    /// Equilibrating pressure and entropy Jacobian in (T,rho) coordinates.
    #[must_use]
    pub const fn derivatives(&self) -> EquilibriumThermodynamicDerivatives {
        self.derivatives
    }
    /// Complete equilibrium T*(ds/dT)_p [J/(kg K)], including redistribution.
    #[must_use]
    pub const fn isobaric_heat_capacity_j_kg_k(&self) -> f64 {
        self.isobaric_heat_capacity_j_kg_k
    }
    /// Largest normalized centered-step error estimate before extrapolation.
    #[must_use]
    pub const fn relative_refinement_error(&self) -> f64 {
        self.refinement_error
    }
    /// Independently differentiated ds/dp + dv/dT, normalized by its scale.
    #[must_use]
    pub const fn relative_maxwell_residual(&self) -> f64 {
        self.maxwell_residual
    }
    /// Exact number of original PT phase selections for an accepted response.
    #[must_use]
    pub const fn pt_evaluations(&self) -> usize {
        PHASE_ACOUSTIC_PT_EVALUATIONS
    }
    /// Sum of observed screen work plus the original split-work ceilings.
    /// This is an upper bound, not an observed total cubic-EOS count.
    #[must_use]
    pub const fn eos_evaluation_upper_bound(&self) -> usize {
        self.eos_evaluation_upper_bound
    }
}
impl FluidAcoustics for PhaseEquilibriumAcoustics {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

#[derive(Clone, Copy)]
struct PtDerivatives {
    volume_temperature: f64,
    volume_pressure: f64,
    entropy_temperature: f64,
    entropy_pressure: f64,
}
impl PtDerivatives {
    fn values(self) -> [f64; 4] {
        [
            self.volume_temperature,
            self.volume_pressure,
            self.entropy_temperature,
            self.entropy_pressure,
        ]
    }
    fn extrapolate(self, coarse: Self) -> Self {
        let r = |fine: f64, coarse: f64| fine + (fine - coarse) / 3.0;
        Self {
            volume_temperature: r(self.volume_temperature, coarse.volume_temperature),
            volume_pressure: r(self.volume_pressure, coarse.volume_pressure),
            entropy_temperature: r(self.entropy_temperature, coarse.entropy_temperature),
            entropy_pressure: r(self.entropy_pressure, coarse.entropy_pressure),
        }
    }
}

fn normalized_difference(a: f64, b: f64, scale: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(scale)
}

fn response(
    t: f64,
    p: f64,
    v: f64,
    d: PtDerivatives,
) -> Result<(f64, EquilibriumThermodynamicDerivatives), E> {
    // At fixed overall composition: ds = Cp/T dT - v_T dp.
    // After independently checking Maxwell reciprocity, use the symmetric
    // thermodynamic Jacobian; never average frozen phase heat capacities.
    let cp = t * d.entropy_temperature;
    let cv = cp + t * d.volume_temperature * (d.volume_temperature / d.volume_pressure);
    let pt = -d.volume_temperature / d.volume_pressure;
    let pr = -(v / d.volume_pressure) * v;
    let tr = (t / cv) * pt * v * v;
    let c2 = pr + pt * tr;
    let pu = pt / cv;
    let pru = pr + pu * (t * pt - p) * v * v;
    let sr = -pt * v * v;
    if d.values().into_iter().any(|x| !x.is_finite())
        || [cp, cv, pr, c2].into_iter().any(|x| !positive(x))
        || [pt, tr, pu, pru, sr, cv / t]
            .into_iter()
            .any(|x| !x.is_finite())
        || d.volume_pressure >= 0.0
    {
        return Err(E::Invalid(
            "unstable or unrepresentable equilibrium acoustic derivatives",
        ));
    }
    Ok((
        cp,
        EquilibriumThermodynamicDerivatives {
            isochoric_heat_capacity_j_kg_k: cv,
            pressure_temperature_at_density_pa_k: pt,
            pressure_density_at_temperature_m2_s2: pr,
            entropy_temperature_at_density_j_kg_k2: cv / t,
            entropy_density_at_temperature_j_m3_kg2_k: sr,
            temperature_density_at_entropy_k_m3_kg: tr,
            pressure_density_at_entropy_m2_s2: c2,
            pressure_energy_at_density_kg_m3: pu,
            pressure_density_at_energy_m2_s2: pru,
        },
    ))
}

impl ThermochemicalPengRobinsonModel {
    /// Construct equilibrium acoustics at fixed total component inventory.
    ///
    /// Eight phase-re-equilibrating PT probes provide a converged local
    /// thermodynamic surface; c²=(dp/drho)_(s,z) follows from its Jacobian.
    /// No nested PS/SV inversion, phase-weighted Cp or frozen Wood law is used.
    /// All probes must retain the center's phase count. This selects the
    /// instantaneous phase/thermal/velocity-equilibrium limit, not its physical
    /// frequency range; slip, capillarity and finite-rate relaxation are absent.
    pub fn phase_equilibrium_acoustics(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        controls: PhaseEquilibriumAcousticControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEquilibriumAcoustics, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let t = temperature.value();
        let p = pressure.value();
        if !positive(t) || !positive(p) {
            return Err(E::Invalid(
                "positive equilibrium acoustic temperature and pressure required",
            ));
        }
        let dt = t * controls.relative_temperature_step;
        let dp = p * controls.relative_pressure_step;
        for (center, step) in [(t, dt), (p, dp)] {
            if !positive(center - step)
                || !(center - step < center - 0.5 * step
                    && center - 0.5 * step < center
                    && center < center + 0.5 * step
                    && center + 0.5 * step < center + step)
                || !(center + step).is_finite()
            {
                return Err(E::Invalid(
                    "equilibrium acoustic stencil is not representable",
                ));
            }
        }
        self.admit_caloric_flash(
            pressure,
            [Temperature::new(t - dt), Temperature::new(t + dt)],
            None,
            &mut cancelled,
        )?;
        let center =
            self.caloric_phase_sample(temperature, pressure, controls.phase, &mut cancelled)?;
        let mut eos_work = center.eos_evaluation_upper_bound;
        let mass = self.mechanics().molar_mass().value();
        let v = center.molar_volume / mass;
        let rs = R / mass;
        if !positive(v) || !positive(rs) {
            return Err(E::Invalid(
                "equilibrium acoustic reference scales are unrepresentable",
            ));
        }
        let split = center.phase.two_phase().is_some();
        let mut stencil = |factor: f64| -> Result<PtDerivatives, E> {
            let points = [
                (t - factor * dt, p),
                (t + factor * dt, p),
                (t, p - factor * dp),
                (t, p + factor * dp),
            ];
            let mut sample =
                |t: f64, p: f64| -> Result<CaloricSample<PengRobinsonPhaseScreenV1>, E> {
                    let state = self.caloric_phase_sample(
                        Temperature::new(t),
                        Pressure::new(p),
                        controls.phase,
                        &mut cancelled,
                    )?;
                    eos_work += state.eos_evaluation_upper_bound;
                    if state.phase.two_phase().is_some() != split {
                        return Err(E::Invalid(
                            "equilibrium acoustic stencil crosses a phase boundary",
                        ));
                    }
                    Ok(state)
                };
            let t_low = sample(points[0].0, points[0].1)?;
            let t_high = sample(points[1].0, points[1].1)?;
            let p_low = sample(points[2].0, points[2].1)?;
            let p_high = sample(points[3].0, points[3].1)?;
            let t_span = points[1].0 - points[0].0;
            let p_span = points[3].1 - points[2].1;
            Ok(PtDerivatives {
                volume_temperature: ((t_high.molar_volume - t_low.molar_volume) / mass) / t_span,
                volume_pressure: ((p_high.molar_volume - p_low.molar_volume) / mass) / p_span,
                entropy_temperature: self.caloric_entropy_change(&t_low, &t_high)? / t_span,
                entropy_pressure: self.caloric_entropy_change(&p_low, &p_high)? / p_span,
            })
        };
        let coarse = stencil(1.0)?;
        let fine = stencil(0.5)?;
        let (coarse_cp, coarse_response) = response(t, p, v, coarse)?;
        let (fine_cp, fine_response) = response(t, p, v, fine)?;
        let mut error: f64 = 0.0;
        for ((a, b), scale) in
            coarse
                .values()
                .into_iter()
                .zip(fine.values())
                .zip([v / t, v / p, rs / t, rs / p])
        {
            error = error.max(normalized_difference(a, b, scale) / 3.0);
        }
        let derived = |cp: f64, d: EquilibriumThermodynamicDerivatives| {
            [
                cp,
                d.isochoric_heat_capacity_j_kg_k,
                d.pressure_temperature_at_density_pa_k,
                d.pressure_density_at_temperature_m2_s2,
                d.pressure_density_at_entropy_m2_s2,
                d.pressure_energy_at_density_kg_m3,
                d.pressure_density_at_energy_m2_s2,
            ]
        };
        for ((a, b), scale) in derived(coarse_cp, coarse_response)
            .into_iter()
            .zip(derived(fine_cp, fine_response))
            .zip([rs, rs, p / t, p * v, p * v, 1.0 / v, p * v])
        {
            error = error.max(normalized_difference(a, b, scale) / 3.0);
        }
        if !error.is_finite() || error > controls.relative_tolerance {
            return Err(E::Invalid(
                "equilibrium acoustic derivatives did not converge",
            ));
        }
        let accepted = fine.extrapolate(coarse);
        let maxwell = normalized_difference(
            accepted.entropy_pressure,
            -accepted.volume_temperature,
            v / t,
        );
        if !maxwell.is_finite() || maxwell > controls.relative_tolerance {
            return Err(E::Invalid(
                "equilibrium acoustic Maxwell reciprocity is unresolved",
            ));
        }
        let (cp, derivatives) = response(t, p, v, accepted)?;
        let acoustic = FluidAcousticState::try_new(
            mass / center.molar_volume,
            det::sqrt(derivatives.pressure_density_at_entropy_m2_s2),
        )?;
        let source = self.publish_caloric_phase(center, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(PhaseEquilibriumAcoustics {
            source,
            acoustic,
            derivatives,
            isobaric_heat_capacity_j_kg_k: cp,
            refinement_error: error,
            maxwell_residual: maxwell,
            eos_evaluation_upper_bound: eos_work,
        })
    }
}
