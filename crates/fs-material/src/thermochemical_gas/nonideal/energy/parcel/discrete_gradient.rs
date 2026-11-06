//! Symmetric discrete Gibbs work on the original homogeneous PR/NASA parcel.
//!
//! At the four corners (T0,nu0), (T1,nu0), (T0,nu1), (T1,nu1), average the
//! two fixed-volume energy/entropy changes into E_T,S_T and the two isothermal
//! changes into E_V,S_V. theta=E_T/S_T and p_bar=(theta*S_V-E_V)/(nu1-nu0)
//! give Delta u+p_bar*Delta nu=theta*Delta s. Original finite source increments
//! cancel formation/entropy references before evaluation; the ideal limit is
//! the existing logarithmic-volume rule. No ideal EOS is imposed on PR.
//!
//! Work pressure is separate from the original material sampled at theta and
//! arithmetic mean specific volume. Its pressure is NEVER replaced by p_bar.
//! Native PR rectangle admission excludes spinodals and attraction-branch
//! crossings; a conservative complete-Cv bound includes the NASA source.
//! Homogeneous metastability is still possible: this is not a phase-equilibrium
//! or latent-heat model. A positive thermal tangent is required at every trial.

use super::{
    IsochoricRealGasParcel, ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel,
};
use crate::thermochemical_gas::energy::parcel::discrete_gradient::{
    GasDiscreteGradientCandidate, GasDiscreteGradientStep,
};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, map_root, poll, root, updated_point,
};
use crate::thermochemical_gas::{failure, positive};
use fs_qty::{Density, Temperature};
use fs_thermochem::eos::peng_robinson::{
    PengRobinsonDensityCaloricEvaluationV1 as EosState,
    caloric::{PengRobinsonCaloricChangeV1, PengRobinsonCaloricReferenceV1},
};

/// The common proposal retains this real inventory and its complete source state.
pub type RealGasDiscreteGradientCandidate =
    GasDiscreteGradientCandidate<IsochoricRealGasParcel, ResolvedThermochemicalRealGasState>;

fn invalid(message: &'static str) -> GasEnergyError {
    GasEnergyError::Invalid(message)
}

const NODES: [f64; 4] = [
    0.06943184420297371,
    0.33000947820757187,
    0.6699905217924281,
    0.9305681557970262,
];
const WEIGHTS: [f64; 4] = [
    0.17392742256872692,
    0.32607257743127305,
    0.32607257743127305,
    0.17392742256872692,
];

struct Grid<'a> {
    model: &'a ThermochemicalPengRobinsonModel,
    reference: PengRobinsonCaloricReferenceV1<'a>,
    initial_volume_leg: PengRobinsonCaloricChangeV1,
    rho0: f64,
    rho1: f64,
    t0: f64,
    dnu: f64,
    molar_mass: f64,
}

struct Trial {
    end: EosState,
    fixed_volume: EosState,
    theta: f64,
    energy: f64,
    pressure: f64,
    pressure_t: f64,
    pressure_nu: f64,
    cv: f64,
}

impl Trial {
    fn thermal_tangent(&self, dnu: f64) -> f64 {
        self.cv + dnu * self.pressure_t
    }
}

impl Grid<'_> {
    fn scalar(&self, temperature: f64, density: f64) -> Result<EosState, GasEnergyError> {
        self.model
            .mechanics()
            .evaluate_density_caloric(Temperature::new(temperature), Density::new(density))
            .map_err(|error| failure("real-gas Gibbs scalar state", error).into())
    }

    fn evaluate(
        &self,
        t: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial, GasEnergyError> {
        poll(cancelled)?;
        let change = if t == self.t0 {
            self.initial_volume_leg
        } else {
            self.reference
                .evaluate(Temperature::new(t))
                .map_err(|error| failure("real-gas Gibbs finite coordinate changes", error))?
        };
        let end = *change.state();
        let fixed_volume = if self.rho0 == self.rho1 {
            end
        } else if t == self.t0 {
            *self.reference.initial_state()
        } else {
            self.scalar(t, self.rho0)?
        };
        let cv = self.model.caloric_from_eos(&end)?.cv_j_kg_k;
        let cv_fixed = self.model.caloric_from_eos(&fixed_volume)?.cv_j_kg_k;
        let source = self.model.source().thermodynamics();
        let ideal_e = source
            .specific_internal_energy_change(Temperature::new(self.t0), Temperature::new(t))
            .map_err(|error| failure("real-gas Gibbs original thermal energy", error))?
            .value();
        poll(cancelled)?;
        let ideal_s = source
            .specific_isochoric_entropy_change(Temperature::new(self.t0), Temperature::new(t))
            .map_err(|error| failure("real-gas Gibbs original thermal entropy", error))?
            .value();
        let (et_r, ev_r) = change.coordinate_internal_energy_changes_j_mol();
        let (st_r, sv_r) = change.coordinate_entropy_changes_j_mol_k();
        let et = ideal_e + et_r / self.molar_mass;
        let st = ideal_s + st_r / self.molar_mass;
        let ev = ev_r / self.molar_mass;
        let sv =
            (change.ideal_isochoric_entropy_density_correction_j_mol_k() + sv_r) / self.molar_mass;
        let theta = if t == self.t0 { self.t0 } else { et / st };
        if ![et, st, ev, sv].iter().all(|v| v.is_finite())
            || !positive(theta)
            || theta < self.t0.min(t) * (1.0 - 32.0 * f64::EPSILON)
            || theta > self.t0.max(t) * (1.0 + 32.0 * f64::EPSILON)
            || (t != self.t0 && (et == 0.0 || st == 0.0))
        {
            return Err(invalid(
                "real-gas Gibbs thermal mean is not finite and represented",
            ));
        }
        let (theta_t, theta_nu) = if t == self.t0 {
            (0.5, 0.0)
        } else if (t - self.t0).abs() <= 1e-3 * self.t0.min(t) {
            self.close_temperature_derivatives(t, theta, st, 0.5 * (cv + cv_fixed), cancelled)?
        } else {
            let initial = self.initial_volume_leg.state();
            let a_nu = 0.5
                * (end.internal_energy_volume_derivative_pa()
                    - initial.internal_energy_volume_derivative_pa());
            let b_nu = 0.5 * (end.pressure_derivatives().0 - initial.pressure_derivatives().0);
            (
                (0.5 * (cv + cv_fixed) / st) * ((t - theta) / t),
                (a_nu - theta * b_nu) / st,
            )
        };
        if !positive(theta_t) || !theta_nu.is_finite() {
            return Err(invalid(
                "real-gas Gibbs temperature sensitivities are unrepresentable",
            ));
        }
        let old = self.reference.initial_state();
        let free_volume = (old
            .molar_volume_m3_per_mol()
            .min(end.molar_volume_m3_per_mol())
            - old.covolume_m3_per_mol())
            / self.molar_mass;
        let (pressure, pressure_t, pressure_nu) = if self.dnu.abs() <= 1e-3 * free_volume {
            self.close_volume_pressure(t, theta, theta_t, theta_nu, &fixed_volume, cancelled)?
        } else {
            let b = sv / self.dnu;
            let a = ev / self.dnu;
            let b_t = 0.5
                * ((end.residual_cv_j_mol_k() - fixed_volume.residual_cv_j_mol_k())
                    / self.molar_mass)
                / (t * self.dnu);
            let b_nu = (0.5
                * (self.initial_volume_leg.state().pressure_derivatives().0
                    + end.pressure_derivatives().0)
                - b)
                / self.dnu;
            let a_nu = (0.5
                * (self
                    .initial_volume_leg
                    .state()
                    .internal_energy_volume_derivative_pa()
                    + end.internal_energy_volume_derivative_pa())
                - a)
                / self.dnu;
            (
                theta * b - a,
                theta_t * b + (theta - t) * b_t,
                theta_nu * b + theta * b_nu - a_nu,
            )
        };
        let trial = Trial {
            end,
            fixed_volume,
            theta,
            energy: et + ev,
            pressure,
            pressure_t,
            pressure_nu,
            cv,
        };
        if !positive(pressure)
            || !pressure_t.is_finite()
            || !pressure_nu.is_finite()
            || !trial.energy.is_finite()
            || !positive(trial.thermal_tangent(self.dnu))
        {
            return Err(invalid(
                "real-gas discrete work pressure or increasing thermal branch refused",
            ));
        }
        poll(cancelled)?;
        Ok(trial)
    }

    // Positive Cv/T moment avoids cancellation in theta_T near equal T.
    // The signed p_TT moment gives theta_nu directly. Four-node Gauss only
    // spans relative temperature <=1e-3 and splits every original NASA seam;
    // residual energies and entropies always use the exact source integrals.
    fn close_temperature_derivatives(
        &self,
        t1: f64,
        theta: f64,
        entropy: f64,
        cv_end: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), GasEnergyError> {
        let delta = t1 - self.t0;
        let mut left = 0.0_f64;
        let mut moments = [0.0; 2];
        for _ in 0..128 {
            poll(cancelled)?;
            let mut right = 1.0_f64;
            for species in self.model.source().thermodynamics().components() {
                for region in species.regions() {
                    let z = (region.temperature_min().value() - self.t0) / delta;
                    if z > left && z < right {
                        right = z;
                    }
                }
            }
            let a = self.t0 + left * delta;
            let b = self.t0 + right * delta;
            let lo = a.min(b);
            let hi = a.max(b);
            for (node, weight) in NODES.into_iter().zip(WEIGHTS) {
                poll(cancelled)?;
                let z = left + (right - left) * node;
                let t = self.t0 + z * delta;
                let sample_t = if lo.next_up() >= hi {
                    lo
                } else {
                    t.max(lo.next_up()).min(hi.next_down())
                };
                let first = self.scalar(sample_t, self.rho0)?;
                let last = if self.rho0 == self.rho1 {
                    first
                } else {
                    self.scalar(sample_t, self.rho1)?
                };
                let cv = 0.5
                    * (self.model.caloric_from_eos(&first)?.cv_j_kg_k
                        + self.model.caloric_from_eos(&last)?.cv_j_kg_k);
                let scale = (right - left) * weight;
                moments[0] += scale * (cv / t) * (1.0 - z);
                moments[1] += scale
                    * ((self.t0 - theta) + z * delta)
                    * last.pressure_temperature_second_derivatives().0;
            }
            if right == 1.0 {
                let divided_entropy = entropy / delta;
                return Ok((
                    (cv_end / t1) * (moments[0] / divided_entropy) / divided_entropy,
                    0.5 * moments[1] / divided_entropy,
                ));
            }
            left = right;
        }
        Err(invalid(
            "too many original thermal seams in real-gas Gibbs tangent; narrow bracket",
        ))
    }

    // Stable limits of the volume divided differences and their derivatives.
    // The four-node rule spans <=1e-3 of the minimum covolume clearance.
    // It resolves the ORIGINAL pressure/Jacobian, never finite differences of
    // rounded material states. Final energy/Gibbs residual gates remain exact.
    fn close_volume_pressure(
        &self,
        t1: f64,
        theta: f64,
        theta_t: f64,
        theta_nu: f64,
        fixed: &EosState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64, f64), GasEnergyError> {
        let initial = self.reference.initial_state();
        let mut pressure = 0.0;
        let mut b = 0.0;
        let mut b_t = 0.0;
        let mut volume_tangent = 0.0;
        for (z, weight) in NODES.into_iter().zip(WEIGHTS) {
            poll(cancelled)?;
            let (first, last) = if self.dnu == 0.0 {
                (*initial, *fixed)
            } else {
                let molar_volume =
                    initial.molar_volume_m3_per_mol() + z * (self.dnu * self.molar_mass);
                let density = self.molar_mass / molar_volume;
                let first = self.scalar(self.t0, density)?;
                let last = if t1 == self.t0 {
                    first
                } else {
                    self.scalar(t1, density)?
                };
                (first, last)
            };
            let (pt0, pv0) = first.pressure_derivatives();
            let (pt1, pv1) = last.pressure_derivatives();
            let (_, ptv0) = first.pressure_temperature_second_derivatives();
            let (ptt1, ptv1) = last.pressure_temperature_second_derivatives();
            b += weight * 0.5 * (pt0 + pt1);
            b_t += weight * 0.5 * ptt1;
            pressure += weight
                * 0.5
                * ((theta - self.t0) * pt0
                    + first.pressure().value()
                    + (theta - t1) * pt1
                    + last.pressure().value());
            volume_tangent += weight
                * 0.5
                * z
                * self.molar_mass
                * ((theta - self.t0) * ptv0 + pv0 + (theta - t1) * ptv1 + pv1);
        }
        Ok((
            pressure,
            theta_t * b + (theta - t1) * b_t,
            theta_nu * b + volume_tangent,
        ))
    }
}

impl IsochoricRealGasParcel {
    /// Propose simultaneous heat and volume work with the original real-gas
    /// discrete Gibbs identity. Positive heat enters this conserved inventory.
    ///
    /// Use step().work_pressure_pa in the enclosing momentum/work equation.
    /// thermodynamic_mean() retains a complete original PR state at theta and
    /// arithmetic mean specific volume; its ordinary EOS pressure is a material
    /// sampling coordinate, not a substituted mechanical work effort.
    ///
    /// The native source screen admits a mechanically stable positive-pressure
    /// temperature/density rectangle and excludes attraction-branch crossings.
    /// The original continuous NASA branch plus the minimum PR residual Cv
    /// must bound COMPLETE Cv positively. This conservative screen can refuse
    /// a larger rectangle even when the sought endpoint is valid: narrow the
    /// temperature bracket or refine volume change. Metastability is permitted;
    /// no global phase stability, flash, reaction or latent heat is inferred.
    ///
    /// Four caloric corners are checked during thermal trials. Complete source
    /// transport resolves at accepted corners and the representative mean only;
    /// there is no cubic solve, source clone or transport query in Newton trials.
    /// Near-coincident derivatives use bounded original-source Gauss moments.
    /// Actual energy, entropy-transfer and Gibbs residuals independently meet
    /// the caller's energy allowance before any proposal is returned.
    pub fn discrete_gradient_heat_volume(
        &self,
        volume_m3: f64,
        heat_j: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasDiscreteGradientCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let t0 = self.state.temperature_k();
        let dv = volume_m3 - self.volume_m3;
        let density = if dv == 0.0 {
            self.density_kg_m3
        } else {
            self.mass_kg / volume_m3
        };
        let dnu = dv / self.mass_kg;
        let q = heat_j / self.mass_kg;
        let tolerance = self.mass_kg * controls.energy_tolerance_j_kg;
        if !positive(volume_m3)
            || !heat_j.is_finite()
            || !positive(density)
            || !dnu.is_finite()
            || !q.is_finite()
            || !positive(tolerance)
            || (heat_j != 0.0 && q == 0.0)
            || (dv != 0.0 && (dnu == 0.0 || density == self.density_kg_m3))
            || t0 < controls.minimum_temperature_k
            || t0 > controls.maximum_temperature_k
        {
            return Err(invalid(
                "real-gas Gibbs step needs represented volume/heat and bracketed initial temperature",
            ));
        }
        let molar_mass = self.model.mechanics().molar_mass().value();
        if dv == 0.0 && heat_j == 0.0 {
            let cv = self.state.specific_heat_cv_j_kg_k();
            let (pt, pv) = self.state.root().pressure_derivatives();
            let dpdv = 0.5 * (molar_mass * pv - t0 * pt * pt / cv) / self.mass_kg;
            let dpdq = 0.5 * (pt / cv) / self.mass_kg;
            if !dpdv.is_finite() || dpdv >= 0.0 || !dpdq.is_finite() {
                return Err(invalid(
                    "stationary real-gas Gibbs pressure tangent is unrepresentable",
                ));
            }
            poll(&mut cancelled)?;
            return Ok(GasDiscreteGradientCandidate {
                endpoint: self.clone(),
                thermodynamic_mean: self.state.clone(),
                step: GasDiscreteGradientStep {
                    initial_volume_m3: self.volume_m3,
                    final_volume_m3: self.volume_m3,
                    supplied_heat_j: 0.0,
                    work_pressure_pa: self.state.pressure_pa(),
                    work_by_gas_j: 0.0,
                    stored_energy_change_j: 0.0,
                    energy_residual_j: 0.0,
                    entropy_change_j_k: 0.0,
                    entropy_transfer_j_k: 0.0,
                    entropy_residual_j_k: 0.0,
                    gibbs_residual_j: 0.0,
                    pressure_final_volume_derivative_pa_m3: dpdv,
                    pressure_heat_derivative_pa_j: dpdq,
                    iterations: 0,
                },
            });
        }
        let branch = self
            .model
            .source()
            .prepare_caloric_branch(controls, &mut cancelled)?;
        let reference = self
            .model
            .caloric_reference(self.density_kg_m3, t0, density)?;
        let residual_cv = reference
            .admit_stable_rectangle([
                Temperature::new(controls.minimum_temperature_k),
                Temperature::new(controls.maximum_temperature_k),
            ])
            .map_err(|error| failure("real-gas Gibbs source rectangle", error))?;
        if !positive((branch.minimum_cv_j_kg_k() + residual_cv / molar_mass).next_down()) {
            return Err(invalid(
                "complete real-gas Cv is not bounded positively on the source rectangle",
            ));
        }
        poll(&mut cancelled)?;
        let initial_volume_leg = reference
            .evaluate(Temperature::new(t0))
            .map_err(|error| failure("real-gas Gibbs initial volume leg", error))?;
        self.model.caloric_from_eos(reference.initial_state())?;
        self.model.caloric_from_eos(initial_volume_leg.state())?;
        let grid = Grid {
            model: &self.model,
            reference,
            initial_volume_leg,
            rho0: self.density_kg_m3,
            rho1: density,
            t0,
            dnu,
            molar_mass,
        };
        let initial = grid.evaluate(t0, &mut cancelled)?;
        let f0 = (initial.energy + initial.pressure * dnu) - q;
        let isothermal_tolerance = (8.0
            * f64::EPSILON
            * initial
                .energy
                .abs()
                .max((initial.pressure * dnu).abs())
                .max(q.abs()))
        .min(0.25 * controls.energy_tolerance_j_kg);
        let (t1, iterations) = if f0.abs() <= isothermal_tolerance {
            (t0, 0)
        } else {
            let solve_tolerance = (0.25 * controls.energy_tolerance_j_kg).min(0.5 * f0.abs());
            if !positive(solve_tolerance) {
                return Err(invalid(
                    "real-gas Gibbs increment is below represented solve resolution",
                ));
            }
            let (lo, hi) = if f0 < 0.0 {
                (t0, controls.maximum_temperature_k)
            } else {
                (controls.minimum_temperature_k, t0)
            };
            let evaluate = |t| {
                let trial = grid.evaluate(t, &mut cancelled)?;
                Ok::<_, GasEnergyError>(root::Sample {
                    residual: (trial.energy + trial.pressure * dnu) - q,
                    tangent: trial.thermal_tangent(dnu),
                    pressure_pa: trial.end.pressure().value(),
                })
            };
            // The physical evaluator polls entry, every source seam/node and
            // return; do not alias its callback through the root driver's poll.
            let solved = root::solve(
                lo,
                hi,
                solve_tolerance,
                controls.maximum_iterations,
                evaluate,
                &mut || false,
            )
            .map_err(map_root)?;
            (solved.temperature_k, solved.iterations)
        };
        let trial = grid.evaluate(t1, &mut cancelled)?;
        // Check the complete original material at both off-diagonal corners;
        // local EOS stability alone is not a transport/source-domain receipt.
        for (t, rho, pressure) in [
            (
                t1,
                self.density_kg_m3,
                trial.fixed_volume.pressure().value(),
            ),
            (
                t0,
                density,
                grid.initial_volume_leg.state().pressure().value(),
            ),
        ] {
            poll(&mut cancelled)?;
            let point = updated_point(&self.point, t, pressure)?;
            self.model.evaluate_density_temperature(&point, rho)?;
        }
        poll(&mut cancelled)?;
        let point = updated_point(&self.point, t1, trial.end.pressure().value())?;
        let state = self.model.evaluate_density_temperature(&point, density)?;
        let mean_density = if dv == 0.0 {
            self.density_kg_m3
        } else {
            self.mass_kg / f64::midpoint(self.volume_m3, volume_m3)
        };
        let mean_point = updated_point(&self.point, trial.theta, self.state.pressure_pa())?;
        let mean = self
            .model
            .evaluate_density_temperature(&mean_point, mean_density)?;
        if state.density_kg_m3().to_bits() != density.to_bits()
            || mean.density_kg_m3().to_bits() != mean_density.to_bits()
            || mean.temperature_k().to_bits() != trial.theta.to_bits()
        {
            return Err(invalid(
                "real-gas Gibbs reconstruction changed its conserved sampling coordinates",
            ));
        }
        poll(&mut cancelled)?;
        let stored = self.mass_kg
            * self
                .model
                .specific_internal_energy_change_between_densities(
                    self.density_kg_m3,
                    t0,
                    density,
                    t1,
                )?;
        let entropy = self.mass_kg
            * self.model.specific_entropy_change_between_densities(
                self.density_kg_m3,
                t0,
                density,
                t1,
            )?;
        let work = trial.pressure * dv;
        let transferred = heat_j / trial.theta;
        let energy_residual = (stored + work) - heat_j;
        let gibbs_residual = (stored + work) - trial.theta * entropy;
        let entropy_residual = entropy - transferred;
        let thermal_tangent = trial.thermal_tangent(dnu);
        let dpdv = ((trial.pressure_nu * trial.cv
            - trial.pressure_t
                * (trial.end.internal_energy_volume_derivative_pa() + trial.pressure))
            / thermal_tangent)
            / self.mass_kg;
        let dpdq = (trial.pressure_t / thermal_tangent) / self.mass_kg;
        if ![
            stored,
            entropy,
            work,
            transferred,
            energy_residual,
            gibbs_residual,
            entropy_residual,
            dpdv,
            dpdq,
        ]
        .iter()
        .all(|v| v.is_finite())
            || (dv != 0.0 && work == 0.0)
            || dpdv >= 0.0
            || energy_residual.abs() > tolerance
            || gibbs_residual.abs() > tolerance
            || entropy_residual.abs() > tolerance / trial.theta
        {
            return Err(invalid(
                "real-gas Gibbs energy, entropy or pressure sensitivity budget refused",
            ));
        }
        let mut endpoint = self.clone();
        endpoint.point = point;
        endpoint.state = state;
        endpoint.density_kg_m3 = density;
        endpoint.volume_m3 = volume_m3;
        let step = GasDiscreteGradientStep {
            initial_volume_m3: self.volume_m3,
            final_volume_m3: volume_m3,
            supplied_heat_j: heat_j,
            work_pressure_pa: trial.pressure,
            work_by_gas_j: work,
            stored_energy_change_j: stored,
            energy_residual_j: energy_residual,
            entropy_change_j_k: entropy,
            entropy_transfer_j_k: transferred,
            entropy_residual_j_k: entropy_residual,
            gibbs_residual_j: gibbs_residual,
            pressure_final_volume_derivative_pa_m3: dpdv,
            pressure_heat_derivative_pa_j: dpdq,
            iterations,
        };
        poll(&mut cancelled)?;
        Ok(GasDiscreteGradientCandidate {
            endpoint,
            thermodynamic_mean: mean,
            step,
        })
    }
}
