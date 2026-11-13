//! Bounded tabulated hardening with native associative Perzyna kinetics.
//!
//! For the existing effective J2 stress q and Prager modulus H_k, one return
//! solves `(3 mu + H_k) d + Y(alpha+d)-Y(alpha) + f_v(d) = q_trial-Y(alpha)`.
//! Here `f_v(d)=s0[eta*d/(dt*s0)]^(1/n)`, with n>=1 and the actual physical dt.
//! The original tensor return receives the combined secant, while its tangent
//! receives the endpoint hardening slope and `f_v/(n*d)`. No second tensor flow
//! or extra material history is introduced.
//!
//! The linear-rate case uses the table's exact piecewise-linear return. The
//! power-law case first locates the root's measured segment and then solves one
//! monotone scalar equation in logarithmic increment. Bracket trials never
//! extrapolate the supplied curve. Point operations allocate no memory.
//!
//! Stored hardening is the original integral of `Y-Y(0)`. Physical heat is
//! `(Y(0)+f_v)*d`; the reduced kinetic potential is `n/(n+1)*f_v*d`, and the
//! convex hardening remainder remains numerical loss. This energetic partition
//! is explicitly selected, not inferred from measured flow stress. The model
//! retains the native small-strain Cauchy-stress convention and excludes
//! sub-yield creep, evolving thermal coefficients and nonlinear backstress.

use super::{TabulatedError, TabulatedJ2, TabulatedJ2Response};
use crate::plastic::energy::{EnergyError, EnergyIncrement, EnergyTolerances, StoredEnergy};
use crate::plastic::rate::PowerLawOverstressJ2;
use crate::plastic::{J2Plasticity, J2State, Update};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard};
use fs_math::det;

/// Same-card measured hardening and independently selected flow kinetics.
pub mod source;
/// Exact free-sided reduction with the original native axial J2 history.
pub mod axial;

/// One immutable table and one explicitly selected physical flow interval.
/// Preparing a later interval does not reset the original native J2 history.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedOverstressJ2 {
    backbone: TabulatedJ2,
    kinetics: PowerLawOverstressJ2,
    linear_selection: bool,
    linear_overstress_modulus_pa: f64,
    log_kinetic_coefficient: f64,
}

impl TabulatedOverstressJ2 {
    /// Select `alpha_dot=max(f,0)/eta` at the supplied actual duration.
    /// No arbitrary stress scale is exposed by this linear selection.
    pub fn new_linear(
        backbone: TabulatedJ2,
        viscosity_pa_s: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        let mut law = Self::new_power_law(backbone, viscosity_pa_s, 1.0, 1.0, duration_s)?;
        law.linear_selection = true;
        Ok(law)
    }

    /// Select `alpha_dot=(s0/eta)*max(f/s0,0)^n`, with positive eta, s0 and dt
    /// and finite n>=1. The existing native kinetics owns parameter admission.
    pub fn new_power_law(
        backbone: TabulatedJ2,
        viscosity_pa_s: f64,
        stress_scale_pa: f64,
        exponent: f64,
        duration_s: f64,
    ) -> Result<Self, MaterialError> {
        let kinetics = PowerLawOverstressJ2::new(
            backbone.backbone().clone(),
            viscosity_pa_s,
            stress_scale_pa,
            exponent,
            duration_s,
        )?;
        let linear_overstress_modulus_pa = viscosity_pa_s / duration_s;
        let (_, mu) = backbone.backbone().elastic.lame();
        if !(3.0 * mu
            + backbone.backbone().kinematic_hardening
            + backbone.hardening().maximum_tangent_pa()
            + linear_overstress_modulus_pa)
            .is_finite()
        {
            return Err(MaterialError::Parameters {
                what: "tabulated hardening and rate return stiffness must be representable".into(),
            });
        }
        let log_kinetic_coefficient = det::ln(stress_scale_pa)
            + (det::ln(viscosity_pa_s) - det::ln(duration_s) - det::ln(stress_scale_pa)) / exponent;
        Ok(Self {
            backbone,
            kinetics,
            linear_selection: false,
            linear_overstress_modulus_pa,
            log_kinetic_coefficient,
        })
    }

    /// Original bounded hardening curve, elasticity and Prager translation.
    pub const fn backbone(&self) -> &TabulatedJ2 {
        &self.backbone
    }
    /// Original viscosity in Pa s.
    pub const fn viscosity_pa_s(&self) -> f64 {
        self.kinetics.viscosity_pa_s()
    }
    /// Actual backward-Euler interval in seconds.
    pub const fn duration_s(&self) -> f64 {
        self.kinetics.duration_s()
    }
    /// Original power exponent, including n=1 for the linear selection.
    pub const fn exponent(&self) -> f64 {
        self.kinetics.exponent()
    }
    /// Original stress normalization; absent for the linear constructor.
    pub const fn stress_scale_pa(&self) -> Option<f64> {
        if self.linear_selection {
            None
        } else {
            Some(self.kinetics.stress_scale_pa())
        }
    }

    fn inverse_rate(&self, delta: f64) -> f64 {
        if delta == 0.0 {
            0.0
        } else if self.exponent() == 1.0 {
            self.linear_overstress_modulus_pa * delta
        } else {
            det::exp(self.log_kinetic_coefficient + det::ln(delta) / self.exponent())
        }
    }

    /// Inverse physical rate law, evaluated without raising stress to a power.
    /// This is overstress alone, not the table's yield resistance.
    pub fn endpoint_overstress_pa(&self, delta: f64) -> Result<f64, TabulatedError> {
        if !delta.is_finite() || delta < 0.0 {
            return Err(TabulatedError::Parameters(
                "rate increment must be finite and nonnegative",
            ));
        }
        let value = self.inverse_rate(delta);
        if !value.is_finite() || delta > 0.0 && value <= 0.0 {
            return Err(TabulatedError::Numerical(
                "endpoint viscous stress is not representable",
            ));
        }
        Ok(value)
    }

    /// Exact scalar restriction for native tensor or free-sided axial returns.
    /// Supply 3 mu or E respectively; the original Prager modulus is added here.
    /// The measured domain is enforced even if an unconstrained root exists.
    pub fn plastic_increment(
        &self,
        trial_pa: f64,
        alpha: f64,
        elastic_correction_pa: f64,
    ) -> Result<f64, TabulatedError> {
        let table = self.backbone.hardening();
        let first = table.segment(alpha)?;
        let radius = table.value_on_segment(alpha, first);
        let correction = elastic_correction_pa + self.backbone.backbone().kinematic_hardening;
        if !trial_pa.is_finite()
            || trial_pa < 0.0
            || !elastic_correction_pa.is_finite()
            || elastic_correction_pa <= 0.0
            || !correction.is_finite()
        {
            return Err(TabulatedError::Parameters(
                "trial stress and elastic return modulus must be finite and nonnegative/positive",
            ));
        }
        if self.exponent() == 1.0 {
            return table.solve_increment(
                trial_pa,
                alpha,
                correction + self.linear_overstress_modulus_pa,
            );
        }
        let excess = trial_pa - radius;
        if excess <= 0.0 {
            return Ok(0.0);
        }
        let last = table.points.len() - 1;
        let maximum = table.maximum_plastic_strain();
        if alpha == maximum {
            return Err(TabulatedError::Domain {
                plastic_strain: alpha + excess / correction,
                maximum,
            });
        }
        // Direct strength increments preserve the partially consumed first
        // segment even when its rise is below the absolute yield's resolution.
        let rise_to_knot = |index: usize| {
            table.slopes[first] * (table.points[first + 1].plastic_strain - alpha)
                + (table.points[index].yield_stress_pa - table.points[first + 1].yield_stress_pa)
        };
        let required = |index: usize| {
            let delta = table.points[index].plastic_strain - alpha;
            correction * delta + rise_to_knot(index) + self.inverse_rate(delta)
        };
        let tolerance = 32.0 * f64::EPSILON * excess;
        let final_residual = required(last) - excess;
        if final_residual < -tolerance {
            return Err(TabulatedError::Domain {
                plastic_strain: alpha + excess / correction,
                maximum,
            });
        }
        // Monotone flow and nondecreasing strength select one measured cell.
        // An overflowing positive required stress remains a valid upper bound.
        let mut low = first + 1;
        let mut high = last;
        while low < high {
            let middle = low + (high - low) / 2;
            if required(middle) < excess {
                low = middle + 1;
            } else {
                high = middle;
            }
        }
        let segment = low - 1;
        let crossed = (table.points[segment].plastic_strain - alpha).max(0.0);
        let rise = if crossed == 0.0 {
            0.0
        } else {
            rise_to_knot(segment)
        };
        let slope = table.slopes[segment];
        let stiffness = correction + slope;
        if !stiffness.is_finite() {
            return Err(TabulatedError::Numerical(
                "combined rate/hardening derivative is not representable",
            ));
        }
        let measured_end = table.points[segment + 1].plastic_strain - alpha;
        // Both omitted terms are nonnegative: elastic-only and kinetic-only
        // roots are upper bounds. Starting at their minimum avoids spending
        // the fixed iteration budget traversing decades of negligible flow.
        let mut end = measured_end.min(excess / correction);
        if end <= 0.0 {
            return Err(TabulatedError::Numerical(
                "positive rate root is below floating-point resolution",
            ));
        }
        let kinetic_log_bound = self.exponent() * (det::ln(excess) - self.log_kinetic_coefficient);
        let smallest = f64::from_bits(1);
        if kinetic_log_bound < det::ln(smallest) {
            return Err(TabulatedError::Numerical(
                "positive rate root is below floating-point resolution",
            ));
        }
        if kinetic_log_bound < det::ln(end) {
            end = end.min(det::exp(kinetic_log_bound).next_up());
        }
        end = end.max(crossed);
        let evaluate = |delta: f64| {
            let viscous = self.inverse_rate(delta);
            let residual = correction * delta + rise + slope * (delta - crossed) + viscous - excess;
            let log_derivative = stiffness * delta + viscous / self.exponent();
            (residual, log_derivative)
        };
        let finish = |delta: f64| {
            let mut bounded = delta.min(end).max(crossed);
            // Subtraction/addition at a table endpoint can round in opposite
            // directions. Stay on the measured side of that rounding boundary.
            if alpha + bounded > table.points[segment + 1].plastic_strain {
                bounded = bounded.next_down();
            }
            if !bounded.is_finite() || bounded <= 0.0 || alpha + bounded == alpha {
                return Err(TabulatedError::Numerical(
                    "rate increment is below retained-history resolution",
                ));
            }
            table.yield_stress_pa(alpha + bounded)?;
            self.endpoint_overstress_pa(bounded)?;
            Ok(bounded)
        };
        if evaluate(end).0.abs() <= tolerance {
            return finish(end);
        }
        if crossed > 0.0 && evaluate(crossed).0.abs() <= tolerance {
            return finish(crossed);
        }
        let lower_delta = crossed.max(smallest);
        if lower_delta >= end || evaluate(lower_delta).0 > tolerance {
            return Err(TabulatedError::Numerical(
                "positive rate root is not representable",
            ));
        }
        let mut lower = det::ln(lower_delta);
        let end_log = det::ln(end);
        let mut upper = end_log;
        let mut position = upper;
        for _ in 0..96 {
            let delta = if position == end_log {
                end
            } else {
                det::exp(position).max(lower_delta).min(end)
            };
            let (residual, derivative) = evaluate(delta);
            if residual.is_finite() && residual.abs() <= tolerance {
                return finish(delta);
            }
            if residual.is_nan() {
                break;
            }
            if residual > 0.0 {
                upper = position;
            } else {
                lower = position;
            }
            let middle = 0.5 * lower + 0.5 * upper;
            if middle == lower || middle == upper {
                break;
            }
            let newton = position - residual / derivative;
            position = if derivative.is_finite()
                && derivative > 0.0
                && newton.is_finite()
                && newton > lower
                && newton < upper
            {
                newton
            } else {
                middle
            };
        }
        Err(TabulatedError::Numerical(
            "bounded tabulated power-law return did not converge",
        ))
    }

    /// Inspect accepted native history without advancing it or imposing an
    /// instantaneous yield surface. Viscous endpoints can retain overstress.
    pub fn validate_history(&self, strain: &Voigt, state: &J2State) -> Result<(), EnergyError> {
        self.backbone
            .backbone()
            .validate_overstress_history(strain, state)?;
        self.backbone
            .hardening()
            .yield_stress_pa(state.alpha)
            .map_err(|_| {
                EnergyError::Invalid("viscous history leaves the measured hardening domain")
            })?;
        Ok(())
    }

    fn returned(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<(J2Plasticity, f64, Update), TabulatedError> {
        self.validate_history(strain, state)
            .map_err(|_| TabulatedError::State("invalid native tabulated rate history"))?;
        let mut shifted = self.backbone.backbone().clone();
        shifted.yield_stress = self.backbone.hardening().yield_stress_pa(state.alpha)?;
        let trial = shifted.radial_return(strain, state, 0.0);
        if !trial.plastic {
            return Ok((shifted, 0.0, trial));
        }
        let (_, mu) = shifted.elastic.lame();
        let delta = self.plastic_increment(trial.q_trial, state.alpha, 3.0 * mu)?;
        if delta == 0.0 {
            // A roundoff excess with no representable hardening increment.
            return Ok((shifted, 0.0, trial));
        }
        let hardening = self
            .backbone
            .hardening()
            .strength_increment_pa(state.alpha, delta)?;
        let secant = hardening / delta + self.endpoint_overstress_pa(delta)? / delta;
        let secant = self.backbone.bounded_return_secant(
            trial.q_trial - shifted.yield_stress,
            state.alpha,
            secant,
        )?;
        let returned = shifted.radial_return(strain, state, secant);
        self.backbone
            .hardening()
            .yield_stress_pa(state.alpha + returned.dgamma)?;
        if returned.stress.iter().any(|x| !x.is_finite()) {
            return Err(TabulatedError::Numerical(
                "tabulated viscous stress is not representable",
            ));
        }
        Ok((shifted, secant, returned))
    }

    fn tangent_of(&self, state: &J2State, returned: &Update) -> Result<Tangent6, TabulatedError> {
        let delta = returned.dgamma;
        let slope = self.backbone.hardening().tangent_pa(state.alpha + delta)?;
        let viscous = if delta > 0.0 {
            (self.endpoint_overstress_pa(delta)? / delta) / self.exponent()
        } else {
            0.0
        };
        let tangent = self
            .backbone
            .backbone()
            .tangent_from_return(returned, slope + viscous);
        if tangent.iter().flatten().any(|x| !x.is_finite()) {
            return Err(TabulatedError::Numerical(
                "tabulated rate tangent is not representable",
            ));
        }
        Ok(tangent)
    }

    /// Checked response on the original tensor flow, tangent and J2 history.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<TabulatedJ2Response, TabulatedError> {
        let (shifted, secant, returned) = self.returned(strain, state)?;
        let next = shifted.update_state_with_overstress(strain, state, secant);
        self.validate_history(strain, &next)
            .map_err(|_| TabulatedError::State("returned rate history is inadmissible"))?;
        Ok(TabulatedJ2Response {
            stress_pa: returned.stress,
            tangent_pa: self.tangent_of(state, &returned)?,
            state: next,
            plastic_increment: returned.dgamma,
        })
    }

    /// Fixed-history storage; viscosity contributes no recoverable energy.
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<StoredEnergy, EnergyError> {
        self.backbone.stored_energy(strain, state)
    }

    /// Stress and the discrete reduced potential at frozen incoming history.
    /// Its tensor-Voigt gradient has twice the shear stress in shear slots.
    /// Physical viscous heat is f_v*d, distinct from the kinetic potential below.
    pub fn stress_and_potential(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<(Voigt, f64), TabulatedError> {
        let (shifted, secant, returned) = self.returned(strain, state)?;
        let next = shifted.update_state_with_overstress(strain, state, secant);
        let stored = self
            .stored_energy(strain, &next)
            .map_err(|_| TabulatedError::State("rate potential history is inadmissible"))?;
        let delta = returned.dgamma;
        let viscous = self.endpoint_overstress_pa(delta)?;
        let fraction = 1.0 / (1.0 + 1.0 / self.exponent());
        let potential = stored.total_j_m3()
            + self.backbone.backbone().yield_stress * delta
            + fraction * viscous * delta;
        if !potential.is_finite() {
            return Err(TabulatedError::Numerical(
                "tabulated rate potential is not representable",
            ));
        }
        Ok((returned.stress, potential))
    }

    /// Original endpoint-work audit with exact table storage and separately
    /// retained physical viscosity, yield heat and numerical integration loss.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.validate_history(old_strain, old_state)?;
        let (shifted, secant, _) = self.returned(strain, old_state).map_err(|_| {
            EnergyError::Invalid("tabulated rate return leaves its domain or is not representable")
        })?;
        let mut next = shifted
            .advance_with_overstress_energy(old_strain, old_state, strain, tolerances, secant)?;
        self.validate_history(strain, &next.state)?;
        let delta = next.state.alpha - old_state.alpha;
        let viscous = self
            .endpoint_overstress_pa(delta)
            .map_err(|_| EnergyError::Invalid("tabulated endpoint rate is not representable"))?;
        let residual = self
            .backbone
            .yield_function_with_state(&next.stress_pa, &next.state)
            - viscous;
        if !residual.is_finite() {
            return Err(EnergyError::Invalid(
                "tabulated rate consistency is not representable",
            ));
        }
        if residual > tolerances.yield_pa || delta > 0.0 && residual.abs() > tolerances.yield_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            });
        }
        let table = self.backbone.hardening();
        next.energy.initial.hardening_j_m3 = table
            .hardening_energy_j_m3(old_state.alpha)
            .map_err(|_| EnergyError::Invalid("invalid initial table storage"))?;
        next.energy.final_state.hardening_j_m3 = table
            .hardening_energy_j_m3(next.state.alpha)
            .map_err(|_| EnergyError::Invalid("invalid final table storage"))?;
        next.energy.physical_dissipation_j_m3 =
            (self.backbone.backbone().yield_stress + viscous) * delta;
        next.energy.numerical_dissipation_j_m3 += table
            .numerical_loss_j_m3(old_state.alpha, delta)
            .map_err(|_| EnergyError::Invalid("invalid table integration loss"))?;
        next.energy.residual_j_m3 = next.energy.mechanical_work_j_m3
            - next.energy.stored_change_j_m3()
            - next.energy.physical_dissipation_j_m3
            - next.energy.numerical_dissipation_j_m3;
        if !next.energy.initial.total_j_m3().is_finite()
            || !next.energy.final_state.total_j_m3().is_finite()
            || !next.energy.physical_dissipation_j_m3.is_finite()
            || !next.energy.numerical_dissipation_j_m3.is_finite()
            || !next.energy.residual_j_m3.is_finite()
        {
            return Err(EnergyError::Invalid(
                "tabulated rate work is not representable",
            ));
        }
        if next.energy.residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(EnergyError::EnergyResidual {
                residual_j_m3: next.energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(next)
    }
}

impl SmallStrainLaw for TabulatedOverstressJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State {
        self.backbone.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        self.returned(strain, state)
            .map_or([f64::NAN; 6], |(_, _, value)| value.stress)
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        let Ok((_, _, returned)) = self.returned(strain, state) else {
            return [[f64::NAN; 6]; 6];
        };
        self.tangent_of(state, &returned)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        let Ok((shifted, secant, _)) = self.returned(strain, state) else {
            return J2State {
                plastic_strain: [f64::NAN; 6],
                alpha: f64::NAN,
            };
        };
        shifted.update_state_with_overstress(strain, state, secant)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "bounded tabulated J2 with associative Perzyna flow and original Prager backstress",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.j2-tabulated-overstress", "0.1.0", Ambition::Solid,
            vec!["small-strain Cauchy J2".into(),
                "bounded tabulated energetic isotropic/Prager hardening".into(),
                "associative Perzyna n>=1 at explicit physical duration".into()],
            self.backbone.card().validity,
            vec!["no extrapolation, sub-yield creep, finite strain, temperature evolution or empirical heat partition".into()],
            0.05,
        )
    }
}
