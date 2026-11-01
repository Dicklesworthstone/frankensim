//! Simultaneous finite Norton creep and energetic isotropic J2 plasticity.
//!
//! One total isochoric inelastic metric is transported by the original finite
//! exponential return. `alpha` retains equivalent PLASTIC accumulation, while
//! `creep_alpha` independently retains equivalent CREEP accumulation. Their sum
//! bounds metric reachability; only plastic alpha stores `H*alpha²/2` and
//! expands the yield radius. This isotropic model needs no separate orientation
//! memories for the two coaxial endpoint flows; it does not claim anisotropic,
//! kinematic-hardening, damage or general noncoaxial mechanism factorization.
//!
//! Backward Euler jointly satisfies `q+3mu*(dc+dp)=q_trial`,
//! `dc=dt*(s0/eta)*(q/s0)^n`, and plastic complementarity against
//! `q <= Y+H*(alpha_old+dp)`. Creep remains active below yield. At active yield,
//! the scalar problem reduces to the same native monotone rate solve with
//! effective stiffness `3mu*H/(3mu+H)`. No sequential mechanism splitting occurs.
//!
//! Physical losses are separately `Y*dp` and `q*dc`. The reduced potential is
//! `psi+Y*dp+n/(n+1)*q*dc`. Its full derivative is the returned Piola stress;
//! the common native tensor tangent differentiates the actual joint scalar
//! return. All stresses are Kirchhoff and energies use original reference
//! volume. Parameter calibration, primary/tertiary creep and rupture are not
//! supplied by this explicitly bounded constitutive family.

use super::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use super::{FiniteJ2, FiniteJ2Error, FiniteJ2Response, FiniteJ2State, RadialReturn, ReturnMap};
use fs_math::det;

/// One immutable combined material and one positive physical interval.
/// This has no `FiniteJ2Material` adapter that could discard its creep kinetics.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteNortonJ2 {
    native: FiniteOverstressJ2,
}

impl FiniteNortonJ2 {
    pub fn new(
        backbone: FiniteJ2,
        kinetics: FiniteOverstressParameters,
        duration_s: f64,
        maximum_accumulated_creep_strain: f64,
    ) -> Result<Self, FiniteJ2Error> {
        let backbone = backbone.with_creep_history_limit(maximum_accumulated_creep_strain)?;
        Ok(Self {
            native: FiniteOverstressJ2::new(backbone, kinetics, duration_s)?,
        })
    }

    /// Elastic/plastic energy and combined-history inspection owner.
    pub fn backbone(&self) -> &FiniteJ2 {
        self.native.backbone()
    }
    /// Explicit Norton viscosity, stress scale and stress exponent.
    pub fn parameters(&self) -> &FiniteOverstressParameters {
        self.native.parameters()
    }
    pub fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }
    pub fn maximum_accumulated_creep_strain(&self) -> f64 {
        self.backbone()
            .maximum_creep_alpha
            .expect("combined constructor retains its bound")
    }
    pub fn initial_state(&self) -> FiniteJ2State {
        self.backbone().initial_state()
    }
    pub fn validate_history(&self, state: &FiniteJ2State) -> Result<(), FiniteJ2Error> {
        self.backbone().validate_history(state)
    }
    pub fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        self.backbone().stored_energy(f, state)
    }
    /// Accepted endpoints obey the rate-independent PLASTIC yield surface;
    /// creep history remains retained and is never re-integrated by inspection.
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        self.backbone()
            .validate_committed(f, state, yield_tolerance_pa)
    }
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<FiniteJ2Response, FiniteJ2Error> {
        let returned = self.return_map(f, old)?;
        self.backbone().response_from_return(f, old, returned)
    }
    /// Same coupled return without full tangent assembly during line search.
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        let returned = self.return_map(f, old)?;
        Ok((returned.first_piola, returned.potential))
    }

    fn return_map(&self, f: &[f64; 9], old: &FiniteJ2State) -> Result<ReturnMap, FiniteJ2Error> {
        let backbone = self.backbone();
        let trial = backbone.elastic_endpoint(f, old)?;
        let stress = trial
            .deviator
            .map(|x| 2.0 * backbone.parameters.shear_modulus_pa * x);
        let q_trial = det::sqrt(1.5) * super::norm(&stress);
        if !q_trial.is_finite() {
            return Err(FiniteJ2Error::Numerical(
                "combined trial stress is not finite",
            ));
        }
        let radial = self.radial_return(q_trial, old.alpha)?;
        backbone.finish_return(f, old, trial, q_trial, radial)
    }

    fn radial_return(&self, trial: f64, alpha: f64) -> Result<RadialReturn, FiniteJ2Error> {
        let p = self.backbone().parameters();
        let h = 3.0 * p.shear_modulus_pa;
        let hardening = p.hardening_modulus_pa;
        let yield_old = self.backbone().yield_radius(alpha);
        let exponent = self.parameters().exponent;
        let potential_fraction = 1.0 / (1.0 + 1.0 / exponent);
        if trial == 0.0 {
            let scale = if exponent == 1.0 {
                let modulus = self.parameters().viscosity_pa_s / self.duration_s();
                modulus / (h + modulus)
            } else {
                1.0
            };
            return Ok(RadialReturn {
                plastic_increment: 0.0,
                creep_increment: 0.0,
                scale,
                tangent_correction: 0.0,
                plastic_dissipation: 0.0,
                creep_dissipation: 0.0,
                primary_potential: 0.0,
                rate_potential: 0.0,
            });
        }
        // Monotonicity identifies the active set with one log-rate evaluation.
        // A floating-point tie takes the inactive side; no increment is clipped.
        let inactive = trial <= yield_old || {
            let log_h = det::ln(h);
            let log_creep = self.native.log_increment_at_stress(yield_old);
            let log_flow = log_h + log_creep;
            let log_gap = det::ln(trial - yield_old);
            let tie = 8.0 * f64::EPSILON * (log_h.abs() + log_creep.abs() + log_gap.abs());
            log_flow >= log_gap || log_flow.is_finite() && log_gap - log_flow <= tie
        };
        let (q, creep, plastic, correction) = if inactive {
            let flow = self.native.solve(trial, h)?;
            let scale = flow.overstress_pa / trial;
            let correction =
                -(h / (h + flow.tangent_modulus_pa)) * ((1.0 - 1.0 / exponent) * scale);
            (flow.overstress_pa, flow.increment, 0.0, correction)
        } else if hardening == 0.0 {
            let creep = self.native.increment_at_stress(yield_old)?;
            let plastic = ((trial - yield_old) - h * creep) / h;
            (yield_old, creep, plastic, -yield_old / trial)
        } else {
            let denominator = h + hardening;
            let weight = hardening / denominator;
            let effective = h * weight;
            if weight <= 0.0 || effective <= 0.0 {
                return Err(FiniteJ2Error::Numerical(
                    "combined effective hardening is not representable",
                ));
            }
            let effective_trial = yield_old + weight * (trial - yield_old);
            let flow = self.native.solve(effective_trial, effective)?;
            let creep = flow.increment;
            let plastic = ((trial - yield_old) - h * creep) / denominator;
            let q = yield_old + hardening * plastic;
            let v = flow.tangent_modulus_pa;
            // dq/dqtrial - q/qtrial, as a sum of nonpositive terms. This
            // avoids subtracting similar radial moduli near linear creep.
            let correction = -(effective / (effective + v)) * ((1.0 - 1.0 / exponent) * q / trial)
                - (h / denominator) * (v / (effective + v)) * (yield_old / trial);
            (q, creep, plastic, correction)
        };
        if !q.is_finite()
            || q <= 0.0
            || !creep.is_finite()
            || creep <= 0.0
            || !plastic.is_finite()
            || plastic < 0.0
        {
            return Err(FiniteJ2Error::Numerical(
                "joint creep/plastic increments must be representable and nonnegative",
            ));
        }
        // Independently audit the actual kinetic, stress and complementarity
        // equations before constructing any new metric or retained state.
        let kinetic_q = self.native.endpoint_overstress_pa(creep)?;
        let balance = h.mul_add(creep + plastic, q) - trial;
        let yield_residual = q - (yield_old + hardening * plastic);
        let scale = trial.max(yield_old).max(q);
        let tolerance = 512.0 * f64::EPSILON * scale;
        if (kinetic_q - q).abs() > 512.0 * f64::EPSILON * q
            || !balance.is_finite()
            || balance.abs() > tolerance
            || (if plastic > 0.0 {
                yield_residual.abs()
            } else {
                yield_residual.max(0.0)
            }) > tolerance
        {
            return Err(FiniteJ2Error::Numerical(
                "joint finite creep/plastic return violates its physical equations",
            ));
        }
        let plastic_dissipation = p.yield_stress_pa * plastic;
        let creep_dissipation = q * creep;
        Ok(RadialReturn {
            plastic_increment: plastic,
            creep_increment: creep,
            scale: q / trial,
            tangent_correction: correction,
            plastic_dissipation,
            creep_dissipation,
            primary_potential: plastic_dissipation,
            rate_potential: potential_fraction * creep_dissipation,
        })
    }
}
