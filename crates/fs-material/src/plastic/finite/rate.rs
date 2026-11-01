//! Explicit-time Perzyna flow on the original multiplicative finite J2 owner.
//!
//! The declared Kirchhoff overstress is `f=q-sigma_y-H*alpha`, with
//! `alpha_dot=(s0/eta)*positive(f/s0)^n`, `eta` in Pa s, `s0` in Pa and
//! `n>=1`. This normalization reduces exactly to `alpha_dot=f/eta` at n=1.
//! Backward Euler solves `(3mu+H)*dgamma+f(dgamma)=positive(f_trial)`.
//! The scalar kinetics changes neither stored hardening nor the native
//! exponential metric update. The consistent tangent differentiates this
//! actual implicit flow, including `df/dgamma=f/(n*dgamma)`.
//!
//! Physical dissipation is `(sigma_y+f)*dgamma`. The reduced incremental
//! potential instead adds `sigma_y*dgamma+n/(n+1)*f*dgamma` to stored energy.
//! Its gradient is Piola stress; replacing that integral with the full viscous
//! work breaks the Newton potential. Endpoint quadrature error remains separate.
//! A held overstressed deformation relaxes over physical time; a committed
//! viscous state need not lie inside the rate-independent yield surface.
//!
//! This is a threshold law, not sub-yield creep, damage or measured material
//! calibration. All original finite strain/history bounds still apply. The
//! [associative Perzyna definition](https://neml.readthedocs.io/en/latest/vp_flow/perzyna.html)
//! provides the flow-family reference; its alternative fluidity normalization
//! is not silently substituted for the explicit Pa-s convention used here.

use super::{FiniteJ2, FiniteJ2Error, FiniteJ2Response, FiniteJ2State};
use fs_math::det;

/// Physical kinetics at one fixed temperature; there are no inferred defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteOverstressParameters {
    /// Positive physical viscosity in Pa s at this fixed temperature.
    pub viscosity_pa_s: f64,
    /// Positive Kirchhoff stress scale in the declared power-law normalization.
    pub stress_scale_pa: f64,
    /// Dimensionless Perzyna exponent, finite and at least one.
    pub exponent: f64,
}

impl FiniteOverstressParameters {
    /// Admit the coefficients independently of any chosen physical duration.
    pub fn validate(&self) -> Result<(), FiniteJ2Error> {
        if !self.viscosity_pa_s.is_finite()
            || self.viscosity_pa_s <= 0.0
            || !self.stress_scale_pa.is_finite()
            || self.stress_scale_pa <= 0.0
            || !self.exponent.is_finite()
            || self.exponent < 1.0
        {
            return Err(FiniteJ2Error::Parameters(
                "finite Perzyna flow requires positive finite eta and stress scale, and finite exponent >= 1",
            ));
        }
        Ok(())
    }
}

/// Native finite backbone plus physical kinetics and one immutable interval.
/// This deliberately does not implement `FiniteJ2Material`: exposing only the
/// backbone to a static solver would silently discard the selected rate law.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteOverstressJ2 {
    backbone: FiniteJ2,
    parameters: FiniteOverstressParameters,
    duration_s: f64,
    log_viscosity_over_duration: f64,
}

pub(super) struct FlowReturn {
    pub(super) increment: f64,
    pub(super) overstress_pa: f64,
    pub(super) tangent_modulus_pa: f64,
    pub(super) transverse_overstress_pa: f64,
    pub(super) potential_fraction: f64,
}

impl FiniteOverstressJ2 {
    /// Admit a native backbone, explicit kinetics and one physical interval.
    pub fn new(
        backbone: FiniteJ2,
        parameters: FiniteOverstressParameters,
        duration_s: f64,
    ) -> Result<Self, FiniteJ2Error> {
        parameters.validate()?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FiniteJ2Error::Parameters(
                "finite Perzyna flow requires a positive finite physical duration",
            ));
        }
        if parameters.exponent == 1.0 {
            let modulus = parameters.viscosity_pa_s / duration_s;
            let p = backbone.parameters();
            if !modulus.is_finite()
                || modulus <= 0.0
                || !(3.0 * p.shear_modulus_pa + p.hardening_modulus_pa + modulus).is_finite()
            {
                return Err(FiniteJ2Error::Parameters(
                    "linear finite overstress modulus and return denominator must be representable",
                ));
            }
        }
        let log_viscosity_over_duration = det::ln(parameters.viscosity_pa_s) - det::ln(duration_s);
        Ok(Self {
            backbone,
            parameters,
            duration_s,
            log_viscosity_over_duration,
        })
    }

    /// Native elastic, energetic-hardening and applicability owner.
    pub const fn backbone(&self) -> &FiniteJ2 {
        &self.backbone
    }
    /// Explicit, fixed-temperature physical kinetics.
    pub const fn parameters(&self) -> &FiniteOverstressParameters {
        &self.parameters
    }
    /// Physical duration used by this backward-Euler return.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Unstrained, virgin native history.
    pub fn initial_state(&self) -> FiniteJ2State {
        self.backbone.initial_state()
    }
    /// Check the same complete metric/accumulation bounds as the native owner.
    pub fn validate_history(&self, state: &FiniteJ2State) -> Result<(), FiniteJ2Error> {
        self.backbone.validate_history(state)
    }
    /// Inspect stored energy at fixed history, without an additional flow step.
    pub fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        self.backbone.stored_energy(f, state)
    }

    /// Inspect a valid accepted viscous state without performing another flow
    /// step or projecting its legitimate overstress onto equilibrium yield.
    /// The tolerance retains the common endpoint API and must be finite and
    /// nonnegative; the discrete flow residual is checked in `response` itself.
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        if !yield_tolerance_pa.is_finite() || yield_tolerance_pa < 0.0 {
            return Err(FiniteJ2Error::Parameters(
                "yield tolerance must be finite and nonnegative",
            ));
        }
        let endpoint = self.backbone.elastic_endpoint(f, state)?;
        self.backbone
            .check_elastic_domain(endpoint.log_j, &endpoint.deviator)?;
        self.backbone
            .energy(endpoint.log_j, &endpoint.deviator, state.alpha)?;
        let p = self.backbone.parameters();
        let mut tau = endpoint
            .deviator
            .map(|value| 2.0 * p.shear_modulus_pa * value);
        for i in [0, 4, 8] {
            tau[i] += p.bulk_modulus_pa * endpoint.log_j;
        }
        let inverse_transpose = super::transpose(&super::inverse(f, endpoint.jacobian)?);
        let piola = super::multiply(&tau, &inverse_transpose);
        if tau.iter().chain(&piola).any(|x| !x.is_finite()) {
            return Err(FiniteJ2Error::Numerical(
                "accepted viscous stress is not representable",
            ));
        }
        Ok(())
    }

    /// Native implicit return and full analytic derivative of its Piola stress.
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<FiniteJ2Response, FiniteJ2Error> {
        let returned = self.backbone.return_map_with_rate(f, old, Some(self))?;
        self.backbone.response_from_return(f, old, returned)
    }

    /// Same return as `response`, with no 9-by-9 tangent assembly in line search.
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        let returned = self.backbone.return_map_with_rate(f, old, Some(self))?;
        Ok((returned.first_piola, returned.potential))
    }

    /// Evaluate the declared kinetic relation at a nonnegative flow increment.
    pub fn endpoint_overstress_pa(&self, increment: f64) -> Result<f64, FiniteJ2Error> {
        if !increment.is_finite() || increment < 0.0 {
            return Err(FiniteJ2Error::Parameters(
                "finite nonnegative flow increment required",
            ));
        }
        if increment == 0.0 {
            return Ok(0.0);
        }
        let p = &self.parameters;
        let stress = if p.exponent == 1.0 {
            (p.viscosity_pa_s / self.duration_s) * increment
        } else {
            let log_scale = det::ln(p.stress_scale_pa);
            det::exp(
                log_scale
                    + (self.log_viscosity_over_duration + det::ln(increment) - log_scale)
                        / p.exponent,
            )
        };
        if !stress.is_finite() || stress <= 0.0 {
            return Err(FiniteJ2Error::Numerical(
                "positive endpoint viscous stress is not representable",
            ));
        }
        Ok(stress)
    }

    pub(super) fn solve(&self, excess: f64, a: f64) -> Result<FlowReturn, FiniteJ2Error> {
        if self.parameters.exponent == 1.0 {
            return self.from_secant(self.parameters.viscosity_pa_s / self.duration_s, excess, a);
        }
        // Solve for x=ln(f/dgamma), rather than forming overflowing powers.
        // r(x)=x/n-(1-1/n)softplus(ln(A)-x)-rhs is strictly increasing.
        // Newton stays in a verified bracket; bounded bisection is the fallback.
        let inverse_n = 1.0 / self.parameters.exponent;
        let complement = 1.0 - inverse_n;
        let log_a = det::ln(a);
        let rhs = self.log_viscosity_over_duration * inverse_n
            + complement * (det::ln(self.parameters.stress_scale_pa) - det::ln(excess));
        let residual = |x: f64| x * inverse_n - complement * softplus(log_a - x) - rhs;
        let mut lower = det::ln(f64::from_bits(1));
        let mut upper = det::ln(f64::MAX);
        if !rhs.is_finite() || residual(lower) > 0.0 || residual(upper) < 0.0 {
            return Err(FiniteJ2Error::Numerical(
                "finite Perzyna secant root is outside representable bounds",
            ));
        }
        let mut x = self.log_viscosity_over_duration.clamp(lower, upper);
        for _ in 0..96 {
            let r = residual(x);
            if !r.is_finite() {
                break;
            }
            if r.abs() <= 32.0 * f64::EPSILON * (1.0 + x.abs() + rhs.abs()) {
                if let Ok(flow) = self.from_secant(det::exp(x), excess, a) {
                    return Ok(flow);
                }
            }
            if r < 0.0 {
                lower = x;
            } else {
                upper = x;
            }
            let difference = log_a - x;
            let sigmoid = if difference >= 0.0 {
                1.0 / (1.0 + det::exp(-difference))
            } else {
                let e = det::exp(difference);
                e / (1.0 + e)
            };
            let derivative = inverse_n + complement * sigmoid;
            let newton = x - r / derivative;
            let middle = 0.5 * lower + 0.5 * upper;
            if middle == lower || middle == upper {
                return self.from_secant(det::exp(middle), excess, a);
            }
            x = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                middle
            };
        }
        Err(FiniteJ2Error::Numerical(
            "bounded finite Perzyna return did not converge",
        ))
    }

    fn from_secant(&self, modulus: f64, excess: f64, a: f64) -> Result<FlowReturn, FiniteJ2Error> {
        let denominator = a + modulus;
        if !modulus.is_finite() || modulus <= 0.0 || !denominator.is_finite() {
            return Err(FiniteJ2Error::Numerical(
                "finite Perzyna return denominator is not representable",
            ));
        }
        let increment = excess / denominator;
        let overstress_pa = self.endpoint_overstress_pa(increment)?;
        let residual = a.mul_add(increment, overstress_pa) - excess;
        if increment <= 0.0
            || !residual.is_finite()
            || residual.abs() > 128.0 * f64::EPSILON * excess
        {
            return Err(FiniteJ2Error::Numerical(
                "finite Perzyna endpoint violates its physical flow equation",
            ));
        }
        let inverse_n = 1.0 / self.parameters.exponent;
        let tangent_modulus_pa = (overstress_pa / increment) * inverse_n;
        if !tangent_modulus_pa.is_finite() || !(a + tangent_modulus_pa).is_finite() {
            return Err(FiniteJ2Error::Numerical(
                "finite Perzyna derivative is not representable",
            ));
        }
        Ok(FlowReturn {
            increment,
            overstress_pa,
            tangent_modulus_pa,
            transverse_overstress_pa: (1.0 - inverse_n) * overstress_pa,
            potential_fraction: 1.0 / (1.0 + inverse_n),
        })
    }
}

fn softplus(x: f64) -> f64 {
    if x > 0.0 {
        x + det::ln(1.0 + det::exp(-x))
    } else {
        det::ln(1.0 + det::exp(x))
    }
}
