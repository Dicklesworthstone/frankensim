//! Calibrated SI shear-thinning viscosity shared by material-driven flow owners.
//!
//! Carreau-Yasuda: mu(g) = mu_inf + (mu_0-mu_inf)
//! * [1 + (lambda*g)^a]^((n-1)/a), with g = sqrt(2 D:D) in 1/s.
//! See Yasuda, Armstrong & Cohen (1981), doi:10.1007/BF01513059, and
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.polymer/polymer_ug_fluidflow.06.06.html>.
//!
//! Temperature dependence explicitly selects mu(T,g)=a_T*mu_ref(g) or
//! mu(T,g)=a_T*mu_ref(a_T*g). These are different calibrations; neither is
//! inferred from a material name. The original WLF/Arrhenius declaration in
//! TemperatureShift supplies a_T and its derivative, without a second fit.
//! See ANSYS Polyflow, "Temperature Dependence of the Viscosity", equations
//! 2-159/160: <https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug_workspaces/flu_ug_matpro_physics_equations.html>.
//!
//! This is generalized Newtonian shear thinning, not elastic memory, normal
//! stress differences, a yield surface, freezing, or an extensional-flow fit.
//! No calibration data, pressure dependence or experimental authority is added.

use super::{FluidTransport, FluidTransportState};
use crate::MaterialError;
use crate::visco::temperature::{TemperatureShift, TemperatureShiftState};
use fs_math::det;

#[cfg(test)]
mod tests;

fn parameters(what: &'static str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn state_error(what: &'static str) -> MaterialError {
    MaterialError::State { what: what.into() }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// Checked reference-temperature viscosity curve, with a finite calibrated
/// shear-rate interval. Viscosities are dynamic [Pa s], never kinematic.
/// n is restricted to (0,1]; a may be any positive finite transition exponent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauYasuda {
    viscosity_zero_pa_s: f64,
    viscosity_infinite_pa_s: f64,
    time_s: f64,
    power_index: f64,
    transition: f64,
    rate_support_per_s: [f64; 2],
    log_time_s: f64,
    log_excess_pa_s: f64,
    log_thinning: f64,
}

struct CurveSample {
    viscosity: f64,
    scaled_derivative: f64,
    stress_tangent: f64,
    // None means an exactly zero derivative. A nonfinite logarithm instead
    // means a nonrepresentable derivative, including the singular origin.
    log_derivative_magnitude: Option<f64>,
}

impl CarreauYasuda {
    /// Admit a declared reference curve and inclusive shear-rate support.
    /// a=2 recovers Carreau; n=1, lambda=0 or equal plateaus is Newtonian.
    /// mu_inf=0 is allowed, but every evaluated positive viscosity and
    /// differential stress must remain representable. No cutoff is invented.
    ///
    /// # Errors
    /// Refuses invalid parameters, reversed plateaus or invalid support.
    pub fn new(
        viscosity_zero_pa_s: f64,
        viscosity_infinite_pa_s: f64,
        time_s: f64,
        power_index: f64,
        transition: f64,
        rate_support_per_s: [f64; 2],
    ) -> Result<Self, MaterialError> {
        if !positive(viscosity_zero_pa_s)
            || !viscosity_infinite_pa_s.is_finite()
            || viscosity_infinite_pa_s < 0.0
            || viscosity_infinite_pa_s > viscosity_zero_pa_s
            || !time_s.is_finite()
            || time_s < 0.0
            || !positive(power_index)
            || power_index > 1.0
            || !positive(transition)
            || !rate_support_per_s
                .into_iter()
                .all(|r| r.is_finite() && r >= 0.0)
            || rate_support_per_s[0] > rate_support_per_s[1]
        {
            return Err(parameters(
                "invalid Carreau-Yasuda parameters or reference shear-rate support",
            ));
        }
        let excess = viscosity_zero_pa_s - viscosity_infinite_pa_s;
        Ok(Self {
            viscosity_zero_pa_s,
            viscosity_infinite_pa_s,
            time_s,
            power_index,
            transition,
            rate_support_per_s,
            log_time_s: if time_s > 0.0 { det::ln(time_s) } else { 0.0 },
            log_excess_pa_s: if excess > 0.0 { det::ln(excess) } else { 0.0 },
            log_thinning: if power_index < 1.0 {
                det::ln(1.0 - power_index)
            } else {
                0.0
            },
        })
    }

    /// Original (mu_0 [Pa s], mu_inf [Pa s], lambda [s], n, a) declaration.
    pub const fn parameters(&self) -> [f64; 5] {
        [
            self.viscosity_zero_pa_s,
            self.viscosity_infinite_pa_s,
            self.time_s,
            self.power_index,
            self.transition,
        ]
    }

    /// Inclusive support of the reference curve, before temperature shifting.
    pub const fn rate_support_per_s(&self) -> [f64; 2] {
        self.rate_support_per_s
    }

    /// Evaluate the reference-temperature dynamic viscosity [Pa s].
    ///
    /// # Errors
    /// Refuses unsupported rates or unrepresentable positive responses.
    pub fn viscosity_pa_s(&self, rate_per_s: f64) -> Result<f64, MaterialError> {
        Ok(self.sample(rate_per_s)?.viscosity)
    }

    fn sample(&self, rate: f64) -> Result<CurveSample, MaterialError> {
        if !rate.is_finite()
            || !(self.rate_support_per_s[0]..=self.rate_support_per_s[1]).contains(&rate)
        {
            return Err(state_error(
                "shear rate leaves the calibrated reference curve",
            ));
        }
        let delta = self.viscosity_zero_pa_s - self.viscosity_infinite_pa_s;
        let constant = self.time_s == 0.0 || self.power_index == 1.0 || delta == 0.0;
        if constant || rate == 0.0 {
            let derivative = if constant || self.transition > 1.0 {
                None
            } else if self.transition == 1.0 {
                Some(self.log_excess_pa_s + self.log_time_s + self.log_thinning)
            } else {
                Some(f64::INFINITY)
            };
            return Ok(CurveSample {
                viscosity: self.viscosity_zero_pa_s,
                scaled_derivative: 0.0,
                stress_tangent: self.viscosity_zero_pa_s,
                log_derivative_magnitude: derivative,
            });
        }

        let x = self.time_s * rate;
        let log_x = if x.is_normal() {
            det::ln(x)
        } else {
            // lambda*rate may overflow/underflow while the actual viscosity
            // remains representable. The declared input itself is unchanged.
            self.log_time_s + det::ln(rate)
        };
        let z = self.transition * log_x;
        let (softened, fraction, complement, log_fraction) = if log_x >= 0.0 {
            let tail = det::exp(-z);
            let correction = log1p_small_positive(tail);
            (
                log_x + correction / self.transition,
                1.0 / (1.0 + tail),
                tail / (1.0 + tail),
                -correction,
            )
        } else {
            let tail = det::exp(z);
            let correction = log1p_small_positive(tail);
            (
                correction / self.transition,
                tail / (1.0 + tail),
                1.0 / (1.0 + tail),
                z - correction,
            )
        };
        // softplus(a*ln(x))/a is formed without infinity/infinity, even
        // when the transition exponent or lambda*rate is extreme.
        let log_shape = (self.power_index - 1.0) * softened;
        let shape = det::exp(log_shape);
        let excess = if shape > 0.0 {
            delta * shape
        } else {
            // An underflowed shape can still carry a finite viscosity when
            // multiplied by a large declared zero-shear plateau.
            det::exp(self.log_excess_pa_s + log_shape)
        };
        let viscosity = if log_shape > -core::f64::consts::LN_2 {
            self.viscosity_zero_pa_s - delta * -det::expm1(log_shape)
        } else {
            self.viscosity_infinite_pa_s + excess
        };
        let scaled_derivative = -(1.0 - self.power_index) * excess * fraction;
        // Do not obtain this by subtracting two nearly equal viscosities:
        // its strictly positive n*excess limit matters for strong thinning.
        let stress_tangent =
            self.viscosity_infinite_pa_s + excess * (complement + self.power_index * fraction);
        if !positive(viscosity) || !positive(stress_tangent) || !scaled_derivative.is_finite() {
            return Err(state_error(
                "Carreau-Yasuda viscosity or differential stress is unrepresentable",
            ));
        }
        Ok(CurveSample {
            viscosity,
            scaled_derivative,
            stress_tangent,
            log_derivative_magnitude: Some(
                self.log_excess_pa_s + log_shape + self.log_thinning + log_fraction - det::ln(rate),
            ),
        })
    }
}

// Argument is exp(-|z|) in [0,1]. The short series retains a small departure
// from the plateau that forming ln(1+x) directly would round away.
fn log1p_small_positive(x: f64) -> f64 {
    if x < 1e-4 {
        x * (1.0 + x * (-0.5 + x * (1.0 / 3.0 + x * (-0.25 + x / 5.0))))
    } else {
        det::ln(1.0 + x)
    }
}

/// The two distinct thermal calibrations use the SAME declared a_T.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ViscosityShiftMode {
    /// mu(T,g) = a_T * mu_ref(g).
    Vertical,
    /// mu(T,g) = a_T * mu_ref(a_T*g).
    VerticalAndRate,
}

/// One checked reference curve and the original WLF/Arrhenius calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalCarreauYasuda {
    curve: CarreauYasuda,
    shift: TemperatureShift,
    mode: ViscosityShiftMode,
}

impl ThermalCarreauYasuda {
    /// Combine admitted declarations without refitting or extending support.
    pub const fn new(
        curve: CarreauYasuda,
        shift: TemperatureShift,
        mode: ViscosityShiftMode,
    ) -> Self {
        Self { curve, shift, mode }
    }

    pub const fn curve(&self) -> &CarreauYasuda {
        &self.curve
    }
    pub const fn temperature_shift(&self) -> &TemperatureShift {
        &self.shift
    }
    pub const fn mode(&self) -> ViscosityShiftMode {
        self.mode
    }

    /// Resolve the actual physical rate and temperature against both supports.
    /// In VerticalAndRate, reference support applies to a_T*g, so the physical
    /// admitted shear interval moves with temperature. No numerical clipping,
    /// rescaling of the authored interval, or viscosity floor is applied.
    ///
    /// # Errors
    /// Refuses unsupported temperature/rate and unrepresentable constitutive
    /// values. A singular scalar zero-rate derivative alone is reported as
    /// None, since viscosity and the stress tangent are still well defined.
    pub fn evaluate(
        &self,
        temperature_k: f64,
        rate_per_s: f64,
    ) -> Result<ThermalRheologyState, MaterialError> {
        if !rate_per_s.is_finite() || rate_per_s < 0.0 {
            return Err(state_error(
                "physical shear rate must be finite and nonnegative",
            ));
        }
        let temperature =
            self.shift
                .evaluate(temperature_k)
                .map_err(|error| MaterialError::State {
                    what: format!("rheology temperature shift: {error}"),
                })?;
        let reduced_rate = match self.mode {
            ViscosityShiftMode::Vertical => rate_per_s,
            ViscosityShiftMode::VerticalAndRate => temperature.factor() * rate_per_s,
        };
        if !reduced_rate.is_finite() || rate_per_s > 0.0 && reduced_rate == 0.0 {
            return Err(state_error(
                "reduced reference shear rate is unrepresentable",
            ));
        }
        let reference = self.curve.sample(reduced_rate)?;
        let derivative = reference.log_derivative_magnitude.map(|log| {
            log + temperature.log_factor()
                * match self.mode {
                    ViscosityShiftMode::Vertical => 1.0,
                    ViscosityShiftMode::VerticalAndRate => 2.0,
                }
        });
        let viscosity = temperature.factor() * reference.viscosity;
        let scaled_derivative = match derivative {
            Some(log) if rate_per_s > 0.0 => {
                let direct = temperature.factor() * reference.scaled_derivative;
                if direct != 0.0 && direct.is_finite() {
                    direct
                } else {
                    // A vanishing reference slope can become representable
                    // after shifting. Restore the actual g*dmu/dg in log space.
                    -det::exp(log + det::ln(rate_per_s))
                }
            }
            _ => 0.0,
        };
        let stress_tangent = temperature.factor() * reference.stress_tangent;
        let temperature_derivative = match self.mode {
            ViscosityShiftMode::Vertical => viscosity,
            ViscosityShiftMode::VerticalAndRate => stress_tangent,
        } * temperature.log_derivative_per_k();
        if !positive(viscosity)
            || !positive(stress_tangent)
            || !scaled_derivative.is_finite()
            || !temperature_derivative.is_finite()
            || temperature_derivative >= 0.0
        {
            return Err(state_error(
                "shifted viscosity or constitutive tangent is unrepresentable",
            ));
        }
        let shear_derivative = match derivative {
            None => Some(0.0),
            Some(log) => {
                let magnitude = det::exp(log);
                positive(magnitude).then_some(-magnitude)
            }
        };
        Ok(ThermalRheologyState {
            law: *self,
            temperature,
            rate_per_s,
            reduced_rate_per_s: reduced_rate,
            viscosity_pa_s: viscosity,
            shear_derivative_pa_s2: shear_derivative,
            scaled_shear_derivative_pa_s: scaled_derivative,
            stress_tangent_pa_s: stress_tangent,
            temperature_derivative_pa_s_per_k: temperature_derivative,
        })
    }
}

/// Original declaration and complete response at one actual constitutive query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalRheologyState {
    law: ThermalCarreauYasuda,
    temperature: TemperatureShiftState,
    rate_per_s: f64,
    reduced_rate_per_s: f64,
    viscosity_pa_s: f64,
    shear_derivative_pa_s2: Option<f64>,
    scaled_shear_derivative_pa_s: f64,
    stress_tangent_pa_s: f64,
    temperature_derivative_pa_s_per_k: f64,
}

impl ThermalRheologyState {
    pub const fn law(&self) -> &ThermalCarreauYasuda {
        &self.law
    }
    pub const fn temperature(&self) -> &TemperatureShiftState {
        &self.temperature
    }
    pub const fn shear_rate_per_s(&self) -> f64 {
        self.rate_per_s
    }
    pub const fn reduced_rate_per_s(&self) -> f64 {
        self.reduced_rate_per_s
    }
    pub const fn dynamic_viscosity_pa_s(&self) -> f64 {
        self.viscosity_pa_s
    }
    /// dmu/dg [Pa s^2]. None identifies an unbounded zero-rate right derivative
    /// (a<1), overflow, or nonzero underflow; it is never a fabricated zero.
    pub const fn shear_derivative_pa_s2(&self) -> Option<f64> {
        self.shear_derivative_pa_s2
    }
    /// g*dmu/dg [Pa s]; exactly zero at rest even for a<1.
    pub const fn scaled_shear_derivative_pa_s(&self) -> f64 {
        self.scaled_shear_derivative_pa_s
    }
    /// d(mu*g)/dg [Pa s], evaluated positively without cancellation.
    pub const fn stress_tangent_pa_s(&self) -> f64 {
        self.stress_tangent_pa_s
    }
    /// Analytic dmu/dT at fixed physical shear rate [Pa s/K].
    pub const fn temperature_derivative_pa_s_per_k(&self) -> f64 {
        self.temperature_derivative_pa_s_per_k
    }

    /// Replace viscosity once in an original current-temperature transport
    /// response, retaining its density, heat capacity and conductivity.
    ///
    /// The source must already describe this SAME temperature. FluidTransport
    /// does not expose a temperature coordinate, so matching that state is the
    /// caller's responsibility. Its original viscosity is retained for audit
    /// but is not multiplied into the declared Carreau-Yasuda viscosity.
    /// Owned sources and borrowed sources (via FluidTransport for &T) both work.
    ///
    /// # Errors
    /// Preserves the original material refusal, or refuses an unrepresentable
    /// coherent transport combination. No source or accepted solver state mutates.
    pub fn bind<F: FluidTransport>(
        &self,
        source: F,
    ) -> Result<RheologicalTransport<F>, MaterialError> {
        let original = source.transport_state()?;
        let transport = FluidTransportState::try_new(
            original.density_kg_m3(),
            original.specific_heat_j_kg_k(),
            original.thermal_conductivity_w_m_k(),
            self.viscosity_pa_s,
        )?;
        Ok(RheologicalTransport {
            source,
            original,
            rheology: *self,
            transport,
        })
    }
}

/// Exact source response plus the explicitly selected constitutive replacement.
/// The source keeps its own receipts; this model does not mint source authority.
#[derive(Debug)]
pub struct RheologicalTransport<F> {
    source: F,
    original: FluidTransportState,
    rheology: ThermalRheologyState,
    transport: FluidTransportState,
}

impl<F> RheologicalTransport<F> {
    pub const fn source(&self) -> &F {
        &self.source
    }
    pub const fn original_transport(&self) -> &FluidTransportState {
        &self.original
    }
    pub const fn rheology(&self) -> &ThermalRheologyState {
        &self.rheology
    }
}

impl<F> FluidTransport for RheologicalTransport<F> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}
