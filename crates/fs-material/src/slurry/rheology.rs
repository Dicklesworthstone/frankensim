//! Source-specific, steady Herschel–Bulkley rheology and laminar pipe hydraulics.
//!
//! The signed scalar law is `|tau| = tau_y + tau_ref (|gamma|/gamma_ref)^n`
//! when flowing. At rest it is set-valued: `|tau| <= tau_y`. Bingham behavior
//! is the `n=1` case; Newtonian behavior additionally requires `tau_y=0`.
//! Referencing an excess stress at a shear rate avoids fractional SI dimensions
//! for the usual consistency coefficient K. All parameters remain caller supplied.
//!
//! Formula references: Hackley & Ferraris, NIST SP 946 (2001), §6, p.10,
//! <https://www.nist.gov/system/files/documents/2017/07/24/b01014.pdf>;
//! Ashena et al. (2023), §§Fluid velocity / Flow rate, equations 1–4 and 26,
//! <https://doi.org/10.1007/s13202-022-01566-2>.
//! Constant parameters apply only inside the declared source domain. No mineral
//! name, density, or solids fraction supplies rheology. There is no settling,
//! thixotropy, wall slip, turbulence, viscoelasticity, or LBM simulation claim.

use fs_matdb::Provenance;
use fs_math::det;
use fs_qty::{Dimensionless, Frequency, Qty, Stress, Temperature};

#[path = "rheology/pipe.rs"]
mod pipe;
pub use pipe::{LaminarCircularPipe, PipeFlow, PressureSolution};

/// Local irreversible heating in watts per cubic metre.
pub type DissipationDensity = Qty<-1, 1, -3, 0, 0, 0>;

/// Invalid input, unsupported source state, or bounded inverse refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SlurryError {
    /// Invalid declaration or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Actual query outside the original inclusive support, in coherent SI.
    OutsideValidity {
        /// Physical coordinate whose support was exceeded.
        axis: &'static str,
        /// Requested value in coherent SI, or a dimensionless fraction.
        value: f64,
        /// Inclusive lower limit in the same units as `value`.
        minimum: f64,
        /// Inclusive upper limit in the same units as `value`.
        maximum: f64,
    },
    /// Requested flow tolerance was not achieved within the iteration budget.
    NoConvergence {
        /// Number of bisection iterations attempted before refusal.
        iterations: u32,
    },
}

impl core::fmt::Display for SlurryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "slurry: {self:?}")
    }
}
impl std::error::Error for SlurryError {}

fn require(ok: bool, message: &'static str) -> Result<(), SlurryError> {
    if ok {
        Ok(())
    } else {
        Err(SlurryError::Invalid(message))
    }
}
fn positive(x: f64, message: &'static str) -> Result<f64, SlurryError> {
    require(x.is_finite() && x > 0.0, message)?;
    Ok(x)
}
fn admit(axis: &'static str, x: f64, bounds: [f64; 2]) -> Result<(), SlurryError> {
    require(x.is_finite(), "query must be finite")?;
    if x < bounds[0] || x > bounds[1] {
        return Err(SlurryError::OutsideValidity {
            axis,
            value: x,
            minimum: bounds[0],
            maximum: bounds[1],
        });
    }
    Ok(())
}

/// Actual mixture state; solids fraction is dry-solids mass / total wet mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlurryConditions {
    /// Absolute temperature in kelvin.
    pub temperature: Temperature,
    /// Dimensionless mass fraction, distinct from a volume fraction.
    pub solids_mass_fraction: Dimensionless,
}

impl SlurryConditions {
    /// Transfer the stream's dry-solid mass fraction with an explicit actual
    /// temperature. This does not infer rheology from density, equilibrate inlet
    /// temperatures, or assert that a fit applies to the stream's constituents.
    /// Zero-flow and entirely dry streams cannot supply slurry conditions.
    pub fn from_stream(
        stream: &super::SlurryStream,
        temperature: Temperature,
    ) -> Result<Self, SlurryError> {
        positive(
            temperature.value(),
            "temperature must be absolute and positive",
        )?;
        let fraction = stream
            .state()
            .solids_mass_fraction
            .ok_or(SlurryError::Invalid("zero-flow stream has no composition"))?;
        require(
            fraction.value() < 1.0,
            "slurry rheology requires a liquid carrier",
        )?;
        Ok(Self {
            temperature,
            solids_mass_fraction: fraction.quantity(),
        })
    }
}

/// Caller-declared applicability of this one constant-parameter fit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlurryDomain {
    shear_rate: [Frequency; 2],
    temperature: [Temperature; 2],
    solids_mass_fraction: [Dimensionless; 2],
}

impl SlurryDomain {
    /// Inclusive shear-rate magnitudes, temperature, and wet-basis mass fraction.
    /// Declaring support down to zero is an explicit model assertion, not a
    /// claim that a rheometer measured every arbitrarily small positive rate.
    pub fn new(
        shear_rate: [Frequency; 2],
        temperature: [Temperature; 2],
        solids_mass_fraction: [Dimensionless; 2],
    ) -> Result<Self, SlurryError> {
        for bounds in [
            shear_rate.map(|v| v.value()),
            temperature.map(|v| v.value()),
            solids_mass_fraction.map(|v| v.value()),
        ] {
            require(
                bounds.iter().all(|v| v.is_finite()) && bounds[0] <= bounds[1],
                "source bounds must be finite and ordered",
            )?;
        }
        require(
            shear_rate[0].value() >= 0.0 && shear_rate[1].value() > 0.0,
            "shear support must be nonnegative with positive upper bound",
        )?;
        require(
            temperature[0].value() > 0.0,
            "temperature must be absolute and positive",
        )?;
        require(
            solids_mass_fraction[0].value() >= 0.0 && solids_mass_fraction[1].value() < 1.0,
            "slurry mass fractions require 0 <= solids < 1",
        )?;
        Ok(Self {
            shear_rate,
            temperature,
            solids_mass_fraction,
        })
    }

    /// Inclusive supported magnitudes in inverse seconds.
    #[must_use]
    pub const fn shear_rate_bounds(self) -> [Frequency; 2] {
        self.shear_rate
    }

    fn admit_conditions(self, state: SlurryConditions) -> Result<(), SlurryError> {
        admit(
            "temperature",
            state.temperature.value(),
            self.temperature.map(|v| v.value()),
        )?;
        admit(
            "solids mass fraction",
            state.solids_mass_fraction.value(),
            self.solids_mass_fraction.map(|v| v.value()),
        )
    }
    fn admit_rate(self, magnitude: f64) -> Result<(), SlurryError> {
        admit(
            "shear-rate magnitude",
            magnitude,
            self.shear_rate.map(|v| v.value()),
        )
    }
}

/// Original sample, preparation/rheometry method, and source usage terms.
#[derive(Clone, Debug, PartialEq)]
pub struct SlurryTestContext {
    /// Original sample identity, or an explicitly labelled authored fixture.
    pub specimen: String,
    /// Rheometry procedure, sample preparation and parameter-fitting method.
    pub method: String,
    /// Original source, usage terms and optional source artifact hash.
    pub provenance: Provenance,
}

/// A dimensionally explicit parameterization of the Herschel–Bulkley law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HerschelBulkleyParameters {
    /// Nonnegative limiting stress magnitude for an unsheared material, in pascals.
    pub yield_stress: Stress,
    /// Excess above yield at `reference_shear_rate`, in pascals.
    pub reference_excess_stress: Stress,
    /// Positive supported shear-rate magnitude in inverse seconds.
    pub reference_shear_rate: Frequency,
    /// Positive flow index; 1 is Bingham, below 1 is shear thinning.
    pub flow_index: Dimensionless,
}

/// Stress is not uniquely determined by a zero shear rate in a yield-stress fluid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ShearResponse {
    /// Zero shear rate, with no dissipative heating and a set-valued stress.
    Unyielded {
        /// Inclusive signed static stress interval `[-tau_y, tau_y]`, in pascals.
        stress_bounds: [Stress; 2],
    },
    /// Nonzero supported shear rate with uniquely determined signed stress.
    Yielded {
        /// Stress in pascals, with the same sign as the imposed shear rate.
        stress: Stress,
        /// Positive irreversible heating `tau * gamma`, in watts per cubic metre.
        dissipation: DissipationDensity,
    },
}

/// Immutable source-qualified scalar constitutive law, not a calibrated plant model.
#[derive(Clone, Debug, PartialEq)]
pub struct HerschelBulkleyCard {
    parameters: HerschelBulkleyParameters,
    domain: SlurryDomain,
    test: SlurryTestContext,
}

impl HerschelBulkleyCard {
    /// Admit a finite constitutive fit, supported reference rate and complete source context.
    pub fn new(
        parameters: HerschelBulkleyParameters,
        domain: SlurryDomain,
        test: SlurryTestContext,
    ) -> Result<Self, SlurryError> {
        require(
            parameters.yield_stress.value().is_finite() && parameters.yield_stress.value() >= 0.0,
            "yield stress must be finite and nonnegative",
        )?;
        positive(
            parameters.reference_excess_stress.value(),
            "reference excess stress must be positive",
        )?;
        positive(
            parameters.reference_shear_rate.value(),
            "reference shear rate must be positive",
        )?;
        let n = positive(parameters.flow_index.value(), "flow index must be positive")?;
        positive(1.0 / n, "reciprocal flow index is unrepresentable")?;
        domain.admit_rate(parameters.reference_shear_rate.value())?;
        for text in [
            &test.specimen,
            &test.method,
            &test.provenance.source,
            &test.provenance.license,
        ] {
            require(
                !text.trim().is_empty(),
                "rheology requires specimen, method, provenance and usage terms",
            )?;
        }
        Ok(Self {
            parameters,
            domain,
            test,
        })
    }

    /// Original yield, reference stress/rate and flow-index declaration.
    #[must_use]
    pub const fn parameters(&self) -> HerschelBulkleyParameters {
        self.parameters
    }
    /// Original applicability domain; construction does not establish calibration.
    #[must_use]
    pub const fn domain(&self) -> SlurryDomain {
        self.domain
    }
    /// Original specimen, method and provenance retained by this fit.
    #[must_use]
    pub fn test_context(&self) -> &SlurryTestContext {
        &self.test
    }

    /// Signed stress and nonnegative dissipation, or the complete static yield interval.
    pub fn stress_at_shear_rate(
        &self,
        rate: Frequency,
        state: SlurryConditions,
    ) -> Result<ShearResponse, SlurryError> {
        self.domain.admit_conditions(state)?;
        self.domain.admit_rate(rate.value().abs())?;
        if rate.value() == 0.0 {
            let y = self.parameters.yield_stress.value();
            return Ok(ShearResponse::Unyielded {
                stress_bounds: [Stress::new(-y), Stress::new(y)],
            });
        }
        let magnitude = self.stress_magnitude(rate.value().abs())?;
        let dissipation = positive(
            magnitude * rate.value().abs(),
            "dissipation is unrepresentable",
        )?;
        Ok(ShearResponse::Yielded {
            stress: Stress::new(magnitude.copysign(rate.value())),
            dissipation: DissipationDensity::new(dissipation),
        })
    }

    /// Invert the original signed law; stresses in the yield interval give zero rate.
    pub fn shear_rate_at_stress(
        &self,
        stress: Stress,
        state: SlurryConditions,
    ) -> Result<Frequency, SlurryError> {
        self.domain.admit_conditions(state)?;
        require(stress.value().is_finite(), "stress must be finite")?;
        let s = stress.value().abs();
        let y = self.parameters.yield_stress.value();
        let rate = if s <= y {
            0.0
        } else {
            let lower = self.domain.shear_rate[0].value();
            let upper = self.domain.shear_rate[1].value();
            let lower_stress = if lower == 0.0 {
                y
            } else {
                self.stress_magnitude(lower)?
            };
            let upper_stress = self.stress_magnitude(upper)?;
            admit("stress magnitude", s, [lower_stress, upper_stress])?;
            if s == upper_stress {
                upper
            } else if s == lower_stress {
                lower
            } else {
                positive(
                    self.parameters.reference_shear_rate.value()
                        * det::pow(
                            (s - y) / self.parameters.reference_excess_stress.value(),
                            1.0 / self.parameters.flow_index.value(),
                        ),
                    "shear rate is unrepresentable",
                )?
            }
        };
        self.domain.admit_rate(rate)?;
        Ok(Frequency::new(rate.copysign(stress.value())))
    }

    fn stress_magnitude(&self, rate: f64) -> Result<f64, SlurryError> {
        let p = self.parameters;
        let excess = positive(
            p.reference_excess_stress.value()
                * det::pow(rate / p.reference_shear_rate.value(), p.flow_index.value()),
            "excess stress is unrepresentable",
        )?;
        let stress = positive(p.yield_stress.value() + excess, "stress is unrepresentable")?;
        require(
            stress > p.yield_stress.value(),
            "flowing stress is indistinguishable from yield",
        )?;
        Ok(stress)
    }
}

#[cfg(test)]
#[path = "rheology/tests.rs"]
mod tests;
