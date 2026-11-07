//! Source-bound electrical heating on the original lumped enthalpy owner.
//!
//! The prepared Ohmic material owns resistance and Joule power; fs-conduction
//! owns sensible/latent storage, convection, radiation and the implicit solve.
//! An accepted temperature can drive fs-render thermal emission independently.
//! Observing that light must not debit radiative cooling a second time.

use fs_conduction::ConductionError;
use fs_conduction::lumped::LumpedThermalEnvironment;
use fs_conduction::lumped::endpoint::{DrivenLumpedEnthalpyStep, LumpedEnthalpyEndpoint};
use fs_material::conductor::graph::MaterialOhmicNode;
use fs_material::conductor::{ConductorError, OhmicDrive, ResolvedConductor};
use fs_material::phase::UniformEnthalpyStepInput;
use fs_thermal::enthalpy::ThermalLaw;

/// One held ideal drive on an existing, continuously supported material law.
///
/// Preparation bounds positive electrical feedback once over the complete
/// source interval. Reuse this object for successive thermal steps; changing
/// the drive prepares a new object. There is no fitted hot resistance,
/// constant heat-capacity substitute, circuit storage or material-name preset.
/// Mass and geometry are prescribed specimen inputs; density-dependent free
/// expansion and spatial current/temperature gradients need their own owners.
pub struct PreparedIncandescentSource<'a> {
    source: &'a MaterialOhmicNode,
    drive: OhmicDrive,
    maximum_positive_power_slope_w_k: f64,
}

/// Original endpoint electrical state, including its complete source receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct IncandescentElectricalState {
    pub conductor: ResolvedConductor,
    pub current_a: f64,
    pub voltage_v: f64,
    pub joule_power_w: f64,
    /// Independently checked terminal work rate minus native Joule power [W].
    pub electrical_power_residual_w: f64,
}

/// One proposed physical step; neither original input is mutated on failure.
#[derive(Clone, Debug, PartialEq)]
pub struct IncandescentStep {
    /// Original caloric state and integrated storage/boundary/source balance.
    pub thermal: DrivenLumpedEnthalpyStep,
    /// The SAME final temperature used by the implicit thermal source.
    pub electrical: IncandescentElectricalState,
}

/// Original thermal and electrical refusals stay distinguishable.
#[derive(Debug)]
pub enum IncandescentError {
    Thermal(ConductionError),
    Electrical(ConductorError),
    Invalid(&'static str),
}

impl core::fmt::Display for IncandescentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "incandescent thermal owner: {error}"),
            Self::Electrical(error) => write!(f, "incandescent electrical owner: {error}"),
            Self::Invalid(what) => write!(f, "incandescent coupling: {what}"),
        }
    }
}
impl std::error::Error for IncandescentError {}
impl From<ConductionError> for IncandescentError {
    fn from(error: ConductionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<ConductorError> for IncandescentError {
    fn from(error: ConductorError) -> Self {
        Self::Electrical(error)
    }
}

impl<'a> PreparedIncandescentSource<'a> {
    /// Retain a fixed-geometry, temperature-dependent source and signed drive.
    ///
    /// Current and voltage control retain their distinct feedback. Electrical
    /// phase-fraction dependence and free expansion require additional tangent
    /// bounds and are refused here; native caloric latent plateaus remain valid
    /// when the electrical source has no independent phase coordinate.
    pub fn try_new(
        source: &'a MaterialOhmicNode,
        drive: OhmicDrive,
    ) -> Result<Self, IncandescentError> {
        let maximum_positive_power_slope_w_k =
            source.maximum_positive_joule_temperature_slope_w_k(drive)?;
        Ok(Self {
            source,
            drive,
            maximum_positive_power_slope_w_k,
        })
    }

    /// Original prepared source; its selected claims and conditions are intact.
    #[must_use]
    pub const fn source(&self) -> &'a MaterialOhmicNode {
        self.source
    }

    /// Held signed terminal drive, with no hidden hot-resistance correction.
    #[must_use]
    pub const fn drive(&self) -> OhmicDrive {
        self.drive
    }

    /// Sufficient interval for a single-valued implicit thermal/source step.
    ///
    /// Infinity means no positive-feedback restriction, not unlimited temporal
    /// accuracy. Choose a smaller duration when resolving turn-on/off dynamics
    /// and check temporal convergence independently of energy closure.
    pub fn maximum_step_s(
        &self,
        environment: &LumpedThermalEnvironment,
        input: UniformEnthalpyStepInput<'_>,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<f64, IncandescentError> {
        checkpoint()?;
        let (endpoint, feedback) = self.prepare(environment, input)?;
        Ok(endpoint.maximum_step_s_for_internal_power(feedback, checkpoint)?)
    }

    /// Advance one declared duration with resistance evaluated at the unknown
    /// final temperature. This reuses the original bounded enthalpy bisection.
    ///
    /// `input.internal_heat_j` is an ADDITIONAL prescribed source, separate from
    /// terminal electrical work. The returned thermal report counts Joule heat
    /// once and reports signed convection/radiation into the body. The native
    /// Biot gate, source support, enthalpy tolerance and joule residual budget
    /// all remain active. Overlong positive-feedback steps refuse explicitly;
    /// use `maximum_step_s` when choosing substeps. Cancellation publishes no
    /// partial state and never alters the caller's accepted state.
    pub fn advance(
        &self,
        environment: &LumpedThermalEnvironment,
        input: UniformEnthalpyStepInput<'_>,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<IncandescentStep, IncandescentError> {
        checkpoint()?;
        let (endpoint, feedback) = self.prepare(environment, input)?;
        let mut electrical_error = None;
        let result = endpoint.solve_with_internal_power(
            feedback,
            |_, temperature| match self.electrical_state(temperature) {
                Ok(sample) => Ok(sample.joule_power_w),
                Err(error) => {
                    let what = error.to_string();
                    electrical_error = Some(error);
                    Err(ConductionError::ScenarioRow {
                        region: "incandescent source".to_owned(),
                        what,
                        fix: "retain an admitted original electrical source and drive".to_owned(),
                    })
                }
            },
            &mut checkpoint,
        );
        if let Some(error) = electrical_error {
            return Err(error.into());
        }
        let thermal = result?;
        let electrical = self.electrical_state(thermal.state.temperature_k())?;
        if electrical.joule_power_w * input.duration_s != thermal.additional_internal_heat_j {
            return Err(IncandescentError::Invalid(
                "accepted electrical work must equal the native thermal deposition",
            ));
        }
        checkpoint()?;
        Ok(IncandescentStep {
            thermal,
            electrical,
        })
    }

    fn prepare<'body>(
        &self,
        environment: &LumpedThermalEnvironment,
        input: UniformEnthalpyStepInput<'body>,
    ) -> Result<(LumpedEnthalpyEndpoint<'body>, f64), IncandescentError> {
        let volume = self.source.length_m() * self.source.area_m2();
        if !volume.is_finite()
            || volume <= 0.0
            || !input.volume_m3.is_finite()
            || input.volume_m3 <= 0.0
            || (volume / input.volume_m3 - 1.0).abs() > 64.0 * f64::EPSILON
            || input.initial.material_card_identity()
                != self.source.source().original().material().card_identity()
        {
            return Err(IncandescentError::Invalid(
                "electrical geometry must span the whole body on the same original material card",
            ));
        }
        let knots = input.curve.knots();
        let [low, high] = self.source.temperature_range_k();
        if knots[0].temperature_k < low || knots[knots.len() - 1].temperature_k > high {
            return Err(IncandescentError::Invalid(
                "electrical source support must cover the complete caloric chart",
            ));
        }
        let material = self.source.source().original().material();
        // A card may contain measurements at many pressures or compositions.
        // Sharing its hash alone does not bind those physical conditions.
        // Authored charts carry no receipts and retain their declared authority.
        for receipt in input
            .curve
            .source_receipts()
            .iter()
            .chain(environment.transport.conductivity_receipts())
        {
            for (axis, value) in &receipt.query_point {
                if receipt.axis_quantities.get(axis) != material.axis_quantities().get(axis)
                    || (axis != self.source.temperature_axis()
                        && material.state_coordinate(axis) != Some(*value))
                {
                    return Err(IncandescentError::Invalid(
                        "caloric and electrical source conditions or quantity conventions differ",
                    ));
                }
            }
        }
        let endpoint = environment.prepare_endpoint(input)?;
        // The material thermal owner bounds its ORIGINAL Cp law, including
        // nonlinear sensible storage; latent plateaus contribute zero dT/dh.
        let temperature_slope = input.curve.maximum_temperature_slope();
        let slope = self.maximum_positive_power_slope_w_k;
        let feedback = if slope == 0.0 || temperature_slope == 0.0 {
            0.0
        } else {
            let direct = slope * (temperature_slope / input.mass_kg);
            let feedback = if direct.is_finite() && direct > 0.0 {
                direct
            } else {
                (slope / input.mass_kg) * temperature_slope
            };
            feedback.next_up()
        };
        if !feedback.is_finite() || feedback < 0.0 {
            return Err(IncandescentError::Invalid(
                "positive electrothermal feedback is unrepresentable",
            ));
        }
        Ok((endpoint, feedback))
    }

    fn electrical_state(
        &self,
        temperature_k: f64,
    ) -> Result<IncandescentElectricalState, ConductorError> {
        let conductor = self.source.resolve_at_temperature_k(temperature_k)?;
        let current_a = conductor.current_a_for_drive(self.drive)?;
        let joule_power_w = conductor.joule_power_w(current_a)?;
        let voltage_v = match self.drive {
            OhmicDrive::Voltage(value) => value,
            OhmicDrive::Current(_) => current_a * conductor.resistance_ohm(),
        };
        let delivered = voltage_v * current_a;
        let electrical_power_residual_w = delivered - joule_power_w;
        if !voltage_v.is_finite()
            || !delivered.is_finite()
            || !electrical_power_residual_w.is_finite()
            || (current_a != 0.0 && (voltage_v == 0.0 || delivered == 0.0))
            || electrical_power_residual_w.abs() > 64.0 * f64::EPSILON * joule_power_w
        {
            return Err(ConductorError::InvalidDerived {
                quantity: "incandescent terminal work and native Joule power",
            });
        }
        Ok(IncandescentElectricalState {
            conductor,
            current_a,
            voltage_v,
            joule_power_w,
            electrical_power_residual_w,
        })
    }
}
