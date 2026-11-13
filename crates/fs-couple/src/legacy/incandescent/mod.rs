//! Source-bound electrical heating on the original lumped enthalpy owner.
//!
//! The prepared Ohmic material owns resistance and Joule power; fs-conduction
//! owns sensible/latent storage, convection, radiation and the implicit solve.
//! An accepted temperature can drive fs-render thermal emission independently.
//! Observing that light must not debit radiative cooling a second time.

/// Series electrical assembly on original spatial thermal cells.
pub mod filament;
/// Held-current and held-voltage histories on the original physical clock.
pub mod history;

use fs_blake3::ContentHash;
use fs_conduction::ConductionError;
use fs_conduction::lumped::endpoint::{
    DrivenLumpedEnthalpyStep, LumpedEndpointEnvironment, LumpedEnthalpyEndpoint,
};
use fs_conduction::lumped::{LumpedEnthalpyBody, LumpedEnthalpyTransport};
use fs_conduction::radiation::MaterialRadiationNodeState;
use fs_material::conductor::graph::MaterialOhmicNode;
use fs_material::conductor::{ConductorError, OhmicDrive, ResolvedConductor};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, PhaseStateError, UniformEnthalpyStepInput,
};
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
    feedback: PreparedFeedback,
}

#[derive(Clone, Copy)]
enum PreparedFeedback {
    Temperature(f64),
    SpecificEnthalpy {
        curve: ContentHash,
        maximum_positive_power_slope: f64,
    },
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
    /// The SAME final temperature and phase used by the implicit heat source.
    pub electrical: IncandescentElectricalState,
    /// Actual accepted material/emissivity receipt for a phase transport body.
    /// Rendering can bind this state through `MaterialThermalBinding` without
    /// selecting another material or subtracting another radiative heat loss.
    /// Temperature-table/declared transport has no fresh material query here.
    pub thermal_surface: Option<MaterialRadiationNodeState>,
}

/// Original thermal and electrical refusals stay distinguishable.
#[derive(Debug)]
pub enum IncandescentError {
    Thermal(ConductionError),
    Electrical(ConductorError),
    Caloric(PhaseStateError),
    Invalid(&'static str),
}

impl core::fmt::Display for IncandescentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "incandescent thermal owner: {error}"),
            Self::Electrical(error) => write!(f, "incandescent electrical owner: {error}"),
            Self::Caloric(error) => write!(f, "incandescent caloric owner: {error}"),
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
impl From<PhaseStateError> for IncandescentError {
    fn from(error: PhaseStateError) -> Self {
        Self::Caloric(error)
    }
}

impl<'a> PreparedIncandescentSource<'a> {
    /// Retain a fixed-geometry, temperature-dependent source and signed drive.
    ///
    /// Current and voltage control retain their distinct feedback. Electrical
    /// phase dependence uses [`Self::try_new_with_caloric_curve`] instead;
    /// native latent storage remains valid here when the electrical source has
    /// no independent phase coordinate. Free expansion requires its own owner.
    pub fn try_new(
        source: &'a MaterialOhmicNode,
        drive: OhmicDrive,
    ) -> Result<Self, IncandescentError> {
        let maximum_positive_power_slope_w_k =
            source.maximum_positive_joule_temperature_slope_w_k(drive)?;
        Ok(Self {
            source,
            drive,
            feedback: PreparedFeedback::Temperature(maximum_positive_power_slope_w_k),
        })
    }

    /// Bind electrical feedback to the original temperature/phase trajectory.
    ///
    /// Resistance may change as liquid mass fraction evolves at constant
    /// melting temperature. Preparation follows every original caloric span
    /// through the pinned electrical source, retaining source knots, surface
    /// crossings and the original sensible Cp law. It caches a positive dP/dh
    /// bound, including latent feedback, for this exact caloric chart and drive.
    /// Successive steps reuse that bound; a different chart must be prepared
    /// separately. Temperature-only electrical sources are also supported.
    pub fn try_new_with_caloric_curve(
        source: &'a MaterialOhmicNode,
        drive: OhmicDrive,
        curve: &EquilibriumEnthalpyPhaseCurve,
    ) -> Result<Self, IncandescentError> {
        let maximum_positive_power_slope =
            source.maximum_positive_joule_specific_enthalpy_slope(drive, curve)?;
        Ok(Self {
            source,
            drive,
            feedback: PreparedFeedback::SpecificEnthalpy {
                curve: curve.identity(),
                maximum_positive_power_slope,
            },
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
        environment: &impl LumpedEndpointEnvironment,
        input: UniformEnthalpyStepInput<'_>,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<f64, IncandescentError> {
        checkpoint()?;
        let endpoint = environment.prepare_endpoint(input)?;
        let feedback = self.feedback_for_body(endpoint.body())?;
        Ok(endpoint.maximum_step_s_for_internal_power(feedback, checkpoint)?)
    }

    /// Advance one declared duration with resistance evaluated at the unknown
    /// final temperature and liquid fraction. This reuses the original bounded
    /// enthalpy bisection, including isothermal melting and solidification.
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
        environment: &impl LumpedEndpointEnvironment,
        input: UniformEnthalpyStepInput<'_>,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<IncandescentStep, IncandescentError> {
        checkpoint()?;
        let endpoint = environment.prepare_endpoint(input)?;
        self.advance_prepared_endpoint(&endpoint, checkpoint)
    }

    /// Drive an already admitted native endpoint, including material phase
    /// conductivity and emissivity. Radiation stays a boundary heat term and
    /// uses the same final enthalpy as electrical work and thermal storage.
    pub fn advance_prepared_endpoint(
        &self,
        endpoint: &LumpedEnthalpyEndpoint<'_>,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<IncandescentStep, IncandescentError> {
        checkpoint()?;
        let feedback = self.feedback_for_body(endpoint.body())?;
        let curve = endpoint.body().phase_curve();
        let mut source_error = None;
        let result = endpoint.solve_with_internal_power(
            feedback,
            |enthalpy, temperature| match self.trial_electrical_state(curve, enthalpy, temperature)
            {
                Ok(sample) => Ok(sample.joule_power_w),
                Err(error) => {
                    let what = error.to_string();
                    source_error = Some(error);
                    Err(ConductionError::ScenarioRow {
                        region: "incandescent source".to_owned(),
                        what,
                        fix: "retain an admitted original electrical source and drive".to_owned(),
                    })
                }
            },
            &mut checkpoint,
        );
        if let Some(error) = source_error {
            return Err(error);
        }
        let thermal = result?;
        let electrical = self.electrical_state(
            thermal.state.temperature_k(),
            thermal.state.liquid_mass_fraction(),
        )?;
        if electrical.joule_power_w * endpoint.duration_s() != thermal.additional_internal_heat_j {
            return Err(IncandescentError::Invalid(
                "accepted electrical work must equal the native thermal deposition",
            ));
        }
        let thermal_surface = match endpoint.body().transport() {
            LumpedEnthalpyTransport::MaterialPhase(transport) => {
                let (_, sample) = transport
                    .radiation()
                    .evaluate(thermal.state.specific_enthalpy_j_kg())?;
                if sample.phase != thermal.state {
                    return Err(IncandescentError::Invalid(
                        "accepted emissivity and stored enthalpy must describe the same physical state",
                    ));
                }
                Some(sample)
            }
            LumpedEnthalpyTransport::Temperature(_) => None,
        };
        checkpoint()?;
        Ok(IncandescentStep {
            thermal,
            electrical,
            thermal_surface,
        })
    }

    fn feedback_for_body(&self, body: &LumpedEnthalpyBody<'_>) -> Result<f64, IncandescentError> {
        let caloric = body.phase_curve();
        if let PreparedFeedback::SpecificEnthalpy { curve, .. } = &self.feedback
            && *curve != caloric.identity()
        {
            return Err(IncandescentError::Invalid(
                "prepared electrical feedback belongs to a different caloric chart",
            ));
        }
        let volume = self.source.length_m() * self.source.area_m2();
        let thermal_volume = body.characteristic_length_m() * body.surface_area_m2();
        if !volume.is_finite()
            || volume <= 0.0
            || !thermal_volume.is_finite()
            || thermal_volume <= 0.0
            || (volume / thermal_volume - 1.0).abs() > 64.0 * f64::EPSILON
            || caloric.material_card_identity()
                != self.source.source().original().material().card_identity()
        {
            return Err(IncandescentError::Invalid(
                "electrical geometry must span the whole body on the same original material card",
            ));
        }
        let knots = caloric.knots();
        let [low, high] = self.source.temperature_range_k();
        if knots[0].temperature_k < low || knots[knots.len() - 1].temperature_k > high {
            return Err(IncandescentError::Invalid(
                "electrical source support must cover the complete caloric chart",
            ));
        }
        let material = self.source.source().original().material();
        let phase_axis = self.source.liquid_mass_fraction_axis();
        // A card may contain measurements at many pressures or compositions.
        // Sharing its hash alone does not bind those physical conditions.
        // Authored charts carry no receipts and retain their declared authority.
        for receipt in caloric.source_receipts() {
            for (axis, value) in &receipt.query_point {
                if receipt.axis_quantities.get(axis) != material.axis_quantities().get(axis)
                    || (axis != self.source.temperature_axis()
                        && Some(axis.as_str()) != phase_axis
                        && material.state_coordinate(axis) != Some(*value))
                {
                    return Err(IncandescentError::Invalid(
                        "caloric and electrical source conditions or quantity conventions differ",
                    ));
                }
            }
        }
        let varying_phase =
            knots[0].liquid_mass_fraction != knots[knots.len() - 1].liquid_mass_fraction;
        let (receipts, evolving_phase) = match body.transport() {
            LumpedEnthalpyTransport::Temperature(transport) => {
                (transport.conductivity_receipts(), false)
            }
            LumpedEnthalpyTransport::MaterialPhase(transport) => {
                if transport.conductivity().temperature_axis() != self.source.temperature_axis()
                    || transport.conductivity().liquid_fraction_axis() != phase_axis
                {
                    return Err(IncandescentError::Invalid(
                        "electrical and thermal phase sources need the same declared temperature and liquid-fraction coordinates",
                    ));
                }
                // Phase preparation already admits equal conductivity/emissivity
                // context across all paths. Its original receipt retains that
                // context without revisiting every path knot on every time step.
                (
                    std::slice::from_ref(transport.radiation().reference_emissivity().receipt()),
                    true,
                )
            }
        };
        for receipt in receipts {
            for (axis, value) in &receipt.query_point {
                let fixed_coordinate_matches = if Some(axis.as_str()) == phase_axis {
                    evolving_phase || (!varying_phase && *value == knots[0].liquid_mass_fraction)
                } else {
                    axis == self.source.temperature_axis()
                        || material.state_coordinate(axis) == Some(*value)
                };
                if receipt.axis_quantities.get(axis) != material.axis_quantities().get(axis)
                    || !fixed_coordinate_matches
                {
                    return Err(IncandescentError::Invalid(
                        "thermal transport must share fixed conditions and cannot freeze an evolving phase coordinate",
                    ));
                }
            }
        }
        // The material thermal owner bounds its ORIGINAL Cp law, including
        // nonlinear sensible storage; latent plateaus contribute zero dT/dh.
        let feedback = match self.feedback {
            PreparedFeedback::SpecificEnthalpy {
                maximum_positive_power_slope: slope,
                ..
            } => {
                if slope == 0.0 {
                    0.0
                } else {
                    (slope / body.mass_kg()).next_up()
                }
            }
            PreparedFeedback::Temperature(slope) => {
                let temperature_slope = caloric.maximum_temperature_slope();
                if slope == 0.0 || temperature_slope == 0.0 {
                    0.0
                } else {
                    let direct = slope * (temperature_slope / body.mass_kg());
                    let feedback = if direct.is_finite() && direct > 0.0 {
                        direct
                    } else {
                        (slope / body.mass_kg()) * temperature_slope
                    };
                    feedback.next_up()
                }
            }
        };
        if !feedback.is_finite() || feedback < 0.0 {
            return Err(IncandescentError::Invalid(
                "positive electrothermal feedback is unrepresentable",
            ));
        }
        Ok(feedback)
    }

    fn trial_electrical_state(
        &self,
        curve: &EquilibriumEnthalpyPhaseCurve,
        enthalpy_j_kg: f64,
        temperature_k: f64,
    ) -> Result<IncandescentElectricalState, IncandescentError> {
        // Reuse the native numeric chart without hashing a trial state at each
        // bisection. The native endpoint already supplies its exact temperature.
        let fraction = if self.source.liquid_mass_fraction_axis().is_some() {
            curve
                .values_at_specific_enthalpy(enthalpy_j_kg)?
                .liquid_mass_fraction
        } else {
            0.0
        };
        Ok(self.electrical_state(temperature_k, fraction)?)
    }

    fn electrical_state(
        &self,
        temperature_k: f64,
        liquid_mass_fraction: f64,
    ) -> Result<IncandescentElectricalState, ConductorError> {
        let conductor = if self.source.liquid_mass_fraction_axis().is_some() {
            self.source
                .resolve_at_temperature_and_liquid_mass_fraction(
                    temperature_k,
                    liquid_mass_fraction,
                )?
        } else {
            self.source.resolve_at_temperature_k(temperature_k)?
        };
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
