//! State-dependent kinetic scale from the original material-property owner.
//!
//! The scale multiplies dimensionless activities p_i/p0, so its units are
//! mol/(m3 s), NOT the order-dependent units of a concentration rate constant.
//! The caller binds a law to this reaction and explicitly declares its p0.
//! The original thermodynamic reverse closure supplies detailed balance.
//! No extra species inventory, caloric law, heat credit or time integration.

use super::{ClosedGasParcel, GasReaction, GasReactionError, Temperature, invalid, poll, positive};
use fs_matdb::MaterialCard;
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Dims, Pressure};
use fs_thermochem::equilibrium::{
    DeclaredStoichiometricMassActionEvaluationV1, DeclaredStoichiometricMassActionV1,
    IdealGasSpeciesActivityV1, ReactionProgressRateQuantityV1, ReactionProgressRateScaleV1,
    ThermodynamicReverseRateClosureV1,
};

/// Activity-based forward progress-rate scale, in coherent SI mol/(m3 s).
pub const ACTIVITY_RATE_SCALE_PROPERTY: &str = "reaction_progress_rate_scale";
/// Use the existing quantity owner's dimensions rather than a second unit law.
pub const ACTIVITY_RATE_SCALE_DIMS: Dims = ReactionProgressRateQuantityV1::DIMS;

/// Immutable source binding for ONE caller-selected reaction's forward scale.
/// Queries use the actual reacting parcel's complete condition point. Thus
/// temperature/pressure cannot be overridden with a stale independent context.
/// All additional condition axes (for example catalyst/loading data) must be
/// present explicitly in that same parcel point; none are inferred from names.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialRateScale {
    material: MaterialCard,
    selection: MaterialPropertySelection,
    reference_pressure: Pressure,
}
impl MaterialRateScale {
    /// Bind source data expressed for dimensionless p_i/p0 activities. A
    /// concentration-based rate constant must first be converted by its owner;
    /// a 1/s claim is not interchangeable with this mol/(m3 s) property.
    /// Source availability and validity are checked at each actual gas state.
    pub fn new(
        material: MaterialCard, selection: MaterialPropertySelection,
        reference_pressure: Pressure,
    ) -> Result<Self, GasReactionError> {
        if !positive(reference_pressure.value()) {
            return Err(invalid("material kinetic scale needs an explicit positive activity reference pressure"));
        }
        Ok(Self { material, selection, reference_pressure })
    }
    /// Original immutable rate law and its provenance, not cached coefficients.
    #[must_use]
    pub const fn material(&self) -> &MaterialCard { &self.material }
    /// Original explicit claim-selection policy.
    #[must_use]
    pub const fn selection(&self) -> &MaterialPropertySelection { &self.selection }
    /// Reference pressure used when authoring this activity-based scale.
    #[must_use]
    pub const fn reference_pressure(&self) -> Pressure { self.reference_pressure }
}
impl From<MaterialStatePointError> for GasReactionError {
    fn from(error: MaterialStatePointError) -> Self { Self::RateMaterial(error) }
}

impl GasReaction {
    /// Opt in to material-driven kinetics across ALL existing parcel, piston,
    /// nozzle/diffusion and electrical/thermal reaction entry points. The
    /// stoichiometry, species sources and reverse-equilibrium closure are kept.
    /// An already bound source cannot be silently replaced by another source.
    /// This is an initialization choice; it does not mutate an accepted state.
    pub fn with_material_rate(mut self, source: MaterialRateScale) -> Result<Self, GasReactionError> {
        if self.rate_material.is_some() {
            return Err(invalid("reaction already has an explicit material rate source"));
        }
        let reference = self.kinetics.reverse_rate_closure().equilibrium()
            .convention().reference_pressure();
        if source.reference_pressure.value() != reference.value() {
            return Err(invalid("material kinetic scale and reaction use different activity reference pressures"));
        }
        self.rate_material = Some(source);
        Ok(self)
    }
    /// Optional original material source. None retains the original fixed scale.
    #[must_use]
    pub const fn rate_material(&self) -> Option<&MaterialRateScale> { self.rate_material.as_ref() }

    /// Evaluate once from the chemical stage's actual input, after transport
    /// when applicable. Later nonlinear probes never reselect a kinetic rate.
    /// The same rate evaluation and source receipt travel with the stage report.
    pub(super) fn evaluate_rate(
        &self, gas: &ClosedGasParcel, activities: Vec<IdealGasSpeciesActivityV1>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(DeclaredStoichiometricMassActionEvaluationV1, Option<ResolvedMaterialStatePoint>), GasReactionError> {
        let temperature = Temperature::new(gas.state().state().temperature);
        let Some(source) = &self.rate_material else {
            return Ok((self.kinetics.evaluate(temperature, activities)?, None));
        };
        poll(cancelled)?;
        let required = [ScalarPropertyRequirement::try_new(
            ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?];
        let resolved = resolve_material_state_point(
            &source.material, gas.point(), &required, source.selection.clone(),
        )?;
        let scale = resolved.property(ACTIVITY_RATE_SCALE_PROPERTY)
            .ok_or_else(|| invalid("resolved material omitted the required kinetic scale"))?.value_si();
        let scale = ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(scale))
            .map_err(|_| invalid("resolved material scale is outside the original kinetic quantity domain"))?;
        // Reuse the certified equilibrium and declared activity orders. Reverse
        // rates remain thermodynamically tied to this very same forward scale.
        let kinetics = DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(
                self.kinetics.reverse_rate_closure().equilibrium().clone(), scale,
            ), self.kinetics.kinetic_order_convention(),
        );
        poll(cancelled)?;
        let evaluated = kinetics.evaluate(temperature, activities)?;
        poll(cancelled)?;
        Ok((evaluated, Some(resolved)))
    }
}
