//! Thermal coefficient changes at fixed total strain and fixed crack geometry.
//!
//! The constitutive ordering is a coefficient change at held strain/damage,
//! followed by irreversible relaxation onto the new fracture envelope. The
//! native cohesive temperature owner retains the original physical onset mode,
//! damage, and accumulated fracture work. Its area accounts are divided by the
//! original band width; the positive in-plane elastic remainder contributes its
//! own signed parameter work. No mechanical strain work occurs in this step.
//!
//! Young's modulus, both strengths, both fracture energies and the BK exponent
//! may vary. Poisson's ratio, the plane and band width must remain identical:
//! changing them changes the work-conjugate modal coordinates and needs another
//! history-transfer model. Cooling and increased toughness never heal damage.
//! This endpoint ordering is not an integration of an unseen thermal path,
//! thermal expansion, a heat/entropy law, crushing or thermal fatigue. Fracture
//! dissipation and parameter work are not automatically converted into heat.

use super::source::{BrittleSourceError, ResolvedBrittleCrackBand};
use super::{BrittleCrackBand, BrittleError, BrittleState, add_outer, finite, positive};
use crate::state_point::MaterialStatePointError;
use fs_matdb::MaterialCard;
use fs_qty::QuantitySpec;

/// Accepted coefficient/history transfer, with energy densities in J/m³.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittleTemperatureStep<'a> {
    /// Same total strain and irreversible history, bound to the new owner.
    pub state: BrittleState<'a>,
    /// Signed work supplied by changing coefficients at held strain/damage.
    pub parameter_work_j_m3: f64,
    /// New nonnegative fracture release during relaxation at fixed strain.
    pub damage_dissipated_j_m3: f64,
    /// Complete recoverable storage change, including the in-plane remainder.
    pub stored_change_j_m3: f64,
    /// Parameter work minus storage change minus new fracture release.
    pub defect_j_m3: f64,
}

impl BrittleState<'_> {
    /// Transfer this accepted state to explicitly supplied thermal coefficients.
    ///
    /// The source adapter below additionally verifies the original card, claims
    /// and temperature coordinate. This numerical operation alone does not
    /// establish source applicability. The immutable input remains retryable on
    /// refusal. Only the new state may subsequently be evolved by `target`.
    pub fn transfer_temperature<'b>(
        &self,
        target: &'b BrittleCrackBand,
        energy_tolerance_j_m3: f64,
    ) -> Result<BrittleTemperatureStep<'b>, BrittleError> {
        positive(energy_tolerance_j_m3)?;
        self.owner.check_strain(&self.strain)?;
        target.check_strain(&self.strain)?;
        if !core::ptr::eq(self.point.law(), &self.owner.cohesive) {
            return Err(BrittleError::StateMismatch);
        }
        require_same_projection(self.owner, target)?;
        let width = target.settings.band_width_m;
        let step = self
            .point
            .transfer_temperature(&target.cohesive, positive(energy_tolerance_j_m3 * width)?)?;
        let state = BrittleState {
            owner: target,
            strain: self.strain,
            point: step.point,
        };
        let parallel_change = finite(
            target.parallel_energy(&self.strain) - self.owner.parallel_energy(&self.strain),
        )?;
        let parameter_work = finite(parallel_change + step.parameter_work_j_m2 / width)?;
        let stored_change = finite(state.stored_energy_j_m3() - self.stored_energy_j_m3())?;
        let dissipated = finite(step.damage_dissipated_j_m2 / width)?;
        finite(state.dissipated_j_m3())?;
        let defect = finite((parameter_work - stored_change) - dissipated)?;
        if defect.abs() > energy_tolerance_j_m3 {
            return Err(BrittleError::EnergyDefect(defect));
        }
        Ok(BrittleTemperatureStep {
            state,
            parameter_work_j_m3: parameter_work,
            damage_dissipated_j_m3: dissipated,
            stored_change_j_m3: stored_change,
            defect_j_m3: defect,
        })
    }
}

impl ResolvedBrittleCrackBand {
    /// Resolve the SAME material card and original selected claims at a new
    /// positive absolute temperature. All extra properties are refreshed too.
    /// Other coordinates, geometry, Poisson ratio and model controls stay fixed.
    /// This prepares coefficients; use [`Self::transfer_temperature`] to move
    /// existing damage/history, rather than constructing a virgin state.
    pub fn at_temperature(
        &self,
        card: &MaterialCard,
        temperature_axis: &str,
        temperature_k: f64,
    ) -> Result<Self, BrittleSourceError> {
        require_temperature_axis(self, temperature_axis)?;
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(BrittleSourceError::Invalid(
                "brittle temperature must be finite positive kelvin",
            ));
        }
        let query = self
            .material()
            .query()
            .with(temperature_axis, temperature_k)
            .map_err(|source| MaterialStatePointError::Query {
                property: format!("brittle temperature:{temperature_axis}"),
                source,
            })?;
        let mut next = Self::try_from_material(
            self.material().at_query_point(card, &query)?,
            self.law().settings(),
        )?;
        if self.law.elastic.poisson != next.law.elastic.poisson {
            return Err(BrittleSourceError::Invalid(
                "brittle thermal transfer cannot change the Poisson modal projection",
            ));
        }
        // Reuse the original prepared geometry exactly. Re-normalizing even a
        // unit normal can change an oblique mode by an ulp; that would move the
        // held cohesive jump without any physical strain increment.
        next.law.settings.normal = self.law.settings.normal;
        next.law.axes = self.law.axes;
        next.law.modes = self.law.modes;
        next.law.parallel_modes = self.law.parallel_modes;
        next.law.parallel_tangent = [[0.0; 6]; 6];
        for (mode, modulus) in next.law.parallel_modes.iter().zip(next.law.parallel_moduli) {
            add_outer(&mut next.law.parallel_tangent, mode, mode, modulus);
        }
        self.require_temperature_successor(&next, temperature_axis)?;
        Ok(next)
    }

    /// Admit a separately prepared thermal successor without changing source
    /// selection, fixed external conditions, model bounds or modal geometry.
    pub fn require_temperature_successor(
        &self,
        next: &Self,
        temperature_axis: &str,
    ) -> Result<(), BrittleSourceError> {
        require_temperature_axis(self, temperature_axis)?;
        require_temperature_axis(next, temperature_axis)?;
        let a = self.material();
        let b = next.material();
        if a.card_identity() != b.card_identity()
            || a.material() != b.material()
            || a.axis_quantities() != b.axis_quantities()
            || self.law().settings() != next.law().settings()
        {
            return Err(BrittleSourceError::Invalid(
                "brittle temperature successor changed original source, geometry or model controls",
            ));
        }
        require_same_projection(self.law(), next.law())?;
        if a.query_point().len() != b.query_point().len()
            || a.query_point()
                .iter()
                .zip(b.query_point())
                .any(|((name, value), (other, new))| {
                    name != other || (name != temperature_axis && value != new)
                })
        {
            return Err(BrittleSourceError::Invalid(
                "brittle temperature successor changed a held external coordinate",
            ));
        }
        if a.properties().len() != b.properties().len()
            || a.properties().iter().zip(b.properties()).any(|(old, new)| {
                old.requirement() != new.requirement()
                    || old.answer().receipt.selected != new.answer().receipt.selected
            })
        {
            return Err(BrittleSourceError::Invalid(
                "brittle temperature successor changed original selected claims",
            ));
        }
        Ok(())
    }

    /// Transfer the original accepted history through a source-admitted thermal
    /// change. Even numerically identical coefficients require matching sources.
    pub fn transfer_temperature<'b>(
        &self,
        state: &BrittleState<'_>,
        target: &'b Self,
        temperature_axis: &str,
        energy_tolerance_j_m3: f64,
    ) -> Result<BrittleTemperatureStep<'b>, BrittleSourceError> {
        if !core::ptr::eq(state.owner, self.law()) {
            return Err(BrittleError::StateMismatch.into());
        }
        self.require_temperature_successor(target, temperature_axis)?;
        Ok(state.transfer_temperature(target.law(), energy_tolerance_j_m3)?)
    }
}

fn require_same_projection(
    original: &BrittleCrackBand,
    target: &BrittleCrackBand,
) -> Result<(), BrittleError> {
    if original.settings.band_width_m != target.settings.band_width_m
        || original.axes != target.axes
        || original.elastic.poisson != target.elastic.poisson
        || original.modes != target.modes
    {
        return Err(BrittleError::Invalid(
            "brittle thermal transfer requires the same width, plane and Poisson ratio",
        ));
    }
    Ok(())
}

fn require_temperature_axis(
    source: &ResolvedBrittleCrackBand,
    axis: &str,
) -> Result<(), BrittleSourceError> {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    let material = source.material();
    let quantity = material.axis_quantities().get(axis).copied();
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let legacy = axis == "T"
        && (quantity.is_none()
            || quantity == Some(QuantitySpec::dimensional(fs_qty::Temperature::DIMS)));
    if !(legacy || quantity == Some(absolute))
        || !material
            .state_coordinate(axis)
            .is_some_and(|value| value.is_finite() && value > 0.0)
        || material
            .properties()
            .iter()
            .any(|property| !property.has_only_fixed_external_conditions(Some(axis), false))
    {
        return Err(BrittleSourceError::Invalid(
            "brittle thermal transfer requires one original absolute-temperature axis",
        ));
    }
    Ok(())
}
