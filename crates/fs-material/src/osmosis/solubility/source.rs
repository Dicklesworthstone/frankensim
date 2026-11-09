//! Material-card admission for the original neutral finite-solubility owners.
//!
//! A complete, explicitly selected van't Hoff model supplies equilibrium;
//! an independent, explicitly selected film model supplies kinetics. The full
//! immutable cards, model identities, calibration hashes and provenance remain
//! inspectable. Admission verifies the declared source convention and support;
//! it neither calibrates a material nor upgrades a citation to measured evidence.
//!
//! Every model declares exactly three validity coordinates: semantic absolute
//! `temperature`, and dimensionless singleton-1 flags `source-solvent:<id>` and
//! `source-solute:<id>`. The exact species IDs are part of the original model
//! content identity, so a water/solute fit cannot be relabelled for another pair.
//! The equilibrium model means mol of neutral solute per m³ of incompressible
//! solvent, its pure solid, ideal dilute activities and constant dissolution
//! enthalpy. It does not supply ionic speciation or phase heat-capacity curves.
//!
//! Film conductance kA is an apparatus/geometry coefficient, not an intrinsic
//! material diffusivity. The fixed-film source declares its prescribed area and
//! precipitation policy; the particle source declares kA at its original solid
//! reference stock and fixed count/shape. Neither bath volume nor material name
//! rescales it. Changing area, mixing or film thickness requires a new declared
//! source. Only the original native particle law evolves its admitted area.
//!
//! Source construction prepares native owners once; runtime phase partition,
//! particle kinetics, osmotic particles, transport and caloric response keep
//! using the existing solvers. Cloning a numeric law into another owner does
//! not itself copy these receipts: retain this source binding alongside it.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, MaterialCard};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, AmountConcentration, Dims, QuantitySpec, Temperature, VolumetricFlowRate};
use fs_thermochem::{
    MolarEnergyQuantityV1, SpeciesId,
    solution::solubility::{
        IdealSolubilityV1, SolubilityErrorV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticsV1},
    },
};

/// Version 1: constant-enthalpy van't Hoff saturation on solvent-volume basis.
/// State schema 0, `ZeroInternalState`, and exactly the three reference
/// parameters below are required.
pub const NEUTRAL_SOLUBILITY_LAW: &str = "neutral-solubility-vant-hoff";
/// Version 1: prescribed fixed kA, with an explicit dimensionless
/// `precipitation_policy` (0 = existing solid; 1 = persistent seeded surface).
/// State schema 1 requires caller-declared actual dissolved and solid amounts.
pub const FIXED_FILM_KINETICS_LAW: &str = "neutral-solubility-fixed-film";
/// Version 1: kA at `reference_solid_amount`, self-similar fixed-count particles,
/// and existing solid only. State schema 1 requires declared actual phases.
pub const PARTICLE_FILM_KINETICS_LAW: &str = "neutral-solubility-particle-film";
/// The sole evolving source coordinate, semantic absolute kelvin.
pub const TEMPERATURE_AXIS: &str = "temperature";
/// Absolute reference temperature, K.
pub const REFERENCE_TEMPERATURE_PARAMETER: &str = "reference_temperature";
/// Reference saturation, mol/m³ of solvent, not mass fraction or molality.
pub const REFERENCE_SOLUBILITY_PARAMETER: &str = "reference_solubility";
/// Dissolved-minus-solid enthalpy, J/mol; positive dissolution absorbs heat.
pub const DISSOLUTION_ENTHALPY_PARAMETER: &str = "dissolution_enthalpy";
/// Prescribed fixed conductance, or conductance at the declared particle stock,
/// m³/s. A diffusivity [m²/s], mass-transfer velocity [m/s] or kLa [1/s] refuses.
pub const TRANSFER_CONDUCTANCE_PARAMETER: &str = "transfer_conductance";
/// Original particle stock at the declared kA, mol; never reset after transport.
pub const REFERENCE_SOLID_PARAMETER: &str = "reference_solid_amount";
/// Dimensionless discrete fixed-film policy: exactly 0 or 1 as documented above.
pub const PRECIPITATION_POLICY_PARAMETER: &str = "precipitation_policy";

/// An original material card, its explicitly pinned model, and the admitted
/// neutral saturation law. All source storage is immutable and shared on clone.
#[derive(Clone, Debug, PartialEq)]
pub struct NeutralSolubilitySource {
    material: Arc<MaterialCard>,
    model_index: usize,
    law: IdealSolubilityV1,
}

impl NeutralSolubilitySource {
    /// Admit a full model identity from this exact card for the explicitly
    /// source-qualified solvent/solute pair. The source must declare finite,
    /// positive temperature support containing the reference temperature.
    ///
    /// # Errors
    /// Missing/ambiguous full model identity, foreign pair, absent calibration
    /// hashes, unsupported law/state schema, extra/missing parameters, wrong SI
    /// dimensions, unsupported conditions and original numeric-law refusals.
    pub fn from_material(
        material: Arc<MaterialCard>,
        model: ContentHash,
        solvent: SpeciesId,
        solute: SpeciesId,
    ) -> Result<Self, NeutralSolubilitySourceError> {
        let model_index = selected_model(&material, model)?;
        let source = &material.models()[model_index];
        source.validate()?;
        require(
            source.law.0 == NEUTRAL_SOLUBILITY_LAW
                && source.law_version == 1
                && source.state_schema_version == 0
                && source.initial_state == InitialStatePolicy::ZeroInternalState,
            "require version-one neutral van't Hoff law with stateless schema zero",
        )?;
        let domain = source_domain(source, &solvent, &solute)?;
        let parameters = exact_parameters(
            source,
            &[
                (REFERENCE_TEMPERATURE_PARAMETER, Temperature::DIMS),
                (REFERENCE_SOLUBILITY_PARAMETER, AmountConcentration::DIMS),
                (DISSOLUTION_ENTHALPY_PARAMETER, MolarEnergyQuantityV1::DIMS),
            ],
        )?;
        let law = IdealSolubilityV1::new(
            solvent,
            solute,
            Temperature::new(parameters[0]),
            AmountConcentration::new(parameters[1]),
            parameters[2],
            domain,
        )?;
        Ok(Self {
            material,
            model_index,
            law,
        })
    }

    /// Complete original immutable material card, including unrelated claims.
    pub const fn material(&self) -> &Arc<MaterialCard> {
        &self.material
    }

    /// Exactly selected model, including original applicability and sources.
    pub fn model_card(&self) -> &ConstitutiveModelCard {
        &self.material.models()[self.model_index]
    }

    /// Original native saturation/phase partition owner, ready for the existing
    /// finite-solution, caloric and transport adapters without another model.
    pub const fn law(&self) -> &IdealSolubilityV1 {
        &self.law
    }

    /// Bind an independently sourced film coefficient for the same exact pair.
    /// The two original temperature domains are intersected, never extended.
    /// A narrower domain that excludes the thermodynamic reference temperature
    /// rebases the same constant-enthalpy curve at its nearest admitted endpoint
    /// using the original saturation evaluator. This incurs ordinary floating
    /// point evaluation; the complete original thermodynamic fit is retained.
    ///
    /// # Errors
    /// The source/schema/unit refusals of [`Self::from_material`], unsupported
    /// fixed-film policy, disjoint support, or native kinetic-parameter refusal.
    pub fn with_kinetics(
        &self,
        material: Arc<MaterialCard>,
        model: ContentHash,
    ) -> Result<NeutralSolubilityKineticsSource, NeutralSolubilitySourceError> {
        let model_index = selected_model(&material, model)?;
        let source = &material.models()[model_index];
        source.validate()?;
        require(
            source.law_version == 1
                && source.state_schema_version == 1
                && source.initial_state == InitialStatePolicy::RequiresDeclaredState,
            "film kinetics require version-one law/state and declared actual phase inventories",
        )?;
        let domain = source_domain(source, self.law.solvent(), self.law.solute())?;
        let (conductance, policy, reference_solid) = match source.law.0.as_str() {
            FIXED_FILM_KINETICS_LAW => {
                let parameters = exact_parameters(
                    source,
                    &[
                        (TRANSFER_CONDUCTANCE_PARAMETER, VolumetricFlowRate::DIMS),
                        (PRECIPITATION_POLICY_PARAMETER, Dims::NONE),
                    ],
                )?;
                let policy = match parameters[1] {
                    0.0 => PrecipitationPolicyV1::ExistingSolidOnly,
                    1.0 => PrecipitationPolicyV1::SeededSurface,
                    _ => {
                        return Err(NeutralSolubilitySourceError::Invalid(
                            "fixed-film precipitation policy must be exactly zero or one",
                        ));
                    }
                };
                (parameters[0], policy, None)
            }
            PARTICLE_FILM_KINETICS_LAW => {
                let parameters = exact_parameters(
                    source,
                    &[
                        (TRANSFER_CONDUCTANCE_PARAMETER, VolumetricFlowRate::DIMS),
                        (REFERENCE_SOLID_PARAMETER, Amount::DIMS),
                    ],
                )?;
                (
                    parameters[0],
                    PrecipitationPolicyV1::ExistingSolidOnly,
                    Some(Amount::new(parameters[1])),
                )
            }
            _ => {
                return Err(NeutralSolubilitySourceError::Invalid(
                    "require explicitly fixed-film or self-similar-particle kinetic source",
                ));
            }
        };
        let original_domain = self.law.temperature_domain();
        let low = domain[0].value().max(original_domain[0].value());
        let high = domain[1].value().min(original_domain[1].value());
        require(
            low <= high,
            "equilibrium and kinetic source supports do not overlap",
        )?;
        let reference = Temperature::new(self.law.reference_temperature().value().clamp(low, high));
        let law = IdealSolubilityV1::new(
            self.law.solvent().clone(),
            self.law.solute().clone(),
            reference,
            self.law.saturation_concentration(reference)?,
            self.law.dissolution_enthalpy_j_mol(),
            [Temperature::new(low), Temperature::new(high)],
        )?;
        let mut kinetics =
            SolubilityKineticsV1::new(law, VolumetricFlowRate::new(conductance), policy)?;
        if let Some(reference) = reference_solid {
            kinetics = kinetics.with_particle_surface(reference)?;
        }
        Ok(NeutralSolubilityKineticsSource {
            equilibrium: self.clone(),
            material,
            model_index,
            kinetics,
        })
    }
}

/// Independently sourced film kinetics retaining both complete original fits.
#[derive(Clone, Debug, PartialEq)]
pub struct NeutralSolubilityKineticsSource {
    equilibrium: NeutralSolubilitySource,
    material: Arc<MaterialCard>,
    model_index: usize,
    kinetics: SolubilityKineticsV1,
}

impl NeutralSolubilityKineticsSource {
    /// Original equilibrium source before the kinetic support intersection.
    pub const fn equilibrium_source(&self) -> &NeutralSolubilitySource {
        &self.equilibrium
    }

    /// Complete original film source card and its apparatus/reference-area basis.
    pub const fn material(&self) -> &Arc<MaterialCard> {
        &self.material
    }

    /// Exact selected kinetic model, with original source hashes and support.
    pub fn model_card(&self) -> &ConstitutiveModelCard {
        &self.material.models()[self.model_index]
    }

    /// Native finite-bath kinetics. Its pair law is valid only on the common
    /// temperature support; use that same law when preparing caloric consumers.
    pub const fn kinetics(&self) -> &SolubilityKineticsV1 {
        &self.kinetics
    }
}

/// Original source/physics refusal; no partial binding is published.
#[derive(Clone, Debug, PartialEq)]
pub enum NeutralSolubilitySourceError {
    Material(MatDbError),
    Solubility(SolubilityErrorV1),
    Invalid(&'static str),
}
impl core::fmt::Display for NeutralSolubilitySourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(source) => write!(f, "neutral solubility material source: {source}"),
            Self::Solubility(source) => write!(f, "neutral solubility source parameters: {source}"),
            Self::Invalid(what) => write!(f, "neutral solubility source refusal: {what}"),
        }
    }
}
impl std::error::Error for NeutralSolubilitySourceError {}
impl From<MatDbError> for NeutralSolubilitySourceError {
    fn from(value: MatDbError) -> Self {
        Self::Material(value)
    }
}
impl From<SolubilityErrorV1> for NeutralSolubilitySourceError {
    fn from(value: SolubilityErrorV1) -> Self {
        Self::Solubility(value)
    }
}

fn selected_model(
    material: &MaterialCard,
    identity: ContentHash,
) -> Result<usize, NeutralSolubilitySourceError> {
    let mut selected = None;
    for (index, candidate) in material.models().iter().enumerate() {
        if candidate.content_hash() == identity {
            require(
                selected.is_none(),
                "selected neutral solubility model is ambiguous",
            )?;
            selected = Some(index);
        }
    }
    selected.ok_or(NeutralSolubilitySourceError::Invalid(
        "selected neutral solubility model is absent from the original material card",
    ))
}

fn source_domain(
    source: &ConstitutiveModelCard,
    solvent: &SpeciesId,
    solute: &SpeciesId,
) -> Result<[Temperature; 2], NeutralSolubilitySourceError> {
    require(
        !source.sources.is_empty(),
        "require original solubility calibration artifact hashes",
    )?;
    let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    require(
        solvent != solute
            && source.validity.bounds().len() == 3
            && source.validity.axis_quantities().len() == 3
            && source.validity.axis_quantities().get(TEMPERATURE_AXIS)
                == Some(&absolute_temperature),
        "source needs exactly semantic absolute temperature and two distinct exact species qualifiers",
    )?;
    for axis in [
        format!("source-solvent:{}", solvent.as_str()),
        format!("source-solute:{}", solute.as_str()),
    ] {
        require(
            source.validity.bound(&axis) == Some((1.0, 1.0))
                && source.validity.axis_quantities().get(&axis)
                    == Some(&QuantitySpec::dimensional(Dims::NONE)),
            "source does not qualify this exact solvent/solute pair",
        )?;
    }
    let (low, high) =
        source
            .validity
            .bound(TEMPERATURE_AXIS)
            .ok_or(NeutralSolubilitySourceError::Invalid(
                "source is missing temperature support",
            ))?;
    require(
        low.is_finite() && high.is_finite() && low > 0.0 && low <= high,
        "source temperature support must be finite, positive and ordered",
    )?;
    Ok([Temperature::new(low), Temperature::new(high)])
}

fn exact_parameters<const N: usize>(
    source: &ConstitutiveModelCard,
    expected: &[(&str, Dims); N],
) -> Result<[f64; N], NeutralSolubilitySourceError> {
    require(
        source.parameters.len() == N,
        "source requires exactly the declared solubility or film parameter block",
    )?;
    let mut values = [0.0; N];
    for (i, &(name, dims)) in expected.iter().enumerate() {
        let parameter =
            source
                .parameters
                .get(name)
                .ok_or(NeutralSolubilitySourceError::Invalid(
                    "missing required solubility source parameter",
                ))?;
        require(
            parameter.dims == dims && parameter.value.is_finite(),
            "solubility source parameter has wrong SI dimensions or is not finite",
        )?;
        values[i] = parameter.value;
    }
    Ok(values)
}

fn require(admitted: bool, what: &'static str) -> Result<(), NeutralSolubilitySourceError> {
    if admitted {
        Ok(())
    } else {
        Err(NeutralSolubilitySourceError::Invalid(what))
    }
}

#[cfg(test)]
mod tests;
