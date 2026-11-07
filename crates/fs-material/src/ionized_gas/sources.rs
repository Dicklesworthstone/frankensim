//! Momentum laws retained independently of which carriers currently exist.
//!
//! Finite-rate ionization may create the first electron, exhaust an ionic
//! stage, or expose a previously absent neutral target. An evaluated transport
//! snapshot only retains present rows. This bank instead admits every charged
//! species and, for kinetic electron drag, every neutral in the original
//! thermodynamic bank. Each query selects the actual populated rows and invokes
//! the existing material/collision owners at the actual Th, Te and pressure.
//! No population floor, LTE projection or guessed collision source is used.

use super::refresh::{poll, refresh_material_source, transport_error};
use super::{
    CarrierMomentumRelaxation, IonizedGasConductivity, IonizedGasDcLimits, IonizedGasError,
    IonizedGasRefreshError, SpeciesId,
};
use fs_matdb::MaterialCard;
use fs_thermochem::coulomb::{CoulombErrorV1, CoulombModelV1};
use fs_thermochem::electron_collision::{
    ElectronCollisionErrorV1, ElectronCollisionKindV1, ElectronCrossSectionV1,
    MaxwellianRateControlsV1, MaxwellianRateWeightingV1,
};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;

/// One authored momentum law, independent of its carrier's current abundance.
#[derive(Clone, Debug, PartialEq)]
pub enum IonizedGasMomentumSource {
    /// Original resolved rate and immutable card. The template must be material
    /// backed; all selected claims, requirements and fixed coordinates survive
    /// each query. Its source must support the mixture/background being modeled.
    Material {
        template: CarrierMomentumRelaxation,
        original_card: MaterialCard,
    },
    /// Electron drag from one momentum-transfer law per original neutral row,
    /// including zero-population targets. The same integration controls and
    /// stationary-heavy domain apply after every population change.
    ElectronNeutral {
        cross_sections: Vec<ElectronCrossSectionV1>,
        controls: MaxwellianRateControlsV1,
        maximum_target_thermal_speed_ratio: f64,
    },
}

#[derive(Clone, Debug, PartialEq)]
struct SourceRow {
    carrier: usize,
    source: IonizedGasMomentumSource,
    // Canonical target indices, prepared once instead of searched each query.
    neutral_indices: Vec<usize>,
}

/// Complete immutable source support for a reacting rigid gas parcel.
///
/// Construction covers every charged row in the original thermodynamic bank,
/// including absent carriers. Source selection and neutral-target indices are
/// prepared once. Evaluation integrates only present targets, constructs a
/// bounded transport snapshot, and reuses the native DC, Coulomb,
/// harmonic and magnetized response owners without duplicating their physics.
#[derive(Clone, Debug, PartialEq)]
pub struct IonizedGasTransportBank {
    reference: FrozenIonizedGasStateV1,
    rows: Vec<SourceRow>,
    coulomb: Option<CoulombModelV1>,
    limits: IonizedGasDcLimits,
}

impl IonizedGasTransportBank {
    /// Admit complete sources even when the reference parcel is exactly neutral.
    ///
    /// Every original charged species needs exactly one law. Electron-neutral
    /// tables cover every original neutral species exactly once, with no ionic
    /// or foreign targets. Material laws and all cross sections are queried at
    /// the reference temperatures to check their initial source support. Later
    /// endpoints are rechecked independently; no intervening path is certified.
    /// An optional Coulomb model must bind this exact reference population and
    /// retains its screening, cutoff and limits for all subsequent queries.
    pub fn new(
        reference: &FrozenIonizedGasStateV1,
        sources: Vec<IonizedGasMomentumSource>,
        coulomb: Option<CoulombModelV1>,
        limits: IonizedGasDcLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasRefreshError> {
        use IonizedGasRefreshError::Invalid;
        poll(&mut cancelled)?;
        limits.validate().map_err(transport_error)?;
        let chemistry = reference.thermodynamics();
        let models = chemistry.models();
        let charged = (0..models.len())
            .filter(|&i| chemistry.charge_vector().get(i) != Some(0))
            .count();
        if sources.len() != charged || charged == 0 {
            return Err(Invalid(
                "one momentum law per original charged species is required",
            ));
        }
        if coulomb
            .as_ref()
            .is_some_and(|model| !model.is_bound_to_populations(reference))
        {
            return Err(Invalid("Coulomb source must bind the reference gas"));
        }
        let mut rows = Vec::with_capacity(charged);
        for mut source in sources {
            poll(&mut cancelled)?;
            let mut neutral_indices = Vec::new();
            let carrier = match &mut source {
                IonizedGasMomentumSource::Material {
                    template,
                    original_card,
                } => {
                    let i = models
                        .binary_search_by(|model| model.species().cmp(template.species()))
                        .map_err(|_| Invalid("momentum source names a foreign species"))?;
                    if chemistry.charge_vector().get(i) == Some(0) {
                        return Err(Invalid("a neutral cannot own a charged-carrier law"));
                    }
                    refresh_material_source(template, original_card, reference)?;
                    i
                }
                IonizedGasMomentumSource::ElectronNeutral {
                    cross_sections,
                    controls,
                    maximum_target_thermal_speed_ratio,
                } => {
                    let i = reference
                        .electron_index()
                        .ok_or(Invalid("structural electron source required"))?;
                    if !maximum_target_thermal_speed_ratio.is_finite()
                        || *maximum_target_thermal_speed_ratio <= 0.0
                        || *maximum_target_thermal_speed_ratio >= 1.0
                    {
                        return Err(Invalid("stationary-target speed bound must lie in (0,1)"));
                    }
                    neutral_indices = (0..models.len())
                        .filter(|&i| chemistry.charge_vector().get(i) == Some(0))
                        .collect();
                    if cross_sections.len() != neutral_indices.len() || cross_sections.is_empty() {
                        return Err(Invalid(
                            "momentum tables must cover every original neutral species",
                        ));
                    }
                    cross_sections.sort_by(|a, b| a.target_species().cmp(b.target_species()));
                    for (&target, table) in neutral_indices.iter().zip(cross_sections) {
                        poll(&mut cancelled)?;
                        if table.target_species() != models[target].species()
                            || table.kind() != ElectronCollisionKindV1::ElasticMomentumTransfer
                        {
                            return Err(Invalid(
                                "neutral momentum tables must be unique original targets",
                            ));
                        }
                        table
                            .evaluate_properties(
                                reference.electron_temperature(),
                                MaxwellianRateWeightingV1::WeakDriftMomentum,
                                *controls,
                                &mut cancelled,
                            )
                            .map_err(collision_error)?;
                    }
                    i
                }
            };
            rows.push(SourceRow {
                carrier,
                source,
                neutral_indices,
            });
        }
        rows.sort_by_key(|row| row.carrier);
        if rows
            .windows(2)
            .any(|pair| pair[0].carrier == pair[1].carrier)
        {
            return Err(Invalid("duplicate charged-carrier momentum law"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            reference: reference.clone(),
            rows,
            coulomb,
            limits,
        })
    }

    /// Original bank/volume binding and reference source-query point.
    pub const fn reference_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.reference
    }
    /// Canonical carrier identities and their complete immutable laws.
    pub fn sources(
        &self,
    ) -> impl ExactSizeIterator<Item = (&SpeciesId, &IonizedGasMomentumSource)> {
        self.rows.iter().map(|row| {
            (
                self.reference.thermodynamics().models()[row.carrier].species(),
                &row.source,
            )
        })
    }
    /// Original optional Coulomb prescription, never a stale evaluated pair set.
    pub const fn coulomb_source(&self) -> Option<&CoulombModelV1> {
        self.coulomb.as_ref()
    }
    pub const fn limits(&self) -> IonizedGasDcLimits {
        self.limits
    }

    /// Evaluate actual populations, including carrier birth and extinction.
    ///
    /// `None` means exactly no charged particles: the modeled free-carrier
    /// current and heating vanish. It does not infer neutral dielectric loss or
    /// supply a positive conductivity to a PDE requiring ellipticity. Every
    /// positive carrier amount, however small, is evaluated or explicitly
    /// refused; no cutoff hides trace carriers or substitutes a conductivity
    /// floor. A nonempty result feeds the existing DC/harmonic/magnetic owners.
    ///
    /// Volume and thermodynamic sources remain fixed. Actual temperatures,
    /// pressure, species support and Coulomb pairs are refreshed atomically.
    /// A missing neutral background, source-domain exit or cancellation refuses
    /// the complete query without changing this bank or any previous snapshot.
    pub fn evaluate(
        &self,
        state: &FrozenIonizedGasStateV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Option<IonizedGasConductivity>, IonizedGasRefreshError> {
        poll(&mut cancelled)?;
        if state.thermodynamics() != self.reference.thermodynamics()
            || state.volume() != self.reference.volume()
        {
            return Err(IonizedGasRefreshError::Invalid(
                "transport bank requires the original thermodynamics and volume",
            ));
        }
        let mut rates = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            poll(&mut cancelled)?;
            if state.amounts()[row.carrier].value() == 0.0 {
                continue;
            }
            let rate = match &row.source {
                IonizedGasMomentumSource::Material {
                    template,
                    original_card,
                } => refresh_material_source(template, original_card, state)?,
                IonizedGasMomentumSource::ElectronNeutral {
                    cross_sections,
                    controls,
                    maximum_target_thermal_speed_ratio,
                } => {
                    let mut evaluations = Vec::with_capacity(cross_sections.len());
                    for (&target, table) in row.neutral_indices.iter().zip(cross_sections) {
                        poll(&mut cancelled)?;
                        if state.amounts()[target].value() == 0.0 {
                            continue;
                        }
                        evaluations.push(
                            table
                                .evaluate(
                                    state.electron_temperature(),
                                    MaxwellianRateWeightingV1::WeakDriftMomentum,
                                    *controls,
                                    &mut cancelled,
                                )
                                .map_err(collision_error)?,
                        );
                    }
                    CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
                        state,
                        evaluations,
                        *maximum_target_thermal_speed_ratio,
                        &mut cancelled,
                    )
                    .map_err(transport_error)?
                }
            };
            rates.push(rate);
        }
        if rates.is_empty() {
            poll(&mut cancelled)?;
            return Ok(None);
        }
        let mut result =
            IonizedGasConductivity::from_state(state.clone(), rates, self.limits, &mut cancelled)
                .map_err(transport_error)?;
        if let Some(source) = &self.coulomb {
            let model = CoulombModelV1::new_with_cutoff(
                state,
                source.screening(),
                source.cutoff(),
                source.limits(),
                &mut cancelled,
            )
            .map_err(|error| match error {
                CoulombErrorV1::Cancelled => IonizedGasRefreshError::Cancelled,
                error => transport_error(IonizedGasError::Coulomb(error)),
            })?;
            result = result
                .with_coulomb_drag(model, &mut cancelled)
                .map_err(transport_error)?;
        }
        poll(&mut cancelled)?;
        Ok(Some(result))
    }
}

fn collision_error(error: ElectronCollisionErrorV1) -> IonizedGasRefreshError {
    match error {
        ElectronCollisionErrorV1::Cancelled => IonizedGasRefreshError::Cancelled,
        error => IonizedGasRefreshError::Collision(error),
    }
}
