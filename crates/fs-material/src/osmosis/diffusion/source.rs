//! State-dependent binary friction from original, immutable material claims.
//!
//! This owner supplies kinetics to the existing solvent-frame solver. It owns
//! actual absolute temperature and all solvent/solute amount concentrations;
//! it does not freeze pressure or another unowned physical coordinate. Source
//! axes must have the matching explicit static quantity declarations. A curve
//! or surface is evaluated by the original material resolver, including its
//! joint interpolation support. Every call qualifies both endpoints AND the
//! arithmetic concentration face against the same originally selected claims.
//!
//! Source resolution allocates queries, evidence and receipts. These are
//! returned to the caller, never cached in the prepared law. The separate
//! constant-data `SolutionMaxwellStefan` path remains allocation-free.

use super::{
    AqueousSolution, MAX_OSMOTIC_SOLUTES, MolarMass, OsmosisError, Result, SolutionDiffusionError,
    SolutionDiffusionReport, SolutionDiffusionSolve, SolutionDiffusionWorkspace,
    SolutionMaxwellStefan, SpeciesId, VirialSolutionModel, poll, positive,
};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{ClaimId, MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::sync::Arc;

/// Maxwell–Stefan binary diffusivity [m²/s], not a tracer/Fick coefficient.
pub const BINARY_SOLUTION_DIFFUSIVITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);

/// Explicit physical assignment for one unordered binary pair. Every pair,
/// including solute/solute pairs, is required; names never infer a coefficient.
#[derive(Clone, Debug)]
pub struct BinarySolutionDiffusionSource {
    pub species: [SpeciesId; 2],
    pub material: Arc<MaterialCard>,
    pub property: PropertyKey,
    /// Applied once at the actual reference solution, then pinned permanently.
    pub selection: MaterialPropertySelection,
}

/// Source coordinate names mapped to their actual physical owners. Exactly one
/// concentration axis per species, INCLUDING solvent, must be supplied.
#[derive(Clone, Debug)]
pub struct SolutionDiffusionSourceAxes {
    pub temperature: String,
    pub concentrations: Vec<(SpeciesId, String)>,
}

/// Immutable selected pair identity. The full original card remains available
/// for receipt replay; unselected claims are not silently discarded.
#[derive(Clone, Debug)]
pub struct PinnedBinarySolutionDiffusionSource {
    species: [SpeciesId; 2],
    material: Arc<MaterialCard>,
    property: PropertyKey,
    selected: ClaimId,
    original_indices: [usize; 2],
}
impl PinnedBinarySolutionDiffusionSource {
    pub fn species(&self) -> &[SpeciesId; 2] {
        &self.species
    }
    pub fn material(&self) -> &MaterialCard {
        &self.material
    }
    pub fn property(&self) -> &PropertyKey {
        &self.property
    }
    pub const fn selected_claim(&self) -> ClaimId {
        self.selected
    }
}

/// Complete original source receipts at the three actual physical samples.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinarySolutionDiffusivity {
    pub species: [SpeciesId; 2],
    pub left: ResolvedMaterialStatePoint,
    pub face: ResolvedMaterialStatePoint,
    pub right: ResolvedMaterialStatePoint,
}

/// One complete source-qualified matrix. The face values enter the shared
/// solver; endpoint values qualify applicability without arithmetic averaging
/// of a nonlinear material law.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedSolutionDiffusivities {
    canonical_species: Box<[SpeciesId]>,
    binary_diffusivities_m2_s: Box<[f64]>,
    pairs: Vec<ResolvedBinarySolutionDiffusivity>,
}
impl ResolvedSolutionDiffusivities {
    pub fn canonical_species(&self) -> &[SpeciesId] {
        &self.canonical_species
    }
    pub fn binary_diffusivities_m2_s(&self) -> &[f64] {
        &self.binary_diffusivities_m2_s
    }
    pub fn pairs(&self) -> &[ResolvedBinarySolutionDiffusivity] {
        &self.pairs
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MaterialSolutionDiffusionReport {
    pub constitutive: SolutionDiffusionReport,
    pub sources: ResolvedSolutionDiffusivities,
}

/// Original nonideal solution law with complete material-sourced binary pairs.
/// Original thermodynamic species order and the shared MS canonical order are
/// mapped explicitly. No kinetic datum is inferred from virial free energy.
#[derive(Clone, Debug)]
pub struct MaterialSolutionMaxwellStefan {
    reference_law: SolutionMaxwellStefan,
    sources: Box<[PinnedBinarySolutionDiffusionSource]>,
    temperature_axis: String,
    /// Original [solvent, thermodynamic solutes...] order.
    concentration_axes: Box<[String]>,
}
impl MaterialSolutionMaxwellStefan {
    /// Select original claims at an actual reference solution. Pair and axis
    /// assignment order is immaterial; solutes retain the virial model's order.
    /// Only typed absolute temperature and typed amount concentration axes are
    /// admitted. Pressure-dependent data need a separate liquid pressure owner.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermodynamics: Arc<VirialSolutionModel>,
        solvent: (SpeciesId, MolarMass),
        solutes: Vec<(SpeciesId, MolarMass)>,
        mut sources: Vec<BinarySolutionDiffusionSource>,
        axes: SolutionDiffusionSourceAxes,
        reference: &AqueousSolution,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let n = solutes.len() + 1;
        if solutes.len() != thermodynamics.species_count()
            || solutes.len() > MAX_OSMOTIC_SOLUTES
            || sources.len() != n * (n - 1) / 2
            || axes.concentrations.len() != n
            || reference.nonideal_model() != Some(thermodynamics.as_ref())
        {
            return Err(SolutionDiffusionError::Invalid(
                "same original solution model and complete bounded binary/coordinate assignments required",
            ));
        }
        let original: Vec<_> = core::iter::once(&solvent.0)
            .chain(solutes.iter().map(|(id, _)| id))
            .collect();
        if original
            .iter()
            .enumerate()
            .any(|(i, id)| original[..i].contains(id))
            || !valid_axis(&axes.temperature)
            || axes
                .concentrations
                .iter()
                .enumerate()
                .any(|(i, (id, axis))| {
                    !original.contains(&id)
                        || !valid_axis(axis)
                        || axis == &axes.temperature
                        || axes.concentrations[..i]
                            .iter()
                            .any(|(prior, name)| prior == id || name == axis)
                })
        {
            return Err(SolutionDiffusionError::Invalid(
                "distinct original species and complete unique owned source axes required",
            ));
        }
        let concentration_axes: Box<[String]> = original
            .iter()
            .map(|id| {
                axes.concentrations
                    .iter()
                    .find(|(candidate, _)| candidate == *id)
                    .expect("admitted complete coordinate assignments")
                    .1
                    .clone()
            })
            .collect();
        for pair in &mut sources {
            if pair.species[0] > pair.species[1] {
                pair.species.swap(0, 1);
            }
        }
        sources.sort_by(|a, b| a.species.cmp(&b.species));
        if sources.iter().enumerate().any(|(i, pair)| {
            pair.species[0] == pair.species[1]
                || !original.contains(&&pair.species[0])
                || !original.contains(&&pair.species[1])
                || (i > 0 && sources[i - 1].species == pair.species)
                || pair.property.quantity()
                    != QuantitySpec::dimensional(BINARY_SOLUTION_DIFFUSIVITY_DIMS)
                || pair.property.hardness_test().is_some()
                || pair.property.elastic_component().is_some()
                || pair.property.strain_component().is_some()
                || pair.property.stress_component().is_some()
        }) {
            return Err(SolutionDiffusionError::Invalid(
                "one explicit scalar m2/s property for every distinct original species pair required",
            ));
        }
        let c = concentrations(
            &thermodynamics,
            reference.solvent_volume_m3(),
            reference.solute_moles(),
        )?;
        let query = query(
            &thermodynamics,
            &axes.temperature,
            &concentration_axes,
            &c[..n - 1],
        )?;
        let mut matrix = vec![0.0; n * n];
        let mut pinned = Vec::with_capacity(sources.len());
        for source in sources {
            poll(&mut cancelled)?;
            let requirement = ScalarPropertyRequirement::try_with_key(
                &source.property,
                ScalarAdmissibility::StrictlyPositive,
            )?;
            let resolved = resolve_material_state_point(
                &source.material,
                &query,
                &[requirement],
                source.selection,
            )?;
            let property = resolved
                .property_by_key(&source.property)
                .expect("original resolver retains the exact required property");
            let selected = property.answer().receipt.selected;
            let claim = source
                .material
                .claims()
                .claim(selected)
                .expect("selected claim belongs to the original card");
            if claim
                .validity
                .bounds()
                .keys()
                .any(|axis| !query.axis_quantities().contains_key(axis))
                || claim.value.axes().any(|(axis, dims)| {
                    query
                        .axis_quantities()
                        .get(axis)
                        .is_none_or(|quantity| quantity.dims() != dims)
                })
            {
                return Err(SolutionDiffusionError::Invalid(
                    "binary source has an unowned physical axis",
                ));
            }
            let indices = source.species.each_ref().map(|id| {
                original
                    .iter()
                    .position(|candidate| *candidate == id)
                    .expect("admitted original binary species")
            });
            matrix[indices[0] * n + indices[1]] = property.value_si();
            matrix[indices[1] * n + indices[0]] = property.value_si();
            // Rebuild the admitted context-free key so its sole heap payload
            // has the exact retained name length used in the payload charge.
            pinned.push(PinnedBinarySolutionDiffusionSource {
                species: source.species,
                material: source.material,
                property: PropertyKey::new(source.property.name(), source.property.dims()),
                selected,
                original_indices: indices,
            });
        }
        let reference_law = SolutionMaxwellStefan::new(thermodynamics, solvent, solutes, matrix)?;
        reference_law.admit_solution(reference)?;
        poll(&mut cancelled)?;
        Ok(Self {
            reference_law,
            sources: pinned.into_boxed_slice(),
            temperature_axis: axes.temperature,
            concentration_axes,
        })
    }

    /// Original thermodynamic/friction owner with coefficients at the reference
    /// state only. Use this adapter's `flux_into` for evolving material data.
    pub fn reference_law(&self) -> &SolutionMaxwellStefan {
        &self.reference_law
    }
    pub fn sources(&self) -> &[PinnedBinarySolutionDiffusionSource] {
        &self.sources
    }
    pub fn workspace(&self) -> Result<SolutionDiffusionWorkspace> {
        self.reference_law.workspace()
    }

    /// Estimated retained payload excluding this fixed-size owner. Counts the
    /// full original card per pair even when Arcs share it. Includes unselected
    /// source data and the original thermodynamic model; excludes allocator
    /// overhead and transient query/receipt allocations, which are NOT bounded
    /// by this prepared-storage charge. No runtime receipts are retained here.
    pub fn estimated_heap_payload_bytes(&self) -> usize {
        self.reference_law.heap_bytes()
            + core::mem::size_of_val(self.sources.as_ref())
            + core::mem::size_of_val(self.concentration_axes.as_ref())
            + self.temperature_axis.capacity()
            + self
                .concentration_axes
                .iter()
                .map(String::capacity)
                .sum::<usize>()
            + self
                .sources
                .iter()
                .map(|source| {
                    source.material.estimated_retained_payload_bytes()
                        + 2 * core::mem::size_of::<usize>()
                        + source
                            .species
                            .iter()
                            .map(|id| id.as_str().len())
                            .sum::<usize>()
                        + source.property.name().len()
                })
                .sum::<usize>()
    }

    /// Qualify actual endpoint inventories AND their arithmetic concentration
    /// face against every pinned source. Useful for checking a proposed finite
    /// update before publication. Fixed model T and solvent concentration are
    /// read from their original owner; no reference-state coordinate is frozen.
    pub fn resolve_inventories(
        &self,
        volumes_m3: [f64; 2],
        amounts_mol: [&[f64]; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedSolutionDiffusivities> {
        poll(&mut cancelled)?;
        let model = self.reference_law.thermodynamics();
        let n = model.species_count() + 1;
        let a = concentrations(model, volumes_m3[0], amounts_mol[0])?;
        let b = concentrations(model, volumes_m3[1], amounts_mol[1])?;
        let mut face = [0.0; MAX_OSMOTIC_SOLUTES];
        for i in 0..n - 1 {
            face[i] = 0.5 * a[i] + 0.5 * b[i];
        }
        model.evaluate(&face[..n - 1]).map_err(OsmosisError::from)?;
        let queries = [&a, &face, &b].map(|c| {
            query(
                model,
                &self.temperature_axis,
                &self.concentration_axes,
                &c[..n - 1],
            )
        });
        let [left_query, face_query, right_query] = queries;
        let queries = [left_query?, face_query?, right_query?];
        let mut binary = vec![0.0; n * n];
        let mut pairs = Vec::with_capacity(self.sources.len());
        for source in &self.sources {
            poll(&mut cancelled)?;
            let requirement = ScalarPropertyRequirement::try_with_key(
                &source.property,
                ScalarAdmissibility::StrictlyPositive,
            )?;
            let mut at = |query: &QueryPoint| -> Result<ResolvedMaterialStatePoint> {
                poll(&mut cancelled)?;
                Ok(resolve_material_state_point(
                    &source.material,
                    query,
                    core::slice::from_ref(&requirement),
                    MaterialPropertySelection::PinnedByProperty(vec![(
                        source.property.name().into(),
                        source.selected,
                    )]),
                )?)
            };
            let left = at(&queries[0])?;
            let face = at(&queries[1])?;
            let right = at(&queries[2])?;
            let value = face
                .property_by_key(&source.property)
                .expect("original resolver retains the exact required property")
                .value_si();
            let [i, j] = source
                .original_indices
                .map(|i| self.reference_law.original_to_canonical[i]);
            binary[i * n + j] = value;
            binary[j * n + i] = value;
            pairs.push(ResolvedBinarySolutionDiffusivity {
                species: source.species.clone(),
                left,
                face,
                right,
            });
        }
        poll(&mut cancelled)?;
        Ok(ResolvedSolutionDiffusivities {
            canonical_species: self.reference_law.owner.species().into(),
            binary_diffusivities_m2_s: binary.into_boxed_slice(),
            pairs,
        })
    }

    /// Evaluate actual state-dependent pairs and then the SAME shared
    /// solvent-frame physical law. Failure/cancellation preserves `output`.
    /// Successful calls allocate material queries and returned source receipts.
    #[allow(clippy::too_many_arguments)]
    pub fn flux_into(
        &self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        distance_m: f64,
        controls: SolutionDiffusionSolve,
        workspace: &mut SolutionDiffusionWorkspace,
        output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialSolutionDiffusionReport> {
        poll(&mut cancelled)?;
        self.reference_law.admit_solution(left)?;
        self.reference_law.admit_solution(right)?;
        let sources = self.resolve_inventories(
            [left.solvent_volume_m3(), right.solvent_volume_m3()],
            [left.solute_moles(), right.solute_moles()],
            &mut cancelled,
        )?;
        let constitutive = self.reference_law.flux_with_binary_into(
            left,
            right,
            distance_m,
            controls,
            &sources.binary_diffusivities_m2_s,
            workspace,
            output,
            cancelled,
        )?;
        Ok(MaterialSolutionDiffusionReport {
            constitutive,
            sources,
        })
    }
}

fn valid_axis(axis: &str) -> bool {
    !axis.trim().is_empty() && axis.len() <= 128 && !axis.chars().any(char::is_control)
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn concentrations(
    model: &VirialSolutionModel,
    volume: f64,
    amounts: &[f64],
) -> Result<[f64; MAX_OSMOTIC_SOLUTES]> {
    model
        .evaluate_inventories(volume, amounts)
        .map_err(OsmosisError::from)?;
    let mut c = [0.0; MAX_OSMOTIC_SOLUTES];
    for (slot, amount) in c.iter_mut().zip(amounts) {
        *slot = amount / volume;
        if !positive(*slot) {
            return Err(SolutionDiffusionError::Invalid(
                "positive representable source concentrations required",
            ));
        }
    }
    Ok(c)
}
fn query(
    model: &VirialSolutionModel,
    temperature_axis: &str,
    axes: &[String],
    c: &[f64],
) -> Result<QueryPoint> {
    let build = || {
        let mut query = QueryPoint::new().with_quantity(
            temperature_axis,
            typed(QuantityKind::AbsoluteTemperature),
            model.temperature_k(),
        )?;
        for (axis, value) in axes.iter().zip(
            core::iter::once(1.0 / model.solvent_molar_volume_m3_mol()).chain(c.iter().copied()),
        ) {
            query = query.with_quantity(axis, typed(QuantityKind::AmountConcentration), value)?;
        }
        Ok::<_, fs_matdb::MatDbError>(query)
    };
    build().map_err(|source| {
        MaterialStatePointError::Query {
            property: "solution diffusion coordinates".into(),
            source,
        }
        .into()
    })
}
