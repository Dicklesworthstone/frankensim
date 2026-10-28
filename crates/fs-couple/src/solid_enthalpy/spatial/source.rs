//! Source selection, continuous support and native constitutive lowering.

use fs_blake3::{ContentHash, DomainHasher};
use fs_conduction::material::{ConductivityModel, ConductivityTable};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, QueryPoint};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_solid::linear3::TetElasticMaterial;

use super::{
    SpatialMaterialSource, SpatialSolidEnthalpyError as Error, ThermalExpansionSource, poll,
    same_measure,
};

pub(super) struct BoundSource<'a> {
    pub source: SpatialMaterialSource<'a>,
    pub reference_elastic: ResolvedMaterialStatePoint,
    pub reference_density_kg_m3: f64,
    pub conductivity: ConductivityModel,
    reference_phase: EquilibriumPhaseState,
    reference_point: QueryPoint,
}

impl<'a> BoundSource<'a> {
    pub fn try_new(
        mut source: SpatialMaterialSource<'a>,
        native: ReferenceEnthalpyMaterial<'_>,
        mesh_frame: ContentHash,
        maximum_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, Error> {
        poll(cx)?;
        let [low, high] = source.temperature_interval_k;
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || low >= high
            || source.point.temperature_axis != "T"
            || source.curve.identity() != native.curve.identity()
            || source.curve.material_card_identity() != source.point.material.content_hash()
        {
            return Err(Error::Invalid(
                "exact native caloric chart/card, canonical T and positive increasing source interval required",
            ));
        }
        let reference_phase = source
            .curve
            .state_at_specific_enthalpy(source.reference_specific_enthalpy_j_kg)?;
        if !reference_phase.is_fully_solid()
            || reference_phase.temperature_k() < low
            || reference_phase.temperature_k() > high
        {
            return Err(Error::Invalid(
                "stress-free reference must be solid and inside source interval",
            ));
        }
        let reference_point = source.point.at(reference_phase)?;
        let elastic = source.stiffness_model.resolve_reference(
            source.point.material,
            &reference_point,
            source.mechanical_selection.clone(),
        )?;
        for property in elastic.properties() {
            source.point.admit_static(property)?;
        }
        let density = elastic
            .property(DENSITY_PROPERTY)
            .ok_or(Error::Invalid("reference density source required"))?
            .value_si();
        if !same_measure(density, native.reference_density_kg_m3)
            || !same_measure(density, reference_phase.bulk_density_kg_m3())
        {
            return Err(Error::Invalid(
                "original caloric, prepared heat and mechanical reference densities must agree",
            ));
        }
        let material = source.stiffness_model.material(
            &elastic,
            native.reference_density_kg_m3,
            elastic.identity(),
            maximum_strain,
        )?;
        let expansion = source.expansion_model.resolve(
            source.point.material,
            &reference_point,
            &reference_point,
            source.expansion_selection.clone(),
            &material,
            mesh_frame,
            maximum_strain,
        )?;
        for property in expansion.current().properties() {
            source.point.admit_static(property)?;
        }
        source.mechanical_selection = MaterialPropertySelection::PinnedByProperty(
            elastic
                .properties()
                .iter()
                .filter(|p| p.requirement().name() != DENSITY_PROPERTY)
                .map(|p| {
                    (
                        p.requirement().name().to_owned(),
                        p.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        source.expansion_selection = MaterialPropertySelection::PinnedByProperty(
            expansion
                .reference()
                .properties()
                .iter()
                .map(|p| {
                    (
                        p.requirement().name().to_owned(),
                        p.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let lower = at_temperature(&reference_point, low)?;
        let upper = at_temperature(&reference_point, high)?;
        let queries = [lower.clone(), upper.clone()];
        let table = match source.conductivity_selection {
            ClaimSelection::Policy(policy) => ConductivityTable::from_claims_at_query_points(
                source.point.material.claims(),
                &source.conductivity,
                "T",
                &queries,
                policy,
            )?,
            ClaimSelection::Pinned(pin) => ConductivityTable::from_claims_pinned_at_query_points(
                source.point.material.claims(),
                &source.conductivity,
                "T",
                &queries,
                pin,
            )?,
        };
        // The conductivity table proves continuous support and retains every
        // knot. Admit its exact selected source axes as well: a fixed-looking
        // query must not silently freeze an unsupported dynamic coordinate.
        let conductivity_pin = table
            .receipts()
            .first()
            .ok_or(Error::Invalid("source conductivity receipts required"))?
            .selected;
        let requirements = [ScalarPropertyRequirement::try_with_key(
            &source.conductivity,
            ScalarAdmissibility::StrictlyPositive,
        )?];
        let conductivity = resolve_material_state_point(
            source.point.material,
            &reference_point,
            &requirements,
            MaterialPropertySelection::PinnedByProperty(vec![(
                source.conductivity.name().to_owned(),
                conductivity_pin,
            )]),
        )?;
        for property in conductivity.properties() {
            source.point.admit_static(property)?;
        }
        let mut grid = vec![low, high, reference_phase.temperature_k()];
        for property in elastic
            .properties()
            .iter()
            .filter(|p| p.requirement().name() != DENSITY_PROPERTY)
            .chain(expansion.reference().properties())
        {
            poll(cx)?;
            let selected = property.answer().receipt.selected;
            let claim = source
                .point
                .material
                .claims()
                .claim(selected)
                .ok_or(Error::Invalid("pinned source claim is absent"))?;
            let support = source
                .point
                .material
                .claims()
                .query_line_typed(
                    &claim.key,
                    "T",
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|error| Error::Support {
                    property: claim.key.name().to_owned(),
                    reason: error.to_string(),
                })?;
            if grid
                .len()
                .checked_add(support.knots.len())
                .is_none_or(|n| n > 4096)
            {
                return Err(Error::Invalid(
                    "continuous mechanical source grid exceeds 4096 knots",
                ));
            }
            grid.extend(support.knots.into_iter().map(|(t, _)| t));
        }
        grid.sort_by(f64::total_cmp);
        grid.dedup();
        let bound = Self {
            source,
            reference_elastic: elastic,
            reference_density_kg_m3: native.reference_density_kg_m3,
            conductivity: ConductivityModel::isotropic(table),
            reference_phase,
            reference_point,
        };
        // Validate every original source knot, not just the two endpoint
        // values. This is a source-support admission, not a substitute for
        // native constitutive checks at each later corner and quadrature point.
        for t in grid {
            poll(cx)?;
            let point = at_temperature(&bound.reference_point, t)?;
            bound.resolve_point(&point, mesh_frame, maximum_strain)?;
        }
        poll(cx)?;
        Ok(bound)
    }

    pub fn contains(&self, t: f64) -> bool {
        let [low, high] = self.source.temperature_interval_k;
        t.is_finite() && t >= low && t <= high
    }

    pub fn resolve(
        &self,
        phase: EquilibriumPhaseState,
        frame: ContentHash,
        maximum_strain: f64,
    ) -> Result<
        (
            ResolvedMaterialStatePoint,
            ThermalExpansionSource,
            TetElasticMaterial,
        ),
        Error,
    > {
        if !phase.is_fully_solid() || !self.contains(phase.temperature_k()) {
            return Err(Error::Invalid(
                "material state must remain solid inside its admitted source interval",
            ));
        }
        let point = self.source.point.at(phase)?;
        self.resolve_point(&point, frame, maximum_strain)
    }

    fn resolve_point(
        &self,
        point: &QueryPoint,
        frame: ContentHash,
        maximum_strain: f64,
    ) -> Result<
        (
            ResolvedMaterialStatePoint,
            ThermalExpansionSource,
            TetElasticMaterial,
        ),
        Error,
    > {
        let elastic = self.source.stiffness_model.resolve(
            self.source.point.material,
            point,
            self.source.mechanical_selection.clone(),
        )?;
        for property in elastic.properties() {
            self.source.point.admit_static(property)?;
        }
        let mut identity = DomainHasher::new("org.frankensim.fs-couple.spatial-solid-material.v1");
        identity.update(elastic.identity().as_bytes());
        identity.update(self.reference_elastic.identity().as_bytes());
        identity.update(self.reference_phase.identity().as_bytes());
        identity.update(&self.reference_density_kg_m3.to_bits().to_le_bytes());
        identity.update(frame.as_bytes());
        let material = self.source.stiffness_model.material(
            &elastic,
            self.reference_density_kg_m3,
            identity.finalize(),
            maximum_strain,
        )?;
        let expansion = self.source.expansion_model.resolve(
            self.source.point.material,
            &self.reference_point,
            point,
            self.source.expansion_selection.clone(),
            &material,
            frame,
            maximum_strain,
        )?;
        for property in expansion.current().properties() {
            self.source.point.admit_static(property)?;
        }
        Ok((elastic, expansion, material))
    }
}

fn at_temperature(point: &QueryPoint, t: f64) -> Result<QueryPoint, Error> {
    point.clone().with("T", t).map_err(|error| Error::Support {
        property: "temperature".into(),
        reason: error.to_string(),
    })
}
