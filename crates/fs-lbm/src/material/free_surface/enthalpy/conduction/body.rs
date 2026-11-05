//! Finite, source-bound solid thermal cells on the liquid's existing timeline.
//!
//! A contact names an actual permanent native wall face and its physical wetted
//! area. Solid storage is finite; equal and opposite pair heat is handled by
//! the same thermal kernel as liquid conduction. Geometry stays fixed and a
//! source/phase refusal keeps BOTH solid and liquid histories unchanged.

use super::*;
use fs_matdb::PropertyKey;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_material::thermal::source::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ThermalMaterialPoint,
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// One stationary, uniform-temperature solid cell with retained finite storage.
///
/// Its original caloric chart supplies mass density, temperature and storage;
/// its material supplies scalar k(T). The initial selected conductivity claim
/// is pinned for the cell's lifetime. A different material is new DATA supplied
/// through this constructor, never a dispatch by material or object name.
/// Spatial resolution and uniform-cell/Biot validity remain caller obligations.
#[derive(Clone, Debug)]
pub struct FreeSurfaceThermalBody<'m> {
    id: u64,
    source: ThermalMaterialPoint<'m>,
    caloric: &'m EquilibriumEnthalpyPhaseCurve,
    geometry: CellGeometry,
    mass_kg: f64,
    specific_h_j_kg: f64,
    conductivity_key: PropertyKey,
    requirements: Vec<ScalarPropertyRequirement>,
    selection: MaterialPropertySelection,
    material: ResolvedMaterialStatePoint,
    conduction_links: Vec<FreeSurfaceBodyLink>,
    gas_faces: Vec<atmosphere::FreeSurfaceBodyGasFace>,
    gas_convection: Option<atmosphere::PreparedBodyGasConvection<'m>>,
    radiation_faces: Vec<radiation::FreeSurfaceBodyRadiationFace>,
}

/// One declared solid/solid face, owned by the cell on its first side.
///
/// Declare each unordered body pair ONCE, on either participant. The two
/// distances follow owner/neighbor order. Both cells retain finite storage;
/// scalar half-cell resistances are in series with zero contact resistance.
/// Actual face area and centroid geometry are explicit caller data. The
/// complete timeline checks both surface-coverage budgets before advancing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceBodyLink {
    /// Stable identity of the other finite solid cell in the same step.
    pub other_body_id: u64,
    pub area_m2: f64,
    /// Positive owner and neighbor centroid-to-face distances [m].
    pub distances_m: [f64; 2],
}

/// Source data actually used by one finite solid/solid conduction link.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceBodyConductionSample {
    pub body_id: u64,
    pub link: FreeSurfaceBodyLink,
    /// Original pinned conductivity sources, in owner/neighbor order.
    pub materials: [ResolvedMaterialStatePoint; 2],
    pub conductance_w_k: f64,
}

/// One explicitly resolved wetted patch on a native stationary wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceBodyContact {
    /// Stable finite solid cell identity.
    pub body_id: u64,
    /// Occupied liquid cell at the beginning of this flow step.
    pub liquid: [usize; 3],
    /// Face-adjacent permanent wall cell. Signed coordinates allow the six
    /// outer domain walls; inactive gas tiles are not admitted as solids.
    pub wall: [i64; 3],
    /// Actual wetted contact area [m2], no greater than the Cartesian face.
    pub area_m2: f64,
    /// Liquid and solid thermal-centroid distances to the contact [m].
    /// The liquid distance is at most dx; the solid may be a coarser cell.
    pub distances_m: [f64; 2],
}

/// Original source evaluations used by a solid/liquid contact.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceBodyContactSample {
    /// Physical patch actually used for this evaluation.
    pub contact: FreeSurfaceBodyContact,
    /// Original liquid density/viscosity/conductivity tuple at its current h.
    pub liquid: ResolvedLiquidState,
    /// Source and effective liquid conductivity; the solid keeps its own law.
    pub liquid_conductivity: FreeSurfaceConductivity,
    /// Original pinned solid conductivity source at its current h.
    pub body: ResolvedMaterialStatePoint,
    /// Original ordered interface at the SAME two caloric states. None is the
    /// explicitly selected zero-contact-resistance approximation.
    pub interface: Option<fs_material::state_point::ResolvedInterfaceStatePoint>,
    /// Total series conductance of the bulk layers and any interface [W/K].
    pub conductance_w_k: f64,
}

/// Finite solid storage and the heat exchanged inside the combined system.
#[derive(Clone, Debug)]
pub struct FreeSurfaceBodyExchange {
    /// Sum of original finite solid enthalpies before this interval [J].
    pub enthalpy_before_j: f64,
    /// Sum of accepted finite solid enthalpies after this interval [J].
    pub enthalpy_after_j: f64,
    /// Signed energy from solids INTO liquid. Never external supply or loss.
    pub heat_to_liquid_j: f64,
    /// Direct signed heat from the SAME finite atmosphere into these solids.
    /// This is distinct from `heat_to_liquid_j`; never debit the gas twice.
    pub heat_from_gas_j: f64,
    /// External radiant heat deposited in the solid, before any contact transfer.
    pub heat_from_radiation_j: f64,
    /// `H_liquid_delta + H_solid_delta - liquid_external_heat - viscous_heat
    /// - heat_from_gas - heat_from_radiation` [J]. Gas exchange is internal
    /// when gas is included; prescribed radiant surroundings remain external.
    pub combined_balance_defect_j: f64,
    /// Source/conductance samples used by the final thermal substep.
    pub last_applied: Vec<FreeSurfaceBodyContactSample>,
    /// Supported source/conductance samples at the thermal endpoint.
    pub thermal_endpoint: Vec<FreeSurfaceBodyContactSample>,
    /// Finite solid/solid links used by the last thermal substep. Their heat is
    /// internal to solid storage and does not enter `heat_to_liquid_j`.
    pub solid_last_applied: Vec<FreeSurfaceBodyConductionSample>,
    /// The same links reevaluated at the accepted thermal endpoint.
    pub solid_thermal_endpoint: Vec<FreeSurfaceBodyConductionSample>,
}

impl<'m> FreeSurfaceThermalBody<'m> {
    /// Bind finite solid storage and conductivity to the same material card.
    /// The chart must have constant density; a mixed/liquid current state
    /// refuses because no moving-solid or phase-boundary handoff is supplied.
    /// Later source or phase-domain crossings likewise refuse atomically.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: u64,
        source: ThermalMaterialPoint<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        geometry: CellGeometry,
        specific_h_j_kg: f64,
        conductivity: &PropertyKey,
        selection: MaterialPropertySelection,
    ) -> Result<Self, Error> {
        let quantity = conductivity.quantity();
        if quantity != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && quantity
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                ))
        {
            return Err(Error::Invalid(
                "solid contact needs scalar static thermal conductivity",
            ));
        }
        let phase = caloric.state_at_specific_enthalpy(specific_h_j_kg)?;
        let density = phase.bulk_density_kg_m3();
        if !phase.is_fully_solid()
            || caloric
                .knots()
                .iter()
                .any(|k| k.bulk_density_kg_m3.to_bits() != density.to_bits())
        {
            return Err(Error::Invalid(
                "finite wall needs solid phase and fixed caloric density",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_with_key(
            conductivity,
            ScalarAdmissibility::NonNegative,
        )?;
        let point = source.at(phase)?;
        let material = resolve_material_state_point(
            source.material,
            &point,
            std::slice::from_ref(&requirement),
            selection,
        )?;
        source.admit_static(&material.properties()[0])?;
        let selection = MaterialPropertySelection::PinnedByProperty(vec![(
            requirement.name().to_owned(),
            material.properties()[0].answer().receipt.selected,
        )]);
        let body = Self {
            id,
            source,
            caloric,
            geometry,
            mass_kg: density * geometry.volume_m3,
            specific_h_j_kg,
            conductivity_key: conductivity.clone(),
            requirements: vec![requirement],
            selection,
            material,
            conduction_links: Vec::new(),
            gas_faces: Vec::new(),
            gas_convection: None,
            radiation_faces: Vec::new(),
        };
        // The common owner admits geometry, extensive storage and the actual
        // borrowed caloric law; there is no separate wall thermal integrator.
        DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(body.cell(CellId(0), MaterialId(0))),
            }],
            &[(MaterialId(0), caloric as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 1,
                max_links: 0,
                max_events: 1,
                max_substeps: 1,
            },
        )?
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 0.0,
            },
            || Ok(()),
        )?;
        Ok(body)
    }

    /// Stable identity used by solid/liquid contact data.
    pub const fn id(&self) -> u64 {
        self.id
    }
    /// Attach explicit solid/solid geometry without changing accepted storage.
    ///
    /// Each pair must be declared only once, on either of its two bodies.
    /// Every neighbor must be supplied alongside this body on EACH subsequent
    /// thermal step. Full graph, geometry, surface coverage and work-budget
    /// admission happens on that common timeline; invalid data cannot advance
    /// liquid, solid or gas history. Empty links retain isolated-cell behavior.
    /// All existing finite-wall and surface-exchange entry points use these
    /// links, including finite-gas and thermocapillary coupling.
    #[must_use]
    pub fn with_conduction_links(mut self, links: Vec<FreeSurfaceBodyLink>) -> Self {
        self.conduction_links = links;
        self
    }

    /// Explicit outgoing geometry, with the owning cell on the first side.
    pub fn conduction_links(&self) -> &[FreeSurfaceBodyLink] {
        &self.conduction_links
    }

    /// Attach exposed patches to the same finite atmosphere as the liquid.
    ///
    /// Every subsequent thermal step must supply that finite atmosphere.
    /// Current native wall/gas geometry and the combined solid contact/film
    /// area budget are admitted on each step, before any participant advances.
    /// No independent bath temperature or gas transport coefficient is stored.
    #[must_use]
    pub fn with_gas_faces(mut self, faces: Vec<atmosphere::FreeSurfaceBodyGasFace>) -> Self {
        self.gas_faces = faces;
        self
    }

    pub fn gas_faces(&self) -> &[atmosphere::FreeSurfaceBodyGasFace] {
        &self.gas_faces
    }

    /// Resolve the original solid k(T) line once and apply gas convection at
    /// the conducting surface, in series with each solid centroid half-path.
    ///
    /// The original selected conductivity claim and pure-solid caloric branch
    /// remain pinned. Every participating face must have positive area and solid
    /// distance, with ZERO gas thickness: the provider replaces that film.
    /// The shared finite gas owns bulk temperature, EOS pressure and all heat.
    /// Body clones and thermal trials share this prepared source line.
    ///
    /// Radiation on the same native face must belong to this body and have the
    /// exact same area. Both then use one solved surface and original emissivity
    /// at that surface; radiation on other faces retains its cell-temperature law.
    pub fn with_gas_convection(
        mut self,
        provider: &'m dyn atmosphere::FreeSurfaceBodyGasConvection,
        temperature_range_k: [f64; 2],
        config: fs_material::thermal::conduction::LayerConductionConfig,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        use fs_material::thermal::conduction::source::ResolvedConductionLayer;
        if self.gas_convection.is_some() {
            return Err(Error::Invalid("finite wall convection is already source-bound"));
        }
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(Error::Invalid("invalid finite wall convection solve controls"));
        }
        let selected = self.material.property_by_key(&self.conductivity_key)
            .ok_or(Error::Invalid("finite wall lost its original conductivity key"))?
            .answer().receipt.selected;
        let layer = ResolvedConductionLayer::new(
            self.caloric,
            self.source,
            &self.conductivity_key,
            fs_matdb::ClaimSelection::Pinned(selected),
            fs_material::phase::SolidLiquidPhase::Solid,
            temperature_range_k,
            || checkpoint(gate),
        )?;
        layer.temperature_at_enthalpy(self.specific_h_j_kg)?;
        self.gas_convection = Some(atmosphere::PreparedBodyGasConvection {
            provider,
            layer: std::sync::Arc::new(layer),
            config,
        });
        Ok(self)
    }

    pub(super) fn gas_convection(&self) -> Option<&atmosphere::PreparedBodyGasConvection<'m>> {
        self.gas_convection.as_ref()
    }

    /// Attach exposed diffuse-gray patches and select emissivity from this
    /// body's ORIGINAL material at its actual accepted caloric state.
    ///
    /// The canonical total hemispherical emissivity joins the conductivity
    /// bundle. Both selected claims are pinned for every later substep and
    /// endpoint, retaining all original nonthermal query coordinates. A known
    /// zero is a reflecting surface; absent data never becomes a guessed value.
    /// This one-time attachment cannot silently replace an already pinned law.
    /// Geometry and aggregate exposed-area/work budgets are checked by the
    /// shared timeline. An equal-area convection patch on the same body/face
    /// shares its solved surface; other patches use the original cell temperature.
    pub fn with_radiation_faces(
        mut self,
        faces: Vec<radiation::FreeSurfaceBodyRadiationFace>,
        selection: MaterialPropertySelection,
    ) -> Result<Self, Error> {
        if self
            .requirements
            .iter()
            .any(|r| r.name() == SURFACE_EMISSIVITY_PROPERTY)
        {
            return Err(Error::Invalid(
                "finite wall emissivity is already source-bound",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_new(
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        let point = self.source.at(self
            .caloric
            .state_at_specific_enthalpy(self.specific_h_j_kg)?)?;
        let emissivity = resolve_material_state_point(
            self.source.material,
            &point,
            std::slice::from_ref(&requirement),
            selection,
        )?;
        let mut pinned: Vec<_> = self
            .material
            .properties()
            .iter()
            .map(|p| {
                (
                    p.requirement().name().to_owned(),
                    p.answer().receipt.selected,
                )
            })
            .collect();
        pinned.push((
            SURFACE_EMISSIVITY_PROPERTY.into(),
            emissivity.properties()[0].answer().receipt.selected,
        ));
        self.requirements.push(requirement);
        self.selection = MaterialPropertySelection::PinnedByProperty(pinned);
        self.material = self.resolve(self.specific_h_j_kg)?;
        self.radiation_faces = faces;
        Ok(self)
    }

    /// Explicit exposed geometry on the body's retained emissivity law.
    pub fn radiation_faces(&self) -> &[radiation::FreeSurfaceBodyRadiationFace] {
        &self.radiation_faces
    }

    /// Fixed mass from the original chart density and cell volume [kg].
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Accepted caloric state [J/kg], without an independent temperature cache.
    pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.specific_h_j_kg
    }
    /// Accepted extensive storage, mass times specific enthalpy [J].
    pub fn enthalpy_j(&self) -> f64 {
        self.mass_kg * self.specific_h_j_kg
    }
    /// Original conductivity and any attached emissivity at accepted temperature.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub(super) const fn caloric(&self) -> &'m EquilibriumEnthalpyPhaseCurve {
        self.caloric
    }
    pub(super) fn accept(&mut self, h: f64, material: ResolvedMaterialStatePoint) {
        self.specific_h_j_kg = h;
        self.material = material;
    }

    pub(super) fn cell(&self, id: CellId, material_id: MaterialId) -> DepositedCell {
        DepositedCell {
            id,
            segment_id: self.id,
            body_id: self.id,
            material_id,
            mass_kg: self.mass_kg,
            specific_h_j_kg: self.specific_h_j_kg,
            geometry: self.geometry,
            boundary: Boundary {
                ambient_k: self
                    .caloric
                    .state_at_specific_enthalpy(self.specific_h_j_kg)
                    .expect("admitted thermal body")
                    .temperature_k(),
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        }
    }

    pub(super) fn resolve(&self, h: f64) -> Result<ResolvedMaterialStatePoint, Error> {
        let phase = self.caloric.state_at_specific_enthalpy(h)?;
        if !phase.is_fully_solid() {
            return Err(Error::Invalid(
                "finite wall melted without a mechanical phase handoff",
            ));
        }
        let point = self.source.at(phase)?;
        let material = resolve_material_state_point(
            self.source.material,
            &point,
            &self.requirements,
            self.selection.clone(),
        )?;
        for property in material.properties() {
            self.source.admit_static(property)?;
        }
        if let Some(emissivity) = material.property(SURFACE_EMISSIVITY_PROPERTY) {
            radiation::admit_emissivity(emissivity.value_si())?;
        }
        Ok(material)
    }

    pub(super) fn conductivity(&self, material: &ResolvedMaterialStatePoint) -> f64 {
        material
            .property(self.requirements[0].name())
            .expect("required conductivity")
            .value_si()
    }

    /// Query only the original pinned emissivity at a surface already admitted
    /// by the prepared pure-solid conductivity span. Do not rebuild k(T) or
    /// replace the cell's accepted caloric/material state during root trials.
    pub(super) fn emissivity_at_surface(
        &self,
        specific_h_j_kg: f64,
        temperature_k: f64,
    ) -> Result<ResolvedScalarProperty, Error> {
        let phase = self.caloric.state_at_specific_enthalpy(specific_h_j_kg)?;
        if !phase.is_fully_solid() {
            return Err(Error::Invalid("conducting radiative wall left its solid phase"));
        }
        let requirement = self.requirements.iter()
            .find(|r| r.name() == SURFACE_EMISSIVITY_PROPERTY)
            .ok_or(Error::Invalid("conducting radiative wall has no emissivity source"))?;
        let selected = self.material.property(SURFACE_EMISSIVITY_PROPERTY)
            .ok_or(Error::Invalid("conducting radiative wall lost its emissivity claim"))?
            .answer().receipt.selected;
        let point = self.source.at(phase)?
            .with(self.source.temperature_axis, temperature_k)
            .map_err(|_| Error::Invalid("invalid conducting wall surface temperature"))?;
        let material = resolve_material_state_point(
            self.source.material,
            &point,
            std::slice::from_ref(requirement),
            MaterialPropertySelection::PinnedByProperty(vec![(
                requirement.name().to_owned(), selected,
            )]),
        )?;
        let emissivity = material.property(SURFACE_EMISSIVITY_PROPERTY)
            .ok_or(Error::Invalid("conducting radiative wall lost its surface emissivity"))?;
        self.source.admit_static(emissivity)?;
        radiation::admit_emissivity(emissivity.value_si())?;
        Ok(emissivity.clone())
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Evolve finite solid cells, liquid conduction and actual VOF motion on
    /// one accepted interval. Contacts must be supplied from current wetted
    /// geometry at each call. A wall cell can serve several disjoint patches;
    /// each Cartesian liquid-wall face appears once. Contact resistance is
    /// explicitly neglected; material half-cell resistances are in series.
    ///
    /// All source/conduction checks and the native flow candidate must succeed
    /// before either participant changes. Pair heat is internal. Attached
    /// radiant wall patches supply separately reported external heat on this
    /// same timeline. Fixed solid
    /// geometry, pure liquid, neglected capillarity and the original first-order
    /// thermal/flow split remain. This is not a finite-strain or melting handoff.
    pub fn step_with_body_conduction(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        self.advance_with_thermal_bodies(faces, heat, controls, &[], bodies, contacts, &[], gate)
    }

    /// The same finite-body exchange with the original molecular heating owner.
    /// Only physical molecular loss is deposited; numerical loss is not heat.
    pub fn step_with_body_conduction_and_viscous_heating(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        let applied = self.viscous_sources(gate)?;
        let conduction = self.advance_with_thermal_bodies(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            &[],
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}
