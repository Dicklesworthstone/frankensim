//! Data-only bindings from the SAME caloric cards to circuit and heat transport.
//! This module owns no material or stepping law. It queries the existing
//! ThermalMaterialPoint and MaterialOhmicBody, then lowers scalar or explicitly
//! oriented principal conduction across declared uniform half-cell distances.

mod prepared;
pub mod conductivity;
mod directional;
pub use conductivity::CellConductionResponse;
pub use directional::{DirectionalConductionLink, DirectionalConductionSide};
pub use prepared::{
    ConductorFreeExpansion, ExpandingBoundary, ExpandingBoundarySample,
    ExpandingMaterialNetwork, ExpandingNetworkSample, MaterialTransportObserver,
    PreparedMaterialNetwork,
};
use conductivity::ElectricalAxes;

use super::*;
use crate::cavity_phs::piston::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::cavity_phs::piston::thermal::power::joule::{JouleHeatingError, MaterialOhmicBody};
use crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point};
use fs_matdb::PropertyKey;
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_thermal::enthalpy::{MaterialId, ThermalLaw};

/// A physical thermal cell and its original immutable caloric/material binding.
pub struct MaterialCell<'a> {
    /// Active thermal cell.
    pub cell: CellId,
    /// Exact law-table entry used by that cell.
    pub material_id: MaterialId,
    /// The SAME borrowed curve instance used to initialize the thermal timeline.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    /// Card and physical coordinates; T and declared phase fraction evolve.
    pub point: ThermalMaterialPoint<'a>,
    /// Scalar conductivity key and explicit claim selection. A linked side needs
    /// this binding or an explicit prepared directional binding. None does not
    /// imply insulation, nor an inferred material conductivity.
    pub conductivity: Option<(&'a PropertyKey, MaterialPropertySelection)>,
}
/// One whole-cell uniform electrical branch; no independent ideal drive.
pub struct CellResistor {
    /// Original circuit branch and its sole receiving cell.
    pub destination: ResistorCell,
    /// Original electrical source selection; legacy methods require resistivity.
    pub selection: MaterialPropertySelection,
    /// Fixed current path length [m].
    pub length_m: f64,
    /// Fixed current-carrying area [m2]; their product must equal cell volume.
    pub area_m2: f64,
}
/// One existing thermal link with explicit homogeneous layer lengths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductionLink {
    /// Existing stable thermal link ID; its endpoints/area stay authoritative.
    pub link: u64,
    /// Positive distance [m] through the link's `a` cell.
    pub distance_a_m: f64,
    /// Positive distance [m] through the link's `b` cell.
    pub distance_b_m: f64,
}
/// Complete source binding, not a parallel electrical or thermal state owner.
pub struct MaterialNetwork<'a> {
    /// Exactly one binding per active cell, including passive thermal cells.
    pub cells: &'a [MaterialCell<'a>],
    /// One whole-cell binding per circuit resistor. Distinct resistor cells.
    pub resistors: &'a [CellResistor],
    /// Exactly one material-conduction binding per active thermal link, except
    /// links explicitly supplied to the additional-conductance methods.
    /// Those methods require complete, disjoint ownership at every trial.
    pub links: &'a [ConductionLink],
}
/// Original failures remain distinguishable from missing bindings.
#[derive(Debug)]
pub enum MaterialNetworkError {
    /// Geometry, source/law identity, missing or duplicate binding.
    Invalid(&'static str),
    /// Original enthalpy curve refusal.
    Phase(PhaseStateError),
    /// Original source-point or static-property refusal.
    Material(ThermalMaterialError),
    /// Original Ohmic geometry/source refusal.
    Conductor(JouleHeatingError),
    /// Original prepared expansion/surface geometry refusal.
    Expansion(fs_material::thermal::transport::ThermalTransportError),
    /// Original thermal state or cancellation refusal.
    Thermal(DepositionError),
}
impl From<DepositionError> for MaterialNetworkError { fn from(e: DepositionError) -> Self { Self::Thermal(e) } }
impl core::fmt::Display for MaterialNetworkError {
    fn fmt(&self, f:&mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f,"material circuit network: {self:?}") }
}
impl std::error::Error for MaterialNetworkError {}
/// Exact accepted cell phase and optional original transport source bundle.
#[derive(Clone, Debug, PartialEq)]
pub struct CellMaterialSample {
    /// Stable receiving cell ID.
    pub cell: CellId,
    /// Actual endpoint enthalpy/phase, bound to its original chart/card.
    pub phase: EquilibriumPhaseState,
    /// Complete same-point transport bundle for a conductive thermal cell,
    /// including its electrical property when it also binds a resistor.
    pub conductivity: Option<ResolvedMaterialStatePoint>,
}
/// Actual link law used by the shared thermal residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinkMaterialSample {
    /// Original thermal link ID.
    pub link: u64,
    /// Endpoints in original link order, not sorted circuit branch order.
    pub cells: [CellId;2],
    /// Actual scalar or n^T K n pair [W/(m K)] at those endpoint phase states.
    pub conductivity_w_m_k: [f64;2],
    /// Resulting total conductance [W/K].
    pub conductance_w_k: f64,
}
/// Exact complete material evaluation attached to the accepted physical interval.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialNetworkSample {
    /// Stable CellId ordering, independent of input binding order.
    pub cells: Vec<CellMaterialSample>,
    /// Stable thermal-link ordering.
    pub links: Vec<LinkMaterialSample>,
}

impl MaterialNetwork<'_> {
    /// Resolve all endpoint material properties without changing any history.
    /// Source cards, phase charts, law-table entries and the actual borrowed
    /// caloric objects must agree. Pointer identity is an in-process binding
    /// check, not measured-material evidence or portable content authentication.
    /// Both linked conductivities are queried even at zero transfer area. Only
    /// explicit zero conductivity/area means isolation. Missing data refuses.
    pub fn resolve(&self, thermal:&DepositionSimulation<'_>, mut cancelled:impl FnMut()->bool)
        -> Result<CircuitThermalCoefficients<MaterialNetworkSample>,MaterialNetworkError> {
        if cancelled() {return Err(MaterialNetworkError::Thermal(DepositionError::Cancelled));}
        if self.resistors.len()>MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(MaterialNetworkError::Invalid("material resistor bindings exceed circuit capacity"));
        }
        let bases=self.resistors.iter().map(|r|(r.destination.branch,ElectricalPropertyBasis::Resistivity)).collect::<Vec<_>>();
        self.resolve_with_bases(thermal,&bases,&mut cancelled)
    }

    /// Resolve an explicit basis for every resistor branch. Thermal conductivity
    /// and electrical response share one original material bundle per cell.
    /// Missing/duplicate branch choices and missing selected sources refuse.
    pub fn resolve_with_bases(&self, thermal:&DepositionSimulation<'_>,
        electrical_bases:&[(usize,ElectricalPropertyBasis)], mut cancelled:impl FnMut()->bool)
        -> Result<CircuitThermalCoefficients<MaterialNetworkSample>,MaterialNetworkError> {
        if cancelled() {return Err(MaterialNetworkError::Thermal(DepositionError::Cancelled));}
        let bases=self.basis_plan(electrical_bases)?;
        self.resolve_bound(thermal,&bases,None,None,None,cancelled)
    }

    fn basis_plan(&self, electrical_bases:&[(usize,ElectricalPropertyBasis)])
        -> Result<BTreeMap<usize,ElectricalPropertyBasis>,MaterialNetworkError> {
        if self.resistors.len()>MAX_MATERIAL_CIRCUIT_BRANCHES
            || electrical_bases.len()!=self.resistors.len() {
            return Err(MaterialNetworkError::Invalid("complete bounded electrical basis plan required"));
        }
        let bases:BTreeMap<_,_>=electrical_bases.iter().copied().collect();
        if bases.len()!=electrical_bases.len()
            || self.resistors.iter().any(|r|!bases.contains_key(&r.destination.branch)) {
            return Err(MaterialNetworkError::Invalid("exactly one explicit electrical basis per resistor is required"));
        }
        Ok(bases)
    }

    fn resolve_bound(
        &self,
        thermal: &DepositionSimulation<'_>,
        bases: &BTreeMap<usize, ElectricalPropertyBasis>,
        pinned: Option<&BTreeMap<CellId, MaterialPropertySelection>>,
        electrical_axes: Option<&ElectricalAxes>,
        directional: Option<&directional::DirectionalPlan>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        self.resolve_bound_with_links(thermal, bases, pinned, electrical_axes, None, directional, cancelled)
    }

    // Additional physical link owners supply their actual step-local G here.
    // They do not replace Fourier bindings or silently inherit authored G.
    fn resolve_bound_with_links(
        &self,
        thermal: &DepositionSimulation<'_>,
        bases: &BTreeMap<usize, ElectricalPropertyBasis>,
        pinned: Option<&BTreeMap<CellId, MaterialPropertySelection>>,
        electrical_axes: Option<&ElectricalAxes>,
        additional_links: Option<&BTreeMap<u64, f64>>,
        directional: Option<&directional::DirectionalPlan>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        let invalid=MaterialNetworkError::Invalid;
        let mut check=|| if cancelled() {Err(DepositionError::Cancelled)} else {Ok(())};
        check()?;
        if let Some(directional) = directional { directional.admit_topology(thermal, &mut check)?; }
        if self.cells.len()!=thermal.cells().len()
            || self.links.len().checked_add(additional_links.map_or(0, BTreeMap::len))
                != Some(thermal.links().len())
            || self.resistors.len()>MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(invalid("complete bounded cell/link bindings required"));
        }
        if let Some(additional) = additional_links {
            let owned = self.links.iter().map(|link| link.link)
                .collect::<std::collections::BTreeSet<_>>();
            for (&id, &conductance) in additional {
                check()?;
                let link = thermal.links().get(&id)
                    .ok_or(invalid("additional thermal link is inactive"))?;
                if owned.contains(&id) {
                    return Err(invalid("additional thermal links cannot replace material conduction"));
                }
                if link.conductance_w_k != 0.0
                    || !conductance.is_finite() || conductance < 0.0 {
                    return Err(invalid(
                        "additional thermal links need zero authored and finite nonnegative resolved conductance",
                    ));
                }
            }
        }
        let mut resistor_cells=BTreeMap::new(); let mut seen_branches=BTreeMap::new();
        for resistor in self.resistors {
            if seen_branches.insert(resistor.destination.branch,()).is_some()
                || resistor_cells.insert(resistor.destination.cell,resistor).is_some() {
                return Err(invalid("duplicate whole-cell resistor or branch binding"));
            }
        }
        let mut bindings=BTreeMap::new(); let mut samples=BTreeMap::new();
        let mut materials=BTreeMap::new();
        for binding in self.cells {
            check()?;
            let cell=thermal.cells().get(&binding.cell).ok_or(invalid("inactive material cell"))?;
            if cell.material_id!=binding.material_id || bindings.insert(binding.cell,binding).is_some()
                || !thermal.uses_material_law(binding.material_id,binding.curve as &dyn ThermalLaw) {
                return Err(invalid("cell binding must retain the actual thermal law instance and material ID"));
            }
            let phase=binding.curve.state_at_specific_enthalpy(cell.specific_h_j_kg).map_err(MaterialNetworkError::Phase)?;
            let point=binding.point.at(phase).map_err(MaterialNetworkError::Material)?;
            let mut requested=Vec::new();
            let operating_axes = resistor_cells
                .get(&binding.cell)
                .map(|resistor| {
                    conductivity::admitted_axes(electrical_axes, resistor.destination.branch)
                })
                .unwrap_or_default();
            if let Some((key,selection))=&binding.conductivity {
                let q=key.quantity();
                if q!=QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
                    && q!=QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity,ValueForm::Static)) {
                    return Err(invalid("scalar static conductivity key required"));
                }
                let requirement=ScalarPropertyRequirement::try_with_key(key,ScalarAdmissibility::NonNegative)
                    .map_err(|e|MaterialNetworkError::Material(ThermalMaterialError::Material(e)))?;
                requested.push((requirement,selection.clone()));
            }
            if let Some(resistor)=resistor_cells.get(&binding.cell) {
                let basis=bases[&resistor.destination.branch];
                let requirement=basis.requirement()
                    .map_err(|e|MaterialNetworkError::Material(ThermalMaterialError::Material(e)))?;
                requested.push((requirement,resistor.selection.clone()));
            }
            let has_directional = if let Some(directional) = directional {
                directional.requests(binding.cell, &mut requested)?
            } else { false };
            let conductivity=if requested.is_empty() {None} else {
                let selection=if let Some(pins)=pinned {
                    pins.get(&binding.cell).ok_or(invalid("initial transport source pins are absent"))?.clone()
                } else {
                    let mut pins=Vec::with_capacity(requested.len());
                    for (requirement,selection) in &requested {
                        let selected=resolve_material_state_point(binding.point.material,&point,
                            std::slice::from_ref(requirement),selection.clone())
                            .map_err(|e|MaterialNetworkError::Material(ThermalMaterialError::Material(e)))?;
                        let property=&selected.properties()[0];
                        pins.push((requirement.name().to_owned(),property.answer().receipt.selected));
                    }
                    MaterialPropertySelection::PinnedByProperty(pins)
                };
                let requirements=requested.into_iter().map(|(r,_)|r).collect::<Vec<_>>();
                let material=resolve_material_state_point(binding.point.material,&point,&requirements,selection)
                    .map_err(|e|MaterialNetworkError::Material(ThermalMaterialError::Material(e)))?;
                for property in material.properties() {
                    // Only the explicit conductivity route supplies E/J axes.
                    // Its complete bundle is refreshed at the solved electrical
                    // point before the thermal residual lowers these links.
                    binding
                        .point
                        .admit_static_with_axes(property, &operating_axes)
                        .map_err(MaterialNetworkError::Material)?;
                }
                let conductivity=(binding.conductivity.is_some() || has_directional).then(|| material.clone());
                materials.insert(binding.cell,material);
                conductivity
            };
            samples.insert(binding.cell,CellMaterialSample {cell:binding.cell,phase,conductivity});
        }
        let mut resistors=Vec::with_capacity(self.resistors.len());
        for resistor in self.resistors {
            check()?;
            let d=resistor.destination;
            let binding=bindings.get(&d.cell).ok_or(invalid("resistor cell lacks material binding"))?;
            let cell=&thermal.cells()[&d.cell];
            let body=MaterialOhmicBody {material:binding.point,selection:resistor.selection.clone(),
                length_m:resistor.length_m,area_m2:resistor.area_m2};
            let material=materials.remove(&d.cell).ok_or(invalid("resistor material bundle is absent"))?;
            let operating_axes = conductivity::admitted_axes(electrical_axes, d.branch);
            let conductor = body
                .conductor_from_material_with_axes(
                    samples[&d.cell].phase,
                    cell.geometry.volume_m3,
                    material,
                    bases[&d.branch],
                    &operating_axes,
                )
                .map_err(MaterialNetworkError::Conductor)?;
            resistors.push(MaterialResistor {branch:d.branch,conductor});
        }
        let (mut conductances, links) = self.resolve_links(thermal, &bindings, &samples, directional, &mut check)?;
        if let Some(additional) = additional_links {
            for (&id, &conductance) in additional {
                check()?;
                conductances.insert(id, conductance);
            }
        }
        check()?;
        Ok(CircuitThermalCoefficients {
            resistors,
            sources_w: BTreeMap::new(),
            conductances_w_k: conductances,
            source: MaterialNetworkSample {
                cells: samples.into_values().collect(),
                links,
            },
        })
    }

    // One canonical half-cell lowering serves seed resolution and actual-field
    // transport. The returned observations describe the coefficients supplied
    // to the same backward-Euler residual.
    fn resolve_links(
        &self,
        thermal: &DepositionSimulation<'_>,
        bindings: &BTreeMap<CellId, &MaterialCell<'_>>,
        samples: &BTreeMap<CellId, CellMaterialSample>,
        directional: Option<&directional::DirectionalPlan>,
        check: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(BTreeMap<u64, f64>, Vec<LinkMaterialSample>), MaterialNetworkError> {
        let invalid = MaterialNetworkError::Invalid;
        let mut conductances = BTreeMap::new();
        let mut links = BTreeMap::new();
        for source in self.links {
            check()?;
            let link = thermal
                .links()
                .get(&source.link)
                .ok_or(invalid("inactive material thermal link"))?;
            if !source.distance_a_m.is_finite()
                || source.distance_a_m <= 0.0
                || !source.distance_b_m.is_finite()
                || source.distance_b_m <= 0.0
                || links.contains_key(&source.link)
            {
                return Err(invalid(
                    "unique links and positive finite conduction distances required",
                ));
            }
            let mut k = [0.0; 2];
            for (i, id) in [link.a, link.b].into_iter().enumerate() {
                let material = samples.get(&id)
                    .ok_or(invalid("linked material sample is absent"))?.conductivity.as_ref()
                    .ok_or(invalid("linked conductivity is absent"))?;
                if let Some(directional) = directional {
                    if let Some(projected) = directional.project(source.link, i, material)? {
                        k[i] = projected;
                        continue;
                    }
                }
                let (key, _) = bindings
                    .get(&id)
                    .ok_or(invalid("linked cell has no material binding"))?
                    .conductivity
                    .as_ref()
                    .ok_or(invalid("linked cell has no conductivity binding"))?;
                k[i] = material.property(key.name())
                    .ok_or(invalid("linked conductivity key is absent"))?
                    .value_si();
            }
            let g = fs_thermal::enthalpy::series_conductance_w_k(
                link.area_m2,
                [source.distance_a_m, source.distance_b_m],
                k,
            )?;
            conductances.insert(source.link, g);
            links.insert(
                source.link,
                LinkMaterialSample {
                    link: source.link,
                    cells: [link.a, link.b],
                    conductivity_w_m_k: k,
                    conductance_w_k: g,
                },
            );
        }
        check()?;
        Ok((conductances, links.into_values().collect()))
    }

    /// Data-only network advance. External authored cell boundaries remain active;
    /// use `resolve` with the generic driver to add a separate source callback.
    /// Every Ohmic resistor gets one destination from this SAME binding list.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'m>(&self, thermal:&DepositionSimulation<'m>, target_s:f64,
        circuit:&MaterialCircuit, ports:&[f64], controls:CircuitEnthalpySolve,
        guess_j:&[f64], mut cancelled:impl FnMut()->bool)
        -> Result<CircuitEnthalpyCandidate<'m,MaterialNetworkSample>,CircuitEnthalpyError<MaterialNetworkError>> {
        if cancelled() {return Err(CircuitEnthalpyError::Cancelled);}
        if self.resistors.len()>MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(CircuitEnthalpyError::Invalid("material resistor bindings exceed circuit capacity"));
        }
        let bases=self.resistors.iter().map(|r|(r.destination.branch,ElectricalPropertyBasis::Resistivity)).collect::<Vec<_>>();
        self.step_candidate_with_bases(thermal,target_s,circuit,ports,&bases,controls,guess_j,cancelled)
    }

    /// Run the existing coupled solve with an explicit electrical basis for
    /// every resistor. Initially selected transport claims stay pinned through
    /// all nonlinear trials and final admission; crossing their support refuses
    /// the whole candidate. Supplied source pins also remain authoritative.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_bases<'m>(&self, thermal:&DepositionSimulation<'m>, target_s:f64,
        circuit:&MaterialCircuit, ports:&[f64], electrical_bases:&[(usize,ElectricalPropertyBasis)],
        controls:CircuitEnthalpySolve, guess_j:&[f64], cancelled:impl FnMut()->bool)
        -> Result<CircuitEnthalpyCandidate<'m,MaterialNetworkSample>,CircuitEnthalpyError<MaterialNetworkError>> {
        let cancel=RefCell::new(cancelled);
        if (*cancel.borrow_mut())() { return Err(CircuitEnthalpyError::Cancelled); }
        if self.resistors.len()>MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(CircuitEnthalpyError::Invalid("material resistor bindings exceed circuit capacity"));
        }
        let prepared=match self.prepare_with_bases(thermal,electrical_bases,||(*cancel.borrow_mut())()) {
            Err(MaterialNetworkError::Thermal(DepositionError::Cancelled)) => return Err(CircuitEnthalpyError::Cancelled),
            other => other.map_err(CircuitEnthalpyError::Source)?,
        };
        let destinations:Vec<_>=self.resistors.iter().map(|r|r.destination).collect();
        let result=super::step_candidate(thermal,target_s,CircuitEnthalpyInput {circuit,ports,destinations:&destinations},
            controls,guess_j,|s|prepared.resolve(s,||(*cancel.borrow_mut())()),||(*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(MaterialNetworkError::Thermal(DepositionError::Cancelled))) => Err(CircuitEnthalpyError::Cancelled),
            other=>other,
        }
    }

    /// Data-only Joule heating and reciprocal field-gap motion. Conductivity or
    /// resistivity and thermal transport are resolved from these SAME caloric
    /// cells at every proposed endpoint; gap/current/momentum come from the
    /// original actuator. One electrical basis is explicit for each resistor.
    /// Initial source choices stay pinned through the whole nonlinear interval.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_candidate<'m>(&self, thermal:&DepositionSimulation<'m>, target_s:f64,
        motion:&crate::material_circuit::motion::MaterialFieldMotion, ports:&[f64], external_force_n:f64,
        electrical_bases:&[(usize,ElectricalPropertyBasis)], controls:CircuitEnthalpySolve,
        guess_j:&[f64], cancelled:impl FnMut()->bool)
        -> Result<super::motion::FieldMotionEnthalpyCandidate<'m,MaterialNetworkSample>,CircuitEnthalpyError<MaterialNetworkError>> {
        let cancel=RefCell::new(cancelled);
        if (*cancel.borrow_mut())() {return Err(CircuitEnthalpyError::Cancelled);}
        let prepared=match self.prepare_with_bases(thermal,electrical_bases,||(*cancel.borrow_mut())()) {
            Err(MaterialNetworkError::Thermal(DepositionError::Cancelled)) => return Err(CircuitEnthalpyError::Cancelled),
            other=>other.map_err(CircuitEnthalpyError::Source)?,
        };
        let destinations=self.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
        let result=super::motion::step_candidate(thermal,target_s,super::motion::FieldMotionEnthalpyInput {
            motion,ports,external_force_n,destinations:&destinations},controls,guess_j,
            |state|prepared.resolve(state,||(*cancel.borrow_mut())()),||(*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(MaterialNetworkError::Thermal(DepositionError::Cancelled))) => Err(CircuitEnthalpyError::Cancelled),
            other=>other,
        }
    }
}
