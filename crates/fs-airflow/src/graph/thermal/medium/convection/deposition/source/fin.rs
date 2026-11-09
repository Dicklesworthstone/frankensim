//! Transient fins on the original spatial enthalpy and material heat owners.
//!
//! Axial cells own explicit mass and the original source calorics. Exact
//! source-integrated half-cell conduction connects them to each other and the
//! finite base. Each side has its own original finite-fluid film. An exposed
//! tip additionally resolves the last axial half-cell and its actual surface
//! temperature. This module supplies coefficients, never a second heat solver.

mod observer;
pub use observer::SourceTransientFinObserver;

use super::links::conjugate::{
    PreparedSourceConjugateLinks, ResolvedSourceConjugateHeat, SourceConjugateError,
    SourceConjugateLinkBinding,
};
use super::links::{ResolvedLinkConvection, SourceConvectiveLink, evaluate_source_link_convection};
use super::observer::with_checkpoint;
use super::{SourceCellConvectionError, ThermalFluidTransport};
use crate::graph::thermal::medium::convection::deposition::fin::TransientFinTopology;
use fs_material::phase::SolidLiquidPhase;
use fs_material::thermal::conduction::{
    LayerConductionConfig, ReservoirLayerGeometry,
    network::{LayerLinkBinding, LayerLinkSample, PreparedLayerConduction},
    source::{LayerSourceError, ResolvedConductionLayer},
};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, ThermalLaw, backward_euler::StepCoefficients,
};
use std::collections::{BTreeMap, BTreeSet};

/// Explicit spatial-approximation and source-work controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientFinSourceControls {
    /// All root, interior, side and tip links together.
    pub maximum_links: usize,
    /// Sampled lateral h A_c/(k P) at each actual fin-cell temperature.
    /// This admits transverse lumping, not the independent axial mesh error.
    pub maximum_transverse_biot: f64,
    /// Original interface-power tolerances apply to each axial/root/tip path.
    /// They do not replace the native integrator's joule and timestep controls.
    pub layer: LayerConductionConfig,
}

#[derive(Debug)]
pub enum SourceTransientFinError<E = core::convert::Infallible> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Layer(LayerSourceError),
    Convection(SourceCellConvectionError<E>),
    Tip(SourceConjugateError<E>),
    TransverseBiot {
        cell: CellId,
        actual: f64,
        maximum: f64,
    },
}
impl<E> From<DepositionError> for SourceTransientFinError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<LayerSourceError> for SourceTransientFinError<E> {
    fn from(error: LayerSourceError) -> Self {
        Self::Layer(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SourceTransientFinError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source transient fin: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SourceTransientFinError<E> {}

/// Original material and film results at one actual native thermal state.
#[derive(Debug)]
pub struct SourceTransientFinSample<F> {
    pub conduction: Vec<LayerLinkSample>,
    pub sides: ResolvedLinkConvection<F>,
    /// Only the combined tip coefficient is applied. Its bare film is evidence.
    pub tip: Option<ResolvedSourceConjugateHeat<F>>,
    pub maximum_sampled_transverse_biot: f64,
    /// Sum of native m*h for the fin cells only, on their original datum [J].
    /// Already part of DepositionSimulation::energy_j, never additional storage.
    pub fin_energy_j: f64,
    /// Actual native G*(T_base-T_first_cell) [W], including root resistance.
    pub root_heat_into_fin_w: f64,
    /// Sum of actual lateral and combined tip powers into the finite fluid [W].
    pub heat_into_fluid_w: f64,
    /// Native authored fin source watts; other composed observers are separate.
    pub authored_fin_source_w: f64,
    /// Root heat + authored source - fluid heat. Fin storage permits this to
    /// be nonzero; it is a physical instantaneous power, not a heat defect.
    pub net_power_into_fin_w: f64,
    /// Sum of absolute axial/root native-vs-integral and layer heat defects [W].
    pub absolute_conduction_power_defect_w: f64,
}

/// Immutable source and topology binding for one finite, fixed-shape solid fin.
///
/// The same original all-solid fin caloric object must own every fin cell.
/// The base may have a different original solid material. Neither current
/// density nor a steady efficiency changes the explicitly declared masses.
/// Arbitrary solid Cp(T) is retained; entering a melting branch refuses until
/// a separately admitted phase/geometry model is selected.
///
/// Side films are local closure assumptions with their authored correlation
/// lengths and property conventions. No fin-array interference, resolved
/// boundary layer, transverse temperature field or momentum law is inferred.
/// The fluid callback must retain the finite fluid cell's original material,
/// pressure/composition and source choices at every trial.
pub struct PreparedSourceTransientFin<'t, 'a, 'm> {
    topology: &'t TransientFinTopology,
    fin_cells: BTreeSet<CellId>,
    owned_links: BTreeSet<u64>,
    fin_material: &'a ResolvedConductionLayer<'m>,
    base_material: &'a ResolvedConductionLayer<'m>,
    conduction: PreparedLayerConduction<'a, 'm>,
    sides: Vec<SourceConvectiveLink>,
    tip: Option<PreparedSourceConjugateLinks<'a, 'm>>,
    controls: TransientFinSourceControls,
}

impl<'t, 'a, 'm> PreparedSourceTransientFin<'t, 'a, 'm> {
    /// Bind an already activated native fin; every side and selected tip needs
    /// exactly one original film declaration. Order is canonicalized by link ID.
    /// Base bare-area cooling stays a separate disjoint physical surface.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        topology: &'t TransientFinTopology,
        base_material: &'a ResolvedConductionLayer<'m>,
        fin_material: &'a ResolvedConductionLayer<'m>,
        convection: &[SourceConvectiveLink],
        controls: TransientFinSourceControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SourceTransientFinError> {
        checkpoint()?;
        if topology.links().len() > controls.maximum_links {
            return Err(DepositionError::Budget("transient fin source links").into());
        }
        if !controls.maximum_transverse_biot.is_finite()
            || controls.maximum_transverse_biot < 0.0
            || controls.layer.max_iterations == 0
            || !controls.layer.absolute_heat_tolerance_w.is_finite()
            || controls.layer.absolute_heat_tolerance_w < 0.0
            || !controls.layer.relative_heat_tolerance.is_finite()
            || controls.layer.relative_heat_tolerance < 0.0
            || base_material.phase() != SolidLiquidPhase::Solid
            || fin_material.phase() != SolidLiquidPhase::Solid
        {
            return Err(SourceTransientFinError::Invalid(
                "fixed fins need original all-solid materials and finite bounded source controls",
            ));
        }
        topology
            .links()
            .len()
            .checked_mul(controls.layer.max_iterations)
            .ok_or(DepositionError::Budget(
                "transient fin source work overflow",
            ))?;
        let fin_cells = topology.cells().iter().map(|cell| cell.id).collect();
        let owned_links = topology.links().iter().map(|link| link.id).collect();
        admit_topology(thermal, topology, &fin_cells, &owned_links, &mut checkpoint)?;
        let expected: BTreeSet<_> = topology
            .side_links()
            .iter()
            .copied()
            .chain(topology.tip_link())
            .collect();
        if convection.len() != expected.len() {
            return Err(SourceTransientFinError::Invalid(
                "every fin side and selected tip needs exactly one film law",
            ));
        }
        let mut mapped = BTreeMap::new();
        for &binding in convection {
            checkpoint()?;
            if binding.fluid_cell != topology.fluid_cell()
                || !expected.contains(&binding.link)
                || mapped.insert(binding.link, binding).is_some()
            {
                return Err(SourceTransientFinError::Invalid(
                    "fin films must own the original distinct side/tip links and finite fluid",
                ));
            }
        }
        let mut layers = Vec::with_capacity(topology.conductive_links().len());
        for binding in topology.conductive_links() {
            checkpoint()?;
            layers.push(LayerLinkBinding {
                link: binding.link_id,
                layers: [
                    if binding.link_id == topology.root_link() {
                        base_material
                    } else {
                        fin_material
                    },
                    fin_material,
                ],
                distances_m: binding.half_distances_m,
                contact_area_resistance_m2_k_w: binding.contact_resistance_m2_k_w,
            });
        }
        let conduction =
            PreparedLayerConduction::new(thermal, layers, controls.maximum_links, &mut checkpoint)?;
        let tip = if let Some(id) = topology.tip_link() {
            let edges = topology.axial_edges_m();
            let half = 0.5 * (edges[edges.len() - 1] - edges[edges.len() - 2]);
            Some(
                PreparedSourceConjugateLinks::new(
                    thermal,
                    &[SourceConjugateLinkBinding {
                        convection: mapped
                            .remove(&id)
                            .ok_or(SourceTransientFinError::Invalid("missing fin tip film"))?,
                        wall_layer: fin_material,
                        geometry: ReservoirLayerGeometry::planar(
                            topology.geometry().cross_section_m2(),
                            half,
                        )?,
                    }],
                    1,
                    &mut checkpoint,
                )
                .map_err(SourceTransientFinError::Tip)?,
            )
        } else {
            None
        };
        let result = Self {
            topology,
            fin_cells,
            owned_links,
            fin_material,
            base_material,
            conduction,
            sides: mapped.into_values().collect(),
            tip,
            controls,
        };
        result.admit_materials::<core::convert::Infallible>(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }

    pub fn topology(&self) -> &TransientFinTopology {
        self.topology
    }
    pub const fn controls(&self) -> TransientFinSourceControls {
        self.controls
    }
    pub fn side_bindings(&self) -> &[SourceConvectiveLink] {
        &self.sides
    }

    fn admit_materials<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SourceTransientFinError<E>> {
        admit_topology(
            thermal,
            self.topology,
            &self.fin_cells,
            &self.owned_links,
            checkpoint,
        )?;
        for (id, material) in std::iter::once((self.topology.base_cell(), self.base_material))
            .chain(
                self.topology
                    .cells()
                    .iter()
                    .map(|cell| (cell.id, self.fin_material)),
            )
        {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(SourceTransientFinError::Invalid("inactive fin/base cell"))?;
            if !thermal.uses_material_law(
                cell.material_id,
                material.caloric_curve() as &dyn ThermalLaw,
            ) {
                return Err(SourceTransientFinError::Invalid(
                    "fin storage must retain its original source caloric object",
                ));
            }
            material.temperature_at_enthalpy(cell.specific_h_j_kg)?;
            // The original native Kelvin value is authoritative, even at an
            // enthalpy support endpoint whose inverse snaps to a source knot.
            material
                .conductivity_curve()
                .conductivity_w_m_k(thermal.temperature_k(id)?)?;
        }
        Ok(())
    }

    /// Resolve original k(T), surface films and reciprocal heat at this state.
    /// Use these coefficients inside every native implicit trial, or use the
    /// observer below for the existing explicit enthalpy timeline. Sample watts
    /// and stored joules are diagnostics; do not deposit or add them again.
    pub fn resolve<F: ThermalFluidTransport, E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut resolve_fluid: impl FnMut(
            &SourceConvectiveLink,
            &DepositionSimulation<'_>,
            Temperature,
        ) -> Result<F, E>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<StepCoefficients<SourceTransientFinSample<F>>, SourceTransientFinError<E>> {
        self.admit_materials(thermal, &mut checkpoint)?;
        let mut axial = self
            .conduction
            .resolve(thermal, self.controls.layer, &mut checkpoint)?;
        let sides = with_checkpoint(&mut checkpoint, |cancelled| {
            evaluate_source_link_convection(
                thermal,
                &self.sides,
                self.controls.maximum_links,
                &mut resolve_fluid,
                cancelled,
            )
        })
        .map_err(SourceTransientFinError::Convection)?;
        let mut maximum_biot: f64 = 0.0;
        let geometry = self.topology.geometry();
        for sample in sides.samples() {
            checkpoint()?;
            let k = self
                .fin_material
                .conductivity_curve()
                .conductivity_w_m_k(sample.wall_temperature_k)?;
            let h = sample.conductance_w_k / sample.area_m2;
            let biot = (h / k) * (geometry.cross_section_m2() / geometry.exposed_perimeter_m());
            if !biot.is_finite() || biot < 0.0 || biot > self.controls.maximum_transverse_biot {
                return Err(SourceTransientFinError::TransverseBiot {
                    cell: sample.wall_cell,
                    actual: biot,
                    maximum: self.controls.maximum_transverse_biot,
                });
            }
            maximum_biot = maximum_biot.max(biot);
        }
        let tip = self
            .tip
            .as_ref()
            .map(|tip| {
                tip.resolve(
                    thermal,
                    self.controls.layer,
                    &mut resolve_fluid,
                    &mut checkpoint,
                )
                .map_err(SourceTransientFinError::Tip)
            })
            .transpose()?;
        merge_disjoint(
            &mut axial.conductances_w_k,
            sides.conductances_w_k(),
            &mut checkpoint,
        )?;
        if let Some(tip) = &tip {
            merge_disjoint(
                &mut axial.conductances_w_k,
                tip.conductances_w_k(),
                &mut checkpoint,
            )?;
        }
        let mut defect = Sum::default();
        let mut root_heat = None;
        for sample in &axial.payload {
            checkpoint()?;
            let a = thermal.temperature_k(sample.cells[0])?;
            let b = thermal.temperature_k(sample.cells[1])?;
            let native = sample.response.conductance_w_k * (a - b);
            let error = native - sample.response.power_a_to_b_w;
            let allowance = self.controls.layer.absolute_heat_tolerance_w
                + self.controls.layer.relative_heat_tolerance
                    * native.abs().max(sample.response.power_a_to_b_w.abs());
            if !native.is_finite() || !allowance.is_finite() || error.abs() > allowance {
                return Err(DepositionError::Arithmetic.into());
            }
            defect.add(error.abs() + sample.response.layer_heat_residual_w.abs())?;
            if sample.link == self.topology.root_link() {
                root_heat = Some(native);
            }
        }
        let mut fluid_heat = Sum::default();
        for sample in sides.samples() {
            checkpoint()?;
            fluid_heat.add(sample.heat_into_fluid_w)?;
        }
        if let Some(tip) = &tip {
            for sample in tip.samples() {
                checkpoint()?;
                fluid_heat.add(sample.applied_heat_into_fluid_w)?;
            }
        }
        let mut energy = Sum::default();
        let mut sources = Sum::default();
        for expected in self.topology.cells() {
            checkpoint()?;
            let cell = &thermal.cells()[&expected.id];
            energy.add(cell.mass_kg * cell.specific_h_j_kg)?;
            sources.add(cell.boundary.source_w)?;
        }
        let root_heat = root_heat.ok_or(SourceTransientFinError::Invalid(
            "fin has no root heat owner",
        ))?;
        let net = root_heat + sources.value() - fluid_heat.value();
        if !net.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        checkpoint()?;
        Ok(StepCoefficients {
            sources_w: axial.sources_w,
            conductances_w_k: axial.conductances_w_k,
            payload: SourceTransientFinSample {
                conduction: axial.payload,
                sides,
                tip,
                maximum_sampled_transverse_biot: maximum_biot,
                fin_energy_j: energy.value(),
                root_heat_into_fin_w: root_heat,
                heat_into_fluid_w: fluid_heat.value(),
                authored_fin_source_w: sources.value(),
                net_power_into_fin_w: net,
                absolute_conduction_power_defect_w: defect.value(),
            },
        })
    }
}

// Geometry and fixed fin inventory remain authoritative at every actual state.
// Authored watts and ambient temperature may change through ordinary events.
fn admit_topology(
    thermal: &DepositionSimulation<'_>,
    topology: &TransientFinTopology,
    fin_cells: &BTreeSet<CellId>,
    owned_links: &BTreeSet<u64>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    checkpoint()?;
    for expected in topology.cells() {
        checkpoint()?;
        let actual = thermal
            .cells()
            .get(&expected.id)
            .ok_or(DepositionError::InvalidInput("inactive transient fin cell"))?;
        if actual.material_id != expected.material_id
            || actual.mass_kg != expected.mass_kg
            || actual.geometry != expected.geometry
            || actual.body_id != expected.body_id
            || actual.segment_id != expected.segment_id
            || actual.boundary.htc_w_m2_k != 0.0
            || actual.boundary.reservoir.is_some()
        {
            return Err(DepositionError::InvalidInput(
                "fin inventory or boundary changed; rebuild its physical binding",
            ));
        }
    }
    for expected in topology.links() {
        checkpoint()?;
        let patch = thermal
            .contact_patch(expected.id)
            .ok_or(DepositionError::InvalidInput("inactive transient fin link"))?;
        if patch.link != *expected || patch.covered_area_m2 != [expected.area_m2; 2] {
            return Err(DepositionError::InvalidInput(
                "fin link or surface coverage changed",
            ));
        }
    }
    // The fin owns all selected faces, including an explicitly insulated tip.
    // Base and fluid can have other disjoint physical connections. This native
    // graph sweep is bounded by its existing link budget; sets are prepared
    // once, not allocated on each nonlinear or explicit source query.
    for (index, link) in thermal.links().values().enumerate() {
        if index % 64 == 0 {
            checkpoint()?;
        }
        if (fin_cells.contains(&link.a) || fin_cells.contains(&link.b))
            && !owned_links.contains(&link.id)
        {
            return Err(DepositionError::InvalidInput(
                "unowned contact changes the prepared fin boundary",
            ));
        }
    }
    checkpoint()?;
    Ok(())
}

fn merge_disjoint(
    destination: &mut BTreeMap<u64, f64>,
    source: &BTreeMap<u64, f64>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    for (&id, &value) in source {
        checkpoint()?;
        if destination.insert(id, value).is_some() {
            return Err(DepositionError::InvalidInput(
                "fin heat link has multiple owners",
            ));
        }
    }
    Ok(())
}

#[derive(Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        let adjusted = value - self.correction;
        let next = self.total + adjusted;
        self.correction = (next - self.total) - adjusted;
        self.total = next;
        if !value.is_finite() || !self.total.is_finite() || !self.correction.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        Ok(())
    }
    fn value(&self) -> f64 {
        self.total
    }
}
