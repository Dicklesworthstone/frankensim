//! Original source-owned conduction/convection on the native explicit heat history.
use super::*;
use crate::solid_enthalpy::{SolidEnthalpy, TetCellMaterial};
use crate::solid_enthalpy::gas::{PreparedContact, SolidGasConvectionSample, StructuralSource};
use crate::solid_enthalpy::gas::convection::endpoint_law;
use fs_material::gas_species::input::temperature::IsobaricGasLaw;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_material::thermal::conduction::source::ResolvedConductionLayer;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use std::sync::Arc;

/// One source film's accepted integral and last applied native substep sample.
#[derive(Clone, Debug)]
pub struct CavityContactHeat {
    pub id: u64,
    pub cell: CellId,
    pub heat_into_solid_j: f64,
    /// Actual cell temperature when the last explicit substep's power was
    /// sampled; this can differ from the returned wall endpoint temperature.
    pub last_cell_temperature_k: f64,
    pub last_power_into_solid_w: f64,
    /// Held-cell-k series diagnostic, or the integrated profile's actual flux
    /// secant. Actual power is authoritative; the explicit timestep retains
    /// its separate flux secant with an explicit zero at equilibrium.
    pub last_conductance_w_k: f64,
    /// Original conduction budget, or actual convective power divided by
    /// gas/cell contrast [W/K]. A convective equilibrium uses zero, not an
    /// inferred natural-convection tangent.
    pub last_explicit_decay_conductance_w_k: f64,
    pub source: ResolvedMaterialStatePoint,
    /// Last applied substep's actual surface, source queries and power residual.
    pub convection: Option<SolidGasConvectionSample>,
}

#[derive(Clone)]
pub(super) struct Film<'m> {
    contact: SolidGasContact,
    element: usize,
    integrated_solid: Option<Arc<ResolvedConductionLayer<'m>>>,
}
impl<'m> Film<'m> {
    fn prepared(&self) -> PreparedContact<'_, 'm> {
        PreparedContact::new(&self.contact, self.element, self.integrated_solid.clone())
    }
}

fn binding<'a>(source: &TetCellMaterial<'a>, reference_h_j_kg: f64) -> StructuralSource<'a> {
    StructuralSource { cell: source.cell, curve: source.curve,
        point: source.point, reference_h_j_kg }
}

// This is the current explicit flux coefficient, not a derivative or a
// time-error estimate. Holding it for the native Euler substep preserves the
// original nonnegative-decay budget even when the surface balance is tolerant.
fn flux_decay<E>(power_w: f64, difference_k: f64) -> Result<f64, CavityError<E>> {
    if difference_k == 0.0 {
        return if power_w == 0.0 { Ok(0.0) }
            else { Err(CavityError::Invalid("equilibrium gas film must carry zero heat")) };
    }
    let conductance = number(power_w / difference_k)?;
    if conductance < 0.0 || (power_w != 0.0 && conductance == 0.0) {
        return Err(CavityError::Invalid("actual gas-film decay is negative or underflows"));
    }
    Ok(conductance)
}

fn phase<E>(
    source: &TetCellMaterial<'_>, thermal: &DepositionSimulation<'_>,
) -> Result<EquilibriumPhaseState, CavityError<E>> {
    let cell = thermal.cells().get(&source.cell)
        .ok_or(CavityError::Invalid("inactive source-bound gas film cell"))?;
    Ok(source.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
        .map_err(SolidEnthalpyError::from)?)
}

pub(super) fn prepare<'m>(
    source: &SolidEnthalpy<'m>, thermal: &DepositionSimulation<'_>, gas: &ClosedGasParcel,
    contacts: &[SolidGasContact], cx: &Cx<'_>,
) -> Result<Vec<Film<'m>>, CavityError> {
    let mut ordered = contacts.to_vec();
    ordered.sort_by_key(|contact| contact.id);
    let mut areas = BTreeMap::<CellId, f64>::new();
    let mut films = Vec::with_capacity(ordered.len());
    let law = if ordered.iter().any(|contact| contact.convection.is_some()) {
        Some(endpoint_law(gas.model(), gas.point(), gas.state())?)
    } else { None };
    for (i, contact) in ordered.iter().enumerate() {
        poll(cx)?;
        if i > 0 && ordered[i - 1].id == contact.id {
            return Err(CavityError::Invalid("unique gas-film patch IDs required"));
        }
        let element = source.cells.iter().position(|bound| bound.source.cell == contact.cell)
            .ok_or(CavityError::Invalid("gas film must name an actual source-bound structural cell"))?;
        let cell = thermal.cells().get(&contact.cell)
            .ok_or(CavityError::Invalid("inactive gas-film cell"))?;
        let area = areas.entry(contact.cell).or_insert(0.0);
        *area = number(*area + contact.geometry.area_m2)?;
        if !contact.geometry.area_m2.is_finite() || contact.geometry.area_m2 < 0.0
            || *area > thermal.exposed_area_m2(contact.cell)?
            || cell.boundary.htc_w_m2_k != 0.0
        {
            return Err(CavityError::Invalid("gas films exceed exposed area or duplicate ambient convection"));
        }
        let bound = &source.cells[element];
        let reference = bound.source.curve.state_at_specific_enthalpy(
            bound.reference_cell.specific_h_j_kg,
        ).map_err(SolidEnthalpyError::from)?;
        let source = binding(&bound.source, bound.reference_cell.specific_h_j_kg);
        let mut film = Film { contact: contact.clone(), element, integrated_solid: None };
        // Pin only the material at its original reference. A convection regime
        // at that irrelevant stress-free temperature is not an admission gate.
        let material = film.prepared().reference_material(&source, gas.state(), reference)?;
        let selected = material.property_by_key(&contact.conductivity)
            .ok_or(CavityError::Invalid("gas film lacks its selected conductivity"))?
            .answer().receipt.selected;
        film.contact.selection = MaterialPropertySelection::PinnedByProperty(vec![
            (contact.conductivity.name().to_owned(), selected),
        ]);
        film.integrated_solid = film.prepared().prepare_solid_profile(
            &source, &mut || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
        )?;
        film.prepared().resolve_with_law(&source, gas.state(), phase(&bound.source, thermal)?,
            None, law.as_ref(), &mut || cx.checkpoint().map_err(|_| DepositionError::Cancelled))?;
        films.push(film);
    }
    poll(cx)?;
    Ok(films)
}

#[derive(Clone)]
pub(super) struct FilmObserver<'a, 's, 'm, O> {
    cavity: &'a ThermoelasticCavity<'s, 'm>,
    gas: &'a ResolvedThermochemicalGasState,
    // One full source law per physical gas endpoint. Cloned observers share
    // it while their own thermal substep histories remain independent.
    law: Option<Arc<IsobaricGasLaw>>,
    pub(super) inner: O,
    applied_w: BTreeMap<CellId, f64>,
    conductances_w_k: BTreeMap<CellId, f64>,
    pub(super) samples: Vec<CavityContactHeat>,
    pub(super) heat_j: f64,
}
impl<'a, 's, 'm, O> FilmObserver<'a, 's, 'm, O> {
    pub(super) fn new<E>(
        cavity: &'a ThermoelasticCavity<'s, 'm>, gas: &'a ClosedGasParcel, inner: O,
    ) -> Result<Self, CavityError<E>> {
        let law = if cavity.films.iter().any(|film| film.contact.convection.is_some()) {
            Some(Arc::new(endpoint_law(gas.model(), gas.point(), gas.state())?))
        } else { None };
        Ok(Self { cavity, gas: gas.state(), law, inner, applied_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(), samples: Vec::new(), heat_j: 0.0 })
    }
}
impl<O: DepositionObserver> DepositionObserver for FilmObserver<'_, '_, '_, O> {
    type Error = CavityError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.heat_j = 0.0;
        self.samples.clear();
        self.applied_w.clear();
        self.conductances_w_k.clear();
        self.inner.begin().map_err(CavityError::Observer)
    }
    fn sources(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let mut sources = self.inner.sources(thermal, checkpoint).map_err(CavityError::Observer)?;
        let mut watts = BTreeMap::<CellId, f64>::new();
        let mut conductances = BTreeMap::<CellId, f64>::new();
        let mut areas = BTreeMap::<CellId, f64>::new();
        let mut samples = Vec::with_capacity(self.cavity.films.len());
        for (i, film) in self.cavity.films.iter().enumerate() {
            checkpoint()?;
            let bound = &self.cavity.solid.source.cells[film.element];
            let area = areas.entry(film.contact.cell).or_insert(0.0);
            *area = number(*area + film.contact.geometry.area_m2)?;
            let cell = thermal.cells().get(&film.contact.cell)
                .ok_or(CavityError::Invalid("inactive gas-film cell"))?;
            if *area > thermal.exposed_area_m2(film.contact.cell)?
                || cell.boundary.htc_w_m2_k != 0.0 {
                return Err(CavityError::Invalid("accepted gas-film coverage changed or duplicates ambient convection"));
            }
            let local = phase(&bound.source, thermal)?;
            let source = binding(&bound.source, bound.reference_cell.specific_h_j_kg);
            let sample = film.prepared().resolve_with_law(&source, self.gas, local,
                None, self.law.as_deref(), &mut || checkpoint())?;
            let power = number(-sample.power_to_gas_w)?;
            let decay = if sample.convection.is_some() {
                flux_decay(power, self.gas.state().temperature - local.temperature_k())?
            } else { sample.conductance_w_k };
            let old = watts.entry(film.contact.cell).or_insert(0.0);
            *old = number(*old + power)?;
            let g = conductances.entry(film.contact.cell).or_insert(0.0);
            *g = number(*g + decay)?;
            samples.push(CavityContactHeat {
                id: film.contact.id, cell: film.contact.cell,
                heat_into_solid_j: self.samples.get(i).map_or(0.0, |row| row.heat_into_solid_j),
                last_cell_temperature_k: local.temperature_k(),
                last_power_into_solid_w: power, last_conductance_w_k: sample.conductance_w_k,
                last_explicit_decay_conductance_w_k: decay,
                source: sample.source,
                convection: sample.convection,
            });
        }
        let mut applied = BTreeMap::new();
        for (&cell, &power) in &watts {
            checkpoint()?;
            let previous = sources.get(&cell).copied().unwrap_or(0.0);
            let next = number(previous + power)?;
            let actual = number(next - previous)?;
            let authored = thermal.cells().get(&cell)
                .ok_or(CavityError::Invalid("inactive gas-film source cell"))?.boundary.source_w;
            let before = number(authored + previous)?;
            let combined = number(authored + next)?;
            if (power != 0.0 && (actual == 0.0 || combined == before))
                || (previous != 0.0 && next == power)
                || (authored != 0.0 && combined == next)
            {
                return Err(CavityError::Invalid("gas-film heat is not represented in the native cell source"));
            }
            sources.insert(cell, next);
            applied.insert(cell, actual);
        }
        self.applied_w = applied;
        self.conductances_w_k = conductances;
        self.samples = samples;
        Ok(sources)
    }
    fn conductances(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner.conductances(thermal, checkpoint).map_err(CavityError::Observer)
    }
    fn reservoir_conductances(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.reservoir_conductances(thermal, checkpoint).map_err(CavityError::Observer)
    }

    fn ambient_conductances(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.ambient_conductances(thermal, checkpoint).map_err(CavityError::Observer)
    }
    fn heat_transfers(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner.heat_transfers(thermal, checkpoint).map_err(CavityError::Observer)
    }
    fn limit_timestep_s(
        &mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        let proposed = self.inner.limit_timestep_s(thermal, proposed_s, checkpoint)
            .map_err(CavityError::Observer)?;
        // Add the actual film decay to the native conductive budget. The
        // original caloric slope/mass owner supplies the bound; no Cp copy.
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed, &self.conductances_w_k, checkpoint,
        )?)
    }
    fn requires_step_interval(&self) -> bool { self.inner.requires_step_interval() }
    fn next_transition_time_s(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.inner.next_transition_time_s(thermal, checkpoint).map_err(CavityError::Observer)
    }
    fn step_sources(
        &mut self, thermal: &DepositionSimulation<'_>, start: f64, dt: f64, end: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.step_sources(thermal, start, dt, end, checkpoint).map_err(CavityError::Observer)
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.inner.event(event, thermal).map_err(CavityError::Observer)
    }
    fn step(
        &mut self, start: f64, dt: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner.step(start, dt, old, thermal, checkpoint).map_err(CavityError::Observer)?;
        for &power in self.applied_w.values() {
            checkpoint()?;
            self.heat_j = number(self.heat_j + dt * power)?;
        }
        for row in &mut self.samples {
            checkpoint()?;
            row.heat_into_solid_j = number(row.heat_into_solid_j + dt * row.last_power_into_solid_w)?;
        }
        Ok(())
    }
}
