//! Conducting moving liquid through the original finite thermal timeline.
//!
//! Geometry is supplied at the accepted VOF configuration. Each declared
//! Cartesian face has a wetted area and two liquid-centroid distances. Nothing
//! infers a contact from fill fraction. Undeclared faces are thermally insulated.
//! The source law supplies k(T) at every thermal substep, and the existing
//! timeline owns pair heat, caloric storage, CFL subcycling and energy admission.
//! Its accepted h feeds the SAME native momentum and mass-transport step.

use super::*;
use fs_conduction::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, EnergyAudit, Event, MaterialId, StepConfig,
    ThermalLaw, ThermalLink, TimedEvent, series_conductance_w_k, series_contact_conductance_w_k,
};

/// Finite source-bound solid cells exchanging heat with native stationary walls.
pub mod body;
pub use body::{
    FreeSurfaceBodyContact, FreeSurfaceBodyContactSample, FreeSurfaceBodyExchange,
    FreeSurfaceThermalBody,
};
/// Source-resolved gas-film exchange on explicitly exposed liquid faces.
pub mod environment;
pub use environment::{
    FreeSurfaceGasContact, FreeSurfaceGasExchange, FreeSurfaceGasFace, FreeSurfaceGasSample,
};
/// Finite atmospheric storage on the same thermal and native VOF interval.
pub mod atmosphere;
pub use atmosphere::FreeSurfaceFiniteGasExchange;
/// Native source-driven current and Joule heat in a declared fully wet core.
pub mod electrical;
pub use electrical::{
    FreeSurfaceElectricalExchange, FreeSurfaceElectricalInput, FreeSurfaceElectricalState,
};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Source-bound resistance at an ordered liquid/solid thermal contact.
pub mod interface;
pub use interface::FreeSurfaceBodyInterface;

/// Actual wet geometry of one face at the beginning of this flow step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceConductionFace {
    /// Two distinct face-adjacent Cartesian cells, in source order.
    pub cells: [[usize; 3]; 2],
    /// Shared liquid-contact area [m2], bounded by the Cartesian face area.
    /// Zero explicitly insulates, but does not waive endpoint/source admission.
    pub area_m2: f64,
    /// Each liquid thermal centroid's normal distance to that face [m].
    /// Both must be positive and no greater than the grid cell width.
    pub distances_m: [f64; 2],
}

/// Thermal work and accuracy bounds, independent of native flow/VOF budgets.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceConductionControls {
    pub maximum_faces: usize,
    pub maximum_substeps: usize,
    /// The original explicit thermal owner may shorten this via its CFL bound.
    /// Nonlinear k(T) and the thermal/flow split still require time refinement.
    pub maximum_dt_s: f64,
    pub energy_tolerance_j: f64,
}

/// The actual source evaluations used by one two-layer face conductance.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceConductionSample {
    pub face: FreeSurfaceConductionFace,
    pub materials: [ResolvedLiquidState; 2],
    pub conductance_w_k: f64,
}

#[derive(Clone, Debug)]
pub struct FreeSurfaceConductionReport {
    /// Original thermal-stage audit. Its pair heat is internal, not a source.
    pub thermal: EnergyAudit,
    /// Combined original-to-final VOF/enthalpy account over the SAME interval.
    pub flow: FreeSurfaceEnthalpyReport,
    /// Last conductances actually applied by the native thermal subcycling.
    pub last_applied: Vec<FreeSurfaceConductionSample>,
    /// Source-supported thermal endpoint BEFORE the subsequent VOF advection.
    pub thermal_endpoint: Vec<FreeSurfaceConductionSample>,
    /// Present only when finite solid cells participate in this timeline.
    pub bodies: Option<FreeSurfaceBodyExchange>,
    /// Environmental heat entering the liquid subsystem's external heat once.
    /// `finite_gas` supplies the opposite gas storage when that bath is finite.
    pub environment: Option<FreeSurfaceGasExchange>,
    /// Additional finite gas storage account when the environment has an owner.
    /// Gas/liquid exchange is internal to this combined thermal system.
    pub finite_gas: Option<FreeSurfaceFiniteGasExchange>,
    /// Actual native electrode supply and Joule heat, counted once in external heat.
    pub electrical: Option<FreeSurfaceElectricalExchange>,
}

fn checkpoint(gate: &CancelGate) -> Result<(), DepositionError> {
    if gate.is_requested() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

#[derive(Clone)]
struct Conduction<'a, 'm> {
    gate: &'a CancelGate,
    law: &'a LiquidTemperatureLaw<'m>,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    density: f64,
    buoyancy: Option<BoussinesqLiquid>,
    rheology: Option<rheology::RheologyModel<'m>>,
    rheology_stress: BTreeMap<CellId, rheology::RheologyStress>,
    faces: &'a [FreeSurfaceConductionFace],
    endpoints: Vec<[CellId; 2]>,
    liquid_cells: Vec<CellId>,
    liquid_coordinates: BTreeMap<CellId, [usize; 3]>,
    wetting: Option<WettingBinding<'m>>,
    bodies: &'a [FreeSurfaceThermalBody<'m>],
    body_ids: Vec<CellId>,
    contacts: &'a [FreeSurfaceBodyContact],
    interfaces: BTreeMap<usize, &'a FreeSurfaceBodyInterface<'a>>,
    contact_endpoints: Vec<(CellId, usize)>,
    powers: BTreeMap<CellId, f64>,
    last: Vec<FreeSurfaceConductionSample>,
    endpoint: Vec<FreeSurfaceConductionSample>,
    body_last: Vec<FreeSurfaceBodyContactSample>,
    body_endpoint: Vec<FreeSurfaceBodyContactSample>,
    body_states_endpoint: Vec<ResolvedMaterialStatePoint>,
    environment: &'a [FreeSurfaceGasContact<'a>],
    environment_endpoints: Vec<CellId>,
    environment_last: Vec<FreeSurfaceGasSample>,
    environment_endpoint: Vec<FreeSurfaceGasSample>,
    environment_reference: Vec<FreeSurfaceGasSample>,
    environment_heat_j: f64,
    environment_absolute_heat_j: f64,
    finite_atmosphere: Option<atmosphere::FiniteAtmosphere>,
    electrical: Option<electrical::ElectricalTimeline>,
}

struct Samples {
    liquid: Vec<FreeSurfaceConductionSample>,
    contacts: Vec<FreeSurfaceBodyContactSample>,
    bodies: Vec<ResolvedMaterialStatePoint>,
}

impl Conduction<'_, '_> {
    fn sample_liquid(
        &self,
        id: CellId,
        h: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedLiquidState, Error> {
        if let Some(model) = &self.rheology {
            Ok(model
                .sample(h, self.rheology_stress[&id], || {
                    check().map_err(Error::from)
                })?
                .material)
        } else {
            Ok(resolve_liquid(self.law, self.caloric, h, self.density, self.buoyancy)?.0)
        }
    }

    fn sample(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Samples, Error> {
        let mut states = BTreeMap::new();
        for &id in &self.liquid_cells {
            check()?;
            let cell = &thermal.cells()[&id];
            let state = self.sample_liquid(id, cell.specific_h_j_kg, check)?;
            if let Some(wetting) = &self.wetting {
                wetting.admit_temperature(
                    self.liquid_coordinates[&id], self.law.temperature_axis(),
                    super::wetting::temperature(&state, self.law.temperature_axis())?,
                )?;
            }
            states.insert(id, state);
        }
        let mut samples = Vec::with_capacity(self.faces.len());
        for (&face, &[a, b]) in self.faces.iter().zip(&self.endpoints) {
            check()?;
            let materials = [states[&a].clone(), states[&b].clone()];
            let conductance_w_k = series_conductance_w_k(
                face.area_m2,
                face.distances_m,
                materials.each_ref().map(|s| s.thermal_conductivity_w_m_k()),
            )?;
            samples.push(FreeSurfaceConductionSample {
                face,
                materials,
                conductance_w_k,
            });
        }
        let mut bodies = Vec::with_capacity(self.bodies.len());
        for (body, &id) in self.bodies.iter().zip(&self.body_ids) {
            check()?;
            bodies.push(body.resolve(thermal.cells()[&id].specific_h_j_kg)?);
        }
        let mut contacts = Vec::with_capacity(self.contacts.len());
        for (index, (&contact, &(liquid_id, body_index))) in self
            .contacts
            .iter()
            .zip(&self.contact_endpoints)
            .enumerate()
        {
            check()?;
            let liquid = states[&liquid_id].clone();
            let body = bodies[body_index].clone();
            let interface = self
                .interfaces
                .get(&index)
                .map(|binding| {
                    binding.at([
                        self.caloric.state_at_specific_enthalpy(
                            thermal.cells()[&liquid_id].specific_h_j_kg,
                        )?,
                        self.bodies[body_index]
                            .caloric()
                            .state_at_specific_enthalpy(
                                thermal.cells()[&self.body_ids[body_index]].specific_h_j_kg,
                            )?,
                    ])
                })
                .transpose()?;
            let resistance = interface.as_ref().map_or(0.0, |source| {
                source
                    .property(interface::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY)
                    .expect("required interface resistance")
                    .value_si()
            });
            let conductance_w_k = series_contact_conductance_w_k(
                contact.area_m2,
                contact.distances_m,
                [
                    liquid.thermal_conductivity_w_m_k(),
                    self.bodies[body_index].conductivity(&body),
                ],
                resistance,
            )?;
            contacts.push(FreeSurfaceBodyContactSample {
                contact,
                liquid,
                body,
                interface,
                conductance_w_k,
            });
        }
        Ok(Samples {
            liquid: samples,
            contacts,
            bodies,
        })
    }
}

impl DepositionObserver for Conduction<'_, '_> {
    type Error = Error;

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Error> {
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Error> {
        check()?;
        let mut powers = self.powers.clone();
        if let Some(current) = self.sample_current(thermal, check)? {
            let electrical = self.electrical.as_mut().expect("bound electrical core");
            for (&id, &watts) in electrical.cells.iter().zip(&current.solution.joule_w) {
                check()?;
                let power = powers.entry(id).or_insert(0.0);
                *power = electrical::finite(*power + watts)?;
            }
            electrical.last = Some(current);
        }
        self.environment_last = self.sample_environment(thermal, check)?;
        if self.environment_reference.is_empty() {
            self.environment_reference = self.environment_last.clone();
        }
        for (&id, sample) in self
            .environment_endpoints
            .iter()
            .zip(&self.environment_last)
        {
            check()?;
            if sample.heat_into_liquid_w == 0.0 {
                continue;
            }
            let watts = powers.entry(id).or_insert(0.0);
            *watts = environment::finite(*watts + sample.heat_into_liquid_w)?;
        }
        Ok(powers)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Error> {
        let samples = self.sample(thermal, check)?;
        self.last = samples.liquid;
        self.body_last = samples.contacts;
        let mut values: BTreeMap<_, _> = self
            .last
            .iter()
            .enumerate()
            .map(|(i, s)| (i as u64, s.conductance_w_k))
            .collect();
        for (i, sample) in self.body_last.iter().enumerate() {
            values.insert((self.faces.len() + i) as u64, sample.conductance_w_k);
        }
        Ok(values)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Error> {
        self.environment_timestep(thermal, proposed_s, check)
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        // Admit all accepted thermal substeps, including isolated heated cells.
        let samples = self.sample(thermal, check)?;
        self.endpoint = samples.liquid;
        self.body_endpoint = samples.contacts;
        self.body_states_endpoint = samples.bodies;
        if let Some(endpoint) = self.sample_current(thermal, check)? {
            self.electrical.as_mut().expect("bound electrical core")
                .accept(start_s, duration_s, endpoint)?;
        }
        let mut gas_heat = 0.0;
        for sample in &self.environment_last {
            check()?;
            let heat = environment::finite(duration_s * sample.heat_into_liquid_w)?;
            if self.finite_atmosphere.is_some() {
                gas_heat = environment::finite(gas_heat + heat)?;
            }
            self.environment_heat_j = environment::finite(self.environment_heat_j + heat)?;
            self.environment_absolute_heat_j =
                environment::finite(self.environment_absolute_heat_j + heat.abs())?;
        }
        if let Some(gas) = &mut self.finite_atmosphere {
            gas.accept_heat(-gas_heat, check)?;
        }
        // Both film endpoints use the actually accepted gas T and EOS pressure.
        self.environment_endpoint = self.sample_environment(thermal, check)?;
        Ok(())
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Conduct between explicitly declared wetted faces, then advance the
    /// original native flow and extensive-enthalpy transport over one interval.
    ///
    /// The same card, selected claims and caloric chart own every substep.
    /// Face geometry is held at the accepted beginning-of-step configuration;
    /// submit the actual new wetted geometry at the next call. Each face may
    /// occur once, its area cannot exceed dx^2, and both endpoints must contain
    /// positive native liquid mass. Undeclared faces are explicitly insulated.
    /// No ambient or gas-interface exchange is inferred from unused face area.
    ///
    /// `heat` is external energy supplied uniformly over this interval, unlike
    /// the instantaneous beginning-of-step heat in `step_with_heat`. Internal
    /// conduction never enters its external-work account. Zero-mass interface
    /// ghosts and global conversion carry cannot receive face heat.
    ///
    /// This is a first-order thermal/flow split with explicit, CFL-subcycled
    /// thermal conduction. It retains the original constant-density, single-
    /// liquid and source-constant capillarity limits. It does not add surface energy,
    /// viscous heating, phase conversion or pressure work. All histories remain
    /// unchanged on thermal, source, VOF, transport, budget or cancellation error.
    pub fn step_with_conduction(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        self.advance_with_conduction(faces, heat, controls, &[], gate)
    }

    pub(super) fn advance_with_conduction(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        viscous: &[viscous::FreeSurfaceViscousCell],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        self.advance_with_thermal_bodies(faces, heat, controls, viscous, &mut [], &[], &[], gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_thermal_bodies(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        viscous: &[viscous::FreeSurfaceViscousCell],
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        self.advance_with_interfaces_and_environment(
            faces,
            heat,
            controls,
            viscous,
            bodies,
            contacts,
            interfaces,
            &[],
            None,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_environment(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        viscous: &[viscous::FreeSurfaceViscousCell],
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        environment: &[FreeSurfaceGasContact<'_>],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        self.advance_with_interfaces_and_environment(
            faces,
            heat,
            controls,
            viscous,
            bodies,
            contacts,
            &[],
            environment,
            None,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_interfaces_and_environment(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        viscous: &[viscous::FreeSurfaceViscousCell],
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        environment: &[FreeSurfaceGasContact<'_>],
        finite_atmosphere: Option<&mut atmosphere::FiniteAtmosphere>,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        poll(gate)?;
        if controls.maximum_faces > self.controls.max_mass_routes
            || faces
                .len()
                .saturating_add(contacts.len())
                .saturating_add(environment.len())
                > controls.maximum_faces
            || bodies.len() > self.controls.max_cells
            || heat.len() > self.controls.max_cells
            || controls.maximum_substeps == 0
            || controls.maximum_substeps > 1_000_000
            || !controls.maximum_dt_s.is_finite()
            || controls.maximum_dt_s <= 0.0
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j < 0.0
        {
            return Err(Error::Invalid(
                "invalid liquid-conduction work or time/energy bounds",
            ));
        }
        let mut interface_contacts = BTreeMap::new();
        if interfaces.len() > contacts.len() {
            return Err(Error::Invalid(
                "more interface laws than liquid/solid contacts",
            ));
        }
        for binding in interfaces {
            poll(gate)?;
            if binding.contact_index() >= contacts.len()
                || interface_contacts
                    .insert(binding.contact_index(), binding)
                    .is_some()
            {
                return Err(Error::Invalid(
                    "each interface law needs one unique existing contact",
                ));
            }
        }
        let dx = self.flow.scaling.cell_size_m();
        let area = dx * dx;
        let surface = 6.0 * area;
        if !area.is_finite() || area <= 0.0 || !surface.is_finite() {
            return Err(Error::Invalid(
                "liquid conduction face geometry is unrepresentable",
            ));
        }
        let duration = self.flow.scaling.timestep_s();
        let rho = self.flow.scaling.material().density_kg_m3();
        let mass_unit = self.flow.scaling.mass_unit_kg();
        let masses = self.flow.lattice.mass_snapshot();
        admit_masses(&masses, self.controls.max_cells)?;
        let mut ids = BTreeMap::new();
        let mut events = Vec::new();
        for (&location, &native_mass) in &masses {
            poll(gate)?;
            let Some(p) = FreeSurface3::transport_coordinates(location) else {
                continue;
            };
            if native_mass == 0.0 {
                continue;
            }
            let id = CellId(ids.len() as u64);
            ids.insert(location, id);
            let center = p.map(|c| (c as f64 + 0.5) * dx);
            let mass = native_mass * mass_unit;
            let mut start = center;
            let mut end = center;
            start[0] -= 0.5 * dx;
            end[0] += 0.5 * dx;
            let temperature = self
                .caloric
                .state_at_specific_enthalpy(self.h[&location])?
                .temperature_k();
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: 0,
                    material_id: MaterialId(0),
                    mass_kg: mass,
                    specific_h_j_kg: self.h[&location],
                    geometry: CellGeometry {
                        start_m: start,
                        end_m: end,
                        volume_m3: mass / rho,
                        surface_area_m2: surface,
                    },
                    boundary: Boundary {
                        ambient_k: temperature,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
        }
        if ids.is_empty() {
            return Err(Error::Invalid(
                "liquid conduction needs positive located cell mass",
            ));
        }
        if ids.len().saturating_add(bodies.len()) > self.controls.max_cells {
            return Err(Error::Invalid(
                "combined solid/liquid thermal cell budget exhausted",
            ));
        }
        let liquid_cells: Vec<_> = ids.values().copied().collect();
        let mut materials = vec![(MaterialId(0), self.caloric as &dyn ThermalLaw)];
        let mut body_by_id = BTreeMap::new();
        let mut body_ids = Vec::with_capacity(bodies.len());
        let mut body_before = 0.0;
        for (i, body) in bodies.iter().enumerate() {
            poll(gate)?;
            if body_by_id.insert(body.id(), i).is_some() {
                return Err(Error::Invalid("finite thermal body ids must be unique"));
            }
            body.resolve(body.specific_enthalpy_j_kg())?;
            let id = CellId((liquid_cells.len() + i) as u64);
            let material = MaterialId((i + 1) as u64);
            body_ids.push(id);
            materials.push((material, body.caloric() as &dyn ThermalLaw));
            body_before += body.enthalpy_j();
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(body.cell(id, material)),
            });
        }
        if !body_before.is_finite() {
            return Err(Error::Invalid("finite solid enthalpy sum overflows"));
        }
        let cell_id = |p| -> Result<CellId, Error> {
            let location = self
                .location_at(p)
                .ok_or(Error::Invalid("conduction cell is not occupied"))?;
            ids.get(&location).copied().ok_or(Error::Invalid(
                "conduction cell has no positive liquid mass",
            ))
        };
        let mut pairs = BTreeSet::new();
        let mut endpoints = Vec::with_capacity(faces.len());
        for (index, &face) in faces.iter().enumerate() {
            poll(gate)?;
            let a = cell_id(face.cells[0])?;
            let b = cell_id(face.cells[1])?;
            let separation = (0..3).fold(0usize, |s, i| {
                s.saturating_add(face.cells[0][i].abs_diff(face.cells[1][i]))
            });
            if separation != 1
                || !face.area_m2.is_finite()
                || face.area_m2 < 0.0
                || face.area_m2 > area
                || face
                    .distances_m
                    .iter()
                    .any(|d| !d.is_finite() || *d <= 0.0 || *d > dx)
                || !pairs.insert((a.min(b), a.max(b)))
            {
                return Err(Error::Invalid(
                    "wetted face must be unique, adjacent and fit its Cartesian cells",
                ));
            }
            endpoints.push([a, b]);
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: index as u64,
                    a,
                    b,
                    area_m2: face.area_m2,
                    conductance_w_k: 0.0,
                }),
            });
        }
        let mut wall_faces = BTreeSet::new();
        let mut contact_endpoints = Vec::with_capacity(contacts.len());
        for (index, &contact) in contacts.iter().enumerate() {
            poll(gate)?;
            let a = cell_id(contact.liquid)?;
            let body = *body_by_id
                .get(&contact.body_id)
                .ok_or(Error::Invalid("liquid contact names an absent finite body"))?;
            let separation = (0..3).fold(0u64, |sum, axis| {
                sum.saturating_add((contact.liquid[axis] as i64).abs_diff(contact.wall[axis]))
            });
            if separation != 1
                || !self.flow.lattice.transport_wall_at(contact.wall)
                || !wall_faces.insert((a, contact.wall))
                || !contact.area_m2.is_finite()
                || contact.area_m2 < 0.0
                || contact.area_m2 > area
                || contact
                    .distances_m
                    .iter()
                    .any(|d| !d.is_finite() || *d <= 0.0)
                || contact.distances_m[0] > dx
            {
                return Err(Error::Invalid(
                    "solid contact needs one unique adjacent permanent wetted wall face",
                ));
            }
            contact_endpoints.push((a, body));
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: (faces.len() + index) as u64,
                    a,
                    b: body_ids[body],
                    area_m2: contact.area_m2,
                    conductance_w_k: 0.0,
                }),
            });
        }
        let environment_endpoints = self.admit_environment(environment, &cell_id, dx, gate)?;
        let mut powers = BTreeMap::new();
        let mut absolute_heat = 0.0;
        let mut external_heat = 0.0;
        for source in heat {
            poll(gate)?;
            let id = cell_id(source.coordinates)?;
            let power = source.energy_j / duration;
            if !source.energy_j.is_finite()
                || !power.is_finite()
                || (source.energy_j != 0.0 && power == 0.0)
                || powers.insert(id, power).is_some()
            {
                return Err(Error::Invalid(
                    "conduction heat must be finite, unique and representable as power",
                ));
            }
            absolute_heat += source.energy_j.abs();
            external_heat += source.energy_j;
        }
        let mut viscous_heat = 0.0;
        for source in viscous {
            poll(gate)?;
            let id = cell_id(source.coordinates)?;
            if source.heat_j != 0.0 {
                let power = source.heat_j / duration;
                let combined = powers.entry(id).or_insert(0.0);
                *combined += power;
                if !power.is_finite() || power == 0.0 || !combined.is_finite() {
                    return Err(Error::Invalid("liquid molecular heat power overflows"));
                }
            }
            viscous_heat += source.heat_j;
        }
        absolute_heat += viscous_heat;
        if !absolute_heat.is_finite() {
            return Err(Error::Invalid("conduction heat magnitude overflows"));
        }
        let count = events.len();
        let config = StepConfig {
            max_dt_s: controls.maximum_dt_s,
            energy_tolerance_j: controls.energy_tolerance_j,
        };
        let thermal = DepositionSimulation::new(
            events,
            &materials,
            DepositionBudget {
                max_cells: self.controls.max_cells,
                max_links: controls.maximum_faces,
                max_events: count,
                max_substeps: controls.maximum_substeps,
            },
        )?
        .advance_to(0.0, config, || checkpoint(gate))?
        .0;
        let electrical = self.electrical.as_ref().map(|core| {
            core.admit(&self.flow.lattice, &self.states, gate)?;
            electrical::ElectricalTimeline::new(core, &cell_id, &self.flow)
        }).transpose()?;
        let rheology = self.rheology_model();
        let mut rheology_stress = BTreeMap::new();
        if rheology.is_some() {
            for (&location, &id) in &ids {
                poll(gate)?;
                rheology_stress.insert(id, self.rheology_stress(&self.flow.lattice, location)?);
            }
        }
        let observer = Conduction {
            gate,
            law: self.law,
            caloric: self.caloric,
            density: rho,
            buoyancy: self.buoyancy,
            rheology,
            rheology_stress,
            faces,
            endpoints,
            liquid_cells,
            liquid_coordinates: ids.iter().filter_map(|(&location, &id)| {
                FreeSurface3::transport_coordinates(location).map(|position| (id, position))
            }).collect(),
            wetting: self.wetting.clone(),
            bodies,
            body_ids: body_ids.clone(),
            contacts,
            interfaces: interface_contacts,
            contact_endpoints,
            powers,
            last: Vec::new(),
            endpoint: Vec::new(),
            body_last: Vec::new(),
            body_endpoint: Vec::new(),
            body_states_endpoint: Vec::new(),
            environment,
            environment_endpoints,
            environment_last: Vec::new(),
            environment_endpoint: Vec::new(),
            environment_reference: Vec::new(),
            environment_heat_j: 0.0,
            environment_absolute_heat_j: 0.0,
            finite_atmosphere: finite_atmosphere.as_deref().cloned(),
            electrical,
        };
        let (thermal, mut observer, audit) =
            thermal.advance_observed(&observer, duration, config, || checkpoint(gate))?;
        let electrical_report = observer.electrical.as_ref()
            .map(|current| current.report(controls.energy_tolerance_j)).transpose()?;
        let electrical_forces = observer.electrical.as_ref()
            .map(|current| current.average_forces(&self.flow.scaling)).transpose()?.flatten();
        let joule_heat = electrical_report.as_ref().map_or(0.0, |current| current.joule_heat_j);
        let mut finite_gas_report = observer.finite_atmosphere.as_ref()
            .map(|gas| gas.balance(observer.environment_heat_j, &audit, controls.energy_tolerance_j))
            .transpose()?;
        let mut heated = self.h.clone();
        for (location, id) in ids {
            poll(gate)?;
            heated.insert(location, thermal.cells()[&id].specific_h_j_kg);
        }
        let external = if viscous.is_empty() && environment.is_empty() && electrical_report.is_none() {
            audit.external_heat_j
        } else {
            let defect =
                audit.external_heat_j - external_heat - viscous_heat - observer.environment_heat_j
                    - joule_heat;
            if !defect.is_finite() || defect.abs() > controls.energy_tolerance_j {
                return Err(Error::Invalid(
                    "authored, molecular, electrical and environmental liquid heat do not close",
                ));
            }
            external_heat + observer.environment_heat_j + joule_heat
        };
        let mut body_after = 0.0;
        let mut body_change = 0.0;
        let mut body_enthalpies = Vec::with_capacity(bodies.len());
        for (body, &id) in bodies.iter().zip(&body_ids) {
            poll(gate)?;
            let h = thermal.cells()[&id].specific_h_j_kg;
            body_after += body.mass_kg() * h;
            body_change += body.mass_kg() * (h - body.specific_enthalpy_j_kg());
            body_enthalpies.push(h);
        }
        if !body_after.is_finite() || !body_change.is_finite() {
            return Err(Error::Invalid("finite solid endpoint enthalpy overflows"));
        }
        let body_heat = if bodies.is_empty() { 0.0 } else { -body_change };
        let body_sources = observer.body_states_endpoint;
        if body_sources.len() != bodies.len() {
            return Err(Error::Invalid(
                "finite solid endpoint source state is incomplete",
            ));
        }
        let body_last = observer.body_last;
        let body_endpoint = observer.body_endpoint;
        let last_applied = observer.last;
        let thermal_endpoint = observer.endpoint;
        let environment_report = (!environment.is_empty()).then_some(FreeSurfaceGasExchange {
            heat_to_liquid_j: observer.environment_heat_j,
            last_applied: observer.environment_last,
            thermal_endpoint: observer.environment_endpoint,
        });
        let flow = self.advance_prepared_enthalpy(
            heated,
            external,
            viscous_heat,
            body_heat,
            absolute_heat + body_heat.abs() + observer.environment_absolute_heat_j + joule_heat,
            electrical_forces.as_ref(),
            gate,
        )?;
        if let Some(report) = &mut finite_gas_report {
            // Native advection already checked its own extensive-enthalpy
            // defect. This is the combined actual liquid/wall/gas account;
            // internal gas and wall transfers cancel, molecular heat stays
            // in its original separate account.
            report.combined_thermal_balance_defect_j =
                flow.balance_defect_j + report.gas_balance_defect_j;
        }
        let accepted_finite_atmosphere = observer.finite_atmosphere.take();
        // Flow has passed every endpoint/advection gate; accepting the already
        // admitted solid candidate is now infallible. No partial thermal commit.
        for ((body, h), material) in bodies.iter_mut().zip(body_enthalpies).zip(body_sources) {
            body.accept(h, material);
        }
        if let (Some(target), Some(accepted)) = (finite_atmosphere, accepted_finite_atmosphere) {
            *target = accepted;
        }
        let body_report = (!bodies.is_empty()).then_some(FreeSurfaceBodyExchange {
            enthalpy_before_j: body_before,
            enthalpy_after_j: body_after,
            heat_to_liquid_j: body_heat,
            combined_balance_defect_j: flow.balance_defect_j,
            last_applied: body_last,
            thermal_endpoint: body_endpoint,
        });
        Ok(FreeSurfaceConductionReport {
            thermal: audit,
            flow,
            last_applied,
            thermal_endpoint,
            bodies: body_report,
            environment: environment_report,
            finite_gas: finite_gas_report,
            electrical: electrical_report,
        })
    }
}
