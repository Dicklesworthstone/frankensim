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
    ThermalLaw, ThermalLink, TimedEvent, series_conductance_w_k,
};

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
    law: &'a LiquidTemperatureLaw<'m>,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    density: f64,
    buoyancy: Option<BoussinesqLiquid>,
    faces: &'a [FreeSurfaceConductionFace],
    endpoints: Vec<[CellId; 2]>,
    powers: BTreeMap<CellId, f64>,
    last: Vec<FreeSurfaceConductionSample>,
    endpoint: Vec<FreeSurfaceConductionSample>,
}

impl Conduction<'_, '_> {
    fn sample(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreeSurfaceConductionSample>, Error> {
        let mut states = BTreeMap::new();
        for (&id, cell) in thermal.cells() {
            check()?;
            let (state, _) = resolve_liquid(
                self.law,
                self.caloric,
                cell.specific_h_j_kg,
                self.density,
                self.buoyancy,
            )?;
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
        Ok(samples)
    }
}

impl DepositionObserver for Conduction<'_, '_> {
    type Error = Error;

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Error> {
        Ok(())
    }

    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Error> {
        check()?;
        Ok(self.powers.clone())
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Error> {
        self.last = self.sample(thermal, check)?;
        Ok(self
            .last
            .iter()
            .enumerate()
            .map(|(i, s)| (i as u64, s.conductance_w_k))
            .collect())
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        // Admit all accepted thermal substeps, including isolated heated cells.
        self.endpoint = self.sample(thermal, check)?;
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
    /// liquid and neglected-capillarity limits. It does not add surface energy,
    /// viscous heating, phase conversion or pressure work. All histories remain
    /// unchanged on thermal, source, VOF, transport, budget or cancellation error.
    pub fn step_with_conduction(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        poll(gate)?;
        if controls.maximum_faces > self.controls.max_mass_routes
            || faces.len() > controls.maximum_faces
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
        let mut powers = BTreeMap::new();
        let mut absolute_heat = 0.0;
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
        }
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
            &[(MaterialId(0), self.caloric as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: self.controls.max_cells,
                max_links: controls.maximum_faces,
                max_events: count,
                max_substeps: controls.maximum_substeps,
            },
        )?
        .advance_to(0.0, config, || checkpoint(gate))?
        .0;
        let observer = Conduction {
            law: self.law,
            caloric: self.caloric,
            density: rho,
            buoyancy: self.buoyancy,
            faces,
            endpoints,
            powers,
            last: Vec::new(),
            endpoint: Vec::new(),
        };
        let (thermal, observer, audit) =
            thermal.advance_observed(&observer, duration, config, || checkpoint(gate))?;
        let mut heated = self.h.clone();
        for (location, id) in ids {
            poll(gate)?;
            heated.insert(location, thermal.cells()[&id].specific_h_j_kg);
        }
        let flow =
            self.advance_prepared_enthalpy(heated, audit.external_heat_j, absolute_heat, gate)?;
        Ok(FreeSurfaceConductionReport {
            thermal: audit,
            flow,
            last_applied: observer.last,
            thermal_endpoint: observer.endpoint,
        })
    }
}
