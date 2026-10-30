//! File-supplied material states lowered into the existing coupled solvers.
//! This is CLI orchestration, not a new constitutive law or discretization.

pub(super) mod transient;
mod conductivity;
mod fluid;
pub(super) use fluid::{Fluid, load_fluid};
pub(super) use fluid::NO_CLAIM as FLUID_NO_CLAIM;

use fs_exec::CancelGate;
use fs_lbm::conjugate::material_coupling::{
    MaterialCoupledSolution, MaterialCouplingConfig,
};
use fs_lbm::conjugate::material_coupling::interfaces::{
    MaterialContact, solve as solve_material_interfaces,
};
use fs_lbm::conjugate::{
    ChtError, ElectricSetup, EnergyConfig, EnergySolution, FlowField,
    FluidProperties, SolidMaterial, ThermalSetup, VoxelDomain,
};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::conductor::{ResolvedConductorState, resolve_conductor_state};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::interfaces::{self, InterfaceSource, ThermalInterfaceState};
use super::radiation::{FaceReport, Radiation};
use fs_lbm::conjugate::boundary_exchange::PreparedBoundaryExchange;
use conductivity::Conductivity;

pub(super) struct Source {
    pub pack: NormalizedMaterialCardPack,
    point: QueryPoint,
    pub(super) temperature_axis: String,
    temperature_quantity: Option<QuantitySpec>,
    pub conducting: bool,
    pub(super) conductivity: Conductivity,
}

pub(super) struct CellState {
    pub thermal: ResolvedMaterialStatePoint,
    pub conductor: Option<ResolvedConductorState>,
    pub k: [f64; 3],
}

impl Source {
    #[cfg(test)]
    pub fn new(
        pack: NormalizedMaterialCardPack,
        point: QueryPoint,
        temperature_axis: String,
        temperature_quantity: Option<QuantitySpec>,
        conducting: bool,
    ) -> Result<Self, ChtError> {
        Self::with_conductivity(pack, point, temperature_axis, temperature_quantity, conducting, None)
    }

    pub fn with_conductivity(
        pack: NormalizedMaterialCardPack,
        point: QueryPoint,
        temperature_axis: String,
        temperature_quantity: Option<QuantitySpec>,
        conducting: bool,
        declaration: Option<&super::J>,
    ) -> Result<Self, ChtError> {
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        if temperature_axis.trim().is_empty()
            || point.axes().contains_key(&temperature_axis)
            || temperature_quantity.is_some_and(|q| {
                q != absolute && q != QuantitySpec::dimensional(Temperature::DIMS)
            })
        {
            return Err(invalid("material.temperature_axis", "declare one separate absolute kelvin axis"));
        }
        let conductivity = Conductivity::load(&pack, declaration)?;
        Ok(Self { pack, point, temperature_axis, temperature_quantity, conducting, conductivity })
    }

    pub(super) fn point_at(&self, temperature: f64) -> Result<QueryPoint, ChtError> {
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(invalid("material.temperature", "temperature must be finite positive kelvin"));
        }
        match self.temperature_quantity {
            Some(quantity) => self.point.clone().with_quantity(
                &self.temperature_axis, quantity, temperature,
            ),
            None => self.point.clone().with(&self.temperature_axis, temperature),
        }.map_err(|e| invalid("material.state", e))
    }

    pub fn at(&self, temperature: f64) -> Result<CellState, ChtError> {
        let point = self.point_at(temperature)?;
        let thermal = resolve_material_state_point(
            self.pack.card(), &point, self.conductivity.requirements(),
            MaterialPropertySelection::SingleClaimOnly,
        ).map_err(|e| invalid("material.thermal_conductivity", e))?;
        let k = self.conductivity.grid_values(&thermal);
        let conductor = if self.conducting {
            Some(resolve_conductor_state(
                self.pack.card(), &point, MaterialPropertySelection::SingleClaimOnly,
            ).map_err(|e| invalid("material.electrical_resistivity", e))?)
        } else {
            None
        };
        Ok(CellState { thermal, conductor, k })
    }
}

pub(super) struct Problem {
    pub domain: VoxelDomain,
    pub material_indices: Vec<usize>,
    pub materials: Vec<Source>,
    pub electrical_interfaces: Vec<InterfaceSource>,
    pub thermal_interfaces: Vec<InterfaceSource>,
    pub thermal_contact_tolerance: f64,
    pub electric: ElectricSetup,
    pub thermal: ThermalSetup,
    pub initial_temperature: Vec<f64>,
    pub coupling: MaterialCouplingConfig,
    pub energy: EnergyConfig,
    pub max_thermal_passes: usize,
    pub conductivity_tolerance: f64,
    pub fluid: Option<Fluid>,
    pub radiation: Option<Radiation>,
}

pub(super) struct ResultField {
    pub coupled: MaterialCoupledSolution,
    pub sources: Vec<CellState>,
    pub interfaces: Vec<MaterialContact>,
    pub thermal_interfaces: Vec<ThermalInterfaceState>,
    pub thermal_contact_relative_defect: f64,
    pub thermal_solves: usize,
    pub conductivity_relative_defect: f64,
    pub radiation: Vec<FaceReport>,
}

impl Problem {
    fn states(&self, temperatures: &[f64], gate: &CancelGate) -> Result<Vec<CellState>, ChtError> {
        if temperatures.len() != self.domain.cell_count() {
            return Err(invalid("material.temperature", "one temperature per cell is required"));
        }
        let mut states = Vec::with_capacity(temperatures.len());
        for (cell, &t) in temperatures.iter().enumerate() {
            if cell.is_multiple_of(256) { poll(gate)?; }
            states.push(self.materials[self.material_indices[cell]].at(t)?);
        }
        poll(gate)?;
        Ok(states)
    }

    pub fn solve(&mut self, gate: &CancelGate) -> Result<ResultField, ChtError> {
        self.solve_at(&self.initial_temperature, &self.electric, &self.thermal, None, gate)
    }

    fn solve_at(
        &self, start: &[f64], electric: &ElectricSetup, thermal: &ThermalSetup,
        storage: Option<&transient::Storage<'_>>, gate: &CancelGate,
    ) -> Result<ResultField, ChtError> {
        poll(gate)?;
        let initial = self.states(start, gate)?;
        let mut solids = Vec::with_capacity(self.materials.len());
        for (m, source) in self.materials.iter().enumerate() {
            let cell = self.material_indices.iter().position(|&i| i == m)
                .ok_or_else(|| invalid("material_indices", "each declared material must be assigned to a cell"))?;
            solids.push(SolidMaterial::new(source.pack.pack_id(), initial[cell].k[0])
                .with_orthotropic(initial[cell].k));
        }
        // Only all-solid scenes use this unconsumed placeholder. Fluid scenes
        // use their exact source-resolved reference and prescribed face flux.
        let unused_fluid = FluidProperties {
            density_kg_m3: 1.0, specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: 1.0, kinematic_viscosity_m2_s: 1.0,
        };
        let quiescent = self.fluid.is_none().then(|| FlowField::quiescent(&self.domain));
        let (fluid, flow) = match &self.fluid {
            Some(binding) => (&binding.properties, &binding.flow),
            None => (&unused_fluid, quiescent.as_ref().expect("all-solid quiescent field")),
        };
        let mut thermal_solves = 0usize;
        let mut thermal_interfaces = Vec::new();
        let mut thermal_contact_relative_defect = 0.0;
        let mut radiation = Vec::new();
        let mut surface_temperatures = self.radiation.as_ref()
            .map_or_else(Vec::new, Radiation::initial_temperatures);
        // Electrical and thermal laws share physical faces/state coordinates,
        // but neither response is inferred from the other. Each owner retains
        // the exact final source data and its own convergence defect.
        let interface_solution = solve_material_interfaces(
            &self.domain, electric, thermal, start, self.coupling,
            |cell, _, temperature| {
                Ok(self.materials[self.material_indices[cell]].at(temperature)?.conductor
                    .map(|c| c.resistivity_ohm_m()))
            },
            |temperature| interfaces::resolve(&self.electrical_interfaces, temperature, gate),
            |heated, temperature| {
                let (energy, joints, defect, surfaces) = self.solve_heat(
                    heated, temperature, &solids, fluid, flow, storage, &mut thermal_solves,
                    &mut surface_temperatures, gate,
                )?;
                thermal_interfaces = joints;
                radiation = surfaces;
                thermal_contact_relative_defect = defect;
                Ok(energy)
            },
            gate,
        )?;
        let coupled = interface_solution.coupled;
        let interfaces = interface_solution.contacts;
        let sources = self.states(&coupled.thermal.temperature, gate)?;
        let conductivity_relative_defect = conductivity_defect(&sources, &coupled.thermal.conductivity);
        if conductivity_relative_defect > self.conductivity_tolerance {
            return Err(invalid("material.conductivity_consistency", "final conductivity differs from the thermal solve beyond tolerance"));
        }
        poll(gate)?;
        Ok(ResultField { coupled, sources, interfaces, thermal_interfaces,
            thermal_contact_relative_defect, thermal_solves, conductivity_relative_defect, radiation })
    }

    #[allow(clippy::too_many_arguments)] // existing physical owners and the accumulated work count
    fn solve_heat(
        &self, setup: &ThermalSetup, start: &[f64], solids: &[SolidMaterial],
        fluid: &FluidProperties, flow: &FlowField, storage: Option<&transient::Storage<'_>>,
        solves: &mut usize, surface_temperatures: &mut [f64], gate: &CancelGate,
    ) -> Result<(EnergySolution, Vec<ThermalInterfaceState>, f64, Vec<FaceReport>), ChtError> {
        let mut previous = start.to_vec();
        let mut states = self.states(&previous, gate)?;
        let mut configured = setup.clone();
        let mut defect = f64::INFINITY;
        let mut radiation_defect = 0.0;
        let passes = self.radiation.as_ref().map_or(self.max_thermal_passes,
            |r| self.max_thermal_passes.min(r.config.max_iterations));
        for _ in 0..passes {
            poll(gate)?;
            let field: Vec<[f64; 3]> = states.iter().map(|state| state.k).collect();
            let joints = interfaces::resolve_thermal(&self.thermal_interfaces, &previous, gate)?;
            configured.face_contacts.clone_from(&setup.face_contacts);
            configured.face_contacts.extend(joints.iter().map(ThermalInterfaceState::face));
            let boundary = self.radiation.as_ref().map(|r|
                r.prepare(self, &configured, &field, surface_temperatures, gate)).transpose()?;
            let with_boundary = boundary.as_ref().map_or(&configured, PreparedBoundaryExchange::setup);
            let capacity = storage.map(|s| s.capacities(&previous, gate)).transpose()?;
            let advection = storage.map(|s| s.advection(self, &previous, with_boundary, gate))
                .transpose()?.flatten();
            let solution = match (storage, &capacity, &advection) {
                (Some(s), Some(capacity), Some(advection)) => {
                    // Variable-Cp transport is conservative upwind enthalpy;
                    // the reference-Cp Patankar blend is not its constitutive law.
                    let energy = EnergyConfig { scheme: fs_lbm::conjugate::ConvectionScheme::Upwind, ..self.energy };
                    energy.step_with_enthalpy(
                        &self.domain, fluid, solids, flow, with_boundary, &field,
                        advection.field(), capacity, s.previous, s.dt, gate,
                    )?.0
                }
                (Some(s), Some(capacity), None) => self.energy.step_with_conductivity(
                    &self.domain, fluid, solids, flow, with_boundary, &field,
                    capacity, s.previous, s.dt, gate,
                )?.0,
                _ => self.energy.solve_with_conductivity(
                    &self.domain, fluid, solids, flow, with_boundary, &field, gate,
                )?,
            };
            *solves += 1;
            let (radiation, surface_defect) = match (&self.radiation, &boundary) {
                (Some(r), Some(b)) => r.audit(self, b, &solution, surface_temperatures, gate)?,
                _ => (Vec::new(), 0.0),
            };
            radiation_defect = surface_defect;
            // Re-query every thermal response at the proposed endpoint. Small
            // temperature change cannot mask a large change in joint resistance.
            let next = self.states(&solution.temperature, gate)?;
            let next_joints = interfaces::resolve_thermal(&self.thermal_interfaces, &solution.temperature, gate)?;
            let contact_defect = interfaces::thermal_defect(&next_joints, &solution)?;
            let storage_consistent = match (storage, &capacity) {
                (Some(s), Some(capacity)) => s.consistent(&solution.temperature, capacity, gate)?,
                _ => true,
            };
            let advection_consistent = match (storage, &advection) {
                (Some(s), Some(a)) => s.advection_consistent(self, &solution.temperature, a, gate)?,
                _ => true,
            };
            let conductivity_defect = conductivity_defect(&next, &solution.conductivity);
            defect = conductivity_defect.max(contact_defect);
            let change = previous.iter().zip(&solution.temperature)
                .map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            if conductivity_defect <= self.conductivity_tolerance
                && contact_defect <= self.thermal_contact_tolerance
                && change <= self.coupling.temperature_tolerance_k && storage_consistent
                && advection_consistent && radiation_defect <= 1.0
            {
                for (temperature, face) in surface_temperatures.iter_mut().zip(&radiation) {
                    *temperature = face.surface_temperature_k;
                }
                return Ok((solution, next_joints, contact_defect, radiation));
            }
            if let Some(r) = &self.radiation {
                for (temperature, face) in surface_temperatures.iter_mut().zip(&radiation) {
                    *temperature = r.config.relaxation.mul_add(face.surface_temperature_k - *temperature, *temperature);
                }
            }
            previous = solution.temperature;
            states = next;
        }
        if radiation_defect > 1.0 {
            return Err(ChtError::SolverNotConverged {
                system: "source-resolved surface radiation", iterations: passes,
                relative_residual: radiation_defect, tolerance: 1.0,
            });
        }
        Err(ChtError::SolverNotConverged {
            system: "source-resolved thermal conductivity/contact/storage/advection", iterations: passes,
            relative_residual: defect, tolerance: self.conductivity_tolerance.min(self.thermal_contact_tolerance),
        })
    }
}

fn conductivity_defect(states: &[CellState], field: &[[f64; 3]]) -> f64 {
    if states.len() != field.len() { return f64::INFINITY; }
    states.iter().zip(field).flat_map(|(state, axes)| {
        state.k.iter().zip(axes).map(|(&resolved, &solved)|
            (resolved - solved).abs() / resolved.max(solved))
    }).fold(0.0, f64::max)
}

pub(super) fn invalid(field: &'static str, reason: impl core::fmt::Display) -> ChtError {
    ChtError::InvalidInput { field, reason: reason.to_string() }
}

pub(super) fn poll(gate: &CancelGate) -> Result<(), ChtError> {
    if gate.is_requested() { Err(ChtError::Cancelled) } else { Ok(()) }
}
