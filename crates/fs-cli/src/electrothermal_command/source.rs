//! File-supplied material states lowered into the existing coupled solvers.
//! This is CLI orchestration, not a new constitutive law or discretization.

use fs_exec::CancelGate;
use fs_lbm::conjugate::material_coupling::{
    MaterialCoupledSolution, MaterialCouplingConfig, solve_material_electrothermal,
};
use fs_lbm::conjugate::{
    ChtError, Conductor, ElectricSetup, EnergyConfig, EnergySolution, FlowField,
    FluidProperties, SolidMaterial, ThermalSetup, VoxelDomain,
};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::conductor::{ResolvedConductorState, resolve_conductor_state};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

pub(super) struct Source {
    pub pack: NormalizedMaterialCardPack,
    point: QueryPoint,
    temperature_axis: String,
    temperature_quantity: Option<QuantitySpec>,
    pub conducting: bool,
    conductivity: ScalarPropertyRequirement,
}

pub(super) struct CellState {
    pub thermal: ResolvedMaterialStatePoint,
    pub conductor: Option<ResolvedConductorState>,
    pub k: f64,
}

impl Source {
    pub fn new(
        pack: NormalizedMaterialCardPack,
        point: QueryPoint,
        temperature_axis: String,
        temperature_quantity: Option<QuantitySpec>,
        conducting: bool,
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
        // The property name selects the physical response. Preserve its exact
        // dimensional/semantic schema; do not strip a source quantity kind.
        let thermal_kind = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity, ValueForm::Static,
        ));
        let candidates = pack.card().claims_for(THERMAL_CONDUCTIVITY_PROPERTY);
        let (_, first) = candidates.first().ok_or_else(|| {
            invalid("material.thermal_conductivity", "card has no thermal-conductivity claim")
        })?;
        let key = &first.key;
        if (key.quantity() != thermal_kind
            && key.quantity() != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS))
            || candidates.iter().any(|(_, claim)| claim.key != *key)
        {
            return Err(invalid("material.thermal_conductivity", "unsupported or conflicting conductivity quantity schemas"));
        }
        let conductivity = ScalarPropertyRequirement::try_with_key(
            key, ScalarAdmissibility::StrictlyPositive,
        ).map_err(|e| invalid("material.thermal_conductivity", e))?;
        Ok(Self { pack, point, temperature_axis, temperature_quantity, conducting, conductivity })
    }

    pub fn at(&self, temperature: f64) -> Result<CellState, ChtError> {
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(invalid("material.temperature", "temperature must be finite positive kelvin"));
        }
        let point = match self.temperature_quantity {
            Some(quantity) => self.point.clone().with_quantity(
                &self.temperature_axis, quantity, temperature,
            ),
            None => self.point.clone().with(&self.temperature_axis, temperature),
        }.map_err(|e| invalid("material.state", e))?;
        let thermal = resolve_material_state_point(
            self.pack.card(), &point, core::slice::from_ref(&self.conductivity),
            MaterialPropertySelection::SingleClaimOnly,
        ).map_err(|e| invalid("material.thermal_conductivity", e))?;
        let k = thermal.property(THERMAL_CONDUCTIVITY_PROPERTY)
            .expect("required conductivity was resolved").value_si();
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
    pub electric: ElectricSetup,
    pub thermal: ThermalSetup,
    pub initial_temperature: Vec<f64>,
    pub coupling: MaterialCouplingConfig,
    pub energy: EnergyConfig,
    pub max_thermal_passes: usize,
    pub conductivity_tolerance: f64,
}

pub(super) struct ResultField {
    pub coupled: MaterialCoupledSolution,
    pub sources: Vec<CellState>,
    pub thermal_solves: usize,
    pub conductivity_relative_defect: f64,
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
        poll(gate)?;
        let initial = self.states(&self.initial_temperature, gate)?;
        let mut solids = Vec::with_capacity(self.materials.len());
        self.electric.conductors.clear();
        for (m, source) in self.materials.iter().enumerate() {
            let cell = self.material_indices.iter().position(|&i| i == m)
                .ok_or_else(|| invalid("material_indices", "each declared material must be assigned to a cell"))?;
            solids.push(SolidMaterial::new(source.pack.pack_id(), initial[cell].k));
            self.electric.conductors.push(initial[cell].conductor.as_ref()
                .map(|c| Conductor::new(c.resistivity_ohm_m())));
        }
        // This command admits all-solid domains. No fluid cell consumes these
        // reference coefficients; they do not describe invented air properties.
        let unused_fluid = FluidProperties {
            density_kg_m3: 1.0, specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: 1.0, kinematic_viscosity_m2_s: 1.0,
        };
        let flow = FlowField::quiescent(&self.domain);
        let mut thermal_solves = 0usize;
        let coupled = solve_material_electrothermal(
            &self.domain, &self.electric, &self.thermal, &self.initial_temperature,
            self.coupling,
            |_, material, temperature| {
                self.materials[usize::from(material)].at(temperature)?.conductor
                    .map(|c| c.resistivity_ohm_m())
                    .ok_or_else(|| invalid("material.electrical_resistivity", "insulator requested as conductor"))
            },
            |heated, temperature| {
                self.solve_heat(heated, temperature, &solids, &unused_fluid, &flow,
                    &mut thermal_solves, gate)
            },
            gate,
        )?;
        let sources = self.states(&coupled.thermal.temperature, gate)?;
        let conductivity_relative_defect = conductivity_defect(&sources, &coupled.thermal.conductivity);
        if conductivity_relative_defect > self.conductivity_tolerance {
            return Err(invalid("material.conductivity_consistency", "final conductivity differs from the thermal solve beyond tolerance"));
        }
        poll(gate)?;
        Ok(ResultField { coupled, sources, thermal_solves, conductivity_relative_defect })
    }

    #[allow(clippy::too_many_arguments)] // existing physical owners and the accumulated work count
    fn solve_heat(
        &self, setup: &ThermalSetup, start: &[f64], solids: &[SolidMaterial],
        fluid: &FluidProperties, flow: &FlowField, solves: &mut usize, gate: &CancelGate,
    ) -> Result<EnergySolution, ChtError> {
        let mut previous = start.to_vec();
        let mut states = self.states(&previous, gate)?;
        let mut defect = f64::INFINITY;
        for _ in 0..self.max_thermal_passes {
            poll(gate)?;
            let field: Vec<[f64; 3]> = states.iter().map(|state| [state.k; 3]).collect();
            let solution = self.energy.solve_with_conductivity(
                &self.domain, fluid, solids, flow, setup, &field, gate,
            )?;
            *solves += 1;
            // This re-query is essential: the last temperature must be supported
            // by the actual source, not by a clamped/extended conductivity table.
            let next = self.states(&solution.temperature, gate)?;
            defect = conductivity_defect(&next, &solution.conductivity);
            let change = previous.iter().zip(&solution.temperature)
                .map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            if defect <= self.conductivity_tolerance && change <= self.coupling.temperature_tolerance_k {
                return Ok(solution);
            }
            previous = solution.temperature;
            states = next;
        }
        Err(ChtError::SolverNotConverged {
            system: "source-resolved thermal conductivity", iterations: self.max_thermal_passes,
            relative_residual: defect, tolerance: self.conductivity_tolerance,
        })
    }
}

fn conductivity_defect(states: &[CellState], field: &[[f64; 3]]) -> f64 {
    if states.len() != field.len() { return f64::INFINITY; }
    states.iter().zip(field).flat_map(|(state, axes)| {
        axes.iter().map(move |&k| (state.k - k).abs() / state.k.max(k))
    }).fold(0.0, f64::max)
}

pub(super) fn invalid(field: &'static str, reason: impl core::fmt::Display) -> ChtError {
    ChtError::InvalidInput { field, reason: reason.to_string() }
}

pub(super) fn poll(gate: &CancelGate) -> Result<(), ChtError> {
    if gate.is_requested() { Err(ChtError::Cancelled) } else { Ok(()) }
}
