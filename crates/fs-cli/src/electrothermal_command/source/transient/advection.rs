//! Transport and storage use the same admitted source chart and enthalpy datum.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    AdvectiveEnthalpyField, AffineSpecificEnthalpy, ChtError, EnergySolution, ThermalFace,
    ThermalSetup,
};
use fs_material::phase::ResolvedHeatCapacityCurve;

use super::{
    super::{Problem, invalid, poll},
    Storage,
};

pub(crate) struct Advection {
    cells: Vec<Option<AffineSpecificEnthalpy>>,
    boundary: [Option<f64>; 6],
}

impl Advection {
    pub(crate) fn field(&self) -> AdvectiveEnthalpyField<'_> {
        AdvectiveEnthalpyField {
            cells: &self.cells,
            boundary_j_kg: self.boundary,
        }
    }

    /// The steady and transient owners share one original caloric field.
    /// Solids receive no fluid Cp, and reservoirs use their exact enthalpy.
    pub(crate) fn resolve(
        curve: &ResolvedHeatCapacityCurve,
        problem: &Problem,
        temperature: &[f64],
        setup: &ThermalSetup,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        if temperature.len() != problem.domain.cell_count() {
            return Err(invalid("material.advection", "temperature shape mismatch"));
        }
        let datum_k = curve.curve().knots()[0].temperature_k;
        let mut cells = Vec::with_capacity(temperature.len());
        for (cell, &t) in temperature.iter().enumerate() {
            if cell.is_multiple_of(256) {
                poll(gate)?;
            }
            cells.push(if problem.domain.is_fluid(cell) {
                Some(AffineSpecificEnthalpy {
                    temperature_k: t,
                    enthalpy_j_kg: curve
                        .specific_enthalpy_change(datum_k, t)
                        .map_err(|e| invalid("material.advection", e))?,
                    heat_capacity_j_kg_k: curve
                        .secant_heat_capacity(t, t)
                        .map_err(|e| invalid("material.advection", e))?,
                })
            } else {
                None
            });
        }
        let mut boundary = [None; 6];
        for (face, rule) in setup.faces.iter().enumerate() {
            let t = match rule {
                ThermalFace::Inflow { temperature } => *temperature,
                ThermalFace::Outflow {
                    backflow_temperature,
                } => *backflow_temperature,
                _ => continue,
            };
            boundary[face] = Some(
                curve
                    .specific_enthalpy_change(datum_k, t)
                    .map_err(|e| invalid("material.advection", e))?,
            );
        }
        poll(gate)?;
        Ok(Self { cells, boundary })
    }
}

impl Storage<'_> {
    pub(crate) fn advection(
        &self,
        problem: &Problem,
        temperature: &[f64],
        setup: &ThermalSetup,
        gate: &CancelGate,
    ) -> Result<Option<Advection>, ChtError> {
        problem
            .fluid
            .as_ref()
            .map(|fluid| {
                Advection::resolve(
                    &self.curves[fluid.material_index],
                    problem,
                    temperature,
                    setup,
                    gate,
                )
            })
            .transpose()
    }

    /// Sum absolute *face-local* linearization errors, in joules, on both
    /// incident control volumes. Opposite caloric errors cannot hide each other.
    pub(crate) fn advection_consistent(
        &self,
        problem: &Problem,
        temperature: &[f64],
        frozen: &Advection,
        gate: &CancelGate,
    ) -> Result<bool, ChtError> {
        let fluid = problem.fluid.as_ref().expect("admitted advective storage");
        let curve = &self.curves[fluid.material_index];
        let datum_k = curve.curve().knots()[0].temperature_k;
        let mut actual = vec![0.0; temperature.len()];
        let mut linear = vec![0.0; temperature.len()];
        for (cell, &t) in temperature.iter().enumerate() {
            if cell.is_multiple_of(256) {
                poll(gate)?;
            }
            if let Some(h) = frozen.cells[cell] {
                actual[cell] = curve
                    .specific_enthalpy_change(datum_k, t)
                    .map_err(|e| invalid("material.advection", e))?;
                linear[cell] = h.at_temperature(t);
            }
        }
        let (mut defect, mut scale) = (0.0, 0.0);
        for cell in 0..temperature.len() {
            if cell.is_multiple_of(256) {
                poll(gate)?;
            }
            if !problem.domain.is_fluid(cell) {
                continue;
            }
            for face in 0..6 {
                let mass = fluid.properties.density_kg_m3
                    * fluid.flow.outward(&problem.domain, cell, face);
                if mass == 0.0 {
                    continue;
                }
                let donor = if mass > 0.0 {
                    Some(cell)
                } else {
                    problem.domain.neighbor(cell, face)
                };
                if let Some(donor) = donor {
                    let weight = self.dt * mass.abs();
                    defect += weight * (actual[donor] - linear[donor]).abs();
                    scale += weight * (actual[donor].abs() + linear[donor].abs());
                }
                // Incoming reservoir enthalpy is already evaluated exactly.
            }
        }
        let limit = self.absolute_j + self.relative * scale;
        if !defect.is_finite() || !limit.is_finite() {
            return Err(invalid(
                "material.advection",
                "local advected enthalpy defect is unrepresentable",
            ));
        }
        Ok(defect <= limit)
    }

    /// Replace only the frozen boundary advective contribution with true chart
    /// enthalpy. Native conduction and physical sinks retain their own accounts.
    pub(super) fn physical_boundary_outflow_w(
        &self,
        problem: &Problem,
        setup: &ThermalSetup,
        solution: &EnergySolution,
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        let Some(fluid) = &problem.fluid else {
            return Ok(solution.report.balance.boundary_outflow_w);
        };
        let actual = self
            .advection(problem, &solution.temperature, setup, gate)?
            .expect("fluid caloric field");
        let mut advected = 0.0;
        for cell in 0..problem.domain.cell_count() {
            if cell.is_multiple_of(256) {
                poll(gate)?;
            }
            if !problem.domain.is_fluid(cell) {
                continue;
            }
            for face in 0..6 {
                if problem.domain.neighbor(cell, face).is_some() {
                    continue;
                }
                let mass = fluid.properties.density_kg_m3
                    * fluid.flow.outward(&problem.domain, cell, face);
                if mass == 0.0 {
                    continue;
                }
                let h = if mass > 0.0 {
                    actual.cells[cell].expect("fluid cell").enthalpy_j_kg
                } else {
                    actual.boundary[face].expect("open boundary")
                };
                advected += mass * h;
            }
        }
        let balance = &solution.report.balance;
        let outflow = balance.boundary_outflow_w + (advected - balance.advective_outflow_w);
        if !outflow.is_finite() {
            return Err(invalid(
                "material.advection",
                "physical boundary outflow overflow",
            ));
        }
        Ok(outflow)
    }
}
