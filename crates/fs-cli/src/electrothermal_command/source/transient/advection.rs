//! Transport and storage use the same admitted source chart and enthalpy datum.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    AdvectiveEnthalpyField, AffineSpecificEnthalpy, ChtError, ThermalFace,
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

}
