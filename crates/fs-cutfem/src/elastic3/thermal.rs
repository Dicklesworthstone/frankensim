//! Source-owned isothermal elasticity on the retained reference geometry.
//!
//! Each active cell carries one original thermoelastic response: current C(T),
//! integrated free strain, thermal derivatives, domain and reference inertia.
//! Temperature is constant within a cell, as is its material tensor. The native
//! potential owns the thermodynamics; this module integrates that response and
//! uses the existing eigenstrain weak form, B-bar map and support terms.

use std::ops::ControlFlow;

use fs_material::elastic::thermal::IsothermalThermoelasticState;
use fs_material::tensor::from_mandel;
use fs_solver::op::LinearOp;

use super::eigenstrain::{EigenstrainResponse3, PreparedEigenstrain3, physical_response};
use super::{Cell3, CutElasticity3, ElasticityError3};

mod adjoint;
mod entropy;
pub(in crate::elastic3) mod equilibrium;
mod exchange;
mod jacobian;
mod linear_response;
mod sensitivity;
mod stress;
pub(in crate::elastic3) mod state;
pub use state::{ThermoelasticState3, ThermoelasticStateSettings3};

pub use exchange::{ThermoelasticCellExchange3, ThermoelasticExchangeLimits3};
pub use equilibrium::{ThermoelasticEquilibrium3, ThermoelasticEquilibriumSettings3};
pub use entropy::{
    PreparedThermoelasticEntropy3, ThermoelasticCellPotential3, ThermoelasticEntropyLimits3,
};
pub use jacobian::{PreparedThermalResidualJacobian3, ThermalJacobianLimits3};
pub use linear_response::{
    PreparedThermoelasticLinearResponse3, ThermoelasticLinearResponseLimits3,
    ThermoelasticLinearSolveSettings3, ThermoelasticTemperaturePullback3,
    ThermoelasticTemperatureResponse3,
};
pub(in crate::elastic3) use jacobian::ThermalCoordinateMap3;
pub use stress::ThermoelasticStressPullback3;
pub use adjoint::{
    ThermoelasticAdjointSettings3, ThermoelasticDesignObjective3,
    ThermoelasticDesignSensitivity3, ThermoelasticMotion3, ThermoelasticObjective3,
    ThermoelasticSensitivity3,
};

/// Mechanical response and elastic contribution to the native thermal potential.
///
/// Extensive quantities include this point's reference-volume weight and the
/// current cell stiffness scale. They exclude caloric storage, heat transport,
/// ghost stabilization and support terms. In particular, heat_capacity_j_k
/// is the elastic contribution at fixed total strain and may be negative; the
/// coupled thermal owner must combine it with its original caloric chart.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermoelasticResponse3 {
    pub mechanical: EigenstrainResponse3,
    pub temperature_k: f64,
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    pub heat_capacity_j_k: f64,
    /// Local partial derivative at fixed constitutive TOTAL strain, in physical
    /// tensor shear [xx, yy, zz, xy, yz, zx], including the stiffness scale.
    /// This is not the derivative through the solved equilibrium.
    pub stress_temperature_derivative: [f64; 6],
}

fn checkpoint(poll: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if poll().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

impl CutElasticity3 {
    /// Original fixed-temperature potential retained by this active cell.
    /// Cold/legacy material fields and out-of-range indices return None.
    #[must_use]
    pub fn cell_thermoelastic_state(&self, cell: usize) -> Option<&IsothermalThermoelasticState> {
        self.cells.get(cell)?.material.thermal_state()
    }

    /// Prepare the native integrated expansion with the SAME current C(T) as K.
    ///
    /// No alpha-times-temperature approximation, adiabatic correction or density
    /// update occurs. Every cell must retain a complete thermal source. The
    /// returned preparation works with existing equilibrium, stress, dynamics,
    /// design and prestress APIs. After refresh_material_field succeeds it
    /// must be prepared again, since both C(T) and free strain may have changed.
    /// The point allowance includes bulk and selected embedded support samples.
    pub fn prepare_thermoelastic(
        &self,
        max_points: usize,
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedEigenstrain3, ElasticityError3> {
        checkpoint(&mut poll)?;
        for cell in &self.cells {
            checkpoint(&mut poll)?;
            if cell.material.thermal_state().is_none() {
                return Err(ElasticityError3::Invalid(
                    "thermoelastic preparation requires a native thermal source in every cell",
                ));
            }
        }
        self.prepare_eigenstrain(
            &|id, _| {
                let source = self.cells[id]
                    .material
                    .thermal_state()
                    .expect("thermal sources checked before preparation")
                    .source();
                Ok(from_mandel(&source.free_strain_mandel()))
            },
            max_points,
            poll,
        )
    }

    /// Recover the native thermal potential and physical mechanical response.
    ///
    /// This uses the source's own integrated free strain, with no user-supplied
    /// preparation that could refer to a different thermal experiment. Actual
    /// and B-bar total/elastic strains satisfy the source's original component
    /// box as well as the caller's total-strain Frobenius budget in (0,1).
    /// Potentials are evaluated at the constitutive total strain. Only retained
    /// bulk points are observed; no between-point extrema or whole-system
    /// thermodynamic balance is asserted. max_points bounds the output first.
    pub fn thermoelastic_response(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<ThermoelasticResponse3>, ElasticityError3> {
        self.thermoelastic_response_mapped::<true>(
            u,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            |cell, field| {
                std::array::from_fn(|i| {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        0.0
                    } else {
                        field[3 * node + i % 3]
                    }
                })
            },
            poll,
        )
    }

    /// Admit the same native thermal/mechanical state without allocating output.
    pub fn admit_thermoelastic(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.thermoelastic_response_mapped::<false>(
            u,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            |cell, field| {
                std::array::from_fn(|i| {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        0.0
                    } else {
                        field[3 * node + i % 3]
                    }
                })
            },
            poll,
        )
        .map(|_| ())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn thermoelastic_response_mapped<const COLLECT: bool>(
        &self,
        u: &[f64],
        n: usize,
        maximum_kinematic_strain: f64,
        max_points: usize,
        mut local: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<ThermoelasticResponse3>, ElasticityError3> {
        checkpoint(&mut poll)?;
        if u.len() != n
            || !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(ElasticityError3::Invalid(
                "invalid thermoelastic displacement or kinematic bound",
            ));
        }
        for chunk in u.chunks(256) {
            checkpoint(&mut poll)?;
            if chunk.iter().any(|v| !v.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "nonfinite thermoelastic displacement",
                ));
            }
        }
        let mut points = 0usize;
        for cell in &self.cells {
            checkpoint(&mut poll)?;
            if cell.material.thermal_state().is_none() {
                return Err(ElasticityError3::Invalid(
                    "thermoelastic response requires a native thermal source in every cell",
                ));
            }
            points =
                points
                    .checked_add(cell.rules.bulk().len())
                    .ok_or(ElasticityError3::Invalid(
                        "thermoelastic point count overflow",
                    ))?;
            if points > max_points {
                return Err(ElasticityError3::Invalid(
                    "thermoelastic point allowance exhausted",
                ));
            }
        }
        let mut response = if COLLECT {
            Vec::with_capacity(points)
        } else {
            Vec::new()
        };
        for (id, cell) in self.cells.iter().enumerate() {
            checkpoint(&mut poll)?;
            let source = cell
                .material
                .thermal_state()
                .expect("preflight checked source")
                .source();
            let free_mandel = source.free_strain_mandel();
            let free_strain = from_mandel(&free_mandel);
            let mean_trace = free_mandel[..3].iter().sum::<f64>();
            let local = local(cell, u);
            let scale = self.scales[id];
            for &(position, weight) in cell.rules.bulk() {
                checkpoint(&mut poll)?;
                let mechanical = physical_response(
                    cell,
                    &local,
                    position,
                    free_strain,
                    mean_trace,
                    scale,
                    maximum_kinematic_strain,
                )?;
                let potential = source
                    .state(mechanical.constitutive_total_mandel, weight)
                    .map_err(|_| {
                        ElasticityError3::Invalid(
                            "native thermoelastic potential is outside its admitted domain",
                        )
                    })?;
                let free_energy_j = scale * potential.free_energy_j;
                let entropy_j_k = scale * potential.entropy_j_k;
                let internal_energy_j = scale * potential.internal_energy_j;
                let heat_capacity_j_k = scale * potential.heat_capacity_j_k;
                let stress_temperature_derivative =
                    from_mandel(&potential.stress_temperature_derivative_mandel_pa_k)
                        .map(|value| scale * value);
                if [
                    free_energy_j,
                    entropy_j_k,
                    internal_energy_j,
                    heat_capacity_j_k,
                ]
                .iter()
                .chain(&stress_temperature_derivative)
                .any(|v| !v.is_finite())
                {
                    return Err(ElasticityError3::Invalid(
                        "scaled thermoelastic potential overflow",
                    ));
                }
                if COLLECT {
                    response.push(ThermoelasticResponse3 {
                        mechanical: EigenstrainResponse3 {
                            cell: id,
                            position,
                            weight,
                            strain: mechanical.strain,
                            free_strain,
                            elastic_strain: mechanical.elastic_strain,
                            constitutive_elastic_strain: mechanical.constitutive_elastic_strain,
                            reference_stress: mechanical.reference_stress,
                            stress: mechanical.stress,
                            energy_density: mechanical.energy_density,
                        },
                        temperature_k: source.temperature_k(),
                        free_energy_j,
                        entropy_j_k,
                        internal_energy_j,
                        heat_capacity_j_k,
                        stress_temperature_derivative,
                    });
                }
            }
        }
        checkpoint(&mut poll)?;
        Ok(response)
    }
}
