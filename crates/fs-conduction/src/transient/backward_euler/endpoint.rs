//! Evaluate a proposed physical endpoint without advancing or solving it.
//!
//! Partitioned temperature-dependent boundary laws must pass the complete
//! free nodal balance as well as total energy closure. Equal and opposite
//! nodal defects can disappear from a whole-body heat balance.

use super::radiation::boundary_at_temperature;
use super::{BackwardEuler, StepConfig, energy_balance, finite, invalid, poll, sum};
use crate::{
    AmbientRadiationPatch, ConductionError, ConductionProblem, RobinFlux,
    ThermalBoundary, ThermalInterfaces,
};
use fs_exec::Cx;

/// Physical residual and boundary exchanges at one supplied endpoint.
///
/// These values do not establish convergence: the caller owns its fixed
/// residual and energy gates. No temperature, matrix or solver state is
/// returned, and neither the candidate nor its physical history is changed.
#[derive(Debug, Clone)]
pub struct EndpointEvaluation {
    /// Euclidean norm of C(T - T_old) + dt (K(T)T - b), on free nodes, J.
    pub residual_norm_j: f64,
    /// Number of free temperature degrees of freedom in that norm.
    pub free_dofs: usize,
    /// Full 1^T C(T - T_old), including prescribed nodes, J.
    pub stored_energy_change_j: f64,
    /// Integrated endpoint volumetric generation, W.
    pub source_w: f64,
    /// Outward endpoint Neumann heat, W.
    pub neumann_out_w: f64,
    /// Physical prescribed-temperature reaction, including fixed-node storage, W inward.
    pub dirichlet_in_w: f64,
    /// Original Robin convection only, W outward.
    pub convective_out_w: f64,
    /// Nonlinear ambient gray-patch heat, W outward; negative heats the solid.
    pub radiation_out_w: f64,
    /// Storage minus dt times the complete physical net input, J.
    pub energy_residual_j: f64,
    /// Original convective traces, with radiation kept out of air exchange.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Boundary re-evaluated at the candidate's radiative patch temperatures.
    pub combined_boundary: ThermalBoundary,
}

impl BackwardEuler<'_> {
    /// Reassemble and evaluate a supplied endpoint without a linear or nonlinear solve.
    ///
    /// The original boundary supplies the caller's current convection laws.
    /// When patches are present, the shared area-mean gray-radiation law
    /// supplies their physical coefficients at this candidate. Conductivity
    /// and matching contact use the same endpoint assembly as nonlinear
    /// backward Euler. The full old field remains the capacity history.
    ///
    /// Candidate prescribed temperatures must already equal the supplied
    /// endpoint boundary. The reaction includes their change from history,
    /// including a first-step prescribed-temperature jump exactly once.
    ///
    /// A finite off-equilibrium candidate is returned with its actual residual,
    /// even when the residual or energy mismatch exceeds the supplied budget.
    /// This permits an outer coupling owner to apply a fixed physical gate.
    /// StepConfig supplies the same input admission as a prescribed step; this
    /// evaluation consumes no Krylov or Newton iterations.
    ///
    /// # Errors
    ///
    /// Refuses an invalid step, mesh/field mismatch, nonfinite data or output,
    /// a candidate that violates prescribed values, invalid material/contact
    /// data, invalid or empty supplied radiation patches, or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_prescribed_endpoint(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old: &[f64],
        candidate: &[f64],
        dt_s: f64,
        step_config: StepConfig,
        patches: Option<&[AmbientRadiationPatch]>,
    ) -> Result<EndpointEvaluation, ConductionError> {
        poll(cx, 0)?;
        let n = self.mesh.vertex_count();
        // Boundary objects may be reused by callers; validate their nodal
        // fields before any assembly or prescribed-index access.
        for (index, condition) in problem.boundary.conditions().iter().enumerate() {
            if index % 512 == 0 {
                poll(cx, index)?;
            }
            condition.validate(n)?;
        }
        for &(vertex, temperature) in problem.boundary.dirichlet() {
            if vertex >= n {
                return Err(invalid("prescribed endpoint vertex is outside the capacity mesh"));
            }
            finite(temperature)?;
        }
        let dofs = self.admit_endpoint_step(cx, problem, old, dt_s, step_config)?;
        if candidate.len() != n {
            return Err(ConductionError::FieldLength {
                field: "endpoint temperature",
                expected: n,
                found: candidate.len(),
            });
        }
        for (index, &temperature) in candidate.iter().enumerate() {
            if index % 512 == 0 {
                poll(cx, index)?;
            }
            finite(temperature)?;
        }
        for (index, &vertex) in dofs.fixed().iter().enumerate() {
            if index % 512 == 0 {
                poll(cx, index)?;
            }
            if candidate[vertex] != dofs.prescribed()[vertex] {
                return Err(invalid(
                    "evaluated candidate must match prescribed endpoint temperatures",
                ));
            }
        }
        let combined_boundary = match patches {
            Some(patches) => boundary_at_temperature(cx, problem, patches, candidate)?,
            None => problem.boundary.clone(),
        };
        let physical_problem = ConductionProblem {
            boundary: &combined_boundary,
            ..problem
        };
        let physical = self.evaluate_endpoint(
            cx, physical_problem, interfaces, old, candidate, dt_s, &dofs,
        )?;
        // The physical system supplies the reaction. Original boundary
        // integrals supply convection, without folding radiation into it.
        let (energy, convective_robin_fluxes) = energy_balance(
            self.mesh, problem.boundary, problem.source, &physical.system, &dofs, candidate,
        );
        poll(cx, n)?;
        let source_w = finite(energy.source_w)?;
        let neumann_out_w = finite(energy.neumann_out_w)?;
        let convective_out_w = finite(energy.robin_out_w)?;
        for (index, flux) in convective_robin_fluxes.iter().enumerate() {
            if index % 512 == 0 {
                poll(cx, index)?;
            }
            for value in [
                flux.area_m2,
                flux.mean_htc_w_per_m2_k,
                flux.mean_wall_temperature_k,
                flux.mean_reference_temperature_k,
                flux.heat_rate_w,
            ] {
                finite(value)?;
            }
        }
        let mut radiation_out_w = 0.0;
        if let Some(patches) = patches {
            for (index, patch) in patches.iter().enumerate() {
                poll(cx, index)?;
                let flux = convective_robin_fluxes.iter()
                    .find(|flux| flux.region == patch.region())
                    .ok_or_else(|| invalid("evaluated radiation patch has no convective trace"))?;
                let heat_w = finite(
                    patch.heat_flux_w_m2(flux.mean_wall_temperature_k)? * flux.area_m2,
                )?;
                radiation_out_w = finite(radiation_out_w + heat_w)?;
            }
        }
        let delta = candidate.iter().zip(old)
            .map(|(temperature, previous)| finite(temperature - previous))
            .collect::<Result<Vec<_>, _>>()?;
        let mut storage = vec![0.0; n];
        self.capacity.spmv(&delta, &mut storage);
        poll(cx, n)?;
        let stored_energy_change_j = sum(storage.iter().copied())?;
        let dirichlet_in_w = finite(
            finite(energy.dirichlet_in_w)?
                + sum(dofs.fixed().iter().map(|&vertex| storage[vertex]))? / dt_s,
        )?;
        let net_input_w = finite(
            source_w + dirichlet_in_w - neumann_out_w - convective_out_w - radiation_out_w,
        )?;
        let energy_residual_j = finite(stored_energy_change_j - dt_s * net_input_w)?;
        poll(cx, n)?;
        Ok(EndpointEvaluation {
            residual_norm_j: physical.norm,
            free_dofs: dofs.n(),
            stored_energy_change_j,
            source_w,
            neumann_out_w,
            dirichlet_in_w,
            convective_out_w,
            radiation_out_w,
            energy_residual_j,
            convective_robin_fluxes,
            combined_boundary,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ConductivityModel, ConductionMesh, LinearConfig, ScalarField, ThermalBc,
        ThermalBoundaryBuilder};
    use crate::transient::VolumetricHeatCapacity;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_rep_mesh::TetComplex;

    fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(
            gate, arena,
            StreamKey { seed: 83, kernel_id: 719, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        )))
    }

    fn mesh() -> ConductionMesh {
        ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ).unwrap()
    }

    fn config() -> StepConfig {
        StepConfig {
            linear: LinearConfig { tolerance: 1e-12, ..LinearConfig::default() },
            energy_tolerance_j: 1e-8,
        }
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!((actual - expected).abs() <= tolerance, "{actual:e} != {expected:e}");
    }

    #[test]
    fn equal_opposite_nodal_defects_do_not_disappear_with_closed_total_energy() {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(6.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [300.0; 4];
        let candidate = [300.0, 301.0, 299.0, 300.0];
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let engine = BackwardEuler::uniform(
                cx, &mesh, VolumetricHeatCapacity::declared(24.0).unwrap(),
            ).unwrap();
            let problem = ConductionProblem {
                mesh: &mesh, boundary: &boundary, material: &material,
                source: &source, element_materials: None,
            };
            let evaluation = engine.evaluate_prescribed_endpoint(
                cx, problem, None, &old, &candidate, 1.0, config(), None,
            ).unwrap();
            assert_eq!(evaluation.free_dofs, 4);
            close(evaluation.stored_energy_change_j, 0.0, 1e-12);
            close(evaluation.energy_residual_j, 0.0, 1e-12);
            // Unit right tetra: C_ii=24/24=1. The opposite changes at
            // nodes 1 and 2 give K*T=[0,1,-1,0] for k=6, hence the
            // complete free residual is [0,2,-2,0] J at dt=1.
            close(evaluation.residual_norm_j, 8.0_f64.sqrt(), 1e-12);
            assert!(evaluation.convective_robin_fluxes.is_empty());
            close(evaluation.radiation_out_w, 0.0, 0.0);
            assert_eq!(old, [300.0; 4]);
            assert_eq!(candidate, [300.0, 301.0, 299.0, 300.0]);
        });
    }

    #[test]
    fn prescribed_endpoint_audit_retains_boundary_storage_and_rejects_bad_fields() {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("driven", |face| (face.centroid.iter().sum::<f64>() - 1.0).abs() < 1e-12,
                ThermalBc::dirichlet(330.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [300.0; 4];
        let candidate = [310.0, 330.0, 330.0, 330.0];
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let engine = BackwardEuler::uniform(
                cx, &mesh, VolumetricHeatCapacity::declared(12.0).unwrap(),
            ).unwrap();
            let problem = ConductionProblem {
                mesh: &mesh, boundary: &boundary, material: &material,
                source: &source, element_materials: None,
            };
            let evaluation = engine.evaluate_prescribed_endpoint(
                cx, problem, None, &old, &candidate, 0.25, config(), None,
            ).unwrap();
            assert_eq!(evaluation.free_dofs, 1);
            close(evaluation.residual_norm_j, 0.0, 1e-11);
            close(evaluation.stored_energy_change_j, 50.0, 1e-11);
            close(evaluation.dirichlet_in_w, 200.0, 1e-9);
            close(evaluation.energy_residual_j, 0.0, 1e-9);
            assert!(engine.evaluate_prescribed_endpoint(
                cx, problem, None, &old, &candidate[..3], 0.25, config(), None,
            ).is_err());
            for bad in [
                [310.0, 329.0, 330.0, 330.0],
                [f64::NAN, 330.0, 330.0, 330.0],
            ] {
                assert!(engine.evaluate_prescribed_endpoint(
                    cx, problem, None, &old, &bad, 0.25, config(), None,
                ).is_err());
            }
            gate.request();
            assert!(matches!(engine.evaluate_prescribed_endpoint(
                cx, problem, None, &old, &candidate, 0.25, config(), None,
            ), Err(ConductionError::Cancelled { .. })));
        });
        assert_eq!(old, [300.0; 4]);
    }
}
