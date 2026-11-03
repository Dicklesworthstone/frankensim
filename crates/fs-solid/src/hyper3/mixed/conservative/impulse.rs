//! Same-position impulses preserve the complete conservative history.

use super::{
    ConservativeFormulation, ConservativeTetError, ConservativeTetSolid, MixedUpError,
    SolidPointKinematics, SolidPointMobility, SolidPointMobilityMatrix, Tet10ImpulsePoint,
    checkpoint, gate, inertia_error, invalid,
};
use crate::hyper3::mixed::impulse::{
    ConstrainedImpulseResponse, MixedImpulseControls, MixedImpulseReport,
};
use crate::inertia::point::tet10;
use crate::inertia::{SolidMomentum, admit_reference_impulses};
use fs_exec::Cx;
use std::sync::Arc;

#[derive(Debug)]
pub struct ConservativeTetImpulse<'a> {
    pub solid: ConservativeTetSolid<'a>,
    pub impulse: MixedImpulseReport,
}

/// Contact response at one accepted geometry, retaining the original inertia,
/// supports and constitutive formulation. Exact incompressibility prepares one
/// native KKT factor for every mobility column and the resulting event.
/// Other formulations retain their original supported consistent-mass solves.
#[derive(Debug)]
pub struct ConservativePointImpulseResponse<'s, 'a> {
    solid: &'s ConservativeTetSolid<'a>,
    held_velocity_dofs: &'s [usize],
    controls: MixedImpulseControls,
    constrained: Option<ConstrainedImpulseResponse<'s>>,
}

impl<'a> ConservativeTetSolid<'a> {
    /// Current isoparametric position and work-conjugate quadratic velocity.
    /// Every signed P2 shape weight contributes; no point mass is substituted.
    pub fn point_kinematics(
        &self,
        point: Tet10ImpulsePoint,
        cx: &Cx<'_>,
    ) -> Result<SolidPointKinematics, ConservativeTetError> {
        tet10::kinematics(
            point,
            self.definition.elements,
            self.definition.nodes_m,
            &self.displacement_m,
            &self.velocity_mps,
            cx,
        )
        .map_err(inertia_error)
    }

    /// Prepare supported point response using this solid's actual formulation.
    /// The borrow keeps geometry, material history and supports fixed while
    /// contact mobility and its accepted event share the original operator.
    ///
    /// FullMaterial and OriginalMaterial require no pressure gauge and use
    /// only `controls.balance`. ExactIncompressible retains the native pressure
    /// factor and checks the incoming weak volume rate without projection.
    pub fn point_impulse_response<'s>(
        &'s self,
        held_velocity_dofs: &'s [usize],
        pressure_gauge_node: Option<usize>,
        controls: MixedImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<ConservativePointImpulseResponse<'s, 'a>, ConservativeTetError> {
        checkpoint(cx)?;
        let constrained = if let ConservativeFormulation::ExactIncompressible {
            maximum_point_volume_work_defect,
            ..
        } = self.definition.formulation
        {
            let response = ConstrainedImpulseResponse::prepare(
                &self.mesh,
                &self.mass,
                held_velocity_dofs,
                pressure_gauge_node,
                controls,
                |e, q| {
                    let point = &self.points[e * self.mesh.points_per_element + q];
                    point
                        .volume_state(
                            point.accepted_deformation(),
                            maximum_point_volume_work_defect,
                        )
                        .map(|volume| volume.endpoint_cofactor)
                        .map_err(|error| MixedUpError::Material {
                            element: e,
                            point: q,
                            error,
                        })
                },
                cx,
            )?;
            response.validate_velocity(&self.velocity_mps, cx)?;
            Some(response)
        } else {
            if pressure_gauge_node.is_some() {
                return Err(invalid("a compressible impulse has no pressure gauge"));
            }
            admit_reference_impulses(
                self.mass.ndofs(),
                held_velocity_dofs,
                &[],
                controls.balance,
                cx,
            )
            .map_err(inertia_error)?;
            None
        };
        Ok(ConservativePointImpulseResponse {
            solid: self,
            held_velocity_dofs,
            controls,
            constrained,
        })
    }

    /// Apply an instantaneous load through original consistent mass, including
    /// exact weak volume constraints when required. Geometry, material points,
    /// reference density, stored pressure and time are retained unchanged.
    ///
    /// FullMaterial and OriginalMaterial use the original checked mass solve;
    /// only `controls.balance` applies, no pressure impulse is introduced, and
    /// pressure/rate arrays and `kkt_inertia` are empty. A pressure gauge
    /// must be None for these compressible formulations. Physical event work
    /// enters the same cumulative account as subsequent native time intervals.
    ///
    /// Held components keep their current velocity. In ExactIncompressible,
    /// a pressure gauge is zero impulse on one original corner and is admitted
    /// only for a true connected
    /// constant null mode. Other rank deficiencies refuse in the native solver.
    ///
    /// ExactIncompressible pre-event velocity must satisfy the stated weak
    /// volume-rate allowance.
    /// The time integrator enforces endpoint volume, so its displacement-chord
    /// velocity is not automatically an admissible tangent velocity. This event
    /// does not conceal a large pre-existing defect by a dissipative projection.
    ///
    /// # Errors
    /// Invalid input, incompatible velocity/gauge,
    /// insufficient sparse/workspace budgets, native factorization or physical
    /// balance failure, or cancellation. No accepted state is mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        pressure_gauge_node: Option<usize>,
        origin_m: [f64; 3],
        controls: MixedImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetImpulse<'a>, ConservativeTetError> {
        self.point_impulse_response(held_velocity_dofs, pressure_gauge_node, controls, cx)?
            .nodal_impulse_candidate(nodal_impulses_kg_m_s, origin_m, cx)
    }

    /// One publication owner for nodal and material-point events, including
    /// the original compensated accumulated work and retained bulk history.
    fn finish_impulse_candidate(
        &self,
        velocity_mps: Vec<[f64; 3]>,
        impulse: MixedImpulseReport,
        work_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetImpulse<'a>, ConservativeTetError> {
        // Only external/support work enters the conservative trajectory ledger.
        // Pressure work is independently zero-gated; it is not supplied energy.
        let work = impulse.external_work_j + impulse.momentum.support_work_j;
        let increment = work - self.work_correction_j;
        let cumulative_work_j = self.cumulative_work_j + increment;
        let work_correction_j = (cumulative_work_j - self.cumulative_work_j) - increment;
        gate(
            "impulsive accumulated total work",
            cumulative_work_j
                - (self.stored_energy_j + impulse.momentum.kinetic_after_j
                    - self.initial_total_energy_j),
            work_tolerance_j,
        )?;
        let solid = Self {
            definition: self.definition,
            mesh: Arc::clone(&self.mesh),
            mass: Arc::clone(&self.mass),
            corner_reference_density: self.corner_reference_density.clone(),
            points: self.points.clone(),
            bulk_points: self.bulk_points.clone(),
            bulk_volume_ratios: self.bulk_volume_ratios.clone(),
            displacement_m: self.displacement_m.clone(),
            velocity_mps,
            algorithmic_pressure_pa: self.algorithmic_pressure_pa.clone(),
            pressure_residual_m3: self.pressure_residual_m3.clone(),
            stored_energy_j: self.stored_energy_j,
            bulk_stored_energy_j: self.bulk_stored_energy_j,
            kinetic_energy_j: impulse.momentum.kinetic_after_j,
            initial_total_energy_j: self.initial_total_energy_j,
            cumulative_work_j,
            work_correction_j,
            time_s: self.time_s,
        };
        checkpoint(cx)?;
        Ok(ConservativeTetImpulse { solid, impulse })
    }
}

impl<'a> ConservativePointImpulseResponse<'_, 'a> {
    /// Native pressure factor inertia, absent for finite-compliance mass solves.
    pub fn kkt_inertia(&self) -> Option<fs_sparse::Inertia> {
        self.constrained
            .as_ref()
            .map(|response| response.kkt_inertia())
    }

    /// Complete cross-point/direction mobility through the original supported
    /// mass and, in exact mode, its volume constraint. Unit world directions
    /// are required. Columns retain their actual checked solves; reciprocity
    /// is never imposed by symmetrization or clipping.
    pub fn point_mobility_matrix_per_kg(
        &self,
        rows: &[(Tet10ImpulsePoint, [f64; 3])],
        maximum_rows: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobilityMatrix, ConservativeTetError> {
        tet10::mobility_matrix_with(
            self.solid.velocity_mps.len(),
            self.solid.definition.elements,
            rows,
            maximum_rows,
            self.controls.balance,
            |loads| {
                if let Some(response) = &self.constrained {
                    response
                        .free_velocity_response(loads, cx)
                        .map_err(Into::into)
                } else {
                    self.solid
                        .mass
                        .free_velocity_response(
                            self.held_velocity_dofs,
                            loads,
                            self.controls.balance,
                            cx,
                        )
                        .map_err(inertia_error)
                }
            },
            inertia_error,
            cx,
        )
    }

    pub fn point_directional_inverse_mass_per_kg(
        &self,
        point: Tet10ImpulsePoint,
        direction: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobility, ConservativeTetError> {
        let mut matrix = self.point_mobility_matrix_per_kg(&[(point, direction)], 1, cx)?;
        let inverse_mass_per_kg = matrix.matrix_per_kg[0];
        if inverse_mass_per_kg < 0.0 {
            return Err(invalid("negative conservative material-point mobility"));
        }
        Ok(SolidPointMobility {
            inverse_mass_per_kg,
            solve: matrix.solves.pop().expect("one admitted mobility row"),
        })
    }

    /// Simultaneous point loads use the transpose of their original signed
    /// quadratic traces. Their nodal contributions accumulate with compensation
    /// before one event solve and one original cumulative-work/history transfer.
    pub fn point_impulses_candidate(
        &self,
        impulses: &[(Tet10ImpulsePoint, [f64; 3])],
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetImpulse<'a>, ConservativeTetError> {
        let loads = tet10::distributions(
            self.solid.definition.elements,
            self.solid.velocity_mps.len(),
            impulses,
            self.controls.balance,
            cx,
        )
        .map_err(inertia_error)?;
        self.nodal_impulse_candidate(&loads, origin_m, cx)
    }

    fn nodal_impulse_candidate(
        &self,
        nodal_impulses_kg_m_s: &[(usize, f64)],
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetImpulse<'a>, ConservativeTetError> {
        checkpoint(cx)?;
        let solid = self.solid;
        let mut positions = Vec::with_capacity(solid.velocity_mps.len());
        for (node, (x, u)) in solid
            .definition
            .nodes_m
            .iter()
            .zip(&solid.displacement_m)
            .enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            positions.push(core::array::from_fn(|c| x[c] + u[c]));
        }
        let (velocity_mps, impulse) = if let Some(response) = &self.constrained {
            response.impulse_candidate(
                &positions,
                &solid.velocity_mps,
                nodal_impulses_kg_m_s,
                origin_m,
                cx,
            )?
        } else {
            let (velocity, momentum) = solid
                .mass
                .impulse_candidate(
                    &positions,
                    &solid.velocity_mps,
                    self.held_velocity_dofs,
                    nodal_impulses_kg_m_s,
                    origin_m,
                    self.controls.balance,
                    cx,
                )
                .map_err(inertia_error)?;
            let report = MixedImpulseReport {
                external: momentum.applied,
                external_work_j: momentum.applied_work_j,
                pressure: SolidMomentum {
                    origin_m,
                    linear_kg_m_s: [0.0; 3],
                    angular_kg_m2_s: [0.0; 3],
                },
                momentum,
                pressure_work_j: 0.0,
                pressure_nodes: Vec::new(),
                pressure_impulses_pa_s: Vec::new(),
                volume_rate_before_m3_s: Vec::new(),
                volume_rate_after_m3_s: Vec::new(),
                kkt_inertia: None,
            };
            (velocity, report)
        };
        solid.finish_impulse_candidate(
            velocity_mps,
            impulse,
            self.controls.balance.work_tolerance_j,
            cx,
        )
    }
}
