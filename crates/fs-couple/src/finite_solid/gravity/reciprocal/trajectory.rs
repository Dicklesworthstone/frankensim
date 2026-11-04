//! Symmetric gravity kick / native material-and-rigid motion / gravity kick.
//!
//! Both participants move before the second original field evaluation. Native
//! objective solid motion, consistent mass, rigid free rotation and endpoint
//! impulse owners remain authoritative. The symmetric composition has a
//! second-order splitting error in smooth, adequately resolved regimes; its
//! measured energy defect is not physical heating or a trajectory error bound.

use super::{
    QuadrupoleSourceExtent, ReciprocalGravityKick, ReciprocalGravityKickControls, components, dd,
    vector,
};
use crate::finite_solid::contact::step::free::{
    FreeFlightControls, FreeFlightError, FreeFlightReport,
};
use crate::finite_solid::{Cx, FiniteSolidBody, FiniteSolidError, checkpoint};
use fs_mbd::{Gravity, MassProperties, RigidBodyState, Vec3, Wrench};

#[derive(Clone, Copy, Debug)]
pub struct ReciprocalGravityStepControls {
    /// `duration_s` is the COMPLETE interval; each gravity kick uses half of it.
    pub gravity: ReciprocalGravityKickControls,
    /// Original native nonlinear/step limits and complete-system P/L tolerances.
    pub middle: FreeFlightControls,
    pub maximum_duration_s: f64,
    /// Bounds additional O(nodes) endpoint geometry and native candidate storage.
    pub maximum_nodes: usize,
    /// Strict clearance from every swept relative tetrahedral bounding box [m].
    pub minimum_separation_m: f64,
    /// Absolute gravity kick-work minus interaction-potential decrease [J].
    pub splitting_work_tolerance_j: f64,
    /// Bounds the total defect and the sum of absolute subflow defects [J].
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Debug)]
pub struct ReciprocalGravityStep {
    pub first_kick: ReciprocalGravityKick,
    pub middle: FreeFlightReport,
    pub last_kick: ReciprocalGravityKick,
    /// Conservative largest-axis-gap bound, minimized over original tetrahedra.
    pub path_separation_lower_bound_m: f64,
    pub kinetic_change_j: f64,
    pub interaction_potential_change_j: f64,
    pub gravitational_splitting_defect_j: f64,
    /// Actual kinetic + strain + pair-potential change minus external work.
    pub integration_defect_j: f64,
    pub absolute_integration_defect_j: f64,
    pub work_energy_residual_j: f64,
    /// Recovered from original measured residual identities, without another
    /// full mass action or subtraction of large endpoint momenta.
    pub linear_momentum_change_kg_m_s: Vec3,
    pub angular_momentum_change_kg_m2_s: Vec3,
    /// Sum of both gravity-kick and middle-flow stationary support impulses.
    pub support_impulse_kg_m_s: Vec3,
    pub support_angular_impulse_kg_m2_s: Vec3,
    /// External applied loads are retained separately in `middle`.
    pub linear_residual_kg_m_s: Vec3,
    pub angular_residual_kg_m2_s: Vec3,
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
}

#[derive(Debug)]
pub enum ReciprocalGravityStepError {
    Invalid(&'static str),
    Solid(FiniteSolidError),
    Middle(FreeFlightError),
    Separation { lower_bound_m: f64, minimum_m: f64 },
    Balance(Box<ReciprocalGravityStep>),
}
impl From<FiniteSolidError> for ReciprocalGravityStepError {
    fn from(error: FiniteSolidError) -> Self {
        Self::Solid(error)
    }
}
impl From<FreeFlightError> for ReciprocalGravityStepError {
    fn from(error: FreeFlightError) -> Self {
        Self::Middle(error)
    }
}
impl core::fmt::Display for ReciprocalGravityStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reciprocal solid/source trajectory: {self:?}")
    }
}
impl std::error::Error for ReciprocalGravityStepError {}

impl<'a> FiniteSolidBody<'a> {
    /// Advance both original participants over one shared physical interval.
    ///
    /// The native objective solid step and symmetric rigid free drift sit
    /// between reciprocal half-kicks. Material strain, held solid loads,
    /// optional uniform pressure and pre-existing external fields retain their
    /// original work/support accounts. Do not bind this same source again as
    /// stationary external gravity. Source spin evolves only under its native
    /// torque-free motion; the monopole supplies no central spin impulse.
    ///
    /// The native solid's acceleration/history describe the MATERIAL MIDDLE
    /// flow, not a sampled continuous acceleration including the split kicks.
    /// There is one accepted native timestep. Both final states escape together.
    ///
    /// Per-element relative bounding boxes enclose each original tetrahedron
    /// along linear interpolation of the two configurations and the source COM.
    /// Their clearance rejects source tunneling on that declared numerical
    /// path. It is conservative and is not CCD for the curved physical motion,
    /// finite-source contact, or a quadrature-accuracy guarantee near the source.
    #[allow(clippy::too_many_arguments)]
    pub fn reciprocal_gravity_step_candidate(
        &self,
        source: RigidBodyState,
        properties: MassProperties,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        uniform_pressure_pa: Option<f64>,
        controls: ReciprocalGravityStepControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, RigidBodyState, ReciprocalGravityStep), ReciprocalGravityStepError> {
        self.gravity_source_step_candidate(
            source,
            properties,
            None,
            nodal_forces_n,
            body_force_n_m3,
            uniform_pressure_pa,
            controls,
            cx,
        )
    }

    /// Advance an extended, rotating source and the original deformable solid.
    ///
    /// Each half-kick samples the source's actual COM pose and full inertia and
    /// exchanges both linear impulse and tidal COM spin. Native free rotation
    /// therefore changes the endpoint field. The original material middle flow,
    /// supports, work accounts and separate splitting/total-energy gates remain
    /// unchanged. The accepted swept relative tetrahedra must stay farther than
    /// both `minimum_separation_m` and the declared source radius/ratio; this
    /// conservative linear-path admission is not physical continuous CCD or a
    /// multipole/solid-quadrature accuracy certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn reciprocal_quadrupole_gravity_step_candidate(
        &self,
        source: RigidBodyState,
        properties: MassProperties,
        extent: QuadrupoleSourceExtent,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        uniform_pressure_pa: Option<f64>,
        controls: ReciprocalGravityStepControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, RigidBodyState, ReciprocalGravityStep), ReciprocalGravityStepError> {
        self.gravity_source_step_candidate(
            source,
            properties,
            Some(extent),
            nodal_forces_n,
            body_force_n_m3,
            uniform_pressure_pa,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn gravity_source_step_candidate(
        &self,
        source: RigidBodyState,
        properties: MassProperties,
        extent: Option<QuadrupoleSourceExtent>,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        uniform_pressure_pa: Option<f64>,
        controls: ReciprocalGravityStepControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, RigidBodyState, ReciprocalGravityStep), ReciprocalGravityStepError> {
        poll(cx)?;
        let h = controls.gravity.duration_s;
        if !h.is_finite()
            || h <= 0.0
            || h > controls.maximum_duration_s
            || !controls.maximum_duration_s.is_finite()
            || controls.maximum_duration_s <= 0.0
            || controls.maximum_nodes == 0
            || self.problem.nodes_m.len() > controls.maximum_nodes
            || !controls.minimum_separation_m.is_finite()
            || controls.minimum_separation_m < 0.0
            || [
                controls.splitting_work_tolerance_j,
                controls.energy_tolerance_j,
                controls.middle.linear_tolerance_kg_m_s,
                controls.middle.angular_tolerance_kg_m2_s,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(ReciprocalGravityStepError::Invalid(
                "finite bounded trajectory controls required",
            ));
        }
        let minimum_separation_m = match extent {
            Some(extent) => controls.minimum_separation_m.max(
                extent
                    .minimum_sample_radius_m()
                    .map_err(|error| FiniteSolidError::Source(error.to_string()))?
                    .next_up(),
            ),
            None => controls.minimum_separation_m,
        };
        let before_x = self.current_positions(cx)?;
        self.gravity_path_clearance(
            &before_x,
            &before_x,
            source,
            source,
            minimum_separation_m,
            cx,
        )?;
        let half = ReciprocalGravityKickControls {
            duration_s: 0.5 * h,
            ..controls.gravity
        };
        let (kicked, kicked_source, first_kick) =
            self.gravity_source_kick_candidate(source, properties, extent, half, cx)?;
        let middle_result = match uniform_pressure_pa {
            None => kicked.rigid_free_flight_candidate(
                &[kicked_source],
                &[properties],
                nodal_forces_n,
                body_force_n_m3,
                &[Wrench::ZERO],
                Gravity::ZERO,
                h,
                half.balance_origin_m,
                controls.middle,
                cx,
            ),
            Some(pressure) => kicked.rigid_free_flight_candidate_with_pressure(
                &[kicked_source],
                &[properties],
                nodal_forces_n,
                body_force_n_m3,
                pressure,
                &[Wrench::ZERO],
                Gravity::ZERO,
                h,
                half.balance_origin_m,
                controls.middle,
                cx,
            ),
        };
        let (moved, moved_sources, middle) = middle_result?;
        let moved_source = moved_sources[0];
        let after_x = moved.current_positions(cx)?;
        let path_separation_lower_bound_m = self.gravity_path_clearance(
            &before_x,
            &after_x,
            source,
            moved_source,
            minimum_separation_m,
            cx,
        )?;
        let (candidate, final_source, last_kick) =
            moved.gravity_source_kick_candidate(moved_source, properties, extent, half, cx)?;
        let potential = dd(last_kick.potential_energy_j) - dd(first_kick.potential_energy_j);
        let split = dd(first_kick.impulse_work_j) + dd(last_kick.impulse_work_j) + potential;
        let kinetic = dd(first_kick.kinetic_energy_change_j)
            + dd(middle.kinetic_change_j)
            + dd(last_kick.kinetic_energy_change_j);
        let defect = kinetic + dd(middle.strain_change_j) + potential - dd(middle.source_work_j);
        let work_residual = defect
            - (dd(middle.integration_defect_j)
                + split
                + dd(first_kick.work_residual_j)
                + dd(last_kick.work_residual_j));
        let absolute_defect = middle.absolute_integration_defect_j
            + split.to_f64().abs()
            + first_kick.work_residual_j.abs()
            + last_kick.work_residual_j.abs()
            + work_residual.to_f64().abs();
        let support = sum_vectors([
            first_kick.solid.support.linear_kg_m_s,
            components(middle.support_impulse_kg_m_s),
            last_kick.solid.support.linear_kg_m_s,
        ]);
        let support_angular = sum_vectors([
            first_kick.solid.support.angular_kg_m2_s,
            components(middle.support_angular_impulse_kg_m2_s),
            last_kick.solid.support.angular_kg_m2_s,
        ]);
        let linear = sum_vectors([
            first_kick.linear_residual_kg_m_s,
            components(middle.linear_momentum_residual_kg_m_s),
            last_kick.linear_residual_kg_m_s,
        ]);
        let angular = sum_vectors([
            first_kick.angular_residual_kg_m2_s,
            components(middle.angular_momentum_residual_kg_m2_s),
            last_kick.angular_residual_kg_m2_s,
        ]);
        let absolute_linear = norm(first_kick.linear_residual_kg_m_s)
            + middle.absolute_linear_residual_kg_m_s
            + norm(last_kick.linear_residual_kg_m_s);
        let absolute_angular = norm(first_kick.angular_residual_kg_m2_s)
            + middle.absolute_angular_residual_kg_m2_s
            + norm(last_kick.angular_residual_kg_m2_s);
        let report = ReciprocalGravityStep {
            path_separation_lower_bound_m,
            kinetic_change_j: kinetic.to_f64(),
            interaction_potential_change_j: potential.to_f64(),
            gravitational_splitting_defect_j: split.to_f64(),
            integration_defect_j: defect.to_f64(),
            absolute_integration_defect_j: absolute_defect,
            work_energy_residual_j: work_residual.to_f64(),
            linear_momentum_change_kg_m_s: sum_vectors([
                components(support),
                components(linear),
                components(middle.external_impulse_kg_m_s),
            ]),
            angular_momentum_change_kg_m2_s: sum_vectors([
                components(support_angular),
                components(angular),
                components(middle.external_angular_impulse_kg_m2_s),
            ]),
            support_impulse_kg_m_s: support,
            support_angular_impulse_kg_m2_s: support_angular,
            linear_residual_kg_m_s: linear,
            angular_residual_kg_m2_s: angular,
            absolute_linear_residual_kg_m_s: absolute_linear,
            absolute_angular_residual_kg_m2_s: absolute_angular,
            first_kick,
            middle,
            last_kick,
        };
        if [
            report.kinetic_change_j,
            report.interaction_potential_change_j,
            report.integration_defect_j,
            absolute_defect,
            absolute_linear,
            absolute_angular,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || [
                report.linear_momentum_change_kg_m_s,
                report.angular_momentum_change_kg_m2_s,
                support,
                support_angular,
                linear,
                angular,
            ]
            .iter()
            .any(|v| !v.is_finite())
            || split.to_f64().abs() > controls.splitting_work_tolerance_j
            || defect.to_f64().abs() > controls.energy_tolerance_j
            || absolute_defect > controls.energy_tolerance_j
            || absolute_linear > controls.middle.linear_tolerance_kg_m_s
            || absolute_angular > controls.middle.angular_tolerance_kg_m2_s
        {
            return Err(ReciprocalGravityStepError::Balance(Box::new(report)));
        }
        poll(cx)?;
        Ok((candidate, final_source, report))
    }

    fn gravity_path_clearance(
        &self,
        before: &[[f64; 3]],
        after: &[[f64; 3]],
        source_before: RigidBodyState,
        source_after: RigidBodyState,
        minimum_m: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, ReciprocalGravityStepError> {
        let mut result = f64::INFINITY;
        for tetrahedron in self.problem.tetrahedra {
            poll(cx)?;
            let mut low = [f64::INFINITY; 3];
            let mut high = [f64::NEG_INFINITY; 3];
            for (positions, source) in [(before, source_before), (after, source_after)] {
                let center = components(source.pose().position_world());
                for &node in tetrahedron {
                    for c in 0..3 {
                        let relative = positions[node][c] - center[c];
                        if !relative.is_finite() {
                            return Err(ReciprocalGravityStepError::Invalid(
                                "unrepresentable relative gravity geometry",
                            ));
                        }
                        // Outward neighbors enclose subtraction roundoff. The
                        // largest separating coordinate is a Euclidean lower bound
                        // without a rounded square root or squared-distance overflow.
                        low[c] = low[c].min(relative.next_down());
                        high[c] = high[c].max(relative.next_up());
                    }
                }
            }
            let distance = (0..3).fold(0.0_f64, |gap, c| gap.max(low[c]).max(-high[c]));
            result = result.min(distance);
            if distance <= minimum_m {
                return Err(ReciprocalGravityStepError::Separation {
                    lower_bound_m: distance,
                    minimum_m,
                });
            }
        }
        Ok(result)
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), ReciprocalGravityStepError> {
    checkpoint(cx).map_err(|error| ReciprocalGravityStepError::Solid(error.into()))
}
fn norm(value: [f64; 3]) -> f64 {
    value[0].hypot(value[1]).hypot(value[2])
}
fn sum_vectors(values: [[f64; 3]; 3]) -> Vec3 {
    vector(std::array::from_fn(|c| {
        (dd(values[0][c]) + dd(values[1][c]) + dd(values[2][c])).to_f64()
    }))
}
