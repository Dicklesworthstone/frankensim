//! One original nonlinear pressure solve for both complete fields and mouths.
//! Only accepted controls warm-start the next original physical interval.

use super::*;
use crate::dg_port::DgIndependentMouthsResponse;
use crate::nsf_dg::NsfDgRadiationSolve;
use crate::nsf_dg::reservoir::{exchange, multiple::solve_joint_interval};

pub(super) struct JointPressureSolve {
    policy: NsfDgRadiationSolve,
    previous_coordinates: Vec<f64>,
}

impl<'a, 'c, L: Lattice3<27>> NsfThermoviscousDg<'a, 'c, L> {
    /// Shared explicit nonlinear policy for an upstream field and/or mouths.
    pub fn joint_solve(&self) -> Option<&NsfDgRadiationSolve> {
        self.joint.as_ref().map(|owner| &owner.policy)
    }

    pub(super) fn admit_joint_solve(&self, policy: NsfDgRadiationSolve) -> Result<(), NsfDgError> {
        policy.validate()?;
        if self.joint.as_ref().is_some_and(|owner| owner.policy != policy) {
            return Err(NsfDgError::Input(
                "both complete fields and their mouths require one identical explicit nonlinear policy"));
        }
        Ok(())
    }

    /// Called only after the attachment and shared policy were admitted.
    pub(super) fn bind_joint_solve(&mut self, policy: NsfDgRadiationSolve) {
        if self.joint.is_none() {
            self.joint = Some(JointPressureSolve { policy, previous_coordinates: Vec::new() });
        }
    }

    pub(super) fn step_joint(&mut self, drive: NsfDgDrive, external: NsfBoundaryDrive,
        rise: f64, breath: Option<(usize, NsfOpenState)>, gate: &CancelGate,
    ) -> Result<NsfThermoviscousDgFrame, NsfDgError> {
        let joint = self.joint.as_ref()
            .ok_or(NsfDgError::Input("complete fields lost their shared nonlinear policy"))?;
        let solve = joint.policy;
        let previous_coordinates = joint.previous_coordinates.clone();
        let (inputs, temperatures) = self.reservoir_radiations.as_ref().map_or_else(
            || (Vec::new(), Vec::new()),
            |owner| (owner.inputs.clone(), owner.temperatures_k.clone()),
        );
        let relation = if inputs.is_empty() {
            let scalar = self.load.port_response(gate).map_err(NsfDgError::Thermoviscous)?;
            DgIndependentMouthsResponse { history_pressure_pa: vec![scalar.history_pressure_pa],
                impedance_pa_s_m3: vec![scalar.impedance_pa_s_m3] }
        } else {
            self.load.independent_mouths_response(gate).map_err(NsfDgError::Thermoviscous)?
        };
        let upstream_relation = self.upstream.as_ref().map(|owner|
            owner.load.port_response(gate).map_err(NsfDgError::Thermoviscous)).transpose()?;
        let width = inputs.len() + 1;
        let limits = self.limits;
        let side = usize::from(self.high);
        let reference = self.reference;
        let step = self.steps + 1;
        let dt = self.load.duration_s();
        let load = &mut self.load;
        let upstream = &mut self.upstream;
        let finite = &mut self.finite_reservoir;
        let (frame, normalized) = solve_joint_interval(
            &mut self.zone, &inputs, &temperatures, external, breath, side,
            self.patch_background, reference.pressure, limits, solve,
            &previous_coordinates, &relation, upstream_relation.map(|relation| (relation, rise)), gate,
            |fluid, physical, evaluations| {
                let field = coupled_exchange(&fluid, side)?;
                let mouths = inputs.iter().map(|(aperture, _)| exchange(&fluid, *aperture))
                    .collect::<Result<Vec<_>, _>>()?;
                let mouth_flows: Vec<_> = mouths.iter()
                    .map(|mouth| mouth.normal_volume_flow_m3_s).collect();
                let mut flows = Vec::with_capacity(width);
                flows.push(field.normal_volume_flow_m3_s);
                flows.extend_from_slice(&mouth_flows);
                let scalar_relation = |row: usize| {
                    let mut history = relation.history_pressure_pa[row];
                    for (column, flow) in flows.iter().enumerate() {
                        if column != row {
                            history += relation.impedance_pa_s_m3[row * width + column] * flow;
                        }
                    }
                    DgPortResponse { history_pressure_pa: history,
                        impedance_pa_s_m3: relation.impedance_pa_s_m3[row * width + row] }
                };
                let acoustic = if mouths.is_empty() {
                    load.preview_step(field.normal_volume_flow_m3_s, gate)
                } else {
                    load.preview_step_with_mouth_flows(field.normal_volume_flow_m3_s, &mouth_flows, gate)
                }.map_err(NsfDgError::Thermoviscous)?;
                let mut reservoir_radiations = Vec::with_capacity(inputs.len());
                if !mouths.is_empty() {
                    let radiation = acoustic.independent_mouths_frame()
                        .ok_or(NsfDgError::Input("complete field lost its original mouth receipts"))?;
                    if radiation.mouths.len() != mouths.len() {
                        return Err(NsfDgError::Input("complete-field exterior receipt count changed"));
                    }
                    for (index, (&mouth, receipt)) in mouths.iter().zip(&radiation.mouths).enumerate() {
                        let audit = acoustic_interface::audit(mouth, receipt.mean_pressure_pa, receipt.load_work_j,
                            0.0, scalar_relation(index + 1), reference, dt, limits)?;
                        let resolved_exchange = *fluid.reservoir_exchange_at(inputs[index].0)
                            .ok_or(NsfDgError::Input("complete-field mouth lost its conservative NSF exchange"))?;
                        reservoir_radiations.push(NsfDgReservoirFrame {
                            resolved_exchange, ghost_pressure_offset_pa: physical[index + 1],
                            mean_acoustic_pressure_pa: mouth.mean_acoustic_pressure_pa,
                            outward_flow_m3_s: mouth.normal_volume_flow_m3_s,
                            radiation_pressure_pa: receipt.mean_pressure_pa,
                            radiation_work_j: receipt.load_work_j,
                            pressure_residual_pa: audit.pressure, projection_covariance_j: audit.covariance,
                            acoustic_work_residual_j: audit.work, mass_projection_difference_kg: audit.mass,
                        });
                    }
                }
                let source = upstream.as_mut().map(|owner| {
                    let face = coupled_exchange(&fluid, 1 - side)?;
                    owner.preview_step(1 - side, face, physical[width], rise,
                        upstream_relation.expect("original upstream field response"), limits, gate)
                }).transpose()?;
                accept_trial(fluid, acoustic, source, finite, side, scalar_relation(0), reference, limits,
                    dt, physical[0], drive, evaluations, step, reservoir_radiations, gate)
            },
        )?;
        self.joint.as_mut().expect("original complete-field nonlinear owner")
            .previous_coordinates = normalized;
        self.steps = step;
        Ok(frame)
    }
}
