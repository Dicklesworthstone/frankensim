//! Motional induction in bulk conductors feeding the original winding circuit.
//!
//! Stationary stranded profiles retain their reciprocal terminal relation.
//! Solid sections additionally supply the back-EMF integral of J_unit*v.grad A
//! and impose their prescribed net current on the SAME moving Ohmic field.
//! Every current trial advances the actual moving
//! nonlinear field from accepted history, so its induced flux and screening
//! change the solved load current. No independently fitted speed voltage is
//! appended to the circuit. The field's conductor-frame current supplies both
//! Joule heat and Lorentz work, and private winding-port work cancels once.

use super::{
    FieldAccount, FieldTrial, NonlinearWindingCandidate, NonlinearWindingCircuit,
    NonlinearWindingControls, NonlinearWindingError, Result, admit_profiles, finite, poll, solid,
};
use crate::planar::eddy::{
    EddyConductor,
    nonlinear::{
        NonlinearEddyError,
        moving::{NonlinearMovingEddyControls, NonlinearMovingEddyStep},
    },
};
use crate::winding::ohmic::WindingConductor;

impl<'m> NonlinearWindingCircuit<'m> {
    /// Solve the original native circuit against moving, source-resolved bulk
    /// conductors and reversible nonlinear or stationary fixed-recoil magnets.
    ///
    /// Velocity is a complete held P1 nodal field [m/s]. The original field owner
    /// checks connectivity, prescribed regional current, magnetic Peclet
    /// and true nonlinear residual limits. Mesh, source laws, conductivity and
    /// material occupancy stay fixed; this operation does not transport an
    /// interface or permanent magnetization. Every stranded winding cell and
    /// every recoil-magnet cell must have zero velocity. Otherwise a fixed
    /// winding profile would omit its own motion-induced terminal work.
    ///
    /// The endpoint-current, backward-Euler field scheme is first order. The
    /// returned account closes delta(storage) + heat + numerical loss + signed
    /// mechanical work = external source work. Negative mechanical work is
    /// energy supplied by the imposed-motion owner, not negative Joule heat.
    /// Rejected trials publish no field, circuit, linkage, current or clock.
    pub fn step_with_eulerian_velocity_candidate(
        &self,
        target_s: f64,
        velocity_m_s: &[[f64; 2]],
        external_inputs: &[f64],
        controls: NonlinearWindingControls<NonlinearMovingEddyControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m, NonlinearMovingEddyStep>> {
        self.step_with_eulerian_velocity_and_sources_candidate(
            target_s,
            velocity_m_s,
            self.field.conductors().to_vec(),
            self.conductors.clone(),
            external_inputs,
            controls,
            cancelled,
        )
    }

    /// Refresh original conductivity/resistance sources for the SAME moving
    /// interval. Magnetic laws, geometry, material occupancy and end storage
    /// remain fixed. The field retains A; solid profiles follow n*sigma/S and
    /// terminal induction uses F_new^T delta(A), never delta(F^T A).
    ///
    /// E_conductor=nI/S-P_sigma(delta(A)/dt+v.grad A). The native voltage
    /// includes the conductivity-weighted motional term, even at zero current.
    /// Total local J supplies both Joule heat and Lorentz work exactly once.
    /// Supported Eulerian flow must leave material coefficients/occupancy
    /// invariant; finite interfaces, magnet motion and thermal advection need
    /// their corresponding transport owners. Every stranded cell stays still.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_eulerian_velocity_and_sources_candidate(
        &self,
        target_s: f64,
        velocity_m_s: &[[f64; 2]],
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls<NonlinearMovingEddyControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m, NonlinearMovingEddyStep>> {
        poll(cancelled)?;
        let mesh = self.field.problem().mesh();
        if velocity_m_s.len() != mesh.vertices().len() {
            return Err(NonlinearWindingError::Invalid(
                "complete finite nodal winding/bulk velocity required",
            ));
        }
        for velocity in velocity_m_s {
            poll(cancelled)?;
            if velocity.iter().any(|v| !v.is_finite()) {
                return Err(NonlinearWindingError::Invalid(
                    "complete finite nodal winding/bulk velocity required",
                ));
            }
        }
        for source in &self.conductors {
            poll(cancelled)?;
            if mesh.faces()[source.triangle]
                .iter()
                .any(|&vertex| velocity_m_s[vertex as usize] != [0.0, 0.0])
            {
                return Err(NonlinearWindingError::Invalid(
                    "fixed stranded winding cells require zero material velocity",
                ));
            }
        }
        let solid = if self.solid_regions.is_empty() {
            None
        } else {
            let solid = solid::prepare(
                &self.field,
                &core_conductors,
                self.profiles_per_m2(),
                &self.solid_regions,
                false,
                cancelled,
            )?;
            admit_profiles(
                mesh,
                &solid.profiles,
                self.field.problem().thickness_m(),
                self.profile_balance_tolerance_turns,
                cancelled,
            )?;
            Some(solid)
        };
        self.step_coupled(
            target_s,
            None,
            solid.as_ref().map(|s| s.profiles.clone()),
            solid.as_ref().map(|s| s.layout.resistance_ohm.as_slice()),
            winding_conductors,
            external_inputs,
            controls,
            |currents, j, cancelled| {
                let mut check = || cancelled();
                let candidate = if let Some(solid) = &solid {
                    let (drives, targets) = solid
                        .layout
                        .drives(currents, &mut check)
                        .map_err(NonlinearEddyError::from)?;
                    self.field
                        .step_with_eulerian_velocity_and_region_drives_candidate(
                            target_s,
                            velocity_m_s,
                            core_conductors.clone(),
                            j,
                            &drives,
                            &targets,
                            controls.field,
                            &mut check,
                        )?
                } else {
                    self.field
                        .step_with_eulerian_velocity_and_sources_candidate(
                            target_s,
                            velocity_m_s,
                            core_conductors.clone(),
                            j,
                            controls.field,
                            &mut check,
                        )?
                };
                let (state, step) = candidate.into_parts();
                let motional_voltage_v = if let Some(solid) = &solid {
                    let mut voltages = vec![0.0; currents.len()];
                    for binding in &solid.layout.bindings {
                        poll(&mut check)?;
                        let group = step
                            .region_motional_electric_field_v_m
                            .binary_search_by_key(&binding.region, |(region, _)| *region)
                            .map_err(|_| {
                                NonlinearWindingError::Invalid(
                                    "missing original solid motional field",
                                )
                            })?;
                        let motional = step.region_motional_electric_field_v_m[group].1;
                        let contribution = finite(
                            -self.field.problem().thickness_m()
                                * finite(binding.signed_turns * motional)?,
                        )?;
                        voltages[binding.winding] =
                            finite(voltages[binding.winding] + contribution)?;
                    }
                    Some(voltages)
                } else {
                    None
                };
                let account = FieldAccount {
                    stored_before_j: step.stored_before_j,
                    stored_after_j: step.stored_after_j,
                    heat_j: step.heat_j,
                    supplied_work_j: step.supplied_work_j,
                    uniform_drive_work_j: step.uniform_drive_work_j,
                    motional_drive_work_j: step.motional_drive_work_j,
                    numerical_dissipation_j: step.numerical_dissipation_j,
                    parameter_work_j: 0.0,
                    mechanical_work_j: step.mechanical_work_j,
                };
                Ok(FieldTrial {
                    state,
                    step,
                    account,
                    motional_voltage_v,
                })
            },
            cancelled,
        )
    }
}
