//! Original spatial implicit conduction and finite-solid storage pay native work.
//!
//! This uses the existing heat/work closure and original conjugate heat/flow
//! trial. Only accepted physical heat, fluid routes and finite solids publish.

use super::*;
use crate::conjugate::{
    ConductiveThermalVolumeMomentumControls, ImplicitConjugateThermalVolumeMomentumStep,
    ImplicitThermalVolumeSolveControls, ImplicitThermalVolumeWorkspace, ThermalVolumeBoundaries,
    ThermalVolumeContact, ThermalVolumeSolidCell,
};
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;

/// The original aggregate native thermal-work counter, with no added drag law.
pub use crate::conjugate::ImplicitResistanceHeatWork as ImplicitThermalVolumeHeatWork;

#[derive(Clone, Copy, Debug)]
pub struct ImplicitHeatWorkThermalVolumeControls {
    /// The same source-work iteration, original flow controls and physical
    /// gates. Its local budget here covers fluid AND finite-solid residuals.
    pub heat_work: SuppliedHeatThermalVolumeControls,
    pub maximum_conductive_faces: usize,
}

/// One accepted fluid/solid endpoint with independent physical heat closure.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitHeatWorkThermalVolumeStep {
    /// Original native conjugate result. Its thermal graph's authored source
    /// includes the mechanical debit: that NET graph input is not a second
    /// external physical heat supply. Source/solid/bath receipts are unchanged.
    pub evolution: ImplicitConjugateThermalVolumeMomentumStep,
    /// Accepted authored fluid increments BEFORE adding native Fourier heat.
    /// These are the external fluid heat minus the converged mechanical debit.
    pub authored_net_enthalpy_increment_j: Vec<f64>,
    pub cell_pressure_work_j: Vec<f64>,
    /// Fluid-only closure; supplied_heat_j includes accepted Fourier/contact
    /// heat into fluid. Use external_heat_j for the complete physical boundary.
    pub heat_work: SuppliedHeatThermalVolumeReport,
    pub represented_solid_enthalpy_change_j: f64,
    /// Authored fluid + authored solid + actual signed bath-channel heat [J].
    pub external_heat_j: f64,
    /// Noncancelling represented fluid and finite-solid heat residuals [J].
    pub absolute_local_heat_defect_j: f64,
    /// Actual fluid Delta H - carried boundary H + solid Delta H + native
    /// source-pressure work - physical external heat, in joules.
    pub total_heat_defect_j: f64,
    /// Includes rejected outer trials and every native nonlinear/linear trial.
    pub work: ImplicitThermalVolumeHeatWork,
}

impl FlowField {
    /// Close original spatial conduction, finite storage and expansion work.
    ///
    /// `supplied_heat_j` is physical authored heat into original fluid masses.
    /// `input.volume.net_enthalpy_increment_j` is an INITIAL authored NET guess
    /// before conduction, with the same meaning as the supplied-heat entry.
    /// Each source-work trial sends that guess to the ORIGINAL implicit
    /// conjugate graph, which solves k(T,f), contacts, surface laws and finite
    /// solid storage from their unchanged original states. Its accepted face
    /// heat then enters the SAME original thermal-volume/pressure/M-Mh-P step.
    ///
    /// The shared local closure is represented Delta H_i + dt*p_i*S_i =
    /// supplied Q_i + native conductive Q_i. Work therefore changes the
    /// temperature used by conduction during this same backward-Euler thermal
    /// interval. S includes the native mixing-volume response. Each trial owns
    /// one native heat solve and one complete native flow proposal; discarded
    /// trials contribute neither physical heat nor an advanced wall endpoint.
    ///
    /// The caller explicitly declares the prepared scalar conductivity
    /// isotropic. Original caloric/transport/source identities, latent storage,
    /// geometry and source support remain authoritative. Finite solids retain
    /// fixed mass and geometry, with no invented wall expansion or momentum.
    /// All thermal boundaries retain their actual separate bath channels.
    /// Only the accepted native work enters the physical heat debit.
    ///
    /// As in the supplied-heat entry, every fluid component needs an original
    /// zero-correction-pressure Free opening at the chart reference pressure.
    /// Only correction work beyond reference p/rho is charged; ambient p0*dV
    /// is not added again. No pressure-profile/EOS closure, numerical-loss
    /// heating, second pressure or thermal solver, or coupled entropy claim
    /// is introduced. Heat is implicit at the pre-advection endpoint; the
    /// original first-order conduction/flow split still requires refinement.
    ///
    /// The supplied residual-evaluation and total Krylov ceilings cover ALL
    /// outer trials. Newton/backtrack limits retain their per-trial meaning.
    /// Scratch may change after refusal; original fluid and solid inventories
    /// do not. Source, work, coupling, representation, native budget or
    /// cancellation refusal publishes no participant.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_isotropic_implicit_heat_work_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        supplied_heat_j: &[f64],
        controls: ImplicitHeatWorkThermalVolumeControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitHeatWorkThermalVolumeStep, ChtError> {
        let mut work = ImplicitThermalVolumeHeatWork::default();
        let accepted = close_native_heat_work(
            domain,
            conductivity.caloric_curve(),
            input,
            supplied_heat_j,
            controls.heat_work,
            gate,
            |authored| {
                let mut remaining = solve;
                remaining.maximum_residual_evaluations = solve
                    .maximum_residual_evaluations
                    .saturating_sub(work.residual_evaluations);
                remaining.linear.maximum_iterations = solve
                    .linear
                    .maximum_iterations
                    .saturating_sub(work.linear_iterations);
                if remaining.maximum_residual_evaluations == 0
                    || remaining.linear.maximum_iterations == 0
                {
                    return Err(refuse(
                        "implicit heat/work whole-call thermal work exhausted",
                    ));
                }
                let native = self
                    .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_baths(
                        domain,
                        faces,
                        None,
                        conductivity,
                        baths,
                        solids,
                        contacts,
                        ThermalVolumeMomentumInput {
                            volume: crate::conjugate::CoexistenceVolumeInput {
                                net_enthalpy_increment_j: authored,
                                ..input.volume
                            },
                            ..input
                        },
                        ConductiveThermalVolumeMomentumControls {
                            evolution: controls.heat_work.native,
                            max_conductive_faces: controls.maximum_conductive_faces,
                        },
                        remaining,
                        workspace,
                        gate,
                    )?;
                let solved = &native.report.conduction.solve;
                for (count, added) in [
                    (&mut work.heat_trials, 1),
                    (&mut work.nonlinear_iterations, solved.nonlinear_iterations),
                    (&mut work.linear_iterations, solved.linear_iterations),
                    (&mut work.residual_evaluations, solved.residual_evaluations),
                ] {
                    *count = count
                        .checked_add(added)
                        .ok_or_else(|| refuse("implicit heat/work count overflow"))?;
                }
                if work.residual_evaluations > solve.maximum_residual_evaluations
                    || work.linear_iterations > solve.linear.maximum_iterations
                {
                    return Err(refuse(
                        "implicit heat/work exceeded its aggregate thermal work",
                    ));
                }
                Ok(NativeHeatWorkTrial {
                    evolution: native.evolution,
                    conductive_heat_j: Some(native.conductive_enthalpy_increment_j),
                    receipt: (native.solid_endpoints, native.report),
                })
            },
        )?;
        let NativeHeatWorkTrial {
            evolution,
            conductive_heat_j,
            receipt: (solid_endpoints, report),
        } = accepted.trial;

        // Match actual returned finite-solid states to original fixed masses.
        // Sorting references happens only after acceptance, never per trial.
        let mut ordered: Vec<_> = solids.iter().collect();
        ordered.sort_unstable_by_key(|solid| solid.cell);
        if ordered.len() != solid_endpoints.len() {
            return Err(refuse("implicit heat/work lost a finite-solid endpoint"));
        }
        let mut solid_change = Dd::ZERO;
        let mut solid_authored = Dd::ZERO;
        let mut solid_local_defect = Dd::ZERO;
        for (index, (original, endpoint)) in ordered.iter().zip(&solid_endpoints).enumerate() {
            poll_cell(gate, index)?;
            if original.cell != endpoint.cell
                || original.mass_kg.to_bits() != endpoint.mass_kg.to_bits()
            {
                return Err(refuse(
                    "implicit heat/work changed a finite-solid inventory",
                ));
            }
            let change = Dd::from_f64(original.mass_kg)
                * (Dd::from_f64(endpoint.phase.specific_enthalpy_j_kg())
                    - Dd::from_f64(original.phase.specific_enthalpy_j_kg()));
            let authored = Dd::from_f64(original.net_enthalpy_increment_j);
            let residual =
                change - authored - Dd::from_f64(endpoint.conductive_enthalpy_increment_j);
            solid_change = solid_change + change;
            solid_authored = solid_authored + authored;
            solid_local_defect = solid_local_defect + magnitude(residual);
        }
        // Each accepted native bath channel appears once, including opposing
        // baths on a shared surface. Internal fluid/wall heat cancels globally.
        let mut external = solid_authored;
        for (cell, &heat) in supplied_heat_j.iter().enumerate() {
            poll_cell(gate, cell)?;
            external = external + Dd::from_f64(heat);
        }
        for (index, bath) in report.conduction.bath_exchanges.iter().enumerate() {
            poll_cell(gate, index)?;
            external = external + Dd::from_f64(bath.heat_j);
        }
        let volume = &evolution.thermal.volume;
        let total = Dd::from_f64(volume.report.material_enthalpy_change_j)
            - Dd::from_f64(volume.transport.report.carried_enthalpy_j)
            + Dd::from_f64(volume.transport.projection.report.pressure_source_work_j)
            + solid_change
            - external;
        let local = accepted.absolute_local_heat_defect + solid_local_defect;
        if [solid_change, external, total, local]
            .iter()
            .any(|v| !number(*v).is_finite())
            || exceeds(
                local,
                controls.heat_work.maximum_absolute_local_heat_defect_j,
            )
            || exceeds(total, controls.heat_work.maximum_total_heat_defect_j)
        {
            return Err(refuse(
                "actual fluid/solid storage, boundary enthalpy and native work exceed the heat gates",
            ));
        }
        poll(gate)?;
        Ok(ImplicitHeatWorkThermalVolumeStep {
            evolution: ImplicitConjugateThermalVolumeMomentumStep {
                evolution,
                conductive_enthalpy_increment_j: conductive_heat_j
                    .ok_or_else(|| refuse("implicit heat/work lost its native Fourier heat"))?,
                solid_endpoints,
                report,
            },
            authored_net_enthalpy_increment_j: accepted.authored_net_enthalpy_increment_j,
            cell_pressure_work_j: accepted.cell_pressure_work_j,
            heat_work: accepted.report,
            represented_solid_enthalpy_change_j: number(solid_change),
            external_heat_j: number(external),
            absolute_local_heat_defect_j: number(local),
            total_heat_defect_j: number(total),
            work,
        })
    }
}
