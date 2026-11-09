//! Native thermal/endpoint-rate molecular stress with independent pore drag.

use super::super::super::super::viscous::feedback::{
    add_molecular_work, add_work, remaining_molecular_controls,
};
use super::*;
use crate::conjugate::flow::phase_volume::momentum::resolve_endpoint_rheology_sources;

/// One accepted joint thermal/rate molecular endpoint and phase-resistance step.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitThermorheologicalResistiveThermalVolumeMomentumStep {
    /// Retains exact molecular and Darcy source receipts, independent bulk,
    /// used/actual physical heats, finite-solid endpoints and aggregate work.
    pub evolution: ImplicitThermoviscousResistiveThermalVolumeMomentumStep,
    /// Native cell RMS sqrt(2 dev(D):dev(D)) [1/s], before Darcy/pressure,
    /// scattering and material advection. Fixed-solid entries are zero.
    pub endpoint_shear_rate_per_s: Vec<f64>,
    /// Whole coefficient-field updates over every molecular predictor and
    /// every rejected thermal/Darcy trial, not just the accepted endpoint.
    pub rheology_source_updates: usize,
}

impl FlowField {
    /// Close temperature and endpoint-rate molecular shear/bulk response on
    /// evolving mobile mushy density, implicit heat and explicit pore drag.
    ///
    /// Each native heat trial retains ORIGINAL fluid M/P/h and held-solid
    /// histories. Its exact caloric phases remain fixed while the existing
    /// native nonlinear stress solver refreshes `resolve_molecular` at every
    /// current physical RMS shear rate. Acceptance uses the true momentum-row
    /// residual and physical stress heat at those endpoint coefficients.
    /// The same cellwise physical heat then closes the implicit thermal loop.
    /// Every kick retains original M/V inertia even when source density changes.
    ///
    /// Molecular sources must preserve exact phase, density and native mobile
    /// law evidence. Independent bulk response [Pa s] supplies
    /// 2 mu dev(D) + zeta div(u) I. No rate cutoff, extrapolation or substitute
    /// viscosity is invented when the supplied material source refuses.
    ///
    /// `resolve_resistance` explicitly supplies the pore-scale Darcy response
    /// at transported phases. Resolved molecular shear does not determine a
    /// pore shear rate; neither zero nor a frozen molecular rate is fabricated
    /// for that callback. Both source families retain the same chart, density
    /// and mobility admission, with separate exact accepted receipts.
    ///
    /// The existing Darcy/pressure/volume owner consumes the converged molecular
    /// momentum and advances conservative M/Mh/P once. Molecular constitutive
    /// response closes before transport, Darcy after transport. This retains
    /// the declared first-order mechanical split, not unsplit endpoint flow.
    /// Physical molecular/Darcy heating, boundary work and numerical losses
    /// retain separate ledgers and the original shared non-cancelling gates.
    ///
    /// Molecular max_sweeps, feedback maximum_trials and thermal residual/
    /// Krylov limits are whole-call ceilings, including discarded nested work.
    /// No source is re-queried to construct the accepted receipt. Exact solids
    /// still require an immobile owner; this opt-in adds no moving skeleton or
    /// SGS response. Both resolvers own their fixed composition/pressure context.
    ///
    /// # Errors
    /// Source, phase, nonlinear convergence, work, conservation or cancellation
    /// refusal returns no fluid or finite-solid endpoint and preserves history.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_isotropic_implicit_conjugate_thermorheological_phase_resistance_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        mut resolve_molecular: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        )
            -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
        resolve_resistance: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ImplicitThermorheologicalResistiveThermalVolumeMomentumStep, ChtError> {
        let mut source_updates = 0;
        let (evolution, endpoint_shear_rate_per_s) = self
            .advance_implicit_resistive_thermal_feedback(
                domain,
                reservoir,
                conductivity,
                baths,
                solids,
                contacts,
                law,
                input,
                boundaries,
                controls,
                viscosity_controls,
                solve,
                feedback,
                workspace,
                gate,
                resolve_resistance,
                |faces, mouth, heat, previous, work| {
                    // A new thermal phase must close its own exact endpoint-rate
                    // momentum equations; a previous coefficient field is no proof.
                    drop(previous);
                    let remaining = remaining_molecular_controls(viscosity_controls, work)?;
                    let mut phases = Vec::with_capacity(domain.cell_count());
                    for cell in 0..domain.cell_count() {
                        poll_every(gate, cell)?;
                        phases.push(if domain.is_fluid(cell) {
                            Some(feedback_phase_at(law.curve(), input, heat, cell)?)
                        } else {
                            None
                        });
                    }
                    let mut endpoint = None;
                    let prepared = self.prepare_material_rheology_with_source(
                        domain,
                        faces,
                        mouth,
                        input.volume.mass_kg,
                        input.momentum_kg_m_s,
                        boundaries,
                        controls.resistance.thermal.mechanics.reconstruction,
                        controls.resistance.thermal.mechanics.projection.duration_s,
                        remaining,
                        gate,
                        |rates| {
                            let sources = resolve_endpoint_rheology_sources(
                                domain,
                                law.curve(),
                                Some(law),
                                |cell| {
                                    phases[cell].ok_or_else(|| {
                                        refuse("mushy endpoint rheology needs a native fluid phase")
                                    })
                                },
                                rates,
                                &mut resolve_molecular,
                                &mut work.source_evaluations,
                                gate,
                            )?;
                            endpoint = Some((sources.receipts, sources.bulk_pa_s.clone()));
                            Ok((sources.shear_pa_s, sources.bulk_pa_s))
                        },
                    )?;
                    add_work(&mut source_updates, prepared.source_updates)?;
                    add_molecular_work(work, &prepared.viscosity)?;
                    let (receipts, bulk_pa_s) = endpoint.ok_or_else(|| {
                        refuse("mushy endpoint rheology returned no original source receipts")
                    })?;
                    // Read immutable accepted receipts once, without another
                    // source query or cloning every shear field during each sweep.
                    let mut shear_pa_s = vec![0.0; domain.cell_count()];
                    for (cell, source) in receipts.iter().enumerate() {
                        poll_every(gate, cell)?;
                        if let Some(source) = source {
                            shear_pa_s[cell] = source
                                .mobile_transport()
                                .ok_or_else(|| {
                                    refuse("accepted molecular source lost mobile transport")
                                })?
                                .dynamic_viscosity_pa_s();
                        }
                    }
                    Ok(PreparedThermalFeedback {
                        viscosity: prepared.viscosity,
                        sources: (
                            ThermalViscositySources {
                                original: receipts,
                                shear_pa_s,
                                bulk_pa_s,
                            },
                            prepared.endpoint_shear_rate_per_s,
                        ),
                    })
                },
            )?;
        Ok(
            ImplicitThermorheologicalResistiveThermalVolumeMomentumStep {
                evolution,
                endpoint_shear_rate_per_s,
                rheology_source_updates: source_updates,
            },
        )
    }
}
