//! Joint thermal and native endpoint-rate feedback on one original inventory.

use super::*;
use crate::conjugate::flow::phase_volume::momentum::resolve_endpoint_rheology_sources;

/// One accepted implicit thermal/rheological endpoint and one material advance.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitThermorheologicalThermalVolumeMomentumStep {
    /// Exact thermal source receipts, independent bulk response, used/actual
    /// molecular heat and aggregate thermal/native work retain their owners.
    pub evolution: ImplicitThermoviscousThermalVolumeMomentumStep,
    /// Cell RMS sqrt(2 dev(D):dev(D)) [1/s] at the accepted native molecular
    /// endpoint, before force/pressure and advection. Fixed solids are zero.
    pub endpoint_shear_rate_per_s: Vec<f64>,
    /// Whole coefficient-field updates across the original predictor and ALL
    /// thermal trials. The nested feedback source_evaluations counts individual
    /// mobile-cell queries; molecular_sweeps and Krylov work also include every
    /// inner nonlinear solve, including rejected thermal trials.
    pub rheology_source_updates: usize,
}

struct ThermorheologicalEndpoint {
    receipts: Vec<Option<ResolvedPhaseTransportState>>,
    bulk_pa_s: Vec<f64>,
    rates_per_s: Vec<f64>,
}

impl FlowField {
    /// Close source shear/bulk rheology and molecular heat at the same native
    /// implicit fluid/finite-solid thermal endpoint.
    ///
    /// Each heat trial starts from ORIGINAL fluid and solid histories. At that
    /// trial's exact pretransport phase, the native coupled momentum rows solve
    /// for endpoint velocity while `resolve(cell, phase, rate)` supplies the
    /// original shear source and independent bulk response [Pa s] at its
    /// current RMS sqrt(2 dev(D):dev(D)) [1/s]. The inner solve accepts only the
    /// true nonlinear row residual evaluated with those endpoint coefficients.
    /// Matching physical stress heat is then compared cellwise with the heat
    /// actually used by the implicit graph, through the same non-cancelling
    /// shared energy gates as the temperature-only entry.
    ///
    /// Original M/P and actual M/V inertia remain fixed across both iterations.
    /// Each thermal trial resolves its own rate-dependent response; coefficients
    /// from another thermal state are never accepted without a fresh solve.
    /// Pressure, force and conservative material transport run once after both
    /// endpoints converge. Temporal, scatter and projection losses remain
    /// numerical and never become molecular heat.
    ///
    /// The source must retain the exact supplied phase and chart density, with
    /// fully liquid mobility and no unrepresented Darcy resistance. The caller
    /// binds its own fixed pressure/composition and must consume BOTH supplied
    /// temperature/phase and physical shear coordinates. Incoming, heated and
    /// transported fluid must remain liquid; held solids keep separate storage.
    ///
    /// `max_sweeps` in `viscosity_controls` and the heat residual/linear limits
    /// in `solve` are WHOLE-CALL bounds, including the predictor and every inner
    /// nonlinear solve. Per-trial Newton/backtrack bounds and native Krylov
    /// bounds are unchanged. This remains a first-order thermal/mechanical
    /// split; no hidden substeps or nonlinear convergence guarantee are added.
    ///
    /// # Errors
    /// Unsupported sources or phases, bounded-work exhaustion, either nonlinear
    /// nonconvergence, conservation refusal or cancellation returns no fluid or
    /// finite-solid endpoint. Original input inventories remain unchanged.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<ImplicitThermorheologicalThermalVolumeMomentumStep, ChtError> {
        let curve = conductivity.caloric_curve();
        let mut source_updates = 0;
        let accepted = self.advance_implicit_thermal_feedback(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            boundaries,
            controls,
            viscosity_controls,
            solve,
            feedback,
            workspace,
            gate,
            |faces, mouth, heat, previous, work| {
                // Only the temperature-only API has a coefficient cache. A rate
                // law at another thermal endpoint must close its own native rows.
                drop(previous);
                let remaining = remaining_molecular_controls(viscosity_controls, work)?;
                // The phase is fixed throughout this inner momentum solve.
                // Recover each native caloric sample once, not at every sweep.
                let mut phases = Vec::with_capacity(domain.cell_count());
                for cell in 0..domain.cell_count() {
                    poll_every(gate, cell)?;
                    phases.push(if domain.is_fluid(cell) {
                        Some(feedback_phase_at(curve, input, heat, cell)?)
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
                    controls.evolution.mechanics.reconstruction,
                    controls.evolution.mechanics.projection.duration_s,
                    remaining,
                    gate,
                    |rates| {
                        let sources = resolve_endpoint_rheology_sources(
                            domain,
                            curve,
                            |cell| {
                                phases[cell].ok_or_else(|| {
                                    refuse("endpoint rheology needs a native fluid phase")
                                })
                            },
                            rates,
                            &mut resolve,
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
                    refuse("endpoint rheology returned no original source receipts")
                })?;
                Ok(PreparedThermalFeedback {
                    viscosity: prepared.viscosity,
                    sources: ThermorheologicalEndpoint {
                        receipts,
                        bulk_pa_s,
                        rates_per_s: prepared.endpoint_shear_rate_per_s,
                    },
                })
            },
        )?;
        Ok(ImplicitThermorheologicalThermalVolumeMomentumStep {
            evolution: ImplicitThermoviscousThermalVolumeMomentumStep {
                evolution: accepted.evolution,
                endpoint_sources: accepted.sources.receipts,
                bulk_viscosity_pa_s: accepted.sources.bulk_pa_s,
                used_molecular_heat_j: accepted.used_molecular_heat_j,
                feedback: accepted.feedback,
            },
            endpoint_shear_rate_per_s: accepted.sources.rates_per_s,
            rheology_source_updates: source_updates,
        })
    }
}
