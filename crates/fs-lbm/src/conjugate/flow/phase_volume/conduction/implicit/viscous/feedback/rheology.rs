//! Joint thermal and native endpoint-rate feedback with optional material SGS.

use super::*;
use crate::conjugate::FvSgsConfig;
use crate::conjugate::flow::phase_volume::momentum::resolve_endpoint_rheology_sources;

/// One accepted implicit thermal/rheological endpoint and one material advance.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitThermorheologicalThermalVolumeMomentumStep {
    /// Exact thermal source receipts, independent bulk response, used/actual
    /// molecular heat and aggregate thermal/native work retain their owners.
    /// Modeled SGS transfer remains in the separate material-viscosity report.
    pub evolution: ImplicitThermoviscousThermalVolumeMomentumStep,
    /// Cell RMS sqrt(2 dev(D):dev(D)) [1/s] at the accepted native stress
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
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<ImplicitThermorheologicalThermalVolumeMomentumStep, ChtError> {
        self.advance_implicit_thermorheological_feedback(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            boundaries,
            None,
            None,
            controls,
            viscosity_controls,
            solve,
            feedback,
            workspace,
            gate,
            resolve,
        )
    }

    /// Close original temperature/rate rheology and resolved SGS at the same
    /// implicit fluid/finite-solid thermal endpoint.
    ///
    /// Every thermal trial starts from the original inventories. Its native
    /// momentum solve resolves mu and independent zeta at the trial's exact
    /// phase and unweighted physical endpoint shear rate. The SGS closure then
    /// refreshes cell, interior-edge and reservoir-mouth stresses and caps with
    /// those molecular coefficients before the common nonlinear residual.
    /// Original `rho = M/V` supplies SGS dynamic stress and local shear caps;
    /// source density does not replace native inertia or the SGS density field.
    ///
    /// Only molecular dissipation supplies heat to the implicit graph. Modeled
    /// transfer remains separate in the accepted material-viscosity report;
    /// it is neither immediate heat nor an evolved subgrid-energy store.
    /// Conductivity retains the original molecular source. Shared finite-solid
    /// storage, baths, contacts and thermal work keep their existing owners.
    ///
    /// The original whole-call thermal and momentum work bounds include every
    /// rejected trial and its inner source/SGS iterations. Once both endpoints
    /// converge, force, pressure and conservative material transport advance
    /// once. Optional opposed reservoir mouths retain their zero tangential
    /// trace and subsequent pressure load. SGS retains its exploratory policy;
    /// zero model coefficients follow the exact original thermorheological path.
    ///
    /// # Errors
    /// Model, source, phase, local cap, thermal or stress nonconvergence, work,
    /// conservation and cancellation refusals publish no partial endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        sgs: FvSgsConfig,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<ImplicitThermorheologicalThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        sgs.validate()?;
        self.advance_implicit_thermorheological_feedback(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            boundaries,
            Some(sgs),
            None,
            controls,
            viscosity_controls,
            solve,
            feedback,
            workspace,
            gate,
            resolve,
        )
    }

    /// Close SGS momentum and modeled fluid heat transport with the original
    /// temperature/rate rheology and finite-solid thermal inventories.
    ///
    /// The caller supplies a finite positive turbulent Prandtl number. On
    /// fluid-fluid faces, the original molecular conductivity is augmented by
    /// `k_t = (M/V) * Cp * nu_t / Pr_t` BEFORE combining the two half-cell
    /// resistances. The original caloric source supplies sensible Cp and its
    /// analytic enthalpy tangent; the existing sparse implicit heat graph owns
    /// flux, Jacobian, nonlinear work and conservative storage.
    ///
    /// Each thermal solve freezes one SGS field. The outer feedback refreshes
    /// native momentum and both molecular coefficients at that trial's exact
    /// phase/rate, and accepts only when molecular heating AND the change in
    /// every modeled face heat fit the shared non-cancelling energy bounds.
    /// Original histories never advance during retries. The feedback report
    /// retains used/accepted eddy fields, Pr_t and their explicit heat defect.
    ///
    /// Baths, source contacts and fluid-solid interfaces retain molecular heat
    /// transfer. No thermal wall function is inferred; near-wall resolution
    /// remains the caller's responsibility. This explicit isotropic gradient
    /// diffusion closure does not supply scalar-variance transport or a
    /// universal turbulent Prandtl law. Modeled momentum-energy transfer stays
    /// separate from caloric heat, exactly as in the momentum-only SGS entry.
    ///
    /// The same whole-call work, original M/V inertia, pressure-mouth ownership,
    /// phase admission and one final pressure/advection split apply. Validated
    /// zero model coefficients follow the exact original thermorheological path
    /// with no thermal SGS field allocation.
    ///
    /// # Errors
    /// Invalid Pr_t, unsupported caloric response, any source/cap/conservation
    /// refusal, bounded nonlinear exhaustion or cancellation returns no partial
    /// fluid or finite-solid endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs_transport(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        sgs: FvSgsConfig,
        turbulent_prandtl_number: f64,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<ImplicitThermorheologicalThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        sgs.validate()?;
        super::super::super::sgs::admit_prandtl_number(turbulent_prandtl_number)?;
        let thermal_prandtl = (sgs.model.coefficient() != 0.0)
            .then_some(turbulent_prandtl_number);
        self.advance_implicit_thermorheological_feedback(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            boundaries,
            Some(sgs),
            thermal_prandtl,
            controls,
            viscosity_controls,
            solve,
            feedback,
            workspace,
            gate,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_implicit_thermorheological_feedback(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        sgs: Option<FvSgsConfig>,
        turbulent_prandtl_number: Option<f64>,
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
            turbulent_prandtl_number,
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
                let resolve_coefficients = |rates: &[f64]| {
                    let sources = resolve_endpoint_rheology_sources(
                        domain,
                        curve,
                        None,
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
                };
                let prepared = match sgs {
                    Some(sgs) => self.prepare_material_rheology_with_source_and_sgs(
                        domain,
                        faces,
                        mouth,
                        input.volume.mass_kg,
                        input.momentum_kg_m_s,
                        boundaries,
                        controls.evolution.mechanics.reconstruction,
                        controls.evolution.mechanics.projection.duration_s,
                        sgs,
                        remaining,
                        gate,
                        resolve_coefficients,
                    ),
                    None => self.prepare_material_rheology_with_source(
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
                        resolve_coefficients,
                    ),
                }?;
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
