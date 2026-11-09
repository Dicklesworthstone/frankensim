//! Endpoint shear rheology on original material inertia and caloric history.

use super::*;
use crate::conjugate::FvSgsConfig;
use crate::conjugate::flow::material_momentum::PreparedMaterialRheology;
use fs_material::fluid_transport::viscous::NewtonianViscosity;

/// One source-resolved nonlinear stress step before pressure and advection.
#[derive(Clone, Debug, PartialEq)]
pub struct EndpointRheologyThermalVolumeMomentumStep {
    pub evolution: ViscousThermalVolumeMomentumStep,
    /// Native RMS sqrt(2 dev(D):dev(D)) at the accepted viscous endpoint [1/s].
    /// These are sampled before scattering impulses back to cell momentum.
    /// Fixed-solid slots are zero.
    pub endpoint_shear_rate_per_s: Vec<f64>,
    /// Original material receipts at that exact rate and the supplied ORIGINAL
    /// thermal phase. The separate thermal advance subsequently inserts heat.
    pub endpoint_sources: Vec<Option<ResolvedPhaseTransportState>>,
    /// Independent bulk response at the same phase and rate [Pa s].
    pub bulk_viscosity_pa_s: Vec<f64>,
    /// One original coefficient field plus one refresh after every sweep.
    pub source_updates: usize,
    /// Actual calls to the original material resolver across all updates.
    pub source_evaluations: usize,
}

/// Same original phase and rate admission for thermal feedback consumers.
pub(in crate::conjugate::flow) struct EndpointRheologySources {
    pub receipts: Vec<Option<ResolvedPhaseTransportState>>,
    pub shear_pa_s: Vec<f64>,
    pub bulk_pa_s: Vec<f64>,
}

impl FlowField {
    /// Resolve nonlinear shear and independent bulk stress at the native
    /// backward-Euler momentum endpoint, then advance physical heat once.
    ///
    /// The native solver keeps the ORIGINAL mass, momentum, face inertia and
    /// time history. Every coupled velocity sweep is followed by material
    /// resolution at its exact physical shear rate; all three momentum
    /// residuals and the final energy account use those refreshed coefficients.
    /// Acceptance therefore bounds the nonlinear stress balance, rather than
    /// merely a coefficient change or a frozen-coefficient solve.
    ///
    /// The resolver receives the original thermal phase and the native RMS
    /// sqrt(2 dev(D):dev(D)) in 1/s. It must consume that rate in the original
    /// source query and return its phase receipt plus independent bulk
    /// viscosity [Pa s]. Stress is 2 mu dev(D) + zeta div(u) I. Material
    /// validity, exact phase/density agreement and fully liquid mobility with
    /// no unrepresented Darcy support remain mandatory. No cutoff or
    /// extrapolation is supplied for an unsupported constitutive response.
    ///
    /// Physical molecular heat and wall work use the SAME accepted velocity
    /// and coefficients. Backward-Euler, scatter and later projection losses
    /// remain separate numerical channels. The subsequent original caloric,
    /// pressure and conservative transport advance runs once. Incoming,
    /// heated and transported fluid states must remain on the liquid branch.
    ///
    /// This resolves endpoint shear dependence at the original temperature;
    /// coupled temperature/rate feedback additionally needs the implicit thermal
    /// owner. Both retain the original first-order pressure/advection split.
    /// The existing viscosity sweep limit bounds the whole nonlinear solve;
    /// convergence is required within that budget, never assumed.
    ///
    /// Optional reservoirs retain the original opposed normal-traction pair
    /// and zero tangential trace. The separate pressure profile supplies
    /// mechanical pressure once, after the zero-pressure molecular predictor.
    ///
    /// # Errors
    /// Source, work, nonlinear residual, heat, phase, conservation or
    /// cancellation refusal returns no endpoint and changes no input history.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_piecewise_thermal_volume_momentum_with_endpoint_rheology(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<EndpointRheologyThermalVolumeMomentumStep, ChtError> {
        self.advance_endpoint_rheology_thermal_volume(
            domain,
            reservoir,
            curve,
            input,
            boundaries,
            None,
            controls,
            viscosity_controls,
            gate,
            resolve,
        )
    }

    /// Resolve original nonlinear molecular shear and independent bulk response
    /// together with SGS on the same native material-momentum endpoint.
    ///
    /// Each coupled sweep first observes the physical molecular shear rate
    /// without eddy-viscosity weighting, then resolves the original phase source
    /// mu(rate) and zeta. The selected closure refreshes its actual-M/V cell and
    /// reciprocal-edge stress before the common residual and work readback.
    /// Its molecular shear cap uses those refreshed source coefficients;
    /// independent bulk viscosity supplies no modeled bulk stress or cap relief.
    ///
    /// Only physical molecular dissipation heats the original phase chart.
    /// SGS transfer remains in the separate material-viscosity report. The
    /// original thermal state is held during the nonlinear predictor; caloric
    /// state, compatible volume, pressure and accepted mass-route transport
    /// subsequently advance once. This does not add implicit temperature/rate
    /// feedback. Optional opposed reservoirs retain the same zero tangential
    /// trace and single subsequent pressure load.
    ///
    /// The original phase/source support and sweep budgets remain mandatory.
    /// SGS retains its explicit exploratory policy and filter assumptions;
    /// a zero model coefficient follows the exact molecular endpoint path.
    ///
    /// # Errors
    /// Source, phase, model, local cap, nonlinear stress, physical-work, caloric,
    /// conservation or cancellation refusal publishes no partial endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_endpoint_rheology_and_sgs(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        sgs: FvSgsConfig,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<EndpointRheologyThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        sgs.validate()?;
        self.advance_endpoint_rheology_thermal_volume(
            domain,
            reservoir,
            curve,
            input,
            boundaries,
            Some(sgs),
            controls,
            viscosity_controls,
            gate,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_endpoint_rheology_thermal_volume(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        sgs: Option<FvSgsConfig>,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<EndpointRheologyThermalVolumeMomentumStep, ChtError> {
        let heat_budget = admit_viscous_thermal_input(domain, input, controls, gate)?;
        for phase in input.volume.boundary_inflow.iter().flatten() {
            admit_endpoint_rheology_phase(curve, *phase)?;
        }
        let (faces, profile, mouth) = match reservoir {
            Some((faces, profile)) => (
                faces,
                Some(profile),
                Some(Self::admit_material_viscosity_pressure_mouth(
                    faces,
                    boundaries,
                    input.boundary_velocity_m_s,
                )?),
            ),
            None => ([FlowFace::Wall; 6], None, None),
        };
        let mut source_evaluations = 0;
        let mut accepted_sources = None;
        let resolve_coefficients = |rates: &[f64]| {
            let source = resolve_endpoint_rheology_sources(
                domain,
                curve,
                |cell| {
                    input.volume.phase_states[cell]
                        .ok_or_else(|| refuse("endpoint rheology needs an original fluid phase"))
                },
                rates,
                &mut resolve,
                &mut source_evaluations,
                gate,
            )?;
            accepted_sources = Some((source.receipts, source.bulk_pa_s.clone()));
            Ok((source.shear_pa_s, source.bulk_pa_s))
        };
        let prepared: PreparedMaterialRheology = match sgs {
            Some(sgs) => self.prepare_material_rheology_with_source_and_sgs(
                domain,
                faces,
                mouth,
                input.volume.mass_kg,
                input.momentum_kg_m_s,
                boundaries,
                controls.mechanics.reconstruction,
                controls.mechanics.projection.duration_s,
                sgs,
                viscosity_controls,
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
                controls.mechanics.reconstruction,
                controls.mechanics.projection.duration_s,
                viscosity_controls,
                gate,
                resolve_coefficients,
            ),
        }?;
        let (endpoint_sources, bulk_viscosity_pa_s) = accepted_sources
            .ok_or_else(|| refuse("endpoint rheology requires accepted source receipts"))?;
        // The same rounded authored-plus-molecular heat used by the existing
        // caloric owner must stay liquid before advection can conceal a change.
        for cell in 0..domain.cell_count() {
            poll_every(gate, cell)?;
            if domain.is_fluid(cell) {
                let phase = input.volume.phase_states[cell]
                    .ok_or_else(|| refuse("missing original liquid rheology state"))?;
                let heat = input.volume.net_enthalpy_increment_j[cell]
                    + prepared.viscosity.report.cell_molecular_heat_j[cell];
                finite("flow.endpoint_rheology_heat", heat)?;
                let heated = curve
                    .advance_specific_energy(phase, heat / input.volume.mass_kg[cell])
                    .map_err(|error| refuse(error.to_string()))?;
                admit_endpoint_rheology_phase(curve, heated)?;
            }
        }
        let evolution = self.finish_viscous_thermal_volume(
            domain,
            faces,
            profile,
            curve,
            input,
            controls,
            viscosity_controls,
            prepared.viscosity,
            heat_budget,
            gate,
        )?;
        for (cell, phase) in evolution
            .evolution
            .thermal
            .volume
            .phase_states
            .iter()
            .enumerate()
        {
            poll_every(gate, cell)?;
            if let Some(phase) = phase {
                admit_endpoint_rheology_phase(curve, *phase)?;
            }
        }
        poll(gate)?;
        Ok(EndpointRheologyThermalVolumeMomentumStep {
            evolution,
            endpoint_shear_rate_per_s: prepared.endpoint_shear_rate_per_s,
            endpoint_sources,
            bulk_viscosity_pa_s,
            source_updates: prepared.source_updates,
            source_evaluations,
        })
    }
}

pub(in crate::conjugate::flow) fn admit_endpoint_rheology_phase(
    curve: &EquilibriumEnthalpyPhaseCurve,
    phase: EquilibriumPhaseState,
) -> Result<(), ChtError> {
    if phase.phase_curve_identity() != curve.identity()
        || phase.material_card_identity() != curve.material_card_identity()
        || phase.liquid_mass_fraction() != 1.0
    {
        return Err(refuse(
            "endpoint rheology needs the original fully liquid chart state",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(in crate::conjugate::flow) fn resolve_endpoint_rheology_sources(
    domain: &VoxelDomain,
    curve: &EquilibriumEnthalpyPhaseCurve,
    phase_at: impl Fn(usize) -> Result<EquilibriumPhaseState, ChtError>,
    rates: &[f64],
    resolve: &mut impl FnMut(
        usize,
        &EquilibriumPhaseState,
        f64,
    ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    evaluations: &mut usize,
    gate: &CancelGate,
) -> Result<EndpointRheologySources, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    if rates.len() != n {
        return Err(refuse(
            "endpoint rheology rates must match the original grid",
        ));
    }
    let mut receipts = Vec::with_capacity(n);
    let mut shear_pa_s = vec![0.0; n];
    let mut bulk_pa_s = vec![0.0; n];
    for cell in 0..n {
        poll_every(gate, cell)?;
        if !domain.is_fluid(cell) {
            receipts.push(None);
            continue;
        }
        let phase = phase_at(cell)?;
        admit_endpoint_rheology_phase(curve, phase)?;
        let rate = rates[cell];
        if !rate.is_finite() || rate < 0.0 {
            return Err(refuse("invalid native endpoint shear rate"));
        }
        poll(gate)?;
        *evaluations = evaluations
            .checked_add(1)
            .ok_or_else(|| refuse("endpoint rheology source count overflow"))?;
        let resolved = resolve(cell, &phase, rate);
        // A cancellation raised by the resolver takes precedence over its
        // simultaneous constitutive error, as for the original thermal owner.
        poll(gate)?;
        let (source, bulk) = resolved?;
        if *source.phase_state() != phase
            || source.source_density_kg_m3().to_bits() != phase.bulk_density_kg_m3().to_bits()
        {
            return Err(refuse(
                "rheology source must match the exact endpoint phase and density",
            ));
        }
        if source.mobility().mobile().is_none_or(|state| {
            state.inverse_permeability_per_m2 != 0.0 || state.volumes.liquid != 1.0
        }) {
            return Err(refuse(
                "endpoint rheology requires liquid without unrepresented Darcy support",
            ));
        }
        let mu = source
            .mobile_transport()
            .ok_or_else(|| refuse("liquid rheology source has no mobile transport"))?
            .dynamic_viscosity_pa_s();
        let response =
            NewtonianViscosity::new(mu, bulk).map_err(|error| refuse(error.to_string()))?;
        shear_pa_s[cell] = response.shear_viscosity_pa_s();
        bulk_pa_s[cell] = response.bulk_viscosity_pa_s();
        receipts.push(Some(source));
    }
    Ok(EndpointRheologySources {
        receipts,
        shear_pa_s,
        bulk_pa_s,
    })
}

#[cfg(test)]
mod tests;
