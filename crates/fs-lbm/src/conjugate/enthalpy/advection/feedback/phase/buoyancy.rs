//! Source buoyancy and stationary support on the existing phase/momentum owner.

use super::*;
use fs_material::phase::{EquilibriumPhaseState, transport::ResolvedPhaseTransportState};

/// Original phase source plus its explicit reference-mass transport projection.
///
/// Exactly solid material retains rho and k without inventing liquid Cp or mu.
/// Mobile material uses the shared Boussinesq projection; the original caloric
/// chart still supplies every temperature, phase fraction and stored enthalpy.
#[derive(Clone, Debug, PartialEq)]
pub struct BoussinesqPhaseMaterial {
    material: ResolvedPhaseTransportState,
    transport: Option<FluidTransportState>,
    reference: FluidTransportState,
    config: BoussinesqConfig<3>,
}

impl BoussinesqPhaseMaterial {
    /// Original unprojected phase/source response and selected source receipts.
    #[must_use]
    pub const fn material(&self) -> &ResolvedPhaseTransportState {
        &self.material
    }

    /// Reference density/Cp retained by native transport, never a replacement h law.
    #[must_use]
    pub const fn reference(&self) -> FluidTransportState {
        self.reference
    }

    #[must_use]
    pub const fn config(&self) -> BoussinesqConfig<3> {
        self.config
    }

    /// Refresh an original phase source while retaining this exact projection.
    ///
    /// Native rheology callbacks can resolve h/T/shear/source density on the
    /// original material, then call this method on the initial phase sample.
    /// Reference transport and bounds remain fixed; they cannot be accidentally
    /// replaced by the current source density or a reference shear value.
    /// The phase curve, mobility law and material identity must also survive.
    ///
    /// Claim selection and fixed nonthermal query context remain the original
    /// resolver's responsibility. Retain its ThermalMaterialPoint and selection;
    /// this projection does not infer which source coordinates may evolve.
    ///
    /// Exact solids still require no liquid Cp or viscosity. The shared
    /// Boussinesq owner bounds mobile rho/Cp while retaining actual mu and k;
    /// caloric mass and latent storage remain owned by the original h chart.
    pub fn at_source(
        &self,
        material: ResolvedPhaseTransportState,
    ) -> Result<Self, MaterialError> {
        if material.material().card_identity() != self.material.material().card_identity()
            || material.phase_state().phase_curve_identity()
                != self.material.phase_state().phase_curve_identity()
            || material.mobility().law_identity() != self.material.mobility().law_identity()
        {
            return Err(MaterialError::State {
                what: "phase buoyancy refresh must retain the original material, chart and mobility law"
                    .into(),
            });
        }
        Self::try_new(material, self.reference, self.config)
    }

    fn try_new(
        material: ResolvedPhaseTransportState,
        reference: FluidTransportState,
        config: BoussinesqConfig<3>,
    ) -> Result<Self, MaterialError> {
        let fail = |what: &str| MaterialError::State { what: what.into() };
        config
            .validate()
            .map_err(|error| fail(&error.to_string()))?;
        let density = reference.density_kg_m3();
        if material.phase_state().bulk_density_kg_m3().to_bits() != density.to_bits() {
            return Err(fail(
                "phase buoyancy must retain the original caloric mass density",
            ));
        }
        let departure = ((material.source_density_kg_m3() - density) / density).abs();
        if !departure.is_finite() || departure > config.max_relative_density_departure {
            return Err(fail(
                "phase buoyancy source density exceeds its admitted departure",
            ));
        }
        let transport = material
            .mobile_transport()
            .map(|actual| {
                BoussinesqMaterial::try_new(*actual, &reference, config)
                    .map(|bounded| *bounded.reference_storage_transport())
            })
            .transpose()?;
        Ok(Self {
            material,
            transport,
            reference,
            config,
        })
    }
}

impl FluidTransport for BoussinesqPhaseMaterial {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.transport.ok_or_else(|| MaterialError::State {
            what: "exactly solid phase has no liquid Boussinesq transport".into(),
        })
    }
}

impl PhaseTransport for BoussinesqPhaseMaterial {
    fn phase_state(&self) -> &EquilibriumPhaseState {
        self.material.phase_state()
    }

    fn transport_density_kg_m3(&self) -> f64 {
        self.reference.density_kg_m3()
    }

    fn thermal_conductivity_w_m_k(&self) -> Result<f64, MaterialError> {
        PhaseTransport::thermal_conductivity_w_m_k(&self.material)
    }
}

impl<'a, R, O: DepositionObserver>
    MaterialFlowFeedback<'a, ResolvedPhaseTransportState, R, O, PhaseFlowTransport>
where
    R: Clone
        + FnMut(usize, &DepositedCell, f64) -> Result<ResolvedPhaseTransportState, MaterialError>,
{
    /// Drive native natural convection as the SAME material melts and flows.
    ///
    /// `resolve_phase_transport_state_with_density_departure` retains source
    /// rho(T,phi), while this explicit projection holds reference mass, volume
    /// and the registered h chart fixed. Mobile rho/Cp departures use the shared
    /// Boussinesq bounds. Local mu and k remain unprojected. Phase classification
    /// precedes liquid-property resolution, including all-solid cold startup.
    ///
    /// Mobile cells receive `(rho_source-rho0)*g` in the ORIGINAL native momentum
    /// RHS alongside their declared mushy Darcy resistance. Exactly solid cells
    /// retain zero velocity; the opposite perturbation load on their stationary
    /// support-on-material reaction is reported in
    /// `buoyancy_support_reaction_n_m3`, with zero work.
    /// The reference hydrostatic load is already removed in this closed-domain
    /// pressure convention. Melting opens native faces at their retained zero
    /// speed without restarting the other velocities or thermal/material history.
    ///
    /// Select `with_unsteady_viscous_heating` for the existing native signed
    /// buoyancy-work, kinetic, physical-dissipation and numerical-loss account.
    /// Work is never a second heat input. This is a bounded reference-density,
    /// stationary-skeleton approximation: no shrinkage, solid elasticity,
    /// independent phase velocities or gravitational-potential store is added.
    /// Every endpoint source refusal uses the original clone-and-publish advance.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn try_new_phase_boussinesq(
        thermal: &DepositionSimulation<'_>,
        domain: &'a VoxelDomain,
        input: FlowEnthalpyInput<'a>,
        config: &'a SimpleConfig,
        inner: O,
        gate: &'a CancelGate,
        max_solves: usize,
        law: &'a MushyPhaseLaw<'a>,
        reference: FluidTransportState,
        buoyancy: BoussinesqConfig<3>,
        mut resolve: R,
    ) -> Result<
        MaterialFlowFeedback<
            'a,
            BoussinesqPhaseMaterial,
            impl Clone
            + FnMut(usize, &DepositedCell, f64) -> Result<BoussinesqPhaseMaterial, MaterialError>
            + use<R, O>,
            O,
            PhaseFlowTransport,
        >,
        MaterialFlowError<O::Error>,
    > {
        buoyancy
            .validate()
            .map_err(|_| invalid("invalid phase Boussinesq bounds"))?;
        if law.fixed_density_kg_m3().map(f64::to_bits) != Some(reference.density_kg_m3().to_bits())
        {
            return Err(
                invalid("phase buoyancy reference must equal the chart mass density").into(),
            );
        }
        let bounded = move |voxel, cell: &DepositedCell, temperature| {
            BoussinesqPhaseMaterial::try_new(
                resolve(voxel, cell, temperature)?,
                reference,
                buoyancy,
            )
        };
        let mut coupled = MaterialFlowFeedback::try_new_phase(
            thermal, domain, input, config, inner, gate, max_solves, law, bounded,
        )?;
        coupled.buoyancy = Some((buoyancy.gravity_m_s2, |sample| {
            sample.material.source_density_kg_m3()
        }));
        Ok(coupled)
    }
}
