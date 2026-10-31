//! Local material feedback on the original conservative enthalpy timeline.
//!
//! Each actual thermal substep resolves material transport at its cells' own
//! enthalpy-derived temperatures. The default quasisteady mode solves native
//! MAC momentum before advecting h. Explicit retained-momentum mode instead
//! advects h with the accepted old velocity and advances the same native
//! backward-Euler momentum owner over the actual accepted thermal interval.
//! Material viscosity, Boussinesq force and declared mushy drag share one
//! starting source epoch; reference mass and geometry remain fixed.

use std::collections::BTreeMap;

use fs_conduction::transient::deposition::advection::{Advected, AdvectionAudit};
use fs_conduction::transient::deposition::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer, series_conductance_w_k,
};
use fs_exec::CancelGate;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::phase::mushy::{MushyPhaseLaw, PhaseMobility};

use super::{FlowEnthalpyInput, admit_map, invalid, mass_graph};
use crate::conjugate::{
    ChtError, FluidProperties, FvFlow, SimpleConfig, Turbulence, ViscousHeatControls, ViscousPower,
    VoxelDomain, creeping_flow_with_body_force_heat, creeping_flow_with_darcy_heat,
    creeping_flow_with_viscous_heat, simple_flow_with_material_resistance,
    simple_flow_with_viscosity, simple_flow_with_viscosity_and_force,
};
use crate::material::transport::boussinesq::{BoussinesqConfig, BoussinesqMaterial};

/// Current and irreversible heat from these exact live-flow material samples.
pub mod electrical;

/// Reciprocal motional current, Lorentz momentum and physical heat on one timeline.
pub mod mhd;

/// Finite source-backed gas storage exchanging heat with the same phase/flow cells.
pub mod environment;

/// Material response to the accepted native symmetric-strain invariant.
pub mod rheology;
use rheology::RheologyBinding;
/// Original porous-skeleton sources on this same flow and heat timeline.
pub mod porous;
use porous::{ForchheimerResistance, PorousFlowSample};
use fs_material::porous::permeability::temperature::PermeabilityTemperatureLaw;

mod viscous;
use viscous::ViscousHeating;
pub use viscous::ViscousHeatingAudit;
mod phase_conduction;
use phase_conduction::PhaseConduction;
mod unsteady;
use crate::conjugate::unsteady::MaterialMomentum;

/// Opt-in inelastic freezing into an explicitly stationary support.
pub mod locking;

/// Conservative complete-melt release into the existing moving free surface.
pub mod release;

/// Exact stationary solid constraints on the same material/enthalpy owner.
pub mod phase;
pub use phase::{FlowTransportSnapshot, PhaseFlowTransport};
use phase::{fluid_snapshot, immobile_flow, immobile_power};

/// A coupled refusal retains the originating material, flow or history error.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialFlowError<E> {
    Thermal(DepositionError),
    Material { cell: CellId, source: MaterialError },
    Momentum(ChtError),
    Inner(E),
}

impl<E> From<DepositionError> for MaterialFlowError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<ChtError> for MaterialFlowError<E> {
    fn from(error: ChtError) -> Self {
        if error == ChtError::Cancelled {
            Self::Thermal(DepositionError::Cancelled)
        } else {
            Self::Momentum(error)
        }
    }
}
impl<E: core::fmt::Display> core::fmt::Display for MaterialFlowError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "material flow thermal owner: {error}"),
            Self::Material { cell, source } => write!(f, "material flow cell {cell:?}: {source}"),
            Self::Momentum(error) => write!(f, "material flow momentum owner: {error}"),
            Self::Inner(error) => write!(f, "material flow coupled history: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MaterialFlowError<E> {}

/// One concrete source result and the exact local query/transport snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowMaterialSample<F, T = FluidTransportState> {
    pub voxel: usize,
    pub cell: CellId,
    pub specific_enthalpy_j_kg: f64,
    pub temperature_k: f64,
    pub material: F,
    pub transport: T,
}

/// Source samples and the face field used by the last thermal substep.
/// Quasisteady mode solves this field from these samples. Retained-momentum
/// mode uses its accepted old field; these samples then drive the endpoint
/// momentum update exposed separately by MaterialFlowFeedback::unsteady_flow.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialFlowEpoch<F, T = FluidTransportState> {
    pub time_s: f64,
    pub samples: Vec<FlowMaterialSample<F, T>>,
    /// Cell-volume RMS sqrt(2 D:D) [1/s] of the accepted starting MAC state.
    /// Present only with explicit native rheology; these rates queried the
    /// SAME material coefficients used by momentum, conductivity and heat.
    pub shear_rate_per_s: Option<Vec<f64>>,
    pub flow: FvFlow,
    /// Applied `(rho(T)-rho0)*g` [N/m³] in native voxel order, zero on solids.
    /// None on the original fixed-density path; never an additional heat input.
    /// With unsteady heating, its signed native-face work is retained as a
    /// subset of total body work in unsteady_mechanical_step.
    pub buoyancy_force_n_m3: Option<Vec<[f64; 3]>>,
    /// Stationary support force ON exact-solid material, opposing its source-
    /// density perturbation load [N/m3]; zero on mobile/nonmaterial cells.
    /// Entries follow native voxel order; the load on the anchor has the opposite sign.
    /// Multiply by the original cell volume for force. This anchor has exactly
    /// zero velocity/work; the reference hydrostatic load is already removed.
    /// None without material buoyancy. It is never fed back as a fluid force.
    pub buoyancy_support_reaction_n_m3: Option<Vec<[f64; 3]>>,
    /// Explicit externally powered force density at this epoch, separate from
    /// buoyancy. Its signed work is in viscous_power for quasisteady heating;
    /// in unsteady_mechanical_step it is total body work minus the buoyancy
    /// subset. Work is never added as heat again.
    /// The MHD adapter adds its native Lorentz force to this applied RHS only;
    /// its report separates electrical transfer from original actuator work.
    pub external_force_n_m3: Option<Vec<[f64; 3]>>,
    /// Applied isotropic `mu/K(phi(h))` [Pa s/m2], in native voxel order.
    /// None without an explicitly attached stationary mushy closure.
    pub mushy_drag_pa_s_m2: Option<Vec<f64>>,
    /// Original skeleton sources sampled at the SAME local caloric temperature.
    pub porous_samples: Vec<PorousFlowSample>,
    /// Full symmetric world/grid-frame Darcy resistance [Pa s/m2], including
    /// any scalar mushy diagonal, when a sourced skeleton is bound. These
    /// coefficients drive both native force and local heat; do not add the
    /// mushy_drag_pa_s_m2 account a second time.
    pub darcy_resistance_pa_s_m2: Option<Vec<[[f64; 3]; 3]>>,
    /// Explicit material-principal rho*beta [kg/m4] and its original frame.
    /// No factor 1/2. Outer None when every source beta is zero; cell None for
    /// clear/immobile voxels. This nonlinear law is not a constant tensor.
    pub forchheimer_resistance: Option<Vec<Option<ForchheimerResistance>>>,
    /// Molecular heat and independent actuator work at this same source epoch.
    /// Present only after explicit quasisteady closed creeping-flow admission.
    /// Transient heat has a separate endpoint unsteady_mechanical_step receipt.
    pub viscous_power: Option<ViscousPower>,
}

/// Live local viscosity/momentum/enthalpy coupling around any existing observer.
///
/// `resolve(voxel, actual_cell, actual_kelvin)` returns an existing concrete
/// `FluidTransport` owner, including its source receipts. This receives h as
/// well as T, so a caller can explicitly admit a phase-dependent closure rather
/// than caching properties at a latent plateau's constant temperature. The
/// callback owns its material selection, fixed pressure/composition and phase
/// mobility admission. It must not mutate shared external physical history.
///
/// Ordinary entry density stays bitwise equal to the initial snapshot; the
/// Boussinesq entry explicitly projects bounded source density to rho0. T(h) and mass
/// stay with the actual registered caloric law; source Cp never replaces it.
/// Conductivity/contact evolution stays with the wrapped observer except for
/// explicitly assigned full fluid faces (see `Self::with_fluid_conduction`).
/// `Self::with_phase_conductivity` optionally scales those same half-cell
/// coefficients with the shared material law at the actual liquid mass fraction.
/// A solid/liquid fraction supplies mobility only through an explicitly
/// attached `Self::with_mushy_mobility` closure; it is never inferred from T.
///
/// By default each source epoch performs one laminar, quasisteady full-stress
/// momentum solve. `Self::with_unsteady_momentum` retains native face momentum
/// and explicitly selects a first-order transient split instead. The shared
/// `Advected` owner contributes its actual face mass rates,
/// boundary h flux, internal transfers and combined conductive/advective CFL.
/// Endpoint source checks run before publication. The applied field belongs to
/// its explicit starting epoch; transient endpoint momentum has a separate
/// receipt. Quasisteady mode requires rapid viscous equilibration relative to
/// thermal evolution. Both modes need timestep refinement and omit pressure-work
/// heat and free-surface motion. The optional unsteady Stokes heating mode
/// separately closes native kinetic storage, physical dissipation and actuator
/// work; ordinary retained momentum exposes only diagnostic kinetic energy.
/// Explicit
/// `Self::with_viscous_heating` admits wall-supplied molecular/Darcy heat in a
/// closed creeping-flow box; `Self::with_body_force_heating` additionally admits
/// explicit actuator work. The ordinary entry adds no viscous heat.
#[derive(Clone)]
pub struct MaterialFlowFeedback<'a, F, R, O, T = FluidTransportState> {
    domain: &'a VoxelDomain,
    input: FlowEnthalpyInput<'a>,
    config: &'a SimpleConfig,
    gate: &'a CancelGate,
    resolve: R,
    // The original liquid entry retains its initial reference exactly. A phase
    // entry has only the chart's fixed density until a mobile source exists.
    reference: Option<FluidTransportState>,
    density_kg_m3: f64,
    snapshot: fn(&F, Option<PhaseMobility>) -> Result<T, MaterialError>,
    exact_phase: bool,
    // Private endpoint-only deferral; never survives a published candidate.
    support_capture_endpoint_s: Option<f64>,
    buoyancy: Option<([f64; 3], fn(&F) -> f64)>,
    external_force_n_m3: Option<Vec<[f64; 3]>>,
    mushy: Option<&'a MushyPhaseLaw<'a>>,
    porous: Option<Vec<Option<PermeabilityTemperatureLaw<'a>>>>,
    porous_endpoint: Vec<PorousFlowSample>,
    phase_conduction: Option<PhaseConduction<'a>>,
    viscous_heating: Option<ViscousHeating>,
    unsteady: Option<MaterialMomentum>,
    rheology: Option<RheologyBinding<'a, F>>,
    bound_time_s: f64,
    fluid_conduction: BTreeMap<(CellId, CellId), FluidConductionFace>,
    advected: Option<Advected<O>>,
    applied: Option<MaterialFlowEpoch<F, T>>,
    endpoint: Vec<FlowMaterialSample<F, T>>,
    audit: AdvectionAudit,
    solves: usize,
    max_solves: usize,
}

impl<'a, F: FluidTransport + Clone, R, O: DepositionObserver> MaterialFlowFeedback<'a, F, R, O>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Bind one existing timeline/domain and bound momentum solves per advance.
    /// The map and reservoirs use the same admission as frozen solved flow.
    /// Each thermal substep costs one native solve; `max_solves` is separate
    /// from the original SIMPLEC iterations and thermal substep budgets.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        domain: &'a VoxelDomain,
        input: FlowEnthalpyInput<'a>,
        config: &'a SimpleConfig,
        inner: O,
        gate: &'a CancelGate,
        max_solves: usize,
        mut resolve: R,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(gate, &mut || Ok(()))?;
        if max_solves == 0 || config.turbulence != Turbulence::Laminar {
            return Err(invalid(
                "live material flow needs laminar momentum and a positive solve budget",
            )
            .into());
        }
        let endpoint = sample(
            thermal,
            domain,
            input,
            gate,
            &mut resolve,
            None,
            fluid_snapshot::<F>,
            &mut || Ok(()),
        )?;
        let reference = endpoint
            .first()
            .ok_or(invalid("live material flow needs an active fluid cell"))?
            .transport;
        admit_density(&endpoint, reference.density_kg_m3())?;
        admit_map(
            thermal,
            domain,
            reference.density_kg_m3(),
            input,
            &mut || Ok(()),
        )?;
        let advected = Advected::new(inner, Vec::new(), input.max_connections)?;
        Ok(Self {
            domain,
            input,
            config,
            gate,
            resolve,
            reference: Some(reference),
            density_kg_m3: reference.density_kg_m3(),
            snapshot: fluid_snapshot::<F>,
            exact_phase: false,
            support_capture_endpoint_s: None,
            buoyancy: None,
            external_force_n_m3: None,
            mushy: None,
            porous: None,
            porous_endpoint: Vec::new(),
            phase_conduction: None,
            viscous_heating: None,
            unsteady: None,
            rheology: None,
            bound_time_s: thermal.time_s(),
            fluid_conduction: BTreeMap::new(),
            advected: Some(advected),
            applied: None,
            endpoint,
            audit: AdvectionAudit::default(),
            solves: 0,
            max_solves,
        })
    }

    /// Add source-density natural/mixed convection to this same heat timeline.
    ///
    /// `reference` is the admitted hydrostatic material state: its rho0 stays
    /// in native inertia, the fixed cell masses and every advected mass flux.
    /// Each local source is evaluated once. The existing `BoussinesqMaterial`
    /// projection bounds rho/Cp departures, retains the original source and
    /// uses its mu and k. Its unprojected density supplies `(rho(T)-rho0)*g`
    /// to the SAME full-stress momentum solve. No independent beta load is
    /// also added. Pressure openings are relative to reference hydrostatics.
    ///
    /// The registered caloric object alone owns h/T/storage; its reference
    /// mass, fixed geometry, phase mobility and datum must be appropriate to
    /// this incompressible approximation. In particular this entry does not
    /// make latent-heat cells mobile or infer a melt permeability. Source
    /// departures are checked again at thermal endpoints before publication.
    /// The caller must justify the small-density approximation and the selected
    /// quasisteady or first-order transient momentum approximation. By default,
    /// mechanical dissipation does not enter the thermal audit. Selecting
    /// `Self::with_unsteady_viscous_heating` deposits physical molecular/Darcy
    /// losses once and separately accounts signed buoyancy work, native kinetic
    /// storage and numerical time loss. Gravity and pressure work are never
    /// additional heat sources; no gravitational-potential store is resolved.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn try_new_boussinesq(
        thermal: &DepositionSimulation<'_>,
        domain: &'a VoxelDomain,
        input: FlowEnthalpyInput<'a>,
        config: &'a SimpleConfig,
        inner: O,
        gate: &'a CancelGate,
        max_solves: usize,
        reference: FluidTransportState,
        buoyancy: BoussinesqConfig<3>,
        mut resolve: R,
    ) -> Result<
        MaterialFlowFeedback<
            'a,
            BoussinesqMaterial<F>,
            impl Clone
            + FnMut(usize, &DepositedCell, f64) -> Result<BoussinesqMaterial<F>, MaterialError>
            + use<F, R, O>,
            O,
        >,
        MaterialFlowError<O::Error>,
    > {
        buoyancy.validate().map_err(|_| {
            invalid("invalid Boussinesq gravity or relative material-departure bounds")
        })?;
        let bounded = move |voxel, cell: &DepositedCell, temperature| {
            BoussinesqMaterial::try_new(resolve(voxel, cell, temperature)?, &reference, buoyancy)
        };
        let mut coupled = MaterialFlowFeedback::try_new(
            thermal, domain, input, config, inner, gate, max_solves, bounded,
        )?;
        coupled.buoyancy = Some((buoyancy.gravity_m_s2, |material| {
            material.resolved_transport().density_kg_m3()
        }));
        Ok(coupled)
    }
}

impl<'a, F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'a, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Attach declared h-dependent stationary mixture mobility to this owner.
    ///
    /// The exact registered caloric chart supplies both T(h) and phase mass
    /// fraction. Equal constant phase densities are mandatory: this fixed-mass
    /// incompressible graph cannot absorb shrinkage into an unreported source.
    /// Actual local source viscosity multiplies the admitted inverse
    /// permeability, implicitly in native momentum. The same native face
    /// fluxes transport total h, including latent energy, exactly once.
    ///
    /// This is the declared homogeneous single-velocity enthalpy/porosity
    /// approximation, not separate moving-liquid/stationary-solid transport.
    /// It retains the underlying callback's source receipts and context; it
    /// does not promote raw permeability/model data to measured authority.
    /// Below the model's minimum liquid volume fraction (including fully
    /// solid) refuses at source and endpoint epochs. No finite artificial
    /// leakage is called solid locking. With explicit `Self::with_viscous_heating`,
    /// molecular and Darcy losses heat the same mixed-phase enthalpy cell at
    /// local thermal equilibrium. Otherwise drag dissipation is not thermalized.
    pub fn with_mushy_mobility(
        mut self,
        law: &'a MushyPhaseLaw<'a>,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        if self.exact_phase {
            return Err(invalid("exact phase mobility is bound by try_new_phase").into());
        }
        let rho = law.fixed_density_kg_m3().ok_or(invalid(
            "fixed-volume mushy flow requires equal constant phase densities",
        ))?;
        if rho.to_bits() != self.density_kg_m3.to_bits()
            || !core::ptr::addr_eq(
                self.input.caloric as *const dyn fs_conduction::transient::deposition::ThermalLaw,
                law.curve() as *const fs_material::phase::EquilibriumEnthalpyPhaseCurve,
            )
        {
            return Err(invalid(
                "mushy mobility must use the original caloric object and fluid density",
            )
            .into());
        }
        for h in self.input.inlet_specific_h_j_kg.into_iter().flatten() {
            law.at_specific_enthalpy(h).map_err(|_| {
                invalid("inlet enthalpy lies outside the admitted mushy mobility domain")
            })?;
        }
        self.mushy = Some(law);
        self.mushy_drag(&self.endpoint, &mut || Ok(()))?;
        Ok(self)
    }

    fn mushy_drag(
        &self,
        samples: &[FlowMaterialSample<F, T>],
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<Vec<f64>>, MaterialFlowError<O::Error>> {
        let Some(law) = self.mushy else {
            return Ok(None);
        };
        let mut drag = vec![0.0; self.domain.cell_count()];
        for sample in samples {
            check_gate(self.gate, check)?;
            let state = match sample.transport.phase_mobility() {
                Some(PhaseMobility::ImmobileSolid { phase, .. }) => {
                    if phase.temperature_k().to_bits() != sample.temperature_k.to_bits() {
                        return Err(invalid(
                            "solid constraint temperature differs from the original caloric owner",
                        )
                        .into());
                    }
                    continue;
                }
                Some(PhaseMobility::Mobile(state)) => *state,
                None => law
                    .at_specific_enthalpy(sample.specific_enthalpy_j_kg)
                    .map_err(|source| MaterialFlowError::Material {
                        cell: sample.cell,
                        source,
                    })?,
            };
            if state.phase.temperature_k().to_bits() != sample.temperature_k.to_bits() {
                return Err(invalid(
                    "mushy mobility temperature differs from the original caloric owner",
                )
                .into());
            }
            let transport = sample
                .transport
                .mobile_transport()
                .ok_or(invalid("mobile mushy state has no material transport"))?;
            let value = transport.dynamic_viscosity_pa_s() * state.inverse_permeability_per_m2;
            if !value.is_finite() || (state.inverse_permeability_per_m2 > 0.0 && value <= 0.0) {
                return Err(invalid("mushy Darcy resistance is not representable").into());
            }
            drag[sample.voxel] = value;
        }
        Ok(Some(drag))
    }

    /// Thermalize molecular and admitted stationary Darcy power in a Stokes box.
    ///
    /// This opt-in explicitly removes momentum convection from each native
    /// solve. Only tangential moving/stationary walls and symmetry are admitted;
    /// fixed no-slip solid voxels are supported. Openings, buoyancy, fans,
    /// turbulence and configured resistances refuse. `Self::with_mushy_mobility`
    /// supplies mu/K(phi(h)); `Self::with_porous_materials` supplies original
    /// permeability and optional calibrated Forchheimer drag. Their actual loss
    /// becomes heat in the ORIGINAL equilibrium cell. The fixed skeleton performs zero reaction work;
    /// its heat and latent energy belong to that same fixed-reference-mass h.
    /// Separate phase temperatures/velocities and moving skeletons are excluded.
    /// Source viscosity and the actual returned velocity must satisfy the declared
    /// global creeping-flow bound and independent actuator-work power check.
    /// `Self::with_body_force_heating` also accounts for prescribed body work.
    ///
    /// The original registered h/T law, source receipts and cell masses remain
    /// the sole thermal owners. Each actual substep adds molecular, Darcy and
    /// explicitly sourced Forchheimer heat
    /// once to the existing source map, alongside any wrapped source. The audit uses
    /// that step's accepted duration. Heat already appears in the original
    /// thermal audit's external_heat_j; never add actuator supply a second time.
    ///
    /// This is a fixed-volume, incompressible caloric approximation. There is no
    /// kinetic store or pressure-work thermalization. Wall drives are prescribed
    /// actuators, not finite mechanical bodies. Establish rapid viscous settling
    /// relative to thermal evolution and refine the original max_dt_s for the
    /// nonlinear viscosity/heat feedback; its thermal CFL alone is insufficient.
    pub fn with_viscous_heating(
        mut self,
        controls: ViscousHeatControls,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if self.buoyancy.is_some() {
            return Err(invalid("viscous heating excludes unaccounted buoyancy work").into());
        }
        if self.unsteady.is_some() {
            return Err(invalid(
                "steady viscous heating cannot account for retained kinetic energy; use a transient mechanical-energy owner",
            ).into());
        }
        crate::conjugate::simple::viscous_heat::admit_heating(self.domain, self.config, controls)?;
        self.viscous_heating = Some(ViscousHeating::new(controls));
        Ok(self)
    }

    /// Drive this same Stokes/material/enthalpy loop by external body forces.
    ///
    /// One finite SI force density [N/m3] per voxel is frozen over the advance.
    /// The original staggered momentum operator computes its work against the
    /// actual returned flow; molecular and stationary Darcy dissipation alone
    /// enter enthalpy. Material viscosity, phase fraction and permeability still
    /// refresh at every actual thermal substep. Signed actuator supply and heat
    /// are accepted together, including negative work by an opposing actuator.
    ///
    /// Closed creeping-flow restrictions and exact stationary solids remain.
    /// Buoyancy is excluded:
    /// no gravitational potential, finite electrical reservoir or kinetic store
    /// is debited by this prescribed actuator. Change a force only between
    /// accepted advances and account for its reported work in the enclosing
    /// source owner. No acceleration-to-force conversion invents another density.
    pub fn with_body_force_heating(
        self,
        force_density_n_m3: Vec<[f64; 3]>,
        controls: ViscousHeatControls,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if force_density_n_m3.len() != self.domain.cell_count() {
            return Err(invalid("body-force heating needs one explicit vector per voxel").into());
        }
        for force in &force_density_n_m3 {
            check_gate(self.gate, &mut || Ok(()))?;
            if force.iter().any(|v| !v.is_finite()) {
                return Err(invalid("body-force heating requires finite force density").into());
            }
        }
        let mut coupled = self.with_viscous_heating(controls)?;
        coupled.external_force_n_m3 = Some(force_density_n_m3);
        let domain = coupled.momentum_domain(&coupled.endpoint, &mut || Ok(()))?;
        coupled.check_phase_drives(domain.as_ref().unwrap_or(coupled.domain), &mut || Ok(()))?;
        Ok(coupled)
    }

    /// Assign existing full voxel-face links to molecular fluid conduction.
    ///
    /// The same applied transport samples supply both k(T) and momentum's
    /// viscosity. The original thermal owner uses the resulting series
    /// conductances for heat exchange and its combined conductive/advective CFL.
    /// No material is queried again, and no graph, contact area or heat source
    /// is created. Other links remain with their authored/wrapped owners.
    ///
    /// Each selected link must cover one entire face of two adjacent mapped
    /// fluid cells, in both contact-patch measures. Their total surface areas
    /// must match the native voxels. Parallel links on an owned face and an
    /// inner observer claiming the same conductance refuse. This models bulk
    /// isotropic conduction with no additional interface resistance; do not
    /// assign authored solid contacts, partial faces or contact-resistance laws.
    /// With `MaterialFlowFeedback::try_new_phase`, the original mapped material faces retain
    /// their sourced conductivity even when either cell becomes immobile.
    pub fn with_fluid_conduction(
        mut self,
        thermal: &DepositionSimulation<'_>,
        links: &[u64],
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        self.check_map(thermal, &mut || Ok(()))?;
        if links.len() > self.input.max_connections {
            return Err(DepositionError::Budget("material-flow conductive faces").into());
        }
        let voxels: BTreeMap<_, _> = self
            .input
            .cell_ids
            .iter()
            .enumerate()
            .filter_map(|(voxel, id)| id.map(|id| (id, voxel)))
            .collect();
        let mut faces = BTreeMap::new();
        for &id in links {
            check_gate(self.gate, &mut || Ok(()))?;
            let link = thermal
                .links()
                .get(&id)
                .ok_or(invalid("inactive material-flow conductive link"))?;
            let cells = ordered_cells(link.a, link.b);
            let a = *voxels
                .get(&cells.0)
                .ok_or(invalid("conductive link needs two mapped fluid cells"))?;
            let b = *voxels
                .get(&cells.1)
                .ok_or(invalid("conductive link needs two mapped fluid cells"))?;
            let axis = (0..6)
                .find(|&face| self.domain.neighbor(a, face) == Some(b))
                .ok_or(invalid("fluid conduction needs adjacent voxel faces"))?
                / 2;
            if faces
                .insert(
                    cells,
                    FluidConductionFace {
                        link: id,
                        voxels: [a, b],
                        axis,
                    },
                )
                .is_some()
            {
                return Err(invalid("duplicate material-flow conductive face").into());
            }
        }
        self.fluid_conduction = faces;
        self.check_fluid_conduction(thermal, &mut || Ok(()))?;
        Ok(self)
    }

    fn check_fluid_conduction(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        if self.fluid_conduction.is_empty() {
            return Ok(());
        }
        for (&cells, face) in &self.fluid_conduction {
            check_gate(self.gate, check)?;
            let patch = thermal
                .contact_patch(face.link)
                .ok_or(invalid("inactive material-flow conductive link"))?;
            let area = self.domain.face_area(face.voxels[0], face.axis);
            if ordered_cells(patch.link.a, patch.link.b) != cells
                || [
                    patch.link.area_m2,
                    patch.covered_area_m2[0],
                    patch.covered_area_m2[1],
                ]
                .into_iter()
                .any(|value| !same_measure(value, area))
            {
                return Err(invalid(
                    "fluid conduction must retain both full-face coverages",
                ));
            }
            for (id, voxel) in [cells.0, cells.1].into_iter().zip(face.voxels) {
                let surface = 2.0
                    * (0..3)
                        .map(|axis| self.domain.face_area(voxel, axis))
                        .sum::<f64>();
                let cell = thermal
                    .cells()
                    .get(&id)
                    .ok_or(invalid("inactive conductive cell"))?;
                if !same_measure(cell.geometry.surface_area_m2, surface) {
                    return Err(invalid(
                        "conducting fluid surface must match the native voxel",
                    ));
                }
            }
        }
        for link in thermal.links().values() {
            check_gate(self.gate, check)?;
            if self
                .fluid_conduction
                .get(&ordered_cells(link.a, link.b))
                .is_some_and(|face| face.link != link.id)
            {
                return Err(invalid(
                    "parallel thermal links duplicate an owned fluid face",
                ));
            }
        }
        Ok(())
    }

    /// Source-backed inputs and the actual field used for the last heat-advection
    /// step. In transient mode this is the starting field, not the new endpoint
    /// momentum; see `Self::unsteady_flow` and `Self::unsteady_step`.
    pub const fn applied(&self) -> Option<&MaterialFlowEpoch<F, T>> {
        self.applied.as_ref()
    }
    /// Material evaluations at the accepted thermal endpoint, including h.
    pub fn endpoint(&self) -> &[FlowMaterialSample<F, T>] {
        &self.endpoint
    }
    /// Summed original advection accounting over the last successful advance.
    /// Its boundary energy is already in the thermal audit; never add it twice.
    pub const fn audit(&self) -> AdvectionAudit {
        self.audit
    }
    /// Wall work and the same heat already included in the thermal owner.
    /// Resets per advance; failed candidates cannot publish these counters.
    pub fn viscous_heating_audit(&self) -> Option<ViscousHeatingAudit> {
        self.viscous_heating.as_ref().map(|owner| owner.audit)
    }
    /// Number of actual momentum solves in the last successful advance.
    pub const fn momentum_solves(&self) -> usize {
        self.solves
    }
    /// The original accepted source, contact and observation history.
    pub fn inner(&self) -> &O {
        self.advected.as_ref().expect("accepted observer").inner()
    }
    pub fn into_inner(self) -> O {
        self.advected.expect("accepted observer").into_inner()
    }

    fn check_map(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        check_gate(self.gate, check)?;
        admit_map(thermal, self.domain, self.density_kg_m3, self.input, check)?;
        self.check_fluid_conduction(thermal, check)?;
        Ok(())
    }
}

impl<F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot> DepositionObserver
    for MaterialFlowFeedback<'_, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    type Error = MaterialFlowError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("unfinished material-flow source epoch"))?
            .begin()
            .map_err(MaterialFlowError::Inner)?;
        self.audit = AdvectionAudit::default();
        if let Some(owner) = &mut self.viscous_heating {
            owner.begin();
        }
        if let Some(owner) = &mut self.unsteady {
            owner.begin_advance()?;
        }
        self.solves = 0;
        Ok(())
    }

    fn requires_step_interval(&self) -> bool {
        self.unsteady
            .as_ref()
            .is_some_and(MaterialMomentum::has_heating)
            || self
                .advected
                .as_ref()
                .is_some_and(|inner| inner.requires_step_interval())
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check_gate(self.gate, check)?;
        // Always forward duration-dependent inner sources, including when
        // this observer has no physical momentum-heating attachment.
        let mut sources = self
            .advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .step_sources(thermal, start_s, dt_s, end_s, check)
            .map_err(MaterialFlowError::Inner)?;
        let heat = self.prepare_unsteady_heat(thermal, start_s, dt_s, end_s, check)?;
        for (cell, watts) in heat {
            check_gate(self.gate, check)?;
            let value = sources.entry(cell).or_insert(0.0);
            *value += watts;
            if !value.is_finite() {
                return Err(DepositionError::Arithmetic.into());
            }
        }
        check_gate(self.gate, check)?;
        Ok(sources)
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_map(thermal, check)?;
        self.check_unsteady_start(thermal, check)?;
        let samples = self.sample_materials(thermal, check)?;
        admit_density(&samples, self.density_kg_m3)?;
        self.check_phase_conduction(&samples, check)?;
        let active_domain = self.momentum_domain(&samples, check)?;
        let domain = active_domain.as_ref().unwrap_or(self.domain);
        self.check_phase_drives(domain, check)?;
        let native_solve = domain.fluid_count() != 0;
        if native_solve && self.solves >= self.max_solves {
            return Err(DepositionError::Budget("live material momentum solves").into());
        }
        // No liquid reference is invented at an all-solid epoch. Once mobile,
        // use an actual retained mobile source for scalar parameter admission.
        let reference = self.reference.or_else(|| {
            samples
                .iter()
                .find_map(|sample| sample.transport.mobile_transport())
        });
        let mut viscosities = Vec::new();
        let fluid = if let Some(reference) = reference {
            viscosities = vec![reference.dynamic_viscosity_pa_s(); self.domain.cell_count()];
            for value in &samples {
                if let Some(transport) = value.transport.mobile_transport() {
                    viscosities[value.voxel] = transport.dynamic_viscosity_pa_s();
                }
            }
            // Caloric evolution remains in the original thermal law. Unused
            // solid slots satisfy native array admission only, never solid μ.
            Some(FluidProperties {
                density_kg_m3: self.density_kg_m3,
                specific_heat_j_kg_k: reference.specific_heat_j_kg_k(),
                conductivity_w_m_k: reference.thermal_conductivity_w_m_k(),
                kinematic_viscosity_m2_s: reference.kinematic_viscosity_m2_s(),
            })
        } else {
            None
        };
        let mut buoyancy_support_reaction_n_m3 = None;
        let buoyancy_force_n_m3 = if let Some((gravity, density)) = self.buoyancy {
            let mut forces = vec![[0.0; 3]; self.domain.cell_count()];
            let mut support = vec![[0.0; 3]; self.domain.cell_count()];
            for sample in &samples {
                check_gate(self.gate, check)?;
                let contrast = density(&sample.material) - self.density_kg_m3;
                let force = gravity.map(|g| contrast * g);
                if force.iter().any(|value| !value.is_finite()) {
                    return Err(invalid("material buoyancy force is not representable").into());
                }
                if sample.transport.mobile_transport().is_some() {
                    forces[sample.voxel] = force;
                } else {
                    support[sample.voxel] = force.map(|value| -value);
                }
            }
            buoyancy_support_reaction_n_m3 = Some(support);
            Some(forces)
        } else {
            None
        };
        let mushy_drag_pa_s_m2 = self.mushy_drag(&samples, check)?;
        let (porous_samples, darcy_resistance_pa_s_m2, forchheimer_resistance) =
            self.porous_drag(&samples, mushy_drag_pa_s_m2.as_deref(), check)?;
        let (flow, viscous_power) = if let Some(owner) = &self.unsteady {
            // Thermal CFL and advection use the accepted old momentum. Its
            // native update waits for the accepted-duration hook when heating
            // is enabled, otherwise for step(). Both use the actual thermal dt.
            // The accepted phase endpoint already admitted every newly fixed
            // face against the retained native velocities, with no clipping.
            (owner.flow.clone(), None)
        } else if !native_solve {
            // All constraints have zero velocity and every admitted drive does
            // zero work. This is an exact empty momentum system, not a solve.
            (
                immobile_flow(domain),
                self.viscous_heating.as_ref().map(|_| {
                    immobile_power(
                        domain,
                        mushy_drag_pa_s_m2.is_some() || darcy_resistance_pa_s_m2.is_some(),
                        self.external_force_n_m3.is_some(),
                    )
                }),
            )
        } else {
            let fluid = fluid
                .as_ref()
                .ok_or(invalid("mobile cells lack a transport source"))?;
            if let Some(owner) = &self.viscous_heating {
                let (flow, power) = if let Some(drag) = &darcy_resistance_pa_s_m2 {
                    crate::conjugate::simple::viscous_heat::creeping_material_heat_with_darcy_tensor(
                        domain, fluid, &viscosities, None, Some(drag), forchheimer_resistance.as_deref(),
                        self.external_force_n_m3.as_deref(), self.config, owner.controls, self.gate,
                    )?
                } else if let Some(force) = &self.external_force_n_m3 {
                    creeping_flow_with_body_force_heat(
                        domain,
                        fluid,
                        &viscosities,
                        force,
                        mushy_drag_pa_s_m2.as_deref(),
                        self.config,
                        owner.controls,
                        self.gate,
                    )?
                } else if let Some(drag) = &mushy_drag_pa_s_m2 {
                    creeping_flow_with_darcy_heat(
                        domain,
                        fluid,
                        &viscosities,
                        drag,
                        self.config,
                        owner.controls,
                        self.gate,
                    )?
                } else {
                    creeping_flow_with_viscous_heat(
                        domain,
                        fluid,
                        &viscosities,
                        self.config,
                        owner.controls,
                        self.gate,
                    )?
                };
                (flow, Some(power))
            } else if let Some(drag) = &darcy_resistance_pa_s_m2 {
                (
                    crate::conjugate::simple::simple_flow_with_darcy_tensor(
                        domain, fluid, &viscosities, drag, forchheimer_resistance.as_deref(),
                        buoyancy_force_n_m3.as_deref(), self.config, self.gate,
                    )?,
                    None,
                )
            } else if let Some(drag) = &mushy_drag_pa_s_m2 {
                (
                    simple_flow_with_material_resistance(
                        domain,
                        fluid,
                        &viscosities,
                        drag,
                        buoyancy_force_n_m3.as_deref(),
                        self.config,
                        self.gate,
                    )?,
                    None,
                )
            } else if let Some(forces) = &buoyancy_force_n_m3 {
                (
                    simple_flow_with_viscosity_and_force(
                        domain,
                        fluid,
                        &viscosities,
                        forces,
                        self.config,
                        self.gate,
                    )?,
                    None,
                )
            } else {
                (
                    simple_flow_with_viscosity(
                        domain,
                        fluid,
                        &viscosities,
                        self.config,
                        self.gate,
                    )?,
                    None,
                )
            }
        };
        check_gate(self.gate, check)?;
        self.admit_forchheimer_speed(&porous_samples, &flow, check)?;
        let (flows, boundaries) =
            mass_graph(domain, &flow.field, self.density_kg_m3, self.input, check)?;
        let old = self
            .advected
            .take()
            .ok_or(invalid("unfinished material-flow source epoch"))?;
        // Rebuilding the graph keeps the accepted wrapped physical history.
        // It does not call begin() again or reset another owner's work counter.
        let mut advected = Advected::new(old.into_inner(), flows, self.input.max_connections)?
            .with_boundaries(boundaries, self.input.max_connections)?;
        let mut sources = advected
            .sources(thermal, check)
            .map_err(MaterialFlowError::Inner)?;
        if let Some(power) = &viscous_power {
            for sample in &samples {
                check_gate(self.gate, check)?;
                let heat = power.cell_heat_w[sample.voxel];
                if heat != 0.0 {
                    let value = sources.entry(sample.cell).or_insert(0.0);
                    *value += heat;
                    if !value.is_finite() {
                        return Err(DepositionError::Arithmetic.into());
                    }
                }
            }
        }
        self.advected = Some(advected);
        self.reference = reference;
        self.applied = Some(MaterialFlowEpoch {
            time_s: thermal.time_s(),
            samples,
            shear_rate_per_s: self.rheology.as_ref().map(|source| source.rates.clone()),
            flow,
            buoyancy_force_n_m3,
            buoyancy_support_reaction_n_m3,
            external_force_n_m3: self.external_force_n_m3.clone(),
            mushy_drag_pa_s_m2,
            porous_samples,
            darcy_resistance_pa_s_m2,
            forchheimer_resistance,
            viscous_power,
        });
        self.solves += usize::from(native_solve);
        Ok(sources)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        let mut values = self
            .advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .conductances(thermal, check)
            .map_err(MaterialFlowError::Inner)?;
        if self.fluid_conduction.is_empty() {
            return Ok(values);
        }
        self.check_map(thermal, check)?;
        let applied = self.applied.as_ref().ok_or(invalid(
            "fluid conduction requires the applied material source epoch",
        ))?;
        if applied.time_s.to_bits() != thermal.time_s().to_bits() {
            return Err(invalid("stale fluid-conduction source epoch").into());
        }
        let mut samples = BTreeMap::new();
        for sample in &applied.samples {
            check_gate(self.gate, check)?;
            let cell = thermal
                .cells()
                .get(&sample.cell)
                .ok_or(invalid("inactive conductive cell"))?;
            if cell.specific_h_j_kg.to_bits() != sample.specific_enthalpy_j_kg.to_bits() {
                return Err(
                    invalid("fluid conduction must use the actual applied enthalpy state").into(),
                );
            }
            samples.insert(sample.cell, sample);
        }
        for (&cells, face) in &self.fluid_conduction {
            check_gate(self.gate, check)?;
            if values.contains_key(&face.link) {
                return Err(invalid("multiple owners supplied one fluid-face conductance").into());
            }
            let mut k = [0.0; 2];
            for (i, id) in [cells.0, cells.1].into_iter().enumerate() {
                let sample = samples.get(&id).ok_or(invalid(
                    "fluid conduction lacks an applied transport sample",
                ))?;
                k[i] = self.effective_conductivity(sample)?;
            }
            let distances = face
                .voxels
                .map(|voxel| 0.5 * self.domain.widths(voxel)[face.axis]);
            let g = series_conductance_w_k(
                self.domain.face_area(face.voxels[0], face.axis),
                distances,
                k,
            )?;
            values.insert(face.link, g);
        }
        check_gate(self.gate, check)?;
        Ok(values)
    }

    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .reservoir_conductances(thermal, check)
            .map_err(MaterialFlowError::Inner)
    }

    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .heat_transfers(thermal, check)
            .map_err(MaterialFlowError::Inner)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?
            .limit_timestep_s(thermal, proposed_s, check)
            .map_err(MaterialFlowError::Inner)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.check_map(thermal, &mut || Ok(()))?;
        self.check_unsteady_start(thermal, &mut || Ok(()))?;
        self.endpoint = self.sample_materials(thermal, &mut || Ok(()))?;
        admit_density(&self.endpoint, self.density_kg_m3)?;
        self.mushy_drag(&self.endpoint, &mut || Ok(()))?;
        self.check_porous_endpoint(&mut || Ok(()))?;
        let domain = self.momentum_domain(&self.endpoint, &mut || Ok(()))?;
        self.check_phase_drives(domain.as_ref().unwrap_or(self.domain), &mut || Ok(()))?;
        self.check_phase_conduction(&self.endpoint, &mut || Ok(()))?;
        self.advected
            .as_mut()
            .ok_or(invalid("unfinished material-flow source epoch"))?
            .event(event, thermal)
            .map_err(MaterialFlowError::Inner)
    }

    fn step(
        &mut self,
        start_s: f64,
        dt: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.check_map(thermal, check)?;
        self.endpoint = self.sample_materials(thermal, check)?;
        admit_density(&self.endpoint, self.density_kg_m3)?;
        self.mushy_drag(&self.endpoint, check)?;
        self.check_porous_endpoint(check)?;
        let domain = self.momentum_domain(&self.endpoint, check)?;
        self.check_phase_drives(domain.as_ref().unwrap_or(self.domain), check)?;
        self.check_phase_conduction(&self.endpoint, check)?;
        let advected = self
            .advected
            .as_mut()
            .ok_or(invalid("missing material-flow source epoch"))?;
        advected
            .step(start_s, dt, old, thermal, check)
            .map_err(MaterialFlowError::Inner)?;
        add_audit(&mut self.audit, advected.audit())?;
        self.advance_unsteady(start_s, dt, old, thermal, check)?;
        if self.exact_phase {
            if let Some(owner) = &mut self.unsteady {
                // Admit the ENDPOINT constraints against the actual native
                // endpoint velocities, never against a cell-centre average
                // or flux rounded to zero. Melting adds zero-speed faces;
                // freezing must remove only zero-energy momentum entries.
                if self.support_capture_endpoint_s != Some(thermal.time_s()) {
                    owner.accept_phase_domain(domain.as_ref().unwrap_or(self.domain), self.gate)?;
                }
            }
        }
        if self.rheology.is_some() {
            // The initial sampling above admits thermal/phase changes before
            // native momentum. Now re-admit the ACTUAL endpoint velocity and T;
            // these coefficients drive only the next accepted interval.
            self.endpoint = self.sample_materials(thermal, check)?;
            admit_density(&self.endpoint, self.density_kg_m3)?;
            self.mushy_drag(&self.endpoint, check)?;
            self.check_porous_endpoint(check)?;
            self.check_phase_conduction(&self.endpoint, check)?;
        }
        if let Some(owner) = &mut self.viscous_heating {
            let applied = self
                .applied
                .as_ref()
                .ok_or(invalid("missing viscous-power source epoch"))?;
            if applied.time_s.to_bits() != start_s.to_bits() {
                return Err(invalid("stale viscous-power source epoch").into());
            }
            owner.accept(
                applied
                    .viscous_power
                    .as_ref()
                    .ok_or(invalid("missing applied molecular power"))?,
                start_s,
                dt,
                thermal.time_s(),
            )?;
        }
        check_gate(self.gate, check)?;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct FluidConductionFace {
    link: u64,
    voxels: [usize; 2],
    axis: usize,
}

fn ordered_cells(a: CellId, b: CellId) -> (CellId, CellId) {
    if a < b { (a, b) } else { (b, a) }
}

fn same_measure(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

fn check_gate(
    gate: &CancelGate,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    check()?;
    if gate.is_requested() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn sample<F, R, E, T: FlowTransportSnapshot>(
    thermal: &DepositionSimulation<'_>,
    domain: &VoxelDomain,
    input: FlowEnthalpyInput<'_>,
    gate: &CancelGate,
    resolve: &mut R,
    phase_law: Option<&MushyPhaseLaw<'_>>,
    snapshot: fn(&F, Option<PhaseMobility>) -> Result<T, MaterialError>,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<FlowMaterialSample<F, T>>, MaterialFlowError<E>>
where
    R: FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    check_gate(gate, check)?;
    if domain.cell_count() > input.max_cells {
        return Err(DepositionError::Budget("live material flow voxels").into());
    }
    if input.cell_ids.len() != domain.cell_count() {
        return Err(invalid("live material flow cell map differs from its domain").into());
    }
    let mut samples = Vec::with_capacity(domain.fluid_count());
    for (voxel, id) in input.cell_ids.iter().enumerate() {
        if !domain.is_fluid(voxel) {
            continue;
        }
        check_gate(gate, check)?;
        let cell = id
            .and_then(|id| thermal.cells().get(&id))
            .ok_or(invalid("every flowing voxel needs an active enthalpy cell"))?;
        let temperature_k = thermal.temperature_k(cell.id)?;
        // Classify the original h before requesting any liquid transport.
        // An unsupported phase therefore cannot trigger a source extrapolation.
        let mobility = phase_law
            .map(|law| law.classify_at_specific_enthalpy(cell.specific_h_j_kg))
            .transpose()
            .map_err(|source| MaterialFlowError::Material {
                cell: cell.id,
                source,
            })?;
        let material =
            resolve(voxel, cell, temperature_k).map_err(|source| MaterialFlowError::Material {
                cell: cell.id,
                source,
            })?;
        let transport =
            snapshot(&material, mobility).map_err(|source| MaterialFlowError::Material {
                cell: cell.id,
                source,
            })?;
        samples.push(FlowMaterialSample {
            voxel,
            cell: cell.id,
            specific_enthalpy_j_kg: cell.specific_h_j_kg,
            temperature_k,
            material,
            transport,
        });
    }
    check_gate(gate, check)?;
    Ok(samples)
}

fn admit_density<F, T: FlowTransportSnapshot>(
    samples: &[FlowMaterialSample<F, T>],
    rho: f64,
) -> Result<(), DepositionError> {
    if samples.iter().any(|value| {
        value
            .transport
            .mobile_transport()
            .is_some_and(|transport| transport.density_kg_m3().to_bits() != rho.to_bits())
    }) {
        Err(invalid(
            "live flow density changes require compatible mass/volume and momentum evolution",
        ))
    } else {
        Ok(())
    }
}

fn add_audit(total: &mut AdvectionAudit, step: AdvectionAudit) -> Result<(), DepositionError> {
    for (sum, value) in [
        (&mut total.transported_mass_kg, step.transported_mass_kg),
        (
            &mut total.absolute_enthalpy_transfer_j,
            step.absolute_enthalpy_transfer_j,
        ),
        (&mut total.inlet_mass_kg, step.inlet_mass_kg),
        (&mut total.outlet_mass_kg, step.outlet_mass_kg),
        (&mut total.external_enthalpy_j, step.external_enthalpy_j),
    ] {
        *sum += value;
        if !sum.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
    }
    total.maximum_courant = total.maximum_courant.max(step.maximum_courant);
    total.substeps = total
        .substeps
        .checked_add(step.substeps)
        .ok_or(DepositionError::Arithmetic)?;
    Ok(())
}
