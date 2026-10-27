//! Pinned liquid/caloric sources on the original solved-flow enthalpy owner.
//!
//! The native face graph, upwind transport, CFL and energy audit remain in
//! `super`. These paths admit fully liquid material. The ordinary entries fix
//! density; `LiquidFlowEnthalpyInput::feedback_boussinesq` admits bounded source
//! density departures only for buoyancy while preserving reference cell mass.
//! Frozen advection additionally fixes viscosity; `LiquidFlowEnthalpyInput::feedback`
//! re-solves momentum from the evolving local viscosity on each thermal substep.
//! The borrowed sourced caloric curve remains authoritative H(T).

use core::fmt;
use std::collections::BTreeMap;

use fs_conduction::transient::deposition::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer, MaterialId,
};
use fs_exec::CancelGate;
use fs_matdb::EvaluationDecision;
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, require_frozen_coefficient_support,
};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_PROPERTY, ResolvedLiquidState,
};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{
    EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialStatePointError, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::Temperature;

use super::feedback::{MaterialFlowError, MaterialFlowFeedback};
use crate::conjugate::{SimpleConfig, VoxelDomain};
use crate::material::transport::boussinesq::{BoussinesqConfig, BoussinesqMaterial};
use fs_material::MaterialError;
use fs_material::fluid_transport::FluidTransport;

use super::{FlowEnthalpyInput, MaterialFlowAdvection, same_positive};
use crate::material::conjugate::MaterialChtFlow;

/// Pinned temperature/shear-rate liquid sources on retained native momentum.
pub mod rheology;

/// Native fixed-cell assignment with one explicitly sourced liquid/caloric law.
#[derive(Clone, Copy)]
pub struct LiquidFlowEnthalpyInput<'a> {
    /// Original voxel-to-thermal map; native geometry and mass admission applies.
    pub cell_ids: &'a [Option<CellId>],
    /// The material ID under which `caloric.curve()` is actually registered.
    pub material_id: MaterialId,
    /// Exact liquid law whose reference state supplied the momentum solve.
    pub liquid: &'a LiquidTemperatureLaw<'a>,
    /// Original borrowed H(T), including its Cp/density selections and context.
    pub caloric: &'a ResolvedHeatCapacityCurve,
    /// Incoming h [J/kg] in Face3::ALL order, in this same caloric datum.
    /// Every native inward/backflow face requires an explicit supported value.
    pub inlet_specific_h_j_kg: [Option<f64>; 6],
    /// Native bound on all voxels, including solids.
    pub max_cells: usize,
    /// Native bound on nonzero interior and boundary mass connections.
    pub max_connections: usize,
}

impl<'a> LiquidFlowEnthalpyInput<'a> {
    fn native(self) -> FlowEnthalpyInput<'a> {
        FlowEnthalpyInput {
            cell_ids: self.cell_ids,
            material_id: self.material_id,
            caloric: self.caloric.curve(),
            inlet_specific_h_j_kg: self.inlet_specific_h_j_kg,
            max_cells: self.max_cells,
            max_connections: self.max_connections,
        }
    }

    /// Continue sourced liquid heating with local viscosity feedback into the
    /// original full-stress momentum solve and its conservative h advection.
    ///
    /// The exact Cp/rho claim pins, pressure/composition coordinates, caloric
    /// law and fully-liquid support use the same admission as the frozen path.
    /// Density remains a validity-wide constant; viscosity may now vary along
    /// the original pinned temperature law. No new property table, name-based
    /// selection, phase mobility or second caloric law is introduced.
    ///
    /// Each actual thermal source epoch resolves its own h-derived T, so its
    /// returned `ResolvedLiquidState` drives the local momentum viscosity.
    /// The source's k remains available on the same samples; conduction is
    /// owned by the wrapped thermal observer. `max_solves` bounds steady
    /// momentum solves per advance. This is quasisteady one-step explicit
    /// thermal feedback, with the limits documented on `MaterialFlowFeedback`.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn feedback<O: DepositionObserver>(
        self,
        thermal: &DepositionSimulation<'_>,
        domain: &'a VoxelDomain,
        config: &'a SimpleConfig,
        inner: O,
        gate: &'a CancelGate,
        max_solves: usize,
    ) -> Result<
        MaterialFlowFeedback<
            'a,
            ResolvedLiquidState,
            impl Clone
            + FnMut(usize, &DepositedCell, f64) -> Result<ResolvedLiquidState, MaterialError>
            + use<'a, O>,
            O,
        >,
        LiquidEnthalpyObserverError<MaterialFlowError<O::Error>>,
    > {
        admit_law_mode(self, false, true).map_err(LiquidEnthalpyObserverError::Binding)?;
        for h in self.inlet_specific_h_j_kg.into_iter().flatten() {
            state_at_mode(self, h, false, true).map_err(LiquidEnthalpyObserverError::Binding)?;
        }
        let resolve = move |_: usize, cell: &DepositedCell, temperature: f64| {
            feedback_state(self, cell.specific_h_j_kg, temperature, true)
        };
        MaterialFlowFeedback::try_new(
            thermal,
            domain,
            self.native(),
            config,
            inner,
            gate,
            max_solves,
            resolve,
        )
        .map_err(LiquidEnthalpyObserverError::Observer)
    }

    /// Natural/mixed convection from this SAME pinned liquid/caloric source.
    ///
    /// The liquid's reference density supplies fixed mass and inertia, while
    /// its actual local density supplies the bounded Boussinesq gravity load.
    /// Its same-state mu and optional full-face k drive native momentum and
    /// thermal transport. Original Cp/density claims, non-temperature context,
    /// fully-liquid support and borrowed caloric identity remain mandatory.
    /// The caloric curve's density is source metadata in this fixed-volume
    /// approximation; it cannot change cell mass or create a volume source.
    /// See `MaterialFlowFeedback::try_new_boussinesq` for the energy boundary.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn feedback_boussinesq<O: DepositionObserver>(
        self,
        thermal: &DepositionSimulation<'_>,
        domain: &'a VoxelDomain,
        config: &'a SimpleConfig,
        inner: O,
        gate: &'a CancelGate,
        max_solves: usize,
        buoyancy: BoussinesqConfig<3>,
    ) -> Result<
        MaterialFlowFeedback<
            'a,
            BoussinesqMaterial<ResolvedLiquidState>,
            impl Clone
            + FnMut(
                usize,
                &DepositedCell,
                f64,
            ) -> Result<BoussinesqMaterial<ResolvedLiquidState>, MaterialError>
            + use<'a, O>,
            O,
        >,
        LiquidEnthalpyObserverError<MaterialFlowError<O::Error>>,
    > {
        admit_law_mode(self, false, false).map_err(LiquidEnthalpyObserverError::Binding)?;
        let reference = self
            .liquid
            .reference_state()
            .transport_state()
            .map_err(|error| {
                LiquidEnthalpyObserverError::Binding(LiquidEnthalpyError::Transport(error))
            })?;
        for h in self.inlet_specific_h_j_kg.into_iter().flatten() {
            let (_, liquid) = state_at_mode(self, h, false, false)
                .map_err(LiquidEnthalpyObserverError::Binding)?;
            BoussinesqMaterial::try_new(liquid, &reference, buoyancy).map_err(|error| {
                LiquidEnthalpyObserverError::Binding(LiquidEnthalpyError::Transport(error))
            })?;
        }
        let resolve = move |_: usize, cell: &DepositedCell, temperature: f64| {
            feedback_state(self, cell.specific_h_j_kg, temperature, false)
        };
        MaterialFlowFeedback::try_new_boussinesq(
            thermal,
            domain,
            self.native(),
            config,
            inner,
            gate,
            max_solves,
            reference,
            buoyancy,
            resolve,
        )
        .map_err(LiquidEnthalpyObserverError::Observer)
    }
}

fn feedback_state(
    input: LiquidFlowEnthalpyInput<'_>,
    h: f64,
    temperature: f64,
    fixed_density: bool,
) -> Result<ResolvedLiquidState, MaterialError> {
    let (phase, liquid) =
        state_at_mode(input, h, false, fixed_density).map_err(|error| MaterialError::State {
            what: error.to_string(),
        })?;
    if phase.temperature_k().to_bits() != temperature.to_bits() {
        return Err(MaterialError::State {
            what: "liquid feedback temperature differs from its original caloric owner".into(),
        });
    }
    Ok(liquid)
}

/// Original selected liquid state at an accepted thermal endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidFlowSample {
    pub voxel: usize,
    pub cell: CellId,
    pub phase: EquilibriumPhaseState,
    pub liquid: ResolvedLiquidState,
}

/// Binding/source refusal, retaining the original numerical-owner error.
#[derive(Debug)]
pub enum LiquidEnthalpyError {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Transport(MaterialError),
    Phase(PhaseStateError),
    Thermal(DepositionError),
}
impl fmt::Display for LiquidEnthalpyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for LiquidEnthalpyError {}
impl From<MaterialStatePointError> for LiquidEnthalpyError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<PhaseStateError> for LiquidEnthalpyError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<DepositionError> for LiquidEnthalpyError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// The user's physical-history error is preserved without type erasure.
#[derive(Debug)]
pub enum LiquidEnthalpyObserverError<E> {
    Binding(LiquidEnthalpyError),
    Observer(E),
}
impl<E: fmt::Debug> fmt::Display for LiquidEnthalpyObserverError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl<E: fmt::Debug> std::error::Error for LiquidEnthalpyObserverError<E> {}
impl<E> From<DepositionError> for LiquidEnthalpyObserverError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Binding(error.into())
    }
}

/// Source checks and the caller's history INSIDE the native Advected observer.
/// Clone deeply copies all mutable endpoint/history state. The original thermal
/// owner discards the entire candidate on any source or inner-history refusal.
#[derive(Clone)]
pub struct LiquidEnthalpyObserver<'a, O> {
    input: LiquidFlowEnthalpyInput<'a>,
    inner: O,
    samples: Vec<LiquidFlowSample>,
}

impl<'a, O> LiquidEnthalpyObserver<'a, O> {
    #[must_use]
    pub const fn inner(&self) -> &O {
        &self.inner
    }
    /// Most recent successful step/event samples in this advance. Empty before
    /// its first endpoint; construction separately admits the initial state.
    #[must_use]
    pub fn samples(&self) -> &[LiquidFlowSample] {
        &self.samples
    }
    #[must_use]
    pub const fn material_law(&self) -> &'a LiquidTemperatureLaw<'a> {
        self.input.liquid
    }
    #[must_use]
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.input.caloric
    }

    fn samples_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<LiquidFlowSample>, LiquidEnthalpyError> {
        check()?;
        if !thermal.uses_material_law(self.input.material_id, self.input.caloric.curve()) {
            return Err(LiquidEnthalpyError::Invalid(
                "the original sourced caloric law changed",
            ));
        }
        let mut samples = Vec::new();
        for (voxel, id) in self.input.cell_ids.iter().enumerate() {
            if voxel % 256 == 0 {
                check()?;
            }
            let Some(id) = id else {
                continue;
            };
            let cell = thermal.cells().get(id).ok_or(LiquidEnthalpyError::Invalid(
                "source-bound flowing cell is inactive",
            ))?;
            if cell.material_id != self.input.material_id {
                return Err(LiquidEnthalpyError::Invalid(
                    "flowing cell caloric assignment changed",
                ));
            }
            let (phase, liquid) = state_at(self.input, cell.specific_h_j_kg)?;
            samples.push(LiquidFlowSample {
                voxel,
                cell: *id,
                phase,
                liquid,
            });
        }
        check()?;
        Ok(samples)
    }
}

impl MaterialChtFlow<'_, ResolvedLiquidState> {
    /// Attach exact source support to the SAME native solved-flow advection.
    ///
    /// Density and viscosity must be validity-wide scalar sources; every caloric
    /// knot must be fully liquid at that density. Cp/density claims and every
    /// non-temperature coordinate/descriptor must agree with the retained liquid
    /// law. A matching material-card hash alone is insufficient. Temperature
    /// updates retain the original pins; leaving their support refuses instead
    /// of switching sources, keeping a stale coefficient, or clipping h.
    ///
    /// The native owner still checks actual geometry, reference mass, the
    /// registered borrowed caloric object, inlet h and local face continuity.
    /// Its roundoff mass-balance gate can require a tighter momentum solve.
    ///
    /// This is frozen laminar momentum with passive sensible-enthalpy transport.
    /// It does not infer phase mobility, pressure/viscous work, buoyancy,
    /// density-change volume flux, turbulent heat transport or surface motion.
    /// Thermal conduction and other physical histories remain with `inner`.
    pub fn advect_liquid_enthalpy<'s, O>(
        &'s self,
        thermal: &DepositionSimulation<'_>,
        input: LiquidFlowEnthalpyInput<'s>,
        inner: O,
        gate: &CancelGate,
    ) -> Result<
        MaterialFlowAdvection<'s, ResolvedLiquidState, LiquidEnthalpyObserver<'s, O>>,
        LiquidEnthalpyError,
    > {
        let mut check = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        if self.material() != input.liquid.reference_state() {
            return Err(LiquidEnthalpyError::Invalid(
                "momentum and thermal binding must retain the exact liquid reference state",
            ));
        }
        let flow = self.flow();
        if [&flow.eddy_viscosity_m2_s, &flow.heat_eddy_viscosity_m2_s]
            .iter()
            .any(|values| {
                values.len() != self.domain().cell_count()
                    || values.iter().any(|&value| value != 0.0)
            })
        {
            return Err(LiquidEnthalpyError::Invalid(
                "source-bound liquid enthalpy requires zero momentum and thermal eddy viscosity",
            ));
        }
        admit_law(input)?;
        for h in input.inlet_specific_h_j_kg.into_iter().flatten() {
            state_at(input, h)?;
        }
        let observer = LiquidEnthalpyObserver {
            input,
            inner,
            samples: Vec::new(),
        };
        // Native admission bounds the map before the source-sample allocation.
        let coupled = self.advect_enthalpy(thermal, input.native(), observer, gate)?;
        coupled
            .advection()
            .inner()
            .samples_at(thermal, &mut check)?;
        Ok(coupled)
    }
}

impl<O: DepositionObserver> DepositionObserver for LiquidEnthalpyObserver<'_, O> {
    type Error = LiquidEnthalpyObserverError<O::Error>;

    fn requires_step_interval(&self) -> bool { self.inner.requires_step_interval() }

    fn step_sources(&mut self, thermal_start: &DepositionSimulation<'_>,
        start_s: f64, dt_s: f64, end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.step_sources(thermal_start, start_s, dt_s, end_s, check)
            .map_err(LiquidEnthalpyObserverError::Observer)
    }

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner
            .begin()
            .map_err(LiquidEnthalpyObserverError::Observer)?;
        self.samples.clear();
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.samples_at(thermal, check)
            .map_err(LiquidEnthalpyObserverError::Binding)?;
        self.inner
            .sources(thermal, check)
            .map_err(LiquidEnthalpyObserverError::Observer)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner
            .conductances(thermal, check)
            .map_err(LiquidEnthalpyObserverError::Observer)
    }

    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner
            .heat_transfers(thermal, check)
            .map_err(LiquidEnthalpyObserverError::Observer)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.inner
            .limit_timestep_s(thermal, proposed_s, check)
            .map_err(LiquidEnthalpyObserverError::Observer)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        let samples = self
            .samples_at(thermal, &mut || Ok(()))
            .map_err(LiquidEnthalpyObserverError::Binding)?;
        self.inner
            .event(event, thermal)
            .map_err(LiquidEnthalpyObserverError::Observer)?;
        self.samples = samples;
        Ok(())
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        let samples = self
            .samples_at(thermal, check)
            .map_err(LiquidEnthalpyObserverError::Binding)?;
        self.inner
            .step(start_s, duration_s, old_temperatures_k, thermal, check)
            .map_err(LiquidEnthalpyObserverError::Observer)?;
        self.samples = samples;
        Ok(())
    }
}

fn admit_law(input: LiquidFlowEnthalpyInput<'_>) -> Result<(), LiquidEnthalpyError> {
    admit_law_mode(input, true, true)
}

fn admit_law_mode(
    input: LiquidFlowEnthalpyInput<'_>,
    frozen_viscosity: bool,
    fixed_density: bool,
) -> Result<(), LiquidEnthalpyError> {
    admit_law_with_rheology(input, frozen_viscosity, fixed_density, None)
}

fn admit_law_with_rheology(
    input: LiquidFlowEnthalpyInput<'_>,
    frozen_viscosity: bool,
    fixed_density: bool,
    shear_rate_axis: Option<&str>,
) -> Result<(), LiquidEnthalpyError> {
    let law = input.liquid;
    let reference = law.reference_state();
    let curve = input.caloric.curve();
    if reference.material().card_identity() != curve.material_card_identity() {
        return Err(LiquidEnthalpyError::Invalid(
            "liquid and caloric material cards differ",
        ));
    }
    if let Some(axis) = shear_rate_axis {
        // Reuse the shared explicit T/shear source admission used by VOF;
        // this marker cannot authorize any frozen-viscosity entry point.
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            reference.material().clone(),
            law.viscosity_basis(),
            law.temperature_axis(),
            axis,
        )?;
    } else {
        ResolvedLiquidSurfaceState::try_from_material(
            reference.material().clone(),
            law.viscosity_basis(),
            LiquidCapillarity::Neglected,
        )?;
    }
    // The momentum admission above checks rho/mu, not the retained thermal
    // roles. Their response axes also need explicit runtime owners.
    for name in [
        SPECIFIC_HEAT_CAPACITY_PROPERTY,
        THERMAL_CONDUCTIVITY_PROPERTY,
    ] {
        let property = reference
            .material()
            .property(name)
            .ok_or(LiquidEnthalpyError::Invalid(
                "missing liquid thermal source property",
            ))?;
        // Conductivity may share the explicitly evolved rate. Caloric storage
        // cannot: the original h(T) chart contains no shear-state energy.
        if shear_rate_axis.is_none() || name == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            require_frozen_coefficient_support(property)?;
        }
    }
    let viscosity = match law.viscosity_basis() {
        LiquidViscosityBasis::Dynamic => DYNAMIC_VISCOSITY_PROPERTY,
        LiquidViscosityBasis::Kinematic => KINEMATIC_VISCOSITY_PROPERTY,
    };
    for name in [DENSITY_PROPERTY, viscosity] {
        let property = reference
            .material()
            .property(name)
            .ok_or(LiquidEnthalpyError::Invalid(
                "missing liquid source property",
            ))?;
        if ((name == DENSITY_PROPERTY && fixed_density) || (name == viscosity && frozen_viscosity))
            && property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
        {
            return Err(LiquidEnthalpyError::Invalid(
                "fixed flow requires validity-wide scalar density and viscosity sources",
            ));
        }
    }
    if curve.knots().iter().any(|knot| {
        knot.liquid_mass_fraction != 1.0
            || (fixed_density && !same_positive(knot.bulk_density_kg_m3, reference.density_kg_m3()))
    }) {
        return Err(LiquidEnthalpyError::Invalid(
            "the complete caloric span must be fully liquid at the fixed source density",
        ));
    }
    let point = law.reference_point();
    for (index, receipt) in input.caloric.receipts().iter().enumerate() {
        // This immutable producer records Cp then rho at every source knot.
        let name = [SPECIFIC_HEAT_CAPACITY_PROPERTY, DENSITY_PROPERTY][index % 2];
        let property = reference
            .material()
            .property(name)
            .ok_or(LiquidEnthalpyError::Invalid(
                "missing liquid caloric source",
            ))?;
        if receipt.property != name || receipt.selected != property.answer().receipt.selected {
            return Err(LiquidEnthalpyError::Invalid(
                "caloric Cp/density source pins differ from the liquid law",
            ));
        }
        if receipt.axis_quantities != *point.axis_quantities()
            || !receipt
                .query_point
                .iter()
                .filter(|(axis, _)| axis.as_str() != law.temperature_axis())
                .map(|(axis, value)| (axis.as_str(), value.to_bits()))
                .eq(point
                    .axes()
                    .iter()
                    .filter(|(axis, _)| axis.as_str() != law.temperature_axis())
                    .map(|(axis, value)| (axis.as_str(), value.to_bits())))
        {
            return Err(LiquidEnthalpyError::Invalid(
                "caloric source context or temperature axis differs from the liquid law",
            ));
        }
    }
    Ok(())
}

fn state_at(
    input: LiquidFlowEnthalpyInput<'_>,
    h: f64,
) -> Result<(EquilibriumPhaseState, ResolvedLiquidState), LiquidEnthalpyError> {
    state_at_mode(input, h, true, true)
}

fn state_at_mode(
    input: LiquidFlowEnthalpyInput<'_>,
    h: f64,
    frozen_viscosity: bool,
    fixed_density: bool,
) -> Result<(EquilibriumPhaseState, ResolvedLiquidState), LiquidEnthalpyError> {
    let phase = input.caloric.curve().state_at_specific_enthalpy(h)?;
    let reference = input.liquid.reference_state();
    if phase.phase() != SolidLiquidPhase::Liquid
        || phase.liquid_mass_fraction() != 1.0
        || (fixed_density && !same_positive(phase.bulk_density_kg_m3(), reference.density_kg_m3()))
    {
        return Err(LiquidEnthalpyError::Invalid(
            "flowing enthalpy state must remain fully liquid at the fixed density",
        ));
    }
    let liquid = input
        .liquid
        .at_temperature(Temperature::new(phase.temperature_k()))?;
    if (fixed_density && liquid.density_kg_m3().to_bits() != reference.density_kg_m3().to_bits())
        || (frozen_viscosity
            && liquid.dynamic_viscosity_pa_s().to_bits()
                != reference.dynamic_viscosity_pa_s().to_bits())
    {
        return Err(LiquidEnthalpyError::Invalid(
            "the solved liquid momentum coefficients changed",
        ));
    }
    Ok((phase, liquid))
}
