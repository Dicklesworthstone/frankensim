//! Finite liquid temperature -> sourced viscosity -> joint Reynolds force/heat.
//!
//! The film has negligible thermal storage and equilibrates to an explicitly
//! bound finite, well-mixed liquid reservoir. That reservoir's original caloric
//! chart remains the energy owner. Flow through the pressure graph does not
//! create another fluid mass inventory. This is a lubrication approximation,
//! not a spatial thermal film, cavitation model or phase-changing flow solver.

use std::collections::BTreeSet;

use fs_contact::normal_patch::source::SourcePointContact;
use fs_exec::CancelGate;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

use super::super::super::{ImpactError, invalid};
use super::super::MaterialContact;
use super::feedback::MaterialFeedbackFrame;
use super::network::{self, NetworkImpactFrame, ReservoirCell};
use super::transport::{CoupledThermalStep, ImpactThermalTrial};
use super::{CaloricReservoir, HeatedImpactSystem};

/// One existing pressure island equilibrated to one existing finite liquid store.
/// Heat destinations remain the separately authored `SqueezeFilmHeatRoute`s:
/// there is no automatic fluid/wall partition or additional heat deposition.
#[derive(Clone, Debug)]
pub struct LiquidFilmBinding<'m> {
    /// Original pressure-island insertion index.
    pub film: usize,
    /// Stable key of the liquid reservoir setting this film's temperature.
    pub reservoir: String,
    /// Shared liquid law, retaining the card, claims and nonthermal context.
    pub law: LiquidTemperatureLaw<'m>,
    /// Maximum `abs(mu_next - mu_old) / mu_old` in one mechanical sample.
    /// Must lie in (0, 1]. A refusal asks for temporal refinement; this bound
    /// is not an a posteriori error estimate for the lagged constitutive law.
    pub maximum_relative_viscosity_change: f64,
}

/// Original prepared mechanics and caloric stores with a live material binding.
///
/// At each accepted time, the retained sources and film operators match those
/// stores' temperature. The original mechanical solve holds viscosity during
/// its sample and computes pressure, tangent and dissipated work together. The
/// thermal proposal must admit the NEXT sourced liquid state before any owner
/// publishes. Only then are the stateless film coefficients replaced.
///
/// This is first-order constitutive splitting; refine the mechanical period.
/// Constant density, fixed reservoir geometry, fully liquid phase and the
/// original thin-film pressure/gap bounds remain required. Temperature-dependent
/// material viscosity is allowed under load because this pressure-eliminated
/// resistance carries no stored elastic or compressible energy.
pub struct LiquidFilmImpactSystem<'m> {
    inner: HeatedImpactSystem<'m>,
    bindings: Vec<LiquidFilmBinding<'m>>,
    reservoir_indices: Vec<usize>,
    sources: Vec<ResolvedLiquidState>,
}

impl<'m> LiquidFilmImpactSystem<'m> {
    /// Bind at the original accepted time without resetting motion or heat.
    /// The exact pressure graphs/limits are retained; only viscosity is compiled.
    /// A store's caloric chart and transport law must name the same material card
    /// and fixed physical context. The caloric chart owns h(T), including any
    /// approximation in its source integration; transport Cp is not reapplied as
    /// a second temperature update. The caller owns spatial lumping validity.
    pub fn new(
        mut inner: HeatedImpactSystem<'m>,
        mut bindings: Vec<LiquidFilmBinding<'m>>,
    ) -> Result<Self, ImpactError> {
        if bindings.is_empty() || bindings.len() > 4 {
            return Err(invalid(
                "liquid film feedback needs one to four existing pressure islands",
            ));
        }
        bindings.sort_by_key(|binding| binding.film);
        let mut reservoir_indices = Vec::with_capacity(bindings.len());
        let mut sources = Vec::with_capacity(bindings.len());
        let mut updates = Vec::with_capacity(bindings.len());
        for (i, binding) in bindings.iter().enumerate() {
            if (i > 0 && bindings[i - 1].film == binding.film)
                || !binding.maximum_relative_viscosity_change.is_finite()
                || binding.maximum_relative_viscosity_change <= 0.0
                || binding.maximum_relative_viscosity_change > 1.0
            {
                return Err(invalid(
                    "duplicate film or invalid constitutive step-change allowance",
                ));
            }
            let film = inner
                .mechanics()
                .squeeze_film(binding.film)
                .ok_or_else(|| invalid("liquid feedback names an absent pressure island"))?;
            let index = inner
                .reservoirs()
                .iter()
                .position(|store| store.key() == binding.reservoir)
                .ok_or_else(|| invalid("liquid feedback names an absent caloric reservoir"))?;
            let store = &inner.reservoirs()[index];
            let source = resolve(binding, store, store.state())?;
            updates.push((
                binding.film,
                film.with_viscosity(source.dynamic_viscosity_pa_s())
                    .map_err(|error| invalid(error.0))?,
            ));
            reservoir_indices.push(index);
            sources.push(source);
        }
        // No fallible operation follows this coefficient publication.
        inner.mechanics.swap_squeeze_films(&mut updates);
        Ok(Self {
            inner,
            bindings,
            reservoir_indices,
            sources,
        })
    }

    /// Read-only accepted mechanics and finite thermal state.
    pub const fn system(&self) -> &HeatedImpactSystem<'m> {
        &self.inner
    }
    /// Retained data prescriptions in pressure-island order.
    pub fn bindings(&self) -> &[LiquidFilmBinding<'m>] {
        &self.bindings
    }
    /// Actual material states driving the operators at the accepted temperature.
    /// Each entry corresponds to the same entry of `bindings`.
    pub fn sources(&self) -> &[ResolvedLiquidState] {
        &self.sources
    }
    /// Select the existing analytic/finite-difference Newton image.
    pub fn set_analytic_newton(&mut self, enabled: bool) {
        self.inner.set_analytic_newton(enabled);
    }
    /// Retain the numerical owner's original iteration ceiling semantics.
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), ImpactError> {
        self.inner.set_iteration_limit(limit)
    }
    /// Release the accepted state and its constitutive data together.
    pub fn into_parts(self) -> (HeatedImpactSystem<'m>, Vec<LiquidFilmBinding<'m>>) {
        (self.inner, self.bindings)
    }

    /// Advance with adiabatic finite stores and frozen contact sources.
    /// Only the original selected mechanical work enters caloric storage.
    pub fn step(
        &mut self,
        forces: &[f64],
        gate: &CancelGate,
    ) -> Result<MaterialFeedbackFrame<()>, ImpactError> {
        let tolerance = self.inner.heat_tolerance_j;
        self.step_with_coupled_transport(
            forces,
            gate,
            |_, _| Ok(None),
            |trial| {
                let transfers = trial
                    .reservoirs()
                    .iter()
                    .map(|input| {
                        let input = input.thermal;
                        let state = input
                            .curve
                            .advance_specific_energy(
                                input.initial,
                                input.internal_heat_j / input.mass_kg,
                            )
                            .map_err(|error| ImpactError::Owner(error.to_string()))?;
                        Ok(UniformEnthalpyStep {
                            state,
                            external_heat_j: 0.0,
                            energy_residual_tolerance_j: tolerance,
                            report: (),
                        })
                    })
                    .collect::<Result<Vec<_>, ImpactError>>()?;
                Ok(CoupledThermalStep {
                    transfers,
                    boundary_heat_j: 0.0,
                })
            },
        )
    }

    /// Use the original joint finite-store transport and unloaded-contact refresh.
    ///
    /// Every proposed liquid temperature must resolve the original pinned claims
    /// and retain constant density/fully liquid state. Failed source support,
    /// excessive viscosity variation, heat/force solve, or late cancellation
    /// preserves motion, heat, contacts, film coefficients and source receipts.
    /// Callbacks must propose isolated state; external side effects and panics
    /// cannot be rolled back. No film pressure or loss is recomputed at the new
    /// temperature and charged retroactively to the completed mechanical sample.
    pub fn step_with_coupled_transport<R>(
        &mut self,
        forces: &[f64],
        gate: &CancelGate,
        contact: impl FnMut(
            &MaterialContact,
            &[CaloricReservoir<'m>],
        ) -> Result<Option<SourcePointContact>, ImpactError>,
        mut transport: impl FnMut(ImpactThermalTrial<'_>) -> Result<CoupledThermalStep<R>, ImpactError>,
    ) -> Result<MaterialFeedbackFrame<R>, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        let Self {
            inner,
            bindings,
            reservoir_indices,
            sources,
        } = self;
        let stores = inner.reservoirs().to_vec();
        let geometry: Vec<_> = bindings
            .iter()
            .map(|binding| {
                inner
                    .mechanics()
                    .squeeze_film(binding.film)
                    .expect("bound film")
                    .clone()
            })
            .collect();
        let mut next_sources = Vec::with_capacity(bindings.len());
        let mut updates = Vec::with_capacity(bindings.len());
        let frame = inner.step_with_material_feedback(forces, gate, contact, |trial| {
            let proposed = transport(trial)?;
            if proposed.transfers.len() != stores.len() {
                return Err(invalid(
                    "liquid feedback requires the complete original thermal proposal",
                ));
            }
            for (i, binding) in bindings.iter().enumerate() {
                if gate.is_requested() {
                    return Err(ImpactError::Cancelled);
                }
                let index = reservoir_indices[i];
                let source = resolve(binding, &stores[index], proposed.transfers[index].state)?;
                let old_mu = sources[i].dynamic_viscosity_pa_s();
                let new_mu = source.dynamic_viscosity_pa_s();
                let relative = (new_mu - old_mu).abs() / old_mu;
                if !relative.is_finite() || relative > binding.maximum_relative_viscosity_change {
                    return Err(invalid(
                        "liquid viscosity changed too far in one mechanical sample; refine dt",
                    ));
                }
                updates.push((
                    binding.film,
                    geometry[i]
                        .with_viscosity(new_mu)
                        .map_err(|error| invalid(error.0))?,
                ));
                next_sources.push(source);
            }
            Ok(proposed)
        })?;
        // All thermal/source/contact/mechanical gates have passed. These swaps
        // cannot fail and carry no change of stored energy or time.
        inner.mechanics.swap_squeeze_films(&mut updates);
        *sources = next_sources;
        Ok(frame)
    }

    /// Couple live film material to the SAME multi-cell conduction/source owner.
    /// The network's actual substeps receive selected mechanical work once.
    /// The observer, network, motion, film coefficients and source states publish
    /// together; a rejected endpoint material state publishes none of them.
    pub fn step_with_thermal_network<'n, O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'n>,
        observer: &mut O,
        bindings: &[ReservoirCell],
        config: StepConfig,
        forces: &[f64],
        gate: &CancelGate,
        contact: impl FnMut(
            &MaterialContact,
            &[CaloricReservoir<'m>],
        ) -> Result<Option<SourcePointContact>, ImpactError>,
    ) -> Result<NetworkImpactFrame, ImpactError>
    where
        O::Error: core::fmt::Display,
    {
        if !config.energy_tolerance_j.is_finite()
            || config.energy_tolerance_j < 0.0
            || config.energy_tolerance_j > self.inner.heat_tolerance_j
        {
            return Err(invalid(
                "network transport cannot enlarge the mechanical heat allowance",
            ));
        }
        let mut staged = None;
        let impact = self.step_with_coupled_transport(forces, gate, contact, |trial| {
            let (next, next_observer, audit, proposals) = network::propose(
                thermal,
                observer,
                bindings,
                config,
                &trial,
                gate,
                BTreeSet::new(),
            )?;
            staged = Some((next, next_observer, audit));
            Ok(proposals)
        })?;
        let (next, next_observer, audit) = staged.expect("accepted network proposal");
        let frame = network::finish(impact, audit);
        *thermal = next;
        *observer = next_observer;
        Ok(frame)
    }
}

fn resolve(
    binding: &LiquidFilmBinding<'_>,
    store: &CaloricReservoir<'_>,
    state: EquilibriumPhaseState,
) -> Result<ResolvedLiquidState, ImpactError> {
    let card = binding.law.card().content_hash();
    if store.curve().material_card_identity() != card
        || state.material_card_identity() != card
        || state.phase_curve_identity() != store.curve().identity()
        || state.liquid_mass_fraction() != 1.0
    {
        return Err(invalid(
            "film temperature requires its original fully liquid caloric material",
        ));
    }
    let source = binding
        .law
        .at_temperature(Temperature::new(state.temperature_k()))
        .map_err(|error| ImpactError::Owner(error.to_string()))?;
    let density = source.density_kg_m3();
    let mass = density * store.volume_m3();
    if density != binding.law.reference_state().density_kg_m3()
        || density != state.bulk_density_kg_m3()
        || !mass.is_finite()
        || (mass - store.mass_kg()).abs() > 64.0 * f64::EPSILON * store.mass_kg()
    {
        return Err(invalid(
            "liquid film feedback requires unchanged consistent density, mass and volume",
        ));
    }
    Ok(source)
}
