//! Finite liquid temperature -> sourced viscosity -> joint Reynolds force/heat.
//!
//! The film has negligible thermal storage. A whole pressure island, a uniform
//! original passage, or each retained through-gap layer equilibrates to its
//! explicitly bound finite liquid reservoir. Its original caloric chart remains
//! the energy owner. The existing thermal network supplies its actual local
//! temperature; no island-average temperature or viscosity is substituted.
//! Flow through the pressure graph creates no second fluid mass inventory.
//! This is a lumped lubrication approximation, not continuum gap heat transport,
//! cavitation model or phase-changing flow solver. Heat destinations remain
//! explicit routes; a temperature binding never deposits additional heat.

mod layers;
mod sliding;
pub use layers::LiquidFilmLayerBinding;
pub use sliding::LiquidSlidingFilmLayerBinding;

use std::collections::BTreeSet;

use fs_contact::normal_patch::source::SourcePointContact;
use fs_exec::CancelGate;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

use super::super::super::squeeze::ResistiveFilm;
use super::super::super::{ImpactError, invalid};
use super::super::MaterialContact;
use super::feedback::MaterialFeedbackFrame;
use super::layered_film::LayerHeatRoute;
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

/// One original hydraulic passage at one finite liquid store's temperature.
/// Channel bindings cover a complete pressure island and retain original
/// channel indices. They cannot overlap a whole-island `LiquidFilmBinding`.
#[derive(Clone, Debug)]
pub struct LiquidFilmChannelBinding<'m> {
    /// Original pressure-island insertion index.
    pub film: usize,
    /// Original hydraulic channel index, including a currently closed passage.
    pub channel: usize,
    /// Stable key of the finite liquid reservoir setting this passage temperature.
    pub reservoir: String,
    /// Original liquid source, pinned claims and fixed nonthermal context.
    pub law: LiquidTemperatureLaw<'m>,
    /// Maximum `abs(mu_next - mu_old) / mu_old` in one accepted mechanical sample.
    /// Must lie in (0, 1]. This is a temporal-refinement gate, not an error bound.
    pub maximum_relative_viscosity_change: f64,
}

/// Original prepared mechanics and caloric stores with a live material binding.
///
/// At each accepted time, the retained sources and film operators match those
/// stores' temperature. The original mechanical solve holds viscosity during
/// its sample and computes pressure, tangent and dissipated work together. The
/// thermal proposal must admit the NEXT sourced liquid state before any owner
/// publishes. Stateless film coefficients and layer-work partitions are then
/// replaced together without changing stored energy.
///
/// This is first-order constitutive splitting; refine the mechanical period.
/// Constant density, fixed reservoir geometry, fully liquid phase and the
/// original thin-film pressure/gap bounds remain required. Temperature-dependent
/// material viscosity is allowed under load because this pressure-eliminated
/// resistance carries no stored elastic or compressible energy.
pub struct LiquidFilmImpactSystem<'m> {
    inner: HeatedImpactSystem<'m>,
    bindings: Vec<LiquidFilmBinding<'m>>,
    channel_bindings: Vec<LiquidFilmChannelBinding<'m>>,
    layer_bindings: Vec<LiquidFilmLayerBinding<'m>>,
    sliding_layer_bindings: Vec<LiquidSlidingFilmLayerBinding<'m>>,
    reservoir_indices: Vec<usize>,
    sources: Vec<ResolvedLiquidState>,
    groups: Vec<FilmBindingGroup>,
    sliding_groups: Vec<sliding::BindingGroup>,
}

// Normal uniform/channel/layer rows precede the sliding-layer source segment.
// Each group covers one complete island in its original mechanical owner.
struct FilmBindingGroup {
    film: usize,
    rows: std::ops::Range<usize>,
    kind: FilmBindingKind,
}

enum FilmBindingKind {
    Uniform,
    Channels,
    Layers(Vec<layers::ChannelBinding>),
}

struct CompiledFilm {
    operator: ResistiveFilm,
    layer_routes: Vec<(usize, LayerHeatRoute)>,
}

impl FilmBindingGroup {
    fn compile(
        &self,
        geometry: &ResistiveFilm,
        sources: &[ResolvedLiquidState],
    ) -> Result<CompiledFilm, ImpactError> {
        let operator = match &self.kind {
            FilmBindingKind::Uniform => {
                geometry.with_viscosity(sources[self.rows.start].dynamic_viscosity_pa_s())
            }
            FilmBindingKind::Channels => {
                let viscosities: Vec<_> = sources[self.rows.clone()]
                    .iter()
                    .map(ResolvedLiquidState::dynamic_viscosity_pa_s)
                    .collect();
                geometry.with_channel_viscosities(&viscosities)
            }
            FilmBindingKind::Layers(channels) => {
                return layers::compile(geometry, sources, channels);
            }
        }
        .map_err(|error| invalid(error.0))?;
        Ok(CompiledFilm {
            operator,
            layer_routes: Vec::new(),
        })
    }
}

impl<'m> LiquidFilmImpactSystem<'m> {
    /// Bind at the original accepted time without resetting motion or heat.
    /// The exact pressure graphs/limits are retained; only viscosity is compiled.
    /// A store's caloric chart and transport law must name the same material card
    /// and fixed physical context. The caloric chart owns h(T), including any
    /// approximation in its source integration; transport Cp is not reapplied as
    /// a second temperature update. The caller owns spatial lumping validity.
    pub fn new(
        inner: HeatedImpactSystem<'m>,
        bindings: Vec<LiquidFilmBinding<'m>>,
    ) -> Result<Self, ImpactError> {
        Self::new_with_channels(inner, bindings, Vec::new())
    }

    /// Bind uniform islands and/or complete channel-resolved islands on the
    /// SAME mechanical and enthalpy clock. Uniform rows retain their legacy
    /// meaning; channel rows resolve each passage at its own reservoir's actual
    /// accepted temperature. One store may control several passages.
    ///
    /// Channel rows are canonicalized by `(film, channel)` and must cover each
    /// selected island exactly once, including closed passages. Missing,
    /// duplicate, unknown and whole-island/channel-overlapping assignments
    /// refuse. Work is bounded by the existing admitted islands/channel counts.
    /// Each passage is uniform across its gap and length at that temperature;
    /// this does not infer a gapwise thermal profile or a mixing law.
    ///
    /// Existing whole-film or explicitly selected per-channel heat routes own
    /// work deposition. A viscosity binding neither chooses those destinations
    /// nor deposits another copy. All source and density/phase gates of `new`
    /// apply independently to every row before any coefficients are published.
    pub fn new_with_channels(
        inner: HeatedImpactSystem<'m>,
        bindings: Vec<LiquidFilmBinding<'m>>,
        channel_bindings: Vec<LiquidFilmChannelBinding<'m>>,
    ) -> Result<Self, ImpactError> {
        Self::new_with_layers(inner, bindings, channel_bindings, Vec::new())
    }

    /// Bind complete source-backed viscosity layers as well as uniform/channel
    /// islands on the original mechanical and thermal clock.
    ///
    /// Every channel of a layer-selected island must already have a retained
    /// layer heat route, and every original layer requires one unique binding.
    /// Its temperature reservoir is its explicit heat destination. Normalized
    /// interfaces, reservoir inventories and source context stay fixed; the
    /// original source laws update viscosity at the accepted local temperatures.
    ///
    /// The hydraulic coefficient and physical layer-work fractions are rebuilt
    /// together, then published with source states only after the existing
    /// mechanical/enthalpy proposal accepts. Every individual layer obeys its
    /// step-change allowance, even if changes cancel in the hydraulic equivalent.
    /// Completed work uses the old profile; no heat is repriced at the endpoint.
    /// This is first-order constitutive splitting with prescribed interfaces,
    /// not layer advection, mixing or a moving-interface thermal solver.
    pub fn new_with_layers(
        inner: HeatedImpactSystem<'m>,
        bindings: Vec<LiquidFilmBinding<'m>>,
        channel_bindings: Vec<LiquidFilmChannelBinding<'m>>,
        layer_bindings: Vec<LiquidFilmLayerBinding<'m>>,
    ) -> Result<Self, ImpactError> {
        Self::new_with_sliding(
            inner,
            bindings,
            channel_bindings,
            layer_bindings,
            Vec::new(),
        )
    }

    /// Bind source-temperature feedback for sliding films on the SAME thermal
    /// proposal as any existing normal-pressure film bindings.
    ///
    /// Sliding indices refer to `sliding_film_count`, independently of the
    /// normal-pressure indices in the other three binding lists. Each selected
    /// sliding island needs every original channel and layer exactly once. A
    /// uniform channel has one implicit layer; channels may share a liquid
    /// reservoir or use different source laws and accepted temperatures.
    /// Each layer's temperature reservoir must be its existing heat destination.
    ///
    /// The original native film recompiles pressure mobility, Couette
    /// entrainment, shear drag and complete layer-work partition together.
    /// Graph geometry, wall rows, layer interfaces and Barus pressure law stay
    /// fixed. Completed mechanical work uses the old coefficients; all updated
    /// coefficients and source states publish with the accepted thermal state.
    /// This is the existing first-order constitutive split, with no extra heat,
    /// material inventory or temperature inferred from pressure or slip speed.
    pub fn new_with_sliding(
        mut inner: HeatedImpactSystem<'m>,
        mut bindings: Vec<LiquidFilmBinding<'m>>,
        mut channel_bindings: Vec<LiquidFilmChannelBinding<'m>>,
        mut layer_bindings: Vec<LiquidFilmLayerBinding<'m>>,
        mut sliding_layer_bindings: Vec<LiquidSlidingFilmLayerBinding<'m>>,
    ) -> Result<Self, ImpactError> {
        let maximum_bindings: usize = (0..inner.mechanics().squeeze_film_count())
            .map(|index| {
                inner
                    .mechanics()
                    .squeeze_film(index)
                    .expect("existing film")
                    .channel_count()
            })
            .sum();
        let ordinary_count = bindings.len().saturating_add(channel_bindings.len());
        let maximum_layers: usize = inner
            .film_routes
            .iter()
            .filter(|route| route.sliding_layers.is_none())
            .filter_map(|route| route.layers.as_ref())
            .map(|route| route.profile.layers().len())
            .sum();
        let maximum_sliding_layers: usize = (0..inner.mechanics().sliding_film_count())
            .map(|index| {
                inner
                    .mechanics()
                    .sliding_film(index)
                    .expect("existing sliding film")
                    .layer_count()
            })
            .sum();
        let normal_count = ordinary_count.saturating_add(layer_bindings.len());
        let count = normal_count.saturating_add(sliding_layer_bindings.len());
        if count == 0
            || ordinary_count > maximum_bindings
            || layer_bindings.len() > maximum_layers
            || sliding_layer_bindings.len() > maximum_sliding_layers
        {
            return Err(invalid(
                "liquid film feedback needs bounded bindings to existing pressure channels",
            ));
        }
        bindings.sort_by_key(|binding| binding.film);
        channel_bindings.sort_by_key(|binding| (binding.film, binding.channel));
        layer_bindings.sort_by_key(|binding| (binding.film, binding.channel, binding.layer));
        sliding_layer_bindings
            .sort_by_key(|binding| (binding.film, binding.channel, binding.layer));
        let mut groups = Vec::new();
        for (i, binding) in bindings.iter().enumerate() {
            if i > 0 && bindings[i - 1].film == binding.film {
                return Err(invalid("duplicate uniform film temperature binding"));
            }
            inner
                .mechanics()
                .squeeze_film(binding.film)
                .ok_or_else(|| invalid("liquid feedback names an absent pressure island"))?;
            groups.push(FilmBindingGroup {
                film: binding.film,
                rows: i..i + 1,
                kind: FilmBindingKind::Uniform,
            });
        }
        let mut start = 0;
        while start < channel_bindings.len() {
            let film_index = channel_bindings[start].film;
            let end = start
                + channel_bindings[start..].partition_point(|binding| binding.film == film_index);
            let film = inner
                .mechanics()
                .squeeze_film(film_index)
                .ok_or_else(|| invalid("liquid feedback names an absent pressure island"))?;
            if groups.iter().any(|group| group.film == film_index)
                || end - start != film.channel_count()
                || channel_bindings[start..end]
                    .iter()
                    .enumerate()
                    .any(|(channel, binding)| binding.channel != channel)
            {
                return Err(invalid(
                    "film temperature needs complete unique channel coverage without a uniform binding",
                ));
            }
            groups.push(FilmBindingGroup {
                film: film_index,
                rows: bindings.len() + start..bindings.len() + end,
                kind: FilmBindingKind::Channels,
            });
            start = end;
        }
        if groups.iter().any(|group| {
            inner.film_routes.iter().any(|route| {
                route.sliding_layers.is_none() && route.film == group.film && route.layers.is_some()
            })
        }) {
            return Err(invalid(
                "uniform/channel liquid feedback cannot replace a retained through-gap layer profile",
            ));
        }
        layers::append_groups(&inner, &layer_bindings, ordinary_count, &mut groups)?;
        groups.sort_by_key(|group| group.film);
        let sliding_groups = sliding::groups(&inner, &sliding_layer_bindings, normal_count)?;
        let mut reservoir_indices = Vec::with_capacity(count);
        let mut sources = Vec::with_capacity(count);
        for (reservoir, law, maximum_change) in source_bindings(
            &bindings,
            &channel_bindings,
            &layer_bindings,
            &sliding_layer_bindings,
        ) {
            if !maximum_change.is_finite() || maximum_change <= 0.0 || maximum_change > 1.0 {
                return Err(invalid("invalid constitutive step-change allowance"));
            }
            let index = inner
                .reservoirs()
                .iter()
                .position(|store| store.key() == reservoir)
                .ok_or_else(|| invalid("liquid feedback names an absent caloric reservoir"))?;
            let store = &inner.reservoirs()[index];
            let source = resolve(law, store, store.state())?;
            reservoir_indices.push(index);
            sources.push(source);
        }
        let mut updates = Vec::with_capacity(groups.len());
        let mut layer_updates = Vec::new();
        for group in &groups {
            let geometry = inner
                .mechanics()
                .squeeze_film(group.film)
                .expect("bound film");
            let compiled = group.compile(geometry, &sources)?;
            updates.push((group.film, compiled.operator));
            layer_updates.extend(compiled.layer_routes);
        }
        let mut sliding_updates = Vec::with_capacity(sliding_groups.len());
        for group in &sliding_groups {
            let film = inner
                .mechanics()
                .sliding_film(group.film)
                .expect("bound sliding film");
            sliding_updates.push((group.film, group.compile(film, &sources)?));
        }
        // No fallible operation follows this coefficient publication.
        inner.mechanics.swap_squeeze_films(&mut updates);
        inner.mechanics.swap_sliding_films(&mut sliding_updates);
        layers::publish(&mut inner, &mut layer_updates);
        Ok(Self {
            inner,
            bindings,
            channel_bindings,
            layer_bindings,
            sliding_layer_bindings,
            reservoir_indices,
            sources,
            groups,
            sliding_groups,
        })
    }

    /// Read-only accepted mechanics and finite thermal state.
    pub const fn system(&self) -> &HeatedImpactSystem<'m> {
        &self.inner
    }
    /// Whole-island prescriptions in pressure-island order. Channel-resolved
    /// islands are inspected separately through `channel_bindings`.
    pub fn bindings(&self) -> &[LiquidFilmBinding<'m>] {
        &self.bindings
    }
    /// Channel prescriptions in `(film, channel)` order, retaining each original
    /// passage index and its own reservoir/source assignment.
    pub fn channel_bindings(&self) -> &[LiquidFilmChannelBinding<'m>] {
        &self.channel_bindings
    }
    /// Layer prescriptions in original `(film, channel, layer)` order.
    pub fn layer_bindings(&self) -> &[LiquidFilmLayerBinding<'m>] {
        &self.layer_bindings
    }
    /// Sliding-layer prescriptions in original `(film, channel, layer)` order.
    /// Uniform sliding passages each retain one implicit layer numbered zero.
    pub fn sliding_layer_bindings(&self) -> &[LiquidSlidingFilmLayerBinding<'m>] {
        &self.sliding_layer_bindings
    }
    /// Actual uniform-island material states at the accepted temperature.
    /// Each entry corresponds to the same entry of `bindings`; channel source
    /// states are available through `channel_sources`.
    pub fn sources(&self) -> &[ResolvedLiquidState] {
        &self.sources[..self.bindings.len()]
    }
    /// Actual per-passage material states at their respective accepted local
    /// temperatures, in the same order as `channel_bindings`.
    pub fn channel_sources(&self) -> &[ResolvedLiquidState] {
        &self.sources[self.bindings.len()..self.bindings.len() + self.channel_bindings.len()]
    }
    /// Actual per-layer liquid states at the accepted temperatures, in the same
    /// canonical order as `layer_bindings`. Active profiles remain available
    /// through `system().layer_heat_profile(film, channel)`.
    pub fn layer_sources(&self) -> &[ResolvedLiquidState] {
        let start = self.bindings.len() + self.channel_bindings.len();
        &self.sources[start..start + self.layer_bindings.len()]
    }
    /// Actual sliding-layer sources at accepted local temperatures, in the
    /// same order as `sliding_layer_bindings`. Native layer profiles and wall
    /// mappings remain available through `system().mechanics().sliding_film`.
    pub fn sliding_layer_sources(&self) -> &[ResolvedLiquidState] {
        &self.sources
            [self.bindings.len() + self.channel_bindings.len() + self.layer_bindings.len()..]
    }
    /// Select the existing analytic/finite-difference Newton image.
    pub fn set_analytic_newton(&mut self, enabled: bool) {
        self.inner.set_analytic_newton(enabled);
    }
    /// Retain the numerical owner's original iteration ceiling semantics.
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), ImpactError> {
        self.inner.set_iteration_limit(limit)
    }
    /// Release accepted mechanics/heat and the uniform-island prescriptions.
    /// This extraction omits channel, layer and sliding prescriptions. Use
    /// `into_all_parts` for uniform/channel sources or `into_layered_parts`
    /// when normal-pressure layer prescriptions must also be retained, or
    /// `into_sliding_parts` for every prescription.
    pub fn into_parts(self) -> (HeatedImpactSystem<'m>, Vec<LiquidFilmBinding<'m>>) {
        (self.inner, self.bindings)
    }
    /// Release accepted mechanics/heat and the uniform/channel prescriptions.
    /// This extraction omits layer and sliding bindings; `into_sliding_parts`
    /// retains every prescription.
    /// Rebinding uses accepted local temperatures without resetting the clock.
    pub fn into_all_parts(
        self,
    ) -> (
        HeatedImpactSystem<'m>,
        Vec<LiquidFilmBinding<'m>>,
        Vec<LiquidFilmChannelBinding<'m>>,
    ) {
        (self.inner, self.bindings, self.channel_bindings)
    }

    /// Release accepted mechanics/heat and the normal-pressure prescriptions.
    /// Sliding prescriptions are retained by `into_sliding_parts`.
    /// Rebinding retains motion, caloric history, layer interfaces and time.
    pub fn into_layered_parts(
        self,
    ) -> (
        HeatedImpactSystem<'m>,
        Vec<LiquidFilmBinding<'m>>,
        Vec<LiquidFilmChannelBinding<'m>>,
        Vec<LiquidFilmLayerBinding<'m>>,
    ) {
        (
            self.inner,
            self.bindings,
            self.channel_bindings,
            self.layer_bindings,
        )
    }

    /// Release accepted mechanics/heat and all normal/sliding source bindings.
    /// Rebinding keeps accepted thermal and mechanical histories and time.
    pub fn into_sliding_parts(
        self,
    ) -> (
        HeatedImpactSystem<'m>,
        Vec<LiquidFilmBinding<'m>>,
        Vec<LiquidFilmChannelBinding<'m>>,
        Vec<LiquidFilmLayerBinding<'m>>,
        Vec<LiquidSlidingFilmLayerBinding<'m>>,
    ) {
        (
            self.inner,
            self.bindings,
            self.channel_bindings,
            self.layer_bindings,
            self.sliding_layer_bindings,
        )
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
            channel_bindings,
            layer_bindings,
            sliding_layer_bindings,
            reservoir_indices,
            sources,
            groups,
            sliding_groups,
        } = self;
        let stores = inner.reservoirs().to_vec();
        let geometry: Vec<_> = groups
            .iter()
            .map(|binding| {
                inner
                    .mechanics()
                    .squeeze_film(binding.film)
                    .expect("bound film")
                    .clone()
            })
            .collect();
        let sliding_geometry: Vec<_> = sliding_groups
            .iter()
            .map(|group| {
                inner
                    .mechanics()
                    .sliding_film(group.film)
                    .expect("bound sliding film")
                    .clone()
            })
            .collect();
        let mut next_sources = Vec::with_capacity(sources.len());
        let mut updates = Vec::with_capacity(groups.len());
        let mut layer_updates = Vec::new();
        let mut sliding_updates = Vec::with_capacity(sliding_groups.len());
        let frame = inner.step_with_material_feedback(forces, gate, contact, |trial| {
            let proposed = transport(trial)?;
            if proposed.transfers.len() != stores.len() {
                return Err(invalid(
                    "liquid feedback requires the complete original thermal proposal",
                ));
            }
            for (i, (_, law, maximum_change)) in source_bindings(
                bindings,
                channel_bindings,
                layer_bindings,
                sliding_layer_bindings,
            )
            .enumerate()
            {
                if gate.is_requested() {
                    return Err(ImpactError::Cancelled);
                }
                let index = reservoir_indices[i];
                let source = resolve(law, &stores[index], proposed.transfers[index].state)?;
                let old_mu = sources[i].dynamic_viscosity_pa_s();
                let new_mu = source.dynamic_viscosity_pa_s();
                let relative = (new_mu - old_mu).abs() / old_mu;
                if !relative.is_finite() || relative > maximum_change {
                    return Err(invalid(
                        "liquid viscosity changed too far in one mechanical sample; refine dt",
                    ));
                }
                next_sources.push(source);
            }
            for (group, geometry) in groups.iter().zip(&geometry) {
                if gate.is_requested() {
                    return Err(ImpactError::Cancelled);
                }
                let compiled = group.compile(geometry, &next_sources)?;
                updates.push((group.film, compiled.operator));
                layer_updates.extend(compiled.layer_routes);
            }
            for (group, film) in sliding_groups.iter().zip(&sliding_geometry) {
                if gate.is_requested() {
                    return Err(ImpactError::Cancelled);
                }
                sliding_updates.push((group.film, group.compile(film, &next_sources)?));
            }
            Ok(proposed)
        })?;
        // All thermal/source/contact/mechanical gates have passed. These swaps
        // cannot fail and carry no change of stored energy or time.
        inner.mechanics.swap_squeeze_films(&mut updates);
        inner.mechanics.swap_sliding_films(&mut sliding_updates);
        layers::publish(inner, &mut layer_updates);
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

fn source_bindings<'s, 'm>(
    uniform: &'s [LiquidFilmBinding<'m>],
    channels: &'s [LiquidFilmChannelBinding<'m>],
    layers: &'s [LiquidFilmLayerBinding<'m>],
    sliding_layers: &'s [LiquidSlidingFilmLayerBinding<'m>],
) -> impl Iterator<Item = (&'s str, &'s LiquidTemperatureLaw<'m>, f64)> {
    uniform
        .iter()
        .map(|binding| {
            (
                binding.reservoir.as_str(),
                &binding.law,
                binding.maximum_relative_viscosity_change,
            )
        })
        .chain(channels.iter().map(|binding| {
            (
                binding.reservoir.as_str(),
                &binding.law,
                binding.maximum_relative_viscosity_change,
            )
        }))
        .chain(layers.iter().map(|binding| {
            (
                binding.reservoir.as_str(),
                &binding.law,
                binding.maximum_relative_viscosity_change,
            )
        }))
        .chain(sliding_layers.iter().map(|binding| {
            (
                binding.reservoir.as_str(),
                &binding.law,
                binding.maximum_relative_viscosity_change,
            )
        }))
}

fn resolve(
    law: &LiquidTemperatureLaw<'_>,
    store: &CaloricReservoir<'_>,
    state: EquilibriumPhaseState,
) -> Result<ResolvedLiquidState, ImpactError> {
    let card = law.card().content_hash();
    if store.curve().material_card_identity() != card
        || state.material_card_identity() != card
        || state.phase_curve_identity() != store.curve().identity()
        || state.liquid_mass_fraction() != 1.0
    {
        return Err(invalid(
            "film temperature requires its original fully liquid caloric material",
        ));
    }
    let source = law
        .at_temperature(Temperature::new(state.temperature_k()))
        .map_err(|error| ImpactError::Owner(error.to_string()))?;
    let density = source.density_kg_m3();
    let mass = density * store.volume_m3();
    if density != law.reference_state().density_kg_m3()
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
