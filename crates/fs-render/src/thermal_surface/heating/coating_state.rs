//! Independent native film/substrate states drive one coherent optical solve.
//!
//! Each constitutive law reads its own receiving cell. The original coherent
//! owner then recomputes the entire stack, since interference couples all films.
//! Receiver/law pairs resolve once per interval; identical stack/binding/angle
//! packets are shared across sources. Only changed optical constants invalidate
//! packets, and dark intervals require domain validation without a Fresnel solve.
//!
//! Geometry is fixed. Whole-film thermal cells and prompt thermalization retain
//! the assumptions of `try_new_coated`. Coefficients use the accepted interval
//! start; temporal refinement is necessary for nonlinear optical/thermal feedback.

use super::coating::{AbsorptionReceiver, PreparedCoatingSource};
use super::feedback::SurfaceConductorLaw;
use super::*;
use crate::conductor::{
    CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack, CoherentDielectricCoating,
    ConductorOptics, MAX_COHERENT_STACK_LAYERS,
};
use std::collections::BTreeSet;

/// A finite film keeps its fixed optical law or follows its own native state.
#[derive(Clone, Debug)]
pub enum CoatingOpticsModel {
    Fixed(CoherentAbsorbingCoating),
    /// Retain the original Cauchy/Sellmeier law of an unchanging clear film.
    Lossless(CoherentDielectricCoating),
    Law(Arc<SurfaceCoatingLaw>),
}

/// Opaque substrate constants are fixed or follow its own caloric state.
/// The explicit stack replaces any coating attached to this substrate owner.
#[derive(Clone, Debug)]
pub enum SubstrateOpticsModel {
    Fixed(ConductorOptics),
    Law(Arc<SurfaceConductorLaw>),
}

/// One authored physical film and its original native heat destination.
#[derive(Clone, Debug)]
pub struct StatefulCoatingLayer {
    pub receiver: SurfaceHeatingReceiver,
    pub optics: CoatingOpticsModel,
}

/// One illumination of a coherent stack with independently evolving layers.
/// A law's native cell is always its heat receiver; material/enthalpy identities
/// and temperature domains remain binding in darkness and at accepted endpoints.
#[derive(Clone, Debug)]
pub struct StatefulCoatedSurfaceIlluminationSource<const N: usize> {
    pub source_id: u64,
    pub films: Vec<StatefulCoatingLayer>,
    pub substrate_receiver: SurfaceHeatingReceiver,
    pub substrate: SubstrateOpticsModel,
    pub incidence: SurfaceIncidence<N>,
    pub power: SpectralPowerSchedule<N>,
}

/// Finite films between fixed, real-index surrounding media. Every film reads
/// its own native heat receiver; no thermal cell is assigned to transmitted
/// radiation. Film and receiver order remains authored outside-to-inside even
/// when reverse incidence changes the physical traversal.
#[derive(Clone, Debug)]
pub struct StatefulTransmittingSurfaceIlluminationSource<const N: usize> {
    pub source_id: u64,
    pub films: Vec<StatefulCoatingLayer>,
    /// Incident-side indices and angle for the requested traversal direction.
    pub incidence: SurfaceIncidence<N>,
    pub transmitted_eta: [f64; N],
    pub reverse: bool,
    /// Source for the declared surrounding media and incidence geometry.
    pub source_identity: ContentHash,
    pub power: SpectralPowerSchedule<N>,
}

/// Borrow inputs for one mixed opaque/transmitting coating observer. All
/// retained laws, schedules and packets are owned by the resulting observer.
#[derive(Clone, Copy, Debug)]
pub enum CoatingSurfaceSource<'a, const N: usize> {
    FixedOpaque(&'a CoatedSurfaceIlluminationSource<N>),
    FixedTransmitting(&'a TransmittingSurfaceIlluminationSource<N>),
    StatefulOpaque(&'a StatefulCoatedSurfaceIlluminationSource<N>),
    StatefulTransmitting(&'a StatefulTransmittingSurfaceIlluminationSource<N>),
}

#[derive(Clone, Copy)]
enum SourceExit<'a, const N: usize> {
    Opaque {
        receiver: SurfaceHeatingReceiver,
        optics: &'a SubstrateOpticsModel,
    },
    Transmitting {
        transmitted_eta: [f64; N],
        reverse: bool,
        source_identity: ContentHash,
    },
}

impl<const N: usize> SourceExit<'_, N> {
    fn receiver(self) -> Option<SurfaceHeatingReceiver> {
        match self {
            Self::Opaque { receiver, .. } => Some(receiver),
            Self::Transmitting { .. } => None,
        }
    }
}

#[derive(Clone, Copy)]
struct StatefulSource<'a, const N: usize> {
    source_id: u64,
    films: &'a [StatefulCoatingLayer],
    exit: SourceExit<'a, N>,
    incidence: SurfaceIncidence<N>,
    power: &'a SpectralPowerSchedule<N>,
}

impl<'a, const N: usize> StatefulSource<'a, N> {
    fn receivers(self) -> impl Iterator<Item = SurfaceHeatingReceiver> + 'a {
        self.films
            .iter()
            .map(|film| film.receiver)
            .chain(self.exit.receiver())
    }
}

impl<'a, const N: usize> CoatingSurfaceSource<'a, N> {
    fn id(self) -> u64 {
        match self {
            Self::FixedOpaque(source) => source.source_id,
            Self::FixedTransmitting(source) => source.source_id,
            Self::StatefulOpaque(source) => source.source_id,
            Self::StatefulTransmitting(source) => source.source_id,
        }
    }
    fn power(self) -> &'a SpectralPowerSchedule<N> {
        match self {
            Self::FixedOpaque(source) => &source.power,
            Self::FixedTransmitting(source) => &source.power,
            Self::StatefulOpaque(source) => &source.power,
            Self::StatefulTransmitting(source) => &source.power,
        }
    }
    fn stateful(self) -> Option<StatefulSource<'a, N>> {
        match self {
            Self::StatefulOpaque(source) => Some(StatefulSource {
                source_id: source.source_id,
                films: &source.films,
                exit: SourceExit::Opaque {
                    receiver: source.substrate_receiver,
                    optics: &source.substrate,
                },
                incidence: source.incidence,
                power: &source.power,
            }),
            Self::StatefulTransmitting(source) => Some(StatefulSource {
                source_id: source.source_id,
                films: &source.films,
                exit: SourceExit::Transmitting {
                    transmitted_eta: source.transmitted_eta,
                    reverse: source.reverse,
                    source_identity: source.source_identity,
                },
                incidence: source.incidence,
                power: &source.power,
            }),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
enum Law {
    Film(Arc<SurfaceCoatingLaw>),
    Substrate(Arc<SurfaceConductorLaw>),
}
#[derive(Clone, Debug, PartialEq)]
enum Resolved {
    Film(CoherentAbsorbingCoating),
    Substrate(ConductorOptics),
}
impl Law {
    fn key(&self) -> (u8, ContentHash) {
        match self {
            Self::Film(law) => (0, law.identity()),
            Self::Substrate(law) => (1, law.identity()),
        }
    }
    fn rows(&self) -> usize {
        match self {
            Self::Film(law) => law.row_count(),
            Self::Substrate(law) => law.budget_size(),
        }
    }
    fn coordinate(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<f64, SurfaceHeatingError> {
        match self {
            Self::Film(law) => law.coordinate_at(thermal, cell, material),
            Self::Substrate(law) => law.coordinate_at(thermal, cell, material),
        }
    }
    fn resolve(&self, coordinate: f64) -> Result<Resolved, SurfaceHeatingError> {
        match self {
            Self::Film(law) => Ok(Resolved::Film(law.coating_at_coordinate(coordinate)?)),
            Self::Substrate(law) => Ok(Resolved::Substrate(law.optics_at_coordinate(coordinate)?)),
        }
    }
    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<Option<ConductorOptics>, SurfaceHeatingError> {
        match self {
            Self::Film(law) => law.coordinate_at(thermal, cell, material).map(|_| None),
            Self::Substrate(law) => law.validate_native(thermal, cell, material),
        }
    }
}

#[derive(Debug)]
struct BoundLaw {
    law: Law,
    receiver: SurfaceHeatingReceiver,
    groups: Vec<usize>,
}
#[derive(Clone, Debug)]
enum Film {
    Fixed(CoherentAbsorbingCoating),
    Lossless(CoherentDielectricCoating),
    Law(usize),
}
#[derive(Clone, Debug)]
enum Substrate {
    Fixed(ConductorOptics),
    Law(usize),
}

#[derive(Clone, Debug)]
enum Exit<const N: usize> {
    Opaque(Substrate),
    Transmitting {
        transmitted_eta: [f64; N],
        reverse: bool,
        source_identity: ContentHash,
    },
}

#[derive(Debug)]
struct Stack<const N: usize> {
    films: Box<[Film]>,
    exit: Exit<N>,
    incidence: SurfaceIncidence<N>,
    cells: Box<[CellId]>,
}
impl<const N: usize> Stack<N> {
    fn response(
        &self,
        resolved: &[Resolved],
    ) -> Result<Arc<CoatingPacket<N>>, SurfaceHeatingError> {
        let layers: Vec<_> = self
            .films
            .iter()
            .map(|film| match film {
                Film::Fixed(optics) => CoatingLayer::Absorbing(optics),
                Film::Lossless(optics) => CoatingLayer::Dielectric(optics),
                Film::Law(index) => match &resolved[*index] {
                    Resolved::Film(optics) => CoatingLayer::Absorbing(optics),
                    _ => unreachable!("film binding retains its native optical type"),
                },
            })
            .collect();
        match &self.exit {
            Exit::Opaque(substrate) => {
                let substrate = match substrate {
                    Substrate::Fixed(optics) => optics,
                    Substrate::Law(index) => match &resolved[*index] {
                        Resolved::Substrate(optics) => optics,
                        _ => unreachable!("substrate binding retains its native optical type"),
                    },
                };
                let response = CoatedSurfacePowerResponse::from_coated_conductor(
                    self.incidence.wavelengths_nm,
                    substrate,
                    CoherentCoatingStack::try_new(&layers)?,
                    self.incidence.incident_eta,
                    self.incidence.incidence_cosine,
                    self.incidence.s_power_fraction,
                )?;
                Ok(self.packet(
                    *response.total(),
                    response
                        .film_absorbed_fractions()
                        .iter()
                        .copied()
                        .chain(core::iter::once(*response.substrate_absorbed_fraction())),
                ))
            }
            Exit::Transmitting {
                transmitted_eta,
                reverse,
                source_identity,
            } => {
                let response = TransmittingSurfacePowerResponse::from_dielectric_stack(
                    self.incidence.wavelengths_nm,
                    &layers,
                    self.incidence.incident_eta,
                    *transmitted_eta,
                    self.incidence.incidence_cosine,
                    *reverse,
                    self.incidence.s_power_fraction,
                    *source_identity,
                )?;
                Ok(self.packet(
                    *response.total(),
                    response.film_absorbed_fractions().iter().copied(),
                ))
            }
        }
    }

    fn packet(
        &self,
        total: SurfacePowerResponse<N>,
        fractions: impl Iterator<Item = [f64; N]>,
    ) -> Arc<CoatingPacket<N>> {
        Arc::new(CoatingPacket {
            total,
            receivers: self
                .cells
                .iter()
                .zip(fractions)
                .map(|(&cell, absorbed)| AbsorptionReceiver { cell, absorbed })
                .collect(),
        })
    }
}

#[derive(Debug)]
pub(super) struct CoatingPacket<const N: usize> {
    pub(super) total: SurfacePowerResponse<N>,
    pub(super) receivers: Box<[AbsorptionReceiver<N>]>,
}
#[derive(Debug)]
struct FeedbackModel<const N: usize> {
    laws: Vec<BoundLaw>,
    stacks: Vec<Stack<N>>,
    sources: Vec<Option<usize>>,
}
#[derive(Clone, Debug)]
struct Cache<const N: usize> {
    resolved: Vec<Resolved>,
    packets: Vec<Option<Arc<CoatingPacket<N>>>>,
}
#[derive(Clone, Debug)]
pub(super) struct CoatingFeedbackState<const N: usize> {
    model: Arc<FeedbackModel<N>>,
    cache: Arc<Cache<N>>,
}

#[derive(Default)]
struct BindingBuilder {
    slots: BTreeMap<(u8, CellId, MaterialId, [u8; 32]), usize>,
    identities: BTreeSet<(u8, [u8; 32])>,
    rows: usize,
    laws: Vec<BoundLaw>,
    resolved: Vec<Resolved>,
}
impl BindingBuilder {
    fn bind(
        &mut self,
        law: Law,
        receiver: SurfaceHeatingReceiver,
        thermal: &DepositionSimulation<'_>,
        max_rows: usize,
    ) -> Result<usize, SurfaceHeatingError> {
        let (kind, identity) = law.key();
        let key = (
            kind,
            receiver.thermal_cell,
            receiver.material_id,
            *identity.as_bytes(),
        );
        if let Some(&index) = self.slots.get(&key) {
            return Ok(index);
        }
        if self.identities.insert((kind, *identity.as_bytes())) {
            self.rows = self
                .rows
                .checked_add(law.rows())
                .ok_or(SurfaceHeatingError::Budget)?;
            if self.rows > max_rows {
                return Err(SurfaceHeatingError::Budget);
            }
        }
        let coordinate = law.coordinate(thermal, receiver.thermal_cell, receiver.material_id)?;
        let resolved = law.resolve(coordinate)?;
        let index = self.laws.len();
        self.laws.push(BoundLaw {
            law,
            receiver,
            groups: Vec::new(),
        });
        self.resolved.push(resolved);
        self.slots.insert(key, index);
        Ok(index)
    }
}

fn receiver_hash(hash: &mut DomainHasher, receiver: SurfaceHeatingReceiver, bounds: bool) {
    hash.update(&receiver.thermal_cell.0.to_le_bytes());
    hash.update(&receiver.material_id.0.to_le_bytes());
    if bounds {
        hash.update(&receiver.temperature_bounds_k.0.to_bits().to_le_bytes());
        hash.update(&receiver.temperature_bounds_k.1.to_bits().to_le_bytes());
    }
}
fn stack_identity<const N: usize>(source: StatefulSource<'_, N>) -> ContentHash {
    let mut hash = DomainHasher::new(match source.exit {
        SourceExit::Opaque { .. } => "org.frankensim.render.coating-feedback-stack.v1",
        SourceExit::Transmitting { .. } => "org.frankensim.render.transmitting-feedback-stack.v1",
    });
    hash.update(&(N as u64).to_le_bytes());
    hash.update(&source.incidence.incidence_cosine.to_bits().to_le_bytes());
    for values in [
        &source.incidence.wavelengths_nm,
        &source.incidence.incident_eta,
        &source.incidence.s_power_fraction,
    ] {
        for value in values {
            hash.update(&value.to_bits().to_le_bytes());
        }
    }
    hash.update(&(source.films.len() as u64).to_le_bytes());
    for film in source.films {
        receiver_hash(&mut hash, film.receiver, false);
        match &film.optics {
            CoatingOpticsModel::Fixed(optics) => {
                hash.update(&[0]);
                let layers = [CoatingLayer::Absorbing(optics)];
                hash.update(
                    CoherentCoatingStack::try_new(&layers)
                        .expect("one admitted film fits the native bound")
                        .content_identity()
                        .as_bytes(),
                );
            }
            CoatingOpticsModel::Lossless(optics) => {
                hash.update(&[2]);
                let layers = [CoatingLayer::Dielectric(optics)];
                hash.update(
                    CoherentCoatingStack::try_new(&layers)
                        .expect("one admitted film fits the native bound")
                        .content_identity()
                        .as_bytes(),
                );
            }
            CoatingOpticsModel::Law(law) => {
                hash.update(&[1]);
                hash.update(law.identity().as_bytes());
            }
        }
    }
    match source.exit {
        SourceExit::Opaque { receiver, optics } => {
            receiver_hash(&mut hash, receiver, false);
            match optics {
                SubstrateOpticsModel::Fixed(optics) => {
                    hash.update(&[0]);
                    hash.update(optics.content_identity().as_bytes());
                }
                SubstrateOpticsModel::Law(law) => {
                    hash.update(&[1]);
                    hash.update(law.identity().as_bytes());
                }
            }
        }
        SourceExit::Transmitting {
            transmitted_eta,
            reverse,
            source_identity,
        } => {
            hash.update(source_identity.as_bytes());
            hash.update(&[u8::from(reverse)]);
            for value in transmitted_eta {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
    }
    hash.finalize()
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Couple fixed and state-dependent coherent stacks through the same native
    /// surface observer. Each film and substrate reads its own receiving cell;
    /// fixed entries keep their explicitly declared temperature support.
    ///
    /// `max_total_optical_knots` counts each distinct film/substrate law's rows
    /// once. Source/segment/evaluation budgets are the existing native budgets.
    /// Repeated identical stacks share response packets, including after a law
    /// changes. Accepted-start coefficients require temporal refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_coated_with_stateful(
        fixed: &[CoatedSurfaceIlluminationSource<N>],
        stateful: &[StatefulCoatedSurfaceIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        max_total_optical_knots: usize,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        let count = fixed
            .len()
            .checked_add(stateful.len())
            .ok_or(SurfaceHeatingError::Budget)?;
        if count > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        if stateful.is_empty() {
            return Self::try_new_coated(fixed, thermal, budget, tolerance, checkpoint);
        }
        let mut sources = Vec::new();
        sources
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        sources.extend(fixed.iter().map(CoatingSurfaceSource::FixedOpaque));
        sources.extend(stateful.iter().map(CoatingSurfaceSource::StatefulOpaque));
        Self::try_new_coating_stacks(
            &sources,
            thermal,
            budget,
            max_total_optical_knots,
            tolerance,
            checkpoint,
        )
    }

    /// Compose fixed and independently state-dependent opaque/transmitting
    /// stacks through one native observer and shared response cache. Every
    /// finite film retains its authored heat receiver under reverse incidence;
    /// a transmitting exit has no thermal destination. Distinct films may use
    /// temperature or exact-caloric-enthalpy laws, including latent plateaus.
    ///
    /// Row budgets count each distinct law once. Changing a shared receiver/law
    /// invalidates only its dependent stack packets; identical physical bindings
    /// and geometry share their next solve. Dark domains remain binding without
    /// rebuilding optical packets. Coefficients use accepted interval starts,
    /// so nonlinear feedback requires temporal refinement.
    pub fn try_new_coating_stacks(
        sources: &[CoatingSurfaceSource<'_, N>],
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        max_total_optical_knots: usize,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        let count = sources.len();
        if count > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut source_ids = BTreeSet::new();
        let mut segments = 0usize;
        for source in sources {
            checkpoint()?;
            if !source_ids.insert(source.id()) {
                return Err(invalid("surface illumination source IDs must be unique"));
            }
            segments = segments
                .checked_add(source.power().segment_count())
                .ok_or(SurfaceHeatingError::Budget)?;
            if segments > budget.max_total_power_segments {
                return Err(SurfaceHeatingError::Budget);
            }
        }
        if sources
            .iter()
            .filter_map(|source| source.stateful())
            .any(|source| source.films.len() > MAX_COHERENT_STACK_LAYERS)
        {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut ordered: Vec<_> = sources.iter().copied().collect();
        ordered.sort_unstable_by_key(|source| source.id());
        let mut bindings = BindingBuilder::default();
        let mut stacks = Vec::<Stack<N>>::new();
        let mut packets = Vec::new();
        let mut groups = BTreeMap::new();
        let mut by_source = BTreeMap::new();
        let mut source_hashes = BTreeMap::new();
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for input in ordered {
            checkpoint()?;
            let fixed = match input {
                CoatingSurfaceSource::FixedOpaque(source) => Some((
                    PreparedCoatingSource::opaque(source)?,
                    source.response.identity(),
                    "org.frankensim.render.fixed-coating-input.v1",
                )),
                CoatingSurfaceSource::FixedTransmitting(source) => Some((
                    PreparedCoatingSource::transmitting(source)?,
                    source.response.identity(),
                    "org.frankensim.render.fixed-transmitting-input.v1",
                )),
                _ => None,
            };
            if let Some((source, identity, domain)) = fixed {
                let mut hash = DomainHasher::new(domain);
                hash.update(identity.as_bytes());
                hash.update(source.power.identity().as_bytes());
                for (receiver, _) in &source.bindings {
                    receiver_hash(&mut hash, *receiver, true);
                }
                source_hashes.insert(source.source_id, hash.finalize());
                prepared.push(source);
                continue;
            }
            let source = input.stateful().expect("fixed inputs already prepared");
            let identity = stack_identity(source);
            let group = if let Some(&index) = groups.get(identity.as_bytes()) {
                index
            } else {
                let mut films = Vec::with_capacity(source.films.len());
                let mut cells = Vec::with_capacity(source.films.len() + 1);
                let mut law_slots = BTreeSet::new();
                for film in source.films {
                    checkpoint()?;
                    films.push(match &film.optics {
                        CoatingOpticsModel::Fixed(optics) => Film::Fixed(*optics),
                        CoatingOpticsModel::Lossless(optics) => Film::Lossless(*optics),
                        CoatingOpticsModel::Law(law) => {
                            let slot = bindings.bind(
                                Law::Film(Arc::clone(law)),
                                film.receiver,
                                thermal,
                                max_total_optical_knots,
                            )?;
                            law_slots.insert(slot);
                            Film::Law(slot)
                        }
                    });
                    cells.push(film.receiver.thermal_cell);
                }
                let exit = match source.exit {
                    SourceExit::Opaque { receiver, optics } => {
                        let substrate = match optics {
                            SubstrateOpticsModel::Fixed(optics) => Substrate::Fixed(*optics),
                            SubstrateOpticsModel::Law(law) => {
                                let slot = bindings.bind(
                                    Law::Substrate(Arc::clone(law)),
                                    receiver,
                                    thermal,
                                    max_total_optical_knots,
                                )?;
                                law_slots.insert(slot);
                                Substrate::Law(slot)
                            }
                        };
                        cells.push(receiver.thermal_cell);
                        Exit::Opaque(substrate)
                    }
                    SourceExit::Transmitting {
                        transmitted_eta,
                        reverse,
                        source_identity,
                    } => Exit::Transmitting {
                        transmitted_eta,
                        reverse,
                        source_identity,
                    },
                };
                let stack = Stack {
                    films: films.into_boxed_slice(),
                    exit,
                    incidence: source.incidence,
                    cells: cells.into_boxed_slice(),
                };
                let packet = stack.response(&bindings.resolved)?;
                let group = stacks.len();
                packets.push(Some(packet));
                for slot in law_slots {
                    bindings.laws[slot].groups.push(group);
                }
                stacks.push(stack);
                groups.insert(*identity.as_bytes(), group);
                group
            };
            by_source.insert(source.source_id, group);
            let packet = packets[group]
                .as_ref()
                .expect("admitted stack has its initial packet");
            prepared.push(PreparedCoatingSource {
                source_id: source.source_id,
                receiver: source.exit.receiver(),
                bindings: source
                    .receivers()
                    .zip(packet.receivers.iter().map(|receiver| receiver.absorbed))
                    .collect(),
                response: packet.total,
                power: source.power.clone(),
            });
            let mut hash = DomainHasher::new("org.frankensim.render.stateful-coating-input.v1");
            hash.update(identity.as_bytes());
            hash.update(source.power.identity().as_bytes());
            for receiver in source.receivers() {
                receiver_hash(&mut hash, receiver, true);
            }
            source_hashes.insert(source.source_id, hash.finalize());
        }
        let mut result =
            Self::try_new_partitioned(prepared, thermal, budget, tolerance, &mut checkpoint)?;
        if stacks.is_empty() {
            return Ok(result);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.stateful-coating-heating.v1");
        hash.update(&(N as u64).to_le_bytes());
        hash.update(&(count as u64).to_le_bytes());
        let sources = result
            .model
            .sources
            .iter()
            .map(|source| {
                hash.update(&source.source_id.to_le_bytes());
                hash.update(source_hashes[&source.source_id].as_bytes());
                by_source.get(&source.source_id).copied()
            })
            .collect();
        Arc::get_mut(&mut result.model)
            .expect("new surface model is uniquely owned")
            .identity = hash.finalize();
        result.coating_feedback = Some(CoatingFeedbackState {
            model: Arc::new(FeedbackModel {
                laws: bindings.laws,
                stacks,
                sources,
            }),
            cache: Arc::new(Cache {
                resolved: bindings.resolved,
                packets,
            }),
        });
        result.validate_receivers(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }
}

impl<const N: usize> CoatingFeedbackState<N> {
    pub(super) fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SurfaceHeatingError> {
        for bound in &self.model.laws {
            checkpoint()?;
            if let Some(optics) = bound.law.validate(
                thermal,
                bound.receiver.thermal_cell,
                bound.receiver.material_id,
            )? {
                // Film/row wavelength support is fixed and was admitted with
                // each stack. Original material wavelengths may move with T;
                // recheck every bound mode against this fresh substrate once.
                for &group in &bound.groups {
                    checkpoint()?;
                    for wavelength in self.model.stacks[group].incidence.wavelengths_nm {
                        optics.complex_ior(wavelength)?;
                    }
                }
            }
        }
        Ok(())
    }
    pub(super) fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SurfaceHeatingError> {
        for (index, bound) in self.model.laws.iter().enumerate() {
            checkpoint()?;
            let coordinate = bound.law.coordinate(
                thermal,
                bound.receiver.thermal_cell,
                bound.receiver.material_id,
            )?;
            let resolved = bound.law.resolve(coordinate)?;
            if resolved != self.cache.resolved[index] {
                let cache = Arc::make_mut(&mut self.cache);
                cache.resolved[index] = resolved;
                for &group in &bound.groups {
                    cache.packets[group] = None;
                }
            }
        }
        Ok(())
    }
    pub(super) fn response(
        &mut self,
        source: usize,
        illuminated: bool,
        evaluations: &mut usize,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<Arc<CoatingPacket<N>>>, SurfaceHeatingError> {
        let Some(group) = self.model.sources[source] else {
            return Ok(None);
        };
        if let Some(packet) = &self.cache.packets[group] {
            return Ok(Some(Arc::clone(packet)));
        }
        if !illuminated {
            return Ok(None);
        }
        checkpoint()?;
        let next = evaluations
            .checked_add(1)
            .ok_or(SurfaceHeatingError::Budget)?;
        let stack = &self.model.stacks[group];
        let packet = stack.response(&self.cache.resolved)?;
        checkpoint()?;
        Arc::make_mut(&mut self.cache).packets[group] = Some(Arc::clone(&packet));
        *evaluations = next;
        Ok(Some(packet))
    }
}
