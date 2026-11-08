//! Original interacting bubble histories feeding coherent compact acoustics.
//!
//! Every source uses its jointly solved volume acceleration and its own fixed
//! center. Summing accelerations at one arbitrary distance would discard the
//! distinct retarded paths. Feed the exported histories and stationary motions
//! to the existing moving compact renderer to sum pressures at one receiver time.
//!
//! The native cluster retains its incompressible, separated-sphere, long-wave
//! interaction domain. Propagation is one-way and adds no radiation loss or
//! feedback. Temperature and gas inventories remain those of the original
//! members; no translating bubbles, cloud shocks or erosion are inferred.
//! The volume-source pressure relation is Brennen, Cavitation and Bubble
//! Dynamics, chapter 3, section 3.8, equation 3.10:
//! <https://media.library.caltech.edu/CaltechBOOK:1995.001/chap3.htm>.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::{CompactSourceStrength, TimedCompactSource, UniformSourceMotion};
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits, MonopoleSample};
use fs_flux::cavitation::{
    BubbleCluster, BubbleClusterError, BubbleClusterNumerics, BubbleClusterObservation,
    BubbleClusterStep,
};
use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState};
use fs_material::fluid_acoustics::FluidAcoustics;

/// An original owner refused the candidate; no member state has been accepted.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialBubbleClusterError {
    Material(MaterialError),
    Dynamics(BubbleClusterError),
    Acoustic(AeroacError),
    /// Bind the complete liquid snapshot, including viscosity and capillarity.
    MaterialMismatch,
    /// An observer must retain the same material, centers and source bandwidth.
    ObserverMismatch,
}

impl core::fmt::Display for MaterialBubbleClusterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "bubble cluster material: {error}"),
            Self::Dynamics(error) => write!(f, "bubble cluster dynamics: {error}"),
            Self::Acoustic(error) => write!(f, "bubble cluster acoustics: {error}"),
            Self::MaterialMismatch => {
                write!(
                    f,
                    "bubble cluster must retain the original complete liquid state"
                )
            }
            Self::ObserverMismatch => {
                write!(
                    f,
                    "bubble cluster observer material or source geometry changed"
                )
            }
        }
    }
}
impl std::error::Error for MaterialBubbleClusterError {}
impl From<MaterialError> for MaterialBubbleClusterError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<BubbleClusterError> for MaterialBubbleClusterError {
    fn from(error: BubbleClusterError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<AeroacError> for MaterialBubbleClusterError {
    fn from(error: AeroacError) -> Self {
        Self::Acoustic(error)
    }
}

/// One jointly accepted interval and its exact native per-member volume sources.
#[derive(Clone, Debug, PartialEq)]
pub struct CompactClusterStep {
    pub cluster: BubbleClusterStep,
    /// Original member order, at the common accepted source time.
    pub sources: Vec<TimedCompactSource>,
}

/// Joint mechanics and individually delayed acoustic samples.
///
/// Members generally have DIFFERENT arrival times. Do not sum these endpoints
/// as though they shared a receiver clock; render their retained histories.
#[derive(Clone, Debug, PartialEq)]
pub struct ObservedClusterStep {
    pub cluster: BubbleClusterStep,
    pub sources: Vec<TimedCompactSource>,
    pub acoustic: Vec<MonopoleSample>,
}

/// Prepared per-member propagation from restrained centers to one receiver.
pub struct ClusterMonopoleObserver {
    material: CavitationState,
    centers_m: Vec<[f64; 3]>,
    receiver_m: [f64; 3],
    maximum_source_angular_frequency_rad_s: f64,
    paths: Vec<CompactMonopoleObserver>,
}

impl ClusterMonopoleObserver {
    pub const fn receiver_position_m(&self) -> [f64; 3] {
        self.receiver_m
    }

    fn matches(&self, cluster: &BubbleCluster) -> bool {
        self.material == cluster.material()
            && self.maximum_source_angular_frequency_rad_s
                == cluster.limits().maximum_source_angular_frequency_rad_s
            && self.centers_m.len() == cluster.members().len()
            && self
                .centers_m
                .iter()
                .zip(cluster.members())
                .all(|(center, member)| *center == member.center_m)
    }
}

/// Source-retaining binding to the original joint nonlinear physical owner.
pub struct MaterialBubbleCluster<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    cluster: BubbleCluster,
}

impl<'a, F: CavitationProperties + ?Sized> MaterialBubbleCluster<'a, F> {
    /// Bind an already admitted cluster without rebuilding any gas inventory,
    /// member equation, radius, velocity or common physical clock.
    pub fn new(source: &'a F, cluster: BubbleCluster) -> Result<Self, MaterialBubbleClusterError> {
        if source.cavitation_state()? != cluster.material() {
            return Err(MaterialBubbleClusterError::MaterialMismatch);
        }
        Ok(Self { source, cluster })
    }

    pub const fn source(&self) -> &'a F {
        self.source
    }
    pub const fn material(&self) -> CavitationState {
        self.cluster.material()
    }
    pub const fn cluster(&self) -> &BubbleCluster {
        &self.cluster
    }

    /// Caller-declared native source bandwidth, including nonlinear harmonics.
    pub fn source_bandwidth_hz(&self) -> f64 {
        self.cluster.limits().maximum_source_angular_frequency_rad_s / core::f64::consts::TAU
    }

    /// Exact restrained centers, in the same order as every exported history.
    pub fn source_motions(&self) -> Result<Vec<UniformSourceMotion>, MaterialBubbleClusterError> {
        self.cluster
            .members()
            .iter()
            .map(|member| {
                Ok(UniformSourceMotion::new(
                    self.cluster.time_s(),
                    member.center_m,
                    [0.0; 3],
                )?)
            })
            .collect()
    }

    /// Prepare original compact monopoles using the same liquid rho and c.
    /// The complete admitted radius envelope must satisfy distance/compactness;
    /// each actual endpoint additionally passes the native source-strength gate.
    pub fn monopole_observer(
        &self,
        receiver_m: [f64; 3],
        limits: MonopoleLimits,
    ) -> Result<ClusterMonopoleObserver, MaterialBubbleClusterError> {
        if receiver_m.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(AeroacError::NonFinite {
                what: "cluster receiver position",
            }
            .into());
        }
        let medium = self.material().acoustic_state()?;
        let mut paths = Vec::with_capacity(self.cluster.members().len());
        let mut centers_m = Vec::with_capacity(self.cluster.members().len());
        for member in self.cluster.members() {
            let center = member.center_m;
            let distance = (receiver_m[0] - center[0])
                .hypot(receiver_m[1] - center[1])
                .hypot(receiver_m[2] - center[2]);
            let path =
                CompactMonopoleObserver::new(medium, distance, self.source_bandwidth_hz(), limits)?;
            path.sample(
                self.cluster.time_s(),
                member.bubble.limits().maximum_radius_m,
                0.0,
            )?;
            paths.push(path);
            centers_m.push(center);
        }
        Ok(ClusterMonopoleObserver {
            material: self.material(),
            centers_m,
            receiver_m,
            maximum_source_angular_frequency_rad_s: self
                .cluster
                .limits()
                .maximum_source_angular_frequency_rad_s,
            paths,
        })
    }

    /// Read the current joint source without advancing any member.
    pub fn compact_source_samples(
        &self,
        absolute_pressure_pa: f64,
    ) -> Result<Vec<TimedCompactSource>, MaterialBubbleClusterError> {
        source_samples(&self.cluster.observe(absolute_pressure_pa)?, &mut || false)
    }

    /// Advance once, then export those exact coupled sources. As with the
    /// native cluster, no state changes on failure or cancellation. The caller
    /// owns history sampling and the eventual observer's propagation admission.
    pub fn advance_compact_sources<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &BubbleClusterNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<CompactClusterStep, MaterialBubbleClusterError>
    where
        Pressure: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        let proposal = self
            .cluster
            .propose_controlled(end_time_s, numerics, pressure, cancelled)?;
        let sources = source_samples(&proposal.step().observation, cancelled)?;
        poll(cancelled)?;
        Ok(CompactClusterStep {
            cluster: proposal.commit(),
            sources,
        })
    }

    /// Admit every individual acoustic endpoint before accepting any member.
    /// This neither reintegrates the trajectory nor clones a second solver.
    pub fn advance_observed<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &BubbleClusterNumerics,
        pressure: &Pressure,
        observer: &ClusterMonopoleObserver,
        cancelled: &mut Cancel,
    ) -> Result<ObservedClusterStep, MaterialBubbleClusterError>
    where
        Pressure: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        if !observer.matches(&self.cluster) {
            return Err(MaterialBubbleClusterError::ObserverMismatch);
        }
        let proposal = self
            .cluster
            .propose_controlled(end_time_s, numerics, pressure, cancelled)?;
        let observation = &proposal.step().observation;
        let sources = source_samples(observation, cancelled)?;
        let mut acoustic = Vec::with_capacity(observation.bubbles.len());
        for (path, bubble) in observer.paths.iter().zip(&observation.bubbles) {
            poll(cancelled)?;
            acoustic.push(path.sample(
                observation.time_s,
                bubble.state.radius_m(),
                bubble.volume_acceleration_m3_s2,
            )?);
        }
        poll(cancelled)?;
        Ok(ObservedClusterStep {
            cluster: proposal.commit(),
            sources,
            acoustic,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialBubbleClusterError> {
    if cancelled() {
        Err(BubbleClusterError::Cancelled.into())
    } else {
        Ok(())
    }
}

fn source_samples(
    observation: &BubbleClusterObservation,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<TimedCompactSource>, MaterialBubbleClusterError> {
    let mut sources = Vec::with_capacity(observation.bubbles.len());
    for bubble in &observation.bubbles {
        poll(cancelled)?;
        sources.push(TimedCompactSource {
            source_time_s: observation.time_s,
            strength: CompactSourceStrength {
                volume_flow_m3_s: bubble.volume_rate_m3_s,
                volume_flow_rate_m3_s2: bubble.volume_acceleration_m3_s2,
                force_on_fluid_n: [0.0; 3],
                force_rate_on_fluid_n_s: [0.0; 3],
            },
        });
    }
    Ok(sources)
}
