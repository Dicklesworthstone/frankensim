//! Transient storage in an ordered heterogeneous pore-body/throat network.
//!
//! Each well-mixed pore body owns a ClosedGasParcel of volume epsilon*A*dx.
//! Separately authored homogeneous throats sample the existing dusty-gas law;
//! the shared transport owner updates all bodies and both finite reservoirs
//! together. Gas cannot be forwarded through another throat in the same snapshot.
//!
//! Body storage and throat transport have distinct physical geometry. Throat
//! length and area are explicit; no center-spacing or half-cell construction
//! replaces an interface law. Different throats retain different material and
//! pore sources. The resolved transient profile is the ordered body inventory,
//! not a continuous composition field inside a throat or across a laminate.
//!
//! Quasi-steady throats neglect their own gas storage. The caller bounds their
//! total gas-filled void volume relative to stored pore-body volume; the ratio
//! is reported. Represent non-negligible connecting storage as further bodies.
//! Do not count one physical void twice in independently authored geometries.
//! No scalar averaging or nonlinear layered interface solution is inferred.
//!
//! This first-order fixed-volume model supplies no temporal/spatial error
//! certificate. Reservoirs are finite and evolve. [`adsorption`] couples the
//! same time-window driver to the existing finite gas/adsorbent owner. Moving
//! porosity, capillary work and liquid storage remain separate physical models.

/// Finite competitive adsorption and skeleton heat on the same pore timeline.
pub mod adsorption;

use std::collections::BTreeMap;

use fs_matdb::QueryPoint;
use fs_thermochem::SpeciesId;

use super::super::super::diffusion::{
    BinaryDiffusionMaterial,
    dusty_gas::{DustyGasFaceSolve, DustyGasMedium},
    network::admit_entropy,
};
use super::super::super::{Sum, inventory, stored_change};
use super::super::{MAX_SNAPSHOT_PARCELS, SpeciesNetworkError, SpeciesNetworkStep};
use super::{
    GasTransportError, GasTransportPaths, GasTransportSolve, GasTransportStep,
    dusty_gas::DustyGasNetworkLink, step_candidate_with_paths,
};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergyError, GasEnergySolve, parcel::ClosedGasParcel},
};

pub const MAX_PORE_CELLS: usize = MAX_SNAPSHOT_PARCELS - 2;
pub const MAX_PORE_SUBSTEPS: usize = 4096;
pub const MAX_PORE_ATTEMPTS: usize = 16384;

/// Immutable storage geometry and its original pore-medium provenance.
#[derive(Clone, Debug, PartialEq)]
pub struct PoreCellGeometry {
    medium: DustyGasMedium,
}
impl PoreCellGeometry {
    pub fn try_new(medium: DustyGasMedium) -> Result<Self, PoreDomainError> {
        if !positive(medium.thickness_m())
            || !positive(medium.gas_filled_porosity())
            || medium.gas_filled_porosity() > 1.0
        {
            return Err(PoreDomainError::Invalid(
                "pore storage needs positive thickness and gas-filled porosity in (0,1]",
            ));
        }
        Ok(Self { medium })
    }
    pub fn medium(&self) -> &DustyGasMedium {
        &self.medium
    }
    pub fn thickness_m(&self) -> f64 {
        self.medium.thickness_m()
    }
    pub fn gas_filled_porosity(&self) -> f64 {
        self.medium.gas_filled_porosity()
    }
    pub fn volume_m3(&self, area_m2: f64) -> Result<f64, PoreDomainError> {
        let volume = void_volume(self.gas_filled_porosity(), area_m2, self.thickness_m());
        if !positive(area_m2) || !positive(volume) {
            return Err(PoreDomainError::Invalid(
                "positive representable gas-filled pore volume required",
            ));
        }
        Ok(volume)
    }
    /// Initialize through the existing material/EOS owner at the physical void volume.
    pub fn initialize_gas(
        &self,
        area_m2: f64,
        model: ThermochemicalGasModel,
        point: QueryPoint,
    ) -> Result<ClosedGasParcel, PoreDomainError> {
        ClosedGasParcel::new(model, point, self.volume_m3(area_m2)?).map_err(gas_error)
    }
}

/// One homogeneous physical throat. Body geometry does not infer its length or area.
#[derive(Clone, Copy, Debug)]
pub struct PoreFace<'a> {
    pub medium: &'a DustyGasMedium,
    /// Actual throat exposure [m²]; zero closes this connection.
    pub area_m2: f64,
    pub pairs: &'a [BinaryDiffusionMaterial<'a>],
    pub gas_context: &'a QueryPoint,
    pub solve: DustyGasFaceSolve,
}

/// Work bounds and physical admission. No local-error estimate is implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PoreDomainSolve {
    pub maximum_step_s: f64,
    /// Retry floor; a final remainder may be shorter.
    pub minimum_step_s: f64,
    pub maximum_substeps: usize,
    pub maximum_attempts: usize,
    /// Maximum neglected throat void volume / stored pore-body void volume.
    /// An explicit geometric approximation limit in [0,1), not an error estimate.
    pub maximum_throat_volume_fraction: f64,
    /// The same mass/energy/entropy limits also apply to the complete window.
    pub transport: GasTransportSolve,
}

/// Both finite reservoirs and all source-backed pore inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct PoreDomain {
    area_m2: f64,
    cells: Vec<PoreCellGeometry>,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    time_s: f64,
}

/// Gross directional throughput survives reversals within the window.
#[derive(Clone, Debug, PartialEq)]
pub struct PoreFaceTransfer {
    pub face: usize,
    pub species: SpeciesId,
    pub mass_left_to_right_kg: f64,
    pub mass_right_to_left_kg: f64,
    /// Signed integral of original donor partial enthalpy; formation references
    /// are already included, with no additional chemical heat.
    pub enthalpy_left_to_right_j: f64,
}

/// Bounded summary: source reports are retained for the last substep only.
#[derive(Clone, Debug, PartialEq)]
pub struct PoreDomainStep<Step = GasTransportStep> {
    pub duration_s: f64,
    pub accepted_substeps: usize,
    pub attempts: usize,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub faces: Vec<PoreFaceTransfer>,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub total_energy_residual_j: f64,
    pub entropy_change_j_k: f64,
    pub pore_mass_change_kg: f64,
    pub neglected_throat_volume_m3: f64,
    pub throat_volume_fraction: f64,
    pub last_substep: Step,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PoreDomainCandidate {
    domain: PoreDomain,
    step: PoreDomainStep,
}
impl PoreDomainCandidate {
    pub fn domain(&self) -> &PoreDomain {
        &self.domain
    }
    pub fn step(&self) -> &PoreDomainStep {
        &self.step
    }
    pub fn into_parts(self) -> (PoreDomain, PoreDomainStep) {
        (self.domain, self.step)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PoreDomainError {
    Invalid(&'static str),
    Gas(GasEnergyError),
    WindowBalance(GasTransportError),
    Transport {
        elapsed_s: f64,
        attempted_step_s: f64,
        error: GasTransportError,
    },
    /// Original joint gas/solid refusal at the actual attempted pore substep.
    Adsorption {
        elapsed_s: f64,
        attempted_step_s: f64,
        error: crate::adsorption::competitive::network::AdsorptionNetworkError,
    },
    WorkLimit {
        elapsed_s: f64,
        accepted_substeps: usize,
        attempts: usize,
    },
    Cancelled,
}
impl core::fmt::Display for PoreDomainError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient porous gas domain: {self:?}")
    }
}
impl std::error::Error for PoreDomainError {}
use PoreDomainError as Error;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
// Apply porosity to the larger geometric factor first to avoid intermediate
// underflow/overflow when the final physical void volume is representable.
fn void_volume(porosity: f64, area: f64, length: f64) -> f64 {
    (porosity * area.max(length)) * area.min(length)
}
fn same_volume(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 8.0 * f64::EPSILON * a.max(b)
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn gas_error(error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Gas(error)
    }
}
fn retryable(error: &GasTransportError) -> bool {
    matches!(
        error,
        GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { .. })
            | GasTransportError::Exchange(SpeciesNetworkError::Parcel {
                error: GasEnergyError::OutsideBracket { .. } | GasEnergyError::NotConverged { .. },
                ..
            })
    )
}

// The two native physical owners share only bounded stepping and face-throughput
// accounting. Their caloric, inventory and entropy admission remains their own.
trait PoreSubstep {
    fn exchange(&self) -> &SpeciesNetworkStep;
}
impl PoreSubstep for GasTransportStep {
    fn exchange(&self) -> &SpeciesNetworkStep {
        &self.exchange
    }
}
enum PoreAdvanceError {
    Transport(GasTransportError),
    Adsorption(crate::adsorption::competitive::network::AdsorptionNetworkError),
}
impl PoreAdvanceError {
    fn cancelled(&self) -> bool {
        match self {
            Self::Transport(error) => *error == GasTransportError::Cancelled,
            Self::Adsorption(error) => {
                *error == crate::adsorption::competitive::network::AdsorptionNetworkError::Cancelled
            }
        }
    }
    fn retryable(&self) -> bool {
        match self {
            Self::Transport(error) => retryable(error),
            Self::Adsorption(error) => adsorption::retryable(error),
        }
    }
    fn at(self, elapsed_s: f64, attempted_step_s: f64) -> Error {
        match self {
            Self::Transport(error) => Error::Transport {
                elapsed_s,
                attempted_step_s,
                error,
            },
            Self::Adsorption(error) => Error::Adsorption {
                elapsed_s,
                attempted_step_s,
                error,
            },
        }
    }
}
struct PoreWindow<Step> {
    domain: PoreDomain,
    duration_s: f64,
    accepted_substeps: usize,
    attempts: usize,
    minimum_accepted_step_s: f64,
    maximum_accepted_step_s: f64,
    faces: Vec<PoreFaceTransfer>,
    neglected_throat_volume_m3: f64,
    throat_volume_fraction: f64,
    last_substep: Step,
}
impl<Step> PoreWindow<Step> {
    fn finish(self, balances: (f64, f64, f64, f64, f64)) -> (PoreDomain, PoreDomainStep<Step>) {
        let (mass, species, energy, entropy, pore_mass) = balances;
        (
            self.domain,
            PoreDomainStep {
                duration_s: self.duration_s,
                accepted_substeps: self.accepted_substeps,
                attempts: self.attempts,
                minimum_accepted_step_s: self.minimum_accepted_step_s,
                maximum_accepted_step_s: self.maximum_accepted_step_s,
                faces: self.faces,
                total_mass_residual_kg: mass,
                maximum_species_mass_residual_kg: species,
                total_energy_residual_j: energy,
                entropy_change_j_k: entropy,
                pore_mass_change_kg: pore_mass,
                neglected_throat_volume_m3: self.neglected_throat_volume_m3,
                throat_volume_fraction: self.throat_volume_fraction,
                last_substep: self.last_substep,
            },
        )
    }
}

impl PoreDomain {
    /// Pore cells follow their physical left-to-right order. Each provided gas
    /// must already have epsilon*A*dx volume; initialization never rescales gas
    /// mass or erases the work of a supplied volume change.
    ///
    /// Internal IDs are immutable: left=0, pores=1..=n, right=n+1. An enclosing
    /// network must namespace this column's identities when composing domains.
    pub fn try_new(
        area_m2: f64,
        cells: Vec<PoreCellGeometry>,
        left: ClosedGasParcel,
        pore_gases: Vec<ClosedGasParcel>,
        right: ClosedGasParcel,
    ) -> Result<Self, Error> {
        if !positive(area_m2)
            || cells.is_empty()
            || cells.len() > MAX_PORE_CELLS
            || cells.len() != pore_gases.len()
        {
            return Err(Error::Invalid(
                "positive column area and 1..254 matched pore cells required",
            ));
        }
        let mut parcels = BTreeMap::new();
        parcels.insert(0, left);
        for (i, (geometry, gas)) in cells.iter().zip(pore_gases).enumerate() {
            if !same_volume(gas.volume_m3(), geometry.volume_m3(area_m2)?) {
                return Err(Error::Invalid(
                    "pore gas volume must equal the authored gas-filled epsilon*A*dx",
                ));
            }
            parcels.insert((i + 1) as u64, gas);
        }
        parcels.insert((cells.len() + 1) as u64, right);
        Ok(Self {
            area_m2,
            cells,
            parcels,
            time_s: 0.0,
        })
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn cells(&self) -> &[PoreCellGeometry] {
        &self.cells
    }
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn left_reservoir(&self) -> &ClosedGasParcel {
        &self.parcels[&0]
    }
    pub fn right_reservoir(&self) -> &ClosedGasParcel {
        &self.parcels[&((self.cells.len() + 1) as u64)]
    }
    pub fn pore_gas(&self, cell: usize) -> Option<&ClosedGasParcel> {
        (cell < self.cells.len()).then(|| &self.parcels[&((cell + 1) as u64)])
    }
    /// Propose a complete window without modifying any accepted gas.
    ///
    /// Only typed aggregate-withdrawal and caloric bracket/convergence refusals
    /// trigger bounded halving. Source/domain refusals propagate immediately.
    /// A failed later substep discards the entire proposed window. Maximum dt
    /// controls caller refinement; successful admission is not an error estimate.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PoreDomainCandidate, Error> {
        let run = self.run_window(
            duration_s,
            faces,
            endpoints,
            controls,
            &mut cancelled,
            |parcels, paths, dt, endpoints, controls, cancelled| {
                step_candidate_with_paths(parcels, paths, dt, endpoints, controls, cancelled)
                    .map(|candidate| candidate.into_parts())
                    .map_err(PoreAdvanceError::Transport)
            },
        )?;
        let balances = window_balances(
            &self.parcels,
            &run.domain.parcels,
            self.cells.len(),
            controls.transport,
            &mut cancelled,
        )?;
        check(&mut cancelled)?;
        let (domain, step) = run.finish(balances);
        Ok(PoreDomainCandidate { domain, step })
    }

    fn run_window<Step: PoreSubstep, C: FnMut() -> bool>(
        &self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        mut cancelled: C,
        mut advance: impl FnMut(
            &BTreeMap<u64, ClosedGasParcel>,
            GasTransportPaths<'_>,
            f64,
            &BTreeMap<u64, GasEnergySolve>,
            GasTransportSolve,
            &mut C,
        )
            -> Result<(BTreeMap<u64, ClosedGasParcel>, Step), PoreAdvanceError>,
    ) -> Result<PoreWindow<Step>, Error> {
        check(&mut cancelled)?;
        let end_time = self.time_s + duration_s;
        if !positive(duration_s)
            || !end_time.is_finite()
            || end_time <= self.time_s
            || !positive(controls.minimum_step_s)
            || !positive(controls.maximum_step_s)
            || controls.minimum_step_s > controls.maximum_step_s
            || !(1..=MAX_PORE_SUBSTEPS).contains(&controls.maximum_substeps)
            || !(1..=MAX_PORE_ATTEMPTS).contains(&controls.maximum_attempts)
            || !controls.maximum_throat_volume_fraction.is_finite()
            || !(0.0..1.0).contains(&controls.maximum_throat_volume_fraction)
            || faces.len() != self.cells.len() + 1
            || endpoints.len() != self.parcels.len()
        {
            return Err(Error::Invalid(
                "bounded pore window, exact faces/endpoints and positive ordered step limits required",
            ));
        }
        let endpoints: BTreeMap<_, _> = endpoints
            .iter()
            .enumerate()
            .map(|(i, &solve)| (i as u64, solve))
            .collect();
        for (&id, solve) in &endpoints {
            check(&mut cancelled)?;
            solve.validate().map_err(gas_error)?;
            let t = self.parcels[&id].state().state().temperature;
            if t < solve.minimum_temperature_k || t > solve.maximum_temperature_k {
                return Err(Error::Invalid(
                    "each caloric bracket must contain its original pore/reservoir temperature",
                ));
            }
        }
        let mut stored_volume = Sum::default();
        for geometry in &self.cells {
            check(&mut cancelled)?;
            stored_volume
                .add(geometry.volume_m3(self.area_m2)?)
                .map_err(gas_error)?;
        }
        let mut throat_volume = Sum::default();
        let mut links = Vec::with_capacity(faces.len());
        for (i, face) in faces.iter().enumerate() {
            check(&mut cancelled)?;
            if !face.area_m2.is_finite() || face.area_m2 < 0.0 {
                return Err(Error::Invalid(
                    "pore throat area must be finite and nonnegative",
                ));
            }
            let void = void_volume(
                face.medium.gas_filled_porosity(),
                face.area_m2,
                face.medium.thickness_m(),
            );
            if !void.is_finite()
                || void < 0.0
                || (face.medium.gas_filled_porosity() > 0.0 && face.area_m2 > 0.0 && void == 0.0)
            {
                return Err(Error::Invalid("pore throat void volume is unrepresentable"));
            }
            throat_volume.add(void).map_err(gas_error)?;
            links.push(DustyGasNetworkLink {
                id: i as u64,
                parcels: [i as u64, (i + 1) as u64],
                medium: face.medium,
                area_m2: face.area_m2,
                pairs: face.pairs,
                gas_context: face.gas_context,
                solve: face.solve,
            });
        }
        let throat_volume_fraction = throat_volume.value / stored_volume.value;
        if !throat_volume_fraction.is_finite()
            || (throat_volume.value > 0.0 && throat_volume_fraction == 0.0)
            || throat_volume_fraction > controls.maximum_throat_volume_fraction
        {
            return Err(Error::Invalid(
                "neglected throat storage exceeds its declared pore-body volume fraction",
            ));
        }
        let paths = GasTransportPaths {
            dusty_gas: &links,
            ..Default::default()
        };
        let mut parcels = self.parcels.clone();
        let mut elapsed = 0.0;
        let mut ceiling = controls.maximum_step_s;
        let (mut accepted, mut attempts) = (0usize, 0usize);
        let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
        let mut totals: BTreeMap<(usize, SpeciesId), [Sum; 3]> = BTreeMap::new();
        let mut last = None;
        while elapsed < duration_s {
            check(&mut cancelled)?;
            if accepted == controls.maximum_substeps || attempts == controls.maximum_attempts {
                return Err(Error::WorkLimit {
                    elapsed_s: elapsed,
                    accepted_substeps: accepted,
                    attempts,
                });
            }
            let remaining = duration_s - elapsed;
            let dt =
                if remaining <= ceiling || remaining - ceiling <= 8.0 * f64::EPSILON * duration_s {
                    remaining
                } else {
                    ceiling
                };
            if !positive(dt) || elapsed + dt <= elapsed {
                return Err(Error::Invalid(
                    "pore substep no longer advances representable physical time",
                ));
            }
            attempts += 1;
            match advance(
                &parcels,
                paths,
                dt,
                &endpoints,
                controls.transport,
                &mut cancelled,
            ) {
                Ok((next, step)) => {
                    for row in &step.exchange().fluxes {
                        check(&mut cancelled)?;
                        let edge = row.transfer.donor.min(row.transfer.receiver);
                        if row.transfer.donor.max(row.transfer.receiver) != edge + 1
                            || edge as usize >= faces.len()
                        {
                            return Err(Error::Invalid(
                                "pore transport returned a nonlocal face transfer",
                            ));
                        }
                        let forward = row.transfer.donor == edge;
                        let sums = totals
                            .entry((edge as usize, row.transfer.species.clone()))
                            .or_insert_with(|| [Sum::default(), Sum::default(), Sum::default()]);
                        sums[if forward { 0 } else { 1 }]
                            .add(row.transfer.mass_kg)
                            .map_err(gas_error)?;
                        sums[2]
                            .add(if forward {
                                row.carried_enthalpy_j
                            } else {
                                -row.carried_enthalpy_j
                            })
                            .map_err(gas_error)?;
                    }
                    parcels = next;
                    elapsed = if dt == remaining {
                        duration_s
                    } else {
                        elapsed + dt
                    };
                    accepted += 1;
                    minimum = minimum.min(dt);
                    maximum = maximum.max(dt);
                    ceiling = if dt < 0.5 * controls.maximum_step_s {
                        2.0 * dt
                    } else {
                        controls.maximum_step_s
                    };
                    last = Some(step);
                }
                Err(error) if error.cancelled() => return Err(Error::Cancelled),
                Err(error) => {
                    let smaller = 0.5 * dt;
                    if error.retryable()
                        && positive(smaller)
                        && smaller >= controls.minimum_step_s
                        && elapsed + smaller > elapsed
                    {
                        ceiling = smaller;
                    } else {
                        return Err(error.at(elapsed, dt));
                    }
                }
            }
        }
        let faces = totals
            .into_iter()
            .map(|((face, species), sums)| PoreFaceTransfer {
                face,
                species,
                mass_left_to_right_kg: sums[0].value,
                mass_right_to_left_kg: sums[1].value,
                enthalpy_left_to_right_j: sums[2].value,
            })
            .collect();
        check(&mut cancelled)?;
        Ok(PoreWindow {
            domain: Self {
                area_m2: self.area_m2,
                cells: self.cells.clone(),
                parcels,
                time_s: end_time,
            },
            duration_s,
            accepted_substeps: accepted,
            attempts,
            minimum_accepted_step_s: minimum,
            maximum_accepted_step_s: maximum,
            faces,
            neglected_throat_volume_m3: throat_volume.value,
            throat_volume_fraction,
            last_substep: last.ok_or(Error::Invalid("pore window accepted no physical substep"))?,
        })
    }

    /// Publish all pore cells, reservoirs and the clock only after the full window.
    pub fn advance(
        &mut self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PoreDomainStep, Error> {
        let candidate =
            self.step_candidate(duration_s, faces, endpoints, controls, &mut cancelled)?;
        check(&mut cancelled)?;
        let (next, step) = candidate.into_parts();
        *self = next;
        Ok(step)
    }
}

fn window_balances(
    old: &BTreeMap<u64, ClosedGasParcel>,
    next: &BTreeMap<u64, ClosedGasParcel>,
    pore_count: usize,
    controls: GasTransportSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64, f64, f64, f64), Error> {
    let (mut mass, mut energy, mut pore_mass) = (Sum::default(), Sum::default(), Sum::default());
    let mut species: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    for (&id, before) in old {
        check(cancelled)?;
        let after = &next[&id];
        let dm = after.mass_kg() - before.mass_kg();
        mass.add(dm).map_err(gas_error)?;
        energy
            .add(stored_change(before, after))
            .map_err(gas_error)?;
        if id > 0 && id <= pore_count as u64 {
            pore_mass.add(dm).map_err(gas_error)?;
        }
        let a = inventory(before.model(), before.mass_kg()).map_err(gas_error)?;
        let b = inventory(after.model(), after.mass_kg()).map_err(gas_error)?;
        for (name, value) in &b {
            check(cancelled)?;
            species
                .entry(name.clone())
                .or_default()
                .add(value - a.get(name).copied().unwrap_or(0.0))
                .map_err(gas_error)?;
        }
        for (name, value) in &a {
            if !b.contains_key(name) {
                species
                    .entry(name.clone())
                    .or_default()
                    .add(-value)
                    .map_err(gas_error)?;
            }
        }
    }
    let species_error = species
        .values()
        .map(|s| s.value.abs())
        .fold(0.0_f64, f64::max);
    if mass.value.abs() > controls.exchange.mass_tolerance_kg
        || species_error > controls.exchange.mass_tolerance_kg
        || energy.value.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "complete pore window exceeded its species/mass/energy balance budget",
        ));
    }
    let entropy = admit_entropy(
        old,
        next,
        controls.entropy_tolerance_j_k,
        controls.exchange.mass_tolerance_kg,
        &mut *cancelled,
    )
    .map_err(|error| match error {
        super::super::super::diffusion::network::DiffusionNetworkError::Cancelled => {
            Error::Cancelled
        }
        other => Error::WindowBalance(GasTransportError::Diffusion(other)),
    })?;
    check(cancelled)?;
    Ok((
        mass.value,
        species_error,
        energy.value,
        entropy,
        pore_mass.value,
    ))
}

#[cfg(test)]
mod tests;
