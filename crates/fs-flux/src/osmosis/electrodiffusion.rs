//! Finite-rate ionic membrane transport coupled to osmotic water and mechanics.
//!
//! Every original mobile species retains its own finite-rate equation. Two
//! opposite-valence ions use a scalar salt fast path; larger layouts use a
//! prepared charge-constrained Newton solve. Both share the original scalar
//! water solve and backward-Euler endpoint. Initial inventories are retained;
//! no instantaneous Donnan partition replaces finite kinetics.
//! Independent selective water pores and ion channels share the original
//! solution potential. Fixed charges and separately declared trapped particles
//! stay with their compartments. No solvent drag or electrical storage is inferred.

pub mod adaptive;
pub use adaptive::{IonicAdvanceControls, IonicAdvanceReport, IonicAdvanceWorkspace};
mod multicomponent;

use super::{AqueousSolution, OsmoticCompartment, OsmoticNetworkError, compensated_add};
use fs_material::osmosis::{
    OsmosisError,
    electrodiffusion::{
        ElectrodiffusionControls, ElectrodiffusionError, ElectrodiffusionProblem,
        ElectrodiffusionReport, ElectrodiffusiveMembrane,
    },
};
use fs_math::det;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

#[derive(Debug)]
pub enum IonicOsmoticError {
    Network(OsmoticNetworkError),
    Membrane(ElectrodiffusionError),
    Charge {
        side: usize,
        residual_mol: f64,
        maximum_mol: f64,
    },
    TemporalAccuracy {
        step_s: f64,
        scaled_error: f64,
    },
}
impl core::fmt::Display for IonicOsmoticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ionic osmosis: {self:?}")
    }
}
impl std::error::Error for IonicOsmoticError {}
impl From<OsmoticNetworkError> for IonicOsmoticError {
    fn from(value: OsmoticNetworkError) -> Self {
        Self::Network(value)
    }
}
impl From<OsmosisError> for IonicOsmoticError {
    fn from(value: OsmosisError) -> Self {
        Self::Network(value.into())
    }
}
impl From<ElectrodiffusionError> for IonicOsmoticError {
    fn from(value: ElectrodiffusionError) -> Self {
        match value {
            ElectrodiffusionError::Cancelled => Self::Network(OsmoticNetworkError::Cancelled),
            other => Self::Membrane(other),
        }
    }
}
type Result<T> = core::result::Result<T, IonicOsmoticError>;
const CAPACITY: usize = super::MAX_OSMOTIC_SPECIES;

/// Allocation-free values in the original species ordering. Indexing, slices
/// and iteration expose exactly the active species, never unused capacity.
#[derive(Clone, Copy, PartialEq)]
pub struct IonicSpeciesValues {
    values: [f64; CAPACITY],
    count: usize,
}
impl core::fmt::Debug for IonicSpeciesValues {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Debug::fmt(self.as_slice(), f)
    }
}
impl IonicSpeciesValues {
    fn zeroed(count: usize) -> Self {
        Self {
            values: [0.0; CAPACITY],
            count,
        }
    }
    fn from_fn(count: usize, mut value: impl FnMut(usize) -> f64) -> Self {
        let mut result = Self::zeroed(count);
        for i in 0..count {
            result.values[i] = value(i);
        }
        result
    }
    pub const fn len(&self) -> usize {
        self.count
    }
    pub const fn is_empty(&self) -> bool {
        self.count == 0
    }
    pub fn as_slice(&self) -> &[f64] {
        &self.values[..self.count]
    }
    fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.values[..self.count]
    }
}
impl core::ops::Deref for IonicSpeciesValues {
    type Target = [f64];
    fn deref(&self) -> &Self::Target {
        self.as_slice()
    }
}
impl AsRef<[f64]> for IonicSpeciesValues {
    fn as_ref(&self) -> &[f64] {
        self.as_slice()
    }
}
impl IntoIterator for IonicSpeciesValues {
    type Item = f64;
    type IntoIter = std::iter::Take<std::array::IntoIter<f64, CAPACITY>>;
    fn into_iter(self) -> Self::IntoIter {
        self.values.into_iter().take(self.count)
    }
}
impl<'a> IntoIterator for &'a IonicSpeciesValues {
    type Item = &'a f64;
    type IntoIter = std::slice::Iter<'a, f64>;
    fn into_iter(self) -> Self::IntoIter {
        self.as_slice().iter()
    }
}
impl<const N: usize> PartialEq<[f64; N]> for IonicSpeciesValues {
    fn eq(&self, rhs: &[f64; N]) -> bool {
        self.as_slice() == rhs.as_slice()
    }
}
impl<const N: usize> PartialEq<IonicSpeciesValues> for [f64; N] {
    fn eq(&self, rhs: &IonicSpeciesValues) -> bool {
        self.as_slice() == rhs.as_slice()
    }
}

/// Endpoint accuracy and finite nested work, independently of temporal error.
#[derive(Clone, Copy, Debug)]
pub struct IonicStepControls {
    pub absolute_water_tolerance_m3: f64,
    pub absolute_amount_tolerance_mol: f64,
    /// Relative to each actual transfer, never its background inventory.
    pub relative_transfer_tolerance: f64,
    /// Absolute original-to-final total-potential allowance [J].
    pub free_energy_tolerance_j: f64,
    /// Per water solve and per nested ion solve, including Newton damping,
    /// analytic-response replays and refused candidates.
    pub maximum_root_evaluations: usize,
    /// COMPLETE requested step/window, including coarse and rejected attempts.
    pub maximum_potential_evaluations: usize,
    pub electrical: ElectrodiffusionControls,
}
impl Default for IonicStepControls {
    fn default() -> Self {
        Self {
            absolute_water_tolerance_m3: 1e-18,
            absolute_amount_tolerance_mol: 1e-18,
            relative_transfer_tolerance: 1e-10,
            free_energy_tolerance_j: 1e-10,
            maximum_root_evaluations: 128,
            maximum_potential_evaluations: 1_000_000,
            electrical: ElectrodiffusionControls::default(),
        }
    }
}
impl IonicStepControls {
    fn validate(self) -> Result<()> {
        let e = self.electrical;
        if !positive(self.absolute_water_tolerance_m3)
            || !positive(self.absolute_amount_tolerance_mol)
            || !nonnegative(self.relative_transfer_tolerance)
            || self.relative_transfer_tolerance >= 1e-3
            || !nonnegative(self.free_energy_tolerance_j)
            || !(4..=4096).contains(&self.maximum_root_evaluations)
            || self.maximum_potential_evaluations == 0
            || !positive(e.absolute_current_tolerance_mol_m2_s)
            || !positive(e.relative_current_tolerance)
            || e.relative_current_tolerance >= 1.0
            || !(1..=4096).contains(&e.maximum_potential_evaluations)
        {
            return Err(invalid(
                "finite ionic transport needs explicit positive accuracy and bounded work",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct IonicStepReport {
    pub interval_s: [f64; 2],
    pub solvent_transfer_m3: f64,
    /// Actual left-to-right transfer of EACH original ion [mol].
    pub solute_transfer_mol: IonicSpeciesValues,
    pub endpoint: ElectrodiffusionReport,
    pub pressure_work_into_solution_j: f64,
    /// Mobile mixing + trapped-particle mixing + original mechanical potential.
    pub total_potential_change_j: f64,
    pub water_balance_residual_m3: f64,
    pub solute_balance_residual_mol: IonicSpeciesValues,
    pub water_root_residual_m3: f64,
    pub solute_root_residual_mol: IonicSpeciesValues,
    pub water_evaluations: usize,
    /// Original ion-state/response evaluations, including multicomponent work.
    pub salt_evaluations: usize,
    pub potential_evaluations: usize,
}

#[derive(Clone, Debug)]
pub struct IonicOsmoticPair {
    compartments: [OsmoticCompartment; 2],
    candidates: [AqueousSolution; 2],
    membrane: ElectrodiffusiveMembrane,
    fixed_charge_mol: [f64; 2],
    trapped_particles_mol: [f64; 2],
    area_m2: f64,
    maximum_charge_residual_mol: f64,
    totals: IonicSpeciesValues,
    stoichiometry: [f64; 2],
    multicomponent: Option<Box<multicomponent::Workspace>>,
    endpoint: ElectrodiffusionReport,
    time_s: f64,
    workspace_bytes: usize,
}
#[derive(Clone, Copy, Default)]
struct Work {
    water: usize,
    salt: usize,
    potential: usize,
}
#[derive(Clone, Copy)]
struct Probe {
    dt: f64,
    water: f64,
    volumes: [f64; 2],
    pressure_drop: f64,
}
#[derive(Clone, Copy)]
struct Trial {
    water: f64,
    salt: f64,
    volumes: [f64; 2],
    amounts: [IonicSpeciesValues; 2],
    flux: IonicSpeciesValues,
    transport: ElectrodiffusionReport,
    salt_residual: f64,
    water_residual: f64,
    ion_residual: IonicSpeciesValues,
    water_error: f64,
    ion_error: f64,
}

impl IonicOsmoticPair {
    /// Preserve supplied finite inventories and evaluate their initial flux.
    /// Both cells must use the same original ideal/virial law and solvent, with
    /// up to 64 explicitly ordered species; zero valence denotes a neutral
    /// co-solute. The absolute charge budget checks mobile PLUS fixed charge,
    /// without inventing counterions. Two opposite ions retain the scalar path.
    pub fn new(
        compartments: [OsmoticCompartment; 2],
        membrane: ElectrodiffusiveMembrane,
        fixed_charge_mol: [f64; 2],
        trapped_particles_mol: [f64; 2],
        area_m2: f64,
        maximum_charge_residual_mol: f64,
        maximum_workspace_bytes: usize,
    ) -> Result<Self> {
        Self::new_with_electrical_controls(
            compartments,
            membrane,
            fixed_charge_mol,
            trapped_particles_mol,
            area_m2,
            maximum_charge_residual_mol,
            ElectrodiffusionControls::default(),
            maximum_workspace_bytes,
        )
    }

    /// Explicit initial electrical accuracy for the supplied permeance scale.
    /// Later steps retain their own caller-selected electrical controls.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_electrical_controls(
        compartments: [OsmoticCompartment; 2],
        membrane: ElectrodiffusiveMembrane,
        fixed_charge_mol: [f64; 2],
        trapped_particles_mol: [f64; 2],
        area_m2: f64,
        maximum_charge_residual_mol: f64,
        initial_electrical: ElectrodiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self> {
        let z = membrane.valences();
        let reference = compartments[0].solution();
        let count = z.len();
        if count > CAPACITY
            || !positive(area_m2)
            || !positive(maximum_charge_residual_mol)
            || fixed_charge_mol.iter().any(|q| !q.is_finite())
            || trapped_particles_mol.iter().any(|&n| !nonnegative(n))
            || compartments.iter().any(|cell| {
                let s = cell.solution();
                cell.has_phase_inventory()
                    || s.solute_moles().len() != count
                    || s.temperature_k().to_bits() != reference.temperature_k().to_bits()
                    || s.solvent_molar_volume_m3_mol().to_bits()
                        != reference.solvent_molar_volume_m3_mol().to_bits()
                    || s.nonideal_model() != reference.nonideal_model()
            })
        {
            return Err(invalid(
                "finite ionic pairs require bounded original species, common solution physics and explicit charge budgets",
            ));
        }
        for side in 0..2 {
            admit_charge(
                compartments[side].solution(),
                z,
                fixed_charge_mol[side],
                side,
                maximum_charge_residual_mol,
            )?;
        }
        let binary =
            count == 2 && z[0] != 0 && z[1] != 0 && z[0].is_positive() != z[1].is_positive();
        let stoichiometry = if binary {
            let magnitudes = [
                i64::from(z[0]).unsigned_abs(),
                i64::from(z[1]).unsigned_abs(),
            ];
            let (mut a, mut b) = (magnitudes[0], magnitudes[1]);
            while b != 0 {
                (a, b) = (b, a % b);
            }
            [(magnitudes[1] / a) as f64, (magnitudes[0] / a) as f64]
        } else {
            [0.0; 2]
        };
        let totals = IonicSpeciesValues::from_fn(count, |i| {
            compartments[0].solution().solute_moles()[i]
                + compartments[1].solution().solute_moles()[i]
        });
        if totals.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite ionic totals overflow"));
        }
        finite(
            compartments[0].solution().solvent_volume_m3()
                + compartments[1].solution().solvent_volume_m3(),
        )?;
        let mut workspace_bytes = core::mem::size_of::<Self>()
            + membrane.heap_bytes()
            + 2 * compartments
                .iter()
                .map(|c| c.solution().heap_bytes())
                .sum::<usize>();
        if workspace_bytes > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: workspace_bytes,
                budget_bytes: maximum_workspace_bytes,
            }
            .into());
        }
        let multicomponent = if binary {
            None
        } else {
            let prepared = multicomponent::Workspace::new(
                &membrane,
                &totals,
                workspace_bytes,
                maximum_workspace_bytes,
            )?;
            workspace_bytes += prepared.storage_bytes();
            Some(Box::new(prepared))
        };
        let pressure_drop = compartments[0].pressure.pressure_difference_pa(
            compartments[0].solution().solvent_volume_m3(),
            compartments[1].pressure,
            compartments[1].solution().solvent_volume_m3(),
        )?;
        let endpoint = membrane.flux_into(
            ElectrodiffusionProblem {
                left: compartments[0].solution(),
                right: compartments[1].solution(),
                pressure_drop_pa: pressure_drop,
                trapped_particles_mol,
            },
            initial_electrical,
            &mut [0.0; CAPACITY][..count],
            || false,
        )?;
        Ok(Self {
            candidates: compartments.each_ref().map(|c| c.solution().clone()),
            compartments,
            membrane,
            fixed_charge_mol,
            trapped_particles_mol,
            area_m2,
            maximum_charge_residual_mol,
            totals,
            stoichiometry,
            multicomponent,
            endpoint,
            time_s: 0.0,
            workspace_bytes,
        })
    }
    pub const fn compartments(&self) -> &[OsmoticCompartment; 2] {
        &self.compartments
    }
    pub const fn endpoint(&self) -> &ElectrodiffusionReport {
        &self.endpoint
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn mobile_totals_mol(&self) -> &[f64] {
        self.totals.as_slice()
    }
    pub const fn species_count(&self) -> usize {
        self.totals.len()
    }
    pub const fn fixed_charge_mol(&self) -> [f64; 2] {
        self.fixed_charge_mol
    }
    pub const fn trapped_particles_mol(&self) -> [f64; 2] {
        self.trapped_particles_mol
    }

    fn salt_bounds(&self, volumes: [f64; 2]) -> Result<[f64; 2]> {
        let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
        for i in 0..2 {
            let (l, r, nu) = (
                self.compartments[0].solution().solute_moles()[i],
                self.compartments[1].solution().solute_moles()[i],
                self.stoichiometry[i],
            );
            lo = lo.max(-r / nu);
            hi = hi.min(l / nu);
            if let Some(model) = self.compartments[0].solution().nonideal_model() {
                let cap = model.maximum_concentrations_mol_m3()[i];
                lo = lo.max((l - cap * volumes[0]) / nu);
                hi = hi.min((cap * volumes[1] - r) / nu);
            }
        }
        // Adjacent representable extents avoid rounded over-withdrawal at an
        // open finite-rate endpoint; no concentration or amount floor is added.
        if lo < hi {
            if lo != 0.0 {
                lo = lo.next_up();
            }
            if hi != 0.0 {
                hi = hi.next_down();
            }
        }
        if lo > hi || !lo.is_finite() || !hi.is_finite() {
            return Err(OsmoticNetworkError::Unbracketed.into());
        }
        Ok([lo, hi])
    }

    fn sample_salt(
        &mut self,
        probe: Probe,
        salt: f64,
        controls: IonicStepControls,
        work: &mut Work,
        start: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        let transfers = IonicSpeciesValues::from_fn(2, |i| self.stoichiometry[i] * salt);
        let mut trial =
            self.sample_transfers(probe, &transfers, controls, work, start, cancelled)?;
        let pivot = usize::from(self.stoichiometry[1] < self.stoichiometry[0]);
        trial.salt = salt;
        trial.salt_residual = salt
            - integrated(probe.dt, self.area_m2, trial.flux[pivot]) / self.stoichiometry[pivot];
        Ok(trial)
    }

    fn sample_transfers(
        &mut self,
        probe: Probe,
        transfers: &IonicSpeciesValues,
        controls: IonicStepControls,
        work: &mut Work,
        start: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        poll(cancelled)?;
        if work.salt - start >= controls.maximum_root_evaluations {
            return Err(OsmoticNetworkError::NonlinearConvergence.into());
        }
        bump(&mut work.salt)?;
        let count = self.species_count();
        let amounts = [0, 1].map(|side| {
            IonicSpeciesValues::from_fn(count, |i| {
                let old = self.compartments[side].solution().solute_moles()[i];
                let transfer = transfers[i];
                if side == 0 {
                    old - transfer
                } else {
                    old + transfer
                }
            })
        });
        for side in 0..2 {
            self.candidates[side].replace_inventories(probe.volumes[side], &amounts[side])?;
        }
        let remaining = controls
            .maximum_potential_evaluations
            .checked_sub(work.potential)
            .filter(|n| *n > 0)
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        let mut electrical = controls.electrical;
        electrical.maximum_potential_evaluations =
            electrical.maximum_potential_evaluations.min(remaining);
        let mut flux = IonicSpeciesValues::zeroed(count);
        let result = self.membrane.flux_into_counted(
            ElectrodiffusionProblem {
                left: &self.candidates[0],
                right: &self.candidates[1],
                pressure_drop_pa: probe.pressure_drop,
                trapped_particles_mol: self.trapped_particles_mol,
            },
            electrical,
            flux.as_mut_slice(),
            &mut work.potential,
            &mut *cancelled,
        );
        let transport = match result {
            Err(ElectrodiffusionError::Cancelled) => {
                return Err(OsmoticNetworkError::Cancelled.into());
            }
            Err(_) if work.potential >= controls.maximum_potential_evaluations => {
                return Err(OsmoticNetworkError::WorkBudget.into());
            }
            other => other?,
        };
        let (mut ion_error, mut ion_residual) = (0.0_f64, IonicSpeciesValues::zeroed(count));
        for i in 0..count {
            let predicted = integrated(probe.dt, self.area_m2, flux[i]);
            let allowance = finite(
                controls.absolute_amount_tolerance_mol
                    + controls.relative_transfer_tolerance * transfers[i].abs(),
            )?;
            for side in 0..2 {
                let actual = (amounts[side][i]
                    - self.compartments[side].solution().solute_moles()[i])
                    * if side == 0 { -1.0 } else { 1.0 };
                let residual = actual - predicted;
                if residual.is_nan() {
                    return Err(invalid("unrepresentable original ionic endpoint equation"));
                }
                ion_error = ion_error.max(residual.abs() / allowance);
                if side == 0 {
                    ion_residual.values[i] = residual;
                }
            }
        }
        let predicted = integrated(probe.dt, self.area_m2, transport.volume_flux_m_s);
        let allowance = finite(
            controls.absolute_water_tolerance_m3
                + controls.relative_transfer_tolerance * probe.water.abs(),
        )?;
        let mut water_error = 0.0_f64;
        for side in 0..2 {
            let actual = (probe.volumes[side]
                - self.compartments[side].solution().solvent_volume_m3())
                * if side == 0 { -1.0 } else { 1.0 };
            if (actual - predicted).is_nan() {
                return Err(invalid("unrepresentable original water endpoint equation"));
            }
            water_error = water_error.max((actual - predicted).abs() / allowance);
        }
        Ok(Trial {
            water: probe.water,
            salt: 0.0,
            volumes: probe.volumes,
            amounts,
            flux,
            transport,
            salt_residual: 0.0,
            water_residual: probe.water - predicted,
            ion_residual,
            water_error,
            ion_error,
        })
    }

    fn solve_salt(
        &mut self,
        probe: Probe,
        controls: IonicStepControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        let start = work.salt;
        let bounds = self.salt_bounds(probe.volumes)?;
        let seed = 0.0_f64.clamp(bounds[0], bounds[1]);
        let initial = self.sample_salt(probe, seed, controls, work, start, cancelled)?;
        if initial.ion_error <= 1.0 {
            return Ok(initial);
        }
        let mut lo = if initial.salt_residual < 0.0 {
            initial
        } else {
            self.sample_salt(probe, bounds[0], controls, work, start, cancelled)?
        };
        if lo.ion_error <= 1.0 {
            return Ok(lo);
        }
        let mut hi = if initial.salt_residual > 0.0 {
            initial
        } else {
            self.sample_salt(probe, bounds[1], controls, work, start, cancelled)?
        };
        if hi.ion_error <= 1.0 {
            return Ok(hi);
        }
        if lo.salt_residual >= 0.0 || hi.salt_residual <= 0.0 || lo.salt >= hi.salt {
            return Err(OsmoticNetworkError::Unbracketed.into());
        }
        let mut midpoint = false;
        loop {
            let q = interpolate(
                lo.salt,
                hi.salt,
                lo.salt_residual,
                hi.salt_residual,
                midpoint,
            );
            if q <= lo.salt || q >= hi.salt {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            let trial = self.sample_salt(probe, q, controls, work, start, cancelled)?;
            if trial.ion_error <= 1.0 {
                return Ok(trial);
            }
            if trial.salt_residual < 0.0 {
                lo = trial;
            } else {
                hi = trial;
            }
            midpoint = !midpoint;
        }
    }

    fn sample_water(
        &mut self,
        dt: f64,
        water: f64,
        controls: IonicStepControls,
        work: &mut Work,
        start: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        poll(cancelled)?;
        if work.water - start >= controls.maximum_root_evaluations {
            return Err(OsmoticNetworkError::NonlinearConvergence.into());
        }
        bump(&mut work.water)?;
        let volumes = [
            self.compartments[0].solution().solvent_volume_m3() - water,
            self.compartments[1].solution().solvent_volume_m3() + water,
        ];
        if volumes
            .iter()
            .zip(&self.compartments)
            .any(|(&v, c)| !positive(v) || v < c.minimum_volume_m3)
        {
            return Err(OsmoticNetworkError::Unbracketed.into());
        }
        let pressure_drop = self.compartments[0].pressure.pressure_difference_pa(
            volumes[0],
            self.compartments[1].pressure,
            volumes[1],
        )?;
        let probe = Probe {
            dt,
            water,
            volumes,
            pressure_drop,
        };
        if let Some(mut prepared) = self.multicomponent.take() {
            // Numerical scratch has no physical history. Always restore its
            // owner, including after cancellation or a refused Newton trial.
            let result = prepared.solve(self, probe, controls, work, cancelled);
            self.multicomponent = Some(prepared);
            result
        } else {
            self.solve_salt(probe, controls, work, cancelled)
        }
    }

    fn solve_water(
        &mut self,
        dt: f64,
        controls: IonicStepControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        let start = work.water;
        let bounds = [
            -(self.compartments[1].solution().solvent_volume_m3()
                - self.compartments[1].minimum_volume_m3)
                .next_down(),
            (self.compartments[0].solution().solvent_volume_m3()
                - self.compartments[0].minimum_volume_m3)
                .next_down(),
        ];
        // A capped solution may require simultaneous dilution to admit its
        // finite salt transfer. Zero water is then NOT a valid nested seed,
        // even though the original state and coupled endpoint are supported.
        let initial = match self.sample_water(dt, 0.0, controls, work, start, cancelled) {
            Ok(trial) => trial,
            Err(error) if fatal(&error) => return Err(error),
            Err(error) => {
                let predictor = integrated(dt, self.area_m2, self.endpoint.volume_flux_m_s);
                let preferred = if predictor < 0.0 { -1.0 } else { 1.0 };
                let extent = bounds[1].max(-bounds[0]);
                let base = if predictor.is_finite() && predictor != 0.0 {
                    predictor.abs().min(0.875 * extent)
                } else {
                    0.125 * extent
                };
                let (mut small, mut large) = (base, base);
                let mut extreme = [f64::NAN; 2];
                let mut supported = None;
                'search: for k in 0..controls.maximum_root_evaluations {
                    let distance = if k % 2 == 0 { small } else { large };
                    for direction in [preferred, -preferred] {
                        if work.water - start >= controls.maximum_root_evaluations / 2 {
                            break 'search;
                        }
                        let side = usize::from(direction > 0.0);
                        let available = if side == 1 { bounds[1] } else { -bounds[0] };
                        let water = direction * distance.min(0.875 * available);
                        if water == 0.0 || water == extreme[side] {
                            continue;
                        }
                        if k % 2 != 0 {
                            extreme[side] = water;
                        }
                        match self.sample_water(dt, water, controls, work, start, cancelled) {
                            Ok(trial) => {
                                supported = Some(trial);
                                break 'search;
                            }
                            Err(error) if fatal(&error) => return Err(error),
                            Err(_) => {}
                        }
                    }
                    if k % 2 == 0 {
                        small *= 0.5;
                    } else {
                        large = (2.0 * large).min(extent);
                    }
                }
                supported.ok_or(error)?
            }
        };
        if initial.water_error <= 1.0 {
            return Ok(initial);
        }
        let direction = if initial.water_residual < 0.0 {
            1.0
        } else {
            -1.0
        };
        let available = if direction > 0.0 {
            bounds[1] - initial.water
        } else {
            initial.water - bounds[0]
        };
        let mut distance = if initial.water_residual.abs() < available {
            initial.water_residual.abs()
        } else {
            0.125 * available
        };
        let mut previous = initial;
        let mut refused = None;
        let (mut lo, mut hi) = loop {
            let previous_distance = direction * (previous.water - initial.water);
            if distance <= previous_distance
                || distance > available
                || work.water - start >= controls.maximum_root_evaluations
            {
                return Err(OsmoticNetworkError::Unbracketed.into());
            }
            match self.sample_water(
                dt,
                initial.water + direction * distance,
                controls,
                work,
                start,
                cancelled,
            ) {
                Ok(trial) => {
                    if trial.water_error <= 1.0 {
                        return Ok(trial);
                    }
                    if direction * trial.water_residual >= 0.0 {
                        break if direction > 0.0 {
                            (previous, trial)
                        } else {
                            (trial, previous)
                        };
                    }
                    previous = trial;
                    distance = refused.map_or((2.0 * distance).min(available), |upper| {
                        distance + 0.5 * (upper - distance)
                    });
                }
                Err(error) if fatal(&error) => return Err(error),
                Err(_) => {
                    refused = Some(distance);
                    distance = previous_distance + 0.5 * (distance - previous_distance);
                }
            }
        };
        let mut midpoint = false;
        loop {
            let water = interpolate(
                lo.water,
                hi.water,
                lo.water_residual,
                hi.water_residual,
                midpoint,
            );
            if water <= lo.water || water >= hi.water {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            let trial = self.sample_water(dt, water, controls, work, start, cancelled)?;
            if trial.water_error <= 1.0 {
                return Ok(trial);
            }
            if trial.water_residual < 0.0 {
                lo = trial;
            } else {
                hi = trial;
            }
            midpoint = !midpoint;
        }
    }

    pub fn step(
        &mut self,
        duration_s: f64,
        controls: IonicStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonicStepReport> {
        self.step_counted(duration_s, controls, &mut Work::default(), &mut cancelled)
    }
    fn step_counted(
        &mut self,
        duration_s: f64,
        controls: IonicStepControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IonicStepReport> {
        poll(cancelled)?;
        controls.validate()?;
        if !positive(duration_s) {
            return Err(invalid("positive finite ionic step required"));
        }
        let end = finite(self.time_s + duration_s)?;
        if end <= self.time_s {
            return Err(OsmoticNetworkError::StepUnderflow.into());
        }
        let start = *work;
        // The represented clock interval owns the integrated fluxes. At a
        // large clock, adding the requested duration can round to another dt.
        let trial = self.solve_water(end - self.time_s, controls, work, cancelled)?;
        for side in 0..2 {
            self.candidates[side].replace_inventories(trial.volumes[side], &trial.amounts[side])?;
        }
        let before = self.compartments.each_ref().map(|c| c.solution());
        let after = self.candidates.each_ref();
        let (water_balance, solute_balance) = self.balance(before, after, controls)?;
        let (pressure_work, potential_change) = self.potential_change(before, after)?;
        if potential_change > controls.free_energy_tolerance_j {
            return Err(OsmoticNetworkError::CoupledStepFreeEnergy {
                change_j: potential_change,
                maximum_j: controls.free_energy_tolerance_j,
            }
            .into());
        }
        let report = IonicStepReport {
            interval_s: [self.time_s, end],
            solvent_transfer_m3: before[0].solvent_volume_m3() - after[0].solvent_volume_m3(),
            solute_transfer_mol: IonicSpeciesValues::from_fn(self.species_count(), |i| {
                before[0].solute_moles()[i] - after[0].solute_moles()[i]
            }),
            endpoint: trial.transport,
            pressure_work_into_solution_j: pressure_work,
            total_potential_change_j: potential_change,
            water_balance_residual_m3: water_balance,
            solute_balance_residual_mol: solute_balance,
            water_root_residual_m3: before[0].solvent_volume_m3()
                - after[0].solvent_volume_m3()
                - integrated(
                    end - self.time_s,
                    self.area_m2,
                    trial.transport.volume_flux_m_s,
                ),
            solute_root_residual_mol: trial.ion_residual,
            water_evaluations: work.water - start.water,
            salt_evaluations: work.salt - start.salt,
            potential_evaluations: work.potential - start.potential,
        };
        poll(cancelled)?;
        for side in 0..2 {
            core::mem::swap(
                &mut self.compartments[side].solution,
                &mut self.candidates[side],
            );
        }
        self.endpoint = trial.transport;
        self.time_s = end;
        Ok(report)
    }

    fn balance(
        &self,
        before: [&AqueousSolution; 2],
        after: [&AqueousSolution; 2],
        controls: IonicStepControls,
    ) -> Result<(f64, IonicSpeciesValues)> {
        let water = finite(
            (after[0].solvent_volume_m3() - before[0].solvent_volume_m3())
                + (after[1].solvent_volume_m3() - before[1].solvent_volume_m3()),
        )?;
        // One original-inventory rounding allowance, never a per-substep sum.
        let maximum = controls.absolute_water_tolerance_m3
            + 64.0 * f64::EPSILON * (before[0].solvent_volume_m3() + before[1].solvent_volume_m3());
        if water.abs() > maximum {
            return Err(OsmoticNetworkError::CoupledWaterBalance {
                residual_m3: water,
                maximum_m3: maximum,
            }
            .into());
        }
        let mut solute = IonicSpeciesValues::zeroed(self.species_count());
        for i in 0..self.species_count() {
            solute.values[i] = finite(
                (after[0].solute_moles()[i] - before[0].solute_moles()[i])
                    + (after[1].solute_moles()[i] - before[1].solute_moles()[i]),
            )?;
            if solute[i].abs()
                > controls.absolute_amount_tolerance_mol + 64.0 * f64::EPSILON * self.totals[i]
            {
                return Err(OsmoticNetworkError::DiffusionBalance {
                    species: i,
                    residual_mol: solute[i],
                }
                .into());
            }
        }
        for side in 0..2 {
            admit_charge(
                after[side],
                self.membrane.valences(),
                self.fixed_charge_mol[side],
                side,
                self.maximum_charge_residual_mol,
            )?;
        }
        Ok((water, solute))
    }

    fn potential_change(
        &self,
        before: [&AqueousSolution; 2],
        after: [&AqueousSolution; 2],
    ) -> Result<(f64, f64)> {
        let (mut mechanical, mut mc, mut total, mut tc) = (0.0, 0.0, 0.0, 0.0);
        let reference = self.compartments[0].pressure.reference_pressure_pa;
        for side in 0..2 {
            let (a, b) = (before[side], after[side]);
            let mixing = if let Some(model) = a.nonideal_model() {
                model
                    .reference_cancelled_inventory_free_energy_change_j(
                        a.solvent_volume_m3(),
                        a.solute_moles(),
                        b.solvent_volume_m3(),
                        b.solute_moles(),
                    )
                    .map_err(OsmosisError::from)?
            } else {
                ideal_change(a, b)?
            };
            let trapped = R
                * a.temperature_k()
                * self.trapped_particles_mol[side]
                * log_ratio(a.solvent_volume_m3(), b.solvent_volume_m3());
            let change = self.compartments[side]
                .pressure
                .reference_cancelled_potential_change_j(
                    a.solvent_volume_m3(),
                    b.solvent_volume_m3(),
                    reference,
                )?;
            compensated_add(&mut mechanical, &mut mc, change);
            for value in [mixing, trapped, change] {
                compensated_add(&mut total, &mut tc, value);
            }
        }
        Ok((-finite(mechanical)?, finite(total)?))
    }
}

fn admit_charge(
    solution: &AqueousSolution,
    z: &[i32],
    fixed: f64,
    side: usize,
    maximum: f64,
) -> Result<()> {
    let (mut charge, mut correction) = (fixed, 0.0);
    for (&valence, &amount) in z.iter().zip(solution.solute_moles()) {
        compensated_add(&mut charge, &mut correction, f64::from(valence) * amount);
    }
    if !charge.is_finite() || charge.abs() > maximum {
        return Err(IonicOsmoticError::Charge {
            side,
            residual_mol: charge,
            maximum_mol: maximum,
        });
    }
    Ok(())
}
fn ideal_change(a: &AqueousSolution, b: &AqueousSolution) -> Result<f64> {
    let (mut sum, mut correction) = (0.0, 0.0);
    let volumes = log_ratio(a.solvent_volume_m3(), b.solvent_volume_m3());
    for (&old, &new) in a.solute_moles().iter().zip(b.solute_moles()) {
        let change = if old == 0.0 && new == 0.0 {
            0.0
        } else if old == 0.0 {
            new * (det::ln(new / b.solvent_volume_m3()) - 1.0)
        } else if new == 0.0 {
            -old * (det::ln(old / a.solvent_volume_m3()) - 1.0)
        } else {
            (new - old) * (det::ln(old / a.solvent_volume_m3()) - 1.0)
                + new * (log_ratio(new, old) + volumes)
        };
        compensated_add(&mut sum, &mut correction, finite(change)?);
    }
    finite(R * a.temperature_k() * sum)
}
fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let relative = (a - b) / b;
    if (-0.5..=1.0).contains(&relative) {
        let s = relative / (2.0 + relative);
        let (mut term, mut result) = (s, s);
        for k in 1..24 {
            term *= s * s;
            result += term / f64::from(2 * k + 1);
        }
        2.0 * result
    } else {
        det::ln(a) - det::ln(b)
    }
}
fn interpolate(lo: f64, hi: f64, flo: f64, fhi: f64, midpoint: bool) -> f64 {
    let width = hi - lo;
    let scale = flo.abs().max(fhi.abs());
    let proposal = lo + width * (-flo / scale) / (fhi / scale - flo / scale);
    // Forced alternating bisection already bounds bracket shrinkage. Keeping
    // any strictly interior secant avoids resolving an entire bath-sized
    // bracket by bisection when a small time step transfers only trace salt.
    if !midpoint && proposal.is_finite() && proposal > lo && proposal < hi {
        proposal
    } else {
        lo + 0.5 * width
    }
}
fn integrated(dt: f64, area: f64, flux: f64) -> f64 {
    if flux == 0.0 {
        return 0.0;
    }
    let mut factors = [dt, area, flux.abs()];
    factors.sort_by(f64::total_cmp);
    ((factors[0] * factors[2]) * factors[1]).copysign(flux)
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn nonnegative(v: f64) -> bool {
    v.is_finite() && v >= 0.0
}
fn invalid(reason: &'static str) -> IonicOsmoticError {
    OsmoticNetworkError::Invalid(reason).into()
}
fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("nonfinite ionic transport arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(OsmoticNetworkError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn bump(count: &mut usize) -> Result<()> {
    *count = count
        .checked_add(1)
        .ok_or(OsmoticNetworkError::WorkBudget)?;
    Ok(())
}
fn fatal(error: &IonicOsmoticError) -> bool {
    matches!(
        error,
        IonicOsmoticError::Network(
            OsmoticNetworkError::Cancelled | OsmoticNetworkError::WorkBudget
        )
    )
}

#[cfg(test)]
mod tests;
