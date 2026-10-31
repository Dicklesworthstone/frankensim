//! Finite water transport with fast, charge-conserving Donnan ion equilibrium.
//!
//! Mobile ions re-equilibrate at each trial pair of solvent volumes, retaining
//! every species total and both fixed charge amounts. Water responds to the
//! FULL ionic plus explicitly trapped osmotic pressure and the original
//! mechanical pressure/volume laws. There are no independent ionic reflection
//! coefficients or unbalanced ionic leaks. This is an ideal dilute fast-ion
//! limit, not transient Nernst--Planck transport or a nonideal electrolyte law.
//!
//! Each step is backward Euler with a safeguarded scalar water solve. Root
//! tolerances control the accepted balance/transfer, not temporal truncation
//! error; callers refine the step to resolve the transient. All mutable
//! workspace is prepared once. Failure or cancellation publishes neither
//! inventories, equilibrium diagnostics, nor time.

use super::{OsmoticCompartment, OsmoticNetworkError};
use fs_material::osmosis::donnan::{
    DonnanControls, DonnanError, DonnanProblem, DonnanReport, IdealDonnan, MAX_DONNAN_SPECIES,
};
use fs_material::osmosis::{HydraulicPermeance, OsmosisError};

/// Original water/ion admission or a failed conservation check.
#[derive(Debug)]
pub enum DonnanWaterError {
    Network(OsmoticNetworkError),
    Ions(DonnanError),
    InventoryBalance { species: usize, residual_mol: f64 },
}
impl core::fmt::Display for DonnanWaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Donnan water transport: {self:?}")
    }
}
impl std::error::Error for DonnanWaterError {}
impl From<OsmoticNetworkError> for DonnanWaterError {
    fn from(value: OsmoticNetworkError) -> Self {
        Self::Network(value)
    }
}
impl From<OsmosisError> for DonnanWaterError {
    fn from(value: OsmosisError) -> Self {
        Self::Network(value.into())
    }
}
impl From<DonnanError> for DonnanWaterError {
    fn from(value: DonnanError) -> Self {
        if value == DonnanError::Cancelled {
            Self::Network(OsmoticNetworkError::Cancelled)
        } else {
            Self::Ions(value)
        }
    }
}
type Result<T> = core::result::Result<T, DonnanWaterError>;

/// Explicit scalar solve accuracy and finite nested work limits.
#[derive(Clone, Copy, Debug)]
pub struct DonnanStepControls {
    /// Absolute water-transfer/root tolerance [m³].
    pub absolute_water_tolerance_m3: f64,
    /// Relative water-transfer/root tolerance; not a time-error estimate.
    pub relative_water_tolerance: f64,
    /// Includes bracketing, interpolation, bisection and accepted re-evaluation.
    pub maximum_root_evaluations: usize,
    /// Bounded charge solve performed at each candidate water transfer.
    pub ion_equilibrium: DonnanControls,
}
impl Default for DonnanStepControls {
    fn default() -> Self {
        Self {
            absolute_water_tolerance_m3: 1e-18,
            relative_water_tolerance: 1e-10,
            maximum_root_evaluations: 128,
            ion_equilibrium: DonnanControls::default(),
        }
    }
}

/// Accepted physical endpoint, conservation residuals and numerical work counts.
#[derive(Clone, Debug)]
pub struct DonnanStepReport {
    pub interval_s: [f64; 2],
    /// Positive transfer goes from left compartment to right [m³].
    pub solvent_transfer_m3: f64,
    /// Interval-mean solvent flux per original membrane area [m/s].
    pub volume_flux_m_s: f64,
    /// Final electrochemical state and FULL mobile-plus-trapped osmotic pressures.
    pub equilibrium: DonnanReport,
    /// Mechanical boundary work into the solution [J]. This is not a heat balance.
    pub pressure_work_into_solution_j: f64,
    pub water_balance_residual_m3: f64,
    pub solute_balance_residual_mol: [f64; MAX_DONNAN_SPECIES],
    pub species_count: usize,
    /// Actual backward-Euler residual, also reported after bracket-position acceptance.
    pub water_root_residual_m3: f64,
    /// All attempted water-residual evaluations, including refused candidates.
    pub root_evaluations: usize,
    /// Charge-residual evaluations from completed ion equilibria. Refused ion
    /// trials do not return a partial count; each is separately bounded by
    /// `ion_equilibrium.maximum_iterations`, with at most `root_evaluations`
    /// total inner solve attempts.
    pub ion_root_evaluations: usize,
}

#[derive(Clone, Debug)]
pub struct DonnanOsmoticPair {
    compartments: [OsmoticCompartment; 2],
    candidates: [OsmoticCompartment; 2],
    ions: IdealDonnan,
    fixed_charge_mol: [f64; 2],
    trapped_particles_mol: [f64; 2],
    membrane_area_m2: f64,
    hydraulic: HydraulicPermeance,
    mobile_totals_mol: Box<[f64]>,
    trial_left: Box<[f64]>,
    trial_right: Box<[f64]>,
    equilibrium: DonnanReport,
    workspace_bytes: usize,
    time_s: f64,
}

#[derive(Clone, Copy)]
struct Trial {
    water: f64,
    volumes: [f64; 2],
    residual: f64,
    derivative: Option<f64>,
    equilibrium: DonnanReport,
}
#[derive(Default)]
struct Work {
    water: usize,
    ions: usize,
}

impl DonnanOsmoticPair {
    /// Prepare a closed pair and canonicalize the initial mobile ion partition
    /// at its supplied volumes. Input species totals are retained exactly as
    /// the conserved basis; no additional ions or solvent are supplied.
    /// Initial fast ionic relaxation has no simulated duration or heat report.
    pub fn new(
        mut compartments: [OsmoticCompartment; 2],
        ions: IdealDonnan,
        fixed_charge_mol: [f64; 2],
        trapped_particles_mol: [f64; 2],
        membrane_area_m2: f64,
        hydraulic: HydraulicPermeance,
        maximum_workspace_bytes: usize,
    ) -> Result<Self> {
        let count = ions.valences().len();
        let reference = compartments[0].solution();
        if !positive(membrane_area_m2)
            || compartments.iter().any(|cell| {
                let s = cell.solution();
                !s.is_ideal()
                    || s.solute_moles().len() != count
                    || s.temperature_k().to_bits() != reference.temperature_k().to_bits()
                    || s.solvent_molar_volume_m3_mol().to_bits()
                        != reference.solvent_molar_volume_m3_mol().to_bits()
            })
        {
            return Err(invalid(
                "Donnan water transport requires common ideal ionic layout, solvent, temperature and positive area",
            ));
        }
        let workspace_bytes = core::mem::size_of::<Self>()
            + 2 * compartments
                .iter()
                .map(|cell| cell.solution().heap_bytes())
                .sum::<usize>()
            + ions.heap_bytes()
            + 3 * count * core::mem::size_of::<f64>();
        if workspace_bytes > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: workspace_bytes,
                budget_bytes: maximum_workspace_bytes,
            }
            .into());
        }
        let volumes = compartments
            .each_ref()
            .map(|cell| cell.solution().solvent_volume_m3());
        finite(
            volumes[0] + volumes[1],
            "total Donnan water volume is not representable",
        )?;
        let mut totals = vec![0.0; count].into_boxed_slice();
        for (i, total) in totals.iter_mut().enumerate() {
            *total = finite(
                compartments[0].solution().solute_moles()[i]
                    + compartments[1].solution().solute_moles()[i],
                "total mobile ion amount is not representable",
            )?;
        }
        let mut trial_left = vec![0.0; count].into_boxed_slice();
        let mut trial_right = vec![0.0; count].into_boxed_slice();
        let equilibrium = ions.equilibrate_into(
            DonnanProblem {
                temperature_k: reference.temperature_k(),
                volumes_m3: volumes,
                mobile_totals_mol: &totals,
                fixed_charge_mol,
                trapped_particles_mol,
            },
            DonnanControls::default(),
            &mut trial_left,
            &mut trial_right,
            || false,
        )?;
        compartments[0]
            .solution
            .replace_inventories(volumes[0], &trial_left)?;
        compartments[1]
            .solution
            .replace_inventories(volumes[1], &trial_right)?;
        Ok(Self {
            candidates: compartments.clone(),
            compartments,
            ions,
            fixed_charge_mol,
            trapped_particles_mol,
            membrane_area_m2,
            hydraulic,
            mobile_totals_mol: totals,
            trial_left,
            trial_right,
            equilibrium,
            workspace_bytes,
            time_s: 0.0,
        })
    }

    /// Accepted solvent and MOBILE ion inventories. These solution owners do
    /// not contain the separate trapped osmoles; their pressure getters are
    /// mobile-only. Use [`Self::equilibrium`] for the full osmotic pressures.
    pub const fn compartments(&self) -> &[OsmoticCompartment; 2] {
        &self.compartments
    }
    pub const fn equilibrium(&self) -> &DonnanReport {
        &self.equilibrium
    }
    pub fn mobile_totals_mol(&self) -> &[f64] {
        &self.mobile_totals_mol
    }
    pub const fn fixed_charge_mol(&self) -> [f64; 2] {
        self.fixed_charge_mol
    }
    pub const fn trapped_particles_mol(&self) -> [f64; 2] {
        self.trapped_particles_mol
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    fn evaluate(
        &mut self,
        mobility: f64,
        water: f64,
        controls: DonnanStepControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        poll(cancelled)?;
        if work.water >= controls.maximum_root_evaluations {
            return Err(OsmoticNetworkError::NonlinearConvergence.into());
        }
        work.water += 1;
        let volumes = [
            self.compartments[0].solution().solvent_volume_m3() - water,
            self.compartments[1].solution().solvent_volume_m3() + water,
        ];
        if volumes
            .iter()
            .zip(&self.compartments)
            .any(|(&v, cell)| !positive(v) || v < cell.minimum_volume_m3)
        {
            return Err(OsmoticNetworkError::Unbracketed.into());
        }
        let equilibrium = self.ions.equilibrate_into(
            DonnanProblem {
                temperature_k: self.compartments[0].solution().temperature_k(),
                volumes_m3: volumes,
                mobile_totals_mol: &self.mobile_totals_mol,
                fixed_charge_mol: self.fixed_charge_mol,
                trapped_particles_mol: self.trapped_particles_mol,
            },
            controls.ion_equilibrium,
            &mut self.trial_left,
            &mut self.trial_right,
            &mut *cancelled,
        )?;
        work.ions += equilibrium.iterations;
        let pressure = [
            self.compartments[0].pressure.pressure_pa(volumes[0])?,
            self.compartments[1].pressure.pressure_pa(volumes[1])?,
        ];
        let residual = finite(
            water - mobility * equilibrium.solvent_driving_pressure_pa(pressure)?,
            "nonfinite Donnan water residual",
        )?;
        let derivative = equilibrium
            .osmotic_pressure_drop_slope_pa_per_m3
            .and_then(|slope| {
                let value = 1.0
                    + mobility
                        * (self.compartments[0].pressure.elastance_pa_per_m3
                            + self.compartments[1].pressure.elastance_pa_per_m3
                            + slope);
                positive(value).then_some(value)
            });
        Ok(Trial {
            water,
            volumes,
            residual,
            derivative,
            equilibrium,
        })
    }

    fn solve_water(
        &mut self,
        mobility: f64,
        controls: DonnanStepControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        let tolerance = |water: f64| {
            controls.absolute_water_tolerance_m3 + controls.relative_water_tolerance * water.abs()
        };
        let initial = self.evaluate(mobility, 0.0, controls, work, cancelled)?;
        if initial.residual.abs() <= tolerance(0.0) {
            return Ok(initial);
        }
        let direction = if initial.residual < 0.0 { 1.0 } else { -1.0 };
        let donor = usize::from(direction < 0.0);
        let available = (self.compartments[donor].solution().solvent_volume_m3()
            - self.compartments[donor].minimum_volume_m3)
            .next_down();
        let estimate = initial.derivative.map_or(initial.residual.abs(), |slope| {
            initial.residual.abs() / slope
        });
        // The equilibrated ideal free energy makes this residual monotone.
        // Keep an admissible Newton predictor instead of truncating an exact
        // linear-pressure root; start conservatively only outside the domain.
        let mut distance = if estimate < available {
            estimate
        } else {
            0.125 * available
        };
        let mut previous = initial;
        let mut refused_bound = None;
        let (mut lo, mut hi) = loop {
            poll(cancelled)?;
            if work.water >= controls.maximum_root_evaluations {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            if distance <= previous.water.abs() || distance > available {
                return Err(OsmoticNetworkError::Unbracketed.into());
            }
            match self.evaluate(mobility, direction * distance, controls, work, cancelled) {
                Ok(trial) => {
                    if trial.residual.abs() <= tolerance(trial.water) {
                        return Ok(trial);
                    }
                    if direction * trial.residual >= 0.0 {
                        break if direction > 0.0 {
                            (previous, trial)
                        } else {
                            (trial, previous)
                        };
                    }
                    previous = trial;
                    distance = refused_bound.map_or_else(
                        || (2.0 * distance).min(available),
                        |upper| distance + 0.5 * (upper - distance),
                    );
                }
                Err(DonnanWaterError::Network(OsmoticNetworkError::Cancelled)) => {
                    return Err(OsmoticNetworkError::Cancelled.into());
                }
                Err(_) => {
                    refused_bound = Some(distance);
                    distance = previous.water.abs() + 0.5 * (distance - previous.water.abs());
                }
            }
        };
        let mut bisect_next = false;
        loop {
            let width = hi.water - lo.water;
            let midpoint = lo.water + 0.5 * width;
            if width <= 2.0 * tolerance(midpoint) {
                // Re-evaluate the accepted point so both trial ion buffers
                // describe precisely the water transfer being published.
                return self.evaluate(mobility, midpoint, controls, work, cancelled);
            }
            let best = if lo.residual.abs() < hi.residual.abs() {
                lo
            } else {
                hi
            };
            let proposed = best
                .derivative
                .map_or(f64::NAN, |slope| best.water - best.residual / slope);
            let interpolate =
                !bisect_next && proposed.is_finite() && proposed > lo.water && proposed < hi.water;
            let water = if interpolate { proposed } else { midpoint };
            if water == lo.water || water == hi.water {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            let trial = self.evaluate(mobility, water, controls, work, cancelled)?;
            if trial.residual.abs() <= tolerance(water) {
                return Ok(trial);
            }
            if trial.residual < 0.0 {
                lo = trial;
            } else {
                hi = trial;
            }
            bisect_next = interpolate && hi.water - lo.water > 0.5 * width;
        }
    }

    /// Advance one first-order fully implicit water step. The charge solve and
    /// full osmotic force are evaluated at the accepted finite inventories.
    /// A narrow sign-changing water bracket may satisfy position accuracy even
    /// when a stiff residual is limited by floating-point precision.
    pub fn step(
        &mut self,
        duration_s: f64,
        controls: DonnanStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DonnanStepReport> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || !positive(controls.absolute_water_tolerance_m3)
            || !positive(controls.relative_water_tolerance)
            || controls.relative_water_tolerance >= 1.0
            || !(2..=4096).contains(&controls.maximum_root_evaluations)
        {
            return Err(invalid(
                "invalid Donnan water step or scalar solve controls",
            ));
        }
        let end = finite(self.time_s + duration_s, "Donnan time is not representable")?;
        let area_time = duration_s * self.membrane_area_m2;
        let mobility = area_time * self.hydraulic.value();
        if end <= self.time_s
            || !positive(area_time)
            || !mobility.is_finite()
            || self.hydraulic.value() > 0.0 && mobility == 0.0
        {
            return Err(OsmoticNetworkError::StepUnderflow.into());
        }
        let mut work = Work::default();
        let trial = self.solve_water(mobility, controls, &mut work, &mut cancelled)?;
        self.candidates[0]
            .solution
            .replace_inventories(trial.volumes[0], &self.trial_left)?;
        self.candidates[1]
            .solution
            .replace_inventories(trial.volumes[1], &self.trial_right)?;
        let before_water = self
            .compartments
            .iter()
            .map(|cell| cell.solution().solvent_volume_m3())
            .sum::<f64>();
        let water_balance = finite(
            (trial.volumes[0] + trial.volumes[1]) - before_water,
            "nonfinite Donnan water balance",
        )?;
        if water_balance.abs()
            > controls.absolute_water_tolerance_m3 + 32.0 * f64::EPSILON * before_water
        {
            return Err(invalid("Donnan water inventory balance failed"));
        }
        let mut balances = [0.0; MAX_DONNAN_SPECIES];
        for (i, (&left, &right)) in self.trial_left.iter().zip(&self.trial_right).enumerate() {
            balances[i] = (left + right) - self.mobile_totals_mol[i];
            if !balances[i].is_finite()
                || balances[i].abs() > 32.0 * f64::EPSILON * self.mobile_totals_mol[i]
            {
                return Err(DonnanWaterError::InventoryBalance {
                    species: i,
                    residual_mol: balances[i],
                });
            }
        }
        let pressure_before =
            self.compartments[0].pressure_pa()? - self.compartments[1].pressure_pa()?;
        let pressure_after =
            self.candidates[0].pressure_pa()? - self.candidates[1].pressure_pa()?;
        let pressure_work = finite(
            (0.5 * pressure_before + 0.5 * pressure_after) * trial.water,
            "Donnan mechanical work is not representable",
        )?;
        let report = DonnanStepReport {
            interval_s: [self.time_s, end],
            solvent_transfer_m3: trial.water,
            volume_flux_m_s: finite(
                trial.water / area_time,
                "Donnan volume flux is not representable",
            )?,
            equilibrium: trial.equilibrium,
            pressure_work_into_solution_j: pressure_work,
            water_balance_residual_m3: water_balance,
            solute_balance_residual_mol: balances,
            species_count: self.mobile_totals_mol.len(),
            water_root_residual_m3: trial.residual,
            root_evaluations: work.water,
            ion_root_evaluations: work.ions,
        };
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.compartments, &mut self.candidates);
        self.equilibrium = trial.equilibrium;
        self.time_s = end;
        Ok(report)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(reason: &'static str) -> DonnanWaterError {
    OsmoticNetworkError::Invalid(reason).into()
}
fn finite(value: f64, reason: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(reason))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(OsmoticNetworkError::Cancelled.into())
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "donnan/tests.rs"]
mod tests;
