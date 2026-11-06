//! Prepared chemical equilibrium at prescribed volume and internal-energy change.
//!
//! This uses the original TV/Gibbs owner, source-relative NASA calorics and the
//! same safeguarded UV inversion as absolute-energy queries. The source and
//! elemental inventory are immutable; only numerical warm starts change.
//! Heat and mechanical work belong to the consuming balance, not to this L1
//! closure. No reaction rate, ignition or entropy-production law is inferred.

use super::super::expansion::path::Reference;
use super::*;

/// One accepted scalar query; source receipts are allocated only on publication.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumEnergyVolumeSampleV1 {
    /// Actual absolute temperature [K].
    pub temperature: Temperature,
    /// Actual chemical-equilibrium pressure [Pa].
    pub pressure: Pressure,
    /// Caller-declared volume [m3].
    pub volume: Volume,
    /// Original reacting derivatives at this same endpoint.
    pub response: EquilibriumThermodynamicResponseV1,
    /// Requested U-U_reference [J].
    pub target_energy_change_j: f64,
    /// Source-relative actual U-U_reference [J].
    pub internal_energy_change_j: f64,
    /// Actual minus requested internal-energy change [J].
    pub energy_residual_j: f64,
    /// Source-relative actual S-S_reference [J/K], including chemical mixing.
    pub entropy_change_j_k: f64,
    /// Actual (NRT/p-V)/V.
    pub relative_volume_residual: f64,
    /// Actual (product mass-feed mass)/feed mass.
    pub relative_mass_residual: f64,
    /// All chemical-equilibrium iterations, including bracket trials.
    pub composition_iterations: usize,
    /// Interior temperature trials, including a preferred predictor and excluding bracket endpoints.
    pub temperature_iterations: usize,
    /// All fixed-volume pressure iterations, including bracket trials.
    pub pressure_iterations: usize,
}

/// Reusable UV closure over one original chemical inventory and energy reference.
///
/// Failed/cancelled queries may change numerical guesses but invalidate the
/// publishable trial. The immutable reference is never rebased or accumulated.
/// Native species/element caps bound workspace storage. No allocation is made
/// for source receipts during trials; the existing full-KKT fallback can allocate.
pub struct EquilibriumEnergyVolumePathV1<'a> {
    owner: &'a IdealGasEquilibriumV1,
    work: Workspace,
    reference: Reference,
    reference_state: IdealGasIsochoricEquilibriumV1,
    controls: EquilibriumInternalEnergySolveV1,
    current: Option<EquilibriumEnergyVolumeSampleV1>,
}

impl<'a> EquilibriumEnergyVolumePathV1<'a> {
    /// Establish an equilibrium reference at explicit T,V and feed elements.
    /// The entire temperature interval must have the original source support.
    pub fn new(
        owner: &'a IdealGasEquilibriumV1,
        temperature: Temperature,
        volume: Volume,
        feed: &[(SpeciesId, Amount)],
        controls: EquilibriumInternalEnergySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, EquilibriumVolumeErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
            .contains(&temperature.value())
        {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "energy-volume reference temperature is outside the source bracket",
            ));
        }
        let mut work = owner.volume_workspace(volume, feed)?;
        work.admit_temperature_range(
            owner,
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        )?;
        let mut counts = Counts::default();
        let residual = work.solve_tv(
            owner,
            temperature,
            volume,
            controls.volume,
            &mut counts,
            &mut cancelled,
        )?;
        let reference = Reference::new(owner, &work, temperature, &mut cancelled)?;
        let reference_state = work.publish_volume(
            owner,
            temperature,
            volume,
            controls.volume,
            counts,
            0,
            None,
            residual,
            &mut cancelled,
        )?;
        let current = Some(EquilibriumEnergyVolumeSampleV1 {
            temperature,
            pressure: reference_state.equilibrium().pressure(),
            volume,
            response: reference_state.thermodynamic().response(),
            target_energy_change_j: 0.0,
            internal_energy_change_j: 0.0,
            energy_residual_j: 0.0,
            entropy_change_j_k: 0.0,
            relative_volume_residual: residual,
            relative_mass_residual: reference_state.equilibrium().relative_mass_residual(),
            composition_iterations: 0,
            temperature_iterations: 0,
            pressure_iterations: 0,
        });
        let mass = work.initial_mass * work.scale;
        if !mass.is_finite() || mass <= 0.0 {
            return Err(E::Arithmetic("original energy-volume feed mass").into());
        }
        poll(&mut cancelled)?;
        Ok(Self {
            owner,
            work,
            reference,
            reference_state,
            controls,
            current,
        })
    }

    /// Immutable original source-bearing reference.
    pub const fn reference_state(&self) -> &IdealGasIsochoricEquilibriumV1 {
        &self.reference_state
    }
    /// Original conserved feed mass, independent of solver product residuals.
    pub fn feed_mass_kg(&self) -> f64 {
        self.work.initial_mass * self.work.scale
    }

    /// Solve at a prescribed volume and total energy change from the reference.
    ///
    /// Heat/work contributions must already be combined with consistent signs.
    /// The native reacting Cv supplies the Newton derivative; chemistry, entropy
    /// and molecular count are allowed to change. Queries reuse Gibbs storage.
    pub fn sample(
        &mut self,
        volume: Volume,
        energy_change: Energy,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnergyVolumeSampleV1, EquilibriumVolumeErrorV1> {
        let previous = self.current.take();
        poll(&mut cancelled)?;
        if !volume.value().is_finite()
            || volume.value() <= 0.0
            || !energy_change.value().is_finite()
        {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "finite positive volume and finite energy change required",
            ));
        }
        if let Some(mut same) = previous
            .filter(|s| s.volume == volume && s.target_energy_change_j == energy_change.value())
        {
            same.composition_iterations = 0;
            same.temperature_iterations = 0;
            same.pressure_iterations = 0;
            self.current = Some(same);
            return Ok(same);
        }
        let preferred = previous.map(|s| {
            s.temperature.value()
                + (energy_change.value() - s.internal_energy_change_j)
                    / s.response.heat_capacity_v_j_k
        });
        let owner = self.owner;
        let reference = &self.reference;
        let (temperature, counts, temperature_iterations, volume_residual) =
            self.work.solve_internal_energy(
                owner,
                volume,
                energy_change.value(),
                self.controls,
                preferred,
                &mut cancelled,
                |work, t, cancel| Ok(reference.energy_change(owner, work, t, cancel)?),
            )?;
        let response = self.work.response(
            owner,
            temperature,
            self.controls.volume.equilibrium,
            &mut cancelled,
        )?;
        let du = reference.energy_change(owner, &self.work, temperature, &mut cancelled)?;
        let ds = reference.entropy_change(owner, &self.work, temperature, &mut cancelled)?;
        let mut mass = Sum::default();
        for &i in &self.work.active {
            mass.add(self.work.n[i] * owner.models[i].molar_mass().value());
        }
        let mass_residual = (mass.value() - self.work.initial_mass) / self.work.initial_mass;
        if !mass_residual.is_finite()
            || mass_residual.abs() > self.controls.volume.equilibrium.relative_balance_tolerance
        {
            return Err(E::MassBalance {
                relative_residual: mass_residual,
            }
            .into());
        }
        let residual = du - energy_change.value();
        if !residual.is_finite() || residual.abs() > self.controls.energy_tolerance_j {
            return Err(EquilibriumVolumeErrorV1::InternalEnergyNoConvergence {
                iterations: temperature_iterations,
                residual_j: residual,
            });
        }
        let sample = EquilibriumEnergyVolumeSampleV1 {
            temperature,
            pressure: self.work.pressure,
            volume,
            response,
            target_energy_change_j: energy_change.value(),
            internal_energy_change_j: du,
            energy_residual_j: residual,
            entropy_change_j_k: ds,
            relative_volume_residual: volume_residual,
            relative_mass_residual: mass_residual,
            composition_iterations: counts.composition,
            temperature_iterations,
            pressure_iterations: counts.pressure,
        };
        poll(&mut cancelled)?;
        self.current = Some(sample);
        Ok(sample)
    }

    /// Publish the last successful endpoint with original per-species receipts.
    /// The relative target/residual stays in the sample: adding a small change
    /// to a large absolute formation reference must not redefine the UV target.
    pub fn publish_current(
        &mut self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        poll(&mut cancelled)?;
        let sample = self.current.ok_or(EquilibriumVolumeErrorV1::Invalid(
            "no accepted energy-volume query to publish",
        ))?;
        self.work.publish_volume(
            self.owner,
            sample.temperature,
            sample.volume,
            self.controls.volume,
            Counts {
                composition: sample.composition_iterations,
                pressure: sample.pressure_iterations,
            },
            sample.temperature_iterations,
            None,
            sample.relative_volume_residual,
            &mut cancelled,
        )
    }
}

impl Workspace {
    // Shared absolute/reference-relative UV root. Each energy sample is evaluated
    // after TV closure. The optional previous-state guess is numerical only.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn solve_internal_energy<C: FnMut() -> bool>(
        &mut self,
        owner: &IdealGasEquilibriumV1,
        volume: Volume,
        target: f64,
        controls: EquilibriumInternalEnergySolveV1,
        preferred: Option<f64>,
        cancelled: &mut C,
        mut energy: impl FnMut(&Workspace, Temperature, &mut C) -> Result<f64, EquilibriumVolumeErrorV1>,
    ) -> Result<(Temperature, Counts, usize, f64), EquilibriumVolumeErrorV1> {
        let mut counts = Counts::default();
        let mut steps = 0;
        let mut low = controls.minimum_temperature_k;
        let mut high = controls.maximum_temperature_k;
        if let Some(t) = preferred.filter(|t| t.is_finite() && *t >= low && *t <= high) {
            // A predictor is an actual interior trial, even when it succeeds.
            // Mandatory bracket endpoints remain outside this trial allowance.
            steps = 1;
            let temperature = Temperature::new(t);
            let volume_residual = self.solve_tv(
                owner,
                temperature,
                volume,
                controls.volume,
                &mut counts,
                cancelled,
            )?;
            if (energy(self, temperature, cancelled)? - target).abs() <= controls.energy_tolerance_j
            {
                return Ok((temperature, counts, steps, volume_residual));
            }
        }
        let mut endpoints = [0.0; 2];
        for (index, t) in [low, high].into_iter().enumerate() {
            let temperature = Temperature::new(t);
            let volume_residual = self.solve_tv(
                owner,
                temperature,
                volume,
                controls.volume,
                &mut counts,
                cancelled,
            )?;
            endpoints[index] = energy(self, temperature, cancelled)?;
            if (endpoints[index] - target).abs() <= controls.energy_tolerance_j {
                return Ok((temperature, counts, steps, volume_residual));
            }
        }
        if !(endpoints[0] < target && target < endpoints[1]) {
            return Err(EquilibriumVolumeErrorV1::InternalEnergyNotBracketed {
                lower_energy_j: endpoints[0],
                upper_energy_j: endpoints[1],
            });
        }
        let mut t = preferred
            .filter(|t| t.is_finite() && *t > low && *t < high)
            .unwrap_or(low + 0.5 * (high - low));
        let mut residual = endpoints[1] - target;
        for step in (steps + 1)..=controls.maximum_temperature_iterations {
            steps = step;
            let temperature = Temperature::new(t);
            let volume_residual = self.solve_tv(
                owner,
                temperature,
                volume,
                controls.volume,
                &mut counts,
                cancelled,
            )?;
            residual = energy(self, temperature, cancelled)? - target;
            if residual.abs() <= controls.energy_tolerance_j {
                return Ok((temperature, counts, step, volume_residual));
            }
            if residual < 0.0 {
                low = t;
            } else {
                high = t;
            }
            let response =
                self.response(owner, temperature, controls.volume.equilibrium, cancelled)?;
            let next = t - residual / response.heat_capacity_v_j_k;
            let margin = 0.01 * (high - low);
            t = if next.is_finite() && next > low + margin && next < high - margin {
                next
            } else {
                low + 0.5 * (high - low)
            };
            if t == low || t == high {
                break;
            }
        }
        Err(EquilibriumVolumeErrorV1::InternalEnergyNoConvergence {
            iterations: steps,
            residual_j: residual,
        })
    }
}
