//! Statistically stationary porous absorption across independent frequencies.
//!
//! Each bin retains its coherent left/right peak-pressure phasors and the
//! original native porous stack. Different frequencies contribute independent
//! cycle-mean powers, not complex pressure amplitudes. This is an ensemble or
//! sufficiently long-time average: finite-time beats, acoustic startup/storage,
//! frequency conversion and fluctuating thermal response are outside its scope.
//! Distinct bins describe independent frequency content, not independent emitters
//! at the same frequency. Combine same-frequency coherent phasors explicitly.
//!
//! Each bin's pressure is already integrated over that bin. For a continuous
//! band, [`angular_pressure_density_bin`] applies a positive authored quadrature
//! weight to a one-sided angular-frequency pressure amplitude density. Spectral
//! discretization accuracy remains the caller's responsibility. Two-sided
//! partial coherence within a bin is supplied by [`super::covariance`]; the
//! two phasors here share one coherent statistical component.

use super::*;

/// Hard ceiling for independently resolved frequencies, before allocation.
pub const MAX_POROUS_SPECTRUM_BINS: usize = 4096;
/// Hard ceiling for the complete cell × frequency coefficient/source cache.
pub const MAX_POROUS_SPECTRUM_CELL_BINS: usize = 65_536;

/// Explicit spectral storage/work limits plus the original stack admission.
#[derive(Clone, Copy, Debug)]
pub struct PorousSpectrumOptions {
    /// Maximum number of distinct, positive angular frequencies.
    pub maximum_bins: usize,
    /// Maximum thermal-cell count multiplied by frequency-bin count.
    pub maximum_cell_bins: usize,
    /// Original per-frequency scattering condition and composition limits.
    pub composition: ComposeOptions,
}
impl Default for PorousSpectrumOptions {
    fn default() -> Self {
        Self {
            maximum_bins: 128,
            maximum_cell_bins: 16_384,
            composition: ComposeOptions::default(),
        }
    }
}

fn admit_drive(drive: PorousHarmonicDrive) -> Result<(), PorousHeatingError> {
    if ![
        drive.angular_frequency_rad_s,
        drive.area_m2,
        drive.reference_impedance_pa_s_m3,
    ]
    .into_iter()
    .all(|x| x.is_finite() && x > 0.0)
        || drive
            .incident_pressure_pa
            .into_iter()
            .any(|p| !p.re.is_finite() || !p.im.is_finite())
    {
        return Err(PorousHeatingError::Invalid(
            "spectral bins need positive frequency/area/reference and finite peak pressure",
        ));
    }
    Ok(())
}

/// Convert one quadrature node of a ONE-SIDED angular-frequency spectrum.
/// `amplitude_density_pa_per_sqrt_rad_s[j]` is the complex pressure amplitude
/// density at exterior port j: its squared magnitude is S_pp(omega), with
/// mean-square pressure = integral_0^infinity S_pp(omega) d omega. The two
/// densities retain their relative phase and are fully coherent within this bin.
/// The returned PEAK pressure is density * sqrt(2 * quadrature_weight_rad_s),
/// so incident watts are weight * sum_j |density[j]|² / reference_impedance.
/// The weight already includes any authored quadrature coefficient. No other
/// bandwidth, factor 1/2, or 2 pi is applied by the heating adapter.
/// For a per-Hz density, first use density_omega = density_Hz / sqrt(2 pi)
/// and weight_omega = 2 pi * weight_Hz. Both weights must be strictly positive.
pub fn angular_pressure_density_bin(
    angular_frequency_rad_s: f64,
    area_m2: f64,
    reference_impedance_pa_s_m3: f64,
    amplitude_density_pa_per_sqrt_rad_s: [C64; 2],
    quadrature_weight_rad_s: f64,
) -> Result<PorousHarmonicDrive, PorousHeatingError> {
    if !quadrature_weight_rad_s.is_finite() || quadrature_weight_rad_s <= 0.0 {
        return Err(PorousHeatingError::Invalid(
            "pressure spectral density needs a positive finite angular weight",
        ));
    }
    // Avoid overflowing 2*weight when its square root is representable.
    let scale = det::sqrt(quadrature_weight_rad_s) * core::f64::consts::SQRT_2;
    let drive = PorousHarmonicDrive {
        angular_frequency_rad_s,
        area_m2,
        reference_impedance_pa_s_m3,
        incident_pressure_pa: amplitude_density_pa_per_sqrt_rad_s.map(|p| p.scale(scale)),
    };
    admit_drive(drive)?;
    Ok(drive)
}

/// Current mean power summed across the statistically independent frequencies.
/// There is no aggregate complex pressure or aggregate reflection coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousSpectrumPower {
    /// Prescribed incident power through the two exterior faces [W].
    pub incident_power_w: f64,
    /// Outgoing power through those same exterior faces [W].
    pub outgoing_power_w: f64,
    /// Native positive viscous and thermal work across all original cells [W].
    pub dissipation: PorousLayerDissipation,
    /// Incident minus outgoing minus deposited power; never used as cell heat.
    pub residual_w: f64,
}

/// One native thermal specimen driven by a bounded stationary acoustic spectrum.
/// Bins are stored in increasing frequency order for deterministic summation.
/// Their material objects, cell geometry and gas specification are common;
/// pressure phases and fixed exterior references may differ across frequencies.
/// All isobaric, caloric-storage and quasisteady assumptions of the original
/// [`TemperaturePorousHeating`] also apply to this shared feedback trajectory.
#[derive(Clone)]
pub struct PorousSpectrumHeating<'a> {
    bins: Vec<TemperaturePorousHeating<'a>>,
    heat: PreparedHeat,
}

/// Jointly accepted spectrum cache, native thermal state and physical accounts.
pub struct PorousSpectrumStep<'t, 'a> {
    /// Reusable per-frequency caches at the returned native temperatures.
    pub heating: PorousSpectrumHeating<'a>,
    /// Shared thermal advance and total acoustic heat, credited exactly once.
    pub step: PorousHeatingStep<'t>,
    /// Per-frequency acoustic energy [J], in increasing frequency order.
    /// These diagnose the total account; they are not additional heat sources.
    pub bin_energy: Vec<PorousHeatingEnergy>,
    /// Gas evaluations are shared; stack solves count each affected frequency.
    pub work: PorousFeedbackWork,
}

impl<'a> PorousSpectrumHeating<'a> {
    /// Admit one common physical cell partition and a nonempty set of distinct
    /// positive frequencies. Each peak-pressure drive already represents its
    /// bin's integrated power. Duplicate frequencies refuse before any acoustic
    /// preparation; the caller must combine same-frequency coherent waves.
    /// Frequency and total cell×bin limits are checked before cache allocation.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: IsobaricPorousGas<'a>,
        bins: &[PorousHarmonicDrive],
        options: PorousSpectrumOptions,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        checkpoint()?;
        if options.maximum_bins == 0
            || options.maximum_bins > MAX_POROUS_SPECTRUM_BINS
            || options.maximum_cell_bins == 0
            || options.maximum_cell_bins > MAX_POROUS_SPECTRUM_CELL_BINS
            || bins.is_empty()
            || bins.len() > options.maximum_bins
            || cells.is_empty()
            || cells.len() > MAX_POROUS_HEATING_CELLS
            || cells.len() > options.composition.maximum_blocks
            || cells
                .len()
                .checked_mul(bins.len())
                .is_none_or(|count| count > options.maximum_cell_bins)
        {
            return Err(DepositionError::Budget("porous spectrum cells and frequencies").into());
        }
        for &drive in bins {
            checkpoint()?;
            admit_drive(drive)?;
            if drive.area_m2 != bins[0].area_m2 {
                return Err(PorousHeatingError::Invalid(
                    "all frequency bins must illuminate the same full specimen area",
                ));
            }
        }
        let mut drives = bins.to_vec();
        drives.sort_by(|a, b| {
            a.angular_frequency_rad_s
                .total_cmp(&b.angular_frequency_rad_s)
        });
        if drives
            .windows(2)
            .any(|pair| pair[0].angular_frequency_rad_s == pair[1].angular_frequency_rad_s)
        {
            return Err(PorousHeatingError::Invalid(
                "duplicate frequency: combine coherent phasors in one authored bin",
            ));
        }
        let mut states = Vec::with_capacity(cells.len());
        for cell in cells {
            checkpoint()?;
            states.push(gas.evaluate(thermal.temperature_k(cell.cell)?)?);
        }
        let mut prepared = Vec::with_capacity(drives.len());
        for drive in drives {
            checkpoint()?;
            prepared.push(TemperaturePorousHeating::prepare(
                thermal,
                cells,
                gas,
                drive,
                options.composition,
                Some(&states),
                &mut checkpoint,
            )?);
        }
        let mut result = Self {
            heat: prepared[0].heat.clone(),
            bins: prepared,
        };
        result.sum_power(&mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }

    /// Native single-frequency diagnostics in canonical increasing order.
    pub fn bins(&self) -> &[TemperaturePorousHeating<'a>] {
        &self.bins
    }
    /// Common physical cells with the sum of their positive per-bin heat rates.
    pub fn cells(&self) -> &[PorousCellHeating] {
        &self.heat.cells
    }
    /// One native gas state per cell, evaluated once and shared across bins.
    pub fn gas_states(&self) -> &[GasState] {
        self.bins[0].gas_states()
    }
    /// Original fixed gas composition, pressure and admitted temperature range.
    pub fn gas(&self) -> IsobaricPorousGas<'a> {
        self.bins[0].gas()
    }
    /// Total mean exterior power and direct physical dissipation.
    pub fn power(&self) -> PorousSpectrumPower {
        PorousSpectrumPower {
            incident_power_w: self.heat.boundary_power_w[0],
            outgoing_power_w: self.heat.boundary_power_w[1],
            dissipation: self.heat.dissipation,
            residual_w: self.heat.power_residual_w,
        }
    }

    fn sum_power(
        &mut self,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        sum_heat_sources(
            &mut self.heat,
            self.bins.iter().map(|bin| &bin.heat),
            checkpoint,
        )
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        work: &mut PorousFeedbackWork,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        self.heat.validate(thermal, checkpoint)?;
        let mut changed = Vec::new();
        for (index, cell) in self.heat.cells.iter().enumerate() {
            checkpoint()?;
            let temperature = thermal.temperature_k(cell.cell)?;
            if temperature.to_bits() == self.bins[0].states[index].temperature.to_bits() {
                continue;
            }
            changed.push((index, self.gas().evaluate(temperature)?));
            work.gas_evaluations = work
                .gas_evaluations
                .checked_add(1)
                .ok_or(DepositionError::Budget("porous spectral gas evaluations"))?;
        }
        if changed.is_empty() {
            return Ok(());
        }
        for bin in &mut self.bins {
            checkpoint()?;
            bin.refresh_changed(&changed, work, checkpoint)?;
        }
        self.sum_power(checkpoint)
    }

    /// One accepted-start explicit feedback trajectory for the ENTIRE spectrum.
    /// All bins see the same evolving native temperatures, authored heat and
    /// conduction. Their positive heat is summed before the native step, not
    /// applied through separate thermal advances. Endpoint source admission and
    /// cancellation are transactional. Refine both spectral quadrature and
    /// thermal timesteps; native conduction CFL does not control source feedback.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousSpectrumStep<'t, 'a>, PorousHeatingError> {
        checkpoint()?;
        let mut observer = SpectrumObserver {
            source: self.clone(),
            account: HeatAccount::new(self.heat.cells.len()),
            bin_accounts: self
                .bins
                .iter()
                .map(|bin| HeatAccount::new(bin.heat.cells.len()))
                .collect(),
            work: PorousFeedbackWork::default(),
        };
        observer
            .source
            .refresh(initial, &mut observer.work, &mut checkpoint)?;
        let (thermal, observer, heat_report) =
            initial.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        let mut bin_energy = Vec::with_capacity(self.bins.len());
        for account in &observer.bin_accounts {
            checkpoint()?;
            bin_energy.push(account.acoustic_energy()?);
        }
        let step = observer.account.finish(
            &observer.source.heat,
            thermal,
            heat_report,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(PorousSpectrumStep {
            heating: observer.source,
            step,
            bin_energy,
            work: observer.work,
        })
    }
}

// One positive source sum serves coherent and covariance spectra.
pub(super) fn sum_heat_sources<'h>(
    heat: &mut PreparedHeat,
    bins: impl ExactSizeIterator<Item = &'h PreparedHeat> + Clone,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), PorousHeatingError> {
    // Retain exact original single-bin arithmetic and diagnostics.
    if bins.len() == 1 {
        heat.clone_from(bins.clone().next().expect("nonempty spectrum"));
        return Ok(());
    }
    let mut total = [Dd::ZERO; 3];
    for (index, cell) in heat.cells.iter_mut().enumerate() {
        let mut sum = [Dd::ZERO; 3];
        for bin in bins.clone() {
            checkpoint()?;
            let loss = bin.cells[index].dissipation;
            for (sum, value) in sum.iter_mut().zip([
                loss.viscous_power_w,
                loss.thermal_power_w,
                loss.absorbed_power_w,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        cell.dissipation = PorousLayerDissipation {
            viscous_power_w: finite(sum[0].to_f64())?,
            thermal_power_w: finite(sum[1].to_f64())?,
            absorbed_power_w: finite(sum[2].to_f64())?,
        };
        *heat
            .sources_w
            .get_mut(&cell.cell)
            .expect("original admitted spectral cell") = cell.dissipation.absorbed_power_w;
        for (total, value) in total.iter_mut().zip(sum) {
            *total = *total + value;
        }
    }
    let mut boundary = [Dd::ZERO; 2];
    for bin in bins.clone() {
        checkpoint()?;
        for (sum, value) in boundary.iter_mut().zip(bin.boundary_power_w) {
            *sum = *sum + Dd::from_f64(value);
        }
    }
    heat.dissipation = PorousLayerDissipation {
        viscous_power_w: finite(total[0].to_f64())?,
        thermal_power_w: finite(total[1].to_f64())?,
        absorbed_power_w: finite(total[2].to_f64())?,
    };
    heat.boundary_power_w = [finite(boundary[0].to_f64())?, finite(boundary[1].to_f64())?];
    // Check boundary closure independently of the positive local sources.
    // Do not recover absorption by subtracting nearly equal boundary watts.
    let count = heat.cells.len() * bins.len();
    let residual = finite((boundary[0] - boundary[1] - total[2]).to_f64())?;
    let tolerance = finite(
        512.0
            * f64::EPSILON
            * (count + 1) as f64
            * heat.boundary_power_w[0]
                .max(heat.boundary_power_w[1])
                .max(heat.dissipation.absorbed_power_w),
    )?;
    if residual.abs() > tolerance {
        return Err(PorousHeatingError::PowerBalance {
            residual_w: residual,
            tolerance_w: tolerance,
        });
    }
    heat.power_residual_w = residual;
    Ok(())
}

#[derive(Clone)]
struct SpectrumObserver<'a> {
    source: PorousSpectrumHeating<'a>,
    account: HeatAccount,
    bin_accounts: Vec<HeatAccount>,
    work: PorousFeedbackWork,
}
impl DepositionObserver for SpectrumObserver<'_> {
    type Error = PorousHeatingError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.source.refresh(thermal, &mut self.work, checkpoint)?;
        Ok(self.source.heat.sources_w.clone())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        // All accounts use the rates that actually drove this accepted step.
        self.account.record(&self.source.heat, dt, checkpoint)?;
        for (account, bin) in self.bin_accounts.iter_mut().zip(&self.source.bins) {
            account.record(&bin.heat, dt, checkpoint)?;
        }
        self.source.refresh(thermal, &mut self.work, checkpoint)
    }
}
