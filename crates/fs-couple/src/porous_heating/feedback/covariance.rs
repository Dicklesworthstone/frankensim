//! Partially coherent stationary spectra on the original thermal timeline.
//! Each physical frequency owns one Hermitian PSD RMS pressure covariance.
//! Its at-most-two positive components share ONE native composed porous stack.
//! Gas state is evaluated once per changed cell across all frequencies; direct
//! viscous/thermal work sums before one accepted-start thermal advance.
//!
//! Frequencies remain statistically independent; same-frequency independent
//! emitters must first add their covariance matrices. Coherent emitters must
//! instead combine phasors or include their cross correlation. This is a
//! stationary ensemble/long-time mean, not finite-time beating, noise-driven
//! thermal fluctuations or acoustic transients. All original isobaric gas,
//! rigid-frame, parallel plane-wave and thermally lumped storage assumptions apply.

use super::spectrum::{
    MAX_POROUS_SPECTRUM_BINS, MAX_POROUS_SPECTRUM_CELL_BINS, PorousSpectrumOptions,
    PorousSpectrumPower, sum_heat_sources,
};
use super::*;
use fs_duct::porous::{PorousStackExcitation, PressureCovariance};

/// One bin's integrated RMS pressure covariance at the two EXTERIOR faces.
#[derive(Clone, Copy, Debug)]
pub struct PorousCovarianceBin {
    /// Distinct positive angular frequency [rad/s], inside native source support.
    pub angular_frequency_rad_s: f64,
    /// Common full physical specimen area [m²], never the pore area.
    pub area_m2: f64,
    /// Fixed positive real reference [Pa s/m³] for both faces at this frequency.
    pub reference_impedance_pa_s_m3: f64,
    /// C_ij=E[u_i conj(u_j)] [Pa²] for RMS amplitudes; peak p gives p p*/2.
    /// `PressureCovariance::from_angular_density` includes bin quadrature weight.
    pub covariance: PressureCovariance,
}

#[derive(Clone)]
struct BinResponse<'a> {
    layers: Vec<PorousLayer<'a>>,
    tangential_wavenumber_rad_m: f64,
    heat: PreparedHeat,
}

/// A common native thermal specimen under a partially coherent spectrum.
/// Native coefficients/geometry are prepared once per physical frequency,
/// not once per covariance component. The source retains the original skeleton,
/// caloric law, mass and explicit isobaric gas binding at every thermal cell.
#[derive(Clone)]
pub struct PorousCrossSpectrumHeating<'a> {
    bins: Vec<PorousCovarianceBin>,
    responses: Vec<BinResponse<'a>>,
    states: Vec<GasState>,
    gas: IsobaricPorousGas<'a>,
    options: ComposeOptions,
    heat: PreparedHeat,
}

/// One jointly accepted native state, spectrum cache and physical energy account.
pub struct PorousCrossSpectrumStep<'t, 'a> {
    /// Cache evaluated at the returned native thermal temperatures.
    pub heating: PorousCrossSpectrumHeating<'a>,
    /// Native thermal candidate and total acoustic heat, credited once.
    pub step: PorousHeatingStep<'t>,
    /// Diagnostic per-frequency energy [J], in increasing frequency order.
    pub bin_energy: Vec<PorousHeatingEnergy>,
    /// Gas counts are per changed cell; stack counts are per physical frequency.
    pub work: PorousFeedbackWork,
}

fn power(heat: &PreparedHeat) -> PorousSpectrumPower {
    PorousSpectrumPower {
        incident_power_w: heat.boundary_power_w[0],
        outgoing_power_w: heat.boundary_power_w[1],
        dissipation: heat.dissipation,
        residual_w: heat.power_residual_w,
    }
}

// Both components use the original analytic per-layer loss; no boundary deficit
// supplies cell heat. With at most two terms, each positive addition rounds only
// once. The first component retains the original coherent arithmetic exactly.
fn component_heat(
    heat: &mut PreparedHeat,
    fields: &[PorousStackExcitation<'_, '_, '_>],
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), PorousHeatingError> {
    let first = &fields[0];
    for (i, cell) in heat.cells.iter_mut().enumerate() {
        checkpoint()?;
        let local = &first.layer_excitations()[i];
        cell.dissipation = local.dissipation_between(0.0, local.layer().thickness_m())?;
        if fields.len() == 2 {
            let local = &fields[1].layer_excitations()[i];
            let next = local.dissipation_between(0.0, local.layer().thickness_m())?;
            cell.dissipation.viscous_power_w =
                finite(cell.dissipation.viscous_power_w + next.viscous_power_w)?;
            cell.dissipation.thermal_power_w =
                finite(cell.dissipation.thermal_power_w + next.thermal_power_w)?;
            cell.dissipation.absorbed_power_w =
                finite(cell.dissipation.absorbed_power_w + next.absorbed_power_w)?;
        }
        *heat
            .sources_w
            .get_mut(&cell.cell)
            .expect("original covariance heat cell") = cell.dissipation.absorbed_power_w;
    }
    heat.dissipation = first.dissipation();
    heat.boundary_power_w = [
        first.boundary_power().incident_power_w,
        first.boundary_power().outgoing_power_w,
    ];
    heat.power_residual_w = first.power_residual_w();
    if fields.len() == 2 {
        let next = fields[1].dissipation();
        let b = fields[1].boundary_power();
        let incoming = Dd::from_f64(heat.boundary_power_w[0]) + Dd::from_f64(b.incident_power_w);
        let outgoing = Dd::from_f64(heat.boundary_power_w[1]) + Dd::from_f64(b.outgoing_power_w);
        let absorbed =
            Dd::from_f64(heat.dissipation.absorbed_power_w) + Dd::from_f64(next.absorbed_power_w);
        heat.dissipation.viscous_power_w =
            finite(heat.dissipation.viscous_power_w + next.viscous_power_w)?;
        heat.dissipation.thermal_power_w =
            finite(heat.dissipation.thermal_power_w + next.thermal_power_w)?;
        heat.dissipation.absorbed_power_w = finite(absorbed.to_f64())?;
        heat.boundary_power_w = [finite(incoming.to_f64())?, finite(outgoing.to_f64())?];
        heat.power_residual_w = finite((incoming - outgoing - absorbed).to_f64())?;
    }
    checkpoint()?;
    Ok(())
}

impl<'a> PorousCrossSpectrumHeating<'a> {
    /// Admit the complete physical thermal partition and a nonempty covariance
    /// spectrum. Frequencies are canonicalized; duplicates refuse rather than
    /// being interpreted as uncorrelated same-frequency drives. Bin/cell-product
    /// limits precede allocations. One native gas state per cell is shared by
    /// all bins, each composed once before its at-most-two coherent excitations.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: IsobaricPorousGas<'a>,
        bins: &[PorousCovarianceBin],
        options: PorousSpectrumOptions,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        Self::prepare(thermal, cells, gas, bins, None, options, checkpoint)
    }

    /// Attach one conserved real tangential wavenumber [rad/m] to each INPUT
    /// frequency bin before canonical sorting. Both covariance components use
    /// that same spatial mode. The fixed real reference is the exterior normal
    /// pressure/volume-flow wave impedance; all original area/storage rules apply.
    pub fn try_new_oblique(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: IsobaricPorousGas<'a>,
        bins: &[PorousCovarianceBin],
        tangential_wavenumbers_rad_m: &[f64],
        options: PorousSpectrumOptions,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        Self::prepare(
            thermal,
            cells,
            gas,
            bins,
            Some(tangential_wavenumbers_rad_m),
            options,
            checkpoint,
        )
    }

    fn prepare(
        thermal: &DepositionSimulation<'_>,
        cells: &[PorousTemperatureCell<'a>],
        gas: IsobaricPorousGas<'a>,
        bins: &[PorousCovarianceBin],
        tangential: Option<&[f64]>,
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
                .is_none_or(|n| n > options.maximum_cell_bins)
        {
            return Err(
                DepositionError::Budget("porous cross-spectrum cells and frequencies").into(),
            );
        }
        if let Some(values) = tangential {
            if values.len() != bins.len() {
                return Err(PorousHeatingError::Invalid(
                    "one tangential wavenumber per covariance frequency required",
                ));
            }
            for &value in values {
                checkpoint()?;
                admit_tangential_wavenumber(value)?;
            }
        }
        for bin in bins {
            checkpoint()?;
            if ![
                bin.angular_frequency_rad_s,
                bin.area_m2,
                bin.reference_impedance_pa_s_m3,
            ]
            .into_iter()
            .all(|v| v.is_finite() && v > 0.0)
                || bin.area_m2 != bins[0].area_m2
            {
                return Err(PorousHeatingError::Invalid(
                    "cross-spectrum needs positive frequency/reference and one full area",
                ));
            }
        }
        let mut ordered: Vec<_> = bins
            .iter()
            .enumerate()
            .map(|(i, &bin)| (bin, tangential.map_or(0.0, |values| values[i])))
            .collect();
        ordered.sort_by(|a, b| {
            a.0.angular_frequency_rad_s
                .total_cmp(&b.0.angular_frequency_rad_s)
        });
        if ordered
            .windows(2)
            .any(|p| p[0].0.angular_frequency_rad_s == p[1].0.angular_frequency_rad_s)
        {
            return Err(PorousHeatingError::Invalid(
                "duplicate frequency: author one total covariance per frequency",
            ));
        }
        let mut states = Vec::with_capacity(cells.len());
        for cell in cells {
            checkpoint()?;
            states.push(gas.evaluate(thermal.temperature_k(cell.cell)?)?);
        }
        let mut responses = Vec::with_capacity(bins.len());
        for &(bin, tangential) in &ordered {
            let mut layers = Vec::with_capacity(cells.len());
            let mut partitions = Vec::with_capacity(cells.len());
            for (cell, state) in cells.iter().zip(&states) {
                checkpoint()?;
                let original = thermal
                    .cells()
                    .get(&cell.cell)
                    .ok_or(PorousHeatingError::Invalid("inactive covariance heat cell"))?;
                let length = finite(norm(std::array::from_fn(|i| {
                    original.geometry.end_m[i] - original.geometry.start_m[i]
                })))?;
                layers.push(PorousLayer::try_new_oblique(
                    cell.material,
                    state,
                    bin.angular_frequency_rad_s,
                    length,
                    bin.area_m2,
                    bin.reference_impedance_pa_s_m3,
                    tangential,
                )?);
                partitions.push([PorousHeatCell {
                    cell: cell.cell,
                    end_depth_m: length,
                }]);
            }
            let heat = {
                let stack = acoustic(&mut checkpoint, |cancelled| {
                    PorousLayerStack::try_new(&layers, options.composition, cancelled)
                })?;
                let fields = acoustic(&mut checkpoint, |cancelled| {
                    stack.prepare_covariance(bin.covariance, cancelled)
                })?;
                let partition_refs: Vec<&[PorousHeatCell]> =
                    partitions.iter().map(|p| p.as_slice()).collect();
                let bound = PorousStackHeating::try_new(
                    &fields[0],
                    thermal,
                    &partition_refs,
                    cells.len(),
                    &mut checkpoint,
                )?;
                let mut heat = PreparedHeat {
                    cells: bound.cells().to_vec(),
                    sources_w: bound
                        .cells()
                        .iter()
                        .map(|c| (c.cell, c.dissipation.absorbed_power_w))
                        .collect(),
                    dissipation: bound.dissipation(),
                    power_residual_w: bound.power_residual_w(),
                    boundary_power_w: [
                        fields[0].boundary_power().incident_power_w,
                        fields[0].boundary_power().outgoing_power_w,
                    ],
                };
                component_heat(&mut heat, &fields, &mut checkpoint)?;
                heat
            };
            responses.push(BinResponse {
                layers,
                heat,
                tangential_wavenumber_rad_m: if tangential == 0.0 { 0.0 } else { tangential },
            });
        }
        let mut heat = responses[0].heat.clone();
        sum_heat_sources(
            &mut heat,
            responses.iter().map(|b| &b.heat),
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(Self {
            bins: ordered.into_iter().map(|(bin, _)| bin).collect(),
            responses,
            states,
            gas,
            options: options.composition,
            heat,
        })
    }

    /// Original covariance bins in increasing physical frequency order.
    pub fn bins(&self) -> &[PorousCovarianceBin] {
        &self.bins
    }
    /// Conserved tangential wavenumber of a canonically ordered bin [rad/m].
    pub fn tangential_wavenumber_rad_m(&self, index: usize) -> Option<f64> {
        self.responses
            .get(index)
            .map(|response| response.tangential_wavenumber_rad_m)
    }
    /// Common native physical heat cells with all positive bin/component watts.
    pub fn cells(&self) -> &[PorousCellHeating] {
        &self.heat.cells
    }
    /// One cached original gas state per native thermal cell.
    pub fn gas_states(&self) -> &[GasState] {
        &self.states
    }
    /// Prescribed native gas composition, pressure and explicit support.
    pub const fn gas(&self) -> IsobaricPorousGas<'a> {
        self.gas
    }
    /// Exterior mean powers and direct native cell dissipation across all bins.
    pub fn power(&self) -> PorousSpectrumPower {
        power(&self.heat)
    }
    /// One physical frequency's powers, including all its covariance components.
    pub fn bin_power(&self, index: usize) -> Option<PorousSpectrumPower> {
        self.responses.get(index).map(|bin| power(&bin.heat))
    }
    /// Per-cell heat for one physical frequency in original depth order.
    pub fn bin_cells(&self, index: usize) -> Option<&[PorousCellHeating]> {
        self.responses
            .get(index)
            .map(|bin| bin.heat.cells.as_slice())
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        work: &mut PorousFeedbackWork,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        self.heat.validate(thermal, checkpoint)?;
        let mut changed = Vec::new();
        for (i, cell) in self.heat.cells.iter().enumerate() {
            checkpoint()?;
            let temperature = thermal.temperature_k(cell.cell)?;
            if temperature.to_bits() == self.states[i].temperature.to_bits() {
                continue;
            }
            self.states[i] = self.gas.evaluate(temperature)?;
            changed.push(i);
            work.gas_evaluations = work
                .gas_evaluations
                .checked_add(1)
                .ok_or(DepositionError::Budget("cross-spectrum gas evaluations"))?;
        }
        if changed.is_empty() {
            return Ok(());
        }
        for (bin, response) in self.bins.iter().zip(&mut self.responses) {
            for &i in &changed {
                checkpoint()?;
                response.layers[i] = PorousLayer::try_new_oblique(
                    response.layers[i].material(),
                    &self.states[i],
                    bin.angular_frequency_rad_s,
                    response.layers[i].thickness_m(),
                    bin.area_m2,
                    bin.reference_impedance_pa_s_m3,
                    response.tangential_wavenumber_rad_m,
                )?;
            }
            let stack = acoustic(checkpoint, |cancelled| {
                PorousLayerStack::try_new(&response.layers, self.options, cancelled)
            })?;
            let fields = acoustic(checkpoint, |cancelled| {
                stack.prepare_covariance(bin.covariance, cancelled)
            })?;
            component_heat(&mut response.heat, &fields, checkpoint)?;
            work.stack_solves = work
                .stack_solves
                .checked_add(1)
                .ok_or(DepositionError::Budget("cross-spectrum stack evaluations"))?;
        }
        sum_heat_sources(
            &mut self.heat,
            self.responses.iter().map(|b| &b.heat),
            checkpoint,
        )
    }

    /// Advance all covariance bins through ONE accepted-start native observer.
    /// Source energy is charged before endpoint refresh; returned caches match
    /// the accepted thermal temperatures. Original caloric storage, conduction,
    /// heat/events and binding identities remain native. Domain/cancellation
    /// failure publishes neither candidate. Refine thermal steps and frequency
    /// quadrature independently of the unchanged native conduction CFL.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousCrossSpectrumStep<'t, 'a>, PorousHeatingError> {
        checkpoint()?;
        let mut observer = CrossSpectrumObserver {
            source: self.clone(),
            account: HeatAccount::new(self.heat.cells.len()),
            bin_accounts: self
                .responses
                .iter()
                .map(|b| HeatAccount::new(b.heat.cells.len()))
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
        Ok(PorousCrossSpectrumStep {
            heating: observer.source,
            step,
            bin_energy,
            work: observer.work,
        })
    }
}

#[derive(Clone)]
struct CrossSpectrumObserver<'a> {
    source: PorousCrossSpectrumHeating<'a>,
    account: HeatAccount,
    bin_accounts: Vec<HeatAccount>,
    work: PorousFeedbackWork,
}
impl DepositionObserver for CrossSpectrumObserver<'_> {
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
        self.account.record(&self.source.heat, dt, checkpoint)?;
        for (account, bin) in self.bin_accounts.iter_mut().zip(&self.source.responses) {
            account.record(&bin.heat, dt, checkpoint)?;
        }
        self.source.refresh(thermal, &mut self.work, checkpoint)
    }
}
