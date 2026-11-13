//! Mutually incoherent spectral-line heating on the original thermal cells.
//!
//! Each authored line keeps its own frequency, incidence, polarization and
//! two-face coherence. The existing Maxwell owner integrates its standing-wave
//! loss over each cell. Only those physical watts are added between lines;
//! one native thermal step receives their sum. A frequency difference alone
//! does not establish incoherence: selecting this owner declares zero mutual
//! cross-line second moments. It supplies no finite-time beat or sum-frequency
//! reservoir work, thermal emission, or transient electromagnetic storage.

use super::*;

/// One stationary optical line. Frequencies need not be sorted or distinct:
/// all lines, including equal-frequency entries, are explicitly uncorrelated.
/// Coherence between the two faces remains specified by each line's drive.
#[derive(Clone, Copy, Debug)]
pub struct ThermoOpticalSpectralLine {
    pub angular_frequency_rad_s: f64,
    pub incidence_cosine: f64,
    pub polarization: Polarization,
}

/// Cardinality bounds admitted before source refresh or field propagation.
/// A fixed number of O(cells * lines) passes admits/prepares each spectrum.
/// Factor contraction costs cells times the sum of per-line coherent modes:
/// one for Coherent, two for Independent, or the authored factor count.
/// Original per-line IlluminationBudget limits remain independently binding.
#[derive(Clone, Copy, Debug)]
pub struct BroadbandThermoOpticalLimits {
    /// At most 4096 lines are admitted, even if this bound is larger.
    pub maximum_lines: usize,
    pub maximum_cell_line_pairs: usize,
    pub maximum_cell_factor_pairs: usize,
}

/// Source-pinned broadband optical feedback with one shared thermal timeline.
/// At most 64 native cells and 4096 authored lines are retained. Operators are
/// prepared/evaluated one line at a time, using O(cells + lines) storage rather
/// than storing every line's field or a line-by-cell heat matrix.
pub struct BroadbandThermoOpticalWaveHeating<
    'b,
    'a,
    'm,
    S: ThermoOpticalSource = VisibleConductorStatePoint,
> {
    thermal: &'b BoundLayeredThermalBody<'a>,
    incident: DielectricResponse<'m>,
    exit: DielectricResponse<'m>,
    sources: Vec<ThermoOpticalLayer<'m, S>>,
    lines: Vec<ThermoOpticalSpectralLine>,
    limits: BroadbandThermoOpticalLimits,
}

/// One jointly admitted broadband thermal endpoint. The native step contains
/// summed physical cell energies and face powers; line_power retains each
/// original Maxwell account in authored line order, in W/m². Source vectors
/// retain original CellId/thickness order and are shared by every wavelength.
pub struct BroadbandThermoOpticalWaveStep<'t, S: ThermoOpticalSource = VisibleConductorStatePoint> {
    pub step: LayeredWaveStep<'t>,
    pub line_power: Vec<LayerIlluminationPower>,
    pub applied_optical_states: Vec<S>,
    pub endpoint_optical_states: Vec<S>,
}

impl<'b, 'a, 'm, S: ThermoOpticalSource> BroadbandThermoOpticalWaveHeating<'b, 'a, 'm, S> {
    /// Bind an explicitly mutually incoherent spectrum to the original body.
    /// Every source/line and both propagating exterior ports are admitted by
    /// the existing monochromatic owner. No spectrum or interpolation beyond
    /// a source's actual wavelength band and gap policy is manufactured.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: &'m ResolvedMaterialStatePoint,
        exit: DielectricResponse<'m>,
        lines: &[ThermoOpticalSpectralLine],
        sources: &[ThermoOpticalLayer<'m, S>],
        limits: BroadbandThermoOpticalLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        Self::try_new_with_incident_response(
            thermal,
            DielectricResponse::Lossless(incident),
            exit,
            lines,
            sources,
            limits,
            cancelled,
        )
    }

    /// Retain one original incident law across all authored spectral lines.
    /// Every line resolves its own source epsilon(omega), conserved transverse
    /// wave number and TE/TM port admittance. Its actual source band and exact
    /// zero-loss requirement are admitted per line by the original wave owner.
    /// The exterior material state is held during thermal feedback; no body
    /// temperature, averaged refractive index or additional bath heat replaces it.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_incident_response(
        thermal: &'b BoundLayeredThermalBody<'a>,
        incident: DielectricResponse<'m>,
        exit: DielectricResponse<'m>,
        lines: &[ThermoOpticalSpectralLine],
        sources: &[ThermoOpticalLayer<'m, S>],
        limits: BroadbandThermoOpticalLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let cells = thermal.prepared_body().cell_count();
        if lines.is_empty()
            || lines.len() > 4096
            || lines.len() > limits.maximum_lines
            || cells > MAX_LAYERS
            || cells
                .checked_mul(lines.len())
                .is_none_or(|work| work > limits.maximum_cell_line_pairs)
        {
            return Err(LayeredWaveHeatingError::Invalid(
                "broadband thermo-optical cell/line work exceeds its limits",
            ));
        }
        admit_sources(thermal, sources, &mut cancelled)?;
        for line in lines {
            poll(&mut cancelled)?;
            let mut responses = Vec::with_capacity(sources.len());
            for source in sources {
                poll(&mut cancelled)?;
                responses.push(source.source.wave_response(line.angular_frequency_rad_s)?);
            }
            let owner = LayeredWaveHeating::try_new_with_incident_response(
                thermal,
                incident,
                exit,
                line.angular_frequency_rad_s,
                &responses,
                &mut cancelled,
            )?;
            owner.prepare(line.incidence_cosine, line.polarization, &mut cancelled)?;
        }
        poll(&mut cancelled)?;
        Ok(Self {
            thermal,
            incident,
            exit,
            sources: sources.to_vec(),
            lines: lines.to_vec(),
            limits,
        })
    }

    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.thermal
    }
    /// Original held incident source, evaluated independently at each line.
    pub const fn incident_response(&self) -> DielectricResponse<'m> {
        self.incident
    }
    pub fn sources(&self) -> &[ThermoOpticalLayer<'m, S>] {
        &self.sources
    }
    pub fn lines(&self) -> &[ThermoOpticalSpectralLine] {
        &self.lines
    }
    pub const fn limits(&self) -> BroadbandThermoOpticalLimits {
        self.limits
    }

    fn admit_drives(
        &self,
        drives: &[LayeredWaveDrive<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), LayeredWaveHeatingError> {
        if drives.len() != self.lines.len() {
            return Err(LayeredWaveHeatingError::Invalid(
                "one explicitly incoherent drive per authored spectral line is required",
            ));
        }
        let cells = self.thermal.prepared_body().cell_count();
        let mut total = 0usize;
        for drive in drives {
            poll(cancelled)?;
            let modes = match drive {
                LayeredWaveDrive::Coherent(_) => 1,
                LayeredWaveDrive::Independent(_) => 2,
                LayeredWaveDrive::Coherence { factors, budget } => {
                    if factors.len() > 4096
                        || factors.len() > budget.maximum_modes
                        || factors
                            .len()
                            .checked_mul(cells)
                            .is_none_or(|work| work > budget.maximum_layer_mode_pairs)
                    {
                        return Err(LayeredWaveHeatingError::Invalid(
                            "broadband line coherence exceeds its original illumination budget",
                        ));
                    }
                    factors.len()
                }
            };
            total = modes
                .checked_mul(cells)
                .and_then(|work| total.checked_add(work))
                .filter(|&work| work <= self.limits.maximum_cell_factor_pairs)
                .ok_or(LayeredWaveHeatingError::Invalid(
                    "broadband thermo-optical factor work exceeds its limit",
                ))?;
        }
        Ok(())
    }

    /// Advance every line over the SAME physical interval, then advance heat
    /// once. Drives are peak tangential E [V/m], or spectral-line second-moment
    /// factors, in authored line order. They are not densities per Hz/metre;
    /// continuum quadrature weights must already be represented in their
    /// physical line powers (amplitudes scale by the square root of weight).
    ///
    /// Each original cell source is refreshed once at its accepted temperature.
    /// Coefficients and line statistics are held over this explicit interval;
    /// caller-owned step refinement remains necessary for thermo-optic accuracy.
    /// All endpoint source bands/gaps/frequencies are readmitted before return.
    /// Any refusal or cancellation leaves the accepted initial state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        drives: &[LayeredWaveDrive<'_>],
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BroadbandThermoOpticalWaveStep<'t, S>, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        self.admit_drives(drives, &mut cancelled)?;
        // Original timeline/geometry admission precedes spectral work. There
        // is exactly one held thermal step, independent of the line count.
        let thermal_step =
            LayeredThermalStep::new(self.thermal, initial, target_s, conduction, || {
                checkpoint(&mut cancelled)
            })?;
        let applied_optical_states =
            refresh_cell_sources(self.thermal, &self.sources, initial, &mut cancelled)?;
        let cells = applied_optical_states.len();
        let mut cell_sums = [Dd::ZERO; MAX_LAYERS];
        let mut power_sums = [Dd::ZERO; 5];
        let mut cell_power = [0.0; MAX_LAYERS];
        let cell_power = &mut cell_power[..cells];
        let mut line_power = Vec::with_capacity(self.lines.len());
        for (line, &drive) in self.lines.iter().zip(drives) {
            poll(&mut cancelled)?;
            let mut responses = Vec::with_capacity(cells);
            for source in &applied_optical_states {
                poll(&mut cancelled)?;
                responses.push(source.wave_response(line.angular_frequency_rad_s)?);
            }
            let owner = LayeredWaveHeating::try_new_cellwise_with_incident_response(
                self.thermal,
                self.incident,
                self.exit,
                line.angular_frequency_rad_s,
                &responses,
                &mut cancelled,
            )?;
            let power = owner
                .prepare(line.incidence_cosine, line.polarization, &mut cancelled)?
                .cell_power_into(drive, cell_power, &mut cancelled)?;
            // Accumulate positive per-cell dissipation, never incoming minus
            // outgoing fluxes and never fields from unrelated frequencies.
            for (sum, &value) in cell_sums.iter_mut().zip(cell_power.iter()) {
                *sum = *sum + Dd::from_f64(value);
            }
            for (sum, value) in power_sums.iter_mut().zip([
                power.incoming_w_m2[0],
                power.incoming_w_m2[1],
                power.outgoing_w_m2[0],
                power.outgoing_w_m2[1],
                power.absorbed_w_m2,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            line_power.push(power);
        }
        let totals = power_sums.map(Dd::to_f64);
        for (power, sum) in cell_power.iter_mut().zip(cell_sums) {
            *power = sum.to_f64();
        }
        if totals
            .iter()
            .chain(cell_power.iter())
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(LayeredWaveHeatingError::Invalid(
                "summed broadband optical power is unrepresentable",
            ));
        }
        let step = solve_held_cell_power(
            thermal_step,
            self.thermal.prepared_body().conduction().area_m2(),
            LayerIlluminationPower {
                incoming_w_m2: [totals[0], totals[1]],
                outgoing_w_m2: [totals[2], totals[3]],
                absorbed_w_m2: totals[4],
            },
            cell_power,
            controls,
            workspace,
            &mut cancelled,
        )?;
        let endpoint_optical_states =
            refresh_cell_sources(self.thermal, &self.sources, &step.thermal, &mut cancelled)?;
        for source in &endpoint_optical_states {
            for line in &self.lines {
                poll(&mut cancelled)?;
                source.wave_response(line.angular_frequency_rad_s)?;
            }
        }
        poll(&mut cancelled)?;
        Ok(BroadbandThermoOpticalWaveStep {
            step,
            line_power,
            applied_optical_states,
            endpoint_optical_states,
        })
    }
}
