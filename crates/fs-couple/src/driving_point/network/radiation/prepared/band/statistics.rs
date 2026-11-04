//! Full stationary port moments, receiver fields and original signed work.
use super::*;
use fs_duct::network::NetworkSourceWork;

/// Reusable output and per-bin workspaces, bound to the admitted band. All
/// getters containing observations become unavailable after ANY refusal.
pub struct RadiatingWelchBandStatistics<'b, 'p, 'a> {
    band: PreparedRadiatingWelchBand<'b, 'p, 'a>,
    valid: bool,
    per_bin: Vec<RadiatingSourceStatistics<'p, 'a>>,
    pressure: Vec<C64>,
    flow: Vec<C64>,
    pressure_flow: Vec<C64>,
    fields: Vec<HarmonicFieldStatistics>,
    work: NetworkSourceWork,
    exterior_surface_power_w: f64,
    workspace_bytes: usize,
}

impl<'b, 'p, 'a> PreparedRadiatingWelchBand<'b, 'p, 'a> {
    /// Admit all retained original per-bin statistics, their copied port-label
    /// payloads and the aggregate arrays before allocation. The logical byte
    /// count excludes borrowed prepared bins and allocator overhead.
    pub fn statistics(
        &self,
        maximum_workspace_bytes: usize,
        gate: &CancelGate,
    ) -> Result<RadiatingWelchBandStatistics<'b, 'p, 'a>, ExteriorNetworkError> {
        poll(gate)?;
        let n = self.ports().len();
        let r = self.receivers().len();
        let matrix = mul(n, n)?;
        // Original duct retains value/correction arrays for three matrices,
        // two per-port work vectors and one dissipation value.
        let original_moments = mul(2, add(add(mul(3, matrix)?, mul(2, n)?)?, 1)?)?;
        let mut labels = bytes::<String>(n)?;
        for label in self.ports() {
            labels = add(labels, label.len())?;
        }
        let per_bin = add(
            add(bytes::<RadiatingSourceStatistics<'_, '_>>(1)?, labels)?,
            add(
                bytes::<C64>(original_moments)?,
                bytes::<HarmonicFieldStatistics>(r)?,
            )?,
        )?;
        let workspace_bytes = add(
            bytes::<RadiatingWelchBandStatistics<'_, '_, '_>>(1)?,
            add(
                mul(self.bins.len(), per_bin)?,
                add(
                    bytes::<C64>(mul(3, matrix)?)?,
                    bytes::<HarmonicFieldStatistics>(r)?,
                )?,
            )?,
        )?;
        storage(workspace_bytes, maximum_workspace_bytes)?;
        let mut per_bin = reserve(self.bins.len())?;
        for bin in self.bins {
            poll(gate)?;
            per_bin.push(bin.owner.statistics_buffer());
        }
        let result = RadiatingWelchBandStatistics {
            band: *self,
            valid: false,
            per_bin,
            pressure: zeros(matrix)?,
            flow: zeros(matrix)?,
            pressure_flow: zeros(matrix)?,
            fields: zeros(r)?,
            work: NetworkSourceWork::default(),
            exterior_surface_power_w: 0.0,
            workspace_bytes,
        };
        poll(gate)?;
        Ok(result)
    }
}

impl RadiatingWelchBandStatistics<'_, '_, '_> {
    pub fn band(&self) -> &PreparedRadiatingWelchBand<'_, '_, '_> {
        &self.band
    }
    pub fn is_valid(&self) -> bool {
        self.valid
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    /// Original per-frequency statistics, in the explicitly selected bin order.
    pub fn bin_statistics(&self) -> Option<&[RadiatingSourceStatistics<'_, '_>]> {
        self.valid.then_some(&self.per_bin)
    }
    /// Summed bin p/p moments [Pa²], in original graph-port order.
    pub fn pressure_second_moments(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.pressure)
    }
    /// Summed bin U/U moments [m⁶/s²], retaining graph-inflow signs/phases.
    pub fn flow_second_moments(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.flow)
    }
    /// Summed bin p/U moments [W]; real diagonal is signed port work.
    pub fn pressure_flow_moments(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.pressure_flow)
    }
    pub fn receiver_fields(&self) -> Option<&[HarmonicFieldStatistics]> {
        self.valid.then_some(&self.fields)
    }
    pub fn network_work(&self) -> Option<NetworkSourceWork> {
        self.valid.then_some(self.work)
    }
    /// Already included in external_load_power_w; never an additional debit.
    pub fn exterior_surface_power_w(&self) -> Option<f64> {
        self.valid.then_some(self.exterior_surface_power_w)
    }

    /// Allocation-free original per-bin evaluations followed by compensated
    /// sums. `band_budget` bounds bins*(network allowance + field allowance +
    /// 3*ports² + 13*receivers + 6), not an average or pooled per-bin allowance.
    /// Every bin must independently pass its original work and physical gates.
    pub fn evaluate(
        &mut self,
        network_budget: NetworkCoherenceBudget,
        field_budget: HarmonicWorkBudget,
        band_budget: HarmonicWorkBudget,
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        self.valid = false;
        poll(gate)?;
        let reduced = add(
            add(mul(3, self.pressure.len())?, mul(13, self.fields.len())?)?,
            6,
        )?;
        work(
            self.band.bins.len(),
            add(
                network_budget.maximum_complex_products,
                field_budget.maximum_products,
            )?,
            reduced,
            band_budget,
        )?;
        for (bin, result) in self.band.bins.iter().zip(&mut self.per_bin) {
            poll(gate)?;
            bin.evaluate_into(network_budget, field_budget, result, gate)?;
        }
        for index in 0..self.pressure.len() {
            self.pressure[index] =
                reduce_checked(self.per_bin.iter().map(|bin| moment(bin, 0, index)), gate)?;
            self.flow[index] =
                reduce_checked(self.per_bin.iter().map(|bin| moment(bin, 1, index)), gate)?;
            self.pressure_flow[index] =
                reduce_checked(self.per_bin.iter().map(|bin| moment(bin, 2, index)), gate)?;
        }
        for (r, field) in self.fields.iter_mut().enumerate() {
            field.mean_square_pressure_pa2 = reduce(
                self.per_bin
                    .iter()
                    .map(|bin| C64::from_re(bin.fields[r].mean_square_pressure_pa2)),
                gate,
            )?
            .re;
            for i in 0..3 {
                field.complex_intensity_w_m2[i] = reduce(
                    self.per_bin
                        .iter()
                        .map(|bin| bin.fields[r].complex_intensity_w_m2[i]),
                    gate,
                )?;
                for j in 0..3 {
                    field.velocity_cross_m2_s2[i][j] = reduce(
                        self.per_bin
                            .iter()
                            .map(|bin| bin.fields[r].velocity_cross_m2_s2[i][j]),
                        gate,
                    )?;
                }
            }
        }
        let value = |project: fn(NetworkSourceWork) -> f64| {
            reduce_checked(
                self.per_bin.iter().map(|bin| {
                    bin.network
                        .work()
                        .map(|w| C64::from_re(project(w)))
                        .ok_or_else(|| bad("incomplete original bin work"))
                }),
                gate,
            )
            .map(|z| z.re)
        };
        self.work = NetworkSourceWork {
            source_power_w: value(|w| w.source_power_w)?,
            external_load_power_w: value(|w| w.external_load_power_w)?,
            network_power_w: value(|w| w.network_power_w)?,
            generator_power_w: value(|w| w.generator_power_w)?,
            source_dissipation_w: value(|w| w.source_dissipation_w)?,
        };
        self.exterior_surface_power_w = reduce(
            self.per_bin
                .iter()
                .map(|bin| C64::from_re(bin.exterior_surface_power_w)),
            gate,
        )?
        .re;
        poll(gate)?;
        self.valid = true;
        Ok(())
    }
}

fn moment(
    bin: &RadiatingSourceStatistics<'_, '_>,
    kind: usize,
    index: usize,
) -> Result<C64, ExteriorNetworkError> {
    let values = match kind {
        0 => bin.network.pressure_second_moments(),
        1 => bin.network.flow_second_moments(),
        _ => bin.network.pressure_flow_moments(),
    };
    values
        .and_then(|v| v.get(index))
        .copied()
        .ok_or_else(|| bad("incomplete original bin port moments"))
}
fn reduce_checked(
    values: impl Iterator<Item = Result<C64, ExteriorNetworkError>>,
    gate: &CancelGate,
) -> Result<C64, ExteriorNetworkError> {
    let mut error = None;
    let value = reduce(
        values.scan((), |_, item| match item {
            Ok(value) => Some(value),
            Err(failure) => {
                error = Some(failure);
                None
            }
        }),
        gate,
    )?;
    error.map_or(Ok(value), Err)
}
