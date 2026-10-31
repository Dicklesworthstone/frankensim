//! Fixed-boundary response and explicitly supplied source second moments.
//!
//! A source factor is a vector of peak pressure/flow coordinates, not a power
//! vector. Its complete field is combined before any quadratic observation.
//! Different factors have unit, mutually uncorrelated latent coefficients:
//! E[x x^H] = B B^H. No phases, correlations, bandwidth or RMS conversion are
//! inferred. The existing harmonic time-average factor 1/2 appears once.

use super::*;

/// Physical coordinate of one prepared source. Mixed-coordinate covariances
/// have the product of the two indicated SI units; they are not dimensionless.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NetworkSourceQuantity {
    /// Prescribed peak flow into the graph [m³/s].
    VolumeFlow,
    /// Peak open-circuit pressure, before the entire source impedance [Pa].
    OpenCircuitPressure,
}

/// Stable labeled source coordinate, sorted lexicographically by label.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NetworkSourceCoordinate {
    /// Exact external source-port label.
    pub label: String,
    /// Physical quantity and units of this input coordinate.
    pub quantity: NetworkSourceQuantity,
}

/// Bounds on one statistical application. Source count and rank are also
/// independently capped at 32; storage never depends on an unbounded rank.
#[derive(Clone, Copy, Debug)]
pub struct NetworkCoherenceBudget {
    /// Largest admitted number of supplied factors (or independent sources).
    pub maximum_rank: usize,
    /// Maximum complex products in field, moment and work contractions.
    /// Dense factors cost rank*(4*ports*sources + 3*ports² + 6*ports);
    /// independent sources replace 4*ports*sources by 4*ports per source.
    pub maximum_complex_products: usize,
}

impl Default for NetworkCoherenceBudget {
    fn default() -> Self {
        Self {
            maximum_rank: MAX_NETWORK_PORTS,
            maximum_complex_products: 262_144,
        }
    }
}

/// Mean signed work for a coherent vector or source second moment [W].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkSourceWork {
    /// Work delivered at all source/network interfaces, excluding source loss.
    pub source_power_w: f64,
    /// Work delivered to the complete scalar/coupled exterior loads.
    pub external_load_power_w: f64,
    /// Net work entering the graph (wall, gas and internal nodal loads).
    pub network_power_w: f64,
    /// Signed generator work before internal source impedances.
    pub generator_power_w: f64,
    /// Complete Hermitian source loss; passive up to admitted roundoff.
    pub source_dissipation_w: f64,
}

/// Complete source-to-port response for one fixed frequency and boundary.
/// Preparation shares one loaded boundary inverse across every source RHS.
/// Application performs bounded contractions only: no solve, eigenanalysis,
/// impedance admission, graph traversal or heap allocation is repeated.
#[derive(Clone, Debug)]
pub struct PreparedNetworkSources {
    omega: f64,
    ports: Vec<String>,
    coordinates: Vec<NetworkSourceCoordinate>,
    defaults: Vec<C64>,
    source_ports: Vec<usize>,
    source_index: Vec<Option<usize>>,
    canonical: Vec<usize>,
    pressure: Vec<C64>,
    flow: Vec<C64>,
    internal: Vec<C64>,
    dissipative: Vec<C64>,
}

impl NetworkResponse {
    /// Prepare all source columns from one loaded boundary factorization.
    /// The supplied excitations become the default vector in canonical source
    /// order. Changing an impedance, medium, topology or frequency requires a
    /// new preparation; changing only the complex source vector does not.
    pub fn prepare_sources(
        &self,
        sources: &[PortSource],
        bank: Option<&CoupledSourceBank>,
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<PreparedNetworkSources, DuctError> {
        let labels: Vec<&str> = self.ports.iter().map(String::as_str).collect();
        let sources = source_conditions(&labels, sources, bank, self.omega)?;
        let system = self.source_system(&sources, loads, coupled, &self.gas)?;
        let driven: Vec<usize> = system
            .canonical
            .iter()
            .copied()
            .filter(|&port| sources.conditions[port].is_some())
            .collect();
        let (n, count, m) = (self.ports.len(), driven.len(), system.live.len());
        let mut flow = vec![C64::ZERO; n * count];
        let mut pressure = flow.clone();
        for (column, &source) in driven.iter().enumerate() {
            if sources.conditions[source]
                .and_then(PortExcitation::flow)
                .is_some()
            {
                flow[source * count + column] = C64::ONE;
                for (row, &port) in system.live.iter().enumerate() {
                    flow[port * count + column] =
                        -coherent_sum(
                            system.live.iter().enumerate().map(|(j, &q)| {
                                system.inverse[row * m + j] * self.z[q * n + source]
                            }),
                        );
                }
            } else {
                let rhs = system
                    .live
                    .iter()
                    .position(|&port| port == source)
                    .expect("pressure-source current retained");
                for (row, &port) in system.live.iter().enumerate() {
                    flow[port * count + column] = system.inverse[row * m + rhs];
                }
            }
            for port in 0..n {
                pressure[port * count + column] = coherent_sum(
                    system
                        .canonical
                        .iter()
                        .map(|&other| self.z[port * n + other] * flow[other * count + column]),
                );
            }
        }
        PreparedNetworkSources::from_transfers(self.omega, &self.ports, &sources, pressure, flow)
    }
}

impl PreparedNetworkSources {
    pub(super) fn from_transfers(
        omega: f64,
        ports: &[String],
        sources: &ResolvedSources<'_>,
        pressure: Vec<C64>,
        flow: Vec<C64>,
    ) -> Result<Self, DuctError> {
        let mut canonical: Vec<usize> = (0..ports.len()).collect();
        canonical.sort_unstable_by(|&a, &b| ports[a].cmp(&ports[b]));
        let source_ports: Vec<usize> = canonical
            .iter()
            .copied()
            .filter(|&port| sources.conditions[port].is_some())
            .collect();
        let count = source_ports.len();
        let mut coordinates = Vec::with_capacity(count);
        let mut defaults = Vec::with_capacity(count);
        let mut source_index = vec![None; ports.len()];
        for (index, &port) in source_ports.iter().enumerate() {
            let (quantity, amplitude) = match sources.conditions[port].expect("source selected") {
                PortExcitation::VolumeFlow(value) => (NetworkSourceQuantity::VolumeFlow, value),
                PortExcitation::Pressure {
                    open_circuit_pressure_pa,
                    ..
                } => (
                    NetworkSourceQuantity::OpenCircuitPressure,
                    open_circuit_pressure_pa,
                ),
            };
            coordinates.push(NetworkSourceCoordinate {
                label: ports[port].clone(),
                quantity,
            });
            defaults.push(amplitude);
            source_index[port] = Some(index);
        }
        // Retain physical source effort and its Hermitian part as transfers.
        // These contractions happen once, not for every realization/factor.
        let mut internal = vec![C64::ZERO; ports.len() * count];
        let mut dissipative = internal.clone();
        for &port in &source_ports {
            for column in 0..count {
                internal[port * count + column] =
                    coherent_sum(canonical.iter().filter_map(|&other| {
                        sources
                            .impedance(port, other)
                            .map(|z| z * flow[other * count + column])
                    }));
                dissipative[port * count + column] =
                    coherent_sum(canonical.iter().filter_map(|&other| {
                        let forward = sources.impedance(port, other)?;
                        let reverse = sources
                            .impedance(other, port)
                            .expect("same source operator");
                        let hermitian = forward.scale(0.5) + reverse.conj().scale(0.5);
                        Some(hermitian * flow[other * count + column])
                    }));
            }
        }
        if pressure
            .iter()
            .chain(&flow)
            .chain(&internal)
            .chain(&dissipative)
            .any(|&v| !finite(v))
        {
            return Err(DuctError::Singular);
        }
        Ok(Self {
            omega,
            ports: ports.to_vec(),
            coordinates,
            defaults,
            source_ports,
            source_index,
            canonical,
            pressure,
            flow,
            internal,
            dissipative,
        })
    }

    /// Fixed angular frequency [rad/s]; every input uses exp(-i omega t).
    #[must_use]
    pub const fn omega(&self) -> f64 {
        self.omega
    }

    /// Physical output-port order, as in the original graph.
    #[must_use]
    pub fn ports(&self) -> &[String] {
        &self.ports
    }

    /// Canonical input order and explicit SI quantity for each coordinate.
    #[must_use]
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        &self.coordinates
    }

    /// Graph-port index belonging to each canonical source coordinate.
    #[must_use]
    pub fn source_port_indices(&self) -> &[usize] {
        &self.source_ports
    }

    /// Original source amplitudes, reordered into [`Self::source_coordinates`].
    #[must_use]
    pub fn default_excitations(&self) -> &[C64] {
        &self.defaults
    }

    /// Row-major pressure transfer: Pa per input coordinate's SI unit.
    #[must_use]
    pub fn pressure_transfer(&self) -> &[C64] {
        &self.pressure
    }

    /// Row-major graph-inflow transfer: m³/s per input coordinate's SI unit.
    #[must_use]
    pub fn flow_transfer(&self) -> &[C64] {
        &self.flow
    }

    /// Allocate one reusable coherent result. Application writes it only after
    /// all inputs, fields and work have passed finite-value checks.
    #[must_use]
    pub fn response_buffer(&self) -> NetworkSourceResponse {
        let n = self.ports.len();
        NetworkSourceResponse {
            fields: NetworkCoherentResponse {
                omega: self.omega,
                ports: self.ports.clone(),
                port_pressures: vec![C64::ZERO; n],
                port_inflows: vec![C64::ZERO; n],
                port_powers_w: vec![0.0; n],
                source_power_w: 0.0,
                external_load_power_w: 0.0,
                network_power_w: 0.0,
            },
            generator_powers_w: vec![0.0; n],
            source_dissipations_w: vec![0.0; n],
            generator_power_w: 0.0,
            source_dissipation_w: 0.0,
        }
    }

    /// Apply one coherent peak vector without allocation or refactorization.
    /// Coordinates are in `source_coordinates` order, never source-list order.
    /// Refusal preserves the previous coherent result.
    pub fn drive_into(
        &self,
        excitation: &[C64],
        output: &mut NetworkSourceResponse,
    ) -> Result<(), DuctError> {
        let n = self.ports.len();
        if output.fields.omega.to_bits() != self.omega.to_bits()
            || output.fields.ports != self.ports
            || output.fields.port_pressures.len() != n
            || output.fields.port_inflows.len() != n
            || output.fields.port_powers_w.len() != n
            || output.generator_powers_w.len() != n
            || output.source_dissipations_w.len() != n
        {
            return Err(bad(
                "prepared coherent result has different ports or storage",
            ));
        }
        self.validate_excitation(excitation)?;
        let sample = self.sample(excitation, None)?;
        let mut port_work = [0.0; MAX_NETWORK_PORTS];
        let mut generator = port_work;
        let mut internal = port_work;
        for port in 0..n {
            port_work[port] = mean_work(sample.pressure[port], sample.flow[port]);
            generator[port] = self.generator_work(port, &sample, excitation);
            internal[port] = mean_work(sample.internal[port], sample.flow[port]);
        }
        let dissipation = coherent_sum(
            self.canonical
                .iter()
                .map(|&port| C64::from_re(mean_work(sample.dissipative[port], sample.flow[port]))),
        )
        .re;
        let work = self.summarize(&port_work, &generator, dissipation)?;
        if internal[..n].iter().any(|value| !value.is_finite()) {
            return Err(DuctError::Singular);
        }
        output
            .fields
            .port_pressures
            .copy_from_slice(&sample.pressure[..n]);
        output
            .fields
            .port_inflows
            .copy_from_slice(&sample.flow[..n]);
        output.fields.port_powers_w.copy_from_slice(&port_work[..n]);
        output.generator_powers_w.copy_from_slice(&generator[..n]);
        output.source_dissipations_w.copy_from_slice(&internal[..n]);
        output.fields.source_power_w = work.source_power_w;
        output.fields.external_load_power_w = work.external_load_power_w;
        output.fields.network_power_w = work.network_power_w;
        output.generator_power_w = work.generator_power_w;
        output.source_dissipation_w = work.source_dissipation_w;
        Ok(())
    }

    /// Allocate reusable storage for all output-port cross moments and work.
    #[must_use]
    pub fn statistics_buffer(&self) -> NetworkSourceStatistics {
        let n = self.ports.len();
        NetworkSourceStatistics {
            omega: self.omega,
            ports: self.ports.clone(),
            valid: false,
            pressure: Moments::new(n * n),
            flow: Moments::new(n * n),
            pressure_flow: Moments::new(n * n),
            generator: Moments::new(n),
            internal: Moments::new(n),
            dissipation: Moments::new(1),
            work: NetworkSourceWork::default(),
        }
    }

    /// Contract C=B B^H for caller-supplied source factors. Storage is
    /// factor-major: each consecutive `source_coordinates().len()` entries
    /// is one full coherent source vector. Distinct factors are uncorrelated;
    /// within each factor every phase and cross term is retained. A single
    /// factor recovers deterministic coherent statistics. All returned moments
    /// include the peak-phasor time average 1/2 and are spectral-line quantities,
    /// not power spectral densities per Hz. No allocation occurs here.
    /// Any refusal invalidates the statistics buffer and its getters.
    pub fn coherence_into(
        &self,
        factors: &[C64],
        rank: usize,
        budget: NetworkCoherenceBudget,
        output: &mut NetworkSourceStatistics,
    ) -> Result<(), DuctError> {
        output.valid = false;
        self.admit_statistics(output, rank, false, budget)?;
        let count = self.coordinates.len();
        if factors.len() != count * rank || factors.iter().any(|&v| !finite(v)) {
            return Err(bad("coherence factors need finite complete source vectors"));
        }
        output.clear();
        for excitation in factors.chunks_exact(count) {
            let sample = self.sample(excitation, None)?;
            self.accumulate(&sample, excitation, output)?;
        }
        self.finish_statistics(output)
    }

    /// Explicitly independent sources with E[x_i conj(x_j)]=0 for i!=j.
    /// Each nonnegative input is E[|x_i|²] in that coordinate's squared SI
    /// units, before the common time-average factor 1/2. This direct-column
    /// path avoids multiplying dense one-hot factors. It does NOT remove
    /// source/load coupling: each response column includes the entire boundary.
    pub fn independent_into(
        &self,
        second_moments: &[f64],
        budget: NetworkCoherenceBudget,
        output: &mut NetworkSourceStatistics,
    ) -> Result<(), DuctError> {
        output.valid = false;
        let count = self.coordinates.len();
        self.admit_statistics(output, count, true, budget)?;
        if second_moments.len() != count
            || second_moments.iter().any(|&v| !v.is_finite() || v < 0.0)
        {
            return Err(bad(
                "independent source moments must be finite and nonnegative",
            ));
        }
        output.clear();
        let mut excitation = [C64::ZERO; MAX_NETWORK_PORTS];
        for (column, &moment) in second_moments.iter().enumerate() {
            if moment == 0.0 {
                continue;
            }
            excitation[column] = C64::from_re(moment.sqrt());
            let sample = self.sample(&excitation[..count], Some(column))?;
            self.accumulate(&sample, &excitation[..count], output)?;
            excitation[column] = C64::ZERO;
        }
        self.finish_statistics(output)
    }

    fn validate_excitation(&self, excitation: &[C64]) -> Result<(), DuctError> {
        if excitation.len() != self.coordinates.len() || excitation.iter().any(|&v| !finite(v)) {
            Err(bad(
                "prepared excitation needs every finite source coordinate",
            ))
        } else {
            Ok(())
        }
    }

    fn sample(&self, excitation: &[C64], independent: Option<usize>) -> Result<Sample, DuctError> {
        let mut sample = Sample::default();
        let count = self.coordinates.len();
        let contract = |matrix: &[C64], port: usize| {
            independent.map_or_else(
                || coherent_sum((0..count).map(|j| matrix[port * count + j] * excitation[j])),
                |j| matrix[port * count + j] * excitation[j],
            )
        };
        for port in 0..self.ports.len() {
            sample.pressure[port] = contract(&self.pressure, port);
            sample.flow[port] = contract(&self.flow, port);
            sample.internal[port] = contract(&self.internal, port);
            sample.dissipative[port] = contract(&self.dissipative, port);
        }
        if sample
            .pressure
            .iter()
            .chain(&sample.flow)
            .chain(&sample.internal)
            .chain(&sample.dissipative)
            .any(|&v| !finite(v))
        {
            return Err(DuctError::Singular);
        }
        Ok(sample)
    }

    fn generator_work(&self, port: usize, sample: &Sample, excitation: &[C64]) -> f64 {
        self.source_index[port].map_or(0.0, |source| {
            let pressure = match self.coordinates[source].quantity {
                NetworkSourceQuantity::VolumeFlow => sample.pressure[port],
                NetworkSourceQuantity::OpenCircuitPressure => excitation[source],
            };
            mean_work(pressure, sample.flow[port])
        })
    }

    fn summarize(
        &self,
        port_work: &[f64],
        generator: &[f64],
        dissipation: f64,
    ) -> Result<NetworkSourceWork, DuctError> {
        let sum = |values: &[f64], driven: Option<bool>| {
            coherent_sum(self.canonical.iter().filter_map(|&port| {
                driven
                    .is_none_or(|value| self.source_index[port].is_some() == value)
                    .then(|| C64::from_re(values[port]))
            }))
            .re
        };
        let work = NetworkSourceWork {
            source_power_w: sum(port_work, Some(true)),
            external_load_power_w: -sum(port_work, Some(false)),
            network_power_w: sum(port_work, None),
            generator_power_w: sum(generator, None),
            source_dissipation_w: dissipation,
        };
        if port_work[..self.ports.len()]
            .iter()
            .chain(&generator[..self.ports.len()])
            .any(|v| !v.is_finite())
            || [
                work.source_power_w,
                work.external_load_power_w,
                work.network_power_w,
                work.generator_power_w,
                work.source_dissipation_w,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(DuctError::Singular);
        }
        Ok(work)
    }

    fn admit_statistics(
        &self,
        output: &NetworkSourceStatistics,
        rank: usize,
        independent: bool,
        budget: NetworkCoherenceBudget,
    ) -> Result<(), DuctError> {
        if output.omega.to_bits() != self.omega.to_bits() || output.ports != self.ports {
            return Err(bad("prepared statistics have different frequency or ports"));
        }
        if rank == 0 || rank > MAX_NETWORK_PORTS || rank > budget.maximum_rank {
            return Err(bad(
                "source coherence rank exceeds the explicit 1..=32 budget",
            ));
        }
        let n = self.ports.len();
        let columns = if independent {
            1
        } else {
            self.coordinates.len()
        };
        let products = rank * (4 * n * columns + 3 * n * n + 6 * n);
        if products > budget.maximum_complex_products || products > 262_144 {
            return Err(bad("source coherence contraction exceeds its work budget"));
        }
        Ok(())
    }

    fn accumulate(
        &self,
        sample: &Sample,
        excitation: &[C64],
        output: &mut NetworkSourceStatistics,
    ) -> Result<(), DuctError> {
        let n = self.ports.len();
        for row in 0..n {
            for column in 0..n {
                let index = row * n + column;
                output.pressure.add(
                    index,
                    (sample.pressure[row].scale(0.5)) * sample.pressure[column].conj(),
                )?;
                output.flow.add(
                    index,
                    sample.flow[row].scale(0.5) * sample.flow[column].conj(),
                )?;
                output.pressure_flow.add(
                    index,
                    sample.pressure[row].scale(0.5) * sample.flow[column].conj(),
                )?;
            }
            output.generator.add(
                row,
                C64::from_re(self.generator_work(row, sample, excitation)),
            )?;
            output.internal.add(
                row,
                C64::from_re(mean_work(sample.internal[row], sample.flow[row])),
            )?;
        }
        let dissipated = coherent_sum(
            self.canonical
                .iter()
                .map(|&port| C64::from_re(mean_work(sample.dissipative[port], sample.flow[port]))),
        );
        output.dissipation.add(0, dissipated)
    }

    fn finish_statistics(&self, output: &mut NetworkSourceStatistics) -> Result<(), DuctError> {
        for moments in [
            &mut output.pressure,
            &mut output.flow,
            &mut output.pressure_flow,
            &mut output.generator,
            &mut output.internal,
            &mut output.dissipation,
        ] {
            moments.finish()?;
        }
        let mut port_work = [0.0; MAX_NETWORK_PORTS];
        let mut generator = port_work;
        for port in 0..self.ports.len() {
            port_work[port] = output.pressure_flow.values[port * self.ports.len() + port].re;
            generator[port] = output.generator.values[port].re;
        }
        output.work = self.summarize(&port_work, &generator, output.dissipation.values[0].re)?;
        output.valid = true;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Sample {
    pressure: [C64; MAX_NETWORK_PORTS],
    flow: [C64; MAX_NETWORK_PORTS],
    internal: [C64; MAX_NETWORK_PORTS],
    dissipative: [C64; MAX_NETWORK_PORTS],
}

impl Default for Sample {
    fn default() -> Self {
        Self {
            pressure: [C64::ZERO; MAX_NETWORK_PORTS],
            flow: [C64::ZERO; MAX_NETWORK_PORTS],
            internal: [C64::ZERO; MAX_NETWORK_PORTS],
            dissipative: [C64::ZERO; MAX_NETWORK_PORTS],
        }
    }
}

fn mean_work(pressure: C64, flow: C64) -> f64 {
    (pressure.scale(0.5) * flow.conj()).re
}

#[derive(Clone, Debug)]
struct Moments {
    values: Vec<C64>,
    correction: Vec<C64>,
}

impl Moments {
    fn new(count: usize) -> Self {
        Self {
            values: vec![C64::ZERO; count],
            correction: vec![C64::ZERO; count],
        }
    }
    fn clear(&mut self) {
        self.values.fill(C64::ZERO);
        self.correction.fill(C64::ZERO);
    }
    fn add(&mut self, index: usize, value: C64) -> Result<(), DuctError> {
        if !finite(value) {
            return Err(DuctError::Singular);
        }
        let sum = self.values[index];
        let next = sum + value;
        let correction = &mut self.correction[index];
        correction.re += if sum.re.abs() >= value.re.abs() {
            (sum.re - next.re) + value.re
        } else {
            (value.re - next.re) + sum.re
        };
        correction.im += if sum.im.abs() >= value.im.abs() {
            (sum.im - next.im) + value.im
        } else {
            (value.im - next.im) + sum.im
        };
        self.values[index] = next;
        Ok(())
    }
    fn finish(&mut self) -> Result<(), DuctError> {
        for (value, correction) in self.values.iter_mut().zip(&self.correction) {
            *value = *value + *correction;
            if !finite(*value) {
                return Err(DuctError::Singular);
            }
        }
        Ok(())
    }
}

/// Reusable full port second moments and signed work. Every matrix is row-major
/// in `ports()` order and includes the harmonic time-average factor 1/2.
/// Getters return `None` until a successful application, and after any refusal;
/// an incomplete contraction can never be observed as a valid statistic.
#[derive(Clone, Debug)]
pub struct NetworkSourceStatistics {
    omega: f64,
    ports: Vec<String>,
    valid: bool,
    pressure: Moments,
    flow: Moments,
    pressure_flow: Moments,
    generator: Moments,
    internal: Moments,
    dissipation: Moments,
    work: NetworkSourceWork,
}

impl NetworkSourceStatistics {
    /// Original graph's physical port order.
    #[must_use]
    pub fn ports(&self) -> &[String] {
        &self.ports
    }
    /// Common angular frequency [rad/s].
    #[must_use]
    pub const fn omega(&self) -> f64 {
        self.omega
    }
    /// Whether the most recent statistical contraction completed successfully.
    #[must_use]
    pub const fn is_valid(&self) -> bool {
        self.valid
    }
    /// 0.5 E[p_i conj(p_j)] [Pa²], including all cross-port/source terms.
    #[must_use]
    pub fn pressure_second_moments(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.pressure.values)
    }
    /// 0.5 E[U_i conj(U_j)] [m⁶/s²], with U into the graph.
    #[must_use]
    pub fn flow_second_moments(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.flow.values)
    }
    /// 0.5 E[p_i conj(U_j)] [W]; the real diagonal is signed graph-port work.
    #[must_use]
    pub fn pressure_flow_moments(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.pressure_flow.values)
    }
    /// Signed mean generator work at one port [W], zero at passive-load ports.
    #[must_use]
    pub fn generator_power_w(&self, port: usize) -> Option<f64> {
        self.valid
            .then(|| self.generator.values.get(port).map(|value| value.re))
            .flatten()
    }
    /// Signed internal-source work at one port [W], including mutual exchange.
    #[must_use]
    pub fn source_internal_power_w(&self, port: usize) -> Option<f64> {
        self.valid
            .then(|| self.internal.values.get(port).map(|value| value.re))
            .flatten()
    }
    /// Complete mean generator, source-loss, interface, exterior and graph work.
    #[must_use]
    pub fn work(&self) -> Option<NetworkSourceWork> {
        self.valid.then_some(self.work)
    }
    fn clear(&mut self) {
        for moments in [
            &mut self.pressure,
            &mut self.flow,
            &mut self.pressure_flow,
            &mut self.generator,
            &mut self.internal,
            &mut self.dissipation,
        ] {
            moments.clear();
        }
    }
}
