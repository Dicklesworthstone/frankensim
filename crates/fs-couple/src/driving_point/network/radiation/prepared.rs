//! Fixed loaded source response and observations of its actual opening flows.
//! Source pressure/flow coordinates never become independent exterior powers.

use super::*;
use fs_bem::near_field::FirstOrder;
use fs_bem::observation_transfers::{
    BemObservation, HarmonicField, HarmonicFieldStatistics, HarmonicWorkBudget,
    PreparedHarmonicObservation, PreparedMicrophoneArray, PreparedRadiationObservation,
    PreparedRadiationOptions, ReceiverSpec, SourceId,
};
use fs_duct::network::{
    CoupledSourceBank, MAX_NETWORK_PORTS, NetworkCoherenceBudget, NetworkSourceCoordinate,
    NetworkSourceResponse, NetworkSourceStatistics, PortSource, PreparedNetworkSources,
    network_prepare_sources, network_prepare_sources_regional,
};

const FACTOR_ENTRIES: usize = MAX_NETWORK_PORTS * MAX_NETWORK_PORTS;

mod angular;
pub use angular::PreparedRadiatingFarField;

/// One immutable source/network/exterior binding at a common positive frequency.
/// Both the original source laws and the exterior geometry remain borrowed.
/// Preparation shares one BEM factorization across openings and one loaded duct
/// factorization across sources, including its retained-wave pole fallback.
/// Subsequent evaluations allocate nothing and perform no boundary solves.
pub struct PreparedRadiatingNetwork<'a> {
    model: &'a RadiatingNetwork,
    sources: &'a [PortSource],
    bank: Option<&'a CoupledSourceBank>,
    scalar_loads: &'a [PortTermination],
    loss: LossModel,
    boundary: CoupledTermination,
    network: PreparedNetworkSources,
    opening_indices: Vec<usize>,
    observer: PreparedHarmonicObservation,
    radiation: Option<PreparedRadiationObservation>,
    // Source-major, receiver-minor, per the canonical source coordinate's SI unit.
    receiver_columns: Vec<HarmonicField>,
}

impl RadiatingNetwork {
    /// Prepare an explicit replacement partition of the original graph's
    /// NON-opening ports. `sources`, `bank` and `scalar_loads` must cover every
    /// such port exactly once. No load from the legacy single-source driving
    /// point is silently retained or discarded; the complete replacement is
    /// supplied here. Opening labels, orientation, geometry and exterior gas
    /// remain those of this owner, with their full mutual radiation matrix.
    ///
    /// The existing field allowance additionally bounds 5*panels*openings +
    /// 8*openings*receivers + 4*sources*(ports+receivers) complex field entries
    /// during preparation. The duct owner retains its own fixed matrix limits.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_sources<'a>(
        &'a self,
        omega: f64,
        loss: LossModel,
        sources: &'a [PortSource],
        bank: Option<&'a CoupledSourceBank>,
        scalar_loads: &'a [PortTermination],
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingNetwork<'a>, ExteriorNetworkError> {
        self.prepare_sources_inner(omega, loss, sources, bank, scalar_loads, None, gate)
    }

    /// Additionally prepare far-field amplitudes and full-sphere radiated power
    /// from the SAME opening solutions that supply the complete radiation load.
    /// Requested unit directions are separate from the positive quadrature
    /// ladder. `radiation` independently bounds additional angular kernel work,
    /// storage and drive-specific power/moment convergence; no new solve occurs.
    /// The ordinary [`Self::prepare_sources`] path pays no angular cost.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_sources_with_radiation<'a>(
        &'a self,
        omega: f64,
        loss: LossModel,
        sources: &'a [PortSource],
        bank: Option<&'a CoupledSourceBank>,
        scalar_loads: &'a [PortTermination],
        directions: &[[f64; 3]],
        radiation: PreparedRadiationOptions,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingNetwork<'a>, ExteriorNetworkError> {
        self.prepare_sources_inner(
            omega,
            loss,
            sources,
            bank,
            scalar_loads,
            Some((directions, radiation)),
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_sources_inner<'a>(
        &'a self,
        omega: f64,
        loss: LossModel,
        sources: &'a [PortSource],
        bank: Option<&'a CoupledSourceBank>,
        scalar_loads: &'a [PortTermination],
        angular: Option<(&[[f64; 3]], PreparedRadiationOptions)>,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingNetwork<'a>, ExteriorNetworkError> {
        poll(gate)?;
        let count = sources
            .len()
            .checked_add(bank.map_or(0, |b| b.impedance.labels.len()))
            .ok_or_else(|| bad("source count overflow"))?;
        if !omega.is_finite()
            || omega <= 0.0
            || count == 0
            || count > MAX_NETWORK_PORTS
            || scalar_loads.len() > MAX_NETWORK_PORTS
        {
            return Err(bad(
                "prepared radiation needs bounded sources and a positive frequency",
            ));
        }
        let (n, p, r) = (
            self.surface.areas().len(),
            self.patches.len(),
            self.receivers.len(),
        );
        let ports = self
            .driving
            .network
            .nodes
            .iter()
            .filter(|node| matches!(node.boundary, NodeBoundary::Port { .. }))
            .count();
        // All dimensions have already passed their original fixed owner caps.
        let entries = 5 * n * p + 8 * p * r + 4 * count * (ports + r);
        let terms = n * n * n + 3 * p * n * n + 4 * p * n * r + 4 * count * p * r;
        if entries > self.budget.max_field_entries || terms > self.budget.max_solve_terms {
            return Err(bad(
                "complete prepared source/exterior basis exceeds its declared budget",
            ));
        }
        let (boundary, solutions) = self.exterior_basis(omega, gate)?;
        let network = match &self.media {
            Some(media) => network_prepare_sources_regional(
                self.driving.network(),
                media,
                omega,
                loss,
                sources,
                bank,
                scalar_loads,
                Some(&boundary),
            )?,
            None => network_prepare_sources(
                self.driving.network(),
                &self.gas,
                omega,
                loss,
                sources,
                bank,
                scalar_loads,
                Some(&boundary),
            )?,
        };
        poll(gate)?;
        let opening_indices = self
            .patches
            .iter()
            .map(|patch| {
                network
                    .ports()
                    .iter()
                    .position(|label| label == &patch.port)
                    .ok_or_else(|| bad("prepared opening label lost"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let source_ids: Vec<_> = self
            .patches
            .iter()
            .map(|patch| SourceId::ApertureVolumeFlow(patch.port.clone()))
            .collect();
        let traces: Vec<Vec<f64>> = self
            .traces
            .iter()
            .map(|trace| trace.iter().map(|value| value.re).collect())
            .collect();
        let receivers: Vec<_> = self
            .receivers
            .iter()
            .enumerate()
            .map(|(i, &position_m)| ReceiverSpec {
                id: format!("network-receiver-{i}"),
                position_m,
            })
            .collect();
        let owner = BemObservation::new(
            &self.surface,
            Medium {
                density: self.gas.density,
                sound_speed: self.gas.sound_speed,
            },
            &source_ids,
            &traces,
            &receivers,
            self.clearance_m,
        )?;
        // The owner verifies exact source traces, surface, medium and frequency.
        // It observes the SAME solutions that produced the full boundary above.
        let observer =
            owner.prepare_harmonic_from_solutions(omega, &solutions, || gate.is_requested())?;
        let radiation = angular
            .map(|(directions, options)| {
                owner.prepare_radiation_from_solutions(omega, &solutions, directions, options, gate)
            })
            .transpose()?;
        let mut prepared = PreparedRadiatingNetwork {
            model: self,
            sources,
            bank,
            scalar_loads,
            loss,
            boundary,
            network,
            opening_indices,
            observer,
            radiation,
            receiver_columns: Vec::with_capacity(count * r),
        };
        let mut field = vec![HarmonicField::default(); r];
        let mut openings = [C64::ZERO; MAX_NETWORK_PORTS];
        for source in 0..count {
            poll(gate)?;
            for (opening, &port) in prepared.opening_indices.iter().enumerate() {
                openings[opening] = -prepared.network.flow_transfer()[port * count + source];
            }
            prepared
                .observer
                .evaluate_coherent(&openings[..p], &mut field, || gate.is_requested())?;
            prepared.receiver_columns.extend_from_slice(&field);
        }
        poll(gate)?;
        Ok(prepared)
    }
}

impl<'a> PreparedRadiatingNetwork<'a> {
    /// Original immutable geometry, gas records and exterior source identity.
    pub fn model(&self) -> &'a RadiatingNetwork {
        self.model
    }
    /// Common angular frequency [rad/s], exp(-i omega t).
    pub fn omega(&self) -> f64 {
        self.network.omega()
    }
    /// Original interior propagation/loss law.
    pub fn loss_model(&self) -> LossModel {
        self.loss
    }
    /// Original scalar source laws and their default phasors.
    pub fn sources(&self) -> &'a [PortSource] {
        self.sources
    }
    /// Original full reciprocal generator impedance, including mutual terms.
    pub fn source_bank(&self) -> Option<&'a CoupledSourceBank> {
        self.bank
    }
    /// Complete explicit scalar-load part of this prepared partition.
    pub fn scalar_loads(&self) -> &'a [PortTermination] {
        self.scalar_loads
    }
    /// Exact opening order and full work-conjugate radiation load [Pa s/m³].
    pub fn boundary(&self) -> &CoupledTermination {
        &self.boundary
    }
    /// Canonical source labels with pressure [Pa] or graph-inflow [m³/s] units.
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        self.network.source_coordinates()
    }
    /// Original phasors reordered into the canonical source coordinates.
    pub fn default_excitations(&self) -> &[C64] {
        self.network.default_excitations()
    }
    /// Original graph-port order, shared by every returned network field/moment.
    pub fn ports(&self) -> &[String] {
        self.network.ports()
    }

    /// Allocate a reusable buffer bound to this exact preparation.
    pub fn response_buffer(&self) -> RadiatingSourceResponse<'_, 'a> {
        RadiatingSourceResponse {
            owner: self,
            valid: false,
            network: self.network.response_buffer(),
            outward_flows: vec![C64::ZERO; self.opening_indices.len()],
            fields: vec![HarmonicField::default(); self.model.receivers.len()],
            exterior_surface_power_w: 0.0,
        }
    }

    /// Apply one coherent peak source vector. Cold source-to-receiver columns
    /// avoid repeating even the opening-to-receiver matrix composition here.
    /// Partial network or receiver updates stay inaccessible after any refusal.
    pub fn drive_into(
        &self,
        excitation: &[C64],
        output: &mut RadiatingSourceResponse<'_, '_>,
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        output.valid = false;
        if !std::ptr::eq(self, output.owner) {
            return Err(bad("coherent buffer belongs to another preparation"));
        }
        poll(gate)?;
        self.network.drive_into(excitation, &mut output.network)?;
        let r = output.fields.len();
        for (receiver, field) in output.fields.iter_mut().enumerate() {
            poll(gate)?;
            field.pressure_pa = sum(excitation
                .iter()
                .enumerate()
                .map(|(s, &x)| self.receiver_columns[s * r + receiver].pressure_pa * x))?;
            for axis in 0..3 {
                field.particle_velocity_m_s[axis] =
                    sum(excitation.iter().enumerate().map(|(s, &x)| {
                        self.receiver_columns[s * r + receiver].particle_velocity_m_s[axis] * x
                    }))?;
            }
        }
        for (outward, &port) in output.outward_flows.iter_mut().zip(&self.opening_indices) {
            *outward = -output.network.fields.port_inflows[port];
        }
        output.exterior_surface_power_w = sum(self
            .opening_indices
            .iter()
            .map(|&port| C64::from_re(-output.network.fields.port_powers_w[port])))?
        .re;
        poll(gate)?;
        output.valid = true;
        Ok(())
    }

    /// Allocate reusable physical field and complete network statistics.
    pub fn statistics_buffer(&self) -> RadiatingSourceStatistics<'_, 'a> {
        RadiatingSourceStatistics {
            owner: self,
            valid: false,
            network: self.network.statistics_buffer(),
            fields: vec![HarmonicFieldStatistics::default(); self.model.receivers.len()],
            exterior_surface_power_w: 0.0,
        }
    }

    /// Source second moment C=B Bᴴ, factor-major in canonical mixed SI units.
    /// Each coherent factor drives the SAME loaded opening basis before any
    /// receiver product. All moments contain the peak time-average 1/2 once.
    /// Rank is at most 32; mapping additionally costs rank*openings*sources
    /// products (at most 32³), with cancellation checked between factors.
    /// The supplied budgets independently bound original duct and field work.
    #[allow(clippy::too_many_arguments)]
    pub fn coherence_into(
        &self,
        factors: &[C64],
        rank: usize,
        network_budget: NetworkCoherenceBudget,
        field_budget: HarmonicWorkBudget,
        output: &mut RadiatingSourceStatistics<'_, '_>,
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        output.valid = false;
        if !std::ptr::eq(self, output.owner) {
            return Err(bad("statistics buffer belongs to another preparation"));
        }
        poll(gate)?;
        self.network
            .coherence_into(factors, rank, network_budget, &mut output.network)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self.opening_factors(factors, rank, &mut scratch, gate)?;
        self.observer.evaluate_coherence_modes(
            openings,
            rank,
            field_budget,
            &mut output.fields,
            || gate.is_requested(),
        )?;
        self.finish_statistics(output, gate)
    }

    /// Explicit E[|x_i|²] in each source coordinate's squared SI units.
    /// This uses direct loaded response columns, not dense one-hot sources.
    /// Independent generators generally produce CORRELATED opening flows:
    /// their mapped factors, never diagonal opening powers, reach the observer.
    pub fn independent_into(
        &self,
        second_moments: &[f64],
        network_budget: NetworkCoherenceBudget,
        field_budget: HarmonicWorkBudget,
        output: &mut RadiatingSourceStatistics<'_, '_>,
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        output.valid = false;
        if !std::ptr::eq(self, output.owner) {
            return Err(bad("statistics buffer belongs to another preparation"));
        }
        poll(gate)?;
        self.network
            .independent_into(second_moments, network_budget, &mut output.network)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self.independent_openings(second_moments, &mut scratch, gate)?;
        self.observer.evaluate_coherence_modes(
            openings,
            second_moments.len(),
            field_budget,
            &mut output.fields,
            || gate.is_requested(),
        )?;
        self.finish_statistics(output, gate)
    }

    fn finish_statistics(
        &self,
        output: &mut RadiatingSourceStatistics<'_, '_>,
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        let moments = output
            .network
            .pressure_flow_moments()
            .ok_or_else(|| bad("incomplete network statistics"))?;
        let n = self.network.ports().len();
        output.exterior_surface_power_w = sum(self
            .opening_indices
            .iter()
            .map(|&port| C64::from_re(-moments[port * n + port].re)))?
        .re;
        poll(gate)?;
        output.valid = true;
        Ok(())
    }

    /// Prepare ideal first-order microphones of the original receiver p/u.
    /// Projection retains near-field pressure/velocity cross terms, direction,
    /// and all complex inter-microphone phases; no plane-wave substitution.
    pub fn microphones(
        &self,
        patterns: &[FirstOrder],
    ) -> Result<PreparedRadiatingMicrophones<'_, 'a>, ExteriorNetworkError> {
        Ok(PreparedRadiatingMicrophones {
            owner: self,
            array: self.observer.microphones(patterns)?,
        })
    }

    fn opening_factors<'s>(
        &self,
        factors: &[C64],
        rank: usize,
        scratch: &'s mut [C64; FACTOR_ENTRIES],
        gate: &CancelGate,
    ) -> Result<&'s [C64], ExteriorNetworkError> {
        let (m, p) = (self.source_coordinates().len(), self.opening_indices.len());
        if rank > MAX_NETWORK_PORTS
            || factors.len() != rank * m
            || factors.iter().any(|&z| !finite(z))
        {
            return Err(bad(
                "source coherence needs bounded finite complete factors",
            ));
        }
        poll(gate)?;
        for (mode, factor) in factors.chunks_exact(m).enumerate() {
            poll(gate)?;
            for (opening, &port) in self.opening_indices.iter().enumerate() {
                scratch[mode * p + opening] = -sum(factor
                    .iter()
                    .enumerate()
                    .map(|(s, &x)| self.network.flow_transfer()[port * m + s] * x))?;
            }
        }
        Ok(&scratch[..rank * p])
    }

    fn independent_openings<'s>(
        &self,
        moments: &[f64],
        scratch: &'s mut [C64; FACTOR_ENTRIES],
        gate: &CancelGate,
    ) -> Result<&'s [C64], ExteriorNetworkError> {
        let (m, p) = (self.source_coordinates().len(), self.opening_indices.len());
        if moments.len() != m || moments.iter().any(|&x| !x.is_finite() || x < 0.0) {
            return Err(bad(
                "independent source moments must be finite and nonnegative",
            ));
        }
        poll(gate)?;
        for (s, &moment) in moments.iter().enumerate() {
            poll(gate)?;
            let amplitude = moment.sqrt();
            for (opening, &port) in self.opening_indices.iter().enumerate() {
                let flow = -self.network.flow_transfer()[port * m + s].scale(amplitude);
                if !finite(flow) {
                    return Err(bad("independent opening flow overflow"));
                }
                scratch[s * p + opening] = flow;
            }
        }
        Ok(&scratch[..m * p])
    }
}

/// Coherent output bound to its exact prepared source/geometry/medium owner.
/// Getters are unavailable until a complete successful application.
pub struct RadiatingSourceResponse<'p, 'a> {
    owner: &'p PreparedRadiatingNetwork<'a>,
    valid: bool,
    network: NetworkSourceResponse,
    outward_flows: Vec<C64>,
    fields: Vec<HarmonicField>,
    exterior_surface_power_w: f64,
}
impl RadiatingSourceResponse<'_, '_> {
    pub fn is_valid(&self) -> bool {
        self.valid
    }
    /// Original signed generator/interface/loss accounting and graph fields.
    pub fn network(&self) -> Option<&NetworkSourceResponse> {
        self.valid.then_some(&self.network)
    }
    /// Actual physical opening flows, positive out of the graph [m³/s].
    pub fn outward_flows_m3_s(&self) -> Option<&[C64]> {
        self.valid.then_some(&self.outward_flows)
    }
    /// Exterior pressure [Pa] and particle velocity [m/s] at original receivers.
    pub fn receiver_fields(&self) -> Option<&[HarmonicField]> {
        self.valid.then_some(&self.fields)
    }
    /// Signed work delivered to just the full exterior opening boundary [W].
    /// Already included in network external-load power; never debit it again.
    pub fn exterior_surface_power_w(&self) -> Option<f64> {
        self.valid.then_some(self.exterior_surface_power_w)
    }
}

/// Reusable time/ensemble statistics. A refused or cancelled application
/// invalidates every getter, including already-computed network quantities.
pub struct RadiatingSourceStatistics<'p, 'a> {
    owner: &'p PreparedRadiatingNetwork<'a>,
    valid: bool,
    network: NetworkSourceStatistics,
    fields: Vec<HarmonicFieldStatistics>,
    exterior_surface_power_w: f64,
}
impl RadiatingSourceStatistics<'_, '_> {
    pub fn is_valid(&self) -> bool {
        self.valid
    }
    /// Full p-p, U-U and p-U moments and signed physical power ownership.
    pub fn network(&self) -> Option<&NetworkSourceStatistics> {
        self.valid.then_some(&self.network)
    }
    /// Mean-square pressure, velocity cross moments and complex intensity.
    pub fn receiver_fields(&self) -> Option<&[HarmonicFieldStatistics]> {
        self.valid.then_some(&self.fields)
    }
    /// Mean exterior opening work [W], already included in network load power.
    pub fn exterior_surface_power_w(&self) -> Option<f64> {
        self.valid.then_some(self.exterior_surface_power_w)
    }
}

/// First-order microphone array bound to the same loaded source response.
/// Returned row-major spectra are 0.5 E[m_i conj(m_j)] [Pa²]. The underlying
/// BEM projection and spectrum owner is reused without independent source powers.
pub struct PreparedRadiatingMicrophones<'p, 'a> {
    owner: &'p PreparedRadiatingNetwork<'a>,
    array: PreparedMicrophoneArray<'p>,
}
impl PreparedRadiatingMicrophones<'_, '_> {
    pub fn coherence(
        &mut self,
        factors: &[C64],
        rank: usize,
        budget: HarmonicWorkBudget,
        gate: &CancelGate,
    ) -> Result<&[C64], ExteriorNetworkError> {
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .opening_factors(factors, rank, &mut scratch, gate)?;
        Ok(self
            .array
            .evaluate_coherence_modes(openings, rank, budget, || gate.is_requested())?)
    }
    /// Independent generators retain the induced cross-opening coherence.
    pub fn independent(
        &mut self,
        moments: &[f64],
        budget: HarmonicWorkBudget,
        gate: &CancelGate,
    ) -> Result<&[C64], ExteriorNetworkError> {
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .independent_openings(moments, &mut scratch, gate)?;
        Ok(self
            .array
            .evaluate_coherence_modes(openings, moments.len(), budget, || gate.is_requested())?)
    }
}

// Compensated complex reduction preserves weak residuals after phase cancellation.
fn sum(values: impl Iterator<Item = C64>) -> Result<C64, ExteriorNetworkError> {
    fn add(sum: &mut f64, correction: &mut f64, value: f64) {
        let next = *sum + value;
        *correction += if sum.abs() >= value.abs() {
            (*sum - next) + value
        } else {
            (value - next) + *sum
        };
        *sum = next;
    }
    let (mut value, mut correction) = (C64::ZERO, C64::ZERO);
    for item in values {
        add(&mut value.re, &mut correction.re, item.re);
        add(&mut value.im, &mut correction.im, item.im);
    }
    let result = value + correction;
    if finite(result) {
        Ok(result)
    } else {
        Err(bad("prepared physical field/work overflow"))
    }
}
