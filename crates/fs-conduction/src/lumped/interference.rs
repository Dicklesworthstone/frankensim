//! Causal linear thermal response to explicitly bound interference power.
//!
//! Each original [LumpedNode](super::LumpedNode) obeys
//! C dT/dt = P(t) - G (T - T_ambient). The complete finite line spectrum
//! drives that equation, including beat and carrier-sum terms; no sampling,
//! carrier averaging, time stepping or thermal state mutation occurs here.
//!
//! A binding is an explicit modelling declaration of APPLIED THERMAL POWER.
//! A general quadratic spectrum can be signed. Acoustic intensity, reversible
//! field storage and mechanical work are not automatically dissipative heat.
//! When a spectrum describes an ensemble mean, this linear solve returns the
//! mean temperature, not temperature fluctuations or a noise realization.
//!
//! This is the original constant-capacity, constant-conductance, independent-
//! node model. The original Biot gate remains mandatory. Radiation, nonlinear
//! material response and exchange between nodes require a different model.

/// Original phase/transport enthalpy evolution driven by prescribed interference power.
pub mod enthalpy;

use fs_fft::harmonic_window::{PreparedHarmonicProductSpectrum, ProductSpectrumError};
use fs_math::dd::Dd;

use super::{BiotGate, LumpedNetwork, LumpedNode, MAX_LUMPED_NODES, ValidityVerdict, lumped_error};
use crate::ConductionError;

/// Hard bound on total per-query spectral work, counting DC once per node.
///
/// Reusing a spectrum at several nodes counts it several times because each
/// node can have a different thermal decay rate.
pub const MAX_LUMPED_INTERFERENCE_TERMS: usize = 1_048_576;

/// One explicitly authored conversion from a product spectrum to node watts.
///
/// Labels and spectra are borrowed. The input_unit must describe the actual
/// spectrum output units, and the finite positive scale must convert one of
/// those units to watts. For a spatially uniform W/m³ spectrum, for example,
/// the scale is the body's volume in m³. No spatial integration, absorption
/// fraction or unit conversion is inferred.
///
/// The source_clock names the absolute time coordinate of the spectrum's
/// authored phasors. All bindings in one preparation must declare the same
/// clock. Labels are assertions by the caller, not proofs of common phase
/// origin or dimensional compatibility; those remain the source owner's job.
#[derive(Clone, Copy, Debug)]
pub struct LumpedPowerSpectrumBinding<'a> {
    pub node_name: &'a str,
    pub spectrum: &'a PreparedHarmonicProductSpectrum,
    pub input_unit: &'a str,
    pub watts_per_spectrum_unit: f64,
    pub source_clock: &'a str,
}

/// Bounds admitted before the response owner or any response is returned.
#[derive(Clone, Copy, Debug)]
pub struct LumpedInterferenceLimits {
    pub maximum_nodes: usize,
    pub maximum_spectral_terms: usize,
    /// Combined logical vector payload of the prepared node bindings and ONE
    /// response. Borrowed networks, spectra, labels, caller-owned inputs,
    /// previously returned responses and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

#[derive(Clone, Copy, Debug)]
struct PreparedNodePower<'a> {
    node: &'a LumpedNode,
    binding: LumpedPowerSpectrumBinding<'a>,
    decay_rate_per_s: f64,
}

/// Biot-gated response owner tied to one original immutable thermal network.
///
/// Preparation maps exactly one named power binding to every original node.
/// Queries retain the shared absolute source clock, use the spectrum owner's
/// allocation-free causal integrals, and allocate just one node-result vector.
/// Source spectra are never cloned, re-solved, flattened to DC or modified.
#[derive(Debug)]
pub struct PreparedLumpedInterference<'a> {
    network: &'a LumpedNetwork,
    nodes: Vec<PreparedNodePower<'a>>,
    source_clock: &'a str,
    gate: BiotGate,
    verdict: ValidityVerdict,
    spectral_terms: usize,
    limits: LumpedInterferenceLimits,
    prepared_bytes: usize,
}

impl<'a> PreparedLumpedInterference<'a> {
    pub fn new(
        network: &'a LumpedNetwork,
        gate: BiotGate,
        bindings: &[LumpedPowerSpectrumBinding<'a>],
        limits: LumpedInterferenceLimits,
    ) -> Result<Self, ConductionError> {
        Self::new_cancellable(network, gate, bindings, limits, || false)
    }

    /// Cancellation or any refusal discards preparation atomically.
    ///
    /// The original gate is applied before admitting a response owner. In
    /// addition, this physical thermal consumer requires a positive absolute
    /// ambient temperature and finite derived Biot/decay values.
    pub fn new_cancellable(
        network: &'a LumpedNetwork,
        gate: BiotGate,
        bindings: &[LumpedPowerSpectrumBinding<'a>],
        limits: LumpedInterferenceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        poll(&mut cancelled, "lumped interference preparation", 0)?;
        let count = network.nodes().len();
        if count > limits.maximum_nodes || count > MAX_LUMPED_NODES {
            return Err(lumped_error(
                "lumped interference",
                "node budget exceeded".to_string(),
            ));
        }
        if bindings.len() != count {
            return Err(lumped_error(
                "lumped interference",
                format!(
                    "{} power bindings for {count} original nodes",
                    bindings.len()
                ),
            ));
        }
        if network.ambient_k() <= 0.0 || !network.ambient_k().is_finite() {
            return Err(lumped_error(
                "lumped interference",
                "ambient must be a finite positive absolute temperature in K".to_string(),
            ));
        }
        for (index, node) in network.nodes().iter().enumerate() {
            poll(&mut cancelled, "lumped interference preparation", index)?;
            if !node.biot().is_finite() {
                return Err(lumped_error(
                    node.name(),
                    "derived Biot number is not finite".to_string(),
                ));
            }
        }
        let verdict = gate.adjudicate(network);
        if let ValidityVerdict::Refused {
            node,
            biot,
            ceiling,
        } = &verdict
        {
            return Err(lumped_error(
                node,
                format!(
                    "Biot {biot} exceeds the declared ceiling {ceiling}; \
                     escalate to the full transient rung"
                ),
            ));
        }
        let result_bytes = vector_bytes::<LumpedInterferenceNodeResponse<'a>>(count)?;
        let expected_prepared_bytes = vector_bytes::<PreparedNodePower<'a>>(count)?;
        admit_bytes(
            expected_prepared_bytes,
            result_bytes,
            limits.maximum_workspace_bytes,
        )?;
        let mut nodes = Vec::new();
        nodes.try_reserve_exact(count).map_err(|_| {
            lumped_error(
                "lumped interference",
                "binding allocation failed".to_string(),
            )
        })?;
        let prepared_bytes = vector_bytes::<PreparedNodePower<'a>>(nodes.capacity())?;
        admit_bytes(prepared_bytes, result_bytes, limits.maximum_workspace_bytes)?;
        // The original network constructor guarantees a nonempty node set.
        let source_clock = bindings[0].source_clock;
        let mut spectral_terms = 0usize;
        for (index, binding) in bindings.iter().copied().enumerate() {
            poll(&mut cancelled, "lumped interference preparation", index)?;
            if binding.node_name.trim().is_empty()
                || binding.input_unit.trim().is_empty()
                || binding.source_clock.trim().is_empty()
            {
                return Err(lumped_error(
                    binding.node_name,
                    "node name, spectrum unit and source clock must be explicit".to_string(),
                ));
            }
            if binding.source_clock != source_clock {
                return Err(lumped_error(
                    binding.node_name,
                    "power bindings declare different source clocks".to_string(),
                ));
            }
            if !binding.watts_per_spectrum_unit.is_finite()
                || binding.watts_per_spectrum_unit <= 0.0
            {
                return Err(lumped_error(
                    binding.node_name,
                    "watts-per-spectrum-unit scale must be finite and positive".to_string(),
                ));
            }
            let node_index = network
                .nodes()
                .binary_search_by(|node| node.name().cmp(binding.node_name))
                .map_err(|_| {
                    lumped_error(
                        binding.node_name,
                        "power binding does not name an original network node".to_string(),
                    )
                })?;
            let node = &network.nodes()[node_index];
            let decay_rate_per_s = node.conductance_w_per_k() / node.capacitance_j_per_k();
            if !decay_rate_per_s.is_finite() || decay_rate_per_s <= 0.0 {
                return Err(lumped_error(
                    node.name(),
                    "original G/C decay rate must be representable and positive".to_string(),
                ));
            }
            spectral_terms = spectral_terms
                .checked_add(binding.spectrum.tones().len())
                .and_then(|sum| sum.checked_add(1))
                .ok_or_else(|| {
                    lumped_error("lumped interference", "spectral work overflow".to_string())
                })?;
            if spectral_terms
                > limits
                    .maximum_spectral_terms
                    .min(MAX_LUMPED_INTERFERENCE_TERMS)
            {
                return Err(lumped_error(
                    "lumped interference",
                    "spectral work budget exceeded".to_string(),
                ));
            }
            nodes.push(PreparedNodePower {
                node,
                binding,
                decay_rate_per_s,
            });
        }
        // The fixed original network bound also bounds this synchronous sort.
        nodes.sort_unstable_by(|a, b| a.node.name().cmp(b.node.name()));
        for (index, (node, prepared)) in network.nodes().iter().zip(&nodes).enumerate() {
            poll(&mut cancelled, "lumped interference preparation", index)?;
            if node.name() != prepared.node.name() {
                return Err(lumped_error(
                    node.name(),
                    "power bindings duplicate a node and omit another original node".to_string(),
                ));
            }
        }
        poll(&mut cancelled, "lumped interference preparation", count)?;
        Ok(Self {
            network,
            nodes,
            source_clock,
            gate,
            verdict,
            spectral_terms,
            limits,
            prepared_bytes,
        })
    }

    /// The actual network that owns the capacities, conductances and ambient.
    pub const fn network(&self) -> &'a LumpedNetwork {
        self.network
    }

    pub const fn source_clock(&self) -> &'a str {
        self.source_clock
    }

    pub const fn spectral_terms(&self) -> usize {
        self.spectral_terms
    }

    pub const fn gate(&self) -> BiotGate {
        self.gate
    }

    pub fn verdict(&self) -> &ValidityVerdict {
        &self.verdict
    }

    /// Solve one physical interval in the declared absolute source clock.
    ///
    /// Initial temperatures follow network.nodes()' deterministic name
    /// order, irrespective of the supplied power-binding order. Duration may
    /// be zero. Input power is signed and can represent heating or extraction.
    ///
    /// Initial and endpoint absolute temperatures must be finite and positive.
    /// This endpoint query does NOT prove that the entire intervening path
    /// remains positive, meets a temperature limit, or stays within a material
    /// model's validity range. No new validity is inferred from a finite result.
    pub fn response_at(
        &self,
        start_s: f64,
        duration_s: f64,
        initial_k: &[f64],
    ) -> Result<LumpedInterferenceResponse<'a>, ConductionError> {
        self.response_at_cancellable(start_s, duration_s, initial_k, || false)
    }

    /// Refusal/cancellation drops all partial results and never mutates the
    /// original network, prepared spectra or any previously returned response.
    pub fn response_at_cancellable(
        &self,
        start_s: f64,
        duration_s: f64,
        initial_k: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LumpedInterferenceResponse<'a>, ConductionError> {
        self.response_impl(
            start_s,
            duration_s,
            initial_k.iter().copied().map(Dd::from_f64),
            cancelled,
        )
    }

    /// Accept DD initial temperatures so successive responses retain weak
    /// heating even when each individual change is below an f64 endpoint ulp.
    pub fn response_at_dd(
        &self,
        start_s: f64,
        duration_s: f64,
        initial_k: &[Dd],
    ) -> Result<LumpedInterferenceResponse<'a>, ConductionError> {
        self.response_at_dd_cancellable(start_s, duration_s, initial_k, || false)
    }

    pub fn response_at_dd_cancellable(
        &self,
        start_s: f64,
        duration_s: f64,
        initial_k: &[Dd],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LumpedInterferenceResponse<'a>, ConductionError> {
        self.response_impl(
            start_s,
            duration_s,
            initial_k
                .iter()
                .map(|value| Dd::from_f64(value.hi) + Dd::from_f64(value.lo)),
            cancelled,
        )
    }

    fn response_impl(
        &self,
        start_s: f64,
        duration_s: f64,
        initial_k: impl ExactSizeIterator<Item = Dd> + Clone,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LumpedInterferenceResponse<'a>, ConductionError> {
        poll(&mut cancelled, "lumped interference response", 0)?;
        if !start_s.is_finite() || !duration_s.is_finite() || duration_s < 0.0 {
            return Err(lumped_error(
                "lumped interference",
                "source-clock start must be finite and duration finite and non-negative"
                    .to_string(),
            ));
        }
        if initial_k.len() != self.nodes.len() {
            return Err(lumped_error(
                "lumped interference",
                "initial temperatures must cover every original node in name order".to_string(),
            ));
        }
        for (index, (prepared, temperature)) in self.nodes.iter().zip(initial_k.clone()).enumerate()
        {
            poll(&mut cancelled, "lumped interference response", index)?;
            let temperature = finite(prepared.node.name(), temperature, "initial temperature")?;
            if temperature.to_f64() <= 0.0 {
                return Err(lumped_error(
                    prepared.node.name(),
                    "initial absolute temperature must be finite and positive in K".to_string(),
                ));
            }
        }
        let mut nodes = Vec::new();
        nodes.try_reserve_exact(self.nodes.len()).map_err(|_| {
            lumped_error(
                "lumped interference",
                "response allocation failed".to_string(),
            )
        })?;
        admit_bytes(
            self.prepared_bytes,
            vector_bytes::<LumpedInterferenceNodeResponse<'a>>(nodes.capacity())?,
            self.limits.maximum_workspace_bytes,
        )?;
        let ambient = Dd::from_f64(self.network.ambient_k());
        for (index, (prepared, initial)) in self.nodes.iter().zip(initial_k).enumerate() {
            poll(&mut cancelled, "lumped interference response", index)?;
            let original = prepared.node;
            let name = original.name();
            let integrals = prepared
                .binding
                .spectrum
                .causal_integral_cancellable(
                    start_s,
                    duration_s,
                    prepared.decay_rate_per_s,
                    &mut cancelled,
                )
                .map_err(|error| match error {
                    ProductSpectrumError::Cancelled => ConductionError::Cancelled {
                        stage: "lumped interference spectrum",
                        at: index,
                    },
                    error => lumped_error(name, format!("thermal power spectrum refused: {error}")),
                })?;
            let capacity = Dd::from_f64(original.capacitance_j_per_k());
            let scale = Dd::from_f64(prepared.binding.watts_per_spectrum_unit);
            let initial_dd = initial;
            let initial_offset = finite(name, initial_dd - ambient, "initial offset")?;
            // Use 1-exp(-lambda*dt) supplied by the causal kernel, not a
            // subtraction of rounded temperatures or rounded exponentials.
            let relaxation_k = finite(
                name,
                initial_offset * integrals.one_minus_decay_dd(),
                "ambient relaxation",
            )?;
            let relaxation_j = finite(name, relaxation_k * capacity, "relaxation energy")?;
            let retained_j = finite(
                name,
                scale * integrals.weighted_integral_dd(),
                "retained source energy",
            )?;
            let input_j = finite(name, scale * integrals.integral_dd(), "input energy")?;
            let stored_j = finite(name, retained_j - relaxation_j, "storage change")?;
            let rejected_source = finite(
                name,
                integrals.integral_dd() - integrals.weighted_integral_dd(),
                "source exchange integral",
            )?;
            let ambient_j = finite(
                name,
                relaxation_j + finite(name, scale * rejected_source, "source exchange energy")?,
                "ambient exchange",
            )?;
            let delta_k = finite(name, stored_j / capacity, "temperature change")?;
            let temperature_k = finite(name, initial_dd + delta_k, "endpoint temperature")?;
            if temperature_k.to_f64() <= 0.0 {
                return Err(lumped_error(
                    name,
                    "endpoint absolute temperature is not positive in K".to_string(),
                ));
            }
            nodes.push(LumpedInterferenceNodeResponse {
                node: original,
                input_unit: prepared.binding.input_unit,
                watts_per_spectrum_unit: prepared.binding.watts_per_spectrum_unit,
                initial_temperature_k: initial,
                temperature_k,
                temperature_change_k: delta_k,
                input_energy_j: input_j,
                stored_energy_change_j: stored_j,
                ambient_exchange_j: ambient_j,
            });
        }
        poll(
            &mut cancelled,
            "lumped interference response",
            self.nodes.len(),
        )?;
        Ok(LumpedInterferenceResponse {
            network: self.network,
            source_clock: self.source_clock,
            start_s,
            duration_s,
            gate: self.gate,
            verdict: self.verdict.clone(),
            nodes,
        })
    }
}

/// An admitted endpoint and signed energy accounting for one original node.
#[derive(Debug)]
pub struct LumpedInterferenceNodeResponse<'a> {
    node: &'a LumpedNode,
    input_unit: &'a str,
    watts_per_spectrum_unit: f64,
    initial_temperature_k: Dd,
    temperature_k: Dd,
    temperature_change_k: Dd,
    input_energy_j: Dd,
    stored_energy_change_j: Dd,
    ambient_exchange_j: Dd,
}

impl<'a> LumpedInterferenceNodeResponse<'a> {
    pub const fn node(&self) -> &'a LumpedNode {
        self.node
    }

    pub const fn input_unit(&self) -> &'a str {
        self.input_unit
    }

    pub const fn watts_per_spectrum_unit(&self) -> f64 {
        self.watts_per_spectrum_unit
    }

    pub const fn initial_temperature_k_dd(&self) -> Dd {
        self.initial_temperature_k
    }

    pub const fn temperature_k_dd(&self) -> Dd {
        self.temperature_k
    }

    /// Retains sub-ulp changes even when the rounded endpoint equals its start.
    pub const fn temperature_change_k_dd(&self) -> Dd {
        self.temperature_change_k
    }

    /// Unweighted physical input integral, positive into the node.
    pub const fn input_energy_j_dd(&self) -> Dd {
        self.input_energy_j
    }

    /// C times the temperature change, including the initial ambient relaxation.
    pub const fn stored_energy_change_j_dd(&self) -> Dd {
        self.stored_energy_change_j
    }

    /// Integrated exchange through the original G path, positive OUTWARD.
    /// A colder-than-ambient node can have negative exchange.
    pub const fn ambient_exchange_j_dd(&self) -> Dd {
        self.ambient_exchange_j
    }

    /// Algebraic balance diagnostic, not an independent error certificate.
    pub fn energy_balance_residual_j_dd(&self) -> Dd {
        self.input_energy_j - self.stored_energy_change_j - self.ambient_exchange_j
    }
}

/// Complete atomic response, retaining original network/gate/clock provenance.
#[derive(Debug)]
pub struct LumpedInterferenceResponse<'a> {
    network: &'a LumpedNetwork,
    source_clock: &'a str,
    start_s: f64,
    duration_s: f64,
    gate: BiotGate,
    verdict: ValidityVerdict,
    nodes: Vec<LumpedInterferenceNodeResponse<'a>>,
}

impl<'a> LumpedInterferenceResponse<'a> {
    pub const fn network(&self) -> &'a LumpedNetwork {
        self.network
    }

    pub const fn source_clock(&self) -> &'a str {
        self.source_clock
    }

    pub const fn start_s(&self) -> f64 {
        self.start_s
    }

    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    pub const fn gate(&self) -> BiotGate {
        self.gate
    }

    pub fn verdict(&self) -> &ValidityVerdict {
        &self.verdict
    }

    /// Responses follow the original network's deterministic node-name order.
    pub fn nodes(&self) -> &[LumpedInterferenceNodeResponse<'a>] {
        &self.nodes
    }
}

fn finite(node: &str, value: Dd, quantity: &str) -> Result<Dd, ConductionError> {
    if value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite() {
        Ok(value)
    } else {
        Err(lumped_error(
            node,
            format!("interference-driven {quantity} is not representable"),
        ))
    }
}

fn poll(
    cancelled: &mut impl FnMut() -> bool,
    stage: &'static str,
    at: usize,
) -> Result<(), ConductionError> {
    if cancelled() {
        Err(ConductionError::Cancelled { stage, at })
    } else {
        Ok(())
    }
}

fn vector_bytes<T>(count: usize) -> Result<usize, ConductionError> {
    count
        .checked_mul(core::mem::size_of::<T>())
        .ok_or_else(|| lumped_error("lumped interference", "workspace size overflow".to_string()))
}

fn admit_bytes(prepared: usize, response: usize, maximum: usize) -> Result<(), ConductionError> {
    if prepared
        .checked_add(response)
        .is_none_or(|sum| sum > maximum)
    {
        Err(lumped_error(
            "lumped interference",
            "combined preparation/response workspace budget exceeded".to_string(),
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
