//! The lumped-network reduced transient rung: the cheap tier of the fidelity
//! graph, with an explicit validity gate.
//!
//! Bead `frankensim-extreal-program-f85xj.5.13`, the third item staged by
//! `f85xj.5.9`. The method-of-lines march in [`crate::transient`] is the
//! expensive answer; most early design questions ("does this die exceed its
//! limit during a thirty-second burst?") are answered by a small RC network in
//! microseconds.
//!
//! # Cost and authority are separate axes
//!
//! The point of a fidelity graph is that a cheaper model is not merely a
//! faster version of an expensive one — it is a WEAKER one, and both facts
//! have to be explicit. This rung is cheaper AND it is only valid in a stated
//! regime. It therefore refuses outside that regime rather than returning a
//! cheap wrong number, because a cheap wrong number is worse than no number:
//! it looks like an answer.
//!
//! # The validity gate is the Biot number
//!
//! Lumping asserts that a body is isothermal — that internal conduction is
//! fast compared with surface transfer. The Biot number `Bi = h Lc / k`
//! measures exactly that, and the `fs-vvreg` Level-A lumped row states the
//! admitted context as `Bi <= 0.1`. [`LumpedNode`] therefore carries the
//! quantities that determine Biot and [`BiotGate`] adjudicates them, so a
//! caller cannot construct a network that silently lies about its own
//! applicability.
//!
//! # One power declaration, two rungs
//!
//! Node power comes from the same [`crate::power::PowerMap`] the full rung
//! consumes — the caller supplies the per-node watts that map already
//! validated and audited. Two rungs describing the same hardware differently
//! is the failure mode a fidelity graph exists to prevent, so the reduced
//! rung deliberately has no power vocabulary of its own.

use crate::ConductionError;

/// The admitted Biot ceiling for lumped treatment.
///
/// Matches the `fs-vvreg` Level-A `thermal-a-lumped-transient` row's declared
/// context. Duplicated as a constant rather than read at runtime because this
/// crate must not depend on the corpus registry; `tests/lumped.rs` asserts the
/// two agree, so the duplication cannot drift silently.
pub const LUMPED_BIOT_CEILING: f64 = 0.1;

/// Maximum nodes admitted in one reduced network.
pub const MAX_LUMPED_NODES: usize = 1_024;

/// One lumped node: an isothermal body with a surface path to ambient.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedNode {
    name: String,
    capacitance_j_per_k: f64,
    conductance_w_per_k: f64,
    characteristic_length_m: f64,
    conductivity_w_per_m_k: f64,
    surface_area_m2: f64,
}

impl LumpedNode {
    /// Declare a node from the quantities that determine both its dynamics
    /// and its validity.
    ///
    /// `conductance_w_per_k` is the surface path to ambient, `h·A`. The
    /// characteristic length, conductivity and area are not redundant with
    /// it: they are what make the Biot number computable, and a node that
    /// could not state its own Biot could not be gated.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for a blank name or any non-positive
    /// or non-finite quantity.
    pub fn new(
        name: impl Into<String>,
        capacitance_j_per_k: f64,
        conductance_w_per_k: f64,
        characteristic_length_m: f64,
        conductivity_w_per_m_k: f64,
        surface_area_m2: f64,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(lumped_error(
                "<unnamed>",
                "lumped node has a blank name".to_string(),
            ));
        }
        for (value, field) in [
            (capacitance_j_per_k, "capacitance"),
            (conductance_w_per_k, "surface conductance"),
            (characteristic_length_m, "characteristic length"),
            (conductivity_w_per_m_k, "conductivity"),
            (surface_area_m2, "surface area"),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(lumped_error(
                    &name,
                    format!("{field} {value} is not finite and positive"),
                ));
            }
        }
        Ok(Self {
            name,
            capacitance_j_per_k,
            conductance_w_per_k,
            characteristic_length_m,
            conductivity_w_per_m_k,
            surface_area_m2,
        })
    }

    /// Node name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Thermal capacitance, J/K.
    #[must_use]
    pub const fn capacitance_j_per_k(&self) -> f64 {
        self.capacitance_j_per_k
    }

    /// Surface conductance to ambient, W/K.
    #[must_use]
    pub const fn conductance_w_per_k(&self) -> f64 {
        self.conductance_w_per_k
    }

    /// Effective surface transfer coefficient `h = (hA) / A`, W/(m²·K).
    #[must_use]
    pub fn transfer_coefficient_w_per_m2_k(&self) -> f64 {
        self.conductance_w_per_k / self.surface_area_m2
    }

    /// The node's Biot number, `h Lc / k`.
    ///
    /// This is the quantity that decides whether lumping is admissible at
    /// all, so it is derived from declared inputs rather than accepted as one.
    #[must_use]
    pub fn biot(&self) -> f64 {
        self.transfer_coefficient_w_per_m2_k() * self.characteristic_length_m
            / self.conductivity_w_per_m_k
    }

    /// The node's time constant `C / (hA)`, s.
    #[must_use]
    pub fn time_constant_s(&self) -> f64 {
        self.capacitance_j_per_k / self.conductance_w_per_k
    }
}

/// The verdict of the lumped validity gate.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidityVerdict {
    /// Every node is inside the declared Biot ceiling.
    Admitted {
        /// Largest Biot number in the network.
        worst_biot: f64,
    },
    /// At least one node is outside it. The rung REFUSES rather than
    /// returning a cheap wrong number.
    Refused {
        /// The offending node.
        node: String,
        /// Its Biot number.
        biot: f64,
        /// The ceiling it exceeded.
        ceiling: f64,
    },
}

impl ValidityVerdict {
    /// Whether the network may be solved on this rung.
    #[must_use]
    pub const fn admitted(&self) -> bool {
        matches!(self, Self::Admitted { .. })
    }
}

/// The Biot gate, with its ceiling declared rather than assumed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BiotGate {
    ceiling: f64,
}

impl BiotGate {
    /// The gate at the corpus-declared ceiling.
    #[must_use]
    pub const fn corpus_default() -> Self {
        Self {
            ceiling: LUMPED_BIOT_CEILING,
        }
    }

    /// A gate at a caller-declared ceiling.
    ///
    /// Loosening this is a modelling decision the caller owns; the returned
    /// verdict always reports the ceiling actually applied, so a loosened
    /// gate cannot be mistaken for the corpus one.
    ///
    /// # Errors
    /// A non-positive or non-finite ceiling.
    pub fn at(ceiling: f64) -> Result<Self, ConductionError> {
        if !ceiling.is_finite() || ceiling <= 0.0 {
            return Err(lumped_error(
                "biot gate",
                format!("ceiling {ceiling} is not finite and positive"),
            ));
        }
        Ok(Self { ceiling })
    }

    /// The applied ceiling.
    #[must_use]
    pub const fn ceiling(self) -> f64 {
        self.ceiling
    }

    /// Adjudicate a network.
    #[must_use]
    pub fn adjudicate(self, network: &LumpedNetwork) -> ValidityVerdict {
        let mut worst = 0.0f64;
        for node in &network.nodes {
            let biot = node.biot();
            if biot > self.ceiling {
                return ValidityVerdict::Refused {
                    node: node.name.clone(),
                    biot,
                    ceiling: self.ceiling,
                };
            }
            worst = worst.max(biot);
        }
        ValidityVerdict::Admitted { worst_biot: worst }
    }
}

/// A reduced thermal network: isothermal nodes coupled to one ambient.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedNetwork {
    nodes: Vec<LumpedNode>,
    ambient_k: f64,
}

impl LumpedNetwork {
    /// Admit a network.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for an empty or oversized node set, a
    /// duplicated node name, or a non-finite ambient.
    pub fn new(nodes: Vec<LumpedNode>, ambient_k: f64) -> Result<Self, ConductionError> {
        if nodes.is_empty() {
            return Err(lumped_error(
                "lumped network",
                "declares no nodes".to_string(),
            ));
        }
        if nodes.len() > MAX_LUMPED_NODES {
            return Err(lumped_error(
                "lumped network",
                format!(
                    "declares {} nodes, above the admitted maximum {MAX_LUMPED_NODES}",
                    nodes.len()
                ),
            ));
        }
        if !ambient_k.is_finite() {
            return Err(lumped_error(
                "lumped network",
                format!("ambient temperature {ambient_k} K is not finite"),
            ));
        }
        let mut nodes = nodes;
        nodes.sort_by(|a, b| a.name.cmp(&b.name));
        for pair in nodes.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(lumped_error(
                    &pair[0].name,
                    "lumped node name is duplicated".to_string(),
                ));
            }
        }
        Ok(Self { nodes, ambient_k })
    }

    /// Nodes in deterministic name order.
    #[must_use]
    pub fn nodes(&self) -> &[LumpedNode] {
        &self.nodes
    }

    /// Declared ambient temperature, K.
    #[must_use]
    pub const fn ambient_k(&self) -> f64 {
        self.ambient_k
    }

    /// Steady temperature of every node under a constant power vector, K.
    ///
    /// Nodes are uncoupled here (each has its own path to one ambient), so
    /// the steady state is `T_amb + P/(hA)` node by node.
    ///
    /// # Errors
    /// A power vector of the wrong length, or a non-finite entry.
    pub fn steady(&self, power_w: &[f64]) -> Result<Vec<f64>, ConductionError> {
        self.admit_power(power_w)?;
        Ok(self
            .nodes
            .iter()
            .zip(power_w.iter())
            .map(|(node, power)| self.ambient_k + power / node.conductance_w_per_k)
            .collect())
    }

    /// The closed-form response at time `t` from an initial temperature.
    ///
    /// Each node is a first-order system, so the reduced rung has an ANALYTIC
    /// solution — no time stepping and no step-size choice. That is exactly
    /// what makes it the cheap tier, and it also means the rung introduces no
    /// integration error of its own to confuse with model error.
    ///
    /// # Errors
    /// A mismatched power or initial vector, a non-finite entry, or a
    /// negative time.
    pub fn response_at(
        &self,
        power_w: &[f64],
        initial_k: &[f64],
        time_s: f64,
    ) -> Result<Vec<f64>, ConductionError> {
        self.admit_power(power_w)?;
        if initial_k.len() != self.nodes.len() {
            return Err(lumped_error(
                "lumped network",
                format!(
                    "initial vector has {} entries for {} nodes",
                    initial_k.len(),
                    self.nodes.len()
                ),
            ));
        }
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(lumped_error(
                "lumped network",
                format!("time {time_s} s is not finite and non-negative"),
            ));
        }
        let steady = self.steady(power_w)?;
        let mut out = Vec::with_capacity(self.nodes.len());
        for ((node, start), settle) in self.nodes.iter().zip(initial_k.iter()).zip(steady.iter()) {
            let decay = (-time_s / node.time_constant_s()).exp();
            let value = settle + (start - settle) * decay;
            if !value.is_finite() {
                return Err(lumped_error(
                    &node.name,
                    format!("response at {time_s} s left the finite range"),
                ));
            }
            out.push(value);
        }
        Ok(out)
    }

    fn admit_power(&self, power_w: &[f64]) -> Result<(), ConductionError> {
        if power_w.len() != self.nodes.len() {
            return Err(lumped_error(
                "lumped network",
                format!(
                    "power vector has {} entries for {} nodes",
                    power_w.len(),
                    self.nodes.len()
                ),
            ));
        }
        for (node, power) in self.nodes.iter().zip(power_w.iter()) {
            if !power.is_finite() || *power < 0.0 {
                return Err(lumped_error(
                    &node.name,
                    format!("node power {power} W is not finite and non-negative"),
                ));
            }
        }
        Ok(())
    }
}

/// A solved reduced-rung response, with the verdict that admitted it.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedSolution {
    /// Nodal temperatures at the requested time, K.
    pub temperature_k: Vec<f64>,
    /// Steady temperatures the response is approaching, K.
    pub steady_k: Vec<f64>,
    /// The validity verdict under which this was produced.
    pub verdict: ValidityVerdict,
}

/// Solve the reduced rung, GATED.
///
/// The gate runs first and a refusal short-circuits: outside the declared
/// Biot regime this returns an error rather than a number, because the whole
/// value of a cheap rung is destroyed if it answers questions it cannot
/// answer.
///
/// # Errors
/// [`ConductionError::ScenarioRow`] when the gate refuses or any input is
/// inadmissible.
pub fn solve_gated(
    network: &LumpedNetwork,
    gate: BiotGate,
    power_w: &[f64],
    initial_k: &[f64],
    time_s: f64,
) -> Result<LumpedSolution, ConductionError> {
    let verdict = gate.adjudicate(network);
    match &verdict {
        ValidityVerdict::Refused {
            node,
            biot,
            ceiling,
        } => Err(lumped_error(
            node,
            format!(
                "Biot {biot} exceeds the admitted ceiling {ceiling}: the body is not isothermal, so the lumped rung cannot answer this and escalating to the full transient is the correct move rather than loosening the gate"
            ),
        )),
        ValidityVerdict::Admitted { .. } => Ok(LumpedSolution {
            temperature_k: network.response_at(power_w, initial_k, time_s)?,
            steady_k: network.steady(power_w)?,
            verdict,
        }),
    }
}

/// Build a single-node reduced model from a measured steady response.
///
/// Given the steady temperature rise a body reaches under a known power, the
/// surface conductance is `hA = P / dT`. That is an EXTRACTION, not a
/// derivation: it assumes the steady rise is dominated by the surface path
/// the reduced model represents, and it inherits every modelling assumption
/// of the run it was extracted from. Callers should record it as model-form
/// evidence, not as a measured property.
///
/// # Errors
/// [`ConductionError::ScenarioRow`] for a non-positive rise or power, or any
/// non-finite input.
#[allow(clippy::too_many_arguments)]
pub fn extract_node_from_steady_rise(
    name: impl Into<String>,
    power_w: f64,
    steady_rise_k: f64,
    capacitance_j_per_k: f64,
    characteristic_length_m: f64,
    conductivity_w_per_m_k: f64,
    surface_area_m2: f64,
) -> Result<LumpedNode, ConductionError> {
    let name = name.into();
    for (value, field) in [
        (power_w, "extraction power"),
        (steady_rise_k, "steady rise"),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(lumped_error(
                &name,
                format!("{field} {value} is not finite and positive"),
            ));
        }
    }
    let conductance = power_w / steady_rise_k;
    LumpedNode::new(
        name,
        capacitance_j_per_k,
        conductance,
        characteristic_length_m,
        conductivity_w_per_m_k,
        surface_area_m2,
    )
}

fn lumped_error(node: &str, what: String) -> ConductionError {
    ConductionError::ScenarioRow {
        region: node.to_string(),
        what,
        fix: "correct the declared lumped network, or escalate to the full transient rung"
            .to_string(),
    }
}
