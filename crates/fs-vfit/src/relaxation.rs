//! Positive-real relaxation terms in series with the existing R-L-C load.
//!
//! Z(s) = Z_base(s) + sum R_j s/(s + w_j). Each term is a parallel
//! resistance/inertance pair, not an output equalizer. Its internal flow q
//! obeys q' = w (U-q), pressure is R (U-q), storage is R q^2/(2w),
//! and dissipation is R (U-q)^2. The same midpoint port work closes every
//! branch. Elimination reduces the complete step to the existing series-load
//! solver at an effective port impedance and a history-dependent incident wave.
//! Bilinear frequency warping remains explicit. Positive coefficients establish
//! passivity of this model, not agreement with any particular material or wall.

use crate::impedance::{ImpedanceFrame, SeriesImpedance, SeriesImpedanceSpec};
use crate::waveguide::WaveguideError;

/// Bounded inline branch storage, matching the existing fs-phs Foster producer.
pub const MAX_RELAXATION_BRANCHES: usize = 8;

/// One passive term R s/(s+w), with retained inertance R/w.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct RelaxationTerm {
    /// High-frequency resistance [Pa s/m^3], strictly positive; omit zero terms.
    pub resistance_pa_s_m3: f64,
    /// Positive relaxation rate [1/s] (the pole is at -rate).
    pub rate_per_s: f64,
}

/// Admitted, immutable passive load. Inactive inline slots are inaccessible.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelaxationImpedanceSpec {
    base: SeriesImpedanceSpec,
    terms: [RelaxationTerm; MAX_RELAXATION_BRANCHES],
    count: usize,
}

impl RelaxationImpedanceSpec {
    /// Keep the supplied branch order; do not fit, repair or merge coefficients.
    /// An empty branch list is exactly the existing series R-L-C load.
    ///
    /// # Errors
    /// Invalid base, nonpositive/nonfinite branch data, or too many branches.
    pub fn new(base: SeriesImpedanceSpec, terms: &[RelaxationTerm]) -> Result<Self, WaveguideError> {
        base.validate()?;
        if terms.len() > MAX_RELAXATION_BRANCHES {
            return Err(WaveguideError("relaxation load exceeds the eight-branch limit"));
        }
        for term in terms {
            if !term.resistance_pa_s_m3.is_finite() || term.resistance_pa_s_m3 <= 0.0
                || !term.rate_per_s.is_finite() || term.rate_per_s <= 0.0
            {
                return Err(WaveguideError("relaxation branches require positive finite resistance and rate"));
            }
        }
        let mut out = Self { base, terms: [RelaxationTerm::default(); MAX_RELAXATION_BRANCHES], count: terms.len() };
        out.terms[..terms.len()].copy_from_slice(terms);
        Ok(out)
    }

    /// The declared R-L-C contribution, excluding every relaxation branch.
    #[must_use]
    pub const fn base(&self) -> SeriesImpedanceSpec { self.base }

    /// The actual, ordered branch coefficients; no hidden fitted values.
    #[must_use]
    pub fn terms(&self) -> &[RelaxationTerm] { &self.terms[..self.count] }
}

/// Immutable candidate containing the complete port observation and all states.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct RelaxationFrame {
    /// Complete load pressure, flow, energy, dissipation and work. The state
    /// inside this record is the base R-L-C state only; branch states follow.
    pub port: ImpedanceFrame,
    /// End-of-step energy in the relaxation branches alone [J].
    pub relaxation_stored_energy_j: f64,
    flows: [f64; MAX_RELAXATION_BRANCHES],
    count: usize,
}

impl RelaxationFrame {
    /// End-of-step internal branch flows [m^3/s], in declared order.
    #[must_use]
    pub fn branch_flows_m3_s(&self) -> &[f64] { &self.flows[..self.count] }
}

/// Finite-state boundary loss; construction and stepping allocate no heap storage.
/// All branch histories share the same solved port flow, not lagged copies.
#[derive(Debug, Clone, Copy)]
pub struct RelaxationImpedance {
    spec: RelaxationImpedanceSpec,
    base: SeriesImpedance,
    port_impedance: f64,
    dt: f64,
    flows: [f64; MAX_RELAXATION_BRANCHES],
    retained: [f64; MAX_RELAXATION_BRANCHES],
    driven: [f64; MAX_RELAXATION_BRANCHES],
    resistance: [f64; MAX_RELAXATION_BRANCHES],
    sqrt_inertance: [f64; MAX_RELAXATION_BRANCHES],
}

impl RelaxationImpedance {
    /// Bind to a characteristic port at zero stored energy. The effective port
    /// resistance is derived from midpoint elimination; it is not physical loss.
    ///
    /// # Errors
    /// Invalid port/time, unrepresentable branch storage or discrete coefficients.
    pub fn new(spec: RelaxationImpedanceSpec, port_impedance: f64, dt: f64) -> Result<Self, WaveguideError> {
        // Validate the physical port before doing branch arithmetic.
        let mut base = SeriesImpedance::new(spec.base, port_impedance, dt)?;
        let mut retained = [0.0; MAX_RELAXATION_BRANCHES];
        let mut driven = [0.0; MAX_RELAXATION_BRANCHES];
        let mut resistance = [0.0; MAX_RELAXATION_BRANCHES];
        let mut sqrt_inertance = [0.0; MAX_RELAXATION_BRANCHES];
        let mut effective_port = port_impedance;
        for (i, term) in spec.terms().iter().enumerate() {
            let k = term.rate_per_s * (0.5 * dt);
            let inertia = term.resistance_pa_s_m3 / term.rate_per_s;
            retained[i] = 1.0 / (1.0 + k);
            driven[i] = k / (1.0 + k);
            resistance[i] = term.resistance_pa_s_m3 * retained[i];
            if ![k, inertia, retained[i], driven[i], resistance[i]].iter().all(|v| v.is_finite() && *v > 0.0) {
                return Err(WaveguideError("relaxation storage or midpoint coefficient is not representable"));
            }
            sqrt_inertance[i] = fs_math::det::sqrt(inertia);
            effective_port += resistance[i];
        }
        if spec.count != 0 {
            base = SeriesImpedance::new(spec.base, effective_port, dt)?;
        }
        Ok(Self { spec, base, port_impedance, dt, flows: [0.0; MAX_RELAXATION_BRANCHES],
            retained, driven, resistance, sqrt_inertance })
    }

    /// Immutable physical coefficient declaration.
    #[must_use]
    pub const fn spec(&self) -> &RelaxationImpedanceSpec { &self.spec }

    /// Accepted branch histories, not a fitted response or accumulated energy.
    #[must_use]
    pub fn branch_flows_m3_s(&self) -> &[f64] { &self.flows[..self.spec.count] }

    fn branch_energy(&self, flows: &[f64; MAX_RELAXATION_BRANCHES]) -> f64 {
        let mut sum = 0.0;
        for (i, flow) in flows.iter().enumerate().take(self.spec.count) {
            let scaled = self.sqrt_inertance[i] * flow;
            sum += 0.5 * scaled * scaled;
        }
        sum
    }

    /// Independently evaluated base plus branch energy [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.base.stored_energy_j() + self.branch_energy(&self.flows)
    }

    /// Preview a simultaneous base/branch/wave-port step without changing state.
    ///
    /// # Errors
    /// Nonfinite input, history drive, state, reflected wave, storage or work.
    pub fn preview_step(&self, incident: f64) -> Result<RelaxationFrame, WaveguideError> {
        if self.spec.count == 0 {
            return Ok(RelaxationFrame { port: self.base.preview_step(incident)?, ..RelaxationFrame::default() });
        }
        let mut effective_incident = incident;
        for i in 0..self.spec.count {
            effective_incident += (0.5 * self.resistance[i]) * self.flows[i];
        }
        let base = self.base.preview_step(effective_incident)?;
        let flow = base.flow_m3_s;
        let mut flows = self.flows;
        let mut dissipated = base.dissipated_energy_j;
        for (i, term) in self.spec.terms().iter().enumerate() {
            let difference = flow - self.flows[i];
            // Compute k/(1+k) directly, not 1-retained (which loses tiny k).
            flows[i] += (2.0 * self.driven[i]) * difference;
            let resistor_flow = self.retained[i] * difference;
            dissipated += (term.resistance_pa_s_m3 * resistor_flow) * (resistor_flow * self.dt);
            if !flows[i].is_finite() {
                return Err(WaveguideError("relaxation branch state left the finite set"));
            }
        }
        let relaxed_energy = self.branch_energy(&flows);
        let reflected = incident - self.port_impedance * flow;
        let pressure = incident + reflected;
        let stored = base.stored_energy_j + relaxed_energy;
        let port = ImpedanceFrame {
            state: base.state, reflected_pressure_pa: reflected, pressure_pa: pressure,
            flow_m3_s: flow, stored_energy_j: stored,
            storage_change_j: stored - self.stored_energy_j(),
            dissipated_energy_j: dissipated, supplied_work_j: pressure * (flow * self.dt),
        };
        if ![incident, effective_incident, reflected, pressure, relaxed_energy, stored,
            port.storage_change_j, dissipated, port.supplied_work_j, port.balance_residual_j()]
            .iter().all(|v| v.is_finite())
        {
            return Err(WaveguideError("relaxation load observation left the finite set"));
        }
        Ok(RelaxationFrame { port, relaxation_stored_energy_j: relaxed_energy, flows, count: self.spec.count })
    }

    /// Only the owner can commit a preview; callers cannot inject forged history.
    pub(crate) fn accept_frame(&mut self, frame: RelaxationFrame) {
        self.base.accept_frame(frame.port);
        self.flows = frame.flows;
    }

    /// Advance only after a complete immutable preview succeeds.
    ///
    /// # Errors
    /// Same refusal boundary as `preview_step`; failed calls preserve all history.
    pub fn step(&mut self, incident: f64) -> Result<RelaxationFrame, WaveguideError> {
        let frame = self.preview_step(incident)?;
        self.accept_frame(frame);
        Ok(frame)
    }
}
