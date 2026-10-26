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
//!
//! Complementary parallel resistance/compliance arms add R*w/(s+w). Their
//! pressure y obeys y'=w*(R*U-y), storage is y^2/(2*R*w), and loss is y^2/R.
//! They represent stress relaxation (for example a Maxwell material wall), not
//! the high-pass inertive branches above. Both families share ONE midpoint port
//! solve and publish together; no branch is stepped against a lagged flow.
//!
//! [`ParallelCellLoad`] is a separate, flow-driven series chain of parallel
//! R, L and C cells, Z = sum 1/(1/R + 1/(sL) + sC), for a caller that has
//! already solved the step-mean flow with the chain's affine pressure in its
//! own equation. A cell with C = 0 is exactly the high-pass term above. Every
//! R-L network is a Foster sum of such terms, so its Re Z / w^2 can only fall
//! with frequency. A load whose Re Z / w^2 rises over some band (an opening's
//! radiation, as diffraction turns the free monopole into a baffled one) needs
//! cells with C > 0. Each cell's storage is (L i^2 + C v^2)/2 and its loss is
//! v^2/R. The midpoint rule closes every cell's work exactly. Its own
//! sixteen-cell budget leaves the eight-branch load above unchanged in size.

use crate::impedance::{ImpedanceFrame, SeriesImpedance, SeriesImpedanceSpec};
use crate::waveguide::WaveguideError;
use fs_math::c64::C64;

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

/// One passive parallel-RC term R*w/(s+w), with compliance 1/(R*w).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ComplianceRelaxationTerm {
    /// Static resistance [Pa s/m^3], strictly positive; omit zero arms.
    pub resistance_pa_s_m3: f64,
    /// Positive relaxation rate [1/s].
    pub rate_per_s: f64,
}

/// Admitted, immutable passive load. Inactive inline slots are inaccessible.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelaxationImpedanceSpec {
    base: SeriesImpedanceSpec,
    terms: [RelaxationTerm; MAX_RELAXATION_BRANCHES],
    count: usize,
    compliance_terms: [ComplianceRelaxationTerm; MAX_RELAXATION_BRANCHES],
    compliance_count: usize,
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
        let mut out = Self { base, terms: [RelaxationTerm::default(); MAX_RELAXATION_BRANCHES], count: terms.len(),
            compliance_terms: [ComplianceRelaxationTerm::default(); MAX_RELAXATION_BRANCHES], compliance_count: 0 };
        out.terms[..terms.len()].copy_from_slice(terms);
        Ok(out)
    }

    /// Add the complementary low-pass compliance arms without changing RL arms.
    /// The combined family has at most eight branches. Repeated attachment
    /// refuses rather than discarding any previously declared material spectrum.
    ///
    /// # Errors
    /// Already attached arms, exceeded combined budget, or nonpositive data.
    pub fn with_compliance_terms(mut self, terms: &[ComplianceRelaxationTerm]) -> Result<Self, WaveguideError> {
        if self.compliance_count != 0 || terms.len() > MAX_RELAXATION_BRANCHES - self.count {
            return Err(WaveguideError("compliance relaxation requires one attachment within the combined eight-branch limit"));
        }
        for term in terms {
            if !term.resistance_pa_s_m3.is_finite() || term.resistance_pa_s_m3 <= 0.0
                || !term.rate_per_s.is_finite() || term.rate_per_s <= 0.0 {
                return Err(WaveguideError("compliance relaxation requires positive finite resistance and rate"));
            }
        }
        self.compliance_terms[..terms.len()].copy_from_slice(terms);
        self.compliance_count = terms.len();
        Ok(self)
    }

    /// Actual ordered R*w/(s+w) arms, distinct from the R*s/(s+w) terms.
    #[must_use]
    pub fn compliance_terms(&self) -> &[ComplianceRelaxationTerm] {
        &self.compliance_terms[..self.compliance_count]
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
    pressures: [f64; MAX_RELAXATION_BRANCHES],
    pressure_count: usize,
}

impl RelaxationFrame {
    /// End-of-step internal branch flows [m^3/s], in declared order.
    #[must_use]
    pub fn branch_flows_m3_s(&self) -> &[f64] { &self.flows[..self.count] }
    /// End-of-step internal RC pressures [Pa], in compliance-term order.
    #[must_use]
    pub fn branch_pressures_pa(&self) -> &[f64] { &self.pressures[..self.pressure_count] }
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
    pressures: [f64; MAX_RELAXATION_BRANCHES],
    pressure_retained: [f64; MAX_RELAXATION_BRANCHES],
    pressure_driven: [f64; MAX_RELAXATION_BRANCHES],
    sqrt_compliance: [f64; MAX_RELAXATION_BRANCHES],
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
        let mut pressure_retained = [0.0; MAX_RELAXATION_BRANCHES];
        let mut pressure_driven = [0.0; MAX_RELAXATION_BRANCHES];
        let mut sqrt_compliance = [0.0; MAX_RELAXATION_BRANCHES];
        for (i, term) in spec.compliance_terms().iter().enumerate() {
            let k = term.rate_per_s * (0.5 * dt);
            // Dividing sequentially avoids overflowing R*w merely to invert it.
            let compliance = (1.0 / term.resistance_pa_s_m3) / term.rate_per_s;
            pressure_retained[i] = 1.0 / (1.0 + k);
            pressure_driven[i] = k / (1.0 + k);
            let r = term.resistance_pa_s_m3 * pressure_driven[i];
            if ![k, compliance, pressure_retained[i], pressure_driven[i], r]
                .iter().all(|v| v.is_finite() && *v > 0.0) {
                return Err(WaveguideError("compliance relaxation storage or midpoint coefficient is not representable"));
            }
            sqrt_compliance[i] = fs_math::det::sqrt(compliance);
            effective_port += r;
        }
        if spec.count != 0 || spec.compliance_count != 0 {
            base = SeriesImpedance::new(spec.base, effective_port, dt)?;
        }
        Ok(Self { spec, base, port_impedance, dt, flows: [0.0; MAX_RELAXATION_BRANCHES],
            retained, driven, resistance, sqrt_inertance,
            pressures: [0.0; MAX_RELAXATION_BRANCHES], pressure_retained, pressure_driven, sqrt_compliance })
    }

    /// Immutable physical coefficient declaration.
    #[must_use]
    pub const fn spec(&self) -> &RelaxationImpedanceSpec { &self.spec }

    /// Accepted coordinates of the base R-L-C contribution, excluding branch flows.
    #[must_use]
    pub const fn base_state(&self) -> crate::impedance::ImpedanceState { self.base.state() }

    /// Accepted branch histories, not a fitted response or accumulated energy.
    #[must_use]
    pub fn branch_flows_m3_s(&self) -> &[f64] { &self.flows[..self.spec.count] }

    /// Accepted RC pressure histories [Pa], not the RL internal flows.
    #[must_use]
    pub fn branch_pressures_pa(&self) -> &[f64] { &self.pressures[..self.spec.compliance_count] }

    fn compliance_energy(&self, pressures: &[f64; MAX_RELAXATION_BRANCHES]) -> f64 {
        let mut sum = 0.0;
        for (i, pressure) in pressures.iter().enumerate().take(self.spec.compliance_count) {
            let scaled = self.sqrt_compliance[i] * pressure;
            sum += 0.5 * scaled * scaled;
        }
        sum
    }

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
        self.base.stored_energy_j() + self.branch_energy(&self.flows) + self.compliance_energy(&self.pressures)
    }

    /// Preview a simultaneous base/branch/wave-port step without changing state.
    ///
    /// # Errors
    /// Nonfinite input, history drive, state, reflected wave, storage or work.
    pub fn preview_step(&self, incident: f64) -> Result<RelaxationFrame, WaveguideError> {
        self.preview_with_loss(incident, None)
    }

    pub(crate) fn preview_with_loss(&self, incident: f64,
        loss: Option<crate::quadratic_loss::QuadraticLossSpec>)
        -> Result<RelaxationFrame, WaveguideError> {
        if self.spec.count == 0 && self.spec.compliance_count == 0 {
            return Ok(RelaxationFrame { port: self.base.preview_with_loss(incident, loss)?, ..RelaxationFrame::default() });
        }
        let mut effective_incident = incident;
        for i in 0..self.spec.count {
            effective_incident += (0.5 * self.resistance[i]) * self.flows[i];
        }
        for i in 0..self.spec.compliance_count {
            effective_incident -= (0.5 * self.pressure_retained[i]) * self.pressures[i];
        }
        let base = self.base.preview_with_loss(effective_incident, loss)?;
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
        let mut pressures = self.pressures;
        for (i, term) in self.spec.compliance_terms().iter().enumerate() {
            let drive = term.resistance_pa_s_m3 * flow;
            pressures[i] += (2.0 * self.pressure_driven[i]) * (drive - self.pressures[i]);
            // Convex form avoids subtracting two almost equal endpoint pressures.
            let midpoint = self.pressure_retained[i] * self.pressures[i] + self.pressure_driven[i] * drive;
            dissipated += (midpoint / term.resistance_pa_s_m3) * (midpoint * self.dt);
            if !pressures[i].is_finite() || !midpoint.is_finite() {
                return Err(WaveguideError("compliance relaxation branch state left the finite set"));
            }
        }
        let relaxed_energy = self.branch_energy(&flows) + self.compliance_energy(&pressures);
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
        Ok(RelaxationFrame { port, relaxation_stored_energy_j: relaxed_energy, flows, count: self.spec.count, pressures, pressure_count: self.spec.compliance_count })
    }

    /// Only the owner can commit a preview; callers cannot inject forged history.
    pub(crate) fn accept_frame(&mut self, frame: RelaxationFrame) {
        self.base.accept_frame(frame.port);
        self.flows = frame.flows;
        self.pressures = frame.pressures;
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

/// Bounded inline cell storage of a [`ParallelCellLoad`].
pub const MAX_PARALLEL_CELLS: usize = 16;

/// One parallel R, L and C carrying the series flow: 1/(1/R + 1/(sL) + sC).
/// C = 0 is exactly the high-pass term R s/(s + R/L).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ParallelCell {
    /// Parallel resistance [Pa s/m^3], positive.
    pub resistance_pa_s_m3: f64,
    /// Parallel inertance [Pa s^2/m^3], positive: the cell's static mass.
    pub inertance_pa_s2_m3: f64,
    /// Parallel compliance [m^3/Pa], zero or positive.
    pub compliance_m3_pa: f64,
}

impl ParallelCell {
    /// Continuous impedance at a nonzero Laplace variable `s`.
    #[must_use]
    pub fn impedance(&self, s: C64) -> C64 {
        (C64::from_re(1.0 / self.resistance_pa_s_m3) + s.scale(self.inertance_pa_s2_m3).recip()
            + s.scale(self.compliance_m3_pa)).recip()
    }
}

/// Admitted series chain of parallel cells, in declared order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParallelCellSpec {
    cells: [ParallelCell; MAX_PARALLEL_CELLS],
    count: usize,
}

impl ParallelCellSpec {
    /// Keep the supplied cells and their order; nothing is fitted or merged.
    ///
    /// # Errors
    /// No cells, more than sixteen, or a nonpositive R or L, or a negative C.
    pub fn new(cells: &[ParallelCell]) -> Result<Self, WaveguideError> {
        if cells.is_empty() || cells.len() > MAX_PARALLEL_CELLS {
            return Err(WaveguideError("a parallel-cell load requires one to sixteen cells"));
        }
        for cell in cells {
            if ![cell.resistance_pa_s_m3, cell.inertance_pa_s2_m3].iter().all(|v| v.is_finite() && *v > 0.0)
                || !cell.compliance_m3_pa.is_finite() || cell.compliance_m3_pa < 0.0 {
                return Err(WaveguideError("parallel cells require positive finite R and L and a finite C >= 0"));
            }
        }
        let mut out = Self { cells: [ParallelCell::default(); MAX_PARALLEL_CELLS], count: cells.len() };
        out.cells[..cells.len()].copy_from_slice(cells);
        Ok(out)
    }

    /// The actual, ordered cells.
    #[must_use]
    pub fn cells(&self) -> &[ParallelCell] { &self.cells[..self.count] }

    /// The chain's static mass, the sum of its inertances [Pa s^2/m^3].
    #[must_use]
    pub fn static_inertance_pa_s2_m3(&self) -> f64 {
        self.cells().iter().map(|cell| cell.inertance_pa_s2_m3).sum()
    }

    /// Continuous series impedance at a nonzero Laplace variable `s`.
    #[must_use]
    pub fn impedance(&self, s: C64) -> C64 {
        self.cells().iter().fold(C64::ZERO, |sum, cell| sum + cell.impedance(s))
    }
}

/// Immutable candidate of one flow-driven step.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ParallelCellFrame {
    /// Step-mean pressure across the chain [Pa].
    pub pressure_pa: f64,
    /// Step-mean flow through the chain [m^3/s].
    pub flow_m3_s: f64,
    /// Storage after the step [J].
    pub stored_energy_j: f64,
    /// Change of storage [J].
    pub storage_change_j: f64,
    /// Resistor loss over the step [J].
    pub dissipated_energy_j: f64,
    /// Pressure times flow times dt [J]; negative on release.
    pub supplied_work_j: f64,
    flows: [f64; MAX_PARALLEL_CELLS],
    pressures: [f64; MAX_PARALLEL_CELLS],
    count: usize,
}

impl ParallelCellFrame {
    /// Uncorrected local energy-balance residual [J].
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.dissipated_energy_j - self.supplied_work_j
    }
    /// End-of-step inductor flows [m^3/s], in cell order.
    #[must_use]
    pub fn inductor_flows_m3_s(&self) -> &[f64] { &self.flows[..self.count] }
    /// End-of-step cell pressures [Pa]: the capacitor's, or for C = 0 the
    /// step mean, in cell order.
    #[must_use]
    pub fn pressures_pa(&self) -> &[f64] { &self.pressures[..self.count] }
}

/// Midpoint coefficients of one cell. With G = 1/R + dt/(2L) + 2C/dt, the
/// step-mean pressure is v_mid = (U + (2C/dt) v - i) / G.
#[derive(Debug, Clone, Copy, Default)]
struct CellCoefficients {
    direct: f64,
    charge: f64,
    advance: f64,
    sqrt_inertance: f64,
    sqrt_compliance: f64,
}

/// A flow-driven series chain of parallel cells at zero initial storage. Per
/// cell, v_mid = (U + (2C/dt) v - i)/G, i' = i + (dt/L) v_mid and
/// v' = 2 v_mid - v. The work U v_mid dt equals the change of
/// (L i^2 + C v^2)/2 plus v_mid^2 dt / R exactly. No heap storage.
#[derive(Debug, Clone, Copy)]
pub struct ParallelCellLoad {
    spec: ParallelCellSpec,
    dt: f64,
    direct: f64,
    flows: [f64; MAX_PARALLEL_CELLS],
    pressures: [f64; MAX_PARALLEL_CELLS],
    coefficients: [CellCoefficients; MAX_PARALLEL_CELLS],
}

impl ParallelCellLoad {
    /// Bind the chain to a time step.
    ///
    /// # Errors
    /// Invalid step or an unrepresentable storage or midpoint coefficient.
    pub fn new(spec: ParallelCellSpec, dt: f64) -> Result<Self, WaveguideError> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(WaveguideError("a parallel-cell load requires a positive finite time step"));
        }
        let mut coefficients = [CellCoefficients::default(); MAX_PARALLEL_CELLS];
        let mut direct = 0.0;
        for (slot, cell) in coefficients.iter_mut().zip(spec.cells()) {
            let charge = 2.0 * cell.compliance_m3_pa / dt;
            let advance = dt / cell.inertance_pa_s2_m3;
            let conductance = 1.0 / cell.resistance_pa_s_m3 + 0.5 * advance + charge;
            *slot = CellCoefficients { direct: 1.0 / conductance, charge, advance,
                sqrt_inertance: fs_math::det::sqrt(cell.inertance_pa_s2_m3),
                sqrt_compliance: fs_math::det::sqrt(cell.compliance_m3_pa) };
            if ![advance, conductance, slot.direct, slot.sqrt_inertance].iter().all(|v| v.is_finite() && *v > 0.0)
                || !charge.is_finite() || !slot.sqrt_compliance.is_finite() {
                return Err(WaveguideError("parallel-cell storage or midpoint coefficient is not representable"));
            }
            direct += slot.direct;
        }
        Ok(Self { spec, dt, direct, flows: [0.0; MAX_PARALLEL_CELLS], pressures: [0.0; MAX_PARALLEL_CELLS],
            coefficients })
    }

    /// Immutable cell declaration.
    #[must_use]
    pub const fn spec(&self) -> &ParallelCellSpec { &self.spec }

    /// Step-mean pressure = `direct * U + history` from the accepted states.
    #[must_use]
    pub fn affine_pressure(&self) -> (f64, f64) {
        let mut history = 0.0;
        for (i, cell) in self.coefficients.iter().enumerate().take(self.spec.count) {
            history += cell.direct * (cell.charge * self.pressures[i] - self.flows[i]);
        }
        (self.direct, history)
    }

    /// Accepted inductor flows [m^3/s].
    #[must_use]
    pub fn inductor_flows_m3_s(&self) -> &[f64] { &self.flows[..self.spec.count] }

    /// Accepted cell pressures [Pa].
    #[must_use]
    pub fn pressures_pa(&self) -> &[f64] { &self.pressures[..self.spec.count] }

    fn energy(&self, flows: &[f64; MAX_PARALLEL_CELLS], pressures: &[f64; MAX_PARALLEL_CELLS]) -> f64 {
        let mut sum = 0.0;
        for (i, cell) in self.coefficients.iter().enumerate().take(self.spec.count) {
            let inductive = cell.sqrt_inertance * flows[i];
            let capacitive = cell.sqrt_compliance * pressures[i];
            sum += 0.5 * inductive * inductive + 0.5 * capacitive * capacitive;
        }
        sum
    }

    /// Stored energy of the accepted states [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 { self.energy(&self.flows, &self.pressures) }

    /// Preview one step at the solved step-mean flow without changing state.
    ///
    /// # Errors
    /// Nonfinite flow, state, storage or work.
    pub fn preview_step(&self, flow: f64) -> Result<ParallelCellFrame, WaveguideError> {
        if !flow.is_finite() {
            return Err(WaveguideError("parallel-cell load requires a finite flow"));
        }
        let (mut flows, mut pressures) = (self.flows, self.pressures);
        let (mut pressure, mut dissipated) = (0.0, 0.0);
        for (i, cell) in self.spec.cells().iter().enumerate() {
            let c = self.coefficients[i];
            let midpoint = c.direct * (flow + c.charge * self.pressures[i] - self.flows[i]);
            flows[i] += c.advance * midpoint;
            pressures[i] = if cell.compliance_m3_pa > 0.0 { 2.0 * midpoint - self.pressures[i] } else { midpoint };
            pressure += midpoint;
            dissipated += (midpoint / cell.resistance_pa_s_m3) * (midpoint * self.dt);
            if ![midpoint, flows[i], pressures[i]].iter().all(|v| v.is_finite()) {
                return Err(WaveguideError("parallel-cell state left the finite set"));
            }
        }
        let stored = self.energy(&flows, &pressures);
        let frame = ParallelCellFrame { pressure_pa: pressure, flow_m3_s: flow, stored_energy_j: stored,
            storage_change_j: stored - self.stored_energy_j(), dissipated_energy_j: dissipated,
            supplied_work_j: pressure * (flow * self.dt), flows, pressures, count: self.spec.count };
        if ![pressure, dissipated, stored, frame.storage_change_j, frame.supplied_work_j].iter().all(|v| v.is_finite()) {
            return Err(WaveguideError("parallel-cell observation left the finite set"));
        }
        Ok(frame)
    }

    /// Advance only after a complete preview succeeds; a refusal keeps all history.
    ///
    /// # Errors
    /// Same refusal boundary as `preview_step`.
    pub fn step(&mut self, flow: f64) -> Result<ParallelCellFrame, WaveguideError> {
        let frame = self.preview_step(flow)?;
        self.flows = frame.flows;
        self.pressures = frame.pressures;
        Ok(frame)
    }
}
