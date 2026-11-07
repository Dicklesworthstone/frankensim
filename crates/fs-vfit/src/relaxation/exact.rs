//! Exact continuous cell dynamics under an explicitly held series flow.

use super::{MAX_PARALLEL_CELLS, ParallelCellFrame, ParallelCellSpec};
use crate::waveguide::WaveguideError;
use fs_math::c64::C64;
use fs_modal::time::{LinearForceTransition, ModalAcousticState};

#[derive(Clone, Copy, Debug)]
enum ExactCell {
    Rl {
        retention: f64,
        decrement: f64,
        mean_pressure_gain: f64,
        sqrt_heat_gain: f64,
    },
    Rlc {
        transition: LinearForceTransition,
        coordinate_per_flow: f64,
        velocity_per_pressure: f64,
    },
}

/// Flow-driven series chain of parallel RLC cells with exact held-flow decay.
///
/// The caller holds the solved series flow constant throughout each interval.
/// The returned pressure is its continuous interval mean, so `U * p_mean * dt`
/// is the physical input work. The resistor loss is independently integrated
/// as a nonnegative quadratic form; it is not obtained by subtracting storage
/// from work. Critical damping, strong overdamping and absent resistors share
/// the same oscillator kernel as modal acoustics.
///
/// An RL cell's decay is `exp(-R*dt/L)`. Consequently a fast physical relaxation
/// vanishes within a coarse step, rather than retaining the alternating tail
/// `(1-R*dt/(2L))/(1+R*dt/(2L))` of implicit midpoint. This does not recover
/// variation in the input between its declared held values.
///
/// Preparation caches all transitions. Stepping uses bounded inline storage,
/// allocates nothing, and evaluates no transcendental functions. This owner
/// is separate from [`super::ParallelCellLoad`] so existing bilinear fits keep
/// their sampling convention and their smaller midpoint state.
#[derive(Clone, Copy, Debug)]
pub struct ExactParallelCellLoad {
    spec: ParallelCellSpec,
    dt: f64,
    coefficients: [Option<ExactCell>; MAX_PARALLEL_CELLS],
    sqrt_inertances: [f64; MAX_PARALLEL_CELLS],
    sqrt_compliances: [f64; MAX_PARALLEL_CELLS],
    flows: [f64; MAX_PARALLEL_CELLS],
    pressures: [f64; MAX_PARALLEL_CELLS],
    direct: f64,
}

impl ExactParallelCellLoad {
    /// Prepare a zero-storage chain for positive finite held-flow intervals.
    ///
    /// # Errors
    /// An invalid interval or an unrepresentable physical transition/scaling.
    pub fn new(spec: ParallelCellSpec, dt: f64) -> Result<Self, WaveguideError> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(WaveguideError(
                "an exact parallel-cell load requires a positive finite time step",
            ));
        }
        let mut out = Self {
            spec,
            dt,
            coefficients: [None; MAX_PARALLEL_CELLS],
            sqrt_inertances: [0.0; MAX_PARALLEL_CELLS],
            sqrt_compliances: [0.0; MAX_PARALLEL_CELLS],
            flows: [0.0; MAX_PARALLEL_CELLS],
            pressures: [0.0; MAX_PARALLEL_CELLS],
            direct: 0.0,
        };
        for (i, cell) in spec.cells().iter().enumerate() {
            out.sqrt_inertances[i] = fs_math::det::sqrt(cell.inertance_pa_s2_m3);
            out.sqrt_compliances[i] = fs_math::det::sqrt(cell.compliance_m3_pa);
            let (coefficient, direct) = if cell.compliance_m3_pa == 0.0 {
                let drag = (cell.resistance_pa_s_m3 / cell.inertance_pa_s2_m3) * dt;
                if !drag.is_finite() || drag <= 0.0 {
                    return Err(WaveguideError(
                        "exact RL cell drag must be positive and representable",
                    ));
                }
                let retention = fs_math::det::exp(-drag);
                let decrement = -fs_math::det::expm1(-drag);
                let mean_pressure_gain = cell.resistance_pa_s_m3 * (decrement / drag);
                let sqrt_heat_gain = out.sqrt_inertances[i]
                    * fs_math::det::sqrt(0.5 * decrement * (1.0 + retention));
                if ![mean_pressure_gain, sqrt_heat_gain]
                    .into_iter()
                    .all(|v| v.is_finite() && v > 0.0)
                {
                    return Err(WaveguideError(
                        "exact RL cell pressure and heat coefficients must be representable",
                    ));
                }
                (
                    ExactCell::Rl {
                        retention,
                        decrement,
                        mean_pressure_gain,
                        sqrt_heat_gain,
                    },
                    mean_pressure_gain,
                )
            } else {
                let sqrt_l = out.sqrt_inertances[i];
                let sqrt_c = out.sqrt_compliances[i];
                let frequency = (1.0 / sqrt_l) / sqrt_c;
                let damping = (1.0 / cell.resistance_pa_s_m3) / cell.compliance_m3_pa;
                let coordinate_per_flow = cell.inertance_pa_s2_m3 * sqrt_c;
                if !frequency.is_finite()
                    || frequency <= 0.0
                    || !coordinate_per_flow.is_finite()
                    || coordinate_per_flow <= 0.0
                    || !(1.0 / coordinate_per_flow).is_finite()
                {
                    return Err(WaveguideError(
                        "exact RLC cell frequency and coordinates must be representable",
                    ));
                }
                let transition =
                    LinearForceTransition::try_new_with_damping_rate(frequency, damping, dt)
                        .map_err(|_| {
                            WaveguideError("exact RLC cell transition is not representable")
                        })?;
                let (_, mobility) = transition
                    .held_force_mean_velocity_relation(ModalAcousticState::default(), 0.0)
                    .map_err(|_| WaveguideError("exact RLC cell mobility is not representable"))?;
                let direct = (mobility / sqrt_c) / sqrt_c;
                if !direct.is_finite() || direct < 0.0 {
                    return Err(WaveguideError(
                        "exact RLC cell mean-pressure response must be finite and passive",
                    ));
                }
                (
                    ExactCell::Rlc {
                        transition,
                        coordinate_per_flow,
                        velocity_per_pressure: sqrt_c,
                    },
                    direct,
                )
            };
            out.coefficients[i] = Some(coefficient);
            out.direct += direct;
        }
        if !out.direct.is_finite() {
            return Err(WaveguideError(
                "exact parallel-cell chain response is not representable",
            ));
        }
        Ok(out)
    }

    /// Immutable continuous cell declaration, in the caller's order.
    #[must_use]
    pub const fn spec(&self) -> &ParallelCellSpec {
        &self.spec
    }

    /// Steady sinusoidal impedance from held flow to interval-mean pressure.
    ///
    /// The input sequence is `exp(i*2*pi*f*n*dt)` and each value is held over
    /// its interval. This is the sampled transfer of this exact owner, not a
    /// bilinear frequency warp or a sample of endpoint pressure. It is
    /// independent of the currently accepted history. DC is exactly zero:
    /// every cell has an inductive shunt.
    ///
    /// # Errors
    /// A frequency outside `[0, Nyquist]`, a lossless sampled pole, or an
    /// unrepresentable response. Construction does not impose a Nyquist guard
    /// on the internal poles; this query only bounds the input frequency.
    pub fn sampled_impedance(&self, frequency_hz: f64) -> Result<C64, WaveguideError> {
        if !frequency_hz.is_finite() || frequency_hz < 0.0 || frequency_hz > 0.5 / self.dt {
            return Err(WaveguideError(
                "exact parallel-cell frequency must lie between DC and Nyquist",
            ));
        }
        if frequency_hz == 0.0 {
            return Ok(C64::ZERO);
        }
        let angle = core::f64::consts::TAU * (frequency_hz * self.dt);
        let half_sine = fs_math::det::sin(0.5 * angle);
        // z - 1 without subtracting nearly equal cosines near DC.
        let s = C64::new(-2.0 * half_sine * half_sine, fs_math::det::sin(angle));
        let mut impedance = C64::ZERO;
        for (i, cell) in self.spec.cells().iter().enumerate() {
            impedance = impedance
                + match self.coefficients[i] {
                    Some(ExactCell::Rl {
                        decrement,
                        mean_pressure_gain,
                        ..
                    }) => (s / (s + C64::from_re(decrement))).scale(mean_pressure_gain),
                    Some(ExactCell::Rlc { transition, .. }) => {
                        let omega = (1.0 / self.sqrt_inertances[i]) / self.sqrt_compliances[i];
                        let damping = (1.0 / cell.resistance_pa_s_m3) / cell.compliance_m3_pa;
                        let probe_force = (1.0 / self.dt).min(1.0);
                        let response = transition
                            .advance(ModalAcousticState::default(), probe_force, probe_force)
                            .map_err(|_| {
                                WaveguideError("exact RLC sampled transition is not representable")
                            })?;
                        let (_, mobility) = transition
                            .held_force_mean_velocity_relation(ModalAcousticState::default(), 0.0)
                            .map_err(|_| {
                                WaveguideError("exact RLC sampled mobility is not representable")
                            })?;
                        // E in energy coordinates has diagonal (1-dq, 1-dv)
                        // and off-diagonal (+b, -b). Keep the unscaled impulse
                        // response B: omega*B may underflow even when c*B and
                        // the physical pressure response remain representable.
                        let impulse_response = response.velocity_m_sqrt_kg_per_s / probe_force;
                        let direct =
                            (mobility / self.sqrt_compliances[i]) / self.sqrt_compliances[i];
                        let dq = (omega * self.dt) * (omega * mobility);
                        let b = omega * impulse_response;
                        let dv = dq + damping * impulse_response;
                        let determinant = dq * dv + b * b;
                        // L*Delta(i) = dt*p_mean. Factor the exact DC zero
                        // and keep its numerator in pressure units: multiplying
                        // tiny dimensionless decrements by L/dt afterwards can
                        // erase a representable resistance through underflow.
                        let pressure_impulse = impulse_response / self.sqrt_compliances[i];
                        let pressure_determinant =
                            direct * dv + pressure_impulse * (pressure_impulse / self.dt);
                        let numerator = s * (C64::from_re(pressure_determinant) + s.scale(direct));
                        let denominator = C64::from_re(determinant) + s.scale(dq + dv) + s * s;
                        if cell.resistance_pa_s_m3.is_infinite() {
                            let scale =
                                determinant.abs() + s.abs() * (dq.abs() + dv.abs() + s.abs());
                            if denominator.abs() <= 64.0 * f64::EPSILON * scale {
                                return Err(WaveguideError("exact LC sampled pole is unresolved"));
                            }
                        }
                        numerator / denominator
                    }
                    None => unreachable!("all declared cell coefficients are prepared"),
                };
        }
        if !impedance.re.is_finite() || !impedance.im.is_finite() {
            return Err(WaveguideError(
                "exact parallel-cell sampled impedance is not representable",
            ));
        }
        Ok(impedance)
    }

    /// Accepted inductor flows [m^3/s].
    #[must_use]
    pub fn inductor_flows_m3_s(&self) -> &[f64] {
        &self.flows[..self.spec.count]
    }

    /// Accepted endpoint pressures [Pa], including an RL cell's algebraic pressure.
    #[must_use]
    pub fn pressures_pa(&self) -> &[f64] {
        &self.pressures[..self.spec.count]
    }

    fn energy(
        &self,
        flows: &[f64; MAX_PARALLEL_CELLS],
        pressures: &[f64; MAX_PARALLEL_CELLS],
    ) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.spec.count {
            let inductive = self.sqrt_inertances[i] * flows[i];
            let capacitive = self.sqrt_compliances[i] * pressures[i];
            energy += 0.5 * inductive * inductive + 0.5 * capacitive * capacitive;
        }
        energy
    }

    /// Physical accepted storage `sum (L*i^2+C*p^2)/2` [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.energy(&self.flows, &self.pressures)
    }

    /// Exact same-window mean pressure `direct * U + history` [Pa].
    ///
    /// Use this relation in the caller's implicit flow solve. No state advances.
    ///
    /// # Errors
    /// An accepted state has an unrepresentable mean-pressure response.
    pub fn affine_pressure(&self) -> Result<(f64, f64), WaveguideError> {
        let mut history = 0.0;
        for (i, coefficient) in self.coefficients.iter().enumerate().take(self.spec.count) {
            history += match coefficient {
                Some(ExactCell::Rl {
                    mean_pressure_gain, ..
                }) => -mean_pressure_gain * self.flows[i],
                Some(ExactCell::Rlc {
                    transition,
                    coordinate_per_flow,
                    velocity_per_pressure,
                }) => {
                    let state = ModalAcousticState {
                        displacement_m_sqrt_kg: coordinate_per_flow * self.flows[i],
                        velocity_m_sqrt_kg_per_s: velocity_per_pressure * self.pressures[i],
                    };
                    let (mean, _) = transition
                        .held_force_mean_velocity_relation(state, 0.0)
                        .map_err(|_| {
                            WaveguideError("exact parallel-cell history is not representable")
                        })?;
                    mean / velocity_per_pressure
                }
                None => unreachable!("all declared cell coefficients are prepared"),
            };
        }
        if !history.is_finite() {
            return Err(WaveguideError(
                "exact parallel-cell history is not representable",
            ));
        }
        Ok((self.direct, history))
    }

    /// Preview one interval at its constant solved flow, preserving all history.
    ///
    /// # Errors
    /// Nonfinite flow or unrepresentable candidate state, pressure, energy or work.
    pub fn preview_step(&self, flow: f64) -> Result<ParallelCellFrame, WaveguideError> {
        if !flow.is_finite() {
            return Err(WaveguideError(
                "exact parallel-cell load requires a finite held flow",
            ));
        }
        let (mut flows, mut pressures) = (self.flows, self.pressures);
        let (mut pressure, mut heat) = (0.0, 0.0);
        for (i, cell) in self.spec.cells().iter().enumerate() {
            match self.coefficients[i] {
                Some(ExactCell::Rl {
                    retention,
                    decrement,
                    mean_pressure_gain,
                    sqrt_heat_gain,
                }) => {
                    let imbalance = flow - self.flows[i];
                    flows[i] = if decrement < 0.5 {
                        self.flows[i] + decrement * imbalance
                    } else {
                        flow - retention * imbalance
                    };
                    pressures[i] = cell.resistance_pa_s_m3 * (retention * imbalance);
                    pressure += mean_pressure_gain * imbalance;
                    let root = sqrt_heat_gain * imbalance;
                    heat += root * root;
                }
                Some(ExactCell::Rlc {
                    transition,
                    coordinate_per_flow,
                    velocity_per_pressure,
                }) => {
                    let state = ModalAcousticState {
                        displacement_m_sqrt_kg: coordinate_per_flow * self.flows[i],
                        velocity_m_sqrt_kg_per_s: velocity_per_pressure * self.pressures[i],
                    };
                    let force = flow / velocity_per_pressure;
                    let result = transition
                        .advance_with_energy(state, force, force)
                        .map_err(|_| {
                            WaveguideError("exact parallel-cell trajectory is not representable")
                        })?;
                    let (mean, _) = transition
                        .held_force_mean_velocity_relation(state, force)
                        .map_err(|_| {
                            WaveguideError("exact parallel-cell mean pressure is not representable")
                        })?;
                    flows[i] = result.state.displacement_m_sqrt_kg / coordinate_per_flow;
                    pressures[i] = result.state.velocity_m_sqrt_kg_per_s / velocity_per_pressure;
                    pressure += mean / velocity_per_pressure;
                    heat += result.viscous_dissipation_j;
                }
                None => unreachable!("all declared cell coefficients are prepared"),
            }
            if ![flows[i], pressures[i]].into_iter().all(f64::is_finite) {
                return Err(WaveguideError(
                    "exact parallel-cell state left the finite set",
                ));
            }
        }
        let stored = self.energy(&flows, &pressures);
        let frame = ParallelCellFrame {
            pressure_pa: pressure,
            flow_m3_s: flow,
            stored_energy_j: stored,
            storage_change_j: stored - self.stored_energy_j(),
            dissipated_energy_j: heat,
            supplied_work_j: pressure * (flow * self.dt),
            flows,
            pressures,
            count: self.spec.count,
        };
        if ![
            pressure,
            heat,
            stored,
            frame.storage_change_j,
            frame.supplied_work_j,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(WaveguideError(
                "exact parallel-cell observation left the finite set",
            ));
        }
        Ok(frame)
    }

    /// Accept the complete candidate only after every cell and observation succeeds.
    ///
    /// # Errors
    /// Same refusal boundary as [`Self::preview_step`]; a refusal retains history.
    pub fn step(&mut self, flow: f64) -> Result<ParallelCellFrame, WaveguideError> {
        let frame = self.preview_step(flow)?;
        self.flows = frame.flows;
        self.pressures = frame.pressures;
        Ok(frame)
    }
}
