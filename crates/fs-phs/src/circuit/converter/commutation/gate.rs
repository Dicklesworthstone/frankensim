//! Original channel-current laws and finite-impedance gate drivers.
//!
//! Command bits select an authored driver branch. Actual drain/source and
//! gate/source voltages select a memoryless channel current and its two
//! analytical partials. The original positive G(Vgs) table remains one such
//! law; complete authored I(Vds,Vgs) responses can retain saturation and
//! asymmetric reverse conduction without dividing current by drain voltage.
//! Channel and body-diode currents have separate ownership. Carrier memory,
//! avalanche and extrapolation beyond a source's applicability are not inferred.
//!
//! Gate drivers are Thevenin sources with finite resistance and explicit current
//! support. Their signed source work, gate-terminal work and positive resistor
//! heat are separate. Finite source impedance does not imply finite rail energy.

use super::{invalid, number, positive};
use crate::PhsError;
use fs_math::dd::Dd;
use std::sync::Arc;

/// A memoryless, immutable channel-only current response on its original domain.
///
/// Drain/source current is positive from drain to source. The physical channel
/// must obey I*Vds >= 0 and I(0,Vgs) = 0 throughout its admitted support. Zero
/// current and zero drain derivative are valid authored cutoff/saturation
/// states. Neither a leakage floor nor a saturation threshold is invented.
///
/// The returned derivatives must be those of the SAME current evaluator. An
/// original piecewise law may select one deterministic incident cell at a knee;
/// that is a semismooth Newton branch, not a classical smoothness claim.
/// Implementations retain and enforce every drain/gate support condition,
/// allocate no source query history, and return the same response on repeated
/// calls. They must exclude any body-diode current already owned by the device.
/// Higher layers own material selection and temperature/phase applicability.
pub trait GateChannelResponse: std::fmt::Debug + Send + Sync {
    fn gate_voltage_range_v(&self) -> [f64; 2];

    fn evaluate(
        &self,
        drain_source_voltage_v: f64,
        gate_source_voltage_v: f64,
    ) -> Result<GateChannelPoint, PhsError>;
}

/// Original channel current and analytic active-cell voltage derivatives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GateChannelPoint {
    pub current_a: f64,
    pub drain_derivative_s: f64,
    pub gate_derivative_s: f64,
}

impl GateChannelPoint {
    pub(super) fn admit(self, drain_source_voltage_v: f64) -> Result<Self, PhsError> {
        if [
            drain_source_voltage_v,
            self.current_a,
            self.drain_derivative_s,
            self.gate_derivative_s,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite original channel current and voltage partials required",
            ));
        }
        if drain_source_voltage_v == 0.0 && self.current_a != 0.0 {
            return Err(invalid(
                "a passive channel cannot supply current at zero drain voltage",
            ));
        }
        let power = number(Dd::from_f64(drain_source_voltage_v) * Dd::from_f64(self.current_a))?;
        if power < 0.0 || (self.current_a != 0.0 && power == 0.0) {
            return Err(invalid(
                "passive physical channel work is unrepresentable or negative",
            ));
        }
        Ok(self)
    }
}

/// Cold-erased immutable source response shared by all native gate intervals.
///
/// Equality identifies the same retained allocation, never an inferred equality
/// of arbitrary constitutive functions. Independently prepared material laws
/// compare their complete original numerical data in their own thermal cache.
#[derive(Clone, Debug)]
pub struct GateChannelLaw {
    response: Arc<dyn GateChannelResponse>,
    gate_voltage_range_v: [f64; 2],
}

impl PartialEq for GateChannelLaw {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.response, &other.response)
            && self.gate_voltage_range_v.map(f64::to_bits)
                == other.gate_voltage_range_v.map(f64::to_bits)
    }
}

impl GateChannelLaw {
    pub fn new<C: GateChannelResponse + 'static>(response: C) -> Result<Self, PhsError> {
        let range = response.gate_voltage_range_v();
        if range.iter().any(|v| !v.is_finite())
            || range[0] >= range[1]
            || !positive(number(Dd::from_f64(range[1]) - Dd::from_f64(range[0]))?)
        {
            return Err(invalid(
                "finite ordered original channel gate support required",
            ));
        }
        Ok(Self {
            response: Arc::new(response),
            gate_voltage_range_v: range,
        })
    }

    pub const fn gate_voltage_range_v(&self) -> [f64; 2] {
        self.gate_voltage_range_v
    }

    pub fn evaluate(
        &self,
        drain_source_voltage_v: f64,
        gate_source_voltage_v: f64,
    ) -> Result<GateChannelPoint, PhsError> {
        let gate = gate_source_voltage_v;
        let range = self.gate_voltage_range_v;
        if !drain_source_voltage_v.is_finite()
            || !gate.is_finite()
            || gate < range[0]
            || gate > range[1]
        {
            return Err(invalid(
                "actual channel voltage leaves original source support",
            ));
        }
        self.response
            .evaluate(drain_source_voltage_v, gate)?
            .admit(drain_source_voltage_v)
    }
}

impl GateChannelResponse for GateChannelLaw {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        GateChannelLaw::gate_voltage_range_v(self)
    }

    fn evaluate(&self, vds: f64, vgs: f64) -> Result<GateChannelPoint, PhsError> {
        GateChannelLaw::evaluate(self, vds, vgs)
    }
}

/// Bounds retained constitutive data and each binary search.
pub const MAX_GATE_CONDUCTANCE_KNOTS: usize = 4096;

#[derive(Clone, Debug, PartialEq)]
struct ConductanceData {
    knots: Vec<(f64, f64)>,
    slopes: Vec<f64>,
}

/// Complete original (gate-source voltage [V], channel conductance [S]) law.
///
/// Conductance is strictly positive throughout every retained affine segment.
/// Its gate derivative may have either sign. Clones share the immutable original
/// table; evaluation allocates nothing and never extrapolates or clips.
#[derive(Clone, Debug, PartialEq)]
pub struct GateConductanceCurve {
    data: Arc<ConductanceData>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GateConductancePoint {
    pub conductance_s: f64,
    /// Exact affine-cell derivative. Interior knots select the right cell;
    /// the final support endpoint selects its left cell.
    pub slope_s_per_v: f64,
}

impl GateConductanceCurve {
    pub fn try_new(knots: Vec<(f64, f64)>) -> Result<Self, PhsError> {
        if !(2..=MAX_GATE_CONDUCTANCE_KNOTS).contains(&knots.len())
            || knots.iter().any(|&(v, g)| !v.is_finite() || !positive(g))
        {
            return Err(invalid(
                "bounded original finite gate-voltage and positive conductance knots required",
            ));
        }
        let mut slopes = Vec::with_capacity(knots.len() - 1);
        for pair in knots.windows(2) {
            if pair[0].0 >= pair[1].0 {
                return Err(invalid(
                    "original gate-voltage knots must be strictly ordered",
                ));
            }
            let width = Dd::from_f64(pair[1].0) - Dd::from_f64(pair[0].0);
            if !positive(number(width)?) {
                return Err(invalid("original gate-voltage span is unrepresentable"));
            }
            let delta = Dd::from_f64(pair[1].1) - Dd::from_f64(pair[0].1);
            let slope = number(delta / width)?;
            if pair[1].1 != pair[0].1 && slope == 0.0 {
                return Err(invalid("original gate-conductance slope underflows"));
            }
            slopes.push(slope);
        }
        Ok(Self {
            data: Arc::new(ConductanceData { knots, slopes }),
        })
    }

    pub fn knots(&self) -> &[(f64, f64)] {
        &self.data.knots
    }

    pub fn gate_voltage_range_v(&self) -> [f64; 2] {
        let k = self.knots();
        [k[0].0, k[k.len() - 1].0]
    }

    pub fn evaluate(&self, gate_source_voltage_v: f64) -> Result<GateConductancePoint, PhsError> {
        let k = self.knots();
        let v = gate_source_voltage_v;
        if !v.is_finite() || v < k[0].0 || v > k[k.len() - 1].0 {
            return Err(invalid(
                "actual gate voltage leaves original channel-law support",
            ));
        }
        let i = k
            .partition_point(|&(x, _)| x <= v)
            .saturating_sub(1)
            .min(k.len() - 2);
        let g = if v == k[i].0 {
            k[i].1
        } else if v == k[i + 1].0 {
            k[i + 1].1
        } else {
            let t = (Dd::from_f64(v) - Dd::from_f64(k[i].0))
                / (Dd::from_f64(k[i + 1].0) - Dd::from_f64(k[i].0));
            number((Dd::from_f64(1.0) - t) * Dd::from_f64(k[i].1) + t * Dd::from_f64(k[i + 1].1))?
        };
        if !positive(g) {
            return Err(invalid(
                "original positive gate conductance is unrepresentable",
            ));
        }
        Ok(GateConductancePoint {
            conductance_s: g,
            slope_s_per_v: self.data.slopes[i],
        })
    }
}

impl GateChannelResponse for GateConductanceCurve {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        GateConductanceCurve::gate_voltage_range_v(self)
    }

    fn evaluate(&self, vds: f64, vgs: f64) -> Result<GateChannelPoint, PhsError> {
        if !vds.is_finite() {
            return Err(invalid("finite original channel drain voltage required"));
        }
        let point = GateConductanceCurve::evaluate(self, vgs)?;
        let current_a = number(Dd::from_f64(vds) * Dd::from_f64(point.conductance_s))?;
        if vds != 0.0 && current_a == 0.0 {
            return Err(invalid(
                "original nonzero gate-controlled channel current underflows",
            ));
        }
        GateChannelPoint {
            current_a,
            drain_derivative_s: point.conductance_s,
            gate_derivative_s: number(Dd::from_f64(vds) * Dd::from_f64(point.slope_s_per_v))?,
        }
        .admit(vds)
    }
}

/// The on/off labels refer to command bits, not numerical voltage ordering.
/// Both branches have explicit finite impedance. The channel responds to actual
/// Vgs, even during a command transition or a Miller-induced gate excursion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GateDriver {
    pub on_voltage_v: f64,
    pub off_voltage_v: f64,
    pub on_resistance_ohm: f64,
    pub off_resistance_ohm: f64,
    /// Source-data operating support; exceeding it refuses rather than clamps.
    pub maximum_current_a: f64,
}

/// An immutable admitted branch of the original driver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GateDriverBranch {
    source_voltage_v: f64,
    resistance_ohm: f64,
    maximum_current_a: f64,
    conductance_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GateDriverPoint {
    pub source_voltage_v: f64,
    pub gate_source_voltage_v: f64,
    /// Signed current delivered into the gate charge coordinate.
    pub current_a: f64,
    /// Signed power from the ideal source behind the finite driver resistance.
    pub source_power_w: f64,
    /// Signed power delivered to the physical gate terminal.
    pub terminal_power_w: f64,
    pub resistor_power_w: f64,
    /// Positive G=1/R; dI/dVgs is -G.
    pub conductance_s: f64,
    pub(super) current: Dd,
    pub(super) supplied: Dd,
    pub(super) delivered: Dd,
    pub(super) loss: Dd,
}

impl GateDriver {
    pub fn validate(&self) -> Result<(), PhsError> {
        if !self.on_voltage_v.is_finite()
            || !self.off_voltage_v.is_finite()
            || !positive(self.maximum_current_a)
            || [self.on_resistance_ohm, self.off_resistance_ohm]
                .into_iter()
                .any(|r| !positive(r) || !positive(1.0 / r))
        {
            return Err(invalid(
                "gate driver needs finite authored voltages, positive representable branch resistances and current support",
            ));
        }
        Ok(())
    }

    pub fn branch(&self, on: bool) -> Result<GateDriverBranch, PhsError> {
        self.validate()?;
        let (source_voltage_v, resistance_ohm) = if on {
            (self.on_voltage_v, self.on_resistance_ohm)
        } else {
            (self.off_voltage_v, self.off_resistance_ohm)
        };
        Ok(GateDriverBranch {
            source_voltage_v,
            resistance_ohm,
            maximum_current_a: self.maximum_current_a,
            conductance_s: 1.0 / resistance_ohm,
        })
    }

    pub fn evaluate(
        &self,
        on: bool,
        gate_source_voltage_v: f64,
    ) -> Result<GateDriverPoint, PhsError> {
        self.branch(on)?.evaluate(gate_source_voltage_v)
    }
}

impl GateDriverBranch {
    pub const fn source_voltage_v(&self) -> f64 {
        self.source_voltage_v
    }
    pub const fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }
    pub const fn maximum_current_a(&self) -> f64 {
        self.maximum_current_a
    }

    pub fn evaluate(&self, gate_source_voltage_v: f64) -> Result<GateDriverPoint, PhsError> {
        if !gate_source_voltage_v.is_finite() {
            return Err(invalid("finite actual gate-source voltage required"));
        }
        let source = Dd::from_f64(self.source_voltage_v);
        let gate = Dd::from_f64(gate_source_voltage_v);
        let drop = source - gate;
        number(drop)?;
        let current_a = number(drop / Dd::from_f64(self.resistance_ohm))?;
        if current_a.abs() > self.maximum_current_a {
            return Err(invalid(
                "actual gate-driver current leaves authored support",
            ));
        }
        if (drop.hi != 0.0 || drop.lo != 0.0) && current_a == 0.0 {
            return Err(invalid("nonzero gate-driver current underflows"));
        }
        // The native state flow is this represented current. Use the same
        // current in every power term, especially near a balanced biased gate.
        let current = Dd::from_f64(current_a);
        let supplied = source * current;
        let delivered = gate * current;
        let loss = drop * current;
        let source_power_w = number(supplied)?;
        let terminal_power_w = number(delivered)?;
        let resistor_power_w = number(loss)?;
        if resistor_power_w < 0.0 || (current_a != 0.0 && resistor_power_w == 0.0) {
            return Err(invalid(
                "positive physical gate-driver heat is unrepresentable",
            ));
        }
        Ok(GateDriverPoint {
            source_voltage_v: self.source_voltage_v,
            gate_source_voltage_v,
            current_a,
            source_power_w,
            terminal_power_w,
            resistor_power_w,
            conductance_s: self.conductance_s,
            current,
            supplied,
            delivered,
            loss,
        })
    }
}
