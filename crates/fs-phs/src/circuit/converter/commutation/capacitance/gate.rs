//! Reciprocal intrinsic drain, gate/source and Miller charge in one storage.
//!
//! The nine coordinates are lambda, two bus electrode charges, two GROUPED
//! node charges, and four gate charges. A grouped node charge includes its
//! high-side gate electrode: its conjugate gate effort is Vgs directly.
//! Intrinsic Qds must be authored separately from Cgs/Cgd. Gate-clamped Coss
//! cannot be augmented with Cgd again without counting the same storage twice.
//!
//! Constant positive Cgs and nonnegative Cgd permit an exact Schur coordinate
//! transformation into the existing nonlinear two-electrode charge owner.
//! Its original inverse, AVF cell tracer and analytic tangent are reused.
//! No gate delay, switching-loss map, carrier memory or extra solver is added.

use super::charge::{
    ChargeCommutationParameters, ChargeCommutationStorage, change_curve, curve_charge_dd,
};
use super::{CommutationParameters, invalid, number, positive};
use crate::circuit::capacitance::CapacitorChargeCurve;
use crate::circuit::converter::commutation::gate::{
    GateChannelLaw, GateChannelResponse, GateConductanceCurve, GateDriver,
};
use crate::{ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, Storage};
use fs_math::dd::Dd;
fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}

#[derive(Clone, Debug, PartialEq)]
pub struct GateCommutationParameters<C = GateConductanceCurve> {
    /// Original diode/bounds cards are retained. Their Ron/Roff do not select
    /// driven-channel current: channel_conductance is the sole gate-channel law.
    pub converter: CommutationParameters,
    pub intrinsic_drain_charge_curves: [CapacitorChargeCurve; 4],
    pub gate_source_capacitance_f: [f64; 4],
    pub gate_drain_capacitance_f: [f64; 4],
    /// Original channel-only responses in input-high/input-low/output-high/
    /// output-low order. The historical field name is retained for existing
    /// conductance-curve callers; C may also be a full signed I(Vds, Vgs) law.
    /// Body-diode current must not be included in this response.
    pub channel_conductance: [C; 4],
    pub gate_drivers: [GateDriver; 4],
    pub maximum_cells_per_chord: usize,
}
impl<C: GateChannelResponse + Clone + 'static> GateCommutationParameters<C> {
    pub fn validate(&self) -> Result<(), PhsError> {
        GateCommutationStorage::new(self).map(|_| ())
    }
    /// Construct physical charges for this exact original terminal decomposition.
    /// The old five-state node charges are not grouped charges and cannot be
    /// copied unchanged. Gate voltages are source referenced, including high side.
    pub fn state(
        &self,
        current_a: f64,
        bus_voltage_v: [f64; 2],
        node_voltage_v: [f64; 2],
        gate_source_voltage_v: [f64; 4],
    ) -> Result<[f64; 9], PhsError> {
        GateCommutationStorage::new(self)?.state(
            current_a,
            bus_voltage_v,
            node_voltage_v,
            gate_source_voltage_v,
        )
    }
    pub fn stored_energy_j(&self, state: &[f64; 9]) -> Result<f64, PhsError> {
        GateCommutationStorage::new(self)?.energy(state)
    }
}

/// Independently retained reversible gate energy, finite-source work and heat.
#[derive(Clone, Debug)]
pub struct GatePowerAccount {
    pub gate_source_energy_change_j: [f64; 4],
    pub gate_drain_energy_change_j: [f64; 4],
    pub source_work_j: [f64; 4],
    pub resistor_heat_j: [f64; 4],
    /// Work-conjugate mean Vgs, not absolute high-side gate potential.
    pub mean_voltage_v: [f64; 4],
    pub mean_current_a: [f64; 4],
    pub last_driver_current_a: [f64; 4],
}

#[derive(Clone, Debug)]
pub struct GateCommutationStorage {
    parameters: GateCommutationParameters<GateChannelLaw>,
    base: ChargeCommutationStorage,
    sum_f: [Dd; 4],
    alpha: [Dd; 4],
    scales: [f64; 9],
}
impl GateCommutationStorage {
    /// Prepare reversible storage once while retaining the same immutable
    /// channel owner for actual-voltage support checks. Channel response changes
    /// neither the charge transform nor the original AVF energy rule.
    pub fn new<C: GateChannelResponse + Clone + 'static>(
        parameters: &GateCommutationParameters<C>,
    ) -> Result<Self, PhsError> {
        let parameters = GateCommutationParameters {
            converter: parameters.converter,
            intrinsic_drain_charge_curves: parameters.intrinsic_drain_charge_curves.clone(),
            gate_source_capacitance_f: parameters.gate_source_capacitance_f,
            gate_drain_capacitance_f: parameters.gate_drain_capacitance_f,
            channel_conductance: super::four(core::array::from_fn(|k| {
                GateChannelLaw::new(parameters.channel_conductance[k].clone())
            }))?,
            gate_drivers: parameters.gate_drivers,
            maximum_cells_per_chord: parameters.maximum_cells_per_chord,
        };
        parameters.converter.validate()?;
        let mut sum_f = [Dd::ZERO; 4];
        let mut alpha = [Dd::ZERO; 4];
        let mut curves = parameters.intrinsic_drain_charge_curves.clone();
        for k in 0..4 {
            parameters.gate_drivers[k].validate()?;
            let gs = parameters.gate_source_capacitance_f[k];
            let gd = parameters.gate_drain_capacitance_f[k];
            if !positive(gs) || !gd.is_finite() || gd < 0.0 {
                return Err(invalid(
                    "authored positive Cgs and nonnegative Cgd required",
                ));
            }
            sum_f[k] = dd(gs) + dd(gd);
            if !positive(number(sum_f[k])?) || !positive(number(dd(1.0) / sum_f[k])?) {
                return Err(invalid("gate storage sum or reciprocal is unrepresentable"));
            }
            alpha[k] = dd(gd) / sum_f[k];
            if gd > 0.0 && (!(0.0..1.0).contains(&number(alpha[k])?) || number(alpha[k])? == 0.0) {
                return Err(invalid("positive Miller coupling is unrepresentable"));
            }
            let parallel = dd(gs) * alpha[k];
            if gd > 0.0 && !positive(number(parallel)?) {
                return Err(invalid("Miller Schur capacitance is unrepresentable"));
            }
            if gd != 0.0 {
                // Exact analytical parallel-capacitor composition at every
                // ORIGINAL voltage knot; no source knot is sampled or removed.
                let knots = curves[k]
                    .knots()
                    .iter()
                    .map(|&(q, v)| {
                        let augmented = number(dd(q) + parallel * dd(v))?;
                        if v != 0.0 && augmented.to_bits() == q.to_bits() {
                            return Err(invalid(
                                "authored Miller augmentation is below charge-knot resolution",
                            ));
                        }
                        Ok((augmented, v))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                curves[k] = CapacitorChargeCurve::try_new(knots)
                    .map_err(|_| invalid("original intrinsic/Miller charge composition refused"))?;
            }
        }
        let base = ChargeCommutationStorage::new(&ChargeCommutationParameters {
            converter: parameters.converter,
            output_charge_curves: curves,
            maximum_cells_per_chord: parameters.maximum_cells_per_chord,
        })?;
        let mut scales = [0.0_f64; 9];
        for k in 0..4 {
            let voltage = parameters.channel_conductance[k].gate_voltage_range_v();
            let curve = parameters.intrinsic_drain_charge_curves[k].knots();
            let drain = [curve[0].1, curve[curve.len() - 1].1];
            for gate in voltage {
                for vds in drain {
                    scales[5 + k] = scales[5 + k].max(
                        number(
                            sum_f[k] * dd(gate)
                                - dd(parameters.gate_drain_capacitance_f[k]) * dd(vds),
                        )?
                        .abs(),
                    );
                }
            }
            if !positive(scales[5 + k]) {
                return Err(invalid("gate charge support magnitude is unrepresentable"));
            }
        }
        for k in 0..2 {
            scales[1 + k] = number(
                dd(base.coordinate_scale(1 + k).unwrap()) + alpha[2 * k] * dd(scales[5 + 2 * k]),
            )?;
            scales[3 + k] = number(
                dd(base.coordinate_scale(3 + k).unwrap())
                    + alpha[2 * k] * dd(scales[5 + 2 * k])
                    + alpha[2 * k + 1] * dd(scales[6 + 2 * k]),
            )?;
        }
        Ok(Self {
            parameters,
            base,
            sum_f,
            alpha,
            scales,
        })
    }
    pub(super) fn parameters(&self) -> &GateCommutationParameters<GateChannelLaw> {
        &self.parameters
    }
    pub(super) fn maximum_cells_per_chord(&self) -> usize {
        self.parameters.maximum_cells_per_chord
    }
    fn transform(&self, x: &[Dd; 9]) -> [Dd; 5] {
        let mut y = core::array::from_fn(|k| x[k]);
        for k in 0..2 {
            let high = self.alpha[2 * k] * x[5 + 2 * k];
            y[1 + k] = y[1 + k] + high;
            y[3 + k] = y[3 + k] - high + self.alpha[2 * k + 1] * x[6 + 2 * k];
        }
        y
    }
    fn pullback(&self, e: [Dd; 5], gate: [Dd; 4]) -> [Dd; 9] {
        let mut out = [Dd::ZERO; 9];
        out[..5].copy_from_slice(&e);
        for k in 0..2 {
            out[5 + 2 * k] = self.alpha[2 * k] * (e[1 + k] - e[3 + k]) + gate[2 * k];
            out[6 + 2 * k] = self.alpha[2 * k + 1] * e[3 + k] + gate[2 * k + 1];
        }
        out
    }
    pub(super) fn state(
        &self,
        current: f64,
        bus: [f64; 2],
        node: [f64; 2],
        gate: [f64; 4],
    ) -> Result<[f64; 9], PhsError> {
        if !current.is_finite() || bus.iter().chain(&node).chain(&gate).any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite physical gate/bus/node initialization required",
            ));
        }
        let p = &self.parameters;
        let mut x = [0.0; 9];
        x[0] = number(dd(p.converter.inductance_h) * dd(current))?;
        for k in 0..2 {
            let v = [dd(bus[k]) - dd(node[k]), dd(node[k])];
            for j in 0..2 {
                p.channel_conductance[2 * k + j].evaluate(number(v[j])?, gate[2 * k + j])?;
            }
            let q = [
                curve_charge_dd(&p.intrinsic_drain_charge_curves[2 * k], v[0])?,
                curve_charge_dd(&p.intrinsic_drain_charge_curves[2 * k + 1], v[1])?,
            ];
            let gd = [
                dd(p.gate_drain_capacitance_f[2 * k]) * (v[0] - dd(gate[2 * k])),
                dd(p.gate_drain_capacitance_f[2 * k + 1]) * (v[1] - dd(gate[2 * k + 1])),
            ];
            x[1 + k] = number(dd(p.converter.bus_capacitance_f[k]) * dd(bus[k]) + q[0] + gd[0])?;
            x[3 + k] = number(q[1] - q[0] + gd[1] - gd[0])?;
            for j in 0..2 {
                x[5 + 2 * k + j] = number(
                    self.sum_f[2 * k + j] * dd(gate[2 * k + j])
                        - dd(p.gate_drain_capacitance_f[2 * k + j]) * v[j],
                )?;
            }
        }
        self.energy(&x)?;
        Ok(x)
    }
    pub(super) fn energy(&self, x: &[f64; 9]) -> Result<f64, PhsError> {
        let x = x.map(dd);
        let mut h = self.base.energy_dd(&self.transform(&x))?;
        for k in 0..4 {
            h = h + dd(0.5) * x[5 + k] * x[5 + k] / self.sum_f[k];
        }
        let value = number(h)?;
        if value < 0.0 || (value == 0.0 && x.iter().any(|v| v.hi != 0.0 || v.lo != 0.0)) {
            return Err(invalid("nonzero original gate storage is unrepresentable"));
        }
        Ok(value)
    }
    fn effort_dd(&self, x: &[Dd; 9]) -> Result<[Dd; 9], PhsError> {
        let e = self.base.efforts_dd(&self.transform(x))?;
        Ok(self.pullback(e, core::array::from_fn(|k| x[5 + k] / self.sum_f[k])))
    }
    fn chord(&self, a: &[Dd; 9], b: &[Dd; 9]) -> Result<([Dd; 9], [[Dd; 3]; 2]), PhsError> {
        let (e, j) = self.base.chord_dd(&self.transform(a), &self.transform(b))?;
        Ok((
            self.pullback(
                e,
                core::array::from_fn(|k| (dd(0.5) * a[5 + k] + dd(0.5) * b[5 + k]) / self.sum_f[k]),
            ),
            j,
        ))
    }
    pub(super) fn parts(&self, a: &[f64; 9], b: &[f64; 9]) -> Result<[f64; 15], PhsError> {
        let a = a.map(dd);
        let b = b.map(dd);
        let ea = self.effort_dd(&a)?;
        let eb = self.effort_dd(&b)?;
        let base = self
            .base
            .parts_dd(&self.transform(&a), &self.transform(&b))?;
        let drain_a = self.base.device_voltages_dd(&self.transform(&a))?;
        let drain_b = self.base.device_voltages_dd(&self.transform(&b))?;
        let mut out = [0.0; 15];
        for k in 0..3 {
            out[k] = number(base[k])?;
        }
        for k in 0..4 {
            let va = drain_a[k];
            let vb = drain_b[k];
            let curve = &self.parameters.intrinsic_drain_charge_curves[k];
            out[3 + k] = number(change_curve(
                curve,
                curve_charge_dd(curve, va)?,
                curve_charge_dd(curve, vb)?,
            )?)?;
            let ga = ea[5 + k];
            let gb = eb[5 + k];
            out[7 + k] = number(
                dd(self.parameters.gate_source_capacitance_f[k])
                    * (gb - ga)
                    * (dd(0.5) * ga + dd(0.5) * gb),
            )?;
            let ma = ga - va;
            let mb = gb - vb;
            out[11 + k] = number(
                dd(self.parameters.gate_drain_capacitance_f[k])
                    * (mb - ma)
                    * (dd(0.5) * ma + dd(0.5) * mb),
            )?;
        }
        Ok(out)
    }
    fn remainder(&self, a: &[f64; 9], b: &[f64; 9], m: &[f64; 9]) -> Result<f64, PhsError> {
        let a = a.map(dd);
        let b = b.map(dd);
        let m = m.map(dd);
        let mut r = self.base.remainder_dd(
            &self.transform(&a),
            &self.transform(&b),
            &self.transform(&m),
        )?;
        for k in 0..4 {
            r = r
                + (b[5 + k] - a[5 + k]) * dd(0.5) / self.sum_f[k]
                    * ((a[5 + k] - m[5 + k]) + (b[5 + k] - m[5 + k]));
        }
        // The public gradient is rounded once. Retain the corresponding tiny
        // constitutive remainder, independently of any accepted port work.
        let e = self.effort_dd(&m)?;
        for k in 0..9 {
            r = r + (e[k] - dd(number(e[k])?)) * (b[k] - a[k]);
        }
        number(r)
    }
    fn action(
        &self,
        a: &[f64; 9],
        b: Option<&[f64; 9]>,
        direction: &[f64; 9],
    ) -> Result<[Dd; 9], PhsError> {
        let a = a.map(dd);
        let d = direction.map(dd);
        let y = self.transform(&a);
        let dy = self.transform(&d);
        let (e, weight) = if let Some(b) = b {
            let (_, j) = self.base.chord_dd(&y, &self.transform(&b.map(dd)))?;
            let mut e = [Dd::ZERO; 5];
            e[0] = dd(0.5) * dy[0] / dd(self.parameters.converter.inductance_h);
            for k in 0..2 {
                e[1 + k] = j[k][0] * dy[1 + k] + j[k][1] * dy[3 + k];
                e[3 + k] = j[k][1] * dy[1 + k] + j[k][2] * dy[3 + k];
            }
            (e, 0.5)
        } else {
            (self.base.smooth_action_dd(&y, &dy)?, 1.0)
        };
        Ok(self.pullback(
            e,
            core::array::from_fn(|k| dd(weight) * d[5 + k] / self.sum_f[k]),
        ))
    }
}
impl Storage for GateCommutationStorage {
    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(self.clone()))
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        x.try_into()
            .ok()
            .and_then(|x| self.energy(x).ok())
            .unwrap_or(f64::NAN)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let e = x
            .try_into()
            .ok()
            .and_then(|x: &[f64; 9]| self.effort_dd(&x.map(dd)).ok());
        if out.len() == 9 {
            if let Some(e) = e {
                for k in 0..9 {
                    out[k] = e[k].to_f64();
                }
                return;
            }
        }
        out.fill(f64::NAN);
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], m: &[f64]) -> Option<f64> {
        Some(
            a.try_into()
                .ok()
                .zip(b.try_into().ok())
                .zip(m.try_into().ok())
                .and_then(|((a, b), m)| self.remainder(a, b, m).ok())
                .unwrap_or(f64::NAN),
        )
    }
    fn coordinate_scale(&self, k: usize) -> Option<f64> {
        self.scales.get(k).copied().filter(|&s| positive(s))
    }
    fn quadratic_part(&self) -> Option<crate::linear_exact::QuadraticPart> {
        let mut matrix = vec![0.0; 81];
        matrix[0] = 1.0 / self.parameters.converter.inductance_h;
        for k in 0..4 {
            matrix[(5 + k) * 9 + 5 + k] = (dd(1.0) / self.sum_f[k]).to_f64();
        }
        Some(crate::linear_exact::QuadraticPart {
            matrix,
            exact: false,
        })
    }
}
impl ConstitutiveDiscreteGradient for GateCommutationStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], _: &[f64], out: &mut [f64]) {
        let e = a
            .try_into()
            .ok()
            .zip(b.try_into().ok())
            .and_then(|(a, b): (&[f64; 9], &[f64; 9])| self.chord(&a.map(dd), &b.map(dd)).ok());
        if out.len() == 9 {
            if let Some((e, _)) = e {
                for k in 0..9 {
                    out[k] = e[k].to_f64();
                }
                return;
            }
        }
        out.fill(f64::NAN);
    }
    fn end_jacobian_action(&self, a: &[f64], b: &[f64], d: &[f64], out: &mut [f64]) -> bool {
        let result = a
            .try_into()
            .ok()
            .zip(b.try_into().ok())
            .zip(d.try_into().ok())
            .and_then(|((a, b), d)| self.action(a, Some(b), d).ok());
        if out.len() != 9 {
            return false;
        }
        let Some(e) = result else {
            return false;
        };
        for k in 0..9 {
            out[k] = e[k].to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }
}
impl ContinuousStorageTangent for GateCommutationStorage {
    fn apply(&self, a: &[f64], d: &[f64], out: &mut [f64]) -> bool {
        let result = a
            .try_into()
            .ok()
            .zip(d.try_into().ok())
            .and_then(|(a, d)| self.action(a, None, d).ok());
        if out.len() != 9 {
            return false;
        }
        let Some(e) = result else {
            return false;
        };
        for k in 0..9 {
            out[k] = e[k].to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }
}
