//! A linear conductance between two heat traces. Ports `a` and `b` are field heat ports: they
//! receive the temperatures at `n_a` and `n_b` points and return the heat load [W] into the
//! element at each point.
//!
//! The interface is a list of points ([`InterfacePoint`]). At each one the temperature of either
//! trace is interpolated from up to three of its points, with weights that are not negative and
//! sum to one, and the heat `G (T_a - T_b)` crosses there from `a` to `b` through the point's
//! conductance `G` [W/K]. Each side's share of that heat goes back to its points by the same
//! weights. The loads are therefore `K T` with `K = D^T W D`: one row of `D` per interface point,
//! `a`'s weights and minus `b`'s, and `W` the diagonal of the conductances.
//!
//! A contact between two P1 traces is the quadrature of `(T_a - T_b) (v_a - v_b) / R''` over
//! their common refinement: one interface point per quadrature point, its weights the two
//! traces' basis values there, its conductance the quadrature weight over `R''`. With one point
//! on each side and one interface point it is [`super::Conductance`].

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use super::{ElementError, temperature};

/// One point of the interface between two traces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfacePoint {
    /// Up to three points of trace `a` and the weights that interpolate its temperature here.
    /// A slot that is not used has the weight zero.
    pub a: [(usize, f64); 3],
    /// The same for trace `b`.
    pub b: [(usize, f64); 3],
    /// The conductance of the interface at this point [W/K].
    pub conductance_w_k: f64,
}

/// The most a side's weights at one interface point may differ from one in their sum.
const UNIT_SUM: f64 = 8.0 * f64::EPSILON;

/// A conductance between the field heat ports `a` and `b` through a list of interface points.
///
/// It stores nothing and books no dissipation: what enters at the hot points leaves at the cold
/// ones. At every interface point heat crosses from the hotter interpolated temperature to the
/// colder, so the entropy it generates, `dt sum_q G_q (T_a - T_b)_q^2 / (T_a T_b)_q`, is never
/// negative.
///
/// The kernel books the entropy through a field heat port point by point, `dt sum Q_i / T_i`.
/// That nodal sum over both ports is not minus what the element generates: the two differ by the
/// difference between a trace's temperature at an interface point and at the points it is
/// interpolated from. The element's own entropy balance closes only to that difference.
#[derive(Debug, Clone)]
pub struct ConductanceMatrix {
    widths: [usize; 2],
    points: Vec<InterfacePoint>,
    /// The most interface points that load one point of a trace.
    depth: usize,
    ports: [PortSpec; 2],
}

impl ConductanceMatrix {
    /// A conductance between a trace of `width_a` points at port `a` and one of `width_b` points
    /// at port `b`, through `points`.
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when a width is zero, there is no interface point, a weight's
    /// point is outside its trace, a weight is negative or not finite, a side's weights do not
    /// sum to one within eight roundoff units, or a conductance is not positive and finite.
    pub fn new(
        width_a: usize,
        width_b: usize,
        points: Vec<InterfacePoint>,
    ) -> Result<Self, ElementError> {
        let ports = [("a", width_a), ("b", width_b)].map(|(name, width)| {
            NonZeroUsize::new(width).map(|width| {
                let id = StableId::new(name).expect("a valid stable id");
                PortSpec::heat_field(id, PortRole::EffortIn, width).expect("a valid heat port")
            })
        });
        let [Some(port_a), Some(port_b)] = ports else {
            return Err(ElementError::Parameter {
                name: "port width (points, at least 1)",
                value: width_a.min(width_b) as f64,
            });
        };
        if points.is_empty() {
            return Err(ElementError::Parameter {
                name: "interface points (at least 1)",
                value: 0.0,
            });
        }
        let mut loaded = vec![0_usize; width_a + width_b];
        for point in &points {
            for (side, width, offset) in [(&point.a, width_a, 0), (&point.b, width_b, width_a)] {
                let mut sum = 0.0;
                for &(index, weight) in side {
                    if index >= width {
                        return Err(ElementError::Parameter {
                            name: "interface weight's point (inside its trace)",
                            value: index as f64,
                        });
                    }
                    if !(weight.is_finite() && weight >= 0.0) {
                        return Err(ElementError::Parameter {
                            name: "interface weight (not negative)",
                            value: weight,
                        });
                    }
                    if weight > 0.0 {
                        loaded[offset + index] += 1;
                    }
                    sum += weight;
                }
                if !((sum - 1.0).abs() <= UNIT_SUM) {
                    return Err(ElementError::Parameter {
                        name: "interface weights' sum (one within 8 eps)",
                        value: sum,
                    });
                }
            }
            if !(point.conductance_w_k.is_finite() && point.conductance_w_k > 0.0) {
                return Err(ElementError::Parameter {
                    name: "interface conductance [W/K]",
                    value: point.conductance_w_k,
                });
            }
        }
        Ok(Self {
            widths: [width_a, width_b],
            depth: loaded.into_iter().max().unwrap_or(0),
            points,
            ports: [port_a, port_b],
        })
    }

    /// The points of ports `a` and `b`.
    #[must_use]
    pub fn widths(&self) -> [usize; 2] {
        self.widths
    }

    /// The interface points.
    #[must_use]
    pub fn points(&self) -> &[InterfacePoint] {
        &self.points
    }

    /// `K`, row-major over the `n_a + n_b` points, `a`'s first. Dense: for tests and small traces.
    #[must_use]
    pub fn matrix(&self) -> Vec<f64> {
        let n = self.widths[0] + self.widths[1];
        let mut column = vec![0.0; n];
        let mut matrix = vec![0.0; n * n];
        let mut unit = vec![0.0; n];
        for j in 0..n {
            unit[j] = 1.0;
            self.apply(&unit, &mut column);
            unit[j] = 0.0;
            for (i, value) in column.iter().enumerate() {
                matrix[i * n + j] = *value;
            }
        }
        matrix
    }

    /// The values `values` take at interface point `point` on trace `a` and on trace `b`.
    fn at(&self, point: &InterfacePoint, values: &[f64]) -> (f64, f64) {
        let mut sides = [0.0_f64; 2];
        for (side, weights, offset) in [(0, &point.a, 0), (1, &point.b, self.widths[0])] {
            for &(index, weight) in weights {
                sides[side] = weight.mul_add(values[offset + index], sides[side]);
            }
        }
        (sides[0], sides[1])
    }

    /// The heat rate `rate` that crosses at `point` from `a` to `b`, added to the loads `out`.
    fn spread(&self, point: &InterfacePoint, rate: f64, out: &mut [f64]) {
        for &(index, weight) in &point.a {
            out[index] = weight.mul_add(rate, out[index]);
        }
        for &(index, weight) in &point.b {
            let load = &mut out[self.widths[0] + index];
            *load = (-weight).mul_add(rate, *load);
        }
    }

    /// `K x` into `out`.
    fn apply(&self, values: &[f64], out: &mut [f64]) {
        let n = self.widths[0] + self.widths[1];
        out[..n].fill(0.0);
        for point in &self.points {
            let (a, b) = self.at(point, values);
            self.spread(point, point.conductance_w_k * (a - b), out);
        }
    }
}

impl Component for ConductanceMatrix {
    type State = ();
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::ConductanceMatrix",
        state_family: "stateless",
        schema_grammar: "(): a memoryless conductance between two traces stores nothing",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        // The widths are part of the system's port description; the interface points are not.
        words.push(self.points.len() as u64);
        for point in &self.points {
            for &(index, weight) in point.a.iter().chain(&point.b) {
                words.push(index as u64);
                words.push(weight.to_bits());
            }
            words.push(point.conductance_w_k.to_bits());
        }
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        let n = self.widths[0] + self.widths[1];
        for &kelvin in &inputs[..n] {
            temperature(kelvin)?;
        }
        outputs[..n].fill(0.0);
        let (mut generated, mut crossing) = (0.0_f64, 0.0_f64);
        for point in &self.points {
            // Both are weighted means of temperatures above zero, by weights that are not
            // negative and sum to one: above zero themselves.
            let (hot, cold) = self.at(point, inputs);
            let rate_w = point.conductance_w_k * (hot - cold);
            self.spread(point, rate_w, outputs);
            generated += dt_s * point.conductance_w_k * (hot - cold) * (hot - cold) / (hot * cold);
            crossing += rate_w.abs();
        }
        // The loads sum to `sum_q rate_q (s_a - s_b)_q`, with `s` the sums of a point's weights,
        // each within 8 eps of one. A load is a sum of at most `depth` products and a port's
        // energy a sum of its `n` loads, each sum within its length in roundoff units of the
        // sizes summed (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., sec.
        // 3.1); the two products with dt round once each.
        let units = (2 * self.depth + n + 18) as f64 * f64::EPSILON;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: generated,
            energy_tolerance_j: units * dt_s * crossing,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for ConductanceMatrix {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.widths[0] + self.widths[1];
        self.apply(&d_inputs[..n], d_outputs);
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // K = D^T W D with W diagonal, so K^T bar is K bar.
        let n = self.widths[0] + self.widths[1];
        self.apply(&bar_outputs[..n], bar_inputs);
        Ok(())
    }
}
