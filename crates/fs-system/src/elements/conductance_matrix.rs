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
//!
//! The conductances may follow the temperatures. With a [`ConductanceLaw`]
//! ([`ConductanceMatrix::following`]) every interface point's conductance is its constant `G`
//! times the law's factor at the two interpolated temperatures there, evaluated in every
//! proposal. The loads are then no longer linear in the temperatures; the tangent carries the
//! factor's two partial derivatives, and its transpose is exact.

use core::fmt;
use core::num::NonZeroUsize;
use std::sync::Arc;

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

/// A place where an element's law is asked and no heat crosses: the two interpolation rows of
/// an interface point, without a conductance ([`ConductanceMatrix::guarded_at`]).
///
/// A law is asked at the interface points, at temperatures interpolated from a trace's points.
/// A trace point outside the law that reaches every interface point through a mean inside it
/// would go unseen. A guard row with the weight one on that point asks the law there as well.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuardRow {
    /// Up to three points of trace `a` and the weights that interpolate its temperature here.
    pub a: [(usize, f64); 3],
    /// The same for trace `b`.
    pub b: [(usize, f64); 3],
}

/// The most a side's weights at one interface point may differ from one in their sum.
const UNIT_SUM: f64 = 8.0 * f64::EPSILON;

/// One side's `(point, weight)` pairs at an interface point or a guard row of a trace of
/// `width` points: every point inside the trace, every weight finite and not negative, and
/// the weights' sum one within [`UNIT_SUM`].
fn admitted(side: &[(usize, f64); 3], width: usize) -> Result<(), ElementError> {
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
        sum += weight;
    }
    if (sum - 1.0).abs() <= UNIT_SUM {
        Ok(())
    } else {
        Err(ElementError::Parameter {
            name: "interface weights' sum (one within 8 eps)",
            value: sum,
        })
    }
}

/// How an interface's conductance follows the temperatures of its two traces.
///
/// The element multiplies every interface point's conductance by [`ConductanceLaw::factor`] at
/// the two temperatures interpolated there. A law that returns one, with no slope, is the
/// constant element.
pub trait ConductanceLaw: fmt::Debug + Send + Sync {
    /// The factor on a point's conductance at the temperatures `kelvin_a` and `kelvin_b` of the
    /// two traces there, and its partial derivatives with respect to them [1/K].
    ///
    /// The factor is positive and finite and the partials are finite. Where a piecewise law has
    /// a knot, the partials are those of the inward segment, the one the state would enter by a
    /// step towards the law's interior: the convention of fs-conduction's Newton linearization
    /// of its temperature contacts (`linearize_temperature_for_newton`).
    ///
    /// # Errors
    /// Why the law has no answer at these temperatures, for example a state outside the
    /// validity of the card behind it. The element then refuses the step with
    /// [`ElementError::Law`]; it never clamps the state or extrapolates the law.
    fn factor(&self, kelvin_a: f64, kelvin_b: f64) -> Result<(f64, [f64; 2]), String>;

    /// Optional direct partials of `factor(T_a, T_b) * (T_a - T_b)`,
    /// with respect to the two temperatures. These partials are dimensionless;
    /// the element multiplies them by each point's reference conductance [W/K].
    ///
    /// The default uses the product rule from [`Self::factor`]. A law may
    /// supply an equivalent form to avoid cancellation, for example the cold
    /// surface derivative of radiation at a large temperature contrast.
    /// The factor's domain and positivity checks still run first, and direct
    /// partials must be finite. This changes neither heat nor its entropy.
    fn rate_partials(&self, _kelvin_a: f64, _kelvin_b: f64) -> Option<[f64; 2]> {
        None
    }

    /// Words that identify the law. They join the element's identity after its interface
    /// points, so two elements that differ in their law differ in identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// Surface-only diffuse-gray exchange. The reference temperature is 1 K, so
/// `reference_conductance_w_k = 4 sigma A epsilon_eff (1 K)^3` and the factor
/// is dimensionless. The rounded SI constant and effective-emissivity reduction
/// match fs-thermal's `enthalpy::conduction::implicit::link::GrayRadiativeLink`.
#[derive(Debug)]
struct GrayRadiation {
    area_m2: f64,
    emissivities: [f64; 2],
    reference_conductance_w_k: f64,
}

impl ConductanceLaw for GrayRadiation {
    fn factor(&self, kelvin_a: f64, kelvin_b: f64) -> Result<(f64, [f64; 2]), String> {
        let a = temperature(kelvin_a).map_err(|error| error.to_string())?;
        let b = temperature(kelvin_b).map_err(|error| error.to_string())?;
        // With a and b measured in units of 1 K, factor(1, 1) = 1.
        // Factoring a^4 - b^4 preserves small temperature differences.
        let sum = a + b;
        let squares = a * a + b * b;
        let factor = (0.25 * sum) * squares;
        let partials = [
            0.25 * squares + (0.5 * a) * sum,
            0.25 * squares + (0.5 * b) * sum,
        ];
        let conductance = self.reference_conductance_w_k * factor;
        if !(factor.is_finite()
            && factor > 0.0
            && partials.iter().all(|value| value.is_finite())
            && conductance.is_finite()
            && conductance > 0.0
            && (conductance * (a - b)).is_finite())
        {
            return Err("gray radiation has an unrepresentable conductance or heat rate".into());
        }
        Ok((factor, partials))
    }

    fn rate_partials(&self, kelvin_a: f64, kelvin_b: f64) -> Option<[f64; 2]> {
        // The normalized flux is (Ta^4 - Tb^4) / (4 (1 K)^3).
        // Its derivative at the colder side must not subtract large hot-side terms.
        Some([
            (kelvin_a * kelvin_a) * kelvin_a,
            -(kelvin_b * kelvin_b) * kelvin_b,
        ])
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            u64::from_le_bytes(*b"grayrad1"),
            self.area_m2.to_bits(),
            self.emissivities[0].to_bits(),
            self.emissivities[1].to_bits(),
        ]);
    }
}

/// `law`'s factor and partials at one interface point, checked, or the step's refusal.
fn checked_factor(
    law: &dyn ConductanceLaw,
    kelvin_a: f64,
    kelvin_b: f64,
) -> Result<(f64, [f64; 2]), ElementError> {
    let (factor, partials) = law
        .factor(kelvin_a, kelvin_b)
        .map_err(|reason| ElementError::Law { reason })?;
    if factor.is_finite() && factor > 0.0 && partials.iter().all(|partial| partial.is_finite()) {
        Ok((factor, partials))
    } else {
        Err(ElementError::Law {
            reason: format!(
                "the factor {factor} with partials {partials:?} at {kelvin_a} K and {kelvin_b} K \
                 is not positive and finite"
            ),
        })
    }
}

/// Partials of the normalized rate, preserving the factor's domain checks.
fn checked_rate_partials(
    law: &dyn ConductanceLaw,
    kelvin_a: f64,
    kelvin_b: f64,
) -> Result<[f64; 2], ElementError> {
    let (factor, [by_a, by_b]) = checked_factor(law, kelvin_a, kelvin_b)?;
    let partials = law.rate_partials(kelvin_a, kelvin_b).unwrap_or_else(|| {
        let jump = kelvin_a - kelvin_b;
        [factor + jump * by_a, jump * by_b - factor]
    });
    if partials.iter().all(|value| value.is_finite()) {
        Ok(partials)
    } else {
        Err(ElementError::Law {
            reason: "the normalized heat-rate partials are not finite".into(),
        })
    }
}

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
    /// What scales the points' conductances with the temperatures there; none for constant ones.
    law: Option<Arc<dyn ConductanceLaw>>,
    /// Where the law is asked besides the interface points. Empty without a law.
    guards: Vec<GuardRow>,
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
                admitted(side, width)?;
                for &(index, weight) in side {
                    if weight > 0.0 {
                        loaded[offset + index] += 1;
                    }
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
            law: None,
            guards: Vec::new(),
        })
    }

    /// Radiation between two opaque diffuse-gray surfaces of equal area [m2]
    /// and full mutual view, at the absolute temperatures [K] of ports `a`
    /// and `b`. Both heat ports have width one. The fixed emissivities are
    /// ordered by port, and heat into `a` is
    /// `sigma A epsilon_eff (T_a^4 - T_b^4)`, with
    /// `epsilon_eff = 1 / (1/epsilon_a + 1/epsilon_b - 1)`.
    ///
    /// This reuses the conductance component's conservative heat exchange,
    /// entropy generation and analytic tangent. Its stored reference
    /// conductance is `4 sigma A epsilon_eff (1 K)^3` [W/K]; the dimensionless
    /// law factor is `(T_a + T_b)(T_a^2 + T_b^2) / (4 (1 K)^3)`.
    /// Thus neither near-equilibrium heat nor its tangent subtracts fourth
    /// powers. The state is `()`, and the class is [`Class::Smooth`].
    ///
    /// These are surface temperatures: no half-cell or coating resistance,
    /// surface storage, participating medium, variable emissivity, view-factor
    /// geometry or enclosure radiosity solve is supplied. The original area
    /// and emissivity bits participate in durable problem identity.
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when area is not positive and finite, an
    /// emissivity is outside `(0, 1]`, or the reference conductance is not
    /// representable as a positive finite number. Perfect reflectors are
    /// outside this positive-conductance constructor. Unrepresentable live
    /// conductance or heat rates refuse through [`ElementError::Law`].
    pub fn gray_radiation(area_m2: f64, emissivities: [f64; 2]) -> Result<Self, ElementError> {
        if !(area_m2.is_finite() && area_m2 > 0.0) {
            return Err(ElementError::Parameter {
                name: "radiating area [m2]",
                value: area_m2,
            });
        }
        for (name, value) in ["emissivity_a", "emissivity_b"]
            .into_iter()
            .zip(emissivities)
        {
            if !(value.is_finite() && value > 0.0 && value <= 1.0) {
                return Err(ElementError::Parameter { name, value });
            }
        }
        // Same stable reduction as the native finite-cell radiation owner:
        // never multiply two tiny emissivities or take their reciprocals.
        let low = emissivities[0].min(emissivities[1]);
        let high = emissivities[0].max(emissivities[1]);
        let effective = low / ((1.0 - low) + low / high);
        let reference_conductance_w_k = (4.0 * 5.670_374_419e-8 * area_m2) * effective;
        let weights = [(0, 1.0), (0, 0.0), (0, 0.0)];
        Ok(Self::new(
            1,
            1,
            vec![InterfacePoint {
                a: weights,
                b: weights,
                conductance_w_k: reference_conductance_w_k,
            }],
        )?
        .following(Arc::new(GrayRadiation {
            area_m2,
            emissivities,
            reference_conductance_w_k,
        })))
    }

    /// The same interface with every point's conductance scaled by `law`'s factor at the two
    /// temperatures there, in every proposal.
    ///
    /// The element is then [`Class::Smooth`], its tangent carries the law's partial derivatives,
    /// and a step at which the law has no answer is refused with [`ElementError::Law`]. The
    /// entropy rule is unchanged: the scaled conductance is positive.
    #[must_use]
    pub fn following(mut self, law: Arc<dyn ConductanceLaw>) -> Self {
        self.law = Some(law);
        self
    }

    /// The same element with its law asked at `rows` too, in every proposal, before the
    /// interface points. No heat crosses at a guard row and it is no part of the tangent: it
    /// only refuses, with [`ElementError::Law`], a step whose temperatures there the law has no
    /// answer for. A contact between matching traces is guarded at each pair of trace points
    /// that face each other, as fs-conduction guards its pointwise contact on one mesh.
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when the element follows no law, since a row would then
    /// guard nothing, or when a row's weights are not an interface point's: a point outside
    /// its trace, a weight that is negative or not finite, or a side's weights that do not sum
    /// to one within eight roundoff units.
    pub fn guarded_at(mut self, rows: Vec<GuardRow>) -> Result<Self, ElementError> {
        if self.law.is_none() {
            return Err(ElementError::Parameter {
                name: "guard rows (of an element that follows a law)",
                value: rows.len() as f64,
            });
        }
        for row in &rows {
            admitted(&row.a, self.widths[0])?;
            admitted(&row.b, self.widths[1])?;
        }
        self.guards = rows;
        Ok(self)
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

    /// `K`, row-major over the `n_a + n_b` points, `a`'s first, at the points' constant
    /// conductances (a law's factor of one). Dense: for tests and small traces.
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
        // A constant element's words end here, as they did before a law could be set. Under a
        // law its words follow, then the number of guard rows and each row's points and weights.
        if let Some(law) = &self.law {
            law.identity_words(words);
            words.push(self.guards.len() as u64);
            for row in &self.guards {
                for &(index, weight) in row.a.iter().chain(&row.b) {
                    words.push(index as u64);
                    words.push(weight.to_bits());
                }
            }
        }
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        // Under a law the slopes move with the temperatures.
        if self.law.is_some() {
            Class::Smooth
        } else {
            Class::ExactLinear
        }
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
        // A law is asked at its guard rows first. No heat crosses there; a refusal there ends
        // the proposal before any load is formed.
        if let Some(law) = &self.law {
            for row in &self.guards {
                let here = InterfacePoint {
                    a: row.a,
                    b: row.b,
                    conductance_w_k: 0.0,
                };
                let (kelvin_a, kelvin_b) = self.at(&here, inputs);
                checked_factor(law.as_ref(), kelvin_a, kelvin_b)?;
            }
        }
        outputs[..n].fill(0.0);
        let (mut generated, mut crossing) = (0.0_f64, 0.0_f64);
        for point in &self.points {
            // Both are weighted means of temperatures above zero, by weights that are not
            // negative and sum to one: above zero themselves.
            let (hot, cold) = self.at(point, inputs);
            // A law's factor is positive, so the scaled conductance is too: heat still crosses
            // from the hotter side and the entropy below is still not negative.
            let conductance_w_k = match &self.law {
                None => point.conductance_w_k,
                Some(law) => point.conductance_w_k * checked_factor(law.as_ref(), hot, cold)?.0,
            };
            let rate_w = conductance_w_k * (hot - cold);
            self.spread(point, rate_w, outputs);
            generated += dt_s * conductance_w_k * (hot - cold) * (hot - cold) / (hot * cold);
            crossing += rate_w.abs();
        }
        // The loads sum to `sum_q rate_q (s_a - s_b)_q`, with `s` the sums of a point's weights,
        // each within 8 eps of one. A load is a sum of at most `depth` products and a port's
        // energy a sum of its `n` loads, each sum within its length in roundoff units of the
        // sizes summed (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., sec.
        // 3.1); the two products with dt round once each. A law changes a point's rate and
        // nothing in this count: both sides are handed the same rate.
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
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.widths[0] + self.widths[1];
        let Some(law) = &self.law else {
            self.apply(&d_inputs[..n], d_outputs);
            return Ok(());
        };
        // The rate at a point is `g f (T_a - T_b)`, so its differential is the block
        // `g [f + (T_a - T_b) f_a, -f + (T_a - T_b) f_b]` on `(dT_a, dT_b)`.
        d_outputs[..n].fill(0.0);
        for point in &self.points {
            let (hot, cold) = self.at(point, inputs);
            let (d_hot, d_cold) = self.at(point, d_inputs);
            let [by_a, by_b] = checked_rate_partials(law.as_ref(), hot, cold)?;
            let d_rate = point.conductance_w_k * (by_a * d_hot + by_b * d_cold);
            self.spread(point, d_rate, d_outputs);
        }
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.widths[0] + self.widths[1];
        let Some(law) = &self.law else {
            // K = D^T W D with W diagonal, so K^T bar is K bar.
            self.apply(&bar_outputs[..n], bar_inputs);
            return Ok(());
        };
        // The transpose of the block in `jvp`, point by point. The loads are the rate times
        // a's weights and minus the rate times b's, so the rate's covector is the interpolated
        // covector of a's loads minus that of b's; it goes back through the two entries of the
        // block, each along its own side's weights. Under a law the block is not symmetric.
        bar_inputs[..n].fill(0.0);
        for point in &self.points {
            let (hot, cold) = self.at(point, inputs);
            let (bar_a, bar_b) = self.at(point, bar_outputs);
            let [by_a, by_b] = checked_rate_partials(law.as_ref(), hot, cold)?;
            let bar_rate = point.conductance_w_k * (bar_a - bar_b);
            let (to_a, to_b) = (bar_rate * by_a, bar_rate * by_b);
            for &(index, weight) in &point.a {
                bar_inputs[index] = weight.mul_add(to_a, bar_inputs[index]);
            }
            for &(index, weight) in &point.b {
                let bar = &mut bar_inputs[self.widths[0] + index];
                *bar = weight.mul_add(to_b, *bar);
            }
        }
        Ok(())
    }
}
