//! Finite reversible exchange between two original narrow-beam snapshots.
//!
//! Mechanical work uses the endpoint-average generalized stresses. This is a
//! symmetric, refinement-consistent work quadrature, not the exact integral
//! along an arbitrary varying-temperature mechanical path. Heat into the
//! original caloric inventory is `W - delta U`; it is signed reversible heat,
//! never a conversion of an energy residual into dissipation.

use super::{LaminateBeamThermoelasticError, ResolvedLaminateBeamThermoelastic, invalid};
use fs_math::dd::Dd;

/// One original layer's finite exchange, per ORIGINAL beam length.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LaminateBeamThermoelasticExchange {
    /// Endpoint-average N/M contracted with the change in eps0/kappa [J/m].
    pub mechanical_work_j_m: f64,
    pub free_energy_change_j_m: f64,
    pub entropy_change_j_m_k: f64,
    /// `delta Psi + delta(T*S)` using this layer's own endpoint temperatures.
    pub internal_energy_change_j_m: f64,
    /// Positive into the original caloric inventory: `W - delta U` [J/m].
    pub reversible_heat_j_m: f64,
    /// Rounded endpoint `Psi1-Psi0` minus the polarized increment [J/m].
    /// Diagnostic only: this subtraction error is never assigned to heat.
    pub free_energy_balance_residual_j_m: f64,
}

impl ResolvedLaminateBeamThermoelastic {
    /// Exchange one original layer between two compatible source snapshots.
    ///
    /// Both endpoints retain the same original source owner, reference section,
    /// width, axis, affine-source interval and full strain envelope. Current
    /// temperatures may differ. Each endpoint's original potential is evaluated
    /// exactly once per Gauss point after checking both layer faces.
    ///
    /// Changes in the original quadratic potential and its entropy are
    /// polarized before summation, avoiding subtraction of large stored energies
    /// for small mechanical increments. Double-double contractions retain the
    /// cancellation in `W-delta U`; they do not increase source-data precision.
    /// No allocation, source query, new material law or mechanical time step is
    /// introduced. The consumer owns length quadrature and its finite-time heat
    /// balance. A fixed-temperature elastic increment can exchange heat through
    /// `-T*delta S` even though its mechanical work equals `delta Psi`.
    pub fn layer_exchange(
        &self,
        end: &Self,
        layer: usize,
        before: [f64; 2],
        after: [f64; 2],
    ) -> Result<LaminateBeamThermoelasticExchange, LaminateBeamThermoelasticError> {
        let start_layer = self.layer(layer)?;
        let end_layer = end.layer(layer)?;
        if !core::ptr::eq(
            self.current.source().layers(),
            end.current.source().layers(),
        ) || !core::ptr::eq(
            self.current.source().reference_section(),
            end.current.source().reference_section(),
        ) || self.layers.len() != end.layers.len()
            || self.beam.width_m.to_bits() != end.beam.width_m.to_bits()
            || self.beam.axis_offset_m.to_bits() != end.beam.axis_offset_m.to_bits()
            || start_layer.temperature_range.map(f64::to_bits)
                != end_layer.temperature_range.map(f64::to_bits)
            || start_layer.limit.to_bits() != end_layer.limit.to_bits()
        {
            return Err(invalid(
                "beam exchange requires the same original source, reference geometry, interval and strain envelope",
            ));
        }
        let start_points = self.layer_points(layer, before[0], before[1])?;
        let end_points = end.layer_points(layer, after[0], after[1])?;
        let q0 = start_layer.resolved.stiffness_derivatives_pa();
        let q1 = end_layer.resolved.stiffness_derivatives_pa();
        let alpha0 = start_layer.resolved.free_strain_temperature_derivatives()[0].map(dd);
        let alpha1 = end_layer.resolved.free_strain_temperature_derivatives()[0].map(dd);
        let delta_alpha = difference(alpha1, alpha0);
        let average_alpha = average(alpha1, alpha0);
        let t0 = dd(start_layer.resolved.temperature_k());
        let t1 = dd(end_layer.resolved.temperature_k());
        let delta_t = t1 - t0;
        let average_t = (t0 + t1) * dd(0.5);
        let delta_eps = dd(after[0]) - dd(before[0]);
        let delta_kappa = dd(after[1]) - dd(before[1]);
        let mut work = Dd::ZERO;
        let mut free = Dd::ZERO;
        let mut entropy = Dd::ZERO;
        let mut internal = Dd::ZERO;
        let mut rounded_free = Dd::ZERO;
        for (point, (start, finish)) in start_points.into_iter().zip(end_points).enumerate() {
            let e0 = core::array::from_fn(|i| dd(start.state.elastic_strain_engineering[i]));
            let e1 = core::array::from_fn(|i| dd(finish.state.elastic_strain_engineering[i]));
            // Preserve the authored increment before total-minus-free strain
            // rounds it away. In particular, an isothermal one-ULP eps change
            // must still exchange -T*delta S when the thermal strain is large.
            let delta_x = delta_eps
                - dd(start.y) * delta_kappa
                - (dd(end.beam.source.samples[layer].free[point][0])
                    - dd(self.beam.source.samples[layer].free[point][0]));
            let average_x = (dd(start.elastic_x) + dd(finish.elastic_x)) * dd(0.5);
            let r0 = self.beam.layers[layer].relaxation.map(dd);
            let r1 = end.beam.layers[layer].relaxation.map(dd);
            let relaxed_delta =
                |i: usize| -(r1[i] + r0[i]) * dd(0.5) * delta_x - (r1[i] - r0[i]) * average_x;
            let de = principal_increment(
                [delta_x, relaxed_delta(0), relaxed_delta(1)],
                start_layer.sin,
                start_layer.cos,
            );
            let e_mean = average(e1, e0);
            let sigma0 = apply(q0[0], e0);
            let sigma1 = apply(q1[0], e1);
            let sigma_mean = average(sigma1, sigma0);
            let delta_sigma = sum(
                apply_average(q1[0], q0[0], de),
                apply_difference(q1[0], q0[0], e_mean),
            );
            let volume = dd(start.volume);
            // Polarization of the SAME original in-plane quadratic potential.
            let dpsi = volume
                * (dot(de, sigma_mean) + dd(0.5) * dot(e0, apply_difference(q1[0], q0[0], e1)));
            // S/V = alpha dot(Q e) - e dot(Q_T e)/2. Differentiate no new law;
            // all matrices and expansion derivatives belong to the source owner.
            let ds = volume
                * (dot(delta_alpha, sigma_mean) + dot(average_alpha, delta_sigma)
                    - dd(0.5)
                        * (dot(de, sum(apply(q1[1], e1), apply(q0[1], e0)))
                            + dot(e0, apply_difference(q1[1], q0[1], e1))));
            let s0 = volume * (dot(alpha0, sigma0) - dd(0.5) * dot(e0, apply(q0[1], e0)));
            let s1 = volume * (dot(alpha1, sigma1) - dd(0.5) * dot(e1, apply(q1[1], e1)));
            let du = dpsi + average_t * ds + (s0 + s1) * dd(0.5) * delta_t;
            let sigma_x0 = dd(self.beam.layers[layer].axial_pa) * dd(start.elastic_x);
            let sigma_x1 = dd(end.beam.layers[layer].axial_pa) * dd(finish.elastic_x);
            work = work
                + volume
                    * dd(0.5)
                    * (sigma_x0 + sigma_x1)
                    * (delta_eps - dd(start.y) * delta_kappa);
            free = free + dpsi;
            entropy = entropy + ds;
            internal = internal + du;
            rounded_free =
                rounded_free + dd(finish.state.free_energy_j) - dd(start.state.free_energy_j);
        }
        let result = LaminateBeamThermoelasticExchange {
            mechanical_work_j_m: work.to_f64(),
            free_energy_change_j_m: free.to_f64(),
            entropy_change_j_m_k: entropy.to_f64(),
            internal_energy_change_j_m: internal.to_f64(),
            reversible_heat_j_m: (work - internal).to_f64(),
            free_energy_balance_residual_j_m: (rounded_free - free).to_f64(),
        };
        if [
            result.mechanical_work_j_m,
            result.free_energy_change_j_m,
            result.entropy_change_j_m_k,
            result.internal_energy_change_j_m,
            result.reversible_heat_j_m,
            result.free_energy_balance_residual_j_m,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            Ok(result)
        } else {
            Err(invalid(
                "beam finite exchange arithmetic is unrepresentable",
            ))
        }
    }
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn sum(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| a[i] + b[i])
}
fn difference(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn average(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    sum(a, b).map(|v| v * dd(0.5))
}
fn dot(a: [Dd; 3], b: [Dd; 3]) -> Dd {
    (0..3).fold(Dd::ZERO, |v, i| v + a[i] * b[i])
}
fn principal_increment(e: [Dd; 3], s: f64, c: f64) -> [Dd; 3] {
    let (ss, cc, sc) = (dd(s * s), dd(c * c), dd(s * c));
    [
        cc * e[0] + ss * e[1] + sc * e[2],
        ss * e[0] + cc * e[1] - sc * e[2],
        dd(2.0) * sc * (e[1] - e[0]) + (cc - ss) * e[2],
    ]
}
fn apply(q: [f64; 9], x: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| (0..3).fold(Dd::ZERO, |v, j| v + dd(q[3 * i + j]) * x[j]))
}
fn apply_difference(q1: [f64; 9], q0: [f64; 9], x: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| {
        (0..3).fold(Dd::ZERO, |v, j| {
            v + (dd(q1[3 * i + j]) - dd(q0[3 * i + j])) * x[j]
        })
    })
}
fn apply_average(q1: [f64; 9], q0: [f64; 9], x: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| {
        (0..3).fold(Dd::ZERO, |v, j| {
            v + (dd(q1[3 * i + j]) + dd(q0[3 * i + j])) * dd(0.5) * x[j]
        })
    })
}
