//! One simultaneous control-volume path, including its separate heat and work.

use super::super::{inflow, outflow};
use fs_math::dd::Dd;
use inflow::IncomingMixture;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct ThroughflowCell {
    pub(crate) streams: [Option<IncomingMixture>; 2],
    pub(crate) mass_kg: [f64; 2],
}

impl ThroughflowCell {
    pub(crate) fn add(&mut self, side: usize, stream: IncomingMixture) -> Option<()> {
        if let Some(existing) = &mut self.streams[side] {
            existing.combine(stream)
        } else {
            self.streams[side] = Some(stream);
            Some(())
        }
    }

    pub(crate) fn requested(&self) -> [f64; 2] {
        self.streams.map(|stream| stream.map_or(0.0, |s| s.mass_kg))
    }

    /// The smaller directional batch is unchanged. The larger batch accepts
    /// the native rounding of their NET mass; no division by a tiny net flux.
    pub(crate) fn represent(&mut self, dm: f64) -> Option<()> {
        let requested = self.requested();
        let net = requested[0] - requested[1];
        let mut represented = requested;
        if net > 0.0 {
            represented[0] = (Dd::from_f64(requested[1]) + Dd::from_f64(dm)).to_f64();
        } else if net < 0.0 {
            represented[1] = (Dd::from_f64(requested[0]) - Dd::from_f64(dm)).to_f64();
        }
        if !dm.is_finite()
            || (net == 0.0 && dm != 0.0)
            || (net != 0.0 && (dm == 0.0 || net.is_sign_positive() != dm.is_sign_positive()))
            || represented.iter().enumerate().any(|(side, m)| {
                !m.is_finite() || *m < 0.0 || (self.streams[side].is_some() && *m == 0.0)
            })
        {
            return None;
        }
        self.mass_kg = represented;
        Some(())
    }

    pub(crate) fn mass_defect(&self, dm: f64) -> f64 {
        let requested = self.requested();
        (self.mass_kg[0] - requested[0]).abs()
            + (self.mass_kg[1] - requested[1]).abs()
            + (Dd::from_f64(self.mass_kg[0]) - Dd::from_f64(self.mass_kg[1]) - Dd::from_f64(dm))
                .to_f64()
                .abs()
    }

    pub(crate) fn heat(&self) -> f64 {
        self.one_way().map_or_else(
            || self.streams.iter().flatten().map(|s| s.heat_j).sum(),
            |stream| stream.heat_j,
        )
    }

    fn one_way(&self) -> Option<IncomingMixture> {
        self.streams[0].xor(self.streams[1])
    }

    fn net(&self, value: impl Fn(IncomingMixture) -> Dd) -> Dd {
        let term = |side: usize| {
            self.streams[side].map_or(Dd::ZERO, |s| Dd::from_f64(self.mass_kg[side]) * value(s))
        };
        term(0) - term(1)
    }

    pub(crate) fn carried_enthalpy(&self, dm: f64) -> f64 {
        self.one_way().map_or_else(
            || {
                self.net(|s| Dd::from_f64(s.specific_enthalpy_j_kg))
                    .to_f64()
            },
            |s| dm * s.specific_enthalpy_j_kg,
        )
    }

    pub(crate) fn relative_enthalpy(&self, dm: f64, old_h: f64) -> f64 {
        self.one_way().map_or_else(
            || {
                (self.net(|s| Dd::from_f64(s.specific_enthalpy_j_kg))
                    - Dd::from_f64(dm) * Dd::from_f64(old_h))
                .to_f64()
            },
            |s| dm * (s.specific_enthalpy_j_kg - old_h),
        )
    }

    pub(crate) fn carried_momentum(&self, dm: f64, axis: usize) -> f64 {
        self.one_way().map_or_else(
            || self.net(|s| Dd::from_f64(s.velocity_m_s[axis])).to_f64(),
            |s| dm * s.velocity_m_s[axis],
        )
    }

    pub(crate) fn carried_kinetic(&self, dm: f64) -> f64 {
        if let Some(stream) = self.one_way() {
            let mut carried = 0.5 * dm * stream.velocity_m_s.iter().map(|v| v * v).sum::<f64>();
            if stream.specific_kinetic_variance_j_kg != 0.0 {
                carried += dm * stream.specific_kinetic_variance_j_kg;
            }
            return carried;
        }
        self.net(|s| {
            Dd::from_f64(0.5)
                * s.velocity_m_s
                    .into_iter()
                    .map(Dd::from_f64)
                    .fold(Dd::ZERO, |sum, v| sum + v * v)
                + Dd::from_f64(s.specific_kinetic_variance_j_kg)
        })
        .to_f64()
    }

    pub(crate) fn state(
        &self,
        old_mass: f64,
        dm: f64,
        new_mass: f64,
        velocity: [f64; 3],
    ) -> Option<([f64; 3], f64, f64)> {
        if let Some(stream) = self.one_way() {
            // Keep the original arithmetic and singleton limits exactly.
            return if self.streams[1].is_none() {
                let (next, mut heat) =
                    inflow::mixing_state(old_mass, dm, new_mass, velocity, stream.velocity_m_s)?;
                if stream.specific_kinetic_variance_j_kg != 0.0 {
                    heat += dm * stream.specific_kinetic_variance_j_kg;
                }
                heat.is_finite().then_some((next, heat, 0.0))
            } else {
                let (next, mut work) =
                    outflow::ejection_state(old_mass, dm, new_mass, velocity, stream.velocity_m_s)?;
                if stream.specific_kinetic_variance_j_kg != 0.0 {
                    work += -dm * stream.specific_kinetic_variance_j_kg;
                }
                (work.is_finite() && work >= 0.0).then_some((next, 0.0, work))
            };
        }
        if !dm.is_finite() || velocity.iter().any(|u| !u.is_finite()) {
            return None;
        }
        let weights = PathWeights::new(old_mass, new_mass)?;
        let delta: [f64; 3] = core::array::from_fn(|axis| {
            ((self.net(|s| Dd::from_f64(s.velocity_m_s[axis]))
                - Dd::from_f64(dm) * Dd::from_f64(velocity[axis]))
                / Dd::from_f64(new_mass))
            .to_f64()
        });
        let next = core::array::from_fn(|axis| velocity[axis] + delta[axis]);
        let energy = |side: usize| {
            self.streams[side].map_or(0.0, |stream| {
                let relative =
                    core::array::from_fn(|axis| stream.velocity_m_s[axis] - velocity[axis]);
                self.mass_kg[side]
                    * (0.5 * weights.integral(relative, delta)
                        + stream.specific_kinetic_variance_j_kg)
            })
        };
        let (heat, work) = (energy(0), energy(1));
        (next.iter().all(|u| u.is_finite())
            && heat.is_finite()
            && work.is_finite()
            && heat >= 0.0
            && work >= 0.0)
            .then_some((next, heat, work))
    }
}

/// Integrate |v_stream-u(s)|² along m(s)=m0+s*(m1-m0) and the SAME
/// simultaneous momentum path. With phi=s*m1/m(s), u=u0+phi*(u1-u0).
/// Store the Gram matrix of (1-phi, phi), pivoted on its larger diagonal.
/// This evaluates a sum of nonnegative squares, including m1=m0, without
/// subtracting large stream energies or dividing by near-zero net transfer.
struct PathWeights {
    increasing: bool,
    large: f64,
    cross_over_large: f64,
    remainder: f64,
}

impl PathWeights {
    fn new(m0: f64, m1: f64) -> Option<Self> {
        if !m0.is_finite() || !m1.is_finite() || m0 <= 0.0 || m1 <= 0.0 {
            return None;
        }
        let t = m0.min(m1) / m0.max(m1);
        if t == 0.0 {
            return None;
        }
        let x = 1.0 - t;
        let (small, cross) = if t >= 0.5 {
            // Positive Taylor sums avoid cubic cancellation near equal mass.
            // x<=1/2; 64 terms leave less than 2^-64 absolute series tail.
            let (mut a, mut b) = (0.0, 0.0);
            for j in (0..64).rev() {
                let denominator = f64::from((j + 2) * (j + 3));
                a = x.mul_add(a, 2.0 / denominator);
                b = x.mul_add(b, f64::from(j + 1) / denominator);
            }
            (t * a, t * b)
        } else {
            let log_t = fs_math::det::ln(t);
            let denominator = x * x * x;
            (
                t * (1.0 - t * t + 2.0 * t * log_t) / denominator,
                t * (-(1.0 + t) * log_t - 2.0 * x) / denominator,
            )
        };
        let large = 1.0 - small - 2.0 * cross;
        let cross_over_large = cross / large;
        let remainder = small - cross * cross_over_large;
        (large.is_finite() && large > 0.0 && cross >= 0.0 && remainder >= 0.0).then_some(Self {
            increasing: m1 >= m0,
            large,
            cross_over_large,
            remainder,
        })
    }

    fn integral(&self, initial: [f64; 3], delta: [f64; 3]) -> f64 {
        (0..3)
            .map(|axis| {
                let final_relative = initial[axis] - delta[axis];
                let (small, large) = if self.increasing {
                    (initial[axis], final_relative)
                } else {
                    (final_relative, initial[axis])
                };
                let combined = self.cross_over_large.mul_add(small, large);
                self.large * combined * combined + self.remainder * small * small
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::{MassLocation3, TransferMomentum};
    use super::*;
    use std::collections::BTreeMap;

    fn near(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    fn norm2(v: [f64; 3]) -> f64 {
        v.iter().map(|x| x * x).sum()
    }

    fn make_cell(a: f64, b: f64, incoming: [f64; 3], outgoing: [f64; 3]) -> ThroughflowCell {
        let mut result = ThroughflowCell::default();
        for (side, mass, velocity, h) in [(0, a, incoming, 120_000.0), (1, b, outgoing, 240_000.0)]
        {
            if mass != 0.0 {
                result
                    .add(side, IncomingMixture::new(mass, velocity, h, 0.0).unwrap())
                    .unwrap();
            }
        }
        result.represent(a - b).unwrap();
        result
    }

    // Independent continuous control-volume equations, integrated by RK4.
    // This reference neither uses the analytic weights nor any energy balance
    // to infer a missing heat/work term: both powers are integrated directly.
    fn reference(cell: &ThroughflowCell, mass: f64, velocity: [f64; 3]) -> [f64; 5] {
        let [a, b] = cell.mass_kg;
        let [incoming, outgoing] = cell.streams.map(Option::unwrap);
        let rhs = |s: f64, y: [f64; 5]| {
            let m = mass + s * (a - b);
            let ri: [f64; 3] = core::array::from_fn(|k| incoming.velocity_m_s[k] - y[k]);
            let ro: [f64; 3] = core::array::from_fn(|k| outgoing.velocity_m_s[k] - y[k]);
            [
                (a * ri[0] - b * ro[0]) / m,
                (a * ri[1] - b * ro[1]) / m,
                (a * ri[2] - b * ro[2]) / m,
                a * (0.5 * norm2(ri) + incoming.specific_kinetic_variance_j_kg),
                b * (0.5 * norm2(ro) + outgoing.specific_kinetic_variance_j_kg),
            ]
        };
        let mut y = [velocity[0], velocity[1], velocity[2], 0.0, 0.0];
        let steps = 4096;
        let h = 1.0 / f64::from(steps);
        for step in 0..steps {
            let s = f64::from(step) * h;
            let k1 = rhs(s, y);
            let k2 = rhs(
                s + 0.5 * h,
                core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]),
            );
            let k3 = rhs(
                s + 0.5 * h,
                core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]),
            );
            let k4 = rhs(s + h, core::array::from_fn(|i| y[i] + h * k3[i]));
            for i in 0..5 {
                y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
            }
        }
        y
    }

    #[test]
    fn g1_simultaneous_path_matches_independent_continuous_momentum_heat_and_work() {
        for (a, b) in [
            (0.3, 0.3),
            (0.8, 0.2),
            (0.2, 0.8),
            (2.0, 2.0),
            (0.2 * (1.0 + 2.0_f64.powi(-40)), 0.2),
            (0.2 * (1.0 - 2.0_f64.powi(-40)), 0.2),
        ] {
            let cell = make_cell(a, b, [0.7, 0.2, -0.2], [-0.3, 0.4, 0.1]);
            let initial = [0.4, -0.3, 0.2];
            let location = MassLocation3::Cell(0, 0);
            let cells = BTreeMap::from([(location, cell)]);
            let model = TransferMomentum::Throughflow {
                cells: &cells,
                maximum_work_j: 100.0,
            };
            let result = model
                .state(&location, 1.0, a - b, 1.0 + a - b, initial)
                .unwrap();
            let expected = reference(&cell, 1.0, initial);
            for k in 0..3 {
                near(result.0[k], expected[k], 2e-11);
                near(
                    (1.0 + a - b) * result.0[k] - initial[k],
                    cell.carried_momentum(a - b, k),
                    4e-15,
                );
            }
            near(result.1, expected[3], 2e-11);
            near(result.2, expected[4], 2e-11);
            let delta = core::array::from_fn(|k| result.0[k] - initial[k]);
            let dk = inflow::kinetic_change(a - b, 1.0 + a - b, initial, delta);
            near(dk + result.1 - result.2, cell.carried_kinetic(a - b), 4e-15);
        }
    }

    #[test]
    fn g1_equal_mass_pump_has_impulse_caloric_transfer_and_separate_positive_heat_work() {
        let cell = make_cell(0.3, 0.3, [1.0, 0.0, 0.0], [0.0; 3]);
        let result = cell.state(1.0, 0.0, 1.0, [0.0; 3]).unwrap();
        assert_eq!(cell.mass_kg, [0.3, 0.3]);
        assert_eq!(result.0, [0.3, 0.0, 0.0]);
        near(result.1, 0.15 * (1.0 - 0.3 + 0.09 / 3.0), 2e-17);
        near(result.2, 0.15 * 0.09 / 3.0, 2e-18);
        assert!(cell.carried_enthalpy(0.0) < -35_000.0);
        assert_eq!(
            cell.relative_enthalpy(0.0, 1e100),
            cell.carried_enthalpy(0.0)
        );
        let sequential = inflow::mixing_state(1.0, 0.3, 1.3, [0.0; 3], [1.0, 0.0, 0.0]).unwrap();
        assert!((sequential.1 - result.1).abs() > 0.005);
        // Identical inlet/outlet speeds still dissipate relative motion at the
        // inlet and require the same work at the outlet, with no net impulse.
        let common = make_cell(0.2, 0.2, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let still = common.state(1.0, 0.0, 1.0, [0.0; 3]).unwrap();
        assert_eq!(still.0, [0.0; 3]);
        near(still.1, 0.1, 3e-17);
        assert_eq!(still.1, still.2);
    }

    #[test]
    fn g0_one_way_cells_preserve_original_motion_and_all_stream_arithmetic() {
        let location = MassLocation3::Cell(0, 0);
        for side in [0, 1] {
            let mut cell = ThroughflowCell::default();
            for velocity in [[0.03, -0.04, 0.02], [-0.01, 0.02, 0.05]] {
                cell.add(
                    side,
                    IncomingMixture::new(0.05, velocity, 150_000.0, 0.01).unwrap(),
                )
                .unwrap();
            }
            let stream = cell.streams[side].unwrap();
            let mixtures = BTreeMap::from([(location, stream)]);
            let original = if side == 0 {
                TransferMomentum::Mixing(&mixtures)
            } else {
                TransferMomentum::Ejection {
                    mixtures: &mixtures,
                    maximum_work_j: 1.0,
                }
            };
            for fraction in [0.37, 1.0, 1.0 + 2.0 * f64::EPSILON] {
                let dm = if side == 0 {
                    0.1 * fraction
                } else {
                    -0.1 * fraction
                };
                cell.represent(dm).unwrap();
                let velocity = [0.02, -0.01, 0.03];
                assert_eq!(
                    cell.state(0.4, dm, 0.4 + dm, velocity),
                    original.state(&location, 0.4, dm, 0.4 + dm, velocity)
                );
                assert_eq!(
                    cell.carried_enthalpy(dm),
                    dm * stream.specific_enthalpy_j_kg
                );
                assert_eq!(
                    cell.relative_enthalpy(dm, 200_000.0),
                    dm * (stream.specific_enthalpy_j_kg - 200_000.0)
                );
                for k in 0..3 {
                    assert_eq!(cell.carried_momentum(dm, k), dm * stream.velocity_m_s[k]);
                }
                let kinetic = 0.5 * dm * norm2(stream.velocity_m_s)
                    + dm * stream.specific_kinetic_variance_j_kg;
                assert_eq!(cell.carried_kinetic(dm), kinetic);
                assert_eq!(cell.mass_defect(dm), (dm.abs() - 0.1).abs());
            }
        }
    }

    #[test]
    fn g3_gross_stream_variance_permutations_density_boost_and_datum_are_retained() {
        let sources = [
            (0, 0.0625, [0.125, -0.25, 0.375], 120_000.0),
            (0, 0.1875, [-0.375, 0.25, 0.125], 240_000.0),
            (1, 0.0625, [0.25, 0.125, -0.125], 180_000.0),
            (1, 0.0625, [-0.125, -0.25, 0.25], 160_000.0),
        ];
        let prepare = |order: [usize; 4], scale: f64, boost: [f64; 3], datum: f64| {
            let mut cell = ThroughflowCell::default();
            for i in order {
                let (side, mass, velocity, h) = sources[i];
                cell.add(
                    side,
                    IncomingMixture::new(
                        scale * mass,
                        core::array::from_fn(|k| velocity[k] + boost[k]),
                        h + datum,
                        scale * 0.01,
                    )
                    .unwrap(),
                )
                .unwrap();
            }
            let dm = 0.125 * scale;
            cell.represent(dm).unwrap();
            let initial = core::array::from_fn(|k| [0.0625, -0.125, 0.25][k] + boost[k]);
            let state = cell.state(scale, dm, 1.125 * scale, initial).unwrap();
            (cell, state)
        };
        let base = prepare([0, 1, 2, 3], 1.0, [0.0; 3], 0.0);
        for scale in [1e-10, 0.5, 1.0, 1e10] {
            for boost in [[0.0; 3], [4.0, -2.0, 1.0]] {
                for datum in [0.0, -1e8, 1e8] {
                    for order in [[0, 1, 2, 3], [3, 1, 2, 0], [2, 3, 1, 0]] {
                        let (cell, state) = prepare(order, scale, boost, datum);
                        for k in 0..3 {
                            near(state.0[k], base.1.0[k] + boost[k], 2e-15);
                        }
                        near(state.1, scale * base.1.1, 2e-16 * scale);
                        near(state.2, scale * base.1.2, 2e-16 * scale);
                        near(
                            cell.relative_enthalpy(0.125 * scale, 200_000.0 + datum),
                            scale * base.0.relative_enthalpy(0.125, 200_000.0),
                            3e-8 * scale,
                        );
                        near(cell.heat(), 0.04 * scale, 1e-17 * scale);
                        let kinetic: f64 = sources
                            .iter()
                            .map(|(side, m, u, _)| {
                                let v = core::array::from_fn(|k| u[k] + boost[k]);
                                let sign = if *side == 0 { 1.0 } else { -1.0 };
                                0.5 * scale * m * sign * norm2(v)
                            })
                            .sum();
                        near(cell.carried_kinetic(0.125 * scale), kinetic, 2e-15 * scale);
                    }
                }
            }
        }
    }

    #[test]
    fn g0_near_balanced_mass_projection_and_analytic_weights_remain_regular() {
        let mut cell = make_cell(1.0 + 2.0_f64.powi(-40), 1.0, [0.5; 3], [-0.5; 3]);
        let requested = cell.requested();
        let net = requested[0] - requested[1];
        cell.represent(net + 2.0_f64.powi(-52)).unwrap();
        assert_eq!(cell.mass_kg[1], requested[1]);
        assert_eq!(cell.mass_kg[0] - cell.mass_kg[1], net + 2.0_f64.powi(-52));
        assert_eq!(cell.mass_defect(net + 2.0_f64.powi(-52)), 2.0_f64.powi(-52));
        // Exact constant relative velocity integrates to itself at every ratio;
        // the endpoint-exchange identity checks the mass-decreasing branch.
        for ratio in [
            1e-200,
            1e-12,
            0.01,
            0.499999,
            0.5,
            0.999999999999,
            1.0,
            2.0,
            100.0,
            1e12,
            1e200,
        ] {
            let weights = PathWeights::new(1.0, ratio).unwrap();
            near(weights.integral([1.0, 2.0, -3.0], [0.0; 3]), 14.0, 2e-14);
            let reverse = PathWeights::new(ratio, 1.0).unwrap();
            near(
                weights.integral([0.2, -0.3, 0.4], [0.1, 0.3, -0.2]),
                reverse.integral([0.1, -0.6, 0.6], [-0.1, -0.3, 0.2]),
                2e-15,
            );
        }
    }

    #[test]
    fn g4_unrepresentable_paths_refuse_without_changing_stream_batches() {
        let mut cell = make_cell(0.3, 0.3, [1.0; 3], [-1.0; 3]);
        let accepted = cell;
        for dm in [f64::NAN, f64::INFINITY, 0.1] {
            assert!(cell.represent(dm).is_none());
            assert_eq!(cell, accepted);
        }
        for (m0, m1) in [
            (0.0, 1.0),
            (1.0, 0.0),
            (1.0, -1.0),
            (f64::INFINITY, 1.0),
            (1e-300, 1e300),
        ] {
            assert!(cell.state(m0, m1 - m0, m1, [0.0; 3]).is_none());
        }
        assert!(cell.state(1.0, 0.0, 1.0, [f64::NAN; 3]).is_none());
        let large = make_cell(1e300, 1e300, [1e100; 3], [-1e100; 3]);
        assert!(large.state(1e300, 0.0, 1e300, [0.0; 3]).is_none());
        assert_eq!(cell, accepted);
    }
}
