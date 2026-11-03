//! Kinetic transport from the same admitted upwind traces as momentum.
//!
//! Keep the original bounded quadrature so every carried velocity can be
//! shifted before products and reduction, without subtracting large
//! frame-dependent kinetic totals. These are floating quadrature samples,
//! not a complete thermodynamic stream or a whole-scheme energy certificate.

use super::{BoundaryWrench, StokesBoundaryError};

/// Outward convective transport per unit depth, divided by source density.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BoundaryTransport {
    /// Physical origin [m] for the angular momentum flux.
    pub origin_m: [f64; 2],
    /// Uniform offset [m/s] added to carried velocities. The retained advector's
    /// normal flow relative to the control surface is unchanged.
    pub velocity_offset_m_s: [f64; 2],
    /// Signed outward volume flow [m²/s].
    pub outward_volume_flux_m2_s: f64,
    /// Outward linear/angular momentum transport divided by density.
    pub outward_momentum: BoundaryWrench,
    /// Signed outward carried kinetic power divided by density [m⁴/s³]:
    /// integral 0.5 * (w.n) * |u_upwind + velocity_offset|².
    pub outward_kinetic_energy_flux: f64,
    /// Nonnegative boundary upwind loss divided by density [m⁴/s³]:
    /// integral 0.5 * max(-w.n, 0) * |u - u_upwind|².
    /// This is numerical boundary dissipation, unchanged by a common boost;
    /// it excludes interior upwinding, viscosity and backflow resistance.
    pub upwind_kinetic_dissipation: f64,
}

impl BoundaryTransport {
    pub(super) fn finite(&self) -> bool {
        self.origin_m
            .iter()
            .chain(&self.velocity_offset_m_s)
            .all(|v| v.is_finite())
            && self.outward_volume_flux_m2_s.is_finite()
            && self.outward_momentum.force.iter().all(|v| v.is_finite())
            && self.outward_momentum.moment.is_finite()
            && self.outward_kinetic_energy_flux.is_finite()
            && self.upwind_kinetic_dissipation.is_finite()
            && self.upwind_kinetic_dissipation >= 0.0
    }
}

/// Fixed-order compensated reductions; no per-face heap storage.
pub(super) struct Sums<const N: usize> {
    values: [f64; N],
    corrections: [f64; N],
}

impl<const N: usize> Sums<N> {
    pub(super) const fn new() -> Self {
        Self {
            values: [0.0; N],
            corrections: [0.0; N],
        }
    }

    pub(super) fn add(&mut self, values: [f64; N]) {
        for (i, value) in values.into_iter().enumerate() {
            let old = self.values[i];
            let next = old + value;
            self.corrections[i] += if old.abs() >= value.abs() {
                (old - next) + value
            } else {
                (value - next) + old
            };
            self.values[i] = next;
        }
    }

    pub(super) fn values(self) -> [f64; N] {
        core::array::from_fn(|i| self.values[i] + self.corrections[i])
    }
}

/// The original upwind quadrature, retained before any frame-dependent sum.
/// At most two Gauss points on each sign branch are needed for affine BDM1.
#[derive(Clone, Copy, Debug)]
pub(super) struct TransportQuadrature {
    samples: [TransportSample; crate::ns::CONVECTION_EDGE_QUADRATURE_CAPACITY],
    count: usize,
    upwind_loss: f64,
}

#[derive(Clone, Copy, Debug)]
struct TransportSample {
    point: [f64; 2],
    velocity: [f64; 2],
    volume_flux: f64,
}

pub(super) struct TransportAccumulator {
    quadrature: TransportQuadrature,
    loss: Sums<1>,
}

impl TransportAccumulator {
    pub(super) const fn new() -> Self {
        Self {
            quadrature: TransportQuadrature {
                samples: [TransportSample {
                    point: [0.0; 2],
                    velocity: [0.0; 2],
                    volume_flux: 0.0,
                }; crate::ns::CONVECTION_EDGE_QUADRATURE_CAPACITY],
                count: 0,
                upwind_loss: 0.0,
            },
            loss: Sums::new(),
        }
    }

    pub(super) fn sample(
        &mut self,
        point: [f64; 2],
        volume_flux: f64,
        upwind: [f64; 2],
        interior: [f64; 2],
    ) -> Result<(), StokesBoundaryError> {
        if !volume_flux.is_finite()
            || !point
                .iter()
                .chain(&upwind)
                .chain(&interior)
                .all(|v| v.is_finite())
        {
            return Err(invalid("nonfinite boundary transport sample"));
        }
        if volume_flux == 0.0 {
            return Ok(());
        }
        let count = self.quadrature.count;
        let sample = self
            .quadrature
            .samples
            .get_mut(count)
            .ok_or_else(|| invalid("boundary transport quadrature capacity exceeded"))?;
        *sample = TransportSample {
            point,
            velocity: upwind,
            volume_flux,
        };
        self.quadrature.count += 1;
        if volume_flux < 0.0 {
            self.loss
                .add([half_weighted_square(-volume_flux, interior[0] - upwind[0])
                    + half_weighted_square(-volume_flux, interior[1] - upwind[1])]);
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<TransportQuadrature, StokesBoundaryError> {
        self.quadrature.upwind_loss = self.loss.values()[0];
        if !self.quadrature.upwind_loss.is_finite() || self.quadrature.upwind_loss < 0.0 {
            return Err(invalid("unrepresentable boundary upwind dissipation"));
        }
        Ok(self.quadrature)
    }
}

impl TransportQuadrature {
    pub(super) fn volume_flux_moments(&self, origin: [f64; 2]) -> (f64, [f64; 2]) {
        let mut sums = Sums::<3>::new();
        for sample in &self.samples[..self.count] {
            sums.add([
                sample.volume_flux,
                (sample.point[0] - origin[0]) * sample.volume_flux,
                (sample.point[1] - origin[1]) * sample.volume_flux,
            ]);
        }
        let values = sums.values();
        (values[0], [values[1], values[2]])
    }

    pub(super) fn evaluate(
        &self,
        origin: [f64; 2],
        offset: [f64; 2],
    ) -> Result<BoundaryTransport, StokesBoundaryError> {
        if !origin.iter().chain(&offset).all(|v| v.is_finite()) {
            return Err(invalid(
                "finite transport origin and velocity offset required",
            ));
        }
        let mut sums = Sums::<5>::new();
        for sample in &self.samples[..self.count] {
            // Shift EACH carried velocity before products and reduction. A
            // tiny quadrature contribution must not vanish inside a large
            // net flow before its shifted momentum/energy can be evaluated.
            let u = [
                sample.velocity[0] + offset[0],
                sample.velocity[1] + offset[1],
            ];
            if !u.iter().all(|v| v.is_finite()) {
                return Err(invalid("unrepresentable boundary transport velocity"));
            }
            let q = sample.volume_flux;
            let force = u.map(|v| q * v);
            sums.add([
                q,
                force[0],
                force[1],
                cross(
                    [sample.point[0] - origin[0], sample.point[1] - origin[1]],
                    force,
                ),
                half_weighted_square(q, u[0]) + half_weighted_square(q, u[1]),
            ]);
        }
        let values = sums.values();
        let result = BoundaryTransport {
            origin_m: origin,
            velocity_offset_m_s: offset,
            outward_volume_flux_m2_s: values[0],
            outward_momentum: BoundaryWrench {
                force: [values[1], values[2]],
                moment: values[3],
            },
            outward_kinetic_energy_flux: values[4],
            upwind_kinetic_dissipation: self.upwind_loss,
        };
        if !result.finite() {
            return Err(invalid("unrepresentable boundary transport frame"));
        }
        Ok(result)
    }
}

// Apply the weight before the second velocity factor: neither |u|² nor a
// reciprocal density is formed just to obtain a representable kinetic flux.
fn half_weighted_square(weight: f64, value: f64) -> f64 {
    (0.5 * value * weight) * value
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}

fn invalid(reason: &'static str) -> StokesBoundaryError {
    StokesBoundaryError::InvalidInput(reason)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn common_velocity_offset_preserves_small_carried_energy_and_upwind_loss() {
        let boost = 1_099_511_627_776.0;
        let mut sums = TransportAccumulator::new();
        sums.sample(
            [2.0, 3.0],
            -2.0,
            [boost + 3.0, -boost - 2.0],
            [boost + 1.0, -boost + 2.0],
        )
        .unwrap();
        let moments = sums.finish().unwrap();
        let moving = moments.evaluate([1.0, 1.0], [-boost, boost]).unwrap();
        assert_eq!(moving.outward_volume_flux_m2_s, -2.0);
        assert_eq!(moving.outward_momentum.force, [-6.0, 4.0]);
        assert_eq!(moving.outward_momentum.moment, 16.0);
        assert_eq!(moving.outward_kinetic_energy_flux, -13.0);
        assert!((moving.upwind_kinetic_dissipation - 20.0).abs() < 1e-13);
        // u=(1,2), incoming=(3,-2), q=-2:
        // q*(u.g-|u|²/2) = 7 = incoming kinetic flux + upwind loss.
        assert!(
            (moving.outward_kinetic_energy_flux + moving.upwind_kinetic_dissipation - 7.0).abs()
                < 1e-13
        );
        let fixed = moments.evaluate([1.0, 1.0], [0.0; 2]).unwrap();
        assert_eq!(
            fixed.upwind_kinetic_dissipation,
            moving.upwind_kinetic_dissipation
        );
    }

    #[test]
    fn signed_flux_with_zero_net_volume_keeps_energy_and_refuses_bad_frames() {
        let mut sums = TransportAccumulator::new();
        sums.sample([0.0, 0.0], 1.0, [2.0, 0.0], [2.0, 0.0])
            .unwrap();
        sums.sample([0.0, 1.0], -1.0, [4.0, 0.0], [4.0, 0.0])
            .unwrap();
        let moments = sums.finish().unwrap();
        let transport = moments.evaluate([0.0; 2], [-2.0, 0.0]).unwrap();
        assert_eq!(transport.outward_volume_flux_m2_s, 0.0);
        assert_eq!(transport.outward_momentum.force, [-2.0, 0.0]);
        assert_eq!(transport.outward_momentum.moment, 2.0);
        assert_eq!(transport.outward_kinetic_energy_flux, -2.0);
        assert_eq!(transport.upwind_kinetic_dissipation, 0.0);
        assert!(moments.evaluate([f64::NAN, 0.0], [0.0; 2]).is_err());
        assert!(moments.evaluate([0.0; 2], [f64::INFINITY, 0.0]).is_err());
        let empty = TransportAccumulator::new().finish().unwrap();
        let zero = empty.evaluate([0.0; 2], [1e300, -1e300]).unwrap();
        assert_eq!(zero.outward_kinetic_energy_flux, 0.0);
    }

    #[test]
    fn tiny_flux_keeps_shifted_energy_for_both_signs_and_sample_orders() {
        let boost = 1_073_741_824.0; // 2^30, with exactly representable oracles.
        for sign in [-1.0, 1.0] {
            for reverse in [false, true] {
                let mut samples = [
                    ([0.0, 0.0], sign / (boost * boost), [0.0; 2]),
                    ([0.0, 1.0], 1.0, [boost, 0.0]),
                ];
                if reverse {
                    samples.reverse();
                }
                let mut sums = TransportAccumulator::new();
                for (point, flow, velocity) in samples {
                    sums.sample(point, flow, velocity, velocity).unwrap();
                }
                let transport = sums
                    .finish()
                    .unwrap()
                    .evaluate([0.0; 2], [-boost, 0.0])
                    .unwrap();
                // The volume sum rounds to 1. Reducing it before frame
                // conversion would erase the small branch's finite energy.
                assert_eq!(transport.outward_volume_flux_m2_s, 1.0);
                assert_eq!(transport.outward_momentum.force, [-sign / boost, 0.0]);
                assert_eq!(transport.outward_momentum.moment, 0.0);
                assert_eq!(transport.outward_kinetic_energy_flux, sign * 0.5);
                assert_eq!(transport.upwind_kinetic_dissipation, 0.0);
            }
        }
    }

    #[test]
    fn retained_quadrature_refuses_excess_samples_and_nonfinite_data() {
        let mut full = TransportAccumulator::new();
        for _ in 0..crate::ns::CONVECTION_EDGE_QUADRATURE_CAPACITY {
            full.sample([0.0; 2], 1.0, [1.0; 2], [1.0; 2]).unwrap();
        }
        assert!(full.sample([0.0; 2], 1.0, [1.0; 2], [1.0; 2]).is_err());
        assert!(
            TransportAccumulator::new()
                .sample([0.0; 2], f64::NAN, [0.0; 2], [0.0; 2])
                .is_err()
        );
        assert!(
            TransportAccumulator::new()
                .sample([f64::NAN, 0.0], 0.0, [0.0; 2], [0.0; 2])
                .is_err()
        );
        assert!(
            TransportAccumulator::new()
                .sample([0.0; 2], 1.0, [f64::INFINITY, 0.0], [0.0; 2])
                .is_err()
        );
        assert!(
            TransportAccumulator::new()
                .sample([0.0; 2], -1.0, [0.0; 2], [f64::NAN, 0.0])
                .is_err()
        );
    }
}
