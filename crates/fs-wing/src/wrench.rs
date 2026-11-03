//! Circulatory force and moment from accepted panel circulation and local flow.
//!
//! The bound-segment Kutta–Joukowsky force is applied at its midpoint. This
//! owner adds no profile drag, separate induced drag, section couple, thrust,
//! weight or apparent mass. It does not solve for circulation or broaden the
//! nonlinear solver's section-plane assumptions.

use crate::{Panel, Refusal};
use fs_ga::{Vec3, Wrench};

fn vec(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn finite(value: Vec3) -> bool {
    value.x.is_finite() && value.y.is_finite() && value.z.is_finite()
}

fn reject(message: &str) -> Refusal {
    Refusal {
        code: "circulatory-wrench-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "supply finite accepted circulation, local velocity, bound endpoints and origin".into(),
        ],
    }
}

/// Allocation-free streaming assembly of a circulatory wrench about one origin.
///
/// Geometry, origin and accepted incident velocities share one Cartesian frame;
/// the output uses that frame, with force [N] and torque about the origin [N m].
/// Uniform flow repeats the same velocity for each panel. Strip-local flow
/// streams each strip's owning panels with its accepted velocity, preserving
/// the caller's strip/panel order. No ownership map or additional solve is built.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CirculatoryWrenchAccumulator {
    density_kg_m3: f64,
    origin_m: Vec3,
    total: Wrench,
}

impl CirculatoryWrenchAccumulator {
    /// Start an empty zero wrench at a fixed origin and positive density.
    ///
    /// # Errors
    /// `circulatory-wrench-invalid` for nonfinite origin or nonpositive density.
    pub fn new(density_kg_m3: f64, origin_m: [f64; 3]) -> Result<Self, Refusal> {
        if !density_kg_m3.is_finite() || density_kg_m3 <= 0.0 || !finite(vec(origin_m)) {
            return Err(reject(
                "positive finite density and a finite origin are required",
            ));
        }
        Ok(Self {
            density_kg_m3,
            origin_m: vec(origin_m),
            total: Wrench::default(),
        })
    }

    /// Add one panel and return its individual wrench about the same origin.
    ///
    /// `gamma_m2_s` is the accepted signed circulation. `incident_mps` includes
    /// the local wash/gust/rigid-motion velocity used by the force model. The
    /// full bound vector B−A supports swept, canted and vertical panels; width,
    /// control point and normal are not substitutes for that vector. For each
    /// panel the circulatory force is perpendicular to its incident velocity.
    ///
    /// The existing native Vec3 arithmetic and accumulation order are retained.
    /// This call allocates nothing on success. A refusal leaves the accumulator
    /// unchanged, including when an individually finite contribution would
    /// overflow the total.
    ///
    /// # Errors
    /// `circulatory-wrench-invalid` for nonfinite inputs, a zero bound segment,
    /// or unrepresentable force, arm, torque or accumulated wrench.
    pub fn add_panel(
        &mut self,
        panel: &Panel,
        gamma_m2_s: f64,
        incident_mps: [f64; 3],
    ) -> Result<Wrench, Refusal> {
        let a = vec(panel.a);
        let b = vec(panel.b);
        let velocity = vec(incident_mps);
        if !finite(a) || !finite(b) || a == b || !gamma_m2_s.is_finite() || !finite(velocity) {
            return Err(reject(
                "finite circulation/velocity and a nonzero finite bound segment are required",
            ));
        }
        let segment = b - a;
        let force = velocity
            .cross(segment)
            .scale(self.density_kg_m3 * gamma_m2_s);
        let arm = (a + b).scale(0.5) - self.origin_m;
        let item = Wrench::new(arm.cross(force), force);
        let total = Wrench::new(
            self.total.torque + item.torque,
            self.total.force + item.force,
        );
        if !finite(segment)
            || !finite(arm)
            || !item.to_array().iter().all(|x| x.is_finite())
            || !total.to_array().iter().all(|x| x.is_finite())
        {
            return Err(reject(
                "circulatory force, moment arm or wrench accumulation overflow",
            ));
        }
        self.total = total;
        Ok(item)
    }

    /// Accumulated wrench, in insertion order, about the declared origin.
    #[must_use]
    pub const fn total(&self) -> Wrench {
        self.total
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SurfaceId;
    use fs_ga::{Se3, So3, So3Tangent, Twist};

    fn array(value: Vec3) -> [f64; 3] {
        [value.x, value.y, value.z]
    }

    fn panels() -> [Panel; 3] {
        core::array::from_fn(|i| Panel {
            surface: if i == 2 {
                SurfaceId::Vertical
            } else {
                SurfaceId::WingLower
            },
            a: [[-0.4, -2.0, 0.1], [0.2, 0.7, -1.0], [-0.5, 1.0, -1.0]][i],
            b: [[0.8, -0.5, 0.4], [0.6, 2.0, -0.4], [-0.5, 1.0, 1.0]][i],
            ctrl: [0.0; 3],
            normal: if i == 2 {
                [0.0, 1.0, 0.0]
            } else {
                [0.0, 0.0, -1.0]
            },
            width_m: 1.0,
        })
    }

    #[test]
    fn g1_spatial_panel_wrenches_are_workless_and_retain_native_accumulation_bits() {
        let panels = panels();
        let gamma = [0.75, -1.25, 0.5];
        let rho = 1.225;
        let origin = [0.2, -0.1, 0.05];
        // A non-contiguous first strip keeps caller ownership/order live.
        let strips: [(&[usize], [f64; 3]); 2] =
            [(&[2, 0], [12.0, 1.5, -0.7]), (&[1], [8.0, -0.5, 1.2])];
        let mut owner = CirculatoryWrenchAccumulator::new(rho, origin).unwrap();
        let mut old_force = Vec3::default();
        let mut old_moment = Vec3::default();
        for (indices, velocity) in strips {
            for &j in indices {
                let panel = &panels[j];
                // Exact previously native Vec3 force/moment loop, retained as
                // a migration oracle alongside independent physical identities.
                let segment = vec(panel.b) - vec(panel.a);
                let force = vec(velocity).cross(segment).scale(rho * gamma[j]);
                let arm = (vec(panel.a) + vec(panel.b)).scale(0.5) - vec(origin);
                old_force = old_force + force;
                old_moment = old_moment + arm.cross(force);
                let item = owner.add_panel(panel, gamma[j], velocity).unwrap();
                let scale = item.force.norm() * vec(velocity).norm();
                assert!(item.force.dot(vec(velocity)).abs() <= 2e-15 * scale);
                let omega = vec([0.3, -0.2, 0.4]);
                let at_origin = vec(velocity) - omega.cross(arm);
                assert!(item.pairing(Twist::new(omega, at_origin)).abs() <= 5e-15 * scale);
                if j == 2 {
                    assert_eq!(item.force.z, 0.0);
                    assert!(item.force.x > 0.0 && item.force.y < -10.0);
                }
            }
        }
        let reference = Wrench::new(old_moment, old_force);
        assert_eq!(
            owner.total().to_array().map(f64::to_bits),
            reference.to_array().map(f64::to_bits)
        );
    }

    #[test]
    fn g3_origin_shift_and_rigid_frame_transport_preserve_the_complete_wrench() {
        let panels = panels();
        let gamma = [0.75, -1.25, 0.5];
        let velocities = [[12.0, 1.5, -0.7], [8.0, -0.5, 1.2], [11.0, 2.0, 0.4]];
        let origin = vec([0.2, -0.1, 0.05]);
        let shift = vec([0.7, -0.3, 0.4]);
        let rotation = So3::exp(So3Tangent::new(vec([0.3, -0.4, 0.2]))).unwrap();
        let translation = vec([3.0, -4.0, 0.7]);
        let mut base = CirculatoryWrenchAccumulator::new(1.225, array(origin)).unwrap();
        let mut shifted = CirculatoryWrenchAccumulator::new(1.225, array(origin + shift)).unwrap();
        let transformed_origin = rotation.rotate(origin).unwrap() + translation;
        let mut transformed =
            CirculatoryWrenchAccumulator::new(1.225, array(transformed_origin)).unwrap();
        let mut at_world_zero = CirculatoryWrenchAccumulator::new(1.225, [0.0; 3]).unwrap();
        for i in 0..panels.len() {
            base.add_panel(&panels[i], gamma[i], velocities[i]).unwrap();
            shifted
                .add_panel(&panels[i], gamma[i], velocities[i])
                .unwrap();
            let moved = Panel {
                a: array(rotation.rotate(vec(panels[i].a)).unwrap() + translation),
                b: array(rotation.rotate(vec(panels[i].b)).unwrap() + translation),
                ..panels[i]
            };
            let velocity = array(rotation.rotate(vec(velocities[i])).unwrap());
            transformed.add_panel(&moved, gamma[i], velocity).unwrap();
            at_world_zero.add_panel(&moved, gamma[i], velocity).unwrap();
        }
        let value = base.total();
        assert!((shifted.total().force - value.force).norm() < 1e-13);
        assert!(
            (shifted.total().torque - (value.torque - shift.cross(value.force))).norm() < 1e-12
        );
        assert!((transformed.total().force - rotation.rotate(value.force).unwrap()).norm() < 1e-12);
        assert!(
            (transformed.total().torque - rotation.rotate(value.torque).unwrap()).norm() < 1e-12
        );
        // The public result is the existing force/moment dual, including its
        // coadjoint transport to a different world moment origin.
        let pose = Se3::from_parts(rotation, transformed_origin).unwrap();
        let expected = value.transform_by(&pose);
        assert!((at_world_zero.total().force - expected.force).norm() < 1e-12);
        assert!((at_world_zero.total().torque - expected.torque).norm() < 1e-12);
    }

    #[test]
    fn g0_invalid_or_overflowing_panel_additions_are_atomic() {
        for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(CirculatoryWrenchAccumulator::new(rho, [0.0; 3]).is_err());
        }
        assert!(CirculatoryWrenchAccumulator::new(1.0, [f64::NAN, 0.0, 0.0]).is_err());
        let panel = Panel {
            a: [0.0; 3],
            b: [0.0, 1.0, 0.0],
            ..panels()[0]
        };
        let mut owner = CirculatoryWrenchAccumulator::new(1.0, [0.0; 3]).unwrap();
        assert_eq!(owner.total(), Wrench::default());
        owner
            .add_panel(&panel, 0.75 * f64::MAX, [1.0, 0.0, 0.0])
            .unwrap();
        let before = owner;
        assert!(
            owner
                .add_panel(&panel, 0.75 * f64::MAX, [1.0, 0.0, 0.0])
                .is_err()
        );
        assert_eq!(owner, before);
        for (gamma, velocity) in [
            (f64::NAN, [1.0, 0.0, 0.0]),
            (1.0, [1.0, f64::INFINITY, 0.0]),
        ] {
            assert!(owner.add_panel(&panel, gamma, velocity).is_err());
            assert_eq!(owner, before);
        }
        let zero_span = Panel {
            b: panel.a,
            ..panel
        };
        assert!(owner.add_panel(&zero_span, 1.0, [1.0, 0.0, 0.0]).is_err());
        assert_eq!(owner, before);
    }
}
