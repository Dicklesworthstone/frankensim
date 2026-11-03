//! Fluid relaxation removed by the original immersed solid collision.
//!
//! Apply relaxation to physical g = f + S/2, never to stored forced f. Every
//! admitted collision conserves momentum there, so replacing its increment
//! leaves the original Guo source outside the solid blend, counted once.

use super::{CollisionModel3, OPP3, Q3};

#[derive(Clone, Copy)]
pub(super) struct ImmersedRelaxation3 {
    /// Even/shear relaxation controls the original saturation weight.
    pub(super) tau: f64,
    odd_rate: Option<f64>,
    physical_moments: bool,
}

impl ImmersedRelaxation3 {
    /// Grid construction already validates rates. The saturation law uses the
    /// second-order viscous time, also for central moments and cumulants.
    /// Extremely small moment rates can have an unrepresentable reciprocal;
    /// those cannot supply this finite saturation weight.
    pub(super) fn new(model: CollisionModel3) -> Option<Self> {
        match model {
            CollisionModel3::Bgk { tau } => Some(Self {
                tau,
                odd_rate: None,
                physical_moments: false,
            }),
            CollisionModel3::Trt { tau, magic } => Some(Self {
                tau,
                odd_rate: Some(1.0 / (0.5 + magic / (tau - 0.5))),
                physical_moments: false,
            }),
            CollisionModel3::CentralMoment {
                second_order_rate, ..
            }
            | CollisionModel3::ReducedCumulant {
                second_order_rate, ..
            } => {
                let tau = 1.0 / second_order_rate;
                (tau.is_finite() && tau > 0.5).then_some(Self {
                    tau,
                    odd_rate: None,
                    physical_moments: true,
                })
            }
        }
    }

    /// One stack array per wet cell, shared by every body in that cell. The
    /// BGK branch retains its original operations. TRT evaluates the same
    /// nine opposite pairs as the native collision, with no extra collision,
    /// moment transform, allocation, or division inside the direction loop.
    ///
    /// Native moment collisions already return C(g) + S/2. Recover their
    /// actual nonlinear increment C(g)-g by removing the same represented
    /// half-source g-f. Unlike substituting BGK for this increment, this keeps
    /// each selected higher-order relaxation and cumulant product. No second
    /// collision/transform or source evaluation is needed. The two rounded
    /// additions of S/2 can leave only ordinary population roundoff.
    pub(super) fn increments(
        self,
        physical: &[f64; Q3],
        equilibrium: &[f64; Q3],
        stored: &[f64; Q3],
        ordinary: &[f64; Q3],
    ) -> [f64; Q3] {
        if self.physical_moments {
            return core::array::from_fn(|q| {
                (ordinary[q] - physical[q]) - (physical[q] - stored[q])
            });
        }
        let Some(odd_rate) = self.odd_rate else {
            return core::array::from_fn(|q| (equilibrium[q] - physical[q]) / self.tau);
        };
        let even_rate = 1.0 / self.tau;
        let mut increments = [0.0; Q3];
        increments[0] = -even_rate * (physical[0] - equilibrium[0]);
        for q in (1..Q3).step_by(2) {
            let opposite = OPP3[q];
            let neq = physical[q] - equilibrium[q];
            let neq_opposite = physical[opposite] - equilibrium[opposite];
            let even = -even_rate * (0.5 * (neq + neq_opposite));
            let odd = -odd_rate * (0.5 * (neq - neq_opposite));
            increments[q] = even + odd;
            increments[opposite] = even - odd;
        }
        increments
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::{E3, W3, collide_cell3, equilibrium3};

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 8e-16,
            "{actual:e} != {expected:e}"
        );
    }

    fn moment(populations: &[f64; Q3], axis: usize) -> f64 {
        (0..Q3)
            .map(|q| populations[q] * f64::from([E3[q].0, E3[q].1, E3[q].2][axis]))
            .sum()
    }

    // Independent Hermite source, not the native source-frame polynomial.
    fn source(velocity: [f64; 3], force: [f64; 3]) -> [f64; Q3] {
        let uf: f64 = velocity.iter().zip(force).map(|(u, f)| u * f).sum();
        core::array::from_fn(|q| {
            let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
            let eu: f64 = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
            let ef: f64 = e.iter().zip(force).map(|(e, f)| e * f).sum();
            W3[q] * (3.0 * ef - 3.0 * uf + 9.0 * eu * ef)
        })
    }

    fn perturb(mut physical: [f64; Q3]) -> [f64; Q3] {
        // Independent diagonal stress, shear stress, third and fourth modes;
        // each has zero mass and first moment. Nonzero stress also activates
        // the reduced cumulant's nonlinear fourth-order product correction.
        physical[1] += 2e-4;
        physical[2] += 2e-4;
        physical[0] -= 4e-4;
        for q in [7, 8] {
            physical[q] += 1.7e-4;
        }
        for q in [9, 10] {
            physical[q] -= 1.7e-4;
        }
        physical[7] += 1e-4;
        physical[8] -= 1e-4;
        physical[9] += 1e-4;
        physical[10] -= 1e-4;
        physical[1] -= 2e-4;
        physical[2] += 2e-4;
        for q in [7, 8, 9, 10] {
            physical[q] += 3e-4;
        }
        for q in [1, 2, 3, 4] {
            physical[q] -= 6e-4;
        }
        physical[0] += 12e-4;
        physical
    }

    #[test]
    fn native_moment_increment_keeps_all_fill_limits_and_original_force() {
        let (rho, u, wall) = (1.07, [0.024, -0.017, 0.011], [-0.011, 0.024, -0.006]);
        let equilibrium = equilibrium3(rho, u);
        let physical = perturb(equilibrium);
        let wall_equilibrium = equilibrium3(rho, wall);
        for rate in [0.4, 1.2, 1.91] {
            for model in [
                CollisionModel3::CentralMoment {
                    second_order_rate: rate,
                    higher_order_rate: 0.73,
                },
                CollisionModel3::ReducedCumulant {
                    second_order_rate: rate,
                    third_order_rate: 0.87,
                    fourth_order_rate: 1.63,
                },
            ] {
                let relaxation = ImmersedRelaxation3::new(model).unwrap();
                close(relaxation.tau, 1.0 / rate);
                // This is a separate call to the original *unforced* owner,
                // rather than repeating the production recovery expression.
                let reference = collide_cell3(physical, model, [0.0; 3]).unwrap();
                let bgk = collide_cell3(
                    physical,
                    CollisionModel3::Bgk {
                        tau: relaxation.tau,
                    },
                    [0.0; 3],
                )
                .unwrap();
                assert!((0..Q3).any(|q| (reference[q] - bgk[q]).abs() > 1e-5));
                for force in [[0.0; 3], [6e-4, -3e-4, 2e-4]] {
                    let source = source(u, force);
                    let stored = core::array::from_fn(|q| physical[q] - 0.5 * source[q]);
                    let ordinary = collide_cell3(stored, model, force).unwrap();
                    let increment =
                        relaxation.increments(&physical, &equilibrium, &stored, &ordinary);
                    for q in 0..Q3 {
                        close(increment[q], reference[q] - physical[q]);
                    }
                    close(increment.iter().sum(), 0.0);
                    for axis in 0..3 {
                        close(moment(&increment, axis), 0.0);
                    }
                    for epsilon in [0.0, 0.37, 0.999, 1.0] {
                        let b = epsilon * (relaxation.tau - 0.5)
                            / (1.0 - epsilon + relaxation.tau - 0.5);
                        let post = core::array::from_fn(|q| {
                            let solid = physical[OPP3[q]] - physical[q] + wall_equilibrium[q]
                                - equilibrium[OPP3[q]];
                            ordinary[q] + b * (solid - increment[q])
                        });
                        if epsilon == 0.0 {
                            assert_eq!(post, ordinary);
                        }
                        let change = core::array::from_fn(|q| post[q] - stored[q]);
                        close(change.iter().sum(), 0.0);
                        for axis in 0..3 {
                            close(
                                moment(&change, axis),
                                force[axis] + b * rho * (wall[axis] - u[axis]),
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn existing_parity_increments_match_the_native_collision_owners() {
        let eq = equilibrium3(1.07, [0.024, -0.017, 0.011]);
        let physical = perturb(eq);
        let unused = [f64::NAN; Q3];
        for tau in [0.5007, 0.83, 4.3] {
            for model in [
                CollisionModel3::Bgk { tau },
                CollisionModel3::Trt { tau, magic: 0.003 },
                CollisionModel3::Trt {
                    tau,
                    magic: 3.0 / 16.0,
                },
                CollisionModel3::Trt { tau, magic: 3.0 },
            ] {
                let increment = ImmersedRelaxation3::new(model)
                    .unwrap()
                    .increments(&physical, &eq, &unused, &unused);
                let reference = collide_cell3(physical, model, [0.0; 3]).unwrap();
                for q in 0..Q3 {
                    close(physical[q] + increment[q], reference[q]);
                }
            }
        }
    }

    #[test]
    fn unrepresentable_moment_relaxation_time_refuses() {
        assert!(
            ImmersedRelaxation3::new(CollisionModel3::ReducedCumulant {
                second_order_rate: f64::from_bits(1),
                third_order_rate: 1.0,
                fourth_order_rate: 1.0,
            })
            .is_none()
        );
    }
}
