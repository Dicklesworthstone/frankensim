//! Geometric passive shape for one original scalar viscous interface law.
//!
//! In material A's axes, z = [v, ell*omega] has uniform velocity units.
//! A dimensionless positive-semidefinite S defines power c(T) z^T S z.
//! Thus C = c(T) diag(I,ell I) S diag(I,ell I). The shape and lever length
//! are explicit fixed geometry/model data; c(T) remains the sourced material
//! coefficient in N s/m. No angular coefficient or empirical scale is inferred.

use fs_ga::So3;
use fs_math::dd::Dd;
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingCoordinateScales, BushingMatrixError, PassiveBushingMatrix,
};

/// Fixed passive geometric resistance in the ORIGINAL first body's material axes.
///
/// The normalized form acts on [relative translation velocity, ell times
/// relative angular velocity]. It may contain reciprocal cross-coupling and
/// exact zero directions. The source's original area is already included in
/// c(T); this dimensionless shape does not multiply the area a second time.
/// Its geometry stays fixed during thermal expansion, matching the original
/// damper's fixed-area interface approximation.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialResistanceShape {
    normalized: PassiveBushingMatrix,
    rotary_length_m: f64,
}

impl SpatialResistanceShape {
    /// Admit an explicitly authored complete passive shape and rotary length.
    ///
    /// The same native strict passive-form owner validates S. Active,
    /// asymmetric, unresolved and nonfinite shapes refuse without projection.
    pub fn new(
        normalized: [[f64; 6]; 6],
        rotary_length_m: f64,
    ) -> Result<Self, BushingMatrixError> {
        if !rotary_length_m.is_finite() || rotary_length_m <= 0.0 {
            return Err(BushingMatrixError("finite positive rotary length required"));
        }
        Ok(Self {
            normalized: PassiveBushingMatrix::try_new(
                normalized,
                BushingCoordinateScales {
                    translation_m: 1.0,
                    rotation_rad: 1.0,
                },
            )?,
            rotary_length_m,
        })
    }

    /// Original dimensionless coefficients, retained without fitting.
    pub const fn normalized_coefficients(&self) -> &[[f64; 6]; 6] {
        self.normalized.native_coefficients()
    }

    /// Declared fixed physical length converting angular rate to speed [m].
    pub const fn rotary_length_m(&self) -> f64 {
        self.rotary_length_m
    }

    /// Rotate the original material form into the actual native principal axes.
    /// Rebinding inertia may relabel principal axes; it must not rotate the
    /// physical material law along with that numerical frame.
    pub(super) fn at_coefficient(
        &self,
        coefficient_n_s_m: f64,
        body_from_principal: So3,
    ) -> Result<PassiveBushingMatrix, BushingMatrixError> {
        if !coefficient_n_s_m.is_finite() || coefficient_n_s_m <= 0.0 {
            return Err(BushingMatrixError(
                "finite positive sourced viscosity required",
            ));
        }
        let rotation = body_from_principal.matrix().m;
        let transform: [[f64; 6]; 6] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i / 3 == j / 3 {
                    rotation[3 * (i % 3) + j % 3]
                } else {
                    0.0
                }
            })
        });
        let shape = self.normalized.native_coefficients();
        let mut native = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..=i {
                let mut value = Dd::ZERO;
                for k in 0..6 {
                    for l in 0..6 {
                        value = value
                            + Dd::from_f64(transform[k][i])
                                * Dd::from_f64(shape[k][l])
                                * Dd::from_f64(transform[l][j]);
                    }
                }
                let mut scaled = value * Dd::from_f64(coefficient_n_s_m);
                if i >= 3 {
                    scaled = scaled * Dd::from_f64(self.rotary_length_m);
                }
                if j >= 3 {
                    scaled = scaled * Dd::from_f64(self.rotary_length_m);
                }
                let entry = scaled.to_f64();
                if !entry.is_finite() || ((value.hi != 0.0 || value.lo != 0.0) && entry == 0.0) {
                    return Err(BushingMatrixError(
                        "physical spatial resistance is unrepresentable",
                    ));
                }
                native[i][j] = entry;
                native[j][i] = entry;
            }
        }
        PassiveBushingMatrix::try_new(
            native,
            BushingCoordinateScales {
                translation_m: self.rotary_length_m,
                rotation_rad: 1.0,
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_ga::Vec3;

    fn action(matrix: &[[f64; 6]; 6], value: [f64; 6]) -> [f64; 6] {
        core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * value[j]).sum())
    }
    fn power(matrix: &[[f64; 6]; 6], value: [f64; 6]) -> f64 {
        action(matrix, value)
            .into_iter()
            .zip(value)
            .map(|(a, b)| a * b)
            .sum()
    }

    #[test]
    fn source_and_length_supply_distinct_translational_rotary_and_cross_units() {
        let mode = [1.0, -0.5, 0.0, 0.25, 0.0, 0.75];
        let s = core::array::from_fn(|i| {
            core::array::from_fn(|j| mode[i] * mode[j] + if i == j { 1.0 } else { 0.0 })
        });
        let shape = SpatialResistanceShape::new(s, 0.2).unwrap();
        let native = shape.at_coefficient(7.0, So3::identity()).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let expected =
                    7.0 * s[i][j] * if i < 3 { 1.0 } else { 0.2 } * if j < 3 { 1.0 } else { 0.2 };
                assert!((native.native_coefficients()[i][j] - expected).abs() < 2e-15);
            }
        }
        let doubled = shape.at_coefficient(14.0, So3::identity()).unwrap();
        let slip = [0.2, 0.1, -0.3, 0.4, -0.2, 0.6];
        assert!(
            (power(doubled.native_coefficients(), slip)
                - 2.0 * power(native.native_coefficients(), slip))
            .abs()
                < 1e-14
        );
        assert_eq!(shape.normalized_coefficients(), &s);
    }

    #[test]
    fn material_law_survives_relabelled_native_principal_axes() {
        let s = core::array::from_fn(|i| {
            core::array::from_fn(|j| if i == j { 1.0 + i as f64 } else { 0.0 })
        });
        let shape = SpatialResistanceShape::new(s, 0.3).unwrap();
        let rotation = So3::exp(fs_ga::So3Tangent::new(Vec3::new(0.27, -0.19, 0.31))).unwrap();
        let native = shape.at_coefficient(2.0, rotation).unwrap();
        let value = [0.4, -0.7, 0.2, 0.5, -0.1, 0.8];
        let linear = rotation
            .rotate(Vec3::new(value[0], value[1], value[2]))
            .unwrap();
        let angular = rotation
            .rotate(Vec3::new(value[3], value[4], value[5]))
            .unwrap();
        let material = [
            linear.x, linear.y, linear.z, angular.x, angular.y, angular.z,
        ];
        let original = shape.at_coefficient(2.0, So3::identity()).unwrap();
        assert!(
            (power(native.native_coefficients(), value)
                - power(original.native_coefficients(), material))
            .abs()
                < 2e-14
        );
        let expected = action(original.native_coefficients(), material);
        let force = rotation
            .inverse()
            .rotate(Vec3::new(expected[0], expected[1], expected[2]))
            .unwrap();
        let moment = rotation
            .inverse()
            .rotate(Vec3::new(expected[3], expected[4], expected[5]))
            .unwrap();
        let transformed = [force.x, force.y, force.z, moment.x, moment.y, moment.z];
        for (actual, expected) in action(native.native_coefficients(), value)
            .into_iter()
            .zip(transformed)
        {
            assert!((actual - expected).abs() < 2e-14);
        }
    }

    #[test]
    fn exact_free_directions_and_invalid_shapes_are_not_repaired() {
        let mut s = [[0.0; 6]; 6];
        s[0][0] = 1.0;
        let shape = SpatialResistanceShape::new(s, 0.4).unwrap();
        let native = shape.at_coefficient(3.0, So3::identity()).unwrap();
        assert_eq!(
            action(native.native_coefficients(), [0.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            [0.0; 6]
        );
        let zero = SpatialResistanceShape::new([[0.0; 6]; 6], 0.4).unwrap();
        assert_eq!(
            zero.at_coefficient(3.0, So3::identity())
                .unwrap()
                .native_coefficients(),
            &[[0.0; 6]; 6]
        );
        s[1][1] = -1.0;
        assert!(SpatialResistanceShape::new(s, 0.4).is_err());
        assert!(SpatialResistanceShape::new([[0.0; 6]; 6], 0.0).is_err());
        assert!(shape.at_coefficient(f64::NAN, So3::identity()).is_err());
        assert!(shape.at_coefficient(0.0, So3::identity()).is_err());
    }

    #[test]
    fn material_shape_drives_native_finite_translation_spin_recoil_and_positive_heat() {
        use fs_mbd::spatial_loads::internal_impulse::InternalWrenchImpulseControls;
        use fs_mbd::spatial_loads::internal_impulse::viscous::{
            PreparedViscousConnector, ViscousConnectorControls,
        };
        use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3 as RigidVec3};

        // The physical material axes are fixed in world space. Only their
        // numerical principal-frame labels change. Unequal resistance in the
        // other directions makes an incorrect frame transform observable.
        let diagonal = [1.0, 2.0, 3.0, 1.0, 4.0, 5.0];
        let shape = SpatialResistanceShape::new(
            core::array::from_fn(|i| {
                core::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 })
            }),
            0.5,
        )
        .unwrap();
        let properties = [MassProperties::new(
            2.0,
            RigidVec3::ZERO,
            RigidVec3::new(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
        )
        .unwrap(); 2];
        let controls = ViscousConnectorControls {
            impulse: InternalWrenchImpulseControls {
                angular_origin_world_m: RigidVec3::ZERO,
                linear_tolerance_ns: 2e-11,
                angular_tolerance_nms: 2e-11,
                work_tolerance_j: 2e-11,
            },
            linear_speed_tolerance_m_s: 2e-11,
            angular_speed_tolerance_rad_s: 2e-11,
        };
        let close = |actual: f64, expected: f64| {
            assert!((actual - expected).abs() < 2e-11, "{actual} != {expected}");
        };
        for rotation in [
            So3::identity(),
            So3::exp(fs_ga::So3Tangent::new(Vec3::new(0.27, -0.19, 0.31))).unwrap(),
        ] {
            let q = rotation.as_quat();
            let orientation = UnitQuaternion::new(q.w, q.x, q.y, q.z).unwrap();
            let states = [1.0, -1.0].map(|sign| {
                RigidBodyState::new(
                    Pose::new(RigidVec3::new(-sign, 0.0, 0.0), orientation).unwrap(),
                    RigidVec3::new(2.0 * sign, 0.0, 0.0),
                    orientation.rotate_world_to_body(RigidVec3::new(sign / 3.0, 0.0, 0.0)),
                )
                .unwrap()
            });
            for coefficient in [1.0, 3.0] {
                let prepared = PreparedViscousConnector::new(
                    states.map(|state| state.pose()),
                    properties,
                    [RigidVec3::ZERO; 2],
                    shape.at_coefficient(coefficient, rotation).unwrap(),
                    0.1,
                    || false,
                )
                .unwrap();
                let receipt = prepared.relax(states, controls, || false).unwrap();
                // Independent finite-body balances: v'=-c*v and
                // omega'=-(2*c*ell^2/I)*omega=-1.5*c*omega.
                let linear = f64::exp(-0.1 * coefficient);
                let spin = f64::exp(-0.15 * coefficient);
                let heat =
                    -2.0 * f64::exp_m1(-0.2 * coefficient) - f64::exp_m1(-0.3 * coefficient) / 3.0;
                for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
                    let actual = receipt.impulse.bodies[i].state_after;
                    let momentum = actual.linear_momentum_world();
                    let angular = actual
                        .pose()
                        .orientation()
                        .rotate_body_to_world(actual.angular_momentum_body());
                    close(momentum.x, sign * 2.0 * linear);
                    close(momentum.y, 0.0);
                    close(momentum.z, 0.0);
                    close(angular.x, sign * spin / 3.0);
                    close(angular.y, 0.0);
                    close(angular.z, 0.0);
                    assert_eq!(actual.pose(), states[i].pose());
                }
                assert!(receipt.dissipated_energy_j > 0.0);
                close(receipt.dissipated_energy_j, heat);
                close(receipt.impulse.kinetic_energy_change_j, -heat);
                close(receipt.kinetic_heat_residual_j, 0.0);
                close(receipt.connector_heat_residual_j, 0.0);
            }
        }
    }
}
