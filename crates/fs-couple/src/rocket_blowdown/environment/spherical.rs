//! Spherical altitude and co-rotating air in one explicit inertial frame.
//!
//! The cached USSA1976 profile supplies local gas properties at the ACTUAL
//! queried point. Its hydrostatic reference gravity remains the standard's
//! convention; the independent supplied central gravitational parameter sets
//! vehicle acceleration and specific potential. No planetary radius or spin
//! is inferred, and no Coriolis or centrifugal force is added to inertial motion.
//!
//! Atmosphere velocity is Omega cross (position-center) plus the declared
//! inertial background wind. This is a spherical reference profile with rigid
//! bulk rotation, not latitude-dependent weather, terrain, oblateness or a
//! self-consistent rotating hydrostatic atmosphere. The source ends at 80 km;
//! unsupported air is refused rather than replaced with vacuum.

use super::{CentralGravity, FlightAirSample, Result, invalid, norm};
use fs_material::gas::atmosphere::StandardAtmosphere1976;
use fs_mbd::Vec3;

/// Local geometric altitude and the original coherent air/wind snapshot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalAirSample {
    /// Radial distance minus the caller's reference radius [m].
    pub geometric_altitude_m: f64,
    /// The same local gas state supplies drag and nozzle ambient pressure.
    pub air: FlightAirSample,
}

/// Cached standard atmosphere, rigid bulk rotation and independent central gravity.
/// All coordinates and velocities are in the identified inertial frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalStandardAtmosphere {
    inertial_frame: u64,
    center_world_m: Vec3,
    reference_radius_m: f64,
    angular_velocity_world_rad_s: Vec3,
    background_wind_world_m_s: Vec3,
    profile: StandardAtmosphere1976,
    gravity: CentralGravity,
}

impl SphericalStandardAtmosphere {
    /// Bind explicit spherical geometry, central mu, rigid atmosphere rotation
    /// and background wind. Angular velocity is an inertial-frame vector.
    pub fn new(
        inertial_frame: u64,
        center_world_m: Vec3,
        reference_radius_m: f64,
        gravitational_parameter_m3_s2: f64,
        angular_velocity_world_rad_s: Vec3,
        background_wind_world_m_s: Vec3,
    ) -> Result<Self> {
        let gravity = CentralGravity::new(
            inertial_frame,
            center_world_m,
            gravitational_parameter_m3_s2,
        )?;
        if !reference_radius_m.is_finite()
            || reference_radius_m <= 0.0
            || !angular_velocity_world_rad_s.is_finite()
            || !background_wind_world_m_s.is_finite()
        {
            return Err(invalid(
                "spherical atmosphere needs a positive finite radius and finite rotation/wind",
            ));
        }
        Ok(Self {
            inertial_frame,
            center_world_m,
            reference_radius_m,
            angular_velocity_world_rad_s,
            background_wind_world_m_s,
            profile: StandardAtmosphere1976::new(),
            gravity,
        })
    }

    /// Caller-owned identity of the common inertial frame.
    #[must_use]
    pub const fn inertial_frame(&self) -> u64 {
        self.inertial_frame
    }
    /// Fixed center of the spherical atmosphere and central gravity source [m].
    #[must_use]
    pub const fn center_world_m(&self) -> Vec3 {
        self.center_world_m
    }
    /// Caller-declared zero-altitude spherical radius [m].
    #[must_use]
    pub const fn reference_radius_m(&self) -> f64 {
        self.reference_radius_m
    }
    /// Rigid atmosphere angular velocity in the declared inertial frame [rad/s].
    #[must_use]
    pub const fn angular_velocity_world_rad_s(&self) -> Vec3 {
        self.angular_velocity_world_rad_s
    }
    /// Explicit velocity added to the rotating atmosphere, in the same frame [m/s].
    #[must_use]
    pub const fn background_wind_world_m_s(&self) -> Vec3 {
        self.background_wind_world_m_s
    }

    /// Sample atmosphere at an actual inertial position, e.g. a nozzle exit or
    /// aerodynamic application point. The original profile enforces geometric
    /// altitude 0..80,000 m; no source extension, surface clipping or vacuum
    /// substitution is performed. Distinct application points need distinct
    /// samples, including their distinct rigid-rotation velocities.
    pub fn sample_at(&self, frame: u64, position_world_m: Vec3) -> Result<SphericalAirSample> {
        if frame != self.inertial_frame {
            return Err(invalid(
                "spherical atmosphere and query must share the declared inertial frame",
            ));
        }
        let displacement = position_world_m.sub(self.center_world_m);
        let radius = norm(displacement);
        if !displacement.is_finite() || !radius.is_finite() || radius <= 0.0 {
            return Err(invalid(
                "spherical atmosphere query requires finite position away from its center",
            ));
        }
        let geometric_altitude_m = radius - self.reference_radius_m;
        let state = self.profile.at_geometric_altitude(geometric_altitude_m)?;
        let wind = self
            .angular_velocity_world_rad_s
            .cross(displacement)
            .add(self.background_wind_world_m_s);
        let air = FlightAirSample::new(self.inertial_frame, state.gas, wind)?;
        Ok(SphericalAirSample {
            geometric_altitude_m,
            air,
        })
    }

    /// Independent central acceleration [m/s²] and specific potential [J/kg]
    /// relative to infinity. This gravity law remains valid outside the air
    /// profile's height domain; it supplies no atmosphere there. No fictitious
    /// rotating-frame forces are included.
    pub fn gravity_at(&self, frame: u64, position_world_m: Vec3) -> Result<(Vec3, f64)> {
        self.gravity.sample(frame, position_world_m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_mbd::UnitQuaternion;

    const FRAME: u64 = 17;
    const RADIUS: f64 = 6_371_000.0;
    const MU: f64 = 3.986_004_418e14;
    const SPIN: f64 = 7.292_115e-5;

    fn environment(center: Vec3, omega: Vec3, wind: Vec3) -> SphericalStandardAtmosphere {
        SphericalStandardAtmosphere::new(FRAME, center, RADIUS, MU, omega, wind).unwrap()
    }
    fn close(a: f64, b: f64, tolerance: f64) {
        assert!((a - b).abs() <= tolerance, "{a} != {b}");
    }
    fn vector_close(a: Vec3, b: Vec3, tolerance: f64) {
        assert!(norm(a.sub(b)) <= tolerance, "{a:?} != {b:?}");
    }

    #[test]
    fn corotating_ground_and_rotating_body_point_have_zero_relative_airflow() {
        let center = Vec3::new(100.0, -200.0, 300.0);
        let omega = Vec3::new(0.0, 0.0, SPIN);
        let model = environment(center, omega, Vec3::ZERO);
        let body_position = center.add(Vec3::new(RADIUS, 0.0, 0.0));
        let body_velocity = omega.cross(body_position.sub(center));
        let ground = model.sample_at(FRAME, body_position).unwrap();
        assert_eq!(ground.geometric_altitude_m, 0.0);
        let stationary = ground
            .air
            .at_body_point(
                FRAME,
                body_velocity,
                UnitQuaternion::IDENTITY,
                omega,
                Vec3::ZERO,
                1.0,
            )
            .unwrap();
        assert_eq!(stationary.speed_m_s(), 0.0);
        assert_eq!(stationary.drag(0.1, 0.5).unwrap().force_world_n, Vec3::ZERO);
        // A point fixed in the rotating body has Omega x arm velocity.
        // Sampling wind at that actual point cancels it; COM wind would not.
        let arm = Vec3::new(0.0, 10.0, 0.0);
        let point = model.sample_at(FRAME, body_position.add(arm)).unwrap();
        let relative = point
            .air
            .at_body_point(
                FRAME,
                body_velocity,
                UnitQuaternion::IDENTITY,
                omega,
                arm,
                1.0,
            )
            .unwrap();
        close(relative.speed_m_s(), 0.0, 1e-12);
        assert!(point.geometric_altitude_m > 0.0);
        assert!(
            ground
                .air
                .at_body_point(
                    FRAME,
                    body_velocity,
                    UnitQuaternion::IDENTITY,
                    omega,
                    arm,
                    1.0,
                )
                .unwrap()
                .speed_m_s()
                > 1e-4
        );
    }

    #[test]
    fn crosswind_height_and_rigid_coordinate_rotation_share_one_physical_sample() {
        let center = Vec3::new(1e6, -2e6, 3e5);
        let omega = Vec3::new(0.0, 0.0, SPIN);
        let wind = Vec3::new(12.0, -3.0, 2.0);
        let model = environment(center, omega, wind);
        let radial = Vec3::new(0.0, RADIUS + 5000.0, 0.0);
        let position = center.add(radial);
        let sample = model.sample_at(FRAME, position).unwrap();
        assert_eq!(sample.geometric_altitude_m, 5000.0);
        assert_eq!(
            *sample.air.gas(),
            StandardAtmosphere1976::new()
                .at_geometric_altitude(5000.0)
                .unwrap()
                .gas
        );
        let velocity = omega.cross(radial);
        let relative = sample
            .air
            .at_body_point(
                FRAME,
                velocity,
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
                Vec3::ZERO,
                1.0,
            )
            .unwrap();
        vector_close(
            relative.relative_velocity_world_m_s(),
            wind.scale(-1.0),
            1e-12,
        );
        let upper = model
            .sample_at(FRAME, center.add(Vec3::new(0.0, RADIUS + 15000.0, 0.0)))
            .unwrap();
        assert!(upper.air.gas().density < sample.air.gas().density);
        let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.72).unwrap();
        let rotate = |v| rotation.rotate_body_to_world(v);
        let rotated = environment(rotate(center), rotate(omega), rotate(wind));
        let rotated_sample = rotated.sample_at(FRAME, rotate(position)).unwrap();
        close(rotated_sample.geometric_altitude_m, 5000.0, 5e-9);
        close(
            rotated_sample.air.pressure_pa() / sample.air.pressure_pa(),
            1.0,
            2e-12,
        );
        vector_close(
            rotated_sample.air.wind_world_m_s(),
            rotate(sample.air.wind_world_m_s()),
            2e-12,
        );
        let load = relative.drag(0.1, 0.4).unwrap();
        let rotated_flow = rotated_sample
            .air
            .at_body_point(
                FRAME,
                rotate(velocity),
                rotation,
                Vec3::ZERO,
                Vec3::ZERO,
                1.0,
            )
            .unwrap();
        vector_close(
            rotated_flow.drag(0.1, 0.4).unwrap().force_world_n,
            rotate(load.force_world_n),
            1e-10,
        );
        let (g, potential) = model.gravity_at(FRAME, position).unwrap();
        let (gr, pr) = rotated.gravity_at(FRAME, rotate(position)).unwrap();
        vector_close(gr, rotate(g), 1e-12);
        close(pr / potential, 1.0, 1e-14);
    }

    #[test]
    fn supplied_inertial_gravity_does_not_rewrite_standard_hydrostatics() {
        let a = environment(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
        let b = SphericalStandardAtmosphere::new(
            FRAME,
            Vec3::ZERO,
            RADIUS,
            2.0 * MU,
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap();
        let position = Vec3::new(RADIUS + 1000.0, 0.0, 0.0);
        assert_eq!(
            a.sample_at(FRAME, position).unwrap(),
            b.sample_at(FRAME, position).unwrap()
        );
        let (ga, pa) = a.gravity_at(FRAME, position).unwrap();
        let (gb, pb) = b.gravity_at(FRAME, position).unwrap();
        assert_eq!(gb, ga.scale(2.0));
        assert_eq!(pb, 2.0 * pa);
        close(ga.x, -MU / position.x / position.x, 1e-14);
        assert_eq!(ga.y, 0.0);
        assert_eq!(ga.z, 0.0);
    }

    #[test]
    fn frame_source_height_and_nonfinite_geometry_are_refused() {
        let model = environment(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
        let ground = Vec3::new(RADIUS, 0.0, 0.0);
        assert!(model.sample_at(FRAME + 1, ground).is_err());
        assert!(model.gravity_at(FRAME + 1, ground).is_err());
        for position in [
            Vec3::ZERO,
            Vec3::new(RADIUS - 0.001, 0.0, 0.0),
            Vec3::new(RADIUS + 80000.001, 0.0, 0.0),
            Vec3::new(f64::NAN, 0.0, 0.0),
            Vec3::new(f64::INFINITY, 0.0, 0.0),
        ] {
            assert!(model.sample_at(FRAME, position).is_err());
        }
        assert!(
            model
                .sample_at(FRAME, Vec3::new(RADIUS + 80000.0, 0.0, 0.0))
                .is_ok()
        );
        assert!(model.gravity_at(FRAME, Vec3::ZERO).is_err());
        assert!(
            model
                .gravity_at(FRAME, Vec3::new(RADIUS + 100000.0, 0.0, 0.0))
                .is_ok()
        );
        for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                SphericalStandardAtmosphere::new(
                    FRAME,
                    Vec3::ZERO,
                    radius,
                    MU,
                    Vec3::ZERO,
                    Vec3::ZERO,
                )
                .is_err()
            );
        }
        assert!(
            SphericalStandardAtmosphere::new(
                FRAME,
                Vec3::ZERO,
                RADIUS,
                MU,
                Vec3::new(f64::NAN, 0.0, 0.0),
                Vec3::ZERO,
            )
            .is_err()
        );
        assert!(
            SphericalStandardAtmosphere::new(
                FRAME,
                Vec3::ZERO,
                RADIUS,
                0.0,
                Vec3::ZERO,
                Vec3::ZERO,
            )
            .is_err()
        );
    }
}
