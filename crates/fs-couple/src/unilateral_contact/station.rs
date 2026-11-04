//! Existing distributed contact laws at physical scalar force/velocity ports.
//! Station measures remain metres; no quadrature normalization or modal basis
//! is inferred. The caller owns its kinematics and accepted time discretization.
use super::{
    DContactError, Obstacle, UnilateralObstacle, friction_error, friction_law,
    validate_span_stations,
};
use fs_dcontact::AffineContactStep;
use fs_tribo::FrictionLaw;

pub(crate) struct StationContact {
    station_frac: f64,
    normal: Obstacle,
    friction: Option<(FrictionLaw, f64)>,
}

impl StationContact {
    /// Preserve original obstacle/station order, gaps, measures and provenance.
    pub(crate) fn prepare(specs: &[UnilateralObstacle]) -> Result<Vec<Self>, DContactError> {
        let mut stations = Vec::new();
        for spec in specs {
            validate_span_stations(spec)?;
            let friction = friction_law(spec)?;
            for ((&station_frac, &gap), &measure) in
                spec.stations.iter().zip(&spec.gaps_m).zip(&spec.measures_m)
            {
                let normal = Obstacle::new(
                    vec![1.0],
                    1,
                    1,
                    vec![gap],
                    vec![measure],
                    spec.stiffness,
                    spec.alpha,
                    spec.provenance.clone(),
                )?
                .with_internal_loss(spec.internal_loss)?;
                stations.push(Self {
                    station_frac,
                    normal,
                    friction: friction.clone(),
                });
            }
        }
        Ok(stations)
    }

    pub(crate) fn station_frac(&self) -> f64 {
        self.station_frac
    }

    /// Instantaneous force, in N, opposing positive penetration. Coincident
    /// endpoints select the original potential derivative, while the actual
    /// supplied port velocity still drives Hunt--Crossley loss. This evaluation
    /// does not claim finite-step work closure for a held-force time owner.
    /// Friction uses the same elastic normal capacity and explicit scalar ramp
    /// as the fixed-end modal adapter, with no inferred tangential frame.
    pub(crate) fn force(
        &self,
        displacement_m: f64,
        velocity_m_s: f64,
    ) -> Result<f64, DContactError> {
        let response = AffineContactStep::new(&self.normal, displacement_m, 1)?
            .response(displacement_m, velocity_m_s)?;
        let mut force = response.force;
        if let Some((law, regularization)) = &self.friction {
            force += law
                .regularized_traction_1d(velocity_m_s, -response.elastic_force, *regularization)
                .map_err(friction_error)?;
        }
        if !force.is_finite() {
            return Err(DContactError::Parameter {
                what: "station contact force overflow",
            });
        }
        Ok(force)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec() -> UnilateralObstacle {
        UnilateralObstacle {
            stations: vec![0.35],
            gaps_m: vec![0.001],
            measures_m: vec![0.25],
            stiffness: 1e6,
            alpha: 2.0,
            mu_kinetic: 0.3,
            mu_static: 0.4,
            friction_regularization_m_s: Some(0.01),
            internal_loss: 2.0,
            provenance: "physical line contact test".into(),
        }
    }

    #[test]
    fn forces_retain_line_units_port_velocity_and_the_existing_modal_laws() {
        let spec = spec();
        let stations = StationContact::prepare(&[spec.clone()]).unwrap();
        let (q, v) = (0.003, 0.02);
        // 0.25 m * 1e6 N/m^3 * (0.002 m)^2 = 1 N. The actual
        // port velocity adds 0.04 N of normal loss and 0.3 N of friction.
        assert!((stations[0].force(q, v).unwrap() + 1.34).abs() < 1e-13);
        let string = fs_scenario::PrestressedString {
            length_m: 0.8,
            tension_n: 20.0,
            lin_density_kg_m: 0.005,
            axial_stiffness_n: 0.0,
            width_m: 0.001,
            n_modes: 1,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
            tension_loss_factor: None,
            bending_loss_factor: None,
        };
        let phi = fs_math::det::sqrt(2.0 / (string.lin_density_kg_m * string.length_m))
            * fs_math::det::sin(core::f64::consts::PI * spec.stations[0]);
        for velocity in [-1.0, -0.002, 0.0, v] {
            let modal = super::super::modal_contact_forces(
                &string,
                &[spec.clone()],
                &[q / phi, velocity / phi],
            )
            .unwrap()[0]
                / phi;
            let actual = stations[0].force(q, velocity).unwrap();
            assert!((actual - modal).abs() < 1e-12 * actual.abs().max(1.0));
        }
    }

    #[test]
    fn quadrature_partition_preserves_force_and_normal_unloading_never_attracts() {
        let mut original = spec();
        let whole = StationContact::prepare(&[original.clone()]).unwrap();
        let mut split = original.clone();
        split.stations = vec![original.stations[0]; 2];
        split.gaps_m = vec![original.gaps_m[0]; 2];
        split.measures_m = vec![original.measures_m[0] / 2.0; 2];
        let split = StationContact::prepare(&[split]).unwrap();
        for velocity in [-1.0, -0.002, 0.0, 0.02] {
            assert_eq!(
                whole[0].force(0.003, velocity).unwrap(),
                split
                    .iter()
                    .map(|s| s.force(0.003, velocity).unwrap())
                    .sum::<f64>()
            );
        }
        original.mu_kinetic = 0.0;
        let normal = StationContact::prepare(&[original]).unwrap();
        assert_eq!(normal[0].force(0.003, -1.0).unwrap(), 0.0);
        assert!((normal[0].force(0.003, 0.02).unwrap() + 1.04).abs() < 1e-13);
        assert!((normal[0].force(0.003, 0.0).unwrap() + 1.0).abs() < 1e-13);
    }

    #[test]
    fn malformed_physical_sites_are_refused_before_force_evaluation() {
        let mutations: [fn(&mut UnilateralObstacle); 10] = [
            |s: &mut UnilateralObstacle| s.measures_m.clear(),
            |s: &mut UnilateralObstacle| s.measures_m[0] = 0.0,
            |s: &mut UnilateralObstacle| s.gaps_m.clear(),
            |s: &mut UnilateralObstacle| s.gaps_m[0] = f64::NAN,
            |s: &mut UnilateralObstacle| s.stations[0] = 0.0,
            |s: &mut UnilateralObstacle| s.stations[0] = 1.1,
            |s: &mut UnilateralObstacle| s.stiffness = -1.0,
            |s: &mut UnilateralObstacle| s.alpha = 0.5,
            |s: &mut UnilateralObstacle| s.internal_loss = -0.1,
            |s: &mut UnilateralObstacle| s.friction_regularization_m_s = None,
        ];
        for mutate in mutations {
            let mut invalid = spec();
            mutate(&mut invalid);
            assert!(StationContact::prepare(&[invalid]).is_err());
        }
    }
}
