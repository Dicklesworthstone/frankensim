//! G1/G3: mixed physical source units, coherent fields and fixed source state.
use super::*;
use crate::helmholtz::bc::{ForceDipole, IncidentSource, Monopole, PlaneWave};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};

fn exterior() -> MolecularExterior {
    let data = "declared 1 synthetic 1 vibration 1 100 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap()
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-18),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn prepared_mixed_sources_keep_original_material_phases_units_and_euler_velocity() {
    let exterior = exterior();
    let plane = PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(0.8, -0.3),
    };
    let monopole = Monopole {
        position_m: [-0.2, 0.1, 0.0],
        volume_velocity_m3_s: C64::new(0.001, -0.0003),
    };
    let force = ForceDipole {
        position_m: [0.1, -0.2, 0.3],
        force_on_fluid_n: [
            C64::new(0.3, 0.2),
            C64::new(-0.1, 0.15),
            C64::new(0.05, -0.1),
        ],
    };
    let sources = [
        IncidentSource::PlaneWave(plane),
        IncidentSource::Monopole(monopole),
        IncidentSource::ForceDipole(force),
    ];
    let omega = 1000.0;
    let prepared = exterior
        .prepare_incident(omega, IncidentField::Sources(&sources), || false)
        .unwrap();
    let k = prepared.wavenumber_per_m();
    let momentum = exterior.source().gas().density * omega;
    assert!(k.im > 0.02 * k.re);
    assert_eq!(prepared.angular_frequency_rad_s(), omega);
    assert_eq!(prepared.source_count(), 3);
    assert_eq!(
        prepared.source_positions(),
        [monopole.position_m, force.position_m]
    );
    assert!(Arc::ptr_eq(&prepared.source().source, &exterior.source));
    let rebound = prepared
        .with_field(IncidentField::PlaneWave(plane), || false)
        .unwrap();
    assert_eq!(rebound.wavenumber_per_m(), k);
    assert_eq!(rebound.angular_frequency_rad_s(), omega);
    assert_eq!(rebound.source_count(), 1);
    assert!(rebound.source_positions().is_empty());
    assert!(Arc::ptr_eq(&rebound.source().source, &exterior.source));
    let phase = |distance: f64| {
        C64::new((k.re * distance).cos(), (k.re * distance).sin()).scale((-k.im * distance).exp())
    };
    // Independent analytic pressure, including the monopole's rho*omega and
    // the force source's absence of an additional density factor.
    let pressure = |point: [f64; 3]| {
        let delta = std::array::from_fn::<_, 3, _>(|a| point[a] - monopole.position_m[a]);
        let r = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
        let monopole_pressure = C64::new(0.0, -momentum)
            * monopole.volume_velocity_m3_s
            * phase(r).scale(1.0 / (4.0 * std::f64::consts::PI * r));
        let delta = std::array::from_fn::<_, 3, _>(|a| point[a] - force.position_m[a]);
        let r = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
        let radial_force = (0..3).fold(C64::ZERO, |sum, a| {
            sum + force.force_on_fluid_n[a].scale(delta[a] / r)
        });
        let force_pressure = phase(r).scale(1.0 / (4.0 * std::f64::consts::PI * r))
            * C64::new(1.0 / r + k.im, -k.re)
            * radial_force;
        plane.amplitude_pa * phase(point[2]) + monopole_pressure + force_pressure
    };
    let points = [[0.7, 0.2, 1.3], [-0.3, 0.4, -0.5]];
    let actual = prepared
        .pressure_velocity_at_many(&points, 6, || false)
        .unwrap();
    for (sample, point) in actual.iter().zip(points) {
        close(sample.pressure_pa, pressure(point), 3e-13);
        let h = 1e-4;
        for axis in 0..3 {
            let displaced = |offset: f64| {
                let mut x = point;
                x[axis] += offset;
                pressure(x)
            };
            let gradient = (displaced(-2.0 * h) - displaced(-h).scale(8.0)
                + displaced(h).scale(8.0)
                - displaced(2.0 * h))
            .scale(1.0 / (12.0 * h));
            close(
                sample.velocity_m_s[axis],
                gradient * C64::new(0.0, -1.0 / momentum),
                2e-8,
            );
        }
        let original = exterior
            .incident_pressure_velocity_at(omega, IncidentField::Sources(&sources), point, || false)
            .unwrap();
        assert_eq!(sample.pressure_pa, original.pressure_pa);
        assert_eq!(sample.velocity_m_s, original.velocity_m_s);
    }
    let again = prepared
        .clone()
        .pressure_velocity_at_many(&points, 6, || false)
        .unwrap();
    assert_eq!(actual[0].pressure_pa, again[0].pressure_pa);
    assert_eq!(actual[1].velocity_m_s, again[1].velocity_m_s);
}

#[test]
fn prepared_incident_admits_whole_batches_and_refuses_partial_or_cancelled_fields() {
    let exterior = exterior();
    let monopoles = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::ONE,
    }];
    let field = IncidentField::Monopoles(&monopoles);
    let prepared = exterior.prepare_incident(1000.0, field, || false).unwrap();
    let points = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(
        prepared
            .pressure_velocity_at_many(&points, 1, || false)
            .is_err()
    );
    assert_eq!(
        prepared
            .pressure_velocity_at_many(&points, 2, || false)
            .unwrap()
            .len(),
        2
    );
    assert!(
        prepared
            .pressure_velocity_at_many(&[], 1, || false)
            .unwrap()
            .is_empty()
    );
    assert!(
        prepared
            .pressure_velocity_at_many(&[], 0, || false)
            .is_err()
    );
    assert!(matches!(
        prepared.pressure_velocity_at_many(&points, 2, || true),
        Err(HelmholtzError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        prepared.pressure_velocity_at_many(&points, 2, || {
            polls += 1;
            polls == 5
        }),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(
        prepared
            .pressure_velocity_at_many(&[[1.0, 0.0, 0.0], [0.0; 3]], 2, || false)
            .is_err()
    );
    assert!(
        prepared
            .pressure_velocity_at([f64::NAN, 0.0, 0.0], || false)
            .is_err()
    );
    assert!(matches!(
        prepared.pressure_velocity_at(points[0], || true),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(exterior.prepare_incident(99.0, field, || false).is_err());
    assert!(exterior.prepare_incident(2001.0, field, || false).is_err());
    assert!(
        exterior
            .prepare_incident(f64::NAN, field, || false)
            .is_err()
    );
    assert!(matches!(
        exterior.prepare_incident(1000.0, field, || true),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    assert!(
        exterior
            .prepare_incident(1000.0, IncidentField::Monopoles(&[]), || false)
            .is_err()
    );
    assert!(
        exterior
            .prepare_incident(
                1000.0,
                IncidentField::PlaneWave(PlaneWave {
                    direction: [2.0, 0.0, 0.0],
                    amplitude_pa: C64::ONE
                }),
                || false
            )
            .is_err()
    );
    assert!(
        prepared
            .pressure_velocity_at_many(&points, 2, || false)
            .is_ok()
    );
    assert!(
        prepared
            .with_field(IncidentField::Monopoles(&[]), || false)
            .is_err()
    );
    assert!(matches!(
        prepared.with_field(field, || true),
        Err(HelmholtzError::Cancelled)
    ));
}
