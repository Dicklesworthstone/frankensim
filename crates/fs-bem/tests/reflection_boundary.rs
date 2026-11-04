//! G0/G1/G3: the original Robin equation across its zero-impedance limit.
use fs_bem::helmholtz::bc::{
    Admittance, IncidentSource, IncidentTraces, PlaneWave, PreparedBoundary,
    PreparedReflectionBoundary,
};
use fs_bem::helmholtz::{Formulation, Medium};
use fs_bem::near_field::{Geometry, Options};
use fs_bem::panel3d::SpherePanels;
use fs_math::{c64::C64, det};

const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 340.0,
};

fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-15),
        "{a:?} != {b:?}"
    );
}

#[test]
fn mixed_reflection_matches_original_finite_robin_drives_and_incidents() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let n = surface.areas().len();
    let beta: Vec<_> = (0..n)
        .map(|i| match i % 3 {
            0 => C64::ZERO,
            1 => C64::new(0.2, -0.8),
            _ => C64::new(0.5, 1.2),
        })
        .collect();
    let reflection: Vec<_> = beta
        .iter()
        .map(|&b| (C64::ONE - b) / (C64::ONE + b))
        .collect();
    let old = PreparedBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        Admittance::Local {
            beta: &beta,
            declared_active: false,
        },
        Formulation::BurtonMiller,
    )
    .unwrap();
    let prepared = PreparedReflectionBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        &reflection,
        false,
        Formulation::BurtonMiller,
    )
    .unwrap();
    let drive: Vec<_> = (0..n)
        .map(|i| C64::new(0.01 * (i + 1) as f64, -0.02))
        .collect();
    let reference = old.solve(&drive).unwrap();
    let actual = prepared.solve(&drive).unwrap();
    for i in 0..n {
        close(
            actual.radiation.pressure[i],
            reference.radiation.pressure[i],
            2e-11,
        );
        close(
            actual.radiation.velocity[i],
            reference.radiation.velocity[i],
            2e-11,
        );
    }
    assert!((actual.absorbed_power_w / reference.absorbed_power_w - 1.0).abs() < 2e-11);
    assert!((actual.driving_power_w / reference.driving_power_w - 1.0).abs() < 2e-11);
    let sources = [IncidentSource::PlaneWave(PlaneWave {
        direction: [0.6, 0.0, 0.8],
        amplitude_pa: C64::new(0.7, -0.2),
    })];
    let incident = IncidentTraces::from_sources(&surface, 3.0, MEDIUM, &sources, 1e-4).unwrap();
    let reference = old.scatter(&incident).unwrap();
    let actual = prepared.scatter(&incident).unwrap();
    for i in 0..n {
        close(
            actual.scattered.pressure[i],
            reference.scattered.pressure[i],
            2e-11,
        );
        close(
            actual.scattered.velocity[i],
            reference.scattered.velocity[i],
            2e-11,
        );
    }
    assert!((actual.absorbed_power_w / reference.absorbed_power_w - 1.0).abs() < 2e-11);
}

#[test]
fn pressure_release_scattering_recovers_the_analytic_spherical_monopole() {
    // A regular radial incident wave j0(kr) excites only the spherical l=0
    // mode. Its exact scattered pressure is -j0(ka)*a/r*exp(ik(r-a)). The
    // polygon/constant-panel approximation is refined independently of r=-1.
    // At ka=0.3, CBIE is below the first interior resonance and retains the
    // established triangle quadrature oracle; BM is exercised above and by
    // the mechanical-resonance observation tests.
    let radius = 0.1;
    let k = 3.0;
    let receiver = [0.0, 0.0, 0.4];
    let (sine, cosine) = det::sin_cos(k * (receiver[2] - radius));
    let expected =
        C64::new(cosine, sine).scale(-det::sin(k * radius) / (k * radius) * radius / receiver[2]);
    let mut errors = Vec::new();
    for refinement in [0, 1, 2] {
        let surface = SpherePanels::icosphere(radius, refinement).unwrap();
        let n = surface.areas().len();
        let mut pressure = Vec::with_capacity(n);
        let mut velocity = Vec::with_capacity(n);
        for (point, normal) in surface.centroids().iter().zip(surface.normals()) {
            let r = det::sqrt(point.iter().map(|v| v * v).sum());
            let kr = k * r;
            let (s, c) = det::sin_cos(kr);
            pressure.push(C64::from_re(s / kr));
            let radial_normal: f64 = point.iter().zip(normal).map(|(x, n)| x * n / r).sum();
            velocity.push(C64::new(
                0.0,
                -(kr * c - s) / (kr * kr) * radial_normal / (MEDIUM.density * MEDIUM.sound_speed),
            ));
        }
        let incident =
            IncidentTraces::from_traces(&surface, k, MEDIUM, &pressure, &velocity).unwrap();
        let prepared = PreparedReflectionBoundary::new(
            &surface,
            k,
            MEDIUM,
            &vec![-C64::ONE; n],
            false,
            Formulation::PlainCbie,
        )
        .unwrap();
        let result = prepared.scatter(&incident).unwrap();
        assert!(result.total_pressure_pa().iter().all(|&p| p == C64::ZERO));
        assert_eq!(result.absorbed_power_w, 0.0);
        assert!(result.scattered.velocity.iter().any(|v| v.abs() > 1e-6));
        let geometry = Geometry::new(&surface, &[receiver], 1e-4).unwrap();
        let received = geometry
            .prepare(k, MEDIUM, Options::default())
            .unwrap()
            .evaluate(&result.scattered)
            .unwrap();
        errors.push((received.pressure[0] - expected).abs() / expected.abs());
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 0.06, "{errors:?}");
}

#[test]
fn rigid_endpoint_retains_original_arithmetic_and_bad_reflections_refuse() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let n = surface.areas().len();
    let drive = vec![C64::new(0.01, 0.02); n];
    let old = PreparedBoundary::new(
        &surface,
        2.0,
        MEDIUM,
        Admittance::Rigid,
        Formulation::BurtonMiller,
    )
    .unwrap()
    .solve(&drive)
    .unwrap();
    let new = PreparedReflectionBoundary::new(
        &surface,
        2.0,
        MEDIUM,
        &vec![C64::ONE; n],
        false,
        Formulation::BurtonMiller,
    )
    .unwrap()
    .solve(&drive)
    .unwrap();
    assert_eq!(old.radiation.pressure, new.radiation.pressure);
    assert_eq!(old.radiation.velocity, new.radiation.velocity);
    for invalid in [C64::from_re(1.1), C64::new(f64::NAN, 0.0)] {
        assert!(
            PreparedReflectionBoundary::new(
                &surface,
                2.0,
                MEDIUM,
                &vec![invalid; n],
                false,
                Formulation::BurtonMiller
            )
            .is_err()
        );
    }
    assert!(
        PreparedReflectionBoundary::new(
            &surface,
            2.0,
            MEDIUM,
            &[C64::ONE],
            false,
            Formulation::BurtonMiller
        )
        .is_err()
    );
}
