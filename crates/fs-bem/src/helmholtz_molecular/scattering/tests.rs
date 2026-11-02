use super::*;
use crate::helmholtz::bc::{
    BoundaryOperator, PlaneWave,
    incident::{self, ForceDipole, IncidentSource, Monopole},
};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

fn model(relaxing: bool) -> MolecularExterior {
    // Manufactured reservoir, deliberately not a measured atmospheric model.
    let text = if relaxing {
        "declared 1 synthetic 1 vibration 1 100 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap()
}

fn options(surface: &SpherePanels) -> MolecularScatteringOptions {
    let n = surface.areas().len();
    MolecularScatteringOptions {
        boundary: ComplexOptions {
            max_panels: n,
            max_bytes: 32 * n * n,
            max_work: 512 * (n * n) as u64,
        },
        minimum_source_clearance_m: 0.02,
        maximum_source_point_products: 32 * n as u64,
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-20),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn molecular_incident_momentum_and_force_gradient_keep_dispersion_and_near_field() {
    let exterior = model(true);
    let omega = 1000.0;
    let rho = exterior.source().gas().density;
    let k = exterior.source().relaxing_wavenumber(omega).unwrap();
    assert!(k.im > 0.02 * k.re);
    let plane = IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(0.8, -0.3),
    });
    let sample = exterior
        .incident_pressure_velocity_at(omega, plane, [0.0, 0.0, 0.7], || false)
        .unwrap();
    let expected = C64::new(0.8, -0.3)
        * C64::new((k.re * 0.7).cos(), (k.re * 0.7).sin()).scale((-k.im * 0.7).exp());
    close(sample.pressure_pa, expected, 2e-13);
    close(
        sample.velocity_m_s[2],
        expected * k.scale(1.0 / (rho * omega)),
        2e-13,
    );
    let lossless_shortcut = expected.scale(1.0 / (rho * exterior.source().gas().sound_speed));
    assert!(
        (sample.velocity_m_s[2] - lossless_shortcut).abs() > 0.02 * sample.velocity_m_s[2].abs()
    );

    let force = ForceDipole {
        position_m: [0.1, -0.2, 0.0],
        force_on_fluid_n: [
            C64::new(0.3, 0.2),
            C64::new(-0.1, 0.15),
            C64::new(0.05, -0.1),
        ],
    };
    // Independent pressure oracle, then a three-step Cartesian difference
    // ladder for the force-source Hessian and physical Euler velocity.
    let pressure = |point: [f64; 3]| {
        let delta = std::array::from_fn::<_, 3, _>(|a| point[a] - force.position_m[a]);
        let r = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
        let g = C64::new((k.re * r).cos(), (k.re * r).sin())
            .scale((-k.im * r).exp() / (4.0 * core::f64::consts::PI * r));
        let axial = (0..3).fold(C64::ZERO, |sum, a| {
            sum + force.force_on_fluid_n[a].scale(delta[a] / r)
        });
        g * C64::new(1.0 / r + k.im, -k.re) * axial
    };
    let sources = [force];
    let point = [0.7, 0.2, -0.3];
    let field = exterior
        .incident_pressure_velocity_at(omega, IncidentField::ForceDipoles(&sources), point, || {
            false
        })
        .unwrap();
    close(field.pressure_pa, pressure(point), 2e-13);
    let mut errors = Vec::new();
    for h in [0.002, 0.001, 0.0005] {
        let mut error = 0.0_f64;
        for axis in 0..3 {
            let mut plus = point;
            plus[axis] += h;
            let mut minus = point;
            minus[axis] -= h;
            let u = (pressure(plus) - pressure(minus)).scale(0.5 / h)
                * C64::new(0.0, -1.0 / (rho * omega));
            error = error.max((u - field.velocity_m_s[axis]).abs());
        }
        errors.push(error);
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    let force = [ForceDipole {
        position_m: [0.0; 3],
        force_on_fluid_n: [C64::ONE, C64::ZERO, C64::ZERO],
    }];
    let null = exterior
        .incident_pressure_velocity_at(
            omega,
            IncidentField::ForceDipoles(&force),
            [0.0, 0.0, 0.3],
            || false,
        )
        .unwrap();
    assert_eq!(null.pressure_pa, C64::ZERO);
    assert!(null.velocity_m_s[0].abs() > 1e-4);
    assert_eq!(null.velocity_m_s[1], C64::ZERO);
    assert_eq!(null.velocity_m_s[2], C64::ZERO);
}

#[test]
fn empty_molecular_sources_and_scattering_recover_original_lossless_owner() {
    let exterior = model(false);
    let omega = 650.0;
    let surface = SpherePanels::icosphere(0.12, 1).unwrap();
    let beta = vec![C64::new(0.25, 0.13); surface.areas().len()];
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.8, 0.1),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [0.7, 0.0, 0.0],
            volume_velocity_m3_s: C64::new(0.002, -0.0004),
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [-0.5, 0.2, 0.0],
            force_on_fluid_n: [C64::new(0.1, 0.04), C64::ZERO, C64::new(-0.07, 0.01)],
        }),
    ];
    let field = IncidentField::Sources(&sources);
    let k = exterior.source().relaxing_wavenumber(omega).unwrap();
    let medium = Medium {
        density: exterior.source().gas().density,
        sound_speed: omega / k.re,
    };
    let point = [0.3, -0.4, 0.5];
    let original = field.pressure_velocity_at(k.re, medium, point).unwrap();
    let molecular = exterior
        .incident_pressure_velocity_at(omega, field, point, || false)
        .unwrap();
    close(molecular.pressure_pa, original.pressure_pa, 2e-13);
    for axis in 0..3 {
        close(
            molecular.velocity_m_s[axis],
            original.velocity_m_s[axis],
            2e-13,
        );
    }
    for law in [
        Admittance::Rigid,
        Admittance::Local {
            beta: &beta,
            declared_active: false,
        },
    ] {
        for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
            let reference =
                incident::solve_scattering(&surface, k.re, medium, law, field, 0.02, formulation)
                    .unwrap();
            let solved = exterior
                .solve_scattering_batch(
                    &surface,
                    omega,
                    &[field],
                    law,
                    formulation,
                    &options(&surface),
                    || false,
                )
                .unwrap()
                .remove(0);
            for (&actual, &expected) in solved
                .solution()
                .scattered
                .pressure_pa()
                .iter()
                .zip(&reference.solution().scattered.pressure)
            {
                close(actual, expected, 3e-11);
            }
            let expected = reference.solution().absorbed_power_w;
            assert!(
                (solved.solution().absorbed_power_w - expected).abs()
                    <= 3e-11 * expected.abs().max(1e-20)
            );
        }
    }
}

struct MeanWall {
    areas: Vec<f64>,
    beta: C64,
}
impl BoundaryOperator for MeanWall {
    fn len(&self) -> usize {
        self.areas.len()
    }
    fn fingerprint(&self) -> u64 {
        0x4d4f4c5343415454
    }
    fn apply(&self, pressure: &[C64], output: &mut [C64]) -> Result<(), HelmholtzError> {
        let area: f64 = self.areas.iter().sum();
        let mean = pressure
            .iter()
            .zip(&self.areas)
            .fold(C64::ZERO, |sum, (&p, &a)| sum + p.scale(a / area));
        output.fill(self.beta * mean);
        Ok(())
    }
}

#[test]
fn molecular_nonlocal_scattering_preserves_total_wall_work_and_source_interference() {
    let exterior = model(true);
    let omega = 1000.0;
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let wall = MeanWall {
        areas: surface.areas().to_vec(),
        beta: C64::new(0.3, 0.12),
    };
    let law = Admittance::NonLocal {
        operator: &wall,
        declared_active: false,
    };
    let first = Monopole {
        position_m: [0.5, 0.0, 0.1],
        volume_velocity_m3_s: C64::new(0.001, 0.0002),
    };
    let second = Monopole {
        position_m: [-0.4, 0.2, 0.3],
        volume_velocity_m3_s: C64::new(-0.0003, 0.0008),
    };
    let a = [first];
    let b = [second];
    let both = [first, second];
    let fields = [
        IncidentField::Monopoles(&a),
        IncidentField::Monopoles(&b),
        IncidentField::Monopoles(&both),
    ];
    let solved = exterior
        .solve_scattering_batch(
            &surface,
            omega,
            &fields,
            law,
            Formulation::PlainCbie,
            &options(&surface),
            || false,
        )
        .unwrap();
    let total = solved[2].solution();
    let p = total.total_pressure_pa();
    let velocity = total.total_normal_velocity_m_s();
    let area: f64 = surface.areas().iter().sum();
    let mean = p
        .iter()
        .zip(surface.areas())
        .fold(C64::ZERO, |sum, (&p, &a)| sum + p.scale(a / area));
    let rho_c = exterior.source().gas().density * exterior.source().gas().sound_speed;
    let expected_velocity = -(wall.beta * mean).scale(1.0 / rho_c);
    for &v in &velocity {
        close(v, expected_velocity, 2e-11);
    }
    let absorbed = 0.5 * area * wall.beta.re * mean.norm_sq() / rho_c;
    assert!((total.absorbed_power_w - absorbed).abs() < 2e-11 * absorbed);
    assert!((total.total_boundary_work_w + absorbed).abs() < 2e-11 * absorbed);
    let independent = solved[0].solution().absorbed_power_w + solved[1].solution().absorbed_power_w;
    assert!(
        (independent - absorbed).abs() > 0.01 * absorbed,
        "coherent wall cross term disappeared"
    );
    for i in 0..surface.areas().len() {
        close(
            total.scattered.pressure_pa()[i],
            solved[0].solution().scattered.pressure_pa()[i]
                + solved[1].solution().scattered.pressure_pa()[i],
            2e-11,
        );
    }
    let points = [[0.1, -0.4, 0.4], [0.0, 0.0, 0.7]];
    let a = solved[0]
        .total_pressure_velocity_at(&points, 0.02, || false)
        .unwrap();
    let b = solved[1]
        .total_pressure_velocity_at(&points, 0.02, || false)
        .unwrap();
    let both = solved[2]
        .total_pressure_velocity_at(&points, 0.02, || false)
        .unwrap();
    for i in 0..points.len() {
        close(
            both[i].pressure_pa,
            a[i].pressure_pa + b[i].pressure_pa,
            3e-11,
        );
        for axis in 0..3 {
            close(
                both[i].velocity_m_s[axis],
                a[i].velocity_m_s[axis] + b[i].velocity_m_s[axis],
                3e-11,
            );
        }
    }
}

#[test]
fn molecular_scattering_refuses_singular_geometry_work_and_cancelled_queries() {
    let exterior = model(true);
    let omega = 800.0;
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let source = [Monopole {
        position_m: [0.5, 0.0, 0.0],
        volume_velocity_m3_s: C64::ONE,
    }];
    let fields = [IncidentField::Monopoles(&source)];
    let mut cap = options(&surface);
    cap.boundary.max_panels -= 1;
    assert!(matches!(
        exterior.solve_scattering_batch(
            &surface,
            omega,
            &fields,
            Admittance::Rigid,
            Formulation::PlainCbie,
            &cap,
            || false
        ),
        Err(MolecularError::Operator(ComplexError::PanelBudget { .. }))
    ));
    for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
        let mut cap = options(&surface);
        cap.boundary.max_bytes -= 1;
        assert!(matches!(
            exterior.solve_scattering_batch(
                &surface,
                omega,
                &fields,
                Admittance::Rigid,
                formulation,
                &cap,
                || false
            ),
            Err(MolecularError::Operator(ComplexError::Budget { .. }))
        ));
        let mut cap = options(&surface);
        let units = if formulation == Formulation::PlainCbie {
            512
        } else {
            8
        };
        cap.boundary.max_work = units * (surface.areas().len().pow(2) as u64) - 1;
        assert!(matches!(
            exterior.solve_scattering_batch(
                &surface,
                omega,
                &fields,
                Admittance::Rigid,
                formulation,
                &cap,
                || false
            ),
            Err(MolecularError::Operator(ComplexError::Budget { .. }))
        ));
    }
    let mut opts = options(&surface);
    opts.maximum_source_point_products = surface.areas().len() as u64 - 1;
    assert!(
        exterior
            .solve_scattering_batch(
                &surface,
                omega,
                &fields,
                Admittance::Rigid,
                Formulation::PlainCbie,
                &opts,
                || false
            )
            .is_err()
    );
    let opts = options(&surface);
    for position in [[0.0; 3], surface.centroids()[0]] {
        let sources = [Monopole {
            position_m: position,
            ..source[0]
        }];
        assert!(
            exterior
                .solve_scattering_batch(
                    &surface,
                    omega,
                    &[IncidentField::Monopoles(&sources)],
                    Admittance::Rigid,
                    Formulation::PlainCbie,
                    &opts,
                    || false
                )
                .is_err()
        );
    }
    for frequency in [0.0, 99.0, 2001.0, f64::NAN] {
        assert!(
            exterior
                .solve_scattering_batch(
                    &surface,
                    frequency,
                    &fields,
                    Admittance::Rigid,
                    Formulation::PlainCbie,
                    &opts,
                    || false
                )
                .is_err()
        );
    }
    let mut calls = 0;
    assert!(matches!(
        exterior.solve_scattering_batch(
            &surface,
            omega,
            &fields,
            Admittance::Rigid,
            Formulation::PlainCbie,
            &opts,
            || {
                calls += 1;
                calls > 4
            }
        ),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    let solved = exterior
        .solve_scattering_batch(
            &surface,
            omega,
            &fields,
            Admittance::Rigid,
            Formulation::PlainCbie,
            &opts,
            || false,
        )
        .unwrap()
        .remove(0);
    assert!(
        solved
            .total_pressure_velocity_at(&[source[0].position_m], 0.02, || false)
            .is_err()
    );
    assert!(
        solved
            .total_pressure_velocity_at(&[[0.0; 3]], 0.02, || false)
            .is_err()
    );
    assert!(
        solved
            .total_pressure_velocity_at(&[[0.0, 0.0, 0.4]], 0.02, || true)
            .is_err()
    );
}

#[test]
fn dispersive_source_sum_keeps_a_weak_field_between_large_opposing_waves() {
    let exterior = model(true);
    let term = |amplitude| {
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::from_re(amplitude),
        })
    };
    let sources = [term(1e16), term(1.0), term(-1e16)];
    let weak = [term(1.0)];
    let point = [0.2, 0.3, 0.7];
    let combined = exterior
        .incident_pressure_velocity_at(900.0, IncidentField::Sources(&sources), point, || false)
        .unwrap();
    let expected = exterior
        .incident_pressure_velocity_at(900.0, IncidentField::Sources(&weak), point, || false)
        .unwrap();
    close(combined.pressure_pa, expected.pressure_pa, 1e-13);
    close(combined.velocity_m_s[2], expected.velocity_m_s[2], 1e-13);
}
