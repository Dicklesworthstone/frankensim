use super::*;
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

fn source(relaxing: bool, time_scale: f64) -> RelaxationState {
    // Synthetic slow internal reservoir isolates constitutive causality. Its
    // constants are deliberately not offered as a measured dry-air source.
    let data = if relaxing {
        format!(
            "declared 1 synthetic 1 vibration 1 900 1 power {} 293.15 0 200 500",
            101325.0 * 0.001 * time_scale
        )
    } else {
        "declared 1 synthetic 0".to_string()
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    law.evaluate(&RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap())
        .unwrap()
}
fn model(relaxing: bool, time_scale: f64) -> MolecularExterior {
    MolecularExterior::new(source(relaxing, time_scale), [100.0, 2000.0]).unwrap()
}
fn limits(n: usize) -> ComplexOptions {
    ComplexOptions {
        max_panels: n,
        max_bytes: 32 * n * n,
        max_work: 512 * (n * n) as u64,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn values_only_owner_retains_complex_operator_bits_with_half_matrix_storage() {
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
        for k in [C64::from_re(1.5), C64::new(1.5, 0.3)] {
            let value =
                complex::assemble_values_with_cancel(&surface, k, formulation, &limits(n), || {
                    false
                })
                .unwrap();
            let mut original_limits = limits(n);
            original_limits.max_bytes *= 2;
            let derivative =
                complex::assemble_with_cancel(&surface, k, formulation, &original_limits, || false)
                    .unwrap();
            for (a, b) in value
                .a
                .iter()
                .chain(&value.b)
                .zip(derivative.a.iter().chain(&derivative.b))
            {
                assert_eq!(
                    (a.re.to_bits(), a.im.to_bits()),
                    (b.re.to_bits(), b.im.to_bits())
                );
            }
            let mut insufficient = limits(n);
            insufficient.max_bytes -= 1;
            assert!(matches!(
                complex::assemble_values_with_cancel(
                    &surface,
                    k,
                    formulation,
                    &insufficient,
                    || false
                ),
                Err(ComplexError::Budget { .. })
            ));
        }
    }
}

#[test]
fn relaxing_sphere_reaction_and_receiver_converge_to_same_material_oracle() {
    let material = model(true, 1.0);
    let omega = 600.0;
    let radius = 0.2;
    let k = material.source().relaxing_wavenumber(omega).unwrap();
    assert!(k.im > 0.0);
    let q = C64::new(0.0, material.source().gas().density * omega);
    let z = q * C64::new(-1.0 / radius - k.im, k.re).recip();
    let distance = 0.8;
    let phase = k.scale(distance - radius);
    let expected_p = z * C64::new(det::cos(phase.re), det::sin(phase.re))
        .scale(det::exp(-phase.im) * radius / distance);
    let expected_gradient = expected_p * C64::new(-1.0 / distance - k.im, k.re);
    let mut errors = Vec::new();
    let mut observed_errors = Vec::new();
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(radius, level).unwrap();
        let n = surface.areas().len();
        let velocity = vec![C64::ONE; n];
        let solution = material
            .solve_batch(
                &surface,
                omega,
                &[&velocity],
                Formulation::BurtonMiller,
                &limits(n),
                || false,
            )
            .unwrap()
            .remove(0);
        let area: f64 = surface.areas().iter().sum();
        let mean = solution
            .pressure_pa()
            .iter()
            .zip(surface.areas())
            .fold(C64::ZERO, |s, (p, a)| s + p.scale(*a / area));
        errors.push((mean - z).abs() / z.abs());
        let observed = solution
            .observe(&surface, &[[0.0, 0.0, distance]], 0.1, || false)
            .unwrap()[0];
        observed_errors.push((observed.pressure_pa - expected_p).abs() / expected_p.abs());
        assert!(solution.boundary_work_w() > 0.0);
        assert_eq!(solution.source(), material.source());
        if level == 2 {
            close(observed.gradient_pa_per_m[2], expected_gradient, 0.08);
            assert!(observed.gradient_pa_per_m[0].abs() < 1e-10 * expected_gradient.abs());
            assert!(observed.gradient_pa_per_m[1].abs() < 1e-10 * expected_gradient.abs());
        }
    }
    println!("molecular sphere impedance error {errors:?}; receiver error {observed_errors:?}");
    assert!(
        errors[1] < errors[0] * 0.7 && errors[2] < errors[1] * 0.7,
        "{errors:?}"
    );
    assert!(
        observed_errors[1] < observed_errors[0] * 0.7
            && observed_errors[2] < observed_errors[1] * 0.7,
        "{observed_errors:?}"
    );
    assert!(errors[2] < 0.05 && observed_errors[2] < 0.08);
}

#[test]
fn empty_reservoirs_recover_original_radiation_and_rates_change_reaction() {
    let surface = SpherePanels::icosphere(0.2, 1).unwrap();
    let n = surface.areas().len();
    let velocity = vec![C64::ONE; n];
    let omega = 700.0;
    let empty = model(false, 1.0);
    let k = empty.source().relaxing_wavenumber(omega).unwrap();
    assert_eq!(k.im, 0.0);
    let medium = Medium {
        density: empty.source().gas().density,
        sound_speed: omega / k.re,
    };
    for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
        let original = solve_radiation(&surface, k.re, medium, &velocity, formulation).unwrap();
        let solved = empty
            .solve_batch(
                &surface,
                omega,
                &[&velocity],
                formulation,
                &limits(n),
                || false,
            )
            .unwrap()
            .remove(0);
        for (&a, &b) in original.pressure.iter().zip(solved.pressure_pa()) {
            close(a, b, 1e-12);
        }
    }
    let a = model(true, 0.1)
        .solve_batch(
            &surface,
            omega,
            &[&velocity],
            Formulation::BurtonMiller,
            &limits(n),
            || false,
        )
        .unwrap()
        .remove(0);
    let b = model(true, 10.0)
        .solve_batch(
            &surface,
            omega,
            &[&velocity],
            Formulation::BurtonMiller,
            &limits(n),
            || false,
        )
        .unwrap()
        .remove(0);
    assert!((a.boundary_work_w() - b.boundary_work_w()).abs() > 1e-3 * a.boundary_work_w().abs());
    let pa = a
        .observe(&surface, &[[0.0, 0.0, 3.0]], 0.1, || false)
        .unwrap()[0]
        .pressure_pa;
    let pb = b
        .observe(&surface, &[[0.0, 0.0, 3.0]], 0.1, || false)
        .unwrap()[0]
        .pressure_pa;
    assert!((pa - pb).abs() > 1e-3 * pa.abs());
}

#[test]
fn molecular_aperture_matrix_closes_existing_duct_and_observes_solved_flows() {
    use fs_duct::compose::ComposeOptions;
    use fs_duct::molecular::MolecularDuct;
    use fs_duct::network::CoupledTermination;
    use fs_duct::{ConeSlices, Duct, HoleState, Segment};
    let material = model(true, 1.0);
    let surface = SpherePanels::icosphere(0.2, 1).unwrap();
    let n = surface.areas().len();
    // Two declared hemispherical generalized radiators exercise a coupled
    // plane-port fixture. They are not a measured instrument mouth geometry.
    let apertures: Vec<_> = [-1.0, 1.0]
        .into_iter()
        .map(|side| Aperture {
            shape: surface
                .centroids()
                .iter()
                .map(|p| if side * p[2] > 0.0 { 1.0 } else { 0.0 })
                .collect(),
            reference_point: [0.0, 0.0, side * 0.2],
            normal: [0.0, 0.0, side],
        })
        .collect();
    let omega = 700.0;
    let exterior = material
        .aperture_response(
            &surface,
            &apertures,
            omega,
            Formulation::BurtonMiller,
            &limits(n),
            || false,
        )
        .unwrap();
    let matrix = exterior.impedance().acoustic(0).unwrap();
    assert!(matrix[1].abs() > 0.01 * matrix[0].abs());
    assert!(exterior.impedance().rows[0].passivity_margin > 0.0);
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.05,
                length: 0.2,
            },
            Segment::ToneHole {
                hole_radius: 0.025,
                chimney_height: 0.01,
                bore_radius: 0.05,
                state: HoleState::Open,
            },
            Segment::Cylinder {
                radius: 0.05,
                length: 0.2,
            },
        ],
    };
    let duct = MolecularDuct::new(
        &geometry,
        material.source().clone(),
        None,
        ConeSlices::Fixed(1),
        5,
    )
    .unwrap();
    let labels = vec!["hole".to_string(), "bell".to_string()];
    let response = duct
        .prepare_exterior_ports(omega, &labels, ComposeOptions::default())
        .unwrap()
        .response(&CoupledTermination {
            omega,
            labels,
            z: matrix.clone(),
            source_id: format!("molecular-bem:{}", exterior.impedance().declaration),
        })
        .unwrap();
    let flows = &response.exterior_flows_over_u_in;
    let mut power = C64::ZERO;
    for a in 0..2 {
        for b in 0..2 {
            power = power + flows[a].conj() * matrix[2 * a + b] * flows[b];
        }
    }
    assert!(
        (0.5 * power.re - response.exterior_power_per_squared_inlet_flow).abs()
            < 1e-11 * power.re.abs()
    );
    let points = [[0.0, 0.0, 0.9], [1.0, 0.0, 0.0]];
    let observed = exterior
        .observe_volume_flows(&surface, flows, &points, 0.1, || false)
        .unwrap();
    let mut direct = vec![
        PointField {
            pressure_pa: C64::ZERO,
            gradient_pa_per_m: [C64::ZERO; 3]
        };
        2
    ];
    for (i, solution) in exterior.boundary_solutions().iter().enumerate() {
        let amplitude = flows[i].scale(1.0 / exterior.impedance().volume_gains[i]);
        for (out, value) in direct
            .iter_mut()
            .zip(solution.observe(&surface, &points, 0.1, || false).unwrap())
        {
            out.pressure_pa = out.pressure_pa + value.pressure_pa * amplitude;
            for axis in 0..3 {
                out.gradient_pa_per_m[axis] =
                    out.gradient_pa_per_m[axis] + value.gradient_pa_per_m[axis] * amplitude;
            }
        }
    }
    for (a, b) in observed.iter().zip(&direct) {
        close(a.pressure_pa, b.pressure_pa, 2e-12);
    }
    assert!(observed.iter().all(|p| p.pressure_pa.abs() > 0.0));
}

#[test]
fn source_band_geometry_budgets_and_late_cancellation_refuse_without_results() {
    let material = model(true, 1.0);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let velocity = vec![C64::ONE; n];
    for omega in [0.0, 99.0, 2001.0, f64::NAN] {
        assert!(
            material
                .solve_batch(
                    &surface,
                    omega,
                    &[&velocity],
                    Formulation::BurtonMiller,
                    &limits(n),
                    || false
                )
                .is_err()
        );
    }
    let mut calls = 0;
    let first = material
        .solve_batch(
            &surface,
            600.0,
            &[&velocity],
            Formulation::BurtonMiller,
            &limits(n),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        material.solve_batch(
            &surface,
            600.0,
            &[&velocity],
            Formulation::BurtonMiller,
            &limits(n),
            || {
                count += 1;
                count == calls
            }
        ),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    let replay = material
        .solve_batch(
            &surface,
            600.0,
            &[&velocity],
            Formulation::BurtonMiller,
            &limits(n),
            || false,
        )
        .unwrap();
    assert_eq!(first[0].pressure_pa(), replay[0].pressure_pa());
    let wrong_surface = SpherePanels::icosphere(0.21, 0).unwrap();
    assert!(
        first[0]
            .observe(&wrong_surface, &[[0.0, 0.0, 1.0]], 0.1, || false)
            .is_err()
    );
    assert!(
        first[0]
            .observe(&surface, &[[0.0, 0.0, 0.0]], 0.01, || false)
            .is_err()
    );
    let mut tight = limits(n);
    tight.max_bytes -= 1;
    assert!(matches!(
        material.solve_batch(
            &surface,
            600.0,
            &[&velocity],
            Formulation::BurtonMiller,
            &tight,
            || false
        ),
        Err(MolecularError::Operator(ComplexError::Budget { .. }))
    ));
}
