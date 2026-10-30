use super::*;

#[test]
fn original_volume_compression_recovers_complex_bulk_work_and_loaded_fields() {
    let size = [0.02; 3];
    let (mesh, points) = box_mesh(1, size);
    let fe = HelmholtzFe::new(&mesh, &points, 1, AIR).unwrap();
    // This explicitly declared constant-pressure Galerkin coordinate verifies
    // compression and work; it is not a full 3-D acoustic accuracy benchmark.
    let basis = HelmholtzModes {
        lambda: vec![0.0],
        frequency_hz: vec![0.0],
        residual: vec![0.0],
        shapes: vec![vec![1.0; fe.ndof()]],
        below_low: 0,
        below_high: 1,
    };
    let model = InteriorMultiport::new(
        &fe,
        &basis,
        &[BoundaryMode::uniform(&end(&fe, 0.0))],
        false,
        &|| false,
    )
    .unwrap();
    let f = 700.0;
    let omega = std::f64::consts::TAU * f;
    let volume = size.iter().product::<f64>();
    let bulk = C64::new(130000.0, 4000.0);
    let reference = 1e6;
    let incident = C64::new(0.3, -0.2);
    let prepared = model
        .prepare_scattering_with_bulk_modulus(f, bulk, &[reference], &|| false)
        .unwrap();
    let response = prepared.response(&[incident], &|| false).unwrap();
    let impedance = C64::new(bulk.im, -bulk.re).scale(1.0 / (omega * volume));
    let flow = incident.scale(2.0 * reference.sqrt()) / (impedance + C64::from_re(reference));
    let pressure = impedance * flow;
    assert!(relative(response.pressure_pa[0], pressure) < 2e-12);
    assert!(relative(response.inward_flow_m3_s[0], flow) < 2e-12);
    let work = prepared
        .bulk_dissipation_w(&response.amplitudes, &|| false)
        .unwrap();
    let expected = 0.5 * bulk.im / (omega * volume) * flow.norm_sq();
    assert!(work > 0.0 && (work - expected).abs() < 3e-12 * expected);
    assert!((response.input_power_w() - expected).abs() < 3e-12 * expected);
    assert_eq!(
        work,
        prepared
            .bulk_dissipation_w(&response.amplitudes, &|| false)
            .unwrap()
    );
    assert!(prepared.bulk_dissipation_w(&[], &|| false).is_err());
    assert!(
        prepared
            .bulk_dissipation_w(&response.amplitudes, &|| true)
            .is_err()
    );
    let outgoing = prepared.scattering(&|| false).unwrap()[0] * incident;
    let wave_work = 0.5 * (incident.norm_sq() - outgoing.norm_sq());
    assert!((wave_work - work).abs() < 5e-12 * work);
}

fn complex_sin_cos(z: C64) -> (C64, C64) {
    // Independent analytic oracle, using ordinary libm only in tests.
    (
        C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh()),
        C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh()),
    )
}

#[test]
fn full_pressure_mesh_converges_to_complex_bulk_duct_with_positive_volume_loss() {
    let dimensions = [0.1, 0.2, 0.15];
    let f = 500.0;
    let omega = std::f64::consts::TAU * f;
    let bulk = C64::new(139000.0, 6000.0);
    let k = (C64::from_re(AIR.density_kg_m3) / bulk).sqrt().scale(omega);
    assert!(k.re > 0.0 && k.im < 0.0);
    let zc = C64::from_re(AIR.density_kg_m3 * omega / (dimensions[1] * dimensions[2])) / k;
    let (sin, cos) = complex_sin_cos(k.scale(dimensions[0]));
    let minus_i = C64::new(0.0, -1.0);
    let z = [
        minus_i * zc * cos / sin,
        minus_i * zc / sin,
        minus_i * zc / sin,
        minus_i * zc * cos / sin,
    ];
    let references = [12000.0, 19000.0];
    let incident = [C64::new(0.7, 0.2), C64::new(-0.3, 0.6)];
    let mut a = z;
    a[0] = a[0] + C64::from_re(references[0]);
    a[3] = a[3] + C64::from_re(references[1]);
    let mut q = [
        incident[0].scale(2.0 * references[0].sqrt()),
        incident[1].scale(2.0 * references[1].sqrt()),
    ];
    lu_complex_with_cancel(&a, 2, || false)
        .unwrap()
        .solve_with_cancel(&mut q, || false)
        .unwrap();
    let expected = [z[0] * q[0] + z[1] * q[1], z[2] * q[0] + z[3] * q[1]];
    let mut errors = Vec::new();
    for n in [1, 2, 3] {
        let (mesh, points) = box_mesh(n, dimensions);
        let fe = HelmholtzFe::new(&mesh, &points, 2, AIR).unwrap();
        let ports = [
            BoundaryMode::uniform(&end(&fe, 0.0)),
            BoundaryMode::uniform(&end(&fe, dimensions[0])),
        ];
        let model = InteriorMultiport::new(&fe, &complete(&fe), &ports, false, &|| false).unwrap();
        let prepared = model
            .prepare_scattering_with_bulk_modulus(f, bulk, &references, &|| false)
            .unwrap();
        let response = prepared.response(&incident, &|| false).unwrap();
        let error = (0..2)
            .map(|i| {
                relative(response.pressure_pa[i], expected[i])
                    .max(relative(response.inward_flow_m3_s[i], q[i]))
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "complex-bulk pressure mesh n={n} dofs={} field_error={error:.15e}",
            model.dim()
        );
        errors.push(error);
        let work = prepared
            .bulk_dissipation_w(&response.amplitudes, &|| false)
            .unwrap();
        assert!(work > 0.0);
        assert!((response.input_power_w() - work).abs() < 2e-10 * work);
        let s = prepared.scattering(&|| false).unwrap();
        assert!((s[1] - s[2]).abs() < 2e-12);
        let b = [
            s[0] * incident[0] + s[1] * incident[1],
            s[2] * incident[0] + s[3] * incident[1],
        ];
        let wave_work = 0.5
            * (incident.iter().map(|z| z.norm_sq()).sum::<f64>()
                - b.iter().map(|z| z.norm_sq()).sum::<f64>());
        assert!((wave_work - work).abs() < 2e-10 * work);
        let ordinary = model.prepare_scattering(f, &references, &|| false).unwrap();
        let real = model
            .prepare_scattering_with_bulk_modulus(
                f,
                C64::from_re(AIR.density_kg_m3 * AIR.sound_speed_m_s.powi(2)),
                &references,
                &|| false,
            )
            .unwrap();
        for (&a, b) in ordinary
            .scattering(&|| false)
            .unwrap()
            .iter()
            .zip(real.scattering(&|| false).unwrap())
        {
            assert!((a - b).abs() < 2e-12);
        }
        let real_field = real.response(&incident, &|| false).unwrap();
        assert_eq!(
            real.bulk_dissipation_w(&real_field.amplitudes, &|| false)
                .unwrap(),
            0.0
        );
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 3e-5, "{errors:?}");
}

#[test]
fn passive_bulk_admission_refuses_wrong_phase_and_equilibrium_wall_overlay() {
    let (mesh, points) = box_mesh(1, [0.3, 0.15, 0.12]);
    let fe = HelmholtzFe::new(&mesh, &points, 2, AIR).unwrap();
    let left = end(&fe, 0.0);
    let right = end(&fe, 0.3);
    let ports = [BoundaryMode::uniform(&left), BoundaryMode::uniform(&right)];
    let basis = complete(&fe);
    let model = InteriorMultiport::new(&fe, &basis, &ports, false, &|| false).unwrap();
    for bad in [
        C64::ZERO,
        C64::new(-1.0, 1.0),
        C64::new(140000.0, -1.0),
        C64::new(f64::INFINITY, 0.0),
        C64::new(140000.0, f64::NAN),
    ] {
        assert!(
            model
                .prepare_scattering_with_bulk_modulus(500.0, bad, &[1e4; 2], &|| false)
                .is_err()
        );
    }
    let bulk = C64::new(140000.0, 2000.0);
    assert!(
        model
            .prepare_scattering_with_bulk_modulus(0.0, bulk, &[1e4; 2], &|| false)
            .is_err()
    );
    assert!(
        model
            .prepare_scattering_with_bulk_modulus(500.0, bulk, &[1e4; 2], &|| true)
            .is_err()
    );
    let walls = fe
        .boundary_faces()
        .into_iter()
        .filter(|f| !left.contains(f) && !right.contains(f))
        .collect::<Vec<_>>();
    let with_wall = InteriorMultiport::new_with_boundary_layers(
        &fe,
        &basis,
        &ports,
        false,
        &walls,
        TRANSPORT,
        GEOMETRY,
        &|| false,
    )
    .unwrap();
    assert!(matches!(
        with_wall.prepare_scattering_with_bulk_modulus(500.0, bulk, &[1e4; 2], &|| false),
        Err(AcousticFeError::InvalidInput {
            what: "molecular bulk modulus cannot overlay Helmholtz boundary layers"
        })
    ));
}
