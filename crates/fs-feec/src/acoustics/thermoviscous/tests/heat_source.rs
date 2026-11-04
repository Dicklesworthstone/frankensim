use super::*;

fn sealed(mesh: &TetComplex, positions: &[[f64; 3]], size: [f64; 3]) -> ThermoviscousBoundary {
    let mut boundary = boundaries(mesh, positions, size, false);
    for port in boundary.ports.drain(..) {
        boundary
            .no_slip_faces
            .extend(port.traces.into_iter().map(|(face, _)| face));
    }
    boundary
}

#[test]
fn g1_closed_volume_heat_release_has_the_original_fluid_eos_and_phase() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let gas = air();
    let thermodynamics = ThermoviscousThermodynamics::ideal_gas(gas).unwrap();
    let frequency = 1300.0;
    let omega = std::f64::consts::TAU * frequency;
    let amplitude = C64::new(0.001, 0.0003);
    let cv_volume = gas.density_kg_m3 * (gas.cp_j_kg_k - gas.gas_constant_j_kg_k());
    let expected_temperature = amplitude / C64::new(0.0, omega * volume * cv_volume);
    let expected_pressure = expected_temperature.scale(gas.pressure_pa / gas.temperature_k);
    // Three original-mesh rungs. This spatially uniform manufactured solution
    // is representable at every rung; no fitted acoustic decay or pressure law.
    for cells in [[1, 1, 1], [2, 1, 1], [3, 1, 1]] {
        let (mesh, positions) = box_mesh(cells, size);
        let boundary = sealed(&mesh, &positions, size);
        let model = ThermoviscousFe::new_with_thermodynamics(
            &mesh,
            &positions,
            gas,
            thermodynamics,
            &boundary,
            limits(),
            &mut || false,
        )
        .unwrap();
        let source = model
            .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
            .unwrap();
        assert!((source.integrated_shape() - 1.0).abs() < 2e-14);
        let prepared = model.prepare(frequency, &mut || false).unwrap();
        let field = prepared
            .response_with_heat_sources(&[], &[], &[(&source, amplitude)], &mut || false)
            .unwrap();
        for pressure in &field.pressure_pa {
            assert!((*pressure - expected_pressure).abs() < 2e-9 * expected_pressure.abs());
        }
        for temperature in &field.temperature_k[..mesh.vertex_count] {
            assert!(
                (*temperature - expected_temperature).abs() < 2e-9 * expected_temperature.abs()
            );
        }
        for temperature in &field.temperature_k[mesh.vertex_count..] {
            assert!(temperature.abs() < 2e-9 * expected_temperature.abs());
        }
        assert!(field.velocity_m_s.iter().flatten().all(|v| v.abs() < 1e-10));
        let scale = amplitude.abs() * expected_temperature.abs() / gas.temperature_k;
        assert!(field.power.heat_source_work_w.abs() < 2e-9 * scale);
        assert!(field.power.balance_residual_w().abs() < 2e-9 * scale);
    }
}

#[test]
fn g1_sealed_heat_release_has_absolute_passive_caloric_work() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let boundary = sealed(&mesh, &positions, size);
    let gas = air();
    let model =
        ThermoviscousFe::new(&mesh, &positions, gas, &boundary, limits(), &mut || false).unwrap();
    // A signed, unnormalized shape must retain both its sign and magnitude.
    let shape_integral = -2.5;
    let source = model
        .volume_heat_source(4, &mut |_, _| Ok(shape_integral / volume), &mut || false)
        .unwrap();
    assert!((source.integrated_shape() - shape_integral).abs() < 5e-14);
    let drive = C64::new(0.001, 0.0003);
    let total_heat = drive.scale(shape_integral);
    let cp = C64::new(970.0, -45.0);
    let cv = cp - C64::from_re(gas.gas_constant_j_kg_k());
    for frequency in [80.0, 1700.0, 25000.0] {
        let omega = std::f64::consts::TAU * frequency;
        let field = model
            .prepare_with_heat_capacity(frequency, cp, &mut || false)
            .unwrap()
            .response_with_heat_sources(&[], &[], &[(&source, drive)], &mut || false)
            .unwrap();
        let temperature = total_heat / (cv * C64::new(0.0, omega * volume * gas.density_kg_m3));
        let pressure = temperature.scale(gas.pressure_pa / gas.temperature_k);
        for actual in &field.temperature_k[..mesh.vertex_count] {
            assert!((*actual - temperature).abs() < 2e-9 * temperature.abs());
        }
        for actual in &field.pressure_pa {
            assert!((*actual - pressure).abs() < 2e-9 * pressure.abs());
        }
        // The uniform sealed-volume solution gives this positive real work
        // directly from the prescribed heat and constitutive coefficients.
        // It independently fixes both the phasor factor 1/2 and the 1/T0 scale.
        let expected_power = total_heat.norm_sq() * -cp.im
            / (2.0 * omega * volume * gas.density_kg_m3 * gas.temperature_k * cv.norm_sq());
        assert!(expected_power > 0.0);
        assert!((field.power.heat_source_work_w - expected_power).abs() < 2e-9 * expected_power);
        assert!((field.power.caloric_w - expected_power).abs() < 2e-9 * expected_power);
        assert_work(&field);
    }
}

#[test]
fn g2_heat_release_and_traction_share_resolved_losses_and_the_same_factor() {
    let size = [0.004, 0.001, 0.0012];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &positions, size, true);
    let model =
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap();
    let source = model
        .volume_heat_source(
            4,
            &mut |_, point| Ok(2.0 * point[0] / size[0] / volume),
            &mut || false,
        )
        .unwrap();
    assert!((source.integrated_shape() - 1.0).abs() < 2e-14);
    let differential = model
        .volume_heat_source(
            4,
            &mut |_, point| Ok((2.0 * point[0] / size[0] - 1.0) / volume),
            &mut || false,
        )
        .unwrap();
    assert!(differential.integrated_shape().abs() < 2e-14);
    let prepared = model
        .prepare_with_heat_capacity(1700.0, C64::new(970.0, -45.0), &mut || false)
        .unwrap();
    let drive = [C64::new(1.0, 0.2), C64::new(-0.3, 0.1)];
    let amplitude = C64::new(0.0003, -0.0002);
    let differential_amplitude = C64::new(-0.0001, 0.0002);
    let combined = prepared
        .response_with_heat_sources(
            &drive,
            &[],
            &[
                (&source, amplitude),
                (&differential, differential_amplitude),
            ],
            &mut || false,
        )
        .unwrap();
    assert_work(&combined);
    assert!(combined.power.conductive_w > 0.0);
    assert!(combined.power.caloric_w > 0.0);
    let mechanical = prepared.response(&drive, &mut || false).unwrap();
    let thermal = prepared
        .response_with_heat_sources(
            &[C64::ZERO; 2],
            &[],
            &[
                (&source, amplitude),
                (&differential, differential_amplitude),
            ],
            &mut || false,
        )
        .unwrap();
    assert_work(&thermal);
    assert!(thermal.power.heat_source_work_w > 0.0);
    for ((actual, a), b) in combined
        .pressure_pa
        .iter()
        .zip(&mechanical.pressure_pa)
        .zip(&thermal.pressure_pa)
    {
        assert!((*actual - *a - *b).abs() < 1e-10 * (a.abs() + b.abs()));
    }
    let split = prepared
        .response_with_heat_sources(
            &drive,
            &[],
            &[
                (&source, amplitude.scale(0.25)),
                (&source, amplitude.scale(0.75)),
                (&differential, differential_amplitude),
            ],
            &mut || false,
        )
        .unwrap();
    assert_work(&split);
    for (actual, expected) in split.pressure_pa.iter().zip(&combined.pressure_pa) {
        assert!((*actual - *expected).abs() < 1e-10 * expected.abs());
    }
    let zero = prepared
        .response_with_heat_sources(&drive, &[], &[], &mut || false)
        .unwrap();
    assert_eq!(zero.pressure_pa, mechanical.pressure_pa);
    assert_eq!(zero.power.heat_source_work_w, 0.0);
}

#[test]
fn g1_liquid_heat_release_uses_signed_expansion_and_constant_volume_capacity() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let boundary = sealed(&mesh, &positions, size);
    // Manufactured stable fluid, not a measured liquid material card.
    let medium = ThermoviscousMedium {
        density_kg_m3: 1000.0,
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        cp_j_kg_k: 4000.0,
        dynamic_viscosity_pa_s: 0.001,
        thermal_conductivity_w_m_k: 0.6,
        bulk_viscosity_pa_s: 0.002,
    };
    let beta = 5.1e-10;
    let frequency = 1300.0;
    let omega = std::f64::consts::TAU * frequency;
    let amplitude = C64::new(0.001, 0.0003);
    for alpha in [2.0e-4, -2.0e-4, 0.0] {
        let model = ThermoviscousFe::new_with_thermodynamics(
            &mesh,
            &positions,
            medium,
            ThermoviscousThermodynamics {
                isothermal_compressibility_pa_inv: beta,
                isobaric_expansion_per_k: alpha,
            },
            &boundary,
            limits(),
            &mut || false,
        )
        .unwrap();
        let source = model
            .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
            .unwrap();
        let field = model
            .prepare(frequency, &mut || false)
            .unwrap()
            .response_with_heat_sources(&[], &[], &[(&source, amplitude)], &mut || false)
            .unwrap();
        let cv_volume =
            medium.density_kg_m3 * medium.cp_j_kg_k - medium.temperature_k * alpha * alpha / beta;
        let temperature = amplitude / C64::new(0.0, omega * volume * cv_volume);
        let pressure = temperature.scale(alpha / beta);
        for actual in &field.temperature_k[..mesh.vertex_count] {
            assert!((*actual - temperature).abs() < 2e-9 * temperature.abs());
        }
        for actual in &field.pressure_pa {
            assert!((*actual - pressure).abs() < 2e-9 * pressure.abs().max(1e-15));
        }
        assert!(field.velocity_m_s.iter().flatten().all(|v| v.abs() < 1e-10));
    }
}

#[test]
fn g4_heat_source_requires_its_original_owner_and_finite_cancel_correct_data() {
    let size = [0.004, 0.002, 0.002];
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let boundary = sealed(&mesh, &positions, size);
    let make = || {
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap()
    };
    let model = make();
    let other = make();
    let source = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0), &mut || false)
        .unwrap();
    let foreign = other
        .volume_heat_source(4, &mut |_, _| Ok(1.0), &mut || false)
        .unwrap();
    let prepared = model.prepare(1000.0, &mut || false).unwrap();
    assert!(
        prepared
            .response_with_heat_sources(&[], &[], &[(&foreign, C64::ONE)], &mut || false)
            .is_err()
    );
    assert!(
        prepared
            .response_with_heat_sources(
                &[],
                &[],
                &[(&source, C64::new(f64::NAN, 0.0))],
                &mut || false
            )
            .is_err()
    );
    for order in [0, 3, 17] {
        assert!(
            model
                .volume_heat_source(order, &mut |_, _| Ok(1.0), &mut || false)
                .is_err()
        );
    }
    assert!(
        model
            .volume_heat_source(4, &mut |_, _| Ok(f64::NAN), &mut || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        model.volume_heat_source(4, &mut |_, _| Ok(1.0), &mut || {
            polls += 1;
            polls == 9
        }),
        Err(ThermoviscousError::Cancelled)
    ));
    assert!(matches!(
        prepared.response_with_heat_sources(&[], &[], &[(&source, C64::ONE)], &mut || true),
        Err(ThermoviscousError::Cancelled)
    ));
}

fn same_resolved_fields(actual: &ThermoviscousField, expected: &ThermoviscousField) {
    for (a, b) in [
        (&actual.pressure_pa, &expected.pressure_pa),
        (&actual.temperature_k, &expected.temperature_k),
        (&actual.inward_flow_m3_s, &expected.inward_flow_m3_s),
    ] {
        let scale = b.iter().map(|z| z.abs()).fold(1e-40, f64::max);
        assert!(a.iter().zip(b).all(|(&x, &y)| (x - y).abs() < 2e-8 * scale));
    }
}

#[test]
fn g2_loaded_heat_source_matches_independent_pressure_flow_closure() {
    let size = [0.004, 0.001, 0.0012];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &positions, size, true);
    let model =
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap();
    let source = model
        .volume_heat_source(
            4,
            &mut |_, point| Ok(2.0 * point[0] / size[0] / volume),
            &mut || false,
        )
        .unwrap();
    let sources = [(&source, C64::new(0.0003, -0.0002))];
    let cp = C64::new(970.0, -45.0);
    let reference: [f64; 2] = [1.3e8, 0.7e8];
    let incident = [C64::new(0.0001, 0.0002), C64::new(-0.0002, 0.00005)];
    let unloaded = model
        .prepare_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    let source_only = unloaded
        .response_with_heat_sources(&[C64::ZERO; 2], &[], &sources, &mut || false)
        .unwrap();
    let y = model
        .admittance_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    // Independent physical-port closure: Q=Y P+Q_source and
    // P+R Q=2 sqrt(R) a. No scattering matrix enters this oracle.
    let mut matrix = vec![C64::ZERO; 4];
    let mut pressure = vec![C64::ZERO; 2];
    for i in 0..2 {
        pressure[i] = incident[i].scale(2.0 * reference[i].sqrt())
            - source_only.inward_flow_m3_s[i].scale(reference[i]);
        for j in 0..2 {
            matrix[2 * i + j] =
                y[2 * i + j].scale(reference[i]) + C64::from_re(if i == j { 1.0 } else { 0.0 });
        }
    }
    lu_complex_with_cancel(&matrix, 2, || false)
        .unwrap()
        .solve_with_cancel(&mut pressure, || false)
        .unwrap();
    let expected = unloaded
        .response_with_heat_sources(&pressure, &[], &sources, &mut || false)
        .unwrap();
    let loaded = model
        .prepare_scattering(1700.0, cp, &reference, &mut || false)
        .unwrap();
    let actual = loaded
        .response_with_sources(&incident, &[], &sources, &mut || false)
        .unwrap();
    same_resolved_fields(&actual, &expected);
    assert_work(&actual);
    let offset = loaded
        .source_outgoing(&[], &sources, &mut || false)
        .unwrap();
    let scattering = loaded.scattering(&mut || false).unwrap();
    let mut boundary_work = 0.0;
    for i in 0..2 {
        let outgoing = (actual.pressure_drive_pa[i]
            - actual.inward_flow_m3_s[i].scale(reference[i]))
        .scale(0.5 / reference[i].sqrt());
        let expected = (0..2).fold(offset[i], |sum, j| {
            sum + scattering[2 * i + j] * incident[j]
        });
        assert!((outgoing - expected).abs() < 2e-9 * (outgoing.abs() + expected.abs()));
        boundary_work += 0.5 * (incident[i].norm_sq() - outgoing.norm_sq());
    }
    let scale = actual.power.dissipated_w().abs()
        + actual.power.heat_source_work_w.abs()
        + boundary_work.abs();
    assert!(
        (actual.power.input_w - boundary_work - actual.power.heat_source_work_w).abs()
            < 2e-9 * scale
    );
    let emitted = loaded
        .response_with_sources(&[C64::ZERO; 2], &[], &sources, &mut || false)
        .unwrap();
    assert_work(&emitted);
    assert!(emitted.power.heat_source_work_w > emitted.power.dissipated_w());
    assert!(offset.iter().any(|z| z.abs() > 0.0));
}

#[test]
fn g4_loaded_sources_preserve_original_drives_owner_and_cancellation() {
    let size = [0.004, 0.001, 0.0012];
    let (mesh, positions) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &positions, size, true);
    let make = || {
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap()
    };
    let model = make();
    let other = make();
    let foreign = other
        .volume_heat_source(4, &mut |_, _| Ok(1.0), &mut || false)
        .unwrap();
    let prepared = model
        .prepare_scattering(
            1700.0,
            C64::from_re(air().cp_j_kg_k),
            &[1e8; 2],
            &mut || false,
        )
        .unwrap();
    let incident = [C64::new(0.0001, 0.0002), C64::new(-0.0002, 0.00005)];
    let original = prepared.response(&incident, &mut || false).unwrap();
    let zero = prepared
        .response_with_sources(&incident, &[], &[], &mut || false)
        .unwrap();
    assert_eq!(original.pressure_pa, zero.pressure_pa);
    assert_eq!(original.temperature_k, zero.temperature_k);
    assert_eq!(zero.power.heat_source_work_w, 0.0);
    assert!(
        prepared
            .response_with_sources(&incident, &[], &[(&foreign, C64::ONE)], &mut || false)
            .is_err()
    );
    assert!(
        prepared
            .response_with_sources(&incident[..1], &[], &[], &mut || false)
            .is_err()
    );
    assert!(
        prepared
            .response_with_sources(&incident, &[C64::ONE], &[], &mut || false)
            .is_err()
    );
    assert!(matches!(
        prepared.response_with_sources(&incident, &[], &[], &mut || true),
        Err(ThermoviscousError::Cancelled)
    ));
}
