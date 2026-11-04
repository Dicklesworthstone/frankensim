//! G1 driven sealed-cavity compliance; G3 coherent work, reciprocity and basis scaling.
use super::*;

fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 2e-11 * scale.max(1e-25), "{a:e} != {b:e}");
}

fn complex_near(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(1e-25),
        "{a:?} != {b:?}"
    );
}

fn model(medium: AcousticMedium, cavity_scale: f64) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![150.0, 280.0],
            shapes: vec![vec![1.0], vec![1.0]],
            loss_factor: 0.02,
        },
        &CavityModes {
            omegas: vec![0.0, 210.0],
            lambdas: vec![0.002 * cavity_scale * cavity_scale, 0.001],
            interface: vec![vec![cavity_scale], vec![1.0]],
            loss_factor: 0.0,
            rho0: medium.rho0,
            c0: medium.c0,
        },
        vec![0.003 * cavity_scale, 0.001, -0.001 * cavity_scale, 0.002],
        Some(vec![
            C64::new(0.3, -0.1),
            C64::ZERO,
            C64::ZERO,
            C64::from_re(0.2),
        ]),
    )
    .unwrap()
}

#[test]
fn boundary_projection_preserves_cancelled_modes_and_physical_surface_work() {
    let cavity = CavityModes {
        omegas: vec![0.0, 100.0],
        lambdas: vec![0.1, 0.05],
        interface: vec![vec![1.0, 1.0], vec![2.0, -1.0]],
        loss_factor: 0.0,
        rho0: 1.2,
        c0: 340.0,
    };
    let areas = [0.01, 0.02];
    let velocity = [C64::new(0.2, -0.1); 2];
    let flow = project_cavity_boundary_flow(&cavity, &areas, &velocity).unwrap();
    complex_near(flow[0], velocity[0].scale(0.03));
    assert_eq!(flow[1], C64::ZERO, "a cancelling mode must remain silent");
    let velocity = [velocity[0], C64::new(-0.13, 0.07)];
    let flow = project_cavity_boundary_flow(&cavity, &areas, &velocity).unwrap();
    let pressure_modes = [C64::new(4.0, -1.0), C64::new(2.0, 3.0)];
    let modal_work: f64 = (0..2)
        .map(|q| 0.5 * (pressure_modes[q] * flow[q].conj()).re)
        .sum();
    let surface_work: f64 = (0..2)
        .map(|i| {
            let pressure = pressure_modes[0] + pressure_modes[1].scale(cavity.interface[1][i]);
            0.5 * areas[i] * (pressure * velocity[i].conj()).re
        })
        .sum();
    near(modal_work, surface_work, surface_work.abs());
    assert!(project_cavity_boundary_flow(&cavity, &areas[..1], &velocity).is_err());
    assert!(project_cavity_boundary_flow(&cavity, &[-0.01, 0.02], &velocity).is_err());
}

#[test]
fn injected_volume_drives_the_wall_and_radiation_with_analytic_sealed_pressure() {
    let medium = AcousticMedium::air();
    let coupled = model(medium, 1.0);
    let bulk = C64::from_re(medium.rho0 * medium.c0 * medium.c0);
    let flow = C64::new(2e-5, -7e-6);
    for omega in [5.0, 80.0, 150.0, 330.0] {
        for force in [C64::ZERO, C64::new(0.3, 0.1)] {
            let result = coupled
                .frf_truncated_with_sources(omega, &[force, C64::ZERO], &[flow, C64::ZERO], 1, 1)
                .unwrap();
            let z = C64::new(0.3, -0.1);
            let d = C64::new(150.0_f64.powi(2) - omega * omega, -0.02 * 150.0_f64.powi(2))
                - C64::new(0.0, omega) * z;
            // Independent continuity and wall force balance:
            // p = K/V (i Q/omega - A b), D b = F + A p.
            let volume_injected = C64::new(0.0, 1.0 / omega) * flow;
            let b = (force + bulk.scale(0.003 / 0.002) * volume_injected)
                / (d + bulk.scale(0.003 * 0.003 / 0.002));
            let p = bulk.scale(1.0 / 0.002) * (volume_injected - b.scale(0.003));
            complex_near(result.b[0], b);
            complex_near(result.a[0], p);
            let v = C64::new(0.0, -omega) * b;
            let acoustic_power = 0.5 * (p * flow.conj()).re;
            near(
                result.power.acoustic_input,
                acoustic_power,
                acoustic_power.abs(),
            );
            near(
                result.power.radiated,
                0.5 * z.re * v.norm_sq(),
                result.power.radiated,
            );
            assert!(result.power.balance_residual.abs() < 1e-11);
            assert!(result.power.interface_residual.abs() < 1e-11);
            if force == C64::ZERO {
                assert!(
                    result.power.interface < 0.0,
                    "the driven cavity powers the wall"
                );
                assert!(result.power.radiated > 0.0);
                assert_eq!(result.power.structural_input, 0.0);
            }
        }
    }
}

#[test]
fn acoustic_drive_retains_molecular_compressibility_and_independent_cycle_work() {
    use fs_material::gas_species::input::read_gas_record;
    use fs_material::relaxation::{RelaxationAbsorption, input::medium_from_record};
    let record = read_gas_record(
        &mut "mixture 293.15 101325 mole 4 N2 0.77319 O2 0.20691 Ar 0.0099 H2O 0.01"
            .split_ascii_whitespace(),
        None,
        None,
    )
    .unwrap();
    let relaxation = RelaxationAbsorption::iso9613_terrestrial_air()
        .unwrap()
        .evaluate(&medium_from_record(&record).unwrap())
        .unwrap();
    let medium = AcousticMedium::try_from_gas(record.state()).unwrap();
    let coupled = model(medium, 1.0)
        .with_cavity_relaxation(relaxation.clone())
        .unwrap();
    let omega = 700.0;
    let flow = C64::new(2e-5, -7e-6);
    let result = coupled
        .frf_truncated_with_sources(omega, &[C64::ZERO; 2], &[flow, C64::ZERO], 1, 1)
        .unwrap();
    let bulk = relaxation.caloric_response(omega).unwrap().bulk_modulus_pa;
    let d = C64::new(150.0_f64.powi(2) - omega * omega, -0.02 * 150.0_f64.powi(2))
        - C64::new(0.0, omega) * C64::new(0.3, -0.1);
    let b = (bulk.scale(0.003 / 0.002) * C64::new(0.0, 1.0 / omega) * flow)
        / (d + bulk.scale(0.003 * 0.003 / 0.002));
    let p = bulk.scale(1.0 / 0.002) * (C64::new(0.0, 1.0 / omega) * flow - b.scale(0.003));
    complex_near(result.b[0], b);
    complex_near(result.a[0], p);
    let mut source_work = 0.0;
    let mut molecular_heat = 0.0;
    let compression_flow = flow + C64::new(0.0, omega * 0.003) * b;
    for sample in 0..128 {
        let phase = core::f64::consts::TAU * sample as f64 / 128.0;
        let (sin, cos) = phase.sin_cos();
        let real = |z: C64| z.re * cos + z.im * sin;
        source_work += real(p) * real(flow) / 128.0;
        molecular_heat += real(p) * real(compression_flow) / 128.0;
    }
    near(result.power.acoustic_input, source_work, source_work.abs());
    near(result.power.cavity, molecular_heat, source_work.abs());
    assert!(result.power.cavity > 0.0);
    assert!(result.power.balance_residual.abs() < 1e-11);
    assert!(result.power.interface_residual.abs() < 1e-11);
}

#[test]
fn coherent_sources_keep_reciprocity_signed_work_and_modal_normalization() {
    let medium = AcousticMedium::air();
    let ordinary = model(medium, 1.0);
    let reversed = model(medium, -3.0);
    let omega = 180.0;
    let forces = [C64::new(-0.08, 0.04), C64::new(0.3, -0.1)];
    let flows = [C64::new(3e-5, 2e-6), C64::new(-1e-5, 3e-6)];
    let full = ordinary.frf_with_sources(omega, &forces, &flows).unwrap();
    assert!(full.power.structural_input > 0.0);
    assert!(
        full.power.acoustic_input < 0.0,
        "a coherent source may absorb work"
    );
    let scaled = reversed
        .frf_with_sources(omega, &forces, &[flows[0].scale(-3.0), flows[1]])
        .unwrap();
    for r in 0..2 {
        complex_near(full.b[r], scaled.b[r]);
    }
    complex_near(full.a[0], scaled.a[0].scale(-3.0));
    complex_near(full.a[1], scaled.a[1]);
    near(full.power.input, scaled.power.input, full.power.input.abs());
    near(
        full.power.acoustic_input,
        scaled.power.acoustic_input,
        full.power.input.abs(),
    );
    let force_only = ordinary.frf(omega, &forces).unwrap();
    let flow_only = ordinary
        .frf_with_sources(omega, &[C64::ZERO; 2], &flows)
        .unwrap();
    for r in 0..2 {
        complex_near(full.b[r], force_only.b[r] + flow_only.b[r]);
        complex_near(full.a[r], force_only.a[r] + flow_only.a[r]);
    }
    let zero_flow = ordinary
        .frf_with_sources(omega, &forces, &[C64::ZERO; 2])
        .unwrap();
    assert_eq!(zero_flow.a, force_only.a);
    assert_eq!(zero_flow.b, force_only.b);
    assert_eq!(zero_flow.power.acoustic_input, 0.0);
    // Mixed force/flow transfer is antisymmetric: the structural velocity
    // points OUT of the cavity, whereas source flow is positive INTO it.
    let mechanical = ordinary.frf(omega, &[C64::ONE, C64::ZERO]).unwrap();
    let acoustic = ordinary
        .frf_with_sources(omega, &[C64::ZERO; 2], &[C64::ONE, C64::ZERO])
        .unwrap();
    complex_near(
        C64::new(0.0, -omega) * acoustic.b[0],
        mechanical.a[0].scale(-1.0),
    );
    assert!(full.power.balance_residual.abs() < 1e-11);
    assert!(full.power.interface_residual.abs() < 1e-11);
}

#[test]
fn acoustic_source_truncation_uses_the_original_full_basis_and_checks_sources() {
    let coupled = model(AcousticMedium::air(), 1.0);
    let forces = [C64::ZERO; 2];
    let flows = [C64::new(2e-5, 1e-6), C64::new(8e-6, -3e-6)];
    let full = coupled.frf_with_sources(190.0, &forces, &flows).unwrap();
    let half = coupled
        .frf_truncated_with_sources(190.0, &forces, &flows, 1, 1)
        .unwrap();
    let delta = (full.power.input - half.power.input).abs() / full.power.input.abs();
    assert!(delta > 1e-5);
    let checked = coupled
        .frf_with_sources_and_convergence(190.0, &forces, &flows, 2.0 * delta)
        .unwrap();
    near(checked.truncation_delta, delta, delta);
    assert_eq!(checked.mode_counts, (2, 2, 1, 1));
    assert!(matches!(
        coupled.frf_with_sources_and_convergence(190.0, &forces, &flows, 0.5 * delta),
        Err(VibroError::TruncationNotConverged { .. })
    ));
    assert!(matches!(
        coupled.frf_with_sources(190.0, &forces, &flows[..1]),
        Err(VibroError::ShapeMismatch { .. })
    ));
    assert!(
        coupled
            .frf_truncated_with_sources(190.0, &forces, &flows, 0, 1)
            .is_err()
    );
    for bad in [f64::NAN, f64::INFINITY] {
        // Even a source outside the retained prefix remains an invalid input.
        assert!(
            coupled
                .frf_truncated_with_sources(190.0, &forces, &[flows[0], C64::new(bad, 0.0)], 1, 1)
                .is_err()
        );
        assert!(
            coupled
                .frf(190.0, &[C64::new(bad, 0.0), C64::ZERO])
                .is_err()
        );
    }
}
