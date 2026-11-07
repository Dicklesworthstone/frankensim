//! G1 independent continuity/force solutions; G3 modal work and normalization.
use super::*;
use crate::vibroacoustic::{AcousticMedium, CavityModes, StructuralModes};

fn close(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 3e-11 * a.abs().max(b.abs()).max(1e-25),
        "{a:?} != {b:?}"
    );
}

fn scalar(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * a.abs().max(b.abs()).max(1e-25),
        "{a:e} != {b:e}"
    );
}

fn one_mode(medium: AcousticMedium) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![150.0],
            shapes: vec![vec![1.0]],
            loss_factor: 0.02,
        },
        &CavityModes {
            omegas: vec![0.0],
            lambdas: vec![0.002],
            interface: vec![vec![1.0]],
            loss_factor: 0.0,
            rho0: medium.rho0,
            c0: medium.c0,
        },
        vec![0.003],
        Some(vec![C64::new(0.3, -0.1)]),
    )
    .unwrap()
}

fn boundary(
    omega: f64,
    shapes: &[Vec<f64>],
    areas: &[f64],
    ys: &[C64],
) -> CavityBoundaryAdmittance {
    CavityBoundaryAdmittance::from_samples(omega, shapes, areas, ys, &CancelGate::new_clock_free())
        .unwrap()
}

#[test]
fn g1_impedance_boundary_matches_continuity_and_wall_reaction() {
    let medium = AcousticMedium::air();
    let model = one_mode(medium);
    let bulk = C64::from_re(medium.rho0 * medium.c0 * medium.c0);
    let force = C64::new(0.2, -0.1);
    let flow = C64::new(2e-5, -7e-6);
    for omega in [5.0, 80.0, 150.0, 330.0, 1500.0] {
        // The same positive resistance plus compliant-wall admittance at
        // each frequency. No normal-incidence rho*c normalization is hidden.
        let y = C64::new(2e-4, -omega * 1e-7);
        let loaded = model
            .clone()
            .with_cavity_boundary_admittance(boundary(omega, &[vec![1.0]], &[0.02], &[y]))
            .unwrap();
        let result = loaded.frf_with_sources(omega, &[force], &[flow]).unwrap();
        let iw = C64::new(0.0, omega);
        let d = C64::new(150.0_f64.powi(2) - omega * omega, -0.02 * 150.0_f64.powi(2))
            - iw * C64::new(0.3, -0.1);
        // Independent physical equations:
        // D b = F + A p; Q = Y_tot p - i omega (V p/K + A b).
        let p = (flow + iw * force.scale(0.003) / d)
            / (y.scale(0.02)
                - iw * (C64::from_re(0.002) / bulk + C64::from_re(0.003_f64.powi(2)) / d));
        let b = (force + p.scale(0.003)) / d;
        close(result.a[0], p);
        close(result.b[0], b);
        scalar(
            result.power.cavity_boundary,
            0.5 * 0.02 * y.re * p.norm_sq(),
        );
        scalar(result.power.cavity, result.power.cavity_boundary);
        assert!(result.power.cavity_boundary > 0.0);
        assert!(result.power.balance_residual.abs() < 3e-11);
        assert!(result.power.interface_residual.abs() < 3e-11);
        let rigid = model
            .clone()
            .with_cavity_boundary_admittance(boundary(omega, &[vec![1.0]], &[0.02], &[C64::ZERO]))
            .unwrap()
            .frf_with_sources(omega, &[force], &[flow])
            .unwrap();
        let original = model.frf_with_sources(omega, &[force], &[flow]).unwrap();
        close(rigid.a[0], original.a[0]);
        close(rigid.b[0], original.b[0]);
        assert_eq!(rigid.power.cavity_boundary, 0.0);
    }
}

fn multimode(scale: [f64; 2], omega: f64) -> VibroacousticModel {
    let shapes = vec![vec![scale[0], scale[0]], vec![2.0 * scale[1], -scale[1]]];
    let cavity = CavityModes {
        omegas: vec![0.0, 210.0],
        lambdas: vec![0.002 * scale[0].powi(2), 0.001 * scale[1].powi(2)],
        interface: shapes.clone(),
        loss_factor: 0.0,
        rho0: 1.2,
        c0: 340.0,
    };
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![150.0],
            shapes: vec![vec![1.0]],
            loss_factor: 0.02,
        },
        &cavity,
        vec![0.003 * scale[0], -0.001 * scale[1]],
        None,
    )
    .unwrap()
    .with_cavity_boundary_admittance(boundary(
        omega,
        &shapes,
        &[0.01, 0.02],
        &[C64::new(2e-4, 1e-4), C64::new(3e-4, -2e-4)],
    ))
    .unwrap()
}

#[test]
fn g3_cross_mode_boundary_work_and_basis_changes_are_physical() {
    let omega = 180.0;
    let original = multimode([1.0; 2], omega);
    let forces = [C64::new(0.3, -0.1)];
    let flows = [C64::new(3e-5, 2e-6), C64::new(-1e-5, 3e-6)];
    for keep in [1, 2] {
        let response = original
            .frf_truncated_with_sources(omega, &forces, &flows, 1, keep)
            .unwrap();
        let mut p = [response.a[0]; 2];
        if keep == 2 {
            p[0] = p[0] + response.a[1].scale(2.0);
            p[1] = p[1] - response.a[1];
        }
        let work = 0.5 * (0.01 * 2e-4 * p[0].norm_sq() + 0.02 * 3e-4 * p[1].norm_sq());
        scalar(response.power.cavity_boundary, work);
        assert!(response.power.balance_residual.abs() < 3e-11);
        assert!(response.power.interface_residual.abs() < 3e-11);
        let scaled = multimode([-3.0, 0.2], omega)
            .frf_truncated_with_sources(
                omega,
                &forces,
                &[flows[0].scale(-3.0), flows[1].scale(0.2)],
                1,
                keep,
            )
            .unwrap();
        close(response.b[0], scaled.b[0]);
        close(response.a[0], scaled.a[0].scale(-3.0));
        if keep == 2 {
            close(response.a[1], scaled.a[1].scale(0.2));
        }
        scalar(response.power.cavity_boundary, scaled.power.cavity_boundary);
    }
    let b = original.cavity_boundary.as_ref().unwrap();
    assert_ne!(
        b.modal_admittance()[1],
        C64::ZERO,
        "nonuniform loss couples modes"
    );
    assert_eq!(b.modal_admittance()[1], b.modal_admittance()[2]);
}

#[test]
fn g0_boundary_pressure_nullspace_stays_lossless() {
    let omega = 130.0;
    let cavity = CavityModes {
        omegas: vec![0.0, 210.0],
        lambdas: vec![0.002, 0.001],
        interface: vec![vec![1.0]; 2],
        loss_factor: 0.0,
        rho0: 1.2,
        c0: 340.0,
    };
    let model = VibroacousticModel::try_new(
        &StructuralModes {
            omegas: vec![150.0],
            shapes: vec![vec![1.0]],
            loss_factor: 0.0,
        },
        &cavity,
        vec![0.0; 2],
        None,
    )
    .unwrap()
    .with_cavity_boundary_admittance(boundary(
        omega,
        &cavity.interface,
        &[0.02],
        &[C64::new(3e-3, -1e-3)],
    ))
    .unwrap();
    let desired = [C64::new(1.0, 0.3), C64::new(-1.0, -0.3)];
    let flows: Vec<_> = (0..2)
        .map(|q| {
            C64::new(
                0.0,
                (cavity.omegas[q].powi(2) - omega * omega) * cavity.lambdas[q]
                    / (omega * cavity.rho0 * cavity.c0 * cavity.c0),
            ) * desired[q]
        })
        .collect();
    let result = model.frf_with_sources(omega, &[C64::ZERO], &flows).unwrap();
    close(result.a[0], desired[0]);
    close(result.a[1], desired[1]);
    assert!(result.power.cavity_boundary < 1e-30);
}

#[test]
fn g1_molecular_bulk_loss_and_boundary_loss_close_separately() {
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
    let omega = 700.0;
    let model = one_mode(AcousticMedium::try_from_gas(record.state()).unwrap())
        .with_cavity_relaxation(relaxation.clone())
        .unwrap()
        .with_cavity_boundary_admittance(boundary(
            omega,
            &[vec![1.0]],
            &[0.02],
            &[C64::new(2e-4, -7e-5)],
        ))
        .unwrap();
    let result = model
        .frf_with_sources(omega, &[C64::ZERO], &[C64::new(2e-5, -7e-6)])
        .unwrap();
    let beta = relaxation
        .caloric_response(omega)
        .unwrap()
        .bulk_modulus_pa
        .recip();
    let bulk_work = 0.5 * omega * 0.002 * beta.im * result.a[0].norm_sq();
    scalar(
        result.power.cavity,
        result.power.cavity_boundary + bulk_work,
    );
    assert!(bulk_work > 0.0 && result.power.cavity_boundary > 0.0);
    assert!(result.power.balance_residual.abs() < 3e-11);
    assert!(result.power.interface_residual.abs() < 3e-11);
}

#[test]
fn invalid_active_duplicate_or_wrong_frequency_boundary_refuses() {
    let gate = CancelGate::new_clock_free();
    let shapes = [vec![1.0]];
    for y in [C64::new(-1e-12, 0.0), C64::new(0.0, f64::NAN)] {
        assert!(
            CavityBoundaryAdmittance::from_samples(100.0, &shapes, &[0.02], &[y], &gate).is_err()
        );
    }
    assert!(
        CavityBoundaryAdmittance::from_samples(100.0, &shapes, &[-0.02], &[C64::ONE], &gate)
            .is_err()
    );
    assert!(
        CavityBoundaryAdmittance::from_samples(100.0, &shapes, &[0.02, 0.03], &[C64::ONE], &gate)
            .is_err()
    );
    let b = boundary(100.0, &shapes, &[0.02], &[C64::from_re(2e-4)]);
    assert_eq!(b.angular_frequency_rad_s(), 100.0);
    let base = one_mode(AcousticMedium::air());
    let loaded = base
        .clone()
        .with_cavity_boundary_admittance(b.clone())
        .unwrap();
    assert!(loaded.frf(101.0, &[C64::ONE]).is_err());
    assert!(loaded.with_cavity_boundary_admittance(b.clone()).is_err());
    let mut scalar_loss = base;
    scalar_loss.eta_a = 0.01;
    assert!(scalar_loss.with_cavity_boundary_admittance(b).is_err());
    gate.request();
    assert!(matches!(
        CavityBoundaryAdmittance::from_samples(100.0, &shapes, &[0.02], &[C64::ONE], &gate),
        Err(VibroError::Cancelled)
    ));
}
