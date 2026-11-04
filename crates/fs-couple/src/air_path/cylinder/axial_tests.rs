//! Axial wall shear, independent complex-Bessel reference, and native work balance.
use super::*;

fn flow(density: f64, viscosity: f64, intervals: usize) -> CylinderFlow<FluidTransportState> {
    CylinderFlow::from_fluid(
        0.2,
        FluidTransportState::try_new(density, 4000.0, 0.5, viscosity).unwrap(),
        CylinderFlowResolution {
            minimum_wavenumber_radius: if intervals > 8 { 1e-5 } else { 0.01 },
            maximum_wavenumber_radius: if intervals > 8 { 3000.0 } else { 2.0 },
            intervals,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 1e-10 * (1.0 + a.abs() + b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn axial_impedance_matches_independent_decaying_bessel_solution() {
    // rho*a²/mu = 1 s. References evaluate 2*pi*z*K1(z)/K0(z) directly
    // with defining complex Bessel series, z=sqrt(i*omega); they do not
    // evaluate the real-Bessel quadrature or subtract its transverse result.
    let flow = flow(2.0, 0.08, 512);
    for (omega, re, im) in [
        (0.01, 2.355263668273521, 0.7180285506223489),
        (0.1, 3.7446710226358345, 1.6869423520690306),
        (1.0, 7.173152271637606, 4.662146819754941),
        (10.0, 17.028218762615095, 14.169253733746251),
        (100.0, 47.51554208582318, 44.47734831005997),
    ] {
        let z = flow.axial_impedance_per_length(omega).unwrap();
        let error = (z.re / 0.08 - re).hypot(z.im / 0.08 - im) / re.hypot(im);
        assert!(
            error < 2.5e-3,
            "omega={omega}, axial complex error={error:e}"
        );
        assert!(z.re > 0.0 && z.im > 0.0);
    }
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(flow.axial_impedance_per_length(omega).is_err());
    }
}

#[test]
fn axial_native_motion_has_dry_inertia_and_reciprocal_memory_work() {
    let flow = flow(2.0, 0.1, 4);
    let gate = CancelGate::new_clock_free();
    let body = flow.axial_body(0.7, 0.6, 4.0, 6, &gate).unwrap();
    assert_eq!(body.translation(), CylinderTranslation::Axial);
    assert_eq!(body.moving_mass_kg(), 0.6);
    assert!(body.memory_masses_kg().iter().all(|m| *m > 0.0));
    assert!(flow.axial_body(0.7, 0.6, 4.0, 5, &gate).is_err());
    let fluid_v = [0.1, -0.2, 0.25, -0.1];
    let mut x = body.state_from_motion(0.02, 0.3, &fluid_v).unwrap();
    let initial = body.system().hamiltonian(&x);
    let mut expected_energy = 0.5 * 4.0 * 0.02_f64.powi(2) + 0.5 * 0.6 * 0.3_f64.powi(2);
    for (k, branch) in flow.axial_branches().enumerate() {
        expected_energy += 0.5 * 0.7 * branch.mass_per_length_kg_m * fluid_v[k].powi(2);
    }
    near(initial, expected_energy);
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..50 {
        let (force, dt) = (0.03, 0.002);
        let record = fs_phs::step(body.system(), &x, &[force], dt).unwrap();
        let mid: Vec<_> = x
            .iter()
            .zip(&record.x)
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        let velocity = mid[1] / 0.6;
        let low = 0.7 * flow.axial_low_tail_resistance_per_length_n_s_m2();
        let mut force_on_body = force - 4.0 * mid[0] - low * velocity;
        let mut power = low * velocity * velocity;
        for (k, branch) in flow.axial_branches().enumerate() {
            let mass = 0.7 * branch.mass_per_length_kg_m;
            let resistance = 0.7 * branch.resistance_per_length_n_s_m2;
            let relative_velocity = velocity - mid[k + 2] / mass;
            let reaction = resistance * relative_velocity;
            near((record.x[k + 2] - x[k + 2]) / dt, reaction);
            force_on_body -= reaction;
            power += resistance * relative_velocity * relative_velocity;
        }
        near((record.x[1] - x[1]) / dt, force_on_body);
        near((record.x[0] - x[0]) / dt, velocity);
        heat += dt * power;
        work += force * (record.x[0] - x[0]);
        x = record.x;
    }
    assert!(heat > 0.0);
    near(body.system().hamiltonian(&x) - initial + heat, work);
    assert!(x[2..].iter().any(|p| *p != 0.0));
}

#[test]
fn axial_drag_changes_actual_motion_without_inventing_an_added_body_mass() {
    let gate = CancelGate::new_clock_free();
    let ambient = flow(2.0, 0.1, 4);
    let viscous = flow(2.0, 0.4, 4);
    let a = ambient.axial_body(0.7, 0.6, 0.0, 6, &gate).unwrap();
    let b = viscous.axial_body(0.7, 0.6, 0.0, 6, &gate).unwrap();
    assert_eq!(a.moving_mass_kg(), b.moving_mass_kg());
    assert_eq!(a.memory_masses_kg(), b.memory_masses_kg());
    let start = [0.0; 6];
    let ra = fs_phs::step(a.system(), &start, &[1.0], 0.01).unwrap();
    let rb = fs_phs::step(b.system(), &start, &[1.0], 0.01).unwrap();
    assert!(rb.x[0] < ra.x[0]);
    assert!(ra.x[0] < 0.5 / 0.6 * 0.01_f64.powi(2));
    assert!(ra.x[2..].iter().all(|p| *p > 0.0));
}
