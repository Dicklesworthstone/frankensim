//! G0/G1/G3: authored Vinen onset, independent direct-time dynamics and energy.
//! The coefficients are numerical fixtures, not helium calibration.

use fs_material::superfluid::vinen::{VinenParameters, VinenStepControls, VinenTangle};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState, TwoFluidVelocity};

fn material() -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: 0.0,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap()
}

fn parameters() -> VinenParameters {
    VinenParameters {
        circulation_quantum_m2_s: 3.0,
        dissipative_b: 1.0,
        production_coefficient: 0.8,
        decay_coefficient_m2_s: 0.5,
        line_energy_j_m: 0.05,
    }
}

fn model(onset: f64, fraction: f64, budget: usize) -> VinenTangle {
    VinenTangle::try_new(
        material(),
        parameters(),
        VinenStepControls {
            maximum_fractional_step: fraction,
            maximum_substeps: budget,
        },
    )
    .unwrap()
    .with_counterflow_onset_m_s(onset)
    .unwrap()
}

fn velocity(relative: f64) -> TwoFluidVelocity {
    material()
        .velocities_from_barycentric([0.25, -0.125, 0.5], [relative, 0.0, 0.0])
        .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}

// Independent physical-time RK4 of the stated nonlinear ODE. It does not use
// exposure, frozen-coefficient exponentials or the production implementation.
fn reference(mut w: f64, mut l: f64, onset: f64, duration: f64) -> (f64, f64) {
    let h = duration / 32768.0;
    let rhs = |w: f64, l: f64| {
        (
            -l * w,
            0.8 * (w - onset).max(0.0) * l * l.sqrt() - 0.5 * l * l,
        )
    };
    for _ in 0..32768 {
        let k1 = rhs(w, l);
        let k2 = rhs(w + h * k1.0 / 2.0, l + h * k1.1 / 2.0);
        let k3 = rhs(w + h * k2.0 / 2.0, l + h * k2.1 / 2.0);
        let k4 = rhs(w + h * k3.0, l + h * k3.1);
        w += h * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0) / 6.0;
        l += h * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1) / 6.0;
    }
    (w, l)
}

#[test]
fn zero_offset_is_the_original_phase_curve_bit_for_bit() {
    let original =
        VinenTangle::try_new(material(), parameters(), model(0.0, 0.05, 512).controls()).unwrap();
    for w in [0.0, 0.1, 2.0, 4.0] {
        let a = original.advance(velocity(w), 1.3, 0.2).unwrap();
        let b = original
            .with_counterflow_onset_m_s(0.0)
            .unwrap()
            .advance(velocity(w), 1.3, 0.2)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(
            a.line_density_per_m2.to_bits(),
            b.line_density_per_m2.to_bits()
        );
        assert_eq!(
            a.dissipated_energy_j_m3.to_bits(),
            b.dissipated_energy_j_m3.to_bits()
        );
    }
}

#[test]
fn subcritical_remanence_uses_exact_decay_and_keeps_drag() {
    let tangle = model(1.0, 1e-12, 1);
    for w in [0.0, 0.5, 1.0] {
        let initial = velocity(w);
        let step = tangle.advance(initial, 4.0, 0.75).unwrap();
        close(step.line_density_per_m2, 1.6, 2e-14);
        close(
            material().counterflow_velocity_m_s(step.velocity).unwrap()[0],
            w / 6.25,
            2e-15,
        );
        assert_eq!(step.substeps, 1);
        close(step.energy_balance_defect_j_m3, 0.0, 2e-15);
        assert!(step.dissipated_energy_j_m3 >= 0.12 - 2e-15);
    }
}

#[test]
fn supercritical_growth_has_the_authored_excess_speed_derivative() {
    let initial = velocity(2.0);
    let dt = 1e-7;
    let step = model(0.75, 0.05, 64).advance(initial, 1.0, dt).unwrap();
    // L'=0.8*(2-0.75)-0.5 = 0.5; W'=-L*W=-2.
    close((step.line_density_per_m2 - 1.0) / dt, 0.5, 3e-7);
    let w = material().counterflow_velocity_m_s(step.velocity).unwrap()[0];
    close((w - 2.0) / dt, -2.0, 3e-7);
    assert!(step.dissipated_energy_j_m3 < step.counterflow_energy_loss_j_m3);
}

#[test]
fn smooth_source_converges_at_second_order_against_direct_time_rk4() {
    let expected = reference(3.0, 1.0, 0.4, 0.5);
    let error = |fraction| {
        let step = model(0.4, fraction, 4096)
            .advance(velocity(3.0), 1.0, 0.5)
            .unwrap();
        let w = material().counterflow_velocity_m_s(step.velocity).unwrap()[0];
        (w - expected.0).abs() + (step.line_density_per_m2 - expected.1).abs()
    };
    let coarse = error(0.08);
    let medium = error(0.04);
    let fine = error(0.02);
    assert!(medium < 0.32 * coarse, "{coarse:e}, {medium:e}, {fine:e}");
    assert!(fine < 0.32 * medium, "{coarse:e}, {medium:e}, {fine:e}");
}

#[test]
fn crossing_onset_retains_momentum_energy_and_converges_to_decay() {
    let expected = reference(2.0, 1.0, 1.1, 2.0);
    assert!(expected.0 < 1.1);
    let initial = velocity(2.0);
    let step = model(1.1, 0.002, 4096).advance(initial, 1.0, 2.0).unwrap();
    let w = material().counterflow_velocity_m_s(step.velocity).unwrap()[0];
    close(w, expected.0, 2e-6);
    close(step.line_density_per_m2, expected.1, 2e-6);
    assert!(step.dissipated_energy_j_m3 > 0.0);
    let before = material()
        .kinetic_energy_density(initial)
        .unwrap()
        .total_j_m3
        + 0.05;
    let after = material()
        .kinetic_energy_density(step.velocity)
        .unwrap()
        .total_j_m3
        + step.vortex_energy_after_j_m3;
    close(before - after, step.dissipated_energy_j_m3, 2e-14);
    close(step.energy_balance_defect_j_m3, 0.0, 2e-14);
    let mean = material().barycentric_velocity_m_s(step.velocity).unwrap();
    for (got, old) in mean.into_iter().zip([0.25, -0.125, 0.5]) {
        close(got, old, 2e-15);
    }
    for i in 0..3 {
        assert_eq!(
            step.normal_impulse_kg_m2_s[i],
            -step.superfluid_impulse_kg_m2_s[i]
        );
    }
}

#[test]
fn imposed_flow_equilibrium_and_zero_seed_are_distinct() {
    let tangle = model(0.75, 0.05, 512);
    assert_eq!(tangle.counterflow_onset_m_s(), 0.75);
    for speed in [0.0, 0.5, 0.75] {
        assert_eq!(
            tangle
                .imposed_counterflow_equilibrium_per_m2(speed)
                .unwrap(),
            0.0
        );
    }
    close(
        tangle.imposed_counterflow_equilibrium_per_m2(2.0).unwrap(),
        4.0,
        1e-15,
    );
    let initial = velocity(2.0);
    let step = tangle.advance(initial, 0.0, 1e6).unwrap();
    assert_eq!(step.velocity, initial);
    assert_eq!(step.line_density_per_m2, 0.0);
    assert_eq!(step.dissipated_energy_j_m3, 0.0);
}

#[test]
fn invalid_onset_and_work_exhaustion_refuse() {
    let tangle = model(0.5, 0.05, 512);
    for bad in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(tangle.with_counterflow_onset_m_s(bad).is_err());
        assert!(tangle.imposed_counterflow_equilibrium_per_m2(bad).is_err());
    }
    assert!(
        model(0.5, 1e-6, 1)
            .advance(velocity(3.0), 1.0, 1.0)
            .is_err()
    );
}
