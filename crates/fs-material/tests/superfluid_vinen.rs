//! G0/G1/G3: independent Vinen source equations, line storage and frame laws.
//! Coefficients below are authored numerical fixtures, not helium calibration.

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

fn controls() -> VinenStepControls {
    VinenStepControls {
        maximum_fractional_step: 0.125,
        maximum_substeps: 4096,
    }
}

fn model(p: VinenParameters) -> VinenTangle {
    VinenTangle::try_new(material(), p, controls()).unwrap()
}

fn velocity(mean: [f64; 3], relative: [f64; 3]) -> TwoFluidVelocity {
    material()
        .velocities_from_barycentric(mean, relative)
        .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tol={tolerance:e}"
    );
}

#[test]
fn exact_decay_advances_drag_with_the_same_changing_density_and_releases_line_energy() {
    let mut p = parameters();
    p.production_coefficient = 0.0;
    let tangle = VinenTangle::try_new(
        material(),
        p,
        VinenStepControls {
            maximum_fractional_step: 1e-12,
            maximum_substeps: 1,
        },
    )
    .unwrap();
    let initial = velocity([0.2, -0.3, 0.1], [0.8, -0.6, 0.0]);
    let step = tangle.advance(initial, 4.0, 0.75).unwrap();
    // L=L0/(1+bL0t); w=w0*(1+bL0t)^(-c/b), here c/b=2.
    close(step.line_density_per_m2, 1.6, 2e-14);
    let relative = material().counterflow_velocity_m_s(step.velocity).unwrap();
    for (got, old) in relative.into_iter().zip([0.8, -0.6, 0.0]) {
        close(got, old / 6.25, 2e-15);
    }
    let released = 0.5 * 0.75 * (1.0 - 1.0 / 6.25_f64.powi(2));
    close(step.counterflow_energy_loss_j_m3, released, 2e-15);
    close(
        step.dissipated_energy_j_m3,
        released + 0.05 * (4.0 - 1.6),
        2e-15,
    );
    close(step.energy_balance_defect_j_m3, 0.0, 2e-15);
    assert_eq!(step.substeps, 1);
}

#[test]
fn coflow_decays_seeded_lines_without_accelerating_either_component() {
    let initial = velocity([0.25, -0.5, 0.75], [0.0; 3]);
    let step = model(parameters()).advance(initial, 8.0, 2.0).unwrap();
    assert_eq!(step.velocity, initial);
    assert_eq!(step.normal_impulse_kg_m2_s, [0.0; 3]);
    assert_eq!(step.counterflow_energy_loss_j_m3, 0.0);
    close(step.line_density_per_m2, 8.0 / 9.0, 2e-14);
    close(step.dissipated_energy_j_m3, 0.05 * 64.0 / 9.0, 2e-15);
    assert_eq!(step.substeps, 1);
}

#[test]
fn growth_is_funded_by_counterflow_and_heat_is_not_double_counted() {
    let initial = velocity([0.2, 0.0, 0.0], [2.0, 0.0, 0.0]);
    let tangle = model(parameters());
    let step = tangle.advance(initial, 1.0, 0.02).unwrap();
    assert!(step.line_density_per_m2 > 1.0);
    assert!(step.dissipated_energy_j_m3 > 0.0);
    assert!(step.dissipated_energy_j_m3 < step.counterflow_energy_loss_j_m3);
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
    close(before, after + step.dissipated_energy_j_m3, 2e-15);
    close(step.energy_balance_defect_j_m3, 0.0, 2e-15);
    close(
        material().barycentric_velocity_m_s(step.velocity).unwrap()[0],
        0.2,
        2e-16,
    );
    for i in 0..3 {
        assert_eq!(
            step.superfluid_impulse_kg_m2_s[i],
            -step.normal_impulse_kg_m2_s[i]
        );
        close(
            step.normal_impulse_kg_m2_s[i],
            material().normal_density_kg_m3()
                * (step.velocity.normal_m_s[i] - initial.normal_m_s[i]),
            3e-16,
        );
    }
}

#[test]
fn zero_seed_and_zero_duration_are_exact_identities() {
    let initial = TwoFluidVelocity {
        normal_m_s: [1.0, -0.0, 0.0],
        superfluid_m_s: [-1.0, 0.0, -0.0],
    };
    for (density, duration) in [(0.0, 1e9), (3.0, 0.0)] {
        let step = model(parameters())
            .advance(initial, density, duration)
            .unwrap();
        assert_eq!(
            step.velocity.normal_m_s.map(f64::to_bits),
            initial.normal_m_s.map(f64::to_bits)
        );
        assert_eq!(
            step.velocity.superfluid_m_s.map(f64::to_bits),
            initial.superfluid_m_s.map(f64::to_bits)
        );
        assert_eq!(step.line_density_per_m2, density);
        assert_eq!(step.dissipated_energy_j_m3, 0.0);
        assert_eq!(step.substeps, 0);
    }
}

// Independent direct-time fourth-order integration of the ORIGINAL nonlinear
// differential equations, not the production exposure transformation.
fn rk4(mut w: f64, mut l: f64, time: f64, n: usize, p: VinenParameters) -> (f64, f64) {
    let h = time / n as f64;
    let rhs = |w: f64, l: f64| {
        (
            -p.dissipative_b * p.circulation_quantum_m2_s / 3.0 * l * w,
            p.production_coefficient * w.abs() * l * l.sqrt() - p.decay_coefficient_m2_s * l * l,
        )
    };
    for _ in 0..n {
        let k1 = rhs(w, l);
        let k2 = rhs(w + 0.5 * h * k1.0, l + 0.5 * h * k1.1);
        let k3 = rhs(w + 0.5 * h * k2.0, l + 0.5 * h * k2.1);
        let k4 = rhs(w + h * k3.0, l + h * k3.1);
        w += h / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
        l += h / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
    }
    (w, l)
}

#[test]
fn source_clock_has_second_order_convergence_against_independent_nonlinear_ode() {
    let p = parameters();
    let tangle = VinenTangle::try_new(
        material(),
        p,
        VinenStepControls {
            maximum_fractional_step: 0.25,
            maximum_substeps: 4096,
        },
    )
    .unwrap();
    let reference = rk4(1.0, 0.3, 1.0, 20000, p);
    let solve = |n: usize| {
        let mut v = velocity([0.0; 3], [1.0, 0.0, 0.0]);
        let mut l = 0.3;
        let mut heat = 0.0;
        for _ in 0..n {
            let step = tangle.advance(v, l, 1.0 / n as f64).unwrap();
            assert_eq!(step.substeps, 1);
            heat += step.dissipated_energy_j_m3;
            v = step.velocity;
            l = step.line_density_per_m2;
        }
        let w = material().counterflow_velocity_m_s(v).unwrap()[0];
        close(
            0.5 * 0.75 + 0.05 * 0.3,
            0.5 * 0.75 * w * w + 0.05 * l + heat,
            2e-14,
        );
        ((w - reference.0).abs() + (l - reference.1).abs(), w, l)
    };
    let coarse = solve(8);
    let medium = solve(16);
    let fine = solve(32);
    assert!(
        (3.5..4.5).contains(&(coarse.0 / medium.0)),
        "{coarse:?} {medium:?}"
    );
    assert!(
        (3.5..4.5).contains(&(medium.0 / fine.0)),
        "{medium:?} {fine:?}"
    );
    close(fine.1, reference.0, 2e-6);
    close(fine.2, reference.1, 2e-6);
}

#[test]
fn equal_and_nearly_equal_exposure_rates_have_the_same_continuous_limit() {
    let mut p = parameters();
    p.circulation_quantum_m2_s = 0.75; // c=b/2=0.25
    let initial = velocity([0.0; 3], [1.5, 0.0, 0.0]);
    let center = model(p).advance(initial, 0.4, 0.5).unwrap();
    let tau = center.integrated_line_density_s_m2;
    let expected =
        (0.4_f64.sqrt() + 0.5 * p.production_coefficient * 1.5 * tau) * (-0.25 * tau).exp();
    close(center.line_density_per_m2, expected * expected, 3e-15);
    for offset in [-1e-12, 1e-12] {
        p.circulation_quantum_m2_s = 0.75 * (1.0 + offset);
        let nearby = model(p).advance(initial, 0.4, 0.5).unwrap();
        close(
            nearby.line_density_per_m2,
            center.line_density_per_m2,
            1e-12,
        );
        close(
            nearby.dissipated_energy_j_m3,
            center.dissipated_energy_j_m3,
            1e-12,
        );
    }
}

#[test]
fn weak_and_stiff_decay_retain_positive_tails_and_report_the_same_energy_store() {
    let mut p = parameters();
    p.production_coefficient = 0.0;
    let tangle = model(p);
    let initial = velocity([0.0; 3], [1.0, 0.0, 0.0]);
    let weak = tangle.advance(initial, 1.0, 1e-18).unwrap();
    assert!(weak.normal_impulse_kg_m2_s[0] < 0.0);
    assert!(weak.dissipated_energy_j_m3 > 0.0);
    close(weak.dissipated_energy_j_m3, 0.775e-18, 1e-32);
    let stiff = tangle.advance(initial, 1e8, 1e4).unwrap();
    assert!(stiff.line_density_per_m2 > 0.0);
    assert!(stiff.velocity.normal_m_s[0] > 0.0);
    close(stiff.line_density_per_m2, 1e8 / (1.0 + 5e11), 1e-17);
    close(
        stiff.dissipated_energy_j_m3,
        0.375 + 0.05 * (1e8 - stiff.line_density_per_m2),
        2e-8,
    );
}

#[test]
fn scaled_decay_and_drag_do_not_erase_representable_weak_or_late_responses() {
    let mut p = parameters();
    p.production_coefficient = 0.0;
    p.decay_coefficient_m2_s = 1e-300;
    let weak = model(p)
        .advance(velocity([0.0; 3], [1.0, 0.0, 0.0]), 1.0, 1e-100)
        .unwrap();
    close(weak.integrated_line_density_s_m2, 1e-100, 1e-115);
    close(weak.normal_impulse_kg_m2_s[0], -0.75e-100, 2e-115);
    close(weak.dissipated_energy_j_m3, 0.75e-100, 2e-115);

    p = parameters();
    p.production_coefficient = 0.0;
    p.circulation_quantum_m2_s = 1200.0; // c=400; exposure=2 -> drag=800.
    let late = model(p)
        .advance(
            velocity([0.0; 3], [1e150, 0.0, 0.0]),
            1.0,
            2.0 * (1.0_f64.exp() - 1.0),
        )
        .unwrap();
    let relative = material().counterflow_velocity_m_s(late.velocity).unwrap()[0];
    let expected = (1e150 * (-400.0_f64).exp()) * (-400.0_f64).exp();
    assert!(relative > 0.0);
    close(relative / expected, 1.0, 5e-12);
}

#[test]
fn boosts_and_rotations_preserve_tangle_history_heat_and_opposite_impulses() {
    let tangle = model(parameters());
    let initial = velocity([0.0; 3], [0.8, -0.6, 0.0]);
    let base = tangle.advance(initial, 0.5, 0.3).unwrap();
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let boost = [0.3, -0.1, 0.4];
    let transform = |v: [f64; 3]| std::array::from_fn(|i| rotate(v)[i] + boost[i]);
    let changed = tangle
        .advance(
            TwoFluidVelocity {
                normal_m_s: transform(initial.normal_m_s),
                superfluid_m_s: transform(initial.superfluid_m_s),
            },
            0.5,
            0.3,
        )
        .unwrap();
    close(changed.line_density_per_m2, base.line_density_per_m2, 1e-15);
    close(
        changed.dissipated_energy_j_m3,
        base.dissipated_energy_j_m3,
        1e-15,
    );
    for i in 0..3 {
        close(
            changed.velocity.normal_m_s[i],
            transform(base.velocity.normal_m_s)[i],
            2e-15,
        );
        close(
            changed.normal_impulse_kg_m2_s[i],
            rotate(base.normal_impulse_kg_m2_s)[i],
            2e-15,
        );
    }
}

#[test]
fn thermodynamically_incompatible_coefficients_and_exhausted_budgets_refuse() {
    let mut p = parameters();
    p.production_coefficient = 100.0;
    assert!(VinenTangle::try_new(material(), p, controls()).is_err());
    p = parameters();
    p.dissipative_b = 0.0;
    assert!(VinenTangle::try_new(material(), p, controls()).is_err());
    p.production_coefficient = 0.0;
    let decaying = VinenTangle::try_new(material(), p, controls()).unwrap();
    let initial = velocity([0.0; 3], [1.0, 0.0, 0.0]);
    let decay = decaying.advance(initial, 1.0, 1.0).unwrap();
    assert_eq!(decay.velocity, initial);
    assert!(decay.dissipated_energy_j_m3 > 0.0);
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(decaying.advance(initial, bad, 1.0).is_err());
        assert!(decaying.advance(initial, 1.0, bad).is_err());
    }
    let limited = VinenTangle::try_new(
        material(),
        parameters(),
        VinenStepControls {
            maximum_fractional_step: 1e-6,
            maximum_substeps: 1,
        },
    )
    .unwrap();
    let saved = limited;
    assert!(
        limited
            .advance(initial, 1.0, 1.0)
            .unwrap_err()
            .to_string()
            .contains("maximum_substeps")
    );
    assert_eq!(limited, saved);
}
