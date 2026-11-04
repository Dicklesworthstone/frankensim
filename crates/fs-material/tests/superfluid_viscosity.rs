//! G1/G3 checks of the frozen normal-wall ODE, its physical accounts and limits.
//! Authored numerical fixtures do not qualify helium transport measurements.

use fs_material::superfluid::{
    NormalViscousWall, TwoFluidProperties, TwoFluidState, TwoFluidVelocity,
};

fn properties() -> TwoFluidProperties {
    TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 12.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: 0.0,
        normal_viscosity_pa_s: 0.2,
    }
}

fn state() -> TwoFluidState {
    TwoFluidState::try_new(properties()).unwrap()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 128.0 * f64::EPSILON * scale,
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

fn close_vector(actual: [f64; 3], expected: [f64; 3], scale: f64) {
    for (a, e) in actual.into_iter().zip(expected) {
        close(a, e, scale);
    }
}

fn velocity() -> TwoFluidVelocity {
    TwoFluidVelocity {
        normal_m_s: [0.6, -0.8, 0.3],
        superfluid_m_s: [1.0, -0.25, -0.0],
    }
}

#[test]
fn permeability_helpers_recover_poiseuille_resistance_with_normal_inertia() {
    let material = state();
    let circle = NormalViscousWall::circular_capillary(2.0).unwrap();
    let slit = NormalViscousWall::parallel_plates(3.0).unwrap();
    assert_eq!(circle.permeability_m2(), 0.5);
    assert_eq!(slit.permeability_m2(), 0.75);
    // rho_n=1, rho=4, rho_s=3. Neither total nor relative inertia belongs here.
    close(
        material.normal_wall_decay_rate_s_inv(circle).unwrap(),
        0.4,
        1.0,
    );
    close(
        material.normal_wall_decay_rate_s_inv(slit).unwrap(),
        4.0 / 15.0,
        1.0,
    );
    // Independently integrate the circular parabolic profile and its wall
    // shear: <u> = G*R²/(8*eta), <force> = -(2/R)*eta*|du/dr|_wall = -G.
    let radius: f64 = 0.006;
    let gradient: f64 = 120.0;
    let eta = material.normal_viscosity_pa_s();
    let mean = gradient * radius.powi(2) / (8.0 * eta);
    let law = NormalViscousWall::circular_capillary(radius).unwrap();
    let dt = 1e-12;
    let result = material
        .relax_normal_wall(
            TwoFluidVelocity {
                normal_m_s: [mean, 0.0, 0.0],
                superfluid_m_s: [0.0; 3],
            },
            law,
            dt,
        )
        .unwrap();
    let observed_force = result.normal_impulse_kg_m2_s[0] / dt;
    assert!((observed_force / -gradient - 1.0).abs() < 3e-8);
}

#[test]
fn exact_normal_decay_routes_external_wall_impulse_and_heat_once() {
    let material = state();
    let wall = NormalViscousWall::try_new(0.1).unwrap();
    let initial = velocity();
    let result = material.relax_normal_wall(initial, wall, 0.3).unwrap();
    let decay = (-0.6_f64).exp();
    close_vector(
        result.velocity.normal_m_s,
        initial.normal_m_s.map(|v| v * decay),
        1.0,
    );
    assert_eq!(
        result.velocity.superfluid_m_s.map(f64::to_bits),
        initial.superfluid_m_s.map(f64::to_bits)
    );
    assert_eq!(
        result.wall_impulse_kg_m2_s,
        result.normal_impulse_kg_m2_s.map(|v| -v)
    );
    let before = material.currents(initial).unwrap().mass_kg_m2_s;
    let after = material.currents(result.velocity).unwrap().mass_kg_m2_s;
    for i in 0..3 {
        close(after[i] - before[i], result.normal_impulse_kg_m2_s[i], 4.0);
        close(
            result.normal_impulse_kg_m2_s[i],
            initial.normal_m_s[i] * (decay - 1.0),
            1.0,
        );
    }
    let normal_energy = 0.5 * (0.6_f64.powi(2) + 0.8_f64.powi(2) + 0.3_f64.powi(2));
    close(
        result.dissipated_energy_j_m3,
        normal_energy * (1.0 - (-1.2_f64).exp()),
        normal_energy,
    );
    close(
        material.kinetic_energy_density(initial).unwrap().total_j_m3,
        material
            .kinetic_energy_density(result.velocity)
            .unwrap()
            .total_j_m3
            + result.dissipated_energy_j_m3,
        4.0,
    );
}

#[test]
fn prepared_source_is_reusable_and_composes_as_a_semigroup() {
    let material = state();
    let wall = NormalViscousWall::try_new(0.1).unwrap();
    let initial = velocity();
    let prepared = material.prepare_normal_wall_relaxation(wall, 0.3).unwrap();
    assert_eq!(prepared.normal_density_kg_m3(), 1.0);
    assert_eq!(prepared.decay_rate_s_inv(), 2.0);
    assert_eq!(prepared.duration_s(), 0.3);
    let whole = prepared.apply(initial).unwrap();
    assert_eq!(
        whole,
        material.relax_normal_wall(initial, wall, 0.3).unwrap()
    );
    assert_eq!(prepared.apply(initial).unwrap(), whole);
    let first = material.relax_normal_wall(initial, wall, 0.1).unwrap();
    let split = material
        .relax_normal_wall(first.velocity, wall, 0.2)
        .unwrap();
    close_vector(split.velocity.normal_m_s, whole.velocity.normal_m_s, 1.0);
    close(
        first.dissipated_energy_j_m3 + split.dissipated_energy_j_m3,
        whole.dissipated_energy_j_m3,
        1.0,
    );
    for i in 0..3 {
        close(
            first.wall_impulse_kg_m2_s[i] + split.wall_impulse_kg_m2_s[i],
            whole.wall_impulse_kg_m2_s[i],
            1.0,
        );
    }
    let rotate = |v: [f64; 3]| [-v[1], v[0], v[2]];
    let rotated = prepared
        .apply(TwoFluidVelocity {
            normal_m_s: rotate(initial.normal_m_s),
            superfluid_m_s: rotate(initial.superfluid_m_s),
        })
        .unwrap();
    close_vector(
        rotated.velocity.normal_m_s,
        rotate(whole.velocity.normal_m_s),
        1.0,
    );
    close(
        rotated.dissipated_energy_j_m3,
        whole.dissipated_energy_j_m3,
        1.0,
    );
}

#[test]
fn walls_damp_coflow_but_leave_pure_superflow_unchanged() {
    let material = state();
    let wall = NormalViscousWall::try_new(0.1).unwrap();
    let coflow = TwoFluidVelocity {
        normal_m_s: [1.0, -2.0, 0.5],
        superfluid_m_s: [1.0, -2.0, 0.5],
    };
    let step = material.relax_normal_wall(coflow, wall, 0.2).unwrap();
    assert!(step.dissipated_energy_j_m3 > 0.0);
    assert_ne!(step.velocity.normal_m_s, step.velocity.superfluid_m_s);
    assert_eq!(step.velocity.superfluid_m_s, coflow.superfluid_m_s);
    let superflow = TwoFluidVelocity {
        normal_m_s: [0.0; 3],
        ..coflow
    };
    let unchanged = material.relax_normal_wall(superflow, wall, 1.0).unwrap();
    assert_eq!(unchanged.velocity, superflow);
    assert_eq!(unchanged.dissipated_energy_j_m3, 0.0);
    assert_eq!(unchanged.normal_impulse_kg_m2_s, [0.0; 3]);
}

#[test]
fn short_and_stiff_intervals_preserve_impulse_heat_and_exponential_tails() {
    let material = state();
    let wall = NormalViscousWall::try_new(0.1).unwrap();
    let initial = TwoFluidVelocity {
        normal_m_s: [1.0, 0.0, 0.0],
        ..velocity()
    };
    let tiny = material.relax_normal_wall(initial, wall, 1e-30).unwrap();
    close(tiny.normal_impulse_kg_m2_s[0], -2e-30, 2e-30);
    close(tiny.dissipated_energy_j_m3, 2e-30, 2e-30);
    let tail = material.relax_normal_wall(initial, wall, 40.0).unwrap();
    close(
        tail.velocity.normal_m_s[0],
        (-80.0_f64).exp(),
        (-80.0_f64).exp(),
    );
    assert_eq!(tail.dissipated_energy_j_m3, 0.5);
    let stiff = material
        .relax_normal_wall(initial, wall, f64::MAX / 2.0)
        .unwrap();
    assert_eq!(stiff.velocity.normal_m_s, [0.0; 3]);
    assert_eq!(stiff.velocity.superfluid_m_s, initial.superfluid_m_s);
    assert_eq!(stiff.dissipated_energy_j_m3, 0.5);
    assert_eq!(stiff.wall_impulse_kg_m2_s, [1.0, 0.0, 0.0]);
    assert_eq!(
        material
            .relax_normal_wall(initial, wall, 0.0)
            .unwrap()
            .velocity,
        initial
    );
    let mut inviscid = properties();
    inviscid.normal_viscosity_pa_s = 0.0;
    let inviscid = TwoFluidState::try_new(inviscid).unwrap();
    let zero = inviscid.relax_normal_wall(initial, wall, f64::MAX).unwrap();
    assert_eq!(zero.velocity, initial);
    assert_eq!(zero.dissipated_energy_j_m3, 0.0);
    assert_eq!(zero.normal_impulse_kg_m2_s, [0.0; 3]);
}

#[test]
fn rate_admission_handles_representable_scales_and_refuses_overflow() {
    let large_circle = NormalViscousWall::circular_capillary(3e154).unwrap();
    close(large_circle.permeability_m2(), 1.125e308, 1.125e308);
    let mut p = properties();
    p.density_kg_m3 = 1e100;
    p.normal_viscosity_pa_s = 1e-300;
    let material = TwoFluidState::try_new(p).unwrap();
    assert_eq!(material.normal_kinematic_viscosity_m2_s(), 0.0);
    close(
        material
            .normal_wall_decay_rate_s_inv(NormalViscousWall::try_new(1e-200).unwrap())
            .unwrap(),
        4e-200,
        4e-200,
    );
    p.normal_viscosity_pa_s = 1e300;
    let material = TwoFluidState::try_new(p).unwrap();
    close(
        material
            .normal_wall_decay_rate_s_inv(NormalViscousWall::try_new(1e300).unwrap())
            .unwrap(),
        4e-100,
        4e-100,
    );
    p.density_kg_m3 = 1e-200;
    p.normal_viscosity_pa_s = 1e-300;
    let material = TwoFluidState::try_new(p).unwrap();
    close(
        material
            .normal_wall_decay_rate_s_inv(NormalViscousWall::try_new(1e-200).unwrap())
            .unwrap(),
        4e100,
        4e100,
    );
    // A directly formed rho_n*k would be subnormal and lose about 1% here.
    close(
        material
            .normal_wall_decay_rate_s_inv(NormalViscousWall::try_new(1e-122).unwrap())
            .unwrap(),
        4e22,
        4e22,
    );
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(NormalViscousWall::try_new(bad).is_err());
        assert!(NormalViscousWall::circular_capillary(bad).is_err());
        assert!(NormalViscousWall::parallel_plates(bad).is_err());
    }
    assert!(NormalViscousWall::circular_capillary(f64::MAX).is_err());
    assert!(NormalViscousWall::parallel_plates(f64::from_bits(1)).is_err());
    let material = state();
    let wall = NormalViscousWall::try_new(0.1).unwrap();
    for bad in [-1.0, f64::INFINITY, f64::NAN, f64::MAX] {
        assert!(material.prepare_normal_wall_relaxation(wall, bad).is_err());
    }
    assert!(
        material
            .normal_wall_decay_rate_s_inv(NormalViscousWall::try_new(f64::from_bits(1)).unwrap())
            .is_err()
    );
    let bad_velocity = TwoFluidVelocity {
        normal_m_s: [f64::NAN; 3],
        ..velocity()
    };
    assert!(material.relax_normal_wall(bad_velocity, wall, 0.0).is_err());
    let mut p = properties();
    p.normal_viscosity_pa_s = 1e-300;
    let material = TwoFluidState::try_new(p).unwrap();
    assert!(
        material
            .normal_wall_decay_rate_s_inv(NormalViscousWall::try_new(1e100).unwrap())
            .is_err()
    );
}
