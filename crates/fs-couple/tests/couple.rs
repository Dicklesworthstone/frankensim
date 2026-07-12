//! Battery for port-Hamiltonian coupling (fs-couple). Covers power-conjugate
//! ports, the Dirac interconnection's exact power conservation, the energy
//! audit (passivity measured, not assumed), the Aitken relaxation factor, and
//! the load-bearing added-mass comparison: naive staggering diverges where
//! Aitken-relaxed coupling stays stable.

use fs_couple::{
    AitkenRelaxation, CoupleError, EnergyAudit, Port, PortKind, interconnect, interface_power,
    iterate_aitken, iterate_fixed_relaxation,
};

#[test]
fn ports_are_power_conjugate_by_physical_type() {
    let force = Port::new(10.0, 2.0, PortKind::MechanicalForceVelocity);
    let force2 = Port::new(5.0, 1.0, PortKind::MechanicalForceVelocity);
    let pressure = Port::new(3.0, 4.0, PortKind::FluidPressureFlux);
    assert!((force.power() - 20.0).abs() < 1e-12); // effort × flow
    assert!(force.conjugate_to(&force2)); // same physical type
    assert!(!force.conjugate_to(&pressure)); // force can't couple to pressure
}

#[test]
fn the_dirac_interconnection_conserves_interface_power_exactly() {
    let c = interconnect(
        PortKind::MechanicalForceVelocity,
        PortKind::MechanicalForceVelocity,
        7.0,
        3.0,
    )
    .unwrap();
    // shared effort, opposite flow -> net interface power is exactly zero (G0).
    assert!(c.interface_power.abs() < 1e-15);
    assert!((c.port_a.effort - c.port_b.effort).abs() < 1e-15);
    assert!((c.port_a.flow + c.port_b.flow).abs() < 1e-15);
    // incompatible ports are refused at composition time.
    assert!(matches!(
        interconnect(
            PortKind::MechanicalForceVelocity,
            PortKind::FluidPressureFlux,
            1.0,
            1.0
        ),
        Err(CoupleError::IncompatiblePorts { .. })
    ));
}

#[test]
fn the_energy_audit_measures_passivity_and_alarms_on_generation() {
    let mut audit = EnergyAudit::new();
    // a correct interconnection conserves power.
    let good = interconnect(
        PortKind::FluidPressureFlux,
        PortKind::FluidPressureFlux,
        4.0,
        2.0,
    )
    .unwrap();
    audit.record(good.interface_power);
    assert!(audit.is_passive(1e-12));
    // a BROKEN coupling (both ports inject power) generates energy -> alarm.
    let broken = interface_power(&[
        Port::new(2.0, 1.0, PortKind::MechanicalForceVelocity),
        Port::new(2.0, 1.0, PortKind::MechanicalForceVelocity),
    ]);
    audit.record(broken);
    assert!(!audit.is_passive(1e-12));
    assert!((audit.max_generation() - 4.0).abs() < 1e-12);
}

#[test]
fn the_energy_audit_fails_closed_on_a_nan_interface_power() {
    // Regression: a NaN interface power is a hard numerical breakdown — the
    // worst thing the passivity audit exists to flag. `f64::max` drops NaN
    // (`f64::max(0.0, NaN) == 0.0`), so the old fold reported ZERO generation
    // and certified the blown-up coupling as passive — a false certificate.
    let mut audit = EnergyAudit::new();
    audit.record(0.0); // a clean, conserved exchange first
    assert!(audit.is_passive(1e-12), "a conserved exchange is passive");
    audit.record(f64::NAN); // then a diverged exchange
    assert!(
        audit.max_generation().is_nan(),
        "a NaN balance must poison the metric, not vanish"
    );
    assert!(
        !audit.is_passive(1e-12),
        "a NaN interface power must never certify as passive"
    );
    // An arbitrarily large tolerance cannot rescue a NaN, either.
    assert!(!audit.is_passive(f64::INFINITY));
}

#[test]
fn the_aitken_factor_follows_the_delta_squared_formula() {
    let mut a = AitkenRelaxation::new(0.5, 2.0);
    // first call returns the initial ω.
    assert!((a.next_omega(3.0) - 0.5).abs() < 1e-12);
    // ω₁ = −ω₀·r₀/(r₁−r₀) = −0.5·3/(−1.5−3) = 1/3.
    assert!((a.next_omega(-1.5) - 1.0 / 3.0).abs() < 1e-9);
}

#[test]
fn naive_staggering_diverges_where_aitken_stays_stable() {
    // dense fluid on a light structure: added-mass ratio μ = 2 (> 1).
    let (mu, c, x0) = (2.0, 3.0, 0.0);
    // naive Gauss-Seidel staggering (ω = 1) DIVERGES.
    let naive = iterate_fixed_relaxation(mu, c, x0, 1.0, 100, 1e-9);
    assert!(!naive.converged, "naive should diverge, got {naive:?}");
    // Aitken-relaxed coupling CONVERGES to the fixed point x* = c/(1+μ) = 1.
    let aitken = iterate_aitken(mu, c, x0, 0.5, 2.0, 100, 1e-9);
    assert!(aitken.converged);
    assert!((aitken.solution - 1.0).abs() < 1e-6);
    assert!(
        aitken.steps <= 5,
        "Aitken should converge fast, took {}",
        aitken.steps
    );
}

#[test]
fn aitken_accelerates_over_a_stable_fixed_relaxation() {
    let (mu, c, x0) = (2.0, 3.0, 0.0);
    // a stable but slower under-relaxation.
    let fixed = iterate_fixed_relaxation(mu, c, x0, 0.3, 200, 1e-12);
    let aitken = iterate_aitken(mu, c, x0, 0.5, 2.0, 200, 1e-12);
    assert!(fixed.converged && aitken.converged);
    assert!(
        aitken.steps <= fixed.steps,
        "Aitken {} !<= fixed {}",
        aitken.steps,
        fixed.steps
    );
}

#[test]
fn light_added_mass_converges_even_naively() {
    // μ < 1 (heavy structure): naive staggering is already stable.
    let r = iterate_fixed_relaxation(0.5, 3.0, 0.0, 1.0, 100, 1e-9);
    assert!(r.converged);
    assert!((r.solution - 2.0).abs() < 1e-6); // x* = 3/(1+0.5) = 2
}

#[test]
fn coupling_is_deterministic() {
    let a = iterate_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9);
    let b = iterate_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9);
    assert_eq!(a, b);
}
