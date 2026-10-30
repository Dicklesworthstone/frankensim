//! G1/G3: one Stokes operator for gas and liquid transport, without an EOS fiction.
use super::*;
use fs_material::{MaterialError, gas::GasSpec};
use std::cell::Cell;

fn resolution() -> CylinderFlowResolution {
    CylinderFlowResolution {
        minimum_wavenumber_radius: 1e-5,
        maximum_wavenumber_radius: 3000.0,
        intervals: 512,
    }
}
fn liquid(density: f64, viscosity: f64) -> FluidTransportState {
    // Explicit numerical test data, not a sourced water/material claim.
    FluidTransportState::try_new(density, 4000.0, 0.5, viscosity).unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 3e-13 * expected.abs(),
        "actual={actual:e}, expected={expected:e}");
}

#[test]
fn gas_entry_point_and_shared_transport_have_identical_numerics() {
    let gate = CancelGate::new_clock_free();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let old = CylinderFlow::from_gas(0.001, gas, resolution(), &gate).unwrap();
    let generic = CylinderFlow::from_fluid(0.001, gas, resolution(), &gate).unwrap();
    let tuple = CylinderFlow::from_fluid(0.001, gas.transport_state().unwrap(), resolution(), &gate).unwrap();
    assert_eq!(old, generic);
    assert_eq!(old.gas(), &gas);
    assert_eq!(old.fluid(), &gas);
    assert_eq!(old.branches(), tuple.branches());
    assert_eq!(old.added_mass_per_length_kg_m(), tuple.added_mass_per_length_kg_m());
    assert_eq!(old.low_tail_resistance_per_length_n_s_m2(), tuple.low_tail_resistance_per_length_n_s_m2());
    for omega in [0.1, 1.0, 100.0, 10000.0] {
        assert_eq!(old.impedance_per_length(omega).unwrap(), tuple.impedance_per_length(omega).unwrap());
    }
}

#[test]
fn liquid_force_matches_independent_complex_bessel_solution() {
    let source = liquid(1000.0, 0.001);
    // rho*a^2/mu = 1 s. Independent complex K0/K1 reference, not quadrature.
    let radius = source.kinematic_viscosity_m2_s().sqrt();
    let flow = CylinderFlow::from_fluid(radius, source, resolution(), &CancelGate::new_clock_free()).unwrap();
    let oracle = [
        (0.01, 4.710527336547042, 1.4674730277805959),
        (0.1, 7.489342045271669, 3.6880439694970404),
        (1.0, 14.346304543275211, 12.465886293099675),
        (10.0, 34.05643752523019, 59.754434003390436),
        (100.0, 95.03108417164636, 403.11396197909926),
    ];
    assert_eq!(flow.fluid(), &source);
    assert_eq!(flow.transport(), &source);
    for (omega, real, imag) in oracle {
        let z = flow.impedance_per_length(omega).unwrap();
        let error = (z.re / 0.001 - real).hypot(z.im / 0.001 - imag);
        assert!(error / real.hypot(imag) < 9e-4);
        assert!(z.re > 0.0 && z.im > flow.added_mass_per_length_kg_m() * omega);
    }
}

#[test]
fn density_and_viscosity_change_every_physical_branch_with_correct_scaling() {
    let gate = CancelGate::new_clock_free();
    let a = CylinderFlow::from_fluid(0.001, liquid(1000.0, 0.001), resolution(), &gate).unwrap();
    let dense = CylinderFlow::from_fluid(0.001, liquid(4000.0, 0.001), resolution(), &gate).unwrap();
    let viscous = CylinderFlow::from_fluid(0.001, liquid(1000.0, 0.004), resolution(), &gate).unwrap();
    close(dense.added_mass_per_length_kg_m(), 4.0 * a.added_mass_per_length_kg_m());
    assert_eq!(viscous.added_mass_per_length_kg_m(), a.added_mass_per_length_kg_m());
    close(viscous.low_tail_resistance_per_length_n_s_m2(), 4.0 * a.low_tail_resistance_per_length_n_s_m2());
    for ((base, density), viscosity) in a.branches().iter().zip(dense.branches()).zip(viscous.branches()) {
        close(density.rate_per_s, base.rate_per_s / 4.0);
        close(density.mass_per_length_kg_m, 4.0 * base.mass_per_length_kg_m);
        assert_eq!(density.resistance_per_length_n_s_m2, base.resistance_per_length_n_s_m2);
        close(viscosity.rate_per_s, 4.0 * base.rate_per_s);
        close(viscosity.mass_per_length_kg_m, base.mass_per_length_kg_m);
        close(viscosity.resistance_per_length_n_s_m2, 4.0 * base.resistance_per_length_n_s_m2);
    }
    for omega in [0.01, 1.0, 100.0] {
        let z = a.impedance_per_length(omega).unwrap();
        let d = dense.impedance_per_length(omega / 4.0).unwrap();
        let v = viscous.impedance_per_length(omega * 4.0).unwrap();
        close(d.re, z.re);
        close(d.im, z.im);
        close(v.re, 4.0 * z.re);
        close(v.im, 4.0 * z.im);
    }
}

#[test]
fn thermal_properties_do_not_invent_incompressible_mechanical_loss() {
    let gate = CancelGate::new_clock_free();
    let a = CylinderFlow::from_fluid(0.001, liquid(1000.0, 0.001), resolution(), &gate).unwrap();
    let b = CylinderFlow::from_fluid(0.001,
        FluidTransportState::try_new(1000.0, 2000.0, 2.0, 0.001).unwrap(), resolution(), &gate).unwrap();
    assert_ne!(a.transport(), b.transport());
    assert_eq!(a.branches(), b.branches());
    assert_eq!(a.impedance_per_length(10.0).unwrap(), b.impedance_per_length(10.0).unwrap());
}

#[test]
fn incoherent_gas_cannot_bypass_admission_through_generic_entry_point() {
    let gate = CancelGate::new_clock_free();
    let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    gas.prandtl = 100.0;
    assert!(CylinderFlow::from_gas(0.001, gas, resolution(), &gate).is_err());
    assert!(CylinderFlow::from_fluid(0.001, gas, resolution(), &gate).is_err());
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE, f64::MAX] {
        assert!(CylinderFlow::from_fluid(radius, liquid(1000.0, 0.001), resolution(), &gate).is_err());
    }
}

#[test]
fn original_source_is_retained_and_evaluated_once_not_per_memory_branch() {
    struct CountedSource {
        identity: &'static str,
        calls: Cell<usize>,
        state: FluidTransportState,
    }
    impl FluidTransport for CountedSource {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.state)
        }
    }
    let source = CountedSource { identity: "explicit-test-transport", calls: Cell::new(0),
        state: liquid(1000.0, 0.001) };
    let flow = CylinderFlow::from_fluid(0.001, source, resolution(), &CancelGate::new_clock_free()).unwrap();
    assert_eq!(flow.fluid().identity, "explicit-test-transport");
    assert_eq!(flow.fluid().calls.get(), 1);
    for omega in [1.0, 10.0, 100.0] { flow.impedance_per_length(omega).unwrap(); }
    assert_eq!(flow.fluid().calls.get(), 1);
}
