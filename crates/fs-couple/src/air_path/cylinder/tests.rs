//! G1/G3: complex Stokes impedance and three independently live refinement axes.
use super::*;

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn resolution(min: f64, max: f64, intervals: usize) -> CylinderFlowResolution {
    CylinderFlowResolution { minimum_wavenumber_radius: min,
        maximum_wavenumber_radius: max, intervals }
}
fn flow(resolution: CylinderFlowResolution) -> CylinderFlow {
    let gas = gas(101_325.0);
    // Set rho*a²/mu = 1 s, so omega itself is the dimensionless frequency.
    let radius = (gas.dynamic_viscosity / gas.density).sqrt();
    CylinderFlow::from_gas(radius, gas, resolution, &CancelGate::new_clock_free()).unwrap()
}
fn error(flow: &CylinderFlow, omega: f64, expected: C64) -> f64 {
    let z = flow.impedance_per_length(omega).unwrap();
    let mu = flow.gas().dynamic_viscosity;
    (z.re / mu - expected.re).hypot(z.im / mu - expected.im)
}

// Independent defining complex I0/I1/K0 power series, differentiated for K1,
// evaluated at z=(1+i)*sqrt(omega/2); NOT the real-Bessel spectral quadrature.
// Formula: pi*i*w + 4*pi*z*K1(z)/K0(z), NASA/CR-2007-215050 equation (7).
// Constants need only 1e-8 relative accuracy here; cutoffs dominate the 9e-4 gate.
const ORACLE: [(f64, C64); 5] = [
    (0.01, C64 { re: 4.710527336547042, im: 1.4674730277805959 }),
    (0.1, C64 { re: 7.489342045271669, im: 3.6880439694970404 }),
    (1.0, C64 { re: 14.346304543275211, im: 12.465886293099675 }),
    (10.0, C64 { re: 34.05643752523019, im: 59.754434003390436 }),
    (100.0, C64 { re: 95.03108417164636, im: 403.11396197909926 }),
];

#[test]
fn complex_cylinder_force_matches_independent_modified_bessel_solution() {
    let flow = flow(resolution(1e-5, 3000.0, 512));
    for &(omega, expected) in &ORACLE {
        let relative = error(&flow, omega, expected) / expected.re.hypot(expected.im);
        assert!(relative < 9e-4, "omega={omega}, relative complex-force error={relative:e}");
        let z = flow.impedance_per_length(omega).unwrap();
        assert!(z.re > 0.0 && z.im > flow.added_mass_per_length_kg_m() * omega);
    }
    for branch in flow.branches() {
        assert!(branch.mass_per_length_kg_m > 0.0 && branch.resistance_per_length_n_s_m2 > 0.0);
        assert!((branch.resistance_per_length_n_s_m2
            / branch.mass_per_length_kg_m / branch.rate_per_s - 1.0).abs() < 3e-15);
    }
}

#[test]
fn grid_count_and_both_spectral_cutoffs_are_live_independent_axes() {
    let mut previous = f64::INFINITY;
    for n in [8, 16, 32] {
        let f = flow(resolution(1e-4, 3000.0, n));
        let worst = ORACLE.iter().map(|&(w, expected)|
            error(&f, w, expected) / expected.re.hypot(expected.im)).fold(0.0_f64, f64::max);
        assert!(worst < previous / 3.0, "quadrature axis n={n}: {worst:e} vs {previous:e}");
        previous = worst;
    }
    assert!(previous < 8e-4);

    let low_oracle = C64 { re: 3.3581195206967602, im: 0.7336873702912944 };
    previous = f64::INFINITY;
    for minimum in [0.5, 0.1, 0.02, 0.004] {
        let f = flow(resolution(minimum, 3000.0, 512));
        let e = error(&f, 0.001, low_oracle);
        assert!(e < previous / 3.0, "slow-tail axis minimum={minimum}: {e:e}");
        previous = e;
    }
    assert!(previous < 0.003);

    previous = f64::INFINITY;
    // After positive fast-tail completion the leading omission is cubic in
    // frequency. Resolve that error before the finite log grid becomes dominant.
    for maximum in [8.0, 16.0, 30.0] {
        let f = flow(resolution(1e-5, maximum, 512));
        let e = error(&f, 100.0, ORACLE[4].1);
        assert!(e < previous / 8.0, "fast-tail axis maximum={maximum}: {e:e}");
        previous = e;
    }
    assert!(previous < 0.04);
}

#[test]
fn fast_tail_completion_recovers_complex_force_without_more_memory_states() {
    // The former omitted tail gives errors of about 8 and 0.267 respectively
    // at omega=100. These are independent complex K-Bessel force references,
    // not a comparison to a finer copy of the real-Bessel spectral algorithm.
    for (maximum, error_limit) in [(100.0, 5e-4), (3000.0, 2.5e-5)] {
        let f = flow(resolution(1e-5, maximum, 512));
        assert_eq!(f.branches().len(), 512);
        for &(omega, expected) in &ORACLE {
            let e = error(&f, omega, expected);
            assert!(e < error_limit, "maximum={maximum}, omega={omega}, error={e:e}");
            let loss = f.impedance_per_length(omega).unwrap().re / f.gas().dynamic_viscosity;
            assert!((loss - expected.re).abs() < 5e-5,
                "maximum={maximum}, omega={omega}, resistive error={:e}", loss - expected.re);
        }
    }
}

#[test]
fn one_completed_tail_state_retains_native_work_and_positive_heat() {
    let gate = CancelGate::new_clock_free();
    let source = FluidTransportState::try_new(2.0, 4000.0, 0.5, 0.1).unwrap();
    let f = CylinderFlow::from_fluid(0.2, source, resolution(0.1, 2.0, 1), &gate).unwrap();
    let body = f.transverse_body(0.7, 0.6, 4.0, 3, &gate).unwrap();
    assert_eq!(body.system().state_dim(), 3);
    let mut state = body.state_from_motion(0.02, 0.3, &[-0.1]).unwrap();
    let initial_energy = body.system().hamiltonian(&state);
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..20 {
        let (force, dt) = (0.03, 0.02);
        let step = fs_phs::step(body.system(), &state, &[force], dt).unwrap();
        let velocity = 0.5 * (state[1] + step.x[1]) / body.moving_mass_kg();
        let fluid_velocity = 0.5 * (state[2] + step.x[2]) / body.memory_masses_kg()[0];
        let relative_velocity = velocity - fluid_velocity;
        heat += dt * body.length_m() * (
            f.low_tail_resistance_per_length_n_s_m2() * velocity * velocity
            + f.branches()[0].resistance_per_length_n_s_m2 * relative_velocity * relative_velocity);
        work += force * (step.x[0] - state[0]);
        state = step.x;
    }
    assert!(heat > 0.0);
    let energy_change = body.system().hamiltonian(&state) - initial_energy;
    assert!((energy_change + heat - work).abs() < 2e-12 * initial_energy);
}

#[test]
fn one_interval_fast_tail_has_independent_physical_inertia_and_resistance() {
    // Independent adaptive integration of 1/[k^3(J0^2+Y0^2)] and
    // 1/[k^5(J0^2+Y0^2)] on [0.001,infinity), with SciPy's real Bessel
    // implementation and QUADPACK, gives Mweight=27954.841894480702 and
    // Rweight=Mweight^2/Nweight=0.061008776802984721. This includes both
    // the slow decades and the infinite fast end of the original spectrum.
    let gate = CancelGate::new_clock_free();
    for (density, viscosity, radius) in [(1.2, 1.8e-5, 0.0003), (1000.0, 0.001, 0.001)] {
        let fluid = FluidTransportState::try_new(density, 1000.0, 0.025, viscosity).unwrap();
        let f = CylinderFlow::from_fluid(radius, fluid, resolution(0.001, 0.1, 1), &gate).unwrap();
        assert_eq!(f.branches().len(), 1);
        let branch = f.branches()[0];
        let mass = (16.0 / core::f64::consts::PI) * density * radius * radius
            * 27954.841894480702;
        let resistance = (16.0 / core::f64::consts::PI) * viscosity * 0.061008776802984721;
        assert!((branch.mass_per_length_kg_m / mass - 1.0).abs() < 2e-9);
        assert!((branch.resistance_per_length_n_s_m2 / resistance - 1.0).abs() < 2e-9);
        assert!(branch.rate_per_s > 0.0);
        assert_eq!(branch.rate_per_s, branch.resistance_per_length_n_s_m2 / branch.mass_per_length_kg_m);
    }
}

#[test]
fn pressure_changes_inertia_and_all_relaxation_rates_without_material_name_branches() {
    let gate = CancelGate::new_clock_free();
    let ambient = gas(101_325.0);
    let pressurized = gas(4.0 * 101_325.0);
    let r = resolution(1e-3, 100.0, 12);
    let a = CylinderFlow::from_gas(0.0003, ambient, r, &gate).unwrap();
    let b = CylinderFlow::from_gas(0.0003, pressurized, r, &gate).unwrap();
    assert_eq!(ambient.dynamic_viscosity, pressurized.dynamic_viscosity);
    assert!((b.added_mass_per_length_kg_m() / a.added_mass_per_length_kg_m() - 4.0).abs() < 1e-14);
    for (x, y) in a.branches().iter().zip(b.branches()) {
        assert_eq!(x.resistance_per_length_n_s_m2, y.resistance_per_length_n_s_m2);
        assert!((y.mass_per_length_kg_m / x.mass_per_length_kg_m - 4.0).abs() < 1e-14);
        assert!((y.rate_per_s / x.rate_per_s - 0.25).abs() < 1e-14);
    }
    for omega in [10.0, 100.0, 1000.0] {
        let za = a.impedance_per_length(omega).unwrap();
        let zb = b.impedance_per_length(omega / 4.0).unwrap();
        assert!((za.re-zb.re).hypot(za.im-zb.im) < 1e-14 * za.re.hypot(za.im));
    }
}

#[test]
fn invalid_or_unrepresentable_resolution_is_refused_without_silent_caps() {
    for r in [resolution(0.0, 1.0, 10), resolution(f64::NAN, 1.0, 10),
        resolution(0.6, 1.0, 10), resolution(0.1, 0.1, 10),
        resolution(0.1, 3000.1, 10), resolution(0.1, 1.0, 0),
        resolution(0.1, 1.0, 4097)] {
        assert!(r.validate().is_err());
    }
    assert!(resolution(0.5, 3000.0, 4096).validate().is_ok());
    let f = flow(resolution(1e-3, 100.0, 12));
    assert!(f.impedance_per_length(0.0).is_err());
    assert!(f.impedance_per_length(f64::INFINITY).is_err());
    assert!(CylinderFlow::from_gas(0.0, gas(101_325.0), f.resolution(),
        &CancelGate::new_clock_free()).is_err());
}
