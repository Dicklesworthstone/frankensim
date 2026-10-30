//! The same material drives diffusion, buoyancy and live coefficient refresh.
//! Authored liquid coefficients below are numerical fixtures, not measurements.

use fs_lbm::material::{MaterialBuoyancyBoundary, plan_scaling_from_material};
use fs_material::MaterialError;
use fs_material::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_qty::semantic::{Composition, CompositionBasis};

fn boundary() -> MaterialBuoyancyBoundary {
    MaterialBuoyancyBoundary {
        bottom_temperature_k: 302.0,
        top_temperature_k: 298.0,
        gravity_m_s2: 9.81,
        max_relative_density_change: 0.02,
    }
}

fn liquid(beta: Option<f64>, conductivity: f64, viscosity: f64) -> ThermalFluidTransportState {
    ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(1000.0, 4000.0, conductivity, viscosity).unwrap(),
        300.0,
        beta,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-20),
        "{actual} != {expected}"
    );
}

#[test]
fn gas_pressure_and_composition_reach_both_solvers_without_separate_buoyancy_data() {
    let air = GasSpec::dry_air_ussa1976();
    let light = GasSpec {
        molar_mass: 0.004,
        gamma: 5.0 / 3.0,
        sutherland_beta: 1e-6,
        sutherland_s: 80.0,
        conductivity: ConductivityModel::Eucken,
    };
    let fractions = Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap();
    let gases = [
        GasState::try_new(&air, 300.0, 100_000.0).unwrap(),
        GasState::try_new(&air, 300.0, 200_000.0).unwrap(),
        GasState::try_new_mixture(&[air, light], &fractions, 300.0, 100_000.0, None).unwrap(),
    ];
    let mut outcomes = Vec::new();
    for gas in &gases {
        let scaling = plan_scaling_from_material(gas, 0.001, 0.001, 0.01).unwrap();
        let mut two = scaling
            .thermal_slab_from_material([4, 8], boundary(), 40)
            .unwrap();
        let mut three = scaling
            .thermal_slab3_from_material([4, 2, 8], boundary(), 80)
            .unwrap();
        // Independent SI Rayleigh number at the actual 8 mm fluid height.
        let transport = gas.transport_state().unwrap();
        let expected_ra = 9.81 * (302.0 - 298.0) * 0.008_f64.powi(3)
            / (300.0 * transport.kinematic_viscosity_m2_s() * transport.thermal_diffusivity_m2_s());
        close(two.lattice.rayleigh(), expected_ra);
        close(two.lattice.gbeta, three.lattice.gbeta);
        close(
            (two.lattice.tau_g - 0.5) / 3.0,
            (three.lattice.tau_g - 0.5) / 4.0,
        );
        let i2 = two.lattice.grid.idx(0, 4);
        let i3 = three.lattice.grid.idx(0, 0, 4);
        two.lattice.gpop[i2][0] += 0.05;
        three.lattice.gpop[i3][0] += 0.05;
        for _ in 0..12 {
            two.lattice.step();
            three.lattice.step();
        }
        let t2 = two.lattice.temperature(0, 4);
        let t3 = three.lattice.temperature(0, 0, 4);
        assert!(t2.is_finite() && t3.is_finite());
        outcomes.push((t2, t3, expected_ra));
    }
    close(outcomes[1].2 / outcomes[0].2, 4.0);
    for outcome in &outcomes[1..] {
        assert!((outcome.0 - outcomes[0].0).abs() > 1e-6);
        assert!((outcome.1 - outcomes[0].1).abs() > 1e-6);
    }
}

#[test]
fn signed_and_zero_material_expansion_drive_the_actual_force() {
    for beta in [-0.001, 0.0, 0.001] {
        let source = liquid(Some(beta), 0.6, 0.001);
        let scaling = plan_scaling_from_material(&source, 0.0001, 0.001, 0.001).unwrap();
        let mut two = scaling
            .thermal_slab_from_material([4, 8], boundary(), 40)
            .unwrap();
        let mut three = scaling
            .thermal_slab3_from_material([4, 2, 8], boundary(), 80)
            .unwrap();
        let theta2 = two.lattice.temperature(0, 2);
        let theta3 = three.lattice.temperature(0, 0, 2);
        two.lattice.step();
        three.lattice.step();
        let expected = 9.81 * beta * 4.0 * 0.001_f64.powi(2) / 0.0001;
        close(
            two.lattice.grid.fext[two.lattice.grid.idx(0, 2)][1],
            expected * (theta2 - 0.5),
        );
        close(
            three.lattice.grid.fext[three.lattice.grid.idx(0, 0, 2)][2],
            expected * (theta3 - 0.5),
        );
    }
}

#[test]
fn off_center_material_temperature_sets_the_actual_buoyancy_reference() {
    let source = liquid(Some(0.001), 0.6, 0.001);
    let scaling = plan_scaling_from_material(&source, 0.0001, 0.001, 0.001).unwrap();
    let walls = MaterialBuoyancyBoundary {
        bottom_temperature_k: 306.0,
        top_temperature_k: 302.0,
        ..boundary()
    };
    let mut two = scaling
        .thermal_slab_from_material([4, 8], walls, 40)
        .unwrap();
    let mut three = scaling
        .thermal_slab3_from_material([4, 2, 8], walls, 80)
        .unwrap();
    assert_eq!(two.lattice.t_ref, -0.5);
    assert_eq!(three.lattice.t_ref, -0.5);
    let t2 = two.temperature_si(two.lattice.temperature(0, 2)).unwrap();
    let t3 = three
        .temperature_si(three.lattice.temperature(0, 0, 2))
        .unwrap();
    two.lattice.step();
    three.lattice.step();
    let force =
        |temperature: f64| 9.81 * 0.001 * (temperature - 300.0) * 0.001_f64.powi(2) / 0.0001;
    close(
        two.lattice.grid.fext[two.lattice.grid.idx(0, 2)][1],
        force(t2),
    );
    close(
        three.lattice.grid.fext[three.lattice.grid.idx(0, 0, 2)][2],
        force(t3),
    );
}

#[test]
fn live_material_refresh_updates_buoyancy_and_both_diffusivities_without_resetting_state() {
    let before = liquid(Some(0.001), 0.6, 0.001);
    let after = ThermalFluidTransportState::try_new(
        liquid(Some(-0.002), 0.9, 0.002).transport(),
        301.0,
        Some(-0.002),
    )
    .unwrap();
    let previous = plan_scaling_from_material(&before, 0.0001, 0.001, 0.001).unwrap();
    let next = plan_scaling_from_material(&after, 0.0001, 0.001, 0.001).unwrap();
    let mut two = previous
        .thermal_slab_from_material([4, 8], boundary(), 40)
        .unwrap();
    let mut three = previous
        .thermal_slab3_from_material([4, 2, 8], boundary(), 80)
        .unwrap();
    two.lattice.step();
    three.lattice.step();
    let old2 = (
        two.lattice.grid.f.clone(),
        two.lattice.gpop.clone(),
        two.lattice.gbeta,
    );
    let old3 = (
        three.lattice.grid.f.clone(),
        three.lattice.gpop.clone(),
        three.lattice.gbeta,
    );
    next.rebind_thermal_slab_from_material(&mut two, &previous, 9.81, 0.02)
        .unwrap();
    next.rebind_thermal_slab3_from_material(&mut three, &previous, 9.81, 0.02)
        .unwrap();
    assert_eq!((&two.lattice.grid.f, &two.lattice.gpop), (&old2.0, &old2.1));
    assert_eq!(
        (&three.lattice.grid.f, &three.lattice.gpop),
        (&old3.0, &old3.1)
    );
    close(two.lattice.gbeta, -2.0 * old2.2);
    close(three.lattice.gbeta, -2.0 * old3.2);
    assert_eq!(two.lattice.t_ref, 0.75);
    assert_eq!(three.lattice.t_ref, 0.75);
    close(two.lattice.grid.tau[8], next.flow().tau);
    close(three.lattice.grid.tau[8], next.flow().tau);
    close(two.lattice.tau_g, next.thermal_relaxation_time());
    close(
        three.lattice.tau_g,
        next.thermal_relaxation_time3().unwrap(),
    );
}

#[test]
fn current_temperature_field_is_checked_before_any_rebind_mutation() {
    let before = liquid(Some(0.001), 0.6, 0.001);
    let after = ThermalFluidTransportState::try_new(
        liquid(Some(0.002), 0.9, 0.002).transport(),
        301.0,
        Some(0.002),
    )
    .unwrap();
    let previous = plan_scaling_from_material(&before, 0.0001, 0.001, 0.001).unwrap();
    let next = plan_scaling_from_material(&after, 0.0001, 0.001, 0.001).unwrap();
    let mut two = previous
        .thermal_slab_from_material([4, 8], boundary(), 40)
        .unwrap();
    let mut three = previous
        .thermal_slab3_from_material([4, 2, 8], boundary(), 80)
        .unwrap();
    two.lattice.gpop[12][0] += 10.0;
    three.lattice.gpop[20][0] += 10.0;
    let old2 = (
        two.lattice.grid.f.clone(),
        two.lattice.gpop.clone(),
        two.lattice.grid.tau.clone(),
        two.lattice.tau_g,
        two.lattice.gbeta,
        two.lattice.t_ref,
    );
    let old3 = (
        three.lattice.grid.f.clone(),
        three.lattice.gpop.clone(),
        three.lattice.grid.tau.clone(),
        three.lattice.tau_g,
        three.lattice.gbeta,
        three.lattice.t_ref,
    );
    assert!(
        next.rebind_thermal_slab_from_material(&mut two, &previous, 9.81, 0.02)
            .is_err()
    );
    assert!(
        next.rebind_thermal_slab3_from_material(&mut three, &previous, 9.81, 0.02)
            .is_err()
    );
    assert_eq!(
        (
            two.lattice.grid.f,
            two.lattice.gpop,
            two.lattice.grid.tau,
            two.lattice.tau_g,
            two.lattice.gbeta,
            two.lattice.t_ref
        ),
        old2
    );
    assert_eq!(
        (
            three.lattice.grid.f,
            three.lattice.gpop,
            three.lattice.grid.tau,
            three.lattice.tau_g,
            three.lattice.gbeta,
            three.lattice.t_ref
        ),
        old3
    );
}

#[test]
fn unknown_expansion_stale_transport_and_out_of_regime_states_refuse() {
    let unknown = liquid(None, 0.6, 0.001);
    let known = liquid(Some(0.001), 0.6, 0.001);
    let missing = plan_scaling_from_material(&unknown, 0.0001, 0.001, 0.001).unwrap();
    assert!(
        missing
            .thermal_slab_from_material([4, 8], boundary(), 40)
            .is_err()
    );
    assert!(
        missing
            .thermal_slab3_from_material([4, 2, 8], boundary(), 80)
            .is_err()
    );

    struct Incoherent(ThermalFluidTransportState);
    impl FluidTransport for Incoherent {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            FluidTransportState::try_new(1000.0, 4000.0, 0.7, 0.001)
        }
    }
    impl ThermalFluidTransport for Incoherent {
        fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
            Ok(self.0)
        }
    }
    let source = Incoherent(known);
    let stale = plan_scaling_from_material(&source, 0.0001, 0.001, 0.001).unwrap();
    assert!(
        stale
            .thermal_slab_from_material([4, 8], boundary(), 40)
            .is_err()
    );

    let scaling = plan_scaling_from_material(&known, 0.0001, 0.001, 0.001).unwrap();
    for limit in [0.0, -0.1, 1.0, f64::NAN, f64::INFINITY, 0.001] {
        assert!(
            scaling
                .thermal_slab_from_material(
                    [4, 8],
                    MaterialBuoyancyBoundary {
                        max_relative_density_change: limit,
                        ..boundary()
                    },
                    40
                )
                .is_err()
        );
    }
    for temperature in [0.0, f64::NAN, f64::INFINITY, 400.0] {
        assert!(
            scaling
                .thermal_slab3_from_material(
                    [4, 2, 8],
                    MaterialBuoyancyBoundary {
                        bottom_temperature_k: temperature,
                        ..boundary()
                    },
                    80
                )
                .is_err()
        );
    }
    assert!(
        scaling
            .thermal_slab_from_material([4, 8], boundary(), 39)
            .is_err()
    );
    assert!(
        scaling
            .thermal_slab3_from_material([4, 2, 8], boundary(), 79)
            .is_err()
    );
}
