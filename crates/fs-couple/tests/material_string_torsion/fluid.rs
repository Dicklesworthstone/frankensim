//! G1/G3/G5: gas inertia and memory change the SAME physical source and receiver.
use super::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{
    StringBodyTranslationConnection, TransverseDirection, TransverseDirections,
};

fn resolution() -> CylinderFlowResolution {
    CylinderFlowResolution { minimum_wavenumber_radius: 1e-3,
        maximum_wavenumber_radius: 30.0, intervals: 12 }
}
fn small_source() -> MaterialStringSource {
    let mut description = template();
    description.n_modes = 1;
    MaterialStringSource::Homogeneous(with_uniform_circular_material_state(
        description, 0.0003, &material(7800.0, 80e9, 1000.0, 200_000.0),
    ).unwrap().with_kelvin_voigt_bending_loss().unwrap())
        .with_torsion("G", "eta_s", 1).unwrap()
        .with_transverse_polarizations(TransverseDirections::try_new([1.0, 0.0], [0.6, 0.8]).unwrap()).unwrap()
}
fn base(cx: &Cx<'_>, gate: &CancelGate, gas: GasState) -> ReciprocalStringBank {
    let mut config = options();
    config.network.max_modes = 64;
    config.network.max_connections = 64;
    config.network.max_setup_terms = 65_536;
    ReciprocalStringBank::from_material_plate(
        cx, vec![small_source()], gas, body(1), law(), config, gate,
    ).unwrap()
}
fn rotation(bank: &ReciprocalStringBank) -> StringBodyRotationConnection {
    let moving = &bank.source().strings()[0].connections[0].body;
    let fixed = ModalAttachment { component: moving.component, shapes: vec![0.0; moving.shapes.len()] };
    rotational_link(ModalRotationAttachment::from_displacement_footprints(moving, &fixed, 0.1).unwrap())
}

#[test]
fn added_mass_changes_both_plane_frequencies_and_force_maps_but_not_solid_or_twist() {
    with_cx(|cx, gate| {
        let a = base(cx, gate, gas()).with_cylinder_gas_loading(resolution(), gate).unwrap();
        let dense = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 4.0*101_325.0).unwrap();
        let b = base(cx, gate, dense).with_cylinder_gas_loading(resolution(), gate).unwrap();
        let loading_a = &a.source().cylinder_gas_loading().unwrap().strings()[0];
        let loading_b = &b.source().cylinder_gas_loading().unwrap().strings()[0];
        assert_eq!(loading_a.solid_card().linear_density_kg_m, loading_b.solid_card().linear_density_kg_m);
        assert_eq!(loading_a.solid_card().tension_n, loading_b.solid_card().tension_n);
        let solid = loading_a.solid_card().linear_density_kg_m;
        let fluid = core::f64::consts::PI * 0.0003_f64.powi(2) * gas().density;
        assert!((a.cards()[0].linear_density_kg_m / (solid + fluid) - 1.0).abs() < 2e-15);
        assert!((b.cards()[0].linear_density_kg_m / (solid + 4.0*fluid) - 1.0).abs() < 1e-14);
        let expected = ((solid + fluid)/(solid + 4.0*fluid)).sqrt();
        assert!((b.cards()[0].fundamental_hz()/a.cards()[0].fundamental_hz() - expected).abs() < 2e-14);
        for bank in [&a, &b] {
            let modes = bank.system().components()[0].modes();
            assert_eq!(modes[0], modes[1], "one physical plane lost the same atmosphere");
            assert_eq!(bank.system().components().len(), 3, "fluid displaced a physical body index");
            let twist = &modes[2];
            assert_eq!(twist, &a.system().components()[0].modes()[2]);
            let card = &bank.cards()[0];
            let omega = card.mode_omega_rad_s(1);
            let k = core::f64::consts::PI / card.length_m;
            let flow = bank.source().cylinder_gas_loading().unwrap().strings()[0].flow();
            let expected_loss = 0.5 * (card.viscous_bending_n_m2_s*k.powi(4)
                + flow.low_tail_resistance_per_length_n_s_m2()) / card.linear_density_kg_m / omega;
            assert_eq!(modes[0].damping_ratio, expected_loss, "sampled gas loss was counted twice");
            assert!(bank.system().acceleration_pressure_weights().unwrap()[5..].iter().all(|w| *w == 0.0));
        }
        let pa = &a.system().port().left.shapes;
        let pb = &b.system().port().left.shapes;
        assert!((pb[0]/pa[0] - expected).abs() < 2e-14);
        assert_eq!(pa[2], pb[2], "gas added mass changed solid polar inertia");
    });
}

#[test]
fn viscous_memory_receives_energy_and_changes_actual_body_pressure_at_equal_density() {
    with_cx(|cx, gate| {
        let mut spec = GasSpec::dry_air_ussa1976();
        spec.sutherland_beta *= 4.0; // synthetic independent transport-property change
        let thick = GasState::try_new(&spec, 293.15, 101_325.0).unwrap();
        assert_eq!(thick.density, gas().density);
        let mut a = base(cx, gate, gas()).with_cylinder_gas_loading(resolution(), gate).unwrap();
        let mut b = base(cx, gate, thick).with_cylinder_gas_loading(resolution(), gate).unwrap();
        assert_eq!(a.cards()[0].fundamental_hz(), b.cards()[0].fundamental_hz());
        assert_eq!(a.system().components()[1].modes(), b.system().components()[1].modes());
        assert_eq!(a.system().acceleration_pressure_weights(), b.system().acceleration_pressure_weights());
        let zero = vec![0.0; a.system().mode_count()];
        let mut pressure_difference = 0.0_f64;
        let mut memory_energy = 0.0_f64;
        let mut work = 0.0;
        let mut loss = 0.0;
        for i in 0..192 {
            let drive = if i < 48 { drive() } else { FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 } };
            let fa = *a.step(&zero, drive).unwrap();
            let fb = *b.step(&zero, drive).unwrap();
            pressure_difference = pressure_difference.max((fa.observer_pressure_pa - fb.observer_pressure_pa).abs());
            work += fa.drive_work_j + fa.external_work_j;
            loss += fa.network_dissipation_j + fa.friction_dissipation_j;
            memory_energy = memory_energy.max(a.system().components()[2].states().iter()
                .map(|s| 0.5*s.velocity_m_sqrt_kg_per_s.powi(2)).sum::<f64>());
            assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
        }
        assert!(memory_energy > 0.0, "gas remained a diagnostic instead of a reciprocal mechanical state");
        assert!(pressure_difference > 1e-15, "equal-density gas viscosity never reached actual radiated pressure");
        assert!(loss > 0.0);
        assert!((a.system().total_energy_j().unwrap() + loss - work).abs() < 2e-11);
        assert!(a.system().components()[0].states()[1].velocity_m_sqrt_kg_per_s != 0.0,
            "body-mediated second-plane motion disappeared");
    });
}

#[test]
fn directed_and_rotational_attachment_order_preserves_complete_gas_memory_and_accepted_history() {
    with_cx(|cx, gate| {
        let plain = base(cx, gate, gas());
        let link = rotation(&plain);
        let translation = StringBodyTranslationConnection {
            string: 0, station_fraction: 0.8, direction: TransverseDirection::Y,
            body: plain.source().strings()[0].connections[0].body.clone(),
            stiffness_n_m: 500.0, damping_n_s_m: 0.02, rest_extension_m: 0.0,
        };
        let mut a = plain.with_translational_body_connections(vec![translation.clone()], gate).unwrap()
            .with_rotational_body_connections(vec![link.clone()], gate).unwrap()
            .with_cylinder_gas_loading(resolution(), gate).unwrap();
        let mut b = base(cx, gate, gas()).with_cylinder_gas_loading(resolution(), gate).unwrap()
            .with_rotational_body_connections(vec![link], gate).unwrap()
            .with_translational_body_connections(vec![translation], gate).unwrap();
        assert!(b.source().cylinder_gas_loading().is_some());
        assert_eq!(a.system().mode_count(), b.system().mode_count());
        assert_eq!(a.system().acceleration_pressure_weights(), b.system().acceleration_pressure_weights());
        assert_eq!(a.source().body_connections(), b.source().body_connections());
        assert_eq!(a.source().rotational_body_connections(), b.source().rotational_body_connections());
        assert_eq!(a.source().translational_body_connections(), b.source().translational_body_connections());
        let zero = vec![0.0; a.system().mode_count()];
        for _ in 0..64 {
            assert_eq!(a.step(&zero, drive()).unwrap(), b.step(&zero, drive()).unwrap());
        }
        for (a, b) in a.system().components().iter().zip(b.system().components()) {
            assert_eq!(a.states(), b.states());
        }
    });
}

#[test]
fn static_preload_retains_zero_fluid_position_gauges_and_releases_into_memory() {
    with_cx(|cx, gate| {
        let bank = base(cx, gate, gas()).with_cylinder_gas_loading(resolution(), gate).unwrap();
        let (_, mut system) = bank.into_parts();
        let mut force = vec![0.0; system.mode_count()];
        for (f, shape) in force.iter_mut().zip(&system.port().left.shapes) { *f = 1e-4*shape; }
        let initial = system.initialize_released_equilibrium(&force, gate).unwrap();
        assert!(initial > 0.0);
        assert!(system.components()[2].states().iter().all(|s|
            s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0));
        let zero = vec![0.0; force.len()];
        let mut dissipation = 0.0;
        for _ in 0..64 {
            let frame = system.step(&zero, FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }).unwrap();
            dissipation += frame.network_dissipation_j;
            assert_eq!(frame.external_work_j, 0.0);
            assert_eq!(frame.drive_work_j, 0.0);
        }
        assert!(system.components()[2].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        assert!((system.total_energy_j().unwrap() + dissipation - initial).abs() < 2e-11);
    });
}

#[test]
fn complete_memory_budget_repeated_binding_and_live_changes_refuse() {
    with_cx(|cx, gate| {
        let over = CylinderFlowResolution { intervals: 64, ..resolution() };
        assert!(base(cx, gate, gas()).with_cylinder_gas_loading(over, gate).is_err());
        let loaded = base(cx, gate, gas()).with_cylinder_gas_loading(resolution(), gate).unwrap();
        assert!(loaded.with_cylinder_gas_loading(resolution(), gate).is_err());
        let mut moving = base(cx, gate, gas());
        moving.step(&vec![0.0; moving.system().mode_count()], drive()).unwrap();
        assert!(moving.with_cylinder_gas_loading(resolution(), gate).is_err());
    });
}
