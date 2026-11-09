//! G1/G3/G4: original conductive flux, nonlinear storage and coupled feedback.
use super::*;
use super::super::{
    self as thermal, FiniteChannelNetworkLink, GasTransportPaths,
    tests::{channel, endpoints, material, material_with_conductivity}, vacuum,
};
use super::super::super::super::tests::{controls, gas, near};
use crate::thermal::transport::volumetric::timeline::FixedVolumeThermalLaw;
use fs_thermal::enthalpy::{CellId, ThermalLink};

fn native<'l, 'p>(
    id: u64,
    ids: [u64; 2],
    laws: [&'l FixedVolumeThermalLaw<'p>; 2],
    area: f64,
) -> FixedVolumeConductionLink<'l, 'p> {
    FixedVolumeConductionLink::try_aligned(
        &ThermalLink {
            id, a: CellId(ids[0]), b: CellId(ids[1]), area_m2: area,
            // The prepared original source replaces this authored placeholder.
            conductance_w_k: 0.0,
        },
        laws, [0.01; 2], 0.0,
    ).unwrap()
}

fn layer_controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-10,
    }
}

#[test]
fn conduction_uses_original_nonlinear_volume_storage_and_material_entropy() {
    let source = material(1e6, 1000.0, [200.0, 1800.0]);
    let gases = BTreeMap::from([(0, gas(600.0, 100_000.0, 0.0))]);
    let mut previous: Option<[f64; 2]> = None;
    for (density, datum) in [(1000.0, 200.0), (8000.0, 1200.0)] {
        let law = FixedVolumeThermalLaw::try_new(&source, density, datum)
            .unwrap().with_sensible_entropy();
        let walls = BTreeMap::from([
            (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap()),
            (10, FixedVolumeChannelWall::try_new_isochoric(&law, 2e-6, 300.0).unwrap()),
        ]);
        let conductor = ChannelWallConduction::Cellwise(native(20, [9, 10], [&law; 2], 0.02));
        let result = thermal::step_candidate_with_conduction(
            &gases, &walls, GasTransportPaths::default(), &[], &[conductor], 0.01,
            &endpoints(&gases), controls(), || false,
        ).unwrap();
        let transfer = &result.step().conduction[0];
        near(transfer.conductance_w_k, 1.0, 2e-15);
        near(transfer.heat_a_to_b_j, 6.0, 2e-14);
        let mut entropy = 0.0;
        for (id, expected_heat) in [(9, -6.0), (10, 6.0)] {
            let old = &walls[&id];
            let next = &result.walls()[&id];
            let (a, b) = (old.temperature_k(), next.temperature_k());
            let energy = old.volume_m3() * (b - a) * (1e6 + 500.0 * (a + b));
            let ds = old.volume_m3() * (1e6 * ((b - a) / a).ln_1p() + 1000.0 * (b - a));
            near(energy, expected_heat, 2e-9);
            let receipt = result.step().walls.iter().find(|row| row.id == id).unwrap();
            near(receipt.energy_change_j, expected_heat, 2e-9);
            near(receipt.entropy_change_j_k, ds, 2e-12);
            assert_eq!(receipt.molecular_heat_into_wall_j, 0.0);
            assert_eq!(receipt.conduction_heat_into_wall_j, receipt.heat_into_wall_j);
            assert_eq!(next.mass_kg(), old.mass_kg());
            assert_eq!(next.law().identity(), old.law().identity());
            entropy += ds;
        }
        near(result.step().total_entropy_change_j_k, entropy, 2e-12);
        assert!(entropy > 0.0);
        assert!(result.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
        assert_eq!(result.parcels()[&0].mass_kg(), gases[&0].mass_kg());
        let temperatures = [result.walls()[&9].temperature_k(), result.walls()[&10].temperature_k()];
        if let Some(old) = previous {
            for i in 0..2 { near(temperatures[i], old[i], 2e-10); }
        }
        previous = Some(temperatures);
    }
}

#[test]
fn integrated_original_k_curves_match_the_independent_kirchhoff_flux() {
    // k(T)=0.01*T, equal source-free layers, no contact jump:
    // Q=A/(d_a+d_b) * integral_(Tb)^(Ta) k(T) dT.
    let source = material_with_conductivity(1e6, 0.0, [200.0, 1800.0], [2.0, 18.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1500.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap()),
        (10, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 300.0).unwrap()),
    ]);
    let gases = BTreeMap::from([(0, gas(500.0, 100_000.0, 0.0))]);
    for ids in [[9, 10], [10, 9]] {
        let path = ChannelWallConduction::Integrated {
            link: native(20, ids, [&law; 2], 0.02).into_integrated_layers().unwrap(),
            controls: layer_controls(),
        };
        let accepted = thermal::step_candidate_with_conduction(
            &gases, &walls, GasTransportPaths::default(), &[], &[path], 1e-4,
            &endpoints(&gases), controls(), || false,
        ).unwrap();
        let q = 0.005 * (walls[&ids[0]].temperature_k().powi(2)
            - walls[&ids[1]].temperature_k().powi(2));
        let step = &accepted.step().conduction[0];
        near(step.power_a_to_b_w, q, 1e-6);
        near(step.heat_a_to_b_j, q * 1e-4, 1e-10);
        let layers = step.layers.unwrap();
        near(layers.interface_temperature_k[0], (0.5_f64 * (900.0_f64.powi(2) + 300.0_f64.powi(2))).sqrt(), 1e-7);
        near(layers.interface_temperature_k[0], layers.interface_temperature_k[1], 1e-10);
        near(accepted.walls()[&9].temperature_k(), 899.64, 2e-8);
        near(accepted.walls()[&10].temperature_k(), 300.36, 2e-8);
    }
    let cellwise = native(20, [9, 10], [&law; 2], 0.02)
        .conductance_at_temperatures([900.0, 300.0]).unwrap();
    near(cellwise * 600.0, 2700.0, 1e-10);
}

#[test]
fn conductive_and_molecular_heat_pool_before_one_inversion_in_both_gas_paths() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 2800.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 500.0).unwrap()),
        (10, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 300.0).unwrap()),
    ]);
    let path = channel(2e-8);
    let conduction = [ChannelWallConduction::Cellwise(native(20, [9, 10], [&law; 2], 0.002))];
    let gases = BTreeMap::from([(0, gas(800.0, 100_000.0, 0.0)), (1, gas(800.0, 100_000.0, 0.0))]);
    let links = [FiniteChannelNetworkLink { id: 7, parcels: [0, 1], channel: &path, area_m2: 1e-6, wall: 9 }];
    let baseline = thermal::step_candidate(
        &gases, &walls, GasTransportPaths::default(), &links, 0.001,
        &endpoints(&gases), controls(), || false,
    ).unwrap();
    let accepted = thermal::step_candidate_with_conduction(
        &gases, &walls, GasTransportPaths::default(), &links, &conduction, 0.001,
        &endpoints(&gases), controls(), || false,
    ).unwrap();
    let wall = &accepted.step().walls[0];
    assert!(wall.molecular_heat_into_wall_j > 0.0);
    near(wall.molecular_heat_into_wall_j, baseline.step().walls[0].heat_into_wall_j, 1e-15);
    near(wall.conduction_heat_into_wall_j, -0.02, 1e-15);
    near(wall.heat_into_wall_j, wall.molecular_heat_into_wall_j - 0.02, 1e-15);
    let (t0, t1) = (500.0, accepted.walls()[&9].temperature_k());
    near(wall.heat_into_wall_j, 1e-6 * (t1 - t0) * (1e6 + 350.0 * (t0 + t1)), 1e-9);
    assert_eq!(accepted.step().channels, baseline.step().channels);
    for id in [0, 1] {
        assert_eq!(accepted.parcels()[&id].state().state().temperature, baseline.parcels()[&id].state().state().temperature);
    }
    let outlets = [vacuum::FiniteVacuumChannelLink { id: 7, parcel: 0, channel: &path, area_m2: 1e-6, wall: 9 }];
    let escaped = vacuum::step_candidate_with_conduction(
        &gases, &walls, GasTransportPaths::default(), &outlets, &conduction, 0.001,
        &endpoints(&gases), controls(), || false,
    ).unwrap();
    let plain = vacuum::step_candidate(
        &gases, &walls, GasTransportPaths::default(), &outlets, 0.001,
        &endpoints(&gases), controls(), || false,
    ).unwrap();
    assert_eq!(escaped.step().gas, plain.step().gas);
    assert!(escaped.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    near(escaped.step().walls.iter().map(|row| row.conduction_heat_into_wall_j).sum(), 0.0, 1e-15);
}

#[test]
fn source_binding_and_integrated_final_support_remain_authoritative() {
    let source_a = material(1e6, 0.0, [200.0, 800.0]);
    let source_b = material(1e6, 0.0, [200.0, 550.0]);
    let a = FixedVolumeThermalLaw::try_new(&source_a, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let b = FixedVolumeThermalLaw::try_new(&source_b, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let wrong = FixedVolumeThermalLaw::try_new(&source_a, 2000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let gases = BTreeMap::from([(0, gas(1000.0, 100_000.0, 0.0)), (1, gas(1000.0, 100_000.0, 0.0))]);
    let mut walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&a, 1e-6, 500.0).unwrap()),
        (10, FixedVolumeChannelWall::try_new_isochoric(&b, 1e-6, 300.0).unwrap()),
    ]);
    let bad = [ChannelWallConduction::Cellwise(native(20, [9, 10], [&wrong, &b], 0.02))];
    assert!(matches!(thermal::step_candidate_with_conduction(
        &gases, &walls, GasTransportPaths::default(), &[], &bad, 1e-4,
        &endpoints(&gases), controls(), || false,
    ), Err(FiniteChannelError::Invalid(_))));
    let path = channel(2e-8);
    let links = [FiniteChannelNetworkLink { id: 7, parcels: [0, 1], channel: &path, area_m2: 1e-6, wall: 9 }];
    let probe = thermal::step_candidate(
        &gases, &walls, GasTransportPaths::default(), &links, 1e-4,
        &endpoints(&gases), controls(), || false,
    ).unwrap();
    // Select a finite capacity for a 100 K molecular-only rise. The endpoint
    // remains inside a's caloric support but leaves the OTHER leg's k support.
    let volume = probe.step().walls[0].heat_into_wall_j / (1e6 * 100.0);
    walls.insert(9, FixedVolumeChannelWall::try_new_isochoric(&a, volume, 500.0).unwrap());
    let conductor = [ChannelWallConduction::Integrated {
        link: native(20, [9, 10], [&a, &b], 2e-10).into_integrated_layers().unwrap(),
        controls: layer_controls(),
    }];
    assert!(matches!(thermal::step_candidate_with_conduction(
        &gases, &walls, GasTransportPaths::default(), &links, &conductor, 1e-4,
        &endpoints(&gases), controls(), || false,
    ), Err(FiniteChannelError::Conduction { id: 20, error: DepositionError::MaterialDomain })));
    assert_eq!(walls[&9].temperature_k(), 500.0);
    assert_eq!(gases[&0].state().state().temperature, 1000.0);
}

fn accuracy(duration: f64, relative: f64) -> thermal::trajectory::EffusionAdvanceControls {
    thermal::trajectory::EffusionAdvanceControls {
        initial_step_s: duration, minimum_step_s: duration * 1e-6, maximum_step_s: duration,
        maximum_attempts: 4096, maximum_accepted_intervals: 1024,
        relative_tolerance: relative, temperature_tolerance_k: 1e-7,
        pressure_tolerance_pa: 1e-5, species_mass_tolerance_kg: 1e-16,
        stream_energy_tolerance_j: 1e-8,
    }
}

#[test]
fn adaptive_conduction_converges_to_analytic_relaxation_and_replays_native_steps() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap()),
        (10, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 300.0).unwrap()),
    ]);
    let gases = BTreeMap::from([(0, gas(600.0, 100_000.0, 0.0))]);
    let conduction = [ChannelWallConduction::Cellwise(native(20, [9, 10], [&law; 2], 0.002))];
    let duration = 0.5;
    let exact = 600.0 + 300.0 * (-0.2_f64 * duration).exp();
    let mut errors = Vec::new();
    for relative in [1e-2, 1e-4] {
        let result = thermal::trajectory::advance_candidate_with_conduction(
            &gases, &walls, GasTransportPaths::default(), &[], &conduction, duration,
            &endpoints(&gases), controls(), accuracy(duration, relative), || false,
        ).unwrap();
        errors.push((result.walls()[&9].temperature_k() - exact).abs());
        assert!(result.step().balance.total_entropy_change_j_k > 0.0);
        assert!(result.step().balance.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
        let (mut g, mut w) = (gases.clone(), walls.clone());
        for step in &result.step().accepted_steps {
            let replay = thermal::step_candidate_with_conduction(
                &g, &w, GasTransportPaths::default(), &[], &conduction, step.duration_s,
                &endpoints(&gases), controls(), || false,
            ).unwrap();
            assert_eq!(replay.step(), step);
            let (next_g, next_w, _) = replay.into_parts();
            (g, w) = (next_g, next_w);
        }
        assert_eq!(w[&9].temperature_k(), result.walls()[&9].temperature_k());
    }
    assert!(errors[1] < errors[0] / 3.0);
    assert!(errors[1] < 0.1);
    // State tolerances alone deliberately permit the coarse step. Conductive
    // heat accuracy must still drive refinement through the shared controller.
    let heat_only = thermal::trajectory::EffusionAdvanceControls {
        relative_tolerance: 0.0, temperature_tolerance_k: 1e6,
        pressure_tolerance_pa: 1e9, species_mass_tolerance_kg: 1.0,
        stream_energy_tolerance_j: 1e-4, ..accuracy(duration, 0.0)
    };
    let result = vacuum::trajectory::advance_candidate_with_conduction(
        &gases, &walls, GasTransportPaths::default(), &[], &conduction, duration,
        &endpoints(&gases), controls(), heat_only, || false,
    ).unwrap();
    assert!(result.step().accepted_intervals > 1);
    assert_eq!(result.step().balance.gas.escaped_mass_kg, 0.0);
    assert!(result.step().balance.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    assert!(result.step().balance.wall_entropy_change_j_k > 0.0);
}

#[test]
fn closed_paths_reordering_budgets_and_late_cancellation_preserve_original_state() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap()),
        (10, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 300.0).unwrap()),
    ]);
    let gases = BTreeMap::from([(0, gas(600.0, 100_000.0, 0.0))]);
    let conductors = [
        ChannelWallConduction::Cellwise(native(21, [9, 10], [&law; 2], 0.001)),
        ChannelWallConduction::Cellwise(native(20, [9, 10], [&law; 2], 0.001)),
    ];
    let advance = |links: &[ChannelWallConduction<'_, '_>], cancel: &mut dyn FnMut() -> bool| {
        thermal::step_candidate_with_conduction(
            &gases, &walls, GasTransportPaths::default(), &[], links, 0.01,
            &endpoints(&gases), controls(), cancel,
        )
    };
    let a = advance(&conductors, &mut || false).unwrap();
    let b = advance(&[conductors[1].clone(), conductors[0].clone()], &mut || false).unwrap();
    assert_eq!(a.step(), b.step());
    let mut polls = 0;
    advance(&conductors, &mut || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(matches!(advance(&conductors, &mut || { count += 1; count == polls }), Err(FiniteChannelError::Cancelled)));
    let closed = [ChannelWallConduction::Cellwise(native(20, [9, 10], [&law; 2], 0.0))];
    let accepted = advance(&closed, &mut || false).unwrap();
    assert_eq!(accepted.walls()[&9].specific_energy_j_kg(), walls[&9].specific_energy_j_kg());
    assert_eq!(accepted.step().conduction[0].heat_a_to_b_j, 0.0);
    assert!(matches!(advance(&[conductors[0].clone(), conductors[0].clone()], &mut || false), Err(FiniteChannelError::Invalid(_))));
    let over = vec![conductors[0].clone(); MAX_SNAPSHOT_TRANSFERS + 1];
    assert!(matches!(advance(&over, &mut || false), Err(FiniteChannelError::Invalid(_))));
    assert_eq!(walls[&9].temperature_k(), 900.0);
    assert_eq!(walls[&10].temperature_k(), 300.0);
}
