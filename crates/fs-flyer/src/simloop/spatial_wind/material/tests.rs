use super::*;
use crate::replayenv::{
    ReplayEnvelope, record_replay_with_material_distributed_spatial_wind,
    record_replay_with_material_spatial_wind,
};
use crate::simloop::{Phase, PilotMode, dec17_scenario};
use fs_atmo::FlatSiteLogLaw;
use fs_material::gas::GasSpec;

fn wind() -> SpatialWindSpec {
    SpatialWindSpec {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: 2.4,
            reference_speed_mps: 11.0,
        },
        sigma_mps: 0.3,
        length_m: 20.0,
        advection_mps: 11.0,
        correlation_times_s: vec![1.8; 8],
    }
}

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, pressure).unwrap()
}

fn initialize(distributed: bool, gas: GasState) -> SimLoop {
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    if distributed {
        SimLoop::init_with_material_distributed_spatial_wind(scenario, wind(), gas).unwrap()
    } else {
        SimLoop::init_with_material_spatial_wind(scenario, wind(), gas).unwrap()
    }
}

#[test]
fn material_density_reaches_original_trim_loads_and_actual_trajectory() {
    for distributed in [false, true] {
        let medium = gas(100_000.0);
        let mut bound = initialize(distributed, medium);
        let mut canonical = dec17_scenario(1903, PilotMode::FixedControls);
        assert_ne!(canonical.rho_kg_m3, medium.density);
        canonical.rho_kg_m3 = medium.density;
        let mut native = if distributed {
            SimLoop::init_with_distributed_spatial_wind(canonical, wind()).unwrap()
        } else {
            SimLoop::init_with_spatial_wind(canonical, wind()).unwrap()
        };
        assert_eq!(bound.scenario().rho_kg_m3, medium.density);
        assert_eq!(bound.tick0().trim, native.tick0().trim);
        assert_eq!(bound.material_gas(), Some(&medium));
        let before = bound.save_checkpoint().unwrap();
        let sample = bound.sample_material_air([4.0, -1.0, 3.0]).unwrap();
        assert_eq!(
            sample.wind,
            bound.sample_spatial_wind([4.0, -1.0, 3.0]).unwrap()
        );
        assert_eq!(sample.gas, medium);
        let observer = [6.0, 0.0, 0.0];
        let speed2: f64 = sample
            .wind
            .u
            .iter()
            .zip(observer)
            .map(|(u, v)| (u - v).powi(2))
            .sum();
        assert_eq!(
            sample.dynamic_pressure_pa(observer).unwrap(),
            0.5 * medium.density * speed2
        );
        let expected = medium.density * speed2.sqrt() * 1.5 / medium.dynamic_viscosity;
        assert!((sample.reynolds(1.5, observer).unwrap() - expected).abs() < 1e-9 * expected);
        assert_eq!(
            before,
            bound.save_checkpoint().unwrap(),
            "sampling cannot alter physics"
        );
        // Native equations are the independent implementation seam: every
        // motion and rotor sample must match a run initialized at the same rho.
        for _ in 0..8 {
            let a = bound.step(None).unwrap();
            let b = native.step(None).unwrap();
            assert_eq!(a, b);
            assert!(!matches!(a.phase, Phase::Ended(_)));
        }
    }
}

#[test]
fn pressure_changes_force_propeller_trim_and_flight_without_hidden_reference_density() {
    let mut low = initialize(true, gas(90_000.0));
    let mut high = initialize(true, gas(105_000.0));
    assert_ne!(low.tick0().trim, high.tick0().trim);
    let load = |sim: &SimLoop| {
        sim.design
            .force_buildup(13.0, 0.06, 0.1, 45.0, 0.0, sim.scenario().rho_kg_m3)
            .unwrap()
    };
    let a = load(&low);
    let b = load(&high);
    let ratio = high.scenario().rho_kg_m3 / low.scenario().rho_kg_m3;
    // At fixed kinematics the same inviscid circulation/BEMT solution scales
    // physical lift and both rotor thrusts with the supplied material density.
    for (left, right) in [
        (a.lift_n, b.lift_n),
        (a.thrust_n[0], b.thrust_n[0]),
        (a.thrust_n[1], b.thrust_n[1]),
    ] {
        assert!(
            (right / left - ratio).abs() < 1e-8,
            "density scaling: {left} {right}"
        );
    }
    assert_ne!(a.torque_imbalance_nm, b.torque_imbalance_nm);
    let mut last = None;
    for _ in 0..8 {
        let a = low.step(None).unwrap();
        let b = high.step(None).unwrap();
        assert!(!matches!(a.phase, Phase::Ended(_)));
        assert!(!matches!(b.phase, Phase::Ended(_)));
        last = Some((a, b));
    }
    let (a, b) = last.unwrap();
    assert_ne!(a.x_m, b.x_m);
    assert_ne!(a.omega_prop_rad_s, b.omega_prop_rad_s);
}

#[test]
fn complete_gas_and_sampling_mode_survive_checkpoint_and_replay() {
    for distributed in [false, true] {
        let medium = GasState::try_new_moist_air(293.15, 101_325.0, 0.6).unwrap();
        let mut original = initialize(distributed, medium);
        for _ in 0..4 {
            original.step(None).unwrap();
        }
        let checkpoint = original.save_checkpoint().unwrap();
        assert_eq!(u32::from_le_bytes(checkpoint[..4].try_into().unwrap()), 5);
        let mut resumed =
            SimLoop::restore_checkpoint(original.scenario().clone(), &checkpoint).unwrap();
        assert_eq!(resumed.material_gas(), Some(&medium));
        assert_eq!(
            resumed.spatial_wind.as_ref().unwrap().distributed,
            distributed
        );
        assert_eq!(resumed.save_checkpoint().unwrap(), checkpoint);
        for _ in 0..4 {
            assert_eq!(original.step(None).unwrap(), resumed.step(None).unwrap());
            assert_eq!(original.digest_hex(), resumed.digest_hex());
        }
        let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
        scenario.max_ticks = 6;
        let envelope = if distributed {
            record_replay_with_material_distributed_spatial_wind(scenario, wind(), medium, 3)
        } else {
            record_replay_with_material_spatial_wind(scenario, wind(), medium, 3)
        }
        .unwrap();
        assert_eq!(envelope.spec.rho_kg_m3, medium.density);
        let decoded = ReplayEnvelope::from_bytes(&envelope.to_bytes()).unwrap();
        let mut direct = if distributed {
            SimLoop::init_with_material_distributed_spatial_wind(
                decoded.spec.clone(),
                wind(),
                medium,
            )
        } else {
            SimLoop::init_with_material_spatial_wind(decoded.spec.clone(), wind(), medium)
        }
        .unwrap();
        assert_eq!(
            decoded.scrub_to_tick(0).unwrap(),
            direct.replay_initial_state().unwrap()
        );
        for tick in 1..=6 {
            assert_eq!(
                decoded.scrub_to_tick(tick).unwrap(),
                direct.step(None).unwrap()
            );
        }
    }
}

#[test]
fn gas_transport_identity_is_retained_even_when_current_polars_do_not_consume_it() {
    let a = initialize(false, gas(100_000.0));
    let mut law = GasSpec::dry_air_ussa1976();
    law.sutherland_beta *= 1.1;
    let changed = GasState::try_new(&law, 288.15, 100_000.0).unwrap();
    let b = initialize(false, changed);
    assert_eq!(a.scenario().rho_kg_m3, b.scenario().rho_kg_m3);
    assert_eq!(a.tick0().trim, b.tick0().trim);
    assert_ne!(a.tick0().digest, b.tick0().digest);
    assert_ne!(a.run_intent_id, b.run_intent_id);
    assert_ne!(a.save_checkpoint().unwrap(), b.save_checkpoint().unwrap());
    let point = [0.0, 0.0, 2.4];
    assert_ne!(
        a.sample_material_air(point)
            .unwrap()
            .reynolds(1.0, [0.0; 3])
            .unwrap(),
        b.sample_material_air(point)
            .unwrap()
            .reynolds(1.0, [0.0; 3])
            .unwrap()
    );
}

#[test]
fn incoherent_gas_and_contradictory_restoration_refuse_before_publication() {
    let medium = gas(100_000.0);
    for field in 0..3 {
        let mut bad = medium;
        match field {
            0 => bad.density *= 1.2,
            1 => bad.dynamic_viscosity = f64::NAN,
            _ => bad.water_mole_fraction = -0.1,
        }
        let result = SimLoop::init_with_material_spatial_wind(
            dec17_scenario(1903, PilotMode::FixedControls),
            wind(),
            bad,
        );
        assert!(matches!(result, Err(e) if e.code == "material-atmosphere-invalid"));
    }
    let result = SimLoop::init_with_material_spatial_wind(
        dec17_scenario(1903, PilotMode::FixedControls),
        wind(),
        gas(10_000.0),
    );
    assert!(matches!(result, Err(e) if e.code == "scenario-invalid"));
    let sim = initialize(false, medium);
    let bytes = sim.save_checkpoint().unwrap();
    assert!(
        matches!(SimLoop::restore_checkpoint(dec17_scenario(1903, PilotMode::FixedControls),
        &bytes), Err(e) if e.code == "checkpoint-malformed")
    );
    let sample = sim.sample_material_air([0.0, 0.0, 2.4]).unwrap();
    assert!(sample.reynolds(0.0, [0.0; 3]).is_err());
    assert!(
        sample
            .dynamic_pressure_pa([f64::INFINITY, 0.0, 0.0])
            .is_err()
    );
}

#[test]
fn material_checkpoint_admits_its_mode_gas_and_identity_before_restoring_state() {
    let sim = initialize(true, gas(100_000.0));
    let checkpoint = sim.save_checkpoint().unwrap();
    let body = &checkpoint[..checkpoint.len() - 32];
    let mut wind = Vec::new();
    sim.spatial_wind
        .as_ref()
        .unwrap()
        .encode(&sim.run_intent_id, &mut wind);
    let gas_start = body.len() - wind.len() - 96;
    let seal = |mut body: Vec<u8>, domain| {
        let digest = hash_domain(domain, &body);
        body.extend_from_slice(digest.as_bytes());
        body
    };
    let restore = |bytes: &[u8]| SimLoop::restore_checkpoint(sim.scenario().clone(), bytes);
    let domain = "org.frankensim.wf.checkpoint-state.v5";

    let mut mode = body.to_vec();
    mode[gas_start - 1] = 2;
    assert!(matches!(restore(&seal(mode, domain)), Err(e) if e.code == "checkpoint-malformed"));

    let mut invalid_gas = body.to_vec();
    invalid_gas[gas_start..gas_start + 8].copy_from_slice(&0.0_f64.to_le_bytes());
    assert!(matches!(restore(&seal(invalid_gas, domain)),
        Err(e) if e.code == "material-atmosphere-invalid"));

    // Coherent but different transport is still a different run, even though
    // the currently used force density has not changed.
    let mut law = GasSpec::dry_air_ussa1976();
    law.sutherland_beta *= 1.1;
    let changed = GasState::try_new(&law, 288.15, 100_000.0).unwrap();
    let mut encoded_gas = Vec::new();
    encode_gas(&changed, &mut encoded_gas);
    let mut changed_gas = body.to_vec();
    changed_gas[gas_start..gas_start + 96].copy_from_slice(&encoded_gas);
    assert!(matches!(restore(&seal(changed_gas, domain)),
        Err(e) if e.code == "checkpoint-malformed"));

    assert!(
        matches!(restore(&seal(body[..gas_start + 80].to_vec(), domain)),
        Err(e) if e.code == "checkpoint-malformed")
    );

    // Removing the material extension cannot turn the same retained run into
    // a historical scalar-density V4 trajectory merely by resealing bytes.
    let mut downgraded = body.to_vec();
    downgraded.drain(gas_start - 1..gas_start + 96);
    downgraded[..4].copy_from_slice(&4_u32.to_le_bytes());
    assert!(
        matches!(restore(&seal(downgraded, "org.frankensim.wf.checkpoint-state.v4")),
        Err(e) if e.code == "checkpoint-malformed")
    );
}
