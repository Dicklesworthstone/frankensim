//! Synthetic fixtures exercise shared mechanical state, not calibrated timbre.
use fs_couple::bowed_string::runtime::reciprocal::{
    InteriorStringConnection, ReciprocalMaterialString, ReciprocalStringOptions,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringMember, MaterialStringSource, ReciprocalStringBank, StringBankOptions, StringBowLocation,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionConfig};
use fs_couple::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use fs_couple::string_specimen::{
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA, ResolvedStringSpecimen,
    StringPrestress, with_uniform_circular_material_state,
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_math::c64::C64;
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{AmbientGas, PrestressedString, gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
use fs_tribo::FrictionLaw;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena, fs_exec::StreamKey {
        seed: 6, kernel_id: 5, tile: 0, iteration: 0,
    }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate))
}
fn material(rho: f64, young: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [("density", Density::DIMS, rho), ("young_modulus", Pressure::DIMS, young),
        ("poisson_ratio", Dims::NONE, 0.3), (ETA, DynViscosity::DIMS, 10.0)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic shared-string test".into(), license: "CC0-1.0".into(), artifact: None },
            observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "fixture".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn template() -> PrestressedString {
    PrestressedString { length_m: 0.5, tension_n: 20.0, lin_density_kg_m: 1.0, axial_stiffness_n: 1.0,
        width_m: 1.0, n_modes: 4, damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None, relaxation_bending: None, polarization_detune: 0.0, moving_end: false,
        tension_loss_factor: None, bending_loss_factor: None }
}
fn specimen(rho: f64, young: f64) -> ResolvedStringSpecimen {
    with_uniform_circular_material_state(template(), 0.0003, &material(rho, young)).unwrap()
        .with_kelvin_voigt_bending_loss().unwrap()
}
fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn connection(body: usize) -> InteriorStringConnection {
    InteriorStringConnection { station_fraction: 0.8,
        body: ModalAttachment { component: body, shapes: vec![1.0 / 0.2_f64.sqrt()] },
        stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0 }
}
fn member(rho: f64) -> MaterialStringMember {
    MaterialStringMember { source: MaterialStringSource::Homogeneous(specimen(rho, 200e9)), connections: vec![connection(0)] }
}
fn layered() -> MaterialStringMember {
    let section = ConcentricElasticSection::try_new(0.0, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0002, material: material(7800.0, 200e9) },
        ConcentricMaterialLayer { outer_radius_m: 0.0003, material: material(15_600.0, 100e9) },
    ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
    MaterialStringMember { source: MaterialStringSource::Concentric { template: template(), section,
        prestress: StringPrestress::FixedTension(20.0) }, connections: vec![connection(0)] }
}
fn body(velocity: f64) -> ModalAcousticTimeModel {
    let mut body = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: (40_000.0 / 0.2_f64).sqrt(), damping_ratio: 0.01,
        pressure_per_modal_velocity: C64::new(1.0, 0.0), // authored synthetic observation
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    body.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: velocity }]).unwrap();
    body
}
fn options(string: usize) -> StringBankOptions {
    StringBankOptions { sample_rate_hz: 48_000, bow: StringBowLocation { string, station_fraction: 0.2 },
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig { max_modes: 32, max_connections: 8, max_setup_terms: 16_384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 100.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8 },
        friction: ModalFrictionConfig { max_iterations: 96, maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0, force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9, sticking_tolerance_m: 1e-12 } }
}
fn law() -> FrictionLaw { FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } }
fn build(cx: &Cx<'_>, gate: &CancelGate, passive_rho: f64) -> ReciprocalStringBank {
    ReciprocalStringBank::new(cx, vec![member(7800.0), member(passive_rho)], gas(101_325.0),
        vec![body(0.0)], law(), options(0), gate).unwrap()
}
fn drive() -> FrictionDrive { FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 } }

#[test]
fn one_string_is_bit_identical_to_the_existing_shared_assembly() {
    with_cx(|cx, gate| {
        let o = options(0);
        let old_options = ReciprocalStringOptions { sample_rate_hz: o.sample_rate_hz,
            bow_station_fraction: o.bow.station_fraction, connections: vec![connection(0)],
            string_budget: o.string_budget, network: o.network, friction: o.friction };
        let mut old = ReciprocalMaterialString::from_specimen(cx, specimen(7800.0, 200e9),
            AmbientGas { temperature_k: 293.15, pressure_pa: 101_325.0, relative_humidity: 0.0 },
            vec![body(0.01)], law(), old_options, gate).unwrap();
        let mut bank = ReciprocalStringBank::new(cx, vec![member(7800.0)], gas(101_325.0),
            vec![body(0.01)], law(), o, gate).unwrap();
        let zero = vec![0.0; bank.system().mode_count()];
        for _ in 0..128 { assert_eq!(old.system_mut().step(&zero, drive()).unwrap(), bank.step(&zero, drive()).unwrap()); }
    });
}

#[test]
fn unplayed_string_vibrates_and_its_material_reacts_back_on_the_played_string() {
    with_cx(|cx, gate| {
        let mut a = build(cx, gate, 7800.0);
        let mut b = build(cx, gate, 15_600.0);
        assert_eq!(a.cards()[0].zetas, b.cards()[0].zetas);
        assert_eq!(a.cards()[0].fundamental_hz(), b.cards()[0].fundamental_hz());
        let zero = vec![0.0; a.system().mode_count()];
        let mut played_difference = 0.0;
        for _ in 0..512 {
            for model in [&mut a, &mut b] {
                let frame = model.step(&zero, drive()).unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                assert!(frame.friction_dissipation_j >= 0.0);
            }
            played_difference += (a.system().components()[0].states()[0].displacement_m_sqrt_kg
                - b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs();
        }
        assert!(played_difference > 1e-10, "passive string did not load the driven string: {played_difference:e}");
        assert!(a.system().components()[1].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-8));
        let before = a.system().total_energy_j().unwrap();
        for _ in 0..256 {
            let f = a.step(&zero, FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }).unwrap();
            assert_eq!(f.drive_work_j, 0.0);
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        }
        assert!(a.system().total_energy_j().unwrap() <= before + 1e-12);
    });
}

#[test]
fn body_indices_shift_by_all_strings_and_existing_body_vibration_is_preserved() {
    with_cx(|cx, gate| {
        let mut members = vec![member(7800.0), member(9000.0), layered()];
        for m in &mut members { m.connections[0] = connection(1); }
        let mut model = ReciprocalStringBank::new(cx, members, gas(101_325.0), vec![body(0.0), body(0.01)],
            law(), options(2), gate).unwrap();
        assert_eq!(model.system().components().len(), 5);
        assert_eq!(model.system().components()[4].states()[0].velocity_m_sqrt_kg_per_s, 0.01);
        assert_eq!(model.system().port().left.component, 2);
        let zero = vec![0.0; model.system().mode_count()];
        for _ in 0..128 { model.step(&zero, FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }).unwrap(); }
        assert_eq!(model.system().components()[3].states()[0].velocity_m_sqrt_kg_per_s, 0.0);
        for component in &model.system().components()[..3] {
            assert!(component.states().iter().any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-8));
            assert!(component.modes().iter().all(|m| m.pressure_per_modal_velocity == C64::ZERO));
        }
    });
}

#[test]
fn gas_and_layer_sources_apply_to_every_member_without_retuning_or_erasing_receipts() {
    with_cx(|cx, gate| {
        let members = vec![member(7800.0), layered()];
        let a = ReciprocalStringBank::new(cx, members.clone(), gas(101_325.0), vec![body(0.0)], law(), options(0), gate).unwrap();
        let b = ReciprocalStringBank::new(cx, members, gas(50_000.0), vec![body(0.0)], law(), options(0), gate).unwrap();
        for (a, b) in a.cards().iter().zip(b.cards()) {
            assert_eq!(a.fundamental_hz(), b.fundamental_hz());
            assert!(a.zetas.iter().zip(&b.zetas).all(|(a,b)| a > b));
        }
        let MaterialStringSource::Concentric { section, .. } = &b.source().strings()[1].source else { panic!("lost layered source") };
        assert_eq!(section.elastic().layers().len(), 2);
        assert_eq!(b.cards()[1].linear_density_kg_m, section.elastic().line_density_kg_m());
        assert_eq!(b.cards()[1].viscous_bending_n_m2_s, section.viscous_stiffness_n_m2_s());
        assert_eq!(*b.source().gas(), gas(50_000.0));
    });
}

#[test]
fn scheduling_is_block_invariant_and_refuses_the_wrong_string() {
    with_cx(|cx, gate| {
        let schedule = |string| GestureSchedule::try_new(1000, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.1, station: 0.2 }, events: vec![],
        }]).unwrap();
        let config = FrictionGestureConfig { sample_rate_hz: 48_000, samples: 129, max_block: 128,
            max_work: 1000, max_events: 16, station_fraction: 0.2 };
        let mut direct = build(cx, gate, 7800.0);
        let zero = vec![0.0; direct.system().mode_count()];
        let expected: Vec<_> = (0..129).map(|_| direct.step(&zero, drive()).unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1, 37, 128] {
            let mut performer = build(cx, gate, 7800.0).into_scheduled(&schedule(0), "bow", config).unwrap();
            assert_eq!(performer.source().strings().len(), 2);
            let mut out = vec![0.0; 129];
            for block in out.chunks_mut(partition) { performer.renderer_mut().block(block).unwrap(); }
            assert_eq!(out.iter().map(|p| p.to_bits()).collect::<Vec<_>>(), expected);
        }
        assert!(build(cx, gate, 7800.0).into_scheduled(&schedule(1), "bow", config).is_err());
    });
}

#[test]
fn detached_members_wrong_geometry_budgets_gas_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        for which in 0..6 {
            let mut members = vec![member(7800.0), layered()];
            let mut o = options(0);
            let mut gas = gas(101_325.0);
            match which {
                0 => members[1].connections.clear(),
                1 => members[1].connections[0].station_fraction = 1.0,
                2 => members[1].connections[0].body.component = 1,
                3 => o.network.max_modes = 8,
                4 => o.bow.string = 2,
                _ => gas.density = 0.0,
            }
            assert!(ReciprocalStringBank::new(cx, members, gas, vec![body(0.0)], law(), o, gate).is_err());
        }
        let mut bank = build(cx, gate, 7800.0);
        let zero = vec![0.0; bank.system().mode_count()];
        let previous = *bank.step(&zero, drive()).unwrap();
        gate.request();
        assert!(bank.step_under_gate(&zero, drive(), gate).is_err());
        assert_eq!(bank.system().last_frame(), Some(&previous));
        assert!(ReciprocalStringBank::new(cx, vec![member(7800.0)], gas(101_325.0), vec![body(0.0)], law(), options(0), gate).is_err());
    });
}

#[path = "material_string_bank/plate.rs"]
mod plate;

#[path = "material_string_bank/shell.rs"]
mod shell;
