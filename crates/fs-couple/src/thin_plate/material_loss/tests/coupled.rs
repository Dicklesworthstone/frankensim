//! Source viscosity must reach actual common-body motion and scheduled pressure.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::{MaterialStringSource,
    ReciprocalStringBank, StringBankOptions, StringBowLocation,
    chart::SharedMaterialPlateChartBody, plate::SharedPlateBridge};
use crate::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity;
use crate::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use crate::modal_acoustic_time::ModalAcousticTimeBudget;
use crate::render::schedule::force::coupled::ModalCouplingConfig;
use crate::render::schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionConfig};
use crate::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use crate::string_specimen::with_uniform_circular_material_state;
use fs_exec::{CancelGate, Cx};
use fs_material::gas::GasSpec;
use fs_scenario::{HelmholtzCavity, PrestressedString,
    gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
use fs_tribo::FrictionLaw;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena, fs_exec::StreamKey {
        seed: 7, kernel_id: 83, tile: 0, iteration: 0,
    }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate));
}
fn build(cx: &Cx<'_>, gate: &CancelGate, eta: f64) -> ReciprocalStringBank {
    let initial = specimen(material(200e9, Some(2e6), 100_000.0), material(100e9, Some(1e6), 100_000.0))
        .with_kelvin_voigt_bending_loss().unwrap();
    let body = SharedMaterialPlateChartBody {
        specimen: initial.with_region_material("right", &material(100e9, Some(eta), 100_000.0)).unwrap(),
        reduction: options(), bridges: vec![SharedPlateBridge { station_fraction: 0.8,
            stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0 }; 2],
        max_plate_dofs: 75, max_plate_elements: 32,
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0,
    };
    let strings = [4, 3].into_iter().map(|n_modes| MaterialStringSource::Homogeneous(
        with_uniform_circular_material_state(PrestressedString {
            length_m: 0.5, tension_n: 20.0, lin_density_kg_m: 1.0, axial_stiffness_n: 0.0,
            width_m: 0.01, n_modes, damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None, relaxation_bending: None, polarization_detune: 0.0, moving_end: false,
            tension_loss_factor: None,
            bending_loss_factor: None,
        }, 0.0003, &material(200e9, Some(10.0), 100_000.0)).unwrap().with_kelvin_voigt_bending_loss().unwrap()
    )).collect();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let cavity = MaterialBodyCavity { cavity: HelmholtzModalCavity::new(HelmholtzCavity {
        volume_m3: 0.004, neck_radius_m: 0.01, neck_length_m: 0.02,
        inner_end: fs_scenario::NeckEnd::Flanged, outer_end: fs_scenario::NeckEnd::Flanged,
    }, gas, CavityNeckLoss::RadiationAtResonance).unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.5 };
    let o = StringBankOptions { sample_rate_hz: 48_000,
        bow: StringBowLocation { string: 0, station_fraction: 0.2 },
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig { max_modes: 32, max_connections: 8, max_setup_terms: 16_384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 100.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8 },
        friction: ModalFrictionConfig { max_iterations: 96, maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0, force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9, sticking_tolerance_m: 1e-12 } };
    ReciprocalStringBank::from_material_plate_chart(cx, strings, gas, body, Some(cavity),
        FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 }, o, gate).unwrap()
}
fn drive() -> FrictionDrive { FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 } }

#[test]
fn regional_viscosity_reacts_back_on_strings_without_changing_gas_or_elastic_ports() {
    with_cx(|cx, gate| {
        let mut a = build(cx, gate, 1e6); let mut b = build(cx, gate, 4e6);
        assert_eq!(a.source().gas(), b.source().gas());
        assert_eq!(a.system().acceleration_pressure_weights(), b.system().acceleration_pressure_weights());
        assert_eq!(a.cards()[0].zetas, b.cards()[0].zetas);
        for (a, b) in a.system().components()[2].modes().iter().zip(b.system().components()[2].modes()) {
            assert_eq!(a.angular_frequency_rad_s, b.angular_frequency_rad_s);
            assert!(b.damping_ratio > a.damping_ratio);
        }
        let zero = vec![0.0; a.system().mode_count()];
        let (mut pressure_difference, mut string_difference) = (0.0, 0.0);
        for _ in 0..512 {
            pressure_difference += (a.step(&zero, drive()).unwrap().observer_pressure_pa
                - b.step(&zero, drive()).unwrap().observer_pressure_pa).abs();
            for model in [&a, &b] {
                let f = model.system().last_frame().unwrap();
                assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
            }
            string_difference += (a.system().components()[0].states()[0].displacement_m_sqrt_kg
                - b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs();
        }
        assert!(pressure_difference > 1e-12);
        assert!(string_difference > 1e-14);
        let before = b.system().total_energy_j().unwrap();
        for _ in 0..128 {
            let f = b.step(&zero, FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 }).unwrap();
            assert_eq!(f.drive_work_j, 0.0);
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        }
        assert!(b.system().total_energy_j().unwrap() <= before + 1e-10);
    });
}

#[test]
fn scheduled_pressure_preserves_regional_loss_sources_and_direct_sample_bits() {
    with_cx(|cx, gate| {
        let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.1, station: 0.2 }, events: vec![],
        }]).unwrap();
        let config = FrictionGestureConfig { sample_rate_hz: 48_000, samples: 129,
            max_block: 128, max_work: 1000, max_events: 16, station_fraction: 0.2 };
        let mut direct = build(cx, gate, 4e6);
        let source = direct.source().material_chart().unwrap().specimen.clone();
        let zero = vec![0.0; direct.system().mode_count()];
        let expected: Vec<_> = (0..129).map(|_| direct.step(&zero, drive()).unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1, 37, 128] {
            let mut p = build(cx, gate, 4e6).into_scheduled(&schedule, "bow", config).unwrap();
            let retained = &p.source().material_chart().unwrap().specimen;
            assert_eq!(retained.specimen_identity(), source.specimen_identity());
            assert_eq!(retained.bending_viscosity(), source.bending_viscosity());
            let mut output = vec![0.0; 129];
            for block in output.chunks_mut(partition) { p.renderer_mut().block(block).unwrap(); }
            assert_eq!(output.iter().map(|p| p.to_bits()).collect::<Vec<_>>(), expected);
        }
    });
}
