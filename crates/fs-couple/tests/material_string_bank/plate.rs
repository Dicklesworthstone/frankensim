use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{SharedMaterialPlateBody, SharedPlateBridge};
use fs_couple::thin_plate::{PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state};

fn source_plate(rho: f64, count: usize) -> SharedMaterialPlateBody {
    let specimen = with_uniform_plate_material_state(fs_scenario::ThinPlate {
        length_m: 0.4, width_m: 0.3, thickness_m: 1.0, density_kg_m3: 1.0,
        e1_pa: 1.0, e2_pa: 1.0, nu12: 0.0, g12_pa: 1.0, material_angle_rad: 0.0,
        damping_ratio: 0.0, thermoelastic: None, kelvin_voigt_bending: None, n_modes: 2,
        geometric_nonlinearity: false, pretension_n_m: 0.0, clamped: false,
    }, &material(rho, 200e9), PlateMaterialModel::Isotropic, PlateThicknessConstraint::FixedThickness(0.001))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap();
    SharedMaterialPlateBody { specimen, bridges: vec![SharedPlateBridge {
        station_fraction: 0.8, stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0,
    }; count], modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0 }
}
fn sources() -> Vec<MaterialStringSource> {
    vec![member(7800.0).source, layered().source]
}
fn bank(cx: &Cx<'_>, gate: &CancelGate, rho: f64, pressure: f64) -> ReciprocalStringBank {
    ReciprocalStringBank::from_material_plate(cx, sources(), gas(pressure), source_plate(rho, 2),
        law(), options(0), gate).unwrap()
}

#[test]
fn one_material_plate_string_matches_the_retained_single_string_image() {
    with_cx(|cx, gate| {
        let plate = source_plate(7800.0, 1);
        let old_body = fs_couple::bowed_string::runtime::reciprocal::plate::MaterialPlateBody {
            specimen: plate.specimen.clone(), bridge_station_fraction: 0.8,
            stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0,
            modal_budget: plate.modal_budget, listener_m: plate.listener_m,
        };
        let o = options(0);
        let old_options = ReciprocalStringOptions { sample_rate_hz: o.sample_rate_hz,
            bow_station_fraction: 0.2, connections: vec![], string_budget: o.string_budget,
            network: o.network, friction: o.friction };
        let mut old = ReciprocalMaterialString::from_material_plate(cx, specimen(7800.0, 200e9),
            AmbientGas { temperature_k: 293.15, pressure_pa: 101_325.0, relative_humidity: 0.0 },
            old_body, law(), old_options, gate).unwrap();
        let mut model = ReciprocalStringBank::from_material_plate(cx, vec![member(7800.0).source],
            gas(101_325.0), plate, law(), o, gate).unwrap();
        let zero = vec![0.0; model.system().mode_count()];
        for _ in 0..128 { assert_eq!(old.system_mut().step(&zero, drive()).unwrap(), model.step(&zero, drive()).unwrap()); }
    });
}

#[test]
fn four_strings_share_one_material_body_and_one_acceleration_pressure_observer() {
    with_cx(|cx, gate| {
        let mut sources = sources(); sources.extend([member(9000.0).source, layered().source]);
        let body = source_plate(7800.0, 4);
        let source_id = body.specimen.material().identity();
        let expected = fs_couple::thin_plate::certified_radiators(body.specimen.plate()).unwrap();
        let mut model = ReciprocalStringBank::from_material_plate(cx, sources, gas(101_325.0), body,
            law(), options(0), gate).unwrap();
        assert_eq!(model.system().components().len(), 5); // NOT four separate body copies
        assert_eq!(model.source().material_plate().unwrap().specimen.material().identity(), source_id);
        let weights = model.system().acceleration_pressure_weights().unwrap();
        assert!(weights[..16].iter().all(|w| *w == 0.0));
        for (k, radiator) in expected.iter().enumerate() {
            assert_eq!(weights[16+k], radiator.radiate_mass_normalized(1.0, gas(101_325.0).density, 1.0));
            assert_eq!(model.system().components()[4].modes()[k].angular_frequency_rad_s, radiator.omega);
            for member in model.source().strings() {
                let b = member.connections[0].body.shapes[k];
                assert!((b-radiator.drive_participation/radiator.mass_kg.sqrt()).abs() < 1e-12);
            }
        }
        let zero = vec![0.0; model.system().mode_count()];
        let mut peak = 0.0_f64;
        for _ in 0..512 {
            let f = model.step(&zero, drive()).unwrap();
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
            peak = peak.max(f.observer_pressure_pa.abs());
        }
        assert!(peak > 1e-12);
        for string in &model.system().components()[1..4] {
            assert!(string.states().iter().any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-8));
        }
    });
}

#[test]
fn shared_body_material_and_gas_changes_propagate_through_all_mechanics_and_pressure() {
    with_cx(|cx, gate| {
        let mut a = bank(cx, gate, 7800.0, 101_325.0);
        let mut b = bank(cx, gate, 15_600.0, 101_325.0);
        let c = bank(cx, gate, 7800.0, 50_000.0);
        let wa = a.system().components()[2].modes()[0].angular_frequency_rad_s;
        let wb = b.system().components()[2].modes()[0].angular_frequency_rad_s;
        assert!((wb/wa - 0.5_f64.sqrt()).abs() < 1e-6);
        assert_eq!(c.system().components()[2].modes()[0].angular_frequency_rad_s, wa);
        for (a, c) in a.cards().iter().zip(c.cards()) {
            assert_eq!(a.fundamental_hz(), c.fundamental_hz());
            assert!(a.zetas.iter().zip(&c.zetas).all(|(a,c)| a > c));
        }
        let pa = a.system().acceleration_pressure_weights().unwrap();
        let pc = c.system().acceleration_pressure_weights().unwrap();
        for (&a, &c) in pa[8..].iter().zip(&pc[8..]) {
            assert!((c-a*gas(50_000.0).density/gas(101_325.0).density).abs() < 1e-12*a.abs().max(1.0));
        }
        let zero = vec![0.0; a.system().mode_count()];
        let mut mechanical_difference = 0.0; let mut pressure_difference = 0.0;
        for _ in 0..256 {
            let pa = a.step(&zero, drive()).unwrap().observer_pressure_pa;
            let pb = b.step(&zero, drive()).unwrap().observer_pressure_pa;
            pressure_difference += (pa-pb).abs();
            mechanical_difference += (a.system().components()[0].states()[0].displacement_m_sqrt_kg
                - b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs();
        }
        assert!(mechanical_difference > 1e-10);
        assert!(pressure_difference > 1e-10);
    });
}

#[test]
fn incomplete_bridge_assignment_invalid_observer_and_budget_never_publish_a_partial_body() {
    with_cx(|cx, gate| {
        for which in 0..5 {
            let mut body = source_plate(7800.0, 2); let mut o = options(0);
            match which {
                0 => { body.bridges.pop(); }
                1 => body.bridges[1].station_fraction = 1.0,
                2 => body.listener_m = f64::NAN,
                3 => body.bridges[1].stiffness_n_m = -1.0,
                _ => o.network.max_modes = 9,
            }
            assert!(ReciprocalStringBank::from_material_plate(cx, sources(), gas(101_325.0), body, law(), o, gate).is_err());
        }
        let body = source_plate(7800.0, 2); gate.request();
        assert!(ReciprocalStringBank::from_material_plate(cx, sources(), gas(101_325.0), body, law(), options(0), gate).is_err());
    });
}

fn cavity(volume: f64, gas: GasState) -> fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity {
    use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
    fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity {
        cavity: HelmholtzModalCavity::new(fs_scenario::HelmholtzCavity {
            volume_m3: volume, neck_radius_m: 0.01, neck_length_m: 0.02,
            inner_end: fs_scenario::NeckEnd::Flanged, outer_end: fs_scenario::NeckEnd::Flanged,
        }, gas, CavityNeckLoss::Resistance(100.0)).unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0,
    }
}

#[test]
fn one_string_cavity_image_is_identical_to_the_existing_cavity_constructor() {
    with_cx(|cx, gate| {
        let plate = source_plate(7800.0, 1);
        let old_body = fs_couple::bowed_string::runtime::reciprocal::plate::MaterialPlateBody {
            specimen: plate.specimen.clone(), bridge_station_fraction: 0.8,
            stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0,
            modal_budget: plate.modal_budget, listener_m: plate.listener_m,
        };
        let o = options(0);
        let old_options = ReciprocalStringOptions { sample_rate_hz: o.sample_rate_hz,
            bow_station_fraction: 0.2, connections: vec![], string_budget: o.string_budget,
            network: o.network, friction: o.friction };
        let gas = gas(101_325.0); let cavity = cavity(0.003, gas);
        let mut old = ReciprocalMaterialString::from_material_plate_with_cavity(cx,
            specimen(7800.0, 200e9), gas, old_body, cavity.clone(), law(), old_options, gate).unwrap();
        let mut bank = ReciprocalStringBank::from_material_plate_with_cavity(cx,
            vec![member(7800.0).source], gas, plate, cavity, law(), o, gate).unwrap();
        let zero = vec![0.0; bank.system().mode_count()];
        for _ in 0..128 { assert_eq!(old.system_mut().step(&zero, drive()).unwrap(), bank.step(&zero, drive()).unwrap()); }
    });
}

#[test]
fn shared_cavity_loads_every_string_and_radiates_the_original_outward_neck_coordinate() {
    with_cx(|cx, gate| {
        let make = |volume| ReciprocalStringBank::from_material_plate_with_cavity(cx,
            sources(), gas(101_325.0), source_plate(7800.0, 2), cavity(volume, gas(101_325.0)),
            law(), options(0), gate).unwrap();
        let mut a = make(0.003); let mut b = make(0.006);
        assert_eq!(a.system().components().len(), 4); // 2 strings, ONE plate, ONE neck
        let c = a.source().material_cavity().unwrap();
        assert_eq!(a.source().body_connections().len(), 2); // spring plus selected resistance
        assert_eq!(c.cavity.geometry().volume_m3, 0.003);
        assert_eq!(c.cavity.gas(), a.source().gas());
        assert_eq!(*a.system().acceleration_pressure_weights().unwrap().last().unwrap(),
            gas(101_325.0).density * c.cavity.neck_volume_shape_m2_per_sqrt_kg()
                / (core::f64::consts::TAU * c.listener_m));
        let zero = vec![0.0; a.system().mode_count()];
        let mut mechanical_difference = 0.0; let mut pressure_difference = 0.0;
        for _ in 0..512 {
            let fa = *a.step(&zero, drive()).unwrap();
            let fb = *b.step(&zero, drive()).unwrap();
            for f in [fa, fb] { assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j); }
            mechanical_difference += (a.system().components()[1].states()[0].displacement_m_sqrt_kg
                - b.system().components()[1].states()[0].displacement_m_sqrt_kg).abs();
            pressure_difference += (fa.observer_pressure_pa - fb.observer_pressure_pa).abs();
        }
        assert!(mechanical_difference > 1e-10, "cavity did not react on the sympathetic string");
        assert!(pressure_difference > 1e-10);
        assert!(a.system().components()[3].states()[0].velocity_m_sqrt_kg_per_s.abs() > 1e-10);
    });
}

#[test]
fn shared_cavity_gas_mismatch_and_total_setup_overflow_refuse() {
    with_cx(|cx, gate| {
        let mut o = options(0); o.network.max_connections = 3; // 2 bridges + spring + resistance = 4
        assert!(ReciprocalStringBank::from_material_plate_with_cavity(cx, sources(), gas(101_325.0),
            source_plate(7800.0, 2), cavity(0.003, gas(101_325.0)), law(), o, gate).is_err());
        let mut o = options(0); o.network.max_modes = 10; // 8 string + 2 plate + 1 neck
        assert!(ReciprocalStringBank::from_material_plate_with_cavity(cx, sources(), gas(101_325.0),
            source_plate(7800.0, 2), cavity(0.003, gas(101_325.0)), law(), o, gate).is_err());
        assert!(ReciprocalStringBank::from_material_plate_with_cavity(cx, sources(), gas(50_000.0),
            source_plate(7800.0, 2), cavity(0.003, gas(101_325.0)), law(), options(0), gate).is_err());
    });
}
