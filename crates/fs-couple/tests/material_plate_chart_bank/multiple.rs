//! Exercise added bows against the real heterogeneous body and shared cavity.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{MaterialBowContact, MultiBowStringBank};
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionRegime;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::{MultiFrictionConfig, MultiFrictionModalSystem};
use fs_couple::modal_acoustic_time::ModalAcousticState;

fn joint() -> MultiFrictionConfig {
    MultiFrictionConfig { max_ports: 4, max_sweeps: 128, max_setup_terms: 16_384 }
}
fn contact() -> MaterialBowContact {
    MaterialBowContact { location: StringBowLocation { string: 1, station_fraction: 0.31 },
        law: FrictionLaw::Coulomb { static_mu: 0.7, kinetic_mu: 0.3 }, friction: options().friction }
}
fn model(cx: &Cx<'_>, gate: &CancelGate) -> MultiBowStringBank {
    build(cx, gate, body(), gas(101_325.0)).with_additional_bows(vec![contact()], joint(), gate).unwrap()
}
fn states(s: &MultiFrictionModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m| m.states().iter().copied()).collect()
}
fn drivers(sample: usize) -> [FrictionDrive; 2] {
    [FrictionDrive { speed_m_s: if sample < 96 { 0.01 } else { -0.01 }, normal_force_n: 0.1 },
        FrictionDrive { speed_m_s: -0.02, normal_force_n: if sample < 144 { 0.07 } else { 0.0 } }]
}

#[test]
fn empty_extension_retains_original_bow_sources_vibration_pressure_and_exact_steps() {
    with_cx(|cx, gate| {
        let mut expected = build(cx, gate, body(), gas(101_325.0));
        let mut original = build(cx, gate, body(), gas(101_325.0));
        let zero = vec![0.0; original.system().mode_count()];
        for _ in 0..19 {
            assert_eq!(original.step(&zero, drive()).unwrap(), expected.step(&zero, drive()).unwrap());
        }
        let identity = original.source().material_chart().unwrap().specimen.specimen_identity();
        let weights = original.system().acceleration_pressure_weights().unwrap().to_vec();
        let mut actual = original.with_additional_bows(vec![], joint(), gate).unwrap();
        assert_eq!(actual.system().samples_rendered(), 19);
        assert_eq!(actual.system().acceleration_pressure_weights().unwrap(), weights);
        assert_eq!(actual.source().bank().material_chart().unwrap().specimen.specimen_identity(), identity);
        assert_eq!(actual.source().bank().gas(), actual.source().bank().material_cavity().unwrap().cavity.gas());
        assert_eq!(actual.source().bows().len(), 1);
        assert_eq!(actual.source().bows()[0].location, options().bow);
        assert_eq!(actual.source().bows()[0].law, law());
        for i in 0..192 {
            let d = drivers(i)[0];
            let a = actual.step(&zero, &[d]).unwrap();
            let e = expected.step(&zero, d).unwrap();
            assert_eq!(a.sample, e.sample);
            assert_eq!(a.observer_pressure_pa.to_bits(), e.observer_pressure_pa.to_bits());
            assert_eq!(a.energy_residual_j, e.energy_residual_j);
            assert_eq!(a.ports[0].traction_n, e.traction_n);
            assert_eq!(a.ports[0].drive_work_j, e.drive_work_j);
            assert_eq!(a.ports[0].dissipation_j, e.friction_dissipation_j);
            for (a, e) in actual.system().components().iter().zip(expected.system().components()) {
                assert_eq!(a.states(), e.states());
            }
        }
    });
}

#[test]
fn additional_bow_uses_its_own_material_basis_and_matches_the_joint_owner_exactly() {
    with_cx(|cx, gate| {
        // Obtain the selected contact independently through the original
        // single-bow constructor. No duplicated shape formula in this oracle.
        let bow = contact(); let mut o = options(); o.bow = bow.location;
        let gas = gas(101_325.0);
        let selected = ReciprocalStringBank::from_material_plate_chart(cx, sources(), gas,
            body(), Some(cavity(gas)), bow.law.clone(), o, gate).unwrap();
        let port = selected.system().port().clone();
        let mut reference = build(cx, gate, body(), gas).into_parts().1
            .with_additional_ports(vec![(port.clone(), bow.friction)], joint(), gate).unwrap();
        let mut actual = model(cx, gate);
        assert_eq!(actual.source().bows()[1], bow);
        assert_eq!(actual.system().port(1).unwrap().left, port.left);
        assert_eq!(actual.system().port(1).unwrap().law, port.law);
        assert_eq!(actual.system().port(1).unwrap().left.shapes.len(), 3);
        assert_eq!(actual.system().port(0).unwrap().left.shapes.len(), 4);
        let zero = vec![0.0; actual.system().mode_count()];
        for i in 0..192 {
            let ds = drivers(i);
            assert_eq!(actual.step(&zero, &ds).unwrap(), reference.step(&zero, &ds).unwrap());
            let f = actual.system().last_frame().unwrap();
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
            assert_eq!(f.drive_work_j, f.ports.iter().map(|p| p.drive_work_j).sum::<f64>());
            assert!(f.ports.iter().all(|p| p.dissipation_j >= 0.0));
            assert_eq!(states(actual.system()), states(&reference));
        }
    });
}

#[test]
fn second_bow_reacts_on_the_first_string_and_release_preserves_sympathetic_motion() {
    with_cx(|cx, gate| {
        let mut released = model(cx, gate); let mut driven = model(cx, gate);
        let zero = vec![0.0; driven.system().mode_count()];
        let off = FrictionDrive { speed_m_s: -0.02, normal_force_n: 0.0 };
        let on = FrictionDrive { speed_m_s: -0.02, normal_force_n: 0.07 };
        let mut difference = 0.0;
        for _ in 0..256 {
            let f = released.step(&zero, &[drive(), off]).unwrap();
            assert_eq!(f.ports[1].regime, ModalFrictionRegime::Released);
            driven.step(&zero, &[drive(), on]).unwrap();
            difference += (released.system().components()[0].states()[0].displacement_m_sqrt_kg
                - driven.system().components()[0].states()[0].displacement_m_sqrt_kg).abs();
        }
        assert!(difference > 1e-14, "second contact did not load the first string: {difference}");
        assert!(released.system().components()[1].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-10));
        let before = driven.system().total_energy_j().unwrap();
        for _ in 0..128 {
            let f = driven.step(&zero, &[off; 2]).unwrap();
            assert_eq!(f.drive_work_j, 0.0);
            assert_eq!(f.friction_dissipation_j, 0.0);
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        }
        assert!(driven.system().total_energy_j().unwrap() <= before + 1e-10);
    });
}

#[test]
fn rejected_joint_samples_retain_every_state_source_and_clock_then_retry_exactly() {
    with_cx(|cx, gate| {
        let mut actual = model(cx, gate); let mut reference = model(cx, gate);
        let zero = vec![0.0; actual.system().mode_count()]; let ds = drivers(0);
        actual.step(&zero, &ds).unwrap(); reference.step(&zero, &ds).unwrap();
        let before = states(actual.system()); let frame = actual.system().last_frame().cloned();
        let identity = actual.source().bank().material_chart().unwrap().specimen.specimen_identity();
        let stop = CancelGate::new_clock_free(); stop.request();
        assert!(actual.step_under_gate(&zero, &ds, &stop).is_err());
        assert!(actual.step(&zero, &ds[..1]).is_err());
        assert!(actual.step(&zero, &[ds[0], FrictionDrive { speed_m_s: f64::NAN, normal_force_n: 0.1 }]).is_err());
        assert!(actual.step(&zero, &[ds[0], FrictionDrive { speed_m_s: 0.0, normal_force_n: -0.1 }]).is_err());
        assert_eq!(states(actual.system()), before);
        assert_eq!(actual.system().last_frame(), frame.as_ref());
        assert_eq!(actual.system().samples_rendered(), 1);
        assert_eq!(actual.source().bank().material_chart().unwrap().specimen.specimen_identity(), identity);
        assert_eq!(actual.step(&zero, &ds).unwrap(), reference.step(&zero, &ds).unwrap());
    });
}

#[test]
fn added_contact_geometry_laws_and_setup_limits_are_admitted_before_use() {
    with_cx(|cx, gate| {
        for defect in 0..9 {
            let original = build(cx, gate, body(), gas(101_325.0));
            let mut bow = contact(); let mut cfg = joint();
            match defect {
                0 => bow.location.string = 2,
                1 => bow.location.station_fraction = 0.0,
                2 => bow.location.station_fraction = 1.0,
                3 => bow.location.station_fraction = f64::NAN,
                4 => cfg.max_ports = 1,
                5 => cfg.max_setup_terms = 0,
                6 => cfg.max_sweeps = 0,
                7 => bow.friction.max_iterations = 0,
                _ => bow.law = FrictionLaw::Coulomb { static_mu: -0.1, kinetic_mu: 0.3 },
            }
            assert!(original.with_additional_bows(vec![bow], cfg, gate).is_err(), "defect {defect}");
        }
        let original = build(cx, gate, body(), gas(101_325.0));
        let stop = CancelGate::new_clock_free(); stop.request();
        assert!(original.with_additional_bows(vec![contact()], joint(), &stop).is_err());
    });
}

#[test]
fn layered_material_changes_additional_bow_shapes_without_rewriting_friction_or_other_sources() {
    use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
    use fs_couple::string_specimen::StringPrestress;
    with_cx(|cx, gate| {
        let make = |density_scale: f64| {
            let mut s = sources();
            let template = match &s[1] { MaterialStringSource::Homogeneous(s) => s.string().clone(), _ => unreachable!() };
            // A section owns its loss, so clear the inherited homogeneous loss.
            let template = PrestressedString { kelvin_voigt_bending: None, ..template };
            let section = ConcentricElasticSection::try_new(0.0, vec![
                ConcentricMaterialLayer { outer_radius_m: 0.0002, material: material(density_scale*7800.0, 200e9) },
                ConcentricMaterialLayer { outer_radius_m: 0.0003, material: material(density_scale*3900.0, 100e9) },
            ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
            s[1] = MaterialStringSource::Concentric { template, section, prestress: StringPrestress::FixedTension(20.0) };
            let gas = gas(101_325.0);
            ReciprocalStringBank::from_material_plate_chart(cx, s, gas, body(), Some(cavity(gas)), law(), options(), gate)
                .unwrap().with_additional_bows(vec![contact()], joint(), gate).unwrap()
        };
        let mut a = make(1.0); let mut b = make(2.0);
        assert_eq!(a.source().bows(), b.source().bows());
        assert_eq!(a.cards()[0].fundamental_hz(), b.cards()[0].fundamental_hz());
        assert_eq!(a.source().bank().material_chart().unwrap().specimen.specimen_identity(),
            b.source().bank().material_chart().unwrap().specimen.specimen_identity());
        for (&a, &b) in a.system().port(1).unwrap().left.shapes.iter().zip(&b.system().port(1).unwrap().left.shapes) {
            assert!((b-a/sqrt(2.0)).abs() < 1e-12);
        }
        let MaterialStringSource::Concentric { section, .. } = &b.source().bank().strings()[1].source else { panic!("lost layered sources") };
        assert_eq!(section.elastic().layers().len(), 2);
        assert_eq!(section.elastic().layers()[0].material.identity(), material(15_600.0, 200e9).identity());
        let zero = vec![0.0; a.system().mode_count()]; let mut difference = 0.0;
        for i in 0..192 {
            difference += (a.step(&zero, &drivers(i)).unwrap().observer_pressure_pa
                - b.step(&zero, &drivers(i)).unwrap().observer_pressure_pa).abs();
        }
        assert!(difference > 1e-12);
    });
}
