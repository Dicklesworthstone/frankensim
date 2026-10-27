//! G0/G3: original bonded spectra use the same constitutive and mechanical owner.
use super::*;
use fs_material::section::relaxation::ConcentricPronySection;

fn section(materials: Vec<ResolvedMaterialStatePoint>) -> ConcentricPronySection {
    let count = materials.len();
    let layers = materials.into_iter().enumerate().map(|(i, material)| {
        let outer_radius_m = if count == 2 {
            if i == 0 { 0.0002 } else { 0.0003 }
        } else { 0.0003 * (i + 1) as f64 / count as f64 };
        ConcentricMaterialLayer { outer_radius_m, material }
    }).collect();
    ConcentricElasticSection::try_new(0.0, layers).unwrap().with_prony_bending(&vec![
        LayerBendingRelaxationProperties {
            equilibrium_modulus: "equilibrium_young_modulus".into(),
            branches: vec![("bending_delta".into(), "bending_time".into())],
        }; count
    ]).unwrap()
}

fn source(section: ConcentricPronySection, modes: usize) -> MaterialStringSource {
    MaterialStringSource::ConcentricProny {
        template: template(modes), section, prestress: StringPrestress::FixedTension(20.0),
    }
}

#[test]
fn exclusive_layer_spectra_match_the_independently_selected_parallel_arm_operator() {
    with_cx(|cx, gate| {
        let materials = [material(40e9, 0.002, 200_000.0), material(80e9, 0.004, 200_000.0)];
        let exclusive = source(section(materials.to_vec()), 2)
            .with_torsion("G", "eta_s", 1).unwrap()
            .with_transverse_polarizations(
                TransverseDirections::try_new([0.6, 0.8], [1.0, 0.0]).unwrap()).unwrap();
        let mut a = bank_with_sources(cx, gate, vec![exclusive]).unwrap();
        let mut b = bank_with_sources(cx, gate, vec![source_from_materials(materials, 2)]).unwrap()
            .with_material_string_bending_relaxation(selection(), None, 16, gate).unwrap();
        assert_eq!(a.system().material_memory().len(), 8);
        assert!(a.source().string_bending_relaxation().is_none());
        for (a, b) in a.system().components().iter().zip(b.system().components()) {
            assert_eq!(a.modes(), b.modes());
        }
        for (a, b) in a.system().material_relaxation_branches().iter()
            .zip(b.system().material_relaxation_branches())
        {
            assert_eq!(a.relaxation_time_s, b.relaxation_time_s);
            for (x, y) in a.projection.iter().zip(&b.projection) {
                let ka = a.stiffness * x * x;
                let kb = b.stiffness * y * y;
                assert!((ka - kb).abs() <= 1e-12 * ka.max(kb));
            }
        }
        let mut forces = vec![0.0; a.system().mode_count()];
        forces[0] = 1e-4;
        let drive = || FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        let mut peak = 0.0_f64;
        let mut error = 0.0_f64;
        for _ in 0..128 {
            let x = a.step(&forces, drive()).unwrap();
            let y = b.step(&forces, drive()).unwrap();
            peak = peak.max(x.observer_pressure_pa.abs());
            error = error.max((x.observer_pressure_pa - y.observer_pressure_pa).abs());
            for frame in [x, y] {
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
        assert!(peak > 1e-15);
        assert!(error <= 1e-9 * peak, "identical source spectra changed the physical trajectory");
        assert!(a.system().coupled_frame().unwrap().material_dissipation_j > 0.0);
    });
}

#[test]
fn independent_layer_count_exceeds_homogeneous_arm_limit_without_merging_times() {
    with_cx(|cx, gate| {
        let section = section((0..10).map(|i|
            material(40e9, 0.002 + 0.001 * i as f64, 200_000.0)).collect());
        let expected = section.branches().to_vec();
        let mut bank = bank_with_sources(cx, gate, vec![source(section, 1)]).unwrap();
        let branches = bank.system().material_relaxation_branches();
        assert_eq!(branches.len(), 10);
        let card = &bank.cards()[0];
        let wave = core::f64::consts::PI / card.length_m;
        for (arm, original) in branches.iter().zip(&expected) {
            assert_eq!(arm.relaxation_time_s, original.relaxation_time_s);
            let stiffness = original.relaxing_stiffness_n_m2 * wave.powi(4)
                / card.linear_density_kg_m;
            assert!((arm.stiffness / stiffness - 1.0).abs() < 1e-12);
        }
        assert!(render(&mut bank).iter().any(|p| p.abs() > 1e-15));
        assert!(bank.system().material_memory().iter().all(|z| z.is_finite()));
    });
}

#[test]
fn exclusive_bonded_zero_arms_preserve_elastic_motion_and_refuse_duplicate_ownership() {
    with_cx(|cx, gate| {
        let materials = [material(0.0, 0.002, 200_000.0), material(0.0, 0.004, 200_000.0)];
        let make = || source(section(materials.to_vec()), 2)
            .with_torsion("G", "eta_s", 1).unwrap()
            .with_transverse_polarizations(
                TransverseDirections::try_new([0.6, 0.8], [1.0, 0.0]).unwrap()).unwrap();
        let mut zero = bank_with_sources(cx, gate, vec![make()]).unwrap();
        let mut elastic = bank_with_sources(cx, gate,
            vec![source_from_materials(materials.clone(), 2)]).unwrap();
        assert!(zero.system().material_memory().is_empty());
        assert_eq!(render(&mut zero), render(&mut elastic));
        let cold = bank_with_sources(cx, gate, vec![make()]).unwrap();
        assert!(cold.with_material_string_bending_relaxation(selection(), None, 16, gate).is_err());
    });
}
