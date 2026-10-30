//! Authored walls stay on the existing DG/plate owner. Synthetic sources are not measurements.
use super::super::{
    AperturePerformance, CoupledAperture, PlateValveInputError, PlateValvePerformance,
};
use crate::bernoulli_aperture::dg::ApertureDg;
use crate::pcm_wav::observation::PressureRenderer;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::{Density, Dims, DynViscosity, Pressure, Temperature};

const DG: &str = "dg_tube 0.05 0.001 open 2 2 8 50 2000 0.01 1e-11 128";
const WALL: &str =
    "dg_wall isotropic wall.matpack 300 0 0.00002 stress-free 20 100000 kelvin-voigt\n";
const REGIONAL: &str = concat!(
    "dg_wall regions 2\n",
    "dg_wall_region 0 0.025 isotropic wall.matpack 300 0 0.00002 stress-free 20 100000 kelvin-voigt\n",
    "dg_wall_region 0.025 0.05 isotropic wall.matpack 300 0 0.00004 strain-fixed 20 100000 kelvin-voigt\n",
);

fn input(wall: &str) -> String {
    let (header, _) = include_str!("../../../../../examples/plate-valve-dg.performance")
        .split_once("\nschedule\n")
        .unwrap();
    format!(
        "{}\nschedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tmouth\tblowing-pressure\ninitial\tpressure-pa\t3e1\nevents\t0\n",
        header
            .replace("audio 96000 9602 20", "audio 192000 64 20")
            .replace("observation inlet\n", &format!("{wall}observation inlet\n"))
    )
}

fn wall_pack(young: f64, omit: Option<&str>) -> NormalizedMaterialCardPack {
    let artifact = hash_domain("org.frankensim.test.authored-dg-wall", &young.to_le_bytes());
    let provenance = Provenance {
        source: "synthetic wall; not measured".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic wall".into(),
            method: "manufactured regression".into(),
            artifact,
            caveats: "not a measured material or instrument".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, dims, values) in [
        ("density", Density::DIMS, [4000.0; 3]),
        (
            "young_modulus",
            Pressure::DIMS,
            [1.1 * young, young, 0.5 * young],
        ),
        ("poisson_ratio", Dims::NONE, [0.3; 3]),
        (
            "kelvin_voigt_young_viscosity",
            DynViscosity::DIMS,
            [10.0; 3],
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [500.0; 3],
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            [5.0, 6.0, 11.0],
        ),
    ] {
        if omit == Some(name) {
            continue;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [200.0, 300.0, 800.0].into_iter().zip(values).collect(),
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 200.0, 800.0)
                    .with("omega", 0.0, 1_000_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        NormalizedPack::new(
            "synthetic-wall",
            "authored-v1",
            artifact,
            "CC0-1.0",
            claims,
            vec![],
            vec![],
        )
        .unwrap(),
    )
    .unwrap()
}

fn loaded(
    text: &str,
    material: &NormalizedMaterialCardPack,
    block: usize,
) -> Result<PlateValvePerformance, PlateValveInputError> {
    PlateValvePerformance::from_bytes_with_all_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "wall.matpack");
            Ok(material.clone())
        },
    )
}

fn dg(renderer: &AperturePerformance) -> &ApertureDg {
    let CoupledAperture::Dg(system) = renderer.system() else {
        panic!("expected original DG owner")
    };
    system
}

fn rendered(performance: PlateValvePerformance, block: usize) -> (Vec<f64>, AperturePerformance) {
    let mut pressure = vec![0.0; performance.info().samples as usize];
    let mut renderer = performance.into_renderer();
    for chunk in pressure.chunks_mut(block) {
        renderer.block(chunk).unwrap();
    }
    assert!(pressure.iter().all(|p| p.is_finite()));
    (pressure, renderer)
}

#[test]
fn omitted_and_explicit_rigid_walls_preserve_pressure_and_original_histories() {
    let material = wall_pack(1e8, None);
    let old = loaded(&input(""), &material, 13).unwrap();
    let rigid = loaded(&input("dg_wall rigid\n"), &material, 1).unwrap();
    assert!(old.dg_wall_layout().is_none());
    assert!(rigid.dg_wall_layout().unwrap().uniform_record().is_none());
    assert!(dg(rigid.renderer()).load().material_wall_source().is_none());
    let (a, left) = rendered(old, 13);
    let (b, right) = rendered(rigid, 1);
    assert_eq!(a, b);
    assert_eq!(dg(&left).duct_state(), dg(&right).duct_state());
    assert_eq!(
        left.system().aperture().state(),
        right.system().aperture().state()
    );
}

#[test]
fn authored_uniform_and_regional_sources_drive_native_wall_motion_and_pressure() {
    let material = wall_pack(1e8, None);
    for declaration in [WALL, REGIONAL] {
        let text = input(declaration);
        let p = loaded(&text, &material, 13).unwrap();
        let layout = p.dg_wall_layout().unwrap();
        let source = dg(p.renderer()).load().material_wall_source().unwrap();
        if declaration == WALL {
            let record = layout.uniform_record().unwrap();
            assert_eq!(record.path(), "wall.matpack");
            assert_eq!(record.pack(), &material);
            assert_eq!(source.uniform_wall().unwrap(), record.wall());
        } else {
            assert_eq!(layout.regions().len(), 2);
            assert_eq!(layout.regions()[0].interval_m(), [0.0, 0.025]);
            assert_eq!(layout.regions()[1].interval_m(), [0.025, 0.05]);
            let walls = source.regional_walls().unwrap();
            for (actual, record) in walls.iter().zip(layout.regions()) {
                assert_eq!(actual, record.source().wall());
                assert_eq!(record.source().pack(), &material);
            }
        }
        let duct = dg(p.renderer()).load().duct();
        let first = duct.compliant_wall_at_element(0).unwrap();
        assert!((first.surface_density / 0.08 - 1.0).abs() < 1e-14);
        assert!((first.stiffness_per_area / 2e9 - 1.0).abs() < 1e-14);
        if declaration == REGIONAL {
            let second = duct.compliant_wall_at_element(1).unwrap();
            assert!((second.surface_density / 0.16 - 1.0).abs() < 1e-14);
            assert!((second.stiffness_per_area / (4e9 / 0.91) - 1.0).abs() < 1e-14);
        }
        let identity = duct.operator_identity();
        let blade = p.section_provenance().to_vec();
        let (pressure, accepted) = rendered(p, 13);
        let wall = dg(&accepted)
            .load()
            .duct()
            .compliant_wall_state(dg(&accepted).duct_state())
            .unwrap();
        assert!(wall.stored_energy_j > 0.0);
        assert!(wall.displacement_m.iter().any(|x| *x != 0.0));
        assert_eq!(dg(&accepted).duct_state().step_count(), 64);
        assert_eq!(
            pressure,
            rendered(loaded(&text, &material, 1).unwrap(), 1).0
        );
        let changed = loaded(&text, &wall_pack(4e8, None), 13).unwrap();
        assert_eq!(changed.section_provenance(), blade.as_slice());
        assert_ne!(
            dg(changed.renderer()).load().duct().operator_identity(),
            identity
        );
        let (changed_pressure, changed) = rendered(changed, 13);
        assert!(
            pressure
                .iter()
                .zip(changed_pressure)
                .any(|(a, b)| (a - b).abs() > 1e-15)
        );
        assert_ne!(
            accepted.system().aperture().state(),
            changed.system().aperture().state()
        );
    }
}

#[test]
fn one_original_pack_is_loaded_once_for_plate_and_all_wall_regions() {
    let material = wall_pack(1e8, None);
    let text = input(REGIONAL).replace(
        "section isotropic 0.0003 900 4000000000 0.3 fixture authored-illustrative-blade-not-measured",
        "section material isotropic 0.0003 wall.matpack manufactured T=300 omega=0",
    ).replace(
        "relaxation 1\ntime_limits 0.1 0.5\nmemory_initial relaxed\nregion 0 4000000000 0.3 0 10000 2\nbranch 2000000000 0.001\nbranch 1000000000 0.01\n",
        "relaxation 0\n",
    );
    let mut loads = 0;
    let performance = PlateValvePerformance::from_bytes_with_all_sources(
        text.as_bytes(),
        13,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "wall.matpack");
            loads += 1;
            Ok(material.clone())
        },
    )
    .unwrap();
    assert_eq!(loads, 1);
    assert!(performance.section_provenance()[0].contains("wall.matpack"));
    for region in performance.dg_wall_layout().unwrap().regions() {
        assert_eq!(region.source().pack(), &material);
        assert_eq!(
            region.source().wall().material().card_identity(),
            material.card().content_hash()
        );
    }
    assert_eq!(dg(performance.renderer()).duct_state().step_count(), 0);
    assert_eq!(
        performance.renderer().system().aperture().accepted_steps(),
        0
    );
}

#[test]
fn incompatible_or_malformed_walls_refuse_before_any_physical_step() {
    let material = wall_pack(1e8, None);
    for text in [
        input(WALL).replace(DG, "tube 0.25 0.007 -0.8 0.002 1048576"),
        input(WALL).replace("wall.matpack", "../wall.matpack"),
        input("dg_wall unsupported\n"),
    ] {
        let mut loads = 0;
        assert!(
            PlateValvePerformance::from_bytes_with_all_sources(
                text.as_bytes(),
                13,
                &CancelGate::new_clock_free(),
                None,
                |_| {
                    loads += 1;
                    Ok(material.clone())
                },
            )
            .is_err()
        );
        assert_eq!(loads, 0);
    }
    for text in [
        input(WALL).replace("300 0 0.00002", "900 0 0.00002"),
        input(WALL).replace("20 100000", "20 1000"),
        input(REGIONAL).replace("0.025", "0.02"),
        input("dg_wall regions 0\n"),
    ] {
        assert!(loaded(&text, &material, 13).is_err());
    }
    assert!(loaded(&input(WALL), &wall_pack(1e8, Some("young_modulus")), 13).is_err());
    let gate = CancelGate::new_clock_free();
    gate.request();
    let mut loads = 0;
    assert!(
        PlateValvePerformance::from_bytes_with_all_sources(
            input(WALL).as_bytes(),
            13,
            &gate,
            None,
            |_| {
                loads += 1;
                Ok(material.clone())
            },
        )
        .is_err()
    );
    assert_eq!(loads, 0);
}
