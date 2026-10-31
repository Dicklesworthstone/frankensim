//! G1/G3: one material-bearing model drives mass, gravity, modes and histories.

use fs_frame::history::StoryParams;
use fs_frame::{BuildingFrame, BuildingSpec, RcMember, building_fragility_from_frame};
use fs_material::fiber::{ManderConcrete, MenegottoPintoSteel};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};
use fs_solid::fiber::FiberLaw;
use fs_solid::fiber::material::{FiberMaterial, RcSection, RcSectionGeometry};
use fs_solid::frame2d::FrameMassFormulation;

fn section(member: RcMember, density_scale: f64, stiffness_scale: f64) -> RcSection {
    let core = FiberMaterial::authored(
        FiberLaw::Concrete(
            ManderConcrete::new(42e6 * stiffness_scale, 0.004, 32e9 * stiffness_scale, 0.015)
                .unwrap(),
        ),
        2400.0 * density_scale,
    )
    .unwrap();
    let cover = FiberMaterial::authored(
        FiberLaw::Concrete(
            ManderConcrete::new(30e6 * stiffness_scale, 0.002, 30e9 * stiffness_scale, 0.005)
                .unwrap(),
        ),
        2400.0 * density_scale,
    )
    .unwrap();
    let steel = FiberMaterial::authored(
        FiberLaw::Steel(
            MenegottoPintoSteel::new(200e9 * stiffness_scale, 450e6 * stiffness_scale, 0.01)
                .unwrap(),
        ),
        7800.0 * density_scale,
    )
    .unwrap();
    RcSection::new(
        RcSectionGeometry {
            depth_m: member.depth,
            width_m: member.width,
            cover_m: 0.04,
            concrete_layers: 16,
            steel_area_per_face_m2: member.steel_area,
            steel_centroid_y_m: 0.5 * member.depth - 0.06,
        },
        core,
        cover,
        steel,
    )
    .unwrap()
}

fn spec() -> BuildingSpec {
    BuildingSpec {
        stories: 1,
        bays: 1,
        floor_mass: 20_000.0,
        sections_per_member: 3,
        ..BuildingSpec::default()
    }
}

fn building(spec: BuildingSpec, density: f64, stiffness: f64) -> BuildingFrame {
    BuildingFrame::with_sections(
        spec,
        section(spec.column, density, stiffness),
        section(spec.beam, density, stiffness),
        FrameMassFormulation::ConsistentReference,
    )
    .unwrap()
}

fn member_mass(member: RcMember) -> f64 {
    // Reinforcement replaces its volume of concrete; it is not added to a
    // gross concrete section a second time.
    2400.0 * (member.depth * member.width - 2.0 * member.steel_area)
        + 7800.0 * 2.0 * member.steel_area
}

#[test]
fn section_materials_supply_member_mass_and_column_scaling_exactly_once() {
    let spec = BuildingSpec {
        scale: 1.4,
        ..spec()
    };
    let model = building(spec, 1.0, 1.0);
    let sections = model.material_sections().unwrap();
    assert_eq!(sections[0].area_scale(), spec.scale);
    assert_eq!(sections[1].area_scale(), 1.0);
    assert!((sections[0].line_mass_kg_m() - member_mass(spec.column) * spec.scale).abs() < 1e-10);
    assert!((sections[1].line_mass_kg_m() - member_mass(spec.beam)).abs() < 1e-10);
    assert!(
        BuildingFrame::new(spec)
            .unwrap()
            .material_sections()
            .is_none()
    );
    let rest = model.clone().run(&[0.0; 21], 0.01).unwrap();
    assert!(rest.peak_drift_ratio < 1e-10, "{rest:?}");
    assert!(rest.internal_work.abs() < 1e-4, "{rest:?}");
}

#[test]
fn independent_density_and_constitutive_changes_reach_the_same_building() {
    let spec = spec();
    let baseline = building(spec, 1.0, 1.0);
    let dense = building(spec, 2.0, 1.0);
    let stiff = building(spec, 1.0, 1.5);
    let m0 = baseline.modal_analysis().base_influence_mass_kg[0];
    let md = dense.modal_analysis().base_influence_mass_kg[0];
    let ms = stiff.modal_analysis().base_influence_mass_kg[0];
    assert!(md > m0);
    assert_eq!(m0, ms);
    assert!(stiff.periods()[0] < baseline.periods()[0]);
    assert_ne!(dense.periods(), baseline.periods());
    let ground: Vec<f64> = (0..=50)
        .map(|i| 0.3 * (4.0 * std::f64::consts::PI * i as f64 * 0.01).sin())
        .collect();
    let a = baseline.clone().run(&ground, 0.01).unwrap();
    let b = dense.clone().run(&ground, 0.01).unwrap();
    let c = stiff.clone().run(&ground, 0.01).unwrap();
    for response in [&a, &b, &c] {
        assert!(response.peak_drift_ratio > 0.0 && response.peak_drift_ratio.is_finite());
        assert!(response.energy_balance_error < 1e-5, "{response:?}");
    }
    assert_ne!(a.roof_drift, b.roof_drift);
    assert_ne!(a.roof_drift, c.roof_drift);
    assert_eq!(a, baseline.clone().run(&ground, 0.01).unwrap());
}

#[test]
fn conflicting_geometry_and_repeated_scaling_refuse_before_assembly() {
    let spec = spec();
    let column = section(spec.column, 1.0, 1.0);
    let beam = section(spec.beam, 1.0, 1.0);
    for wrong_column in [beam.clone(), column.scaled(2.0).unwrap()] {
        assert!(
            BuildingFrame::with_sections(
                spec,
                wrong_column,
                beam.clone(),
                FrameMassFormulation::Lumped,
            )
            .is_err()
        );
    }
    assert!(
        BuildingFrame::with_sections(
            BuildingSpec {
                scale: f64::INFINITY,
                ..spec
            },
            column,
            beam,
            FrameMassFormulation::Lumped,
        )
        .is_err()
    );
}

#[test]
fn fragility_preserves_the_selected_fine_model_and_its_starting_state() {
    let model = building(spec(), 1.5, 0.8);
    let ensemble = StochasticEnsemble {
        name: "material-building".into(),
        seed: 29,
        members: 8,
        duration: QtyAny::new(0.2, Dims([0, 0, 1, 0, 0, 0])),
        dt: QtyAny::new(0.02, Dims([0, 0, 1, 0, 0, 0])),
        model: SpectrumModel::KanaiTajimi {
            s0: 1e-4,
            omega_g: QtyAny::new(12.5, Dims([0, 0, -1, 0, 0, 0])),
            zeta_g: 0.6,
        },
    };
    let motion = ensemble.realize(0).unwrap();
    let before = model.clone().run(&motion.values, 0.02).unwrap();
    assert!(before.peak_drift_ratio > 0.0);
    let report =
        building_fragility_from_frame(&ensemble, &model, StoryParams::default(), 0.01, 0.05, 0.3)
            .unwrap();
    assert_eq!(report.members_used, 8);
    assert_eq!(report.peak_drifts[0], before.peak_drift_ratio);
    assert_eq!(before, model.clone().run(&motion.values, 0.02).unwrap());
    assert_eq!(model.material_sections().unwrap()[0].area_scale(), 1.0);
}

#[test]
fn reference_mass_counts_every_floor_and_full_member_including_supports() {
    let spec = BuildingSpec {
        stories: 2,
        bays: 2,
        scale: 1.4,
        ..spec()
    };
    let expected = spec.stories as f64
        * (spec.floor_mass
            + (spec.bays + 1) as f64
                * spec.story_height
                * member_mass(spec.column)
                * spec.scale
            + spec.bays as f64 * spec.bay_width * member_mass(spec.beam));
    for formulation in [
        FrameMassFormulation::Lumped,
        FrameMassFormulation::ConsistentReference,
    ] {
        let model = BuildingFrame::with_sections(
            spec,
            section(spec.column, 1.0, 1.0),
            section(spec.beam, 1.0, 1.0),
            formulation,
        )
        .unwrap();
        assert!(
            (model.reference_mass_kg() - expected).abs() < 2e-12 * expected,
            "{formulation:?}: {} != {expected}",
            model.reference_mass_kg(),
        );
        if formulation == FrameMassFormulation::Lumped {
            // Supported first-story column mass belongs in physical inventory
            // even though no free mode carries its horizontal kinetic energy.
            let fixed_mass = 0.5
                * (spec.bays + 1) as f64
                * spec.story_height
                * member_mass(spec.column)
                * spec.scale;
            let effective_mass: f64 = model
                .modal_analysis()
                .modes
                .iter()
                .map(|mode| mode.effective_mass_kg[0])
                .sum();
            assert!(
                (model.reference_mass_kg() - effective_mass - fixed_mass).abs()
                    < 1e-8 * expected
            );
        }
    }
}

#[test]
fn reference_mass_uses_actual_authored_line_mass_without_rescaling_it() {
    let spec = BuildingSpec {
        scale: 1.7,
        ..spec()
    };
    let column_kg_m = 131.0;
    let beam_kg_m = 177.0;
    let model = BuildingFrame::with_member_mass(
        spec,
        column_kg_m,
        beam_kg_m,
        FrameMassFormulation::ConsistentReference,
    )
    .unwrap();
    let expected = spec.stories as f64
        * (spec.floor_mass
            + (spec.bays + 1) as f64 * spec.story_height * column_kg_m
            + spec.bays as f64 * spec.bay_width * beam_kg_m);
    assert!((model.reference_mass_kg() - expected).abs() < 2e-12 * expected);
    let floor_only = BuildingFrame::new(spec).unwrap();
    assert!(
        (floor_only.reference_mass_kg() - spec.stories as f64 * spec.floor_mass).abs()
            < 2e-12 * expected
    );
}
