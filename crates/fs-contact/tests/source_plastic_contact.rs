//! Source-preserving pressure saturation through real material-card resolvers.

#[path = "source_plastic_contact/temperature.rs"]
mod temperature;

#[path = "source_plastic_contact/load.rs"]
mod load;

#[path = "source_plastic_contact/drive.rs"]
mod drive;

#[path = "source_plastic_contact/linear.rs"]
mod linear;

#[path = "source_point/common.rs"]
mod common;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_domain;
use fs_contact::interface_binding::BoundNormalInterfaceState;
use fs_contact::normal_patch::FiniteGapGrid;
use fs_contact::normal_patch::plastic_source::{
    PlasticContactControls, PlasticContactEnvelope, SourcePlasticContact, SourcePlasticContactError,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::hardness::{INDENTATION_HARDNESS_QUANTITY, IndentationHardnessStatePoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Dims, Pressure, QtyAny};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn material(name: &str, hardness: f64, test: &str) -> (MaterialCard, PropertyKey) {
    let original = common::solid(name, 1000.0, 2.0e7);
    let mut claims = original.claims().clone();
    let provenance = Provenance {
        source: "synthetic pressure-hardness contact fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("{name}/{test}"),
            method: "synthetic spherical indentation; force/projected area".into(),
            artifact: hash_domain("test", test.as_bytes()),
            caveats: "numerical fixture; no experimental validation".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("indentation-hardness", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic rigid sphere, radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic force/projected-area protocol v1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Scalar {
                value: hardness,
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    (
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap(),
        key,
    )
}

fn sources(
    values: [f64; 2],
    test: &str,
    adhesion: f64,
) -> (
    BoundNormalInterfaceState,
    [IndentationHardnessStatePoint; 2],
) {
    let a = material("A", values[0], test);
    let b = material("B", values[1], test);
    let resolve = |(card, key): &(MaterialCard, PropertyKey)| {
        let state = resolve_material_state_point(
            card,
            &common::query(),
            &[
                ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
                    .unwrap(),
            ],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        IndentationHardnessStatePoint::try_from_resolved(state, key.clone()).unwrap()
    };
    let hardness = [resolve(&a), resolve(&b)];
    let normal = common::bind(&a.0, &b.0, None, adhesion)
        .normal_state()
        .clone();
    (normal, hardness)
}

fn envelope() -> PlasticContactEnvelope {
    PlasticContactEnvelope {
        maximum_approach_m: 1.0e-4,
        maximum_plastic_depth_m: 1.0e-4,
        half_space_depth_m: [1.0e-2; 2],
        maximum_patch_to_depth: [0.1; 2],
    }
}

fn build(
    normal: BoundNormalInterfaceState,
    hardness: [IndentationHardnessStatePoint; 2],
    envelope: PlasticContactEnvelope,
    cx: &Cx<'_>,
) -> Result<SourcePlasticContact, SourcePlasticContactError> {
    let mut gap = vec![1.0e-3; 25];
    gap[12] = 0.0;
    SourcePlasticContact::new(
        normal,
        hardness,
        FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 1.0e-4,
            cell_depth_m: 1.0e-4,
        },
        gap,
        PlasticContactControls {
            maximum_active_set_iterations: 100,
            complementarity_tolerance_m: 1.0e-12,
            boundary_clearance_cells: 1,
        },
        envelope,
        cx,
    )
}

#[test]
fn original_softer_surface_controls_pressure_and_retained_dent() {
    with_cx(|cx, _| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "base", 0.0);
        let law = build(normal.clone(), hardness.clone(), envelope(), cx).unwrap();
        assert_eq!(law.normal_source(), &normal);
        assert_eq!(law.hardness_sources(), &hardness);

        let loaded = law.advance(&law.initial_state(), 5.0e-5, cx).unwrap();
        common::close(loaded.receipt.normal.peak_pressure_pa, 1.0e6);
        common::close(loaded.receipt.normal.normal_force_n, 1.0e6 * 1.0e-8);
        let dent = loaded.candidate.contact_state().plastic_depth_m()[12];
        assert!(dent > 0.0);
        common::close(loaded.receipt.plastic_work_j, 1.0e6 * 1.0e-8 * dent);
        assert_eq!(loaded.source_identity, law.identity());

        let unloaded = law.advance(&loaded.candidate, -1.0e-5, cx).unwrap();
        assert_eq!(unloaded.receipt.normal.normal_force_n, 0.0);
        assert_eq!(unloaded.receipt.plastic_work_j, 0.0);
        assert_eq!(
            unloaded.candidate.contact_state().plastic_depth_m(),
            loaded.candidate.contact_state().plastic_depth_m()
        );
        let reloaded = law.advance(&unloaded.candidate, 5.0e-5, cx).unwrap();
        common::close(
            reloaded.receipt.normal.normal_force_n,
            loaded.receipt.normal.normal_force_n,
        );
        assert!(reloaded.receipt.plastic_work_j <= 1.0e-16);

        let (normal, hardness) = sources([1.0e8; 2], "hard", 0.0);
        let hard = build(normal, hardness, envelope(), cx).unwrap();
        let elastic = hard.advance(&hard.initial_state(), 5.0e-5, cx).unwrap();
        assert!(elastic.receipt.normal.normal_force_n > loaded.receipt.normal.normal_force_n);
        assert_eq!(elastic.receipt.plastic_work_j, 0.0);
        assert_eq!(
            elastic.candidate.contact_state().plastic_depth_m(),
            &[0.0; 25]
        );
    });
}

#[test]
fn equal_numbers_do_not_erase_original_source_or_surface_roles() {
    with_cx(|cx, _| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "base", 0.0);
        let a = build(normal.clone(), hardness.clone(), envelope(), cx).unwrap();
        let (other_normal, other_hardness) = sources([1.0e6, 2.0e6], "different-test", 0.0);
        let b = build(other_normal, other_hardness.clone(), envelope(), cx).unwrap();
        assert_eq!(a.model().identity(), b.model().identity());
        assert_ne!(a.identity(), b.identity());
        assert!(b.advance(&a.initial_state(), 1.0e-6, cx).is_err());
        assert!(build(normal.clone(), other_hardness, envelope(), cx).is_err());
        assert!(
            build(
                normal.clone(),
                [hardness[1].clone(), hardness[0].clone()],
                envelope(),
                cx,
            )
            .is_err()
        );

        let warmer = hardness[0]
            .at_query_point(
                &material("A", 1.0e6, "base").0,
                &QueryPoint::new().with("T", 310.0).unwrap(),
            )
            .unwrap();
        assert!(build(normal, [warmer, hardness[1].clone()], envelope(), cx).is_err());
    });
}

#[test]
fn unsupported_adhesion_and_late_depth_or_domain_failure_do_not_mutate_history() {
    with_cx(|cx, gate| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "base", 0.1);
        assert!(build(normal, hardness, envelope(), cx).is_err());

        let (normal, hardness) = sources([1.0e6, 2.0e6], "base", 0.0);
        let mut bounded = envelope();
        bounded.maximum_plastic_depth_m = 1.0e-6;
        let law = build(normal.clone(), hardness.clone(), bounded, cx).unwrap();
        let before = law.initial_state();
        let saved = before.clone();
        assert!(matches!(
            law.advance(&before, 5.0e-5, cx),
            Err(SourcePlasticContactError::OutsideApplicability {
                quantity: "permanent relative indentation [m]",
                ..
            })
        ));
        assert_eq!(before, saved);
        assert!(law.advance(&before, 1.0e-6, cx).is_ok());
        assert!(law.advance(&before, 1.01e-4, cx).is_err());

        let mut shallow = envelope();
        shallow.half_space_depth_m = [1.0e-5; 2];
        let shallow_law = build(normal.clone(), hardness.clone(), shallow, cx).unwrap();
        assert!(matches!(
            shallow_law.advance(&shallow_law.initial_state(), 1.0e-6, cx),
            Err(SourcePlasticContactError::OutsideApplicability {
                quantity: "patch radius / homogeneous depth",
                ..
            })
        ));

        let law = build(normal, hardness, envelope(), cx).unwrap();
        let before = law.initial_state();
        let saved = before.clone();
        gate.request();
        assert!(law.advance(&before, 1.0e-6, cx).is_err());
        assert_eq!(before, saved);
    });
}
