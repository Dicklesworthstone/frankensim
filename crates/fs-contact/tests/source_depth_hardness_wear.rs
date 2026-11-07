//! Original depth profiles drive spatial contact caps and separate body wear.
//! Synthetic equations/support checks, not measured-material qualification.

#[path = "source_point/common.rs"]
mod common;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_domain;
use fs_contact::interface_binding::{
    ADHESION_ENERGY_DIMS, ADHESION_ENERGY_PROPERTY, BoundNormalInterfaceState,
    bind_dry_interface_system_card, bind_isotropic_elastic_interface, bind_normal_interface_state,
};
use fs_contact::normal_patch::plastic_source::depth_wear::{
    DepthHardnessWearError, SourceDepthHardnessWear,
};
use fs_contact::normal_patch::plastic_source::{PlasticContactControls, PlasticContactEnvelope};
use fs_contact::normal_patch::{FiniteGapGrid, FiniteGapPlasticLinearSolver};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    HardnessLoadStep, HardnessTestContext, InterfaceSystemCard, InterpolationPolicy, MaterialCard,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::depth::{
    DepthHardnessError, EXPOSED_MATERIAL_DEPTH_AXIS as DEPTH, ResolvedDepthHardness,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec, Temperature};
use fs_tribo::InputAuthority;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const MAX_DEPTH: f64 = 2.0e-6;
const AREA: f64 = 1.0e-8;
const K: [f64; 2] = [0.01, 0.02];
const H0: f64 = 1.0e6;
const HB: f64 = 2.0e6;
const GRADIENT: f64 = -4.0e11;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 32,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic original depth-hardness/wear fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}

fn curve(surface: f64, gradient: f64) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: DEPTH.into(),
        abscissa_dims: LENGTH,
        knots: vec![(0.0, surface), (MAX_DEPTH, surface + gradient * MAX_DEPTH)],
        dims: Pressure::DIMS,
    }
}

fn material(name: &str, value: PropertyValue) -> (MaterialCard, PropertyKey) {
    let original = common::solid(name, 1000.0, 2.0e7);
    let mut claims = original.claims().clone();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("{name}/original material depths"),
            method: "synthetic force/projected-area indentation".into(),
            artifact: hash_domain("depth-hardness-wear-test", name.as_bytes()),
            caveats: "numerical fixture, not experimentally calibrated".into(),
            provenance: provenance(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("depth-H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "rigid sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic indentation pressure v1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value,
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0)
                .with_quantity(DEPTH, QuantitySpec::dimensional(LENGTH), 0.0, MAX_DEPTH),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    (
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap(),
        key,
    )
}

struct Fixture {
    cards: [MaterialCard; 2],
    profiles: [ResolvedDepthHardness; 2],
    normal: BoundNormalInterfaceState,
    coefficients: ResolvedInterfaceStatePoint,
    keys: [PropertyKey; 2],
}

impl Fixture {
    fn new(gradient: f64) -> Self {
        let sources = [
            material("A", curve(H0, gradient)),
            material("B", curve(HB, 0.0)),
        ];
        let cards = sources.each_ref().map(|(card, _)| card.clone());
        let point = common::query();
        let select = MaterialPropertySelection::SingleClaimOnly;
        let profiles = sources.each_ref().map(|(card, key)| {
            ResolvedDepthHardness::new(card, &point, key.clone(), select.clone(), MAX_DEPTH)
                .unwrap()
        });
        let base = common::bind(&cards[0], &cards[1], None, 0.0);
        let original = base.normal_state().elastic().interface().card();
        let mut claims = original.claims().clone();
        let keys = ["k_a", "k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
        for (key, value) in keys.clone().into_iter().zip(K) {
            claims
                .insert_claim(PropertyClaim {
                    key,
                    value: PropertyValue::Scalar {
                        value,
                        dims: Dims::NONE,
                    },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: provenance(),
                })
                .unwrap();
        }
        let interface = InterfaceSystemCard::assemble(
            original.surface_a().clone(),
            original.surface_b().clone(),
            original.context().clone(),
            claims,
            original.models().to_vec(),
        )
        .unwrap();
        let solids = cards
            .each_ref()
            .map(|card| resolve_isotropic_solid_state_point(card, &point, select.clone()).unwrap());
        let dry =
            bind_dry_interface_system_card(&interface, InputAuthority::SyntheticFixture).unwrap();
        let elastic = bind_isotropic_elastic_interface(&dry, &solids[0], &solids[1]).unwrap();
        let adhesion = resolve_interface_state_point(
            &interface,
            &point,
            &[ScalarPropertyRequirement::try_new(
                ADHESION_ENERGY_PROPERTY,
                ADHESION_ENERGY_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap()],
            select.clone(),
        )
        .unwrap();
        let normal = bind_normal_interface_state(&elastic, &adhesion).unwrap();
        let coefficients = resolve_interface_state_point(
            &interface,
            &point,
            &keys.each_ref().map(|key| {
                ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative)
                    .unwrap()
            }),
            select,
        )
        .unwrap();
        Self {
            cards,
            profiles,
            normal,
            coefficients,
            keys,
        }
    }

    fn owner(
        &self,
        two_cells: bool,
        solver: FiniteGapPlasticLinearSolver,
        maximum_plastic_depth_m: f64,
        cx: &Cx<'_>,
    ) -> SourceDepthHardnessWear {
        let mut gap = vec![1.0e-3; 25];
        gap[12] = 0.0;
        if two_cells {
            gap[13] = 0.0;
        }
        SourceDepthHardnessWear::new(
            self.normal.clone(),
            self.coefficients.clone(),
            self.keys.clone(),
            self.profiles.clone(),
            FiniteGapGrid {
                cells_x: 5,
                cells_y: 5,
                cell_width_m: 1.0e-4,
                cell_depth_m: 1.0e-4,
            },
            gap,
            PlasticContactControls {
                maximum_active_set_iterations: 512,
                complementarity_tolerance_m: 1.0e-13,
                boundary_clearance_cells: 1,
            },
            PlasticContactEnvelope {
                maximum_approach_m: 1.0e-4,
                maximum_plastic_depth_m,
                half_space_depth_m: [1.0e-2; 2],
                maximum_patch_to_depth: [0.1; 2],
            },
            [1.5e-6; 2],
            solver,
            cx,
        )
        .unwrap()
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-9 * actual.abs().max(expected.abs()).max(1.0e-30),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_local_exposed_hardness_changes_caps_wear_and_native_work() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(GRADIENT);
        let pristine = fixture.owner(true, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx);
        let loaded = pristine.advance_normal(5.0e-5, cx).unwrap();
        let before = loaded.candidate;
        let saved = before.clone();
        assert!(before.contact_state().plastic_depth_m()[12] > MAX_DEPTH);
        assert_eq!(before.volumes_m3(), [0.0; 2]);
        assert!(
            before
                .recession_m()
                .iter()
                .all(|depths| depths.iter().all(|z| *z == 0.0))
        );
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-5;
        slip[13] = 2.0e-5;
        let worn = before.advance_wear(&slip, cx).unwrap();
        assert_eq!(before, saved);
        for i in [12, 13] {
            close(before.contact_state().pressure_pa()[i], H0);
            for (side, h) in [H0, HB].into_iter().enumerate() {
                let expected = K[side] * H0 * slip[i] / h;
                close(worn.recession_increment_m[side][i], expected);
                close(worn.candidate.recession_m()[side][i], expected);
            }
            let depth = worn.candidate.recession_m()[0][i];
            close(
                worn.candidate.model().config().hardness_pa[i],
                H0 + GRADIENT * depth,
            );
            close(
                worn.candidate.contact_state().pressure_pa()[i],
                H0 + GRADIENT * depth,
            );
            let endpoint = fixture.profiles[0]
                .resolve_at_depth(&fixture.cards[0], depth)
                .unwrap();
            assert_eq!(
                endpoint
                    .resolved()
                    .query_point()
                    .iter()
                    .find(|(axis, _)| axis == DEPTH)
                    .map(|(_, value)| *value),
                Some(depth),
            );
            assert_eq!(endpoint.key(), fixture.profiles[0].reference_source().key());
            close(
                endpoint.hardness_pa(),
                worn.candidate.model().config().hardness_pa[i],
            );
        }
        for side in 0..2 {
            close(
                worn.candidate.volumes_m3()[side],
                AREA * worn.candidate.recession_m()[side].iter().sum::<f64>(),
            );
        }
        assert_eq!(worn.recession.normal_work_j, 0.0);
        assert_eq!(worn.hardness.normal_work_j, 0.0);
        assert!(worn.hardness.response.plastic_work_j > 0.0);
        close(
            worn.hardness.elastic_energy_j - loaded.receipt.normal.reversible_energy_j
                + worn.recession.response.plastic_work_j
                + worn.hardness.response.plastic_work_j,
            worn.recession.parameter_work_j + worn.hardness.parameter_work_j,
        );
        assert_eq!(
            worn.recession.predecessor_state_identity,
            before.contact_state().identity()
        );
        assert_eq!(
            worn.hardness.predecessor_state_identity,
            worn.recession.successor_state_identity
        );
        assert_eq!(
            worn.hardness.successor_state_identity,
            worn.candidate.contact_state().identity()
        );
        assert_eq!(worn.candidate.source_identity(), before.source_identity());
    });
}

#[test]
fn g1_refining_slip_matches_graded_hardness_wear_integral() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(GRADIENT);
        let initial = fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap()
            .candidate;
        let total_slip = 1.0e-4;
        // A remains pressure-saturated: dz_A/ds=k_A. The independent continuous
        // B integral is k_B/H_B integral(H0 + gradient*k_A*s) ds.
        let exact_b =
            K[1] / HB * (H0 * total_slip + 0.5 * GRADIENT * K[0] * total_slip * total_slip);
        let mut previous_error = None;
        for intervals in [1, 2, 4, 8] {
            let mut state = initial.clone();
            let mut slip = vec![0.0; 25];
            slip[12] = total_slip / f64::from(intervals);
            for _ in 0..intervals {
                state = state.advance_wear(&slip, cx).unwrap().candidate;
            }
            close(state.recession_m()[0][12], K[0] * total_slip);
            close(
                state.contact_state().pressure_pa()[12],
                H0 + GRADIENT * K[0] * total_slip,
            );
            close(state.volumes_m3()[1], AREA * state.recession_m()[1][12]);
            let error = state.recession_m()[1][12] - exact_b;
            assert!(error > 0.0);
            if let Some(previous) = previous_error {
                close(error * 2.0, previous);
            }
            previous_error = Some(error);
        }
    });
}

#[test]
fn g3_constant_profile_matches_original_native_recession_in_all_lanes() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(0.0);
        for solver in [
            FiniteGapPlasticLinearSolver::Dense,
            FiniteGapPlasticLinearSolver::MatrixFree {
                maximum_iterations: 256,
                maximum_corrections: 3,
            },
            FiniteGapPlasticLinearSolver::Fft {
                maximum_iterations: 256,
                maximum_corrections: 3,
            },
        ] {
            let state = fixture
                .owner(true, solver, 1.0e-4, cx)
                .advance_normal(5.0e-5, cx)
                .unwrap()
                .candidate;
            let mut slip = vec![0.0; 25];
            slip[12] = 1.0e-5;
            slip[13] = 2.0e-5;
            let relative: Vec<_> = state
                .contact_state()
                .pressure_pa()
                .iter()
                .zip(&slip)
                .map(|(&p, &ds)| (K[0] * (p / H0)) * ds + (K[1] * (p / HB)) * ds)
                .collect();
            let native = state
                .model()
                .recede(state.contact_state(), &relative, cx)
                .unwrap();
            let worn = state.advance_wear(&slip, cx).unwrap();
            assert_eq!(worn.candidate.contact_state(), &native.candidate);
            assert_eq!(worn.recession, native.receipt);
            assert_eq!(worn.hardness.parameter_work_j, 0.0);
            assert_eq!(worn.hardness.response.plastic_work_j, 0.0);
            close(
                worn.candidate.volumes_m3()[0],
                AREA * K[0] * (slip[12] + slip[13]),
            );
            close(
                worn.candidate.volumes_m3()[1],
                AREA * K[1] * H0 / HB * (slip[12] + slip[13]),
            );
            let unchanged = worn.candidate.advance_wear(&[0.0; 25], cx).unwrap();
            assert_eq!(unchanged.candidate, worn.candidate);
        }
    });
}

#[test]
fn g0_prepared_profile_requires_complete_source_support_and_exact_tests() {
    let fixture = Fixture::new(GRADIENT);
    let profile = &fixture.profiles[0];
    for depth in [0.0, 0.2e-6, 0.7e-6, MAX_DEPTH] {
        close(
            profile.hardness_pa_at(depth).unwrap(),
            H0 + GRADIENT * depth,
        );
        let source = profile.resolve_at_depth(&fixture.cards[0], depth).unwrap();
        assert_eq!(source.key(), profile.reference_source().key());
        assert_eq!(
            source.property().answer().receipt.selected,
            profile
                .reference_source()
                .property()
                .answer()
                .receipt
                .selected
        );
        close(source.hardness_pa(), profile.hardness_pa_at(depth).unwrap());
    }
    assert!(profile.resolve_at_depth(&fixture.cards[1], 0.0).is_err());
    for depth in [-1.0e-12, MAX_DEPTH * 1.01, f64::NAN, f64::INFINITY] {
        assert!(profile.hardness_pa_at(depth).is_err());
    }
    let mut samples = Vec::new();
    let mut triangles = Vec::new();
    for (lo, hi) in [(0.0, 0.5e-6), (1.5e-6, MAX_DEPTH)] {
        let base = samples.len() as u32;
        for (depth, t) in [(lo, 280.0), (hi, 280.0), (hi, 320.0), (lo, 320.0)] {
            samples.push(([depth, t], H0 + GRADIENT * depth));
        }
        triangles.extend([[base, base + 1, base + 2], [base, base + 2, base + 3]]);
    }
    let (gapped, key) = material(
        "gap",
        PropertyValue::Surface2D {
            axes: [DEPTH.into(), "T".into()],
            axis_dims: [LENGTH, Temperature::DIMS],
            samples,
            triangles,
            dims: Pressure::DIMS,
        },
    );
    assert!(matches!(
        ResolvedDepthHardness::new(
            &gapped,
            &common::query(),
            key,
            MaterialPropertySelection::SingleClaimOnly,
            MAX_DEPTH
        ),
        Err(DepthHardnessError::Support(_))
    ));
    let (negative, key) = material(
        "negative interior",
        PropertyValue::Curve {
            abscissa: DEPTH.into(),
            abscissa_dims: LENGTH,
            knots: vec![(0.0, H0), (1.0e-6, -1.0), (MAX_DEPTH, H0)],
            dims: Pressure::DIMS,
        },
    );
    assert!(
        ResolvedDepthHardness::new(
            &negative,
            &common::query(),
            key,
            MaterialPropertySelection::SingleClaimOnly,
            MAX_DEPTH
        )
        .is_err()
    );
}

#[test]
fn g4_depth_and_late_physical_refusals_keep_the_complete_owner() {
    with_cx(|cx, gate| {
        let fixture = Fixture::new(GRADIENT);
        let state = fixture
            .owner(true, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap()
            .candidate;
        let saved = state.clone();
        assert!(state.advance_wear(&[0.0; 24], cx).is_err());
        for bad in [-1.0, f64::NAN, f64::INFINITY, 1.0] {
            let mut slip = vec![0.0; 25];
            slip[12] = 1.0e-5;
            slip[13] = bad;
            assert!(state.advance_wear(&slip, cx).is_err());
            assert_eq!(state, saved);
        }
        // Initial plastic depth is below this original physical limit. Recession
        // alone passes; exposing a softer layer induces additional plastic flow
        // which fails only at final cap-refresh admission, after private work.
        let limited = fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 4.0e-5, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap()
            .candidate;
        let original = limited.clone();
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-4;
        let error = limited.advance_wear(&slip, cx).unwrap_err();
        assert!(
            matches!(error, DepthHardnessWearError::Contact(_)),
            "{error:?}"
        );
        assert_eq!(limited, original);
        gate.request();
        assert!(state.advance_wear(&[0.0; 25], cx).is_err());
        assert_eq!(state, saved);
    });
}
