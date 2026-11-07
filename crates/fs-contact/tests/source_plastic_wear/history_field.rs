//! Real accepted native histories drive spatial caps, separate wear and refresh.
//! Synthetic material curves exercise equations; no empirical calibration claim.

use std::sync::Arc;

use super::*;
use fs_contact::interface_binding::BoundNormalInterfaceState;
use fs_contact::normal_patch::plastic_source::history_field::{
    ContactHardnessField, HistoryHardnessContactError, HistoryHardnessControl as Control,
    SourceHistoryHardnessWear,
};
use fs_contact::normal_patch::{
    FiniteGapPlasticDrive, FiniteGapPlasticLinearSolver, FiniteGapPlasticLoadControls,
};
use fs_matdb::{ClaimSet, MaterialStateId};
use fs_material::hardness::history_field::AcceptedHardnessField;
use fs_material::hardness::plastic_history::{
    ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA, ResolvedJ2Hardness,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY, resolve_j2_plasticity,
};
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_qty::{Density, Pressure, QuantitySpec};

const AREA: f64 = 1.0e-8;
const K: [f64; 2] = [0.01, 0.02];
const GATES: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1.0e-6,
    yield_pa: 1.0e-7,
};

fn native_card(name: &str, h0: f64, h1: f64) -> (MaterialCard, PropertyKey, PropertyKey) {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", 2.0e7, Pressure::DIMS),
        ("poisson_ratio", 0.25, Dims::NONE),
        ("yield_stress", 1.0e4, Pressure::DIMS),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, 1.0e5, Pressure::DIMS),
    ] {
        claims
            .insert_claim(claim(PropertyKey::new(name, dims), value))
            .unwrap();
    }
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("{name}/native preconditioned history"),
            method: "synthetic force/projected-area indentation".into(),
            artifact: hash_domain("history-hardness-contact-test", name.as_bytes()),
            caveats: "numerical regression only".into(),
            provenance: provenance(),
        })
        .unwrap();
    let context = HardnessTestContext::new(
        "synthetic rigid sphere radius 1 mm",
        vec![
            HardnessLoadStep::new(
                QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
            )
            .unwrap(),
        ],
        "synthetic projected area pressure",
        observation,
    )
    .unwrap();
    let history_key = PropertyKey::with_quantity("history-H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(context.clone())
        .unwrap();
    let fixed_key = PropertyKey::with_quantity("fixed-H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(context)
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: history_key.clone(),
            value: PropertyValue::Curve {
                abscissa: ALPHA.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, h0), (0.04, h1)],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0)
                .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.04),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    claims.insert_claim(claim(fixed_key.clone(), h0)).unwrap();
    (
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "solid".into(),
                process: "synthetic preconditioned history".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        history_key,
        fixed_key,
    )
}

struct Fixture {
    cards: [MaterialCard; 2],
    sources: [ResolvedJ2Hardness; 2],
    fixed: [ContactHardnessField; 2],
    normal: BoundNormalInterfaceState,
    coefficients: ResolvedInterfaceStatePoint,
    keys: [PropertyKey; 2],
}
impl Fixture {
    fn new(softening: bool) -> Self {
        let data = [
            native_card("A", 1.0e6, if softening { 0.4e6 } else { 3.0e6 }),
            native_card("B", 2.0e6, 4.0e6),
        ];
        let point = common::query();
        let select = MaterialPropertySelection::SingleClaimOnly;
        let sources = data.each_ref().map(|(card, key, _)| {
            let plastic = resolve_j2_plasticity(
                card,
                &point,
                select.clone(),
                HardeningPropertyBasis::PlasticStrainModulus,
                0.05,
            )
            .unwrap();
            ResolvedJ2Hardness::new(plastic, card, key.clone(), select.clone(), 0.04).unwrap()
        });
        let fixed = data.each_ref().map(|(card, _, key)| {
            let resolved = resolve_material_state_point(
                card,
                &point,
                &[ScalarPropertyRequirement::try_with_key(
                    key,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap()],
                select.clone(),
            )
            .unwrap();
            ContactHardnessField::Fixed(
                IndentationHardnessStatePoint::try_from_resolved(resolved, key.clone()).unwrap(),
            )
        });
        let cards = data.map(|(card, _, _)| card);
        let base = common::bind(&cards[0], &cards[1], None, 0.0);
        let original = base.normal_state().elastic().interface().card();
        let mut claims = original.claims().clone();
        let keys = ["k_a", "k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
        for (key, value) in keys.clone().into_iter().zip(K) {
            claims.insert_claim(claim(key, value)).unwrap();
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
            sources,
            fixed,
            normal,
            coefficients,
            keys,
        }
    }

    fn field(&self, side: usize, axial: &[f64]) -> ContactHardnessField {
        let source = &self.sources[side];
        let accepted: Vec<_> = axial
            .iter()
            .map(|value| {
                source
                    .advance_candidate(
                        &source.initial_state(),
                        &[*value, 0.0, 0.0, 0.0, 0.0, 0.0],
                        GATES,
                    )
                    .unwrap()
                    .state
            })
            .collect();
        ContactHardnessField::Accepted(Arc::new(
            AcceptedHardnessField::from_j2(
                source,
                accepted
                    .iter()
                    .map(|state| (source.plastic(), state.strain(), state.plastic())),
                GATES.yield_pa,
            )
            .unwrap(),
        ))
    }

    fn owner(
        &self,
        hardness: [ContactHardnessField; 2],
        mapping: Vec<[usize; 2]>,
        two: bool,
        cx: &Cx<'_>,
    ) -> Result<SourceHistoryHardnessWear, HistoryHardnessContactError> {
        let mut gap = vec![1.0e-3; 25];
        gap[12] = 0.0;
        if two {
            gap[13] = 0.0;
        }
        SourceHistoryHardnessWear::new(
            self.normal.clone(),
            self.coefficients.clone(),
            self.keys.clone(),
            hardness,
            mapping,
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
                maximum_plastic_depth_m: 1.0e-4,
                half_space_depth_m: [0.01; 2],
                maximum_patch_to_depth: [0.1; 2],
            },
            [1.0e-4; 2],
            FiniteGapPlasticLinearSolver::Dense,
            cx,
        )
    }
}

fn load_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}
fn local_slip(distance: f64) -> Vec<f64> {
    let mut slip = vec![0.0; 25];
    slip[12] = distance;
    slip
}
fn compliance(owner: &SourceHistoryHardnessWear) -> f64 {
    4.0 * owner.model().config().grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (std::f64::consts::PI * owner.model().config().reduced_modulus_pa)
}

#[test]
fn g1_explicit_native_point_map_drives_caps_and_both_body_wear_denominators() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(false);
        let fields = [
            fixture.field(0, &[0.0, 0.02]),
            fixture.field(1, &[0.0, 0.01]),
        ];
        let mut map = vec![[0, 0]; 25];
        map[12] = [0, 1];
        map[13] = [1, 0];
        let owner = fixture
            .owner(fields.clone(), map.clone(), true, cx)
            .unwrap();
        let before = owner
            .advance_normal(Control::Approach(5.0e-5), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let mut slip = local_slip(1.0e-5);
        slip[13] = 2.0e-5;
        slip[0] = 1.0;
        let step = before
            .advance_wear(&slip, Control::Approach(5.0e-5), cx)
            .unwrap();
        for cell in [12, 13] {
            let h = std::array::from_fn::<_, 2, _>(|side| {
                fields[side].hardness_pa(map[cell][side]).unwrap()
            });
            close(before.contact_state().pressure_pa()[cell], h[0].min(h[1]));
            close(before.model().config().hardness_pa[cell], h[0].min(h[1]));
            for side in 0..2 {
                close(
                    step.recession_increment_m[side][cell],
                    K[side] * before.contact_state().pressure_pa()[cell] * slip[cell] / h[side],
                );
            }
        }
        for side in 0..2 {
            close(
                step.wear[side].increment_m3,
                AREA * (step.recession_increment_m[side][12]
                    + step.recession_increment_m[side][13]),
            );
            close(
                step.candidate.volumes_m3()[side],
                step.wear[side].increment_m3,
            );
            assert_eq!(step.candidate.recession_m()[side][0], 0.0);
        }
        assert_eq!(step.candidate.hardness_fields(), &fields);
        assert_eq!(step.candidate.cell_to_material_point(), map);
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_external_softening_refresh_retains_wear_dent_and_native_parameter_work() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(true);
        let before = fixture
            .owner(
                [fixture.field(0, &[0.0]), fixture.fixed[1].clone()],
                vec![[0, 0]; 25],
                false,
                cx,
            )
            .unwrap()
            .advance_normal(Control::Approach(5.0e-5), cx)
            .unwrap()
            .candidate
            .advance_wear(&local_slip(1.0e-6), Control::Approach(5.0e-5), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let next_fields = [fixture.field(0, &[0.03]), fixture.fixed[1].clone()];
        let refreshed = before
            .refresh_hardness(next_fields.clone(), Control::Approach(5.0e-5), cx)
            .unwrap();
        let next = &refreshed.candidate;
        assert!(next.contact_state().pressure_pa()[12] < before.contact_state().pressure_pa()[12]);
        assert!(
            next.contact_state().plastic_depth_m()[12]
                > before.contact_state().plastic_depth_m()[12]
        );
        assert!(refreshed.receipt.response.plastic_work_j > 0.0);
        close(
            refreshed.receipt.parameter_work_j,
            refreshed.receipt.elastic_energy_change_j + refreshed.receipt.response.plastic_work_j
                - refreshed.receipt.normal_work_j,
        );
        assert_eq!(next.volumes_m3(), before.volumes_m3());
        assert_eq!(next.recession_m(), before.recession_m());
        assert_eq!(next.source_identity(), before.source_identity());
        assert_eq!(next.hardness_fields(), &next_fields);
        assert_eq!(
            refreshed.receipt.predecessor_state_identity,
            before.contact_state().identity()
        );
        assert_eq!(
            refreshed.receipt.successor_state_identity,
            next.contact_state().identity()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_force_and_finite_spring_wear_and_field_refresh_keep_original_boundary() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(false);
        let pristine = fixture
            .owner(
                [fixture.field(0, &[0.0]), fixture.fixed[1].clone()],
                vec![[0, 0]; 25],
                false,
                cx,
            )
            .unwrap();
        let c = compliance(&pristine);
        let force = 0.004;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: AREA / c,
            position_m: 2.0 * c * force / AREA,
        };
        for control in [
            Control::Load {
                force_n: force,
                controls: load_controls(),
            },
            Control::Drive {
                drive,
                controls: load_controls(),
            },
        ] {
            let before = pristine.advance_normal(control, cx).unwrap().candidate;
            let p = before.contact_state().pressure_pa()[12];
            let step = before
                .advance_wear(&local_slip(1.0e-4), control, cx)
                .unwrap();
            let relative = K[0] * p * 1.0e-4 / 1.0e6 + K[1] * p * 1.0e-4 / 2.0e6;
            let expected = match control {
                Control::Load { .. } => p,
                _ => p - relative / (2.0 * c),
            };
            close(step.candidate.contact_state().pressure_pa()[12], expected);
            let before_refresh = step.candidate;
            let refreshed = before_refresh
                .refresh_hardness(
                    [fixture.field(0, &[0.02]), fixture.fixed[1].clone()],
                    control,
                    cx,
                )
                .unwrap();
            close(
                refreshed.candidate.contact_state().pressure_pa()[12],
                expected,
            );
            assert_eq!(
                refreshed.candidate.volumes_m3(),
                before_refresh.volumes_m3()
            );
            assert_eq!(
                refreshed.candidate.recession_m(),
                before_refresh.recession_m()
            );
            if let Control::Load { .. } = control {
                close(step.receipt.normal_work_j, force * relative);
            } else {
                let a = before.contact_state().approach_m();
                let b = before_refresh.contact_state().approach_m();
                close(
                    step.receipt.normal_work_j,
                    0.5 * drive.stiffness_n_per_m
                        * ((drive.position_m - a).powi(2) - (drive.position_m - b).powi(2)),
                );
            }
        }
    });
}

#[test]
fn g0_field_binding_and_refresh_refuse_foreign_sources_mapping_and_alpha_rewind_atomically() {
    with_cx(|cx, gate| {
        let fixture = Fixture::new(false);
        let fields = [fixture.field(0, &[0.01]), fixture.fixed[1].clone()];
        let mut invalid_map = vec![[0, 0]; 25];
        invalid_map[24][0] = 1;
        assert!(
            fixture
                .owner(fields.clone(), invalid_map, false, cx)
                .is_err()
        );
        assert!(
            fixture
                .owner(
                    [fixture.field(1, &[0.01]), fixture.fixed[1].clone()],
                    vec![[0, 0]; 25],
                    false,
                    cx
                )
                .is_err()
        );
        let warm = resolve_j2_plasticity(
            &fixture.cards[0],
            &QueryPoint::new().with("T", 294.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            HardeningPropertyBasis::PlasticStrainModulus,
            0.05,
        )
        .unwrap();
        let warm_source = ResolvedJ2Hardness::new(
            warm,
            &fixture.cards[0],
            fixture.sources[0].hardness_source().key().clone(),
            MaterialPropertySelection::SingleClaimOnly,
            0.04,
        )
        .unwrap();
        let warm_state = warm_source.initial_state();
        let warm_field = ContactHardnessField::Accepted(Arc::new(
            AcceptedHardnessField::from_j2(
                &warm_source,
                [(
                    warm_source.plastic(),
                    warm_state.strain(),
                    warm_state.plastic(),
                )]
                .into_iter(),
                GATES.yield_pa,
            )
            .unwrap(),
        ));
        assert!(
            fixture
                .owner(
                    [warm_field, fixture.fixed[1].clone()],
                    vec![[0, 0]; 25],
                    false,
                    cx
                )
                .is_err()
        );
        let before = fixture
            .owner(fields.clone(), vec![[0, 0]; 25], false, cx)
            .unwrap()
            .advance_normal(Control::Approach(2.0e-6), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        for bad in [
            [fixture.field(0, &[0.0]), fixture.fixed[1].clone()],
            [fixture.field(0, &[0.01, 0.02]), fixture.fixed[1].clone()],
            [fixture.fixed[0].clone(), fixture.fixed[1].clone()],
            [fixture.field(1, &[0.02]), fixture.fixed[1].clone()],
        ] {
            assert!(
                before
                    .refresh_hardness(bad, Control::Approach(2.0e-6), cx)
                    .is_err()
            );
        }
        assert!(
            before
                .advance_wear(&local_slip(1.0e-4), Control::Approach(3.0e-6), cx)
                .is_err()
        );
        let mut bad_slip = local_slip(1.0e-4);
        bad_slip[24] = f64::NAN;
        assert!(
            before
                .advance_wear(&bad_slip, Control::Approach(2.0e-6), cx)
                .is_err()
        );
        assert!(
            before
                .advance_wear(&local_slip(1.0), Control::Approach(2.0e-6), cx)
                .is_err()
        );
        gate.request();
        assert!(
            before
                .refresh_hardness(fields, Control::Approach(2.0e-6), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
