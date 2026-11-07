//! Independent H(T,z), elastic compliance, wear and boundary-work references.
//! Synthetic equations exercise original source adapters, not measured data.

use super::*;
use fs_contact::normal_patch::plastic_source::depth_wear::temperature::DepthHardnessTemperatureError;
use fs_contact::normal_patch::{FiniteGapPlasticDrive, FiniteGapPlasticLoadControls};
use fs_matdb::ClaimSet;

pub(super) const COLD: f64 = 293.15;
pub(super) const HOT: f64 = 313.15;
const DEPTH_SLOPE: [f64; 2] = [-1.0e11, -2.0e11];
const THERMAL_SLOPE: f64 = -2.5e4;

fn h(side: usize, temperature: f64, depth: f64) -> f64 {
    [H0, HB][side] + THERMAL_SLOPE * (temperature - COLD) + DEPTH_SLOPE[side] * depth
}

fn thermal_material(side: usize, support_temperature: f64) -> (MaterialCard, PropertyKey) {
    let surface = PropertyValue::Surface2D {
        axes: [DEPTH.into(), "T".into()],
        axis_dims: [LENGTH, Temperature::DIMS],
        samples: vec![
            ([0.0, COLD], h(side, COLD, 0.0)),
            ([MAX_DEPTH, COLD], h(side, COLD, MAX_DEPTH)),
            (
                [MAX_DEPTH, support_temperature],
                h(side, support_temperature, MAX_DEPTH),
            ),
            (
                [0.0, support_temperature],
                h(side, support_temperature, 0.0),
            ),
        ],
        triangles: vec![[0, 1, 2], [0, 2, 3]],
        dims: Pressure::DIMS,
    };
    let (original, key) = material(["A", "B"][side], surface);
    let mut claims = ClaimSet::new();
    for id in original.claims().observation_ids() {
        claims
            .register_observation(original.claims().observation(id).unwrap().clone())
            .unwrap();
    }
    for (_, old) in original.claims().claims_ordered() {
        let mut claim = old.clone();
        if claim.key.name() == "young_modulus" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(COLD, 2.0e7), (HOT, if side == 0 { 1.0e7 } else { 2.0e7 })],
                dims: Pressure::DIMS,
            };
            claim.validity = ValidityDomain::unconstrained().with("T", COLD, HOT);
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    (
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap(),
        key,
    )
}

pub(super) fn thermal_fixture(hardness_support: f64, coefficient_support: f64) -> Fixture {
    let data = [
        thermal_material(0, hardness_support),
        thermal_material(1, hardness_support),
    ];
    let cards = data.each_ref().map(|(card, _)| card.clone());
    let point = common::query();
    let select = MaterialPropertySelection::SingleClaimOnly;
    let profiles = data.each_ref().map(|(card, key)| {
        ResolvedDepthHardness::new(card, &point, key.clone(), select.clone(), MAX_DEPTH).unwrap()
    });
    let base = common::bind(&cards[0], &cards[1], None, 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    let keys = ["k_a", "k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
    for (key, k) in keys.clone().into_iter().zip(K) {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (COLD, k),
                        (
                            coefficient_support,
                            k * (1.0 + (coefficient_support - COLD) / (HOT - COLD)),
                        ),
                    ],
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained().with("T", COLD, coefficient_support),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
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
    let dry = bind_dry_interface_system_card(&interface, InputAuthority::SyntheticFixture).unwrap();
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
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select,
    )
    .unwrap();
    Fixture {
        cards,
        profiles,
        normal,
        coefficients,
        keys,
    }
}

fn controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}
fn slip(distance: f64) -> Vec<f64> {
    let mut slip = vec![0.0; 25];
    slip[12] = distance;
    slip
}
fn compliance(source: &SourceDepthHardnessWear) -> f64 {
    4.0 * source.model().config().grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (std::f64::consts::PI * source.model().config().reduced_modulus_pa)
}

#[test]
fn g1_independent_temperatures_soften_exposed_profile_and_change_future_wear() {
    with_cx(|cx, _| {
        let fixture = thermal_fixture(HOT, HOT);
        let before = fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap()
            .candidate
            .advance_wear_integrated(&slip(1.0e-5), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let temperatures = [HOT, COLD + 10.0];
        let interface_temperature = COLD + 15.0;
        let hot = before
            .refresh_at_temperatures(
                [&fixture.cards[0], &fixture.cards[1]],
                temperatures,
                interface_temperature,
                cx,
            )
            .unwrap();
        let after = &hot.candidate;
        for side in 0..2 {
            let z = before.recession_m()[side][12];
            close(
                after.hardness_profiles()[side].hardness_pa_at(z).unwrap(),
                h(side, temperatures[side], z),
            );
            assert_eq!(
                after.hardness_profiles()[side]
                    .reference_source()
                    .resolved()
                    .state_coordinate("T"),
                Some(temperatures[side])
            );
            assert_eq!(
                after.hardness_profiles()[side].reference_source().key(),
                before.hardness_profiles()[side].reference_source().key()
            );
            assert_eq!(
                after.hardness_profiles()[side]
                    .reference_source()
                    .property()
                    .answer()
                    .receipt
                    .selected,
                before.hardness_profiles()[side]
                    .reference_source()
                    .property()
                    .answer()
                    .receipt
                    .selected
            );
        }
        assert_eq!(
            after.coefficient_source().state_coordinate("T"),
            Some(interface_temperature)
        );
        close(
            after.model().config().hardness_pa[12],
            h(0, HOT, before.recession_m()[0][12]),
        );
        assert!(
            after.contact_state().plastic_depth_m()[12]
                > before.contact_state().plastic_depth_m()[12]
        );
        assert!(hot.receipt.response.plastic_work_j > 0.0);
        assert_eq!(after.volumes_m3(), before.volumes_m3());
        assert_eq!(after.recession_m(), before.recession_m());
        close(
            hot.receipt.parameter_work_j,
            hot.receipt.elastic_energy_change_j + hot.receipt.response.plastic_work_j,
        );
        let distance = 1.0e-6;
        let p = after.contact_state().pressure_pa()[12];
        let removed = after.advance_wear_integrated(&slip(distance), cx).unwrap();
        for side in 0..2 {
            let old = after.recession_m()[side][12];
            let hs = h(side, temperatures[side], old);
            let measure = 1.75 * K[side] * p * distance;
            let depth = 2.0 * measure / (hs + (hs * hs + 2.0 * DEPTH_SLOPE[side] * measure).sqrt());
            close(removed.recession_increment_m[side][12], depth);
            close(removed.wear[side].increment_m3, AREA * depth);
        }
        let cooled = after
            .refresh_at_temperatures([&fixture.cards[0], &fixture.cards[1]], [COLD; 2], COLD, cx)
            .unwrap()
            .candidate;
        assert_eq!(cooled.recession_m(), after.recession_m());
        assert_eq!(cooled.volumes_m3(), after.volumes_m3());
        assert_eq!(
            cooled.contact_state().plastic_depth_m(),
            after.contact_state().plastic_depth_m()
        );
        let noop = after
            .refresh_at_temperatures(
                [&fixture.cards[0], &fixture.cards[1]],
                temperatures,
                interface_temperature,
                cx,
            )
            .unwrap();
        assert_eq!(noop.candidate.source_identity(), after.source_identity());
        assert_eq!(noop.candidate.model().identity(), after.model().identity());
        assert_eq!(noop.candidate.recession_m(), after.recession_m());
        assert_eq!(noop.candidate.volumes_m3(), after.volumes_m3());
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_depth_thermal_force_and_spring_refresh_match_compliance_and_work() {
    with_cx(|cx, _| {
        let fixture = thermal_fixture(HOT, HOT);
        let pristine = fixture.owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx);
        let cold_c = compliance(&pristine);
        let force = 0.003;
        let p = force / AREA;
        let before = pristine
            .advance_normal_at_load(force, controls(), cx)
            .unwrap()
            .normal
            .candidate
            .advance_wear_integrated_at_load(&slip(1.0e-5), force, controls(), cx)
            .unwrap()
            .wear
            .candidate;
        let hot = before
            .refresh_at_temperatures_at_load(
                [&fixture.cards[0], &fixture.cards[1]],
                [HOT; 2],
                HOT,
                force,
                controls(),
                cx,
            )
            .unwrap();
        let after = &hot.thermal.candidate;
        let hot_c = 1.5 * cold_c;
        close(compliance(after), hot_c);
        close(after.contact_state().pressure_pa()[12], p);
        let da = (hot_c - cold_c) * p;
        close(
            after.contact_state().approach_m() - before.contact_state().approach_m(),
            da,
        );
        close(hot.thermal.receipt.normal_work_j, force * da);
        close(hot.thermal.receipt.parameter_work_j, -0.5 * force * da);
        assert_eq!(after.volumes_m3(), before.volumes_m3());
        assert_eq!(after.recession_m(), before.recession_m());
        assert_eq!(
            after.contact_state().plastic_depth_m(),
            before.contact_state().plastic_depth_m()
        );
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: AREA / cold_c,
            position_m: 2.0 * cold_c * p,
        };
        let before = pristine
            .advance_normal_with_drive(drive, controls(), cx)
            .unwrap()
            .normal
            .candidate
            .advance_wear_integrated_with_drive(&slip(1.0e-5), drive, controls(), cx)
            .unwrap()
            .wear
            .candidate;
        let hot = before
            .refresh_at_temperatures_with_drive(
                [&fixture.cards[0], &fixture.cards[1]],
                [HOT; 2],
                HOT,
                drive,
                controls(),
                cx,
            )
            .unwrap();
        let after = &hot.thermal.candidate;
        let recession = before.recession_m()[0][12] + before.recession_m()[1][12];
        let expected = (drive.position_m - recession) / (hot_c + AREA / drive.stiffness_n_per_m);
        close(after.contact_state().pressure_pa()[12], expected);
        let initial_extension = drive.position_m - before.contact_state().approach_m();
        let final_extension = drive.position_m - after.contact_state().approach_m();
        close(
            hot.thermal.receipt.normal_work_j,
            0.5 * drive.stiffness_n_per_m * (initial_extension.powi(2) - final_extension.powi(2)),
        );
        close(hot.drive.achieved_force_n, expected * AREA);
        assert_eq!(after.volumes_m3(), before.volumes_m3());
        assert_eq!(after.recession_m(), before.recession_m());
    });
}

#[test]
fn g0_depth_thermal_late_support_force_and_cancellation_refuse_atomically() {
    with_cx(|cx, gate| {
        for (hardness_support, coefficient_support) in [(COLD + 10.0, HOT), (HOT, COLD + 10.0)] {
            let fixture = thermal_fixture(hardness_support, coefficient_support);
            let before = fixture
                .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
                .advance_normal(2.0e-6, cx)
                .unwrap()
                .candidate
                .advance_wear_integrated(&slip(1.0e-5), cx)
                .unwrap()
                .candidate;
            let saved = before.clone();
            let error = before
                .refresh_at_temperatures([&fixture.cards[0], &fixture.cards[1]], [HOT; 2], HOT, cx)
                .unwrap_err();
            if hardness_support < HOT {
                assert!(matches!(
                    error,
                    DepthHardnessTemperatureError::Wear(DepthHardnessWearError::Hardness(_))
                ));
            } else {
                assert!(matches!(
                    error,
                    DepthHardnessTemperatureError::CoefficientSource(_)
                ));
            }
            assert_eq!(before, saved);
        }
        let fixture = thermal_fixture(HOT, HOT);
        let force = 0.009;
        let before = fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal_at_load(force, controls(), cx)
            .unwrap()
            .normal
            .candidate;
        let saved = before.clone();
        assert!(
            before
                .refresh_at_temperatures_at_load(
                    [&fixture.cards[0], &fixture.cards[1]],
                    [HOT; 2],
                    HOT,
                    force,
                    controls(),
                    cx
                )
                .is_err()
        );
        let foreign = Fixture::new(GRADIENT);
        assert!(
            before
                .refresh_at_temperatures(
                    [&foreign.cards[0], &fixture.cards[1]],
                    [COLD; 2],
                    COLD,
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(
            before
                .refresh_at_temperatures(
                    [&fixture.cards[0], &fixture.cards[1]],
                    [COLD; 2],
                    COLD,
                    cx
                )
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
