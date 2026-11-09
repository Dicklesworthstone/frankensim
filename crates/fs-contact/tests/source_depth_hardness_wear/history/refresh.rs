//! Further native preconditioning retains removal and irreversible indentation.

use super::*;
use fs_contact::normal_patch::plastic_source::depth_wear::history::DepthHardnessHistoryRefresh;
use fs_contact::normal_patch::{FiniteGapPlasticDrive, FiniteGapPlasticLoadControls};

const AXIAL: [f64; 3] = [0.004, 0.012, 0.024];
const DEPTH_SLOPE: f64 = -1.0e11;

pub(super) struct HistoryFixture {
    pub(super) fixture: Fixture,
    pub(super) stages: [ResolvedDepthHardness; 3],
    key: PropertyKey,
    other_key: PropertyKey,
}

fn profile_material() -> (MaterialCard, PropertyKey, PropertyKey) {
    let (original, key) = native_card(0);
    let original_test = key.hardness_test().unwrap();
    let other_key = PropertyKey::with_quantity(key.name(), INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                original_test.indenter(),
                original_test.loading().to_vec(),
                "different original indentation protocol",
                original_test.observation(),
            )
            .unwrap(),
        )
        .unwrap();
    let mut claims = ClaimSet::new();
    for id in original.claims().observation_ids() {
        claims
            .register_observation(original.claims().observation(id).unwrap().clone())
            .unwrap();
    }
    for (_, old) in original.claims().claims_ordered() {
        let mut claim = old.clone();
        if claim.key == key {
            let mut samples = Vec::new();
            let mut triangles = Vec::new();
            for (alpha, surface) in [(0.0, 1.3e6), (0.008, 0.8e6), (0.016, 1.5e6), (0.04, 1.7e6)] {
                for z in [0.0, MAX_DEPTH] {
                    samples.push(([z, alpha], surface + DEPTH_SLOPE * z));
                }
            }
            for interval in 0..3 {
                let i = 2 * interval;
                triangles.extend([[i, i + 1, i + 3], [i, i + 3, i + 2]]);
            }
            claim.value = PropertyValue::Surface2D {
                axes: [DEPTH.into(), ALPHA.into()],
                axis_dims: [LENGTH, Dims::NONE],
                samples,
                triangles,
                dims: Pressure::DIMS,
            };
            let mut other = claim.clone();
            other.key = other_key.clone();
            claims.insert_claim(other).unwrap();
        }
        claims.insert_claim(claim).unwrap();
    }
    (
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap(),
        key,
        other_key,
    )
}

impl HistoryFixture {
    pub(super) fn new() -> Self {
        let (card, key, other_key) = profile_material();
        let plastic = resolve_j2_plasticity(
            &card,
            &common::query(),
            MaterialPropertySelection::SingleClaimOnly,
            HardeningPropertyBasis::PlasticStrainModulus,
            0.05,
        )
        .unwrap();
        let mut old_strain = [0.0; 6];
        let mut state = plastic.initial_state();
        let stages = std::array::from_fn(|stage| {
            let strain = [AXIAL[stage], 0.0, 0.0, 0.0, 0.0, 0.0];
            let next = plastic
                .law()
                .advance_with_energy(&old_strain, &state, &strain, GATES)
                .unwrap();
            let profile = ResolvedDepthHardness::from_j2_accepted(
                &plastic,
                &strain,
                &next.state,
                GATES.yield_pa,
                &card,
                key.clone(),
                MaterialPropertySelection::SingleClaimOnly,
                MAX_DEPTH,
            )
            .unwrap();
            old_strain = strain;
            state = next.state;
            profile
        });
        let (fixed_card, fixed_key) = material("fixed counterface", curve(HB, 0.0));
        let fixed = ResolvedDepthHardness::new(
            &fixed_card,
            &common::query(),
            fixed_key,
            MaterialPropertySelection::SingleClaimOnly,
            MAX_DEPTH,
        )
        .unwrap();
        let fixture = Fixture::from_profiles([card, fixed_card], [stages[0].clone(), fixed]);
        Self {
            fixture,
            stages,
            key,
            other_key,
        }
    }

    pub(super) fn profiles(&self, stage: usize) -> [ResolvedDepthHardness; 2] {
        [self.stages[stage].clone(), self.fixture.profiles[1].clone()]
    }

    fn variant(
        &self,
        key: &PropertyKey,
        strain_limit: f64,
        maximum_depth: f64,
        temperature: f64,
    ) -> ResolvedDepthHardness {
        let card = &self.fixture.cards[0];
        let plastic = resolve_j2_plasticity(
            card,
            &common::query().with("T", temperature).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            HardeningPropertyBasis::PlasticStrainModulus,
            strain_limit,
        )
        .unwrap();
        let strain = [AXIAL[1], 0.0, 0.0, 0.0, 0.0, 0.0];
        let next = plastic
            .law()
            .advance_with_energy(&[0.0; 6], &plastic.initial_state(), &strain, GATES)
            .unwrap();
        ResolvedDepthHardness::from_j2_accepted(
            &plastic,
            &strain,
            &next.state,
            GATES.yield_pa,
            card,
            key.clone(),
            MaterialPropertySelection::SingleClaimOnly,
            maximum_depth,
        )
        .unwrap()
    }
}

pub(super) fn controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1e-4,
        maximum_approach_m: 1e-4,
        absolute_force_tolerance_n: 1e-14,
        relative_force_tolerance: 1e-12,
        approach_tolerance_m: 1e-18,
        maximum_equilibrium_evaluations: 128,
    }
}
pub(super) fn slip(distance: f64) -> [f64; 25] {
    let mut values = [0.0; 25];
    values[12] = distance;
    values
}
pub(super) fn compliance(owner: &SourceDepthHardnessWear) -> f64 {
    4.0 * owner.model().config().grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * owner.model().config().reduced_modulus_pa)
}
fn retained(before: &SourceDepthHardnessWear, refresh: &DepthHardnessHistoryRefresh) {
    let next = &refresh.candidate;
    assert_eq!(refresh.original_source_identity, before.source_identity());
    assert_eq!(next.recession_m(), before.recession_m());
    assert_eq!(next.volumes_m3(), before.volumes_m3());
    assert_eq!(
        next.contact_state().recession_m(),
        before.contact_state().recession_m()
    );
    assert!(
        next.contact_state()
            .plastic_depth_m()
            .iter()
            .zip(before.contact_state().plastic_depth_m())
            .all(|(new, old)| new >= old)
    );
    close(
        next.contact_state().cumulative_plastic_work_j()
            - before.contact_state().cumulative_plastic_work_j(),
        refresh.receipt.response.plastic_work_j,
    );
    assert_eq!(next.normal_source(), before.normal_source());
    assert_eq!(next.coefficient_source(), before.coefficient_source());
    assert_eq!(next.coefficient_keys(), before.coefficient_keys());
    assert_eq!(
        next.model().config().undeformed_gap_m,
        before.model().config().undeformed_gap_m
    );
    assert_eq!(
        refresh.receipt.predecessor_state_identity,
        before.contact_state().identity()
    );
    assert_eq!(
        refresh.receipt.successor_state_identity,
        next.contact_state().identity()
    );
    let residual = refresh.receipt.elastic_energy_change_j
        + refresh.receipt.response.plastic_work_j
        - refresh.receipt.normal_work_j
        - refresh.receipt.parameter_work_j;
    let scale = refresh.receipt.elastic_energy_j.abs()
        + refresh.receipt.elastic_energy_change_j.abs()
        + refresh.receipt.response.plastic_work_j.abs()
        + refresh.receipt.normal_work_j.abs()
        + refresh.receipt.parameter_work_j.abs();
    assert!(residual.abs() <= 256.0 * f64::EPSILON * scale.max(1e-20));
}

#[test]
fn g1_softened_then_rehardened_history_keeps_wear_and_cannot_heal_the_dent() {
    with_cx(|cx, _| {
        let f = HistoryFixture::new();
        let before = f
            .fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1e-4, cx)
            .advance_normal(5e-5, cx)
            .unwrap()
            .candidate
            .advance_wear_integrated(&slip(1e-5), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let z = before.recession_m()[0][12];
        assert!(z > 0.0);
        let p0 = before.contact_state().pressure_pa()[12];
        let expected = f.stages[1].hardness_pa_at(z).unwrap();
        assert!(expected < p0);
        let soft = before.refresh_hardness_profiles(f.profiles(1), cx).unwrap();
        retained(&before, &soft);
        close(soft.candidate.model().config().hardness_pa[12], expected);
        close(soft.candidate.contact_state().pressure_pa()[12], expected);
        let dent = compliance(&before) * (p0 - expected);
        close(
            soft.candidate.contact_state().plastic_depth_m()[12]
                - before.contact_state().plastic_depth_m()[12],
            dent,
        );
        close(soft.receipt.response.plastic_work_j, AREA * expected * dent);
        close(
            soft.receipt.parameter_work_j,
            -0.5 * AREA * compliance(&before) * (p0 - expected).powi(2),
        );
        assert_eq!(soft.receipt.normal_work_j, 0.0);
        assert_ne!(soft.candidate.source_identity(), before.source_identity());
        let hard = soft
            .candidate
            .refresh_hardness_profiles(f.profiles(2), cx)
            .unwrap();
        retained(&soft.candidate, &hard);
        assert!(hard.candidate.model().config().hardness_pa[12] > expected);
        // Raising the cap only enlarges the admissible pressure range. The same
        // approach and retained permanent dent still give the softened pressure.
        close(hard.candidate.contact_state().pressure_pa()[12], expected);
        assert_eq!(
            hard.candidate.contact_state().plastic_depth_m(),
            soft.candidate.contact_state().plastic_depth_m()
        );
        assert_eq!(hard.receipt.response.plastic_work_j, 0.0);
        assert_eq!(hard.receipt.normal_work_j, 0.0);
        let next = hard
            .candidate
            .advance_wear_integrated(&slip(1e-6), cx)
            .unwrap();
        let hs = hard.candidate.hardness_profiles()[0]
            .hardness_pa_at(z)
            .unwrap();
        let measure = K[0] * expected * 1e-6;
        let removal = 2.0 * measure / (hs + (hs * hs + 2.0 * DEPTH_SLOPE * measure).sqrt());
        close(next.recession_increment_m[0][12], removal);
        let noop = hard
            .candidate
            .refresh_hardness_profiles(f.profiles(2), cx)
            .unwrap();
        assert_eq!(
            noop.candidate.source_identity(),
            hard.candidate.source_identity()
        );
        assert_eq!(
            noop.candidate.model().identity(),
            hard.candidate.model().identity()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_direct_force_and_spring_profile_refresh_keep_their_original_boundary_work() {
    with_cx(|cx, _| {
        let f = HistoryFixture::new();
        let pristine = f
            .fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1e-4, cx);
        let force = 0.003;
        let before = pristine
            .advance_normal_at_load(force, controls(), cx)
            .unwrap()
            .normal
            .candidate
            .advance_wear_integrated_at_load(&slip(1e-5), force, controls(), cx)
            .unwrap()
            .wear
            .candidate;
        let soft = before
            .refresh_hardness_profiles_at_load(f.profiles(1), force, controls(), cx)
            .unwrap();
        retained(&before, &soft.history);
        close(
            soft.history.candidate.contact_state().pressure_pa()[12],
            force / AREA,
        );
        close(soft.load.achieved_force_n, force);
        close(
            soft.history.candidate.contact_state().approach_m(),
            before.contact_state().approach_m(),
        );
        assert_eq!(
            soft.history.candidate.contact_state().plastic_depth_m(),
            before.contact_state().plastic_depth_m()
        );
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: AREA / compliance(&pristine),
            position_m: 7e-5,
        };
        let before = pristine
            .advance_normal_with_drive(drive, controls(), cx)
            .unwrap()
            .normal
            .candidate
            .advance_wear_integrated_with_drive(&slip(1e-5), drive, controls(), cx)
            .unwrap()
            .wear
            .candidate;
        let soft = before
            .refresh_hardness_profiles_with_drive(f.profiles(1), drive, controls(), cx)
            .unwrap();
        retained(&before, &soft.history);
        let after = &soft.history.candidate;
        let pressure = f.stages[1]
            .hardness_pa_at(before.recession_m()[0][12])
            .unwrap();
        close(after.contact_state().pressure_pa()[12], pressure);
        close(
            after.contact_state().approach_m(),
            drive.position_m - AREA * pressure / drive.stiffness_n_per_m,
        );
        close(soft.drive.achieved_force_n, AREA * pressure);
        let old_extension = drive.position_m - before.contact_state().approach_m();
        let new_extension = drive.position_m - after.contact_state().approach_m();
        close(
            soft.history.receipt.normal_work_j,
            0.5 * drive.stiffness_n_per_m * (old_extension.powi(2) - new_extension.powi(2)),
        );
        assert!(soft.history.receipt.response.plastic_work_j > 0.0);
        assert!(soft.history.receipt.normal_work_j > 0.0);
    });
}

#[test]
fn g0_profile_pins_alpha_fixed_side_force_capacity_and_cancellation_refuse_atomically() {
    with_cx(|cx, gate| {
        let f = HistoryFixture::new();
        let before = f
            .fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1e-4, cx)
            .advance_normal(5e-5, cx)
            .unwrap()
            .candidate
            .advance_wear_integrated(&slip(1e-5), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        for profile in [
            f.variant(&f.other_key, 0.05, MAX_DEPTH, 293.15),
            f.variant(&f.key, 0.04, MAX_DEPTH, 293.15),
            f.variant(&f.key, 0.05, 0.9 * MAX_DEPTH, 293.15),
            f.variant(&f.key, 0.05, MAX_DEPTH, 294.15),
        ] {
            assert!(matches!(
                before.refresh_hardness_profiles([profile, f.fixture.profiles[1].clone()], cx),
                Err(DepthHardnessWearError::Invalid(_))
            ));
            assert_eq!(before, saved);
        }
        let changed_fixed = ResolvedDepthHardness::new(
            &f.fixture.cards[1],
            &common::query(),
            f.fixture.profiles[1].reference_source().key().clone(),
            MaterialPropertySelection::SingleClaimOnly,
            0.9 * MAX_DEPTH,
        )
        .unwrap();
        assert!(
            before
                .refresh_hardness_profiles([f.stages[1].clone(), changed_fixed], cx)
                .is_err()
        );
        assert!(
            before
                .refresh_hardness_profiles([f.fixture.profiles[1].clone(), f.stages[1].clone()], cx)
                .is_err()
        );
        let soft = before
            .refresh_hardness_profiles(f.profiles(1), cx)
            .unwrap()
            .candidate;
        let soft_saved = soft.clone();
        assert!(soft.refresh_hardness_profiles(f.profiles(0), cx).is_err());
        assert_eq!(soft, soft_saved);
        let force = 0.01; // Between original and softened single-cell capacities.
        let loaded = f
            .fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1e-4, cx)
            .advance_normal_at_load(force, controls(), cx)
            .unwrap()
            .normal
            .candidate;
        let loaded_saved = loaded.clone();
        assert!(matches!(
            loaded.refresh_hardness_profiles_at_load(f.profiles(1), force, controls(), cx),
            Err(DepthHardnessWearError::Load(_))
        ));
        assert_eq!(loaded, loaded_saved);
        gate.request();
        assert!(before.refresh_hardness_profiles(f.profiles(1), cx).is_err());
        assert_eq!(before, saved);
    });
}
