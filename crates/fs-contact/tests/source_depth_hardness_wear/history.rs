//! Original native preconditioning and graded depth enter actual contact/wear.

#[path = "history/refresh.rs"]
mod refresh;
#[path = "history/sliding_refresh.rs"]
mod sliding_refresh;

use super::*;
use fs_matdb::{ClaimSet, MaterialStateId};
use fs_material::SmallStrainLaw;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY, resolve_j2_plasticity,
};
use fs_qty::Density;

const SLOPES_ALPHA: [f64; 2] = [4.0e7, 2.0e7];
const GRADIENTS: [f64; 2] = [GRADIENT, 1.0e11];
const GATES: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-6,
    yield_pa: 1e-7,
};

fn native_card(side: usize) -> (MaterialCard, PropertyKey) {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", 2.0e7, Pressure::DIMS),
        ("poisson_ratio", 0.25, Dims::NONE),
        ("yield_stress", 1.0e4, Pressure::DIMS),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, 1.0e5, Pressure::DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("synthetic preconditioned graded side {side}"),
            method: "synthetic force/projected-area indentation".into(),
            artifact: hash_domain("native-depth-contact", &[side as u8]),
            caveats: "numerical regression only".into(),
            provenance: provenance(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("native-depth-H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic rigid sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic depth/history indentation",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    let samples = [[0.0, 0.0], [MAX_DEPTH, 0.0], [MAX_DEPTH, 0.04], [0.0, 0.04]]
        .map(|x| {
            (
                x,
                [H0, HB][side] + GRADIENTS[side] * x[0] + SLOPES_ALPHA[side] * x[1],
            )
        })
        .to_vec();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Surface2D {
                axes: [DEPTH.into(), ALPHA.into()],
                axis_dims: [LENGTH, Dims::NONE],
                samples,
                triangles: vec![[0, 1, 2], [0, 2, 3]],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0)
                .with_quantity(DEPTH, QuantitySpec::dimensional(LENGTH), 0.0, MAX_DEPTH)
                .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.04),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    (
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic graded material {side}"),
                phase: "solid".into(),
                process: "uniform native preconditioning".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        key,
    )
}

fn fixture() -> (Fixture, [f64; 2]) {
    let data = [native_card(0), native_card(1)];
    let mut alphas = [0.0; 2];
    let profiles = std::array::from_fn(|side| {
        let (card, key) = &data[side];
        let plastic = resolve_j2_plasticity(
            card,
            &common::query(),
            MaterialPropertySelection::SingleClaimOnly,
            HardeningPropertyBasis::PlasticStrainModulus,
            0.05,
        )
        .unwrap();
        let axial = [0.01, 0.02][side];
        let strain = [axial, 0.0, 0.0, 0.0, 0.0, 0.0];
        let native = plastic
            .law()
            .advance_with_energy(&[0.0; 6], &plastic.initial_state(), &strain, GATES)
            .unwrap();
        alphas[side] = native.state.alpha;
        // Independent radial-return reference, with native tensor convention.
        close(alphas[side], (1.6e7 * axial - 1.0e4) / 2.41e7);
        ResolvedDepthHardness::from_j2_accepted(
            &plastic,
            &strain,
            &native.state,
            GATES.yield_pa,
            card,
            key.clone(),
            MaterialPropertySelection::SingleClaimOnly,
            MAX_DEPTH,
        )
        .unwrap()
    });
    (
        Fixture::from_profiles(data.map(|(card, _)| card), profiles),
        alphas,
    )
}

#[test]
fn g1_native_preconditioning_changes_caps_and_each_bodys_exact_graded_wear() {
    with_cx(|cx, _| {
        let (fixture, alphas) = fixture();
        let surface_h =
            std::array::from_fn::<_, 2, _>(|i| [H0, HB][i] + SLOPES_ALPHA[i] * alphas[i]);
        let loaded = fixture
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap();
        let before = loaded.candidate;
        let saved = before.clone();
        close(before.contact_state().pressure_pa()[12], surface_h[0]);
        assert!(surface_h[0] > H0);
        assert_eq!(before.volumes_m3(), [0.0; 2]);
        assert!(before.contact_state().plastic_depth_m()[12] > MAX_DEPTH);
        let mut slip = [0.0; 25];
        slip[12] = 1.0e-5;
        let worn = before.advance_wear_integrated(&slip, cx).unwrap();
        assert_eq!(before, saved);
        for side in 0..2 {
            // Independent exact root of H(alpha,0)*z + gradient*z²/2 = k*p*ds.
            let measure = K[side] * surface_h[0] * slip[12];
            let expected = 2.0 * measure
                / (surface_h[side]
                    + (surface_h[side] * surface_h[side] + 2.0 * GRADIENTS[side] * measure).sqrt());
            close(worn.recession_increment_m[side][12], expected);
            close(worn.candidate.recession_m()[side][12], expected);
            close(worn.candidate.volumes_m3()[side], AREA * expected);
            let profile = &worn.candidate.hardness_profiles()[side];
            assert_eq!(profile.accepted_history().unwrap().alpha(), alphas[side]);
            assert_eq!(profile.identity(), fixture.profiles[side].identity());
            let endpoint = profile
                .resolve_at_depth(&fixture.cards[side], expected)
                .unwrap();
            close(
                endpoint.hardness_pa(),
                surface_h[side] + GRADIENTS[side] * expected,
            );
            assert_eq!(
                endpoint.resolved().state_coordinate(ALPHA),
                Some(alphas[side])
            );
        }
        let cap = surface_h[0] + GRADIENTS[0] * worn.candidate.recession_m()[0][12];
        close(worn.candidate.model().config().hardness_pa[12], cap);
        close(worn.candidate.contact_state().pressure_pa()[12], cap);
        close(
            worn.hardness.elastic_energy_j - loaded.receipt.normal.reversible_energy_j
                + worn.recession.response.plastic_work_j
                + worn.hardness.response.plastic_work_j,
            worn.recession.parameter_work_j + worn.hardness.parameter_work_j,
        );
        assert_eq!(worn.recession.normal_work_j, 0.0);
        assert_eq!(worn.hardness.normal_work_j, 0.0);
        assert_eq!(worn.candidate.source_identity(), before.source_identity());
    });
}

#[test]
fn g0_history_profile_cannot_be_swapped_between_original_bodies() {
    with_cx(|cx, _| {
        let (fixture, _) = fixture();
        let result = SourceDepthHardnessWear::new(
            fixture.normal.clone(),
            fixture.coefficients.clone(),
            fixture.keys.clone(),
            [fixture.profiles[1].clone(), fixture.profiles[0].clone()],
            FiniteGapGrid {
                cells_x: 5,
                cells_y: 5,
                cell_width_m: 1e-4,
                cell_depth_m: 1e-4,
            },
            vec![0.0; 25],
            PlasticContactControls {
                maximum_active_set_iterations: 512,
                complementarity_tolerance_m: 1e-13,
                boundary_clearance_cells: 1,
            },
            PlasticContactEnvelope {
                maximum_approach_m: 1e-4,
                maximum_plastic_depth_m: 1e-4,
                half_space_depth_m: [1e-2; 2],
                maximum_patch_to_depth: [0.1; 2],
            },
            [1.5e-6; 2],
            FiniteGapPlasticLinearSolver::Dense,
            cx,
        );
        assert!(matches!(result, Err(DepthHardnessWearError::Invalid(_))));
    });
}
