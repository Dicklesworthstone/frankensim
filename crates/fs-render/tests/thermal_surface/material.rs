use super::*;
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, LinearizedSurfaceRadiation, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_render::material_state::{MaterialThermalBinding, MaterialThermalError};

fn card(epsilon: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: epsilon,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 2500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic gray-surface regression".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "declared gray-diffuse finish".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn state(card: &MaterialCard, t: f64) -> fs_material::state_point::ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        &QueryPoint::new().with("T", t).unwrap(),
        &[ScalarPropertyRequirement::try_new(
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn g3_same_card_and_temperature_drive_cooling_and_actual_rendered_radiance() {
    with_cx(|cx| {
        let mut pixels = Vec::new();
        let mut powers = Vec::new();
        for epsilon in [0.2, 0.6] {
            let source = card(epsilon);
            let accepted = state(&source, 1800.0);
            let emission = SurfaceEmissivity::from_card(
                "body",
                &source,
                1800.0,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
            let binding = MaterialThermalBinding::try_gray_diffuse(
                &accepted,
                &emission,
                ContentHash([9; 32]),
            )
            .unwrap();
            assert_eq!(binding.emissivity(), &emission);
            let cooling =
                LinearizedSurfaceRadiation::new("body", emission.clone(), 1790.0, 1780.0, 30.0)
                    .unwrap()
                    .evaluate(1800.0)
                    .unwrap();
            let ambient =
                GrayDiffuseThermalSurface::try_new(1780.0, epsilon, ContentHash([9; 32])).unwrap();
            let net =
                binding.surface().bolometric_exitance_w_m2() - ambient.bolometric_exitance_w_m2();
            assert!((net / cooling.nonlinear_outward_flux_w_m2 - 1.0).abs() < 1.0e-13);
            powers.push(net);
            pixels.push(
                render(
                    &scene(binding.material()),
                    cx,
                    &settings(DirectStrategy::Mis),
                )
                .unwrap()
                .xyz[0],
            );
            assert!(matches!(
                MaterialThermalBinding::try_gray_diffuse(
                    &state(&source, 1801.0),
                    &emission,
                    ContentHash([9; 32])
                ),
                Err(MaterialThermalError::StateMismatch)
            ));
            assert!(matches!(
                MaterialThermalBinding::try_gray_diffuse(
                    &state(&card(0.7), 1800.0),
                    &emission,
                    ContentHash([9; 32])
                ),
                Err(MaterialThermalError::StateMismatch)
            ));
        }
        assert!((powers[1] / powers[0] - 3.0).abs() < 1.0e-12);
        for channel in 0..3 {
            assert!((pixels[1][channel] / pixels[0][channel] - 3.0).abs() < 1.0e-11);
        }
    });
}
