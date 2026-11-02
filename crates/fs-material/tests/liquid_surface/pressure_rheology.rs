//! A pressure response requires its own explicit physical coordinate owner.

use super::*;
use fs_qty::{Frequency, Pressure};

fn pressure_source(kind: QuantityKind) -> (MaterialCard, QueryPoint) {
    let original = card(1000.0, 0.1, LiquidViscosityBasis::Dynamic, None);
    let pressure = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.validity = claim
            .validity
            .with_quantity("pressure", pressure, 90_000.0, 110_000.0);
        if claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "pressure".into(),
                abscissa_dims: Pressure::DIMS,
                knots: vec![(90_000.0, 0.1), (110_000.0, 0.3)],
                dims: DYNAMIC_VISCOSITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let source = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let query = point()
        .with_quantity("pressure", pressure, 100_000.0)
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    (source, query)
}

#[test]
fn g0_pressure_admission_retains_required_owner_and_rejects_stress_semantics() {
    for kind in [QuantityKind::Pressure, QuantityKind::Stress] {
        let (source, query) = pressure_source(kind);
        let requirements = [
            (DENSITY_PROPERTY, Density::DIMS),
            (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        ]
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap()
        });
        let bundle = resolve_material_state_point(
            &source,
            &query,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
                bundle.clone(),
                LiquidViscosityBasis::Dynamic,
                "temperature",
                "rate",
            )
            .is_err()
        );
        let result = ResolvedLiquidSurfaceState::try_for_pressure_shear_rate_evolution(
            bundle.clone(),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected,
            "temperature",
            "rate",
            "pressure",
        );
        if kind == QuantityKind::Stress {
            assert!(
                result.is_err(),
                "pressure dimensions do not turn stress into pressure"
            );
            continue;
        }
        let admitted = result.unwrap();
        assert_eq!(admitted.rheology_axes(), Some(("temperature", "rate")));
        assert_eq!(admitted.rheology_pressure_axis(), Some("pressure"));
        assert_eq!(admitted.material(), &bundle);
        assert!((admitted.dynamic_viscosity_pa_s() - 0.2).abs() < 1e-15);
        assert!(
            ResolvedLiquidSurfaceState::try_for_pressure_shear_rate_evolution(
                bundle.clone(),
                LiquidViscosityBasis::Dynamic,
                LiquidCapillarity::Neglected,
                "temperature",
                "rate",
                "temperature",
            )
            .is_err()
        );
        let outside = query.with("pressure", 120_000.0).unwrap();
        assert!(bundle.at_query_point(&source, &outside).is_err());
    }
}
