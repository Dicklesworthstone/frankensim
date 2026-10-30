//! G0/G1/G3 source data drive native magnetic memory; synthetic fits, not measurements.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
        hysteresis::{UniformHystereticCore, source::*},
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{Dims, QuantitySpec, Temperature};
use std::sync::Arc;

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
}
fn roles() -> HysteresisSources {
    HysteresisSources {
        background_permeability: PropertyKey::new("fitted_background", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("fitted_field_limit", FIELD_STRENGTH_DIMS),
        branches: (0..2)
            .map(|i| HysteresisBranchSources {
                susceptibility: PropertyKey::new(
                    format!("branch_{i}_slope"),
                    MAGNETIC_PERMEABILITY_DIMS,
                ),
                pinning_field: PropertyKey::new(format!("branch_{i}_pinning"), FIELD_STRENGTH_DIMS),
                saturation_polarization: PropertyKey::new(
                    format!("branch_{i}_saturation"),
                    FLUX_DENSITY_DIMS,
                ),
            })
            .collect(),
    }
}
fn card(change: impl Fn(&mut PropertyClaim)) -> Arc<MaterialCard> {
    let roles = roles();
    let mut values = vec![
        (roles.background_permeability, 0.01),
        (roles.field_limit, 30.0),
    ];
    for (i, branch) in roles.branches.into_iter().enumerate() {
        values.extend([
            (branch.susceptibility, if i == 0 { 0.2 } else { 0.3 }),
            (branch.pinning_field, if i == 0 { 2.0 } else { 1.0 }),
            (branch.saturation_polarization, 1.0 + i as f64),
        ]);
    }
    let mut claims = ClaimSet::new();
    for (key, value) in values {
        let dims = key.dims();
        let mut claim = PropertyClaim {
            key,
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                250.0,
                600.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic independent play fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        };
        change(&mut claim);
        claims.insert_claim(claim).unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "solid".into(),
                process: "fitted play test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}
fn resolve(card: Arc<MaterialCard>) -> ResolvedScalarHysteresis {
    ResolvedScalarHysteresis::resolve(
        card,
        &point(300.0),
        roles(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12, "{a:.16e} != {b:.16e}");
}

#[test]
fn source_branches_drive_native_memory_work_and_geometry_without_name_dispatch() {
    let fit = resolve(card(|_| {}));
    let law = fit.law();
    let first = law.candidate_field(&law.initial_state(), 4.0).unwrap();
    // z=(.4,.9), B=.01*4+.4+.9; psi=.01*4^2/2+.4^2/.4+.9^2/.6.
    near(first.state.flux_density_t(), 1.34);
    near(first.state.stored_energy_j_m3(), 1.83);
    near(first.irreversible_heat_j_m3, 1.7);
    near(first.magnetic_work_j_m3, 3.53);
    let reversal = law.candidate_field(&first.state, 0.0).unwrap();
    near(reversal.state.polarizations_t()[0], 0.4);
    near(reversal.state.polarizations_t()[1], 0.3);
    near(reversal.state.flux_density_t(), 0.7);
    near(reversal.irreversible_heat_j_m3, 0.6);
    near(reversal.state.stored_energy_j_m3(), 0.55);

    let core = UniformHystereticCore::from_source(fit, 3.0, 0.5, 2.0).unwrap();
    let step = core
        .candidate_current(&core.initial_state(), 8.0 / 3.0)
        .unwrap();
    near(step.flux_linkage_wb, 2.01);
    near(step.stored_energy_j, 1.83);
    near(step.core_heat_j, 1.7);
    let saved = step.material.state.clone();
    assert!(core.candidate_current(&saved, 21.0).is_err()); // H=31.5 exceeds source fit.
    assert_eq!(step.material.state, saved);
    let mut shared = roles();
    shared.branches = vec![shared.branches[0].clone(); 3];
    let repeated = ResolvedScalarHysteresis::resolve(
        card(|_| {}),
        &point(300.0),
        shared,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    near(
        repeated
            .law()
            .candidate_field(&repeated.law().initial_state(), 4.0)
            .unwrap()
            .state
            .flux_density_t(),
        1.24,
    );
}

#[test]
fn temperature_refresh_preserves_original_coefficients_selection_and_remanence() {
    let fit = resolve(card(|_| {}));
    fit.admit_temperature_interval("T", [250.0, 600.0]).unwrap();
    let core = UniformHystereticCore::from_source(fit, 1.0, 1.0, 1.0).unwrap();
    let prepared = core.candidate_current(&core.initial_state(), 4.0).unwrap();
    let warm = core.at_source_query(&point(500.0)).unwrap();
    let before = core.hold_candidate(&prepared.material.state).unwrap();
    let after = warm.hold_candidate(&prepared.material.state).unwrap();
    assert_eq!(before, after);
    assert_eq!(warm.source().unwrap().material().query(), point(500.0));
    for (a, b) in core
        .source()
        .unwrap()
        .material()
        .properties()
        .iter()
        .zip(warm.source().unwrap().material().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    assert_eq!(
        warm.source().unwrap().material().card_identity(),
        core.source().unwrap().material().card_identity()
    );
    assert!(core.at_source_query(&point(600.1)).is_err());
    assert!(
        core.source()
            .unwrap()
            .admit_temperature_interval("T", [249.0, 600.0])
            .is_err()
    );
    assert_eq!(core.source().unwrap().material().query(), point(300.0));
    let untyped = QueryPoint::new().with("T", 500.0).unwrap();
    assert!(core.at_source_query(&untyped).is_err());
}

#[test]
fn temperature_dependent_fits_dynamic_domains_and_incompatible_roles_refuse() {
    let resolve_result = |card| {
        ResolvedScalarHysteresis::resolve(
            card,
            &point(300.0),
            roles(),
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    let varying = card(|claim| {
        if claim.key.name() == "branch_0_slope" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.2), (600.0, 0.1)],
                dims: claim.key.dims(),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
    });
    assert!(resolve_result(varying).is_err());
    let frequency = card(|claim| {
        claim.validity = claim.validity.clone().with_quantity(
            "omega",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            0.0,
            100.0,
        );
    });
    assert!(
        ResolvedScalarHysteresis::resolve(
            frequency,
            &point(300.0)
                .with_quantity(
                    "omega",
                    QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                    0.0
                )
                .unwrap(),
            roles(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let unconditioned = resolve(card(|claim| {
        claim.validity = ValidityDomain::unconstrained()
    }));
    assert!(
        unconditioned
            .admit_temperature_interval("T", [250.0, 600.0])
            .is_err()
    );
    let mut wrong = roles();
    wrong.branches[0].susceptibility = PropertyKey::new("branch_0_slope", Dims([0; 6]));
    assert!(
        ResolvedScalarHysteresis::resolve(
            card(|_| {}),
            &point(300.0),
            wrong,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}
