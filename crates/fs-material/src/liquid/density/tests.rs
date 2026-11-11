use super::*;
use crate::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY, VOLUMETRIC_EXPANSION_DIMS, VOLUMETRIC_EXPANSION_PROPERTY,
    resolve_liquid_state, resolve_liquid_state_from_kinematic, resolve_liquid_state_with_expansion,
};
use crate::liquid_temperature::LiquidViscosityBasis;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, MaterialStateId, PropertyClaim, PropertyKey, Provenance, UncertaintyModel,
};

const RHO_P: f64 = 1.0 / 1_048_576.0;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", typed(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("p", typed(QuantityKind::Pressure), p)
        .unwrap()
        .with_quantity("salinity", QuantitySpec::dimensional(Dims([0; 6])), 0.25)
        .unwrap()
}

fn validity(t: [f64; 2], p: [f64; 2]) -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity("T", typed(QuantityKind::AbsoluteTemperature), t[0], t[1])
        .with_quantity("p", typed(QuantityKind::Pressure), p[0], p[1])
        .with_quantity(
            "salinity",
            QuantitySpec::dimensional(Dims([0; 6])),
            0.0,
            1.0,
        )
}

fn surface(reverse: bool, dims: Dims, value: impl Fn(f64, f64) -> f64) -> PropertyValue {
    sampled_surface([256.0, 384.0], [65_536.0, 196_608.0], reverse, dims, value)
}

// The same two cells, sampled on the four corners of one T/p box.
fn sampled_surface(
    t_box: [f64; 2],
    p_box: [f64; 2],
    reverse: bool,
    dims: Dims,
    value: impl Fn(f64, f64) -> f64,
) -> PropertyValue {
    let mut axes = ["T".into(), "p".into()];
    let mut axis_dims = [Temperature::DIMS, Pressure::DIMS];
    if reverse {
        axes.swap(0, 1);
        axis_dims.swap(0, 1);
    }
    let samples = [
        [t_box[0], p_box[0]],
        [t_box[1], p_box[0]],
        [t_box[1], p_box[1]],
        [t_box[0], p_box[1]],
    ]
    .map(|[t, p]| (if reverse { [p, t] } else { [t, p] }, value(t, p)))
    .to_vec();
    PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples,
        triangles: vec![[0, 1, 2], [0, 2, 3]],
        dims,
    }
}

fn density_law(rho_t: f64, rho_p: f64) -> impl Fn(f64, f64) -> f64 {
    move |t: f64, p: f64| 1000.0 + rho_t * (t - 256.0) + rho_p * (p - 65_536.0)
}

fn density_surface(rho_t: f64, rho_p: f64, reverse: bool) -> PropertyValue {
    surface(reverse, Density::DIMS, density_law(rho_t, rho_p))
}

fn card_with_density(
    density: PropertyValue,
    interpolation: InterpolationPolicy,
    density_validity: ValidityDomain,
    basis: LiquidViscosityBasis,
    expansion: Option<f64>,
) -> MaterialCard {
    let mut properties = vec![
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            density,
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            surface(false, SPECIFIC_HEAT_CAPACITY_DIMS, |t, p| {
                4000.0 + 2.0 * (t - 256.0) + (p - 65_536.0) / 4096.0
            }),
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            typed(QuantityKind::ThermalConductivity),
            surface(false, THERMAL_CONDUCTIVITY_DIMS, |t, p| {
                0.5 + (t - 256.0) / 1024.0 + (p - 65_536.0) / 1_048_576.0
            }),
        ),
    ];
    properties.push(match basis {
        LiquidViscosityBasis::Dynamic => (
            DYNAMIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(DYNAMIC_VISCOSITY_DIMS),
            surface(false, DYNAMIC_VISCOSITY_DIMS, |t, p| {
                1.0 / 1024.0 + (t - 256.0) / 1_048_576.0 + (p - 65_536.0) / 1_073_741_824.0
            }),
        ),
        LiquidViscosityBasis::Kinematic => (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            surface(false, KINEMATIC_VISCOSITY_DIMS, |t, p| {
                1.0 / 1_048_576.0
                    + (t - 256.0) / 1_073_741_824.0
                    + (p - 65_536.0) / 1_099_511_627_776.0
            }),
        ),
    });
    if let Some(value) = expansion {
        properties.push((
            VOLUMETRIC_EXPANSION_PROPERTY,
            QuantitySpec::dimensional(VOLUMETRIC_EXPANSION_DIMS),
            PropertyValue::Scalar {
                value,
                dims: VOLUMETRIC_EXPANSION_DIMS,
            },
        ));
    }
    let mut claims = ClaimSet::new();
    for (name, quantity, value) in properties {
        let scalar = matches!(value, PropertyValue::Scalar { .. });
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value,
                validity: if name == DENSITY_PROPERTY {
                    density_validity.clone()
                } else {
                    validity([256.0, 384.0], [65_536.0, 196_608.0])
                },
                interpolation: if name == DENSITY_PROPERTY {
                    interpolation
                } else if scalar {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic local density law, not measured liquid data".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic density fixture".into(),
            phase: "liquid".into(),
            process: "affine single-phase test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn card(rho_t: f64, rho_p: f64, reverse: bool, basis: LiquidViscosityBasis) -> MaterialCard {
    card_with_density(
        density_surface(rho_t, rho_p, reverse),
        InterpolationPolicy::LinearInside,
        validity([256.0, 384.0], [65_536.0, 196_608.0]),
        basis,
        None,
    )
}

fn liquid(
    card: &MaterialCard,
    point: &QueryPoint,
    basis: LiquidViscosityBasis,
) -> ResolvedLiquidState {
    match basis {
        LiquidViscosityBasis::Dynamic => {
            resolve_liquid_state(card, point, MaterialPropertySelection::SingleClaimOnly)
        }
        LiquidViscosityBasis::Kinematic => resolve_liquid_state_from_kinematic(
            card,
            point,
            MaterialPropertySelection::SingleClaimOnly,
        ),
    }
    .unwrap()
}

fn state(
    card: &MaterialCard,
    point: &QueryPoint,
    basis: LiquidViscosityBasis,
) -> ResolvedLiquidDensityState {
    ResolvedLiquidDensityState::try_from_liquid(card, liquid(card, point, basis), "T", "p").unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-12),
        "{actual} != {expected}"
    );
}

#[test]
fn affine_density_closes_caloric_acoustic_and_buoyancy_views_in_both_axis_orders() {
    for reverse in [false, true] {
        let card = card(-0.5, RHO_P, reverse, LiquidViscosityBasis::Dynamic);
        let original = liquid(
            &card,
            &point(320.0, 114_688.0),
            LiquidViscosityBasis::Dynamic,
        );
        let state =
            ResolvedLiquidDensityState::try_from_liquid(&card, original.clone(), "T", "p").unwrap();
        assert_eq!(state.liquid(), &original);
        assert_eq!(state.material().properties().len(), 4);
        assert!(state.liquid().volumetric_expansion_per_k().is_none());
        let response = state.response();
        let rho = 1000.0 - 0.5 * 64.0 + RHO_P * 49_152.0;
        let cp = 4000.0 + 2.0 * 64.0 + 49_152.0 / 4096.0;
        let beta = 0.5 / rho;
        let kappa_t = RHO_P / rho;
        let cv = cp - 320.0 * beta * beta / RHO_P;
        let kappa_s = kappa_t * cv / cp;
        assert_eq!(response.density_kg_m3(), original.density_kg_m3());
        close(response.density_kg_m3(), rho);
        assert_eq!(response.density_temperature_derivative().value(), -0.5);
        assert_eq!(response.density_pressure_derivative().value(), RHO_P);
        close(response.volumetric_expansion_per_k(), beta);
        close(response.isothermal_compressibility_per_pa(), kappa_t);
        close(response.specific_heat_cv_j_kg_k(), cv);
        close(response.adiabatic_compressibility_per_pa(), kappa_s);
        close(response.sound_speed_m_s(), (1.0 / (rho * kappa_s)).sqrt());
        assert!(response.sound_speed_m_s() > (1.0 / RHO_P).sqrt());
        let acoustic = state.acoustic_state().unwrap();
        let transport = state.transport_state().unwrap();
        let thermal = state.thermal_transport_state().unwrap();
        let tv = state.thermoviscous_state().unwrap();
        for density in [
            acoustic.density_kg_m3(),
            transport.density_kg_m3(),
            thermal.transport().density_kg_m3(),
            tv.acoustic().density_kg_m3(),
        ] {
            assert_eq!(density.to_bits(), original.density_kg_m3().to_bits());
        }
        assert_eq!(acoustic.sound_speed_m_s(), response.sound_speed_m_s());
        assert_eq!(thermal.volumetric_expansion_per_k(), Some(beta));
        assert_eq!(
            tv.specific_heat_cv_j_kg_k(),
            response.specific_heat_cv_j_kg_k()
        );
        assert_eq!(tv.heat_capacity_ratio(), response.heat_capacity_ratio());
        assert_eq!(transport, original.transport_state().unwrap());
    }
}

#[test]
fn simultaneous_pressure_temperature_refresh_retains_claims_and_shared_preparation() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let card = card(-0.5, RHO_P, true, basis);
        let initial = state(&card, &point(288.0, 98_304.0), basis);
        let next_point = point(352.0, 163_840.0);
        let next = initial.at_query_point(&card, &next_point).unwrap();
        let independently_resolved = liquid(&card, &next_point, basis);
        assert_eq!(
            next.material().query(),
            independently_resolved.material().query()
        );
        // A refresh records pinned-claim selection; an independent resolution
        // records SingleClaimOnly. Their physical tuples, not policy receipts,
        // must agree exactly.
        assert_eq!(
            next.transport_state().unwrap(),
            independently_resolved.transport_state().unwrap()
        );
        assert!(Arc::ptr_eq(&initial.source, &next.source));
        assert_eq!(next.temperature_axis(), "T");
        assert_eq!(next.pressure_axis(), "p");
        for (old, new) in initial
            .material()
            .properties()
            .iter()
            .zip(next.material().properties())
        {
            assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
            assert_ne!(old.value_si(), new.value_si());
        }
        assert_ne!(initial.response(), next.response());
        assert_eq!(
            next.response().specific_heat_cp_j_kg_k(),
            next.liquid().specific_heat_j_kg_k()
        );
        assert_eq!(
            next.thermal_transport_state().unwrap().temperature_k(),
            352.0
        );
        let transport = next.transport_state().unwrap();
        close(
            transport.kinematic_viscosity_m2_s(),
            transport.dynamic_viscosity_pa_s() / next.response().density_kg_m3(),
        );
        if basis == LiquidViscosityBasis::Kinematic {
            let source_nu = next
                .material()
                .property(KINEMATIC_VISCOSITY_PROPERTY)
                .unwrap()
                .value_si();
            assert_eq!(
                transport.dynamic_viscosity_pa_s(),
                next.response().density_kg_m3() * source_nu
            );
        }
    }
}

#[test]
fn negative_and_zero_expansion_retain_their_physics() {
    for rho_t in [0.5, 0.0] {
        let card = card(rho_t, RHO_P, false, LiquidViscosityBasis::Dynamic);
        let state = state(
            &card,
            &point(320.0, 114_688.0),
            LiquidViscosityBasis::Dynamic,
        );
        let response = state.response();
        let beta = state
            .thermal_transport_state()
            .unwrap()
            .volumetric_expansion_per_k()
            .unwrap();
        if rho_t == 0.0 {
            assert_eq!(beta, 0.0);
            assert_eq!(
                response.specific_heat_cv_j_kg_k(),
                response.specific_heat_cp_j_kg_k()
            );
            assert_eq!(
                response.adiabatic_compressibility_per_pa(),
                response.isothermal_compressibility_per_pa()
            );
            assert_eq!(response.heat_capacity_ratio(), 1.0);
        } else {
            assert!(beta < 0.0);
            assert!(response.specific_heat_cv_j_kg_k() < response.specific_heat_cp_j_kg_k());
            assert!(
                state
                    .transport_state()
                    .unwrap()
                    .rayleigh(9.81, beta, 10.0, 0.1)
                    .unwrap()
                    < 0.0
            );
        }
    }
}

#[test]
fn mechanically_or_calorically_unstable_density_slopes_refuse() {
    for (rho_t, rho_p) in [(-0.5, 0.0), (-0.5, -RHO_P), (-0.5, RHO_P / 1024.0)] {
        let card = card(rho_t, rho_p, false, LiquidViscosityBasis::Dynamic);
        let original = liquid(
            &card,
            &point(320.0, 114_688.0),
            LiquidViscosityBasis::Dynamic,
        );
        assert!(matches!(
            ResolvedLiquidDensityState::try_from_liquid(&card, original, "T", "p"),
            Err(LiquidDensityError::Response(
                DensityResponseErrorV1::Unstable(_)
            ))
        ));
    }
}

#[test]
fn source_edges_validity_endpoints_and_unsupported_points_do_not_gain_derivatives() {
    let card = card(-0.5, RHO_P, false, LiquidViscosityBasis::Dynamic);
    let initial = state(
        &card,
        &point(320.0, 114_688.0),
        LiquidViscosityBasis::Dynamic,
    );
    for (t, p) in [
        (256.0, 114_688.0),
        (384.0, 114_688.0),
        (320.0, 65_536.0),
        (320.0, 196_608.0),
        (400.0, 114_688.0),
    ] {
        assert!(initial.at_query_point(&card, &point(t, p)).is_err());
    }
    // Declared T/p validity endpoints refuse
    // (crates/fs-material/contracts/part-04.md:4957-4958). Since b323fed42a,
    // fs-matdb admits a surface claim only when its finite validity bounds
    // enclose every source sample (crates/fs-matdb/CONTRACT.md:348-349), so a
    // validity box can no longer end inside a source cell. The restricted card
    // keeps the same affine density law and the same two cells, sampled on the
    // corners of its own validity box; its validity endpoints are source edges.
    let restricted = card_with_density(
        sampled_surface(
            [280.0, 360.0],
            [90_000.0, 180_000.0],
            false,
            Density::DIMS,
            density_law(-0.5, RHO_P),
        ),
        InterpolationPolicy::LinearInside,
        validity([280.0, 360.0], [90_000.0, 180_000.0]),
        LiquidViscosityBasis::Dynamic,
        None,
    );
    for (t, p) in [
        (280.0, 114_688.0),
        (360.0, 114_688.0),
        (320.0, 90_000.0),
        (320.0, 180_000.0),
    ] {
        let original = liquid(&restricted, &point(t, p), LiquidViscosityBasis::Dynamic);
        assert!(
            ResolvedLiquidDensityState::try_from_liquid(&restricted, original, "T", "p").is_err()
        );
    }
    // A single declared triangle leaves a hole within the claimed rectangle.
    let mut partial = density_surface(-0.5, RHO_P, false);
    if let PropertyValue::Surface2D { triangles, .. } = &mut partial {
        triangles.truncate(1);
    }
    let partial = card_with_density(
        partial,
        InterpolationPolicy::LinearInside,
        validity([256.0, 384.0], [65_536.0, 196_608.0]),
        LiquidViscosityBasis::Dynamic,
        None,
    );
    let initial = state(
        &partial,
        &point(352.0, 98_304.0),
        LiquidViscosityBasis::Dynamic,
    );
    assert!(
        initial
            .at_query_point(&partial, &point(280.0, 163_840.0))
            .is_err()
    );
}

#[test]
fn unequal_cell_gradients_and_source_knots_refuse_but_affine_seams_work() {
    let card = card(-0.5, RHO_P, false, LiquidViscosityBasis::Dynamic);
    // The source diagonal has identical gradients on both sides.
    state(
        &card,
        &point(320.0, 131_072.0),
        LiquidViscosityBasis::Dynamic,
    );
    for knot in [false, true] {
        let mut density = density_surface(-0.5, RHO_P, false);
        if let PropertyValue::Surface2D {
            samples, triangles, ..
        } = &mut density
        {
            if knot {
                samples.push(([320.0, 131_072.0], 969.0625));
                *triangles = vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
            } else {
                samples[1].1 -= 1.0;
            }
        }
        let card = card_with_density(
            density,
            InterpolationPolicy::LinearInside,
            validity([256.0, 384.0], [65_536.0, 196_608.0]),
            LiquidViscosityBasis::Dynamic,
            None,
        );
        let original = liquid(
            &card,
            &point(320.0, 131_072.0),
            LiquidViscosityBasis::Dynamic,
        );
        assert!(matches!(
            ResolvedLiquidDensityState::try_from_liquid(&card, original, "T", "p"),
            Err(LiquidDensityError::Source(
                MaterialStatePointError::Query { .. }
            ))
        ));
    }
}

#[test]
fn original_card_claim_context_and_expansion_authority_are_preserved() {
    let card = card(-0.5, RHO_P, false, LiquidViscosityBasis::Dynamic);
    let initial = state(
        &card,
        &point(320.0, 114_688.0),
        LiquidViscosityBasis::Dynamic,
    );
    let foreign = card_with_density(
        density_surface(-0.25, RHO_P, false),
        InterpolationPolicy::LinearInside,
        validity([256.0, 384.0], [65_536.0, 196_608.0]),
        LiquidViscosityBasis::Dynamic,
        None,
    );
    assert!(
        initial
            .at_query_point(&foreign, &point(320.0, 114_688.0))
            .is_err()
    );
    assert!(
        ResolvedLiquidDensityState::try_from_liquid(&foreign, initial.liquid().clone(), "T", "p")
            .is_err()
    );
    assert!(
        initial
            .at_query_point(
                &card,
                &point(320.0, 114_688.0).with("salinity", 0.3).unwrap()
            )
            .is_err()
    );
    assert!(
        initial
            .at_query_point(
                &card,
                &point(320.0, 114_688.0).with("shear_rate", 1.0).unwrap()
            )
            .is_err()
    );
    let expanded = initial
        .liquid()
        .clone()
        .with_density_expansion(&card, "T", Temperature::new(280.0), Temperature::new(360.0))
        .unwrap();
    assert!(ResolvedLiquidDensityState::try_from_liquid(&card, expanded, "T", "p").is_err());
    let independent = card_with_density(
        density_surface(-0.5, RHO_P, false),
        InterpolationPolicy::LinearInside,
        validity([256.0, 384.0], [65_536.0, 196_608.0]),
        LiquidViscosityBasis::Dynamic,
        Some(1e-4),
    );
    let expanded = resolve_liquid_state_with_expansion(
        &independent,
        &point(320.0, 114_688.0),
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(ResolvedLiquidDensityState::try_from_liquid(&independent, expanded, "T", "p").is_err());
    assert_eq!(initial.material().card_identity(), card.content_hash());
    assert_eq!(initial.response().temperature_k(), 320.0);
}

#[test]
fn isolated_values_and_discrete_only_sources_do_not_become_density_laws() {
    let scalar = PropertyValue::Scalar {
        value: 1000.0,
        dims: Density::DIMS,
    };
    let card = card_with_density(
        scalar,
        InterpolationPolicy::ConstantWithinValidity,
        validity([256.0, 384.0], [65_536.0, 196_608.0]),
        LiquidViscosityBasis::Dynamic,
        None,
    );
    assert!(
        ResolvedLiquidDensityState::try_from_liquid(
            &card,
            liquid(
                &card,
                &point(320.0, 114_688.0),
                LiquidViscosityBasis::Dynamic
            ),
            "T",
            "p"
        )
        .is_err()
    );
    let card = card_with_density(
        density_surface(-0.5, RHO_P, false),
        InterpolationPolicy::TabulatedOnly,
        validity([256.0, 384.0], [65_536.0, 196_608.0]),
        LiquidViscosityBasis::Dynamic,
        None,
    );
    assert!(
        ResolvedLiquidDensityState::try_from_liquid(
            &card,
            liquid(
                &card,
                &point(256.0, 65_536.0),
                LiquidViscosityBasis::Dynamic
            ),
            "T",
            "p"
        )
        .is_err()
    );
}

#[test]
fn density_coordinates_require_explicit_temperature_and_pressure_meanings() {
    let dimensions = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(Temperature::DIMS),
            320.0,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            100_000.0,
        )
        .unwrap();
    assert!(validate_coordinates(&dimensions, "temperature", "pressure").is_ok());
    for temperature in [
        QuantitySpec::dimensional(Pressure::DIMS),
        typed(QuantityKind::TemperatureDifference),
    ] {
        let point = QueryPoint::new()
            .with_quantity("T", temperature, 320.0)
            .unwrap()
            .with_quantity("p", typed(QuantityKind::Pressure), 100_000.0)
            .unwrap();
        assert!(validate_coordinates(&point, "T", "p").is_err());
    }
    let untyped = QueryPoint::new()
        .with("T", 320.0)
        .unwrap()
        .with("p", 100_000.0)
        .unwrap();
    assert!(validate_coordinates(&untyped, "T", "p").is_err());
    assert!(validate_coordinates(&dimensions, "temperature", "temperature").is_err());
}
