//! G1/G3/G4: actual interface temperatures, pinned contact laws and heat balance.
#[path = "contact/generation.rs"]
mod generation;

use super::{close, curve, options, temperature, with_cancelled_cx, with_cx};
use fs_conduction::interface::temperature::TemperatureContactCoordinates;
use fs_conduction::layered::{
    LayeredTemperatureContact, LayeredWall, PlanarThermalLayer, TemperatureContactWallOptions,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionError, ConductivityTable, TEMPERATURE_DIMS,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};

fn contact_options() -> TemperatureContactWallOptions {
    TemperatureContactWallOptions {
        solver: options(),
        max_contact_iterations: 128,
        relaxation: 0.5,
    }
}

fn contact_card(
    axis: &str,
    knots: Vec<(f64, f64)>,
    validity: ValidityDomain,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: TEMPERATURE_DIMS,
                knots,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: validity.with("pressure", 100_000.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured contact law; not empirical material data".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        texture_frame: "declared planar interface".into(),
    };
    InterfaceSystemCard::assemble(
        surface("left"),
        surface("right"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "declared fixture".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn mean_point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
}

fn mean_card(low: f64, high: f64) -> InterfaceSystemCard {
    contact_card(
        "T",
        vec![(low, 0.02 * (low - 290.0)), (high, 0.02 * (high - 290.0))],
        ValidityDomain::unconstrained().with("T", low, high),
    )
}

fn binding(card: &InterfaceSystemCard, index: usize) -> LayeredTemperatureContact<'_> {
    LayeredTemperatureContact::new(
        index,
        "temperature-dependent joint",
        card,
        &mean_point(310.0),
        TemperatureContactCoordinates::Mean("T".into()),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap()
}

fn wall() -> LayeredWall {
    LayeredWall::new(
        vec![
            curve(
                "nonlinear hot solid",
                1.0,
                vec![(300.0, 0.75), (450.0, 1.125)],
            ),
            PlanarThermalLayer::new("middle", 1.0, ConductivityTable::declared(1.0).unwrap())
                .unwrap(),
            PlanarThermalLayer::new("cold", 1.0, ConductivityTable::declared(1.0).unwrap())
                .unwrap(),
        ],
        vec![0.4, 0.5],
        [0.0; 2],
    )
    .unwrap()
}

#[test]
fn g1_source_contacts_close_actual_interface_heat_with_nonlinear_bulk_and_fixed_contacts() {
    let card = mean_card(300.0, 450.0);
    let contacts = [binding(&card, 0)];
    let wall = wall();
    let original = wall.clone();
    // Construct q=20 independently from the cold end:
    // right solid 300->320; fixed contact ->330; middle solid ->350;
    // R''(mean)=0.02*(mean-290) gives jump 30, mean 365 and R''=1.5;
    // integral_380^hot (T/400) dT = q L =20 for the hot nonlinear layer.
    let hot = (380.0_f64.powi(2) + 16_000.0).sqrt();
    let answer = with_cx(|cx| {
        wall.solve_temperature_dependent_with_contacts(
            cx,
            [hot, 300.0],
            2.0,
            &contacts,
            contact_options(),
        )
    })
    .unwrap();
    assert!((answer.conduction.heat_flux_w_per_m2 - 20.0).abs() < 1e-7);
    assert!((answer.conduction.heat_rate_w - 40.0).abs() < 2e-7);
    for (layer, expected) in answer.conduction.layer_temperatures.iter().zip([
        [hot, 380.0],
        [350.0, 330.0],
        [320.0, 300.0],
    ]) {
        close(layer.left_k, expected[0]);
        close(layer.right_k, expected[1]);
    }
    let state = &answer.contacts[0];
    close(state.physical.value_m2_k_per_w(), 1.5);
    close(state.temperatures_k[0], 380.0);
    close(state.temperatures_k[1], 350.0);
    assert!(answer.maximum_flux_residual_w_per_m2 <= 1.1e-8);
    assert!(state.flux_residual_w_per_m2.abs() <= 1.1e-8);
    assert_eq!(
        state.physical.receipt().selected,
        contacts[0].reference_resistance().receipt().selected
    );
    card.claims()
        .verify_receipt(state.applied.receipt())
        .unwrap();
    card.claims()
        .verify_receipt(state.physical.receipt())
        .unwrap();
    let query = &state.physical.receipt().query_point;
    close(query.iter().find(|(axis, _)| axis == "T").unwrap().1, 365.0);
    assert!(query.contains(&("pressure".into(), 100_000.0)));
    close(
        answer.conduction.wall.contact_resistances_m2_k_per_w()[1],
        0.5,
    );
    assert!(answer.contact_iterations > 1);
    assert!(answer.total_shooting_iterations >= answer.conduction.iterations);
    let frozen =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [hot, 300.0], 2.0, options())).unwrap();
    assert!(
        (frozen.heat_flux_w_per_m2 - 20.0).abs() > 4.0,
        "the reference contact value must not stand in for the physical interface law"
    );
    assert_eq!(wall, original);

    let constant_card = contact_card(
        "T",
        vec![(300.0, 0.4), (450.0, 0.4)],
        ValidityDomain::unconstrained().with("T", 300.0, 450.0),
    );
    let constant = with_cx(|cx| {
        wall.solve_temperature_dependent_with_contacts(
            cx,
            [hot, 300.0],
            2.0,
            &[binding(&constant_card, 0)],
            contact_options(),
        )
    })
    .unwrap();
    assert_eq!(constant.conduction, frozen);
    assert_eq!(constant.contact_iterations, 1);
}

#[test]
fn g1_zero_heat_still_resolves_the_actual_local_resistance_and_receipt() {
    let card = mean_card(300.0, 450.0);
    let contacts = [binding(&card, 0)];
    let answer = with_cx(|cx| {
        wall().solve_temperature_dependent_with_contacts(
            cx,
            [355.0; 2],
            3.0,
            &contacts,
            contact_options(),
        )
    })
    .unwrap();
    assert_eq!(answer.conduction.heat_rate_w, 0.0);
    close(answer.contacts[0].physical.value_m2_k_per_w(), 1.3);
    close(answer.contacts[0].applied.value_m2_k_per_w(), 1.3);
    close(
        answer.conduction.effective_resistance_m2_k_per_w,
        400.0 / 355.0 + 2.0 + 0.5 + 1.3,
    );
    assert_eq!(answer.contact_iterations, 1);
    assert_eq!(answer.total_shooting_iterations, 0);
    assert!(
        answer.contacts[0]
            .physical
            .receipt()
            .query_point
            .contains(&("T".into(), 355.0))
    );
}

#[test]
fn g1_g3_two_sided_source_keeps_physical_left_right_order_when_heat_reverses() {
    let card = contact_card(
        "Ta",
        vec![(300.0, 0.4), (400.0, 1.4)],
        ValidityDomain::unconstrained()
            .with_quantity("Ta", temperature(), 300.0, 400.0)
            .with_quantity("Tb", temperature(), 300.0, 400.0),
    );
    let point = QueryPoint::new()
        .with_quantity("Ta", temperature(), 310.0)
        .unwrap()
        .with_quantity("Tb", temperature(), 330.0)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap();
    let contacts = [LayeredTemperatureContact::new(
        0,
        "ordered source",
        &card,
        &point,
        TemperatureContactCoordinates::TwoSided(["Ta".into(), "Tb".into()]),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap()];
    let wall = LayeredWall::new(
        vec![
            PlanarThermalLayer::new("left", 1.0, ConductivityTable::declared(1.0).unwrap())
                .unwrap(),
            PlanarThermalLayer::new("right", 2.0, ConductivityTable::declared(1.0).unwrap())
                .unwrap(),
        ],
        vec![0.5],
        [0.0; 2],
    )
    .unwrap();
    // Forward q=20 gives Ta=360, Tb=340 and R''(Ta)=1.
    // Reverse magnitude m obeys m*(3.4+0.01*m)=80 because Ta=300+m.
    let reverse_magnitude = 160.0 / (3.4 + 14.76_f64.sqrt());
    for (reservoirs, expected) in [([380.0, 300.0], 20.0), ([300.0, 380.0], -reverse_magnitude)] {
        let answer = with_cx(|cx| {
            wall.solve_temperature_dependent_with_contacts(
                cx,
                reservoirs,
                1.0,
                &contacts,
                contact_options(),
            )
        })
        .unwrap();
        assert!((answer.conduction.heat_flux_w_per_m2 - expected).abs() < 1e-7);
        let state = &answer.contacts[0];
        for (axis, value) in ["Ta", "Tb"].into_iter().zip(state.temperatures_k) {
            close(
                state
                    .physical
                    .receipt()
                    .query_point
                    .iter()
                    .find(|(name, _)| name == axis)
                    .unwrap()
                    .1,
                value,
            );
        }
        card.claims()
            .verify_receipt(state.physical.receipt())
            .unwrap();
    }
}

#[test]
fn g4_source_domain_work_budget_binding_and_cancellation_refuse_without_changing_wall() {
    let wall = wall();
    let original = wall.clone();
    let card = mean_card(300.0, 450.0);
    let contacts = [binding(&card, 0)];
    let hot = (380.0_f64.powi(2) + 16_000.0).sqrt();
    let mut exhausted = contact_options();
    exhausted.max_contact_iterations = 1;
    let error = with_cx(|cx| {
        wall.solve_temperature_dependent_with_contacts(cx, [hot, 300.0], 1.0, &contacts, exhausted)
    })
    .unwrap_err();
    assert!(error.to_string().contains("exhausted"));
    let narrow = mean_card(300.0, 330.0);
    assert!(
        with_cx(|cx| wall.solve_temperature_dependent_with_contacts(
            cx,
            [355.0; 2],
            1.0,
            &[binding(&narrow, 0)],
            contact_options(),
        ))
        .is_err()
    );
    assert!(
        with_cx(|cx| wall.solve_temperature_dependent_with_contacts(
            cx,
            [hot, 300.0],
            1.0,
            &[binding(&narrow, 0)],
            contact_options(),
        ))
        .is_err()
    );
    for invalid in [
        vec![binding(&card, 4)],
        vec![binding(&card, 0), binding(&card, 0)],
    ] {
        assert!(
            with_cx(|cx| wall.solve_temperature_dependent_with_contacts(
                cx,
                [hot, 300.0],
                1.0,
                &invalid,
                contact_options(),
            ))
            .is_err()
        );
    }
    for relaxation in [0.0, -1.0, 1.1, f64::NAN] {
        let mut invalid = contact_options();
        invalid.relaxation = relaxation;
        assert!(
            with_cx(|cx| wall.solve_temperature_dependent_with_contacts(
                cx,
                [hot, 300.0],
                1.0,
                &contacts,
                invalid,
            ))
            .is_err()
        );
    }
    assert!(matches!(
        with_cancelled_cx(|cx| wall.solve_temperature_dependent_with_contacts(
            cx,
            [hot, 300.0],
            1.0,
            &contacts,
            contact_options(),
        )),
        Err(ConductionError::Cancelled { .. })
    ));
    assert_eq!(wall, original);
}
