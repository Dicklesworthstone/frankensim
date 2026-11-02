//! G1/G3/G4: real source-temperature continuation retains loaded shear memory.
use super::*;
use fs_contact::interface_binding::mindlin::temperature::SourceMindlinTemperatureError;
use fs_contact::normal_patch::source::state::PointContactMaterialCards;
use fs_matdb::{ClaimSet, MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

struct Fixture {
    bulk: [MaterialCard; 2],
    interface: InterfaceSystemCard,
}
impl Fixture {
    fn new() -> Self {
        let bulk = [
            common::solid("thermal-Mindlin-A", 1000.0, 1e7),
            common::solid("thermal-Mindlin-B", 2000.0, 2e7),
        ];
        let normal = common::bind(&bulk[0], &bulk[1], None, 0.0);
        let original = normal.normal_state().elastic().interface().card();
        let mut claims = original.claims().clone();
        claims.insert_claim(claim(MU, 0.6, Dims::NONE)).unwrap();
        let interface = InterfaceSystemCard::assemble(
            original.surface_a().clone(),
            original.surface_b().clone(),
            original.context().clone(),
            claims,
            original.models().to_vec(),
        )
        .unwrap();
        Self { bulk, interface }
    }
    fn at(&self, temperatures: [f64; 3]) -> Result<SourceMindlinContact, String> {
        let points = temperatures.map(point);
        self.at_points(
            &points,
            2e-5,
            MaterialPropertySelection::PinnedByProperty(vec![(
                MU.into(),
                self.interface.claims().claims_for(MU)[0].0,
            )]),
        )
    }
    fn at_points(
        &self,
        points: &[QueryPoint; 3],
        indentation: f64,
        selection: MaterialPropertySelection,
    ) -> Result<SourceMindlinContact, String> {
        let normal = PointContactMaterialCards {
            surface_a: &self.bulk[0],
            surface_b: &self.bulk[1],
            interface: &self.interface,
            bulk_selection: [
                MaterialPropertySelection::SingleClaimOnly,
                MaterialPropertySelection::SingleClaimOnly,
            ],
            interface_selection: MaterialPropertySelection::SingleClaimOnly,
            model_selection: NormalContactModelSelection::SingleSupported,
            authority: InputAuthority::SyntheticFixture,
        }
        .resolve_at_temperatures(
            [&points[0], &points[1]],
            &points[2],
            PointContactGeometry::SpherePlane {
                effective_radius_m: 0.1,
            },
            common::envelope(),
        )
        .map_err(|e| format!("{e}"))?;
        SourceMindlinContact::resolve(
            normal,
            indentation,
            ScalarPropertyRequirement::try_new(MU, Dims::NONE, ScalarAdmissibility::NonNegative)
                .unwrap(),
            selection,
            64,
        )
        .map_err(|e| format!("{e}"))
    }
}
fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermal Mindlin source".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn point(temperature: f64) -> QueryPoint {
    common::query()
        .with("T", temperature)
        .unwrap()
        .with_quantity(
            "ambient-pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
        )
        .unwrap()
}
fn loaded(
    law: &SourceMindlinContact,
) -> fs_contact::interface_binding::mindlin::SourceMindlinState {
    let critical = law.law().gross_slip_displacement_m();
    let mut state = law.initial_state();
    for fraction in [0.7, -0.3, 0.1] {
        state = law
            .advance(&state, &step(law, fraction * critical))
            .unwrap()
            .state;
    }
    assert!(state.history_intervals() > 1);
    assert!(law.stored_energy_j(&state) > 0.0);
    state
}

#[test]
fn independently_heated_sources_continue_exact_loaded_and_reversal_histories() {
    let f = Fixture::new();
    let reference = f.at([293.15; 3]).unwrap();
    let mut reference_state = loaded(&reference);
    let mut current = reference.clone();
    let mut state = reference_state.clone();
    assert_eq!(
        current
            .continue_at_temperature(&state, &current, &mut || false)
            .unwrap(),
        state
    );
    for (temperatures, fraction) in [
        ([300.0, 310.0, 305.0], -0.4),
        ([315.0, 297.0, 307.0], 0.2),
        ([293.15; 3], -0.1),
    ] {
        let refreshed = f.at(temperatures).unwrap();
        let before = state.clone();
        let transferred = current
            .continue_at_temperature(&state, &refreshed, &mut || false)
            .unwrap();
        assert_eq!(transferred.constitutive_state(), state.constitutive_state());
        assert_eq!(
            refreshed.stored_energy_j(&transferred),
            current.stored_energy_j(&state)
        );
        assert_eq!(transferred.endpoint_force_n(), state.endpoint_force_n());
        assert_eq!(transferred.slip_m(), state.slip_m());
        assert_eq!(transferred.history_intervals(), state.history_intervals());
        assert!(current.stored_energy_j(&transferred).is_nan());
        assert!(refreshed.stored_energy_j(&state).is_nan());
        let solids = refreshed.normal_source().model().normal_state().elastic();
        assert_eq!(
            solids.surface_a_state().resolved().state_coordinate("T"),
            Some(temperatures[0])
        );
        assert_eq!(
            solids.surface_b_state().resolved().state_coordinate("T"),
            Some(temperatures[1])
        );
        assert_eq!(refreshed.temperature_kelvin(), temperatures[2]);
        let du = fraction * reference.law().gross_slip_displacement_m();
        let expected = reference
            .advance(&reference_state, &step(&reference, du))
            .unwrap();
        let continued = refreshed
            .advance(&transferred, &step(&refreshed, du))
            .unwrap();
        assert_eq!(
            continued.state.constitutive_state(),
            expected.state.constitutive_state()
        );
        assert_eq!(continued.traction_n, expected.traction_n);
        assert_eq!(continued.stored_change_j, expected.stored_change_j);
        assert_eq!(continued.dissipated_j, expected.dissipated_j);
        assert_eq!(state, before);
        reference_state = expected.state;
        state = continued.state;
        current = refreshed;
    }
    assert_eq!(current, reference);
    assert_eq!(state, reference_state);
}

#[test]
fn equal_endpoint_values_do_not_hide_temperature_dependent_friction_or_elasticity() {
    for bulk_varies in [false, true] {
        let mut f = Fixture::new();
        let (name, value, dims) = if bulk_varies {
            ("young_modulus", 1e7, Pressure::DIMS)
        } else {
            (MU, 0.6, Dims::NONE)
        };
        let mut varying = claim(name, value, dims);
        varying.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (280.0, value),
                (293.15, value),
                (300.0, 1.3 * value),
                (310.0, value),
                (320.0, value),
            ],
            dims,
        };
        varying.interpolation = InterpolationPolicy::LinearInside;
        let original = if bulk_varies {
            f.bulk[0].claims()
        } else {
            f.interface.claims()
        };
        let mut claims = ClaimSet::new();
        for (_, property) in original
            .claims_ordered()
            .filter(|(_, p)| p.key.name() != name)
        {
            claims.insert_claim(property.clone()).unwrap();
        }
        claims.insert_claim(varying).unwrap();
        if bulk_varies {
            f.bulk[0] = MaterialCard::assemble(f.bulk[0].id().clone(), claims, vec![]).unwrap();
        } else {
            f.interface = InterfaceSystemCard::assemble(
                f.interface.surface_a().clone(),
                f.interface.surface_b().clone(),
                f.interface.context().clone(),
                claims,
                f.interface.models().to_vec(),
            )
            .unwrap();
        }
        let before = f.at([293.15; 3]).unwrap();
        let temperatures = if bulk_varies {
            [310.0, 293.15, 293.15]
        } else {
            [293.15, 293.15, 310.0]
        };
        let after = f.at(temperatures).unwrap();
        assert_eq!(before.law(), after.law());
        let state = loaded(&before);
        // No temperature change is an exact continuation even for a curve.
        assert_eq!(
            before
                .continue_at_temperature(&state, &before, &mut || false)
                .unwrap(),
            state
        );
        assert!(
            before
                .continue_at_temperature(&state, &after, &mut || false)
                .is_err()
        );
        assert_eq!(state, loaded(&before));
    }
}

#[test]
fn original_cards_selected_claims_geometry_and_nonthermal_conditions_cannot_change() {
    let mut f = Fixture::new();
    // Equal competing coefficient claims stay distinct within the SAME full card.
    let mut competing = claim(MU, 0.6, Dims::NONE);
    competing.provenance.source = "different synthetic friction observation".into();
    let mut claims = f.interface.claims().clone();
    let other_id = claims.insert_claim(competing).unwrap();
    f.interface = InterfaceSystemCard::assemble(
        f.interface.surface_a().clone(),
        f.interface.surface_b().clone(),
        f.interface.context().clone(),
        claims,
        f.interface.models().to_vec(),
    )
    .unwrap();
    let original = f.at([293.15; 3]).unwrap();
    let state = loaded(&original);
    let points = [300.0, 301.0, 302.0].map(point);
    let other_claim = f
        .at_points(
            &points,
            2e-5,
            MaterialPropertySelection::PinnedByProperty(vec![(MU.into(), other_id)]),
        )
        .unwrap();
    assert_eq!(original.law(), other_claim.law());
    assert!(
        original
            .continue_at_temperature(&state, &other_claim, &mut || false)
            .is_err()
    );
    let pin = MaterialPropertySelection::PinnedByProperty(vec![(
        MU.into(),
        original.friction_source().properties()[0]
            .answer()
            .receipt
            .selected,
    )]);
    let changed_load = f.at_points(&points, 3e-5, pin.clone()).unwrap();
    assert!(
        original
            .continue_at_temperature(&state, &changed_load, &mut || false)
            .is_err()
    );
    let pressure = points
        .each_ref()
        .map(|p| p.clone().with("ambient-pressure", 2e5).unwrap());
    let pressure = f.at_points(&pressure, 2e-5, pin.clone()).unwrap();
    assert_eq!(original.law(), pressure.law());
    assert!(
        original
            .continue_at_temperature(&state, &pressure, &mut || false)
            .is_err()
    );
    let schema = points.each_ref().map(|p| {
        p.clone()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                p.axes()["T"],
            )
            .unwrap()
    });
    // The original source refuses changing an untyped coordinate to typed T
    // before a candidate contact can exist; continuation does not bypass it.
    assert!(f.at_points(&schema, 2e-5, pin).is_err());
    // An unrelated full-card edit cannot masquerade as temperature continuation.
    let mut revised_bulk = f.bulk[0].claims().clone();
    revised_bulk
        .insert_claim(claim("unrelated-original-datum", 4.0, Dims::NONE))
        .unwrap();
    f.bulk[0] = MaterialCard::assemble(f.bulk[0].id().clone(), revised_bulk, vec![]).unwrap();
    let revised = f.at([300.0; 3]).unwrap();
    assert_eq!(original.law(), revised.law());
    assert!(
        original
            .continue_at_temperature(&state, &revised, &mut || false)
            .is_err()
    );
    assert_eq!(state, loaded(&original));
}

#[test]
fn source_bounds_foreign_history_and_final_cancellation_preserve_exact_retry() {
    let f = Fixture::new();
    let original = f.at([293.15; 3]).unwrap();
    let refreshed = f.at([300.0, 310.0, 305.0]).unwrap();
    let state = loaded(&original);
    let saved = state.clone();
    for temperatures in [
        [321.0, 300.0, 300.0],
        [300.0, 279.0, 300.0],
        [300.0, 300.0, 321.0],
    ] {
        assert!(f.at(temperatures).is_err());
    }
    let mut polls = 0;
    let expected = original
        .continue_at_temperature(&state, &refreshed, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 2);
    let mut calls = 0;
    let stopped = original.continue_at_temperature(&state, &refreshed, &mut || {
        calls += 1;
        calls == polls
    });
    assert_eq!(
        stopped.unwrap_err(),
        SourceMindlinTemperatureError::Cancelled
    );
    assert_eq!(calls, polls);
    assert!(
        refreshed
            .continue_at_temperature(&state, &refreshed, &mut || false)
            .is_err()
    );
    assert!(
        original
            .continue_at_temperature(&expected, &refreshed, &mut || false)
            .is_err()
    );
    assert_eq!(state, saved);
    assert_eq!(
        original
            .continue_at_temperature(&state, &refreshed, &mut || false)
            .unwrap(),
        expected
    );
    assert_eq!(
        refreshed
            .continue_at_temperature(&expected, &refreshed, &mut || false)
            .unwrap(),
        expected
    );
}
