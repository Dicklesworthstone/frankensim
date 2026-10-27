//! Original resistance law with independent thermal states, not measured contact qualification.
#[path = "contact_material_state/common.rs"]
mod common;
use common::*;
mod support;
use fs_conduction::interface::lumped_pair::{MaterialContactLaw, advance_contact_pair};
use fs_conduction::interface::temperature::ContactIterationConfig;
use fs_conduction::lumped::LumpedThermalTransport;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_conduction::InterfaceResistance;
use fs_matdb::{ClaimSet, SelectionPolicy};
use fs_qty::{Dims, QuantitySpec, Temperature};

#[test]
fn phase_changes_resistance_at_fixed_kelvin_and_retains_complete_source_point() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let (source,id)=source(&cards,1.0); let context=point(); let binding=binding(&cards,&source,&context);
    for (h,expected) in [(25.0,1.25),(75.0,1.75)] {
        let states=[curves[0].state_at_specific_enthalpy(h).unwrap(),
            curves[1].state_at_specific_enthalpy(-50.0).unwrap()];
        let point=binding.at(states).unwrap();
        assert_eq!(point.axes()["aT"],400.0); assert_eq!(point.axes()["bT"],350.0);
        assert_eq!(point.axes()["T"],375.0); assert_eq!(point.axes()["pressure"],200000.0);
        assert_eq!(point.axes()["lot"],7.0); assert_eq!(point.axes()["fa"],h/100.0);
        assert_eq!(point.axis_quantities()["aT"],QuantitySpec::dimensional(Temperature::DIMS));
        let resistance=InterfaceResistance::from_card_pinned("contact",&source,&point,id).unwrap();
        assert_eq!(resistance.value_m2_k_per_w(),expected);
        assert_eq!(resistance,InterfaceResistance::from_card_pinned("contact",&source,&point,id).unwrap());
        assert_eq!(resistance.card_identity(),source.content_hash());
    }
    assert_eq!(context.axes()["fa"],0.0);
}

#[test]
fn ordered_surface_and_chart_identity_and_axis_aliasing_are_not_silently_repaired() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let (source,_)=source(&cards,1.0); let context=point(); let mut binding=binding(&cards,&source,&context);
    let states=[curves[0].state_at_specific_enthalpy(25.0).unwrap(),curves[1].state_at_specific_enthalpy(-50.0).unwrap()];
    assert!(binding.at([states[1],states[0]]).is_err());
    binding.materials.swap(0,1); assert!(binding.at([states[1],states[0]]).is_err());
    binding.materials.swap(0,1); binding.temperature_axes=["aT","aT"];
    assert!(binding.at(states).is_err()); binding.temperature_axes=["aT","bT"];
    binding.mean_temperature_axis=Some("aT"); assert!(binding.at(states).is_err());
    binding.mean_temperature_axis=Some("T"); binding.liquid_fraction_axes=[Some("fa"),Some("fa")];
    assert!(binding.at(states).is_err());
    let bad=fs_matdb::QueryPoint::new().with_quantity("aT",QuantitySpec::dimensional(Dims::NONE),400.0).unwrap();
    binding.context=&bad; binding.liquid_fraction_axes=[Some("fa"),Some("fb")];
    assert!(binding.at(states).is_err());
}

#[test]
fn pins_domains_and_static_source_support_remain_load_bearing() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let states=[curves[0].state_at_specific_enthalpy(25.0).unwrap(),curves[1].state_at_specific_enthalpy(-50.0).unwrap()];
    let mut claims=ClaimSet::new(); let first=claims.insert_claim(claim("fa",Dims::NONE,1.0)).unwrap();
    claims.insert_claim(claim("fa",Dims::NONE,2.0)).unwrap();
    let source=interface(&cards,claims); let context=point(); let b=binding(&cards,&source,&context);
    let query=b.at(states).unwrap();
    assert!(InterfaceResistance::from_card("contact",&source,&query,SelectionPolicy::SingleClaimOnly).is_err());
    assert_eq!(InterfaceResistance::from_card_pinned("contact",&source,&query,first)
        .unwrap().value_m2_k_per_w(),1.25);
    let wrong_pressure=context.clone().with("pressure",100000.0).unwrap();
    let wrong_query=binding(&cards,&source,&wrong_pressure).at(states).unwrap();
    assert!(InterfaceResistance::from_card_pinned("contact",&source,&wrong_query,first).is_err());
    assert!(b.is_thermal_axis("fa",Dims::NONE));
    assert!(b.is_thermal_axis("aT",Temperature::DIMS));
    assert!(b.is_thermal_axis("T",Temperature::DIMS));
    assert!(!b.is_thermal_axis("frequency",Dims::NONE));
    assert!(!b.is_thermal_axis("aT",Dims::NONE));
}

fn contact_controls() -> ContactIterationConfig {
    ContactIterationConfig {
        max_iterations: 64,
        relaxation: 1.0,
        residual_tolerance_j: 1e-13,
    }
}

fn caloric(cards: &[fs_matdb::MaterialCard; 2]) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    std::array::from_fn(|side| {
        let cp = [0.02, 0.04][side];
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            cards[side].content_hash(),
            SolidLiquidPhase::Solid,
            [300.0, 320.0]
                .map(|temperature_k| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: cp * (temperature_k - 300.0),
                    temperature_k,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                })
                .to_vec(),
        )
        .unwrap()
    })
}

fn thermal_input(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 0.001,
        volume_m3: 1e-6,
        surface_area_m2: 0.0006,
        internal_heat_j: 0.0,
        duration_s: 0.001,
    }
}

fn temperature_claim(slope: f64, low: f64, high: f64) -> fs_matdb::PropertyClaim {
    let mut source = claim("T", Temperature::DIMS, 1e-4);
    source.value = fs_matdb::PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1e-4), (320.0, 1e-4 * (1.0 + 20.0 * slope))],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    };
    source.validity =
        source
            .validity
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), low, high);
    source
}

fn temperature_interface(
    cards: &[fs_matdb::MaterialCard; 2],
    slope: f64,
) -> (fs_matdb::InterfaceSystemCard, fs_matdb::ClaimId) {
    let mut claims = ClaimSet::new();
    let id = claims
        .insert_claim(temperature_claim(slope, 300.0, 320.0))
        .unwrap();
    (interface(cards, claims), id)
}

#[test]
fn g1_material_contact_closes_the_endpoint_law_against_an_independent_quadratic() {
    let cards = [bulk("A"), bulk("B")];
    let curves = caloric(&cards);
    let inputs = [
        thermal_input(&curves[0], 0.2),
        thermal_input(&curves[1], 0.0),
    ];
    let (source, id) = temperature_interface(&cards, 0.1);
    let context = point();
    let b = binding(&cards, &source, &context);
    let initial = [inputs[0].initial, inputs[1].initial];
    let law =
        MaterialContactLaw::new(b, initial, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let reference = law.reference().clone();
    let frozen = InterfaceResistance::from_card_pinned(
        "material contact",
        &source,
        &b.at(initial).unwrap(),
        id,
    )
    .unwrap();
    let transport = LumpedThermalTransport::try_declared(100.0, 0.0).unwrap();
    support::with_cx(|cx| {
        let result = law
            .advance(cx, inputs, [&transport; 2], 1e-6, 1e-14, contact_controls())
            .unwrap();
        // Independent root of 125000 Q^2 - 225 Q + 0.01 = 0.
        let expected = 0.02 / (225.0 + 45_625.0_f64.sqrt());
        let q = result[1].external_heat_j;
        assert!((q - expected).abs() < 2e-13);
        assert_eq!(result[0].external_heat_j, -q);
        let ta = 310.0 - q / (0.001 * 0.02);
        let tb = 300.0 + q / (0.001 * 0.04);
        let actual_r = 1e-4 * (1.0 + 0.1 * (0.5 * (ta + tb) - 300.0));
        assert!((q - 0.001 * (ta - tb) / (actual_r / 1e-6)).abs() < 1e-13);
        let report = &result[0].report;
        assert!(report.iterations > 1);
        assert!(report.physical_contact_residual_j.abs() <= 1e-13);
        assert_eq!(
            report.endpoint.state_coordinate("aT"),
            Some(result[0].state.temperature_k())
        );
        assert_eq!(
            report.endpoint.state_coordinate("bT"),
            Some(result[1].state.temperature_k())
        );
        assert_eq!(report.endpoint.state_coordinate("pressure"), Some(200000.0));
        assert_eq!(report.endpoint.state_coordinate("lot"), Some(7.0));
        assert_eq!(
            report.endpoint.properties()[0].answer().receipt.selected,
            id
        );
        let held = advance_contact_pair(cx, inputs, [&transport; 2], &frozen, 1e-6, 1e-14).unwrap();
        assert!(q > 1.01 * held[1].external_heat_j);
        let replay = law
            .advance(cx, inputs, [&transport; 2], 1e-6, 1e-14, contact_controls())
            .unwrap();
        for side in 0..2 {
            assert_eq!(result[side].state, replay[side].state);
            assert_eq!(result[side].external_heat_j, replay[side].external_heat_j);
            assert_eq!(result[side].report, replay[side].report);
        }
    });
    assert_eq!(law.reference(), &reference);
}

#[test]
fn g3_constant_material_contact_retains_native_exchange_and_open_link_heating() {
    let cards = [bulk("A"), bulk("B")];
    let curves = caloric(&cards);
    let inputs = [
        thermal_input(&curves[0], 0.2),
        thermal_input(&curves[1], 0.0),
    ];
    let initial = [inputs[0].initial, inputs[1].initial];
    let (source, id) = temperature_interface(&cards, 0.0);
    let context = point();
    let b = binding(&cards, &source, &context);
    let law =
        MaterialContactLaw::new(b, initial, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let frozen = InterfaceResistance::from_card_pinned(
        "material contact",
        &source,
        &b.at(initial).unwrap(),
        id,
    )
    .unwrap();
    let transport = LumpedThermalTransport::try_declared(100.0, 0.0).unwrap();
    support::with_cx(|cx| {
        let result = law
            .advance(cx, inputs, [&transport; 2], 1e-6, 1e-14, contact_controls())
            .unwrap();
        let held = advance_contact_pair(cx, inputs, [&transport; 2], &frozen, 1e-6, 1e-14).unwrap();
        assert_eq!(result[0].report.iterations, 1);
        for side in 0..2 {
            assert_eq!(result[side].state, held[side].state);
            assert_eq!(result[side].external_heat_j, held[side].external_heat_j);
        }
        let mut separated = inputs;
        separated[1].internal_heat_j = 4e-5;
        let open = law
            .advance(
                cx,
                separated,
                [&transport; 2],
                0.0,
                1e-14,
                contact_controls(),
            )
            .unwrap();
        assert_eq!(open[0].state, inputs[0].initial);
        assert!((open[1].state.temperature_k() - 301.0).abs() < 1e-12);
        assert_eq!(open[0].external_heat_j, 0.0);
        assert_eq!(open[1].external_heat_j, 0.0);
        assert_eq!(
            open[0].report.endpoint.state_coordinate("bT"),
            Some(open[1].state.temperature_k())
        );
        assert_eq!(open[0].report.physical_contact_residual_j, 0.0);
        assert_eq!(open[0].report.endpoint_maximum_biot, 0.0);
        assert!(open[0].report.applied.contact.is_none());
    });
}

#[test]
fn g1_latent_fraction_refreshes_contact_at_unchanged_temperature() {
    let cards = [bulk("A"), bulk("B")];
    let curves = [curve(&cards[0]), curve(&cards[1])];
    let mut inputs = [
        thermal_input(&curves[0], 0.0),
        thermal_input(&curves[1], 300.0),
    ];
    // These original phase charts declare density 1 kg/m3.
    for input in &mut inputs {
        input.volume_m3 = 0.001;
        input.surface_area_m2 = 0.06;
    }
    let (source, id) = source(&cards, 0.001);
    let context = point();
    let law = MaterialContactLaw::new(
        binding(&cards, &source, &context),
        [inputs[0].initial, inputs[1].initial],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let transport = LumpedThermalTransport::try_declared(100.0, 0.0).unwrap();
    support::with_cx(|cx| {
        let result = law
            .advance(cx, inputs, [&transport; 2], 1e-4, 1e-14, contact_controls())
            .unwrap();
        // Independent root of 100 Q^2 + 11 Q - 0.2 = 0, Q into A.
        let expected = 0.4 / (11.0 + 201.0_f64.sqrt());
        let heat_into_a = result[0].external_heat_j;
        assert!((heat_into_a - expected).abs() < 2e-13);
        assert_eq!(result[1].external_heat_j, -heat_into_a);
        assert_eq!(result[0].state.temperature_k(), 400.0);
        let fraction = heat_into_a / 0.1;
        let report = &result[0].report;
        assert!((report.endpoint.state_coordinate("fa").unwrap() - fraction).abs() < 2e-12);
        assert!(
            (report.endpoint.properties()[0].value_si() - 0.001 * (1.0 + fraction)).abs() < 1e-15
        );
        assert_eq!(
            report.endpoint.properties()[0].answer().receipt.selected,
            id
        );
        assert!(report.iterations > 1);
        assert!(report.physical_contact_residual_j.abs() <= 1e-13);
    });
}

#[test]
fn g4_material_contact_refuses_source_switches_work_exhaustion_and_cancellation() {
    let cards = [bulk("A"), bulk("B")];
    let curves = caloric(&cards);
    let inputs = [
        thermal_input(&curves[0], 0.2),
        thermal_input(&curves[1], 0.0),
    ];
    let initial = [inputs[0].initial, inputs[1].initial];
    let (source, _) = temperature_interface(&cards, 0.1);
    let context = point();
    let law = MaterialContactLaw::new(
        binding(&cards, &source, &context),
        initial,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let reference = law.reference().clone();
    let transport = LumpedThermalTransport::try_declared(100.0, 0.0).unwrap();
    support::with_cx(|cx| {
        let mut exhausted = contact_controls();
        exhausted.max_iterations = 1;
        assert!(
            law.advance(cx, inputs, [&transport; 2], 1e-6, 1e-14, exhausted)
                .is_err()
        );
        let mut claims = ClaimSet::new();
        let pinned = claims
            .insert_claim(temperature_claim(0.1, 304.9, 320.0))
            .unwrap();
        claims
            .insert_claim(temperature_claim(0.1, 300.0, 304.899))
            .unwrap();
        let limited = interface(&cards, claims);
        let pinned_law = MaterialContactLaw::new(
            binding(&cards, &limited, &context),
            initial,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let original = pinned_law.reference().clone();
        assert_eq!(original.properties()[0].answer().receipt.selected, pinned);
        assert!(
            pinned_law
                .advance(cx, inputs, [&transport; 2], 1e-6, 1e-14, contact_controls())
                .is_err()
        );
        assert_eq!(pinned_law.reference(), &original);
    });
    support::with_cancelled_cx(|cx| {
        assert!(matches!(
            law.advance(cx, inputs, [&transport; 2], 1e-6, 1e-14, contact_controls()),
            Err(fs_conduction::ConductionError::Cancelled { .. })
        ));
    });
    assert_eq!(law.reference(), &reference);
    assert_eq!([inputs[0].initial, inputs[1].initial], initial);
}

#[test]
fn g0_endpoint_contact_biot_cannot_hide_behind_a_valid_starting_resistance() {
    let cards = [bulk("A"), bulk("B")];
    let curves = caloric(&cards);
    let inputs = [
        thermal_input(&curves[0], 0.2),
        thermal_input(&curves[1], 0.0),
    ];
    let initial = [inputs[0].initial, inputs[1].initial];
    let (source, id) = temperature_interface(&cards, 0.1);
    let context = point();
    let b = binding(&cards, &source, &context);
    let law =
        MaterialContactLaw::new(b, initial, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let frozen = InterfaceResistance::from_card_pinned(
        "material contact",
        &source,
        &b.at(initial).unwrap(),
        id,
    )
    .unwrap();
    let transport = LumpedThermalTransport::try_declared(0.19, 0.0).unwrap();
    support::with_cx(|cx| {
        let held = advance_contact_pair(cx, inputs, [&transport; 2], &frozen, 1e-6, 1e-14).unwrap();
        let applied_biot = held[0].report.contact.unwrap().maximum_biot;
        assert!(applied_biot < fs_conduction::lumped::LUMPED_BIOT_CEILING);
        let mean = 0.5 * (held[0].state.temperature_k() + held[1].state.temperature_k());
        let endpoint_r = 1e-4 * (1.0 + 0.1 * (mean - 300.0));
        assert!(
            applied_biot * frozen.value_m2_k_per_w() / endpoint_r
                > fs_conduction::lumped::LUMPED_BIOT_CEILING
        );
        assert!(
            law.advance(cx, inputs, [&transport; 2], 1e-6, 1e-14, contact_controls())
                .is_err()
        );
    });
}
