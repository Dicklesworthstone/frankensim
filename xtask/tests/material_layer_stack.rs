//! Ordered multilayer heat flux, original source evidence and conservative storage.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint,
    SelectionPolicy,
};
use fs_material::thermal::{
    conduction::{
        ConductivityCurve, LayerConductionConfig, LayerStackConduction, MAX_CONDUCTION_LAYERS,
        TwoLayerConduction,
    },
    transport::{
        stack::{LayerStackSourceQuery, PreparedLayerStackConduction},
        storage::{
            HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw,
            PreparedHeatStorage,
        },
    },
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
    backward_euler::{BackwardEulerStep, StepCoefficients},
};
use std::{collections::BTreeMap, path::Path, process::Command};

const CONFIG: LayerConductionConfig = LayerConductionConfig {
    max_iterations: 128,
    absolute_heat_tolerance_w: 1e-10,
    relative_heat_tolerance: 1e-11,
};
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const CVOL: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const CHARTS: [[[f64; 2]; 4]; 4] = [
    [[250.0, 0.2], [310.0, 0.3], [390.0, 0.4], [450.0, 0.25]],
    [[250.0, 0.6], [340.0, 0.4], [370.0, 0.8], [450.0, 0.5]],
    [[250.0, 0.5], [300.0, 0.3], [340.0, 0.6], [450.0, 0.7]],
    [[250.0, 0.15], [315.0, 0.25], [370.0, 0.4], [450.0, 0.5]],
];

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance {tolerance:.3e}"
    );
}
fn constant(k: f64) -> ConductivityCurve {
    ConductivityCurve::new(vec![[250.0, k], [450.0, k]]).unwrap()
}
// Independent affine primitives manufacture geometry from prescribed face
// temperatures. Neither the production integral nor inverse supplies the oracle.
fn integral(knots: &[[f64; 2]], lo: f64, hi: f64) -> f64 {
    let mut value = 0.0;
    for pair in knots.windows(2) {
        let [[a, ka], [b, kb]] = [pair[0], pair[1]];
        let start = lo.max(a);
        let end = hi.min(b);
        if end > start {
            let slope = (kb - ka) / (b - a);
            value += (end - start) * (ka + slope * (0.5 * (start + end) - a));
        }
    }
    value
}

#[test]
fn g1_single_and_constant_stacks_retain_exact_series_resistance_and_curve_inverse() {
    let varying =
        ConductivityCurve::new(vec![[250.0, 0.2], [290.0, 0.4], [330.0, 0.1], [400.0, 0.7]])
            .unwrap();
    let amount = 20.0 + 27.0 / 7.0;
    close(
        varying
            .temperature_after_integral_w_m(270.0, amount)
            .unwrap(),
        360.0,
        2e-12,
    );
    close(
        varying
            .temperature_after_integral_w_m(360.0, -amount)
            .unwrap(),
        270.0,
        2e-12,
    );
    assert_eq!(
        varying.temperature_after_integral_w_m(250.0, 50.0).unwrap(),
        400.0
    );
    assert_eq!(
        varying
            .temperature_after_integral_w_m(400.0, -50.0)
            .unwrap(),
        250.0
    );
    assert_eq!(
        varying.temperature_after_integral_w_m(330.0, 0.0).unwrap(),
        330.0
    );
    assert!(
        varying
            .temperature_after_integral_w_m(250.0, 50.0_f64.next_up())
            .is_err()
    );
    assert!(varying.temperature_after_integral_w_m(249.0, 0.0).is_err());
    assert!(
        varying
            .temperature_after_integral_w_m(300.0, f64::NAN)
            .is_err()
    );

    let one = [varying];
    let single = LayerStackConduction::new(&one, &[0.02], &[], 0.1).unwrap();
    let response = single.evaluate([400.0, 250.0], CONFIG, || Ok(())).unwrap();
    close(response.power_a_to_b_w, 250.0, 1e-12);
    close(response.conductance_w_k, 250.0 / 150.0, 1e-14);
    close(response.endpoint_derivative_w_k[0], 3.5, 1e-14);
    close(response.endpoint_derivative_w_k[1], -1.0, 1e-14);
    assert_eq!(response.layer_temperatures_k(), &[[400.0, 250.0]]);

    let ks = [0.2, 0.7, 2.0];
    let curves = ks.map(constant);
    let lengths = [0.01, 0.025, 0.005];
    let contacts = [0.004, 0.002];
    let area = 0.03;
    let resistance =
        lengths.iter().zip(ks).map(|(d, k)| d / k).sum::<f64>() + contacts.iter().sum::<f64>();
    let flux = 50.0 / resistance;
    let stack = LayerStackConduction::new(&curves, &lengths, &contacts, area).unwrap();
    let response = stack.evaluate([350.0, 300.0], CONFIG, || Ok(())).unwrap();
    close(response.power_a_to_b_w, area * flux, 3e-9);
    close(response.conductance_w_k, area / resistance, 6e-11);
    close(
        response.endpoint_derivative_w_k[0],
        area / resistance,
        1e-12,
    );
    close(
        response.endpoint_derivative_w_k[1],
        -area / resistance,
        1e-12,
    );
    let mut inlet = 350.0;
    for (i, pair) in response.layer_temperatures_k().iter().enumerate() {
        let outlet = inlet - flux * lengths[i] / ks[i];
        close(pair[0], inlet, 1e-8);
        close(pair[1], outlet, 1e-8);
        if i < contacts.len() {
            inlet = outlet - flux * contacts[i];
        }
    }
    let equal = stack.evaluate([315.0; 2], CONFIG, || Ok(())).unwrap();
    assert_eq!(equal.power_a_to_b_w, 0.0);
    close(equal.conductance_w_k, area / resistance, 1e-12);
    assert_eq!(equal.layer_temperatures_k(), &[[315.0; 2]; 3]);
    assert_eq!(equal.iterations, 0);

    let curves = vec![constant(0.5); MAX_CONDUCTION_LAYERS];
    let lengths = vec![0.001; MAX_CONDUCTION_LAYERS];
    let contacts = vec![0.0; MAX_CONDUCTION_LAYERS - 1];
    let maximum = LayerStackConduction::new(&curves, &lengths, &contacts, 1.0).unwrap();
    let response = maximum.evaluate([310.0, 300.0], CONFIG, || Ok(())).unwrap();
    close(
        response.power_a_to_b_w,
        10.0 / (0.002 * MAX_CONDUCTION_LAYERS as f64),
        2e-8,
    );
    assert_eq!(response.layer_temperatures_k().len(), MAX_CONDUCTION_LAYERS);
}

#[test]
fn g1_g3_nonuniform_stacks_resolve_source_breaks_contacts_and_physical_tangents() {
    let faces = [
        [400.0, 380.0],
        [377.0, 350.0],
        [350.0, 330.0],
        [325.0, 300.0],
    ];
    let flux = 7.5;
    let area = 0.04;
    let curves = CHARTS.map(|knots| ConductivityCurve::new(knots.to_vec()).unwrap());
    let lengths: [f64; 4] =
        std::array::from_fn(|i| integral(&CHARTS[i], faces[i][1], faces[i][0]) / flux);
    let contacts = [3.0 / flux, 0.0, 5.0 / flux];
    let stack = LayerStackConduction::new(&curves, &lengths, &contacts, area).unwrap();
    let response = stack.evaluate([400.0, 300.0], CONFIG, || Ok(())).unwrap();
    close(response.power_a_to_b_w, 0.3, 2e-10);
    close(response.conductance_w_k, 0.003, 2e-12);
    for (actual, expected) in response.layer_temperatures_k().iter().zip(faces) {
        close(actual[0], expected[0], 2e-7);
        close(actual[1], expected[1], 2e-7);
    }
    assert!(response.max_layer_heat_residual_w <= 2e-10);
    assert!(
        response.max_contact_temperature_residual_k
            <= response.max_contact_roundoff_allowance_k + 1e-7
    );
    for endpoint in 0..2 {
        let mut plus = [400.0, 300.0];
        let mut minus = plus;
        plus[endpoint] += 0.01;
        minus[endpoint] -= 0.01;
        let q_plus = stack
            .evaluate(plus, CONFIG, || Ok(()))
            .unwrap()
            .power_a_to_b_w;
        let q_minus = stack
            .evaluate(minus, CONFIG, || Ok(()))
            .unwrap()
            .power_a_to_b_w;
        close(
            response.endpoint_derivative_w_k[endpoint],
            (q_plus - q_minus) / 0.02,
            5e-8,
        );
    }
    assert!(response.endpoint_derivative_w_k[0] > 0.0);
    assert!(response.endpoint_derivative_w_k[1] < 0.0);

    // Reverse the physical geometry as well as the endpoints: changing only
    // endpoints on heterogeneous temperature-dependent layers need not give -Q.
    let mut reverse_curves = curves.clone();
    reverse_curves.reverse();
    let mut reverse_lengths = lengths;
    reverse_lengths.reverse();
    let mut reverse_contacts = contacts;
    reverse_contacts.reverse();
    let reversed =
        LayerStackConduction::new(&reverse_curves, &reverse_lengths, &reverse_contacts, area)
            .unwrap()
            .evaluate([300.0, 400.0], CONFIG, || Ok(()))
            .unwrap();
    close(reversed.power_a_to_b_w, -response.power_a_to_b_w, 2e-10);
    for (actual, expected) in reversed
        .layer_temperatures_k()
        .iter()
        .zip(response.layer_temperatures_k().iter().rev())
    {
        close(actual[0], expected[1], 2e-7);
        close(actual[1], expected[0], 2e-7);
    }
    close(
        reversed.endpoint_derivative_w_k[0],
        -response.endpoint_derivative_w_k[1],
        1e-10,
    );
    close(
        reversed.endpoint_derivative_w_k[1],
        -response.endpoint_derivative_w_k[0],
        1e-10,
    );

    let pair =
        LayerStackConduction::new(&curves[..2], &lengths[..2], &contacts[..1], area).unwrap();
    let original = TwoLayerConduction::new(
        [&curves[0], &curves[1]],
        [lengths[0], lengths[1]],
        area,
        contacts[0],
    )
    .unwrap();
    for endpoints in [[400.0, 300.0], [300.0, 400.0]] {
        let stack_response = pair.evaluate(endpoints, CONFIG, || Ok(())).unwrap();
        let pair_response = original.evaluate(endpoints, CONFIG, || Ok(())).unwrap();
        close(
            stack_response.power_a_to_b_w,
            pair_response.power_a_to_b_w,
            3e-10,
        );
        close(
            stack_response.layer_temperatures_k()[0][1],
            pair_response.interface_temperature_k[0],
            2e-7,
        );
        close(
            stack_response.layer_temperatures_k()[1][0],
            pair_response.interface_temperature_k[1],
            2e-7,
        );
        for i in 0..2 {
            close(
                stack_response.endpoint_derivative_w_k[i],
                pair_response.endpoint_derivative_w_k[i],
                1e-9,
            );
        }
    }
    let doubled_lengths = lengths.map(|d| 2.0 * d);
    let doubled_contacts = contacts.map(|r| 2.0 * r);
    let scaled = LayerStackConduction::new(&curves, &doubled_lengths, &doubled_contacts, area)
        .unwrap()
        .evaluate([400.0, 300.0], CONFIG, || Ok(()))
        .unwrap();
    close(scaled.power_a_to_b_w, 0.5 * response.power_a_to_b_w, 2e-10);
}

#[test]
fn g0_g4_stack_geometry_support_and_work_refusals_do_not_change_inputs() {
    let curves = CHARTS.map(|knots| ConductivityCurve::new(knots.to_vec()).unwrap());
    let original = curves.clone();
    let lengths = [0.01, 0.02, 0.03, 0.04];
    let contacts = [0.01, 0.0, 0.03];
    assert!(LayerStackConduction::new(&[], &[], &[], 0.01).is_err());
    assert!(LayerStackConduction::new(&curves, &lengths[..3], &contacts, 0.01).is_err());
    assert!(LayerStackConduction::new(&curves, &lengths, &contacts[..2], 0.01).is_err());
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(LayerStackConduction::new(&curves, &lengths, &contacts, invalid).is_err());
        let mut bad = lengths;
        bad[2] = invalid;
        assert!(LayerStackConduction::new(&curves, &bad, &contacts, 0.01).is_err());
    }
    for invalid in [-1.0, f64::INFINITY, f64::NAN] {
        let mut bad = contacts;
        bad[1] = invalid;
        assert!(LayerStackConduction::new(&curves, &lengths, &bad, 0.01).is_err());
    }
    let excessive = vec![constant(0.5); MAX_CONDUCTION_LAYERS + 1];
    assert!(
        LayerStackConduction::new(
            &excessive,
            &vec![0.01; excessive.len()],
            &vec![0.0; excessive.len() - 1],
            1.0
        )
        .is_err()
    );
    let stack = LayerStackConduction::new(&curves, &lengths, &contacts, 0.01).unwrap();
    for bad in [[249.0, 300.0], [300.0, 451.0], [f64::NAN, 300.0]] {
        assert!(stack.evaluate(bad, CONFIG, || Ok(())).is_err());
    }
    assert!(
        stack
            .evaluate(
                [400.0, 300.0],
                LayerConductionConfig {
                    max_iterations: 0,
                    ..CONFIG
                },
                || Ok(())
            )
            .is_err()
    );
    assert!(matches!(
        stack.evaluate(
            [400.0, 300.0],
            LayerConductionConfig {
                max_iterations: 1,
                ..CONFIG
            },
            || Ok(())
        ),
        Err(DepositionError::Budget(_))
    ));
    assert!(
        stack
            .evaluate(
                [400.0, 300.0],
                LayerConductionConfig {
                    absolute_heat_tolerance_w: -1.0,
                    ..CONFIG
                },
                || Ok(())
            )
            .is_err()
    );
    let mut calls = 0;
    let cancelled = stack.evaluate([400.0, 300.0], CONFIG, || {
        calls += 1;
        if calls == 3 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(cancelled, Err(DepositionError::Cancelled)));
    assert_eq!(calls, 3);
    assert_eq!(curves, original);
    assert!(stack.evaluate([400.0, 300.0], CONFIG, || Ok(())).is_ok());

    let partial = [
        ConductivityCurve::new(vec![[300.0, 0.2], [400.0, 0.3]]).unwrap(),
        ConductivityCurve::new(vec![[250.0, 0.4], [350.0, 0.5]]).unwrap(),
    ];
    let stack = LayerStackConduction::new(&partial, &[0.01; 2], &[0.0], 0.01).unwrap();
    // Supporting each own endpoint is insufficient for a conservative full-span trial.
    assert!(stack.evaluate([390.0, 290.0], CONFIG, || Ok(())).is_err());
}

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn conductivity_key() -> PropertyKey {
    PropertyKey::with_quantity(
        "thermal-conductivity",
        quantity(QuantityKind::ThermalConductivity),
    )
}
fn source_point(card: &MaterialCard, key: &PropertyKey, temperature: f64) -> QueryPoint {
    let (_, claim) = card
        .claims()
        .claims_ordered()
        .find(|(_, c)| &c.key == key)
        .unwrap();
    let mut point = QueryPoint::new();
    for (axis, &(low, high)) in claim.validity.bounds() {
        let value = if axis == "temperature" {
            temperature
        } else {
            assert_eq!(low, high);
            low
        };
        point = point
            .with_quantity(axis, claim.validity.axis_quantities()[axis], value)
            .unwrap();
    }
    point
}
fn compiled(slug: &str) -> MaterialCard {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "layer-stack-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process: "original FSRI population".into(),
            revision: 0,
        },
        pack.claims().clone(),
        vec![],
    )
    .unwrap()
}
fn source_layer<'a>(
    id: u64,
    card: &'a MaterialCard,
    key: &'a PropertyKey,
    thickness: f64,
    lower: &'a QueryPoint,
    upper: &'a QueryPoint,
) -> LayerStackSourceQuery<'a> {
    LayerStackSourceQuery {
        id,
        card,
        conductivity: key,
        thickness_m: thickness,
        temperature_axis: "temperature",
        lower,
        upper,
        selection: SINGLE,
    }
}
fn storage(card: &MaterialCard) -> PreparedHeatStorage {
    let key = PropertyKey::new("volumetric-heat-capacity", CVOL);
    let lower = source_point(card, &key, 288.15);
    let upper = source_point(card, &key, 313.15);
    PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
        id: 1,
        card,
        capacity: &key,
        amount: HeatStorageAmount::ReferenceVolumeM3(0.001),
        temperature_axis: "temperature",
        lower: &lower,
        upper: &upper,
        selection: SINGLE,
    }])
    .unwrap()
}
fn cell(law: &HeatStorageThermalLaw<'_>, id: u64, temperature: f64) -> DepositedCell {
    let cell = DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(id),
        mass_kg: 1.0,
        specific_h_j_kg: law.specific_storage_j_kg(temperature).unwrap(),
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.1, 0.0, 0.0],
            volume_m3: 0.001,
            surface_area_m2: 0.1,
        },
        boundary: Boundary {
            ambient_k: 293.15,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    law.admit_cell(&cell).unwrap();
    cell
}

#[test]
fn g1_g3_sourced_wall_stack_drives_original_implicit_energy_and_entropy_owner() {
    let cards = [
        compiled("fsri-gypsum-wallboard-unconditioned-2021"),
        compiled("fsri-mineral-wool-unconditioned-2021"),
        compiled("fsri-plywood-unconditioned-2021"),
    ];
    let key = conductivity_key();
    let lower = cards
        .each_ref()
        .map(|card| source_point(card, &key, 288.15));
    let upper = cards
        .each_ref()
        .map(|card| source_point(card, &key, 313.15));
    // Literal FSRI conductivity means at 15/45 C. The geometry is an authored
    // homogeneous through-thickness model, not qualification of other panels.
    let source_k = [[0.152, 0.153], [0.033, 0.038], [0.113, 0.120]];
    let faces = [[302.15, 300.15], [299.65, 296.65], [296.15, 294.15]];
    let flux = 10.0;
    let area = 0.01;
    let lengths: [f64; 3] = std::array::from_fn(|i| {
        let [hot, cold] = faces[i];
        let [k15, k45] = source_k[i];
        (hot - cold) * (k15 + (k45 - k15) / 30.0 * (0.5 * (hot + cold) - 288.15)) / flux
    });
    let contacts = [0.05, 0.05];
    // Stable IDs deliberately differ from the physical order.
    let ids = [8, 3, 12];
    let queries: [LayerStackSourceQuery<'_>; 3] = std::array::from_fn(|i| {
        source_layer(ids[i], &cards[i], &key, lengths[i], &lower[i], &upper[i])
    });
    let stack = PreparedLayerStackConduction::try_new(&queries, &contacts, area).unwrap();
    assert_eq!(stack.temperature_span_k(), [288.15, 313.15]);
    assert_eq!(stack.contact_resistances_m2_k_w(), &contacts);
    assert_eq!(stack.area_m2(), area);
    for (i, layer) in stack.layers().iter().enumerate() {
        assert_eq!(layer.id(), ids[i]);
        assert_eq!(layer.material_card_identity(), cards[i].content_hash());
        assert_eq!(layer.thickness_m(), lengths[i]);
        assert_eq!(layer.temperature_axis(), "temperature");
        assert_eq!(
            layer.query_endpoints(),
            &[lower[i].clone(), upper[i].clone()]
        );
        for (_, answer) in &layer.source_line().knots {
            cards[i].claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    let response = stack.evaluate([302.15, 294.15], CONFIG, || Ok(())).unwrap();
    close(response.power_a_to_b_w, 0.1, 2e-10);
    for (actual, expected) in response.layer_temperatures_k().iter().zip(faces) {
        close(actual[0], expected[0], 2e-7);
        close(actual[1], expected[1], 2e-7);
    }
    let mut reordered = queries;
    reordered.reverse();
    assert_ne!(
        stack.identity(),
        PreparedLayerStackConduction::try_new(&reordered, &contacts, area)
            .unwrap()
            .identity()
    );
    assert_ne!(
        stack.identity(),
        PreparedLayerStackConduction::try_new(&queries, &contacts, 2.0 * area)
            .unwrap()
            .identity()
    );
    let mut duplicate = queries;
    duplicate[2].id = ids[0];
    assert!(PreparedLayerStackConduction::try_new(&duplicate, &contacts, area).is_err());
    let bad_context = lower[0]
        .clone()
        .with_quantity(
            "mean-test-temperature-as-local-model",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    let mut wrong_context = queries;
    wrong_context[0].lower = &bad_context;
    assert!(PreparedLayerStackConduction::try_new(&wrong_context, &contacts, area).is_err());
    let capacity_key = PropertyKey::new("volumetric-heat-capacity", CVOL);
    let mut wrong_quantity = queries;
    wrong_quantity[0].conductivity = &capacity_key;
    assert!(PreparedLayerStackConduction::try_new(&wrong_quantity, &contacts, area).is_err());
    assert!(
        stack
            .evaluate([288.15_f64.next_down(), 300.0], CONFIG, || Ok(()))
            .is_err()
    );

    // Two separately declared endpoint stores exchange through a massless,
    // quasi-steady path. This does not add storage inside the conducting layers.
    let hot_storage = storage(&cards[0]);
    let cold_storage = storage(&cards[2]);
    // Explicit reference masses convert native volume storage into the owner's
    // specific coordinate; neither density nor matched source specimens is inferred.
    let hot_law = HeatStorageThermalLaw::try_new(&hot_storage, 1.0, 293.15)
        .unwrap()
        .with_sensible_entropy();
    let cold_law = HeatStorageThermalLaw::try_new(&cold_storage, 1.0, 293.15)
        .unwrap()
        .with_sensible_entropy();
    let capacity = [
        [627851.6666666666, 643764.6666666666],
        [689262.6666666666, 724033.3333333334],
    ];
    let energy = |t: f64, c: [f64; 2]| {
        let x = t - 293.15;
        0.001 * (c[0] * x + 0.5 * (c[1] - c[0]) / 10.0 * x * x)
    };
    let inverse = |e: f64, c: [f64; 2]| {
        let base = 0.001 * c[0];
        let slope = 0.001 * (c[1] - c[0]) / 10.0;
        293.15 + 2.0 * e / (base + (base * base + 2.0 * slope * e).sqrt())
    };
    let duration = 500.0;
    let transfer = duration * area * flux;
    let final_t = [302.15, 294.15];
    let initial_t = [
        inverse(energy(final_t[0], capacity[0]) + transfer, capacity[0]),
        inverse(energy(final_t[1], capacity[1]) - transfer, capacity[1]),
    ];
    let link = ThermalLink {
        id: 7,
        a: CellId(1),
        b: CellId(2),
        area_m2: area,
        conductance_w_k: 0.0,
    };
    let initial = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell(&hot_law, 1, initial_t[0])),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell(&cold_law, 2, initial_t[1])),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(link),
            },
        ],
        &[
            (MaterialId(1), &hot_law as &dyn ThermalLaw),
            (MaterialId(2), &cold_law as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let initial_energy = initial.energy_j();
    let step = BackwardEulerStep::new(&initial, duration, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_entropy(
            &[-transfer, transfer],
            |candidate| {
                let response = stack.evaluate(
                    [
                        candidate.temperature_k(CellId(1))?,
                        candidate.temperature_k(CellId(2))?,
                    ],
                    CONFIG,
                    || Ok(()),
                )?;
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::from([(7, response.conductance_w_k)]),
                    payload: response,
                })
            },
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-7);
    close(trial.payload().power_a_to_b_w, 0.1, 2e-10);
    let entropy = trial.entropy_audit().unwrap();
    let delta_entropy = |a: f64, b: f64, c: [f64; 2]| {
        let slope = (c[1] - c[0]) / 10.0;
        0.001 * ((c[0] - slope * 293.15) * ((b - a) / a).ln_1p() + slope * (b - a))
    };
    let expected_entropy = delta_entropy(initial_t[0], final_t[0], capacity[0])
        + delta_entropy(initial_t[1], final_t[1], capacity[1]);
    close(entropy.thermal_storage_change_j_k, expected_entropy, 2e-10);
    close(
        entropy.conduction_generation_j_k,
        transfer * (1.0 / final_t[1] - 1.0 / final_t[0]),
        2e-10,
    );
    assert!(entropy.thermal_storage_change_j_k > 0.0);
    assert!(entropy.conduction_generation_j_k > 0.0);
    let (accepted, response, audit) = trial.accept(1e-7, 1e-7, || Ok(())).unwrap();
    close(
        accepted.temperature_k(CellId(1)).unwrap(),
        final_t[0],
        2e-12,
    );
    close(
        accepted.temperature_k(CellId(2)).unwrap(),
        final_t[1],
        2e-12,
    );
    close(accepted.energy_j(), initial_energy, 1e-8);
    assert_eq!(audit.external_heat_j, 0.0);
    close(response.power_a_to_b_w, 0.1, 2e-10);
    assert_eq!(accepted.links()[&7], link);
    assert!(
        step.evaluate(
            &[-transfer, transfer],
            |_| Err::<StepCoefficients<()>, _>(DepositionError::Cancelled),
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.energy_j(), initial_energy);
    close(
        initial.temperature_k(CellId(1)).unwrap(),
        initial_t[0],
        2e-12,
    );
}
