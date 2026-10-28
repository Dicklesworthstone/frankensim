//! Complete capacity-storage enclosure reverse history against fresh native solves.
//! This uses the genuine C(T-new - T-old) owner, including smooth k(T), instead
//! of translating the problem into a total-enthalpy chart.
use super::*;
use fs_project::TransientInitialTemperature;

const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";
const SECOND_EPSILON: f64 = 0.6;

fn fixture(nonlinear: bool) -> (ProjectSpec, CardPackSet) {
    let material = if nonlinear {
        transient_product::transient_conductivity_cards("nonlinear")
    } else {
        fixture_cards()
    };
    let (cards, first) = radiation_product::emissivity_cards(&material, EPSILON);
    let (cards, second) = radiation_product::emissivity_cards(&cards, SECOND_EPSILON);
    let mut spec = project(&first);
    if nonlinear {
        rebind_to(&mut spec, &material);
        let binding = &mut spec.materials.as_mut().unwrap()[0];
        binding.temp_lo.value = 280.0;
        binding.temp_hi.value = 360.0;
    }
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap().surfaces[1].card = second;
    let policy = time(&mut spec);
    policy.enthalpy = None;
    policy.capacities = vec![TransientRegionCapacity {
        region: "air".into(),
        volumetric_heat_capacity: QtyAny::new(DENSITY * CP, dims::VOLUMETRIC_HEAT_CAPACITY),
        source: "synthetic positive volumetric capacity for enclosure history derivatives".into(),
    }];
    policy.initial_temperature = TransientInitialTemperature::Uniform(
        QtyAny::new(INITIAL_K, dims::TEMPERATURE));
    (spec, cards)
}

fn evaluate(spec: &ProjectSpec, cards: &CardPackSet, output: Option<&str>) -> (J, Vec<f64>) {
    let mut spec = spec.clone();
    if let Some(name) = output {
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: name.into(), kind: "report".into(), region: None,
        });
    }
    let solved = run(&spec, cards);
    let receipt = J::parse(&solved.conduction).unwrap();
    let field = J::parse(&solved.field).unwrap();
    (receipt, values(&field, "temperature"))
}

fn maximum(field: &[f64]) -> f64 {
    field.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn gradient(receipt: &J) -> &J {
    let row = receipt.get("transient").unwrap().get("adjoint").unwrap();
    assert_eq!(row.str_field("schema"), Some("fs-cli-native-transient-adjoint-v1"));
    row
}

fn surface(row: &J, index: usize) -> &J {
    &row.get("radiation").unwrap().get("surfaces").unwrap().as_array().unwrap()[index]
}

fn derivative_close(actual: f64, expected: f64) {
    close(actual, expected, 5e-6 + 3e-4 * expected.abs());
}

fn set_initial(spec: &mut ProjectSpec, temperature: f64) {
    time(spec).initial_temperature = TransientInitialTemperature::Uniform(
        QtyAny::new(temperature, dims::TEMPERATURE));
}

fn emissivity_probe(spec: &ProjectSpec, cards: &CardPackSet, index: usize, value: f64) -> f64 {
    let (cards, card) = radiation_product::emissivity_cards(cards, value);
    let mut changed = spec.clone();
    changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap().surfaces[index].card = card;
    maximum(&evaluate(&changed, &cards, None).1)
}

#[test]
fn g1_capacity_enclosure_reverse_keeps_storage_k_t_reflections_and_original_robin_controls() {
    // The nonconstant material card makes an omitted endpoint k'(T) term
    // visible in the history derivative, even though the card is held fixed.
    let (spec, cards) = fixture(true);
    let (forward, plain) = evaluate(&spec, &cards, None);
    let (receipt, retained) = evaluate(&spec, &cards, Some(FINAL));
    assert_eq!(plain.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
        retained.iter().map(|value| value.to_bits()).collect::<Vec<_>>());
    assert_eq!(forward.get("transient").unwrap().get("sampled_peak"),
        receipt.get("transient").unwrap().get("sampled_peak"));
    let row = gradient(&receipt);
    close(n(row, "value_k"), maximum(&retained), 0.0);
    close(n(row, "state_index"), 6.0, 0.0);
    close(n(row, "revalidated_solid_endpoints"), 6.0, 0.0);
    close(n(row, "reconstructed_solid_endpoints"), 0.0, 0.0);
    assert!(n(row, "adjoint_krylov_iterations") <= n(row, "max_adjoint_krylov_iterations"));
    assert!(n(row, "maximum_dual_relative_residual") <= 1e-9);
    assert_eq!(row.get("radiation").unwrap().str_field("method"),
        Some("complete-closed-gray-diffuse-capacity-tangent"));
    let epsilon = 0.0005;
    for (index, nominal) in [EPSILON, SECOND_EPSILON].into_iter().enumerate() {
        let expected = (emissivity_probe(&spec, &cards, index, nominal + epsilon)
            - emissivity_probe(&spec, &cards, index, nominal - epsilon)) / (2.0 * epsilon);
        derivative_close(n(surface(row, index), "dtemperature_demissivity_k"), expected);
    }
    let intervals = row.get("intervals").unwrap().as_array().unwrap();
    let first_pulse: f64 = intervals.iter().filter(|entry| n(entry, "start_s") < 0.3)
        .map(|entry| n(&entry.get("regional_power").unwrap().as_array().unwrap()[0],
            "dtemperature_ddelivered_watts_k_per_w")).sum();
    let expected_direction = 1.5 * n(row, "dtemperature_duniform_initial_k_per_k")
        + 0.2 * n(row, "dtemperature_dcapacity_multiplier_k") + 3.0 * first_pulse;
    let probe = |delta: f64| {
        let mut changed = spec.clone();
        set_initial(&mut changed, INITIAL_K + 1.5 * delta);
        let policy = time(&mut changed);
        policy.capacities[0].volumetric_heat_capacity.value *= 1.0 + 0.2 * delta;
        policy.power_schedules[0].steps[0].watts.value += 3.0 * delta;
        maximum(&evaluate(&changed, &cards, None).1)
    };
    let delta = 0.002;
    derivative_close(expected_direction, (probe(delta) - probe(-delta)) / (2.0 * delta));
    for (index, target) in ["wall-x", "wall-z"].into_iter().enumerate() {
        let robin = row.get("robin").unwrap().as_array().unwrap().iter()
            .find(|entry| entry.str_field("region") == Some(target)).unwrap();
        for logarithmic in [false, true] {
            let probe = |delta: f64| {
                let mut changed = spec.clone();
                let ThermalBoundaryCondition::Convection { coefficient, reference_temperature } =
                    &mut changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .boundaries[index].condition
                else { panic!("expected prescribed convection"); };
                if logarithmic { coefficient.value *= delta.exp(); }
                else { reference_temperature.value += delta; }
                maximum(&evaluate(&changed, &cards, None).1)
            };
            let key = if logarithmic { "dtemperature_dlog_htc_k" }
                else { "dtemperature_dreference_k_per_k" };
            derivative_close(n(robin, key), (probe(delta) - probe(-delta)) / (2.0 * delta));
        }
    }
}

#[test]
fn g1_capacity_enclosure_sampled_peak_preserves_initial_state_and_future_zero_controls() {
    let (mut spec, cards) = fixture(false);
    set_initial(&mut spec, 353.15);
    for step in &mut time(&mut spec).power_schedules[0].steps {
        step.watts.value = 0.0;
    }
    let (receipt, _) = evaluate(&spec, &cards, Some(PEAK));
    let row = gradient(&receipt);
    close(n(row, "value_k"), 353.15, 0.0);
    close(n(row, "state_index"), 0.0, 0.0);
    close(n(row, "dtemperature_duniform_initial_k_per_k"), 1.0, 0.0);
    close(n(row, "dtemperature_dcapacity_multiplier_k"), 0.0, 0.0);
    for index in 0..2 {
        close(n(surface(row, index), "dtemperature_demissivity_k"), 0.0, 0.0);
    }
    for interval in row.get("intervals").unwrap().as_array().unwrap() {
        close(n(&interval.get("regional_power").unwrap().as_array().unwrap()[0],
            "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
    }
    set_initial(&mut spec, INITIAL_K);
    time(&mut spec).power_schedules[0].steps[0].watts.value = 200.0;
    for boundary in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        let ThermalBoundaryCondition::Convection { reference_temperature, .. } =
            &mut boundary.condition
        else { panic!("expected prescribed convection"); };
        reference_temperature.value = 290.0;
    }
    let (receipt, _) = evaluate(&spec, &cards, Some(PEAK));
    let row = gradient(&receipt);
    let peak = receipt.get("transient").unwrap().get("sampled_peak").unwrap();
    close(n(row, "value_k"), n(peak, "value_k"), 0.0);
    close(n(row, "state_index"), 2.0, 0.0);
    close(n(row, "time_s"), 0.3, 0.0);
    close(n(row, "revalidated_solid_endpoints"), 2.0, 0.0);
    for interval in row.get("intervals").unwrap().as_array().unwrap() {
        if n(interval, "start_s") >= 0.3 {
            close(n(&interval.get("regional_power").unwrap().as_array().unwrap()[0],
                "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
        }
    }
}

fn mixed_fixture() -> (ProjectSpec, CardPackSet) {
    let (mut spec, cards) = fixture(false);
    set_initial(&mut spec, 330.0);
    let envelope = spec.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 299.0;
    envelope.ambient_hi.value = 351.0;
    let prototype = conjugate_fixture_project(7, &tetra_stl(), "convection.gnielinski");
    let mut upstream = prototype.cooling.as_ref().unwrap().conduction.as_ref().unwrap()
        .boundaries[0].condition.clone();
    let ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } = &mut upstream
    else { panic!("expected airflow fixture"); };
    inlet_temperature.value = 300.0;
    let mut downstream = upstream.clone();
    let ThermalBoundaryCondition::AirflowConvection { order, .. } = &mut downstream
    else { panic!("expected airflow fixture"); };
    *order = 1;
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = upstream;
    setup.boundaries[1].condition = ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(350.0, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    };
    setup.boundaries.push(ThermalBoundary { target: "air".into(), condition: downstream });
    (spec, cards)
}

fn convection_control<'a>(row: &'a J, list: &str, key: &str, name: &str) -> &'a J {
    row.get("convection").unwrap().get(list).unwrap().as_array().unwrap().iter()
        .find(|entry| entry.str_field(key) == Some(name)).unwrap()
}

#[test]
fn g3_capacity_enclosure_reverse_composes_natural_and_streamwise_air_feedback() {
    let (spec, cards) = mixed_fixture();
    let (_, plain) = evaluate(&spec, &cards, None);
    let (receipt, retained) = evaluate(&spec, &cards, Some(FINAL));
    assert_eq!(plain, retained);
    let row = gradient(&receipt);
    close(n(row, "revalidated_solid_endpoints"), 6.0, 0.0);
    let natural = convection_control(row, "natural_reservoirs", "target", "wall-z");
    let inlet = convection_control(row, "air_inlets", "branch", "air");
    let fan = convection_control(row, "fan_banks", "bank", "fixture-bank");
    assert_eq!(fan.str_field("status"), Some("available"));
    let direction = 0.5 * n(natural, "dtemperature_dambient_temperature_k_per_k")
        + 0.3 * n(inlet, "dtemperature_dinlet_temperature_k_per_k")
        + 0.02 * n(fan, "dtemperature_dspeed_ratio_k")
        + 0.015 * n(surface(row, 0), "dtemperature_demissivity_k");
    let probe = |delta: f64| {
        let (changed_cards, card) =
            radiation_product::emissivity_cards(&cards, EPSILON + 0.015 * delta);
        let mut changed = spec.clone();
        let cooling = changed.cooling.as_mut().unwrap();
        cooling.fan_system.as_mut().unwrap().banks[0].speed_ratio = 1.0 + 0.02 * delta;
        let setup = cooling.conduction.as_mut().unwrap();
        setup.enclosure_radiation.as_mut().unwrap().surfaces[0].card = card;
        for boundary in &mut setup.boundaries {
            match &mut boundary.condition {
                ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } =>
                    ambient_temperature.value += 0.5 * delta,
                ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } =>
                    inlet_temperature.value += 0.3 * delta,
                _ => {},
            }
        }
        assert_eq!(changed.envelope, spec.envelope);
        maximum(&evaluate(&changed, &changed_cards, None).1)
    };
    let delta = 0.02;
    derivative_close(direction, (probe(delta) - probe(-delta)) / (2.0 * delta));
}

