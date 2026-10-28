//! G1/G3 capacity baths through native project parsing, imported geometry,
//! accepted time steps and the same complete physical trajectory transpose.

use super::*;

const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";
const DBATH: &str = "dtemperature_dprescribed_temperature_k_per_k";

fn set_bath(spec: &mut ProjectSpec, temperature_k: f64) {
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.iter_mut().find(|row| row.target == "driven-face").unwrap().condition =
        ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(temperature_k, dims::TEMPERATURE),
        };
}

fn set_initial(spec: &mut ProjectSpec, temperature_k: f64) {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap().initial_temperature =
        fs_project::TransientInitialTemperature::Uniform(
            QtyAny::new(temperature_k, dims::TEMPERATURE));
}

fn gradient(solved: &ThermalRun) -> &J {
    solved.transient().get("adjoint").unwrap()
}

fn bath(solved: &ThermalRun) -> &J {
    let result = gradient(solved);
    let family = result.get("prescribed_temperature").unwrap();
    close(n(family, "control_steps"), n(result, "state_index"), 0.0);
    assert_eq!(family.str_field("parameter_scope"),
        Some("independent-uniform-bath-temperatures-held-throughout-selected-trajectory"));
    assert!(family.str_field("scope").unwrap().contains("Initial temperatures"));
    let rows = family.get("boundaries").unwrap().as_array().unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].str_field("target"), Some("driven-face"));
    assert_eq!(rows[0].str_field("status"), Some("available"));
    close(n(&rows[0], "mapped_vertices"), 3.0, 0.0);
    &rows[0]
}

fn requested(spec: &ProjectSpec, cards: &CardPackSet, bytes: &[u8], output: &str) -> ThermalRun {
    let mut changed = spec.clone();
    request_air_adjoint(&mut changed, output);
    run(&changed, cards, bytes)
}

fn simple_fixture() -> (ProjectSpec, CardPackSet, Vec<u8>) {
    let bytes = tetra_stl();
    let mut spec = transient_product::driven_face_project(&bytes);
    let cards = transient_product::transient_conductivity_cards("constant-scalar");
    rebind_to(&mut spec, &cards);
    set_bath(&mut spec, 295.0);
    set_initial(&mut spec, 330.0);
    spec.power.as_mut().unwrap()[0].watts.value = 4.0;
    spec.power.as_mut().unwrap()[0].duty = 1.0;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    // Unequal fine steps prevent an accidental extra or omitted dt from hiding.
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap().horizon.value = 0.7;
    (spec, cards, bytes)
}

#[test]
fn g1_native_capacity_bath_adjoint_matches_independent_scalar_time_balance() {
    let (spec, cards, bytes) = simple_fixture();
    let plain = run(&spec, &cards, &bytes);
    let solved = requested(&spec, &cards, &bytes, FINAL);
    assert_eq!(plain.temperature, solved.temperature);
    assert_eq!(plain.rows("fine"), solved.rows("fine"));
    assert_eq!(plain.rows("coarse"), solved.rows("coarse"));
    let result = gradient(&solved);
    close(n(result, "active_ties"), 1.0, 0.0);
    close(n(bath(&solved), "temperature_k"), 295.0, 0.0);
    // Unit tetrahedron: one free origin, capacity rhoCp*V/4 = 15/24,
    // conductance k/2 = 1 and one quarter of the 4-W source at the free node.
    let (capacity, conductance) = (15.0 / 24.0, 1.0);
    let (mut temperature, mut derivative) = (330.0, 0.0);
    let mut widths = Vec::new();
    for interval in result.get("intervals").unwrap().as_array().unwrap() {
        let dt = n(interval, "end_s") - n(interval, "start_s");
        widths.push(dt);
        let denominator = capacity + dt * conductance;
        temperature = (capacity * temperature + dt * (conductance * 295.0 + 1.0)) / denominator;
        derivative = (capacity * derivative + dt * conductance) / denominator;
    }
    assert!(widths.iter().any(|dt| (dt - widths[0]).abs() > 0.01));
    close(final_k(&solved), temperature, 2e-8);
    close(n(bath(&solved), DBATH), derivative, 2e-9);
    let evaluate = |delta: f64| {
        let mut probe = spec.clone();
        set_bath(&mut probe, 295.0 + delta);
        final_k(&run(&probe, &cards, &bytes))
    };
    close(derivative, (evaluate(0.001) - evaluate(-0.001)) / 0.002, 2e-7);
    let repeated = requested(&spec, &cards, &bytes, FINAL);
    assert_eq!(bath(&solved), bath(&repeated));
}

#[test]
fn g1_native_capacity_bath_direct_goal_and_initial_peak_have_exact_partials() {
    let (mut spec, cards, bytes) = simple_fixture();
    spec.power.as_mut().unwrap()[0].watts.value = 0.0;
    set_initial(&mut spec, 300.0);
    set_bath(&mut spec, 330.0);
    let fixed = requested(&spec, &cards, &bytes, FINAL);
    close(final_k(&fixed), 330.0, 0.0);
    close(n(bath(&fixed), DBATH), 1.0, 1e-12);
    close(n(gradient(&fixed), "dtemperature_duniform_initial_k_per_k"), 0.0, 1e-12);
    set_initial(&mut spec, 330.0);
    set_bath(&mut spec, 295.0);
    let initial = requested(&spec, &cards, &bytes, PEAK);
    assert!(final_k(&initial) < 330.0);
    close(n(gradient(&initial), "state_index"), 0.0, 0.0);
    close(n(gradient(&initial), "value_k"), 330.0, 0.0);
    close(n(bath(&initial), DBATH), 0.0, 0.0);
    close(n(gradient(&initial), "adjoint_krylov_iterations"), 0.0, 0.0);
}

#[test]
fn g3_native_capacity_bath_keeps_nonlinear_air_natural_and_radiation_feedback() {
    let bytes = tetra_stl();
    let mut spec = driven_air_project(&bytes);
    let materials = transient_product::transient_conductivity_cards("nonlinear");
    let (cards, radiation_card) = radiation_product::emissivity_cards(&materials, EMISSIVITY);
    rebind_to(&mut spec, &cards);
    set_bath(&mut spec, 295.0);
    set_initial(&mut spec, 330.0);
    spec.power.as_mut().unwrap()[0].watts.value = 4.0;
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient.as_mut().unwrap().horizon.value = 0.7;
    setup.boundaries.iter_mut().find(|row| row.target == "air-face-2").unwrap().condition =
        ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: QtyAny::new(0.1, dims::LENGTH),
            ambient_temperature: QtyAny::new(350.0, dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    radiation_product::declare(&mut spec, &radiation_card, "air-face-1", 280.0);
    let radiation = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap();
    radiation.max_iterations = 32;
    radiation.relaxation = 1.0;
    let plain = run(&spec, &cards, &bytes);
    let solved = requested(&spec, &cards, &bytes, FINAL);
    assert_eq!(plain.temperature, solved.temperature);
    assert_eq!(plain.rows("fine"), solved.rows("fine"));
    let result = gradient(&solved);
    close(n(result, "active_ties"), 1.0, 0.0);
    close(n(result, "reconstructed_solid_endpoints"), 0.0, 0.0);
    close(n(result, "revalidated_solid_endpoints"), solved.rows("fine").len() as f64, 0.0);
    assert!(result.get("convection").is_some());
    assert!(result.get("radiation").is_some());
    let derivative = n(bath(&solved), DBATH);
    assert!(derivative.abs() > 1e-4);
    let evaluate = |delta: f64| {
        let mut probe = spec.clone();
        set_bath(&mut probe, 295.0 + delta);
        final_k(&run(&probe, &cards, &bytes))
    };
    for epsilon in [0.002, 0.001] {
        let finite_difference = (evaluate(epsilon) - evaluate(-epsilon)) / (2.0 * epsilon);
        close(derivative, finite_difference, 5e-5 * derivative.abs().max(1.0));
    }
}
