//! Source data -> conservative transport -> original reed/duct/microphone.
use super::super::super::{ObservationPolicy, ReedPerformance, RenderVoice};
use super::*;
use fs_material::gas_species::input::ThermochemicalPhase;
#[path = "../tests/support.rs"]
mod support;
use support::packs;

const INPUT_FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/examples/reed-steady-gas.performance"
));
/// The fixture's physics without its sidecar declaration: these library tests
/// add their own observation records, and the library never refuses an
/// undeclared law (ruling 6532 puts that in the sidecar writer).
static INPUT: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| INPUT_FIXTURE.replace("observer_law stokes-kirchhoff\n", ""));

fn load(input: &str, sources: &ThermochemicalSources<'_>, block: usize) -> ReedPerformance {
    ReedPerformance::from_bytes_with_gas_sources(
        input.as_bytes(),
        block,
        ObservationPolicy::Diagnostic,
        Some(sources),
    )
    .unwrap()
}
fn render(mut performance: ReedPerformance, block: usize) -> (Vec<f64>, [u64; 4]) {
    let mut pressure = vec![0.0; performance.info().samples as usize];
    for chunk in pressure.chunks_mut(block) {
        performance.pressure_block(chunk).unwrap();
    }
    let RenderVoice::ReedBore(voice) = &performance.renderer().context().voices[0] else {
        panic!("steady gas must retain the original reed");
    };
    assert!(pressure.iter().all(|p| p.is_finite()));
    (
        pressure,
        [
            voice.reed_y,
            voice.reed_v,
            voice.p_plus_prev,
            voice.line.incoming(),
        ]
        .map(f64::to_bits),
    )
}
fn without_steady(input: &str) -> String {
    let begin = input.find("gas_steady ").unwrap();
    let end = input.find("compile_limits ").unwrap();
    format!("{}{}", &input[..begin], &input[end..])
}

#[test]
fn generated_field_and_original_pressure_path_agree_with_direct_owner_binding() {
    let (a, b) = (
        packs("N2", 0.028, 3.5, 200.0),
        packs("CO2", 0.044, 4.5, 200.0),
    );
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let performance = load(&INPUT, &sources, 37);
    let generated = performance.steady_gas().unwrap();
    let problem = generated.problem();
    let solution = generated.solution();
    let control = generated.controls();
    assert_eq!(problem.cells.len(), 2);
    for cell in &problem.cells {
        assert!(
            (cell.volume_m3 - core::f64::consts::PI * 0.0015_f64.powi(2) * 0.285).abs() < 1e-20
        );
        assert!((cell.wall_area_m2 - core::f64::consts::TAU * 0.0015 * 0.285).abs() < 1e-17);
    }
    assert_eq!(
        problem.inlet.composition.fractions(),
        &[0.0, 1.0],
        "canonical CO2,N2 order"
    );
    assert_eq!(
        problem.injections[0].state.composition.fractions(),
        &[1.0, 0.0]
    );
    assert_eq!(solution.mass_flow_kg_s(), &[1e-6, 1.5e-6]);
    assert!(solution.maximum_power_residual_bound_w() <= control.steady_power_tolerance_w);
    assert!(solution.maximum_energy_residual_j() <= control.energy_tolerance_j);
    assert!(solution.substeps() > 0 && solution.substeps() <= control.maximum_substeps);
    assert!(solution.maximum_mean_mach() <= control.maximum_mean_mach);
    let field = solution.field();
    assert_eq!(performance.gas_field().unwrap().0, field);
    // The upstream state is the enthalpy-balanced supply mixture, inverted
    // from h to T by the caloric solve: 330 K only within its J/kg tolerance.
    // Both synthetic species have cp above 800 J/(kg K) (3.5 R/0.028 and
    // 4.5 R/0.044), so that tolerance bounds the temperature error.
    let inlet_k = field.state_at(0.0).unwrap().state().temperature;
    assert!(
        (inlet_k - 330.0).abs() <= control.caloric.energy_tolerance_j_kg / 800.0,
        "inlet {inlet_k} K"
    );
    let outlet = field.state_at(0.57).unwrap();
    assert!((290.0..330.0).contains(&outlet.state().temperature));
    assert_eq!(*outlet.state(), performance.outlet_receiver().unwrap().gas);
    assert_eq!(field.model().source_packs().len(), 2);
    let direct = ReedPerformance::from_bytes_with_field_resolution(
        without_steady(&INPUT).as_bytes(),
        37,
        ObservationPolicy::Diagnostic,
        Some(&sources),
        field.clone(),
        performance.gas_field_resolution().unwrap(),
    )
    .unwrap();
    assert!(direct.steady_gas().is_none());
    let expected = render(performance, 37);
    assert_eq!(expected, render(direct, 37));
    assert_eq!(expected, render(load(&INPUT, &sources, 113), 113));
}

#[test]
fn species_columns_and_data_only_heat_injection_changes_reach_actual_reed_motion() {
    let (a, b) = (
        packs("N2", 0.028, 3.5, 200.0),
        packs("CO2", 0.044, 4.5, 200.0),
    );
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let permuted = INPUT
        .replace("mole 2 N2 CO2", "mole 2 CO2 N2")
        .replace("330 1 0", "330 0 1")
        .replace("295 0 1", "295 1 0");
    assert_eq!(
        render(load(&INPUT, &sources, 37), 37),
        render(load(&permuted, &sources, 37), 37)
    );
    // Fix the exterior and inlet. Only the downstream injection/wall data change.
    let fixed = INPUT.replace(
        "termination flanged",
        "exterior 300 100000 0\ntermination flanged",
    );
    let changed = fixed
        .replace("295 0 1", "295 1 0")
        .replace("0.0005 290", "0.0005 320");
    let first = load(&fixed, &sources, 37);
    let second = load(&changed, &sources, 37);
    assert_eq!(
        first.outlet_receiver().unwrap().gas,
        second.outlet_receiver().unwrap().gas
    );
    assert_eq!(
        first.gas_field().unwrap().0.state_at(0.0).unwrap(),
        second.gas_field().unwrap().0.state_at(0.0).unwrap()
    );
    assert_ne!(
        first.gas_field().unwrap().0.state_at(0.57).unwrap(),
        second.gas_field().unwrap().0.state_at(0.57).unwrap()
    );
    let (p, motion) = render(first, 37);
    let (q, other_motion) = render(second, 37);
    assert_ne!(
        motion, other_motion,
        "the changed medium must feed back into the reed"
    );
    let energy: f64 = p.iter().map(|x| x * x).sum();
    let difference: f64 = p.iter().zip(&q).map(|(a, b)| (a - b).powi(2)).sum();
    assert!(energy > 0.0 && difference > energy * 1e-12);
}

#[test]
fn thermal_cell_geometry_is_the_exact_authored_frustum_even_across_segments() {
    let segments = [
        Segment::Cylinder {
            radius: 0.002,
            length: 0.25,
        },
        Segment::Cone {
            inlet_radius: 0.002,
            outlet_radius: 0.004,
            length: 0.25,
        },
    ];
    let (v, a) = geometry(&segments, 0.0, 0.5, 1, &mut 0).unwrap();
    let expected_v = core::f64::consts::PI
        * 0.25
        * (0.002_f64.powi(2) + (0.002_f64.powi(2) + 0.002 * 0.004 + 0.004_f64.powi(2)) / 3.0);
    let expected_a = core::f64::consts::TAU * 0.002 * 0.25
        + core::f64::consts::PI * 0.006 * 0.25_f64.hypot(0.002);
    assert!((v - expected_v).abs() < expected_v * 1e-14);
    assert!((a - expected_a).abs() < expected_a * 1e-14);
    let (v0, a0) = geometry(&segments, 0.0, 0.375, 1, &mut 0).unwrap();
    let (v1, a1) = geometry(&segments, 0.375, 0.5, 1, &mut 0).unwrap();
    assert!((v0 + v1 - v).abs() < v * 1e-14);
    assert!((a0 + a1 - a).abs() < a * 1e-14);
}

#[test]
fn smooth_profile_geometry_generates_the_same_field_used_by_pressure_rendering() {
    // G1/G3: all three profile constructors reach the original thermal/species
    // and acoustic owners. An independent analytic exponential-volume check
    // prevents substituting the endpoint cone for the true profile.
    let (a, b) = (packs("N2", 0.028, 3.5, 200.0), packs("CO2", 0.044, 4.5, 200.0));
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    for row in [
        "profile exponential 16 0.0015 0.4 0.57",
        "profile power 16 0.0015 10 1 0.57",
        "profile spline 16 natural 3 0 0.0015 0.285 0.0016 0.57 0.00155",
        "profile spline 16 clamped 0 0 3 0 0.0015 0.285 0.0016 0.57 0.00155",
    ] {
        let input = INPUT.replace("cylinder 0.0015 0.57", row);
        let performance = load(&input, &sources, 37);
        let generated = performance.steady_gas().unwrap();
        let cells = &generated.problem().cells;
        assert_eq!(cells.len(), 2);
        if row.contains("exponential") {
            for (index, cell) in cells.iter().enumerate() {
                let left = index as f64 * 0.285;
                let right = left + 0.285;
                let exact = core::f64::consts::PI * 0.0015_f64.powi(2)
                    * ((0.4 * right).exp() - (0.4 * left).exp()) / 0.4;
                assert!((cell.volume_m3 / exact - 1.0).abs() < 2e-12);
            }
        }
        let direct = ReedPerformance::from_bytes_with_field_resolution(
            without_steady(&input).as_bytes(), 37, ObservationPolicy::Diagnostic,
            Some(&sources), generated.solution().field().clone(),
            performance.gas_field_resolution().unwrap(),
        ).unwrap();
        let expected = render(performance, 37);
        assert!(expected.0.iter().any(|p| p.abs() > 1e-12));
        assert_eq!(expected, render(direct, 37));
        assert_eq!(expected, render(load(&input, &sources, 113), 113));
    }
}

#[test]
fn smooth_profile_cells_keep_exact_junctions_and_share_one_work_budget() {
    // G0/G1: thermal cells may cross profile/ordinary segment boundaries, and
    // the integration cap bounds the complete binding instead of each cell.
    let segments = [
        Segment::Cylinder { radius: 0.002, length: 0.2 },
        Segment::Profile {
            profile: std::sync::Arc::new(fs_duct::RadiusProfile::power(0.002, 0.2, 1.0, 0.3).unwrap()),
            steps: fs_duct::RiccatiSteps(4),
        },
    ];
    let ordinary = [
        segments[0].clone(),
        Segment::Cone { inlet_radius: 0.002, outlet_radius: 0.005, length: 0.3 },
    ];
    let mut work = 0;
    for (start, end) in [(0.1, 0.35), (0.35, 0.5)] {
        let actual = geometry(&segments, start, end, 1, &mut work).unwrap();
        let exact = geometry(&ordinary, start, end, 1, &mut 0).unwrap();
        assert!((actual.0 / exact.0 - 1.0).abs() < 2e-14);
        assert!((actual.1 / exact.1 - 1.0).abs() < 2e-14);
    }
    assert_eq!(work, 8);
    let mut exhausted = fs_duct::profile::MAX_RICCATI_STEPS - 7;
    geometry(&segments, 0.2, 0.35, 1, &mut exhausted).unwrap();
    assert_eq!(exhausted, fs_duct::profile::MAX_RICCATI_STEPS - 3);
    let error = geometry(&segments, 0.35, 0.5, 1, &mut exhausted).unwrap_err();
    assert!(matches!(error, ReedPerformanceError::Input {
        what: "thermal profile geometry exceeds MAX_RICCATI_STEPS", ..
    }));
    assert_eq!(exhausted, fs_duct::profile::MAX_RICCATI_STEPS,
        "the second cell must stop before executing an over-budget panel");
}

#[test]
fn conflicts_unbalanced_flow_and_exhausted_work_refuse_before_rendering() {
    let (a, b) = (
        packs("N2", 0.028, 3.5, 200.0),
        packs("CO2", 0.044, 4.5, 200.0),
    );
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    for input in [
        INPUT.replace("cylinder 0.0015 0.57", "cylinder 0.0015 0.57 gas 300 0"),
        INPUT.replace("gas_cell 0.57", "gas_cell 0.56"),
        INPUT.replace("gas_injection 1", "gas_injection 2"),
        INPUT.replace("mole 2 N2 CO2", "mole 2 N2 N2"),
        INPUT.replace("0.1 2000 100000", "0.1 1 100000"),
        INPUT.replace("gas_inlet 0.000001", "gas_inlet 0.01"),
        INPUT.replace("cylinder 0.0015 0.57", "profile exponential 0 0.0015 0.4 0.57"),
        INPUT.replace("cylinder 0.0015 0.57", "profile spline 16 natural 2 0 0.0015 0 0.0016"),
        INPUT.replace("cylinder 0.0015 0.57", "profile power 16 0.0015 0 1 0.57"),
        INPUT.replace(
            "segments 1\ncylinder 0.0015 0.57",
            "segments 3\ncylinder 0.0015 0.25\nhole 0.001 0.001 0.0015 1\ncylinder 0.0015 0.32",
        ),
    ] {
        assert!(
            ReedPerformance::from_bytes_with_gas_sources(
                input.as_bytes(),
                37,
                ObservationPolicy::Diagnostic,
                Some(&sources)
            )
            .is_err(),
            "{input}"
        );
    }
    assert!(ReedPerformance::from_bytes(INPUT.as_bytes(), 37).is_err());
    let field = load(&INPUT, &sources, 37).gas_field().unwrap().0.clone();
    assert!(
        ReedPerformance::from_bytes_with_gas_field(
            INPUT.as_bytes(),
            37,
            ObservationPolicy::Diagnostic,
            Some(&sources),
            field,
            8
        )
        .is_err()
    );
}
