//! Actual file-to-native-solver checks for source-defined heat anisotropy.
//! The fixture contains no scalar thermal_conductivity claim.
use super::*;
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
};

const KEYS: [&str; 3] = [
    "heat_along_fibre",
    "heat_across_fibre",
    "heat_through_stack",
];
const IDENTITY: &str = "[[1,0,0],[0,1,0],[0,0,1]]";
const SWAP_XZ: &str = "[[0,0,-1],[0,1,0],[1,0,0]]";

fn pack(k: [f64; 3], slopes: [f64; 3], semantic: bool) -> Vec<u8> {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (
            ELECTRICAL_RESISTIVITY_PROPERTY,
            ELECTRICAL_RESISTIVITY_DIMS,
            2.0,
            0.0,
        ),
        (DENSITY_PROPERTY, fs_qty::Density::DIMS, 2.0, 0.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            3.0,
            0.0,
        ),
        (KEYS[0], THERMAL_CONDUCTIVITY_DIMS, k[0], slopes[0]),
        (KEYS[1], THERMAL_CONDUCTIVITY_DIMS, k[1], slopes[1]),
        (KEYS[2], THERMAL_CONDUCTIVITY_DIMS, k[2], slopes[2]),
    ] {
        let key = if semantic && KEYS.contains(&name) {
            PropertyKey::with_quantity(
                name,
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                )),
            )
        } else {
            PropertyKey::new(name, dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (290.0, value * (1.0 - 10.0 * slope)),
                        (350.0, value * (1.0 + 50.0 * slope)),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290.0, 350.0)
                    .with("pressure", 101325.0, 101325.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored orthotropic transport fixture".into(),
                    license: "CC0: authored synthetic test data".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    authored_observation(&mut claims, "authored orthotropic transport fixture");
    let pack = NormalizedPack::new(
        "synthetic-directional",
        "test-v1",
        hash_bytes(b"authored orthotropic transport"),
        "CC0: authored synthetic test data",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "synthetic directional fixture".into(),
            phase: "solid".into(),
            process: "authored analytic law".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
    .to_bytes()
}

fn directional_scene(axes: &str) -> String {
    scene().replace("\"conducting\":true,", &format!(
        "\"conducting\":true,\"thermal_conductivity\":{{\"principal_properties\":[\"{}\",\"{}\",\"{}\"],\"principal_axes\":{axes}}},",
        KEYS[0], KEYS[1], KEYS[2],
    ))
}

fn grid(result: &J, field: &str, cell: usize) -> Vec<f64> {
    result.get(field).unwrap().as_array().unwrap()[cell]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

#[test]
fn principal_source_data_and_frame_change_actual_heat_flow_without_a_scalar_fallback() {
    for semantic in [false, true] {
        for (k, axes, expected_grid, kx) in [
            ([1.0, 4.0, 9.0], IDENTITY, [1.0, 4.0, 9.0], 1.0),
            ([1.0, 4.0, 9.0], SWAP_XZ, [9.0, 4.0, 1.0], 9.0),
            ([2.0, 8.0, 18.0], IDENTITY, [2.0, 8.0, 18.0], 2.0),
        ] {
            let bytes = pack(k, [0.0; 3], semantic);
            let path = fixture(&bytes, &directional_scene(axes));
            let (_, result) = execute(&path);
            // Same 4.5 W source, same geometry. Two 1.5/kx series paths to
            // the boundaries give the independent centre rise 3.375/kx.
            assert!((component(&result, "temperature_k", 1) - (300.0 + 3.375 / kx)).abs() < 1e-8);
            assert!((result.f64_field("total_joule_w").unwrap() - 4.5).abs() < 1e-9);
            assert!(
                result
                    .f64_field("thermal_balance_relative_residual")
                    .unwrap()
                    < 1e-9
            );
            for cell in 0..3 {
                assert_eq!(
                    grid(&result, "resolved_conductivity_grid_w_m_k", cell),
                    expected_grid
                );
                assert_eq!(
                    grid(&result, "solved_conductivity_grid_w_m_k", cell),
                    expected_grid
                );
                assert!(matches!(
                    &result
                        .get("resolved_conductivity_w_m_k")
                        .unwrap()
                        .as_array()
                        .unwrap()[cell],
                    J::Null
                ));
            }
            let declaration = result.get("materials").unwrap().as_array().unwrap()[0]
                .get("thermal_conductivity")
                .unwrap();
            assert_eq!(
                declaration.str_field("kind"),
                Some("grid-aligned-orthotropic")
            );
            assert_eq!(
                declaration.get("principal_axes"),
                Some(&J::parse(axes).unwrap())
            );
            assert_eq!(
                std::fs::read(path.parent().unwrap().join("material.fsmcdpk")).unwrap(),
                bytes
            );
        }
    }
}

#[test]
fn all_directional_laws_are_requeried_at_the_actual_coupled_temperature() {
    let bases = [1.0, 4.0, 9.0];
    let slopes = [0.01, 0.02, 0.03];
    let path = fixture(&pack(bases, slopes, true), &directional_scene(IDENTITY));
    let (_, result) = execute(&path);
    // Independent symmetric network reduction at constant electrical rho:
    // x = 4.5/(4*(1+.01*x));
    // d = 4.5/(2*(1+.01*x)) + 4.5/(4*(1+.01*d)).
    let x = ((1.0_f64 + 0.01 * 4.5).sqrt() - 1.0) / 0.02;
    let left = 4.5 / (2.0 * (1.0 + 0.01 * x));
    let d = (-(1.0 - 0.01 * left) + ((1.0 - 0.01 * left).powi(2) + 0.04 * (left + 1.125)).sqrt())
        / 0.02;
    for (cell, rise) in [x, d, x].into_iter().enumerate() {
        let temperature = component(&result, "temperature_k", cell);
        assert!((temperature - 300.0 - rise).abs() < 1e-7);
        let resolved = grid(&result, "resolved_conductivity_grid_w_m_k", cell);
        for axis in 0..3 {
            assert!(
                (resolved[axis] - bases[axis] * (1.0 + slopes[axis] * (temperature - 300.0))).abs()
                    < 1e-12
            );
        }
    }
    assert!(result.f64_field("conductivity_relative_defect").unwrap() <= 1e-10);
    assert!(
        result
            .f64_field("thermal_balance_relative_residual")
            .unwrap()
            < 1e-9
    );
}

#[test]
fn transient_storage_and_all_three_conductivity_axes_share_the_same_material_state() {
    let mut peaks = Vec::new();
    for (axes, kx) in [(IDENTITY, 1.0), (SWAP_XZ, 9.0)] {
        let scene = directional_scene(axes).replace(
            "{\"temperature_k\":300.0},{\"temperature_k\":300.0}",
            "\"adiabatic\",\"adiabatic\"",
        );
        let bytes = pack([1.0, 4.0, 9.0], [0.0; 3], false);
        let path = fixture(&bytes, &scene);
        let time_path = path.parent().unwrap().join("time.json");
        std::fs::write(
            &time_path,
            r#"{
            "schema":"frankensim.electrothermal.transient.v1",
            "time_step_s":0.1,"steps":1,"temperature_min_k":290.0,"temperature_max_k":350.0,
            "enthalpy_interpolation_error_j_kg":0.00001,
            "energy_absolute_tolerance_j":0.0000001,"energy_relative_tolerance":0.000000001,
            "schedule":[[1,1.0,1.0]]
        }"#,
        )
        .unwrap();
        let output = run(
            &[
                path.as_os_str().to_owned(),
                "--transient".into(),
                time_path.as_os_str().to_owned(),
            ],
            true,
        );
        assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
        let result = J::parse(&output.stdout).unwrap();
        // Independent backward-Euler three-cell system: each mass*Cp/dt=60.
        let d = 4.5 / (60.0 + 2.0 * kx - 2.0 * kx * kx / (60.0 + kx));
        let x = kx / (60.0 + kx) * d;
        for (cell, rise) in [x, d, x].into_iter().enumerate() {
            assert!((component(&result, "temperature_k", cell) - 300.0 - rise).abs() < 1e-8);
            assert_eq!(component(&result, "cell_mass_kg", cell), 2.0);
        }
        let energy: f64 = (0..3)
            .map(|cell| 6.0 * (component(&result, "temperature_k", cell) - 300.0))
            .sum();
        assert!((energy - 0.45).abs() < 1e-7);
        assert_eq!(grid(&result, "resolved_conductivity_grid_w_m_k", 1)[0], kx);
        assert_eq!(grid(&result, "solved_conductivity_grid_w_m_k", 1)[0], kx);
        peaks.push(component(&result, "temperature_k", 1));
    }
    assert!(peaks[0] > peaks[1]);
}

#[test]
fn unsupported_frames_roles_and_missing_sources_refuse_before_publishing_fields() {
    let bytes = pack([1.0, 4.0, 9.0], [0.0; 3], false);
    let valid = directional_scene(IDENTITY);
    for text in [
        directional_scene("[[0.8,0.6,0],[-0.6,0.8,0],[0,0,1]]"),
        directional_scene("[[1,0,0],[1,0,0],[0,0,1]]"),
        directional_scene("[[1,1,0],[0,1,0],[0,0,1]]"),
        directional_scene("[[0,0,0],[0,1,0],[0,0,1]]"),
        valid.replace(KEYS[2], KEYS[0]),
        valid.replace(KEYS[2], "missing_directional_law"),
        valid.replace(KEYS[2], SPECIFIC_HEAT_CAPACITY_PROPERTY),
        scene(), // no isotropic property in this pack; no principal averaging
        valid.replace("\"schema\":", "\"fluid\":{\"material_index\":0,\"reference_temperature_k\":300,\"velocity_m_s\":[0,0,0],\"flow_faces\":[\"wall\",\"wall\",\"wall\",\"wall\",\"wall\",\"wall\"]},\"schema\":"),
    ] {
        let path = fixture(&bytes, &text);
        let output = run(&[path.as_os_str().to_owned()], true);
        assert_eq!(output.exit_code, exit::INPUT, "{}", output.stderr);
        assert!(output.stdout.is_empty());
        assert!(J::parse(&output.stderr).unwrap().str_field("message").unwrap()
            .contains("thermal_conductivity"));
    }
}

#[test]
fn transverse_source_changes_act_on_native_side_face_heat_flux() {
    for k in [[1.0, 4.0, 9.0], [1.0, 8.0, 9.0], [1.0, 4.0, 18.0]] {
        let text = directional_scene(IDENTITY).replace(
            "\"adiabatic\",\"adiabatic\",\"adiabatic\",\"adiabatic\"",
            "{\"temperature_k\":300},{\"temperature_k\":300},{\"temperature_k\":300},{\"temperature_k\":300}",
        );
        let path = fixture(&pack(k, [0.0; 3], false), &text);
        let (_, result) = execute(&path);
        // Four transverse half-cell paths: g = 4*ky + 4*kz. The two
        // endpoint cells also lose heat through their 2*kx end faces.
        let g = 4.0 * (k[1] + k[2]);
        let centre = 4.5 / (2.0 * k[0] + g - 2.0 * k[0] * k[0] / (3.0 * k[0] + g));
        let end = k[0] * centre / (3.0 * k[0] + g);
        for (cell, rise) in [end, centre, end].into_iter().enumerate() {
            assert!((component(&result, "temperature_k", cell) - 300.0 - rise).abs() < 1e-8);
        }
        assert!(
            result
                .f64_field("thermal_balance_relative_residual")
                .unwrap()
                < 1e-9
        );
    }
}
