//! G0/G3: authored material data reach the original thermal and acoustic owners.
use crate::dg_port::input::{DgSourceError, DgSourceInput};
use crate::pcm_wav::observation::ObservationPolicy;
use crate::render::schedule::reed::{ReedPerformance, ReedPerformanceError};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::{Density, Temperature};

#[path = "../../render/schedule/reed/field_input/tests/support.rs"]
mod support;

fn pack(conductivity: f64) -> NormalizedMaterialCardPack {
    let artifact = hash_domain(
        "org.frankensim.tests.steady-wall-input.v1",
        &conductivity.to_le_bytes(),
    );
    let provenance = Provenance {
        source: "manufactured thermal wall; not measured".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic wall".into(),
            method: "manufactured regression".into(),
            artifact,
            caveats: "no measured material or instrument claim".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("cp", SPECIFIC_HEAT_CAPACITY_DIMS, 500.0),
        ("conductivity", THERMAL_CONDUCTIVITY_DIMS, conductivity),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (400.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "synthetic-wall".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        NormalizedPack::new(
            "synthetic-wall",
            "authored-v1",
            artifact,
            "CC0-1.0",
            claims,
            vec![],
            vec![],
        )
        .unwrap(),
    )
    .unwrap()
}

const WALL: &str = concat!(
    "gas_wall_material wall.matpack manufactured\n",
    "gas_wall_source T solid density cp conductivity 200 400 0.000001\n",
    "gas_wall_geometry cylindrical 0.01 0.1",
);
const WALL_SOLVE: &str = "gas_wall_solve 96 0.000000000001 0.0000000000001\n";

fn reed_input() -> String {
    include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/examples/reed-steady-gas.performance"
    ))
    .replace("observer_law stokes-kirchhoff\n", "")
}
fn material_reed_input() -> String {
    reed_input()
        .replace(
            "gas_cell 0.285 0.0005 290 300",
            &format!("gas_cell 0.285 material 290 300\n{WALL}"),
        )
        .replace(
            "gas_cell 0.57 0.0005 290 300",
            &format!("gas_cell 0.57 material 290 300\n{WALL}"),
        )
        .replace("compile_limits", &format!("{WALL_SOLVE}compile_limits"))
}

#[test]
fn g3_material_only_swap_changes_reed_gas_and_original_pressure_renderer() {
    let (a, b) = (
        support::packs("N2", 0.028, 3.5, 200.0),
        support::packs("CO2", 0.044, 4.5, 200.0),
    );
    let bindings = [a.binding(), b.binding()];
    let gases = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let text = material_reed_input();
    let mut results = Vec::new();
    for conductivity in [0.002, 0.02] {
        let material = pack(conductivity);
        let mut loads = 0;
        let mut performance = ReedPerformance::from_bytes_with_sources(
            text.as_bytes(),
            1024,
            ObservationPolicy::Diagnostic,
            Some(&gases),
            |name| {
                assert_eq!(name, "wall.matpack");
                loads += 1;
                Ok(material.clone())
            },
        )
        .unwrap();
        assert_eq!(loads, 1, "one immutable pack for both wall cells");
        let generated = performance.steady_gas().unwrap();
        assert_eq!(generated.material_walls().len(), 2);
        for wall in generated.material_walls() {
            assert_eq!(wall.source.pack().content_hash(), material.content_hash());
            assert_eq!(wall.inner_radius_m, 0.0015);
            assert_eq!(wall.length_m, 0.285);
        }
        assert!(
            generated.solution().maximum_power_residual_bound_w()
                <= generated.controls().steady_power_tolerance_w
        );
        let outlet = generated.solution().field().state_at(0.57).unwrap();
        assert_eq!(outlet.state(), &performance.outlet_receiver().unwrap().gas);
        let mut pressure = [0.0; 1024];
        performance.pressure_block(&mut pressure).unwrap();
        assert!(pressure.iter().all(|value| value.is_finite()));
        results.push((outlet.state().temperature, pressure));
    }
    assert!(results[1].0 < results[0].0);
    assert_ne!(results[1].1, results[0].1);
    // The source-capable API does not invoke any material loader or alter the
    // original arithmetic when every cell keeps its explicit numeric G.
    let mut original = ReedPerformance::from_bytes_with_gas_sources(
        reed_input().as_bytes(),
        1024,
        ObservationPolicy::Diagnostic,
        Some(&gases),
    )
    .unwrap();
    let mut compatible = ReedPerformance::from_bytes_with_sources(
        reed_input().as_bytes(),
        1024,
        ObservationPolicy::Diagnostic,
        Some(&gases),
        |_| panic!("numeric G needs no wall pack"),
    )
    .unwrap();
    let (mut old, mut new) = ([0.0; 1024], [0.0; 1024]);
    original.pressure_block(&mut old).unwrap();
    compatible.pressure_block(&mut new).unwrap();
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));
}

#[test]
fn g3_shared_dg_parser_keeps_wall_sources_after_preparing_the_acoustic_load() {
    let gas = support::packs("N2", 0.028, 3.5, 200.0);
    let bindings = [gas.binding()];
    let gases = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let ambient = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let text = format!(
        concat!(
            "dg_bore 1 open 2 1 8 50 2000 0.01 1e-11 192 4 gas_steady 2 0 mole 1 N2\n",
            "dg_segment cylinder 0.001 0.05\n",
            "gas_inlet 0.000001 330 1\n",
            "gas_cell 0.025 material 290 300\n{WALL}\n",
            "gas_cell 0.05 material 290 300\n{WALL}\n",
            "gas_solve 275 350 0.0000001 64 0.01 5000 100000 0.00000001 0.00000001 0.05\n{WALL_SOLVE}",
        ),
        WALL = WALL,
        WALL_SOLVE = WALL_SOLVE
    );
    let mut results = Vec::new();
    for conductivity in [0.002, 0.02] {
        let material = pack(conductivity);
        let mut rows = text.lines().enumerate();
        let mut fields = rows.next().unwrap().1.split_ascii_whitespace();
        let kind = fields.next().unwrap();
        let source = DgSourceInput::read_with_material_sources(
            kind,
            1,
            &mut fields,
            &mut || {
                let (line, row) = rows.next().expect("complete test input");
                Ok::<_, DgSourceError>((line + 1, row.split_ascii_whitespace()))
            },
            &ambient,
            Some(&gases),
            &mut Vec::new(),
            &CancelGate::new_clock_free(),
            &mut |_| Ok(material.clone()),
        )
        .unwrap();
        assert!(rows.next().is_none());
        let load = source
            .prepare(&ambient, 1.0 / 96000.0, &CancelGate::new_clock_free())
            .unwrap();
        let retained = load
            .steady_gas()
            .expect("prepared load retains the actual thermal source");
        assert_eq!(
            retained.material_walls()[0].source.pack().content_hash(),
            material.content_hash()
        );
        assert_eq!(
            retained.solution().field(),
            load.duct().gas_field().unwrap()
        );
        results.push((
            load.duct().outlet_gas().temperature,
            load.duct().operator_identity(),
        ));
    }
    assert!(results[1].0 < results[0].0);
    assert_ne!(results[1].1, results[0].1);
}

#[test]
fn g0_authored_material_wall_refuses_missing_loader_geometry_and_source_mismatch() {
    let (a, b) = (
        support::packs("N2", 0.028, 3.5, 200.0),
        support::packs("CO2", 0.044, 4.5, 200.0),
    );
    let bindings = [a.binding(), b.binding()];
    let gases = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let text = material_reed_input();
    let missing = ReedPerformance::from_bytes_with_gas_sources(
        text.as_bytes(),
        1,
        ObservationPolicy::Diagnostic,
        Some(&gases),
    );
    assert!(matches!(
        missing,
        Err(ReedPerformanceError::MaterialSource { .. })
    ));
    for invalid in [
        text.replace("cylinder 0.0015 0.57", "cone 0.0015 0.002 0.57"),
        text.replace("wall.matpack manufactured", "wall.matpack other-state"),
        text.replace("T solid density", "T liquid density"),
        text.replace("cp conductivity", "cp absent"),
        text.replace("200 400 0.000001", "300 400 0.000001"),
    ] {
        assert!(
            ReedPerformance::from_bytes_with_sources(
                invalid.as_bytes(),
                1,
                ObservationPolicy::Diagnostic,
                Some(&gases),
                |_| Ok(pack(0.002))
            )
            .is_err()
        );
    }
}
