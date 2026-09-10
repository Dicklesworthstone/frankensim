//! G2/G3 CLI conformance for material, species, and model pack compilation.

#![deny(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    InterpolationPolicy, NormalizedModelPack, NormalizedPack, NormalizedSpeciesPack, PropertyValue,
    SPECIES_MOLAR_MASS_DIMS, SPECIES_PACK_TARGET_BASIS, SPECIES_REFERENCE_PRESSURE_DIMS,
    SpeciesNormalizationTarget, UncertaintyModel,
};
use fs_qty::Dims;

/// Source acquisition checks exercise the real offline compiler and typed
/// evaluator. They establish source transport, not experimental validation.
mod common_material_acquisition {
    use super::*;
    use fs_matdb::{PropertyClaim, QueryPoint, SelectionPolicy};

    pub(super) fn compile(slug: &str) -> (NormalizedPack, PathBuf) {
        let manifest = workspace_path(&format!("data/matdb/seed-v1/{slug}/manifest.tsv"));
        let path = fixture_dir().join(format!("{slug}.fsmatpk"));
        let run = run_compiler(&manifest, &path);
        assert!(
            run.status.success(),
            "{slug}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let bytes = fs::read(&path).unwrap();
        let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
        let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
        assert_eq!(
            pack.schema_version(),
            3,
            "new sources must retain typed axes"
        );
        (pack, path)
    }

    pub(super) fn point(claim: &PropertyClaim, overrides: &[(&str, f64)]) -> QueryPoint {
        for (axis, _) in overrides {
            assert!(
                claim.validity.bound(axis).is_some(),
                "query override {axis} is absent from {} validity",
                claim.key.name()
            );
        }
        let mut point = QueryPoint::new();
        for (axis, &(lo, _)) in claim.validity.bounds() {
            let value = overrides
                .iter()
                .find(|(name, _)| *name == axis)
                .map_or(lo, |(_, value)| *value);
            point = if let Some(quantity) = claim.validity.axis_quantities().get(axis) {
                point.with_quantity(axis, *quantity, value)
            } else {
                point.with(axis, value)
            }
            .unwrap();
        }
        point
    }

    pub(super) fn sample(pack: &NormalizedPack, name: &str, overrides: &[(&str, f64)]) -> f64 {
        let claims = pack.claims().claims_for(name);
        let claim = claims
            .first()
            .unwrap_or_else(|| panic!("missing canonical property {name}"))
            .1;
        pack.claims()
            .query_typed(
                &claim.key,
                &point(claim, overrides),
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap_or_else(|error| panic!("{name} {overrides:?}: {error:?}"))
            .evidence
            .value
            .value
    }

    pub(super) fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-10 * expected.abs().max(1.0),
            "actual {actual}, source-derived expected {expected}"
        );
    }

    /// G1/G3: independently sourced cubic constants reach the actual oriented
    /// tetrahedral operator. This is a cross-source engineering reference at
    /// 25 C, not qualification of a particular silicon wafer or its doping.
    #[test]
    fn g1_g3_sourced_silicon_tensor_reaches_oriented_solid() {
        use fs_alloc::{ArenaConfig, ArenaPool};
        use fs_blake3::ContentHash;
        use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, PropertyKey};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, resolve_elastic_tensor_state_point,
        };
        use fs_solid::{
            TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
        };

        let manifest =
            workspace_path("data/matdb/seed-v1/silicon-cubic-25c-nasa-rp1057/manifest.tsv");
        let scratch = fixture_dir();
        let output = scratch.join("silicon.fsmatpk");
        let run = run_compiler(&manifest, &output);
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
        let bytes = fs::read(output).unwrap();
        let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
        let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
        assert_eq!(
            pack.schema_version(),
            5,
            "tensor coordinates must survive compilation"
        );
        let card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "Si; stiffness dopant/purity unspecified".into(),
                phase: "cubic crystal".into(),
                process: "cross-source engineering reference: handbook stiffness and separate pure-crystal density; not a qualified wafer".into(),
                revision: 0,
            },
            pack,
        ).unwrap();
        let database = scratch.join("silicon.sqlite");
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[CatalogPack::MaterialCard(card.clone())])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(loaded) = store.load_catalog_pack(card.pack_id()).unwrap()
        else {
            panic!("wrong stored family")
        };
        assert_eq!(loaded, card);
        let claims = loaded.card().claims();
        let mut entries: [[Option<PropertyKey>; 6]; 6] =
            core::array::from_fn(|_| core::array::from_fn(|_| None));
        for (_, claim) in claims.claims_ordered() {
            if let Some(component) = claim.key.elastic_component() {
                let (row, column) = component.indices();
                assert!(entries[row][column].replace(claim.key.clone()).is_none());
                assert!(matches!(claim.uncertainty, UncertaintyModel::Unstated));
            }
        }
        let keys = entries
            .map(|row| row.map(|entry| entry.expect("explicit tensor entry, including zeros")));
        let density = claims.claims_for("density")[0].1;
        let at = point(density, &[("T", 298.15)]);
        let resolve = |at: &QueryPoint| {
            resolve_elastic_tensor_state_point(
                loaded.card(),
                at,
                &keys,
                MaterialPropertySelection::SingleClaimOnly,
            )
        };
        let state = resolve(&at).unwrap();
        assert_eq!(state.resolved().properties().len(), 37);
        close(state.density_kg_m3(), 2329.0);
        for property in state.resolved().properties() {
            claims.verify_receipt(&property.answer().receipt).unwrap();
        }
        // Independent printed constants, not values read back from the tested
        // matrix. RP-1057's Mbar convention is checked by its 1.012/Mbar bulk
        // compressibility on p.100; the paired compliances are rounded.
        let (c11, c12, c44) = (165.773e9, 63.924e9, 79.619e9);
        assert!((3.0e11 / (c11 + 2.0 * c12) - 1.012_f64).abs() < 0.011);
        for row in 0..6 {
            for column in 0..6 {
                let expected = if row < 3 && column < 3 {
                    if row == column { c11 } else { c12 }
                } else if row == column {
                    c44
                } else {
                    0.0
                };
                close(state.stiffness_pa()[row][column], expected);
            }
        }
        for rejected in [
            point(density, &[("T", 298.1501)]),
            point(density, &[("source-pressure-known", 1.0)]),
            QueryPoint::new()
                .with_quantity("T", at.axis_quantities()["T"], 298.15)
                .unwrap(),
        ] {
            assert!(
                resolve(&rejected).is_err(),
                "unsupported source state: {rejected:?}"
            );
        }
        let mut misaddressed = keys.clone();
        misaddressed[0][3] = keys[0][4].clone();
        assert!(
            resolve_elastic_tensor_state_point(
                loaded.card(),
                &at,
                &misaddressed,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err(),
            "zero-valued components still need their own coordinates"
        );

        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let epsilon = 1.0e-4;
        let displacement = [
            0.0, 0.0, 0.0, epsilon, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let gradients = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut energies = Vec::new();
        for angle in [0.0_f64, std::f64::consts::FRAC_PI_4] {
            let (s, c) = angle.sin_cos();
            let q = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
            let transformed = TetElasticMaterial::from_resolved_elastic_tensor(
                &state,
                ContentHash([0x51; 32]),
                q,
            )
            .unwrap();
            assert_eq!(
                transformed.source_material_identity(),
                state.resolved().identity()
            );
            let problem = TetLinearElasticProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                materials: TetMaterialField::Uniform(transformed.material()),
                fixed_dofs: &[],
                budget: TetAssemblyBudget::standard(),
            };
            let gate = CancelGate::new();
            let pool = ArenaPool::new(ArenaConfig::default());
            let assembly = pool.scope(|arena| {
                let cx = Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 1,
                        kernel_id: 10,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                );
                let assembled = problem.assemble(&cx).unwrap();
                assert_eq!(
                    assembled.stiffness.to_dense(),
                    problem.assemble(&cx).unwrap().stiffness.to_dense()
                );
                assembled
            });
            close(assembly.total_mass_kg, 2329.0 / 6.0);
            // Independently rotate the imposed physical strain into crystal
            // coordinates, apply the three-constant cubic law, rotate stress
            // back, and integrate the P1 shape gradients over volume 1/6 m^3.
            let local_stress: [[f64; 3]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    let strain = epsilon * q[0][i] * q[0][j];
                    if i == j {
                        c12 * epsilon + (c11 - c12) * strain
                    } else {
                        2.0 * c44 * strain
                    }
                })
            });
            let stress: [[f64; 3]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    (0..3)
                        .flat_map(|a| (0..3).map(move |b| q[i][a] * local_stress[a][b] * q[j][b]))
                        .sum()
                })
            });
            let k = assembly.stiffness.to_dense();
            let force: [f64; 12] =
                core::array::from_fn(|i| (0..12).map(|j| k[i * 12 + j] * displacement[j]).sum());
            for (node, gradient) in gradients.iter().enumerate() {
                for component in 0..3 {
                    let expected: f64 = (0..3)
                        .map(|j| stress[component][j] * gradient[j] / 6.0)
                        .sum();
                    assert!(
                        (force[3 * node + component] - expected).abs() < 1e-7,
                        "angle={angle} node={node} component={component}: force={} expected={expected} N tolerance=1e-7 N",
                        force[3 * node + component]
                    );
                }
            }
            let energy: f64 = displacement
                .iter()
                .zip(force)
                .map(|(u, f)| 0.5 * u * f)
                .sum();
            let expected_energy = epsilon * stress[0][0] / 12.0;
            close(energy, expected_energy);
            assert!(energy > 0.0);
            println!(
                "silicon source={} state={} T=298.15 K pressure=unknown angle={angle} rad mass={} kg energy={energy} J expected={expected_energy} J force_tolerance=1e-7 N selected_claims=37 empirical_qualification=false",
                loaded.pack_id(),
                state.resolved().identity(),
                assembly.total_mass_kg
            );
            energies.push(energy);
        }
        assert!(
            energies[1] > 1.17 * energies[0],
            "the actual oriented operator must change response"
        );
    }

    /// G1/G3: an exact dry-air source condition survives compilation and
    /// persistent storage to drive the actual acoustic cylinder-loss law.
    /// This checks the declared ideal-gas/Sutherland/USSA model, not ambient
    /// measurements, weather, humidity physics, or an experimental loss test.
    #[test]
    fn g1_g3_sourced_dry_air_reaches_acoustic_loss() {
        use fs_couple::air_path::oscillating_cylinder_air_resistance_per_length;
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::gas::{ConductivityModel, resolve_sutherland_gas_state};
        use fs_material::state_point::MaterialPropertySelection;
        use fs_qty::QuantitySpec;
        use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

        let (pack, pack_path) = compile("air-dry-ussa1976");
        assert_eq!(pack.pack_id(), "air-dry-ussa1976-ambient-model");
        let card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "dry air; USSA-1976 reference composition".into(),
                phase: "gas".into(),
                process: "USSA-1976 dry-air constants and transport fits".into(),
                revision: 0,
            },
            pack.clone(),
        )
        .unwrap();
        let database = fixture_dir().join("dry-air-ussa1976.sqlite");
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[CatalogPack::MaterialCard(card.clone())])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(loaded) = store.load_catalog_pack(card.pack_id()).unwrap()
        else {
            panic!("wrong stored family")
        };
        assert_eq!(loaded, card);
        let claims = loaded.card().claims();
        let parameters = [
            "molar_mass",
            "heat_capacity_ratio",
            "sutherland_reference_viscosity",
            "sutherland_reference_temperature",
            "sutherland_temperature",
        ];
        let pins: Vec<_> = parameters
            .iter()
            .map(|name| {
                let available = claims.claims_for(name);
                assert_eq!(
                    available.len(),
                    1,
                    "{name} must have one exact source claim"
                );
                available[0].0
            })
            .collect();
        let discovery = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            workspace_path("examples/material-discovery/dry-air.json")
                .to_str()
                .unwrap()
                .into(),
            pack_path.to_str().unwrap().into(),
        ]);
        assert_eq!(
            discovery.exit_code,
            fs_cli::exit::SUCCESS,
            "{}",
            discovery.stderr
        );
        assert!(
            discovery.stdout.contains("\"status\":\"complete\""),
            "{}",
            discovery.stdout
        );
        println!("dry-air units rho=kg/m3 c=m/s mu=Pa*s k=W/m/K Cp=J/kg/K R_loss=N*s/m2");

        let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        let dimensionless = QuantitySpec::dimensional(Dims::NONE);
        let point = |temperature_k: f64, pressure_pa: f64| {
            QueryPoint::new()
                .with_quantity("temperature", absolute_temperature, temperature_k)
                .unwrap()
                .with_quantity("pressure", pressure, pressure_pa)
                .unwrap()
                .with_quantity("relative-humidity", dimensionless, 0.0)
                .unwrap()
                .with_quantity("source-composition-ussa1976", dimensionless, 1.0)
                .unwrap()
        };
        let resolve = |at: &QueryPoint, selection| {
            resolve_sutherland_gas_state(
                loaded.card(),
                at,
                ConductivityModel::Ussa1976AirFit,
                selection,
            )
        };
        let relative_close = |actual: f64, expected: f64| {
            assert!(
                (actual - expected).abs() <= 1.0e-10 * expected.abs(),
                "actual {actual}, source-derived expected {expected}, relative_tolerance=1e-10"
            );
        };
        let source_observations = claims.observation_ids().collect::<Vec<_>>();
        let reference_temperature = 273.15_f64;
        let sutherland_temperature = 110.4_f64;
        let reference_viscosity = 1.458e-6 * reference_temperature.powf(1.5)
            / (reference_temperature + sutherland_temperature);
        let expected = |temperature_k: f64, pressure_pa: f64| {
            let molar_mass = 0.028_964_4_f64;
            let gamma = 1.4_f64;
            let r = 8.314_32_f64 / molar_mass;
            let beta = reference_viscosity * (reference_temperature + sutherland_temperature)
                / reference_temperature.powf(1.5);
            let temperature_32 = temperature_k.powf(1.5);
            let density = pressure_pa / (r * temperature_k);
            let sound_speed = (gamma * r * temperature_k).sqrt();
            let viscosity = beta * temperature_32 / (temperature_k + sutherland_temperature);
            let cp = gamma * r / (gamma - 1.0);
            let conductivity = 2.646_38e-3 * temperature_32
                / (temperature_k + 245.4 * 10.0_f64.powf(-12.0 / temperature_k));
            (
                density,
                sound_speed,
                viscosity,
                cp,
                conductivity,
                viscosity * cp / conductivity,
            )
        };
        let radius_m = 0.0005_f64;
        let omega_rad_s = core::f64::consts::TAU * 220.0;
        let mut resistance = Vec::new();
        for (temperature_k, pressure_pa) in [
            (273.15, 80_000.0),
            (288.15, 101_325.0),
            (293.15, 101_325.0),
            (313.15, 110_000.0),
        ] {
            let at = point(temperature_k, pressure_pa);
            let resolved = resolve(&at, MaterialPropertySelection::SingleClaimOnly).unwrap();
            let replay = resolve(&at, MaterialPropertySelection::SingleClaimOnly).unwrap();
            assert_eq!(resolved, replay, "deterministic gas replay at {at:?}");
            assert_eq!(
                resolved.conductivity_model(),
                ConductivityModel::Ussa1976AirFit
            );
            assert_eq!(
                resolved.parameters().query_point(),
                replay.parameters().query_point()
            );
            assert_eq!(resolved.parameters().properties().len(), parameters.len());
            for property in resolved.parameters().properties() {
                claims.verify_receipt(&property.answer().receipt).unwrap();
            }
            let state = resolved.state();
            let (rho, c, mu, cp, conductivity, prandtl) = expected(temperature_k, pressure_pa);
            relative_close(state.density, rho);
            relative_close(state.sound_speed, c);
            relative_close(state.dynamic_viscosity, mu);
            relative_close(state.specific_heat_cp, cp);
            relative_close(state.thermal_conductivity, conductivity);
            relative_close(state.prandtl, prandtl);
            relative_close(state.characteristic_impedance, rho * c);
            let cv = state.specific_heat_cp - state.specific_gas_constant;
            relative_close(state.specific_heat_cp / cv, 1.4);
            let actual_resistance =
                oscillating_cylinder_air_resistance_per_length(radius_m, omega_rad_s, state)
                    .unwrap();
            let expected_resistance = core::f64::consts::TAU
                * mu
                * (1.0 + radius_m * (2.0 * rho * omega_rad_s / mu).sqrt());
            relative_close(actual_resistance, expected_resistance);
            assert_eq!(
                actual_resistance.to_bits(),
                oscillating_cylinder_air_resistance_per_length(
                    radius_m,
                    omega_rad_s,
                    replay.state()
                )
                .unwrap()
                .to_bits()
            );
            println!(
                "dry-air source_artifact={} observations={:?} query={:?} T={temperature_k} K p={pressure_pa} Pa rho_expected={rho:.9} rho_actual={:.9} c_expected={c:.9} c_actual={:.9} mu_expected={mu:.12e} mu_actual={:.12e} k_expected={conductivity:.12e} k_actual={:.12e} Pr_expected={prandtl:.9} Pr_actual={:.9} resistance_expected={expected_resistance:.12e} resistance_actual={actual_resistance:.12e} relative_tolerance=1e-10 source_receipts={:?}",
                pack.source_artifact(),
                claims.observation_ids().collect::<Vec<_>>(),
                resolved.parameters().query_point(),
                state.density,
                state.sound_speed,
                state.dynamic_viscosity,
                state.thermal_conductivity,
                state.prandtl,
                resolved
                    .parameters()
                    .properties()
                    .iter()
                    .map(|property| property.answer().receipt.selected)
                    .collect::<Vec<_>>(),
            );
            resistance.push(actual_resistance);
        }
        assert!(
            (resistance[1] - resistance[2]).abs() > 0.0,
            "temperature changes air loss"
        );
        let low_pressure = resolve(
            &point(288.15, 80_000.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let high_pressure = resolve(
            &point(288.15, 110_000.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            oscillating_cylinder_air_resistance_per_length(
                radius_m,
                omega_rad_s,
                high_pressure.state()
            )
            .unwrap()
                > oscillating_cylinder_air_resistance_per_length(
                    radius_m,
                    omega_rad_s,
                    low_pressure.state()
                )
                .unwrap(),
            "pressure changes the actual acoustic-loss result"
        );
        let sea_level = resolve(
            &point(288.15, 101_325.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!((sea_level.state().density - 1.2250).abs() < 5.0e-5);
        assert!((sea_level.state().sound_speed - 340.29).abs() < 5.0e-3);
        assert!((sea_level.state().dynamic_viscosity - 1.7894e-5).abs() < 5.0e-9);
        assert!(
            (sea_level.state().thermal_conductivity - 0.025_325_884_264_263_953).abs() < 1.0e-14
        );

        let compile_synthetic = |suffix: &str, properties: String| {
            let source_dir = workspace_path("data/matdb/seed-v1/air-dry-ussa1976");
            let synthetic_dir = fixture_dir().join(format!("dry-air-{suffix}"));
            fs::create_dir_all(&synthetic_dir).unwrap();
            fs::copy(
                source_dir.join("manifest.tsv"),
                synthetic_dir.join("manifest.tsv"),
            )
            .unwrap();
            fs::write(synthetic_dir.join("properties.tsv"), properties).unwrap();
            let output = synthetic_dir.join("synthetic.fsmatpk");
            let run = run_compiler(&synthetic_dir.join("manifest.tsv"), &output);
            assert!(
                run.status.success(),
                "synthetic {suffix}: {}",
                String::from_utf8_lossy(&run.stderr)
            );
            let bytes = fs::read(output).unwrap();
            let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
            NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap()
        };
        let source_properties = fs::read_to_string(workspace_path(
            "data/matdb/seed-v1/air-dry-ussa1976/properties.tsv",
        ))
        .unwrap();
        let causal_properties = source_properties.replacen(
            "0.000017160792662455268875",
            "0.0000188768719287007957625",
            1,
        );
        assert_ne!(causal_properties, source_properties);
        let causal_pack = compile_synthetic("synthetic-mu-reference", causal_properties);
        let causal_card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "dry air; USSA-1976 reference composition".into(),
                phase: "gas".into(),
                process: "synthetic mu_ref causal control; not a source claim".into(),
                revision: 0,
            },
            causal_pack,
        )
        .unwrap();
        let causal_state = resolve_sutherland_gas_state(
            causal_card.card(),
            &point(288.15, 101_325.0),
            ConductivityModel::Ussa1976AirFit,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let base_resistance = oscillating_cylinder_air_resistance_per_length(
            radius_m,
            omega_rad_s,
            sea_level.state(),
        )
        .unwrap();
        let causal_resistance = oscillating_cylinder_air_resistance_per_length(
            radius_m,
            omega_rad_s,
            causal_state.state(),
        )
        .unwrap();
        assert!(causal_state.state().dynamic_viscosity > sea_level.state().dynamic_viscosity);
        assert!(causal_resistance > base_resistance);
        println!(
            "dry-air synthetic_control=mu_ref_changed source_claim=false baseline_resistance={base_resistance:.12e} changed_resistance={causal_resistance:.12e}"
        );
        let wrong_unit_properties =
            source_properties.replacen("\tkg/mol\tconstant", "\t1\tconstant", 1);
        assert_ne!(wrong_unit_properties, source_properties);
        let wrong_unit_pack = compile_synthetic("wrong-molar-mass-unit", wrong_unit_properties);
        let wrong_unit_card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "dry air; synthetic wrong molar-mass unit".into(),
                phase: "gas".into(),
                process: "synthetic unit refusal control; not a source claim".into(),
                revision: 0,
            },
            wrong_unit_pack,
        )
        .unwrap();
        assert!(
            resolve_sutherland_gas_state(
                wrong_unit_card.card(),
                &point(288.15, 101_325.0),
                ConductivityModel::Ussa1976AirFit,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err()
        );
        assert_eq!(
            claims.observation_ids().collect::<Vec<_>>(),
            source_observations
        );

        let pinned: Vec<_> = parameters
            .iter()
            .zip(&pins)
            .map(|(name, id)| ((*name).to_owned(), *id))
            .collect();
        let pinned_state = resolve(
            &point(288.15, 101_325.0),
            MaterialPropertySelection::PinnedByProperty(pinned.clone()),
        )
        .unwrap();
        assert_eq!(pinned_state.state(), sea_level.state());
        assert_eq!(
            pinned_state
                .parameters()
                .properties()
                .iter()
                .map(|property| property.answer().receipt.selected)
                .collect::<Vec<_>>(),
            sea_level
                .parameters()
                .properties()
                .iter()
                .map(|property| property.answer().receipt.selected)
                .collect::<Vec<_>>()
        );
        assert!(
            resolve(
                &point(288.15, 101_325.0),
                MaterialPropertySelection::PinnedByProperty(pinned[..4].to_vec())
            )
            .is_err()
        );
        let incomplete = QueryPoint::new()
            .with_quantity("pressure", pressure, 101_325.0)
            .unwrap()
            .with_quantity("relative-humidity", dimensionless, 0.0)
            .unwrap()
            .with_quantity("source-composition-ussa1976", dimensionless, 1.0)
            .unwrap();
        assert!(resolve(&incomplete, MaterialPropertySelection::SingleClaimOnly).is_err());
        let wrong_temperature = QueryPoint::new()
            .with("temperature", 288.15)
            .unwrap()
            .with_quantity("pressure", pressure, 101_325.0)
            .unwrap()
            .with_quantity("relative-humidity", dimensionless, 0.0)
            .unwrap()
            .with_quantity("source-composition-ussa1976", dimensionless, 1.0)
            .unwrap();
        assert!(
            resolve(
                &wrong_temperature,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
        let missing_pressure = QueryPoint::new()
            .with_quantity("temperature", absolute_temperature, 288.15)
            .unwrap()
            .with_quantity("relative-humidity", dimensionless, 0.0)
            .unwrap()
            .with_quantity("source-composition-ussa1976", dimensionless, 1.0)
            .unwrap();
        assert!(
            resolve(
                &missing_pressure,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
        let wrong_pressure = QueryPoint::new()
            .with_quantity("temperature", absolute_temperature, 288.15)
            .unwrap()
            .with_quantity("pressure", absolute_temperature, 101_325.0)
            .unwrap()
            .with_quantity("relative-humidity", dimensionless, 0.0)
            .unwrap()
            .with_quantity("source-composition-ussa1976", dimensionless, 1.0)
            .unwrap();
        assert!(resolve(&wrong_pressure, MaterialPropertySelection::SingleClaimOnly).is_err());
        for refused in [
            point(313.151, 101_325.0),
            point(288.15, 110_000.1),
            QueryPoint::new()
                .with_quantity("temperature", absolute_temperature, 288.15)
                .unwrap()
                .with_quantity("pressure", pressure, 101_325.0)
                .unwrap()
                .with_quantity("relative-humidity", dimensionless, 1.0)
                .unwrap()
                .with_quantity("source-composition-ussa1976", dimensionless, 1.0)
                .unwrap(),
            QueryPoint::new()
                .with_quantity("temperature", absolute_temperature, 288.15)
                .unwrap()
                .with_quantity("pressure", pressure, 101_325.0)
                .unwrap()
                .with_quantity("source-composition-ussa1976", dimensionless, 1.0)
                .unwrap(),
            QueryPoint::new()
                .with_quantity("temperature", absolute_temperature, 288.15)
                .unwrap()
                .with_quantity("pressure", pressure, 101_325.0)
                .unwrap()
                .with_quantity("relative-humidity", dimensionless, 0.0)
                .unwrap()
                .with_quantity("source-composition-ussa1976", dimensionless, 0.0)
                .unwrap(),
            QueryPoint::new()
                .with_quantity("temperature", absolute_temperature, 288.15)
                .unwrap()
                .with_quantity("pressure", pressure, 101_325.0)
                .unwrap()
                .with_quantity("relative-humidity", dimensionless, 0.0)
                .unwrap(),
        ] {
            assert!(resolve(&refused, MaterialPropertySelection::SingleClaimOnly).is_err());
        }
    }

    /// G1/G3: two explicitly different frozen source states drive the actual
    /// thermoelastic plate/radiator path. Their reference temperatures differ,
    /// so this is not a controlled same-temperature grade substitution.
    #[test]
    fn g1_g3_sourced_metal_profiles_reach_thermoelastic_plate() {
        use fs_couple::thin_plate::{
            PlateThicknessConstraint, certified_radiators,
            with_uniform_isotropic_thermoelastic_material_state,
        };
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, resolve_isotropic_thermoelastic_state_point,
        };
        use fs_scenario::ThinPlate;

        let (aluminum_pack, aluminum_path) = compile("aluminum-2024-t3-nasa-tn-d6448");
        let (stainless_pack, stainless_path) = compile("stainless-316-20c-engineering-reference");
        let aluminum_card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "Al 2024-T3".into(),
                phase: "solid".into(),
                process: "NASA TN D-6448 panel property set".into(),
                revision: 0,
            },
            aluminum_pack,
        )
        .unwrap();
        let stainless_card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "316 stainless steel / UNS S31600".into(),
                phase: "solid".into(),
                process: "cross-source grade-level 20 C engineering reference".into(),
                revision: 0,
            },
            stainless_pack,
        )
        .unwrap();
        let database = fixture_dir().join("metal-thermoelastic-plate.sqlite");
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[
                    CatalogPack::MaterialCard(aluminum_card.clone()),
                    CatalogPack::MaterialCard(stainless_card.clone()),
                ])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(aluminum) =
            store.load_catalog_pack(aluminum_card.pack_id()).unwrap()
        else {
            panic!("wrong stored aluminum family")
        };
        let CatalogPack::MaterialCard(stainless) =
            store.load_catalog_pack(stainless_card.pack_id()).unwrap()
        else {
            panic!("wrong stored stainless family")
        };
        assert_eq!(aluminum, aluminum_card);
        assert_eq!(stainless, stainless_card);
        let property_names = [
            "density",
            "young_modulus",
            "poisson_ratio",
            "linear_thermal_expansion_coefficient",
            "specific_heat_capacity",
            "thermal_conductivity",
        ];
        let template = ThinPlate {
            length_m: 0.20,
            width_m: 0.15,
            thickness_m: 0.0016,
            density_kg_m3: 1.0,
            e1_pa: 1.0e9,
            e2_pa: 1.0e9,
            nu12: 0.3,
            g12_pa: 1.0e9 / 2.6,
            material_angle_rad: 0.0,
            damping_ratio: 0.0,
            thermoelastic: None,
            kelvin_voigt_bending: None,
            n_modes: 1,
            geometric_nonlinearity: false,
            pretension_n_m: 0.0,
            clamped: false,
        };
        let thickness_m = 0.0016;
        let relative_close = |actual: f64, expected: f64| {
            assert!(
                (actual - expected).abs() <= 1.0e-10 * expected.abs(),
                "actual={actual} expected={expected} relative_tolerance=1e-10"
            );
        };
        let inspect = |label: &str,
                       card: &NormalizedMaterialCardPack,
                       pack_path: &Path,
                       discovery: &Path,
                       temperature_k: f64,
                       expected: [f64; 6],
                       requires_pressure_known_axis: bool| {
            let discovery = fs_cli::run(vec![
                "--json".into(),
                "discover".into(),
                discovery.to_str().unwrap().into(),
                pack_path.to_str().unwrap().into(),
            ]);
            assert_eq!(
                discovery.exit_code,
                fs_cli::exit::SUCCESS,
                "{label}: {}",
                discovery.stderr
            );
            assert!(
                discovery.stdout.contains("\"status\":\"complete\""),
                "{label}: {}",
                discovery.stdout
            );
            let claims = card.card().claims();
            let density = claims.claims_for("density")[0].1;
            let at = point(density, &[("T", temperature_k)]);
            assert_eq!(at.axes().get("T"), Some(&temperature_k));
            if requires_pressure_known_axis {
                assert_eq!(at.axes().get("source-pressure-known"), Some(&0.0));
            }
            let state = resolve_isotropic_thermoelastic_state_point(
                card.card(),
                &at,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            assert_eq!(state.resolved().properties().len(), property_names.len());
            let pins: Vec<_> = property_names
                .iter()
                .zip(expected)
                .map(|(name, expected)| {
                    let claims_for_property = claims.claims_for(name);
                    assert_eq!(claims_for_property.len(), 1, "{label} {name}");
                    let (id, claim) = claims_for_property[0];
                    if requires_pressure_known_axis {
                        assert_eq!(claim.observations.len(), 2, "{label} {name}");
                        assert!(
                            claim.observations.iter().any(|observation| {
                                claims.observation(*observation).is_some_and(|dataset| {
                                    dataset.method.contains("cross-source")
                                        && dataset.caveats.contains("not a common coupon")
                                })
                            }),
                            "{label} {name} must retain its compatibility caveat"
                        );
                    }
                    let actual = state.resolved().property(name).unwrap().value_si();
                    assert!(
                        (actual - expected).abs() <= 1.0e-10 * expected.abs(),
                        "{label} {name}: actual={actual} expected={expected}"
                    );
                    claims
                        .verify_receipt(&state.resolved().property(name).unwrap().answer().receipt)
                        .unwrap();
                    ((*name).to_owned(), id)
                })
                .collect();
            let specimen = with_uniform_isotropic_thermoelastic_material_state(
                template,
                &state,
                PlateThicknessConstraint::FixedThickness(thickness_m),
            )
            .unwrap();
            assert_eq!(specimen.material().identity(), state.resolved().identity());
            let plate = specimen.plate();
            let [rho, e, nu, alpha, cp, conductivity] = expected;
            assert_eq!(plate.density_kg_m3, rho);
            assert_eq!(plate.e1_pa, e);
            assert_eq!(plate.e2_pa, e);
            assert_eq!(plate.nu12, nu);
            assert_eq!(plate.g12_pa, e / (2.0 * (1.0 + nu)));
            let thermal = plate.thermoelastic.unwrap();
            assert_eq!(thermal.temperature_k, temperature_k);
            assert_eq!(thermal.linear_expansion_per_k, alpha);
            assert_eq!(thermal.specific_heat_j_kg_k, cp);
            assert_eq!(thermal.conductivity_w_m_k, conductivity);
            assert_eq!(thermal.state_identity, Some(state.resolved().identity()));
            close(
                specimen.mass_kg(),
                rho * thickness_m * plate.length_m * plate.width_m,
            );
            let bending_rigidity = e * thickness_m.powi(3) / (12.0 * (1.0 - nu * nu));
            let continuum_omega11 = core::f64::consts::PI.powi(2)
                * (bending_rigidity / (rho * thickness_m)).sqrt()
                * (1.0 / plate.length_m.powi(2) + 1.0 / plate.width_m.powi(2));
            let mut radiators = certified_radiators(plate).unwrap();
            let body = radiators.remove(0);
            assert!(
                (body.omega / continuum_omega11 - 1.0).abs() < 0.20,
                "{label}: DKT omega={} continuum omega11={continuum_omega11}",
                body.omega
            );
            let tau =
                thickness_m.powi(2) * rho * cp / (core::f64::consts::PI.powi(2) * conductivity);
            let delta = e * alpha.powi(2) * temperature_k / (rho * cp);
            let expected_zeta = 0.5 * delta * body.omega * tau / (1.0 + (body.omega * tau).powi(2));
            relative_close(body.zeta, expected_zeta);
            let replay = certified_radiators(plate).unwrap().remove(0);
            assert_eq!(body.omega.to_bits(), replay.omega.to_bits());
            assert_eq!(body.zeta.to_bits(), replay.zeta.to_bits());

            let total_force_n = 0.1;
            let dt_s = 1.0e-5;
            let mut damped = body.clone();
            let mut replay_damped = body.clone();
            let mut no_thermal_plate = plate;
            no_thermal_plate.thermoelastic = None;
            let mut undamped = certified_radiators(no_thermal_plate).unwrap().remove(0);
            assert_eq!(undamped.zeta, 0.0);
            let mut pressure_difference = false;
            let mut final_pressure = 0.0;
            let mut max_damped_undamped_delta_pa = 0.0_f64;
            for step in 0..128 {
                let generalized_force = if step == 0 {
                    total_force_n * damped.drive_participation
                } else {
                    0.0
                };
                let pressure = damped
                    .drive_and_radiate(generalized_force, dt_s, 1.2, 1.0)
                    .unwrap();
                let replay_pressure = replay_damped
                    .drive_and_radiate(generalized_force, dt_s, 1.2, 1.0)
                    .unwrap();
                let undamped_pressure = undamped
                    .drive_and_radiate(generalized_force, dt_s, 1.2, 1.0)
                    .unwrap();
                assert!(pressure.is_finite());
                final_pressure = pressure;
                max_damped_undamped_delta_pa =
                    max_damped_undamped_delta_pa.max((pressure - undamped_pressure).abs());
                if step == 0 {
                    let expected_initial_pressure = 1.2 * body.area_m2 * generalized_force
                        / body.mass_kg
                        / (2.0 * core::f64::consts::PI);
                    relative_close(pressure, expected_initial_pressure);
                }
                assert_eq!(pressure.to_bits(), replay_pressure.to_bits());
                pressure_difference |= pressure.to_bits() != undamped_pressure.to_bits();
            }
            assert!(
                pressure_difference,
                "{label}: the resolved thermoelastic zeta must affect evolved pressure"
            );
            assert!(max_damped_undamped_delta_pa > 0.0);
            println!(
                "metal-plate source={} T={temperature_k} K source_pressure_known_axis_present={} source_pressure_known_value={:?} rho={rho} kg/m3 E={e} Pa nu={nu} alpha={alpha} 1/K Cp={cp} J/kg/K k={conductivity} W/m/K mass={} kg omega={} rad/s continuum_omega11={continuum_omega11} rad/s zeta={} final_pressure={final_pressure} Pa max_damped_undamped_delta={max_damped_undamped_delta_pa} Pa source_receipts={:?} different_reference_temperatures=true ambient_density_for_radiation=1.2kg/m3",
                card.pack_id(),
                requires_pressure_known_axis,
                at.axes().get("source-pressure-known"),
                specimen.mass_kg(),
                body.omega,
                body.zeta,
                state
                    .resolved()
                    .properties()
                    .iter()
                    .map(|property| property.answer().receipt.selected)
                    .collect::<Vec<_>>(),
            );
            (state, pins, at, final_pressure)
        };
        let (aluminum_state, aluminum_pins, aluminum_point, aluminum_pressure) = inspect(
            "aluminum",
            &aluminum,
            &aluminum_path,
            &workspace_path(
                "examples/material-discovery/aluminum-2024-t3-thermoelastic-reference.json",
            ),
            300.0,
            [2705.0, 72.4e9, 0.33, 23.0e-6, 840.0, 126.0],
            false,
        );
        let (stainless_state, _stainless_pins, stainless_point, stainless_pressure) = inspect(
            "stainless",
            &stainless,
            &stainless_path,
            &workspace_path(
                "examples/material-discovery/stainless-316-thermoelastic-reference.json",
            ),
            293.15,
            [
                8000.0,
                194_644_525_149.347_6,
                0.30,
                1.538_267_176_297_722_9e-5,
                485.417_304_218_698_3,
                15.127_317_130_127_033,
            ],
            true,
        );
        let t = 293.15_f64;
        let epsilon = 1.0e-5
            * (-295.54 - 0.39811 * t + 0.0092683 * t.powi(2) - 0.000020261 * t.powi(3)
                + 0.000000017127 * t.powi(4));
        let epsilon_prime = 1.0e-5
            * (-0.39811 + 2.0 * 0.0092683 * t - 3.0 * 0.000020261 * t.powi(2)
                + 4.0 * 0.000000017127 * t.powi(3));
        relative_close(epsilon, 0.0000030538065389264904375);
        relative_close(
            stainless_state
                .resolved()
                .property("linear_thermal_expansion_coefficient")
                .unwrap()
                .value_si(),
            epsilon_prime / (1.0 + epsilon),
        );
        relative_close(
            stainless_state
                .resolved()
                .property("young_modulus")
                .unwrap()
                .value_si(),
            1.0e9
                * (207.9488 + 0.07394241 * t - 0.0009627200 * t.powi(2)
                    + 0.000002845560 * t.powi(3)
                    - 0.0000000032408 * t.powi(4)),
        );
        assert_ne!(
            aluminum_state.resolved().identity(),
            stainless_state.resolved().identity()
        );
        assert_ne!(aluminum_pressure.to_bits(), stainless_pressure.to_bits());
        assert!(
            resolve_isotropic_thermoelastic_state_point(
                stainless.card(),
                &stainless_point,
                MaterialPropertySelection::PinnedByProperty(aluminum_pins.clone()),
            )
            .is_err()
        );
        assert!(
            resolve_isotropic_thermoelastic_state_point(
                aluminum.card(),
                &aluminum_point,
                MaterialPropertySelection::PinnedByProperty(aluminum_pins[..5].to_vec()),
            )
            .is_err()
        );
        let aluminum_density = aluminum.card().claims().claims_for("density")[0].1;
        let stainless_density = stainless.card().claims().claims_for("density")[0].1;
        for refused in [
            point(aluminum_density, &[("T", 300.01)]),
            point(stainless_density, &[("T", 293.16)]),
            point(stainless_density, &[("source-pressure-known", 1.0)]),
        ] {
            assert!(
                resolve_isotropic_thermoelastic_state_point(
                    if refused.axes().contains_key("source-pressure-known") {
                        stainless.card()
                    } else {
                        aluminum.card()
                    },
                    &refused,
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .is_err()
            );
        }
    }

    /// G1/G3: real source compilation and persistent transport reach the
    /// nonlinear heat solve. This qualifies the declared interpolant, not a
    /// physical stainless specimen or an unspecified pressure condition.
    #[test]
    fn g1_g3_sourced_stainless_store_to_conduction() {
        use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
        use fs_conduction::field::ScalarField;
        use fs_conduction::material::{ConductivityModel, ConductivityTable, ProvenanceClass};
        use fs_conduction::mesh::ConductionMesh;
        use fs_conduction::solve::{
            ConductionProblem, InitialGuess, SolveConfig, StopRule, element_heat_flux, solve,
        };
        use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point,
        };
        use fs_qty::QuantitySpec;
        use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

        let database = fixture_dir().join("stainless-conduction.sqlite");
        let packs: Vec<_> = [
            (
                "stainless-304-nist-cryogenic",
                "UNS S30400",
                "process unspecified",
            ),
            (
                "stainless-316-nist-cryogenic",
                "UNS S31600",
                "process unspecified",
            ),
            ("aluminum-6061-t6-cryogenic", "UNS A96061", "T6 temper"),
        ]
        .into_iter()
        .map(|(slug, chemistry, process)| {
            let (pack, _) = compile(slug);
            CatalogPack::MaterialCard(
                NormalizedMaterialCardPack::new(
                    MaterialStateId {
                        chemistry: chemistry.into(),
                        phase: "solid".into(),
                        process: format!("NIST cryogenic compilation; {process}"),
                        revision: 0,
                    },
                    pack,
                )
                .unwrap(),
            )
        })
        .collect();
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store.ingest_bundle(&packs).unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let mut tables = Vec::new();
        for original in &packs {
            let CatalogPack::MaterialCard(original) = original else {
                unreachable!()
            };
            let CatalogPack::MaterialCard(loaded) =
                store.load_catalog_pack(original.pack_id()).unwrap()
            else {
                panic!("wrong stored family")
            };
            assert_eq!(loaded, *original);
            let claims = loaded.card().claims();
            let (pin, claim) = claims.claims_for("thermal-conductivity")[0];
            let points = [250.0, 275.0].map(|t| point(claim, &[("temperature", t)]));
            let table = ConductivityTable::from_claims_at_query_points(
                claims,
                &claim.key,
                "temperature",
                &points,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
            let pinned = ConductivityTable::from_claims_pinned_at_query_points(
                claims,
                &claim.key,
                "temperature",
                &points,
                pin,
            )
            .unwrap();
            assert_eq!(pinned.knots(), table.knots());
            assert_eq!(table.knots().len(), table.receipts().len());
            assert_eq!(table.knots().len(), pinned.receipts().len());
            let requirement = ScalarPropertyRequirement::try_with_key(
                &claim.key,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap();
            for (((temperature, _), receipt), pinned_receipt) in table
                .knots()
                .iter()
                .zip(table.receipts())
                .zip(pinned.receipts())
            {
                let at = point(claim, &[("temperature", *temperature)]);
                let resolved = resolve_material_state_point(
                    loaded.card(),
                    &at,
                    std::slice::from_ref(&requirement),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
                let property = resolved.property("thermal-conductivity").unwrap();
                close(
                    table.eval(at.axes()["temperature"]).unwrap(),
                    property.value_si(),
                );
                assert_eq!(*receipt, property.answer().receipt);
                assert_eq!(&receipt.axis_quantities, at.axis_quantities());
                assert!(
                    receipt
                        .query_point
                        .contains(&("source-pressure-known".into(), 0.0))
                );
                assert!(!receipt.source_hashes.is_empty());
                claims.verify_receipt(receipt).unwrap();
                claims.verify_receipt(pinned_receipt).unwrap();
            }

            // The original T-only constructor cannot silently drop the source
            // axis spelling or turn unknown pressure into a declared pressure.
            assert!(
                ConductivityTable::from_claims(
                    claims,
                    "thermal-conductivity",
                    &[250.0, 275.0],
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
            );
            let absolute = claim.validity.axis_quantities()["temperature"];
            let missing_pressure = [250.0, 275.0].map(|t| {
                QueryPoint::new()
                    .with_quantity("temperature", absolute, t)
                    .unwrap()
            });
            let wrong_pressure = [250.0, 275.0]
                .map(|t| point(claim, &[("temperature", t), ("source-pressure-known", 1.0)]));
            let outside = [275.0, 301.0].map(|t| point(claim, &[("temperature", t)]));
            let mut changing_context = points.clone();
            changing_context[1] = wrong_pressure[1].clone();
            let mut wrong_kind = points.clone();
            for at in &mut wrong_kind {
                *at = at
                    .clone()
                    .with_quantity(
                        "temperature",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::TemperatureDifference,
                            ValueForm::Static,
                        )),
                        at.axes()["temperature"],
                    )
                    .unwrap();
            }
            for invalid in [
                &missing_pressure,
                &wrong_pressure,
                &outside,
                &changing_context,
                &wrong_kind,
            ] {
                assert!(
                    ConductivityTable::from_claims_at_query_points(
                        claims,
                        &claim.key,
                        "temperature",
                        invalid,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
                );
            }
            assert!(
                ConductivityTable::from_claims_at_query_points(
                    claims,
                    &claim.key,
                    "T",
                    &points,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
            );
            assert!(table.eval(249.0).is_err());
            assert!(table.eval(276.0).is_err());
            tables.push(table);
        }

        // Source coefficients, independently evaluated at these NIST knots.
        let k250 = 13.9812447520621;
        let k275 = 14.6459875036829;
        close(tables[0].eval(250.0).unwrap(), k250);
        close(tables[0].eval(275.0).unwrap(), k275);
        let slope = (k275 - k250) / 25.0;
        // T = 250 + 25*x, k linear between these two knots, hence
        // -div(k(T) grad T) = -k'(T)*25^2 W/m^3 on the unit cube.
        let (complex, positions) = fs_conduction::fixtures::unit_cube(3);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let reference: Vec<_> = mesh
            .positions()
            .iter()
            .map(|p| 250.0 + 25.0 * p[0])
            .collect();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "manufactured-temperature",
                |_| true,
                ThermalBc::Dirichlet {
                    temperature: ScalarField::Nodal(reference.clone()),
                },
            )
            .unwrap()
            .finish()
            .unwrap();
        let source = ScalarField::Uniform(-slope * 25.0 * 25.0);
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut solutions = Vec::new();
        for (index, table) in tables.into_iter().enumerate() {
            let expected_receipts = table.receipts().len();
            let material = ConductivityModel::isotropic(table);
            let solution = pool.scope(|arena| {
                let cx = Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 0xC071,
                        kernel_id: 51,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                );
                solve(
                    &cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: &material,
                        element_materials: None,
                        source: &source,
                    },
                    SolveConfig {
                        initial: InitialGuess::Uniform(262.5),
                        // The high-conductivity aluminum boundary loads are
                        // much larger than the net source used by the energy
                        // check; its relative closure needs a tighter solve.
                        stop: StopRule {
                            residual_rtol: 1.0e-13,
                            ..StopRule::default()
                        },
                        ..SolveConfig::default()
                    },
                )
                .unwrap()
            });
            assert_eq!(
                solution.report.free_dofs, 8,
                "must actually solve interior temperatures"
            );
            assert_eq!(
                solution.report.material_provenance,
                ProvenanceClass::MatdbReceipts
            );
            assert_eq!(solution.report.material_receipts, expected_receipts);
            assert!(
                solution.report.energy.relative_closure() < 1.0e-8,
                "material {index}: {:?}",
                solution.report.energy
            );
            if index <= 1 {
                for (actual, expected) in solution.temperature.iter().zip(&reference) {
                    assert!(
                        (actual - expected).abs() < 1.0e-7,
                        "temperature {actual}, reference {expected}"
                    );
                }
                let fluxes = element_heat_flux(&mesh, &material, &solution.temperature).unwrap();
                for (tet, flux) in mesh.complex().tets.iter().zip(fluxes) {
                    let mean = tet.iter().map(|&v| reference[v as usize]).sum::<f64>() / 4.0;
                    let expected = -25.0 * (k250 + slope * (mean - 250.0));
                    assert!((flux[0] - expected).abs() < 1.0e-5);
                    assert!(flux[1].abs().max(flux[2].abs()) < 1.0e-5);
                }
            }
            solutions.push(solution.temperature);
        }
        // NIST publishes the same conductivity fit for these two grades.
        // A material-name change alone must not invent a different response.
        assert_eq!(solutions[0], solutions[1]);
        let material_effect = solutions[0]
            .iter()
            .zip(&solutions[2])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(
            material_effect > 1.0e-5,
            "substituting sourced 6061-T6 for 304 must affect the same thermal problem: {material_effect}"
        );
    }

    /// G1/G3: conduction-only liquid water at a fixed pressure and phase.
    /// This steady model needs k(T), not Cp or viscosity; it does not model
    /// buoyancy, fluid motion, transient heating or phase change.
    #[test]
    fn g1_g3_sourced_water_conduction_at_two_states() {
        use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
        use fs_conduction::field::ScalarField;
        use fs_conduction::material::{ConductivityModel, ConductivityTable, ProvenanceClass};
        use fs_conduction::mesh::ConductionMesh;
        use fs_conduction::solve::{
            ConductionProblem, InitialGuess, SolveConfig, StopRule, element_heat_flux, solve,
        };
        use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point,
        };

        let (pack, _) = compile("water-liquid-iapws-sr6-08");
        let original = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "IAPWS ordinary water".into(),
                phase: "liquid".into(),
                process: "SR6-08(2011) correlation at 0.1 MPa".into(),
                revision: 0,
            },
            pack,
        )
        .unwrap();
        let store = MaterialStore::open(":memory:").unwrap();
        store
            .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
            .unwrap();
        store.seal_corpus().unwrap();
        let CatalogPack::MaterialCard(loaded) =
            store.load_catalog_pack(original.pack_id()).unwrap()
        else {
            panic!("wrong stored family")
        };
        assert_eq!(loaded, original);
        let claims = loaded.card().claims();
        let claim = claims.claims_for("thermal-conductivity")[0].1;
        let requirement = ScalarPropertyRequirement::try_with_key(
            &claim.key,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap();
        let (complex, positions) = fs_conduction::fixtures::unit_cube(3);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut mean_fluxes = Vec::new();
        // Independently evaluated SR6-08 equation (8), at adjacent source knots.
        for (lo, k_lo, k_hi) in [
            (293.15, 0.598004651798403, 0.606502307735785),
            (333.15, 0.651015542988153, 0.655588716984780),
        ] {
            let delta = 5.0;
            let points = [lo, lo + delta].map(|t| point(claim, &[("temperature", t)]));
            let table = ConductivityTable::from_claims_at_query_points(
                claims,
                &claim.key,
                "temperature",
                &points,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
            for ((at, expected), receipt) in points.iter().zip([k_lo, k_hi]).zip(table.receipts()) {
                let state = resolve_material_state_point(
                    loaded.card(),
                    at,
                    std::slice::from_ref(&requirement),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
                let selected = state.property("thermal-conductivity").unwrap();
                close(selected.value_si(), expected);
                assert_eq!(&selected.answer().receipt, receipt);
                assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
                assert!(receipt.query_point.contains(&("phase-liquid".into(), 1.0)));
                claims.verify_receipt(receipt).unwrap();
            }
            let material = ConductivityModel::isotropic(table);
            let reference: Vec<_> = mesh.positions().iter().map(|p| lo + delta * p[0]).collect();
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region(
                    "conduction-only-water",
                    |_| true,
                    ThermalBc::Dirichlet {
                        temperature: ScalarField::Nodal(reference.clone()),
                    },
                )
                .unwrap()
                .finish()
                .unwrap();
            // For T=lo+delta*x and linear k(T), -div(k grad T)=-delta*(k_hi-k_lo).
            let source = ScalarField::Uniform(-delta * (k_hi - k_lo));
            let solution = pool.scope(|arena| {
                let cx = Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 0xA91,
                        kernel_id: 52,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                );
                solve(
                    &cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: &material,
                        element_materials: None,
                        source: &source,
                    },
                    SolveConfig {
                        initial: InitialGuess::Uniform(lo + delta / 2.0),
                        stop: StopRule {
                            residual_rtol: 1e-13,
                            ..StopRule::default()
                        },
                        ..SolveConfig::default()
                    },
                )
                .unwrap()
            });
            assert_eq!(solution.report.free_dofs, 8);
            assert_eq!(
                solution.report.material_provenance,
                ProvenanceClass::MatdbReceipts
            );
            assert_eq!(solution.report.material_receipts, 2);
            assert!(
                solution.report.energy.relative_closure() < 1e-8,
                "{:?}",
                solution.report.energy
            );
            for (actual, expected) in solution.temperature.iter().zip(&reference) {
                assert!(
                    (actual - expected).abs() < 1e-7,
                    "T={actual}, expected={expected}"
                );
            }
            let fluxes = element_heat_flux(&mesh, &material, &solution.temperature).unwrap();
            for (tet, flux) in mesh.complex().tets.iter().zip(&fluxes) {
                let mean = tet.iter().map(|&v| reference[v as usize]).sum::<f64>() / 4.0;
                let expected = -delta * (k_lo + (k_hi - k_lo) * (mean - lo) / delta);
                assert!((flux[0] - expected).abs() < 1e-5);
                assert!(flux[1].abs().max(flux[2].abs()) < 1e-5);
            }
            let mean_flux = fluxes.iter().map(|q| q[0]).sum::<f64>() / fluxes.len() as f64;
            println!(
                "water conduction {lo}..{} K, 100000 Pa, liquid: mean x-flux={mean_flux} W/m2; relative energy closure={}",
                lo + delta,
                solution.report.energy.relative_closure()
            );
            mean_fluxes.push(mean_flux);
        }
        // The same 5 K/m gradient transfers more heat in the warmer source state.
        assert!(mean_fluxes[1].abs() > 1.08 * mean_fluxes[0].abs());
        for invalid in [
            ("temperature", 273.15),
            ("pressure", 101325.0),
            ("phase-liquid", 0.0),
        ] {
            // Keep the grid increasing so an unsupported temperature is refused
            // by source coverage rather than by duplicate-grid validation.
            let points = [0.0, 5.0].map(|offset| {
                if invalid.0 == "temperature" {
                    point(claim, &[("temperature", invalid.1 + offset)])
                } else {
                    point(claim, &[("temperature", 293.15 + offset)])
                        .with(invalid.0, invalid.1)
                        .unwrap()
                }
            });
            assert!(
                ConductivityTable::from_claims_at_query_points(
                    claims,
                    &claim.key,
                    "temperature",
                    &points,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
            );
            assert!(
                resolve_material_state_point(
                    loaded.card(),
                    &points[0],
                    std::slice::from_ref(&requirement),
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
            );
        }
    }

    /// G1/G3: an isobaric liquid parcel heats using the sourced enthalpy curve,
    /// without fabricated melting endpoints or a frozen initial heat capacity.
    #[test]
    fn g1_g3_sourced_liquid_water_enthalpy_heating() {
        use fs_conduction::lumped::{
            BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig, solve_lumped_enthalpy,
        };
        use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::phase::{
            EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase,
        };

        let (pack, path) = compile("water-liquid-iapws-sr6-08");
        let original = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "IAPWS ordinary water".into(),
                phase: "liquid".into(),
                process: "SR6-08(2011), isobaric 0.1 MPa".into(),
                revision: 0,
            },
            pack,
        )
        .unwrap();
        let db = fixture_dir().join("liquid-water-enthalpy.sqlite");
        {
            let store = MaterialStore::open(db.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(db.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(loaded) =
            store.load_catalog_pack(original.pack_id()).unwrap()
        else {
            panic!("wrong stored family")
        };
        assert_eq!(loaded, original);
        let card = loaded.card();
        let h = card.claims().claims_for("specific-enthalpy")[0].1;
        let rho = card.claims().claims_for("density")[0].1;
        let k = card.claims().claims_for("thermal-conductivity")[0].1;
        let PropertyValue::Curve { knots: h_knots, .. } = &h.value else {
            panic!("source enthalpy curve required")
        };
        assert_eq!(h_knots.len(), 15);
        let mut phase_knots = Vec::new();
        let mut minimum_k = f64::INFINITY;
        let mut receipts = Vec::new();
        for &(temperature, _) in h_knots {
            let values: Vec<_> = [h, rho, k]
                .iter()
                .map(|claim| {
                    let answer = card
                        .claims()
                        .query_typed(
                            &claim.key,
                            &point(claim, &[("temperature", temperature)]),
                            SelectionPolicy::SingleClaimOnly,
                        )
                        .unwrap();
                    card.claims().verify_receipt(&answer.receipt).unwrap();
                    receipts.push(answer.receipt.clone());
                    answer.evidence.value.value
                })
                .collect();
            phase_knots.push(EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: values[0],
                temperature_k: temperature,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: values[1],
            });
            minimum_k = minimum_k.min(values[2]);
        }
        assert_eq!(receipts.len(), 45);
        assert!(
            EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), phase_knots.clone())
                .is_err()
        );
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Liquid,
            phase_knots,
        )
        .unwrap();
        // No boundary heat exchange: conductivity supplies a conservative Biot
        // input only. No emissivity observation is needed for disabled radiation.
        // Mass is fixed; volume may change at prescribed pressure. This does not
        // solve flow, geometry evolution, or spatial temperature gradients.
        let mass = 0.1;
        let body = LumpedEnthalpyBody::try_new(
            "isobaric-liquid-water",
            mass,
            0.01,
            0.0,
            0.0,
            0.01,
            minimum_k,
            &curve,
        )
        .unwrap();
        // Independent SR6-08 table ordinates at 20 C and 60 C.
        let initial_h = 84005.842699135;
        let final_h = 251246.640203514;
        let config = LumpedEnthalpyMarchConfig {
            initial_specific_enthalpy_j_kg: initial_h,
            ambient_temperature_k: 293.15,
            radiation_temperature_k: 293.15,
            internal_power_w: mass * (final_h - initial_h) / 100.0,
            duration_s: 100.0,
            maximum_step_s: 1.0,
            maximum_steps: 100,
            enthalpy_tolerance_j_kg: 1e-7,
        };
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate, arena, StreamKey { seed: 0xA91, kernel_id: 53, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
            let march = solve_lumped_enthalpy(&cx, &body, BiotGate::corpus_default(), config).unwrap();
            assert_eq!(march, solve_lumped_enthalpy(&cx, &body, BiotGate::corpus_default(), config).unwrap());
            assert_eq!(march.samples().len(), 101);
            for sample in march.samples() {
                let expected_h = initial_h + config.internal_power_w * sample.time_s / mass;
                assert!((sample.phase_state.specific_enthalpy_j_kg() - expected_h).abs() < 1e-5);
                let pair = h_knots.windows(2).find(|pair| pair[0].1 <= expected_h && expected_h <= pair[1].1).unwrap();
                let expected_t = pair[0].0 + (expected_h - pair[0].1) * (pair[1].0 - pair[0].0) / (pair[1].1 - pair[0].1);
                assert!((sample.phase_state.temperature_k() - expected_t).abs() < 1e-7);
                assert_eq!(sample.phase_state.phase(), SolidLiquidPhase::Liquid);
                assert_eq!(sample.phase_state.material_card_identity(), card.content_hash());
                assert_eq!(sample.convection_into_body_w, 0.0);
                assert_eq!(sample.radiation_into_body_w, 0.0);
            }
            let last = march.samples().last().unwrap().phase_state;
            assert!((last.temperature_k() - 333.15).abs() < 1e-7);
            assert!(march.cumulative_absolute_energy_residual_j() < 1e-5);
            let outside = LumpedEnthalpyMarchConfig { internal_power_w: 1000.0, ..config };
            assert!(solve_lumped_enthalpy(&cx, &body, BiotGate::corpus_default(), outside).is_err());
            assert_eq!(march, solve_lumped_enthalpy(&cx, &body, BiotGate::corpus_default(), config).unwrap());
            eprintln!("water isobaric=100000Pa mass={mass}kg samples={} initial=293.15K final={}K heat={}J energy_residual={}J receipts=45 source_span=283.15..353.15K no_phase_transition_or_spatial_flow_claim=1", march.samples().len(), last.temperature_k(), mass*(last.specific_enthalpy_j_kg()-initial_h), march.cumulative_absolute_energy_residual_j());
        });
        for (axis, wrong) in [
            ("pressure", 100001.0),
            ("phase-liquid", 0.0),
            ("temperature", 353.16),
            ("enthalpy-reference-iapws-sr6-08-eq1", 0.0),
        ] {
            assert!(
                card.claims()
                    .query_typed(
                        &h.key,
                        &point(h, &[(axis, wrong)]),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err(),
                "{axis}"
            );
        }
        let discovery = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            workspace_path("examples/material-discovery/liquid-water.json")
                .to_str()
                .unwrap()
                .into(),
            path.to_str().unwrap().into(),
        ]);
        assert_eq!(
            discovery.exit_code,
            fs_cli::exit::SUCCESS,
            "{}",
            discovery.stderr
        );
        assert!(
            discovery.stdout.contains("\"status\":\"complete\""),
            "{}",
            discovery.stdout
        );
    }

    fn check_claims(pack: &NormalizedPack) {
        assert!(pack.claims().claims_ordered().next().is_some());
        for (id, claim) in pack.claims().claims_ordered() {
            assert!(!claim.observations.is_empty());
            let samples = match &claim.value {
                PropertyValue::Curve {
                    abscissa, knots, ..
                } => knots
                    .iter()
                    .map(|&(x, y)| (point(claim, &[(abscissa, x)]), y))
                    .collect::<Vec<_>>(),
                PropertyValue::Scalar { value, .. } => vec![(point(claim, &[]), *value)],
            };
            for (at, expected) in samples {
                let answer = pack
                    .claims()
                    .query_pinned_typed(&claim.key, &at, id)
                    .unwrap();
                close(answer.evidence.value.value, expected);
                pack.claims().verify_receipt(&answer.receipt).unwrap();
            }
            assert!(
                pack.claims()
                    .query_pinned_typed(&claim.key, &QueryPoint::new(), id)
                    .is_err(),
                "source context must not disappear for {}",
                claim.key.name()
            );
            if let Some((_, hi)) = claim.validity.bound("temperature") {
                let outside = point(claim, &[("temperature", hi + 1.0)]);
                assert!(
                    pack.claims()
                        .query_pinned_typed(&claim.key, &outside, id)
                        .is_err()
                );
            }
        }
    }

    #[test]
    fn g0_g3_cryogenic_aluminum_copper_curves_and_discovery() {
        use fs_matdb_store::{CatalogPack, MaterialStore};

        // Published NIST equations are independent of the retained TSV knots.
        // Coefficients are in ascending polynomial order; E is in GPa.
        let cases: [(&str, &str, &[f64], f64); 5] = [
            (
                "aluminum-6061-t6-cryogenic",
                "thermal-conductivity",
                &[
                    0.07918, 1.0957, -0.07277, 0.08084, 0.02803, -0.09464, 0.04179, -0.00571, 0.0,
                ],
                0.00086,
            ),
            (
                "aluminum-6061-t6-cryogenic",
                "specific-heat-capacity",
                &[
                    46.6467, -314.292, 866.662, -1298.3, 1162.27, -637.795, 210.351, -38.3094,
                    2.96344,
                ],
                0.00185,
            ),
            (
                "aluminum-6061-t6-cryogenic",
                "young-modulus",
                &[
                    77.71221,
                    0.01030646,
                    -0.0002924100,
                    0.00000089936,
                    -0.0000000010709,
                ],
                0.000036,
            ),
            (
                "ofhc-copper-rrr100",
                "thermal-conductivity",
                &[
                    2.2154, -0.47461, -0.88068, 0.13871, 0.29505, -0.02043, -0.04831, 0.001281,
                    0.003207,
                ],
                0.00556,
            ),
            (
                "ofhc-copper-rrr100",
                "specific-heat-capacity",
                &[
                    -1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797, 0.0,
                ],
                0.00277,
            ),
        ];
        let compiled = [
            compile("aluminum-6061-t6-cryogenic"),
            compile("ofhc-copper-rrr100"),
        ];
        let store = MaterialStore::open(":memory:").unwrap();
        store
            .ingest_bundle(
                &compiled
                    .iter()
                    .map(|(pack, _)| CatalogPack::Properties(pack.clone()))
                    .collect::<Vec<_>>(),
            )
            .unwrap();
        store.seal_corpus().unwrap();
        for (slug, name, coefficients, sampled_error_limit) in cases {
            let original = &compiled[usize::from(slug.starts_with("ofhc"))].0;
            let CatalogPack::Properties(pack) =
                store.load_catalog_pack(original.pack_id()).unwrap()
            else {
                panic!("wrong property family")
            };
            assert_eq!(&pack, original);
            check_claims(&pack);
            let claim = pack.claims().claims_for(name)[0].1;
            let PropertyValue::Curve { knots, .. } = &claim.value else {
                panic!("missing {name} curve")
            };
            assert_eq!(knots.len(), 24);
            assert_eq!(claim.validity.bound("temperature"), Some((77.0, 293.0)));
            let equation = |t: f64| {
                if slug.starts_with("ofhc") && name == "thermal-conductivity" {
                    let a = coefficients;
                    let s = t.sqrt();
                    10.0f64.powf(
                        (a[0] + a[2] * s + a[4] * t + a[6] * t * s + a[8] * t * t)
                            / (1.0 + a[1] * s + a[3] * t + a[5] * t * s + a[7] * t * t),
                    )
                } else {
                    let x = if name == "young-modulus" {
                        t
                    } else {
                        t.log10()
                    };
                    let p = coefficients.iter().rev().fold(0.0, |p, a| p * x + a);
                    if name == "young-modulus" {
                        p * 1e9
                    } else {
                        10.0f64.powf(p)
                    }
                }
            };
            for &(t, value) in knots {
                let expected = equation(t);
                assert!(
                    (value / expected - 1.0).abs() < 1.0e-8,
                    "{slug} {name} at {t}: {value} != {expected}"
                );
            }
            let mut maximum = 0.0f64;
            for pair in knots.windows(2) {
                for i in 1..100 {
                    let t = pair[0].0 + (pair[1].0 - pair[0].0) * f64::from(i) / 100.0;
                    let value = sample(&pack, name, &[("temperature", t)]);
                    maximum = maximum.max((value / equation(t) - 1.0).abs());
                }
            }
            assert!(
                maximum <= sampled_error_limit,
                "{slug} {name}: sampled interpolation discrepancy {maximum} exceeds {sampled_error_limit}"
            );
            println!(
                "{slug} {name}: 24 source knots, sampled maximum relative interpolation discrepancy {maximum:.9}"
            );
            for overrides in [
                vec![("temperature", 76.0)],
                vec![("temperature", 294.0)],
                vec![("temperature", 373.15)],
                vec![("temperature", 150.0), ("source-pressure-known", 1.0)],
            ] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &claim.key,
                            &point(claim, &overrides),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
        }
        let request = fixture_dir().join("cryogenic-aluminum-copper.json");
        fs::write(
            &request,
            include_str!("../../examples/material-discovery/cryogenic-aluminum-copper.json"),
        )
        .unwrap();
        let run = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            request.to_str().unwrap().into(),
            compiled[0].1.to_str().unwrap().into(),
            compiled[1].1.to_str().unwrap().into(),
        ]);
        assert_eq!(run.exit_code, fs_cli::exit::SUCCESS, "{}", run.stderr);
        assert_eq!(
            run.stdout.matches("\"status\":\"complete\"").count(),
            2,
            "{}",
            run.stdout
        );
        assert_eq!(
            run.stdout.matches("\"status\":\"supported\"").count(),
            4,
            "{}",
            run.stdout
        );
        // Coverage of unbound property packs does not invent a coherent
        // specimen or bind the RRR-unspecified heat capacity to RRR=100.
        assert_eq!(
            run.stdout
                .matches("unbound to a named material state")
                .count(),
            2,
            "{}",
            run.stdout
        );
    }

    /// Source equations and Table 8 are independent of the retained TSV values.
    /// Linear engineering curves have the measured numerical discrepancies
    /// below; this test does not establish a fluid solve or phase transition.
    #[test]
    fn g0_g3_iapws_liquid_water_curves_and_discovery() {
        use fs_matdb_store::{CatalogPack, MaterialStore};

        let (original, path) = compile("water-liquid-iapws-sr6-08");
        let (_, rebuilt) = compile("water-liquid-iapws-sr6-08");
        assert_eq!(fs::read(&path).unwrap(), fs::read(rebuilt).unwrap());
        let store = MaterialStore::open(":memory:").unwrap();
        store
            .ingest_bundle(&[CatalogPack::Properties(original.clone())])
            .unwrap();
        store.seal_corpus().unwrap();
        let CatalogPack::Properties(pack) = store.load_catalog_pack(original.pack_id()).unwrap()
        else {
            panic!("wrong property family")
        };
        assert_eq!(pack, original);
        assert_eq!(pack.claims().claims_ordered().count(), 5);
        check_claims(&pack);

        let equation = |t: f64| {
            let r = 461.51805;
            let tau = t / 10.0;
            let alpha = 10.0 / (593.0 - t);
            let beta = 10.0 / (t - 232.0);
            let a = [(-1.661470539e5, 4), (2.708781640e6, 5), (-1.557191544e8, 7)];
            let b = [
                (-0.8237426256, 2),
                (1.908956353, 3),
                (-2.017597384, 4),
                (0.8546361348, 5),
            ];
            let sum = |terms: &[(f64, i32)], x: f64, derivative: i32| {
                terms
                    .iter()
                    .map(|&(c, n)| {
                        let factor = match derivative {
                            0 => 1.0,
                            1 => f64::from(n),
                            _ => f64::from(n * (n + 1)),
                        };
                        c * factor * x.powi(n + derivative)
                    })
                    .sum::<f64>()
            };
            let g = r
                * 10.0
                * (-245.2093414 + 38.69269598 * tau - 8.983025854 * tau * tau.ln()
                    + sum(&a, alpha, 0)
                    + sum(&b, beta, 0));
            let s = -r
                * (38.69269598 - 8.983025854 * (1.0 + tau.ln()) + sum(&a, alpha, 1)
                    - sum(&b, beta, 1));
            let cp = -r * (-8.983025854 + tau * (sum(&a, alpha, 2) + sum(&b, beta, 2)));
            let volume = r * 10.0 / 100000.0
                * (0.0193763157
                    + sum(
                        &[
                            (6744.58446, 4),
                            (-222521.604, 5),
                            (100231247.0, 7),
                            (-1635521180.0, 8),
                            (8322996580.0, 9),
                        ],
                        alpha,
                        0,
                    )
                    + sum(
                        &[
                            (0.00578545292, 1),
                            (-0.0153195665, 2),
                            (0.0311337859, 3),
                            (-0.0423546241, 4),
                            (0.0338713507, 5),
                            (-0.0119946761, 6),
                        ],
                        beta,
                        0,
                    ));
            let transport = |terms: &[(f64, f64)]| {
                terms
                    .iter()
                    .map(|&(c, n)| c * (t / 300.0).powf(n))
                    .sum::<f64>()
            };
            [
                1.0 / volume,
                cp,
                g + t * s,
                transport(&[
                    (1.6630, -1.15),
                    (-1.7781, -3.4),
                    (1.1567, -6.0),
                    (-0.432115, -7.6),
                ]),
                1e-6 * transport(&[
                    (280.68, -1.9),
                    (511.45, -7.7),
                    (61.131, -19.6),
                    (0.45903, -40.0),
                ]),
            ]
        };
        let cases = [
            ("density", "kg/m3", 997.047013, 0.000039252),
            ("specific-heat-capacity", "J/kg/K", 4181.44618, 0.000078733),
            (
                "specific-enthalpy",
                "J/kg",
                -4561.7537 + 298.15 * 367.20145,
                0.000079557,
            ),
            ("thermal-conductivity", "W/m/K", 0.606502308, 0.000188763),
            ("dynamic-viscosity", "Pa*s", 889.996774e-6, 0.003833467),
        ];
        for (index, (name, unit, table8, interpolation_limit)) in cases.into_iter().enumerate() {
            let claim = pack.claims().claims_for(name)[0].1;
            assert_eq!(
                claim.key.dims(),
                fs_qty::parse::parse_qty(&format!("1 {unit}")).unwrap().dims
            );
            assert_eq!(claim.validity.bound("temperature"), Some((283.15, 353.15)));
            assert_eq!(claim.validity.bound("pressure"), Some((100000.0, 100000.0)));
            let PropertyValue::Curve { knots, .. } = &claim.value else {
                panic!("missing {name} curve")
            };
            assert_eq!(knots.len(), 15);
            // Table 8 rounds printed values; h uses its separately rounded g,s.
            assert!((sample(&pack, name, &[("temperature", 298.15)]) / table8 - 1.0).abs() < 2e-8);
            for &(t, value) in knots {
                assert!(value > 0.0);
                assert!(
                    (value / equation(t)[index] - 1.0).abs() < 1e-11,
                    "{name} at {t}"
                );
            }
            let mut maximum = 0.0_f64;
            for pair in knots.windows(2) {
                for i in 1..100 {
                    let t = pair[0].0 + (pair[1].0 - pair[0].0) * f64::from(i) / 100.0;
                    maximum = maximum.max(
                        (sample(&pack, name, &[("temperature", t)]) / equation(t)[index] - 1.0)
                            .abs(),
                    );
                }
            }
            assert!(maximum <= interpolation_limit, "{name}: {maximum}");
            for override_value in [
                ("temperature", 283.14),
                ("temperature", 353.16),
                ("pressure", 101325.0),
                ("phase-liquid", 0.0),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &claim.key,
                            &point(claim, &[override_value]),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            println!(
                "water {name}: 15 source knots, sampled relative interpolation discrepancy {maximum:.9}"
            );
        }
        let enthalpy = pack.claims().claims_for("specific-enthalpy")[0].1;
        assert!(
            pack.claims()
                .query_typed(
                    &enthalpy.key,
                    &point(enthalpy, &[("enthalpy-reference-iapws-sr6-08-eq1", 0.0)]),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        let PropertyValue::Curve { knots, .. } = &enthalpy.value else {
            unreachable!()
        };
        for pair in knots.windows(2) {
            let dt = pair[1].0 - pair[0].0;
            let dh = pair[1].1 - pair[0].1;
            let cp0 = sample(
                &pack,
                "specific-heat-capacity",
                &[("temperature", pair[0].0)],
            );
            let cp1 = sample(
                &pack,
                "specific-heat-capacity",
                &[("temperature", pair[1].0)],
            );
            assert!((dh - dt * 0.5 * (cp0 + cp1)).abs() <= 1.100);
            assert!((dh / dt - cp0).abs().max((dh / dt - cp1).abs()) <= 3.519);
        }
        let request = workspace_path("examples/material-discovery/liquid-water.json");
        let run = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            request.to_str().unwrap().into(),
            path.to_str().unwrap().into(),
        ]);
        assert_eq!(run.exit_code, fs_cli::exit::SUCCESS, "{}", run.stderr);
        assert_eq!(
            run.stdout.matches("\"status\":\"complete\"").count(),
            1,
            "{}",
            run.stdout
        );
        assert_eq!(
            run.stdout.matches("\"status\":\"supported\"").count(),
            5,
            "{}",
            run.stdout
        );
        assert!(
            run.stdout.contains("unbound to a named material state"),
            "{}",
            run.stdout
        );
    }

    #[test]
    fn g0_g3_liquid_lead_transport_source() {
        let (pack, _) = compile("lead-liquid-transport-nasa-cr144016");
        check_claims(&pack);
        for (name, temperature, expected) in [
            ("density", 773.15, 10390.0),
            ("thermal-conductivity", 873.15, 3.8055 * 4.184),
            ("dynamic-viscosity", 824.15, 1.700 * 0.001),
            ("surface-tension", 673.15, 43800.0 * 0.00001),
        ] {
            close(
                sample(&pack, name, &[("temperature", temperature)]),
                expected,
            );
            let claim = pack.claims().claims_for(name)[0].1;
            assert_eq!(
                claim.validity.bound("source-pressure-known"),
                Some((0.0, 0.0))
            );
            assert!(
                pack.claims()
                    .query_typed(
                        &claim.key,
                        &point(claim, &[("source_phase_liquid", 0.0)]),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
        let density = pack.claims().claims_for("density")[0].1;
        assert!(
            pack.claims()
                .query_typed(
                    &density.key,
                    &point(density, &[("temperature", 650.0)]),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }

    #[test]
    fn g0_g3_lead_phase_tables_and_actual_discovery() {
        let (solid, _) = compile("lead-solid-nasa-tp3287");
        let (liquid, liquid_path) = compile("lead-liquid-nasa-tp3287");
        let (fusion, _) = compile("lead-fusion-nasa-tp3287");
        for pack in [&solid, &liquid, &fusion] {
            check_claims(pack);
        }
        close(
            sample(&solid, "specific-heat-capacity", &[("temperature", 300.0)]),
            26.673 / 0.2072,
        );
        close(
            sample(&liquid, "specific-heat-capacity", &[("temperature", 700.0)]),
            30.313 / 0.2072,
        );
        let at_melt = [("temperature", 600.65)];
        let jump = sample(&liquid, "specific-enthalpy-reference-29815k", &at_melt)
            - sample(&solid, "specific-enthalpy-reference-29815k", &at_melt);
        close(jump, (13.361 - 8.549) * 1000.0 / 0.2072);
        close(sample(&fusion, "latent-heat-fusion", &[]), jump);
        close(sample(&fusion, "melting-point", &[]), 600.65);
        // Independent columns of the printed source: trapezoidal cp integrates
        // to tabulated h within 0.003 kJ/mol on this 50 K solid interval.
        let cp_integral = 25.0
            * (sample(&solid, "specific-heat-capacity", &[("temperature", 400.0)])
                + sample(&solid, "specific-heat-capacity", &[("temperature", 450.0)]));
        let dh = sample(
            &solid,
            "specific-enthalpy-reference-29815k",
            &[("temperature", 450.0)],
        ) - sample(
            &solid,
            "specific-enthalpy-reference-29815k",
            &[("temperature", 400.0)],
        );
        assert!((cp_integral - dh).abs() < 3.0 / 0.2072);

        let request = fixture_dir().join("liquid-heat-capacity.json");
        fs::write(
            &request,
            include_str!("../../examples/material-discovery/lead-liquid-heat-capacity.json"),
        )
        .unwrap();
        let run = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            request.to_str().unwrap().into(),
            liquid_path.to_str().unwrap().into(),
        ]);
        assert_eq!(run.exit_code, fs_cli::exit::SUCCESS, "{}", run.stderr);
        assert!(
            run.stdout.contains("\"status\":\"complete\""),
            "{}",
            run.stdout
        );
        assert!(
            run.stdout.contains("\"unknown_properties\":[]"),
            "{}",
            run.stdout
        );
        assert!(run.stdout.contains("\"claim\":"));
        let cp = liquid.claims().claims_for("specific-heat-capacity")[0].1;
        for overrides in [
            vec![("temperature", 650.0), ("source_phase_liquid", 0.0)],
            vec![("temperature", 650.0), ("pressure", 101325.0)],
        ] {
            assert!(
                liquid
                    .claims()
                    .query_typed(
                        &cp.key,
                        &point(cp, &overrides),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
        println!("{}", run.stdout);
    }

    #[test]
    fn g0_g3_common_element_phase_tables() {
        // Cp anchors and fusion jumps are independent columns of NASA's
        // printed selected-function tables III.2, III.15, III.42 and III.51.
        for (element, mass, cp_molar, tm, fusion_molar) in [
            ("aluminum", 0.026981539, 24.200, 933.61, 10700.0),
            ("copper", 0.063546, 24.440, 1358.0, 13138.0),
            ("tin", 0.118710, 27.112, 505.12, 7195.0),
            ("zinc", 0.06539, 25.390, 692.73, 7300.0),
        ] {
            let (solid, _) = compile(&format!("{element}-solid-nasa-tp3287"));
            let (liquid, _) = compile(&format!("{element}-liquid-nasa-tp3287"));
            let (fusion, _) = compile(&format!("{element}-fusion-nasa-tp3287"));
            for pack in [&solid, &liquid, &fusion] {
                check_claims(pack);
            }
            close(
                sample(&solid, "specific-heat-capacity", &[("temperature", 298.15)]),
                cp_molar / mass,
            );
            close(
                sample(
                    &solid,
                    "specific-enthalpy-reference-29815k",
                    &[("temperature", 298.15)],
                ),
                0.0,
            );
            close(sample(&fusion, "melting-point", &[]), tm);
            close(
                sample(&fusion, "latent-heat-fusion", &[]),
                fusion_molar / mass,
            );
            let jump = sample(
                &liquid,
                "specific-enthalpy-reference-29815k",
                &[("temperature", tm)],
            ) - sample(
                &solid,
                "specific-enthalpy-reference-29815k",
                &[("temperature", tm)],
            );
            close(jump, fusion_molar / mass);
            if element == "tin" {
                // Table III.42: H(200)-H(298.15) = 3.740-6.323 kJ/mol.
                close(
                    sample(
                        &solid,
                        "specific-enthalpy-reference-29815k",
                        &[("temperature", 200.0)],
                    ),
                    -2583.0 / mass,
                );
            }
            println!(
                "{element}: three phase packs, source Cp, reference enthalpy and fusion jump pass"
            );
        }
    }

    #[test]
    fn g0_g3_wood_glass_and_pvc_source_conditions() {
        for (slug, name, temperature, expected) in [
            (
                "plywood-southern-pine-fpl-gtr282",
                "bending-moe",
                None,
                7.70e9,
            ),
            (
                "plywood-douglas-fir-fpl-gtr282",
                "bending-mor",
                None,
                41.37e6,
            ),
            (
                "osb-aspen-pu-1992-mill7-fpl-gtr282",
                "bending-moe-perpendicular",
                None,
                2.03e9,
            ),
            (
                "osb-southern-pine-biblis-1989-mill1-fpl-gtr282",
                "bending-moe-parallel",
                None,
                4.41e9,
            ),
            (
                "waterborne-preservative-treated-lumber-fpl-gtr282",
                "post-treatment-redrying-standard-temperature-limit",
                None,
                347.15,
            ),
            (
                "glass-borosilicate-duran-ntrs-19860021558",
                "specific-heat-capacity",
                Some(293.15),
                0.18 * 4184.0,
            ),
            (
                "glass-soda-lime-srm-1826b-nist",
                "density",
                Some(293.15),
                2548.668,
            ),
            (
                "pvc-cryogenic-nist",
                "specific-heat-capacity",
                Some(300.0),
                1167.15580013972,
            ),
        ] {
            let (pack, _) = compile(slug);
            check_claims(&pack);
            let at = temperature
                .map(|t| ("temperature", t))
                .into_iter()
                .collect::<Vec<_>>();
            close(sample(&pack, name, &at), expected);
            if slug == "glass-borosilicate-duran-ntrs-19860021558" {
                for strength in [
                    "tensile-strength-fire-bright",
                    "tensile-strength-cut-surface",
                ] {
                    assert!(
                        pack.claims().claims_for(strength).is_empty(),
                        "conflicting source unit columns cannot yield a strength claim"
                    );
                }
                let observation = pack
                    .claims()
                    .observation(pack.claims().observation_ids().next().unwrap())
                    .unwrap();
                for retained in [
                    "source unit conflict",
                    "7.8 MPa",
                    "11,400 psia",
                    "3.9 MPa",
                    "5,700 psia",
                ] {
                    assert!(observation.caveats.contains(retained));
                }
                let shock = pack
                    .claims()
                    .claims_for("thermal-shock-resistance-temperature-difference")[0]
                    .1;
                assert_eq!(
                    shock.key.quantity().semantic_type().unwrap().kind(),
                    fs_qty::semantic::QuantityKind::TemperatureDifference
                );
                close(sample(&pack, shock.key.name(), &[]), 250.0);
                let transition = pack.claims().claims_for("glass-transition-temperature")[0].1;
                assert_eq!(
                    transition.key.quantity().semantic_type().unwrap().kind(),
                    fs_qty::semantic::QuantityKind::AbsoluteTemperature
                );
                close(sample(&pack, transition.key.name(), &[]), 803.15);
            }
            if slug == "waterborne-preservative-treated-lumber-fpl-gtr282" {
                let limit = pack.claims().claims_for(name)[0].1;
                assert_eq!(
                    limit.key.quantity().semantic_type().unwrap().kind(),
                    fs_qty::semantic::QuantityKind::AbsoluteTemperature
                );
                assert!(
                    pack.claims().claims_for("young-modulus").is_empty(),
                    "process thresholds cannot manufacture a treated-wood modulus"
                );
            }
            if slug == "glass-soda-lime-srm-1826b-nist" {
                let density = pack.claims().claims_for("density")[0].1;
                assert_eq!(
                    density.uncertainty,
                    UncertaintyModel::HalfWidth {
                        half_width: 0.032,
                        confidence: 0.95
                    }
                );
            }
            if slug == "pvc-cryogenic-nist" {
                for (_, claim) in pack.claims().claims_ordered() {
                    assert_eq!(
                        claim.interpolation,
                        InterpolationPolicy::ConstantWithinValidity
                    );
                    assert_eq!(claim.validity.bound("temperature"), Some((300.0, 300.0)));
                    assert!(
                        pack.claims()
                            .query_typed(
                                &claim.key,
                                &point(claim, &[("temperature", 299.0)]),
                                SelectionPolicy::SingleClaimOnly,
                            )
                            .is_err()
                    );
                }
                let conductivity = pack.claims().claims_for("thermal-conductivity");
                assert_eq!(conductivity.len(), 2);
                // Different foam densities and fill gases must select their
                // own claims without a pin or an arbitrary tie breaker.
                for (_, claim) in conductivity {
                    let answer = pack
                        .claims()
                        .query_typed(
                            &claim.key,
                            &point(claim, &[]),
                            SelectionPolicy::SingleClaimOnly,
                        )
                        .unwrap();
                    let PropertyValue::Scalar { value, .. } = &claim.value else {
                        panic!("NIST fit evaluations are exact points");
                    };
                    close(answer.evidence.value.value, *value);
                }
            }
            println!("source pack {slug}: source value, typed context and query selection pass");
        }
    }

    #[test]
    fn g0_g3_stainless_thermal_envelope() {
        let mut paths = Vec::new();
        for (grade, cp_coefficients, cp_100, cp_300) in [
            (
                "304",
                [
                    22.0061, -127.5528, 303.647, -381.0098, 274.0328, -112.9212, 24.7593,
                    -2.239153, 0.0,
                ],
                275.496445572426,
                469.448840671653,
            ),
            (
                "316",
                [
                    -1879.464, 3643.198, 76.70125, -6176.028, 7437.6247, -4305.7217, 1382.4627,
                    -237.22704, 17.05262,
                ],
                273.023481294007,
                490.213382301636,
            ),
        ] {
            let (pack, path) = compile(&format!("stainless-{grade}-nist-cryogenic"));
            paths.push(path);
            check_claims(&pack);
            close(
                sample(&pack, "specific-heat-capacity", &[("temperature", 100.0)]),
                cp_100,
            );
            close(
                sample(&pack, "specific-heat-capacity", &[("temperature", 300.0)]),
                cp_300,
            );
            for (name, coefficients) in [
                ("specific-heat-capacity", cp_coefficients),
                (
                    "thermal-conductivity",
                    [
                        -1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199,
                    ],
                ),
            ] {
                let claims = pack.claims().claims_for(name);
                assert_eq!(
                    claims.len(),
                    1,
                    "old point claims must not create ambiguous lookup"
                );
                let claim = claims[0].1;
                assert_eq!(claim.interpolation, InterpolationPolicy::LinearInside);
                assert_eq!(claim.validity.bound("temperature"), Some((77.0, 300.0)));
                assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                // Numerical approximation check against the published fit at
                // 223 non-knot temperatures. This is neither a global error
                // bound nor physical validation of the original measurements.
                for integer in 77..300 {
                    let temperature = f64::from(integer) + 0.5;
                    let x = temperature.log10();
                    let polynomial = coefficients.iter().rev().fold(0.0, |sum, a| sum * x + a);
                    let source_fit = 10.0_f64.powf(polynomial);
                    let interpolated = sample(&pack, name, &[("temperature", temperature)]);
                    assert!(
                        (interpolated / source_fit - 1.0).abs() < 0.015,
                        "{grade} {name} at {temperature} K"
                    );
                }
                for overrides in [
                    vec![("temperature", 76.0)],
                    vec![("temperature", 301.0)],
                    vec![("temperature", 150.0), ("source-pressure-known", 1.0)],
                ] {
                    assert!(
                        pack.claims()
                            .query_typed(
                                &claim.key,
                                &point(claim, &overrides),
                                SelectionPolicy::SingleClaimOnly
                            )
                            .is_err()
                    );
                }
            }
            // The longer thermal range cannot extend the source's shorter
            // modulus fit or manufacture density for a heating trajectory.
            let modulus = pack.claims().claims_for("young-modulus")[0].1;
            assert!(
                pack.claims()
                    .query_typed(
                        &modulus.key,
                        &point(modulus, &[("temperature", 295.0)]),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
            assert!(pack.claims().claims_for("density").is_empty());
        }
        let request = fixture_dir().join("stainless-thermal.json");
        fs::write(
            &request,
            include_str!("../../examples/material-discovery/stainless-thermal.json"),
        )
        .unwrap();
        let run = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            request.to_str().unwrap().into(),
            paths[0].to_str().unwrap().into(),
            paths[1].to_str().unwrap().into(),
        ]);
        assert_eq!(run.exit_code, fs_cli::exit::SUCCESS, "{}", run.stderr);
        assert_eq!(
            run.stdout.matches("\"status\":\"complete\"").count(),
            2,
            "{}",
            run.stdout
        );
        assert_eq!(
            run.stdout.matches("\"status\":\"supported\"").count(),
            4,
            "{}",
            run.stdout
        );
        for grade in ["304", "316"] {
            assert!(
                run.stdout
                    .contains(&format!("\"pack\":\"stainless-{grade}-nist-cryogenic\""))
            );
        }
        assert!(run.stdout.contains("\"unknown_properties\":[]"));
        println!("{}", run.stdout);
    }

    #[test]
    fn g0_g3_solid_lead_source_conditions() {
        let (pack, _) = compile("lead-cast-expansion-nbs-rp500");
        check_claims(&pack);
        assert_eq!(pack.claims().claims_ordered().count(), 5);
        let density = pack.claims().claims_for("density")[0].1;
        close(
            sample(&pack, "density", &[("temperature", 298.15)]),
            11310.0,
        );
        assert!(density.validity.bound("source-heating").is_none());
        assert!(
            pack.claims()
                .query_typed(
                    &density.key,
                    &point(density, &[("temperature", 299.15)]),
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
        let name = "mean-linear-expansion-coefficient-from-20c";
        let claims = pack.claims().claims_for(name);
        assert_eq!(claims.len(), 4);
        for (end, mean) in [
            (333.15, 28.3e-6),
            (373.15, 28.6e-6),
            (473.15, 29.5e-6),
            (573.15, 31.2e-6),
        ] {
            close(
                sample(&pack, name, &[("source-range-end-temperature", end)]),
                mean,
            );
        }
        for (_, claim) in claims {
            assert_eq!(claim.key.dims(), Dims([0, 0, 0, -1, 0, 0]));
            assert_eq!(
                claim.interpolation,
                InterpolationPolicy::ConstantWithinValidity
            );
            assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
            assert_eq!(
                claim.validity.bound("source-reference-temperature"),
                Some((293.15, 293.15))
            );
            for axis in [
                "source-reference-temperature",
                "source-range-end-temperature",
            ] {
                assert_eq!(
                    claim.validity.axis_quantities()[axis]
                        .semantic_type()
                        .unwrap()
                        .kind(),
                    fs_qty::semantic::QuantityKind::AbsoluteTemperature
                );
            }
            for overrides in [
                vec![("source-reference-temperature", 298.15)],
                vec![("source-range-end-temperature", 400.0)],
                vec![("source-sample-id", 1215.0)],
                vec![("source-heating", 0.0)],
                vec![("source-pressure-known", 1.0)],
            ] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &claim.key,
                            &point(claim, &overrides),
                            SelectionPolicy::SingleClaimOnly,
                        )
                        .is_err()
                );
            }
        }
        // The source's interval means cannot fill an instantaneous alpha law.
        assert!(
            pack.claims()
                .claims_for("linear-thermal-expansion-coefficient")
                .is_empty()
        );
        let (conductivity, _) = compile("lead-solid-conductivity-nbs-rp668");
        check_claims(&conductivity);
        let claim = conductivity.claims().claims_for("thermal-conductivity")[0].1;
        for (temperature, expected) in [
            (273.15, 35.2),
            (373.15, 33.2),
            (473.15, 31.2),
            (573.15, 29.2),
            (423.15, 32.2),
        ] {
            close(
                sample(
                    &conductivity,
                    "thermal-conductivity",
                    &[("temperature", temperature)],
                ),
                expected,
            );
        }
        assert_eq!(claim.validity.bound("temperature"), Some((273.15, 573.15)));
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        for overrides in [
            vec![("temperature", 273.14)],
            vec![("temperature", 573.16)],
            vec![("source-calibration-assumed", 0.0)],
            vec![("source_phase_liquid", 1.0)],
            vec![("source-pressure-known", 1.0)],
        ] {
            assert!(
                conductivity
                    .claims()
                    .query_typed(
                        &claim.key,
                        &point(claim, &overrides),
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .is_err()
            );
        }
        // A differently prepared conductivity standard supplies no density
        // measurement for the RP500 specimen or for its own temperature range.
        assert!(conductivity.claims().claims_for("density").is_empty());
    }

    #[test]
    fn g0_g3_stainless_thermomechanical_envelope() {
        let mut paths = Vec::new();
        for (grade, upper, coefficients, e_150) in [
            (
                "304",
                293,
                [
                    210.0593,
                    0.1534883,
                    -0.001617390,
                    0.000005117060,
                    -0.0000000061546,
                ],
                210.84558125e9,
            ),
            (
                "316",
                294,
                [
                    207.9488,
                    0.07394241,
                    -0.0009627200,
                    0.000002845560,
                    -0.0000000032408,
                ],
                205.3420715e9,
            ),
        ] {
            let (pack, path) = compile(&format!("stainless-{grade}-nist-cryogenic"));
            paths.push(path);
            check_claims(&pack);
            close(
                sample(&pack, "young-modulus", &[("temperature", 150.0)]),
                e_150,
            );
            for (name, high) in [
                ("young-modulus", f64::from(upper)),
                ("linear-expansion-relative-to-293k", 300.0),
                ("linear-thermal-expansion-coefficient", 300.0),
            ] {
                let claims = pack.claims().claims_for(name);
                assert_eq!(claims.len(), 1, "old points must not make lookup ambiguous");
                let claim = claims[0].1;
                assert_eq!(claim.interpolation, InterpolationPolicy::LinearInside);
                assert_eq!(claim.validity.bound("temperature"), Some((77.0, high)));
                assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                for overrides in [
                    vec![("temperature", 76.0)],
                    vec![("temperature", high + 0.01)],
                    vec![("temperature", 150.0), ("source-pressure-known", 1.0)],
                ] {
                    assert!(
                        pack.claims()
                            .query_typed(
                                &claim.key,
                                &point(claim, &overrides),
                                SelectionPolicy::SingleClaimOnly,
                            )
                            .is_err()
                    );
                }
            }
            // Published polynomial evaluated separately from the tabulated
            // curve. These sampled numerical checks are not physical bounds.
            for integer in 77..upper {
                let t = f64::from(integer) + 0.5;
                let e_fit = coefficients
                    .iter()
                    .enumerate()
                    .map(|(power, coefficient)| coefficient * t.powi(power as i32))
                    .sum::<f64>()
                    * 1e9;
                let interpolated = sample(&pack, "young-modulus", &[("temperature", t)]);
                assert!(
                    (interpolated / e_fit - 1.0).abs() < 0.0004,
                    "{grade} at {t} K"
                );
            }
            for integer in 77..300 {
                let t = f64::from(integer) + 0.5;
                let expansion_fit = (-295.54 - 0.39811 * t + 0.0092683 * t.powi(2)
                    - 0.000020261 * t.powi(3)
                    + 0.000000017127 * t.powi(4))
                    * 1e-5;
                let interpolated = sample(
                    &pack,
                    "linear-expansion-relative-to-293k",
                    &[("temperature", t)],
                );
                assert!(
                    (interpolated - expansion_fit).abs() < 1e-5,
                    "{grade} at {t} K"
                );
            }
            let expansion = pack
                .claims()
                .claims_for("linear-expansion-relative-to-293k")[0]
                .1;
            assert_eq!(expansion.key.dims(), Dims::NONE);
            close(
                sample(&pack, expansion.key.name(), &[("temperature", 150.0)]),
                -0.0020643008125,
            );
            close(
                sample(&pack, expansion.key.name(), &[("temperature", 293.0)]),
                7.4646191727e-7,
            );
            // Alpha is the derivative of log length, not a relabeling of
            // relative strain. NIST does not provide its uncertainty bound.
            let alpha = pack
                .claims()
                .claims_for("linear-thermal-expansion-coefficient")[0]
                .1;
            assert_eq!(alpha.key.dims(), Dims([0, 0, 0, -1, 0, 0]));
            let epsilon = |t: f64| {
                (-295.54 - 0.39811 * t + 0.0092683 * t.powi(2) - 0.000020261 * t.powi(3)
                    + 0.000000017127 * t.powi(4))
                    * 1e-5
            };
            let instantaneous = |t: f64| {
                (-0.39811 + 0.0185366 * t - 0.000060783 * t.powi(2) + 0.000000068508 * t.powi(3))
                    * 1e-5
                    / (1.0 + epsilon(t))
            };
            let PropertyValue::Curve { knots, .. } = &alpha.value else {
                panic!("instantaneous expansion must be a bounded curve")
            };
            assert_eq!(knots.len(), 11);
            for &(t, value) in knots {
                assert!((value / instantaneous(t) - 1.0).abs() < 1e-12);
                if t < 300.0 {
                    let step = 1e-3;
                    let derivative =
                        (epsilon(t + step).ln_1p() - epsilon(t - step).ln_1p()) / (2.0 * step);
                    assert!((value / derivative - 1.0).abs() < 2e-8);
                }
            }
            let mut maximum = 0.0_f64;
            for pair in knots.windows(2) {
                for i in 1..100 {
                    let t = pair[0].0 + (pair[1].0 - pair[0].0) * f64::from(i) / 100.0;
                    let value = sample(&pack, alpha.key.name(), &[("temperature", t)]);
                    maximum = maximum.max((value / instantaneous(t) - 1.0).abs());
                }
            }
            assert!(
                maximum < 0.007,
                "{grade}: sampled alpha discrepancy {maximum}"
            );
            // The real card/store/resolver preserves this derived claim and
            // its source context; it does not invent the missing rho or nu.
            let card = fs_matdb::NormalizedMaterialCardPack::new(
                fs_matdb::MaterialStateId {
                    chemistry: format!("UNS S{grade}00"),
                    phase: "solid".into(),
                    process: "NIST fit; product form and heat treatment unstated".into(),
                    revision: 0,
                },
                pack.clone(),
            )
            .unwrap();
            let store = fs_matdb_store::MaterialStore::open(":memory:").unwrap();
            store
                .ingest_bundle(&[fs_matdb_store::CatalogPack::MaterialCard(card.clone())])
                .unwrap();
            store.seal_corpus().unwrap();
            let fs_matdb_store::CatalogPack::MaterialCard(loaded) =
                store.load_catalog_pack(card.pack_id()).unwrap()
            else {
                panic!("wrong stored family")
            };
            assert_eq!(loaded, card);
            let state = fs_material::state_point::resolve_material_state_point(
                loaded.card(),
                &point(alpha, &[("temperature", 150.0)]),
                &[
                    fs_material::state_point::ScalarPropertyRequirement::try_with_key(
                        &alpha.key,
                        fs_material::state_point::ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap(),
                ],
                fs_material::state_point::MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let selected = state.property(alpha.key.name()).unwrap();
            assert!((selected.value_si() / instantaneous(150.0) - 1.0).abs() < 1e-12);
            loaded
                .card()
                .claims()
                .verify_receipt(&selected.answer().receipt)
                .unwrap();
            println!("{grade}: 11 derived-alpha knots, sampled maximum discrepancy {maximum}");
            assert!(pack.claims().claims_for("density").is_empty());
        }
        let request = fixture_dir().join("stainless-thermomechanical.json");
        let source =
            include_str!("../../examples/material-discovery/stainless-thermomechanical.json");
        for (upper, complete, supported) in [(293, 2, 10), (294, 1, 9), (295, 0, 8)] {
            fs::write(&request, source.replace("293 K", &format!("{upper} K"))).unwrap();
            let run = fs_cli::run(vec![
                "--json".into(),
                "discover".into(),
                request.to_str().unwrap().into(),
                paths[0].to_str().unwrap().into(),
                paths[1].to_str().unwrap().into(),
            ]);
            assert_eq!(run.exit_code, fs_cli::exit::SUCCESS, "{}", run.stderr);
            assert_eq!(
                run.stdout.matches("\"status\":\"complete\"").count(),
                complete,
                "{}",
                run.stdout
            );
            assert_eq!(
                run.stdout.matches("\"status\":\"supported\"").count(),
                supported,
                "{}",
                run.stdout
            );
            assert_eq!(
                run.stdout.matches("\"status\":\"gap\"").count(),
                10 - supported,
                "{}",
                run.stdout
            );
            assert!(
                run.stdout.contains("\"unknown_properties\":[]"),
                "{}",
                run.stdout
            );
            println!("{upper} K envelope: {}", run.stdout);
        }
    }

    #[test]
    fn g0_g3_common_metals_and_construction_sources() {
        // These expected SI values come from the independently reviewed source
        // tables/equations, not a generated golden from the pack under test.
        let cases = [
            (
                "aluminum-pure-nasa-cr-71699",
                "thermal-conductivity",
                Some(300.0),
                237.0,
            ),
            (
                "copper-pure-nasa-cr-71699",
                "thermal-conductivity",
                Some(300.0),
                398.0,
            ),
            (
                "iron-pure-nasa-cr-71699",
                "thermal-conductivity",
                Some(300.0),
                80.0,
            ),
            (
                "nickel-pure-nasa-cr-71699",
                "thermal-conductivity",
                Some(350.0),
                83.0,
            ),
            (
                "titanium-pure-nasa-cr-71699",
                "thermal-conductivity",
                Some(300.0),
                21.9,
            ),
            (
                "aluminum-7075-t6-nasa-cr-71699",
                "thermal-conductivity",
                Some(500.0),
                178.0,
            ),
            (
                "inconel-x750-nasa-cr-71699",
                "thermal-conductivity",
                Some(300.0),
                11.7,
            ),
            (
                "stainless-304a-nasa-cr-71699",
                "thermal-conductivity",
                Some(500.0),
                18.4,
            ),
            (
                "stainless-347-nasa-cr-71699",
                "thermal-conductivity",
                Some(300.0),
                14.8,
            ),
            (
                "titanium-a110at-nasa-cr-71699",
                "thermal-conductivity",
                Some(500.0),
                9.8,
            ),
            (
                "stainless-304-nist-cryogenic",
                "thermal-conductivity",
                Some(293.0),
                15.1233439431416,
            ),
            (
                "stainless-316-nist-cryogenic",
                "specific-heat-capacity",
                Some(293.0),
                485.317885527169,
            ),
            (
                "brass-c26000-nist-cryogenic",
                "thermal-conductivity",
                Some(77.0),
                39.7957752285347,
            ),
            (
                "aluminum-1100-nist-cryogenic",
                "thermal-conductivity",
                Some(293.0),
                211.815871273228,
            ),
            (
                "concrete-normalweight-cfast-sp1041",
                "thermal-conductivity",
                None,
                1.75,
            ),
            (
                "concrete-lightweight-cfast-sp1041",
                "thermal-conductivity",
                None,
                0.125,
            ),
            (
                "cement-mortar-cfast-sp1041",
                "thermal-conductivity",
                None,
                0.72,
            ),
            ("brick-clay-cfast-sp1041", "thermal-conductivity", None, 1.5),
            (
                "brick-common-cfast-sp1041",
                "thermal-conductivity",
                None,
                0.72,
            ),
            (
                "calcium-silicate-board-cfast-sp1041",
                "thermal-conductivity",
                None,
                0.18,
            ),
            (
                "cellulose-insulation-cfast-sp1041",
                "thermal-conductivity",
                None,
                0.039,
            ),
            (
                "glass-fiber-insulation-cfast",
                "thermal-conductivity",
                None,
                0.04,
            ),
            ("gypsum-board-5-8-cfast", "thermal-conductivity", None, 0.16),
            (
                "gypsum-board-type-x-5-8-cfast",
                "thermal-conductivity",
                None,
                0.14,
            ),
            (
                "urethane-rigid-foam-insulation-cfast",
                "thermal-conductivity",
                None,
                0.026,
            ),
            (
                "concrete-nsc-mix-iv-nistir6475",
                "compressive-strength",
                Some(298.15),
                51.9e6,
            ),
        ];
        for (slug, name, temperature, expected) in cases {
            let (pack, _) = compile(slug);
            check_claims(&pack);
            let at = temperature
                .map(|t| ("temperature", t))
                .into_iter()
                .collect::<Vec<_>>();
            close(sample(&pack, name, &at), expected);
            if matches!(
                slug,
                "aluminum-pure-nasa-cr-71699" | "copper-pure-nasa-cr-71699"
            ) {
                let claim = pack.claims().claims_for(name)[0].1;
                assert!(
                    pack.claims()
                        .query_typed(
                            &claim.key,
                            &point(claim, &[("source_phase_liquid", 1.0)]),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
                let lo = claim.validity.bound("temperature").unwrap().0;
                assert!(
                    pack.claims()
                        .query_typed(
                            &claim.key,
                            &point(claim, &[("temperature", lo - 1.0)]),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            println!("source pack {slug}: compiler, typed query, context refusal, receipt pass");
        }
    }

    /// G1/G3: the source-resolved conductor adapter reaches the actual
    /// circuit DAE. These are bulk resistivity source checks at the two
    /// retained wire-table states, not a wire-temperature or contact model.
    #[test]
    fn g1_g3_sourced_conductors_reach_circuit_dissipation() {
        use fs_blake3::ContentHash;
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, PropertyKey, PropertyValue};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::conductor::{
            ConductorError, ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY,
            resolve_uniform_conductor,
        };
        use fs_material::state_point::MaterialPropertySelection;
        use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
        use fs_qty::semantic::{QuantityKind, QuantitySpec, SemanticType, ValueForm};

        const LENGTH_M: f64 = 1.0;
        const AREA_M2: f64 = 1.0e-6;
        const CURRENT_A: f64 = 2.0;
        const COPPER: &str = "copper-annealed-iacs-nbs-hb100";
        const ALUMINUM: &str = "aluminum-ec-h19-nbs-hb109";

        // HB-100 prints copper resistivity as 0.017241 ohm mm²/m at 20 C,
        // with 0.0000681 ohm mm²/m/C slope. These are the retained source
        // table relation's three derived SI knots; no reciprocal conductivity or
        // generic temperature coefficient is substituted.
        let cases: [(&str, &[f64], &[f64]); 2] = [
            (
                COPPER,
                &[283.15, 293.15, 303.15],
                &[1.6560e-8, 1.7241e-8, 1.7922e-8],
            ),
            (
                ALUMINUM,
                &[273.15, 288.15, 293.15, 298.15, 303.15],
                &[2.5986e-8, 2.7695e-8, 2.8264e-8, 2.8834e-8, 2.9403e-8],
            ),
        ];
        let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let resistivity_key = PropertyKey::with_quantity(
            ELECTRICAL_RESISTIVITY_PROPERTY,
            QuantitySpec::dimensional(ELECTRICAL_RESISTIVITY_DIMS),
        );
        let close_relative = |actual: f64, expected: f64| {
            assert!(
                (actual / expected - 1.0).abs() <= 1.0e-10,
                "actual {actual} differs from source-derived expected {expected}"
            );
        };

        for (slug, temperatures, expected_rho) in cases {
            let (pack, pack_path) = compile(slug);
            let claims = pack.claims();
            let resistivity_claims = claims.claims_for(ELECTRICAL_RESISTIVITY_PROPERTY);
            assert_eq!(
                resistivity_claims.len(),
                1,
                "{slug}: exactly one resistivity claim required"
            );
            let (_, claim) = resistivity_claims[0];
            assert_eq!(claim.key, resistivity_key);
            let PropertyValue::Curve {
                abscissa, knots, ..
            } = &claim.value
            else {
                panic!("{slug}: electrical resistivity must be a curve");
            };
            assert_eq!(abscissa, "T");
            assert_eq!(knots.len(), temperatures.len());
            for (&(actual_t, actual_rho), (&expected_t, &expected_value)) in
                knots.iter().zip(temperatures.iter().zip(expected_rho))
            {
                assert_eq!(
                    actual_t, expected_t,
                    "{slug}: source temperature knot moved"
                );
                close_relative(actual_rho, expected_value);
            }
            assert!(
                claim.validity.axis_quantities().get("T") == Some(&absolute_temperature),
                "{slug}: T must retain absolute-temperature semantics"
            );
            assert_eq!(
                claim.validity.bound("source-pressure-known"),
                Some((0.0, 0.0))
            );

            let card = NormalizedMaterialCardPack::new(
                MaterialStateId {
                    chemistry: format!("{slug}; wire-table source state"),
                    phase: "solid conductor".into(),
                    process: "source-table wire temper retained; no thermal evolution".into(),
                    revision: 0,
                },
                pack.clone(),
            )
            .unwrap();
            let card_path = fixture_dir().join(format!("{slug}.fsmcdpk"));
            let card_bytes = card.to_bytes();
            fs::write(&card_path, &card_bytes).unwrap();
            let decoded =
                NormalizedMaterialCardPack::from_bytes_verified(card.content_hash(), &card_bytes)
                    .expect("round-trip card");
            assert_eq!(decoded, card);

            let database = fixture_dir().join(format!("{slug}.sqlite"));
            {
                let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
                store
                    .ingest_bundle(&[CatalogPack::MaterialCard(card.clone())])
                    .unwrap();
                store.seal_corpus().unwrap();
                let stored = store
                    .evaluate_typed(
                        card.pack_id(),
                        &resistivity_key,
                        &query_point(claim, temperatures[0]),
                        fs_matdb::SelectionPolicy::SingleClaimOnly,
                    )
                    .unwrap();
                close_relative(stored.evidence.value.value, expected_rho[0]);
            }
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            let CatalogPack::MaterialCard(loaded) =
                store.load_catalog_pack(card.pack_id()).unwrap()
            else {
                panic!("{slug}: stored family changed");
            };
            assert_eq!(loaded, card);

            let graph = CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::CurrentSource { port: 0 }),
                    (1, 0, Branch::Resistor { ohms: 1.0 }),
                ],
                transformers: vec![],
            };
            let mut operating_states: Vec<_> = temperatures
                .iter()
                .copied()
                .zip(expected_rho.iter().copied())
                .collect();
            // Interior state independently checks the declared linear
            // interpolant; it is not an additional source measurement.
            operating_states.push((
                (temperatures[0] + temperatures[1]) * 0.5,
                (expected_rho[0] + expected_rho[1]) * 0.5,
            ));
            for (temperature, rho) in operating_states {
                let point = query_point(
                    loaded
                        .card()
                        .claims()
                        .claims_for(ELECTRICAL_RESISTIVITY_PROPERTY)[0]
                        .1,
                    temperature,
                );
                let conductor = resolve_uniform_conductor(
                    loaded.card(),
                    &point,
                    MaterialPropertySelection::SingleClaimOnly,
                    LENGTH_M,
                    AREA_M2,
                )
                .unwrap_or_else(|error| panic!("{slug} at {temperature} K: {error}"));
                let expected_resistance = rho * LENGTH_M / AREA_M2;
                close_relative(conductor.resistance_ohm(), expected_resistance);
                close_relative(
                    conductor.joule_power_w(CURRENT_A).unwrap(),
                    CURRENT_A * CURRENT_A * expected_resistance,
                );
                assert_eq!(conductor.joule_power_w(0.0).unwrap(), 0.0);
                assert_eq!(
                    conductor.joule_power_w(-CURRENT_A).unwrap(),
                    conductor.joule_power_w(CURRENT_A).unwrap()
                );
                assert!(conductor.joule_power_w(f64::NAN).is_err());
                let resolved_rho = conductor
                    .material()
                    .property(ELECTRICAL_RESISTIVITY_PROPERTY)
                    .expect("resistivity in resolved conductor");
                loaded
                    .card()
                    .claims()
                    .verify_receipt(&resolved_rho.answer().receipt)
                    .expect("source receipt verifies");
                let scaled = resolve_uniform_conductor(
                    loaded.card(),
                    &point,
                    MaterialPropertySelection::SingleClaimOnly,
                    2.0 * LENGTH_M,
                    0.5 * AREA_M2,
                )
                .expect("scaled uniform geometry");
                close_relative(scaled.resistance_ohm(), 4.0 * expected_resistance);

                let dae = assemble_circuit(&CircuitGraph {
                    branches: vec![
                        (1, 0, Branch::CurrentSource { port: 0 }),
                        (
                            1,
                            0,
                            Branch::Resistor {
                                ohms: conductor.resistance_ohm(),
                            },
                        ),
                    ],
                    ..graph.clone()
                })
                .expect("current-source/resistor circuit admits");
                let zero = dae
                    .consistent_initial_state(&vec![0.0; dae.system.state_dim()], &[0.0])
                    .expect("zero-source algebraic state");
                let (record, defect) = dae
                    .step_audited(&zero, &[CURRENT_A], 1.0)
                    .expect("one-second DAE step");
                let voltage = record.y[0];
                let midpoint_voltage = (zero[dae.node_potential_index[0]]
                    + record.x[dae.node_potential_index[0]])
                    * 0.5;
                close_relative(voltage, midpoint_voltage);
                close_relative(voltage, CURRENT_A * expected_resistance);
                close_relative(
                    record.dissipated,
                    CURRENT_A * CURRENT_A * expected_resistance,
                );
                close_relative(record.supplied, CURRENT_A * voltage);
                assert!(
                    defect <= 1.0e-10,
                    "{slug} at {temperature} K: audit {defect}"
                );
                assert!(record.solver_residual.is_finite());
                assert!(record.solver_residual <= 1.0e-10);
                println!(
                    "conductor source={slug} T={temperature} K rho_expected={rho} Ohm*m R_actual={} Ohm R_expected={expected_resistance} Ohm I={CURRENT_A} A dt=1 s voltage={voltage} V dissipated={} J supplied={} J residual={} supply_defect={defect} J relative_tolerance=1e-10 claim={:?}",
                    conductor.resistance_ohm(),
                    record.dissipated,
                    record.supplied,
                    record.solver_residual,
                    resolved_rho.answer().receipt.selected,
                );
                let (replayed_record, replayed_defect) = dae
                    .step_audited(&zero, &[CURRENT_A], 1.0)
                    .expect("deterministic DAE replay");
                assert_eq!(
                    replayed_record.x, record.x,
                    "{slug}: DAE state replay moved"
                );
                assert_eq!(
                    replayed_record.y, record.y,
                    "{slug}: DAE voltage replay moved"
                );
                assert_eq!(
                    replayed_record.dissipated, record.dissipated,
                    "{slug}: DAE dissipation replay moved"
                );
                assert_eq!(replayed_defect, defect, "{slug}: DAE audit replay moved");
                let replay = resolve_uniform_conductor(
                    loaded.card(),
                    &point,
                    MaterialPropertySelection::SingleClaimOnly,
                    LENGTH_M,
                    AREA_M2,
                )
                .expect("deterministic conductor replay");
                assert_eq!(replay, conductor, "{slug} conductor replay moved");
            }

            let source_claim = loaded
                .card()
                .claims()
                .claims_for(ELECTRICAL_RESISTIVITY_PROPERTY)[0]
                .1;
            let valid_point = query_point(source_claim, temperatures[0]);
            assert!(matches!(
                resolve_uniform_conductor(
                    loaded.card(),
                    &valid_point,
                    MaterialPropertySelection::SingleClaimOnly,
                    0.0,
                    AREA_M2,
                ),
                Err(ConductorError::InvalidInput {
                    quantity: "length_m"
                })
            ));
            assert!(matches!(
                resolve_uniform_conductor(
                    loaded.card(),
                    &valid_point,
                    MaterialPropertySelection::SingleClaimOnly,
                    LENGTH_M,
                    f64::NAN,
                ),
                Err(ConductorError::InvalidInput {
                    quantity: "area_m2"
                })
            ));
            let outside = query_point(source_claim, temperatures[0] - 1.0);
            assert!(
                resolve_uniform_conductor(
                    loaded.card(),
                    &outside,
                    MaterialPropertySelection::SingleClaimOnly,
                    LENGTH_M,
                    AREA_M2,
                )
                .is_err()
            );
            let wrong_unit = PropertyKey::with_quantity(
                ELECTRICAL_RESISTIVITY_PROPERTY,
                fs_qty::QuantitySpec::dimensional(Dims::NONE),
            );
            assert!(
                store
                    .evaluate_typed(
                        card.pack_id(),
                        &wrong_unit,
                        &valid_point,
                        fs_matdb::SelectionPolicy::SingleClaimOnly,
                    )
                    .is_err()
            );
            let missing_card = fs_matdb::MaterialCard::assemble(
                MaterialStateId {
                    chemistry: format!("{slug}; missing resistivity"),
                    phase: "solid conductor".into(),
                    process: "deliberate missing-property refusal".into(),
                    revision: 0,
                },
                fs_matdb::ClaimSet::new(),
                Vec::new(),
            )
            .unwrap();
            assert!(matches!(
                resolve_uniform_conductor(
                    &missing_card,
                    &valid_point,
                    MaterialPropertySelection::SingleClaimOnly,
                    LENGTH_M,
                    AREA_M2,
                ),
                Err(ConductorError::MaterialState(_))
            ));
            let mut wrong_claims = fs_matdb::ClaimSet::new();
            for observation_id in loaded.card().claims().observation_ids() {
                wrong_claims
                    .register_observation(
                        loaded
                            .card()
                            .claims()
                            .observation(observation_id)
                            .expect("source observation")
                            .clone(),
                    )
                    .unwrap();
            }
            let wrong_value = match &source_claim.value {
                PropertyValue::Curve {
                    abscissa,
                    abscissa_dims,
                    knots,
                    ..
                } => PropertyValue::Curve {
                    abscissa: abscissa.clone(),
                    abscissa_dims: *abscissa_dims,
                    knots: knots.clone(),
                    dims: Dims::NONE,
                },
                PropertyValue::Scalar { value, .. } => PropertyValue::Scalar {
                    value: *value,
                    dims: Dims::NONE,
                },
            };
            let wrong_claim = fs_matdb::PropertyClaim {
                key: PropertyKey::with_quantity(
                    ELECTRICAL_RESISTIVITY_PROPERTY,
                    fs_qty::QuantitySpec::dimensional(Dims::NONE),
                ),
                value: wrong_value,
                ..source_claim.clone()
            };
            wrong_claims.insert_claim(wrong_claim).unwrap();
            let wrong_unit_card = fs_matdb::MaterialCard::assemble(
                MaterialStateId {
                    chemistry: format!("{slug}; wrong resistivity unit"),
                    phase: "solid conductor".into(),
                    process: "deliberate dimension refusal".into(),
                    revision: 0,
                },
                wrong_claims,
                Vec::new(),
            )
            .unwrap();
            assert!(matches!(
                resolve_uniform_conductor(
                    &wrong_unit_card,
                    &valid_point,
                    MaterialPropertySelection::SingleClaimOnly,
                    LENGTH_M,
                    AREA_M2,
                ),
                Err(ConductorError::MaterialState(_))
            ));
            assert!(
                resolve_uniform_conductor(
                    loaded.card(),
                    &valid_point,
                    MaterialPropertySelection::PinnedByProperty(vec![(
                        ELECTRICAL_RESISTIVITY_PROPERTY.into(),
                        fs_matdb::ClaimId(ContentHash([0; 32])),
                    )]),
                    LENGTH_M,
                    AREA_M2,
                )
                .is_err()
            );

            let request_path = workspace_path(&format!(
                "examples/material-discovery/{}-dc-conductor.json",
                if slug == COPPER { "copper" } else { "aluminum" }
            ));
            assert!(
                request_path.is_file(),
                "missing discovery request {request_path:?}"
            );
            let discovery = fs_cli::run(vec![
                "--json".into(),
                "discover".into(),
                request_path.to_str().unwrap().into(),
                pack_path.to_str().unwrap().into(),
            ]);
            assert_eq!(
                discovery.exit_code,
                fs_cli::exit::SUCCESS,
                "{}",
                discovery.stderr
            );
            assert!(
                discovery.stdout.contains("\"status\":\"complete\""),
                "{}",
                discovery.stdout
            );
            println!(
                "{slug}: source_rho={expected_rho:?} temperature_K={temperatures:?} geometry_L_m={LENGTH_M} geometry_A_m2={AREA_M2} current_A={CURRENT_A} discovery=complete"
            );
        }
    }

    fn query_point(claim: &PropertyClaim, temperature: f64) -> QueryPoint {
        let mut point = QueryPoint::new();
        for (axis, &(lo, _)) in claim.validity.bounds() {
            let value = if axis == "T" { temperature } else { lo };
            point = if let Some(quantity) = claim.validity.axis_quantities().get(axis) {
                point.with_quantity(axis, *quantity, value)
            } else {
                point.with(axis, value)
            }
            .expect("typed source query point");
        }
        point
    }

    /// G1/G3: two independently sourced component cards reach the actual
    /// humid-gas transport and cylinder-loss consumer. The Buck, ideal-mixture,
    /// Wilke, and WMS calculations below are independent test oracles; this is
    /// a bounded engineering mixture check, not humid-air qualification.
    #[test]
    fn g1_g3_sourced_humid_air_reaches_acoustic_transport() {
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::gas::{ConductivityModel, R_USSA_1976, resolve_sutherland_gas_state};
        use fs_material::moist_air::resolve_moist_air_state;
        use fs_material::state_point::MaterialPropertySelection;
        use fs_qty::semantic::{
            QuantityKind, QuantitySpec as SemanticQuantitySpec, SemanticType, ValueForm,
        };
        use fs_qty::{Dims, QuantitySpec};

        const DRY: &str = "air-dry-ussa1976";
        const VAPOR: &str = "water-vapor-sutherland-ambient";
        const MD: f64 = 28.9644e-3;
        const MV: f64 = 18.01528e-3;
        const GD: f64 = 1.4;
        const GV: f64 = 33.590 / (33.590 - R_USSA_1976);
        const BETA_D: f64 = 1.458e-6;
        const S_D: f64 = 110.4;
        const MU_REF_V: f64 = 1.12e-5;
        const T_REF_V: f64 = 350.0;
        const S_V: f64 = 1064.0;
        const RADIUS_M: f64 = 0.0005;
        const OMEGA_RAD_S: f64 = core::f64::consts::TAU * 220.0;

        let (dry_pack, dry_path) = compile(DRY);
        assert_eq!(dry_pack.pack_id(), "air-dry-ussa1976-ambient-model");
        let (vapor_pack, vapor_path) = compile(VAPOR);
        assert_eq!(vapor_pack.pack_id(), VAPOR);
        let dry_card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "dry air; USSA-1976 reference composition".into(),
                phase: "gas".into(),
                process: "USSA-1976 dry-air constants and transport fits".into(),
                revision: 0,
            },
            dry_pack.clone(),
        )
        .unwrap();
        let vapor_card = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "water vapor; Sutherland ambient reference".into(),
                phase: "gas".into(),
                process: "NIST component source and declared Sutherland/Eucken model".into(),
                revision: 0,
            },
            vapor_pack.clone(),
        )
        .unwrap();
        let database = fixture_dir().join("humid-air-components.sqlite");
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[
                    CatalogPack::MaterialCard(dry_card.clone()),
                    CatalogPack::MaterialCard(vapor_card.clone()),
                ])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(loaded_dry) =
            store.load_catalog_pack(dry_card.pack_id()).unwrap()
        else {
            panic!("dry component changed family");
        };
        let CatalogPack::MaterialCard(loaded_vapor) =
            store.load_catalog_pack(vapor_card.pack_id()).unwrap()
        else {
            panic!("vapor component changed family");
        };
        assert_eq!(loaded_dry, dry_card);
        assert_eq!(loaded_vapor, vapor_card);

        let absolute_temperature = SemanticQuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let pressure = SemanticQuantitySpec::semantic(SemanticType::new(
            QuantityKind::Pressure,
            ValueForm::Static,
        ));
        let dimensionless = QuantitySpec::dimensional(Dims::NONE);
        let component_point = |claim: &PropertyClaim, t: f64, p: f64, flag: &str| {
            let mut point = QueryPoint::new();
            for (axis, &(lo, _)) in claim.validity.bounds() {
                let value = match axis.as_str() {
                    "temperature" => t,
                    "pressure" => p,
                    "relative-humidity" => 0.0,
                    name if name == flag => 1.0,
                    _ => lo,
                };
                let quantity = claim
                    .validity
                    .axis_quantities()
                    .get(axis)
                    .copied()
                    .unwrap_or_else(|| match axis.as_str() {
                        "temperature" => absolute_temperature,
                        "pressure" => pressure,
                        "relative-humidity"
                        | "source-composition-ussa1976"
                        | "source-component-water-vapor" => dimensionless,
                        _ => QuantitySpec::dimensional(Dims::NONE),
                    });
                point = point
                    .with_quantity(axis, quantity, value)
                    .expect("typed humid-air component point");
            }
            point
        };
        let dry_claim = dry_pack.claims().claims_for("molar_mass")[0].1;
        let vapor_claim = vapor_pack.claims().claims_for("molar_mass")[0].1;
        let dry_point =
            |t: f64, p: f64| component_point(dry_claim, t, p, "source-composition-ussa1976");
        let vapor_point =
            |t: f64, p: f64| component_point(vapor_claim, t, p, "source-component-water-vapor");
        let relative_close = |actual: f64, expected: f64| {
            assert!(
                (actual / expected - 1.0).abs() <= 1.0e-10,
                "actual {actual:.16e}, expected {expected:.16e}"
            );
        };
        let saturation_pressure = |temperature_k: f64| {
            let celsius = temperature_k - 273.15;
            611.21 * f64::exp((18.678 - celsius / 234.5) * celsius / (257.14 + celsius))
        };
        let sutherland = |beta: f64, s: f64, t: f64| beta * t.powf(1.5) / (t + s);
        let expected = |t: f64, p: f64, rh: f64| {
            let x = rh * saturation_pressure(t) / p;
            let mm = MD + x * (MV - MD);
            let r_v = R_USSA_1976 / MV;
            let r_mix = R_USSA_1976 / mm;
            let inv_gamma_minus_one = (1.0 - x) / (GD - 1.0) + x / (GV - 1.0);
            let gamma = 1.0 + 1.0 / inv_gamma_minus_one;
            let cp_v = 33.590 / MV;
            let mu_d = sutherland(BETA_D, S_D, t);
            let mu_v = MU_REF_V * (t / T_REF_V).powf(1.5) * (T_REF_V + S_V) / (t + S_V);
            let t32 = t * t.sqrt();
            let k_d =
                2.64638e-3 * t32 / (t + 245.4 * f64::exp((-12.0 / t) * core::f64::consts::LN_10));
            let k_v = mu_v * (cp_v + 1.25 * r_v);
            let phi = |mi: f64, mj: f64, mui: f64, muj: f64| {
                (1.0 + (mui / muj).sqrt() * (mj / mi).powf(0.25)).powi(2)
                    / (8.0 * (1.0 + mi / mj)).sqrt()
            };
            let phi_dv = phi(MD, MV, mu_d, mu_v);
            let phi_vd = phi(MV, MD, mu_v, mu_d);
            let mu =
                (1.0 - x) * mu_d / ((1.0 - x) + x * phi_dv) + x * mu_v / ((1.0 - x) * phi_vd + x);
            let k = (1.0 - x) * k_d / ((1.0 - x) + x * phi_dv) + x * k_v / ((1.0 - x) * phi_vd + x);
            let density = p / (r_mix * t);
            let sound_speed = (gamma * r_mix * t).sqrt();
            let cp = gamma * r_mix / (gamma - 1.0);
            let prandtl = mu * cp / k;
            let loss = core::f64::consts::TAU
                * mu
                * (1.0 + RADIUS_M * (2.0 * density * OMEGA_RAD_S / mu).sqrt());
            (
                x,
                x * MV / mm,
                gamma,
                density,
                sound_speed,
                mu,
                k,
                cp,
                prandtl,
                loss,
            )
        };

        for &(t, p, rh) in &[
            (273.15, 110_000.0, 0.8),
            (293.15, 101_325.0, 0.3),
            (293.15, 101_325.0, 0.8),
            (293.15, 101_325.0, 1.0),
            (313.15, 80_000.0, 0.8),
        ] {
            let dry_at = dry_point(t, p);
            let vapor_at = vapor_point(t, p);
            let resolved = resolve_moist_air_state(
                loaded_dry.card(),
                &dry_at,
                MaterialPropertySelection::SingleClaimOnly,
                loaded_vapor.card(),
                &vapor_at,
                MaterialPropertySelection::SingleClaimOnly,
                rh,
            )
            .unwrap_or_else(|error| panic!("humid state T={t} K p={p} Pa RH={rh}: {error}"));
            let replay = resolve_moist_air_state(
                loaded_dry.card(),
                &dry_at,
                MaterialPropertySelection::SingleClaimOnly,
                loaded_vapor.card(),
                &vapor_at,
                MaterialPropertySelection::SingleClaimOnly,
                rh,
            )
            .expect("humid state replay");
            assert_eq!(resolved, replay, "humid state replay moved");
            let (x, w, gamma, rho, c, mu, k, cp, prandtl, loss) = expected(t, p, rh);
            let state = resolved.state();
            relative_close(state.water_mole_fraction, x);
            relative_close(resolved.water_mass_fraction(), w);
            relative_close(state.gamma, gamma);
            relative_close(state.density, rho);
            relative_close(state.sound_speed, c);
            relative_close(state.dynamic_viscosity, mu);
            relative_close(state.thermal_conductivity, k);
            relative_close(state.specific_heat_cp, cp);
            relative_close(state.prandtl, prandtl);
            let actual_loss = fs_couple::air_path::oscillating_cylinder_air_resistance_per_length(
                RADIUS_M,
                OMEGA_RAD_S,
                state,
            )
            .unwrap();
            relative_close(actual_loss, loss);
            let absorption = OMEGA_RAD_S.powi(2) / (2.0 * rho * c.powi(3))
                * (4.0 * mu / 3.0 + (gamma - 1.0) * k / cp);
            relative_close(state.stokes_kirchhoff_absorption(OMEGA_RAD_S), absorption);
            assert_eq!(resolved.dry_air().parameters().properties().len(), 5);
            assert_eq!(resolved.water_vapor().parameters().properties().len(), 5);
            for (name, expected_value) in [
                ("sutherland_reference_viscosity", MU_REF_V),
                ("sutherland_reference_temperature", T_REF_V),
                ("sutherland_temperature", S_V),
            ] {
                relative_close(
                    resolved
                        .water_vapor()
                        .parameters()
                        .property(name)
                        .unwrap()
                        .value_si(),
                    expected_value,
                );
            }
            for (component, card) in [
                (resolved.dry_air(), loaded_dry.card()),
                (resolved.water_vapor(), loaded_vapor.card()),
            ] {
                for property in component.parameters().properties() {
                    card.claims()
                        .verify_receipt(&property.answer().receipt)
                        .expect("humid component receipt");
                }
            }
            println!(
                "humid-air T={t} K p={p} Pa RH={rh} x_w={x:.9e} w_w={w:.9e} rho={rho:.9e} kg/m3 c={c:.9e} m/s mu_actual={:.9e} mu_expected={mu:.9e} Pa*s k_actual={:.9e} k_expected={k:.9e} W/m/K Cp={cp:.9e} J/kg/K Pr={prandtl:.9e} loss_actual={actual_loss:.9e} loss_expected={loss:.9e} N*s/m2 absorption={absorption:.9e} 1/m relative_tolerance=1e-10 dry_claims={:?} vapor_claims={:?}",
                state.dynamic_viscosity,
                state.thermal_conductivity,
                resolved.dry_air().parameters().identity(),
                resolved.water_vapor().parameters().identity(),
            );
        }

        let dry_at = dry_point(293.15, 101_325.0);
        let vapor_at = vapor_point(293.15, 101_325.0);
        let dry_limit = resolve_moist_air_state(
            loaded_dry.card(),
            &dry_at,
            MaterialPropertySelection::SingleClaimOnly,
            loaded_vapor.card(),
            &vapor_at,
            MaterialPropertySelection::SingleClaimOnly,
            0.0,
        )
        .expect("RH=0 dry limit");
        let dry_state = resolve_sutherland_gas_state(
            loaded_dry.card(),
            &dry_at,
            ConductivityModel::Ussa1976AirFit,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("dry component");
        assert_eq!(
            dry_limit.state(),
            dry_state.state(),
            "RH=0 must be dry limit"
        );
        assert_eq!(dry_limit.relative_humidity(), 0.0);
        assert_eq!(dry_limit.water_mass_fraction(), 0.0);

        let wrong_flag = dry_at
            .clone()
            .with_quantity(
                "source-composition-ussa1976",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        assert!(
            resolve_moist_air_state(
                loaded_dry.card(),
                &wrong_flag,
                MaterialPropertySelection::SingleClaimOnly,
                loaded_vapor.card(),
                &vapor_at,
                MaterialPropertySelection::SingleClaimOnly,
                0.3,
            )
            .is_err()
        );
        assert!(
            resolve_moist_air_state(
                loaded_dry.card(),
                &dry_at,
                MaterialPropertySelection::SingleClaimOnly,
                loaded_vapor.card(),
                &vapor_at,
                MaterialPropertySelection::SingleClaimOnly,
                -0.01,
            )
            .is_err()
        );
        assert!(
            resolve_moist_air_state(
                loaded_dry.card(),
                &dry_at,
                MaterialPropertySelection::SingleClaimOnly,
                loaded_vapor.card(),
                &vapor_at,
                MaterialPropertySelection::SingleClaimOnly,
                1.01,
            )
            .is_err()
        );
        let outside_t = dry_point(250.0, 101_325.0);
        let outside_vapor_t = vapor_point(250.0, 101_325.0);
        assert!(
            resolve_moist_air_state(
                loaded_dry.card(),
                &outside_t,
                MaterialPropertySelection::SingleClaimOnly,
                loaded_vapor.card(),
                &outside_vapor_t,
                MaterialPropertySelection::SingleClaimOnly,
                0.3,
            )
            .is_err()
        );
        for mismatched in [
            vapor_point(294.15, 101_325.0),
            vapor_point(293.15, 100_000.0),
        ] {
            for rh in [0.0, 0.3] {
                assert!(
                    resolve_moist_air_state(
                        loaded_dry.card(),
                        &dry_at,
                        MaterialPropertySelection::SingleClaimOnly,
                        loaded_vapor.card(),
                        &mismatched,
                        MaterialPropertySelection::SingleClaimOnly,
                        rh,
                    )
                    .is_err(),
                    "mismatched components must refuse even at the dry limit"
                );
            }
        }
        for (request, path) in [
            ("examples/material-discovery/dry-air.json", dry_path),
            (
                "examples/material-discovery/water-vapor-component.json",
                vapor_path,
            ),
        ] {
            let discovery = fs_cli::run(vec![
                "--json".into(),
                "discover".into(),
                workspace_path(request).to_str().unwrap().into(),
                path.to_str().unwrap().into(),
            ]);
            assert_eq!(
                discovery.exit_code,
                fs_cli::exit::SUCCESS,
                "{}",
                discovery.stderr
            );
            assert!(
                discovery.stdout.contains("\"status\":\"complete\""),
                "{}",
                discovery.stdout
            );
        }
    }
    #[test]
    fn g1_g3_sourced_glycols_reach_heat_and_flow() {
        use fs_conduction::lumped::{BiotGate, LumpedNetwork, LumpedNode, solve_gated};
        use fs_lbm::Lbm;
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::liquid::resolve_liquid_state;
        use fs_material::state_point::MaterialPropertySelection;

        // Independent literal guide rows: rho[kg/m3], cp[J/kg/K], k[W/m/K], mu[Pa*s].
        // Versions are deliberately not pooled with later, differing Dow TDS values.
        let products = [
            (
                "dowtherm-sr1-eg50-40-60c",
                "ethylene-glycol-volume-fraction",
                "source-formulation-dowtherm-sr1",
                "examples/material-discovery/dowtherm-sr1-eg50.json",
                [
                    [1064.9, 3359.0, 0.378, 0.00226],
                    [1062.1, 3379.0, 0.381, 0.00200],
                    [1059.3, 3398.0, 0.383, 0.00178],
                    [1056.3, 3417.0, 0.385, 0.00159],
                    [1053.2, 3437.0, 0.387, 0.00143],
                ],
            ),
            (
                "dowfrost-pg50-40-60c",
                "propylene-glycol-volume-fraction",
                "source-formulation-dowfrost",
                "examples/material-discovery/dowfrost-pg50.json",
                [
                    [1032.1, 3609.0, 0.353, 0.00310],
                    [1028.8, 3628.0, 0.355, 0.00265],
                    [1025.4, 3648.0, 0.358, 0.00228],
                    [1021.9, 3667.0, 0.360, 0.00199],
                    [1018.2, 3686.0, 0.362, 0.00175],
                ],
            ),
        ];
        let relative_close = |actual: f64, expected: f64| {
            assert!(
                (actual / expected - 1.0).abs() < 1e-10,
                "actual={actual:.16e} expected={expected:.16e} relative_tolerance=1e-10"
            );
        };
        let mut warm_flows = Vec::new();
        for (slug, fraction_axis, formulation_axis, request, rows) in products {
            let (pack, path) = compile(slug);
            let original = NormalizedMaterialCardPack::new(
                MaterialStateId {
                    chemistry: slug.into(),
                    phase: "liquid".into(),
                    process: "versioned inhibited product; 50 volume percent glycol".into(),
                    revision: 0,
                },
                pack,
            )
            .unwrap();
            let database = fixture_dir().join(format!("{slug}.sqlite"));
            {
                let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
                store
                    .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                    .unwrap();
                store.seal_corpus().unwrap();
            }
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            let CatalogPack::MaterialCard(loaded) =
                store.load_catalog_pack(original.pack_id()).unwrap()
            else {
                panic!("wrong stored family")
            };
            assert_eq!(loaded, original);
            let card = loaded.card();
            let claim = card.claims().claims_for("density")[0].1;
            let mut states = rows
                .iter()
                .enumerate()
                .map(|(i, row)| (313.15 + i as f64 * 5.0, *row))
                .collect::<Vec<_>>();
            states.push((
                315.65,
                std::array::from_fn(|j| (rows[0][j] + rows[1][j]) / 2.0),
            ));
            for (t, [rho, cp, k, mu]) in states {
                let at = point(claim, &[("temperature", t)]);
                let state =
                    resolve_liquid_state(card, &at, MaterialPropertySelection::SingleClaimOnly)
                        .unwrap();
                assert_eq!(
                    state,
                    resolve_liquid_state(card, &at, MaterialPropertySelection::SingleClaimOnly)
                        .unwrap()
                );
                for (actual, expected) in [
                    (state.density_kg_m3(), rho),
                    (state.specific_heat_j_kg_k(), cp),
                    (state.thermal_conductivity_w_m_k(), k),
                    (state.dynamic_viscosity_pa_s(), mu),
                    (state.kinematic_viscosity_m2_s(), mu / rho),
                    (state.thermal_diffusivity_m2_s(), k / (rho * cp)),
                    (state.prandtl(), mu * cp / k),
                ] {
                    relative_close(actual, expected);
                }
                for name in [
                    "density",
                    "specific_heat_capacity",
                    "thermal_conductivity",
                    "dynamic_viscosity",
                ] {
                    let selected = state.material().property(name).unwrap();
                    // Receipt validation uses the real original claim set, not a copied expected hash.
                    card.claims()
                        .verify_receipt(&selected.answer().receipt)
                        .unwrap();
                }

                // Frozen-coefficient, small-signal heat response. The 0.1 K excursion
                // stays inside the sourced range; this is not continuously recoupled heat.
                let volume = 1e-6;
                let area = 1e-3;
                let path_length = 0.1;
                let ambient = if t == 333.15 { t - 0.1 } else { t + 0.1 };
                let capacity = state.density_kg_m3() * volume * state.specific_heat_j_kg_k();
                let conductance = state.thermal_conductivity_w_m_k() * area / path_length;
                let node = LumpedNode::new(
                    slug,
                    capacity,
                    conductance,
                    volume / area,
                    state.thermal_conductivity_w_m_k(),
                    area,
                )
                .unwrap();
                relative_close(node.biot(), 0.01);
                let network = LumpedNetwork::new(vec![node], ambient).unwrap();
                let time = 100.0;
                let heat =
                    solve_gated(&network, BiotGate::corpus_default(), &[0.0], &[t], time).unwrap();
                let tau = rho * volume * cp / (k * area / path_length);
                let expected_delta = (ambient - t) * (1.0 - (-time / tau).exp());
                let actual_delta = heat.temperature_k[0] - t;
                assert!((actual_delta - expected_delta).abs() < 1e-11);
                let energy = capacity * actual_delta;
                let supplied =
                    (k * area / path_length) * (ambient - t) * tau * (1.0 - (-time / tau).exp());
                assert!((energy - supplied).abs() < 1e-10);
                assert_eq!(
                    heat,
                    solve_gated(&network, BiotGate::corpus_default(), &[0.0], &[t], time).unwrap()
                );

                // Existing forced D2Q9 channel: fixed physical acceleration and geometry,
                // source viscosity controls tau. Map SI -> lattice and back explicitly.
                let dx: f64 = 1e-4;
                let dt: f64 = 1e-3;
                let acceleration: f64 = 1e-3;
                let height = 16.0 * dx;
                let nu_lattice = state.kinematic_viscosity_m2_s() * dt / (dx * dx);
                let tau_lattice = 0.5 + 3.0 * nu_lattice;
                assert!((0.8..1.5).contains(&tau_lattice));
                let mut flow = Lbm::channel(2, 16, tau_lattice, acceleration * dt * dt / dx);
                let initial_mass = flow.total_mass();
                let mut replay = flow.clone();
                flow.run(12_000);
                replay.run(12_000);
                assert_eq!(flow.x_velocity_profile(), replay.x_velocity_profile());
                let mass_defect = (flow.total_mass() - initial_mass).abs();
                assert!(mass_defect < 1e-8, "{slug} mass defect {mass_defect}");
                let peak_reference = acceleration * height * height / (8.0 * (mu / rho));
                let mut maximum_error = 0.0_f64;
                for (row, u_lattice) in flow.x_velocity_profile().into_iter().enumerate() {
                    let y = (row as f64 + 0.5) * dx;
                    let expected = acceleration * y * (height - y) / (2.0 * (mu / rho));
                    let actual = u_lattice * dx / dt;
                    maximum_error = maximum_error.max((actual - expected).abs() / peak_reference);
                }
                // The established BGK/bounce-back fixture has a spatial wall error.
                // This band covers that error; it is not an experimental accuracy bound.
                assert!(
                    maximum_error < 0.02,
                    "{slug} T={t} profile defect={maximum_error}"
                );
                let midpoint_velocity = flow.velocity(0, 7).0 * dx / dt;
                if t == 313.15 {
                    warm_flows.push(midpoint_velocity);
                }
                eprintln!(
                    "coolant={slug} T={t}K volume_fraction=0.5 pressure=unstated volume_reference_temperature=unstated rho={rho}kg/m3 cp={cp}J/kg/K k={k}W/m/K mu={mu}Pa*s source_relative_tolerance=1e-10 heat_delta_actual={actual_delta:.12e}K heat_delta_expected={expected_delta:.12e}K heat_absolute_tolerance=1e-11K energy={energy:.12e}J energy_defect={:.12e}J energy_tolerance=1e-10J flow_midpoint={midpoint_velocity:.12e}m/s profile_peak_relative_defect={maximum_error:.12e} profile_tolerance=0.02 lattice_mass_defect={mass_defect:.12e} source_bundle={:?}",
                    (energy - supplied).abs(),
                    state.material().identity()
                );
            }
            for (axis, wrong) in [
                (fraction_axis, 0.4),
                (fraction_axis, 0.6),
                (formulation_axis, 0.0),
                ("phase-liquid", 0.0),
                ("source-pressure-known", 1.0),
                ("volume-reference-temperature-known", 1.0),
                ("temperature", 253.15),
                ("temperature", 333.16),
            ] {
                let at = point(claim, &[(axis, wrong)]);
                assert!(
                    resolve_liquid_state(card, &at, MaterialPropertySelection::SingleClaimOnly)
                        .is_err(),
                    "{slug}: must refuse {axis}={wrong}"
                );
            }
            // A mass-fraction label is not a volume fraction. Missing the required
            // volume coordinate must refuse; no unknown-reference conversion is invented.
            let mut mass_basis = QueryPoint::new();
            for (axis, &(lo, _)) in claim.validity.bounds() {
                let name = if axis == fraction_axis {
                    "glycol-mass-fraction"
                } else {
                    axis
                };
                mass_basis = mass_basis
                    .with_quantity(name, claim.validity.axis_quantities()[axis], lo)
                    .unwrap();
            }
            assert!(
                resolve_liquid_state(
                    card,
                    &mass_basis,
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
            );
            let discovery = fs_cli::run(vec![
                "--json".into(),
                "discover".into(),
                workspace_path(request).to_str().unwrap().into(),
                path.to_str().unwrap().into(),
            ]);
            assert_eq!(
                discovery.exit_code,
                fs_cli::exit::SUCCESS,
                "{}",
                discovery.stderr
            );
            assert!(
                discovery.stdout.contains("\"status\":\"complete\""),
                "{}",
                discovery.stdout
            );
        }
        assert!(
            warm_flows[1] < warm_flows[0] * 0.8,
            "source viscosity difference must reach the actual channel"
        );
    }

    #[test]
    fn g1_g3_sourced_mineral_oil_reaches_heat_and_flow() {
        use fs_conduction::lumped::{BiotGate, LumpedNetwork, LumpedNode, solve_gated};
        use fs_lbm::Lbm;
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::liquid::{resolve_liquid_state, resolve_liquid_state_from_kinematic};
        use fs_material::state_point::MaterialPropertySelection::SingleClaimOnly;

        let slug = "shell-heat-transfer-oil-s2-2011";
        let (pack, path) = compile(slug);
        let original = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: slug.into(),
                phase: "liquid".into(),
                process: "May 2011 typical mineral heat-transfer oil; aged oil excluded".into(),
                revision: 0,
            },
            pack,
        )
        .unwrap();
        let database = fixture_dir().join("shell-s2.sqlite");
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(loaded) =
            store.load_catalog_pack(original.pack_id()).unwrap()
        else {
            panic!("wrong stored family")
        };
        assert_eq!(loaded, original);
        let card = loaded.card();
        let claim = card.claims().claims_for("density")[0].1;
        // Literal source T[K], rho[kg/m3], Cp[J/kg/K], k[W/m/K], nu[mm2/s].
        // The fifth state is only the declared linear interpolant, not a measured point.
        let rows = [
            [273.15, 876.0, 1809.0, 0.136, 223.0],
            [313.15, 850.0, 1954.0, 0.133, 25.0],
            [373.15, 811.0, 2173.0, 0.128, 4.7],
            [473.15, 746.0, 2538.0, 0.121, 1.1],
            [343.15, 830.5, 2063.5, 0.1305, 14.85],
        ];
        let close = |actual: f64, expected: f64| {
            assert!(
                (actual / expected - 1.0).abs() < 1e-10,
                "actual={actual:.16e} expected={expected:.16e} relative_tolerance=1e-10"
            );
        };
        let mut knot_flows = Vec::new();
        for [t, rho, cp, k, nu_mm2_s] in rows {
            let at = point(claim, &[("temperature", t)]);
            let state = resolve_liquid_state_from_kinematic(card, &at, SingleClaimOnly).unwrap();
            assert_eq!(
                state,
                resolve_liquid_state_from_kinematic(card, &at, SingleClaimOnly).unwrap()
            );
            // Explicit basis: this pack cannot satisfy the dynamic-source API.
            assert!(resolve_liquid_state(card, &at, SingleClaimOnly).is_err());
            let nu = nu_mm2_s * 1e-6;
            for (actual, expected) in [
                (state.density_kg_m3(), rho),
                (state.specific_heat_j_kg_k(), cp),
                (state.thermal_conductivity_w_m_k(), k),
                (state.kinematic_viscosity_m2_s(), nu),
                (state.dynamic_viscosity_pa_s(), rho * nu),
                (state.thermal_diffusivity_m2_s(), k / (rho * cp)),
                (state.prandtl(), rho * nu * cp / k),
            ] {
                close(actual, expected);
            }
            for name in [
                "density",
                "specific_heat_capacity",
                "thermal_conductivity",
                "kinematic_viscosity",
            ] {
                card.claims()
                    .verify_receipt(&state.material().property(name).unwrap().answer().receipt)
                    .unwrap();
            }
            assert!(state.material().property("dynamic_viscosity").is_none());

            // Existing gated lumped owner, frozen coefficients and 0.1 K excitation.
            let capacity = state.density_kg_m3() * 1e-6 * state.specific_heat_j_kg_k();
            let conductance = state.thermal_conductivity_w_m_k() * 0.01;
            let ambient = if t == 473.15 { t - 0.1 } else { t + 0.1 };
            let node = LumpedNode::new(slug, capacity, conductance, 0.001, k, 0.001).unwrap();
            close(node.biot(), 0.01);
            let network = LumpedNetwork::new(vec![node], ambient).unwrap();
            let heat =
                solve_gated(&network, BiotGate::corpus_default(), &[0.0], &[t], 100.0).unwrap();
            let decay = (-100.0 * k * 0.01 / (rho * 1e-6 * cp)).exp();
            let expected_delta = (ambient - t) * (1.0 - decay);
            let actual_delta = heat.temperature_k[0] - t;
            assert!((actual_delta - expected_delta).abs() < 1e-11);
            let energy_defect = capacity * (actual_delta - expected_delta);
            assert!(energy_defect.abs() < 1e-10);
            assert_eq!(
                heat,
                solve_gated(&network, BiotGate::corpus_default(), &[0.0], &[t], 100.0).unwrap()
            );

            // Hold physical geometry/acceleration fixed; choose dt for tau=1.
            // Adapting numerical time scale accommodates the 200-fold viscosity range.
            let dx: f64 = 1e-4;
            let dt = dx * dx / (6.0 * state.kinematic_viscosity_m2_s());
            let acceleration = 1e-3;
            let mut flow = Lbm::channel(2, 16, 1.0, acceleration * dt * dt / dx);
            let mass = flow.total_mass();
            let mut replay = flow.clone();
            flow.run(12_000);
            replay.run(12_000);
            assert_eq!(flow.x_velocity_profile(), replay.x_velocity_profile());
            let mass_defect = (flow.total_mass() - mass).abs();
            assert!(mass_defect < 1e-8);
            let height = 16.0 * dx;
            let peak = acceleration * height * height / (8.0 * nu);
            let mut defect = 0.0_f64;
            for (row, u) in flow.x_velocity_profile().into_iter().enumerate() {
                let y = (row as f64 + 0.5) * dx;
                let expected = acceleration * y * (height - y) / (2.0 * nu);
                defect = defect.max((u * dx / dt - expected).abs() / peak);
            }
            assert!(defect < 0.02, "T={t} profile defect={defect}");
            let midpoint = flow.velocity(0, 7).0 * dx / dt;
            if t != 343.15 {
                knot_flows.push(midpoint);
            }
            eprintln!(
                "oil={slug} T={t}K pressure=unstated aged_oil=excluded rho={rho}kg/m3 cp={cp}J/kg/K k={k}W/m/K source_nu={nu_mm2_s}mm2/s derived_mu={}Pa*s derived_Pr={} source_Pr_conflict=unresolved heat_delta={actual_delta:.12e}K expected_delta={expected_delta:.12e}K heat_tolerance=1e-11K energy_defect={energy_defect:.12e}J energy_tolerance=1e-10J flow_midpoint={midpoint:.12e}m/s profile_peak_relative_defect={defect:.12e} profile_tolerance=0.02 lattice_mass_defect={mass_defect:.12e} source_bundle={:?}",
                state.dynamic_viscosity_pa_s(),
                state.prandtl(),
                state.material().identity()
            );
        }
        assert!(knot_flows.windows(2).all(|pair| pair[1] > pair[0]));
        // Tiny cold-state lattice forcing loses precision to population subtraction.
        // Compare the physical inverse-viscosity ratio within the same profile band.
        assert!(((knot_flows[3] / knot_flows[0]) / (223.0 / 1.1) - 1.0).abs() < 0.02);
        for (axis, wrong) in [
            ("temperature", 273.14),
            ("temperature", 473.16),
            ("phase-liquid", 0.0),
            ("source-pressure-known", 1.0),
            ("source-formulation-shell-s2-2011", 0.0),
            ("aged-oil", 1.0),
        ] {
            assert!(
                resolve_liquid_state_from_kinematic(
                    card,
                    &point(claim, &[(axis, wrong)]),
                    SingleClaimOnly
                )
                .is_err(),
                "must refuse {axis}={wrong}"
            );
        }
        let discovery = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            workspace_path("examples/material-discovery/shell-heat-transfer-oil-s2.json")
                .to_str()
                .unwrap()
                .into(),
            path.to_str().unwrap().into(),
        ]);
        assert_eq!(
            discovery.exit_code,
            fs_cli::exit::SUCCESS,
            "{}",
            discovery.stderr
        );
        assert!(
            discovery.stdout.contains("\"status\":\"complete\""),
            "{}",
            discovery.stdout
        );
    }

    /// G1/G3: named PP conductivity reaches a reference-state series heat
    /// resistance without inventing heat capacity or a temperature law.
    #[test]
    fn g1_g3_source_card_slab_proteus_pp() {
        use fs_conduction::interface::{
            ResistanceOrigin, ResistanceUncertainty, SeriesThermalResistance, ThermalResistanceTerm,
        };
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, QueryPoint};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point,
        };

        let slug = "mcam-proteus-homopolymer-pp-natural-2023";
        let (pack, path) = compile(slug);
        assert_eq!(pack.claims().claim_count(), 6);
        let original = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "MCAM Proteus Homopolymer PP Natural".into(),
                phase: "dry solid stock shapes".into(),
                process: "2023 producer comparison data; processing schedule unstated".into(),
                revision: 0,
            },
            pack,
        )
        .unwrap();
        let database = fixture_dir().join("proteus-pp.sqlite");
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(loaded) = store.load_catalog_pack(slug).unwrap() else {
            panic!("wrong stored family");
        };
        assert_eq!(loaded, original);
        let card = loaded.card();
        for (name, expected) in [
            ("density", 910.0),
            ("thermal-conductivity", 0.22),
            ("tensile-modulus", 1.8e9),
            ("tensile-strength", 34.0e6),
            ("tensile-yield-strain", 0.06),
            ("mean-linear-thermal-expansion-coefficient", 150e-6),
        ] {
            let claims = card.claims().claims_for(name);
            assert_eq!(claims.len(), 1, "{name}");
            let claim = claims[0].1;
            let answer = card
                .claims()
                .query_typed(
                    &claim.key,
                    &point(claim, &[]),
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected);
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
        let claim = card.claims().claims_for("thermal-conductivity")[0].1;
        assert_eq!(claim.validity.bound("temperature"), Some((296.15, 296.15)));
        let modulus = card.claims().claims_for("tensile-modulus")[0].1;
        assert_eq!(
            modulus.validity.bound("temperature"),
            Some((296.15, 296.15))
        );
        assert_eq!(
            modulus.validity.bound("test-speed"),
            Some((1.0 / 60_000.0, 1.0 / 60_000.0))
        );
        assert!(
            card.claims()
                .query_typed(
                    &modulus.key,
                    &point(modulus, &[("test-speed", 2.0 / 60_000.0)]),
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
        let expansion = card
            .claims()
            .claims_for("mean-linear-thermal-expansion-coefficient")[0]
            .1;
        assert_eq!(
            expansion.validity.bound("expansion-interval-lower"),
            Some((296.15, 296.15))
        );
        assert_eq!(
            expansion.validity.bound("expansion-interval-upper"),
            Some((373.15, 373.15))
        );
        let at = point(claim, &[]);
        let requirement = ScalarPropertyRequirement::try_with_key(
            &claim.key,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap();
        let state = resolve_material_state_point(
            card,
            &at,
            &[requirement],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let slab = |name, length, at: &QueryPoint| {
            ThermalResistanceTerm::slab_from_card(
                name,
                length,
                0.01,
                card,
                &claim.key,
                at,
                SelectionPolicy::SingleClaimOnly,
            )
        };
        let a = slab("layer-a", 0.002, &at).unwrap();
        let b = slab("layer-b", 0.003, &at).unwrap();
        for (term, length) in [(&a, 0.002), (&b, 0.003)] {
            assert_eq!(term.uncertainty(), &ResistanceUncertainty::Unstated);
            let ResistanceOrigin::BulkMaterialCard {
                card_identity,
                receipt,
                length_m,
                area_m2,
                ..
            } = term.origin()
            else {
                panic!("lost bulk source provenance");
            };
            assert_eq!(*card_identity, card.content_hash());
            assert_eq!(*length_m, length);
            assert_eq!(*area_m2, 0.01);
            assert_eq!(
                receipt,
                &state
                    .property("thermal-conductivity")
                    .unwrap()
                    .answer()
                    .receipt
            );
            card.claims().verify_receipt(receipt).unwrap();
        }
        let network = SeriesThermalResistance::new(vec![a.clone(), b.clone()]).unwrap();
        let replay = SeriesThermalResistance::new(vec![b, a]).unwrap();
        assert_eq!(network, replay);
        // Independent Fourier-law conductance for 5 mm total thickness,
        // 100 cm² cross-section and the source conductivity at exactly 23 C.
        close(network.budget().value_k_per_w, 25.0 / 11.0);
        close(1.0 / network.budget().value_k_per_w, 0.44);
        assert_eq!(network.budget().complete_half_width_k_per_w(), None);
        assert_eq!(network.budget().unbounded_terms.len(), 2);
        for (axis, wrong) in [
            ("temperature", 296.14),
            ("temperature", 296.16),
            ("source-grade-proteus-pp-natural-2023", 0.0),
            ("dry-material", 0.0),
        ] {
            assert!(
                slab("invalid", 0.002, &point(claim, &[(axis, wrong)])).is_err(),
                "{axis}"
            );
        }
        assert!(slab("missing-context", 0.002, &QueryPoint::new()).is_err());
        assert!(
            card.claims()
                .claims_for("specific-heat-capacity")
                .is_empty()
        );
        assert!(
            card.claims()
                .claims_for("linear-thermal-expansion-coefficient")
                .is_empty()
        );
        let discovery = fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            workspace_path("examples/material-discovery/proteus-pp-conductivity.json")
                .to_str()
                .unwrap()
                .into(),
            path.to_str().unwrap().into(),
        ]);
        assert_eq!(
            discovery.exit_code,
            fs_cli::exit::SUCCESS,
            "{}",
            discovery.stderr
        );
        assert!(
            discovery.stdout.contains("\"status\":\"complete\""),
            "{}",
            discovery.stdout
        );
        eprintln!(
            "PP={slug} reference_temperature=296.15K k=0.22W/m/K lengths=0.002,0.003m area=0.01m2 resistance={}K/W expected=2.272727272727273K/W conductance=0.44W/K unknown_uncertainty_terms=2 replay=identical no_transient_or_temperature_range_claim=1",
            network.budget().value_k_per_w
        );
    }

    /// G1/G3: exact housing grades survive storage and supply room-state
    /// slab resistances; PC additionally retains its through-plane direction.
    /// Neither grade supplies a transient heating or isotropic elastic law.
    #[test]
    fn g1_g3_housing_polymers_source_card_slab() {
        use fs_conduction::interface::{ResistanceOrigin, ThermalResistanceTerm};
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, QueryPoint};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point,
        };

        for (slug, density, modulus, conductivity) in [
            ("covestro-makrolon-2405", 1200.0, 2.4e9, 0.2),
            ("ineos-terluran-gp-22", 1040.0, 2.3e9, 0.17),
        ] {
            let (pack, path) = compile(slug);
            assert_eq!(pack.claims().claim_count(), 3, "{slug}");
            let original = NormalizedMaterialCardPack::new(
                MaterialStateId {
                    chemistry: slug.into(),
                    phase: "solid polymer".into(),
                    process: "named injection-molding grade; producer typical data".into(),
                    revision: 0,
                },
                pack,
            )
            .unwrap();
            let database = fixture_dir().join(format!("{slug}.sqlite"));
            {
                let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
                store
                    .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                    .unwrap();
                store.seal_corpus().unwrap();
            }
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            let CatalogPack::MaterialCard(loaded) = store.load_catalog_pack(slug).unwrap() else {
                panic!("wrong stored family");
            };
            assert_eq!(loaded, original);
            let card = loaded.card();
            for (name, expected) in [
                ("density", density),
                ("tensile-modulus", modulus),
                ("thermal-conductivity", conductivity),
            ] {
                let claim = card.claims().claims_for(name)[0].1;
                let answer = card
                    .claims()
                    .query_typed(
                        &claim.key,
                        &point(claim, &[]),
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .unwrap();
                close(answer.evidence.value.value, expected);
                card.claims().verify_receipt(&answer.receipt).unwrap();
                assert!(
                    card.claims()
                        .query_typed(
                            &claim.key,
                            &QueryPoint::new(),
                            SelectionPolicy::SingleClaimOnly,
                        )
                        .is_err(),
                    "{slug}: {name} requires source context"
                );
            }
            for absent in [
                "specific-heat-capacity",
                "young-modulus",
                "linear-thermal-expansion-coefficient",
            ] {
                assert!(
                    card.claims().claims_for(absent).is_empty(),
                    "{slug}: {absent}"
                );
            }
            let claim = card.claims().claims_for("thermal-conductivity")[0].1;
            let is_pc = slug == "covestro-makrolon-2405";
            assert_eq!(claim.validity.bound("temperature"), Some((296.15, 296.15)));
            if is_pc {
                let discovery = fs_cli::run(vec![
                    "--json".into(),
                    "discover".into(),
                    workspace_path("examples/material-discovery/makrolon-2405-conductivity.json")
                        .to_str()
                        .unwrap()
                        .into(),
                    path.to_str().unwrap().into(),
                ]);
                assert_eq!(
                    discovery.exit_code,
                    fs_cli::exit::SUCCESS,
                    "{}",
                    discovery.stderr
                );
                assert!(
                    discovery.stdout.contains("\"status\":\"complete\""),
                    "{}",
                    discovery.stdout
                );
                assert_eq!(claim.validity.bound("relative-humidity"), Some((0.5, 0.5)));
                assert_eq!(claim.validity.bound("through-plane"), Some((1.0, 1.0)));
                let tensile = card.claims().claims_for("tensile-modulus")[0].1;
                assert_eq!(
                    tensile.validity.bound("test-speed"),
                    Some((1.0 / 60_000.0, 1.0 / 60_000.0))
                );
                assert!(
                    card.claims()
                        .query_typed(
                            &tensile.key,
                            &point(tensile, &[("test-speed", 2.0 / 60_000.0)]),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
            let at = point(claim, &[]);
            let state = resolve_material_state_point(
                card,
                &at,
                &[ScalarPropertyRequirement::try_with_key(
                    &claim.key,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap()],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let slab = |length, at: &QueryPoint| {
                ThermalResistanceTerm::slab_from_card(
                    "polymer-wall",
                    length,
                    0.01,
                    card,
                    &claim.key,
                    at,
                    SelectionPolicy::SingleClaimOnly,
                )
            };
            let wall = slab(0.002, &at).unwrap();
            let ResistanceOrigin::BulkMaterialCard {
                card_identity,
                receipt,
                length_m,
                area_m2,
                ..
            } = wall.origin()
            else {
                panic!("lost source receipt");
            };
            assert_eq!(*card_identity, card.content_hash());
            assert_eq!(*length_m, 0.002);
            assert_eq!(*area_m2, 0.01);
            assert_eq!(
                receipt,
                &state
                    .property("thermal-conductivity")
                    .unwrap()
                    .answer()
                    .receipt
            );
            card.claims().verify_receipt(receipt).unwrap();
            let network =
                fs_conduction::interface::SeriesThermalResistance::new(vec![wall]).unwrap();
            let expected = if is_pc { 1.0 } else { 20.0 / 17.0 };
            close(network.budget().value_k_per_w, expected);
            let thicker = fs_conduction::interface::SeriesThermalResistance::new(vec![
                slab(0.004, &at).unwrap(),
            ])
            .unwrap();
            close(thicker.budget().value_k_per_w, 2.0 * expected);
            let mut refusals = vec![("temperature", 296.16), ("temperature", 296.14)];
            if is_pc {
                refusals.extend([
                    ("relative-humidity", 0.4),
                    ("through-plane", 0.0),
                    ("source-grade-makrolon-2405", 0.0),
                ]);
            } else {
                refusals.extend([
                    ("source-grade-terluran-gp22", 0.0),
                    ("source-typical-uncolored", 0.0),
                ]);
            }
            for (axis, wrong) in refusals {
                assert!(
                    slab(0.002, &point(claim, &[(axis, wrong)])).is_err(),
                    "{axis}"
                );
            }
            eprintln!(
                "polymer={slug} T=296.15K k={conductivity}W/m/K wall=0.002m area=0.01m2 resistance={}K/W expected={expected}K/W thickness_scaling=2 no_transient_or_temperature_range_claim=1 store_replay=identical",
                network.budget().value_k_per_w
            );
        }
    }

    /// G0/G3: retain the named PE300 source observations and their unresolved
    /// test conditions through compiler, material-card storage and replay.
    /// This is not a thermal, constitutive, or product qualification model.
    #[test]
    fn g0_g3_sourced_ensinger_tecafine_pe300_natural_2017_observations() {
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_qty::QuantitySpec;

        let (pack, _) = compile("ensinger-tecafine-pe300-natural-2017");
        assert_eq!(pack.claims().claim_count(), 4);
        let original = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "PE-HD; Tecafine PE300 natural".into(),
                phase: "solid polymer".into(),
                process: "Ensinger 2017 source observations; test conditions unresolved".into(),
                revision: 0,
            },
            pack,
        )
        .unwrap();
        let database = fixture_dir().join("ensinger-tecafine-pe300-natural-2017.sqlite");
        {
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            store
                .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                .unwrap();
            store.seal_corpus().unwrap();
        }
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let CatalogPack::MaterialCard(loaded) =
            store.load_catalog_pack(original.pack_id()).unwrap()
        else {
            panic!("wrong stored family")
        };
        assert_eq!(loaded, original);
        let claims = loaded.card().claims();
        let axes = [
            ("source-grade-tecafine-pe300-natural-2017", 1.0),
            ("source-test-temperature-known", 0.0),
            ("source-processing-condition-known", 0.0),
            ("source-test-rate-known", 0.0),
        ];
        for (name, expected, dims) in [
            ("density", 960.0, Dims([-3, 1, 0, 0, 0, 0])),
            ("tensile-modulus", 1.1e9, Dims([-1, 1, -2, 0, 0, 0])),
            ("tensile-yield-strength", 23.0e6, Dims([-1, 1, -2, 0, 0, 0])),
            ("tensile-yield-strain", 0.09, Dims::NONE),
        ] {
            let available = claims.claims_for(name);
            assert_eq!(available.len(), 1, "one source observation for {name}");
            let (_, claim) = available[0];
            assert_eq!(claim.key.dims(), dims, "{name} dimensions");
            assert!(
                claim.validity.bound("temperature").is_none(),
                "{name} must not manufacture a source temperature"
            );
            for (axis, value) in axes {
                assert_eq!(claim.validity.bound(axis), Some((value, value)));
                assert_eq!(
                    claim.validity.axis_quantities().get(axis),
                    Some(&QuantitySpec::dimensional(Dims::NONE)),
                    "{name} {axis} must remain a typed dimensionless source condition"
                );
            }
            let answer = claims
                .query_typed(
                    &claim.key,
                    &point(claim, &[]),
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected);
            claims.verify_receipt(&answer.receipt).unwrap();
        }
        let density = claims.claims_for("density")[0].1;
        for (axis, wrong) in [
            ("source-grade-tecafine-pe300-natural-2017", 0.0),
            ("source-test-temperature-known", 1.0),
        ] {
            assert!(
                claims
                    .query_typed(
                        &density.key,
                        &point(density, &[(axis, wrong)]),
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .is_err(),
                "must refuse {axis}={wrong}"
            );
        }
        for absent in [
            "specific-heat-capacity",
            "thermal-conductivity",
            "young-modulus",
            "specific_heat_capacity",
            "thermal_conductivity",
            "young_modulus",
        ] {
            assert!(
                claims.claims_for(absent).is_empty(),
                "source-observed tensile data must not manufacture {absent}"
            );
        }
        println!(
            "HDPE Ensinger 2017 AA: 4 source facts compiled and stored/reopened; rho=960 kg/m3, tensile modulus=1100 MPa, yield=23 MPa, yield strain=0.09; wrong-grade and asserted-known-temperature queries refused; no thermal or isotropic elastic profile"
        );
    }

    /// G1/G3: frozen HDPE/PVC engineering references drive the same geometry
    /// and thermal boundary. The author-declared 25–26 C model range is not
    /// measured property coverage, a creep model, or product qualification.
    #[test]
    fn g1_g3_polymer_reference_comparison_reaches_heat() {
        use fs_conduction::lumped::{BiotGate, LumpedNetwork, LumpedNode, solve_gated};
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point,
        };
        use fs_qty::QuantitySpec;
        use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

        let mut deltas = Vec::new();
        for (slug, identity_flag, process, names, expected, discovery_file) in [
            (
                "iplex-pvc-u-pipe-engineering-reference",
                "reference-iplex-rigid-pvc-u-pipe",
                "extruded pipe-family frozen reference; schedule and source temperatures unknown",
                vec![
                    "density",
                    "tensile-modulus",
                    "poisson-ratio",
                    "specific-heat-capacity",
                    "thermal-conductivity",
                ],
                vec![1470.0, 3.2e9, 0.38, 1045.0, 0.138],
                "examples/material-discovery/iplex-pvc-reference.json",
            ),
            (
                "roechling-polystone-g-natural-reference",
                "reference-roechling-polystone-g-natural",
                "named HDPE stock material frozen reference; specimen history and source temperatures unknown",
                vec!["density", "specific-heat-capacity", "thermal-conductivity"],
                vec![950.0, 1900.0, 0.4],
                "examples/material-discovery/polystone-g-reference.json",
            ),
        ] {
            let (pack, path) = compile(slug);
            assert_eq!(pack.claims().claim_count(), names.len());
            let original = NormalizedMaterialCardPack::new(
                MaterialStateId {
                    chemistry: slug.into(),
                    phase: "solid polymer".into(),
                    process: process.into(),
                    revision: 0,
                },
                pack,
            )
            .unwrap();
            let database = fixture_dir().join(format!("{slug}.sqlite"));
            {
                let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
                store
                    .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                    .unwrap();
                store.seal_corpus().unwrap();
            }
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            let CatalogPack::MaterialCard(loaded) =
                store.load_catalog_pack(original.pack_id()).unwrap()
            else {
                panic!("wrong stored family")
            };
            assert_eq!(loaded, original);
            let card = loaded.card();
            let flags = [
                ("engineering-reference-model", 1.0),
                ("source-temperatures-matched", 0.0),
                (identity_flag, 1.0),
            ];
            let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            ));
            let dimensionless = QuantitySpec::dimensional(Dims::NONE);
            let keys: Vec<_> = names
                .iter()
                .map(|name| {
                    let found = card.claims().claims_for(name);
                    assert_eq!(found.len(), 1, "one source claim for {name}");
                    found[0].1.key.clone()
                })
                .collect();
            let requirements = keys
                .iter()
                .map(|key| {
                    ScalarPropertyRequirement::try_with_key(
                        key,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                })
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            let anchor = card.claims().claims_for("density")[0].1;
            let resolve = |point: &QueryPoint| {
                resolve_material_state_point(
                    card,
                    point,
                    &requirements,
                    MaterialPropertySelection::SingleClaimOnly,
                )
            };
            let lower = point(anchor, &[]);
            let state = resolve(&lower).unwrap();
            assert_eq!(state.properties().len(), names.len());
            for ((name, expected), key) in names
                .iter()
                .copied()
                .zip(expected.iter().copied())
                .zip(&keys)
            {
                close(state.property(name).unwrap().value_si(), expected);
                assert_eq!(state.property_by_key(key).unwrap().value_si(), expected);
            }
            for property in state.properties() {
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
            for (_, claim) in card.claims().claims_ordered() {
                assert_eq!(claim.validity.bound("temperature"), Some((298.15, 299.15)));
                assert_eq!(
                    claim.validity.axis_quantities().get("temperature"),
                    Some(&absolute_temperature)
                );
                for (axis, value) in flags {
                    assert_eq!(claim.validity.bound(axis), Some((value, value)));
                    assert_eq!(
                        claim.validity.axis_quantities().get(axis),
                        Some(&dimensionless)
                    );
                }
            }
            assert!(
                resolve(&QueryPoint::new()).is_err(),
                "opt-in axes are mandatory"
            );
            for (axis, wrong) in [
                ("engineering-reference-model", 0.0),
                ("source-temperatures-matched", 1.0),
                (identity_flag, 0.0),
                ("temperature", 298.14),
                ("temperature", 299.16),
            ] {
                assert!(
                    resolve(&point(anchor, &[(axis, wrong)])).is_err(),
                    "must refuse {axis}={wrong}"
                );
            }
            let upper = resolve(&point(anchor, &[("temperature", 299.15)])).unwrap();
            for (name, expected) in names.iter().copied().zip(expected.iter().copied()) {
                close(upper.property(name).unwrap().value_si(), expected);
            }
            assert!(
                card.claims()
                    .claims_for("linear-thermal-expansion-coefficient")
                    .is_empty(),
                "an observation-only CLTE cannot become an instantaneous alpha law"
            );
            let mut thermal_and_expansion = requirements.clone();
            thermal_and_expansion.push(
                ScalarPropertyRequirement::try_new(
                    "linear-thermal-expansion-coefficient",
                    Dims([0, 0, 0, -1, 0, 0]),
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap(),
            );
            assert!(
                resolve_material_state_point(
                    card,
                    &lower,
                    &thermal_and_expansion,
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .is_err(),
                "thermal inputs must not manufacture an unsupported expansion law"
            );

            let rho = state.property("density").unwrap().value_si();
            let cp = state.property("specific-heat-capacity").unwrap().value_si();
            let conductivity = state.property("thermal-conductivity").unwrap().value_si();
            let volume = 1.0e-6;
            let area = 1.0e-3;
            let conductance = conductivity * area / 0.1;
            let capacity = rho * volume * cp;
            let initial = 298.15;
            let ambient = 299.15;
            let time = 100.0;
            let node = LumpedNode::new(
                slug,
                capacity,
                conductance,
                volume / area,
                conductivity,
                area,
            )
            .unwrap();
            let network = LumpedNetwork::new(vec![node], ambient).unwrap();
            let heat = solve_gated(
                &network,
                BiotGate::corpus_default(),
                &[0.0],
                &[initial],
                time,
            )
            .unwrap();
            assert_eq!(
                heat,
                solve_gated(
                    &network,
                    BiotGate::corpus_default(),
                    &[0.0],
                    &[initial],
                    time,
                )
                .unwrap()
            );
            let tau = capacity / conductance;
            let expected_delta = (ambient - initial) * (1.0 - (-time / tau).exp());
            let actual_delta = heat.temperature_k[0] - initial;
            deltas.push(actual_delta);
            assert!((actual_delta - expected_delta).abs() < 1e-11);
            let energy = capacity * actual_delta;
            let supplied = conductance * (ambient - initial) * tau * (1.0 - (-time / tau).exp());
            assert!((energy - supplied).abs() < 1e-10);

            let discovery = fs_cli::run(vec![
                "--json".into(),
                "discover".into(),
                workspace_path(discovery_file).to_str().unwrap().into(),
                path.to_str().unwrap().into(),
            ]);
            assert_eq!(
                discovery.exit_code,
                fs_cli::exit::SUCCESS,
                "{}",
                discovery.stderr
            );
            assert!(
                discovery.stdout.contains("\"status\":\"complete\""),
                "{}",
                discovery.stdout
            );
            eprintln!(
                "polymer={slug} comparison=fixed_geometry volume=1e-6m3 area=1e-3m2 path=0.1m time=100s author_frozen_reference_temperature=298.15..299.15K source_temperatures_matched=0 rho={rho}kg/m3 Cp={cp}J/kg/K k={conductivity}W/m/K delta={actual_delta:.12e}K tau={tau:.12e}s energy={energy:.12e}J no_measurement_coverage_or_product_qualification=1 source_bundle={:?}",
                state.identity()
            );
        }
        assert!(
            deltas[1] > deltas[0],
            "source-driven HDPE reference warms faster than PVC under identical geometry and boundary"
        );
    }
}

static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

const PACK_BYTES_GOLDEN: usize = 3_177;
const PACK_HASH_GOLDEN: &str = "c1fb2f443708d297423179f4ac6024ee26b1d0c940a229d1d9084726ccbd2bc5";
const NASA9_PACK_BYTES_GOLDEN: usize = 4_940;
const NASA9_PACK_HASH_GOLDEN: &str =
    "006177a7cc6f7b4ae10a9eb4a5bf49faaf21911ef9473190a29ecfc3a818a162";
const MATERIAL_COMPILER_ID: &str = "frankensim-matdb-pack-compiler-v1";
const NASA9_COMPILER_ID: &str = "frankensim-matdb-nasa9-model-pack-compiler-v1";
const KINETICS_COMPILER_ID: &str = "frankensim-matdb-kinetics-model-pack-compiler-v1";
const SPECIES_COMPILER_ID: &str = "frankensim-matdb-species-pack-compiler-v1";
const METHANE_SEED_MANIFEST: &str = "data/matdb/seed-v1/methane/manifest.tsv";
const ALUMINUM_6061_T6_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aluminum-6061-t6-cryogenic/manifest.tsv";
const OFHC_COPPER_SEED_MANIFEST: &str = "data/matdb/seed-v1/ofhc-copper-rrr100/manifest.tsv";
const PTFE_TEFLON_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/ptfe-teflon-nist-cryogenic/manifest.tsv";
const PEEK_THERMIC_SEED_MANIFEST: &str = "data/matdb/seed-v1/peek-nasa-thermic-plate/manifest.tsv";
const NASA_CR_115153_WATER_ETHYLENE_GLYCOL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-cr-115153-water-ethylene-glycol/manifest.tsv";
const N0602_001_NITRILE_JP8_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/n0602-001-nitrile-jp8-compatibility/manifest.tsv";
const NASA_TN_D_8184_M19_MATERIAL_DECK_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-tn-d-8184-m19-material-deck/manifest.tsv";
const NASA_CR_4538_TEMPEL_24N208_M19_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-cr-4538-tempel-24n208-m19/manifest.tsv";
const TORRENT_2018_M19_STEINMETZ_INPUTS_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/torrent-2018-m19-steinmetz-inputs/manifest.tsv";
const NGYC_N42_SINTERED_NICKEL_COATED_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/ngyc-n42-sintered-nickel-coated/manifest.tsv";
const JINSHAN_N42_PRISTINE_TEMPERATURE_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/jinshan-n42-pristine-temperature/manifest.tsv";
const SJOLUND_2020_Y30_CATALOG_MODEL_INPUTS_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/sjolund-2020-y30-catalog-model-inputs/manifest.tsv";
const KIM_BAEK_2026_Y30_AFCP_DEMAGNETIZATION_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/kim-baek-2026-y30-afcp-demagnetization/manifest.tsv";
const NACA_TN_2680_ISOOCTANE_FLAME_SPEED_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/naca-tn-2680-isooctane-flame-speed/manifest.tsv";
const FACE_G_CDTRF_G_2023_V1_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/face-g-cdtrf-g-2023-v1/manifest.tsv";
const NIST_SRM_1720_NORTHERN_CONTINENTAL_AIR_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nist-srm-1720-northern-continental-air/manifest.tsv";
const NIST_SRM_2728_AUTO_EMISSION_REFERENCE_GAS_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nist-srm-2728-auto-emission-reference-gas/manifest.tsv";
const WO2018_125520_FORMULATION_8_5W30_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/wo2018-125520-formulation-8-5w30/manifest.tsv";
const NASA_UAM_MW16C_POLYIMIDE_WIRE_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-uam-mw16c-polyimide-magnet-wire/manifest.tsv";
const NASA_UAM_NOMEX_410_SLOT_LINER_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-uam-nomex-410-slot-liner/manifest.tsv";
const NASA_UAM_COOLTHERM_EP2000_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-uam-cooltherm-ep2000-180c-cure/manifest.tsv";
const AISI_4140_RC33_SEED_MANIFEST: &str = "data/matdb/seed-v1/aisi-4140-rc33/manifest.tsv";
const AISI_1045_COLD_DRAWN_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-1045-cold-drawn/manifest.tsv";
const AISI_52100_CVM_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-52100-cvm-hot-hardness/manifest.tsv";
const AISI_9310_CVM_CARBURIZED_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-9310-cvm-carburized/manifest.tsv";
const NAPC_PE_5_L_1274_GEAR_OIL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/napc-pe-5-l-1274-gear-oil/manifest.tsv";
const NAPC_PE_5_L_1307_1553_GEAR_OIL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/napc-pe-5-l-1307-1553-gear-oil/manifest.tsv";
const RHEOLUBE_2000_PENNZANE_GREASE_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/rheolube-2000-pennzane-grease/manifest.tsv";
const PENNZANE_SHF_X_2000_BEARING_OIL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/pennzane-shf-x-2000-bearing-oil/manifest.tsv";
const GRAY_CAST_IRON_S2_S_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/gray-cast-iron-s2-s/manifest.tsv";
const NASA_CR_195445_OMC_PS200_ROTARY_COATING_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-cr-195445-omc-ps200-rotary-coating/manifest.tsv";
const NASA_SEED_LICENSE: &str = "Work-of-the-US-Government-Public-Use-Permitted";
const PUBLIC_USE_PERMITTED_LICENSE: &str = "Public-Use-Permitted";
const PUBLIC_USE_AND_PATENT_PUBLICATION_LICENSE: &str =
    "Public-Use-Permitted-and-US-Patent-Publication";
const CC_BY_4_0_LICENSE: &str = "CC-BY-4.0";
const NIST_PUBLIC_INFORMATION_LICENSE: &str = "NIST-Public-Information-Attribution-Requested";
const USPTO_PATENT_TEXT_LICENSE: &str = "USPTO-Patent-Text-Typically-No-Copyright-Restrictions";
const NASA_METHANE_MOLAR_MASS_G_PER_MOL: f64 = 16.042_46;
const NIST_SRD69_METHANE_MOLAR_MASS_KG_PER_MOL: f64 = 0.016_042_5;
const NIST_SRD69_DISPLAY_ROUNDING_HALF_WIDTH_KG_PER_MOL: f64 = 0.000_000_05;

#[derive(Clone, Copy)]
struct CommittedSpeciesSeed {
    manifest: &'static str,
    species: &'static str,
    nasa_molar_mass_g_per_mol: f64,
    nist_molar_mass_g_per_mol: f64,
    nist_display_rounding_half_width_g_per_mol: f64,
}

const AIR_EXHAUST_SPECIES_SEEDS: [CommittedSpeciesSeed; 6] = [
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/nitrogen/manifest.tsv",
        species: "N2",
        nasa_molar_mass_g_per_mol: 28.013_40,
        nist_molar_mass_g_per_mol: 28.013_4,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/oxygen/manifest.tsv",
        species: "O2",
        nasa_molar_mass_g_per_mol: 31.998_80,
        nist_molar_mass_g_per_mol: 31.998_8,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/argon/manifest.tsv",
        species: "Ar",
        nasa_molar_mass_g_per_mol: 39.948_00,
        nist_molar_mass_g_per_mol: 39.948,
        nist_display_rounding_half_width_g_per_mol: 0.000_5,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/carbon-dioxide/manifest.tsv",
        species: "CO2",
        nasa_molar_mass_g_per_mol: 44.009_50,
        nist_molar_mass_g_per_mol: 44.009_5,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/water-vapor/manifest.tsv",
        species: "H2O",
        nasa_molar_mass_g_per_mol: 18.015_28,
        nist_molar_mass_g_per_mol: 18.015_3,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/carbon-monoxide/manifest.tsv",
        species: "CO",
        nasa_molar_mass_g_per_mol: 28.010_10,
        nist_molar_mass_g_per_mol: 28.010_1,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
];

const MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tfixture-alloy-x\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture handbook table 7\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tsource.tsv\tmaterial-tsv-v1\n",
);

const SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
    "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant\n",
    "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t1\n",
    "validity\tdensity\ttemperature\t0\t100\tdegC\n",
    "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tlinear\n",
    "uncertainty\tmodulus\trelative\t2\t%\t0.95\t1\n",
    "validity\tmodulus\ttemperature\t0\t100\tdegC\n",
    "frame\tmodulus\tspecimen\tlab\n",
    "joint\tcoupon\tdensity-modulus\tdensity:scalar,modulus:y:0\t0.000025,0,0.000009\t1,0,1\t1\n",
);

const MATERIAL_FAMILIES_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tfixture-material-families\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture material-family extracts\n",
    "license\tCC-BY-4.0\n",
    "source\thandbook\thandbook.tsv\tmaterial-tsv-v1\n",
    "source\tbh-curve\tbh.tsv\tmaterial-tsv-v1\n",
    "source\tsn-curve\tsn.tsv\tmaterial-tsv-v1\n",
    "source\tlubricant\tlubricant.tsv\tmaterial-tsv-v1\n",
);

const HANDBOOK_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\thandbook-coupon\talloy-X-solution-treated\thandbook-table\tdensity extract\n",
    "scalar\thandbook-density\thandbook-coupon\tdensity\t7.85\tg/cm3\tconstant\n",
    "uncertainty\thandbook-density\trelative\t0.5\t%\t0.95\t1\n",
    "validity\thandbook-density\ttemperature\t0\t100\tdegC\n",
);

const BH_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tbh-loop\talloy-X-ring\tquasistatic-hysteresis\tdemagnetized branch\n",
    "curve\tbh-curve\tbh-loop\tmagnetic_flux_density\tmagnetic_field_strength\tA/m\tT\t0:0,100:0.2,1000:1.5\tlinear\n",
    "uncertainty\tbh-curve\trelative\t1\t%\t0.95\t1\n",
    "validity\tbh-curve\tmagnetic_field_strength\t0\t1000\tA/m\n",
    "validity\tbh-curve\ttemperature\t20\t25\tdegC\n",
);

const SN_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tsn-coupons\talloy-X-polished\tconstant-amplitude-fatigue\tfully reversed\n",
    "curve\tsn-curve\tsn-coupons\tfatigue_life\tstress_amplitude\tMPa\t1\t100:10000000,250:500000,400:20000\ttabulated\n",
    "uncertainty\tsn-curve\trelative\t5\t%\t0.90\t1\n",
    "validity\tsn-curve\tstress_amplitude\t100\t400\tMPa\n",
);

const LUBRICANT_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tlubricant-batch\tPAO-4-batch-A\trotational-rheometer\tnew fluid\n",
    "curve\tlubricant-viscosity\tlubricant-batch\tdynamic_viscosity\ttemperature\tdegC\tPa*s\t-20:0.12,40:0.018,100:0.005\tlinear\n",
    "uncertainty\tlubricant-viscosity\trelative\t3\t%\t0.95\t1\n",
    "validity\tlubricant-viscosity\ttemperature\t-20\t100\tdegC\n",
);

const NASA9_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tN2\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture NASA-9 species table\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tnasa9.tsv\tnasa9-v1\n",
);

const NASA9_SOURCE: &str = concat!(
    "frankensim.nasa9-source.v1\n",
    "region\tN2\tlow\t-73.15\t700\tdegC\t100\tkPa\n",
    "coefficient\tN2\tlow\ta0\t0\tK^2\n",
    "coefficient\tN2\tlow\ta1\t0\tK\n",
    "coefficient\tN2\tlow\ta2\t3.5\t1\n",
    "coefficient\tN2\tlow\ta3\t0.001\tK^-1\n",
    "coefficient\tN2\tlow\ta4\t0\tK^-2\n",
    "coefficient\tN2\tlow\ta5\t0\tK^-3\n",
    "coefficient\tN2\tlow\ta6\t0\tK^-4\n",
    "coefficient\tN2\tlow\ta7\t100\tK\n",
    "coefficient\tN2\tlow\ta8\t1\t1\n",
    "region\tN2\thigh\t1000\t6000\tK\t100000\tPa\n",
    "coefficient\tN2\thigh\ta0\t0\tK^2\n",
    "coefficient\tN2\thigh\ta1\t0\tK\n",
    "coefficient\tN2\thigh\ta2\t4\t1\n",
    "coefficient\tN2\thigh\ta3\t0.0001\tK^-1\n",
    "coefficient\tN2\thigh\ta4\t0\tK^-2\n",
    "coefficient\tN2\thigh\ta5\t0\tK^-3\n",
    "coefficient\tN2\thigh\ta6\t0\tK^-4\n",
    "coefficient\tN2\thigh\ta7\t200\tK\n",
    "coefficient\tN2\thigh\ta8\t2\t1\n",
);

const KINETICS_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\twater-formation\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture first-order kinetics table\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tkinetics.tsv\tkinetics-v1\n",
);

const KINETICS_SOURCE: &str = concat!(
    "frankensim.kinetics-source.v1\n",
    "reaction\twater-formation\tfirst-order\t300\t2500\tK\n",
    "parameter\twater-formation\tactivation_temperature\t12000\tK\n",
    "parameter\twater-formation\tpre_exponential\t2.5e7\ts^-1\n",
);

const SPECIES_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tN2\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture licensed species metadata\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tspecies.tsv\tspecies-v1\n",
);

const SPECIES_SOURCE: &str = concat!(
    "frankensim.species-source.v1\n",
    "species\tN2\t28.0134\tg/mol\tgas\tideal-gas\t100\tkPa\tNASA-TP-2002-211556\n",
);

fn fixture_dir() -> PathBuf {
    loop {
        let sequence = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "frankensim-matdb-pack-cli-test-{}-{sequence}",
            std::process::id()
        ));
        match fs::create_dir(&path) {
            Ok(()) => return path,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => panic!("unique fixture directory: {error}"),
        }
    }
}

mod glass_reference {
    use super::common_material_acquisition::{close, compile, point};
    use super::*;

    /// G1/G3: source-reference glass bundles drive a frozen-coefficient heat
    /// operator. The source temperatures differ, so this is not a common-T
    /// material comparison or a temperature-qualified property law.
    #[test]
    fn g1_g3_glass_reference_profiles_reach_heat() {
        use fs_conduction::lumped::{BiotGate, LumpedNetwork, LumpedNode, solve_gated};
        use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, QueryPoint};
        use fs_matdb_store::{CatalogPack, MaterialStore};
        use fs_material::state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point,
        };
        use fs_qty::QuantitySpec;
        use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

        let references = [
            (
                "pilkington-float-glass-reference",
                "Pilkington float glass reference",
                [2500.0, 72.0e9, 0.22, 880.0, 0.937, 8.3e-6],
                (297.038_888_888_888_9, 574.816_666_666_666_7),
            ),
            (
                "schott-borofloat33-reference",
                "SCHOTT BOROFLOAT 33 reference",
                [2230.0, 64.0e9, 0.2, 830.0, 1.2, 3.25e-6],
                (293.15, 573.15),
            ),
        ];
        let names = [
            "density",
            "young-modulus",
            "poisson-ratio",
            "specific-heat-capacity",
            "thermal-conductivity",
            "mean-linear-thermal-expansion-coefficient",
        ];
        let flags = [
            ("engineering-reference-model", 1.0),
            ("source-temperatures-matched", 0.0),
            ("reference-annealed-glass", 1.0),
        ];
        let dimensionless = QuantitySpec::dimensional(Dims::NONE);
        let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let mut temperature_deltas = Vec::new();
        for (slug, chemistry, expected, expansion_interval) in references {
            let (pack, path) = compile(slug);
            assert_eq!(pack.claims().claim_count(), 6);
            let original = NormalizedMaterialCardPack::new(
                MaterialStateId {
                    chemistry: chemistry.into(),
                    phase: "annealed glass reference".into(),
                    process: "cross-source engineering reference; source temperatures unmatched"
                        .into(),
                    revision: 0,
                },
                pack,
            )
            .unwrap();
            let database = fixture_dir().join(format!("{slug}.sqlite"));
            {
                let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
                store
                    .ingest_bundle(&[CatalogPack::MaterialCard(original.clone())])
                    .unwrap();
                store.seal_corpus().unwrap();
            }
            let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
            let CatalogPack::MaterialCard(loaded) =
                store.load_catalog_pack(original.pack_id()).unwrap()
            else {
                panic!("wrong stored family")
            };
            assert_eq!(loaded, original);
            let claims = loaded.card().claims();
            let keys = names.map(|name| {
                let available = claims.claims_for(name);
                assert_eq!(available.len(), 1, "one source claim for {slug} {name}");
                available[0].1.key.clone()
            });
            let requirements = keys
                .iter()
                .map(|key| {
                    ScalarPropertyRequirement::try_with_key(
                        key,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                })
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            let mean = claims.claims_for("mean-linear-thermal-expansion-coefficient")[0].1;
            let at = point(mean, &[]);
            let resolve = |at: &QueryPoint| {
                resolve_material_state_point(
                    loaded.card(),
                    at,
                    &requirements,
                    MaterialPropertySelection::SingleClaimOnly,
                )
            };
            let state = resolve(&at).unwrap();
            assert_eq!(state.properties().len(), 6);
            let without_opt_in = QueryPoint::new()
                .with_quantity("temperature", absolute_temperature, 298.15)
                .unwrap();
            assert!(resolve(&without_opt_in).is_err());
            for ((name, expected), key) in names.into_iter().zip(expected).zip(&keys) {
                close(state.property(name).unwrap().value_si(), expected);
                assert_eq!(state.property_by_key(key).unwrap().value_si(), expected);
            }
            for property in state.properties() {
                claims.verify_receipt(&property.answer().receipt).unwrap();
            }
            for (_, claim) in claims.claims_ordered() {
                assert!(
                    claim.validity.bound("temperature") == Some((298.15, 299.15)),
                    "{slug}: author-selected model range, not measured property coverage"
                );
                assert_eq!(
                    claim.validity.axis_quantities().get("temperature"),
                    Some(&absolute_temperature)
                );
                for (axis, value) in flags {
                    assert_eq!(claim.validity.bound(axis), Some((value, value)));
                    assert_eq!(
                        claim.validity.axis_quantities().get(axis),
                        Some(&dimensionless)
                    );
                }
            }
            assert_eq!(
                mean.validity.bound("expansion-interval-lower"),
                Some((expansion_interval.0, expansion_interval.0))
            );
            assert_eq!(
                mean.validity.bound("expansion-interval-upper"),
                Some((expansion_interval.1, expansion_interval.1))
            );
            for axis in ["expansion-interval-lower", "expansion-interval-upper"] {
                assert_eq!(
                    mean.validity.axis_quantities().get(axis),
                    Some(&absolute_temperature)
                );
            }
            assert!(
                claims
                    .claims_for("linear-thermal-expansion-coefficient")
                    .is_empty(),
                "mean CTE must not manufacture an instantaneous law"
            );
            for (axis, wrong) in [
                ("engineering-reference-model", 0.0),
                ("source-temperatures-matched", 1.0),
                ("reference-annealed-glass", 0.0),
                ("temperature", 297.0),
                ("temperature", 300.0),
                ("expansion-interval-lower", expansion_interval.0 + 1.0),
            ] {
                assert!(resolve(&point(mean, &[(axis, wrong)])).is_err());
            }
            let upper = resolve(&point(mean, &[("temperature", 299.15)])).unwrap();
            for (name, value) in names.into_iter().zip(expected) {
                close(upper.property(name).unwrap().value_si(), value);
            }

            let rho = state.property("density").unwrap().value_si();
            let cp = state.property("specific-heat-capacity").unwrap().value_si();
            let conductivity = state.property("thermal-conductivity").unwrap().value_si();
            let volume = 1.0e-6;
            let area = 1.0e-3;
            let conductance = conductivity * area / 0.1;
            let capacity = rho * volume * cp;
            let node = LumpedNode::new(
                slug,
                capacity,
                conductance,
                volume / area,
                conductivity,
                area,
            )
            .unwrap();
            let ambient = 299.15;
            let initial = 298.15;
            let time = 100.0;
            let network = LumpedNetwork::new(vec![node], ambient).unwrap();
            let heat = solve_gated(
                &network,
                BiotGate::corpus_default(),
                &[0.0],
                &[initial],
                time,
            )
            .unwrap();
            let replay = solve_gated(
                &network,
                BiotGate::corpus_default(),
                &[0.0],
                &[initial],
                time,
            )
            .unwrap();
            assert_eq!(heat, replay);
            let tau = capacity / conductance;
            let expected_delta = (ambient - initial) * (1.0 - (-time / tau).exp());
            let actual_delta = heat.temperature_k[0] - initial;
            assert!((actual_delta - expected_delta).abs() < 1.0e-11);
            let energy = capacity * actual_delta;
            let supplied = conductance * (ambient - initial) * tau * (1.0 - (-time / tau).exp());
            assert!((energy - supplied).abs() < 1.0e-10);
            temperature_deltas.push(actual_delta);
            eprintln!(
                "glass={slug} frozen_nominal_engineering_reference=1 source_temperatures_matched=0 rho={rho}kg/m3 Cp={cp}J/kg/K k={conductivity}W/m/K Tinitial={initial}K Tambient={ambient}K time={time}s delta={actual_delta:.12e}K tau={tau:.12e}s energy={energy:.12e}J source_bundle={:?} no_physical_accuracy_bound=1",
                state.identity()
            );

            let discovery = fs_cli::run(vec![
                "--json".into(),
                "discover".into(),
                workspace_path("examples/material-discovery/glass-reference.json")
                    .to_str()
                    .unwrap()
                    .into(),
                path.to_str().unwrap().into(),
            ]);
            assert_eq!(
                discovery.exit_code,
                fs_cli::exit::SUCCESS,
                "{}",
                discovery.stderr
            );
            assert!(
                discovery.stdout.contains("\"status\":\"complete\""),
                "{}",
                discovery.stdout
            );
        }
        assert!(
            temperature_deltas[1] > temperature_deltas[0],
            "fixed geometry must retain the two sourced frozen-coefficient heat responses"
        );
    }
}

fn write_fixture(source: &str) -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, MANIFEST).expect("write manifest fixture");
    fs::write(directory.join("source.tsv"), source).expect("write source fixture");
    (directory, manifest)
}

fn write_material_families_fixture() -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, MATERIAL_FAMILIES_MANIFEST).expect("write family manifest fixture");
    for (name, source) in [
        ("handbook.tsv", HANDBOOK_SOURCE),
        ("bh.tsv", BH_SOURCE),
        ("sn.tsv", SN_SOURCE),
        ("lubricant.tsv", LUBRICANT_SOURCE),
    ] {
        fs::write(directory.join(name), source).expect("write family source fixture");
    }
    (directory, manifest)
}

fn write_nasa9_fixture() -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, NASA9_MANIFEST).expect("write NASA-9 manifest fixture");
    fs::write(directory.join("nasa9.tsv"), NASA9_SOURCE).expect("write NASA-9 source fixture");
    (directory, manifest)
}

fn write_kinetics_fixture() -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, KINETICS_MANIFEST).expect("write kinetics manifest fixture");
    fs::write(directory.join("kinetics.tsv"), KINETICS_SOURCE)
        .expect("write kinetics source fixture");
    (directory, manifest)
}

fn write_species_fixture(source: &str) -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, SPECIES_MANIFEST).expect("write species manifest fixture");
    fs::write(directory.join("species.tsv"), source).expect("write species source fixture");
    (directory, manifest)
}

fn run_compiler(manifest: &Path, output: &Path) -> Output {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a workspace parent");
    Command::new(env!("CARGO_BIN_EXE_xtask"))
        .arg("matdb-pack")
        .arg("--manifest")
        .arg(manifest)
        .arg("--out")
        .arg(output)
        .env("CARGO_WORKSPACE_DIR", workspace)
        .output()
        .expect("run xtask matdb-pack")
}

fn workspace_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a workspace parent")
        .join(relative)
}

fn assert_decision_compiler(output: &Output, expected: &str) {
    let stdout = std::str::from_utf8(&output.stdout).expect("decision stream is UTF-8");
    assert!(!stdout.is_empty(), "compiler emitted no decision rows");
    let expected_prefix = format!("{{\"check\":\"matdb-pack\",\"compiler\":\"{expected}\",");
    assert!(
        stdout.lines().all(|row| row.starts_with(&expected_prefix)),
        "decision row used the wrong compiler identity:\n{stdout}"
    );
}

#[test]
fn g3_cli_compiles_two_identical_pinned_packs() {
    let (directory, manifest) = write_fixture(SOURCE);
    let first_path = directory.join("first.fsmatpk");
    let second_path = directory.join("second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first compiler run failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second compiler run failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "decision stream moved");
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(&first_path).expect("read first normalized pack");
    let second_bytes = fs::read(&second_path).expect("read second normalized pack");
    assert_eq!(first_bytes, second_bytes, "published pack bytes moved");
    assert_eq!(first_bytes.len(), PACK_BYTES_GOLDEN);
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode compiler output");
    assert_eq!(decoded.content_hash().to_string(), PACK_HASH_GOLDEN);
    NormalizedPack::from_bytes_verified(decoded.content_hash(), &first_bytes)
        .expect("externally pinned bytes re-admit");

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    let rows: Vec<_> = decisions.lines().collect();
    assert!(!rows.is_empty());
    assert!(
        rows.iter()
            .all(|row| row.starts_with("{\"check\":\"matdb-pack\""))
    );
    assert!(
        rows.iter()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{PACK_HASH_GOLDEN}\"")))
    );
    assert!(decisions.contains("\"reason_code\":\"published_new_verified_artifact\""));
    assert!(decisions.contains("\"reason_code\":\"joint_statistics_normalized\""));
}

#[test]
fn g3_cli_compiles_handbook_bh_sn_and_lubricant_material_claims() {
    let (directory, manifest) = write_material_families_fixture();
    let first_path = directory.join("families-first.fsmatpk");
    let second_path = directory.join("families-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first material-family compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second material-family compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "decision stream moved");

    let first_bytes = fs::read(first_path).expect("read first family pack");
    let second_bytes = fs::read(second_path).expect("read second family pack");
    assert_eq!(first_bytes, second_bytes, "material-family pack moved");
    let decoded = NormalizedPack::from_bytes_verified(
        NormalizedPack::from_bytes(&first_bytes)
            .expect("decode material-family pack")
            .content_hash(),
        &first_bytes,
    )
    .expect("verified material-family pack");
    assert_eq!(decoded.claims().claim_count(), 4);

    for (property, source, expects_curve) in [
        ("density", "handbook", false),
        ("magnetic_flux_density", "bh-curve", true),
        ("fatigue_life", "sn-curve", true),
        ("dynamic_viscosity", "lubricant", true),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique {property} claim");
        let claim = claims[0].1;
        assert_eq!(
            matches!(&claim.value, PropertyValue::Curve { .. }),
            expects_curve,
            "unexpected payload kind for {property}"
        );
        assert!(
            claim
                .provenance
                .source
                .contains(&format!("[source:{source}]")),
            "{property} lost source-local provenance: {:?}",
            claim.provenance
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    for source in ["handbook", "bh-curve", "sn-curve", "lubricant"] {
        assert!(
            decisions.contains(&format!("\"subject\":\"source:{source}\"")),
            "missing admission row for {source}"
        );
    }
}

#[test]
fn g3_cli_compiles_committed_aluminum_6061_t6_curve_seed() {
    let compiler_id = "frankensim-matdb-pack-compiler-v3";
    let manifest = workspace_path(ALUMINUM_6061_T6_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Aluminum 6061-T6 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aluminum-6061-t6-first.fsmatpk");
    let second_path = directory.join("aluminum-6061-t6-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Aluminum 6061-T6 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Aluminum 6061-T6 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Aluminum 6061-T6 decision stream moved"
    );
    assert_decision_compiler(&first, compiler_id);

    let first_bytes = fs::read(first_path).expect("read first Aluminum 6061-T6 pack");
    let second_bytes = fs::read(second_path).expect("read second Aluminum 6061-T6 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Aluminum 6061-T6 pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode Aluminum 6061-T6 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Aluminum 6061-T6 pack identity");

    assert_eq!(decoded.pack_id(), "aluminum-6061-t6-cryogenic");
    assert_eq!(decoded.compiler(), compiler_id);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public information")
    );
    assert_eq!(decoded.claims().claim_count(), 3);
    assert!(decoded.joint_statistics().is_empty());

    let expected = [
        (
            "thermal_conductivity",
            77.0,
            83.531_441_947_072_3,
            Dims([1, 1, -3, -1, 0, 0]),
            "0.5 percent curve-fit error",
            "a..i=0.07918,1.0957",
        ),
        (
            "thermal_conductivity",
            293.0,
            154.345_205_650_720,
            Dims([1, 1, -3, -1, 0, 0]),
            "0.5 percent curve-fit error",
            "a..i=0.07918,1.0957",
        ),
        (
            "specific_heat_capacity",
            77.0,
            348.127_924_910_444,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=46.6467,-314.292",
        ),
        (
            "specific_heat_capacity",
            293.0,
            942.911_235_969_257,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=46.6467,-314.292",
        ),
        (
            "young_modulus",
            77.0,
            77.145_050_657_273_1e9,
            Dims([-1, 1, -2, 0, 0, 0]),
            "1 percent curve-fit error",
            "a..e=77.71221,0.01030646",
        ),
        (
            "young_modulus",
            293.0,
            70.358_592_182_729_1e9,
            Dims([-1, 1, -2, 0, 0, 0]),
            "1 percent curve-fit error",
            "a..e=77.71221,0.01030646",
        ),
    ];

    for (property, temperature, expected_value, expected_dims, fit_error_note, coefficient_note) in
        expected
    {
        let canonical = property.replace('_', "-");
        let (_, claim) = decoded.claims().claims_for(&canonical)[0];
        let PropertyValue::Curve { knots, dims, .. } = &claim.value else {
            panic!("{property} must carry a curve");
        };
        assert_eq!(knots.len(), 24);
        assert_eq!(claim.validity.bound("temperature"), Some((77.0, 293.0)));
        let value = common_material_acquisition::sample(
            &decoded,
            &canonical,
            &[("temperature", temperature)],
        );
        assert_eq!(*dims, expected_dims, "{property} dimensions moved");
        let relative_error = (value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "{property} at {temperature} K moved by {relative_error:e} relative"
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.interpolation, InterpolationPolicy::LinearInside);
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NIST_PUBLIC_INFORMATION_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("National Institute of Standards and Technology")
        );
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nist-cryogenic-fit]")
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "aluminum-6061-t6-uns-a96061-temper-t6"
        );
        assert!(observation.method.contains("NIST"));
        assert!(observation.method.contains("linear engineering curve"));
        assert!(observation.caveats.contains(fit_error_note));
        assert!(observation.caveats.contains(coefficient_note));
        assert!(
            observation
                .caveats
                .contains("without a confidence level or degrees of freedom")
        );
    }

    // G3 independent-source evidence only: NASA's 1966 compilation reports
    // 82 W/(m K) at 75 K and 155 W/(m K) at 300 K. These nearby-temperature
    // checks do not replace the exact NIST-derived values stored above.
    for (nist_temperature, nasa_temperature, nasa_value) in
        [(77.0, 75.0, 82.0), (293.0, 300.0, 155.0)]
    {
        let value = common_material_acquisition::sample(
            &decoded,
            "thermal-conductivity",
            &[("temperature", nist_temperature)],
        );
        assert!((nist_temperature - nasa_temperature).abs() <= 7.0);
        let relative_difference = (value - nasa_value).abs() / nasa_value;
        assert!(
            relative_difference <= 0.03,
            "NIST-derived {nist_temperature} K conductivity and NASA {nasa_temperature} K comparison differ by {relative_difference:e}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_ofhc_copper_curve_seed() {
    let compiler_id = "frankensim-matdb-pack-compiler-v3";
    let manifest = workspace_path(OFHC_COPPER_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed OFHC Copper seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("ofhc-copper-first.fsmatpk");
    let second_path = directory.join("ofhc-copper-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first OFHC Copper seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second OFHC Copper seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "OFHC Copper decision stream moved"
    );
    assert_decision_compiler(&first, compiler_id);

    let first_bytes = fs::read(first_path).expect("read first OFHC Copper pack");
    let second_bytes = fs::read(second_path).expect("read second OFHC Copper pack");
    assert_eq!(first_bytes, second_bytes, "OFHC Copper pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode OFHC Copper pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify OFHC Copper pack identity");

    assert_eq!(decoded.pack_id(), "ofhc-copper-cryogenic");
    assert_eq!(decoded.compiler(), compiler_id);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public information")
    );
    assert_eq!(decoded.claims().claim_count(), 2);
    assert!(decoded.joint_statistics().is_empty());

    let expected = [
        (
            "thermal_conductivity",
            77.0,
            547.199_698_079_356,
            Dims([1, 1, -3, -1, 0, 0]),
            "1 percent curve-fit error",
            "a..i=2.2154,-0.47461",
            "ofhc-copper-uns-c10100-c10200-rrr100",
        ),
        (
            "thermal_conductivity",
            293.0,
            396.908_547_137_121,
            Dims([1, 1, -3, -1, 0, 0]),
            "1 percent curve-fit error",
            "a..i=2.2154,-0.47461",
            "ofhc-copper-uns-c10100-c10200-rrr100",
        ),
        (
            "specific_heat_capacity",
            77.0,
            195.920_875_203_320,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=-1.91844,-0.15973",
            "ofhc-copper-uns-c10100-c10200-source-rrr-unspecified",
        ),
        (
            "specific_heat_capacity",
            293.0,
            389.085_653_150_356,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=-1.91844,-0.15973",
            "ofhc-copper-uns-c10100-c10200-source-rrr-unspecified",
        ),
    ];

    for (
        property,
        temperature,
        expected_value,
        expected_dims,
        fit_error_note,
        coefficient_note,
        expected_specimen,
    ) in expected
    {
        let canonical = property.replace('_', "-");
        let (_, claim) = decoded.claims().claims_for(&canonical)[0];
        let PropertyValue::Curve { knots, dims, .. } = &claim.value else {
            panic!("OFHC {property} must carry a curve");
        };
        assert_eq!(knots.len(), 24);
        assert_eq!(claim.validity.bound("temperature"), Some((77.0, 293.0)));
        let value = common_material_acquisition::sample(
            &decoded,
            &canonical,
            &[("temperature", temperature)],
        );
        assert_eq!(*dims, expected_dims, "OFHC {property} dimensions moved");
        let relative_error = (value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "OFHC {property} at {temperature} K moved by {relative_error:e} relative"
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.interpolation, InterpolationPolicy::LinearInside);
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NIST_PUBLIC_INFORMATION_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("Material Properties: OFHC Copper")
        );
        assert!(claim.provenance.source.contains("[source:nist-ofhc-fit]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("OFHC claim observation remains linked");
        assert_eq!(observation.specimen, expected_specimen);
        assert!(observation.method.contains("NIST"));
        assert!(observation.method.contains("linear engineering curve"));
        assert!(observation.caveats.contains(fit_error_note));
        assert!(observation.caveats.contains(coefficient_note));
        assert!(
            observation
                .caveats
                .contains("without a confidence level or degrees of freedom")
        );
    }

    // G3 independent-source evidence only: NASA-CR-134806 reports typical
    // room-temperature OFHC Copper values of 390 W/(m K) and 386 J/(kg K).
    // They remain comparisons and do not replace the NIST-derived claims.
    for (property, nasa_value) in [
        ("thermal_conductivity", 390.0),
        ("specific_heat_capacity", 386.0),
    ] {
        let value = common_material_acquisition::sample(
            &decoded,
            &property.replace('_', "-"),
            &[("temperature", 293.0)],
        );
        let relative_difference = (value - nasa_value).abs() / nasa_value;
        assert!(
            relative_difference <= 0.02,
            "NIST-derived 293 K OFHC {property} and NASA room-temperature comparison differ by {relative_difference:e}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_ptfe_teflon_cryogenic_seed() {
    let manifest = workspace_path(PTFE_TEFLON_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed PTFE/Teflon seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("ptfe-teflon-first.fsmatpk");
    let second_path = directory.join("ptfe-teflon-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first PTFE/Teflon seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second PTFE/Teflon seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "PTFE/Teflon decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first PTFE/Teflon pack");
    let second_bytes = fs::read(second_path).expect("read second PTFE/Teflon pack");
    assert_eq!(first_bytes, second_bytes, "PTFE/Teflon pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode PTFE/Teflon pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify PTFE/Teflon pack identity");

    assert_eq!(decoded.pack_id(), "ptfe-teflon-nist-cryogenic");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public information")
    );
    assert_eq!(decoded.claims().claim_count(), 4);
    assert!(decoded.joint_statistics().is_empty());

    let expected = [
        (
            "thermal_conductivity",
            77.0,
            0.232_391_801_023_681,
            Dims([1, 1, -3, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=2.7380,-30.677",
        ),
        (
            "thermal_conductivity",
            293.0,
            0.272_587_209_362_470,
            Dims([1, 1, -3, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=2.7380,-30.677",
        ),
        (
            "specific_heat_capacity",
            77.0,
            301.115_701_345_352,
            Dims([2, 0, -2, -1, 0, 0]),
            "1.5 percent curve-fit error",
            "a..i=31.88256,-166.51949",
        ),
        (
            "specific_heat_capacity",
            293.0,
            1_015.656_817_896_37,
            Dims([2, 0, -2, -1, 0, 0]),
            "1.5 percent curve-fit error",
            "a..i=31.88256,-166.51949",
        ),
    ];

    for (property, temperature, expected_value, expected_dims, fit_error_note, coefficient_note) in
        expected
    {
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing PTFE/Teflon {property} claim at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PTFE/Teflon {property} at {temperature} K was not an exact-point scalar");
        };
        assert_eq!(
            *dims, expected_dims,
            "PTFE/Teflon {property} dimensions moved"
        );
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "PTFE/Teflon {property} at {temperature} K moved by {relative_error:e} relative"
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NIST_PUBLIC_INFORMATION_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("Material Properties: Teflon")
        );
        assert!(claim.provenance.source.contains("[source:nist-ptfe-fit]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PTFE/Teflon claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "ptfe-teflon-source-grade-and-process-unspecified"
        );
        assert!(observation.method.contains("NIST Teflon"));
        assert!(
            observation
                .method
                .contains("exact-temperature derived scalars")
        );
        assert!(observation.caveats.contains(fit_error_note));
        assert!(observation.caveats.contains(coefficient_note));
        assert!(
            observation
                .caveats
                .contains("data and equation range 4-300 K")
        );
        assert!(
            observation
                .caveats
                .contains("without a confidence level or degrees of freedom")
        );
        assert!(
            observation
                .caveats
                .contains("does not identify resin grade")
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_peek_thermic_plate_seed() {
    let manifest = workspace_path(PEEK_THERMIC_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed PEEK THERMIC seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("peek-thermic-first.fsmatpk");
    let second_path = directory.join("peek-thermic-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first PEEK THERMIC seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second PEEK THERMIC seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "PEEK THERMIC decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first PEEK THERMIC pack");
    let second_bytes = fs::read(second_path).expect("read second PEEK THERMIC pack");
    assert_eq!(first_bytes, second_bytes, "PEEK THERMIC pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode PEEK THERMIC pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify PEEK THERMIC pack identity");

    assert_eq!(decoded.pack_id(), "peek-nasa-thermic-plate-2021");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 9);
    assert!(decoded.joint_statistics().is_empty());

    let conductivity = [
        (300.0, 0.224_458_9),
        (400.0, 0.243_077_8),
        (500.0, 0.265_855_5),
        (525.0, 0.274_943_823_437_5),
    ];
    for (temperature, expected_value) in conductivity {
        let (_, claim) = decoded
            .claims()
            .claims_for("thermal_conductivity")
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing PEEK conductivity at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PEEK conductivity at {temperature} K was not scalar");
        };
        assert_eq!(*dims, Dims([1, 1, -3, -1, 0, 0]));
        assert_eq!((*value).to_bits(), f64::to_bits(expected_value));
        assert_eq!(
            claim.validity.bound("source_pressure_atmospheric"),
            Some((1.0, 1.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA/TM-20210014330"));
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nasa-thermic-peek]")
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PEEK conductivity observation remains linked");
        assert_eq!(
            observation.specimen,
            "nasa-larc-thermic-peek-plate-grade-and-process-unspecified"
        );
        assert!(observation.method.contains("Continuous Genetic Algorithm"));
        assert!(observation.caveats.contains("c0..c3=-4.0607e-2"));
        assert!(
            observation
                .caveats
                .contains("narrower repeated range governs")
        );
        assert!(observation.caveats.contains("differed by about 3 percent"));
        assert!(observation.caveats.contains("does not identify PEEK grade"));
    }

    let specific_heat = [
        (300.0, 1_058.931),
        (400.0, 1_347.916),
        (500.0, 1_765.685),
        (525.0, 1_897.616_390_625),
    ];
    for (temperature, expected_value) in specific_heat {
        let (_, claim) = decoded
            .claims()
            .claims_for("specific_heat_capacity")
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing PEEK specific heat at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PEEK specific heat at {temperature} K was not scalar");
        };
        assert_eq!(*dims, Dims([2, 0, -2, -1, 0, 0]));
        assert_eq!((*value).to_bits(), f64::to_bits(expected_value));
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PEEK specific-heat observation remains linked");
        assert!(
            observation
                .method
                .contains("differential-scanning-calorimeter")
        );
        assert!(observation.caveats.contains("1.0477e-5*T^3"));
        assert!(
            observation
                .caveats
                .contains("Equation 1's dimensional balance")
        );
        assert!(observation.caveats.contains("no residual, dispersion"));
    }

    let density_claims = decoded.claims().claims_for("density");
    assert_eq!(density_claims.len(), 1);
    let (_, density) = density_claims[0];
    let PropertyValue::Scalar { value, dims } = &density.value else {
        panic!("PEEK density was not scalar");
    };
    assert_eq!((*value).to_bits(), 1_264.0f64.to_bits());
    assert_eq!(*dims, Dims([-3, 1, 0, 0, 0, 0]));
    assert_eq!(
        density.validity.bound("source_test_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(density.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(density.provenance.license, NASA_SEED_LICENSE);
    let density_observation = decoded
        .claims()
        .observation(density.observations[0])
        .expect("PEEK density observation remains linked");
    assert!(density_observation.method.contains("Commercial-laboratory"));
    assert!(
        density_observation
            .caveats
            .contains("Netzsch report 621004797")
    );
    assert!(density_observation.caveats.contains("test temperature"));

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_nasa_cr_115153_water_ethylene_glycol_seed() {
    let manifest = workspace_path(NASA_CR_115153_WATER_ETHYLENE_GLYCOL_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NASA-CR-115153 water/ethylene-glycol seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nasa-cr-115153-water-glycol-first.fsmatpk");
    let second_path = directory.join("nasa-cr-115153-water-glycol-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-CR-115153 water/glycol seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-CR-115153 water/glycol seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NASA-CR-115153 water/glycol decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NASA water/glycol pack");
    let second_bytes = fs::read(second_path).expect("read second NASA water/glycol pack");
    assert_eq!(
        first_bytes, second_bytes,
        "NASA-CR-115153 water/glycol pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode NASA water/glycol pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NASA water/glycol pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nasa-cr-115153-inhibited-water-ethylene-glycol-coolant"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 11);
    assert!(decoded.joint_statistics().is_empty());

    let expected_bounds = [
        ("sodium_nitrite_mass_fraction_lower_bound", 0.10),
        ("sodium_nitrite_mass_fraction_upper_bound", 0.25),
        ("sodium_benzoate_mass_fraction_lower_bound", 1.33),
        ("sodium_benzoate_mass_fraction_upper_bound", 1.57),
        ("water_mass_fraction_lower_bound", 36.0),
        ("water_mass_fraction_upper_bound", 38.5),
    ];
    for (property, source_percent) in expected_bounds {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique {property} claim");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("{property} was not scalar");
        };
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "{property} moved by {relative_error:e} relative"
        );
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-CR-115153"));
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nasa-cr-115153-table-5]")
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("composition-bound observation remains linked");
        assert_eq!(
            observation.specimen,
            "nasa-cr-115153-water-ethylene-glycol-sodium-nitrite-sodium-benzoate-solution"
        );
        assert!(observation.method.contains("formulation specification"));
        assert!(observation.caveats.contains("formulation bounds"));
        assert!(observation.caveats.contains("without inventing a midpoint"));
    }
    assert!(
        decoded
            .claims()
            .claims_for("ethylene_glycol_mass_fraction")
            .is_empty(),
        "an unreported ethylene-glycol balance must not be inferred"
    );

    let bulk_properties = [
        (
            "density",
            1_081.246_277_742_31,
            Dims([-3, 1, 0, 0, 0, 0]),
            false,
        ),
        (
            "thermal_conductivity",
            0.380_761_626_601_706,
            Dims([1, 1, -3, -1, 0, 0]),
            true,
        ),
    ];
    for (property, expected_value, expected_dims, btu_conversion) in bulk_properties {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique water/glycol {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("water/glycol {property} was not scalar");
        };
        assert_eq!((*value).to_bits(), f64::to_bits(expected_value));
        assert_eq!(*dims, expected_dims);
        assert_eq!(
            claim.validity.bound("source_test_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_test_pressure_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("bulk-property observation remains linked");
        assert!(observation.caveats.contains("extra SI digits"));
        assert!(observation.caveats.contains("not source precision"));
        if btu_conversion {
            assert_eq!(
                claim.validity.bound("source_btu_convention_known"),
                Some((0.0, 0.0))
            );
            assert!(observation.caveats.contains("Btu_IT=1055.05585262 J"));
        }
    }

    let specific_heat_points = [
        (255.372_222_222_222, 2_805.156),
        (283.15, 4_479.876),
        (310.927_777_777_778, 6_154.596),
    ];
    let specific_heat_claims = decoded.claims().claims_for("specific_heat_capacity");
    assert_eq!(specific_heat_claims.len(), specific_heat_points.len());
    for (temperature, expected_value) in specific_heat_points {
        let (_, claim) = specific_heat_claims
            .iter()
            .copied()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing NASA water/glycol cp at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("NASA water/glycol cp at {temperature} K was not scalar");
        };
        assert_eq!((*value).to_bits(), f64::to_bits(expected_value));
        assert_eq!(*dims, Dims([2, 0, -2, -1, 0, 0]));
        assert_eq!(
            claim.validity.bound("source_test_pressure_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_btu_convention_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("specific-heat observation remains linked");
        assert!(observation.method.contains("three exact temperatures"));
        assert!(
            observation
                .caveats
                .contains("approximately cp=(0.67 + 0.008*T_degF)")
        );
        assert!(
            observation
                .caveats
                .contains("do not expose a continuous law")
        );
    }

    // G3 comparison evidence only: NASA/TM-2019-220019 Table VIII lists a
    // separately sourced, composition-basis-unspecified 50-50 water/ethylene-
    // glycol fluid at 1082 kg/m3 and 0.402 W/(m K). Those condition-mismatched
    // values do not overwrite this NASA-CR-115153 formulation; they only bound
    // a coarse transcription plausibility check.
    let (_, density) = decoded.claims().claims_for("density")[0];
    let PropertyValue::Scalar {
        value: density_value,
        ..
    } = &density.value
    else {
        panic!("NASA water/glycol density was not scalar");
    };
    let density_relative_difference = (*density_value - 1_082.0_f64).abs() / 1_082.0;
    assert!(density_relative_difference <= 0.001);

    let (_, conductivity) = decoded.claims().claims_for("thermal_conductivity")[0];
    let PropertyValue::Scalar {
        value: conductivity_value,
        ..
    } = &conductivity.value
    else {
        panic!("NASA water/glycol conductivity was not scalar");
    };
    let conductivity_relative_difference = (*conductivity_value - 0.402_f64).abs() / 0.402;
    assert!(conductivity_relative_difference <= 0.06);

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_n0602_001_nitrile_jp8_compatibility_seed() {
    let manifest = workspace_path(N0602_001_NITRILE_JP8_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed N0602-001 nitrile seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("n0602-001-nitrile-first.fsmatpk");
    let second_path = directory.join("n0602-001-nitrile-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first N0602-001 nitrile seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second N0602-001 nitrile seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "N0602-001 nitrile decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first N0602-001 nitrile pack");
    let second_bytes = fs::read(second_path).expect("read second N0602-001 nitrile pack");
    assert_eq!(
        first_bytes, second_bytes,
        "N0602-001 nitrile pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode N0602-001 nitrile pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify N0602-001 nitrile pack identity");

    assert_eq!(
        decoded.pack_id(),
        "n0602-001-nitrile-o-ring-jp8-compatibility"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 10);
    assert!(decoded.joint_statistics().is_empty());

    let tga_claims = decoded
        .claims()
        .claims_for("tga_semivolatile_mass_fraction");
    assert_eq!(tga_claims.len(), 1);
    let (_, tga) = tga_claims[0];
    let PropertyValue::Scalar { value, dims } = &tga.value else {
        panic!("N0602-001 TGA semi-volatiles were not scalar");
    };
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    let tga_expected = 10.1_f64 * 0.01;
    assert!((*value - tga_expected).abs() / tga_expected <= 2.0e-15);
    assert_eq!(
        tga.validity.bound("source_tga_temperature_program_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(tga.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(tga.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
    assert!(tga.provenance.source.contains("NTRS 20080003822"));
    assert!(tga.provenance.source.contains("[source:primary]"));
    let tga_observation = decoded
        .claims()
        .observation(tga.observations[0])
        .expect("N0602-001 TGA observation remains linked");
    assert_eq!(
        tga_observation.specimen,
        "n0602-001-nitrile-rubber-o-ring-source-formulation-and-lot-unspecified"
    );
    assert!(tga_observation.method.contains("Thermogravimetric"));
    assert!(tga_observation.caveats.contains("propensity to shrink"));
    assert!(tga_observation.caveats.contains("compound formulation"));

    let absorbed_claims = decoded.claims().claims_for("absorbed_fuel_volume_fraction");
    assert_eq!(absorbed_claims.len(), 2);
    for (source_aromatic_percent, source_absorbed_percent) in [(0.0, 8.7), (25.0, 27.9)] {
        let aromatic_fraction = source_aromatic_percent * 0.01;
        let (_, claim) = absorbed_claims
            .iter()
            .copied()
            .find(|(_, claim)| {
                claim.validity.bound("fuel_aromatic_volume_fraction")
                    == Some((aromatic_fraction, aromatic_fraction))
            })
            .unwrap_or_else(|| {
                panic!("missing N0602-001 fuel absorption at {source_aromatic_percent}% aromatic")
            });
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("N0602-001 fuel absorption was not scalar");
        };
        let expected_value = source_absorbed_percent * 0.01;
        assert!((*value - expected_value).abs() / expected_value <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(
            claim.validity.bound("source_exposure_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_duration_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("N0602-001 fuel-absorption observation remains linked");
        assert!(observation.method.contains("thermal-desorption GC-MS"));
        assert!(
            observation
                .caveats
                .contains("do not define service compatibility")
        );
    }

    for (property, expected_value) in [
        ("jp8_alkane_fuel_polymer_partition_coefficient", 0.120),
        ("jp8_aromatic_fuel_polymer_partition_coefficient", 0.412),
        ("jp8_aromatic_to_alkane_partition_ratio", 3.4),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique N0602-001 {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("N0602-001 {property} was not scalar");
        };
        assert_eq!((*value).to_bits(), f64::to_bits(expected_value));
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(
            claim.validity.bound("fuel_aromatic_volume_fraction"),
            Some((0.0, 25.0_f64 * 0.01))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_duration_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
    }

    let slope_claims = decoded
        .claims()
        .claims_for("jp8_volume_swell_per_aromatic_volume_fraction");
    assert_eq!(slope_claims.len(), 2);
    let mut slopes = slope_claims
        .iter()
        .map(|(_, claim)| match &claim.value {
            PropertyValue::Scalar { value, dims } => {
                assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
                *value
            }
            PropertyValue::Curve { .. } => panic!("N0602-001 slope was not scalar"),
        })
        .collect::<Vec<_>>();
    slopes.sort_by(f64::total_cmp);
    assert_eq!(slopes, vec![0.451, 0.463]);
    assert_ne!(
        slope_claims[0].1.observations[0], slope_claims[1].1.observations[0],
        "conflicting printed slopes must retain distinct observations"
    );
    for (_, claim) in slope_claims {
        assert_eq!(
            claim.validity.bound("fuel_aromatic_volume_fraction"),
            Some((0.0, 25.0_f64 * 0.01))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_duration_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("N0602-001 slope observation remains linked");
        assert!(observation.caveats.contains("conflict"));
    }

    let r_squared = decoded
        .claims()
        .claims_for("jp8_volume_swell_aromatic_fraction_r_squared");
    assert_eq!(r_squared.len(), 1);
    let PropertyValue::Scalar { value, dims } = &r_squared[0].1.value else {
        panic!("N0602-001 R-squared was not scalar");
    };
    assert_eq!((*value).to_bits(), 0.948f64.to_bits());
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    assert_eq!(
        r_squared[0]
            .1
            .validity
            .bound("source_exposure_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        r_squared[0].1.provenance.license,
        PUBLIC_USE_PERMITTED_LICENSE
    );

    let intercepts = decoded
        .claims()
        .claims_for("jp8_volume_swell_zero_aromatic_intercept");
    assert_eq!(intercepts.len(), 1);
    let PropertyValue::Scalar { value, dims } = &intercepts[0].1.value else {
        panic!("N0602-001 regression intercept was not scalar");
    };
    let expected_intercept = -1.167_f64 * 0.01;
    assert!((*value - expected_intercept).abs() / expected_intercept.abs() <= 2.0e-15);
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    assert_eq!(
        intercepts[0]
            .1
            .validity
            .bound("source_exposure_duration_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        intercepts[0].1.provenance.license,
        PUBLIC_USE_PERMITTED_LICENSE
    );
    let intercept_observation = decoded
        .claims()
        .observation(intercepts[0].1.observations[0])
        .expect("N0602-001 intercept observation remains linked");
    assert!(
        intercept_observation
            .caveats
            .contains("not a certified shrinkage value")
    );
    assert!(
        decoded
            .claims()
            .claims_for("jp8_prediction_interval_overlap")
            .is_empty(),
        "the source's approximate 57% overlap must remain observation-only"
    );

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_nasa_tn_d_8184_m19_material_deck_without_inventing_process_state() {
    let manifest = workspace_path(NASA_TN_D_8184_M19_MATERIAL_DECK_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NASA-TN-D-8184 M-19 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nasa-tn-d-8184-m19-first.fsmatpk");
    let second_path = directory.join("nasa-tn-d-8184-m19-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-TN-D-8184 M-19 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-TN-D-8184 M-19 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NASA-TN-D-8184 M-19 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NASA-TN-D-8184 M-19 pack");
    let second_bytes = fs::read(second_path).expect("read second NASA-TN-D-8184 M-19 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "NASA-TN-D-8184 M-19 pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode NASA-TN-D-8184 M-19 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NASA-TN-D-8184 M-19 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nasa-tn-d-8184-m19-silicon-steel-material-deck"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("public use"));
    assert_eq!(decoded.claims().claim_count(), 6);
    assert!(decoded.joint_statistics().is_empty());

    let magnetization = decoded.claims().claims_for("magnetic_flux_density");
    assert_eq!(magnetization.len(), 1);
    let magnetization = magnetization[0].1;
    let PropertyValue::Curve {
        abscissa,
        abscissa_dims,
        knots,
        dims,
    } = &magnetization.value
    else {
        panic!("NASA-TN-D-8184 M-19 magnetization data was not a curve");
    };
    assert_eq!(abscissa, "magnetic_field_strength");
    assert_eq!(*abscissa_dims, Dims([-1, 0, 0, 0, 1, 0]));
    assert_eq!(*dims, Dims([0, 1, -2, 0, -1, 0]));
    assert_eq!(knots.len(), 14);
    assert_eq!(
        magnetization.interpolation,
        InterpolationPolicy::TabulatedOnly
    );
    assert_eq!(magnetization.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(magnetization.provenance.license, NASA_SEED_LICENSE);
    assert!(magnetization.provenance.source.contains("NASA-TN-D-8184"));
    assert!(magnetization.provenance.source.contains("[source:primary]"));
    for missing_axis in [
        "source_manufacturer_known",
        "source_processing_anneal_state_known",
        "source_lamination_thickness_known",
        "source_magnetic_test_method_known",
        "source_test_frequency_known",
        "source_test_temperature_known",
        "source_chemistry_known",
        "source_waveform_known",
        "source_direction_known",
    ] {
        assert_eq!(
            magnetization.validity.bound(missing_axis),
            Some((0.0, 0.0)),
            "M-19 B-H curve must retain missing identity axis {missing_axis}"
        );
    }
    assert_eq!(
        magnetization
            .validity
            .bound("source_curve_points_printed_not_digitized"),
        Some((1.0, 1.0))
    );

    let source_b_kilolines_per_square_inch = [
        26.0_f64, 30.0, 40.0, 50.0, 60.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 110.0, 116.0,
    ];
    let source_h_ampere_turns_per_inch = [
        1.30_f64, 1.45, 1.95, 2.55, 3.50, 5.1, 6.5, 8.8, 13.0, 21.0, 37.0, 60.0, 130.0, 185.0,
    ];
    for ((actual_h, actual_b), (source_h, source_b)) in knots.iter().zip(
        source_h_ampere_turns_per_inch
            .iter()
            .zip(source_b_kilolines_per_square_inch.iter()),
    ) {
        let expected_h = source_h / 0.0254_f64;
        let expected_b = source_b * 1.0e-5_f64 / 0.0254_f64.powi(2);
        assert!((actual_h - expected_h).abs() / expected_h <= 2.0e-14);
        assert!((actual_b - expected_b).abs() / expected_b <= 2.0e-14);
    }

    let scalar_expectations = [
        (
            "specific_core_loss",
            9.4_f64 / 0.453_592_37_f64,
            Dims([2, 0, -3, 0, 0, 0]),
        ),
        (
            "core_loss_frequency_power_law_exponent",
            1.47_f64,
            Dims([0, 0, 0, 0, 0, 0]),
        ),
        (
            "lamination_thickness",
            0.014_f64 * 0.0254_f64,
            Dims([1, 0, 0, 0, 0, 0]),
        ),
        (
            "core_loss_reference_frequency",
            400.0_f64,
            Dims([0, 0, -1, 0, 0, 0]),
        ),
        (
            "core_loss_reference_flux_density",
            64.5_f64 * 1.0e-5_f64 / 0.0254_f64.powi(2),
            Dims([0, 1, -2, 0, -1, 0]),
        ),
    ];
    let mut scalar_observations = Vec::new();
    for (property, expected_value, expected_dims) in scalar_expectations {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing NASA M-19 {property}");
        let claim = claims[0].1;
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("NASA M-19 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let scale = expected_value.abs().max(1.0e-12_f64);
        assert!((*value - expected_value).abs() / scale <= 2.0e-15);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TN-D-8184"));
        assert_eq!(
            claim.validity.bound("source_material_grade_m19"),
            Some((1.0, 1.0))
        );
        for missing_axis in [
            "source_manufacturer_known",
            "source_chemistry_known",
            "source_processing_anneal_state_known",
            "source_magnetic_test_method_known",
            "source_test_temperature_known",
            "source_waveform_known",
            "source_direction_known",
        ] {
            assert_eq!(
                claim.validity.bound(missing_axis),
                Some((0.0, 0.0)),
                "NASA M-19 {property} must retain missing identity axis {missing_axis}"
            );
        }
        scalar_observations.push(claim.observations[0]);
    }
    assert!(
        scalar_observations
            .windows(2)
            .all(|pair| pair[0] == pair[1]),
        "NASA M-19 frequency-loss parameters must share one source observation"
    );

    let curve_observation = decoded
        .claims()
        .observation(magnetization.observations[0])
        .expect("NASA M-19 curve observation remains linked");
    assert!(curve_observation.method.contains("Figure 10"));
    assert!(curve_observation.caveats.contains("tabulated-only"));
    assert!(curve_observation.caveats.contains("anneal"));
    let loss_observation = decoded
        .claims()
        .observation(scalar_observations[0])
        .expect("NASA M-19 frequency-loss observation remains linked");
    assert!(loss_observation.method.contains("WCORE=9.4 W/lb"));
    assert!(
        loss_observation
            .caveats
            .contains("not a complete Steinmetz law")
    );
    assert!(loss_observation.caveats.contains("test method"));

    for refused_property in [
        "steinmetz_coefficient",
        "core_loss_flux_density_power_law_exponent",
        "recoil_relative_permeability",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent NASA M-19 property must remain refused: {refused_property}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_tempel_24n208_m19_rating_without_fusing_material_deck() {
    let manifest = workspace_path(NASA_CR_4538_TEMPEL_24N208_M19_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NASA-CR-4538 Tempel 24N208 M19 manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nasa-cr-4538-tempel-24n208-first.fsmatpk");
    let second_path = directory.join("nasa-cr-4538-tempel-24n208-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Tempel 24N208 M19 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Tempel 24N208 M19 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Tempel 24N208 M19 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Tempel 24N208 M19 pack");
    let second_bytes = fs::read(second_path).expect("read second Tempel 24N208 M19 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Tempel 24N208 M19 pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode Tempel 24N208 M19 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Tempel 24N208 M19 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nasa-cr-4538-tempel-24n208-annealed-m19-rating"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("public use"));
    assert_eq!(decoded.claims().claim_count(), 3);
    assert!(decoded.joint_statistics().is_empty());

    let expectations = [
        (
            "specific_hysteresis_loss_rating",
            2.08_f64 / 0.453_592_37_f64,
            Dims([2, 0, -3, 0, 0, 0]),
        ),
        (
            "lamination_thickness",
            0.025_f64 * 0.0254_f64,
            Dims([1, 0, 0, 0, 0, 0]),
        ),
        (
            "nominal_silicon_mass_fraction",
            3.0_f64 * 0.01_f64,
            Dims::NONE,
        ),
    ];
    let mut observation_ids = Vec::new();
    for (property, expected_value, expected_dims) in expectations {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing Tempel 24N208 {property}");
        let claim = claims[0].1;
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Tempel 24N208 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let scale = expected_value.abs().max(1.0e-12);
        assert!((*value - expected_value).abs() / scale <= 2.0e-15);
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-CR-4538"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        for identity_axis in [
            "source_manufacturer_tempel_steel",
            "source_product_24n208",
            "source_material_grade_m19",
            "source_nonoriented_state",
            "source_annealed_state",
        ] {
            assert_eq!(claim.validity.bound(identity_axis), Some((1.0, 1.0)));
        }
        for missing_axis in ["source_product_lot_known", "source_anneal_schedule_known"] {
            assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
        }
        observation_ids.push(claim.observations[0]);
    }
    assert!(
        observation_ids.windows(2).all(|pair| pair[0] == pair[1]),
        "Tempel identity, thickness, and rating must share one source observation"
    );

    let loss_claim = decoded
        .claims()
        .claims_for("specific_hysteresis_loss_rating")[0]
        .1;
    assert_eq!(loss_claim.validity.bound("frequency"), Some((60.0, 60.0)));
    assert_eq!(
        loss_claim.validity.bound("magnetic_flux_density"),
        Some((1.5, 1.5))
    );
    assert_eq!(
        loss_claim.validity.bound("lamination_thickness"),
        Some((0.000_635, 0.000_635))
    );
    assert_eq!(
        loss_claim.validity.bound("source_with_grain_fraction"),
        Some((0.5, 0.5))
    );
    assert_eq!(
        loss_claim
            .validity
            .bound("source_nominal_silicon_mass_fraction"),
        Some((0.03, 0.03))
    );
    assert_eq!(
        loss_claim.validity.bound("source_manufacturer_rating"),
        Some((1.0, 1.0))
    );
    for missing_axis in [
        "source_report_author_measurement",
        "source_surface_insulation_known",
        "source_magnetic_test_method_known",
        "source_waveform_known",
        "source_test_temperature_known",
        "source_rating_bound_semantics_known",
        "source_loss_includes_eddy_current_known",
        "source_loss_is_hysteresis_only_known",
        "source_repeats_and_dispersion_known",
    ] {
        assert_eq!(loss_claim.validity.bound(missing_axis), Some((0.0, 0.0)));
    }

    let observation = decoded
        .claims()
        .observation(observation_ids[0])
        .expect("Tempel 24N208 observation remains linked");
    assert_eq!(
        observation.specimen,
        "tempel-steel-company-24n208-nonoriented-annealed-nominal-3pct-silicon-steel-aisi-m19-lot-unspecified"
    );
    assert!(observation.method.contains("Hysteresis Loss, Laminations"));
    for retained_source_text in ["2.08 W/lbm", "15 kG", "60 Hz", "50 percent w/ the grain"] {
        assert!(observation.caveats.contains(retained_source_text));
    }
    assert!(
        observation
            .caveats
            .contains("does not say whether the rating")
    );
    assert!(
        observation
            .caveats
            .contains("not fused with NASA-TN-D-8184")
    );

    for refused_property in [
        "specific_core_loss",
        "magnetic_flux_density",
        "steinmetz_coefficient",
        "core_loss_frequency_power_law_exponent",
        "core_loss_flux_density_power_law_exponent",
        "core_loss_curve",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "Tempel point rating crossed the {refused_property} no-claim boundary"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_torrent_2018_m19_steinmetz_inputs_without_cross_state_fusion() {
    let manifest = workspace_path(TORRENT_2018_M19_STEINMETZ_INPUTS_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Torrent 2018 M19 Steinmetz-input manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("torrent-2018-m19-steinmetz-inputs-first.fsmatpk");
    let second_path = directory.join("torrent-2018-m19-steinmetz-inputs-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Torrent 2018 M19 Steinmetz-input compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Torrent 2018 M19 Steinmetz-input compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Torrent 2018 M19 Steinmetz-input decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes =
        fs::read(first_path).expect("read first Torrent 2018 M19 Steinmetz-input pack");
    let second_bytes =
        fs::read(second_path).expect("read second Torrent 2018 M19 Steinmetz-input pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Torrent 2018 M19 Steinmetz-input pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes)
        .expect("decode Torrent 2018 M19 Steinmetz-input pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Torrent 2018 M19 Steinmetz-input pack identity");

    assert_eq!(decoded.pack_id(), "torrent-2018-m19-steinmetz-inputs");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 10);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless_expectations: [(&str, f64); 9] = [
        ("equation_4_reported_k_h_numeric", 4.8),
        ("equation_4_frequency_exponent_a", 1.2),
        ("equation_4_flux_density_exponent_n", 2.0),
        ("equation_4_reported_output_scale_numeric", 0.01),
        ("equation_5_reported_k_f_numeric", 60.0),
        ("equation_5_frequency_exponent_x", 2.05),
        ("equation_5_thickness_exponent_y", 2.0),
        ("equation_5_flux_density_exponent_z", 2.05),
        ("equation_5_reported_output_scale_numeric", 100.0),
    ];
    let mut equation_4_observations = Vec::new();
    let mut equation_5_observations = Vec::new();
    for (property, expected_value) in dimensionless_expectations {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing Torrent M19 input {property}");
        let claim = claims[0].1;
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Torrent M19 input {property} was not scalar");
        };
        assert_eq!(*dims, Dims::NONE);
        let scale = expected_value.abs().max(1.0e-12_f64);
        assert!((*value - expected_value).abs() / scale <= 2.0e-15);
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(claim.provenance.source.contains("10.3390/en11061549"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        assert_eq!(
            claim.validity.bound("source_fit_frequency"),
            Some((50.0, 1000.0))
        );
        assert_eq!(
            claim.validity.bound("source_fit_flux_density"),
            Some((0.1, 1.5))
        );
        for retained_axis in [
            "source_material_nomenclature_is_m19_m290_50a",
            "source_excitation_is_sinusoidal",
            "source_manufacturer_cogent_electrical_steel",
            "source_is_reported_fit_not_executable_pack_model",
        ] {
            assert_eq!(claim.validity.bound(retained_axis), Some((1.0, 1.0)));
        }
        for missing_axis in [
            "source_product_process_anneal_coating_lot_known",
            "source_magnetic_loss_test_method_and_temperature_known",
            "source_fit_uncertainty_dispersion_known",
            "source_coefficients_portable_without_source_equations",
            "source_bh_curve_reported",
        ] {
            assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
        }
        if property.starts_with("equation_4_") {
            equation_4_observations.push(claim.observations[0]);
        } else {
            equation_5_observations.push(claim.observations[0]);
        }
    }

    let thickness = decoded.claims().claims_for("equation_5_sheet_thickness_e");
    assert_eq!(thickness.len(), 1);
    let thickness = thickness[0].1;
    let PropertyValue::Scalar { value, dims } = &thickness.value else {
        panic!("Torrent M19 Equation 5 sheet thickness was not scalar");
    };
    assert_eq!((*value).to_bits(), 0.0005_f64.to_bits());
    assert_eq!(*dims, Dims([1, 0, 0, 0, 0, 0]));
    assert_eq!(thickness.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(thickness.provenance.license, CC_BY_4_0_LICENSE);
    assert!(thickness.provenance.source.contains("10.3390/en11061549"));
    assert!(thickness.provenance.source.contains("[source:primary]"));
    assert_eq!(
        thickness.validity.bound("source_fit_frequency"),
        Some((50.0, 1000.0))
    );
    assert_eq!(
        thickness.validity.bound("source_fit_flux_density"),
        Some((0.1, 1.5))
    );
    for retained_axis in [
        "source_material_nomenclature_is_m19_m290_50a",
        "source_excitation_is_sinusoidal",
        "source_manufacturer_cogent_electrical_steel",
        "source_is_reported_fit_not_executable_pack_model",
    ] {
        assert_eq!(thickness.validity.bound(retained_axis), Some((1.0, 1.0)));
    }
    for missing_axis in [
        "source_product_process_anneal_coating_lot_known",
        "source_magnetic_loss_test_method_and_temperature_known",
        "source_fit_uncertainty_dispersion_known",
        "source_coefficients_portable_without_source_equations",
        "source_bh_curve_reported",
    ] {
        assert_eq!(thickness.validity.bound(missing_axis), Some((0.0, 0.0)));
    }
    equation_5_observations.push(thickness.observations[0]);

    assert!(
        equation_4_observations
            .windows(2)
            .all(|pair| pair[0] == pair[1]),
        "Torrent Equation 4 inputs must share one source observation"
    );
    assert!(
        equation_5_observations
            .windows(2)
            .all(|pair| pair[0] == pair[1]),
        "Torrent Equation 5 inputs must share one source observation"
    );
    assert_ne!(
        equation_4_observations[0], equation_5_observations[0],
        "the two reported source equations must retain distinct observations"
    );

    let equation_4_observation = decoded
        .claims()
        .observation(equation_4_observations[0])
        .expect("Torrent Equation 4 observation remains linked");
    assert_eq!(
        equation_4_observation.specimen,
        "torrent-2018-ave-induction-motor-stator-cogent-electrical-steel-aisi-m19-m290-50a-product-and-process-state-unstated"
    );
    assert!(equation_4_observation.method.contains("Equation 4"));
    assert!(
        equation_4_observation
            .caveats
            .contains("identifies Cogent Electrical Steel")
    );
    assert!(
        equation_4_observation
            .caveats
            .contains("does not identify a Cogent product designation")
    );
    assert!(
        equation_4_observation
            .caveats
            .contains("not a portable executable model")
    );
    let equation_5_observation = decoded
        .claims()
        .observation(equation_5_observations[0])
        .expect("Torrent Equation 5 observation remains linked");
    assert_eq!(
        equation_5_observation.specimen,
        equation_4_observation.specimen
    );
    assert!(equation_5_observation.method.contains("Equation 5"));
    assert!(
        equation_5_observation
            .caveats
            .contains("not a portable executable model")
    );
    assert!(
        equation_5_observation
            .caveats
            .contains("separate NASA and Tempel M-19 states")
    );

    for refused_property in [
        "specific_core_loss",
        "specific_hysteresis_loss",
        "magnetization_curve",
        "magnetic_flux_density",
        "bh_curve",
        "steinmetz_model",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "reported Torrent fit input crossed the {refused_property} no-claim boundary"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_ngyc_n42_sintered_magnet_seed() {
    let manifest = workspace_path(NGYC_N42_SINTERED_NICKEL_COATED_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NGYC N42 magnet seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("ngyc-n42-first.fsmatpk");
    let second_path = directory.join("ngyc-n42-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NGYC N42 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NGYC N42 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NGYC N42 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NGYC N42 pack");
    let second_bytes = fs::read(second_path).expect("read second NGYC N42 pack");
    assert_eq!(first_bytes, second_bytes, "NGYC N42 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode NGYC N42 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NGYC N42 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "ngyc-n42-sintered-ndfeb-nickel-coated-cubes"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 4);
    assert!(decoded.joint_statistics().is_empty());

    let remanence = decoded.claims().claims_for("remanent_flux_density");
    assert_eq!(remanence.len(), 1);
    let PropertyValue::Scalar { value, dims } = &remanence[0].1.value else {
        panic!("NGYC N42 remanence was not scalar");
    };
    let expected_remanence = 1350.0_f64 * 1.0e-3;
    assert!((*value - expected_remanence).abs() / expected_remanence <= 2.0e-15);
    assert_eq!(*dims, Dims([0, 1, -2, 0, -1, 0]));

    let coercivity = decoded.claims().claims_for("coercive_field_strength");
    assert_eq!(coercivity.len(), 1);
    let PropertyValue::Scalar { value, dims } = &coercivity[0].1.value else {
        panic!("NGYC N42 coercivity was not scalar");
    };
    assert_eq!((*value).to_bits(), (923.0_f64 * 1.0e3).to_bits());
    assert_eq!(*dims, Dims([-1, 0, 0, 0, 1, 0]));

    let energy_products = decoded
        .claims()
        .claims_for("maximum_magnetic_energy_product");
    assert_eq!(energy_products.len(), 2);
    assert_ne!(
        energy_products[0].1.observations[0], energy_products[1].1.observations[0],
        "conflicting printed energy products must retain distinct observations"
    );
    let mut energy_values = energy_products
        .iter()
        .map(|(_, claim)| match &claim.value {
            PropertyValue::Scalar { value, dims } => {
                assert_eq!(*dims, Dims([-1, 1, -2, 0, 0, 0]));
                *value
            }
            PropertyValue::Curve { .. } => panic!("NGYC N42 energy product was not scalar"),
        })
        .collect::<Vec<_>>();
    energy_values.sort_by(f64::total_cmp);
    let printed_si = 318.3_f64 * 1.0e3;
    let normalized_42_mgoe = 42.0_f64 * (100_000.0 / (4.0 * std::f64::consts::PI));
    assert!((energy_values[0] - printed_si).abs() / printed_si <= 2.0e-15);
    assert!((energy_values[1] - normalized_42_mgoe).abs() / normalized_42_mgoe <= 2.0e-15);
    assert!(energy_values[1] > energy_values[0]);

    for property in [
        "remanent_flux_density",
        "coercive_field_strength",
        "maximum_magnetic_energy_product",
    ] {
        for (_, claim) in decoded.claims().claims_for(property) {
            assert_eq!(
                claim
                    .validity
                    .bound("source_magnetic_test_temperature_known"),
                Some((0.0, 0.0))
            );
            assert_eq!(
                claim.validity.bound("source_magnetic_test_method_known"),
                Some((0.0, 0.0))
            );
            assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
            assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
            assert!(
                claim
                    .provenance
                    .source
                    .contains("10.1038/s41598-023-47689-2")
            );
            assert!(claim.provenance.source.contains("[source:primary]"));
            let observation = decoded
                .claims()
                .observation(claim.observations[0])
                .expect("NGYC N42 observation remains linked");
            assert_eq!(
                observation.specimen,
                "ngyc-yinxian-ningbo-n42-sintered-ndfeb-nickel-coated-cubes-paper-lot-unspecified"
            );
        }
    }

    let si_observation = decoded
        .claims()
        .observation(remanence[0].1.observations[0])
        .expect("NGYC N42 SI observation remains linked");
    assert!(si_observation.method.contains("Telfah et al. 2023"));
    assert!(si_observation.caveats.contains("supplier nominal values"));
    assert!(si_observation.caveats.contains("not SI-equivalent"));
    let cgs_observation = energy_products
        .iter()
        .find_map(|(_, claim)| {
            let observation = decoded.claims().observation(claim.observations[0])?;
            observation
                .method
                .contains("Exact unit normalization")
                .then_some(observation)
        })
        .expect("NGYC N42 CGS observation remains linked");
    assert!(cgs_observation.method.contains("1 Oe=1000/(4*pi) A/m"));
    assert!(cgs_observation.caveats.contains("source conflict"));

    for refused_property in [
        "intrinsic_coercive_field_strength",
        "recoil_relative_permeability",
        "remanence_temperature_coefficient",
        "coercivity_temperature_coefficient",
        "demagnetization_curve",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent NGYC N42 property must remain refused: {refused_property}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_jinshan_n42_pristine_temperature_endpoints() {
    let manifest = workspace_path(JINSHAN_N42_PRISTINE_TEMPERATURE_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Jinshan N42 temperature seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("jinshan-n42-temperature-first.fsmatpk");
    let second_path = directory.join("jinshan-n42-temperature-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Jinshan N42 temperature seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Jinshan N42 temperature seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Jinshan N42 temperature decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Jinshan N42 temperature pack");
    let second_bytes = fs::read(second_path).expect("read second Jinshan N42 temperature pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Jinshan N42 temperature pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode Jinshan N42 temperature pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Jinshan N42 temperature pack identity");

    assert_eq!(
        decoded.pack_id(),
        "jinshan-n42-pristine-sintered-temperature-endpoints"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 8);
    assert!(decoded.joint_statistics().is_empty());

    let endpoint_expectations = [
        (
            "remanent_flux_density",
            25.0_f64,
            12.75_f64 * 1.0e3 * 1.0e-4,
            Dims([0, 1, -2, 0, -1, 0]),
        ),
        (
            "remanent_flux_density",
            120.0_f64,
            11.18_f64 * 1.0e3 * 1.0e-4,
            Dims([0, 1, -2, 0, -1, 0]),
        ),
        (
            "intrinsic_coercive_field_strength",
            25.0_f64,
            12.07_f64 * 1.0e3 * (1.0e3 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 0, 0, 0, 1, 0]),
        ),
        (
            "intrinsic_coercive_field_strength",
            120.0_f64,
            5.17_f64 * 1.0e3 * (1.0e3 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 0, 0, 0, 1, 0]),
        ),
        (
            "maximum_magnetic_energy_product",
            25.0_f64,
            40.14_f64 * (100_000.0 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 1, -2, 0, 0, 0]),
        ),
        (
            "maximum_magnetic_energy_product",
            120.0_f64,
            29.29_f64 * (100_000.0 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 1, -2, 0, 0, 0]),
        ),
    ];
    let endpoint_observation = endpoint_expectations
        .iter()
        .map(
            |(property, source_temperature_c, expected_value, expected_dims)| {
                let temperature_k = source_temperature_c + 273.15;
                let claims = decoded.claims().claims_for(property);
                let (_, claim) = claims
                    .into_iter()
                    .find(|(_, claim)| {
                        claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
                    })
                    .unwrap_or_else(|| {
                        panic!("missing Jinshan N42 {property} at {source_temperature_c} degC")
                    });
                let PropertyValue::Scalar { value, dims } = &claim.value else {
                    panic!("Jinshan N42 {property} endpoint was not scalar");
                };
                assert_eq!(dims, expected_dims);
                let scale = expected_value.abs().max(1.0e-12);
                assert!((*value - expected_value).abs() / scale <= 2.0e-15);
                assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
                assert!(
                    claim
                        .provenance
                        .source
                        .contains("10.1016/j.jmrt.2024.12.235")
                );
                assert!(claim.provenance.source.contains("[source:primary]"));
                assert_eq!(
                    claim.validity.bound("source_instrument_nim_6500c"),
                    Some((1.0, 1.0))
                );
                assert_eq!(
                    claim
                        .validity
                        .bound("source_authors_heat_treatment_applied"),
                    Some((0.0, 0.0))
                );
                for missing_axis in [
                    "source_supplier_production_lot_known",
                    "source_composition_known",
                    "source_temperature_control_method_known",
                ] {
                    assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                }
                claim.observations[0]
            },
        )
        .collect::<Vec<_>>();
    assert!(
        endpoint_observation
            .windows(2)
            .all(|pair| pair[0] == pair[1]),
        "all Table 1 endpoints must retain one shared observation"
    );
    let endpoint_observation = decoded
        .claims()
        .observation(endpoint_observation[0])
        .expect("Jinshan N42 endpoint observation remains linked");
    assert_eq!(
        endpoint_observation.specimen,
        "jinshan-magnetic-materials-commercial-n42-sintered-ndfeb-pristine-wire-cut-10x10x6-mm-lot-unspecified"
    );
    assert!(endpoint_observation.method.contains("NIM 6500C"));
    assert!(
        endpoint_observation
            .caveats
            .contains("calls the 10 mm x 10 mm x 6 mm pieces cubes")
    );
    assert!(
        endpoint_observation
            .caveats
            .contains("no curve points are digitized")
    );

    for (property, expected_value) in [
        ("remanence_temperature_coefficient", -0.00129_f64),
        ("coercivity_temperature_coefficient", -0.00602_f64),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing Jinshan N42 {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Jinshan N42 {property} was not scalar");
        };
        assert_eq!(*dims, Dims([0, 0, 0, -1, 0, 0]));
        assert_eq!((*value).to_bits(), f64::to_bits(expected_value));
        assert_eq!(
            claim.validity.bound("temperature"),
            Some((25.0 + 273.15, 120.0 + 273.15))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_coefficient_uses_25c_and_120c_endpoints"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_rounded_endpoints_reproduce_printed_coefficient_exactly"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("Jinshan N42 coefficient observation remains linked");
        assert!(
            observation
                .caveats
                .contains("do not reproduce both printed coefficients exactly")
        );
        assert!(
            observation
                .caveats
                .contains("not continuous constitutive laws")
        );
    }

    for refused_property in [
        "coercive_field_strength",
        "recoil_relative_permeability",
        "demagnetization_curve",
        "irreversible_demagnetization_loss_boundary",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent Jinshan N42 property must remain refused: {refused_property}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_y30_catalog_model_inputs_without_recoil_transfer() {
    let manifest = workspace_path(SJOLUND_2020_Y30_CATALOG_MODEL_INPUTS_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Sjolund Y30 catalog-model manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("sjolund-y30-catalog-first.fsmatpk");
    let second_path = directory.join("sjolund-y30-catalog-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Sjolund Y30 catalog compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Sjolund Y30 catalog compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Sjolund Y30 catalog decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Sjolund Y30 catalog pack");
    let second_bytes = fs::read(second_path).expect("read second Sjolund Y30 catalog pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Sjolund Y30 catalog pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode Sjolund Y30 catalog pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Sjolund Y30 catalog pack identity");

    assert_eq!(decoded.pack_id(), "sjolund-2020-y30-catalog-model-inputs");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 5);
    assert!(decoded.joint_statistics().is_empty());

    let mu_0 = 4.0_f64 * std::f64::consts::PI * 1.0e-7;
    let model_relative_permeability = (0.385_f64 * 0.385_f64) / (4.0 * mu_0 * 28_000.0);
    let expectations = [
        (
            "remanent_flux_density",
            385.0_f64 * 1.0e-3,
            Dims([0, 1, -2, 0, -1, 0]),
        ),
        (
            "coercive_field_strength",
            192.5_f64 * 1.0e3,
            Dims([-1, 0, 0, 0, 1, 0]),
        ),
        (
            "intrinsic_coercive_field_strength",
            200.0_f64 * 1.0e3,
            Dims([-1, 0, 0, 0, 1, 0]),
        ),
        (
            "maximum_magnetic_energy_product",
            28.0_f64 * 1.0e3,
            Dims([-1, 1, -2, 0, 0, 0]),
        ),
        (
            "model_relative_permeability",
            model_relative_permeability,
            Dims::NONE,
        ),
    ];
    let mut observation_ids = Vec::new();
    for (property, expected_value, expected_dims) in expectations {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing Sjolund Y30 {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Sjolund Y30 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let scale = expected_value.abs().max(1.0e-12);
        assert!((*value - expected_value).abs() / scale <= 2.0e-15);
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(claim.provenance.source.contains("10.1063/1.5129303"));
        assert!(
            claim
                .provenance
                .source
                .contains("[source:published-table-ii]")
        );
        assert_eq!(
            claim.validity.bound("source_catalog_grade_y30"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim.validity.bound("source_catalog_midpoint_used"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim.validity.bound("source_simulation_temperature"),
            Some((293.15, 293.15))
        );
        for missing_axis in [
            "source_physical_product_identified",
            "source_product_supplier_identified",
            "source_production_lot_known",
            "source_composition_known",
            "source_sinter_process_known",
            "source_magnetic_test_temperature_known",
            "source_magnetic_test_method_known",
        ] {
            assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
        }
        observation_ids.push((property, claim.observations[0]));
    }

    let catalog_observation_id = observation_ids[0].1;
    assert!(
        observation_ids[..4]
            .iter()
            .all(|(_, observation)| *observation == catalog_observation_id),
        "the four Table II midpoint claims must retain one shared observation"
    );
    let catalog_observation = decoded
        .claims()
        .observation(catalog_observation_id)
        .expect("Sjolund Y30 Table II observation remains linked");
    assert_eq!(
        catalog_observation.specimen,
        "e-magnetsuk-y30-online-grade-family-accessed-2019-product-lot-process-unspecified"
    );
    assert!(
        catalog_observation
            .method
            .contains("midpoint plus or minus half-range")
    );
    for printed_range in [
        "Br 385 plus or minus 15 mT",
        "Hcb 192.5 plus or minus 17.5 kA/m",
        "Hcj 200 plus or minus 20 kA/m",
        "BHmax 28 plus or minus 2 kJ/m3",
    ] {
        assert!(catalog_observation.caveats.contains(printed_range));
    }
    assert!(
        catalog_observation
            .caveats
            .contains("not laundered into a material measurement temperature")
    );

    let model_claim = decoded.claims().claims_for("model_relative_permeability");
    let model_claim = model_claim[0].1;
    assert_ne!(
        model_claim.observations[0], catalog_observation_id,
        "the Equation 2 derivation must retain its own observation"
    );
    assert_eq!(
        model_claim
            .validity
            .bound("source_model_mu_equation_2_derived"),
        Some((1.0, 1.0))
    );
    assert_eq!(
        model_claim
            .validity
            .bound("source_model_mu_is_measured_recoil_mu"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        model_claim
            .validity
            .bound("source_minor_loop_recoil_data_known"),
        Some((0.0, 0.0))
    );
    let model_observation = decoded
        .claims()
        .observation(model_claim.observations[0])
        .expect("Sjolund Y30 Equation 2 observation remains linked");
    assert!(model_observation.method.contains("Equation 2"));
    assert!(
        model_observation
            .caveats
            .contains("not measured recoil permeability")
    );
    assert!(
        model_observation
            .caveats
            .contains("no adequate Y30 demagnetization curve")
    );

    for refused_property in [
        "recoil_relative_permeability",
        "demagnetization_curve",
        "irreversible_demagnetization_loss_boundary",
        "remanence_temperature_coefficient",
        "coercivity_temperature_coefficient",
        "continuous_demagnetization_temperature_law",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "catalog/model input crossed the {refused_property} no-claim boundary"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_y30_afcp_application_demagnetization_without_intrinsic_transfer() {
    let manifest = workspace_path(KIM_BAEK_2026_Y30_AFCP_DEMAGNETIZATION_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Kim-Baek Y30 application-demagnetization manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("kim-baek-y30-demagnetization-first.fsmatpk");
    let second_path = directory.join("kim-baek-y30-demagnetization-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Kim-Baek Y30 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Kim-Baek Y30 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Kim-Baek Y30 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Kim-Baek Y30 pack");
    let second_bytes = fs::read(second_path).expect("read second Kim-Baek Y30 pack");
    assert_eq!(first_bytes, second_bytes, "Kim-Baek Y30 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode Kim-Baek Y30 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Kim-Baek Y30 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "kim-baek-2026-y30-afcp-application-demagnetization"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 2);
    assert!(decoded.joint_statistics().is_empty());

    let expected: [(f64, f64); 2] = [(20.0, 1.654 * 0.01), (-40.0, 22.396 * 0.01)];
    let mut observation_ids = Vec::new();
    for (source_temperature_c, expected_fraction) in expected {
        let temperature_k = source_temperature_c + 273.15;
        let claims = decoded
            .claims()
            .claims_for("application_model_maximum_demagnetization_fraction");
        let (id, claim) = claims
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
            })
            .unwrap_or_else(|| {
                panic!("missing Y30 application result at {source_temperature_c} degC")
            });
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Y30 application demagnetization claim {id:?} was not scalar");
        };
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!((*value).to_bits(), expected_fraction.to_bits());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(claim.provenance.source.contains("10.3390/app16021094"));
        assert!(claim.provenance.source.contains("[source:primary]"));

        for (axis, expected_bound) in [
            ("rotational_frequency", (100.0, 100.0)),
            ("source_motor_output_power", (750.0, 750.0)),
            (
                "source_optimal_model_magnet_volume",
                (0.00006016, 0.00006016),
            ),
            ("source_current_multiplier_relative_to_rated", (5.0, 5.0)),
            (
                "source_maximum_stator_magnetic_field_strength",
                (256_490.0, 256_490.0),
            ),
            ("source_grade_label_y30", (1.0, 1.0)),
            ("source_result_is_3d_fea", (1.0, 1.0)),
            ("source_coefficient_is_spatial_maximum", (1.0, 1.0)),
            (
                "source_equation_uses_post_field_recoil_flux_density",
                (1.0, 1.0),
            ),
        ] {
            assert_eq!(claim.validity.bound(axis), Some(expected_bound));
        }
        for missing_axis in [
            "source_magnet_supplier_process_composition_lot_known",
            "source_fea_software_mesh_convergence_known",
            "source_bh_curve_points_tabulated",
            "source_recoil_relative_permeability_known",
            "source_experimental_prototype_validation_performed",
            "source_intrinsic_material_limit_claimed",
            "source_uncertainty_dispersion_known",
        ] {
            assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
        }
        observation_ids.push(claim.observations[0]);
    }
    assert_eq!(observation_ids[0], observation_ids[1]);
    let observation = decoded
        .claims()
        .observation(observation_ids[0])
        .expect("Kim-Baek Y30 observation remains linked");
    assert!(observation.method.contains("five times rated current"));
    assert!(
        observation
            .caveats
            .contains("not intrinsic Y30 material allowables")
    );
    assert!(observation.caveats.contains("no tabulated points"));

    for refused_property in [
        "remanent_flux_density",
        "intrinsic_coercive_field_strength",
        "recoil_relative_permeability",
        "demagnetization_curve",
        "irreversible_demagnetization_loss_boundary",
        "continuous_demagnetization_temperature_law",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "application model crossed the intrinsic {refused_property} no-claim boundary"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_naca_tn_2680_isooctane_flame_speed_seed() {
    let manifest = workspace_path(NACA_TN_2680_ISOOCTANE_FLAME_SPEED_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NACA TN 2680 iso-octane flame-speed seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("naca-tn-2680-isooctane-first.fsmatpk");
    let second_path = directory.join("naca-tn-2680-isooctane-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NACA TN 2680 iso-octane seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NACA TN 2680 iso-octane seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NACA TN 2680 iso-octane decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NACA TN 2680 iso-octane pack");
    let second_bytes = fs::read(second_path).expect("read second NACA TN 2680 iso-octane pack");
    assert_eq!(
        first_bytes, second_bytes,
        "NACA TN 2680 iso-octane pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode NACA TN 2680 iso-octane pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NACA TN 2680 iso-octane pack identity");

    assert_eq!(
        decoded.pack_id(),
        "naca-tn-2680-2-2-4-trimethylpentane-flame-speed"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("Work of the US Gov. Public Use Permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 16);
    assert!(decoded.joint_statistics().is_empty());

    let purity = decoded
        .claims()
        .claims_for("minimum_reported_fuel_mole_fraction_purity");
    assert_eq!(purity.len(), 1);
    let (_, purity_claim) = purity[0];
    let PropertyValue::Scalar { value, dims } = &purity_claim.value else {
        panic!("NACA TN 2680 fuel minimum purity was not scalar");
    };
    let expected_purity = 99.6_f64 * 0.01;
    assert!((*value - expected_purity).abs() / expected_purity <= 2.0e-15);
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    assert_eq!(
        purity_claim
            .validity
            .bound("source_fuel_supplier_identity_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        purity_claim.validity.bound("source_fuel_lot_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        purity_claim.validity.bound("source_exact_assay_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(purity_claim.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(purity_claim.provenance.license, NASA_SEED_LICENSE);
    assert!(
        purity_claim
            .provenance
            .source
            .contains("NACA Technical Note 2680")
    );
    let purity_observation = decoded
        .claims()
        .observation(purity_claim.observations[0])
        .expect("NACA TN 2680 purity observation remains linked");
    assert!(
        purity_observation
            .method
            .contains("minimum-purity statement")
    );
    assert!(purity_observation.caveats.contains("lower-bound statement"));

    let flame_claims = decoded.claims().claims_for("maximum_laminar_flame_speed");
    assert_eq!(flame_claims.len(), 15);
    let expected_rows = [
        (311.0, 0.210, 1000.0, 1.256, 34.6),
        (311.0, 0.250, 1000.0, 0.838, 52.1),
        (311.0, 0.294, 1600.0, 0.838, 72.2),
        (311.0, 0.294, 900.0, 0.617, 67.2),
        (311.0, 0.347, 1200.0, 0.617, 89.1),
        (311.0, 0.496, 1800.0, 0.297, 152.2),
        (367.0, 0.210, 1000.0, 1.256, 44.8),
        (422.0, 0.210, 1000.0, 1.256, 56.1),
        (422.0, 0.210, 1000.0, 1.256, 59.0),
        (422.0, 0.210, 700.0, 0.838, 57.4),
        (422.0, 0.250, 1400.0, 0.838, 83.1),
        (422.0, 0.294, 900.0, 0.617, 108.0),
        (422.0, 0.294, 900.0, 0.617, 102.1),
        (422.0, 0.347, 1400.0, 0.617, 138.0),
        (422.0, 0.496, 1800.0, 0.297, 229.9),
    ];

    for (temperature, oxygen_fraction, reynolds_number, diameter_cm, speed_cm_per_s) in
        expected_rows
    {
        let diameter_m = diameter_cm * 0.01;
        let speed_m_per_s = speed_cm_per_s * 0.01;
        let matching = flame_claims
            .iter()
            .filter(|(_, claim)| {
                claim.validity.bound("initial_mixture_temperature")
                    == Some((temperature, temperature))
                    && claim
                        .validity
                        .bound("oxygen_mole_fraction_in_oxygen_nitrogen")
                        == Some((oxygen_fraction, oxygen_fraction))
                    && claim.validity.bound("stream_flow_reynolds_number")
                        == Some((reynolds_number, reynolds_number))
                    && claim.validity.bound("burner_inside_diameter")
                        == Some((diameter_m, diameter_m))
                    && matches!(
                        &claim.value,
                        PropertyValue::Scalar { value, .. }
                            if (*value).to_bits() == f64::to_bits(speed_m_per_s)
                    )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            matching.len(),
            1,
            "expected one NACA Table I row at T={temperature}, O2={oxygen_fraction}, Re={reynolds_number}, diameter={diameter_cm} cm, speed={speed_cm_per_s} cm/s"
        );
        let (_, claim) = matching[0];
        let PropertyValue::Scalar { dims, .. } = &claim.value else {
            unreachable!("matching predicate admitted only scalar flame speeds");
        };
        assert_eq!(*dims, Dims([1, 0, -1, 0, 0, 0]));
        assert_eq!(
            claim.validity.bound("average_atmospheric_pressure"),
            Some((99.2_f64 * 1.0e3, 99.2_f64 * 1.0e3))
        );
        assert_eq!(
            claim.validity.bound("source_pressure_per_row_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_equivalence_ratio_at_maximum_exact_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_maximum_equivalence_ratio_lower_bound"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_maximum_equivalence_ratio_upper_bound"),
            Some((1.1, 1.1))
        );
        assert_eq!(
            claim.validity.bound("source_fuel_supplier_identity_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_fuel_lot_known"),
            Some((0.0, 0.0))
        );
        if oxygen_fraction.to_bits() == 0.210_f64.to_bits() {
            assert_eq!(
                claim
                    .validity
                    .bound("source_oxidizer_analysis_half_width_known"),
                Some((0.0, 0.0))
            );
            assert_eq!(
                claim
                    .validity
                    .bound("source_oxidizer_analysis_absolute_half_width"),
                None
            );
        } else {
            assert_eq!(
                claim
                    .validity
                    .bound("source_oxidizer_analysis_absolute_half_width"),
                Some((0.1_f64 * 0.01, 0.1_f64 * 0.01))
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NTRS 19930083861"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("NACA TN 2680 flame-speed observation remains linked");
        assert!(observation.method.contains("total-area method"));
        assert!(
            observation
                .caveats
                .contains("not geometry-free bulk-material constants")
        );
        assert!(observation.caveats.contains("Repeated rows"));
    }

    for refused_property in [
        "density",
        "dynamic_viscosity",
        "surface_tension",
        "specific_heat_capacity",
        "heat_of_vaporization",
        "vapor_pressure",
        "research_octane_number",
        "empirical_maximum_flame_speed_fit",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent or model-only NACA TN 2680 property must remain refused: {refused_property}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_face_g_cdtrf_g_2023_v1_surrogate_seed() {
    let manifest = workspace_path(FACE_G_CDTRF_G_2023_V1_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed FACE G CDTRF-G 2023 v1 surrogate seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("face-g-cdtrf-g-2023-v1-first.fsmatpk");
    let second_path = directory.join("face-g-cdtrf-g-2023-v1-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first FACE G CDTRF-G seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second FACE G CDTRF-G seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "FACE G CDTRF-G decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first FACE G CDTRF-G pack");
    let second_bytes = fs::read(second_path).expect("read second FACE G CDTRF-G pack");
    assert_eq!(first_bytes, second_bytes, "FACE G CDTRF-G pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode FACE G CDTRF-G pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify FACE G CDTRF-G pack identity");

    assert_eq!(decoded.pack_id(), "face-g-cdtrf-g-2023-v1");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 7);
    assert!(decoded.joint_statistics().is_empty());

    let expected_components = [
        ("isooctane_component_volume_fraction", 23.75),
        ("n_heptane_component_volume_fraction", 19.0),
        ("toluene_component_volume_fraction", 42.75),
        ("diisobutylene_component_volume_fraction", 9.5),
        ("cyclohexane_component_volume_fraction", 5.0),
    ];
    let mut fraction_sum = 0.0;
    for (property, source_percent) in expected_components {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique CDTRF-G {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("CDTRF-G {property} was not scalar");
        };
        let expected_fraction = source_percent * 0.01;
        assert!((*value - expected_fraction).abs() / expected_fraction <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        fraction_sum += *value;
        assert_eq!(
            claim
                .validity
                .bound("source_composition_basis_is_volume_fraction"),
            Some((1.0, 1.0))
        );
        for missing_axis in [
            "source_component_supplier_known",
            "source_component_lot_known",
            "source_component_purity_known",
            "source_mixing_temperature_known",
            "source_mixing_pressure_known",
            "source_volume_contraction_treatment_known",
        ] {
            assert_eq!(
                claim.validity.bound(missing_axis),
                Some((0.0, 0.0)),
                "CDTRF-G {property} must retain missing axis {missing_axis}"
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("10.3390/molecules28114273")
        );
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("CDTRF-G composition observation remains linked");
        assert_eq!(
            observation.specimen,
            "face-g-targeted-cdtrf-g-2023-v1-source-component-lots-and-purities-unspecified"
        );
        assert!(observation.method.contains("Table 2 CDTRF-G"));
        assert!(observation.caveats.contains("sum to exactly 100 percent"));
        assert!(observation.caveats.contains("molar basis"));
    }
    assert!((fraction_sum - 1.0_f64).abs() <= 2.0e-15);

    let ron_claims = decoded
        .claims()
        .claims_for("reported_calculated_research_octane_number");
    assert_eq!(ron_claims.len(), 2);
    assert_ne!(
        ron_claims[0].1.observations[0], ron_claims[1].1.observations[0],
        "conflicting CDTRF-G RON prints must retain distinct observations"
    );
    let mut ron_values = ron_claims
        .iter()
        .map(|(_, claim)| match &claim.value {
            PropertyValue::Scalar { value, dims } => {
                assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
                *value
            }
            PropertyValue::Curve { .. } => panic!("CDTRF-G RON was not scalar"),
        })
        .collect::<Vec<_>>();
    ron_values.sort_by(f64::total_cmp);
    assert_eq!(ron_values, vec![93.9, 94.0]);
    for (_, claim) in ron_claims {
        assert_eq!(
            claim.validity.bound("source_octane_measurement_performed"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_octane_calculation_basis_unambiguous"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_component_purity_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("CDTRF-G RON observation remains linked");
        assert!(observation.caveats.contains("Table 2"));
        assert!(observation.caveats.contains("Table 7"));
        assert!(observation.caveats.contains("separate"));
    }

    for refused_property in [
        "density",
        "dynamic_viscosity",
        "surface_tension",
        "specific_heat_capacity",
        "heat_of_vaporization",
        "vapor_pressure",
        "motor_octane_number",
        "laminar_flame_speed",
        "ignition_delay_time",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent or model-only CDTRF-G property must remain refused: {refused_property}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_wo2018_formulation_8_5w30_seed() {
    let manifest = workspace_path(WO2018_125520_FORMULATION_8_5W30_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed WO 2018/125520 Formulation 8 5W30 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("wo2018-125520-formulation-8-5w30-first.fsmatpk");
    let second_path = directory.join("wo2018-125520-formulation-8-5w30-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first WO 2018/125520 Formulation 8 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second WO 2018/125520 Formulation 8 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "WO 2018/125520 Formulation 8 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first WO 2018 Formulation 8 pack");
    let second_bytes = fs::read(second_path).expect("read second WO 2018 Formulation 8 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "WO 2018/125520 Formulation 8 pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode WO 2018/125520 Formulation 8 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify WO 2018/125520 Formulation 8 pack identity");

    assert_eq!(decoded.pack_id(), "wo2018-125520-formulation-8-5w30");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("not a patent-practice or trademark license")
    );
    assert_eq!(decoded.claims().claim_count(), 12);
    assert!(decoded.joint_statistics().is_empty());

    let expected_components = [
        ("spectrasyn_4_component_mass_fraction", 60.0),
        ("synesstic_5_component_mass_fraction", 10.0),
        ("spectrasyn_elite_150_component_mass_fraction", 18.0),
        ("infineum_p6003_component_mass_fraction", 12.0),
    ];
    let mut fraction_sum = 0.0;
    let mut composition_observation = None;
    for (property, source_percent) in expected_components {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique Formulation 8 {property}");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Formulation 8 {property} was not scalar");
        };
        let expected_fraction = source_percent * 0.01;
        assert!((*value - expected_fraction).abs() <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        fraction_sum += *value;
        assert_eq!(
            claim
                .validity
                .bound("source_composition_basis_is_mass_fraction"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim.validity.bound("source_formulation_number"),
            Some((8.0, 8.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_component_commercial_identifier_known"),
            Some((1.0, 1.0))
        );
        for missing_axis in [
            "source_component_lot_known",
            "source_component_detailed_chemistry_known",
            "source_final_blend_protocol_known",
            "source_patent_practice_license_granted",
        ] {
            assert_eq!(
                claim.validity.bound(missing_axis),
                Some((0.0, 0.0)),
                "Formulation 8 {property} must retain missing axis {missing_axis}"
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, USPTO_PATENT_TEXT_LICENSE);
        assert!(claim.provenance.source.contains("WO 2018/125520 A1"));
        assert!(claim.provenance.source.contains("US 2018/0179462 A1"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        assert_eq!(id.0, claim.content_hash());
        match composition_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => composition_observation = Some(claim.observations[0]),
        }
    }
    assert!((fraction_sum - 1.0_f64).abs() <= 2.0e-15);

    let composition = decoded
        .claims()
        .observation(composition_observation.expect("composition observation exists"))
        .expect("Formulation 8 composition observation remains linked");
    assert_eq!(
        composition.specimen,
        "wo2018-125520-table-ix-formulation-8-source-products-lots-unspecified"
    );
    assert!(composition.method.contains("Table IX Formulation 8"));
    assert!(composition.caveats.contains("sum to exactly 100.00 wt%"));
    assert!(composition.caveats.contains("not present-day fungible"));
    assert!(composition.caveats.contains("without implying endorsement"));

    let kinematic_viscosity_dims = Dims([2, 0, -1, 0, 0, 0]);
    let viscosity_claims = decoded.claims().claims_for("kinematic_viscosity");
    assert_eq!(viscosity_claims.len(), 2);
    for (source_temperature_c, source_mm2_per_s) in [(40.0, 61.49), (100.0, 10.62)] {
        let temperature_k = source_temperature_c + 273.15;
        let mut matches = viscosity_claims.iter().copied().filter(|(_, claim)| {
            claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
        });
        let (_, claim) = matches.next().unwrap_or_else(|| {
            panic!("missing Formulation 8 viscosity at {source_temperature_c} degC")
        });
        assert!(matches.next().is_none());
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Formulation 8 viscosity was not scalar");
        };
        assert_eq!(*dims, kinematic_viscosity_dims);
        let expected_m2_per_s = source_mm2_per_s * 1.0e-6;
        assert!((*value - expected_m2_per_s).abs() / expected_m2_per_s <= 2.0e-15);
    }

    let dynamic_viscosity_dims = Dims([-1, 1, -1, 0, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let temperature_dims = Dims([0, 0, 0, 1, 0, 0]);
    let expected_unique = [
        ("viscosity_index_scale_reading", 164.0, dimensionless, None),
        (
            "pour_point_temperature",
            -66.0 + 273.15,
            temperature_dims,
            None,
        ),
        (
            "cold_cranking_simulator_dynamic_viscosity",
            4_886.0 * 1.0e-3,
            dynamic_viscosity_dims,
            Some(-30.0 + 273.15),
        ),
        (
            "mini_rotary_viscometer_dynamic_viscosity",
            10_782.0 * 1.0e-3,
            dynamic_viscosity_dims,
            Some(-35.0 + 273.15),
        ),
        (
            "high_temperature_high_shear_dynamic_viscosity",
            3.395 * 1.0e-3,
            dynamic_viscosity_dims,
            Some(150.0 + 273.15),
        ),
        (
            "noack_mass_loss_fraction",
            9.2 * 0.01,
            dimensionless,
            Some(250.0 + 273.15),
        ),
    ];
    let mut performance_observation = None;
    for (property, expected_value, expected_dims, validity_temperature) in expected_unique {
        let expected_value: f64 = expected_value;
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique Formulation 8 {property}");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Formulation 8 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let scale = expected_value.abs().max(1.0e-12);
        assert!((*value - expected_value).abs() / scale <= 2.0e-15);
        match validity_temperature {
            Some(temperature_k) => assert_eq!(
                claim.validity.bound("temperature"),
                Some((temperature_k, temperature_k))
            ),
            None => assert_eq!(claim.validity.bound("temperature"), None),
        }
        assert_eq!(
            claim.validity.bound("source_formulation_number"),
            Some((8.0, 8.0))
        );
        assert_eq!(
            claim.validity.bound("source_viscosity_grade_is_5w30"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_patent_practice_license_granted"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_test_method_edition_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_repeat_count_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, USPTO_PATENT_TEXT_LICENSE);
        assert_eq!(id.0, claim.content_hash());
        match performance_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => performance_observation = Some(claim.observations[0]),
        }
    }
    assert_eq!(
        decoded.claims().claims_for("noack_mass_loss_fraction")[0]
            .1
            .validity
            .bound("source_test_duration_known"),
        Some((0.0, 0.0))
    );

    for (_, claim) in viscosity_claims {
        assert_eq!(
            claim.validity.bound("source_formulation_number"),
            Some((8.0, 8.0))
        );
        assert_eq!(
            claim.validity.bound("source_viscosity_grade_is_5w30"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_patent_practice_license_granted"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, USPTO_PATENT_TEXT_LICENSE);
        match performance_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => performance_observation = Some(claim.observations[0]),
        }
    }

    let performance = decoded
        .claims()
        .observation(performance_observation.expect("performance observation exists"))
        .expect("Formulation 8 performance observation remains linked");
    assert!(performance.method.contains("5W30 property-row"));
    assert!(performance.caveats.contains("absences, not zero-valued"));
    assert!(performance.caveats.contains("do not generalize"));

    assert!(decoded.claims().claims_for("dynamic_viscosity").is_empty());
    assert!(decoded.claims().claims_for("density").is_empty());
    assert!(decoded.claims().claims_for("total_base_number").is_empty());
    assert!(
        decoded
            .claims()
            .claims_for("flash_point_temperature")
            .is_empty()
    );

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_nasa_uam_insulation_stack_constituents() {
    let seeds = [
        (
            NASA_UAM_MW16C_POLYIMIDE_WIRE_SEED_MANIFEST,
            "nasa-uam-mw16c-polyimide-magnet-wire",
            2_usize,
        ),
        (
            NASA_UAM_NOMEX_410_SLOT_LINER_SEED_MANIFEST,
            "nasa-uam-nomex-410-slot-liner",
            1_usize,
        ),
        (
            NASA_UAM_COOLTHERM_EP2000_SEED_MANIFEST,
            "nasa-uam-cooltherm-ep2000-180c-cure",
            2_usize,
        ),
    ];
    let directory = fixture_dir();

    for (manifest_relative, expected_pack_id, expected_claim_count) in seeds {
        let manifest = workspace_path(manifest_relative);
        assert!(
            manifest.is_file(),
            "committed NASA UAM insulation seed manifest is missing: {manifest_relative}"
        );
        let first_path = directory.join(format!("{expected_pack_id}-first.fsmatpk"));
        let second_path = directory.join(format!("{expected_pack_id}-second.fsmatpk"));

        let first = run_compiler(&manifest, &first_path);
        let second = run_compiler(&manifest, &second_path);
        assert!(
            first.status.success(),
            "first {expected_pack_id} compilation failed: {}",
            String::from_utf8_lossy(&first.stderr)
        );
        assert!(
            second.status.success(),
            "second {expected_pack_id} compilation failed: {}",
            String::from_utf8_lossy(&second.stderr)
        );
        assert_eq!(
            first.stdout, second.stdout,
            "{expected_pack_id} decision stream moved"
        );
        assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

        let first_bytes = fs::read(first_path).expect("read first NASA insulation pack");
        let second_bytes = fs::read(second_path).expect("read second NASA insulation pack");
        assert_eq!(
            first_bytes, second_bytes,
            "{expected_pack_id} pack bytes moved"
        );
        let decoded =
            NormalizedPack::from_bytes(&first_bytes).expect("decode NASA insulation pack");
        let pack_hash = decoded.content_hash();
        let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
            .expect("verify NASA insulation pack identity");

        assert_eq!(decoded.pack_id(), expected_pack_id);
        assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
        assert!(
            decoded
                .redistribution_terms()
                .contains("government public use permitted")
        );
        assert_eq!(decoded.claims().claim_count(), expected_claim_count);
        assert!(decoded.joint_statistics().is_empty());

        match expected_pack_id {
            "nasa-uam-mw16c-polyimide-magnet-wire" => {
                let temperature_claims = decoded
                    .claims()
                    .claims_for("thermal_endurance_reference_temperature");
                let duration_claims = decoded
                    .claims()
                    .claims_for("thermal_endurance_reference_duration");
                assert_eq!(temperature_claims.len(), 1);
                assert_eq!(duration_claims.len(), 1);

                let (temperature_id, temperature_claim) = temperature_claims[0];
                let PropertyValue::Scalar {
                    value: temperature,
                    dims: temperature_dims,
                } = &temperature_claim.value
                else {
                    panic!("MW-16C thermal-endurance temperature was not scalar");
                };
                assert_eq!(*temperature_dims, Dims([0, 0, 0, 1, 0, 0]));
                assert!((*temperature - (240.0 + 273.15)).abs() <= 1.0e-12);
                assert_eq!(
                    temperature_claim.validity.bound("reference_duration"),
                    Some((20_000.0 * 3_600.0, 20_000.0 * 3_600.0))
                );

                let (duration_id, duration_claim) = duration_claims[0];
                let PropertyValue::Scalar {
                    value: duration,
                    dims: duration_dims,
                } = &duration_claim.value
                else {
                    panic!("MW-16C thermal-endurance duration was not scalar");
                };
                assert_eq!(*duration_dims, Dims([0, 0, 1, 0, 0, 0]));
                assert_eq!((*duration).to_bits(), (20_000.0_f64 * 3_600.0).to_bits());
                assert_eq!(
                    duration_claim.validity.bound("reference_temperature"),
                    Some((240.0 + 273.15, 240.0 + 273.15))
                );

                for (id, claim) in [
                    (temperature_id, temperature_claim),
                    (duration_id, duration_claim),
                ] {
                    for required_axis in [
                        "source_wire_spec_nema_mw16c",
                        "source_nema_mw1000_2003",
                        "source_test_standard_astm_d2307_2013",
                    ] {
                        assert_eq!(claim.validity.bound(required_axis), Some((1.0, 1.0)));
                    }
                    for missing_axis in [
                        "source_wire_vendor_known",
                        "source_wire_lot_known",
                        "source_thermal_endurance_raw_data_available",
                    ] {
                        assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                    }
                    assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                    assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
                    assert!(claim.provenance.source.contains("NTRS 20240007451"));
                    assert!(claim.provenance.source.contains("[source:primary]"));
                    assert_eq!(id.0, claim.content_hash());
                }
                assert_eq!(temperature_claim.observations, duration_claim.observations);
                let observation = decoded
                    .claims()
                    .observation(temperature_claim.observations[0])
                    .expect("MW-16C observation remains linked");
                assert!(observation.method.contains("ASTM D2307-2013"));
                assert!(
                    observation
                        .caveats
                        .contains("cross-bound classification basis")
                );
                assert!(observation.caveats.contains("not an Arrhenius law"));
                assert!(observation.caveats.contains("different unspecified spool"));
            }
            "nasa-uam-nomex-410-slot-liner" => {
                let claims = decoded.claims().claims_for("selected_slot_liner_thickness");
                assert_eq!(claims.len(), 1);
                let (id, claim) = claims[0];
                let PropertyValue::Scalar { value, dims } = &claim.value else {
                    panic!("Nomex 410 selected thickness was not scalar");
                };
                assert_eq!(*dims, Dims([1, 0, 0, 0, 0, 0]));
                assert!((*value - 0.08e-3).abs() <= 1.0e-15);
                assert_eq!(
                    claim.validity.bound("source_product_is_dupont_nomex_410"),
                    Some((1.0, 1.0))
                );
                assert_eq!(
                    claim.validity.bound("source_material_is_aramid_paper"),
                    Some((1.0, 1.0))
                );
                for missing_axis in [
                    "source_product_lot_known",
                    "source_thickness_measurement_method_known",
                    "source_moisture_condition_known",
                ] {
                    assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                }
                assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
                assert_eq!(id.0, claim.content_hash());
                let observation = decoded
                    .claims()
                    .observation(claim.observations[0])
                    .expect("Nomex 410 observation remains linked");
                assert!(observation.method.contains("component-selection statement"));
                assert!(
                    observation
                        .caveats
                        .contains("not a generic Nomex 410 design allowable")
                );
            }
            "nasa-uam-cooltherm-ep2000-180c-cure" => {
                let completed_claims = decoded
                    .claims()
                    .claims_for("highest_completed_post_cure_temperature");
                let omitted_claims = decoded
                    .claims()
                    .claims_for("manufacturer_recommended_final_cure_temperature");
                assert_eq!(completed_claims.len(), 1);
                assert_eq!(omitted_claims.len(), 1);

                for (claims, source_temperature_c, step_completed) in [
                    (&completed_claims, 180.0, 1.0),
                    (&omitted_claims, 210.0, 0.0),
                ] {
                    let (id, claim) = claims[0];
                    let PropertyValue::Scalar { value, dims } = &claim.value else {
                        panic!("CoolTherm EP-2000 cure temperature was not scalar");
                    };
                    assert_eq!(*dims, Dims([0, 0, 0, 1, 0, 0]));
                    assert!((*value - (source_temperature_c + 273.15)).abs() <= 1.0e-12);
                    assert_eq!(
                        claim
                            .validity
                            .bound("source_product_is_parker_lord_cooltherm_ep2000"),
                        Some((1.0, 1.0))
                    );
                    assert_eq!(
                        claim.validity.bound("source_step_completed"),
                        Some((step_completed, step_completed))
                    );
                    for missing_axis in [
                        "source_epoxy_lot_known",
                        "source_cure_hold_duration_known",
                        "source_degree_of_cure_known",
                    ] {
                        assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                    }
                    assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                    assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
                    assert_eq!(id.0, claim.content_hash());
                }
                assert_eq!(
                    completed_claims[0].1.observations,
                    omitted_claims[0].1.observations
                );
                let observation = decoded
                    .claims()
                    .observation(completed_claims[0].1.observations[0])
                    .expect("CoolTherm EP-2000 observation remains linked");
                assert!(observation.caveats.contains("intentionally not completed"));
                assert!(
                    observation
                        .caveats
                        .contains("deliberately incomplete source process")
                );
            }
            unexpected => panic!("unexpected NASA insulation pack {unexpected}"),
        }

        for refused_property in [
            "partial_discharge_inception_voltage",
            "dielectric_strength",
            "thermal_conductivity",
            "service_life",
            "arrhenius_activation_energy",
        ] {
            assert!(
                decoded.claims().claims_for(refused_property).is_empty(),
                "assembly- or source-absent property must stay refused for {expected_pack_id}: {refused_property}"
            );
        }

        let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
        assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
        assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
        assert!(
            decisions
                .lines()
                .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
        );
    }
}

#[test]
fn g3_cli_compiles_committed_aisi_4140_rc33_exact_condition_seed() {
    let manifest = workspace_path(AISI_4140_RC33_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 4140 Rockwell C33 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-4140-rc33-first.fsmatpk");
    let second_path = directory.join("aisi-4140-rc33-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 4140 Rockwell C33 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 4140 Rockwell C33 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 4140 Rockwell C33 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 4140 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 4140 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 4140 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 4140 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 4140 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-4140-qq-s-624-rc33");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 14);
    assert!(decoded.joint_statistics().is_empty());

    let pressure_dims = Dims([-1, 1, -2, 0, 0, 0]);
    let energy_dims = Dims([2, 1, -2, 0, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let expected = [
        (
            "ultimate_tensile_strength",
            26.7,
            1.074 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "yield_strength_0p2_offset",
            26.7,
            0.985 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_elongation_2in",
            26.7,
            19.4 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_reduction_of_area",
            26.7,
            62.5 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "charpy_v_notch_impact_energy",
            26.7,
            95.2,
            energy_dims,
            "MIL-STD-151 Charpy V-notched impact",
            "four impact tests",
        ),
        (
            "double_shear_ultimate_strength",
            26.7,
            0.66 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
        (
            "double_shear_yield_strength",
            26.7,
            0.56 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
        (
            "ultimate_tensile_strength",
            -73.0,
            1.158 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "yield_strength_0p2_offset",
            -73.0,
            1.060 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_elongation_2in",
            -73.0,
            20.0 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_reduction_of_area",
            -73.0,
            61.0 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "charpy_v_notch_impact_energy",
            -73.0,
            84.6,
            energy_dims,
            "MIL-STD-151 Charpy V-notched impact",
            "four impact tests",
        ),
        (
            "double_shear_ultimate_strength",
            -73.0,
            0.73 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
        (
            "double_shear_yield_strength",
            -73.0,
            0.60 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
    ];

    for (property, temperature_c, expected_value, expected_dims, method_note, sample_note) in
        expected
    {
        let expected_value: f64 = expected_value;
        let temperature_k = temperature_c + 273.15;
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
            })
            .unwrap_or_else(|| {
                panic!("missing AISI 4140 {property} claim at {temperature_c} degC")
            });
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 4140 {property} at {temperature_c} degC was not scalar");
        };
        assert_eq!(
            *dims, expected_dims,
            "AISI 4140 {property} dimensions moved"
        );
        let scale = f64::abs(expected_value).max(1.0);
        let relative_error = (*value - expected_value).abs() / scale;
        assert!(
            relative_error <= 2.0e-15,
            "AISI 4140 {property} at {temperature_c} degC moved by {relative_error:e} relative"
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TM-X-64791"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 4140 claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-4140-QQ-S-624-heat-137M186-1in-bar-Rockwell-C33"
        );
        assert!(observation.method.contains(method_note));
        assert!(observation.caveats.contains(sample_note));
        assert!(observation.caveats.contains("oil quenched"));
        assert!(observation.caveats.contains("tempered 566 degC"));
    }

    // NASA Table IV prints both ksi and GN/m2. This checks transcription and
    // unit normalization against the redundant source columns; it is not an
    // independent-source agreement claim.
    for (property, temperature_c, source_ksi) in [
        ("ultimate_tensile_strength", 26.7, 155.8),
        ("yield_strength_0p2_offset", 26.7, 142.9),
        ("ultimate_tensile_strength", -73.0, 168.0),
        ("yield_strength_0p2_offset", -73.0, 153.7),
    ] {
        let temperature_k = temperature_c + 273.15;
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
            })
            .expect("AISI 4140 redundant-unit comparison point");
        let PropertyValue::Scalar { value, .. } = &claim.value else {
            panic!("AISI 4140 redundant-unit comparison point was not scalar");
        };
        let source_ksi_in_pa = source_ksi * 6_894_757.293_168_361;
        let relative_rounding_difference = (*value - source_ksi_in_pa).abs() / *value;
        assert!(
            relative_rounding_difference <= 5.0e-4,
            "AISI 4140 {property} at {temperature_c} degC disagrees with the source ksi column by {relative_rounding_difference:e}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_aisi_1045_cold_drawn_tensile_seed() {
    let manifest = workspace_path(AISI_1045_COLD_DRAWN_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 1045 cold-drawn seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-1045-cold-drawn-first.fsmatpk");
    let second_path = directory.join("aisi-1045-cold-drawn-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 1045 cold-drawn seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 1045 cold-drawn seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 1045 cold-drawn decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 1045 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 1045 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 1045 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 1045 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 1045 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-1045-cold-drawn-tensile");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("Attribution 4.0 International")
    );
    assert_eq!(decoded.claims().claim_count(), 3);
    assert!(
        decoded.joint_statistics().is_empty(),
        "paired source rows do not authorize an inferred covariance block"
    );

    let pressure_dims = Dims([-1, 1, -2, 0, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let expected = [
        (
            "yield_strength",
            550.51,
            0.005,
            [540.73, 557.59, 553.20],
            1.0e6,
            pressure_dims,
        ),
        (
            "ultimate_tensile_strength",
            695.31,
            0.005,
            [684.58, 707.75, 693.60],
            1.0e6,
            pressure_dims,
        ),
        (
            "tensile_elongation_50mm",
            14.1,
            0.05,
            [14.42, 14.20, 13.68],
            0.01,
            dimensionless,
        ),
    ];
    let student_t_0p975_df2 = 4.302_652_729_911_275;
    let crosshead_speed_m_per_s = 10.0 * 1.0e-3 / 60.0;

    for (
        property,
        reported_mean,
        reported_rounding_half_width,
        samples,
        source_unit_scale,
        expected_dims,
    ) in expected
    {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(
            claims.len(),
            1,
            "expected exactly one AISI 1045 {property} claim"
        );
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 1045 {property} was not scalar");
        };
        assert_eq!(
            *dims, expected_dims,
            "AISI 1045 {property} dimensions moved"
        );
        let expected_value = reported_mean * source_unit_scale;
        let relative_value_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
        assert!(
            relative_value_error <= 2.0e-15,
            "AISI 1045 {property} moved by {relative_value_error:e} relative"
        );

        let sample_mean = samples.iter().copied().sum::<f64>() / 3.0;
        assert!(
            (sample_mean - reported_mean).abs() <= reported_rounding_half_width,
            "AISI 1045 {property} source mean is inconsistent with its printed replicates"
        );
        let sample_variance = samples
            .iter()
            .map(|sample| (sample - sample_mean).powi(2))
            .sum::<f64>()
            / 2.0;
        let expected_half_width =
            student_t_0p975_df2 * sample_variance.sqrt() / 3.0_f64.sqrt() * source_unit_scale;
        let UncertaintyModel::HalfWidth {
            half_width,
            confidence,
        } = &claim.uncertainty
        else {
            panic!("AISI 1045 {property} lost its derived Student-t half-width");
        };
        let relative_half_width_error =
            (*half_width - expected_half_width).abs() / expected_half_width;
        assert!(
            relative_half_width_error <= 2.0e-14,
            "AISI 1045 {property} Student-t half-width moved by {relative_half_width_error:e} relative"
        );
        assert_eq!((*confidence).to_bits(), 0.95f64.to_bits());

        let Some((speed_lo, speed_hi)) = claim.validity.bound("crosshead_speed") else {
            panic!("AISI 1045 {property} lost its crosshead-speed validity point");
        };
        for speed in [speed_lo, speed_hi] {
            let relative_speed_error =
                (speed - crosshead_speed_m_per_s).abs() / crosshead_speed_m_per_s;
            assert!(
                relative_speed_error <= 2.0e-15,
                "AISI 1045 {property} crosshead speed moved by {relative_speed_error:e} relative"
            );
        }
        assert_eq!(
            claim.validity.bound("source_test_temperature_known"),
            Some((0.0, 0.0)),
            "AISI 1045 {property} must require explicit acknowledgement of missing temperature"
        );
        assert_eq!(claim.validity.bounds().len(), 2);
        assert_eq!(claim.validity.bound("temperature"), None);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(claim.provenance.source.contains("doi:10.3390/pr12061171"));
        assert!(claim.provenance.source.contains("[source:primary]"));

        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 1045 claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-1045-cold-drawn-bar-37mm-OD-102mm-length-test-temperature-not-reported"
        );
        assert!(observation.method.contains("ASTM E8"));
        assert!(observation.method.contains("50 mm gauge length"));
        assert!(observation.method.contains("10 mm/min crosshead speed"));
        assert!(observation.caveats.contains("three specimens"));
        assert!(observation.caveats.contains("t(0.975, df=2)"));
        assert!(
            observation
                .caveats
                .contains("source does not report test temperature")
        );
        assert!(
            observation
                .caveats
                .contains("no joint covariance is inferred")
        );
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_aisi_52100_cvm_heat_treatment_states() {
    let manifest = workspace_path(AISI_52100_CVM_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 52100 CVM seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-52100-cvm-first.fsmatpk");
    let second_path = directory.join("aisi-52100-cvm-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 52100 CVM seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 52100 CVM seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 52100 CVM decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 52100 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 52100 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 52100 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 52100 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 52100 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-52100-cvm-nasa-tn-d-6632");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 15);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    for (property, source_percent) in [
        ("carbon_mass_fraction", 0.96),
        ("silicon_mass_fraction", 0.22),
        ("manganese_mass_fraction", 0.36),
        ("sulfur_mass_fraction", 0.012),
        ("phosphorus_mass_fraction", 0.007),
        ("chromium_mass_fraction", 1.36),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one AISI 52100 {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 52100 {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "AISI 52100 {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TN-D-6632"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 52100 chemistry observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-52100-consumable-vacuum-melted-single-ingot-NASA-TN-D-6632"
        );
        assert!(observation.method.contains("Table I actual composition"));
        assert!(observation.caveats.contains("Balance iron"));
        assert!(observation.caveats.contains("no heat identifier"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let hardness_claims = decoded.claims().claims_for("rockwell_c_scale_reading");
    let austenite_claims = decoded
        .claims()
        .claims_for("retained_austenite_volume_fraction");
    assert_eq!(hardness_claims.len(), 5);
    assert_eq!(austenite_claims.len(), 4);
    let states = [
        (Some(505.0), "second-temper-505K", 59.7, None),
        (Some(450.0), "second-temper-450K", 62.3, Some(12.8)),
        (Some(433.0), "second-temper-433K", 63.4, Some(15.6)),
        (Some(394.0), "second-temper-394K", 64.6, Some(18.4)),
        (None, "no-second-temper", 65.1, Some(11.8)),
    ];

    for (second_temper_k, specimen_state, expected_hardness, expected_austenite_percent) in states {
        let matches_state = |claim: &fs_matdb::PropertyClaim| match second_temper_k {
            Some(second_temper_k) => {
                claim.validity.bound("second_temper_temperature")
                    == Some((second_temper_k, second_temper_k))
                    && claim.validity.bound("second_temper_applied").is_none()
            }
            None => {
                claim.validity.bound("second_temper_applied") == Some((0.0, 0.0))
                    && claim.validity.bound("second_temper_temperature").is_none()
            }
        };
        let (hardness_id, hardness_claim) = hardness_claims
            .iter()
            .copied()
            .find(|(_, claim)| matches_state(claim))
            .unwrap_or_else(|| panic!("missing AISI 52100 hardness state {specimen_state}"));
        let PropertyValue::Scalar { value, dims } = &hardness_claim.value else {
            panic!("AISI 52100 hardness state {specimen_state} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        assert_eq!((*value).to_bits(), f64::to_bits(expected_hardness));
        assert_eq!(
            hardness_claim.validity.bound("temperature"),
            Some((294.0, 294.0))
        );
        assert_eq!(hardness_claim.validity.bounds().len(), 2);
        assert_eq!(hardness_claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            hardness_claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(hardness_claim.observations.len(), 1);
        assert_eq!(hardness_claim.provenance.license, NASA_SEED_LICENSE);
        let hardness_observation = decoded
            .claims()
            .observation(hardness_claim.observations[0])
            .expect("AISI 52100 hardness observation remains linked");
        assert!(hardness_observation.specimen.contains(specimen_state));
        assert!(
            hardness_observation
                .specimen
                .contains("austenitize-1116-to-1144K-30min")
        );
        assert!(hardness_observation.specimen.contains("oil-quench-325K"));
        assert!(
            hardness_observation
                .specimen
                .contains("first-temper-394K-60min")
        );
        assert!(hardness_observation.method.contains("150 kg load"));
        assert!(hardness_observation.method.contains("Rockwell C"));
        assert!(hardness_observation.method.contains("294 K reading"));
        assert!(
            hardness_observation
                .caveats
                .contains("Minimum two hardness measurements")
        );
        assert!(
            hardness_observation
                .caveats
                .contains("dispersion not reported")
        );
        assert!(hardness_observation.caveats.contains("ASTM grain size 12"));
        assert!(
            hardness_observation
                .caveats
                .contains("predictive equation is not measurement uncertainty")
        );
        if second_temper_k == Some(505.0) {
            assert!(
                hardness_observation
                    .caveats
                    .contains("censored as less than 2 volume percent")
            );
        }
        assert_eq!(
            hardness_claim.observations[0].0,
            hardness_observation.content_hash()
        );
        assert_eq!(hardness_id.0, hardness_claim.content_hash());

        let matching_austenite = austenite_claims
            .iter()
            .copied()
            .find(|(_, claim)| matches_state(claim));
        match (matching_austenite, expected_austenite_percent) {
            (Some((austenite_id, austenite_claim)), Some(expected_percent)) => {
                let PropertyValue::Scalar { value, dims } = &austenite_claim.value else {
                    panic!("AISI 52100 retained-austenite state {specimen_state} was not scalar");
                };
                assert_eq!(*dims, dimensionless);
                let expected_value = expected_percent * 0.01;
                let relative_error = (*value - expected_value).abs() / expected_value;
                assert!(
                    relative_error <= 2.0e-15,
                    "AISI 52100 retained austenite {specimen_state} moved by {relative_error:e} relative"
                );
                assert_eq!(
                    austenite_claim.validity.bound("temperature"),
                    Some((294.0, 294.0))
                );
                assert_eq!(austenite_claim.validity.bounds().len(), 2);
                assert_eq!(austenite_claim.uncertainty, UncertaintyModel::Unstated);
                assert_eq!(austenite_claim.observations.len(), 1);
                let observation = decoded
                    .claims()
                    .observation(austenite_claim.observations[0])
                    .expect("AISI 52100 retained-austenite observation remains linked");
                assert!(observation.specimen.contains(specimen_state));
                assert!(observation.method.contains("X-ray diffraction"));
                assert!(
                    observation
                        .caveats
                        .contains("uncertainty and replicate count not reported")
                );
                assert!(
                    observation
                        .caveats
                        .contains("no covariance with hardness is inferred")
                );
                assert_ne!(
                    austenite_claim.observations[0],
                    hardness_claim.observations[0]
                );
                assert_eq!(
                    austenite_claim.observations[0].0,
                    observation.content_hash()
                );
                assert_eq!(austenite_id.0, austenite_claim.content_hash());
            }
            (None, None) => {}
            (Some(_), None) => panic!("censored AISI 52100 austenite became an exact scalar"),
            (None, Some(_)) => panic!("missing exact AISI 52100 austenite state {specimen_state}"),
        }
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_aisi_9310_cvm_carburized_gear_seed() {
    let manifest = workspace_path(AISI_9310_CVM_CARBURIZED_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 9310 CVM carburized seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-9310-cvm-carburized-first.fsmatpk");
    let second_path = directory.join("aisi-9310-cvm-carburized-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 9310 CVM carburized seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 9310 CVM carburized seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 9310 CVM carburized decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 9310 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 9310 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 9310 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 9310 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 9310 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-9310-cvm-carburized-nasa-tm-104352");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 13);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    for (property, source_percent) in [
        ("carbon_mass_fraction", 0.10),
        ("manganese_mass_fraction", 0.63),
        ("silicon_mass_fraction", 0.27),
        ("nickel_mass_fraction", 3.22),
        ("chromium_mass_fraction", 1.21),
        ("molybdenum_mass_fraction", 0.12),
        ("copper_mass_fraction", 0.13),
        ("phosphorus_mass_fraction", 0.005),
        ("sulfur_mass_fraction", 0.005),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one AISI 9310 {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 9310 {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "AISI 9310 {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TM-104352"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 9310 chemistry observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-9310-CVM-single-lot-single-heat-28-tooth-spur-gear-NASA-TM-104352"
        );
        assert!(observation.method.contains("Table I nominal composition"));
        assert!(observation.caveats.contains("Nominal grade chemistry"));
        assert!(observation.caveats.contains("balance iron"));
        assert!(observation.caveats.contains("not inferred"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let case_claims = decoded.claims().claims_for("case_rockwell_c_scale_reading");
    assert_eq!(
        case_claims.len(),
        2,
        "the report's conflicting C58 and C60 case statements must both survive"
    );
    let mut case_values = Vec::with_capacity(2);
    for (id, claim) in case_claims {
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 9310 case hardness was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 9310 case-hardness observation remains linked");
        match *value {
            58.0 => {
                assert!(observation.method.contains("Test Materials detailed"));
                assert!(observation.caveats.contains("carburize 1172 K for 8 h"));
                assert!(observation.caveats.contains("austenitize 1117 K for 2.5 h"));
                assert!(
                    observation
                        .caveats
                        .contains("subzero treat 180 K for 3.5 h")
                );
                assert!(observation.caveats.contains("double temper 450 K"));
                assert!(observation.caveats.contains("stress relieve 450 K for 2 h"));
                assert!(observation.caveats.contains("conflicts"));
            }
            60.0 => {
                assert!(observation.method.contains("abstract and summary"));
                assert!(observation.caveats.contains("conflicts"));
                assert!(observation.caveats.contains("not averaged or selected"));
            }
            other => panic!("unexpected AISI 9310 case-hardness value {other}"),
        }
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
        case_values.push(*value);
    }
    case_values.sort_by(f64::total_cmp);
    assert_eq!(case_values, [58.0, 60.0]);

    let core_claims = decoded.claims().claims_for("core_rockwell_c_scale_reading");
    assert_eq!(core_claims.len(), 1);
    let (core_id, core_claim) = core_claims[0];
    let PropertyValue::Scalar {
        value: core_value,
        dims: core_dims,
    } = &core_claim.value
    else {
        panic!("AISI 9310 core hardness was not scalar");
    };
    assert_eq!(*core_dims, dimensionless);
    assert_eq!((*core_value).to_bits(), 40.0f64.to_bits());
    assert!(core_claim.validity.bounds().is_empty());
    assert_eq!(core_claim.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(core_claim.observations.len(), 1);

    let depth_claims = decoded.claims().claims_for("carburized_case_depth");
    assert_eq!(depth_claims.len(), 1);
    let (depth_id, depth_claim) = depth_claims[0];
    let PropertyValue::Scalar {
        value: depth_value,
        dims: depth_dims,
    } = &depth_claim.value
    else {
        panic!("AISI 9310 carburized case depth was not scalar");
    };
    assert_eq!(*depth_dims, Dims([1, 0, 0, 0, 0, 0]));
    let expected_depth_m = 0.97e-3;
    let relative_depth_error = (*depth_value - expected_depth_m).abs() / expected_depth_m;
    assert!(relative_depth_error <= 2.0e-15);
    assert!(depth_claim.validity.bounds().is_empty());
    assert_eq!(depth_claim.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(depth_claim.observations.len(), 1);
    assert_eq!(core_claim.observations, depth_claim.observations);

    let detailed_observation = decoded
        .claims()
        .observation(core_claim.observations[0])
        .expect("AISI 9310 detailed gear observation remains linked");
    assert!(detailed_observation.method.contains("case/core hardness"));
    assert!(detailed_observation.method.contains("case-depth"));
    assert!(
        detailed_observation
            .caveats
            .contains("One lot from one CVM heat")
    );
    assert!(detailed_observation.caveats.contains("replicate count"));
    assert_eq!(
        core_claim.observations[0].0,
        detailed_observation.content_hash()
    );
    assert_eq!(core_id.0, core_claim.content_hash());
    assert_eq!(depth_id.0, depth_claim.content_hash());

    // G3 plausibility only: NASA SP-410 (NTRS 19750018303) reports a different
    // VAR AISI 9310 gear lot at nominal C62 case, C45 core, and 1 mm case depth.
    // These checks bound transcription-scale agreement without fusing the lots.
    let independent_case_hardness: f64 = 62.0;
    let independent_core_hardness: f64 = 45.0;
    let independent_case_depth_m: f64 = 1.0e-3;
    assert!((58.0 - independent_case_hardness).abs() <= 4.0);
    assert!((*core_value - independent_core_hardness).abs() <= 5.0);
    assert!((*depth_value - independent_case_depth_m).abs() / independent_case_depth_m <= 0.031);

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_napc_gear_oils_without_fusing_batches() {
    type ExpectedClaim = (&'static str, &'static str, f64, Dims, Option<f64>);

    let temperature_dims = Dims([0, 0, 0, 1, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let seeds: [(&str, &str, &str, &[ExpectedClaim]); 2] = [
        (
            NAPC_PE_5_L_1274_GEAR_OIL_SEED_MANIFEST,
            "napc-pe-5-l-1274",
            "napc-pe-5-l-1274-polyol-ester-gear-oil",
            &[
                (
                    "flash_point_temperature",
                    "PE-5-L-1274",
                    516.0,
                    temperature_dims,
                    None,
                ),
                (
                    "reported_specific_gravity",
                    "PE-5-L-1274",
                    0.998,
                    dimensionless,
                    Some(298.0),
                ),
                (
                    "total_acid_number_as_koh_mass_per_oil_mass",
                    "PE-5-L-1274",
                    0.07e-3,
                    dimensionless,
                    None,
                ),
            ],
        ),
        (
            NAPC_PE_5_L_1307_1553_GEAR_OIL_SEED_MANIFEST,
            "napc-pe-5-l-1307-1553",
            "napc-pe-5-l-1307-1553-mil-l-23699-gear-oil",
            &[
                (
                    "flash_point_temperature",
                    "PE-5-L-1307-NASA",
                    539.0,
                    temperature_dims,
                    None,
                ),
                (
                    "pour_point_temperature",
                    "PE-5-L-1307-NASA",
                    220.0,
                    temperature_dims,
                    None,
                ),
                (
                    "flash_point_temperature",
                    "PE-5-L-1553-NASA",
                    539.0,
                    temperature_dims,
                    None,
                ),
                (
                    "pour_point_temperature",
                    "PE-5-L-1553-NASA",
                    213.0,
                    temperature_dims,
                    None,
                ),
                (
                    "reported_specific_gravity",
                    "PE-5-L-1307-and-PE-5-L-1553",
                    1.0,
                    dimensionless,
                    Some(289.0),
                ),
                (
                    "total_acid_number_as_koh_mass_per_oil_mass",
                    "PE-5-L-1307-and-PE-5-L-1553",
                    0.03e-3,
                    dimensionless,
                    None,
                ),
            ],
        ),
    ];
    let directory = fixture_dir();

    for (manifest_relative, stem, expected_pack_id, expected_rows) in seeds {
        let manifest = workspace_path(manifest_relative);
        assert!(
            manifest.is_file(),
            "committed NASA/NAPC gear-oil seed manifest is missing: {manifest_relative}"
        );
        let first_path = directory.join(format!("{stem}-first.fsmatpk"));
        let second_path = directory.join(format!("{stem}-second.fsmatpk"));

        let first = run_compiler(&manifest, &first_path);
        let second = run_compiler(&manifest, &second_path);
        assert!(
            first.status.success(),
            "first {stem} seed compilation failed: {}",
            String::from_utf8_lossy(&first.stderr)
        );
        assert!(
            second.status.success(),
            "second {stem} seed compilation failed: {}",
            String::from_utf8_lossy(&second.stderr)
        );
        assert_eq!(first.stdout, second.stdout, "{stem} decision stream moved");
        assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

        let first_bytes = fs::read(first_path).expect("read first NASA/NAPC gear-oil pack");
        let second_bytes = fs::read(second_path).expect("read second NASA/NAPC gear-oil pack");
        assert_eq!(first_bytes, second_bytes, "{stem} pack bytes moved");
        let decoded =
            NormalizedPack::from_bytes(&first_bytes).expect("decode NASA/NAPC gear-oil pack");
        let pack_hash = decoded.content_hash();
        let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
            .expect("verify NASA/NAPC gear-oil pack identity");

        assert_eq!(decoded.pack_id(), expected_pack_id);
        assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
        assert!(
            decoded
                .redistribution_terms()
                .contains("public use is permitted")
        );
        assert_eq!(decoded.claims().claim_count(), expected_rows.len());
        assert!(decoded.joint_statistics().is_empty());

        assert!(
            decoded
                .claims()
                .claims_for("kinematic_viscosity")
                .is_empty(),
            "Table IV omits the viscosity unit, so no viscosity claim is admissible"
        );
        for &(property, observation_token, expected_value, expected_dims, validity_temperature) in
            expected_rows
        {
            let claims = decoded.claims().claims_for(property);
            let mut matches = claims.iter().copied().filter(|(_, claim)| {
                if claim.observations.len() != 1 {
                    return false;
                }
                decoded
                    .claims()
                    .observation(claim.observations[0])
                    .is_some_and(|observation| observation.specimen.contains(observation_token))
            });
            let (id, claim) = matches
                .next()
                .unwrap_or_else(|| panic!("missing {property} for {observation_token}"));
            assert!(
                matches.next().is_none(),
                "duplicate {property} for {observation_token}"
            );
            let PropertyValue::Scalar { value, dims } = &claim.value else {
                panic!("{property} for {observation_token} was not scalar");
            };
            assert_eq!(*dims, expected_dims);
            let relative_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
            assert!(
                relative_error <= 2.0e-15,
                "{property} for {observation_token} moved by {relative_error:e} relative"
            );
            match validity_temperature {
                Some(temperature_k) => {
                    assert_eq!(
                        claim.validity.bound("temperature"),
                        Some((temperature_k, temperature_k))
                    );
                    assert_eq!(claim.validity.bounds().len(), 1);
                }
                None => assert!(claim.validity.bounds().is_empty()),
            }
            assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
            assert_eq!(
                claim.interpolation,
                InterpolationPolicy::ConstantWithinValidity
            );
            assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
            assert!(claim.provenance.source.contains("NASA-TM-104352"));
            assert!(claim.provenance.source.contains("[source:primary]"));

            let observation = decoded
                .claims()
                .observation(claim.observations[0])
                .expect("NASA/NAPC gear-oil observation remains linked");
            assert!(observation.method.contains("Table IV"));
            assert!(observation.caveats.contains("proprietary"));
            if observation_token == "PE-5-L-1274" {
                assert!(observation.caveats.contains("NASA reference lubricant"));
                assert!(observation.caveats.contains("without stating their unit"));
                assert!(observation.caveats.contains("less than 200 K"));
            } else if observation_token.contains("-NASA") {
                assert!(
                    observation
                        .caveats
                        .contains("two batches of the same lubricant")
                );
                assert!(observation.caveats.contains("MIL-L-23699"));
                assert!(
                    observation
                        .caveats
                        .contains("batch-specific values remain separate")
                );
                assert!(observation.caveats.contains("without stating their unit"));
            } else {
                assert!(observation.caveats.contains("MIL-L-23699"));
                assert!(observation.caveats.contains("same MIL-L-23699 lubricant"));
            }
            assert_eq!(claim.observations[0].0, observation.content_hash());
            assert_eq!(id.0, claim.content_hash());
        }

        if expected_rows.len() == 6 {
            for property in ["flash_point_temperature", "pour_point_temperature"] {
                let batch_claims = decoded.claims().claims_for(property);
                assert_eq!(batch_claims.len(), 2);
                assert_ne!(
                    batch_claims[0].1.observations, batch_claims[1].1.observations,
                    "the two NASA/NAPC batches were fused for {property}"
                );
            }
        } else {
            assert!(
                decoded
                    .claims()
                    .claims_for("pour_point_temperature")
                    .is_empty(),
                "the censored PE-5-L-1274 pour point became an exact scalar"
            );
        }

        let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
        assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
        assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
        assert!(
            decisions
                .lines()
                .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
        );
    }
}

#[test]
fn g3_cli_compiles_committed_rheolube_2000_bearing_grease_seed() {
    let manifest = workspace_path(RHEOLUBE_2000_PENNZANE_GREASE_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Rheolube 2000 bearing-grease seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("rheolube-2000-first.fsmatpk");
    let second_path = directory.join("rheolube-2000-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Rheolube 2000 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Rheolube 2000 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Rheolube 2000 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Rheolube 2000 pack");
    let second_bytes = fs::read(second_path).expect("read second Rheolube 2000 pack");
    assert_eq!(first_bytes, second_bytes, "Rheolube 2000 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode Rheolube 2000 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Rheolube 2000 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "rheolube-2000-pennzane-shf-x-2000-bearing-grease"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 3);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let expected = [
        ("nlgi_consistency_grade", 2.0, dimensionless, None, None),
        ("density", 890.0, density_dims, Some(298.15), None),
        (
            "oil_separation_mass_fraction",
            0.033,
            dimensionless,
            Some(373.15),
            Some(86_400.0),
        ),
    ];
    let mut shared_observation = None;

    for (property, expected_value, expected_dims, temperature_k, duration_s) in expected {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(
            claims.len(),
            1,
            "expected one Rheolube 2000 {property} claim"
        );
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Rheolube 2000 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let relative_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
        assert!(
            relative_error <= 2.0e-15,
            "Rheolube 2000 {property} moved by {relative_error:e} relative"
        );
        match temperature_k {
            Some(temperature_k) => assert_eq!(
                claim.validity.bound("temperature"),
                Some((temperature_k, temperature_k))
            ),
            None => assert_eq!(claim.validity.bound("temperature"), None),
        }
        match duration_s {
            Some(duration_s) => assert_eq!(
                claim.validity.bound("duration"),
                Some((duration_s, duration_s))
            ),
            None => assert_eq!(claim.validity.bound("duration"), None),
        }
        assert_eq!(
            claim.validity.bounds().len(),
            usize::from(temperature_k.is_some()) + usize::from(duration_s.is_some())
        );
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-CP-3350"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        assert_eq!(claim.observations.len(), 1);
        match shared_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => shared_observation = Some(claim.observations[0]),
        }

        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("Rheolube 2000 observation remains linked");
        assert_eq!(
            observation.specimen,
            "Rheolube-2000-Pennzane-SHF-X-2000-sodium-octadecylterephthalamate-bearing-grease"
        );
        assert!(observation.method.contains("Bessette Table 7"));
        assert!(observation.method.contains("typical grease properties"));
        assert!(observation.caveats.contains("approximately 20 percent"));
        assert!(observation.caveats.contains("approximately 260 degC"));
        assert!(observation.caveats.contains("labels results typical"));
        assert!(observation.caveats.contains("vacuum-hardening state"));
        assert!(observation.caveats.contains("no printed unit or method"));
        assert!(observation.caveats.contains("not admitted as bulk claims"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    for refused_property in [
        "dropping_point_temperature",
        "penetration_scale_reading",
        "wear_scar_diameter",
        "oxidation_pressure_drop",
        "vapor_pressure",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "Rheolube 2000 {refused_property} crossed the no-claim boundary"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_pennzane_shf_x_2000_bearing_oil_seed() {
    let manifest = workspace_path(PENNZANE_SHF_X_2000_BEARING_OIL_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Pennzane SHF X-2000 bearing-oil seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("pennzane-shf-x-2000-first.fsmatpk");
    let second_path = directory.join("pennzane-shf-x-2000-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Pennzane SHF X-2000 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Pennzane SHF X-2000 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Pennzane SHF X-2000 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Pennzane SHF X-2000 pack");
    let second_bytes = fs::read(second_path).expect("read second Pennzane SHF X-2000 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Pennzane SHF X-2000 pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode Pennzane SHF X-2000 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Pennzane SHF X-2000 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "pennzane-shf-x-2000-mac-aerospace-bearing-oil"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 7);
    assert!(decoded.joint_statistics().is_empty());

    let kinematic_viscosity_dims = Dims([2, 0, -1, 0, 0, 0]);
    let viscosity_claims = decoded.claims().claims_for("kinematic_viscosity");
    assert_eq!(viscosity_claims.len(), 3);
    let mut shared_observation = None;
    for (source_temperature_c, source_mm2_per_s) in
        [(100.0, 14.3), (40.0, 107.0), (-40.0, 80_500.0)]
    {
        let temperature_k = source_temperature_c + 273.15;
        let mut matches = viscosity_claims.iter().copied().filter(|(_, claim)| {
            claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
        });
        let (id, claim) = matches
            .next()
            .unwrap_or_else(|| panic!("missing Pennzane viscosity at {source_temperature_c} degC"));
        assert!(matches.next().is_none());
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Pennzane viscosity at {source_temperature_c} degC was not scalar");
        };
        assert_eq!(*dims, kinematic_viscosity_dims);
        let expected_m2_per_s = source_mm2_per_s * 1.0e-6;
        let relative_error = (*value - expected_m2_per_s).abs() / expected_m2_per_s;
        assert!(relative_error <= 2.0e-15);
        assert_eq!(claim.validity.bounds().len(), 1);
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-CP-3350"));
        assert_eq!(claim.observations.len(), 1);
        match shared_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => shared_observation = Some(claim.observations[0]),
        }
        assert_eq!(id.0, claim.content_hash());
    }

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let temperature_dims = Dims([0, 0, 0, 1, 0, 0]);
    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let expected = [
        ("viscosity_index_scale_reading", 137.0, dimensionless, None),
        (
            "flash_point_temperature",
            300.0 + 273.15,
            temperature_dims,
            None,
        ),
        (
            "pour_point_temperature",
            -55.0 + 273.15,
            temperature_dims,
            None,
        ),
        ("density", 0.84 * 1_000.0, density_dims, Some(298.15)),
    ];

    for (property, expected_value, expected_dims, validity_temperature) in expected {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one Pennzane {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Pennzane {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let relative_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
        assert!(relative_error <= 2.0e-15);
        match validity_temperature {
            Some(temperature_k) => {
                assert_eq!(
                    claim.validity.bound("temperature"),
                    Some((temperature_k, temperature_k))
                );
                assert_eq!(claim.validity.bounds().len(), 1);
            }
            None => assert!(claim.validity.bounds().is_empty()),
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(
            claim.observations[0],
            shared_observation.expect("shared observation")
        );
        assert_eq!(id.0, claim.content_hash());
    }

    let observation_id = shared_observation.expect("Pennzane observation id");
    let observation = decoded
        .claims()
        .observation(observation_id)
        .expect("Pennzane observation remains linked");
    assert_eq!(
        observation.specimen,
        "Pennzane-SHF-X-2000-multiply-alkylated-cyclopentane-aerospace-bearing-oil"
    );
    assert!(observation.method.contains("Bessette Table 6"));
    assert!(observation.method.contains("typical Pennzane properties"));
    assert!(
        observation
            .caveats
            .contains("Tris(2-octyldodecyl) cyclopentane")
    );
    assert!(
        observation
            .caveats
            .contains("approximate molecular weight 910 g/mol")
    );
    assert!(observation.caveats.contains("labels results typical"));
    assert!(
        observation
            .caveats
            .contains("no temperature-interval degree-Celsius token")
    );
    assert!(observation.caveats.contains("tribometer-system result"));
    assert!(observation.caveats.contains("not unambiguous"));
    assert_eq!(observation_id.0, observation.content_hash());

    for refused_property in [
        "volumetric_thermal_expansion_coefficient",
        "wear_scar_diameter",
        "vapor_pressure",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "Pennzane {refused_property} crossed the no-claim boundary"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_s2_s_gray_cast_iron_seed() {
    let manifest = workspace_path(GRAY_CAST_IRON_S2_S_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed S2-S gray-cast-iron seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("gray-cast-iron-s2-s-first.fsmatpk");
    let second_path = directory.join("gray-cast-iron-s2-s-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first S2-S gray-cast-iron seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second S2-S gray-cast-iron seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "S2-S gray-cast-iron decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first S2-S gray-iron pack");
    let second_bytes = fs::read(second_path).expect("read second S2-S gray-iron pack");
    assert_eq!(first_bytes, second_bytes, "S2-S gray-iron pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode S2-S gray-iron pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify S2-S gray-iron pack identity");

    assert_eq!(decoded.pack_id(), "pearlitic-gray-cast-iron-s2-s-sr-fesi");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("Attribution 4.0 International")
    );
    assert_eq!(decoded.claims().claim_count(), 15);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    for (property, source_percent) in [
        ("carbon_mass_fraction", 3.54),
        ("silicon_mass_fraction", 1.62),
        ("manganese_mass_fraction", 0.51),
        ("phosphorus_mass_fraction", 0.025),
        ("sulfur_mass_fraction", 0.028),
        ("molybdenum_mass_fraction", 0.35),
        ("copper_mass_fraction", 0.58),
        ("tin_mass_fraction", 0.060),
        ("carbon_equivalent_ce", 4.05),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one S2-S {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("S2-S {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "S2-S {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(claim.provenance.source.contains("doi:10.3390/ma11101876"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("S2-S composition observation remains linked");
        assert!(
            observation
                .specimen
                .contains("S2-S-pearlitic-gray-cast-iron")
        );
        assert!(observation.specimen.contains("0p4wtpct-SrFeSi-Ino2"));
        assert!(observation.method.contains("Table 1"));
        assert!(observation.caveats.contains("2.0 wt% Sr"));
        assert!(observation.caveats.contains("no balance-iron scalar"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let carbon_equivalent_from_printed_composition: f64 = 3.54 + 0.31 * 1.62 + 0.33 * 0.025;
    assert!(
        (carbon_equivalent_from_printed_composition - 4.05).abs() <= 0.005,
        "S2-S carbon-equivalent transcription exceeds the source's printed rounding"
    );

    for (property, expected_value, expected_dims, caveat_fragment) in [
        (
            "graphite_area_fraction",
            9.0 * 0.01,
            dimensionless,
            "graphite area 9.0 +/- 0.2 percent",
        ),
        (
            "maximum_graphite_flake_length",
            273.0 * 1.0e-6,
            Dims([1, 0, 0, 0, 0, 0]),
            "maximum graphite length 273 +/- 19 um",
        ),
        (
            "primary_dendrite_area_fraction",
            15.6 * 0.01,
            dimensionless,
            "primary-dendrite area 15.6 +/- 0.9 percent",
        ),
        (
            "eutectic_colony_areal_density",
            371.0 * 1.0e4,
            Dims([-2, 0, 0, 0, 0, 0]),
            "eutectic-colony count 371 +/- 19 per cm2",
        ),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one S2-S {property} claim");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("S2-S {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "S2-S {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("S2-S microstructure observation remains linked");
        assert!(observation.method.contains("eight cross-section fields"));
        assert!(observation.caveats.contains("type-A graphite"));
        assert!(observation.caveats.contains(caveat_fragment));
        assert!(observation.caveats.contains("one standard deviation"));
        assert!(
            observation
                .caveats
                .contains("runtime uncertainty remains Unstated")
        );
    }

    let (_, tensile) = decoded
        .claims()
        .claims_for("ultimate_tensile_strength")
        .into_iter()
        .next()
        .expect("S2-S ultimate tensile strength claim");
    let PropertyValue::Scalar {
        value: tensile_value,
        dims: tensile_dims,
    } = &tensile.value
    else {
        panic!("S2-S ultimate tensile strength was not scalar");
    };
    assert_eq!(*tensile_dims, Dims([-1, 1, -2, 0, 0, 0]));
    assert_eq!((*tensile_value).to_bits(), 326.0e6_f64.to_bits());
    assert_eq!(
        tensile.validity.bound("source_test_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(tensile.uncertainty, UncertaintyModel::Unstated);
    let tensile_observation = decoded
        .claims()
        .observation(tensile.observations[0])
        .expect("S2-S tensile observation remains linked");
    assert!(tensile_observation.method.contains("GB/T T228.1-2010"));
    assert!(tensile_observation.method.contains("three tests averaged"));
    assert!(tensile_observation.caveats.contains("nearest 1 MPa"));
    assert!(tensile_observation.caveats.contains("approximately 8 MPa"));
    assert!(
        tensile_observation
            .caveats
            .contains("exact test temperature is not reported")
    );

    let (_, conductivity) = decoded
        .claims()
        .claims_for("thermal_conductivity")
        .into_iter()
        .next()
        .expect("S2-S thermal conductivity claim");
    let PropertyValue::Scalar {
        value: conductivity_value,
        dims: conductivity_dims,
    } = &conductivity.value
    else {
        panic!("S2-S thermal conductivity was not scalar");
    };
    assert_eq!(*conductivity_dims, Dims([1, 1, -3, -1, 0, 0]));
    assert_eq!((*conductivity_value).to_bits(), 58.8f64.to_bits());
    assert_eq!(
        conductivity.validity.bound("source_test_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(conductivity.validity.bounds().len(), 1);
    assert_eq!(conductivity.uncertainty, UncertaintyModel::Unstated);
    let conductivity_observation = decoded
        .claims()
        .observation(conductivity.observations[0])
        .expect("S2-S thermal observation remains linked");
    assert!(conductivity_observation.method.contains("NETZSCH LFA 457"));
    assert!(
        conductivity_observation
            .method
            .contains("Archimedes density")
    );
    assert!(
        conductivity_observation
            .caveats
            .contains("nearest 0.1 W/(m K)")
    );
    assert!(
        conductivity_observation
            .caveats
            .contains("approximately 0.3 W/(m K)")
    );
    assert!(
        conductivity_observation
            .caveats
            .contains("exact room temperature is not reported")
    );

    // G3 independent-source plausibility evidence only: ORNL/TM-2012/506
    // Appendix C gives a broad 42..62 W/(m K) range for generic gray cast
    // iron. It neither identifies S2-S nor overwrites the primary claim.
    assert!((42.0..=62.0).contains(conductivity_value));

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_nasa9_regions_into_identical_verified_model_packs() {
    let (directory, manifest) = write_nasa9_fixture();
    let first_path = directory.join("nasa9-first.fsmodpk");
    let second_path = directory.join("nasa9-second.fsmodpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-9 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-9 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "NASA-9 decision stream moved");
    assert_decision_compiler(&first, NASA9_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NASA-9 pack");
    let second_bytes = fs::read(second_path).expect("read second NASA-9 pack");
    assert_eq!(first_bytes, second_bytes, "NASA-9 pack bytes moved");
    assert_eq!(first_bytes.len(), NASA9_PACK_BYTES_GOLDEN);
    let decoded = NormalizedModelPack::from_bytes(&first_bytes).expect("decode NASA-9 model pack");
    assert_eq!(decoded.content_hash().to_string(), NASA9_PACK_HASH_GOLDEN);
    let decoded = NormalizedModelPack::from_bytes_verified(decoded.content_hash(), &first_bytes)
        .expect("verified NASA-9 model pack");
    assert_eq!(decoded.pack_id(), "N2");
    assert_eq!(
        decoded.compiler(),
        "frankensim-matdb-nasa9-model-pack-compiler-v1"
    );
    assert_eq!(decoded.models().len(), 2);
    assert_eq!(decoded.normalizations().len(), 24);
    assert!(decoded.models().iter().all(|card| {
        card.law.0 == "nasa9-standard-state"
            && card.law_version == 1
            && card.parameters.len() == 10
            && card.validity.bound("T").is_some()
            && card.provenance.source.contains("[species:N2]")
    }));

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"species_pack_id_bound\""));
    assert!(decisions.contains("\"reason_code\":\"nasa9_region_normalized\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_model_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{NASA9_PACK_HASH_GOLDEN}\"")))
    );
}

#[test]
fn g3_cli_compiles_first_order_kinetics_into_an_identical_verified_model_pack() {
    let (directory, manifest) = write_kinetics_fixture();
    let first_path = directory.join("kinetics-first.fsmodpk");
    let second_path = directory.join("kinetics-second.fsmodpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first kinetics compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second kinetics compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "kinetics decision stream moved"
    );
    assert_decision_compiler(&first, KINETICS_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first kinetics pack");
    let second_bytes = fs::read(second_path).expect("read second kinetics pack");
    assert_eq!(first_bytes, second_bytes, "kinetics pack bytes moved");
    let decoded =
        NormalizedModelPack::from_bytes(&first_bytes).expect("decode kinetics model pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedModelPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verified kinetics model pack");
    assert_eq!(decoded.pack_id(), "water-formation");
    assert_eq!(
        decoded.compiler(),
        "frankensim-matdb-kinetics-model-pack-compiler-v1"
    );
    assert_eq!(decoded.models().len(), 1);
    assert_eq!(decoded.normalizations().len(), 4);
    let card = &decoded.models()[0];
    assert_eq!(card.law.0, "arrhenius-first-order-rate");
    assert_eq!(card.law_version, 1);
    assert_eq!(
        card.parameters["activation_temperature"].dims,
        Dims([0, 0, 0, 1, 0, 0])
    );
    assert_eq!(
        card.parameters["pre_exponential"].dims,
        Dims([0, 0, -1, 0, 0, 0])
    );
    assert!(
        card.provenance
            .source
            .contains("[reaction:water-formation]")
    );
    assert!(card.provenance.source.contains("[rate-basis:first-order]"));

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"reaction_pack_id_bound\""));
    assert!(decisions.contains("\"reason_code\":\"kinetics_reaction_normalized\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_model_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
#[allow(clippy::too_many_lines)] // The end-to-end receipt audit is clearer as one ordered assertion path.
fn g3_cli_compiles_species_association_into_identical_verified_species_packs() {
    let (directory, manifest) = write_species_fixture(SPECIES_SOURCE);
    let first_path = directory.join("species-first.fsspcpk");
    let second_path = directory.join("species-second.fsspcpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first species compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second species compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "species decision stream moved");
    assert_decision_compiler(&first, SPECIES_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first species pack");
    let second_bytes = fs::read(second_path).expect("read second species pack");
    assert_eq!(first_bytes, second_bytes, "species pack bytes moved");
    let decoded = NormalizedSpeciesPack::from_bytes(&first_bytes).expect("decode species pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedSpeciesPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verified species pack");
    assert_eq!(decoded.pack_id(), "N2");
    assert_eq!(
        decoded.compiler(),
        "frankensim-matdb-species-pack-compiler-v1"
    );
    let association = decoded.association();
    assert_eq!(association.species().as_str(), "N2");
    assert_eq!(association.molar_mass().to_bits(), 0.028_013_4f64.to_bits());
    assert_eq!(association.standard_state_phase(), "gas");
    assert_eq!(association.reference_eos(), "ideal-gas");
    assert_eq!(
        association.reference_pressure().to_bits(),
        100_000.0f64.to_bits()
    );
    assert_eq!(association.elemental_reference(), "NASA-TP-2002-211556");
    assert_eq!(association.sources().len(), 2);
    assert_eq!(association.provenance().license.as_str(), "CC-BY-4.0");
    assert!(
        association
            .provenance()
            .artifact
            .is_some_and(|artifact| association.sources().contains(&artifact))
    );
    assert_eq!(decoded.normalizations().len(), 2);
    assert_eq!(
        decoded.normalizations()[0].target(),
        SpeciesNormalizationTarget::MolarMass
    );
    assert_eq!(decoded.normalizations()[0].dims(), SPECIES_MOLAR_MASS_DIMS);
    assert_eq!(
        decoded.normalizations()[0].scale().to_bits(),
        0.001f64.to_bits()
    );
    assert_eq!(
        decoded.normalizations()[0].offset().to_bits(),
        0.0f64.to_bits()
    );
    assert_eq!(decoded.normalizations()[0].source_basis(), "g/mol");
    assert_eq!(
        decoded.normalizations()[0].target_basis(),
        SPECIES_PACK_TARGET_BASIS
    );
    assert_eq!(
        decoded.normalizations()[1].target(),
        SpeciesNormalizationTarget::ReferencePressure
    );
    assert_eq!(
        decoded.normalizations()[1].dims(),
        SPECIES_REFERENCE_PRESSURE_DIMS
    );
    assert_eq!(
        decoded.normalizations()[1].scale().to_bits(),
        1_000.0f64.to_bits()
    );

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"species_pack_id_bound\""));
    assert!(decisions.contains("\"reason_code\":\"species_association_normalized\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_species_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_methane_seed_and_records_independent_agreement() {
    let manifest = workspace_path(METHANE_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed methane seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("methane-first.fsspcpk");
    let second_path = directory.join("methane-second.fsspcpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first methane seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second methane seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "methane decision stream moved");
    assert_decision_compiler(&first, SPECIES_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first methane pack");
    let second_bytes = fs::read(second_path).expect("read second methane pack");
    assert_eq!(first_bytes, second_bytes, "methane pack bytes moved");
    let decoded = NormalizedSpeciesPack::from_bytes(&first_bytes).expect("decode methane pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedSpeciesPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify methane pack identity");

    assert_eq!(decoded.pack_id(), "CH4");
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    let association = decoded.association();
    assert_eq!(association.species().as_str(), "CH4");
    assert_eq!(
        association.molar_mass().to_bits(),
        (NASA_METHANE_MOLAR_MASS_G_PER_MOL * 0.001).to_bits()
    );
    assert_eq!(association.standard_state_phase(), "gas");
    assert_eq!(association.reference_eos(), "ideal-gas");
    assert_eq!(
        association.reference_pressure().to_bits(),
        100_000.0f64.to_bits()
    );
    assert_eq!(
        association.elemental_reference(),
        "NASA-TP-2002-211556-reference-elements-298.15K-1bar"
    );
    assert_eq!(association.provenance().license.as_str(), NASA_SEED_LICENSE);
    assert!(
        association
            .provenance()
            .source
            .contains("NASA/TP-2002-211556")
    );
    assert!(association.provenance().source.contains("[source:primary]"));

    let independent_difference =
        (association.molar_mass() - NIST_SRD69_METHANE_MOLAR_MASS_KG_PER_MOL).abs();
    assert!(
        independent_difference <= NIST_SRD69_DISPLAY_ROUNDING_HALF_WIDTH_KG_PER_MOL,
        "NASA seed and NIST SRD 69 display disagree beyond the recorded rounding band: {independent_difference:e} kg/mol"
    );
}

#[test]
fn g3_cli_compiles_committed_air_exhaust_constituents_without_inventing_a_mixture() {
    for seed in AIR_EXHAUST_SPECIES_SEEDS {
        let manifest = workspace_path(seed.manifest);
        assert!(
            manifest.is_file(),
            "committed {} seed manifest is missing",
            seed.species
        );
        let directory = fixture_dir();
        let first_path = directory.join(format!("{}-first.fsspcpk", seed.species));
        let second_path = directory.join(format!("{}-second.fsspcpk", seed.species));

        let first = run_compiler(&manifest, &first_path);
        let second = run_compiler(&manifest, &second_path);
        assert!(
            first.status.success(),
            "first {} seed compilation failed: {}",
            seed.species,
            String::from_utf8_lossy(&first.stderr)
        );
        assert!(
            second.status.success(),
            "second {} seed compilation failed: {}",
            seed.species,
            String::from_utf8_lossy(&second.stderr)
        );
        assert_eq!(
            first.stdout, second.stdout,
            "{} decision stream moved",
            seed.species
        );
        assert_decision_compiler(&first, SPECIES_COMPILER_ID);

        let first_bytes = fs::read(first_path).expect("read first constituent pack");
        let second_bytes = fs::read(second_path).expect("read second constituent pack");
        assert_eq!(
            first_bytes, second_bytes,
            "{} pack bytes moved",
            seed.species
        );
        let decoded =
            NormalizedSpeciesPack::from_bytes(&first_bytes).expect("decode constituent pack");
        let pack_hash = decoded.content_hash();
        let decoded = NormalizedSpeciesPack::from_bytes_verified(pack_hash, &first_bytes)
            .expect("verify constituent pack identity");

        assert_eq!(decoded.pack_id(), seed.species);
        assert!(
            decoded
                .redistribution_terms()
                .contains("public use permitted")
        );
        let association = decoded.association();
        assert_eq!(association.species().as_str(), seed.species);
        assert_eq!(
            association.molar_mass().to_bits(),
            (seed.nasa_molar_mass_g_per_mol * 0.001).to_bits()
        );
        assert_eq!(association.standard_state_phase(), "gas");
        assert_eq!(association.reference_eos(), "ideal-gas");
        assert_eq!(
            association.reference_pressure().to_bits(),
            100_000.0f64.to_bits()
        );
        assert_eq!(
            association.elemental_reference(),
            "NASA-TP-2002-211556-reference-elements-298.15K-1bar"
        );
        assert_eq!(association.provenance().license.as_str(), NASA_SEED_LICENSE);
        assert!(
            association
                .provenance()
                .source
                .contains("NASA/TP-2002-211556")
        );

        let independent_difference_g_per_mol =
            (association.molar_mass() * 1_000.0 - seed.nist_molar_mass_g_per_mol).abs();
        assert!(
            independent_difference_g_per_mol <= seed.nist_display_rounding_half_width_g_per_mol,
            "NASA {} seed and NIST SRD 69 display disagree beyond the recorded rounding band: {independent_difference_g_per_mol:e} g/mol",
            seed.species
        );
    }
}

#[test]
fn g3_cli_compiles_committed_nist_srm_1720_northern_continental_air() {
    let manifest = workspace_path(NIST_SRM_1720_NORTHERN_CONTINENTAL_AIR_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NIST SRM 1720 northern-continental-air manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nist-srm-1720-first.fsmatpk");
    let second_path = directory.join("nist-srm-1720-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NIST SRM 1720 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NIST SRM 1720 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NIST SRM 1720 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NIST SRM 1720 pack");
    let second_bytes = fs::read(second_path).expect("read second NIST SRM 1720 pack");
    assert_eq!(first_bytes, second_bytes, "NIST SRM 1720 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode NIST SRM 1720 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NIST SRM 1720 pack identity");

    assert_eq!(decoded.pack_id(), "nist-srm-1720-northern-continental-air");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("NIST"));
    assert_eq!(decoded.claims().claim_count(), 4);
    assert!(decoded.joint_statistics().is_empty());

    let expected_claims: [(&str, f64); 4] = [
        ("information_oxygen_amount_fraction", 20.93),
        ("information_argon_amount_fraction", 0.935),
        (
            "information_carbon_monoxide_amount_fraction_lower_bound",
            0.000013,
        ),
        (
            "information_carbon_monoxide_amount_fraction_upper_bound",
            0.000018,
        ),
    ];
    let mut carbon_monoxide_bounds = Vec::new();
    for (property, source_percent) in expected_claims {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique SRM 1720 {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("SRM 1720 {property} was not scalar");
        };
        let expected_value: f64 = source_percent * 0.01;
        let comparison_scale: f64 = expected_value.abs().max(1.0e-12);
        assert!((*value - expected_value).abs() / comparison_scale <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NIST_PUBLIC_INFORMATION_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("Standard Reference Material 1720")
        );
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nist-srm-1720-archived-certificate]")
        );

        for (axis, expected) in [
            ("source_value_is_information_not_certified", 1.0),
            ("source_composition_basis_is_amount_fraction", 1.0),
            ("source_nitrogen_is_balance_gas", 1.0),
            ("source_air_was_scrubbed_of_moisture", 1.0),
            ("source_remaining_humidity_quantified", 0.0),
            ("source_cylinder_identity_known", 0.0),
            ("source_certified_greenhouse_values_present", 0.0),
            ("source_use_temperature_known", 0.0),
            ("source_use_pressure_known", 0.0),
            ("source_is_northern_continental_air_lot", 1.0),
            ("source_is_universal_air_composition", 0.0),
            ("source_certificate_is_archived_and_expired", 1.0),
        ] {
            assert_eq!(
                claim.validity.bound(axis),
                Some((expected, expected)),
                "SRM 1720 {property} moved validity axis {axis}"
            );
        }

        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("SRM 1720 information observation remains linked");
        assert!(observation.method.contains("information-value table"));
        assert!(
            observation
                .caveats
                .contains("cannot establish metrological traceability")
        );
        assert!(observation.caveats.contains("SAMPLE placeholders"));
        assert!(
            observation
                .caveats
                .contains("not a universal dry-air composition")
        );

        if property.contains("carbon_monoxide") {
            carbon_monoxide_bounds.push(*value);
        }
    }
    carbon_monoxide_bounds.sort_by(f64::total_cmp);
    assert_eq!(carbon_monoxide_bounds.len(), 2);
    assert!(carbon_monoxide_bounds[0] < carbon_monoxide_bounds[1]);
    for absent in [
        "certified_carbon_dioxide_amount_fraction",
        "certified_methane_amount_fraction",
        "certified_nitrous_oxide_amount_fraction",
        "information_nitrogen_amount_fraction",
    ] {
        assert!(
            decoded.claims().claims_for(absent).is_empty(),
            "SRM 1720 must not invent {absent}"
        );
    }
}

#[test]
fn g3_cli_compiles_committed_nist_srm_2728_auto_emission_reference_gas() {
    let manifest = workspace_path(NIST_SRM_2728_AUTO_EMISSION_REFERENCE_GAS_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NIST SRM 2728 reference-gas manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nist-srm-2728-first.fsmatpk");
    let second_path = directory.join("nist-srm-2728-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NIST SRM 2728 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NIST SRM 2728 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NIST SRM 2728 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NIST SRM 2728 pack");
    let second_bytes = fs::read(second_path).expect("read second NIST SRM 2728 pack");
    assert_eq!(first_bytes, second_bytes, "NIST SRM 2728 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode NIST SRM 2728 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NIST SRM 2728 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nist-srm-2728-auto-emission-reference-gas"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("NIST"));
    assert_eq!(decoded.claims().claim_count(), 4);
    assert!(decoded.joint_statistics().is_empty());

    let expected_claims: [(&str, f64, bool, bool); 4] = [
        ("nominal_carbon_dioxide_amount_fraction", 14.0, true, false),
        ("nominal_carbon_monoxide_amount_fraction", 8.0, true, false),
        ("nominal_propane_amount_fraction", 0.3, true, false),
        (
            "information_total_other_hydrocarbons_propane_equivalent_amount_fraction",
            0.0008,
            false,
            true,
        ),
    ];
    for (property, source_percent, is_nominal, is_information) in expected_claims {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique SRM 2728 {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("SRM 2728 {property} was not scalar");
        };
        let expected_value: f64 = source_percent * 0.01;
        let comparison_scale: f64 = expected_value.abs().max(1.0e-12);
        assert!((*value - expected_value).abs() / comparison_scale <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NIST_PUBLIC_INFORMATION_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("Standard Reference Material 2728")
        );
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nist-srm-2728-archived-certificate]")
        );

        for (axis, expected) in [
            ("source_composition_basis_is_amount_fraction", 1.0),
            ("source_nitrogen_is_balance_gas", 1.0),
            ("source_cylinder_identity_known", 0.0),
            ("source_certified_value_and_95pct_interval_present", 0.0),
            ("source_is_auto_emission_calibration_gas", 1.0),
            ("source_is_engine_generated_exhaust_sample", 0.0),
            ("source_mixture_temperature_known", 0.0),
            ("source_mixture_pressure_known", 0.0),
            ("source_oxygen_water_nox_fractions_known", 0.0),
            ("source_certificate_is_archived_template", 1.0),
        ] {
            assert_eq!(
                claim.validity.bound(axis),
                Some((expected, expected)),
                "SRM 2728 {property} moved validity axis {axis}"
            );
        }
        let nominal_axis = if is_nominal { 1.0 } else { 0.0 };
        let information_axis = if is_information { 1.0 } else { 0.0 };
        assert_eq!(
            claim
                .validity
                .bound("source_value_is_nominal_not_cylinder_certified"),
            Some((nominal_axis, nominal_axis))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_value_is_information_not_certified"),
            Some((information_axis, information_axis))
        );

        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("SRM 2728 composition observation remains linked");
        assert!(observation.method.contains("nominal composition"));
        assert!(
            observation
                .caveats
                .contains("95 percent confidence intervals blank")
        );
        assert!(
            observation
                .caveats
                .contains("not a sampled or equilibrium engine exhaust")
        );
        assert!(
            observation
                .caveats
                .contains("Nitrogen is identified only as the balance gas")
        );
    }

    assert!(
        decoded
            .claims()
            .claims_for("nominal_nitrogen_amount_fraction")
            .is_empty(),
        "SRM 2728 nitrogen balance must not become an inferred scalar"
    );
}

#[test]
fn g3_cli_refuses_malformed_species_without_publishing() {
    let malformed = SPECIES_SOURCE.replacen("28.0134\tg/mol", "28.0134\tkg", 1);
    assert_ne!(malformed, SPECIES_SOURCE);
    let (directory, manifest) = write_species_fixture(&malformed);
    let output = directory.join("refused-species.fsspcpk");

    let refused = run_compiler(&manifest, &output);
    assert!(
        !refused.status.success(),
        "invalid species unexpectedly compiled"
    );
    assert!(!output.exists(), "species refusal published an output");
    assert_decision_compiler(&refused, SPECIES_COMPILER_ID);
    let decisions = String::from_utf8(refused.stdout).expect("decision stream is UTF-8");
    assert_eq!(decisions.matches("\"verdict\":\"refuse\"").count(), 1);
    assert!(decisions.contains("\"reason_code\":\"species_molar_mass_dims_mismatch\""));
    assert!(decisions.contains("\"subject\":\"species:N2\""));
    assert!(
        String::from_utf8_lossy(&refused.stderr)
            .contains("error: matdb pack refused [species_molar_mass_dims_mismatch]")
    );
}

#[test]
fn g3_cli_compiles_committed_nasa_cr_195445_omc_ps200_rotary_coating_system() {
    let manifest = workspace_path(NASA_CR_195445_OMC_PS200_ROTARY_COATING_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NASA-CR-195445 OMC PS-200 coating manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nasa-cr-195445-omc-ps200-first.fsmatpk");
    let second_path = directory.join("nasa-cr-195445-omc-ps200-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-CR-195445 OMC PS-200 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-CR-195445 OMC PS-200 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NASA-CR-195445 OMC PS-200 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first OMC PS-200 pack");
    let second_bytes = fs::read(second_path).expect("read second OMC PS-200 pack");
    assert_eq!(first_bytes, second_bytes, "OMC PS-200 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode OMC PS-200 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify OMC PS-200 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nasa-cr-195445-omc-ps200-rotary-coating-system"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert!(decoded.redistribution_terms().contains("US4728448A"));
    assert_eq!(decoded.claims().claim_count(), 7);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let mut composition_sum = 0.0_f64;
    let mut composition_observation = None;
    for (property, source_percent) in [
        (
            "ps200_bonded_chromium_carbide_feedstock_mass_fraction",
            80.0_f64,
        ),
        ("ps200_silver_feedstock_mass_fraction", 10.0_f64),
        ("ps200_baf2_caf2_eutectic_feedstock_mass_fraction", 10.0_f64),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one PS-200 {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PS-200 {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        assert_eq!(value.to_bits(), f64::to_bits(expected_value));
        composition_sum += *value;
        for required_axis in [
            "source_feedstock_composition_basis_is_mass_fraction",
            "source_ps200_plasma_sprayed_in_engine_report",
        ] {
            assert_eq!(claim.validity.bound(required_axis), Some((1.0, 1.0)));
        }
        for missing_or_refused_axis in [
            "source_post_spray_phase_fractions_known",
            "source_powder_lot_known",
            "source_patent_practice_license_granted",
        ] {
            assert_eq!(
                claim.validity.bound(missing_or_refused_axis),
                Some((0.0, 0.0))
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.provenance.license,
            PUBLIC_USE_AND_PATENT_PUBLICATION_LICENSE
        );
        assert!(claim.provenance.source.contains("NASA-CR-195445"));
        assert!(claim.provenance.source.contains("US Patent 4,728,448"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PS-200 composition observation remains linked");
        assert!(observation.specimen.contains("80wtpct-nickel-bonded-Cr3C2"));
        assert!(observation.method.contains("US4728448A Table II"));
        assert!(observation.caveats.contains("pre-spray feedstock"));
        assert!(observation.caveats.contains("not a patent-practice"));
        match composition_observation {
            Some(expected) => assert_eq!(claim.observations[0], expected),
            None => composition_observation = Some(claim.observations[0]),
        }
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }
    assert!((composition_sum - 1.0).abs() <= f64::EPSILON);

    let roughness_claims = decoded.claims().claims_for("surface_roughness_rms");
    assert_eq!(roughness_claims.len(), 4);
    let mut observed_conditions = Vec::new();
    for (id, claim) in roughness_claims {
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("OMC PS-200 RMS finish was not scalar");
        };
        assert_eq!(*dims, Dims([1, 0, 0, 0, 0, 0]));
        let (test_number, _) = claim
            .validity
            .bound("source_engine_test_number")
            .expect("engine test number remains pinned");
        let (stage_index, _) = claim
            .validity
            .bound("source_surface_stage_index")
            .expect("surface stage remains pinned");
        let condition = (test_number as u8, stage_index as u8);
        let expected_from_microinch: f64 = match condition {
            (3, 0) => 21.0 * 25.4e-9,
            (3, 1) => 7.0 * 25.4e-9,
            (6, 0) => 24.0 * 25.4e-9,
            (6, 1) => 17.0 * 25.4e-9,
            other => panic!("unexpected OMC PS-200 test/stage condition: {other:?}"),
        };
        assert!(
            (*value - expected_from_microinch).abs() <= 1.0e-21,
            "OMC PS-200 microinch-to-metre transcription moved for {condition:?}"
        );
        observed_conditions.push(condition);
        for required_axis in [
            "source_ps200_over_zirconia_or_sx331",
            "source_substrate_is_aluminum_alloy",
        ] {
            assert_eq!(claim.validity.bound(required_axis), Some((1.0, 1.0)));
        }
        for missing_axis in [
            "source_aluminum_alloy_grade_known",
            "source_coating_thickness_known",
            "source_surface_finish_method_known",
        ] {
            assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.provenance.license,
            PUBLIC_USE_AND_PATENT_PUBLICATION_LICENSE
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("OMC PS-200 finish observation remains linked");
        match test_number as u8 {
            3 => {
                assert_eq!(
                    claim.validity.bound("source_narrative_run_duration"),
                    Some((2.5 * 3_600.0, 2.5 * 3_600.0))
                );
                assert!(observation.specimen.contains("Test3-air-cooled-OMC"));
                assert!(observation.caveats.contains("TBC crack"));
                assert!(observation.caveats.contains("scrap that housing"));
                assert!(observation.caveats.contains("no zero-wear"));
            }
            6 => {
                assert_eq!(
                    claim.validity.bound("source_actual_run_duration"),
                    Some((1.5 * 3_600.0, 1.5 * 3_600.0))
                );
                assert_eq!(
                    claim
                        .validity
                        .bound("source_local_ps200_breakthrough_present"),
                    Some((1.0, 1.0))
                );
                assert!(observation.specimen.contains("Test6-air-cooled-OMC"));
                assert!(observation.caveats.contains("500 degF"));
                assert!(observation.caveats.contains("local PS-200 breakthrough"));
                assert!(observation.caveats.contains("do not establish a wear rate"));
            }
            other => panic!("unexpected OMC PS-200 engine test number: {other}"),
        }
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }
    observed_conditions.sort_unstable();
    assert_eq!(observed_conditions, [(3, 0), (3, 1), (6, 0), (6, 1)]);

    for refused_property in [
        "coefficient_of_friction",
        "wear_rate",
        "coating_thickness",
        "thermal_conductivity",
        "specific_fuel_consumption",
        "service_life",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent or system-level coating property must remain refused: {refused_property}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_uses_generic_driver_identity_for_an_unknown_profile() {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    let unsupported = MANIFEST.replace("material-tsv-v1", "future-profile-v1");
    fs::write(&manifest, unsupported).expect("write unsupported-profile manifest");
    fs::write(directory.join("source.tsv"), SOURCE).expect("write source fixture");
    let output = directory.join("unsupported.fsmatpk");

    let refused = run_compiler(&manifest, &output);
    assert!(!refused.status.success());
    assert!(!output.exists(), "unsupported profile published an output");
    assert_decision_compiler(&refused, MATERIAL_COMPILER_ID);
    assert!(
        String::from_utf8_lossy(&refused.stdout)
            .contains("\"reason_code\":\"unsupported_source_profile\"")
    );
}

#[test]
fn g3_cli_uses_generic_driver_identity_before_profile_selection() {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    let incomplete = MANIFEST.replace("license\tCC-BY-4.0\n", "");
    fs::write(&manifest, incomplete).expect("write incomplete manifest");
    fs::write(directory.join("source.tsv"), SOURCE).expect("write source fixture");
    let output = directory.join("incomplete.fsmatpk");

    let refused = run_compiler(&manifest, &output);
    assert!(!refused.status.success());
    assert!(!output.exists(), "incomplete manifest published an output");
    assert_decision_compiler(&refused, MATERIAL_COMPILER_ID);
    assert!(
        String::from_utf8_lossy(&refused.stdout).contains("\"reason_code\":\"missing_license\"")
    );
}

#[test]
fn g3_cli_retains_prior_admissions_when_later_claim_refuses() {
    let rejected_source = SOURCE.replacen(
        "uncertainty\tmodulus\trelative\t2\t%\t0.95\t1\n",
        "uncertainty\tmodulus\tabsolute\t2\tkg\t0.95\t1\n",
        1,
    );
    assert_ne!(rejected_source, SOURCE);
    let (directory, manifest) = write_fixture(&rejected_source);
    let first_path = directory.join("refused-first.fsmatpk");
    let second_path = directory.join("refused-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(!first.status.success());
    assert!(!second.status.success());
    assert_eq!(first.stdout, second.stdout, "refusal transcript moved");
    assert!(!first_path.exists(), "refused compilation published output");
    assert!(
        !second_path.exists(),
        "refused compilation published output"
    );

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"subject\":\"source:primary\""));
    assert!(decisions.contains("\"subject\":\"claim:density\""));
    assert!(decisions.contains("\"reason_code\":\"claim_normalized\""));
    assert_eq!(decisions.matches("\"verdict\":\"refuse\"").count(), 1);
    let refusal = decisions
        .lines()
        .find(|row| row.contains("\"verdict\":\"refuse\""))
        .expect("one terminal refusal row");
    assert!(refusal.contains("\"reason_code\":\"uncertainty_dims_mismatch\""));
    assert!(!refusal.contains("\"source_hash\":\"\""));
    assert!(refusal.contains("\"pack_hash\":\"\""));
    assert!(
        decisions
            .lines()
            .filter(|row| row.contains("\"verdict\":\"admit\""))
            .all(|row| row.contains("\"pack_hash\":\"\""))
    );
    assert!(
        String::from_utf8_lossy(&first.stderr)
            .contains("error: matdb pack refused [uncertainty_dims_mismatch]")
    );
}

const TONEWOOD_SEED_PACKS: [&str; 23] = [
    "ash-white-fpl-gtr282",
    "baldcypress-fpl-gtr282",
    "balsa-fpl-gtr282",
    "basswood-american-fpl-gtr282",
    "birch-yellow-fpl-gtr282",
    "cedar-western-red-fpl-gtr282",
    "cherry-black-fpl-gtr282",
    "cottonwood-eastern-fpl-gtr282",
    "douglas-fir-coast-fpl-gtr282",
    "hemlock-western-fpl-gtr282",
    "larch-western-fpl-gtr282",
    "mahogany-african-fpl-gtr282",
    "mahogany-honduras-fpl-gtr282",
    "maple-red-fpl-gtr282",
    "maple-sugar-fpl-gtr282",
    "redwood-old-growth-fpl-gtr282",
    "rosewood-brazilian-fpl-gtr282",
    "rosewood-indian-fpl-gtr282",
    "spruce-engelmann-fpl-gtr282",
    "spruce-sitka-fpl-gtr282",
    "sweetgum-fpl-gtr282",
    "walnut-black-fpl-gtr282",
    "yellow-poplar-fpl-gtr282",
];

const TONEWOOD_RATIO_PROPERTIES: [&str; 5] = [
    "et_over_el",
    "er_over_el",
    "glr_over_el",
    "glt_over_el",
    "grt_over_el",
];

const TONEWOOD_POISSON_PROPERTIES: [&str; 6] =
    ["nu_lr", "nu_lt", "nu_rt", "nu_tr", "nu_rl", "nu_tl"];

/// The unique scalar SI value for `property`, or `None` when the pack has
/// no claim for it (the source printed no value).
fn tonewood_scalar(pack: &NormalizedPack, property: &str) -> Option<f64> {
    let claims = pack.claims().claims_for(property);
    if claims.is_empty() {
        return None;
    }
    assert_eq!(
        claims.len(),
        1,
        "tonewood packs carry one claim per property, got {} for {property}",
        claims.len()
    );
    match &claims[0].1.value {
        PropertyValue::Scalar { value, .. } => Some(*value),
        PropertyValue::Curve { .. } => panic!("{property} must be a scalar claim"),
    }
}

/// Longitudinal sound speed gate: c_L = sqrt(E_L/rho) must land inside the
/// published clear-wood range of roughly 3,000-6,500 m/s longitudinally.
/// NOTE: FPL-GTR-282 chapter 5 (Vibration Properties) prints a worked
/// example "12.4 GPa / 480 kg/m^3 -> about 3,800 m/s (12,500 ft/s)" that is
/// internally inconsistent — sqrt(12.4e9/480) = 5,082 m/s — so the window is
/// anchored to the species-table spread, not that example. A GPa/MPa or
/// g/cm^3 slip moves c_L by >5x, far outside the window either way.
fn tonewood_sound_speed_in_range(el_pa: f64, rho_si: f64) -> bool {
    let c = (el_pa / rho_si).sqrt();
    (3_000.0..6_500.0).contains(&c)
}

#[test]
#[allow(clippy::too_many_lines)] // one coherent corpus gate
fn g2_cli_compiles_fpl_gtr282_tonewood_seeds_with_derived_quantity_gates() {
    let mut complete_orthotropic = 0usize;
    let mut speed_gated = 0usize;
    for slug in TONEWOOD_SEED_PACKS {
        let manifest = workspace_path(&format!("data/matdb/seed-v1/{slug}/manifest.tsv"));
        let out = fixture_dir().join(format!("{slug}.fsmatpk"));
        let run = run_compiler(&manifest, &out);
        assert!(
            run.status.success(),
            "{slug} refused: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let bytes = fs::read(&out).expect("read compiled tonewood pack");
        let pack = NormalizedPack::from_bytes_verified(
            NormalizedPack::from_bytes(&bytes)
                .expect("decode tonewood pack")
                .content_hash(),
            &bytes,
        )
        .expect("verified tonewood pack");
        assert_eq!(pack.pack_id(), slug);

        let sg = tonewood_scalar(&pack, "specific_gravity").expect("specific gravity present");
        assert!(
            (0.1..1.1).contains(&sg),
            "{slug}: specific gravity {sg} implausible"
        );
        let moe =
            tonewood_scalar(&pack, "modulus_of_elasticity_bending").expect("bending MOE present");
        let el = tonewood_scalar(&pack, "young_modulus_longitudinal").expect("E_L present");
        assert!(
            (el / moe - 1.10).abs() < 1e-9,
            "{slug}: E_L must be exactly the 10% shear-corrected bending MOE"
        );
        assert!(
            (3.0e9..25.0e9).contains(&el),
            "{slug}: E_L = {el:.3e} Pa outside the clear-wood window (unit slip?)"
        );
        if let Some(rho) = tonewood_scalar(&pack, "density") {
            assert!(
                (300.0..1200.0).contains(&rho),
                "{slug}: density {rho} kg/m^3 implausible"
            );
            assert!(
                ((rho / (1000.0 * sg * 1.12)) - 1.0).abs() < 1e-9,
                "{slug}: density must equal the declared 1000*SG*1.12 derivation"
            );
            assert!(
                tonewood_sound_speed_in_range(el, rho),
                "{slug}: c_L = {:.0} m/s outside 3000-6500 (unit slip?)",
                (el / rho).sqrt()
            );
            speed_gated += 1;
        }
        let ratio_count = TONEWOOD_RATIO_PROPERTIES
            .iter()
            .filter_map(|p| tonewood_scalar(&pack, p))
            .inspect(|v| assert!((0.0..0.25).contains(v), "{slug}: ratio {v} implausible"))
            .count();
        let nu_count = TONEWOOD_POISSON_PROPERTIES
            .iter()
            .filter_map(|p| tonewood_scalar(&pack, p))
            .inspect(|v| assert!((0.0..1.0).contains(v), "{slug}: Poisson {v} implausible"))
            .count();
        if ratio_count >= 4 && nu_count >= 4 {
            complete_orthotropic += 1;
        }
    }
    assert!(
        complete_orthotropic >= 21,
        "expected >=21 complete orthotropic tonewood sets, found {complete_orthotropic}"
    );
    assert!(
        speed_gated >= 18,
        "expected >=18 density-bearing packs under the sound-speed gate, found {speed_gated}"
    );

    // The machine-readable axis-convention contract ships beside the data
    // and names every property family the packs use.
    let convention = fs::read_to_string(workspace_path(
        "data/matdb/seed-v1/instrument-axis-convention.tsv",
    ))
    .expect("instrument axis-convention file present");
    assert!(
        convention.starts_with("frankensim.instrument-axis-convention.v1"),
        "axis-convention header moved"
    );
    for property in TONEWOOD_RATIO_PROPERTIES
        .iter()
        .chain(TONEWOOD_POISSON_PROPERTIES.iter())
        .chain(
            [
                "modulus_of_elasticity_bending",
                "young_modulus_longitudinal",
                "density",
                "specific_gravity",
            ]
            .iter(),
        )
    {
        assert!(
            convention.contains(property),
            "axis-convention file must document {property}"
        );
    }

    println!(
        "{{\"suite\":\"xtask-matdb\",\"case\":\"tonewood-seeds\",\"packs\":{},\"complete_orthotropic\":{complete_orthotropic},\"speed_gated\":{speed_gated},\"verdict\":\"pass\"}}",
        TONEWOOD_SEED_PACKS.len()
    );
}

#[test]
fn g3_tonewood_unit_swap_mutation_is_caught_by_the_sound_speed_gate() {
    // MUTATION: a MPa->GPa slip on the modulus rows compiles fine (same
    // dimensions), so the dims check alone CANNOT catch it. The derived
    // sound-speed gate must.
    let source = fs::read_to_string(workspace_path(
        "data/matdb/seed-v1/spruce-sitka-fpl-gtr282/properties.tsv",
    ))
    .expect("read sitka properties");
    let doctored = source.replace("\tMPa\t", "\tGPa\t");
    assert_ne!(source, doctored, "the mutation must actually change units");
    let manifest_text = fs::read_to_string(workspace_path(
        "data/matdb/seed-v1/spruce-sitka-fpl-gtr282/manifest.tsv",
    ))
    .expect("read sitka manifest");
    let directory = fixture_dir();
    fs::write(directory.join("properties.tsv"), doctored).expect("write doctored source");
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, manifest_text).expect("write doctored manifest");
    let out = directory.join("doctored.fsmatpk");
    let run = run_compiler(&manifest, &out);
    assert!(
        run.status.success(),
        "the doctored pack must COMPILE (same dims): {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = fs::read(&out).expect("read doctored pack");
    let pack = NormalizedPack::from_bytes(&bytes).expect("decode doctored pack");
    let el = tonewood_scalar(&pack, "young_modulus_longitudinal").expect("E_L present");
    let rho = tonewood_scalar(&pack, "density").expect("density present");
    assert!(
        !tonewood_sound_speed_in_range(el, rho),
        "the sound-speed gate FAILED to catch a 1000x modulus unit slip"
    );
    println!(
        "{{\"suite\":\"xtask-matdb\",\"case\":\"tonewood-unit-swap-mutation\",\"c_l_m_per_s\":{:.0},\"verdict\":\"caught\"}}",
        (el / rho).sqrt()
    );
}

/// The instrument loss-factor and string/metal tranche
/// (bead frankensim-fsim-instrument-matdb-zwzey, second tranche): every
/// pack compiles through the fail-closed compiler and its values sit
/// inside authored plausibility windows, with derived-quantity gates
/// (sound speed, Poisson consistency, Zener frequency ordering) that
/// catch unit slips dims checks cannot.
const INSTRUMENT_TRANCHE2_PACKS: [&str; 12] = [
    "spruce-sitka-loss-qiu2026",
    "spruce-norway-danihelova2022",
    "maple-sycamore-danihelova2022",
    "aluminum-2024-t4-damping-nasa-tn-d2893",
    "aluminum-2024-t3-nasa-tn-d6448",
    "brass-yellow-half-hard-damping-nasa-tn-d1467",
    "brass-cartridge-c26000-mil-hdbk-698a",
    "phosphor-bronze-c51000-nist-mono177",
    "phosphor-bronze-5a-mil-hdbk-698a",
    "music-wire-nbs-c447",
    "music-wire-spring-moduli-fuchs-1968",
    "bronze-gunmetal-damping-rsic508",
];

#[test]
#[allow(clippy::too_many_lines)] // one coherent corpus gate
fn g2_cli_compiles_instrument_loss_and_string_tranche_with_derived_gates() {
    let mut packs = std::collections::BTreeMap::new();
    for slug in INSTRUMENT_TRANCHE2_PACKS {
        let manifest = workspace_path(&format!("data/matdb/seed-v1/{slug}/manifest.tsv"));
        let out = fixture_dir().join(format!("{slug}.fsmatpk"));
        let run = run_compiler(&manifest, &out);
        assert!(
            run.status.success(),
            "{slug} refused: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let bytes = fs::read(&out).expect("read compiled tranche-2 pack");
        let pack = NormalizedPack::from_bytes_verified(
            NormalizedPack::from_bytes(&bytes)
                .expect("decode tranche-2 pack")
                .content_hash(),
            &bytes,
        )
        .expect("verified tranche-2 pack");
        assert_eq!(pack.pack_id(), slug);
        packs.insert(slug, pack);
    }
    let scalar = |slug: &str, property: &str| -> f64 {
        tonewood_scalar(&packs[slug], property)
            .unwrap_or_else(|| panic!("{slug} must carry {property}"))
    };

    // Wood loss factors: audio-band tan-delta plausibility window.
    for (slug, property, lo, hi) in [
        (
            "spruce-sitka-loss-qiu2026",
            "loss_factor_longitudinal",
            0.005,
            0.03,
        ),
        (
            "spruce-sitka-loss-qiu2026",
            "loss_factor_longitudinal_mode5",
            0.005,
            0.03,
        ),
        (
            "spruce-norway-danihelova2022",
            "loss_factor_longitudinal",
            0.005,
            0.03,
        ),
        (
            "maple-sycamore-danihelova2022",
            "loss_factor_longitudinal",
            0.005,
            0.03,
        ),
    ] {
        let eta = scalar(slug, property);
        assert!(
            (lo..hi).contains(&eta),
            "{slug}: {property} = {eta} implausible"
        );
    }
    // Log decrement / loss factor internal consistency: eta = theta/pi.
    for slug in [
        "spruce-norway-danihelova2022",
        "maple-sycamore-danihelova2022",
    ] {
        let theta = scalar(slug, "log_decrement_longitudinal");
        let eta = scalar(slug, "loss_factor_longitudinal");
        assert!(
            (eta / (theta / core::f64::consts::PI) - 1.0).abs() < 0.02,
            "{slug}: eta {eta} must equal theta/pi = {}",
            theta / core::f64::consts::PI
        );
    }
    // Dynamic-modulus wood packs: the same c_L = sqrt(E/rho) gate the
    // FPL tranche uses (a MPa/GPa slip lands far outside 3000-6500 m/s).
    for slug in [
        "spruce-norway-danihelova2022",
        "maple-sycamore-danihelova2022",
    ] {
        let el = scalar(slug, "dynamic_young_modulus_longitudinal");
        let rho = scalar(slug, "density");
        assert!(
            tonewood_sound_speed_in_range(el, rho),
            "{slug}: c_L = {:.0} m/s outside the clear-wood window",
            (el / rho).sqrt()
        );
    }

    // Aluminum 2024-T4 vacuum damping: Zener thermal relaxation is past
    // its peak here, so the loss factor must DECREASE monotonically with
    // frequency — an ordering gate a single-window check cannot provide.
    let al_band = [
        ("loss_factor_15hz", 15.18),
        ("loss_factor_30hz", 30.30),
        ("loss_factor_70hz", 71.60),
        ("loss_factor_150hz", 153.2),
        ("loss_factor_300hz", 306.5),
        ("loss_factor_700hz", 716.2),
        ("loss_factor_1500hz", 1572.0),
    ];
    let mut previous = f64::INFINITY;
    for (property, _f) in al_band {
        let g = scalar("aluminum-2024-t4-damping-nasa-tn-d2893", property);
        assert!(
            (3.0e-5..5.0e-3).contains(&g),
            "aluminum vacuum loss factor {property} = {g} implausible"
        );
        assert!(
            g < previous,
            "aluminum loss factors must decrease with frequency: {property} = {g}"
        );
        previous = g;
    }
    // Cross-source Zener consistency: the peak height Delta/2 =
    // alpha^2 E T / (2 rho c) is GEOMETRY-INDEPENDENT, so the TN D-6448
    // wire-experiment peak (2.6e-3 at 83 Hz) must agree with the
    // TN D-2893 beam band's low-frequency plateau (2.92e-3 at 15 Hz,
    // near that geometry's own peak) — measured 11% apart; 30% is the
    // authored envelope. This replaces an earlier peak-vs-tail check
    // that a review showed was implied by the window alone.
    let peak = scalar(
        "aluminum-2024-t3-nasa-tn-d6448",
        "loss_factor_thermoelastic_peak",
    );
    assert!(
        (1.0e-3..5.0e-3).contains(&peak),
        "thermoelastic peak {peak} implausible"
    );
    let plateau = scalar("aluminum-2024-t4-damping-nasa-tn-d2893", "loss_factor_15hz");
    assert!(
        (peak / plateau - 1.0).abs() < 0.30,
        "Zener peak height must be geometry-independent across the two \
         NASA sources: wire {peak:.3e} vs beam plateau {plateau:.3e}"
    );

    // Gun-metal bronze specific damping capacity (graph-read survey):
    // window plus the recorded eta = psi/(2 pi) conversion staying in a
    // plausible cast-bronze band.
    // The compiler normalizes the source's percent unit to an SI
    // fraction, so 1 percent arrives as 0.01.
    let psi = scalar(
        "bronze-gunmetal-damping-rsic508",
        "specific_damping_capacity_5ksi_shear",
    );
    assert!(
        (0.001..0.05).contains(&psi),
        "gun metal SDC fraction = {psi} implausible"
    );
    assert!(
        (5.0e-4..1.0e-2).contains(&(psi / (2.0 * core::f64::consts::PI))),
        "gun metal implied loss factor outside the cast-bronze band"
    );

    // Brass specimen damping (graph-read, air, upper bound): window only.
    for property in [
        "specimen_loss_factor_10ksi_50hz",
        "specimen_loss_factor_10ksi_500hz",
    ] {
        let g = scalar("brass-yellow-half-hard-damping-nasa-tn-d1467", property);
        assert!(
            (1.0e-3..3.0e-2).contains(&g),
            "brass specimen loss factor {property} = {g} implausible"
        );
    }

    // Metal elastic windows and cross-source consistency.
    let checks: [(&str, &str, f64, f64); 8] = [
        (
            "brass-cartridge-c26000-mil-hdbk-698a",
            "young_modulus",
            90.0e9,
            130.0e9,
        ),
        (
            "brass-cartridge-c26000-mil-hdbk-698a",
            "density",
            8_300.0,
            8_700.0,
        ),
        (
            "phosphor-bronze-c51000-nist-mono177",
            "young_modulus",
            90.0e9,
            130.0e9,
        ),
        (
            "phosphor-bronze-5a-mil-hdbk-698a",
            "young_modulus",
            90.0e9,
            130.0e9,
        ),
        (
            "phosphor-bronze-5a-mil-hdbk-698a",
            "density",
            8_700.0,
            9_000.0,
        ),
        (
            "aluminum-2024-t3-nasa-tn-d6448",
            "young_modulus",
            60.0e9,
            80.0e9,
        ),
        (
            "music-wire-spring-moduli-fuchs-1968",
            "young_modulus",
            190.0e9,
            220.0e9,
        ),
        ("music-wire-nbs-c447", "density", 7_700.0, 7_900.0),
    ];
    for (slug, property, lo, hi) in checks {
        let value = scalar(slug, property);
        assert!(
            (lo..hi).contains(&value),
            "{slug}: {property} = {value:.3e} outside [{lo:.1e}, {hi:.1e}] (unit slip?)"
        );
    }
    // Poisson consistency: nu = E/(2G) - 1 must be a physical metal value
    // wherever a pack carries both moduli.
    for slug in [
        "brass-cartridge-c26000-mil-hdbk-698a",
        "phosphor-bronze-5a-mil-hdbk-698a",
        "music-wire-spring-moduli-fuchs-1968",
    ] {
        let e = scalar(slug, "young_modulus");
        let g = scalar(slug, "shear_modulus");
        let nu = e / (2.0 * g) - 1.0;
        assert!(
            (0.25..0.40).contains(&nu),
            "{slug}: derived nu = {nu:.3} unphysical (moduli inconsistent)"
        );
    }
    // The two phosphor-bronze sources must agree within stated spread.
    let e_nist = scalar("phosphor-bronze-c51000-nist-mono177", "young_modulus");
    let e_mil = scalar("phosphor-bronze-5a-mil-hdbk-698a", "young_modulus");
    assert!(
        ((e_nist - e_mil).abs() / e_nist) < 0.10,
        "phosphor-bronze sources disagree by more than 10%: {e_nist:.3e} vs {e_mil:.3e}"
    );

    println!(
        "{{\"suite\":\"xtask-matdb\",\"case\":\"instrument-tranche2\",\"packs\":{},\"verdict\":\"pass\"}}",
        INSTRUMENT_TRANCHE2_PACKS.len()
    );
}

/// The MIL-HDBK-5J stainless tranche (bead frankensim-0er85 tranche 1):
/// design-basis stainless data with the DIRECTION-RESOLVED
/// tension/compression asymmetry the owner named, gated by derived
/// physics rather than dims alone.
const STAINLESS_TRANCHE_PACKS: [&str; 2] = [
    "stainless-301-annealed-mil-hdbk-5j",
    "stainless-301-full-hard-mil-hdbk-5j",
];

#[test]
fn g2_cli_compiles_stainless_tranche_with_asymmetry_gates() {
    let mut packs = std::collections::BTreeMap::new();
    for slug in STAINLESS_TRANCHE_PACKS {
        let manifest = workspace_path(&format!("data/matdb/seed-v1/{slug}/manifest.tsv"));
        let out = fixture_dir().join(format!("{slug}.fsmatpk"));
        let run = run_compiler(&manifest, &out);
        assert!(
            run.status.success(),
            "{slug} refused: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let bytes = fs::read(&out).expect("read stainless pack");
        let pack = NormalizedPack::from_bytes_verified(
            NormalizedPack::from_bytes(&bytes)
                .expect("decode stainless pack")
                .content_hash(),
            &bytes,
        )
        .expect("verified stainless pack");
        assert_eq!(pack.pack_id(), slug);
        packs.insert(slug, pack);
    }
    let scalar = |slug: &str, property: &str| -> f64 {
        tonewood_scalar(&packs[slug], property)
            .unwrap_or_else(|| panic!("{slug} must carry {property}"))
    };

    // Full hard, B basis: the cold-rolling texture signature — the
    // LONGITUDINAL compressive yield sits FAR below tensile while the
    // TRANSVERSE compressive yield exceeds tensile. Sign-resolved
    // asymmetry, not a single knockdown factor.
    let fh = "stainless-301-full-hard-mil-hdbk-5j";
    let fty_l = scalar(fh, "tensile_yield_l_b_basis");
    let fcy_l = scalar(fh, "compressive_yield_l_b_basis");
    let fty_lt = scalar(fh, "tensile_yield_lt_b_basis");
    let fcy_lt = scalar(fh, "compressive_yield_lt_b_basis");
    assert!(
        fcy_l < 0.75 * fty_l,
        "longitudinal Fcy must sit far below Fty: {fcy_l:.3e} vs {fty_l:.3e}"
    );
    assert!(
        fcy_lt > 1.05 * fty_lt,
        "transverse Fcy must exceed Fty: {fcy_lt:.3e} vs {fty_lt:.3e}"
    );
    // Annealed: near-isotropic yield (no texture yet) — asymmetry
    // within 15%.
    let ann = "stainless-301-annealed-mil-hdbk-5j";
    let a_fty = scalar(ann, "tensile_yield_l_s_basis");
    let a_fcy = scalar(ann, "compressive_yield_l_s_basis");
    assert!(
        (a_fcy / a_fty - 1.0).abs() < 0.15,
        "annealed asymmetry must be small: {a_fcy:.3e} vs {a_fty:.3e}"
    );
    // Cold work strengthens: full hard ultimate far above annealed.
    assert!(
        scalar(fh, "tensile_ultimate_l_b_basis") > 2.0 * scalar(ann, "tensile_ultimate_l_s_basis"),
        "full hard must be dramatically stronger than annealed"
    );
    // Ordering sanity per condition: ultimate > yield.
    assert!(scalar(fh, "tensile_ultimate_l_b_basis") > fty_l);
    assert!(scalar(ann, "tensile_ultimate_l_s_basis") > a_fty);
    // Elastic consistency: nu = E/2G - 1 within handbook rounding of
    // the printed 0.27 (window 0.20..0.36 — the rounded moduli give
    // 0.29 annealed).
    for slug in [ann, fh] {
        let e = scalar(slug, "young_modulus_l");
        let g = scalar(slug, "shear_modulus");
        let nu = e / (2.0 * g) - 1.0;
        assert!(
            (0.20..0.36).contains(&nu),
            "{slug}: derived nu {nu:.3} inconsistent with handbook moduli"
        );
        let rho = scalar(slug, "density");
        assert!((7_800.0..8_050.0).contains(&rho), "{slug}: density {rho}");
        // Longitudinal sound speed sanity (unit-slip catcher).
        let c = (e / rho).sqrt();
        assert!(
            (4_400.0..5_400.0).contains(&c),
            "{slug}: c_L = {c:.0} m/s outside the steel window"
        );
    }
    // Elevated-temperature yield fractions: below 1 and monotone
    // decreasing (graph-read tier, but the TREND is load-bearing).
    let fractions = [
        scalar(ann, "tensile_yield_fraction_478k"),
        scalar(ann, "tensile_yield_fraction_589k"),
        scalar(ann, "tensile_yield_fraction_811k"),
        scalar(ann, "tensile_yield_fraction_1033k"),
    ];
    let mut previous = 1.0f64;
    for (i, &f) in fractions.iter().enumerate() {
        assert!(
            f > 0.0 && f < 1.0 && f <= previous,
            "yield fraction {i} must decrease below 1: {fractions:?}"
        );
        previous = f;
    }
    println!(
        "{{\"suite\":\"xtask-matdb\",\"case\":\"stainless-tranche\",\"packs\":{},\"fcy_l_over_fty_l\":{:.3},\"fcy_lt_over_fty_lt\":{:.3},\"verdict\":\"pass\"}}",
        STAINLESS_TRANCHE_PACKS.len(),
        fcy_l / fty_l,
        fcy_lt / fty_lt
    );
}

const BREADTH_TRANCHE_PACKS: [&str; 15] = [
    "stainless-17-4ph-h1025-bar-mil-hdbk-5j",
    "steel-4130-sheet-normalized-mil-hdbk-5j",
    "aluminum-2024-t3-sheet-mil-hdbk-5j",
    "lead-pure-nbs-c447",
    "tin-pure-nbs-c447",
    "zinc-pure-nbs-c447",
    "aluminum-pure-nbs-c447",
    "copper-pure-nbs-c447",
    "copper-c11000-mil-hdbk-698a",
    "brass-muntz-c28000-mil-hdbk-698a",
    "phosphor-bronze-c51000-mil-hdbk-698a",
    "brass-70-30-nbs-c447-hardness",
    "carbon-steel-cast-c033-nbs-c447",
    "sapele-fpl-gtr190",
    "teak-fpl-gtr190",
];

/// G0 actual source compiler -> canonical pack -> user-facing discovery.
/// This verifies declared source coverage, not the physical measurements.
#[test]
fn g0_cli_lead_heating_discovery_explains_actual_missing_inputs() {
    let dir = fixture_dir();
    let manifest = workspace_path("data/matdb/seed-v1/lead-pure-nbs-c447/manifest.tsv");
    let pack = dir.join("lead.fsmatpk");
    let compiled = run_compiler(&manifest, &pack);
    assert!(
        compiled.status.success(),
        "{}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    let request = dir.join("lead-heating.json");
    fs::write(
        &request,
        include_str!("../../examples/material-discovery/lead-heating.json"),
    )
    .unwrap();
    let invoke = || {
        fs_cli::run(vec![
            "--json".into(),
            "discover".into(),
            request.to_str().unwrap().into(),
            pack.to_str().unwrap().into(),
        ])
    };
    let heating = invoke();
    assert_eq!(
        heating.exit_code,
        fs_cli::exit::SUCCESS,
        "{}",
        heating.stderr
    );
    assert!(heating.stdout.contains("\"pack\":\"lead-pure-nbs-c447\""));
    assert!(heating.stdout.contains("\"kind\":\"properties\""));
    assert!(heating.stdout.contains("\"status\":\"unavailable\""));
    assert!(heating.stdout.contains("\"unknown_properties\":[\"thermal-conductivity\",\"specific-heat-capacity\",\"latent-heat\"]"));
    assert!(
        heating
            .stdout
            .contains(r#"NoClaimInDomain { property: \"density\", considered: 1 }"#),
        "{}",
        heating.stdout
    );
    assert!(
        heating
            .stdout
            .contains(r#"point: QueryPoint { axes: {\"temperature\": 650.0}"#)
    );
    // Positive existing data stays usable at its actual declared point. The
    // mixed-condition source is still explicitly an unbound property pack.
    fs::write(
        &request,
        r#"{
      "schema":"frankensim.discovery.v1", "target":"properties",
      "properties":[{"name":"density","unit":"kg/m3","kind":"dimensional"}],
      "models":[], "selection":"single-claim-only",
      "domain":{"mode":"local-state","axes":[
        {"name":"temperature","kind":"legacy","value":"293 K"}
      ]}
    }"#,
    )
    .unwrap();
    let local = invoke();
    assert_eq!(local.exit_code, fs_cli::exit::SUCCESS, "{}", local.stderr);
    assert!(local.stdout.contains("\"status\":\"complete\""));
    assert!(local.stdout.contains("\"lower_si\":11340"));
    assert!(local.stdout.contains("\"claim\":"));
    println!("{}{}", heating.stdout, local.stdout);
}

/// G2: the common-materials breadth tranches (bead frankensim-0er85) —
/// MIL-HDBK-5J structural metals, NBS/MIL melting + scale-typed
/// hardness, FPL tropical woods — compile fail-closed and pass their
/// per-tranche derived-quantity gates.
#[test]
fn g2_cli_compiles_breadth_tranches_with_derived_gates() {
    let mut packs = std::collections::BTreeMap::new();
    for slug in BREADTH_TRANCHE_PACKS {
        let manifest = workspace_path(&format!("data/matdb/seed-v1/{slug}/manifest.tsv"));
        let out = fixture_dir().join(format!("{slug}.fsmatpk"));
        let run = run_compiler(&manifest, &out);
        assert!(
            run.status.success(),
            "{slug} refused: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let bytes = fs::read(&out).expect("read breadth pack");
        let pack = NormalizedPack::from_bytes_verified(
            NormalizedPack::from_bytes(&bytes)
                .expect("decode breadth pack")
                .content_hash(),
            &bytes,
        )
        .expect("verified breadth pack");
        assert_eq!(pack.pack_id(), slug);
        packs.insert(slug, pack);
    }
    let scalar = |slug: &str, property: &str| -> f64 {
        tonewood_scalar(&packs[slug], property)
            .unwrap_or_else(|| panic!("{slug} must carry {property}"))
    };

    // TRANCHE 1 (MIL-HDBK-5J metals): E/G/nu consistency — the derived
    // nu = E/2G - 1 must sit within handbook rounding of the printed
    // Poisson ratio.
    for (slug, e_p, g_p, nu_p) in [
        (
            "stainless-17-4ph-h1025-bar-mil-hdbk-5j",
            "young_modulus",
            "shear_modulus",
            "poisson_ratio",
        ),
        (
            "steel-4130-sheet-normalized-mil-hdbk-5j",
            "young_modulus",
            "shear_modulus",
            "poisson_ratio",
        ),
        (
            "aluminum-2024-t3-sheet-mil-hdbk-5j",
            "young_modulus",
            "shear_modulus",
            "poisson_ratio",
        ),
    ] {
        let nu_derived = scalar(slug, e_p) / (2.0 * scalar(slug, g_p)) - 1.0;
        let nu_printed = scalar(slug, nu_p);
        assert!(
            (nu_derived - nu_printed).abs() < 0.06,
            "{slug}: derived nu {nu_derived:.3} vs printed {nu_printed:.3} beyond handbook rounding"
        );
    }
    // Strength ordering: yield below ultimate; 17-4PH Fcy below Fty
    // (bar, longitudinal).
    let ph = "stainless-17-4ph-h1025-bar-mil-hdbk-5j";
    assert!(scalar(ph, "tensile_yield_l_s_basis") < scalar(ph, "tensile_ultimate_l_s_basis"));
    assert!(scalar(ph, "compressive_yield_l_s_basis") < scalar(ph, "tensile_yield_l_s_basis"));
    // 2024-T3: the rolling-texture asymmetry REVERSES with direction
    // (Fcy_L < Fty_L but Fcy_LT > Fty_LT) — same signature class the
    // 301 tranche pinned; a single knockdown factor cannot represent
    // this.
    let al = "aluminum-2024-t3-sheet-mil-hdbk-5j";
    assert!(
        scalar(al, "compressive_yield_l_a_basis") < scalar(al, "tensile_yield_l_a_basis"),
        "2024-T3 longitudinal Fcy must sit below Fty"
    );
    assert!(
        scalar(al, "compressive_yield_lt_a_basis") > scalar(al, "tensile_yield_lt_a_basis"),
        "2024-T3 transverse Fcy must exceed Fty"
    );

    // TRANCHE 2 (melting + hardness): solidus never above liquidus for
    // every alloy pack that carries both.
    for slug in [
        "copper-c11000-mil-hdbk-698a",
        "brass-muntz-c28000-mil-hdbk-698a",
        "phosphor-bronze-c51000-mil-hdbk-698a",
    ] {
        assert!(
            scalar(slug, "melting_point_solidus") <= scalar(slug, "melting_point_liquidus"),
            "{slug}: solidus above liquidus"
        );
    }
    // Cross-source melting pin: C11000 liquidus (MIL-HDBK-698A, deg F)
    // vs elemental copper (NBS C447, deg C) agree within 1 K.
    let cu_mil = scalar("copper-c11000-mil-hdbk-698a", "melting_point_liquidus");
    let cu_nbs = scalar("copper-pure-nbs-c447", "melting_point");
    assert!(
        (cu_mil - cu_nbs).abs() < 1.0,
        "cross-source copper melting disagreement: {cu_mil:.2} vs {cu_nbs:.2} K"
    );
    // Alloying depresses melting: both brasses melt below pure copper.
    assert!(scalar("brass-muntz-c28000-mil-hdbk-698a", "melting_point_liquidus") < cu_nbs);
    // Melting ordering across the pure metals (the extreme-regime
    // doctrine anchor: lead lowest of the structural set).
    let melt = |s: &str| scalar(s, "melting_point");
    assert!(melt("tin-pure-nbs-c447") < melt("lead-pure-nbs-c447"));
    assert!(melt("lead-pure-nbs-c447") < melt("zinc-pure-nbs-c447"));
    assert!(melt("zinc-pure-nbs-c447") < melt("aluminum-pure-nbs-c447"));
    assert!(melt("aluminum-pure-nbs-c447") < melt("copper-pure-nbs-c447"));
    // SCALE TYPING is structural: every hardness property name carries
    // its scale, and no scale-less "hardness" scalar exists anywhere in
    // the tranche.
    for (slug, pack) in &packs {
        for (_id, claim) in pack.claims().claims_ordered() {
            let name = claim.key.name();
            if name.contains("hardness") && !name.contains("side_hardness") {
                assert!(
                    name.contains("brinell")
                        || name.contains("rockwell")
                        || name.contains("vickers"),
                    "{slug}: hardness claim {name} is not scale-typed"
                );
            }
            assert_ne!(name, "hardness", "{slug}: scale-less hardness forbidden");
        }
    }
    // Cold work raises hardness within each lot (copper and brass).
    assert!(
        scalar("copper-pure-nbs-c447", "hardness_brinell_cold_drawn_56pct")
            > scalar("copper-pure-nbs-c447", "hardness_brinell_annealed")
    );
    assert!(
        scalar(
            "brass-70-30-nbs-c447-hardness",
            "hardness_rockwell_b_cold_worked_37pct"
        ) > scalar(
            "brass-70-30-nbs-c447-hardness",
            "hardness_rockwell_b_soft_sheet"
        )
    );
    assert!(
        scalar(
            "brass-70-30-nbs-c447-hardness",
            "hardness_brinell_cold_rolled_11pct"
        ) > scalar("brass-70-30-nbs-c447-hardness", "hardness_brinell_annealed")
    );

    // TRANCHE 3 (FPL woods): plausibility gates — MOR far below MOE,
    // compression parallel below MOR, positive toughness.
    for slug in ["sapele-fpl-gtr190", "teak-fpl-gtr190"] {
        let mor = scalar(slug, "modulus_of_rupture");
        let moe = scalar(slug, "bending_modulus_of_elasticity") * 1.0e3; // GPa -> MPa
        assert!(moe > 50.0 * mor, "{slug}: MOE/MOR ratio implausible");
        assert!(scalar(slug, "compression_parallel_max") < mor);
        assert!(scalar(slug, "work_to_maximum_load") > 0.0);
        assert!(scalar(slug, "side_hardness") > 0.0);
    }

    println!(
        "{{\"suite\":\"matdb-pack\",\"case\":\"breadth-tranches\",\"packs\":{},\"verdict\":\"pass\"}}",
        packs.len()
    );
}

/// G2: the 301-annealed elevated-temperature yield-fraction data is
/// carried as a queryable temperature CURVE (bead frankensim-0er85's
/// curve-claim acceptance), decreasing across the read band.
#[test]
fn g2_stainless_yield_fraction_temperature_curve() {
    let slug = "stainless-301-annealed-mil-hdbk-5j";
    let manifest = workspace_path(&format!("data/matdb/seed-v1/{}/manifest.tsv", slug));
    let out = fixture_dir().join(format!("{}-curve.fsmatpk", slug));
    let run = run_compiler(&manifest, &out);
    assert!(
        run.status.success(),
        "refused: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = fs::read(&out).expect("read pack");
    let pack = NormalizedPack::from_bytes_verified(
        NormalizedPack::from_bytes(&bytes)
            .expect("decode pack")
            .content_hash(),
        &bytes,
    )
    .expect("verified pack");
    let claims = pack.claims().claims_for("tensile_yield_fraction");
    assert_eq!(claims.len(), 1, "one curve claim expected");
    let claim = claims[0].1;
    let PropertyValue::Curve {
        abscissa, knots, ..
    } = &claim.value
    else {
        panic!("yield fraction must be a CURVE claim");
    };
    assert_eq!(abscissa, "temperature");
    assert_eq!(knots.len(), 4);
    // Strength decreases with temperature across the whole read band —
    // the monotonicity gate the bead names.
    for pair in knots.windows(2) {
        assert!(pair[0].0 < pair[1].0, "temperature knots must ascend");
        assert!(
            pair[0].1 > pair[1].1,
            "yield fraction must DECREASE with temperature: {:?}",
            knots
        );
    }
    // Exact endpoints as graph-read: 478 K -> 0.68, 1033 K -> 0.43.
    assert!((knots[0].0 - 478.0).abs() < 1.0e-9 && (knots[0].1 - 0.68).abs() < 1.0e-12);
    assert!((knots[3].0 - 1033.0).abs() < 1.0e-9 && (knots[3].1 - 0.43).abs() < 1.0e-12);
    // QUERIED through the evaluator (the bead's "landed and queried"):
    // tabulated-only answers at an exact knot temperature.
    let point = fs_matdb::QueryPoint::new()
        .with("temperature", 589.0)
        .expect("query point");
    let answer = pack
        .claims()
        .query(
            "tensile_yield_fraction",
            &point,
            fs_matdb::SelectionPolicy::SingleClaimOnly,
        )
        .expect("curve query at a knot");
    assert!((answer.evidence.value.value - 0.62).abs() < 1.0e-12);
    println!("{{\"suite\":\"matdb-pack\",\"case\":\"fty-fraction-curve\",\"verdict\":\"pass\"}}");
}
