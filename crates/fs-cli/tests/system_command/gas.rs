//! Source-backed stream scenes retain native calorics and inventory accounts.

use fs_gas::component::{ClosedGasParcel, GasEnergySolve};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::gas_species::input::packs::{SourcePackBytes, read_source_parcel};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;

use super::{element, example, run_binary, scene_context, scratch, value};

const PACKS: [&str; 3] = ["species.pack", "nasa9.pack", "transport.fsmcdpk"];
const CONDUCTANCE: f64 = 2.0e-9;
const HEAT_W: f64 = 40.0;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual} vs {expected}, tolerance {tolerance}"
    );
}

fn probe(scene: &fs_cli::SystemScene, coordinate: &str) -> f64 {
    scene.probe(coordinate).unwrap().unwrap()
}

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 900.0,
        energy_tolerance_j_kg: 1.0e-7,
        maximum_iterations: 128,
    }
}

fn nitrogen_packs() -> [Vec<u8>; 3] {
    PACKS.map(|name| std::fs::read(example(&format!("n2/{name}"))).unwrap())
}

fn parcel(record: &str, packs: &[[Vec<u8>; 3]]) -> ClosedGasParcel {
    let sources: Vec<_> = packs
        .iter()
        .map(|triple| SourcePackBytes {
            species: &triple[0],
            nasa9: &triple[1],
            transport: &triple[2],
        })
        .collect();
    ClosedGasParcel::from_isochoric(
        read_source_parcel(
            record,
            &sources,
            None,
            ThermochemicalPhase::FrozenGas,
            1.0e-3,
        )
        .unwrap(),
    )
}

fn scene_file(text: &str) -> std::path::PathBuf {
    let path = scratch("n2_stream_exchange.json", text);
    write_packs(&path, "n2", &nitrogen_packs());
    path
}

fn write_packs(scene: &std::path::Path, directory: &str, packs: &[Vec<u8>; 3]) {
    let base = scene.parent().unwrap().join(directory);
    std::fs::create_dir_all(&base).unwrap();
    for (name, bytes) in PACKS.into_iter().zip(packs) {
        std::fs::write(base.join(name), bytes).unwrap();
    }
}

/// An explicitly authored association for numerical tests. The tracer uses
/// nitrogen's original molar calorics and transport with twice its molar mass;
/// it is not a measured second chemical species or an oxygen/air model.
fn authored_species(packs: &[Vec<u8>; 3], id: &str, mass_scale: f64) -> Vec<u8> {
    use fs_matdb::{
        NormalizedSpeciesPack, Provenance, SPECIES_MOLAR_MASS_DIMS, SPECIES_PACK_TARGET_BASIS,
        SPECIES_REFERENCE_PRESSURE_DIMS, SpeciesAssociation, SpeciesNormalizationReceipt,
        SpeciesNormalizationTarget,
    };
    let original = NormalizedSpeciesPack::from_bytes(&packs[0]).unwrap();
    let association = original.association();
    let mass = mass_scale * association.molar_mass();
    let source = fs_blake3::hash_domain(
        "org.frankensim.fs-cli.streaming-gas-test.v1",
        format!("{id};molar_mass_kg_mol={mass:.17e}").as_bytes(),
    );
    let declared = SpeciesAssociation::new(
        fs_qty::chemistry::SpeciesId::new(id).unwrap(),
        mass,
        association.standard_state_phase(),
        association.reference_eos(),
        association.reference_pressure(),
        association.elemental_reference(),
        vec![source],
        Provenance {
            source: "authored numerical stream test association; no measured-data claim".into(),
            license: "CC0: authored numerical fixture".into(),
            artifact: Some(source),
        },
    )
    .unwrap();
    let receipts = [
        (
            SpeciesNormalizationTarget::MolarMass,
            SPECIES_MOLAR_MASS_DIMS,
            "kg/mol",
        ),
        (
            SpeciesNormalizationTarget::ReferencePressure,
            SPECIES_REFERENCE_PRESSURE_DIMS,
            "Pa",
        ),
    ]
    .into_iter()
    .map(|(target, dims, basis)| {
        SpeciesNormalizationReceipt::new(
            target,
            source,
            dims,
            1.0,
            0.0,
            basis,
            SPECIES_PACK_TARGET_BASIS,
        )
    })
    .collect();
    NormalizedSpeciesPack::new(
        id,
        "fs-cli streaming tests",
        source,
        "CC0 authored association",
        declared,
        receipts,
    )
    .unwrap()
    .to_bytes()
}

fn tracer_packs(nitrogen: &[Vec<u8>; 3]) -> [Vec<u8>; 3] {
    let original = fs_matdb::NormalizedModelPack::from_bytes(&nitrogen[1]).unwrap();
    let calorics = fs_matdb::NormalizedModelPack::new(
        "tracer",
        original.compiler(),
        original.source_artifact(),
        original.redistribution_terms(),
        original.models().to_vec(),
        original.normalizations().to_vec(),
    )
    .unwrap();
    [
        authored_species(nitrogen, "tracer", 2.0),
        calorics.to_bytes(),
        nitrogen[2].clone(),
    ]
}

/// The oracle uses the native paired-transfer entry point, independent of
/// GasCavity's species-boundary implementation, then the native heat solve.
fn compare_native(
    text: &str,
    path: &std::path::Path,
    mut native: [ClosedGasParcel; 2],
) -> (fs_cli::SystemScene, [u32; 2]) {
    let species: Vec<_> = native[0]
        .model()
        .thermodynamics()
        .components()
        .iter()
        .map(|model| model.species().clone())
        .collect();
    let mass_before: f64 = native.iter().map(ClosedGasParcel::mass_kg).sum();
    let species_before: Vec<f64> = species
        .iter()
        .map(|id| {
            native
                .iter()
                .map(|gas| gas.species_mass_kg(id).unwrap())
                .sum()
        })
        .collect();
    let energy_before: f64 = native
        .iter()
        .map(|gas| gas.mass_kg() * gas.state().specific_internal_energy_j_kg())
        .sum();
    let solve = GasMassTransferSolve {
        donor: controls(),
        receiver: controls(),
        mass_tolerance_kg: 1.0e-15,
        energy_tolerance_j: 1.0e-9,
    };
    let mut scene = fs_cli::SystemScene::from_text(text, path).unwrap();
    let mut directions = [0; 2];
    scene_context(|cx, _| {
        while !scene.is_finished() {
            let rate = CONDUCTANCE
                * (native[0].state().state().pressure - native[1].state().state().pressure);
            let before = scene.time_s();
            assert!(scene.advance(f64::MAX, cx).unwrap());
            let dt = scene.time_s() - before;
            if rate != 0.0 {
                let from = usize::from(rate < 0.0);
                let to = 1 - from;
                directions[from] += 1;
                let candidate = native[from]
                    .transfer_mass_to(&native[to], dt * rate.abs(), solve, || false)
                    .unwrap();
                let (donor, receiver, _) = candidate.into_parts();
                native[from] = donor;
                native[to] = receiver;
            }
            native[1] = native[1]
                .discrete_gradient_heat_volume(
                    native[1].volume_m3(),
                    HEAT_W * dt,
                    controls(),
                    || false,
                )
                .unwrap()
                .into_parts()
                .0;
            close(probe(&scene, "a_valve.mass_rate_kg_s"), rate, 2.0e-11);
            for (name, expected) in ["left", "right"].into_iter().zip(&native) {
                for (field, value, tolerance) in [
                    ("mass_kg", expected.mass_kg(), 2.0e-15),
                    ("volume_m3", expected.volume_m3(), 0.0),
                    (
                        "temperature_k",
                        expected.state().state().temperature,
                        2.0e-6,
                    ),
                    ("pressure_pa", expected.state().state().pressure, 1.0e-3),
                ] {
                    close(probe(&scene, &format!("{name}.{field}")), value, tolerance);
                }
                for (index, id) in species.iter().enumerate() {
                    let mass = expected.species_mass_kg(id).unwrap();
                    close(
                        probe(&scene, &format!("{name}.species{index}_mass_kg")),
                        mass,
                        2.0e-15,
                    );
                    close(
                        probe(&scene, &format!("{name}.mass_fraction{index}")),
                        mass / expected.mass_kg(),
                        2.0e-12,
                    );
                }
            }
        }
    });
    close(
        probe(&scene, "left.mass_kg") + probe(&scene, "right.mass_kg"),
        mass_before,
        2.0e-15,
    );
    for (index, initial) in species_before.into_iter().enumerate() {
        close(
            probe(&scene, &format!("left.species{index}_mass_kg"))
                + probe(&scene, &format!("right.species{index}_mass_kg")),
            initial,
            2.0e-15,
        );
        close(
            scene.ledger().carried_residual(index).unwrap(),
            0.0,
            2.0e-15,
        );
    }
    let energy_after: f64 = ["left", "right"]
        .into_iter()
        .map(|name| {
            probe(&scene, &format!("{name}.mass_kg"))
                * probe(&scene, &format!("{name}.specific_internal_energy_j_kg"))
        })
        .sum();
    close(
        energy_after - energy_before,
        HEAT_W * scene.time_s(),
        2.0e-6,
    );
    close(scene.ledger().mass_residual_kg().unwrap(), 0.0, 2.0e-15);
    close(scene.ledger().energy_residual_j(), 0.0, 2.0e-6);
    for account in &scene.ledger().nodes {
        assert!(account.has_streams);
        assert_eq!(account.mass.as_ref().unwrap().audited_steps, scene.steps());
        assert_eq!(account.carried.len(), species.len());
        for carried in &account.carried {
            assert_eq!(carried.as_ref().unwrap().audited_steps, scene.steps());
        }
    }
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(json, scene.render_json());
    (scene, directions)
}

#[test]
fn source_backed_exchange_matches_native_transfer_and_reverses_with_heat() {
    let path = std::path::PathBuf::from(example("n2_stream_exchange.json"));
    let text = std::fs::read_to_string(&path).unwrap();
    let packs = [nitrogen_packs()];
    let initial = [
        parcel("thermochemical 300 90000 mole 1 N2 1", &packs),
        parcel("thermochemical 300 80000 mole 1 N2 1", &packs),
    ];
    let (scene, directions) = compare_native(&text, &path, initial);
    assert!(directions.iter().all(|count| *count > 0), "{directions:?}");
    assert!(probe(&scene, "right.pressure_pa") > probe(&scene, "left.pressure_pa"));
    assert!(scene.render_json().contains("\"stream_species\":[\"N2\"]"));
    close(
        value(&scene.render_json(), &element("right"), "mass_fraction0"),
        1.0,
        1.0e-14,
    );
}

#[test]
fn source_backed_mixtures_keep_each_species_mass_and_declared_coordinate() {
    let original = std::fs::read_to_string(example("n2_stream_exchange.json")).unwrap();
    let left = "thermochemical 300 90000 mole 2 N2 0.75 tracer 0.25";
    let right = "thermochemical 300 80000 mole 2 N2 0.25 tracer 0.75";
    let text = original
        .replace("thermochemical 300 90000 mole 1 N2 1", left)
        .replace("thermochemical 300 80000 mole 1 N2 1", right)
        .replace(
            "\"transport\": \"n2/transport.fsmcdpk\"}]",
            "\"transport\": \"n2/transport.fsmcdpk\"}, {\"species\": \"tracer/species.pack\", \
             \"nasa9\": \"tracer/nasa9.pack\", \"transport\": \"tracer/transport.fsmcdpk\"}]",
        );
    let path = scene_file(&text);
    let nitrogen = nitrogen_packs();
    let tracer = tracer_packs(&nitrogen);
    write_packs(&path, "tracer", &tracer);
    let packs = [nitrogen, tracer];
    let initial = [parcel(left, &packs), parcel(right, &packs)];
    let starting_fraction = initial[1]
        .species_mass_kg(initial[1].model().thermodynamics().components()[0].species())
        .unwrap()
        / initial[1].mass_kg();
    let (scene, _) = compare_native(&text, &path, initial);
    assert!(
        scene
            .render_json()
            .contains("\"stream_species\":[\"N2\",\"tracer\"]")
    );
    assert!(probe(&scene, "right.mass_fraction0") > starting_fraction);
}

#[test]
fn stream_scene_admission_preserves_source_identity_and_requires_explicit_policy() {
    let original = std::fs::read_to_string(example("n2_stream_exchange.json")).unwrap();
    let path = scene_file(&original);
    for (changed, reason) in [
        (
            original.replace(
                "\"stream_mass_tolerance_kg\": 1.0e-15",
                "\"stream_mass_tolerance_kg\": 0",
            ),
            "stream_mass_tolerance_kg",
        ),
        (
            original.replace(
                "\"stream_mass_tolerance_kg\":",
                "\"insulated\": true, \"stream_mass_tolerance_kg\":",
            ),
            "insulated",
        ),
        (
            original.replace(
                "\"stream_basis\": \"left.stream\"",
                "\"stream_basis\": \"left.p\"",
            ),
            "stream",
        ),
        (
            original.replace(
                "\"stream_basis\": \"left.stream\"",
                "\"stream_basis\": \"missing.stream\"",
            ),
            "missing",
        ),
        (
            original.replace(
                "\"stream_basis\": \"left.stream\"",
                "\"stream_basis\": \"left.stream\", \"carried\": 1",
            ),
            "carried",
        ),
        (
            original.replace("\"stream_basis\": \"left.stream\"", "\"carried\": 1"),
            "declared medium",
        ),
    ] {
        assert_ne!(changed, original);
        let error = fs_cli::SystemScene::from_text(&changed, &path)
            .err()
            .unwrap();
        assert!(error.message.contains(reason), "expected {reason}: {error}");
    }

    // Same species name and stream width, but a different admitted molar
    // mass: the explicit native basis must refuse instead of erasing it.
    let mut changed_packs = nitrogen_packs();
    changed_packs[0] = authored_species(&changed_packs, "N2", 1.01);
    write_packs(&path, "other", &changed_packs);
    let split = original.find("\"right\":").unwrap();
    let changed = format!(
        "{}{}",
        &original[..split],
        original[split..].replace("n2/", "other/")
    );
    let error = fs_cli::SystemScene::from_text(&changed, &path)
        .err()
        .unwrap();
    assert!(error.message.contains("declared medium"), "{error}");
}

#[test]
fn streaming_checkpoint_continuation_binds_mass_tolerance_and_source() {
    let text = std::fs::read_to_string(example("n2_stream_exchange.json")).unwrap();
    let path = scene_file(&text);
    // A real process resumes evolving native mass and source-owned calorics;
    // the existing helper compares its entire report with an uninterrupted
    // native library run that takes the same physical step at the pause.
    let root = super::checkpoints::checkpoint_continuation_at(&path, &text, 0.375);
    let changed = text.replace(
        "\"stream_mass_tolerance_kg\": 1.0e-15",
        "\"stream_mass_tolerance_kg\": 2.0e-15",
    );
    assert_ne!(changed, text);
    assert!(fs_cli::SystemScene::from_text(&changed, &path).is_ok());
    std::fs::write(&path, &changed).unwrap();
    super::checkpoints::resume_is_refused("stream mass tolerance", &path, &root);

    // Both endpoints now read the same changed source, so this remains a
    // valid connected scene. Its saved state belongs to the old source.
    std::fs::write(&path, &text).unwrap();
    let mut packs = nitrogen_packs();
    packs[0] = authored_species(&packs, "N2", 1.01);
    write_packs(&path, "n2", &packs);
    assert!(fs_cli::SystemScene::from_text(&text, &path).is_ok());
    super::checkpoints::resume_is_refused("stream gas source", &path, &root);
}

fn stream_accounts(species: usize, tight: impl Fn(&str) -> bool) -> Vec<(String, f64)> {
    let mut accounts = Vec::new();
    for cavity in ["left", "right"] {
        let mut coordinates = vec![
            ("volume_m3".to_owned(), 1.0e-3),
            ("temperature_k".to_owned(), 300.0),
            ("mass_kg".to_owned(), 1.0e-3),
            ("energy_j".to_owned(), 100.0),
        ];
        coordinates.extend((0..species).map(|index| (format!("species{index}_mass_kg"), 1.0e-3)));
        for (coordinate, scale) in coordinates {
            let name = format!("{cavity}.{coordinate}");
            let fraction = if tight(&name) { 1.0e-6 } else { 1.0e-2 };
            accounts.push((name, fraction * scale));
        }
    }
    accounts
}

fn adaptive_stream(text: &str, accounts: &[(String, f64)]) -> String {
    let borrowed: Vec<_> = accounts
        .iter()
        .map(|(name, value)| (name.as_str(), *value))
        .collect();
    let settings = super::error_controlled(0.1, &borrowed).replace("\"order\": 2", "\"order\": 1");
    assert!(text.contains("\"step_s\": 0.01,"));
    text.replace("\"step_s\": 0.01,", &settings)
}

/// Fixed-volume mass and enthalpy balances, with NASA-9 coefficients from
/// the original source TSVs. This calls no Component or parcel integrator.
/// The common formation constant cancels from h_d - u_i for pure N2.
fn n2_exchange_reference(steps: u32) -> [f64; 4] {
    let a = super::n2_nasa9_low();
    let r = super::n2_gas_constant();
    let volume = 1.0e-3;
    let mut state = [
        90_000.0 * volume / (r * 300.0),
        300.0,
        80_000.0 * volume / (r * 300.0),
        300.0,
    ];
    let rhs = |y: [f64; 4]| {
        let rate = CONDUCTANCE * r * (y[0] * y[1] - y[2] * y[3]) / volume;
        let donor_temperature = if rate >= 0.0 { y[1] } else { y[3] };
        let enthalpy = r * (super::n2_energy(&a, donor_temperature) + donor_temperature);
        let left_u = r * super::n2_energy(&a, y[1]);
        let right_u = r * super::n2_energy(&a, y[3]);
        [
            -rate,
            rate * (left_u - enthalpy)
                / (y[0] * r * (super::n2_heat_capacity(&a, y[1]) - 1.0)),
            rate,
            (rate * (enthalpy - right_u) + HEAT_W)
                / (y[2] * r * (super::n2_heat_capacity(&a, y[3]) - 1.0)),
        ]
    };
    let dt = 1.0 / f64::from(steps);
    for _ in 0..steps {
        let first = rhs(state);
        let second = rhs(std::array::from_fn(|i| state[i] + 0.5 * dt * first[i]));
        let third = rhs(std::array::from_fn(|i| state[i] + 0.5 * dt * second[i]));
        let fourth = rhs(std::array::from_fn(|i| state[i] + dt * third[i]));
        for i in 0..4 {
            state[i] += dt * (first[i] + 2.0 * second[i] + 2.0 * third[i] + fourth[i]) / 6.0;
        }
    }
    state
}

#[test]
fn adaptive_source_backed_exchange_converges_to_independent_caloric_balances() {
    let original = std::fs::read_to_string(example("n2_stream_exchange.json")).unwrap();
    let path = scene_file(&original);
    let run = |tight| {
        let text = adaptive_stream(&original, &stream_accounts(1, |_| tight));
        let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        scene_context(|cx, _| scene.run(cx).unwrap());
        scene
    };
    let loose = run(false);
    let tight = run(true);
    let reference = n2_exchange_reference(8192);
    let finer = n2_exchange_reference(16384);
    for (index, (actual, expected)) in reference.into_iter().zip(finer).enumerate() {
        close(actual, expected, if index % 2 == 0 { 1.0e-11 } else { 1.0e-5 });
    }
    let distance = |scene: &fs_cli::SystemScene| {
        [
            "left.mass_kg",
            "left.temperature_k",
            "right.mass_kg",
            "right.temperature_k",
        ]
        .into_iter()
        .zip(finer)
        .map(|(name, expected)| (probe(scene, name) / expected - 1.0).abs())
        .fold(0.0, f64::max)
    };
    let (far, near) = (distance(&loose), distance(&tight));
    assert!(near < 0.2 * far && near < 2.0e-4, "loose {far:e}, tight {near:e}");
    assert!(tight.steps() > loose.steps());
    assert!(probe(&tight, "a_valve.mass_rate_kg_s") < 0.0);
    close(
        probe(&tight, "left.mass_kg") + probe(&tight, "right.mass_kg"),
        170_000.0 * 1.0e-3 / (super::n2_gas_constant() * 300.0),
        2.0e-14,
    );
    close(tight.ledger().mass_residual_kg().unwrap(), 0.0, 2.0e-14);
    close(tight.ledger().carried_residual(0).unwrap(), 0.0, 2.0e-14);
    close(tight.ledger().energy_residual_j(), 0.0, 2.0e-6);

    // The committed example takes this same path through the actual binary.
    let example_path = example("adaptive_n2_stream_exchange.json");
    let text = std::fs::read_to_string(&example_path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &example_path).unwrap();
    scene_context(|cx, _| scene.run(cx).unwrap());
    let (code, json, stderr) = run_binary(&["--json", "system", &example_path]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(json, scene.render_json());
}

#[test]
fn adaptive_mixture_species_accounts_control_refinement_and_resume_exactly() {
    let original = std::fs::read_to_string(example("n2_stream_exchange.json")).unwrap();
    let original = original
        .replace(
            "thermochemical 300 90000 mole 1 N2 1",
            "thermochemical 300 90000 mole 2 N2 0.75 tracer 0.25",
        )
        .replace(
            "thermochemical 300 80000 mole 1 N2 1",
            "thermochemical 300 80000 mole 2 N2 0.25 tracer 0.75",
        )
        .replace(
            "\"transport\": \"n2/transport.fsmcdpk\"}]",
            "\"transport\": \"n2/transport.fsmcdpk\"}, {\"species\": \"tracer/species.pack\", \
             \"nasa9\": \"tracer/nasa9.pack\", \"transport\": \"tracer/transport.fsmcdpk\"}]",
        );
    let path = scene_file(&original);
    write_packs(&path, "tracer", &tracer_packs(&nitrogen_packs()));
    let loose_text = adaptive_stream(&original, &stream_accounts(2, |_| false));
    let text = adaptive_stream(&original, &stream_accounts(2, |name| name.contains(".species")));
    let mut loose = fs_cli::SystemScene::from_text(&loose_text, &path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    let totals = |scene: &fs_cli::SystemScene| {
        ["mass_kg", "species0_mass_kg", "species1_mass_kg", "energy_j"].map(|coordinate| {
            probe(scene, &format!("left.{coordinate}")) + probe(scene, &format!("right.{coordinate}"))
        })
    };
    let initial = totals(&scene);
    scene_context(|cx, _| {
        loose.run(cx).unwrap();
        assert!(scene.advance(f64::MAX, cx).unwrap());
        assert!(!scene.is_finished());
        let (bytes, root) = scene.checkpoint().unwrap();
        let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        scene.run(cx).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());
        assert_eq!(resumed.checkpoint().unwrap(), scene.checkpoint().unwrap());
    });
    assert!(
        scene.steps() > loose.steps(),
        "species tolerances did not control the step"
    );
    let final_totals = totals(&scene);
    for index in 0..3 {
        close(final_totals[index], initial[index], 2.0e-14);
    }
    close(final_totals[3] - initial[3], HEAT_W * scene.time_s(), 2.0e-6);
    for index in 0..2 {
        close(scene.ledger().carried_residual(index).unwrap(), 0.0, 2.0e-14);
    }
    close(scene.ledger().mass_residual_kg().unwrap(), 0.0, 2.0e-14);
    close(scene.ledger().energy_residual_j(), 0.0, 2.0e-6);
    for name in ["left", "right"] {
        close(
            probe(&scene, &format!("{name}.energy_j")),
            probe(&scene, &format!("{name}.mass_kg"))
                * probe(&scene, &format!("{name}.specific_internal_energy_j_kg")),
            1.0e-12,
        );
    }
    assert!(scene.render_json().contains("\"stream_species\":[\"N2\",\"tracer\"]"));

    // Every state coordinate is mandatory, including a later catalog species.
    let accounts = stream_accounts(2, |_| true);
    for missing in 0..accounts.len() {
        let mut partial = accounts.clone();
        let (name, _) = partial.remove(missing);
        let error = fs_cli::SystemScene::from_text(&adaptive_stream(&original, &partial), &path)
            .err()
            .unwrap();
        assert!(
            error.message.contains(&format!("adaptive.tolerances.{name}")),
            "{error}"
        );
    }
    let order_two = text.replace("\"order\": 1", "\"order\": 2");
    let error = fs_cli::SystemScene::from_text(&order_two, &path).err().unwrap();
    assert!(
        error.message.contains("adaptive.order") && error.message.contains("order 1"),
        "{error}"
    );
}
