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
            original.replace(
                "\"duration_s\": 1.0,",
                "\"duration_s\": 1.0, \"adaptive\": {\"tolerances\": {}},",
            ),
            "fixed step_s",
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
