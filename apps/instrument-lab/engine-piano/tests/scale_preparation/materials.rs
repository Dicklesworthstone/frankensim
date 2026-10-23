//! Manufactured bonded rods exercise canonical sources through the original
//! section, tuning, modal bank and live receiver. No measured alloy is implied.
use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_live_piano::{
    preparation::materials::{COLUMNS, StringMaterials, TensionPolicy, operating_point},
    refmods::{bank_string_losses, linear::Bank},
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MATDB_PACK_TARGET_BASIS, MaterialStateId, NormalizationReceipt,
    NormalizationTarget, NormalizedMaterialCardPack, NormalizedPack, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, StatisticMember, UncertaintyModel,
};
use fs_material::{
    gas::{GasSpec, GasState},
    state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY},
};
use fs_qty::{
    Density, Dims, DynViscosity, Frequency, Pressure, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use std::{
    f64::consts::{PI, TAU},
    sync::Arc,
};

const ETA: &str = "declared_bending_viscosity";
const LIMIT: &str = "declared_tensile_limit";
// Manufactured limit above every fixture stress unless a test replaces it.
const LIMIT_PA: f64 = 5e9;
const INNER: f64 = 0.00035;
const OUTER: f64 = 0.0005;
const CORE: [f64; 3] = [8000.0, 2e11, 1e8];
const COATING: [f64; 3] = [10000.0, 1e11, 4e8];
const MAX_HZ: f64 = 1e6;
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point() -> fs_matdb::QueryPoint {
    // Use the actual browser/native helper with source-defined axis names.
    operating_point("bath_kelvin", 293.15, "oscillation_hz", 1000.0).unwrap()
}
fn validity(lo: f64, hi: f64) -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity(
            "bath_kelvin",
            typed(QuantityKind::AbsoluteTemperature),
            280.0,
            320.0,
        )
        .with_quantity(
            "oscillation_hz",
            typed(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
            lo,
            hi,
        )
}
fn pack(values: [f64; 3]) -> NormalizedMaterialCardPack {
    pack_with(values, |_| {})
}
fn pack_with(
    values: [f64; 3],
    mut edit: impl FnMut(&mut PropertyClaim),
) -> NormalizedMaterialCardPack {
    let artifact = hash_bytes(format!("manufactured SI coefficients {values:?}").as_bytes());
    let provenance = Provenance {
        source: "authored bonded-rod fixture".into(),
        license: "test fixture; no measured evidence".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "manufactured concentric rod".into(),
            method: "analytic test data".into(),
            artifact,
            caveats: "not an alloy calibration".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let mut receipts = Vec::new();
    for ((name, dims), value) in [
        (DENSITY_PROPERTY, Density::DIMS),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
        (ETA, DynViscosity::DIMS),
        (LIMIT, Pressure::DIMS),
    ]
    .into_iter()
    .zip(values.into_iter().chain([LIMIT_PA]))
    {
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: validity(0.0, MAX_HZ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance.clone(),
        };
        edit(&mut claim);
        let scalar = match &claim.value {
            PropertyValue::Scalar { value, dims } => Some((*value, *dims)),
            _ => None,
        };
        let id = claims.insert_claim(claim).unwrap();
        if let Some((value, dims)) = scalar {
            receipts.push(NormalizationReceipt::new(
                NormalizationTarget::ClaimValue(StatisticMember::scalar(id)),
                hash_bytes(value.to_string().as_bytes()),
                dims,
                1.0,
                0.0,
                MATDB_PACK_TARGET_BASIS,
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ));
        }
    }
    let claims = NormalizedPack::new(
        "manufactured-bonded-rod",
        "test-v1",
        artifact,
        "test redistribution",
        claims,
        vec![],
        receipts,
    )
    .unwrap();
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
    )
    .unwrap()
}
fn small_description() -> Description {
    let mut desc = description();
    desc.courses.retain(|c| [57, 69, 81].contains(&c.midi));
    for c in &mut desc.courses {
        c.unison = 1;
        c.detune_cents = 0.0;
        c.duplex_length_m = 0.0;
    }
    desc
}
fn construction(tension: &str, strength: &str, rows: &str) -> String {
    format!("frankensim-piano-materials-v2\ntension,{tension}\nstrength,{strength}\n{COLUMNS}\n{rows}")
}
fn section_geometry(courses: &[Course]) -> String {
    let mut rows = String::new();
    for course in courses {
        rows.push_str(&format!(
            "{},core,{INNER},,,0\n{},core,{OUTER},,,1\n",
            course.midi, course.midi
        ));
    }
    construction("installed-tension", LIMIT, &rows)
}
fn load(
    courses: &[Course],
    core: &NormalizedMaterialCardPack,
    coating: &NormalizedMaterialCardPack,
) -> Arc<StringMaterials> {
    Arc::new(
        StringMaterials::from_packs(
            &[&core.to_bytes(), &coating.to_bytes()],
            &section_geometry(courses),
            &point(),
            ETA,
        )
        .unwrap(),
    )
}
fn config() -> Config {
    let mut cfg = Config::live(48_000);
    cfg.max_partials = 12;
    cfg.course_mode_budget = 12;
    cfg.max_dynamic_modes = 36;
    cfg.max_courses = 3;
    cfg.duplex = false;
    cfg.partial_band_hz = 6000.0;
    cfg.silent_bank_hz = 0.0;
    cfg.rt0425_string_damping = false;
    cfg
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

#[test]
fn material_geometry_and_the_same_atmosphere_determine_live_cylinder_drag() {
    let mut desc = small_description();
    let source = load(&desc.courses, &pack(CORE), &pack(COATING));
    preparation::prepare_with_materials(&mut desc, None, None, Some(&source)).unwrap();
    let gases = [
        preparation::resolve_atmosphere("293.15 50000 0").unwrap(),
        preparation::resolve_atmosphere("293.15 100000 0").unwrap(),
        preparation::resolve_atmosphere("mixture 293.15 100000 mole 1 CO2 1").unwrap(),
    ];
    let mut prepared = Vec::new();
    for atmosphere in gases {
        let piano = LivePiano::new_with_materials(&desc, config(), atmosphere, source.clone()).unwrap();
        assert_eq!(piano.acoustic_gas(), atmosphere);
        for course in piano.courses() {
            let frequencies = piano.partials_hz(course.midi);
            let rates = piano.string_gas_drag_rates_per_s(course.midi).unwrap();
            assert_eq!(frequencies.len(), rates.len());
            for (&frequency, &rate) in frequencies.iter().zip(rates) {
                // Independent dimensional reference for the existing model:
                // resistance per length / solid mass per length, both SI.
                let expected = TAU * atmosphere.dynamic_viscosity
                    * (1.0 + OUTER * (2.0 * atmosphere.density * TAU * frequency
                        / atmosphere.dynamic_viscosity).sqrt())
                    / course.linear_density_kg_m;
                near(rate, expected);
            }
        }
        assert!(fs_live_piano::describe(&piano).contains(
            "\"stringGasLoading\":\"sampled-cylinder-retained-modes\""));
        prepared.push(piano);
    }
    let a = prepared[0].string_gas_drag_rates_per_s(69).unwrap();
    let b = prepared[1].string_gas_drag_rates_per_s(69).unwrap();
    let c = prepared[2].string_gas_drag_rates_per_s(69).unwrap();
    assert!(a.iter().zip(b).all(|(a, b)| b > a));
    assert_ne!(b, c);
    // This reduced resistance carries no fluid inertia; gas substitutions do
    // not retune the source's elastic frequencies or installed tension.
    assert_eq!(prepared[0].partials_hz(69), prepared[1].partials_hz(69));
    assert_eq!(prepared[0].partials_hz(69), prepared[2].partials_hz(69));
    assert_eq!(prepared[1].courses(), desc.courses.as_slice());
    let effective = LivePiano::new_with_gas(&small_description(), config(), gases[1]).unwrap();
    assert!(effective.string_gas_drag_rates_per_s(69).is_none());
}

#[test]
fn atmosphere_changes_mechanics_and_decay_with_lossless_solid_sources_at_both_clocks() {
    let mut desc = small_description();
    let source = load(&desc.courses,
        &pack([CORE[0], CORE[1], 0.0]), &pack([COATING[0], COATING[1], 0.0]));
    preparation::prepare_with_materials(&mut desc, None, None, Some(&source)).unwrap();
    let low = preparation::resolve_atmosphere("293.15 50000 0").unwrap();
    let high = preparation::resolve_atmosphere("293.15 100000 0").unwrap();
    assert_eq!(low.sound_speed, high.sound_speed);
    for silent_band in [0.0, 800.0] {
        for ratio in [1, 2] {
            let mut cfg = config();
            cfg.silent_bank_hz = silent_band;
            let mut a = PianoAudio::new_with_materials(&desc, cfg.clone(), ratio, low, source.clone()).unwrap();
            let mut b = PianoAudio::new_with_materials(&desc, cfg, ratio, high, source.clone()).unwrap();
            a.note_on(69, 0.25).unwrap();
            b.note_on(69, 0.25).unwrap();
            let (mut pressure, mut difference, mut motion_difference) = (0.0_f64, 0.0_f64, 0.0_f64);
            for _ in 0..4096 {
                // Remove the Rayleigh receiver's known density gain. Any
                // remaining difference comes from the solved source history.
                let x = a.step_pressure().unwrap() / low.density;
                let y = b.step_pressure().unwrap() / high.density;
                pressure = pressure.max(x.abs());
                difference = difference.max((x-y).abs());
                motion_difference = motion_difference.max((a.piano().bridge_velocity(69).unwrap()
                    - b.piano().bridge_velocity(69).unwrap()).abs());
            }
            assert!(pressure > 1e-9);
            assert!(difference > 1e-8 * pressure, "source response must change beyond receiver gain");
            assert!(motion_difference > 1e-14, "actual bridge motion must respond to gas drag");
            let la = a.close_ledger();
            let lb = b.close_ledger();
            // There is no material bending loss in this fixture. Gas loss
            // enters the original string ledger once, with no extra heat credit.
            assert!(la.string_loss_j > 0.0 && lb.string_loss_j > 0.0);
            assert_ne!(la.string_loss_j, lb.string_loss_j);
            assert!(la.residual_j.abs() < 1e-7 * la.input_j.max(1.0));
            assert!(lb.residual_j.abs() < 1e-7 * lb.input_j.max(1.0));
            let mut copied = b.clone();
            assert_eq!(b.piano().string_gas_drag_rates_per_s(69), copied.piano().string_gas_drag_rates_per_s(69));
            let mut whole = [0.0_f32; 257];
            let mut blocks = [0.0_f32; 257];
            b.render(&mut whole).unwrap();
            for block in blocks.chunks_mut(37) { copied.render(block).unwrap(); }
            assert_eq!(whole, blocks);
            assert_eq!(format!("{:?}", b.close_ledger()), format!("{:?}", copied.close_ledger()));
        }
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * actual.abs().max(expected.abs()),
        "{actual} != {expected}"
    );
}

#[test]
fn canonical_layers_drive_mass_axial_bending_and_original_modal_loss_integrals() {
    let desc = small_description();
    let core = pack(CORE);
    let coating = pack(COATING);
    let source = load(&desc.courses, &core, &coating);
    let changed = source.apply(&desc.courses).unwrap();
    // Independent annular integrals; do not use the section's cached geometry.
    let areas = [PI * INNER.powi(2), PI * (OUTER.powi(2) - INNER.powi(2))];
    let moments = [
        PI * INNER.powi(4) / 4.0,
        PI * (OUTER.powi(4) - INNER.powi(4)) / 4.0,
    ];
    let mu = CORE[0] * areas[0] + COATING[0] * areas[1];
    let ea = CORE[1] * areas[0] + COATING[1] * areas[1];
    let ei = CORE[1] * moments[0] + COATING[1] * moments[1];
    let cb = CORE[2] * moments[0] + COATING[2] * moments[1];
    for (old, new) in desc.courses.iter().zip(&changed) {
        near(new.linear_density_kg_m, mu);
        near(new.flexural_rigidity_nm2, ei);
        let bonded = new.bonded_section.unwrap();
        near(bonded.axial_rigidity_n(), ea);
        near(bonded.bending_loss().unwrap().viscous_rigidity_nm2s(), cb);
        assert_eq!(
            *old,
            Course {
                linear_density_kg_m: old.linear_density_kg_m,
                flexural_rigidity_nm2: old.flexural_rigidity_nm2,
                bonded_section: old.bonded_section,
                ..*new
            }
        );
        let section = &source.sections()[&old.midi];
        assert_eq!(section.omega_band_rad_s(), (0.0, TAU * MAX_HZ));
        for (layer, original) in section.elastic().layers().iter().zip([&core, &coating]) {
            assert_eq!(
                layer.material.card_identity(),
                original.card().content_hash()
            );
            for property in layer.material.properties() {
                let receipt = &property.answer().receipt;
                assert!(original.card().claims().claim(receipt.selected).is_some());
                assert_eq!(
                    receipt.axis_quantities["bath_kelvin"],
                    typed(QuantityKind::AbsoluteTemperature)
                );
            }
        }
    }
    assert_eq!(source.packs()[0].to_bytes(), core.to_bytes());
    assert_eq!(source.packs()[1].to_bytes(), coating.to_bytes());
    assert_eq!(source.packs()[0].claims_pack().normalizations().len(), 4);
    let bank =
        Bank::new_with_string_damping(&changed, &desc.board, 48_000, 6000.0, 12, true, None, false)
            .unwrap();
    let (losses, endpoint_drag) = bank_string_losses(&bank);
    assert!(endpoint_drag.iter().all(|drag| *drag == 0.0));
    for string in &bank.strings {
        let course = &changed[string.course];
        for (n, index) in string.modes.clone().enumerate() {
            near(
                losses[index].strain_rate,
                cb * ((n + 1) as f64 * PI / course.length_m).powi(4) / mu,
            );
            assert_eq!(losses[index].mass_rate, 0.0);
        }
    }
    // Explicitly lossless source data must not regain the demonstration's
    // constant velocity drag or its bending spectrum in the original bank.
    let zero = load(
        &desc.courses,
        &pack([CORE[0], CORE[1], 0.0]),
        &pack([COATING[0], COATING[1], 0.0]),
    );
    let lossless = zero.apply(&desc.courses).unwrap();
    let bank = Bank::new_with_string_damping(
        &lossless,
        &desc.board,
        48_000,
        6000.0,
        12,
        true,
        None,
        false,
    )
    .unwrap();
    let (losses, endpoint_drag) = bank_string_losses(&bank);
    assert!(
        losses
            .iter()
            .all(|loss| loss.mass_rate == 0.0 && loss.strain_rate == 0.0)
    );
    assert!(endpoint_drag.iter().all(|drag| *drag == 0.0));
    // Changing one input property cannot quietly change another section integral.
    for (index, value) in [(0, 16000.0), (1, 4e11), (2, 1.6e9)] {
        let mut values = CORE;
        values[index] = value;
        let alternative = load(&desc.courses, &pack(values), &coating);
        let section = &alternative.sections()[&69];
        let expected = [
            mu + (values[0] - CORE[0]) * areas[0],
            ea + (values[1] - CORE[1]) * areas[0],
            ei + (values[1] - CORE[1]) * moments[0],
            cb + (values[2] - CORE[2]) * moments[0],
        ];
        for (a, b) in [
            section.elastic().line_density_kg_m(),
            section.elastic().axial_stiffness_n(),
            section.elastic().bending_stiffness_n_m2(),
            section.viscous_stiffness_n_m2_s(),
        ]
        .into_iter()
        .zip(expected)
        {
            near(a, b);
        }
    }
}

#[test]
fn material_preparation_preserves_installed_tension_then_tunes_the_changed_pde() {
    let mut desc = small_description();
    let original = desc.courses.clone();
    let board = desc.board_bytes();
    let source = load(&original, &pack(CORE), &pack(COATING));
    preparation::prepare_with_materials(&mut desc, None, None, Some(&source)).unwrap();
    for (a, b) in original.iter().zip(&desc.courses) {
        assert_eq!(a.tension_n, b.tension_n);
    }
    assert_eq!(desc.board_bytes(), board);
    let mut expected = source.apply(&original).unwrap();
    stretch_tuning::partials::tune(&mut expected, 442.0, tuning().profile()).unwrap();
    let csv = geometry::write_scale(&original);
    preparation::prepare_with_materials(&mut desc, Some(&csv), Some(tuning()), Some(&source))
        .unwrap();
    assert_eq!(desc.courses, expected);
    assert_eq!(desc.board_bytes(), board);
    let a4 = desc.courses.iter().find(|c| c.midi == 69).unwrap();
    // Independent inversion of the newly integrated stiff-string equation.
    let k = PI / a4.length_m;
    near(
        a4.tension_n,
        (TAU * 442.0 / k).powi(2) * a4.linear_density_kg_m - a4.flexural_rigidity_nm2 * k * k,
    );
    near(a4.partial_hz(1, a4.tension_n), 442.0);
    let native = LivePiano::new_with_materials(&desc, config(), gas(), source.clone()).unwrap();
    for c in &desc.courses {
        near(native.partials_hz(c.midi)[0], c.partial_hz(1, c.tension_n));
    }
    assert!(std::ptr::eq(
        native.string_materials().unwrap(),
        source.as_ref()
    ));
    assert!(std::ptr::eq(
        native.clone().string_materials().unwrap(),
        source.as_ref()
    ));
}

#[test]
fn each_material_only_change_reaches_live_pressure_and_cloned_audio_at_all_clocks() {
    let coating = pack(COATING);
    for ratio in [1, 2, 4] {
        let mut pcm = Vec::new();
        let mut fundamentals = Vec::new();
        for values in [
            CORE,
            [16000.0, CORE[1], CORE[2]],
            [CORE[0], 4e11, CORE[2]],
            [CORE[0], CORE[1], 1.6e9],
        ] {
            let mut desc = small_description();
            let source = load(&desc.courses, &pack(values), &coating);
            preparation::prepare_with_materials(&mut desc, None, None, Some(&source)).unwrap();
            let mut audio =
                PianoAudio::new_with_materials(&desc, config(), ratio, gas(), source.clone())
                    .unwrap();
            fundamentals.push(audio.piano().partials_hz(69)[0]);
            assert_eq!(audio.mechanical_rate(), 48_000 * ratio as u32);
            assert_eq!(audio.piano().courses(), desc.courses.as_slice());
            audio.note_on(69, 1.0).unwrap();
            audio.render(&mut [0.0_f32; 193]).unwrap();
            let mut split = audio.clone();
            assert!(std::ptr::eq(
                split.piano().string_materials().unwrap(),
                source.as_ref()
            ));
            let mut a = vec![0.0_f32; 768];
            let mut b = a.clone();
            audio.render(&mut a).unwrap();
            for chunk in b.chunks_mut(37) {
                split.render(chunk).unwrap();
            }
            assert_eq!(a, b);
            assert!(a.iter().any(|x| *x != 0.0));
            assert!(a.iter().all(|x| x.is_finite()));
            assert_eq!(audio.piano().counters.contact_failures, 0);
            let ledger = audio.close_ledger();
            assert!(ledger.input_j > 0.0 && ledger.string_loss_j > 0.0);
            assert_eq!(format!("{ledger:?}"), format!("{:?}", split.close_ledger()));
            let metadata = fs_live_piano::describe_audio(&audio);
            assert!(metadata.contains("\"stringDamping\":\"source-kelvin-voigt\""));
            assert!(metadata.contains(&format!("\"stringMaterials\":{}", source.json())));
            pcm.push(a);
        }
        assert!(fundamentals[1] < fundamentals[0]); // only rho increased
        assert!(fundamentals[2] > fundamentals[0]); // only E increased
        assert_eq!(fundamentals[3], fundamentals[0]); // eta changes loss, not elastic frequency
        for alternative in &pcm[1..] {
            assert_ne!(&pcm[0], alternative);
        }
    }
}

#[test]
fn unsupported_source_geometry_framing_and_actual_mode_band_refuse() {
    let desc = small_description();
    let core = pack(CORE).to_bytes();
    let coating = pack(COATING).to_bytes();
    let geometry = section_geometry(&desc.courses);
    // A missing loss claim is deferred past the tensile check, then refuses
    // as missing data before anything can play.
    let missing = pack_with(CORE, |claim| {
        if claim.key.name() == ETA {
            claim.key = PropertyKey::new("missing_bending_eta", DynViscosity::DIMS);
        }
    });
    let deferred = StringMaterials::from_packs(&[&missing.to_bytes(), &coating], &geometry, &point(), ETA)
        .unwrap();
    let error = deferred.apply(&desc.courses).unwrap_err();
    assert!(error.contains("missing data") && error.contains(ETA), "{error}");
    let bad_sources = [
        pack_with(CORE, |claim| {
            if claim.key.name() == ETA {
                claim.key = PropertyKey::new(ETA, Pressure::DIMS);
                claim.value = PropertyValue::Scalar {
                    value: CORE[2],
                    dims: Pressure::DIMS,
                };
            }
        }),
        pack_with(CORE, |claim| {
            if claim.key.name() == ETA {
                claim.validity = ValidityDomain::unconstrained();
            }
        }),
        pack_with(CORE, |claim| {
            if claim.key.name() == ETA {
                claim.value = PropertyValue::Curve {
                    abscissa: "oscillation_hz".into(),
                    abscissa_dims: Frequency::DIMS,
                    knots: vec![(0.0, CORE[2]), (MAX_HZ, 2.0 * CORE[2])],
                    dims: DynViscosity::DIMS,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
        }),
        pack([CORE[0], CORE[1], -1.0]),
    ];
    for bad in bad_sources {
        assert!(
            StringMaterials::from_packs(&[&bad.to_bytes(), &coating], &geometry, &point(), ETA)
                .is_err()
        );
    }
    for rows in [
        "69,core,0,,,0",
        "69,core,NaN,,,0",
        "69,core,0.0005,,,0\n69,core,0.0004,,,1",
        "69,core,0.0005,,,2",
    ] {
        let bad = construction("installed-tension", LIMIT, &format!("{rows}\n"));
        assert!(StringMaterials::from_packs(&[&core, &coating], &bad, &point(), ETA).is_err());
    }
    let outside = operating_point("bath_kelvin", 500.0, "oscillation_hz", 1000.0).unwrap();
    assert!(StringMaterials::from_packs(&[&core, &coating], &geometry, &outside, ETA).is_err());
    let bytes = [core.clone(), coating.clone()].concat();
    let lengths = [core.len() as u32, coating.len() as u32];
    let framed =
        StringMaterials::from_packed_bytes(&bytes, &lengths, &geometry, &point(), ETA).unwrap();
    assert_eq!(framed.packs()[0].to_bytes(), core);
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(
        StringMaterials::from_packed_bytes(&trailing, &lengths, &geometry, &point(), ETA).is_err()
    );
    assert!(
        StringMaterials::from_packed_bytes(
            &bytes[..bytes.len() - 1],
            &lengths,
            &geometry,
            &point(),
            ETA
        )
        .is_err()
    );
    let mut bad_pack = core.clone();
    bad_pack.push(0);
    assert!(StringMaterials::from_packs(&[&bad_pack, &coating], &geometry, &point(), ETA).is_err());
    // A valid sampled state is insufficient: actual constructed modes must fit
    // every source coefficient's frequency band, including cyclic-Hz conversion.
    let narrow = pack_with(CORE, |claim| {
        claim.validity = validity(999.0, 1001.0);
    });
    let source = load(&desc.courses, &narrow, &pack(COATING));
    let mut prepared = small_description();
    preparation::prepare_with_materials(&mut prepared, None, None, Some(&source)).unwrap();
    assert!(LivePiano::new_with_materials(&prepared, config(), gas(), source).is_err());
}

#[test]
fn partial_coverage_preset_and_late_tuning_failure_do_not_publish_partial_preparation() {
    let mut desc = small_description();
    let before = desc.courses.clone();
    let board = desc.board_bytes();
    let pointer = desc.courses.as_ptr();
    let partial = load(&before[..2], &pack(CORE), &pack(COATING));
    assert!(preparation::prepare_with_materials(&mut desc, None, None, Some(&partial)).is_err());
    assert_eq!(desc.courses, before);
    assert_eq!(desc.courses.as_ptr(), pointer);
    let complete = load(&before, &pack(CORE), &pack(COATING));
    let band = Tuning::new(442.0, [1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], Some(1.0)).unwrap();
    let mut imported = before.clone();
    imported[0].hammer_mass_kg *= 1.25;
    let csv = geometry::write_scale(&imported);
    assert!(
        preparation::prepare_with_materials(&mut desc, Some(&csv), Some(band), Some(&complete))
            .is_err()
    );
    assert_eq!(desc.courses, before);
    assert_eq!(desc.courses.as_ptr(), pointer);
    assert_eq!(desc.board_bytes(), board);
    assert!(LivePiano::new_with_materials(&desc, config(), gas(), complete.clone()).is_err());
    preparation::prepare_with_materials(&mut desc, None, None, Some(&complete)).unwrap();
    let mut preset = config();
    preset.rt0425_string_damping = true;
    assert!(LivePiano::new_with_materials(&desc, preset, gas(), complete.clone()).is_err());
    let other = load(&before, &pack([16000.0, CORE[1], CORE[2]]), &pack(COATING));
    assert!(PianoAudio::new_with_materials(&desc, config(), 2, gas(), other).is_err());
    assert!(PianoAudio::new_with_materials(&desc, config(), 2, gas(), complete).is_ok());
    assert_eq!(desc.board_bytes(), board);
}

// ---- U2: tension policy, the mandatory tensile check and wound strings ----

const CORE_R: f64 = 0.0005;
fn limited(values: [f64; 3], limit: f64) -> NormalizedMaterialCardPack {
    pack_with(values, |claim| {
        if claim.key.name() == LIMIT {
            claim.value = PropertyValue::Scalar { value: limit, dims: Pressure::DIMS };
        }
    })
}
// A dense, compliant, weak manufactured metal; not a claim about any alloy.
const WEAK: [f64; 3] = [19_300.0, 7.9e10, 1e8];
const WEAK_LIMIT_PA: f64 = 1.2e8;
fn single_cores(courses: &[Course]) -> String {
    courses.iter().map(|c| format!("{},core,{CORE_R},,,0\n", c.midi)).collect()
}
fn weak(tension: &str, strength: &str) -> Result<StringMaterials, String> {
    StringMaterials::from_packs(&[&limited(WEAK, WEAK_LIMIT_PA).to_bytes()],
        &construction(tension, strength, &single_cores(&small_description().courses)), &point(), ETA)
}

// Planted negatives: the retired v1 format and a core pack without the
// declared limit refuse as missing data; a wrap carries no tension and needs none.
#[test]
fn g0_v1_constructions_and_unrated_core_packs_refuse_as_missing_data() {
    let courses = small_description().courses;
    let v1 = format!("frankensim-piano-materials-v1\nkey,outer_radius_m,material_index\n{}",
        courses.iter().map(|c| format!("{},{CORE_R},0\n", c.midi)).collect::<String>());
    let core = pack(CORE).to_bytes();
    let error = StringMaterials::from_packs(&[&core], &v1, &point(), ETA).unwrap_err();
    assert!(error.contains("missing data"), "{error}");
    let unrated = pack_with(CORE, |claim| {
        if claim.key.name() == LIMIT {
            claim.key = PropertyKey::new("some_other_strength", Pressure::DIMS);
        }
    }).to_bytes();
    let geometry = construction("target-pitch", LIMIT, &single_cores(&courses));
    let error = StringMaterials::from_packs(&[&unrated], &geometry, &point(), ETA).unwrap_err();
    assert!(error.contains("missing data") && error.contains(LIMIT), "{error}");
    // The same unrated pack is admitted as a wrap, whose wire carries no tension.
    let wound = construction("target-pitch", LIMIT, &courses.iter()
        .map(|c| format!("{},core,0.0004,,,0\n{},wrap,,0.0003,close,1\n", c.midi, c.midi)).collect::<String>());
    assert!(StringMaterials::from_packs(&[&core, &unrated], &wound, &point(), ETA).is_ok());
    for bad in ["tension,target-pitch\nstrength,", "tension,yield-fraction,0\nstrength,x",
        "tension,yield-fraction,1.5\nstrength,x", "tension,somehow\nstrength,x"] {
        let text = format!("frankensim-piano-materials-v2\n{bad}\n{COLUMNS}\n{}", single_cores(&courses));
        assert!(StringMaterials::from_packs(&[&core], &text, &point(), ETA).is_err(), "{bad}");
    }
}

// G0: target-pitch re-solves tension to the installed fundamental; the
// installed-tension policy keeps the tension and lets the pitch move.
#[test]
fn g0_target_pitch_keeps_installed_fundamentals_and_installed_tension_keeps_tension() {
    let desc = small_description();
    let core = pack(CORE).to_bytes();
    let target = StringMaterials::from_packs(&[&core],
        &construction("target-pitch", LIMIT, &single_cores(&desc.courses)), &point(), ETA).unwrap();
    assert_eq!(target.tension_policy(), TensionPolicy::TargetPitch);
    let kept = StringMaterials::from_packs(&[&core],
        &construction("installed-tension", LIMIT, &single_cores(&desc.courses)), &point(), ETA).unwrap();
    let (a, b) = (target.apply(&desc.courses).unwrap(), kept.apply(&desc.courses).unwrap());
    for ((old, tuned), held) in desc.courses.iter().zip(&a).zip(&b) {
        near(tuned.partial_hz(1, tuned.tension_n), old.partial_hz(1, old.tension_n));
        assert_eq!(held.tension_n, old.tension_n);
        assert_eq!(Course { tension_n: held.tension_n, ..*tuned }, *held);
        assert_ne!(held.partial_hz(1, held.tension_n), old.partial_hz(1, old.tension_n));
    }
}

// G1: at a fixed pitch the core stress is 4 L^2 f^2 rho less the bending
// share, whatever the radius. A dense, weak core cannot reach the scale:
// every key refuses together with its stress, its limit and f_max.
#[test]
fn g1_weak_dense_core_refuses_at_pitch_naming_every_key_stress_limit_and_f_max() {
    let desc = small_description();
    let error = weak("target-pitch", LIMIT).unwrap().apply(&desc.courses).unwrap_err();
    assert!(error.contains("\"refusal\":\"overstressed-strings\""), "{error}");
    for c in &desc.courses { assert!(error.contains(&format!("\"key\":{}", c.midi)), "{error}"); }
    let labelled = weak("target-pitch", &format!("{LIMIT},non-physical-override")).unwrap();
    let played = labelled.apply(&desc.courses).unwrap();
    let report = labelled.strength_report(&played).unwrap();
    assert!(report.non_physical_override);
    assert_eq!(report.overstressed().count(), desc.courses.len());
    let (rho, young, area) = (WEAK[0], WEAK[1], PI * CORE_R * CORE_R);
    for ((old, c), key) in desc.courses.iter().zip(&played).zip(&report.keys) {
        let f = old.partial_hz(1, old.tension_n);
        let k = PI / c.length_m;
        assert_eq!(key.key, c.midi);
        near(key.required_stress_pa, rho * (2.0 * c.length_m * f).powi(2) - young * CORE_R * CORE_R * k * k / 4.0);
        near(key.limit_pa, WEAK_LIMIT_PA);
        near(key.capacity_n, WEAK_LIMIT_PA * area);
        assert!(key.utilization > 1.0);
        // Independent stiff-string pitch at the capacity tension, and the
        // flexible-string limit sqrt(limit / rho) / (2 L) it approaches.
        let ei = young * PI * CORE_R.powi(4) / 4.0;
        near(key.f_max_hz, ((WEAK_LIMIT_PA * area * k * k + ei * k.powi(4)) / (rho * area)).sqrt() / TAU);
        let flexible = (WEAK_LIMIT_PA / rho).sqrt() / (2.0 * c.length_m);
        assert!(key.f_max_hz > flexible && key.f_max_hz / flexible - 1.0 < ei * k * k / (WEAK_LIMIT_PA * area));
        assert!(key.f_max_hz < f);
    }
    assert!(labelled.json().contains("\"nonPhysicalOverride\":true"));
}

// G1: yield-fraction sets tension from the sourced capacity, so the weak
// core plays physically, lower than the scale, and the engine renders it.
#[test]
fn g1_yield_fraction_tension_follows_capacity_and_the_weak_string_plays() {
    let mut desc = small_description();
    let source = Arc::new(weak("yield-fraction,0.5", LIMIT).unwrap());
    let before = desc.courses.clone();
    preparation::prepare_with_materials(&mut desc, None, None, Some(&source)).unwrap();
    let report = source.strength_report(&desc.courses).unwrap();
    for ((old, c), key) in before.iter().zip(&desc.courses).zip(&report.keys) {
        near(c.tension_n, 0.5 * key.capacity_n);
        near(key.utilization, 0.5);
        assert!(c.partial_hz(1, c.tension_n) < old.partial_hz(1, old.tension_n));
        assert!(c.partial_hz(1, c.tension_n) < key.f_max_hz);
    }
    assert_eq!(report.overstressed().count(), 0);
    let mut audio = PianoAudio::new_with_materials(&desc, config(), 1, gas(), source.clone()).unwrap();
    audio.note_on(69, 1.0).unwrap();
    let mut pcm = vec![0.0_f32; 2048];
    audio.render(&mut pcm).unwrap();
    assert!(pcm.iter().all(|x| x.is_finite()) && pcm.iter().any(|x| *x != 0.0));
    assert_eq!(audio.piano().counters.contact_failures, 0);
    assert!(fs_live_piano::describe_audio(&audio).contains("\"yieldFraction\":0.5"));
}

// G0: retuning raises tension, so a string tuned up past its limit refuses
// and the description is not published.
#[test]
fn g0_retuning_past_the_limit_refuses_without_publishing() {
    let mut desc = small_description();
    let rows = single_cores(&desc.courses);
    let probe = StringMaterials::from_packs(&[&pack(CORE).to_bytes()],
        &construction("target-pitch", LIMIT, &rows), &point(), ETA).unwrap();
    let report = probe.strength_report(&probe.apply(&desc.courses).unwrap()).unwrap();
    let stress = report.keys.iter().map(|k| k.required_stress_pa).fold(0.0, f64::max);
    let tight = StringMaterials::from_packs(&[&limited(CORE, 1.001 * stress).to_bytes()],
        &construction("target-pitch", LIMIT, &rows), &point(), ETA).unwrap();
    let before = desc.courses.clone();
    preparation::prepare_with_materials(&mut small_description(), None, None, Some(&tight)).unwrap();
    let error = preparation::prepare_with_materials(&mut desc, None, Some(tuning()), Some(&tight)).unwrap_err();
    assert!(error.contains("overstressed-strings"), "{error}");
    assert_eq!(desc.courses, before);
}

// G1: a wound course: the core carries tension, EA, EI and loss; the wrap adds
// its helix mass only and moves the gas drag to the finished radius.
#[test]
fn g1_wound_course_core_carries_tension_and_the_wrap_adds_only_mass() {
    let mut desc = small_description();
    let (core_r, wire) = (0.0004, 0.0004);
    let rows = desc.courses.iter()
        .map(|c| format!("{},core,{core_r},,,0\n{},wrap,,{wire},close,1\n", c.midi, c.midi)).collect::<String>();
    let wrap_rho = 8960.0;
    let source = Arc::new(StringMaterials::from_packs(
        &[&pack(CORE).to_bytes(), &pack([wrap_rho, 1.2e11, 0.0]).to_bytes()],
        &construction("target-pitch", LIMIT, &rows), &point(), ETA).unwrap());
    let before = desc.courses.clone();
    preparation::prepare_with_materials(&mut desc, None, None, Some(&source)).unwrap();
    let pitch = fs_material::section::wound::close_wound_pitch_m(2.0 * core_r, wire).unwrap();
    let wrap_mass = wrap_rho * PI * wire * wire / 4.0 * (PI * (2.0 * core_r + wire)).hypot(pitch) / pitch;
    let (core_mass, ea, ei) = (CORE[0] * PI * core_r * core_r, CORE[1] * PI * core_r * core_r,
        CORE[1] * PI * core_r.powi(4) / 4.0);
    let report = source.strength_report(&desc.courses).unwrap();
    for ((old, c), key) in before.iter().zip(&desc.courses).zip(&report.keys) {
        near(c.linear_density_kg_m, core_mass + wrap_mass);
        near(c.flexural_rigidity_nm2, ei);
        let section = c.bonded_section.unwrap();
        near(section.axial_rigidity_n(), ea);
        near(section.wrap_line_density_kg_m(), wrap_mass);
        near(section.bending_loss().unwrap().viscous_rigidity_nm2s(), CORE[2] * PI * core_r.powi(4) / 4.0);
        near(c.partial_hz(1, c.tension_n), old.partial_hz(1, old.tension_n));
        // The whole tension sits in the core.
        near(key.required_stress_pa, c.tension_n / (PI * core_r * core_r));
        near(key.capacity_n, LIMIT_PA * PI * core_r * core_r);
    }
    // Only the wrap density changed: mass and pitch follow, core coefficients do not.
    let heavy = StringMaterials::from_packs(
        &[&pack(CORE).to_bytes(), &pack([2.0 * wrap_rho, 1.2e11, 0.0]).to_bytes()],
        &construction("installed-tension", LIMIT, &rows), &point(), ETA).unwrap();
    for (c, h) in desc.courses.iter().zip(heavy.apply(&desc.courses).unwrap()) {
        near(h.linear_density_kg_m, core_mass + 2.0 * wrap_mass);
        assert_eq!(h.flexural_rigidity_nm2, c.flexural_rigidity_nm2);
        assert_eq!(h.tension_n, c.tension_n);
        assert!(h.partial_hz(1, h.tension_n) < c.partial_hz(1, c.tension_n));
    }
    let atmosphere = gas();
    let piano = LivePiano::new_with_materials(&desc, config(), atmosphere, source.clone()).unwrap();
    let outer = core_r + wire;
    assert_eq!(source.outer_radius_m(69), Some(outer));
    for course in piano.courses() {
        let frequencies = piano.partials_hz(course.midi);
        let rates = piano.string_gas_drag_rates_per_s(course.midi).unwrap();
        for (&frequency, &rate) in frequencies.iter().zip(rates) {
            let expected = TAU * atmosphere.dynamic_viscosity
                * (1.0 + outer * (2.0 * atmosphere.density * TAU * frequency
                    / atmosphere.dynamic_viscosity).sqrt())
                / course.linear_density_kg_m;
            near(rate, expected);
        }
    }
    assert!(source.json().contains(&format!("\"outerRadiusM\":{outer}")));
}

// ---- U2: loss-angle cores ----

const TAN_DELTA: f64 = 2e-3;
fn tangent_pack(values: [f64; 3], tangent: f64) -> NormalizedMaterialCardPack {
    pack_with(values, |claim| {
        if claim.key.name() == ETA {
            claim.key = PropertyKey::new(ETA, Dims::NONE);
            claim.value = PropertyValue::Scalar { value: tangent, dims: Dims::NONE };
        }
    })
}
fn loss_angle_construction(courses: &[Course]) -> String {
    format!("frankensim-piano-materials-v2\ntension,target-pitch\nstrength,{LIMIT}\nloss,loss-angle\n{COLUMNS}\n{}",
        single_cores(courses))
}

// G1: a loss tangent gives each mode the equal-energy-per-cycle rate
// H k^4 / (mu omega), H = E I tan(delta): it changes decay, not frequency,
// and it falls one power of omega below a Kelvin-Voigt k^4 law.
#[test]
fn g1_loss_angle_cores_give_each_mode_its_hysteretic_rate_without_moving_frequencies() {
    let desc = small_description();
    let hysteretic = StringMaterials::from_packs(&[&tangent_pack(CORE, TAN_DELTA).to_bytes()],
        &loss_angle_construction(&desc.courses), &point(), ETA).unwrap();
    assert_eq!(hysteretic.loss_angle_sections().len(), 3);
    assert!(hysteretic.sections().is_empty());
    let viscous = StringMaterials::from_packs(&[&pack(CORE).to_bytes()],
        &construction("target-pitch", LIMIT, &single_cores(&desc.courses)), &point(), ETA).unwrap();
    let (a, b) = (hysteretic.apply(&desc.courses).unwrap(), viscous.apply(&desc.courses).unwrap());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(Course { bonded_section: y.bonded_section, ..*x }, *y);
        let loss = x.bonded_section.unwrap().bending_loss().unwrap();
        assert_eq!(loss.viscous_rigidity_nm2s(), 0.0);
        near(loss.hysteretic_rigidity_nm2(), CORE[1] * PI * CORE_R.powi(4) / 4.0 * TAN_DELTA);
    }
    let bank = |courses: &[Course]| Bank::new_with_string_damping(courses, &desc.board, 48_000, 6000.0, 12,
        true, None, false).unwrap();
    let (loss_bank, kv_bank) = (bank(&a), bank(&b));
    let (losses, _) = bank_string_losses(&loss_bank);
    let h = CORE[1] * PI * CORE_R.powi(4) / 4.0 * TAN_DELTA;
    for string in &loss_bank.strings {
        let course = &a[string.course];
        for (n, index) in string.modes.clone().enumerate() {
            let k = (n + 1) as f64 * PI / course.length_m;
            // Independent stiff-string frequency, not the bank's stored omega.
            let omega = TAU * course.partial_hz(n + 1, course.tension_n);
            near(losses[index].strain_rate, h * k.powi(4) / (course.linear_density_kg_m * omega));
            assert_eq!(losses[index].mass_rate, 0.0);
            assert_eq!(loss_bank.modes[index].omega, kv_bank.modes[index].omega);
        }
        let first = losses[string.modes.start].strain_rate;
        let last = string.modes.end - 1;
        let order = (last - string.modes.start + 1) as f64;
        near(losses[last].strain_rate / first, order.powi(4) * loss_bank.modes[string.modes.start].omega
            / loss_bank.modes[last].omega);
    }
    assert!(hysteretic.json().contains("\"bendingLoss\":\"loss-angle\""));
    let mut prepared = small_description();
    let source = Arc::new(hysteretic);
    preparation::prepare_with_materials(&mut prepared, None, None, Some(&source)).unwrap();
    let mut audio = PianoAudio::new_with_materials(&prepared, config(), 1, gas(), source.clone()).unwrap();
    assert!(fs_live_piano::describe_audio(&audio).contains("\"stringDamping\":\"source-loss-angle\""));
    audio.note_on(69, 1.0).unwrap();
    let mut pcm = vec![0.0_f32; 2048];
    audio.render(&mut pcm).unwrap();
    assert!(pcm.iter().all(|x| x.is_finite()) && pcm.iter().any(|x| *x != 0.0));
    assert!(audio.close_ledger().string_loss_j > 0.0);
}

// Planted negatives: a viscosity in Pa s is not a loss tangent, an unknown
// law refuses, and a band that misses a retained mode refuses at build.
#[test]
fn g0_loss_angle_law_refuses_wrong_quantity_unknown_law_and_uncovered_modes() {
    let desc = small_description();
    assert!(StringMaterials::from_packs(&[&pack(CORE).to_bytes()], &loss_angle_construction(&desc.courses),
        &point(), ETA).is_err());
    let unknown = loss_angle_construction(&desc.courses).replace("loss,loss-angle", "loss,quality-factor");
    assert!(StringMaterials::from_packs(&[&tangent_pack(CORE, TAN_DELTA).to_bytes()], &unknown, &point(), ETA)
        .is_err());
    let narrow = pack_with(CORE, |claim| {
        if claim.key.name() == ETA {
            claim.key = PropertyKey::new(ETA, Dims::NONE);
            claim.value = PropertyValue::Scalar { value: TAN_DELTA, dims: Dims::NONE };
            claim.validity = validity(8.0, 3000.0);
        }
    });
    let source = Arc::new(StringMaterials::from_packs(&[&narrow.to_bytes()],
        &loss_angle_construction(&desc.courses), &point(), ETA).unwrap());
    let mut prepared = small_description();
    preparation::prepare_with_materials(&mut prepared, None, None, Some(&source)).unwrap();
    assert!(LivePiano::new_with_materials(&prepared, config(), gas(), source).is_err());
}

// ---- U2: a string that would break is reported before missing loss data ----

// G0: a core pack without any bending-loss claim still resolves its mass,
// stiffness and strength. Overstressed strings refuse first, as strings
// that would break; a physically tensioned set then refuses as missing data.
#[test]
fn g0_overstress_refuses_before_missing_loss_data() {
    let desc = small_description();
    let unrated = |values: [f64; 3], limit: f64| pack_with(values, |claim| {
        if claim.key.name() == ETA {
            claim.key = PropertyKey::new("some_other_viscosity", DynViscosity::DIMS);
        }
        if claim.key.name() == LIMIT {
            claim.value = PropertyValue::Scalar { value: limit, dims: Pressure::DIMS };
        }
    });
    let weak = unrated(WEAK, WEAK_LIMIT_PA).to_bytes();
    let at_pitch = StringMaterials::from_packs(&[&weak],
        &construction("target-pitch", LIMIT, &single_cores(&desc.courses)), &point(), ETA).unwrap();
    assert!(at_pitch.missing_loss().unwrap().contains("missing data"));
    let error = at_pitch.apply(&desc.courses).unwrap_err();
    assert!(error.contains("overstressed-strings") && !error.contains("missing data"), "{error}");
    let report = at_pitch.strength_report(&desc.courses.iter().map(|c| Course {
        linear_density_kg_m: WEAK[0] * PI * CORE_R * CORE_R,
        flexural_rigidity_nm2: WEAK[1] * PI * CORE_R.powi(4) / 4.0, bonded_section: None, ..*c
    }).collect::<Vec<_>>()).unwrap();
    assert_eq!(report.keys.len(), 3);
    let relaxed = StringMaterials::from_packs(&[&weak],
        &construction("yield-fraction,0.5", LIMIT, &single_cores(&desc.courses)), &point(), ETA).unwrap();
    let error = relaxed.apply(&desc.courses).unwrap_err();
    assert!(error.contains("missing data") && error.contains(ETA) && !error.contains("overstressed"), "{error}");
    let json = relaxed.json();
    assert!(json.contains("\"missingLossData\":true") && json.contains("\"viscousRigidityNm2S\":null"), "{json}");
    // A rated pack in the same position plays: only the loss claim was missing.
    assert!(StringMaterials::from_packs(&[&limited(WEAK, WEAK_LIMIT_PA).to_bytes()],
        &construction("yield-fraction,0.5", LIMIT, &single_cores(&desc.courses)), &point(), ETA)
        .unwrap().apply(&desc.courses).is_ok());
}
