//! G0/G1/G3: source capacities compose into one conserved sensible heat account.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedPack, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw, PreparedHeatStorage,
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use std::{path::Path, process::Command};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CVOL: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const A: [(f64, f64); 3] = [(280.0, 100.0), (300.0, 140.0), (340.0, 100.0)];
const B: [(f64, f64); 3] = [(290.0, 1000.0), (310.0, 1600.0), (330.0, 1200.0)];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("specimen", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
}
fn key(dims: Dims) -> PropertyKey {
    PropertyKey::new(
        if dims == CP {
            "specific-heat-capacity"
        } else {
            "volumetric-heat-capacity"
        },
        dims,
    )
}
fn claim(dims: Dims, knots: &[(f64, f64)]) -> PropertyClaim {
    PropertyClaim {
        key: key(dims),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots: knots.to_vec(),
            dims,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("specimen", QuantitySpec::dimensional(Dims::NONE), 1.0, 1.0)
            .with_quantity(
                "temperature",
                semantic(QuantityKind::AbsoluteTemperature),
                knots[0].0,
                knots.last().unwrap().0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "authored independent caloric fixture; no measured material claim".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(name: &str, claims: impl IntoIterator<Item = PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "explicit fixed source context".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn component<'a>(
    id: u64,
    card: &'a MaterialCard,
    capacity: &'a PropertyKey,
    amount: HeatStorageAmount,
    lower: &'a QueryPoint,
    upper: &'a QueryPoint,
) -> HeatStorageComponentQuery<'a> {
    HeatStorageComponentQuery {
        id,
        card,
        capacity,
        amount,
        temperature_axis: "temperature",
        lower,
        upper,
        selection: SINGLE,
    }
}
fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.16e} != {expected:.16e}"
    );
}

// Independent dimensional sum: 2 kg of A plus 0.25 m3 of B. The source
// curves have different breaks and their common domain is 290..330 K.
fn capacity(t: f64) -> f64 {
    if t <= 300.0 {
        490.0 + 11.5 * (t - 290.0)
    } else if t <= 310.0 {
        605.0 + 5.5 * (t - 300.0)
    } else {
        660.0 - 7.0 * (t - 310.0)
    }
}
fn energy(t: f64) -> f64 {
    if t <= 300.0 {
        let x = t - 290.0;
        490.0 * x + 5.75 * x * x
    } else if t <= 310.0 {
        let x = t - 300.0;
        5475.0 + 605.0 * x + 2.75 * x * x
    } else {
        let x = t - 310.0;
        11800.0 + 660.0 * x - 3.5 * x * x
    }
}
fn entropy(t: f64) -> f64 {
    let segment = |lo: f64, hi: f64, c: f64, slope: f64| {
        (c - slope * lo) * (hi / lo).ln() + slope * (hi - lo)
    };
    let mut s = segment(290.0, t.min(300.0), 490.0, 11.5);
    if t > 300.0 {
        s += segment(300.0, t.min(310.0), 605.0, 5.5);
    }
    if t > 310.0 {
        s += segment(310.0, t, 660.0, -7.0);
    }
    s
}

#[test]
fn g1_mixed_mass_and_volume_capacities_keep_breaks_and_exact_sensible_storage() {
    let a = card("mass-capacity", [claim(CP, &A)]);
    let b = card("volume-capacity", [claim(CVOL, &B)]);
    let (ka, kb) = (key(CP), key(CVOL));
    let (a0, a1, b0, b1) = (point(280.0), point(340.0), point(290.0), point(330.0));
    let queries = [
        component(20, &a, &ka, HeatStorageAmount::MassKg(2.0), &a0, &a1),
        component(
            10,
            &b,
            &kb,
            HeatStorageAmount::ReferenceVolumeM3(0.25),
            &b0,
            &b1,
        ),
    ];
    let assembly = PreparedHeatStorage::try_new(&queries).unwrap();
    assert_eq!(assembly.temperature_span_k(), [290.0, 330.0]);
    assert_eq!(
        assembly.temperatures_k().collect::<Vec<_>>(),
        [290.0, 300.0, 310.0, 330.0]
    );
    assert_eq!(
        assembly
            .components()
            .iter()
            .map(|c| c.id())
            .collect::<Vec<_>>(),
        [10, 20]
    );
    for (part, source, low, high, selected_key) in [
        (&assembly.components()[0], &b, &b0, &b1, &kb),
        (&assembly.components()[1], &a, &a0, &a1, &ka),
    ] {
        assert_eq!(part.material_card_identity(), source.content_hash());
        assert_eq!(part.temperature_axis(), "temperature");
        assert_eq!(part.query_endpoints(), &[low.clone(), high.clone()]);
        assert_eq!(
            part.source_line(),
            &source
                .claims()
                .query_line_typed(selected_key, "temperature", low, high, SINGLE)
                .unwrap()
        );
        for (_, answer) in &part.source_line().knots {
            source.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
    }
    for initial in [290.0, 297.25, 300.0, 310.0, 329.0, 330.0] {
        close(
            assembly.heat_capacity_j_k(initial).unwrap(),
            capacity(initial),
            3e-12,
        );
        for end in [290.0, 295.5, 300.0, 305.25, 310.0, 325.75, 330.0] {
            let e = energy(end) - energy(initial);
            close(assembly.energy_change_j(initial, end).unwrap(), e, 3e-11);
            close(
                assembly
                    .temperature_after_energy_change(initial, e)
                    .unwrap(),
                end,
                2e-13,
            );
            close(
                assembly.sensible_entropy_change_j_k(initial, end).unwrap(),
                entropy(end) - entropy(initial),
                2e-12,
            );
            assert_eq!(
                assembly.energy_change_j(initial, end).unwrap(),
                -assembly.energy_change_j(end, initial).unwrap()
            );
        }
    }
    for invalid in [289.999, 330.001, f64::NAN, f64::INFINITY] {
        assert!(assembly.heat_capacity_j_k(invalid).is_err());
        assert!(assembly.energy_change_j(300.0, invalid).is_err());
    }
    assert!(
        assembly
            .temperature_after_energy_change(290.0, -1.0)
            .is_err()
    );
    assert!(
        assembly
            .temperature_after_energy_change(330.0, 1.0)
            .is_err()
    );
    assert!(
        assembly
            .temperature_after_energy_change(300.0, f64::NAN)
            .is_err()
    );
}

#[test]
fn g3_component_permutation_amounts_and_selected_sources_change_only_their_declared_response() {
    let original = claim(CP, &A);
    let mut alternate = original.clone();
    alternate.provenance.source =
        "independent alternative observation with the same numerical model".into();
    let pins = [
        ClaimSelection::Pinned(ClaimId(original.content_hash())),
        ClaimSelection::Pinned(ClaimId(alternate.content_hash())),
    ];
    let a = card("two-selectable-sources", [original, alternate]);
    let b = card("native-volume", [claim(CVOL, &B)]);
    let (ka, kb) = (key(CP), key(CVOL));
    let (lo, hi) = (point(290.0), point(330.0));
    let mut qa = component(2, &a, &ka, HeatStorageAmount::MassKg(2.0), &lo, &hi);
    let qb = component(
        1,
        &b,
        &kb,
        HeatStorageAmount::ReferenceVolumeM3(0.25),
        &lo,
        &hi,
    );
    assert!(PreparedHeatStorage::try_new(&[qa, qb]).is_err());
    qa.selection = pins[0];
    let base = PreparedHeatStorage::try_new(&[qa, qb]).unwrap();
    let reversed = PreparedHeatStorage::try_new(&[qb, qa]).unwrap();
    assert_eq!(base.identity(), reversed.identity());
    assert_eq!(
        base.energy_change_j(292.0, 328.0).unwrap(),
        reversed.energy_change_j(292.0, 328.0).unwrap()
    );
    qa.selection = pins[1];
    let selected = PreparedHeatStorage::try_new(&[qa, qb]).unwrap();
    assert_ne!(base.identity(), selected.identity());
    close(
        base.energy_change_j(292.0, 328.0).unwrap(),
        selected.energy_change_j(292.0, 328.0).unwrap(),
        0.0,
    );
    qa.selection = pins[0];
    qa.amount = HeatStorageAmount::MassKg(4.0);
    let mut doubled_b = qb;
    doubled_b.amount = HeatStorageAmount::ReferenceVolumeM3(0.5);
    let doubled = PreparedHeatStorage::try_new(&[qa, doubled_b]).unwrap();
    assert_ne!(base.identity(), doubled.identity());
    close(
        doubled.energy_change_j(292.0, 328.0).unwrap(),
        2.0 * base.energy_change_j(292.0, 328.0).unwrap(),
        0.0,
    );
    assert!(
        doubled
            .temperature_after_energy_change(300.0, 1000.0)
            .unwrap()
            < base.temperature_after_energy_change(300.0, 1000.0).unwrap()
    );

    // The volume component has positive mass in addition to A's known 2 kg.
    assert!(HeatStorageThermalLaw::try_new(&base, 2.0, 300.0).is_err());
    for reference_mass in [3.0, 20.0] {
        for datum in [290.0, 310.0] {
            let plain = HeatStorageThermalLaw::try_new(&base, reference_mass, datum).unwrap();
            let law = plain.clone().with_sensible_entropy();
            assert_ne!(plain.identity(), law.identity());
            assert_eq!(
                law.clone().with_sensible_entropy().identity(),
                law.identity()
            );
            for t in [290.0, 300.0, 330.0] {
                let h = law.specific_storage_j_kg(t).unwrap();
                close(h * reference_mass, energy(t) - energy(datum), 3e-11);
                close(law.temperature_k(h).unwrap(), t, 2e-13);
            }
            assert!(law.maximum_temperature_slope() >= reference_mass / 490.0);
            assert!(law.maximum_temperature_slope() <= (reference_mass / 490.0).next_up());
            let (h0, h1) = (
                law.specific_storage_j_kg(295.0).unwrap(),
                law.specific_storage_j_kg(325.0).unwrap(),
            );
            assert_eq!(plain.specific_entropy_change_j_kg_k(h0, h1).unwrap(), None);
            close(
                law.specific_entropy_change_j_kg_k(h0, h1).unwrap().unwrap() * reference_mass,
                entropy(325.0) - entropy(295.0),
                2e-12,
            );
        }
    }
}

#[test]
fn g0_unsupported_capacity_basis_context_and_component_domains_refuse() {
    let a = card("mass-capacity", [claim(CP, &A)]);
    let (lo, hi) = (point(280.0), point(340.0));
    let ka = key(CP);
    let q = component(1, &a, &ka, HeatStorageAmount::MassKg(2.0), &lo, &hi);
    assert!(PreparedHeatStorage::try_new(&[]).is_err());
    assert!(PreparedHeatStorage::try_new(&[q, q]).is_err());
    for amount in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut invalid = q;
        invalid.amount = HeatStorageAmount::MassKg(amount);
        assert!(PreparedHeatStorage::try_new(&[invalid]).is_err());
    }
    let mut wrong_basis = q;
    wrong_basis.amount = HeatStorageAmount::ReferenceVolumeM3(1.0);
    assert!(PreparedHeatStorage::try_new(&[wrong_basis]).is_err());
    let changed_context = hi
        .clone()
        .with_quantity("specimen", QuantitySpec::dimensional(Dims::NONE), 2.0)
        .unwrap();
    let mut invalid = q;
    invalid.upper = &changed_context;
    assert!(PreparedHeatStorage::try_new(&[invalid]).is_err());
    invalid.upper = &lo;
    assert!(PreparedHeatStorage::try_new(&[invalid]).is_err());
    let wrong_temperature = hi
        .clone()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::TemperatureDifference),
            340.0,
        )
        .unwrap();
    invalid.upper = &wrong_temperature;
    assert!(PreparedHeatStorage::try_new(&[invalid]).is_err());

    let mut discrete = claim(CP, &A);
    discrete.interpolation = InterpolationPolicy::TabulatedOnly;
    let discrete = card("discrete-capacity", [discrete]);
    assert!(
        PreparedHeatStorage::try_new(&[component(
            1,
            &discrete,
            &ka,
            HeatStorageAmount::MassKg(1.0),
            &lo,
            &hi
        )])
        .is_err()
    );
    let negative = card(
        "nonpositive-capacity",
        [claim(CP, &[(280.0, 10.0), (340.0, -1.0)])],
    );
    assert!(
        PreparedHeatStorage::try_new(&[component(
            1,
            &negative,
            &ka,
            HeatStorageAmount::MassKg(1.0),
            &lo,
            &hi
        )])
        .is_err()
    );
    let distant = card(
        "disjoint-domain",
        [claim(CP, &[(350.0, 10.0), (360.0, 20.0)])],
    );
    let (far_lo, far_hi) = (point(350.0), point(360.0));
    assert!(
        PreparedHeatStorage::try_new(&[
            q,
            component(
                2,
                &distant,
                &ka,
                HeatStorageAmount::MassKg(1.0),
                &far_lo,
                &far_hi
            )
        ])
        .is_err()
    );
    let assembly = PreparedHeatStorage::try_new(&[q]).unwrap();
    for mass in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(HeatStorageThermalLaw::try_new(&assembly, mass, 300.0).is_err());
    }
    assert!(HeatStorageThermalLaw::try_new(&assembly, 2.0, 279.0).is_err());
    assert!(HeatStorageThermalLaw::try_new(&assembly, 3.0, 300.0).is_err());
}

fn compiled(slug: &str) -> MaterialCard {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "heat-storage-{slug}-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process: "original population and fixed caloric conditions".into(),
            revision: 0,
        },
        pack.claims().clone(),
        vec![],
    )
    .unwrap()
}
fn source_point(card: &MaterialCard, key: &PropertyKey, t: f64) -> QueryPoint {
    let claim = card
        .claims()
        .claims_ordered()
        .find(|(_, c)| &c.key == key)
        .unwrap()
        .1;
    let mut point = QueryPoint::new();
    for (axis, &(low, high)) in claim.validity.bounds() {
        let value = if axis == "temperature" {
            t
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
fn heated(storage: &PreparedHeatStorage, mass: f64, start: f64, power: f64, duration: f64) -> f64 {
    let law = HeatStorageThermalLaw::try_new(storage, mass, start)
        .unwrap()
        .with_sensible_entropy();
    let cell = DepositedCell {
        id: CellId(1),
        segment_id: 1,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: law.specific_storage_j_kg(start).unwrap(),
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.1, 0.0, 0.0],
            volume_m3: 0.003,
            surface_area_m2: 0.1,
        },
        boundary: Boundary {
            ambient_k: start,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: power,
        },
    };
    law.admit_cell(&cell).unwrap();
    let mut wrong_mass = cell.clone();
    wrong_mass.mass_kg *= 1.01;
    assert!(law.admit_cell(&wrong_mass).is_err());
    let config = StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-8,
    };
    let initial = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(1), &law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config, || Ok(()))
    .unwrap()
    .0;
    let (final_state, audit) = initial.advance_to(duration, config, || Ok(())).unwrap();
    close(audit.external_heat_j, power * duration, 1e-10);
    close(audit.stored_energy_change_j, power * duration, 1e-10);
    close(audit.residual_j, 0.0, 1e-9);
    assert_eq!(final_state.mass_kg(), mass);
    assert!(
        initial
            .advance_to(duration, config, || Err(DepositionError::Cancelled))
            .is_err()
    );
    close(initial.temperature_k(CellId(1)).unwrap(), start, 0.0);
    final_state.temperature_k(CellId(1)).unwrap()
}

#[test]
fn g1_g3_real_native_and_mass_capacity_sources_drive_the_existing_heat_timeline() {
    let mdf = compiled("fsri-mdf-unconditioned-2021");
    let countertop = compiled("fsri-acrylic-countertop-dried-2021");
    let k = key(CVOL);
    let (m0, m1) = (
        source_point(&mdf, &k, 293.15),
        source_point(&mdf, &k, 303.15),
    );
    let (c0, c1) = (
        source_point(&countertop, &k, 293.15),
        source_point(&countertop, &k, 303.15),
    );
    let assembly = PreparedHeatStorage::try_new(&[
        component(
            1,
            &mdf,
            &k,
            HeatStorageAmount::ReferenceVolumeM3(0.002),
            &m0,
            &m1,
        ),
        component(
            2,
            &countertop,
            &k,
            HeatStorageAmount::ReferenceVolumeM3(0.001),
            &c0,
            &c1,
        ),
    ])
    .unwrap();
    // Original native capacity ordinates, not Cp times an invented density.
    // MDF's local k starts at 303.15 K; sensible storage still supports 293.15 K.
    let c_start = 0.002 * 886846.0 + 0.001 * 2118377.3333333335;
    let slope = (0.002 * (948692.0 - 886846.0) + 0.001 * (2186227.0 - 2118377.3333333335)) / 10.0;
    let energy: f64 = 600.0;
    let expected =
        293.15 + 2.0 * energy / (c_start + (c_start * c_start + 2.0 * slope * energy).sqrt());
    for mass_basis in [1.0, 10.0] {
        close(
            heated(&assembly, mass_basis, 293.15, 30.0, 20.0),
            expected,
            1e-12,
        );
    }
    assert_eq!(assembly.components().len(), 2);

    let beryllium = compiled("beryllium-nist-thermal-reference");
    assert!(
        !beryllium
            .claims()
            .claims_ordered()
            .any(|(_, c)| ["density", "thermal-conductivity"].contains(&c.key.name()))
    );
    let cp = key(CP);
    let (low, high) = (
        source_point(&beryllium, &cp, 200.0),
        source_point(&beryllium, &cp, 201.157448185),
    );
    let metal = PreparedHeatStorage::try_new(&[component(
        1,
        &beryllium,
        &cp,
        HeatStorageAmount::MassKg(0.1),
        &low,
        &high,
    )])
    .unwrap();
    // Two literal sampled source-fit ordinates bound this exact affine segment.
    let cp_start = 1097.364623976881;
    let cp_slope = (1108.949687212369 - cp_start) / (201.157448185 - 200.0);
    let specific_energy: f64 = 500.0;
    let expected = 200.0
        + 2.0 * specific_energy
            / (cp_start + (cp_start * cp_start + 2.0 * cp_slope * specific_energy).sqrt());
    close(heated(&metal, 0.1, 200.0, 10.0, 5.0), expected, 1e-12);
    for (source_card, component) in [
        (&mdf, &assembly.components()[0]),
        (&countertop, &assembly.components()[1]),
        (&beryllium, &metal.components()[0]),
    ] {
        for (_, answer) in &component.source_line().knots {
            source_card
                .claims()
                .verify_receipt(&answer.receipt)
                .unwrap();
        }
    }
}
