//! G1/G3: original volumetric source curves drive the existing thermal timeline.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, NormalizedPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::volumetric::{
    PreparedVolumetricThermalTransport as PathModel, VOLUMETRIC_HEAT_CAPACITY_DIMS as C,
    VolumetricThermalPathQuery,
    timeline::{FixedVolumeConductionLink, FixedVolumeThermalLaw},
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
    backward_euler::{BackwardEulerStep, StepCoefficients},
};
use std::{collections::BTreeMap, path::Path, process::Command};

const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const CONFIG: StepConfig = StepConfig {
    max_dt_s: 1.0,
    energy_tolerance_j: 1e-8,
};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn keys() -> [PropertyKey; 2] {
    [
        PropertyKey::new("volumetric-heat-capacity", C),
        PropertyKey::with_quantity(
            "thermal-conductivity",
            quantity(QuantityKind::ThermalConductivity),
        ),
    ]
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
}
fn state(name: &str) -> MaterialStateId {
    MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "fixed source context".into(),
        revision: 0,
    }
}
fn path(card: &MaterialCard, lower: &QueryPoint, upper: &QueryPoint) -> PathModel {
    PathModel::try_new(
        card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &keys()[0],
            conductivity: &keys()[1],
            temperature_axis: "temperature",
            lower,
            upper,
            selections: [POLICY; 2],
        },
    )
    .unwrap()
}
fn tabulated_card(capacity: &[(f64, f64)], conductivity: &[(f64, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (i, knots) in [capacity, conductivity].into_iter().enumerate() {
        claims
            .insert_claim(PropertyClaim {
                key: keys()[i].clone(),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: knots.to_vec(),
                    dims: [C, Dims([1, 1, -3, -1, 0, 0])][i],
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    knots[0].0,
                    knots.last().unwrap().0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored analytic fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(state("analytic volumetric law"), claims, vec![]).unwrap()
}
fn analytic_card(c0: f64, slope: f64, k0: f64, kslope: f64) -> MaterialCard {
    let knots = |base: f64, derivative: f64| {
        [290.0, 300.0, 320.0, 340.0].map(|t| (t, base + derivative * (t - 300.0)))
    };
    tabulated_card(&knots(c0, slope), &knots(k0, kslope))
}

fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.15e} != {expected:.15e}"
    );
}
fn cell(law: &FixedVolumeThermalLaw<'_>, id: u64, t: f64, volume: f64) -> DepositedCell {
    let cell = DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: MaterialId(id),
        mass_kg: law.reference_mass_kg(volume).unwrap(),
        specific_h_j_kg: law.specific_storage_j_kg(t).unwrap(),
        geometry: CellGeometry {
            start_m: [id as f64 * 0.01, 0.0, 0.0],
            end_m: [id as f64 * 0.01, 0.1, 0.0],
            volume_m3: volume,
            surface_area_m2: 0.1,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    law.admit_cell(&cell).unwrap();
    cell
}
fn timeline<'a>(
    cells: Vec<DepositedCell>,
    laws: &[(MaterialId, &'a dyn ThermalLaw)],
    link: Option<ThermalLink>,
) -> DepositionSimulation<'a> {
    let mut schedule: Vec<_> = cells
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    if let Some(link) = link {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(link),
        });
    }
    DepositionSimulation::new(
        schedule,
        laws,
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, CONFIG, || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_g3_fixed_volume_storage_survives_mass_basis_and_datum_changes() {
    let card = analytic_card(2e6, 20000.0, 0.2, 0.002);
    let p = path(&card, &point(290.0), &point(340.0));
    for rho in [200.0, 1000.0] {
        for datum in [290.0, 300.0, 315.0, 340.0] {
            let law = FixedVolumeThermalLaw::try_new(&p, rho, datum).unwrap();
            assert_eq!(law.path(), &p);
            assert_eq!(law.reference_density_kg_m3(), rho);
            assert_eq!(law.datum_temperature_k(), datum);
            let primitive = |t: f64| {
                let x = t - 300.0;
                2e6 * x + 10000.0 * x * x
            };
            for t in [290.0, 297.25, 300.0, 315.0, 320.0, 337.5, 340.0] {
                let h = law.specific_storage_j_kg(t).unwrap();
                close(h * rho, primitive(t) - primitive(datum), 4e-8);
                close(law.temperature_k(h).unwrap(), t, 6e-14);
                close(
                    law.reference_mass_kg(0.002).unwrap() * h,
                    0.002 * (primitive(t) - primitive(datum)),
                    1e-10,
                );
            }
            assert!(law.maximum_temperature_slope() >= rho / 1.8e6);
            assert!(law.maximum_temperature_slope() <= (rho / 1.8e6).next_up());
            let [lo, hi] = law.specific_storage_span_j_kg();
            assert!(law.temperature_k(lo.next_down()).is_err());
            assert!(law.temperature_k(hi.next_up()).is_err());
            assert!(law.temperature_k(f64::NAN).is_err());
            let mut wrong = cell(&law, 1, 310.0, 0.002);
            wrong.mass_kg *= 1.01;
            assert!(law.admit_cell(&wrong).is_err());
            assert!(law.reference_mass_kg(0.0).is_err());
            assert!(law.reference_mass_kg(f64::MAX).is_err());
        }
    }
    let first = FixedVolumeThermalLaw::try_new(&p, 200.0, 300.0).unwrap();
    for (rho, t) in [(400.0, 300.0), (200.0, 310.0)] {
        assert_ne!(
            first.identity(),
            FixedVolumeThermalLaw::try_new(&p, rho, t)
                .unwrap()
                .identity()
        );
    }
    for rho in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(FixedVolumeThermalLaw::try_new(&p, rho, 300.0).is_err());
    }
    assert!(FixedVolumeThermalLaw::try_new(&p, 200.0, 289.0).is_err());
}

#[test]
fn g1_g3_source_bound_k_and_volumetric_storage_share_one_implicit_heat_balance() {
    let a = analytic_card(2e6, 20000.0, 0.2, 0.002);
    let b = analytic_card(1e6, 10000.0, 0.4, -0.001);
    let pa = path(&a, &point(290.0), &point(340.0));
    let pb = path(&b, &point(290.0), &point(340.0));
    let la = FixedVolumeThermalLaw::try_new(&pa, 200.0, 300.0).unwrap();
    let lb = FixedVolumeThermalLaw::try_new(&pb, 1000.0, 300.0).unwrap();
    let raw_link = ThermalLink {
        id: 7,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.01,
        conductance_w_k: 0.0,
    };
    let link =
        FixedVolumeConductionLink::try_aligned(&raw_link, [&la, &lb], [0.01, 0.01], 0.02).unwrap();
    assert_eq!(link.id(), raw_link.id);
    let initial = timeline(
        vec![cell(&la, 1, 300.0, 0.001), cell(&lb, 2, 330.0, 0.001)],
        &[(MaterialId(1), &la), (MaterialId(2), &lb)],
        Some(raw_link),
    );
    let initial_energy = initial.energy_j();
    let initial_g = link.conductance_w_k(&initial).unwrap();
    close(initial_g, 0.01 / (0.01 / 0.2 + 0.01 / 0.37 + 0.02), 1e-16);

    // Independent analytic inverse of affine C(T), then a bounded scalar solve
    // of the two-cell backward-Euler balance. No production inverse or k sampler
    // supplies this oracle. The existing owner's local and global residuals
    // independently admit the resulting energy increments.
    let inverse =
        |e: f64, c0: f64, slope: f64| 300.0 + 2.0 * e / (c0 + (c0 * c0 + 2.0 * slope * e).sqrt());
    let endpoint = |q: f64| {
        let ta = inverse(q / 0.001, 2e6, 20000.0);
        let tb = inverse((34500.0 - q) / 0.001, 1e6, 10000.0);
        let ka = 0.2 + 0.002 * (ta - 300.0);
        let kb = 0.4 - 0.001 * (tb - 300.0);
        let g = 0.01 / (0.01 / ka + 0.01 / kb + 0.02);
        (ta, tb, g, q - 5000.0 * g * (tb - ta))
    };
    let (mut low, mut high) = (0.0, 34500.0);
    for _ in 0..80 {
        let q = 0.5 * (low + high);
        if endpoint(q).3 > 0.0 {
            high = q;
        } else {
            low = q;
        }
    }
    let q = 0.5 * (low + high);
    let (ta, tb, g, _) = endpoint(q);
    let step = BackwardEulerStep::new(&initial, 5000.0, || Ok(())).unwrap();
    let evaluate = |increments: &[f64]| {
        step.evaluate(
            increments,
            |candidate| {
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::from([(7, link.conductance_w_k(candidate)?)]),
                    payload: (),
                })
            },
            || Ok(()),
        )
    };
    let trial = evaluate(&[q, -q]).unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-8);
    let (accepted, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
    close(accepted.temperature_k(CellId(1)).unwrap(), ta, 6e-14);
    close(accepted.temperature_k(CellId(2)).unwrap(), tb, 6e-14);
    close(link.conductance_w_k(&accepted).unwrap(), g, 1e-16);
    close(accepted.energy_j(), initial_energy, 1e-8);
    close(audit.external_heat_j, 0.0, 0.0);
    close(audit.residual_j, 0.0, 1e-8);
    assert!((q - 5000.0 * initial_g * (tb - ta)).abs() > 10.0);
    assert_eq!(accepted.mass_kg(), initial.mass_kg());
    assert_eq!(accepted.links()[&7], raw_link); // Trial coefficients are not authored state.
    assert!(
        evaluate(&[q + 1.0, -q - 1.0])
            .unwrap()
            .accept(1e-8, 1e-8, || Ok(()))
            .is_err()
    );
    assert!(evaluate(&[1e9, -1e9]).is_err());
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.energy_j(), initial_energy);

    let identical_but_unbound = la.clone();
    let wrong = FixedVolumeConductionLink::try_aligned(
        &raw_link,
        [&identical_but_unbound, &lb],
        [0.01; 2],
        0.02,
    )
    .unwrap();
    assert!(wrong.conductance_w_k(&initial).is_err());
    let changed_geometry = FixedVolumeConductionLink::try_aligned(
        &ThermalLink {
            area_m2: 0.02,
            ..raw_link
        },
        [&la, &lb],
        [0.01; 2],
        0.02,
    )
    .unwrap();
    assert!(changed_geometry.conductance_w_k(&initial).is_err());
    let mut wrong_cell = cell(&la, 1, 300.0, 0.001);
    wrong_cell.mass_kg *= 2.0;
    let wrong_state = timeline(
        vec![wrong_cell, cell(&lb, 2, 330.0, 0.001)],
        &[(MaterialId(1), &la), (MaterialId(2), &lb)],
        Some(raw_link),
    );
    assert!(link.conductance_w_k(&wrong_state).is_err());
    assert!(
        FixedVolumeConductionLink::try_aligned(&raw_link, [&la, &lb], [0.0, 0.01], 0.0).is_err()
    );
    assert!(
        FixedVolumeConductionLink::try_aligned(&raw_link, [&la, &lb], [0.01; 2], -1.0).is_err()
    );
}

fn compiled(slug: &str) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "volumetric-timeline-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let raw = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    let pack = NormalizedMaterialCardPack::new(state(slug), raw).unwrap();
    NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &pack.to_bytes()).unwrap()
}
fn source_point(card: &MaterialCard, t: f64) -> QueryPoint {
    let mut point = QueryPoint::new();
    for key in keys() {
        let (_, claim) = card
            .claims()
            .claims_ordered()
            .find(|(_, c)| c.key == key)
            .unwrap();
        for (name, &(low, high)) in claim.validity.bounds() {
            let value = if name == "temperature" {
                t
            } else {
                assert_eq!(low, high);
                low
            };
            point = point
                .with_quantity(name, claim.validity.axis_quantities()[name], value)
                .unwrap();
        }
    }
    point
}

fn sourced_heating(slug: &str, capacity_at_20_30_c: [f64; 2], reference_ratios: [f64; 2]) {
    let pack = compiled(slug);
    let card = pack.card();
    let p = path(
        card,
        &source_point(card, 288.15),
        &source_point(card, 313.15),
    );
    // This deliberately authored reference ratio is NOT FSRI's normalization
    // density. The fixed-volume result is invariant when m and h use it together.
    let mut outcomes = Vec::new();
    for rho in reference_ratios {
        let law = FixedVolumeThermalLaw::try_new(&p, rho, 293.15).unwrap();
        let mut heated_cell = cell(&law, 1, 293.15, 0.002);
        heated_cell.boundary.source_w = 30.0;
        let initial = timeline(vec![heated_cell], &[(MaterialId(1), &law)], None);
        let (heated, audit) = initial.advance_to(20.0, CONFIG, || Ok(())).unwrap();
        // Literal FSRI capacity means at 20/30 C, integrated independently.
        let c0 = capacity_at_20_30_c[0];
        let slope = (capacity_at_20_30_c[1] - c0) / 10.0;
        let delta = 2.0 * 300000.0 / (c0 + (c0 * c0 + 2.0 * slope * 300000.0).sqrt());
        let t = heated.temperature_k(CellId(1)).unwrap();
        close(t, 293.15 + delta, 1e-12);
        close(audit.external_heat_j, 600.0, 1e-10);
        close(audit.stored_energy_change_j, 600.0, 1e-10);
        close(audit.residual_j, 0.0, 1e-9);
        assert_eq!(heated.time_s(), 20.0);
        assert_eq!(heated.mass_kg(), initial.mass_kg());
        assert!(
            initial
                .advance_to(20.0, CONFIG, || Err(DepositionError::Cancelled))
                .is_err()
        );
        assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 293.15);
        outcomes.push(t);
    }
    assert_eq!(outcomes[0], outcomes[1]);
    for source in p.source_lines() {
        for (_, answer) in &source.knots {
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
}

#[test]
fn g1_g3_real_foam_capacity_drives_original_timeline_heating_without_measured_density() {
    sourced_heating(
        "fsri-pe-foam-pipe-unconditioned-2022",
        [80182.33333333333, 84950.0],
        [35.0, 70.0],
    );
}

#[test]
fn g1_g3_real_frp_capacity_drives_original_timeline_without_a_density_or_cp_claim() {
    sourced_heating(
        "fsri-frp-panel-unconditioned-2022",
        [
            (1546098.0 + 1537905.0 + 1552146.0) / 3.0,
            (1653830.0 + 1628139.0 + 1666615.0) / 3.0,
        ],
        [1000.0, 2000.0],
    );
}

#[test]
fn g1_sensible_entropy_integrates_source_breaks_and_resolves_small_changes() {
    // Independent slope/intercept primitive, including decreasing capacity.
    for slope in [-20000.0, 0.0, 20000.0] {
        let card = analytic_card(2e6, slope, 0.2, 0.001);
        let p = path(&card, &point(290.0), &point(340.0));
        let original = p.clone();
        for (a, b) in [(290.0, 340.0), (297.25, 330.5), (340.0, 292.5)] {
            let expected = slope * (b - a) + (2e6 - 300.0 * slope) * (b / a).ln();
            let actual = p.sensible_entropy_change_j_m3_k(a, b).unwrap();
            close(actual, expected, 2e-9);
            assert_eq!(p.sensible_entropy_change_j_m3_k(b, a).unwrap(), -actual);
        }
        // Adjacent floats straddle an actual retained knot. The integral is
        // bounded independently by endpoint C and reciprocal temperature.
        let a = 300.0_f64.next_down();
        let b = 300.0_f64.next_up();
        let capacity_a = 2e6 + slope * (a - 300.0);
        let capacity_b = 2e6 + slope * (b - 300.0);
        let entropy = p.sensible_entropy_change_j_m3_k(a, b).unwrap();
        let lower = capacity_a.min(capacity_b) / b * (b - a);
        let upper = capacity_a.max(capacity_b) / a * (b - a);
        assert!(entropy > 0.0 && entropy >= lower * (1.0 - 8.0 * f64::EPSILON));
        assert!(entropy <= upper * (1.0 + 8.0 * f64::EPSILON));
        close(
            entropy,
            p.sensible_entropy_change_j_m3_k(a, 300.0).unwrap()
                + p.sensible_entropy_change_j_m3_k(300.0, b).unwrap(),
            1e-23,
        );
        assert_eq!(p.sensible_entropy_change_j_m3_k(300.0, 300.0).unwrap(), 0.0);
        for outside in [0.0, 289.0, 341.0, f64::INFINITY, f64::NAN] {
            assert!(p.sensible_entropy_change_j_m3_k(outside, 300.0).is_err());
            assert!(p.sensible_entropy_change_j_m3_k(300.0, outside).is_err());
        }
        assert_eq!(p, original);
    }
    // Capacity and conductivity have DIFFERENT breaks; cross all of them.
    let card = tabulated_card(
        &[(290.0, 1e6), (310.0, 3e6), (340.0, 1.5e6)],
        &[(290.0, 0.1), (300.0, 0.2), (320.0, 0.3), (340.0, 0.4)],
    );
    let p = path(&card, &point(290.0), &point(340.0));
    let expected = 100000.0 * 20.0 - 28e6 * (310.0_f64 / 290.0).ln() - 50000.0 * 30.0
        + 18.5e6 * (340.0_f64 / 310.0).ln();
    close(
        p.sensible_entropy_change_j_m3_k(290.0, 340.0).unwrap(),
        expected,
        1e-8,
    );
    // The entropy can be representable even when T_high/T_low and the energy
    // integral overflow. This is a source-path query, not an admitted timeline.
    let card = tabulated_card(
        &[(1e-300, 1e100), (1e300, 1e100)],
        &[(1e-300, 1e-100), (1e300, 1e-100)],
    );
    let p = path(&card, &point(1e-300), &point(1e300));
    assert!(p.storage_change_j_m3(1e-300, 1e300).is_err());
    let actual = p.sensible_entropy_change_j_m3_k(1e-300, 1e300).unwrap();
    close(actual / 1e100, 600.0 * 10.0_f64.ln(), 1e-12);
}

#[test]
fn g1_g3_opt_in_entropy_closes_the_existing_implicit_conduction_account() {
    let affine = analytic_card(2e6, 20000.0, 0.2, 0.002);
    let affine_path = path(&affine, &point(290.0), &point(340.0));
    let expected = 20000.0 * (337.5 - 297.25) - 4e6 * (337.5_f64 / 297.25).ln();
    for rho in [200.0, 1000.0] {
        for datum in [290.0, 300.0, 315.0, 340.0] {
            let plain = FixedVolumeThermalLaw::try_new(&affine_path, rho, datum).unwrap();
            let law = plain.clone().with_sensible_entropy();
            assert!(!plain.sensible_entropy_enabled());
            assert!(law.sensible_entropy_enabled());
            assert_ne!(plain.identity(), law.identity());
            assert_eq!(
                law.identity(),
                law.clone().with_sensible_entropy().identity()
            );
            assert_eq!(plain.path(), law.path());
            let hi = law.specific_storage_j_kg(297.25).unwrap();
            let hf = law.specific_storage_j_kg(337.5).unwrap();
            assert_eq!(plain.specific_entropy_change_j_kg_k(hi, hf).unwrap(), None);
            close(
                law.specific_entropy_change_j_kg_k(hi, hf).unwrap().unwrap() * rho,
                expected,
                3e-9,
            );
            assert!(law.specific_entropy_change_j_kg_k(f64::NAN, hf).is_err());
            assert!(
                law.specific_entropy_change_j_kg_k(
                    hi,
                    law.specific_storage_span_j_kg()[1].next_up()
                )
                .is_err()
            );
        }
    }
    let card = analytic_card(1e6, 0.0, 0.5, 0.0);
    let p = path(&card, &point(290.0), &point(340.0));
    let plain = FixedVolumeThermalLaw::try_new(&p, 1000.0, 300.0).unwrap();
    let unavailable = timeline(
        vec![cell(&plain, 1, 300.0, 0.001)],
        &[(MaterialId(1), &plain)],
        None,
    );
    assert!(
        unavailable
            .advance_to_with_entropy(1.0, CONFIG, || Ok(()))
            .is_err()
    );
    assert_eq!(unavailable.time_s(), 0.0);
    let law = plain.with_sensible_entropy();
    // Distinct stored coordinates can round to one T: do not erase their heat.
    let hi = law.specific_storage_j_kg(310.0).unwrap();
    let hf = hi.next_up();
    let t = law.temperature_k(hi).unwrap();
    assert_eq!(t, law.temperature_k(hf).unwrap());
    let ds = law.specific_entropy_change_j_kg_k(hi, hf).unwrap().unwrap();
    assert_eq!(ds, (hf - hi) / t);
    assert!(ds > 0.0);
    assert_eq!(
        law.specific_entropy_change_j_kg_k(hf, hi).unwrap(),
        Some(-ds)
    );

    let raw_link = ThermalLink {
        id: 7,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.01,
        conductance_w_k: 0.0,
    };
    let link =
        FixedVolumeConductionLink::try_aligned(&raw_link, [&law; 2], [0.01; 2], 0.0).unwrap();
    let initial = timeline(
        vec![cell(&law, 1, 300.0, 0.001), cell(&law, 2, 330.0, 0.001)],
        &[(MaterialId(1), &law), (MaterialId(2), &law)],
        Some(raw_link),
    );
    // Closed-form backward Euler for two constant 1000 J/K cells, independently
    // using k=0.5, area=0.01 and a total conduction length of 0.02 metres.
    let g = 0.5 * 0.01 / 0.02;
    let duration = 2000.0;
    let q = duration * g * 30.0 / (1.0 + duration * g * 2.0 / 1000.0);
    let ta = 300.0 + q / 1000.0;
    let tb = 330.0 - q / 1000.0;
    let step = BackwardEulerStep::new(&initial, duration, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_entropy(
            &[q, -q],
            |candidate| {
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::from([(7, link.conductance_w_k(candidate)?)]),
                    payload: (),
                })
            },
            || Ok(()),
        )
        .unwrap();
    let entropy = *trial.entropy_audit().unwrap();
    let state_entropy = 1000.0 * ((ta / 300.0_f64).ln() + (tb / 330.0_f64).ln());
    let production = q * (1.0 / ta - 1.0 / tb);
    close(entropy.thermal_storage_change_j_k, state_entropy, 2e-12);
    close(entropy.conduction_generation_j_k, production, 2e-14);
    close(
        entropy.temporal_defect_j_k,
        state_entropy - production,
        2e-12,
    );
    assert!(entropy.thermal_storage_change_j_k > 0.0);
    assert!(entropy.conduction_generation_j_k > 0.0);
    assert!(entropy.temporal_defect_j_k > 0.0);
    assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
    assert_eq!(entropy.unclassified_source_at_cell_j_k, 0.0);
    assert!(!entropy.has_unclassified_sources);
    let (accepted, (), energy) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
    close(accepted.temperature_k(CellId(1)).unwrap(), ta, 6e-14);
    close(accepted.temperature_k(CellId(2)).unwrap(), tb, 6e-14);
    close(accepted.energy_j(), initial.energy_j(), 1e-8);
    assert_eq!(energy.external_heat_j, 0.0);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_g3_sourced_entropy_retains_native_capacity_and_unclassified_heat() {
    // Literal primary FSRI means at 20/30 C; no query supplies the oracle.
    for (slug, capacities, ratios) in [
        (
            "fsri-pp-unconditioned-2021",
            [1563841.0, 1641098.0],
            [800.0, 1600.0],
        ),
        (
            "fsri-puf-mattress-pad-unconditioned-2022",
            [40175.666666666664, 34740.333333333336],
            [50.0, 100.0],
        ),
    ] {
        let pack = compiled(slug);
        let card = pack.card();
        let p = path(
            card,
            &source_point(card, 288.15),
            &source_point(card, 313.15),
        );
        let original = p.clone();
        let c0 = capacities[0];
        let slope = (capacities[1] - c0) / 10.0;
        let delta_at = |heat_j: f64| {
            let density = heat_j / 0.002;
            2.0 * density / (c0 + (c0 * c0 + 2.0 * slope * density).sqrt())
        };
        let delta = delta_at(600.0);
        assert!((0.0..10.0).contains(&delta));
        let expected = 0.002 * (slope * delta + (c0 - slope * 293.15) * (delta / 293.15).ln_1p());
        let source_at_cell: f64 = (0..20)
            .map(|i| 30.0 / (293.15 + delta_at(30.0 * f64::from(i))))
            .sum();
        for rho in ratios {
            // This ratio is an authored reference mass basis, not a density measurement.
            let law = FixedVolumeThermalLaw::try_new(&p, rho, 293.15)
                .unwrap()
                .with_sensible_entropy();
            let mut heated_cell = cell(&law, 1, 293.15, 0.002);
            heated_cell.boundary.source_w = 30.0;
            let initial = timeline(vec![heated_cell], &[(MaterialId(1), &law)], None);
            let (heated, energy, entropy) = initial
                .advance_to_with_entropy(20.0, CONFIG, || Ok(()))
                .unwrap();
            let (ordinary, ordinary_energy) = initial.advance_to(20.0, CONFIG, || Ok(())).unwrap();
            assert_eq!(
                heated.cells()[&CellId(1)].specific_h_j_kg,
                ordinary.cells()[&CellId(1)].specific_h_j_kg
            );
            assert_eq!(energy, ordinary_energy);
            close(
                heated.temperature_k(CellId(1)).unwrap(),
                293.15 + delta,
                1e-12,
            );
            close(energy.stored_energy_change_j, 600.0, 1e-9);
            close(entropy.thermal_storage_change_j_k, expected, 1e-10);
            close(
                entropy.unclassified_source_at_cell_j_k,
                source_at_cell,
                1e-12,
            );
            close(
                entropy.temporal_defect_j_k,
                expected - source_at_cell,
                1e-10,
            );
            assert!(entropy.has_unclassified_sources);
            assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
            assert_eq!(entropy.conduction_generation_j_k, 0.0);
            assert!(entropy.temporal_defect_j_k < 0.0);
            assert!(
                initial
                    .advance_to_with_entropy(20.0, CONFIG, || Err(DepositionError::Cancelled))
                    .is_err()
            );
            assert_eq!(initial.time_s(), 0.0);
        }
        assert_eq!(p, original);
        for source in p.source_lines() {
            for (_, answer) in &source.knots {
                card.claims().verify_receipt(&answer.receipt).unwrap();
            }
        }
    }
}

const LAYER_CONFIG: fs_material::thermal::conduction::LayerConductionConfig =
    fs_material::thermal::conduction::LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-12,
        relative_heat_tolerance: 1e-12,
    };

#[test]
fn g1_integrated_layers_resolve_affine_conduction_interfaces_and_tangents() {
    let card = analytic_card(2e6, 20000.0, 0.2, 0.01);
    let p = path(&card, &point(290.0), &point(340.0));
    let law = FixedVolumeThermalLaw::try_new(&p, 1000.0, 300.0).unwrap();
    let raw = ThermalLink {
        id: 7,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.01,
        conductance_w_k: 0.0,
    };
    let sampled = FixedVolumeConductionLink::try_aligned(&raw, [&law; 2], [0.01; 2], 0.0).unwrap();
    let integrated = sampled.clone().into_integrated_layers().unwrap();
    assert_eq!(integrated.id(), raw.id);
    assert!(
        integrated
            .source_paths()
            .into_iter()
            .all(|source| std::ptr::eq(source, &p))
    );
    let initial = timeline(
        vec![cell(&law, 1, 300.0, 0.001), cell(&law, 2, 330.0, 0.001)],
        &[(MaterialId(1), &law), (MaterialId(2), &law)],
        Some(raw),
    );
    let result = integrated
        .evaluate(&initial, LAYER_CONFIG, || Ok(()))
        .unwrap();
    // Independent Kirchhoff primitive: integral_300^330 [0.2+0.01*(T-300)] dT = 10.5.
    // Equal layer lengths put half that integral on either side of the interface.
    let interface = 300.0 + 2.0 * 5.25 / (0.2 + (0.2_f64 * 0.2 + 2.0 * 0.01 * 5.25).sqrt());
    close(result.power_a_to_b_w, -5.25, 2e-11);
    close(result.conductance_w_k, 0.175, 1e-12);
    for t in result.interface_temperature_k {
        close(t, interface, 2e-11);
    }
    close(result.endpoint_derivative_w_k[0], 0.1, 2e-13);
    close(result.endpoint_derivative_w_k[1], -0.25, 2e-13);
    close(result.layer_heat_residual_w, 0.0, 2e-11);
    // The original cell-centre approximation remains an explicit, distinct choice.
    let sampled_power = -30.0 * sampled.conductance_w_k(&initial).unwrap();
    close(
        sampled_power,
        -30.0 * 0.01 / (0.01 / 0.2 + 0.01 / 0.5),
        1e-14,
    );
    assert!((sampled_power - result.power_a_to_b_w).abs() > 0.9);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.links()[&7], raw);

    let card = analytic_card(1e6, 0.0, 0.4, 0.0);
    let p = path(&card, &point(290.0), &point(340.0));
    let law = FixedVolumeThermalLaw::try_new(&p, 1000.0, 300.0).unwrap();
    let second_card = analytic_card(1.5e6, 0.0, 0.8, 0.0);
    let second_path = path(&second_card, &point(290.0), &point(340.0));
    let second_law = FixedVolumeThermalLaw::try_new(&second_path, 500.0, 300.0).unwrap();
    let raw = ThermalLink {
        area_m2: 0.02,
        ..raw
    };
    let layers =
        FixedVolumeConductionLink::try_aligned(&raw, [&law, &second_law], [0.01, 0.02], 0.05)
            .unwrap()
            .into_integrated_layers()
            .unwrap();
    let sources = layers.source_paths();
    assert!(std::ptr::eq(sources[0], &p));
    assert!(std::ptr::eq(sources[1], &second_path));
    let initial = timeline(
        vec![
            cell(&law, 1, 330.0, 0.001),
            cell(&second_law, 2, 300.0, 0.001),
        ],
        &[(MaterialId(1), &law), (MaterialId(2), &second_law)],
        Some(raw),
    );
    let result = layers.evaluate(&initial, LAYER_CONFIG, || Ok(())).unwrap();
    // Constant-k series resistance: R_area = .01/.4 + .02/.8 + .05 = .1 m²K/W.
    close(result.power_a_to_b_w, 6.0, 2e-11);
    close(result.interface_temperature_k[0], 322.5, 2e-11);
    close(result.interface_temperature_k[1], 307.5, 2e-11);
    close(result.endpoint_derivative_w_k[0], 0.2, 2e-13);
    close(result.endpoint_derivative_w_k[1], -0.2, 2e-13);
    close(result.contact_temperature_residual_k, 0.0, 2e-11);
}

#[test]
fn g0_g4_integrated_layers_keep_source_binding_geometry_and_work_limits() {
    let card = analytic_card(2e6, 20000.0, 0.2, 0.01);
    let p = path(&card, &point(290.0), &point(340.0));
    let law = FixedVolumeThermalLaw::try_new(&p, 1000.0, 300.0).unwrap();
    let raw = ThermalLink {
        id: 7,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.01,
        conductance_w_k: 0.0,
    };
    let initial = timeline(
        vec![cell(&law, 1, 300.0, 0.001), cell(&law, 2, 330.0, 0.001)],
        &[(MaterialId(1), &law), (MaterialId(2), &law)],
        Some(raw),
    );
    let layers = FixedVolumeConductionLink::try_aligned(&raw, [&law; 2], [0.01; 2], 0.0)
        .unwrap()
        .into_integrated_layers()
        .unwrap();
    let energy = initial.energy_j();
    assert!(matches!(
        layers.evaluate(&initial, LAYER_CONFIG, || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    assert!(matches!(
        layers.evaluate(
            &initial,
            fs_material::thermal::conduction::LayerConductionConfig {
                max_iterations: 1,
                ..LAYER_CONFIG
            },
            || Ok(())
        ),
        Err(DepositionError::Budget(_))
    ));
    let zero_area = ThermalLink {
        area_m2: 0.0,
        ..raw
    };
    assert!(
        FixedVolumeConductionLink::try_aligned(&zero_area, [&law; 2], [0.01; 2], 0.0)
            .unwrap()
            .into_integrated_layers()
            .is_err()
    );
    let unbound_clone = law.clone();
    let wrong =
        FixedVolumeConductionLink::try_aligned(&raw, [&unbound_clone, &law], [0.01; 2], 0.0)
            .unwrap()
            .into_integrated_layers()
            .unwrap();
    assert!(wrong.evaluate(&initial, LAYER_CONFIG, || Ok(())).is_err());
    let altered_link = ThermalLink {
        area_m2: 0.02,
        ..raw
    };
    let wrong = FixedVolumeConductionLink::try_aligned(&altered_link, [&law; 2], [0.01; 2], 0.0)
        .unwrap()
        .into_integrated_layers()
        .unwrap();
    assert!(wrong.evaluate(&initial, LAYER_CONFIG, || Ok(())).is_err());
    let mut wrong_cell = cell(&law, 1, 300.0, 0.001);
    wrong_cell.mass_kg *= 1.1;
    let wrong = timeline(
        vec![wrong_cell, cell(&law, 2, 330.0, 0.001)],
        &[(MaterialId(1), &law), (MaterialId(2), &law)],
        Some(raw),
    );
    assert!(layers.evaluate(&wrong, LAYER_CONFIG, || Ok(())).is_err());

    // Both individual endpoint temperatures are supported, but neither prepared
    // layer covers the WHOLE 300..330 K interval required by the shared kernel.
    let pa = path(&card, &point(290.0), &point(320.0));
    let pb = path(&card, &point(310.0), &point(340.0));
    let la = FixedVolumeThermalLaw::try_new(&pa, 1000.0, 300.0).unwrap();
    let lb = FixedVolumeThermalLaw::try_new(&pb, 1000.0, 330.0).unwrap();
    let narrow = timeline(
        vec![cell(&la, 1, 300.0, 0.001), cell(&lb, 2, 330.0, 0.001)],
        &[(MaterialId(1), &la), (MaterialId(2), &lb)],
        Some(raw),
    );
    let sampled = FixedVolumeConductionLink::try_aligned(&raw, [&la, &lb], [0.01; 2], 0.0).unwrap();
    assert!(sampled.conductance_w_k(&narrow).is_ok());
    assert!(
        sampled
            .into_integrated_layers()
            .unwrap()
            .evaluate(&narrow, LAYER_CONFIG, || Ok(()))
            .is_err()
    );
    assert_eq!(initial.energy_j(), energy);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_g3_fsri_pp_layers_share_the_original_implicit_energy_and_entropy_account() {
    let pack = compiled("fsri-pp-unconditioned-2021");
    let card = pack.card();
    let p = path(
        card,
        &source_point(card, 288.15),
        &source_point(card, 313.15),
    );
    let original = p.clone();
    // Authored reference mass basis; this is not a measured density claim.
    let law = FixedVolumeThermalLaw::try_new(&p, 1000.0, 293.15)
        .unwrap()
        .with_sensible_entropy();
    let raw = ThermalLink {
        id: 7,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.01,
        conductance_w_k: 0.0,
    };
    let sampled = FixedVolumeConductionLink::try_aligned(&raw, [&law; 2], [0.01; 2], 0.0).unwrap();
    let layers = sampled.clone().into_integrated_layers().unwrap();
    let initial = timeline(
        vec![cell(&law, 1, 293.15, 0.001), cell(&law, 2, 303.15, 0.001)],
        &[(MaterialId(1), &law), (MaterialId(2), &law)],
        Some(raw),
    );
    // Literal FSRI Cvol means at 20/30 C and published k at 15/45 C.
    // The oracle integrates these affine laws independently of the prepared
    // storage inverse, layer solver, source sampler and implicit residual.
    let c0 = 1563841.0;
    let capacity_slope = (1641098.0 - c0) / 10.0;
    let k_slope = (0.151 - 0.137) / 30.0;
    let k20 = 0.137 + k_slope * 5.0;
    let initial_hot_energy = 0.001 * (c0 * 10.0 + 0.5 * capacity_slope * 100.0);
    let inverse = |energy: f64| {
        let ev = energy / 0.001;
        293.15 + 2.0 * ev / (c0 + (c0 * c0 + 2.0 * capacity_slope * ev).sqrt())
    };
    let duration = 5000.0;
    let endpoint = |q: f64| {
        let ta = inverse(q);
        let tb = inverse(initial_hot_energy - q);
        let integral_k = (tb - ta) * (k20 + 0.5 * k_slope * ((ta - 293.15) + (tb - 293.15)));
        (ta, tb, q - duration * 0.01 / 0.02 * integral_k)
    };
    let (mut lo, mut hi) = (0.0, 0.5 * initial_hot_energy);
    for _ in 0..80 {
        let q = 0.5 * (lo + hi);
        if endpoint(q).2 > 0.0 {
            hi = q;
        } else {
            lo = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let (ta, tb, _) = endpoint(q);
    let step = BackwardEulerStep::new(&initial, duration, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_entropy(
            &[q, -q],
            |candidate| {
                let response = layers.evaluate(candidate, LAYER_CONFIG, || Ok(()))?;
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::from([(7, response.conductance_w_k)]),
                    payload: response,
                })
            },
            || Ok(()),
        )
        .unwrap();
    close(trial.payload().power_a_to_b_w, -q / duration, 2e-11);
    assert!(trial.maximum_cell_residual_j() < 1e-7);
    let entropy = *trial.entropy_audit().unwrap();
    let ds = |a: f64, b: f64| {
        0.001 * (capacity_slope * (b - a) + (c0 - capacity_slope * 293.15) * ((b - a) / a).ln_1p())
    };
    let expected_entropy = ds(293.15, ta) + ds(303.15, tb);
    let generation = q * (1.0 / ta - 1.0 / tb);
    close(entropy.thermal_storage_change_j_k, expected_entropy, 2e-10);
    close(entropy.conduction_generation_j_k, generation, 2e-10);
    close(
        entropy.temporal_defect_j_k,
        expected_entropy - generation,
        2e-10,
    );
    assert!(entropy.thermal_storage_change_j_k > 0.0 && entropy.conduction_generation_j_k > 0.0);
    assert!(entropy.temporal_defect_j_k > 0.0);
    assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
    assert!(!entropy.has_unclassified_sources);
    let (accepted, response, energy) = trial.accept(1e-7, 1e-7, || Ok(())).unwrap();
    close(accepted.temperature_k(CellId(1)).unwrap(), ta, 1e-12);
    close(accepted.temperature_k(CellId(2)).unwrap(), tb, 1e-12);
    close(accepted.energy_j(), initial.energy_j(), 1e-8);
    close(response.power_a_to_b_w, -q / duration, 2e-11);
    close(energy.external_heat_j, 0.0, 0.0);
    assert!((q - duration * sampled.conductance_w_k(&accepted).unwrap() * (tb - ta)).abs() > 0.01);
    assert_eq!(accepted.links()[&7], raw);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(p, original);
}
