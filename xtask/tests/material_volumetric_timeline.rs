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
fn analytic_card(c0: f64, slope: f64, k0: f64, kslope: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (i, (base, derivative)) in [(c0, slope), (k0, kslope)].into_iter().enumerate() {
        claims
            .insert_claim(PropertyClaim {
                key: keys()[i].clone(),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [290.0, 300.0, 320.0, 340.0]
                        .into_iter()
                        .map(|t| (t, base + derivative * (t - 300.0)))
                        .collect(),
                    dims: [C, Dims([1, 1, -3, -1, 0, 0])][i],
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    290.0,
                    340.0,
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
