//! G2/G3 checks of the three retained figure extractions. The error boxes
//! describe acquisition only; they supply no experimental confidence or
//! constitutive/instrument qualification.

use super::{fixture_dir, run_compiler, workspace_path};
use fs_matdb::{
    MatDbError, NormalizedPack, PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_math::det;
use fs_qty::Dims;
use fs_vvreg::instrument_data::{
    ArtifactRef, Decision, Digitization, INSTRUMENT_DATA_LEDGER_LOCATOR, UncertaintyStatement,
    parse_instrument_data_ledger,
};
use std::{collections::BTreeMap, fs, time::Instant};

const BEAD: &str = "frankensim-irp-v5-instrument-realism-pexfp.354";

const SLUGS: [&str; 3] = [
    "music-wire-nbs-c447",
    "bronze-gunmetal-damping-rsic508",
    "brass-yellow-half-hard-damping-nasa-tn-d1467",
];

struct Read<'a> {
    claim: &'a str,
    pass: &'a str,
    fields: Vec<&'a str>,
}
impl Read<'_> {
    fn number(&self, index: usize) -> f64 {
        self.fields[index].parse::<f64>().unwrap()
    }
    fn value(&self, anchor0: f64, anchor1: f64, pixel: f64) -> f64 {
        let (v0, v1) = (self.number(11), self.number(13));
        (v0 + (pixel - anchor0) * (v1 - v0) / (anchor1 - anchor0)) * self.number(17)
    }
    fn bounds(&self) -> [f64; 2] {
        let (p0, p1, p) = (self.number(10), self.number(12), self.number(14));
        let (reading, anchor) = (self.number(15), self.number(16));
        assert!(reading >= 0.5 && anchor >= 0.5);
        assert!((p1 - p0).abs() > 2.0 * anchor);
        let mut bounds = [f64::INFINITY, f64::NEG_INFINITY];
        // The affine map is monotone in each uncertain coordinate while the
        // anchor intervals are disjoint. Its extrema occur at these corners.
        for a in [p0 - anchor, p0 + anchor] {
            for b in [p1 - anchor, p1 + anchor] {
                for q in [p - reading, p + reading] {
                    let value = self.value(a, b, q);
                    bounds[0] = bounds[0].min(value);
                    bounds[1] = bounds[1].max(value);
                }
            }
        }
        bounds
    }
    fn encloses(&self, nominal: f64, radius: f64) -> bool {
        let [lo, hi] = self.bounds();
        nominal - radius <= lo && hi <= nominal + radius
    }
}

fn property(claim: &str) -> &str {
    match claim {
        "tensile-strength-0p02in" => "tensile_strength_diameter_0p02in",
        "tensile-strength-0p18in" => "tensile_strength_diameter_0p18in",
        "specific-damping-capacity-5ksi" => "specific_damping_capacity_5ksi_shear",
        "specimen-loss-factor-10ksi-50hz" => "specimen_loss_factor_10ksi_50hz",
        "specimen-loss-factor-10ksi-500hz" => "specimen_loss_factor_10ksi_500hz",
        _ => panic!("unregistered source claim {claim}"),
    }
}

fn compile(slug: &str) -> NormalizedPack {
    let path = format!("data/matdb/seed-v1/{slug}/manifest.tsv");
    let output = fixture_dir().join(format!("{slug}.pack"));
    let run = run_compiler(&workspace_path(&path), &output);
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap()
}

fn query_scope(pack: &NormalizedPack, name: &str) -> (f64, QueryPoint) {
    let claims = pack.claims().claims_for(name);
    assert_eq!(claims.len(), 1);
    let claim = claims[0].1;
    let mut point = QueryPoint::new();
    for (axis, &(lo, hi)) in claim.validity.bounds() {
        assert_eq!(lo, hi, "only the source's discrete conditions are supplied");
        point = if let Some(quantity) = claim.validity.axis_quantities().get(axis) {
            point.with_quantity(axis, *quantity, lo)
        } else {
            point.with(axis, lo)
        }
        .unwrap();
    }
    let answer = pack
        .claims()
        .query(name, &point, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    // Each physical coordinate is live, including stress and wire diameter.
    // No value at a neighboring condition is manufactured from these scalars.
    for (axis, &(lo, _)) in claim.validity.bounds() {
        for outside in [lo.next_down(), lo.next_up()] {
            let trial = point.clone().with(axis, outside).unwrap();
            assert!(
                matches!(
                    pack.claims()
                        .query(name, &trial, SelectionPolicy::SingleClaimOnly),
                    Err(MatDbError::NoClaimInDomain { .. })
                ),
                "{name} silently crossed its {axis} scope"
            );
        }
    }
    (answer.evidence.value.value, point)
}

pub(super) fn source_bounds_reach_consumers() {
    let started = Instant::now();
    let ledger = parse_instrument_data_ledger(
        &fs::read_to_string(workspace_path(INSTRUMENT_DATA_LEDGER_LOCATOR)).unwrap(),
    )
    .unwrap();
    let mut reading_count = 0;
    let mut rejection_count = 0;
    for slug in SLUGS {
        let directory = format!("data/matdb/seed-v1/{slug}");
        let text =
            fs::read_to_string(workspace_path(&format!("{directory}/digitization.tsv"))).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines[0], "frankensim.figure-extraction.v1");
        let header = lines
            .iter()
            .position(|line| line.starts_with("claim\t"))
            .unwrap();
        let metadata = lines[1..header]
            .iter()
            .map(|line| line.split_once('\t').unwrap())
            .collect::<BTreeMap<_, _>>();
        assert!(metadata["repeatability"].contains("GoldHarbor passes A/B"));
        assert!(metadata["repeatability"].contains("RusticPlateau"));
        assert!(metadata["authority"].contains("fixture-only"));
        assert!(
            metadata["authority"].contains("experimental scatter and model discrepancy unstated")
        );
        assert!(metadata["bound_rule"].contains("eight"));
        for field in ["source_sha256", "retained_sha256"] {
            assert_eq!(metadata[field].len(), 64);
            assert!(metadata[field].bytes().all(|b| b.is_ascii_hexdigit()));
        }
        let page_path = format!("{directory}/{}", metadata["retained_file"]);
        let page = fs::read(workspace_path(&page_path)).unwrap();
        assert!(page.starts_with(b"%PDF-"));
        let page_hash = fs_blake3::hash_bytes(&page);
        let pack = compile(slug);
        let rows = ledger.for_artifact(&ArtifactRef::Matdb(slug.into()));
        assert_eq!(rows.len(), 1);
        let row = rows[0];
        assert_eq!(row.uncertainty, UncertaintyStatement::Bounded);
        assert_eq!(row.digitization, Digitization::FigureDigitized);
        assert_eq!(row.decision, Decision::FixtureOnly);
        assert!(
            row.notes
                .contains("experimental scatter and model discrepancy remain unstated")
        );
        assert!(
            row.artifacts
                .contains(&ArtifactRef::Table(format!("{directory}/digitization.tsv")))
        );
        assert!(
            row.artifacts
                .contains(&ArtifactRef::Table(page_path.clone()))
        );

        let mut passes = BTreeMap::<&str, Vec<&str>>::new();
        let mut failed_half_error = false;
        for line in &lines[header + 1..] {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 21);
            let read = Read {
                claim: fields[0],
                pass: fields[1],
                fields,
            };
            assert!(matches!(read.pass, "A" | "B" | "C"));
            passes.entry(read.claim).or_default().push(read.pass);
            let (nominal, radius) = (read.number(19), read.number(20));
            let bounds = read.bounds();
            assert!(
                read.encloses(nominal, radius),
                "{slug}/{} {}: {bounds:?}",
                read.claim,
                read.pass
            );
            // A half-sized claim fails the actual selected pixel/anchor box,
            // proving the stated extraction-error term is load-bearing.
            failed_half_error |= !read.encloses(nominal, radius * 0.5);
            let dpi = read.number(2);
            assert!(matches!(dpi as u32, 150 | 240 | 300));
            if read.fields[3] != "none" {
                let (x0, v0, x1, v1, target) = (
                    read.number(4),
                    read.number(5),
                    read.number(6),
                    read.number(7),
                    read.number(8),
                );
                let fraction = match read.fields[3] {
                    "linear" => (target - v0) / (v1 - v0),
                    "log10" => det::ln(target / v0) / det::ln(v1 / v0),
                    other => panic!("unknown axis {other}"),
                };
                let target_pixel = x0 + fraction * (x1 - x0);
                assert!((read.number(9) - target_pixel).abs() <= 0.002);
            }
            let (scale, dims) = match read.fields[18] {
                "MPa" => (1e6, fs_qty::Pressure::DIMS),
                "%" => (0.01, Dims([0; 6])),
                "1" => (1.0, Dims([0; 6])),
                other => panic!("unsupported figure unit {other}"),
            };
            let name = property(read.claim);
            let claim = pack.claims().claims_for(name)[0].1;
            assert_eq!(claim.key.dims(), dims);
            let (axis, expected) = match read.claim {
                "tensile-strength-0p02in" => ("wire_diameter", 0.02 * 0.0254),
                "tensile-strength-0p18in" => ("wire_diameter", 0.18 * 0.0254),
                "specific-damping-capacity-5ksi" => {
                    ("surface_shear_stress", 5000.0 * 6894.757293168)
                }
                _ => ("stress_amplitude", 10000.0 * 6894.757293168),
            };
            let (lo, hi) = claim.validity.bounds()[axis];
            assert!((lo / expected - 1.0).abs() <= 2e-15 && lo == hi);
            if read.fields[18] == "MPa" {
                let ksi_to_mpa = 4.448_221_615_260_5 / (0.0254 * 0.0254) * 1e-3;
                assert!((read.number(17) / ksi_to_mpa - 1.0).abs() <= 1e-12);
            }
            let (value, _) = query_scope(&pack, name);
            assert!((value / (nominal * scale) - 1.0).abs() <= 2e-15);
            let PropertyValue::Scalar { value: scalar, .. } = &claim.value else {
                panic!("expected pointwise scalar")
            };
            assert_eq!(value.to_bits(), scalar.to_bits());
            assert!(value > radius * scale);
            if read.fields[18] == "%" {
                let eta = value / core::f64::consts::TAU;
                let eta_error = radius * scale / core::f64::consts::TAU;
                assert!((eta * core::f64::consts::TAU / value - 1.0).abs() <= 2e-15);
                assert!(eta > eta_error);
                // A percent/percentage-point value cannot masquerade as SI.
                assert!(!read.encloses(value, radius * scale));
            }
            reading_count += 1;
            println!(
                "{{\"suite\":\"matdb-pack\",\"bead\":\"{BEAD}\",\"test\":\"g2_g3_redigitized_metal_source_bounds_reach_consumers\",\"source\":\"{slug}\",\"claim\":\"{}\",\"pass\":\"{}\",\"dpi\":{dpi},\"extracted_bounds_in_source_unit\":{:?},\"nominal\":{nominal},\"extraction_half_width\":{radius},\"unit\":\"{}\",\"SI_value\":{value},\"SI_extraction_half_width\":{},\"total_uncertainty\":\"unstated\",\"decision\":\"fixture-only\",\"pack\":\"{}\",\"source_page_blake3\":\"{page_hash}\",\"verdict\":\"matched\"}}",
                read.claim,
                read.pass,
                bounds,
                read.fields[18],
                radius * scale,
                pack.content_hash()
            );
        }
        assert!(
            failed_half_error,
            "an undersized extraction bound must fail for {slug}"
        );
        rejection_count += 1;
        for actual in passes.values_mut() {
            actual.sort_unstable();
            assert_eq!(actual.as_slice(), &["A", "B", "C"]);
        }
        let observations =
            fs::read_to_string(workspace_path(&format!("{directory}/properties.tsv"))).unwrap();
        let lower = observations.to_ascii_lowercase();
        assert!(
            lower.contains("experimental scatter and model discrepancy remain unstated")
                || lower.contains("experimental scatter or model discrepancy, both unstated")
        );
        if slug == SLUGS[0] {
            let (density, _) = query_scope(&pack, "density");
            assert!((density - 1000.0 * (7.855 - 0.032 * 0.85)).abs() < 1e-10);
            assert_eq!(
                pack.claims().claims_for("density")[0].1.key.dims(),
                fs_qty::Density::DIMS
            );
            for name in [
                "tensile_strength_diameter_0p01in",
                "tensile_strength_diameter_0p20in",
                "young_modulus",
                "loss_factor",
            ] {
                assert!(
                    pack.claims().claims_for(name).is_empty(),
                    "unsupported C447 claim {name}"
                );
            }
            assert!(observations.contains("density is NOT graph-read"));
            assert!(observations.contains("separately authored bending-loss spectrum"));
        } else if slug == SLUGS[1] {
            assert!(observations.contains("88 Cu, 10 Zn, 2 Sn verbatim"));
            assert!(observations.contains("does not identify bell bronze"));
            assert!(observations.contains("nominal room temperature"));
        } else {
            assert!(observations.contains("UPPER BOUND proxy for material loss"));
            assert!(observations.contains("never a separated constitutive law"));
            assert!(pack.claims().claims_for("loss_factor").is_empty());
            let (_, point) = query_scope(&pack, "specimen_loss_factor_10ksi_50hz");
            let point = point.with("frequency", 100.0).unwrap();
            assert!(matches!(
                pack.claims().query(
                    "specimen_loss_factor_10ksi_50hz",
                    &point,
                    SelectionPolicy::SingleClaimOnly
                ),
                Err(MatDbError::NoClaimInDomain { .. })
            ));
        }
    }
    assert_eq!(reading_count, 15);
    println!(
        "{{\"suite\":\"matdb-pack\",\"bead\":\"{BEAD}\",\"test\":\"g2_g3_redigitized_metal_source_bounds_reach_consumers\",\"reads_checked\":{reading_count},\"undersized_bound_mutations_caught\":{rejection_count},\"verdict\":\"matched\",\"scope\":\"figure extraction and scoped pack consumers; no instrument fidelity or acceptance promotion\",\"elapsed_s\":{}}}",
        started.elapsed().as_secs_f64()
    );
}
