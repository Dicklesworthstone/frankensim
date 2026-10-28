//! G0/G3: users can find shipped sources and inspect the real compiled claims.
#![deny(unsafe_code)]

use std::path::Path;
use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    Command::new(env!("CARGO_BIN_EXE_xtask"))
        .arg("matdb-catalog")
        .args(args)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap()
}

fn success(args: &[&str]) -> String {
    let output = run(args);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn catalog_search_and_required_properties_use_compiled_source_contents() {
    let result = success(&[
        "--search",
        "ENSINGER",
        "--property",
        "thermal-conductivity",
        "--property",
        "specific-heat-capacity",
    ]);
    let mut lines = result.lines();
    assert_eq!(
        lines.next(),
        Some("source-id\tpack-id\tfamily\tclaims\tproperties\tmanifest\tstatus\tdiagnostic")
    );
    let rows: Vec<_> = lines.collect();
    assert!(rows.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(
        rows.iter()
            .any(|row| row.starts_with("ensinger-tecamid-66-natural\t"))
    );
    for row in rows {
        let cells: Vec<_> = row.split('\t').collect();
        assert_eq!(cells.len(), 8);
        assert_eq!(cells[6], "compiled");
        let properties: Vec<_> = cells[4].split(',').collect();
        assert!(properties.contains(&"thermal-conductivity"));
        assert!(properties.contains(&"specific-heat-capacity"));
    }
    let melamine = success(&["--search", "melamine"]);
    assert!(melamine.contains("basotect-gplus-melamine-foam\t"));
    assert!(melamine.contains("kronospan-melamine-faced-p2-18mm-2024\t"));
}

#[test]
fn inspection_keeps_directional_curves_surface_topology_and_source_caveats() {
    let g10 = success(&["--pack", "g10cr-nist-cryogenic"]);
    assert!(g10.starts_with("source\tg10cr-nist-cryogenic\npack\tg10cr-nist-cryogenic\n"));
    for record in [
        "\ncontent-hash\t",
        "\ncitation\t",
        "\nlicense\t",
        "\nclaim\t",
        "\nknot\t",
        "\naxis\t",
        "\nobservation\t",
        "\nclaim-observation\t",
    ] {
        assert!(g10.contains(record), "missing {record}");
    }
    assert!(g10.contains("conductivity-direction"));
    assert!(g10.contains("source-pressure-known"));
    let eps = success(&["--pack", "eps-nist-srm1453-2013"]);
    assert!(eps.contains("\nsample\t"));
    assert!(eps.contains("\ntriangle\t"));
    // Inspection is deterministic and does not publish a generated binary pack.
    assert_eq!(g10, success(&["--pack", "g10cr-nist-cryogenic"]));
    // Independently sourced argon bundles legitimately share the internal ID Ar.
    let argon = success(&["--search", "argon"]);
    assert!(argon.contains("argon\tAr\t"));
    assert!(argon.contains("argon-nasa9-mcbride2002\tAr\t"));
    let species = success(&["--pack", "argon"]);
    assert!(species.contains("inspection-scope\tmetadata-only;"));
}

#[test]
fn invalid_options_missing_ids_and_conflicting_filters_refuse_without_partial_output() {
    for args in [
        vec!["--search"],
        vec!["--property", "--pack", "g10cr-nist-cryogenic"],
        vec!["--pack", "missing-pack"],
        vec![
            "--pack",
            "g10cr-nist-cryogenic",
            "--property",
            "nonexistent-property",
        ],
        vec!["--pack", "g10cr-nist-cryogenic", "--pack", "other"],
    ] {
        let output = run(&args);
        assert!(!output.status.success(), "{args:?}");
        assert!(output.stdout.is_empty(), "no partial catalog for {args:?}");
        assert!(!output.stderr.is_empty());
    }
    assert!(success(&["--help"]).contains("common valid state"));
}

#[test]
fn state_search_and_exact_inspection_sample_real_source_curves() {
    let state = [
        "--property",
        "density",
        "--property",
        "thermal-conductivity",
        "--at",
        "temperature=350",
        "--at",
        "source-generic-316-identity=1",
        "--at",
        "source-uns-designation-known=0",
        "--at",
        "source-heat-treatment-known=0",
        "--at",
        "source-product-form-known=0",
        "--at",
        "source-pressure-known=0",
    ];
    let mut args = vec!["--pack", "stainless-316-nasa-tp216435"];
    args.extend_from_slice(&state);
    let exact = success(&args);
    let values = exact
        .lines()
        .find(|line| line.starts_with("evaluated-si\t"))
        .unwrap();
    assert!(values.contains("density=7938.6"), "{values}");
    assert!(values.contains("thermal-conductivity=14.3"), "{values}");
    assert!(values.contains("Estimate; claim="));
    assert!(values.contains("LinearInside"));
    assert!(values.contains("uncertainty=Unstated"));
    let point = exact
        .lines()
        .find(|line| line.starts_with("query-point\t"))
        .unwrap();
    assert!(point.contains("AbsoluteTemperature"));
    assert_eq!(exact, success(&args));

    let mut search = vec!["--search", "stainless-316-nasa-tp216435"];
    search.extend_from_slice(&state);
    let filtered = success(&search);
    assert_eq!(filtered.lines().count(), 2);
    assert!(
        filtered
            .lines()
            .next()
            .unwrap()
            .ends_with("\tevaluated-si\tquery-point")
    );
    assert_eq!(filtered.lines().nth(1).unwrap().split('\t').count(), 10);
    assert!(filtered.contains("density=7938.6"));
    let outside = search
        .iter()
        .map(|value| {
            if *value == "temperature=350" {
                "temperature=500"
            } else {
                *value
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(success(&outside).lines().count(), 1);
}

#[test]
fn state_query_never_infers_missing_conditions_or_ignores_extra_axes() {
    for args in [
        vec![
            "--pack",
            "niobium-espi-reference",
            "--property",
            "density",
            "--at",
            "temperature=293.15",
        ],
        vec![
            "--pack",
            "niobium-espi-reference",
            "--property",
            "density",
            "--at",
            "temperature=293.15",
            "--at",
            "engineering-reference-model=1",
            "--at",
            "source-temperatures-matched=0",
            "--at",
            "pressure=101325",
        ],
        vec![
            "--pack",
            "iridium-espi-reference",
            "--property",
            "density",
            "--property",
            "specific-heat-capacity",
            "--at",
            "temperature=293.15",
            "--at",
            "source-reference-espi-iridium=1",
            "--at",
            "source-temperature-known=0",
            "--at",
            "source-pressure-known=0",
        ],
        vec![
            "--pack",
            "iridium-espi-reference",
            "--property",
            "electrical-resistivity",
            "--at",
            "temperature=283.15",
            "--at",
            "source-reference-espi-iridium=1",
            "--at",
            "source-pressure-known=0",
        ],
    ] {
        let output = run(&args);
        assert!(!output.status.success(), "must refuse {args:?}");
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).contains("inspect --pack"));
    }
}
