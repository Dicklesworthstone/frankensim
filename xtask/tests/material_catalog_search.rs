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
