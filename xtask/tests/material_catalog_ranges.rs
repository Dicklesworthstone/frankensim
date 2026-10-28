//! Real CLI discovery over complete source support, including interior refusals.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn run(root: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtask"))
        .arg("matdb-catalog")
        .args(args)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap()
}

fn success(output: Output) -> String {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

fn refusal(output: Output) -> String {
    assert!(!output.status.success());
    assert!(
        output.stdout.is_empty(),
        "refused exact query published a successful-looking prefix"
    );
    String::from_utf8(output.stderr).unwrap()
}

fn fixture_root() -> PathBuf {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let root = std::env::temp_dir().join(format!("matdb-range-cli-{}-{nonce}", std::process::id()));
    std::fs::create_dir(&root).unwrap();
    root
}

fn source(root: &Path, id: &str, rows: &str) {
    let directory = root.join("data/matdb/seed-v1").join(id);
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(directory.join("manifest.tsv"), format!(
        "frankensim.matdb-manifest.v3\npack_id\t{id}\nredistribution\tpermitted\tOriginal synthetic numerical regression\ncitation\tOriginal range-support regression\nlicense\tinternal-use\nsource\tfixture\tproperties.tsv\tmaterial-tsv-v1\nproperty\tresponse\tresponse\tdimensional\naxis\tstate\tstate\tdimensional\naxis\tfixed\tfixed\tdimensional\n"
    )).unwrap();
    std::fs::write(directory.join("properties.tsv"), format!(
        "frankensim.matdb-source.v1\nobservation\tfixture\tSynthetic material\tAnalytical support fixture\tNot measured material data\n{rows}"
    )).unwrap();
}

#[test]
fn g3_ranges_select_continuous_sources_and_refuse_discrete_gapped_or_ambiguous_interiors() {
    let root = fixture_root();
    let linear = "curve\tlaw\tfixture\tresponse\tstate\t1\t1\t0:10,5:50,10:20\tlinear\nuncertainty\tlaw\tunstated\t-\t-\t-\t-\nvalidity\tlaw\tstate\t0\t10\t1\nvalidity\tlaw\tfixed\t0\t2\t1\n";
    source(&root, "continuous", linear);
    source(
        &root,
        "discrete",
        &linear.replace("\tlinear\n", "\ttabulated\n"),
    );
    source(
        &root,
        "interior-conflict",
        &format!(
            "{linear}scalar\tother\tfixture\tresponse\t30\t1\tconstant\nuncertainty\tother\tunstated\t-\t-\t-\t-\nvalidity\tother\tstate\t4\t6\t1\nvalidity\tother\tfixed\t0\t2\t1\n"
        ),
    );
    source(
        &root,
        "gapped",
        "curve\tleft\tfixture\tresponse\tstate\t1\t1\t0:10,4:30\tlinear\nuncertainty\tleft\tunstated\t-\t-\t-\t-\nvalidity\tleft\tstate\t0\t4\t1\nvalidity\tleft\tfixed\t0\t2\t1\ncurve\tright\tfixture\tresponse\tstate\t1\t1\t6:30,10:20\tlinear\nuncertainty\tright\tunstated\t-\t-\t-\t-\nvalidity\tright\tstate\t6\t10\t1\nvalidity\tright\tfixed\t0\t2\t1\n",
    );

    let listing = success(run(
        &root,
        &[
            "--property",
            "response",
            "--range",
            "state=0:10",
            "--at",
            "fixed=1",
        ],
    ));
    let rows: Vec<_> = listing.lines().collect();
    assert_eq!(rows.len(), 2, "{listing}");
    assert!(rows[0].ends_with("lower-corner-si\tupper-corner-si\tquery-lower\tquery-upper"));
    assert!(rows[1].starts_with("continuous\t"));
    assert_eq!(rows[0].split('\t').count(), rows[1].split('\t').count());
    assert!(rows[1].contains("response=10 ["));
    assert!(rows[1].contains("response=20 ["));

    for id in ["discrete", "gapped", "interior-conflict"] {
        // Both endpoints resolve independently. That does not establish support between them.
        for state in ["state=0", "state=10"] {
            success(run(
                &root,
                &[
                    "--pack",
                    id,
                    "--property",
                    "response",
                    "--at",
                    state,
                    "--at",
                    "fixed=1",
                ],
            ));
        }
        refusal(run(
            &root,
            &[
                "--pack",
                id,
                "--property",
                "response",
                "--range",
                "state=0:10",
                "--at",
                "fixed=1",
            ],
        ));
    }
    let corners = success(run(
        &root,
        &[
            "--pack",
            "continuous",
            "--property",
            "response",
            "--range",
            "state=0:10",
            "--range",
            "fixed=0:2",
        ],
    ));
    assert!(corners.contains("corner values are not extrema or uncertainty bounds"));
    assert!(
        corners
            .lines()
            .find(|line| line.starts_with("lower-corner-si\t"))
            .unwrap()
            .contains("response=10 [")
    );
    assert!(
        corners
            .lines()
            .find(|line| line.starts_with("upper-corner-si\t"))
            .unwrap()
            .contains("response=20 [")
    );
    // The interior maximum is 50, deliberately larger than either displayed corner.
    let middle = success(run(
        &root,
        &[
            "--pack",
            "continuous",
            "--property",
            "response",
            "--at",
            "state=5",
            "--at",
            "fixed=1",
        ],
    ));
    assert!(middle.contains("evaluated-si\tresponse=50 ["));
}

#[test]
fn g0_ranges_require_explicit_finite_typed_compatible_context() {
    let root = fixture_root();
    source(
        &root,
        "scalar",
        "scalar\tlaw\tfixture\tresponse\t7\t1\tconstant\nuncertainty\tlaw\tunstated\t-\t-\t-\t-\nvalidity\tlaw\tstate\t0\t10\t1\nvalidity\tlaw\tfixed\t0\t2\t1\n",
    );
    let base = ["--pack", "scalar", "--property", "response"];
    for extra in [
        vec!["--range", "state=10:0", "--at", "fixed=1"],
        vec!["--range", "state=NaN:10", "--at", "fixed=1"],
        vec!["--range", "state=0:inf", "--at", "fixed=1"],
        vec!["--range", "state=0:10:20", "--at", "fixed=1"],
        vec!["--range", "state=0:10"],
        vec!["--range", "state=0:10", "--at", "fixed=3"],
        vec![
            "--range",
            "state=0:10",
            "--at",
            "state=0",
            "--at",
            "fixed=1",
        ],
        vec![
            "--at",
            "state=0",
            "--range",
            "state=0:10",
            "--at",
            "fixed=1",
        ],
        vec![
            "--range",
            "state=0:10",
            "--range",
            "state=0:10",
            "--at",
            "fixed=1",
        ],
        vec![
            "--range",
            "state=0:10",
            "--range",
            "unrecognized=0:1",
            "--at",
            "fixed=1",
        ],
    ] {
        refusal(run(&root, &[base.as_slice(), extra.as_slice()].concat()));
    }
    refusal(run(&root, &["--range", "state=0:10"]));
    let exact = success(run(
        &root,
        &[
            "--pack",
            "scalar",
            "--property",
            "response",
            "--range",
            "state=5:5",
            "--at",
            "fixed=1",
        ],
    ));
    assert!(exact.contains("lower-corner-si\tresponse=7 ["));
    assert!(exact.contains("upper-corner-si\tresponse=7 ["));
}

#[test]
fn g3_real_nickel_transport_range_keeps_both_claims_and_common_temperature_coverage() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let conditions = [
        "--at",
        "engineering-reference-model=1",
        "--at",
        "source-specimens-matched=0",
        "--at",
        "annealed-material=1",
    ];
    let base = [
        "--pack",
        "nickel-200-annealed-special-metals-2006",
        "--property",
        "thermal-conductivity",
        "--property",
        "electrical-resistivity",
    ];
    let range = ["--range", "temperature=373.15:1273.15"];
    let output = success(run(
        root,
        &[base.as_slice(), conditions.as_slice(), range.as_slice()].concat(),
    ));
    let low = output
        .lines()
        .find(|line| line.starts_with("lower-corner-si\t"))
        .unwrap();
    let high = output
        .lines()
        .find(|line| line.starts_with("upper-corner-si\t"))
        .unwrap();
    assert!(low.contains("thermal-conductivity=66.5 ["));
    assert!(high.contains("thermal-conductivity=68.2 ["));
    assert!(low.contains("electrical-resistivity=0.00000013 ["));
    assert!(high.contains("electrical-resistivity=0.00000051 ["));
    for part in low.split("claim=").skip(1) {
        let claim = part.split(';').next().unwrap();
        assert!(high.contains(&format!("claim={claim};")));
    }
    assert!(output.contains("AbsoluteTemperature"));
    let outside = ["--range", "temperature=373.15:1373.15"];
    let error = refusal(run(
        root,
        &[base.as_slice(), conditions.as_slice(), outside.as_slice()].concat(),
    ));
    assert!(error.contains("thermal-conductivity"));
}
