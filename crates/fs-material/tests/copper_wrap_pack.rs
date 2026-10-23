//! Gate for the string-wrap form of NBS C447's copper density (bead `frankensim-irp-v5-instrument-realism-pexfp.363`, U2): the
//! value equals the original copper-pure-nbs-c447 pack's and is admitted at 20 C on a typed temperature axis.

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).parent().and_then(std::path::Path::parent).expect("root").to_path_buf()
}

fn read(dir: &str, file: &str) -> String {
    std::fs::read_to_string(repo_root().join(format!("data/matdb/seed-v1/{dir}/{file}"))).expect(file)
}

fn density(properties: &str) -> f64 {
    properties.lines().map(|l| l.split('\t').collect::<Vec<_>>())
        .find(|c| c[0] == "scalar" && c[3] == "density").expect("density row")[4].parse().expect("value")
}

#[test]
fn copper_wrap_pack_matches_its_source_row() {
    let dir = "copper-pure-nbs-c447-string-wrap";
    let manifest = read(dir, "manifest.tsv");
    assert!(manifest.starts_with("frankensim.matdb-manifest.v3\n"));
    assert!(manifest.contains("\nlicense\tWork-of-the-US-Government-Public-Use-Permitted\n"));
    assert!(manifest.contains("\naxis\tT\tT\tabsolute-temperature\n"));
    let properties = read(dir, "properties.tsv");
    assert!(properties.starts_with("frankensim.matdb-source.v1\n"));
    assert_eq!(density(&properties), density(&read("copper-pure-nbs-c447", "properties.tsv")));
    assert!(properties.contains("\nvalidity\tdensity\tT\t293.15\t293.15\tK\n"));
    assert!(properties.contains("\nuncertainty\tdensity\tunstated\t-\t-\t-\t-\n"));
    assert_eq!(properties.lines().filter(|l| l.starts_with("scalar\t")).count(), 1);
}
