//! Synthetic source-to-runtime checks; no empirical material accuracy claim.

use super::tests::{manifest, write_fixture};
use super::*;
use fs_matdb::{QueryPoint, SelectionPolicy};

const SOURCE: &str = concat!(
    "frankensim.matdb-source.v3\n",
    "observation\ttest\tsynthetic specimen\tauthored two-axis law\tno empirical claim\n",
    "surface2d\tresponse\ttest\tresponse\ttemp\tdegC\tpressure\tkPa\tGPa\t0:100:10,100:100:20,0:200:40\t0:1:2\tlinear\n",
    "uncertainty\tresponse\tunstated\t-\t-\t-\t-\n",
    "validity\tresponse\ttemp\t0\t100\tdegC\n",
    "validity\tresponse\tpressure\t100\t200\tkPa\n",
    "joint\ttest\tcoordinates-value\tresponse:surface-x:0,response:surface-y:0,response:surface-value:0\t4,0,9,0,0,16\t-\t-\n",
);

fn surface_manifest(profile: &str) -> String {
    manifest(profile, true).replacen(MANIFEST_HEADER, SURFACE_MANIFEST_HEADER, 1)
        + "axis\ttemp\tT\tabsolute-temperature\naxis\tpressure\tp\tdimensional\n"
}

fn point(temperature_k: f64, pressure_pa: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature_k,
        )
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
            pressure_pa,
        )
        .unwrap()
}

#[test]
fn g1_surface_source_normalizes_both_axes_value_statistics_and_shared_pack_families() {
    for profile in [MATERIAL_PROFILE, INTERFACE_PROFILE] {
        let source = if profile == INTERFACE_PROFILE {
            SOURCE.to_owned()
                + "surface_a\tfixture-a\tsolid\tsynthetic\t0\tframe-a\n"
                + "surface_b\tfixture-b\tsolid\tsynthetic\t0\tframe-b\n"
                + "context\tfixture-medium\t-\tfixture-environment\tfixture-history\n"
        } else {
            SOURCE.to_owned()
        };
        let (path, _) = write_fixture(&surface_manifest(profile), &source);
        let first = compile_manifest(&path).unwrap();
        let second = compile_manifest(&path).unwrap();
        assert_eq!(first.bytes, second.bytes);
        assert_eq!(first.decisions, second.decisions);
        let pack = if profile == INTERFACE_PROFILE {
            let interface = NormalizedInterfacePack::from_bytes_verified(
                first.pack.content_hash(),
                &first.bytes,
            )
            .unwrap();
            assert_eq!(interface.compiler(), SURFACE_INTERFACE_COMPILER_ID);
            interface.claims_pack().clone()
        } else {
            let pack = NormalizedPack::from_bytes_verified(first.pack.content_hash(), &first.bytes)
                .unwrap();
            assert_eq!(pack.compiler(), SURFACE_COMPILER_ID);
            pack
        };
        assert_eq!(pack.schema_version(), 8);
        let (id, claim) = pack.claims().claims_for("response")[0];
        let PropertyValue::Surface2D {
            axes,
            axis_dims,
            samples,
            triangles,
            dims,
        } = &claim.value
        else {
            panic!("compiler must preserve the declared two-axis law");
        };
        assert_eq!(axes, &["T", "p"]);
        assert_eq!(
            *axis_dims,
            [fs_qty::Temperature::DIMS, fs_qty::Pressure::DIMS]
        );
        assert_eq!(*dims, fs_qty::Pressure::DIMS);
        assert_eq!(
            samples,
            &vec![
                ([273.15, 100_000.0], 10e9),
                ([373.15, 100_000.0], 20e9),
                ([273.15, 200_000.0], 40e9)
            ]
        );
        assert_eq!(triangles, &vec![[0, 1, 2]]);
        let answer = pack
            .claims()
            .query(
                "response",
                &point(298.15, 125_000.0),
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert!((answer.evidence.value.value - 20e9).abs() < 1e-4);
        pack.claims().verify_receipt(&answer.receipt).unwrap();
        assert!(
            pack.claims()
                .query(
                    "response",
                    &point(363.15, 190_000.0),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err(),
            "the unmeasured rectangle corner must remain outside source support"
        );
        for (member, scale, offset) in [
            (StatisticMember::surface_coordinate(id, 0, 0), 1.0, 273.15),
            (StatisticMember::surface_coordinate(id, 0, 1), 1000.0, 0.0),
            (StatisticMember::surface_ordinate(id, 0), 1e9, 0.0),
        ] {
            let receipt = pack
                .normalizations()
                .iter()
                .find(|receipt| receipt.target() == &NormalizationTarget::ClaimValue(member))
                .unwrap();
            assert_eq!(receipt.scale(), scale);
            assert_eq!(receipt.offset(), offset);
        }
        let joint = &pack.joint_statistics()[0];
        for (row, member) in joint.members().iter().copied().enumerate() {
            let expected = if member == StatisticMember::surface_coordinate(id, 0, 0) {
                4.0
            } else if member == StatisticMember::surface_coordinate(id, 0, 1) {
                9e6
            } else {
                assert_eq!(member, StatisticMember::surface_ordinate(id, 0));
                16e18
            };
            assert_eq!(joint.covariance()[packed_index(row, row)], expected);
            for column in 0..row {
                assert_eq!(joint.covariance()[packed_index(row, column)], 0.0);
            }
        }
    }
}

#[test]
fn g3_surface_data_substitution_and_exact_only_support_reach_runtime() {
    let manifest = surface_manifest(MATERIAL_PROFILE);
    let mut values = Vec::new();
    for source in [SOURCE.to_owned(), SOURCE.replace("0:200:40", "0:200:80")] {
        let (path, _) = write_fixture(&manifest, &source);
        let output = compile_manifest(&path).unwrap();
        let pack = NormalizedPack::from_bytes(&output.bytes).unwrap();
        let answer = pack
            .claims()
            .query(
                "response",
                &point(298.15, 125_000.0),
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        values.push(answer.evidence.value.value);
    }
    assert!((values[0] - 20e9).abs() < 1e-4);
    assert!((values[1] - 30e9).abs() < 1e-4);
    let exact = SOURCE.replace("0:1:2\tlinear", "-\ttabulated");
    let (path, _) = write_fixture(&manifest, &exact);
    let output = compile_manifest(&path).unwrap();
    let pack = NormalizedPack::from_bytes(&output.bytes).unwrap();
    assert_eq!(
        pack.claims()
            .query(
                "response",
                &point(273.15, 100_000.0),
                SelectionPolicy::SingleClaimOnly
            )
            .unwrap()
            .evidence
            .value
            .value,
        10e9
    );
    assert!(
        pack.claims()
            .query(
                "response",
                &point(298.15, 125_000.0),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}

#[test]
fn g4_surface_sources_refuse_old_schemas_missing_units_support_and_bad_connectivity() {
    let manifest = surface_manifest(MATERIAL_PROFILE);
    for (manifest, source, expected) in [
        (
            manifest.replace(SURFACE_MANIFEST_HEADER, STRESS_MANIFEST_HEADER),
            SOURCE.to_owned(),
            "unsupported_source_schema",
        ),
        (
            manifest.clone(),
            SOURCE.replace(SURFACE_SOURCE_HEADER, SOURCE_HEADER),
            "unknown_source_record",
        ),
        (
            manifest.replace("axis\tpressure\tp\tdimensional\n", ""),
            SOURCE.to_owned(),
            "missing_surface_axis_declaration",
        ),
        (
            manifest.clone(),
            SOURCE.replace("validity\tresponse\tpressure\t100\t200\tkPa\n", ""),
            "unused_name_mapping",
        ),
        (
            manifest.clone(),
            SOURCE.replace("100\t200\tkPa", "100\t150\tkPa"),
            "claim_admission_refused",
        ),
        (
            manifest.clone(),
            SOURCE.replace("100\t200\tkPa", "100\t200\tK"),
            "surface_validity_dims_mismatch",
        ),
        (
            manifest.clone(),
            SOURCE.replace("0:1:2\tlinear", "0:1:2\tconstant"),
            "incompatible_interpolation",
        ),
        (
            manifest.clone(),
            SOURCE.replace("0:1:2\tlinear", "0:1:99\tlinear"),
            "claim_admission_refused",
        ),
        (
            manifest.clone(),
            SOURCE.replace("0:1:2\tlinear", "00:1:2\tlinear"),
            "noncanonical_unsigned_integer",
        ),
        (
            manifest.clone(),
            SOURCE.replace("0:1:2\tlinear", "0:1:2,0:1:2\tlinear"),
            "claim_admission_refused",
        ),
        (
            manifest.clone(),
            SOURCE.replace("0:100:10,100:100:20,0:200:40", "0:100:10,0:100:20,0:200:40"),
            "claim_admission_refused",
        ),
        (
            manifest.clone(),
            SOURCE.replace("response:surface-y:0", "response:surface-y:9"),
            "unknown_joint_member",
        ),
    ] {
        let (path, _) = write_fixture(&manifest, &source);
        let error = compile_manifest(&path).unwrap_err();
        assert_eq!(error.code, expected, "{error}");
    }
}
