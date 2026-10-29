//! G0/G1/G3: published varnish TSV -> source card -> passive layered mechanics.
//!
//! These are numerical transport and constitutive checks. A uniaxial modern
//! film fit does not establish isotropy or historical-instrument acoustics.
#![deny(unsafe_code)]

use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, QueryPoint};
use fs_material::{
    coating::ResolvedCoating,
    elastic::OrthotropicElastic,
    laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, resolve_material_state_point,
    },
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use std::{f64::consts::TAU, path::Path, process::Command};

fn compile(slug: &str) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-coating-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let out = dir.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&out)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.pack_id(), slug);
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "source-identified coating; no inferred constitutive symmetry".into(),
            process: "original specimen and conditioning retained in every claim".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
}

fn flag(point: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    point
        .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
        .unwrap()
}

fn point(hz: f64) -> QueryPoint {
    let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let p = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "frequency",
            semantic(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
            hz,
        )
        .unwrap();
    flag(
        flag(
            flag(p, "conditioning_relative_humidity", 0.35),
            "source_fitted_response",
            1.,
        ),
        "source_density_sample_matched",
        0.,
    )
}

fn source(pack: &NormalizedMaterialCardPack) -> ResolvedMaterialStatePoint {
    let state = resolve_material_state_point(
        pack.card(),
        &point(1000.),
        &ResolvedCoating::requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(state.card_identity(), pack.card().content_hash());
    assert_eq!(state.material(), pack.card().id());
    for property in state.properties() {
        pack.card()
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    state
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "actual {actual:e}, expected {expected:e}, relative tolerance {relative:e}"
    );
}

#[test]
fn g0_all_six_source_cards_reach_a_frequency_slice_without_losing_receipts() {
    // Independent complex-arithmetic evaluations of the published Table 2 law
    // at 1 kHz, not additional measured observations. Density remains a group
    // mean; nu and the small-strain envelope below are authored assumptions.
    for (family, sample, storage, loss, density) in [
        ("oil-copal", 1, 1.55486870068e9, 1.44058950841e8, 1097.),
        ("oil-copal", 2, 1.56880652665e9, 1.58642242845e8, 1097.),
        ("spirit-shellac", 1, 1.88845938660e9, 1.41656174937e8, 984.),
        ("spirit-shellac", 2, 1.93275048767e9, 1.53507086324e8, 984.),
        ("mineral-ground", 1, 2.16288588059e9, 1.32446992687e8, 1155.),
        ("mineral-ground", 2, 2.26157869020e9, 1.68120771267e8, 1155.),
    ] {
        let pack = compile(&format!("varnish-{family}-lammlein2021-s{sample}"));
        let state = source(&pack);
        let nu = 0.3;
        let film = ResolvedCoating::from_material(state.clone(), nu, 0.0003).unwrap();
        let thickness = 20e-6; // Authored external coat geometry, not a DMA film measurement.
        let omega = TAU * 1000.;
        let slice = film.at_frequency(thickness, omega).unwrap();
        assert_eq!(film.source(), &state);
        assert_eq!(slice.layer.material.as_ref(), Some(&state));
        assert_eq!(film.poisson_ratio(), nu);
        assert_eq!(film.strain_limit(), 0.0003);
        assert_eq!(slice.layer.density_kg_m3, density);
        assert_eq!(slice.reference_omega_rad_s, omega);
        let LayerLaw::Elastic(elastic) = &slice.layer.law else {
            panic!("a frozen frequency slice must be declared elastic plus viscosity");
        };
        close(elastic.e[0], storage, 1e-10);
        close(elastic.g[0], storage / (2. * (1. + nu)), 1e-10);
        close(
            slice.plane_viscosity_pa_s[0] * omega,
            loss / (1. - nu * nu),
            1e-10,
        );
        close(
            slice.plane_viscosity_pa_s[1] * omega,
            nu * loss / (1. - nu * nu),
            1e-10,
        );
        close(
            slice.transverse_shear_viscosity_pa_s[0] * omega,
            loss / (2. * (1. + nu)),
            1e-10,
        );
        let section =
            LaminateSection::try_new(vec![slice.layer], ReferenceSurface::Midsurface).unwrap();
        close(section.i0_kg_m2(), density * thickness, 1e-14);
        close(
            section.a_pa_m()[0],
            storage * thickness / (1. - nu * nu),
            1e-10,
        );
        let (a_loss, _, d_loss) = section
            .integrate_layer_matrices(&[slice.plane_viscosity_pa_s])
            .unwrap();
        close(a_loss[0] * omega, loss * thickness / (1. - nu * nu), 1e-10);
        close(
            d_loss[0] * omega,
            loss * thickness.powi(3) / (12. * (1. - nu * nu)),
            1e-10,
        );
    }
}

#[test]
fn g1_source_film_prony_retains_offset_mass_and_complex_bending_response() {
    let pack = compile("varnish-oil-copal-lammlein2021-s1");
    let state = source(&pack);
    let nu = 0.3;
    let film = ResolvedCoating::from_material(state.clone(), nu, 0.0003).unwrap();
    // The error gates compare two executable laws; they are not experimental
    // confidence bounds on the source fit or on a historical instrument.
    let (tensor, lowered) = film.lower_to_tensor_prony(16, 0.001, 0.001).unwrap();
    assert!(lowered.sup_rel_storage_err <= 0.001);
    assert!(lowered.sup_rel_loss_err <= 0.001);
    assert_eq!(lowered.verification_points, 512);
    assert!(
        lowered
            .model
            .terms
            .iter()
            .all(|&(e, tau)| e > 0. && tau > 0.)
    );

    // A synthetic support isolates the coating contribution; these constants
    // are not represented as a source-backed violin plate or measured wood.
    let support_thickness = 0.003;
    let support = LaminaLayer {
        thickness_m: support_thickness,
        orientation_rad: 0.,
        density_kg_m3: 400.,
        material: None,
        law: LayerLaw::Elastic(
            OrthotropicElastic::new([12e9; 3], [0.2; 3], [5e9; 3], 0.01).unwrap(),
        ),
    };
    let reference = ReferenceSurface::Offset(support_thickness / 2.);
    let baseline = LaminateSection::try_new(vec![support.clone()], reference).unwrap();
    let thickness = 20e-6;
    let section = LaminateSection::try_new(
        vec![
            support,
            LaminaLayer {
                thickness_m: thickness,
                orientation_rad: 0.,
                density_kg_m3: film.density_kg_m3(),
                material: Some(state.clone()),
                law: LayerLaw::Prony(tensor),
            },
        ],
        reference,
    )
    .unwrap();
    assert_eq!(section.layers()[1].material.as_ref(), Some(&state));
    let z0 = support_thickness / 2.;
    let z1 = z0 + thickness;
    let first_moment = thickness * (z0 + z1) / 2.;
    let second_moment = thickness * (z0 * z0 + z0 * z1 + z1 * z1) / 3.;
    close(
        section.i0_kg_m2() - baseline.i0_kg_m2(),
        film.density_kg_m3() * thickness,
        1e-12,
    );
    close(
        section.i1_kg_m(),
        film.density_kg_m3() * first_moment,
        1e-12,
    );
    close(
        section.i2_kg() - baseline.i2_kg(),
        film.density_kg_m3() * second_moment,
        1e-12,
    );
    for hz in [20., 100., 1000., 10000., 20000.] {
        let omega = TAU * hz;
        let (storage, loss) = film.law().modulus(omega).unwrap();
        let (a, b, d) = section.complex_abd(omega).unwrap();
        let (a_base, _, d_base) = baseline.complex_abd(omega).unwrap();
        let q_storage = storage / (1. - nu * nu);
        let q_loss = loss / (1. - nu * nu);
        close(a[0].re - a_base[0].re, q_storage * thickness, 0.001);
        close(a[0].im, q_loss * thickness, 0.001);
        close(b[0].re, q_storage * first_moment, 0.001);
        close(b[0].im, q_loss * first_moment, 0.001);
        close(d[0].re - d_base[0].re, q_storage * second_moment, 0.001);
        close(d[0].im, q_loss * second_moment, 0.001);
        assert!(d[0].im > 0. && a[8].im > 0.);
    }
    for omega in [TAU * 19.99, TAU * 20000.01, f64::NAN] {
        assert!(section.complex_abd(omega).is_err());
    }
}

#[test]
fn g3_no_source_extrapolation_or_unstated_tensor_and_historical_substitution() {
    let pack = compile("varnish-oil-copal-lammlein2021-s1");
    let state = source(&pack);
    for nu in [-1., 0.5, f64::NAN, f64::INFINITY] {
        assert!(ResolvedCoating::from_material(state.clone(), nu, 0.0003).is_err());
    }
    for strain in [0., -0.001, f64::NAN] {
        assert!(ResolvedCoating::from_material(state.clone(), 0.3, strain).is_err());
    }
    let film = ResolvedCoating::from_material(state.clone(), 0.3, 0.0003).unwrap();
    for hz in [0., 19.99, 20000.01, f64::NAN] {
        assert!(film.at_frequency(20e-6, TAU * hz).is_err());
    }
    for thickness in [0., -20e-6, f64::NAN] {
        assert!(film.at_frequency(thickness, TAU * 1000.).is_err());
    }
    for p in [
        point(19.99),
        point(20000.01),
        flag(point(1000.), "conditioning_relative_humidity", 0.65),
        flag(point(1000.), "source_density_sample_matched", 1.),
    ] {
        assert!(state.at_query_point(pack.card(), &p).is_err());
    }
    let other = compile("varnish-oil-copal-lammlein2021-s2");
    assert!(state.at_query_point(other.card(), &point(1000.)).is_err());
    let historical = compile("varnish-stradivari-echard2010-stratigraphy");
    let p = flag(
        flag(QueryPoint::new(), "specimen_code", 2.),
        "source_stratigraphy_only",
        1.,
    );
    assert!(
        resolve_material_state_point(
            historical.card(),
            &p,
            &ResolvedCoating::requirements().unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    assert_eq!(film.source(), &state);
}
