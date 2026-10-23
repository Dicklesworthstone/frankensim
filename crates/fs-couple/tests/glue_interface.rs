//! Integration tests for `InterfaceSystemCard` -> `GlueLine` resolution and single-owner
//! modal glue damping projection (`M1.2-W2b`, `frankensim-irp-v5-instrument-realism-pexfp.51`).

use fs_couple::plate_viscosity::{
    GLUE_NORMAL_STIFFNESS_PROPERTY, GLUE_NORMAL_VISCOSITY_PROPERTY,
    GLUE_SHEAR_STIFFNESS_1_PROPERTY, GLUE_SHEAR_STIFFNESS_2_PROPERTY,
    GLUE_SHEAR_VISCOSITY_1_PROPERTY, GLUE_SHEAR_VISCOSITY_2_PROPERTY, ModalPlateViscosity,
    ResolvedGlueLine,
};
use fs_couple::shear_shell::connections::{stiffness_per_area_dims, viscosity_per_area_dims};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_modal::{SliceOptions, slice_window};
use fs_plate::PlateSection;
use fs_plate::shell::junction::{
    EdgeTie, RotationCoupling, ShellComponent, TieBudget, TieKind, TieOffset, TieSide, tie_shells,
};
use fs_plate::shell::{ShellMesh, ShellSupport, assemble_shell};
use fs_qty::Pressure;

fn authored_spruce_maple_states() -> (MaterialStateId, MaterialStateId) {
    (
        MaterialStateId {
            chemistry: "Picea-abies-spruce".to_string(),
            phase: "cellular-orthotropic-wood".to_string(),
            process: "quarter-sawn-air-dried-12pct-mc".to_string(),
            revision: 1,
        },
        MaterialStateId {
            chemistry: "Acer-pseudoplatanus-maple".to_string(),
            phase: "cellular-orthotropic-wood".to_string(),
            process: "quarter-sawn-rib-bent-12pct-mc".to_string(),
            revision: 1,
        },
    )
}

fn build_authored_glue_fixture_card(
    omit_property: Option<&str>,
    wrong_dims_property: Option<&str>,
    third_body: Option<String>,
) -> InterfaceSystemCard {
    let (spruce, maple) = authored_spruce_maple_states();
    let k_dims = stiffness_per_area_dims();
    let c_dims = viscosity_per_area_dims();
    let mut claims = ClaimSet::new();

    for (key, value, dims) in [
        (GLUE_NORMAL_STIFFNESS_PROPERTY, 2.4e11, k_dims),
        (GLUE_SHEAR_STIFFNESS_1_PROPERTY, 1.1e11, k_dims),
        (GLUE_SHEAR_STIFFNESS_2_PROPERTY, 1.2e11, k_dims),
        (GLUE_NORMAL_VISCOSITY_PROPERTY, 4.8e5, c_dims),
        (GLUE_SHEAR_VISCOSITY_1_PROPERTY, 2.2e5, c_dims),
        (GLUE_SHEAR_VISCOSITY_2_PROPERTY, 2.4e5, c_dims),
    ] {
        if omit_property == Some(key) {
            continue;
        }
        let actual_dims = if wrong_dims_property == Some(key) {
            Pressure::DIMS
        } else {
            dims
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, actual_dims),
                value: PropertyValue::Scalar {
                    value,
                    dims: actual_dims,
                },
                validity: ValidityDomain::unconstrained().with("temperature", 270.0, 330.0),
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.1,
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "authored-fixture:hide-glue-spruce-maple-v1".to_string(),
                    license: "CC0-1.0".to_string(),
                    artifact: None,
                },
            })
            .unwrap();
    }

    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: spruce,
            texture_frame: "planed-toothed-0.1mm-L".to_string(),
        },
        SurfaceSpec {
            material: maple,
            texture_frame: "scraped-rib-edge-L".to_string(),
        },
        SystemContext {
            medium: "cured-bondline-0.08mm".to_string(),
            third_body,
            environment: "air-293K-50pctRH".to_string(),
            history: "cured-72h-293K-virgin".to_string(),
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

#[test]
fn g0_glue_compliance_is_data() {
    let (spruce, maple) = authored_spruce_maple_states();
    let card = build_authored_glue_fixture_card(
        None,
        None,
        Some("authored-fixture:hot-hide-glue-192g-bloom".to_string()),
    );
    let point = QueryPoint::new().with("temperature", 293.15).unwrap();

    let resolved = ResolvedGlueLine::from_card(
        &card,
        (&spruce, &maple),
        &point,
        SelectionPolicy::SingleClaimOnly,
        0.003,
        true,
    )
    .expect("authored fixture card resolves with receipts");

    assert_eq!(resolved.card_identity, card.content_hash());
    assert!(resolved.adhesive_third_body.starts_with("authored-fixture:"));
    assert_eq!(resolved.history, "cured-72h-293K-virgin");
    assert_eq!(resolved.glue_line.normal_n_m3, 2.4e11);
    assert_eq!(resolved.glue_line.shear_n_m3, [1.1e11, 1.2e11]);
    assert_eq!(
        resolved.glue_line.viscosity_n_s_m3,
        Some([4.8e5, 2.2e5, 2.4e5])
    );
    assert!(resolved.viscosity_receipts.is_some());

    // Reversed surface order (maple, spruce) must refuse.
    assert!(
        ResolvedGlueLine::from_card(
            &card,
            (&maple, &spruce),
            &point,
            SelectionPolicy::SingleClaimOnly,
            0.003,
            true,
        )
        .is_err()
    );

    // Card without SystemContext.third_body adhesive must refuse.
    let no_adhesive_card = build_authored_glue_fixture_card(None, None, None);
    let err_no_adhesive = ResolvedGlueLine::from_card(
        &no_adhesive_card,
        (&spruce, &maple),
        &point,
        SelectionPolicy::SingleClaimOnly,
        0.003,
        true,
    )
    .unwrap_err();
    assert!(err_no_adhesive.to_string().contains("third_body"));
}

#[test]
fn g0_glue_card_missing_property_refuses_by_name() {
    let (spruce, maple) = authored_spruce_maple_states();
    let point = QueryPoint::new().with("temperature", 293.15).unwrap();
    for missing in [
        GLUE_NORMAL_STIFFNESS_PROPERTY,
        GLUE_SHEAR_STIFFNESS_1_PROPERTY,
        GLUE_SHEAR_STIFFNESS_2_PROPERTY,
        GLUE_NORMAL_VISCOSITY_PROPERTY,
        GLUE_SHEAR_VISCOSITY_1_PROPERTY,
        GLUE_SHEAR_VISCOSITY_2_PROPERTY,
    ] {
        let card = build_authored_glue_fixture_card(
            Some(missing),
            None,
            Some("authored-fixture:hot-hide-glue-192g-bloom".to_string()),
        );
        let err = ResolvedGlueLine::from_card(
            &card,
            (&spruce, &maple),
            &point,
            SelectionPolicy::SingleClaimOnly,
            0.003,
            true,
        )
        .unwrap_err();
        assert!(
            err.to_string().contains(missing),
            "missing property {missing} must be named in refusal: {err}"
        );
    }
}

#[test]
fn g0_glue_card_units_checked() {
    let (spruce, maple) = authored_spruce_maple_states();
    let point = QueryPoint::new().with("temperature", 293.15).unwrap();
    for wrong_key in [
        GLUE_NORMAL_STIFFNESS_PROPERTY,
        GLUE_SHEAR_STIFFNESS_2_PROPERTY,
        GLUE_SHEAR_VISCOSITY_1_PROPERTY,
    ] {
        let card = build_authored_glue_fixture_card(
            None,
            Some(wrong_key),
            Some("authored-fixture:hot-hide-glue-192g-bloom".to_string()),
        );
        let err = ResolvedGlueLine::from_card(
            &card,
            (&spruce, &maple),
            &point,
            SelectionPolicy::SingleClaimOnly,
            0.003,
            true,
        )
        .unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("unit mismatch") && msg.contains(wrong_key),
            "wrong dimensions for {wrong_key} must refuse with unit mismatch: {msg}"
        );
    }
}

#[test]
fn g0_glue_loss_nonnegative_and_counted_once() {
    let (spruce, maple) = authored_spruce_maple_states();
    let card = build_authored_glue_fixture_card(
        None,
        None,
        Some("authored-fixture:hot-hide-glue-192g-bloom".to_string()),
    );
    let point = QueryPoint::new().with("temperature", 293.15).unwrap();
    let resolved = ResolvedGlueLine::from_card(
        &card,
        (&spruce, &maple),
        &point,
        SelectionPolicy::SingleClaimOnly,
        0.003,
        true,
    )
    .unwrap();

    // Two rectangular panels joined along x = 0.2 by the resolved compliant glue line.
    let make_rect = |x0: f64, x1: f64| -> ShellMesh {
        let mut nodes = Vec::new();
        for j in 0..=2 {
            let y = 0.1 * (j as f64);
            for i in 0..=2 {
                let x = x0 + 0.5 * (x1 - x0) * (i as f64);
                nodes.push([x, y, 0.0]);
            }
        }
        let idx = |i: usize, j: usize| j * 3 + i;
        let mut tris = Vec::new();
        for j in 0..2 {
            for i in 0..2 {
                tris.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
                tris.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
            }
        }
        ShellMesh::new(nodes, tris).unwrap()
    };

    let m0 = make_rect(0.0, 0.2);
    let m1 = make_rect(0.2, 0.4);
    let sec = PlateSection::isotropic(10.0e9, 0.3, 0.003, 500.0).unwrap();
    let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
    let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
    let h0 = vec![0.003; m0.nodes.len()];
    let h1 = vec![0.003; m1.nodes.len()];
    let comps = [
        ShellComponent {
            id: 0,
            mesh: &m0,
            model: &mod0,
            nodal_thickness_m: &h0,
        },
        ShellComponent {
            id: 1,
            mesh: &m1,
            model: &mod1,
            nodal_thickness_m: &h1,
        },
    ];
    let ties = [EdgeTie {
        a: TieSide {
            component: 0,
            nodes: vec![2, 5, 8],
        },
        b: TieSide {
            component: 1,
            nodes: vec![0, 3, 6],
        },
        offset: TieOffset::Midsurface,
        rotations: RotationCoupling::Full,
        kind: TieKind::Compliant(resolved.glue_line),
    }];
    let budget = TieBudget {
        max_segments: 32,
        max_slave_dofs: 64,
        gap_tolerance_m: 1e-9,
        angle_tolerance_rad: 1e-6,
    };
    let tied = tie_shells(&comps, &ties, budget, || false).unwrap();

    // Solve first 12 modes (6 rigid + 6 flexible) with a positive shift.
    let sigma = 1.0e5_f64;
    let mut k_shift = fs_sparse::Coo::new(tied.model.free, tied.model.free);
    for i in 0..tied.model.free {
        let (k_cols, k_vals) = tied.model.k.row(i);
        for (&c, &v) in k_cols.iter().zip(k_vals.iter()) {
            k_shift.push(i, c, v);
        }
        let (m_cols, m_vals) = tied.model.m.row(i);
        for (&c, &v) in m_cols.iter().zip(m_vals.iter()) {
            k_shift.push(i, c, sigma * v);
        }
    }
    let k_shifted = k_shift.assemble();
    let report = slice_window(
        &k_shifted,
        &tied.model.m,
        (0.1 * sigma, 2.0e8 + sigma),
        &SliceOptions {
            ritz_tol: 1e-13,
            max_lanczos: 1024,
            max_restarts: 64,
            mass_diagonal_equilibration: true,
            ..SliceOptions::default()
        },
    )
    .unwrap();
    let phis: Vec<Vec<f64>> = report.modes[..12].iter().map(|m| m.phi.clone()).collect();

    // Project through the one modal viscous owner (`fs_couple::plate_viscosity::ModalPlateViscosity`).
    let gate = CancelGate::new();
    let modal_visc =
        ModalPlateViscosity::from_tied_shell_glue(&tied, &phis, 100.0, &gate).expect("PSD glue projection");

    let c_mat = modal_visc.matrix_per_s();
    for m in 0..6 {
        assert!(
            c_mat[m * 12 + m].abs() <= 1e-9,
            "rigid mode {m} must have zero glue loss, got {:.3e}",
            c_mat[m * 12 + m]
        );
    }
    let mut active_flexible = 0usize;
    for m in 6..12 {
        let diag = c_mat[m * 12 + m];
        assert!(diag >= 0.0, "per-mode glue loss must be >= 0");
        if diag > 1e-6 {
            active_flexible += 1;
        }
        let direct = tied.glue_dissipation_rate(&phis[m]).unwrap();
        assert!((diag - direct).abs() <= 1e-11 * direct.max(1.0));
    }
    assert!(active_flexible > 0, "flexible modes across glue seam must experience positive glue damping");
    assert!(modal_visc.channel_count() > 0 && modal_visc.channel_count() <= 6);
}
