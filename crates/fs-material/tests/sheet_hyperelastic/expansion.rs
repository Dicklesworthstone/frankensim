//! Source preparation changes natural strain, never the original sheet inventory.
#[path = "expansion/thermoelastic.rs"]
mod thermoelastic;

use super::*;
use fs_blake3::ContentHash;
use fs_matdb::{ClaimSelection, SelectionPolicy};
use fs_material::sheet::hyperelastic::{
    HyperelasticMembrane,
    expansion::{SheetExpansionFrame, SheetExpansionReference},
};
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
};
use fs_material::visco_hyper::membrane::{FrozenViscoMembrane, ViscoMembraneState};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use std::sync::Arc;

const PRINCIPAL: ContentHash = ContentHash([7; 32]);
const SHEET: ContentHash = ContentHash([8; 32]);
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const NAMES: [&str; 3] = [
    "warp_length_change",
    "weft_length_change",
    "normal_length_change",
];

fn thermal_point(axis: &str, t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(axis, QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with("p", 1e5)
        .unwrap()
}

fn thermal_card(axis: &str, model_span: [f64; 2]) -> MaterialCard {
    let original = source(
        vec![model("neo-hookean-log"), model("mooney-rivlin-decoupled")],
        0.0013,
        0.24,
        ValidityDomain::unconstrained(),
        &[("specific_heat_capacity", 1400.0, Dims([2, 0, -2, -1, 0, 0]))],
    );
    let mut claims = original.claims().clone();
    for (name, values) in
        NAMES
            .into_iter()
            .zip([[0.0, 0.02, 0.06], [-0.02, -0.01, -0.03], [0.02, 0.03, 0.08]])
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [290.0, 300.0, 310.0].into_iter().zip(values).collect(),
                    dims: Dims::NONE,
                },
                interpolation: InterpolationPolicy::LinearInside,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        axis,
                        QuantitySpec::dimensional(Temperature::DIMS),
                        290.0,
                        310.0,
                    )
                    .with("p", 1e5, 2e5),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let mut models = original.models().to_vec();
    for model in &mut models {
        model.validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            model_span[0],
            model_span[1],
        );
    }
    MaterialCard::assemble(original.id().clone(), claims, models).unwrap()
}

fn path_at(
    card: &MaterialCard,
    axis: &str,
    reference: QueryPoint,
    isotropic: bool,
) -> PreparedFreeExpansion {
    let keys = NAMES.map(|name| PropertyKey::new(name, Dims::NONE));
    let pick = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let strains = if isotropic {
        FreeExpansionStrainSources::Isotropic {
            strain: &keys[0],
            selection: pick,
        }
    } else {
        FreeExpansionStrainSources::FixedPrincipal {
            frame: PRINCIPAL,
            strains: [&keys[0], &keys[1], &keys[2]],
            selections: [pick; 3],
        }
    };
    let lo = reference
        .clone()
        .with_quantity(axis, QuantitySpec::dimensional(Temperature::DIMS), 290.0)
        .unwrap();
    let hi = reference
        .clone()
        .with_quantity(axis, QuantitySpec::dimensional(Temperature::DIMS), 310.0)
        .unwrap();
    PreparedFreeExpansion::try_new(
        card,
        FreeExpansionGeometryQuery {
            reference: &reference,
            strains,
            temperature_axis: axis,
            lower: &lo,
            upper: &hi,
        },
    )
    .unwrap()
}

fn path(card: &MaterialCard, axis: &str, isotropic: bool) -> PreparedFreeExpansion {
    path_at(card, axis, thermal_point(axis, 300.0), isotropic)
}

fn sheet(card: &MaterialCard, exact: bool) -> ResolvedSheet {
    let requirements = [
        (AREAL_DENSITY, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, Dims([1, 0, 0, 0, 0, 0])),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    let material = resolve_material_state_point(
        card,
        &thermal_point("T", 300.0),
        &requirements,
        Pick::SingleClaimOnly,
    )
    .unwrap();
    let sheet = ResolvedSheet::try_from_material(material).unwrap();
    if exact {
        sheet
            .with_incompressible_hyperelastic_membrane(card, 1, controls())
            .unwrap()
    } else {
        sheet
            .with_hyperelastic_membrane(card, 1, controls())
            .unwrap()
    }
}

fn frame(q: [[f64; 3]; 3]) -> SheetExpansionFrame {
    SheetExpansionFrame::FixedPrincipal {
        principal_frame: PRINCIPAL,
        sheet_frame: SHEET,
        principal_to_sheet: q,
    }
}

fn bound(sheet: &ResolvedSheet) -> Arc<HyperelasticMembrane> {
    let MembraneLaw::Hyperelastic(law) = sheet.membrane_law().unwrap() else {
        panic!("native hyperelastic source required")
    };
    law
}

fn endpoint(law: &HyperelasticMembrane) -> &SheetExpansionReference {
    law.free_expansion_reference().unwrap()
}

#[test]
fn g3_temperature_revisits_keep_original_datum_sources_model_and_inventory() {
    for axis in ["T", "sample_temperature"] {
        let card = thermal_card(axis, [290.0, 310.0]);
        let original = sheet(&card, false);
        let source = path(&card, axis, false);
        let initial = original
            .clone()
            .with_hyperelastic_free_expansion(&card, source.clone(), 300.0, frame(IDENTITY))
            .unwrap();
        assert_eq!(bound(&initial).reference_stretches(), [1.0; 3]);
        let hot = initial.at_hyperelastic_temperature(&card, 308.0).unwrap();
        let revisited = hot.at_hyperelastic_temperature(&card, 304.0).unwrap();
        let direct = initial.at_hyperelastic_temperature(&card, 304.0).unwrap();
        assert_eq!(revisited, direct);
        assert_eq!(
            hot.at_hyperelastic_temperature(&card, 300.0).unwrap(),
            initial
        );
        let law = bound(&direct);
        for (actual, expected) in
            law.reference_stretches()
                .into_iter()
                .zip([1.036 / 1.02, 0.982 / 0.99, 1.05 / 1.03])
        {
            close(actual, expected, 1e-14, 0.0);
        }
        assert_eq!(direct.material(), original.material());
        assert_eq!(law.material(), original.material());
        assert_eq!(direct.areal_density_kg_m2(), original.areal_density_kg_m2());
        assert_eq!(direct.thickness_m(), original.thickness_m());
        assert!(law.material().property("specific_heat_capacity").is_some());
        assert_eq!(law.model().model_card(), &card.models()[1]);
        assert_eq!(law.model().context(), &thermal_point("T", 304.0));
        assert_eq!(
            endpoint(&law).reference_context(),
            &original.material().query()
        );
        assert_eq!(endpoint(&law).current_context(), law.model().context());
        assert_eq!(endpoint(&law).source(), &source);
        assert_eq!(endpoint(&law).geometry().path_identity(), source.identity());
        assert_eq!(endpoint(&law).temperature_k(), 304.0);
        assert_eq!(endpoint(&law).frame(), frame(IDENTITY));
        assert!(law.mean_thermal_expansion().is_none());
    }
}

#[test]
fn g1_aligned_directional_and_isotropic_expansion_recover_their_free_shape() {
    let card = thermal_card("T", [290.0, 310.0]);
    let swap = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]];
    for exact in [false, true] {
        for isotropic in [false, true] {
            let original = sheet(&card, exact);
            let mapping = if isotropic {
                SheetExpansionFrame::Isotropic
            } else {
                frame(swap)
            };
            let target = original
                .clone()
                .with_hyperelastic_free_expansion(
                    &card,
                    path(&card, "T", isotropic),
                    310.0,
                    mapping,
                )
                .unwrap();
            let law = bound(&target);
            let source_g = endpoint(&law).geometry().principal_stretches();
            let g = if isotropic {
                source_g
            } else {
                [source_g[1], source_g[0], source_g[2]]
            };
            assert_eq!(law.reference_stretches(), g);
            assert_eq!(
                *endpoint(&law)
                    .reference_deformation()
                    .deformation_gradient(),
                [g[0], 0.0, 0.0, 0.0, g[1], 0.0, 0.0, 0.0, g[2]]
            );
            let free = [0.5 * (g[0] * g[0] - 1.0), 0.5 * (g[1] * g[1] - 1.0), 0.0];
            let response = law.response(free).unwrap();
            close(response.energy_j_m2, 0.0, 1e-10, 0.0);
            for stress in response.resultants_n_m {
                close(stress, 0.0, 1e-8, 0.0);
            }
            assert!(law.response([0.0; 3]).unwrap().energy_j_m2 > 0.01);
            let explicit = original.with_hyperelastic_reference_stretches(g).unwrap();
            let actual = law.response([0.12, 0.05, 0.08]).unwrap();
            let expected = explicit
                .membrane_law()
                .unwrap()
                .response([0.12, 0.05, 0.08])
                .unwrap();
            // The native oriented reference retains its source-order volume
            // product; a signed permutation can change only roundoff here.
            close(actual.energy_j_m2, expected.energy_j_m2, 1e-11, 1e-12);
            for (a, b) in actual
                .resultants_n_m
                .into_iter()
                .zip(expected.resultants_n_m)
            {
                close(a, b, 1e-10, 1e-12);
            }
            for (a, b) in actual.tangent_n_m.into_iter().zip(expected.tangent_n_m) {
                close(a, b, 1e-9, 1e-12);
            }
        }
    }
}

#[test]
fn g0_unsupported_frames_refuse_even_at_the_unit_reference_temperature() {
    let card = thermal_card("T", [290.0, 310.0]);
    let original = sheet(&card, true);
    let s = 0.5_f64.sqrt();
    for mapping in [
        SheetExpansionFrame::Isotropic,
        frame([[s, -s, 0.0], [s, s, 0.0], [0.0, 0.0, 1.0]]),
        frame([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        frame([[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        SheetExpansionFrame::FixedPrincipal {
            principal_frame: SHEET,
            sheet_frame: SHEET,
            principal_to_sheet: IDENTITY,
        },
        SheetExpansionFrame::FixedPrincipal {
            principal_frame: PRINCIPAL,
            sheet_frame: ContentHash([0; 32]),
            principal_to_sheet: IDENTITY,
        },
        SheetExpansionFrame::FixedPrincipal {
            principal_frame: PRINCIPAL,
            sheet_frame: PRINCIPAL,
            principal_to_sheet: [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
        },
    ] {
        assert!(
            original
                .clone()
                .with_hyperelastic_free_expansion(&card, path(&card, "T", false), 300.0, mapping)
                .is_err()
        );
    }
    assert_eq!(bound(&original).reference_stretches(), [1.0; 3]);
    assert!(bound(&original).free_expansion_reference().is_none());
}

#[test]
fn g0_source_datum_context_model_domain_and_existing_references_are_not_replaced() {
    let card = thermal_card("T", [295.0, 305.0]);
    let original = sheet(&card, false);
    let source = path(&card, "T", false);
    let selected = original
        .clone()
        .with_hyperelastic_free_expansion(&card, source.clone(), 302.0, frame(IDENTITY))
        .unwrap();
    let before = selected.clone();
    for t in [289.0, 308.0, 311.0, f64::NAN, f64::INFINITY] {
        assert!(selected.at_hyperelastic_temperature(&card, t).is_err());
    }
    let foreign = thermal_card("T", [290.0, 310.0]);
    assert!(
        selected
            .at_hyperelastic_temperature(&foreign, 303.0)
            .is_err()
    );
    assert!(
        original
            .clone()
            .with_hyperelastic_free_expansion(&foreign, source.clone(), 300.0, frame(IDENTITY))
            .is_err()
    );
    for query in [
        thermal_point("T", 302.0),
        thermal_point("T", 300.0).with("p", 2e5).unwrap(),
        thermal_point("T", 300.0).with("strain", 0.0).unwrap(),
    ] {
        let wrong = path_at(&card, "T", query, false);
        assert!(
            original
                .clone()
                .with_hyperelastic_free_expansion(&card, wrong, 300.0, frame(IDENTITY))
                .is_err()
        );
    }
    for existing in [
        selected.clone(),
        original
            .clone()
            .with_hyperelastic_reference_stretches([1.0; 3])
            .unwrap(),
    ] {
        assert!(
            existing
                .with_hyperelastic_free_expansion(&card, source.clone(), 300.0, frame(IDENTITY))
                .is_err()
        );
    }
    assert!(original.at_hyperelastic_temperature(&card, 302.0).is_err());
    assert_eq!(selected, before);
}

#[test]
fn g3_source_temperature_changes_reuse_native_membrane_memory_and_work() {
    let card = thermal_card("T", [290.0, 310.0]);
    let strain = [0.16, 0.08, 0.1];
    for exact in [false, true] {
        let original = sheet(&card, exact)
            .with_hyperelastic_free_expansion(
                &card,
                path(&card, "T", false),
                300.0,
                frame(IDENTITY),
            )
            .unwrap();
        let hot = original.at_hyperelastic_temperature(&card, 306.0).unwrap();
        let old = bound(&original);
        let target = bound(&hot);
        let law = ViscoHyperelastic::new(
            old.model().law().clone(),
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 1.6e5 },
                tau_s: 0.2,
            }],
        )
        .unwrap()
        .with_reference_deformation(*endpoint(&old).reference_deformation());
        let target_law = law
            .clone()
            .with_reference_deformation(*endpoint(&target).reference_deformation());
        let state = ViscoMembraneState::initial(&law, old.volumetric_treatment());
        let accepted =
            FrozenViscoMembrane::new(&law, &state, 0.08, original.thickness_m(), controls())
                .unwrap()
                .step_response(strain)
                .unwrap()
                .state;
        let before = accepted.clone();
        let changed = target_law
            .change_membrane_reference_from(
                &law,
                &accepted,
                strain,
                original.thickness_m(),
                controls(),
            )
            .unwrap();
        let metrics = |state: &ViscoMembraneState| match state {
            ViscoMembraneState::Full(s) => s.inverse_viscous_cauchy_green.clone(),
            ViscoMembraneState::Incompressible(s) => s.inverse_viscous_cauchy_green.clone(),
        };
        assert_eq!(metrics(&changed.state), metrics(&accepted));
        assert_eq!(changed.state.dissipated_pa(), accepted.dissipated_pa());
        close(
            changed.reference_work_j_m2 + changed.native_transverse_loading_work_j_m2,
            changed.stored_after_j_m2 - changed.stored_before_j_m2,
            1e-9,
            1e-12,
        );
        assert!(changed.reference_work_j_m2.abs() > 1e-3);
        let continued = FrozenViscoMembrane::new(
            &target_law,
            &changed.state,
            0.08,
            hot.thickness_m(),
            controls(),
        )
        .unwrap()
        .step_response(strain)
        .unwrap();
        assert!(continued.dissipated_j_m2 > 0.0);
        assert_eq!(accepted, before);
    }
}
