//! Original material queries on the same finite conservative network timeline.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::{
    osmosis::diffusion::source::{
        BINARY_SOLUTION_DIFFUSIVITY_DIMS as DIMS, BinarySolutionDiffusionSource,
        MaterialSolutionMaxwellStefan, SolutionDiffusionSourceAxes,
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn material(model: &Arc<VirialSolutionModel>, range: [f64; 2]) -> MaterialSolutionMaxwellStefan {
    let mut claims = ClaimSet::new();
    let key = PropertyKey::new("binary mobility", DIMS);
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "c_A".into(),
                abscissa_dims: typed(QuantityKind::AmountConcentration).dims(),
                knots: range
                    .into_iter()
                    .map(|c| (c, 1e-9 * (1.0 + c / 100.0)))
                    .collect(),
                dims: DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", typed(QuantityKind::AbsoluteTemperature), 290.0, 310.0)
                .with_quantity(
                    "c_A",
                    typed(QuantityKind::AmountConcentration),
                    range[0],
                    range[1],
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured nonlinear solution mobility".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "assigned pair".into(),
                phase: "liquid".into(),
                process: "manufactured reference".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    let solvent = SpeciesId::new("water").unwrap();
    let solute = SpeciesId::new("solute-0").unwrap();
    MaterialSolutionMaxwellStefan::new(
        model.clone(),
        (solvent.clone(), MolarMass::new(0.018)),
        vec![(solute.clone(), MolarMass::new(0.04))],
        vec![BinarySolutionDiffusionSource {
            species: [solvent.clone(), solute.clone()],
            material: card,
            property: key,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        SolutionDiffusionSourceAxes {
            temperature: "T".into(),
            concentrations: vec![(solvent, "c_water".into()), (solute, "c_A".into())],
        },
        cell(model, &[60.0]).solution(),
        || false,
    )
    .unwrap()
}
fn volume_cell(model: &Arc<VirialSolutionModel>, volume: f64, c: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        AqueousSolution::nonideal(model.clone(), volume, vec![volume * c]).unwrap(),
        OsmoticPressureVolume::new(volume, 100000.0, 1e8).unwrap(),
        0.1 * volume,
    )
    .unwrap()
}

fn multicomponent_material(model: &Arc<VirialSolutionModel>) -> MaterialSolutionMaxwellStefan {
    let key = PropertyKey::new("shared manufactured mixture mobility", DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "c_A".into(),
                abscissa_dims: typed(QuantityKind::AmountConcentration).dims(),
                knots: vec![(0.0, 1e-9), (200.0, 3e-9)],
                dims: DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "c_A",
                typed(QuantityKind::AmountConcentration),
                0.0,
                200.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured mixture diffusion".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "assigned liquid pairs".into(),
                phase: "liquid".into(),
                process: "analytic fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    let water = SpeciesId::new("water").unwrap();
    let a = SpeciesId::new("solute-0").unwrap();
    let b = SpeciesId::new("solute-1").unwrap();
    let pairs = [
        [water.clone(), a.clone()],
        [water.clone(), b.clone()],
        [a.clone(), b.clone()],
    ];
    MaterialSolutionMaxwellStefan::new(
        model.clone(),
        (water.clone(), MolarMass::new(0.018)),
        vec![
            (a.clone(), MolarMass::new(0.04)),
            (b.clone(), MolarMass::new(0.06)),
        ],
        pairs
            .into_iter()
            .map(|species| BinarySolutionDiffusionSource {
                species,
                material: card.clone(),
                property: key.clone(),
                selection: MaterialPropertySelection::SingleClaimOnly,
            })
            .collect(),
        SolutionDiffusionSourceAxes {
            temperature: "T".into(),
            concentrations: vec![
                (water, "c_water".into()),
                (a, "c_A".into()),
                (b, "c_B".into()),
            ],
        },
        cell(model, &[60.0, 60.0]).solution(),
        || false,
    )
    .unwrap()
}

#[test]
fn ionic_source_coefficients_refresh_at_the_same_neutral_face_used_for_transport() {
    use crate::osmosis::coupled::CoupledSolutionControls;
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let material = multicomponent_material(&model);
    let policy = ZeroCurrentDiffusion::new(vec![1, -1], 1e-15, 1e-16).unwrap();
    let (left, right) = (cell(&model, &[100.0, 100.0]), cell(&model, &[20.0, 20.0]));
    let mut flux = [0.0; 2];
    let report = material
        .flux_zero_current_into(
            left.solution(),
            right.solution(),
            0.001,
            physical().constitutive,
            &policy,
            &mut material.workspace().unwrap(),
            &mut flux,
            || false,
        )
        .unwrap();
    let mean = 80.0 / 5.0_f64.ln();
    let coefficient = 1e-9 * (1.0 + mean / 100.0);
    let expected = coefficient * 80.0 * (1.0 + 2.0 * VBAR * mean) / 0.001;
    for value in flux {
        assert!((value - expected).abs() < 2e-14);
    }
    for pair in report.sources.pairs() {
        assert!((pair.face.properties()[0].value_si() - coefficient).abs() < 1e-23);
    }
    assert!(
        report
            .constitutive
            .dimensionless_potential_right_minus_left
            .unwrap()
            .abs()
            < 1e-12
    );
    let mut system = OsmoticNetwork::with_zero_current_solution_transport(
        vec!["solute-0".into(), "solute-1".into()],
        vec![left, right],
        Vec::new(),
        vec![SolutionDiffusionLink::from_material_sources(0, 1, 0.01, 0.001, material).unwrap()],
        CoupledSolutionControls {
            diffusion: physical(),
            water_balance_tolerance_m3: 1e-18,
        },
        policy,
        1 << 22,
    )
    .unwrap();
    system.advance(10.0, numerical(1e-4), || false).unwrap();
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    let (a, b) = (
        left.solute_moles()[0] / VOLUME,
        right.solute_moles()[0] / VOLUME,
    );
    let expected = 1e-9 * (1.0 + ((a - b) / (a / b).ln()) / 100.0);
    let link = &system.solution_diffusion_links()[0];
    let resolved = link
        .material_sources()
        .unwrap()
        .resolve_zero_current_inventories(
            [VOLUME, VOLUME],
            [left.solute_moles(), right.solute_moles()],
            link.zero_current_diffusion().unwrap(),
            || false,
        )
        .unwrap();
    assert!((expected - coefficient).abs() > 1e-12);
    for pair in resolved.pairs() {
        assert!((pair.face.properties()[0].value_si() - expected).abs() < 1e-22);
    }
}

#[test]
fn multicomponent_implicit_candidates_refresh_all_original_binary_sources() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let material = multicomponent_material(&model);
    let mut system = OsmoticNetwork::with_solution_diffusion(
        vec!["solute-0".into(), "solute-1".into()],
        vec![cell(&model, &[100.0, 20.0]), cell(&model, &[20.0, 100.0])],
        vec![SolutionDiffusionLink::from_material_sources(0, 1, 0.01, 0.001, material).unwrap()],
        physical(),
        1 << 22,
    )
    .unwrap()
    .with_implicit_multicomponent_diffusion(
        ImplicitMulticomponentDiffusionControls::default(),
        1 << 22,
    )
    .unwrap();
    let saved = bits(&system);
    let mut implicit = system.implicit_multicomponent.take().unwrap();
    let mut transfers = [f64::NAN; 2];
    let mut evaluations = 0;
    implicit
        .transfers_into(
            &system.cells[0].solution,
            &system.cells[1].solution,
            &mut system.diffusion[0],
            100.0,
            physical().constitutive,
            &mut evaluations,
            &mut transfers,
            &mut || false,
        )
        .unwrap();
    let left = [100.0 * VOLUME - transfers[0], 20.0 * VOLUME - transfers[1]];
    let right = [20.0 * VOLUME + transfers[0], 100.0 * VOLUME + transfers[1]];
    let (cl, cr) = (left[0] / VOLUME, right[0] / VOLUME);
    let mean = (cl - cr) / (cl / cr).ln();
    let coefficient = 1e-9 * (1.0 + mean / 100.0);
    let expected = 100.0 * 0.01 * coefficient * (cl - cr) / 0.001;
    assert!((transfers[0] - expected).abs() < 5e-14);
    assert!((transfers[0] + transfers[1]).abs() < 2e-14);
    assert!((coefficient - 1e-9 * (1.0 + (80.0 / 5.0_f64.ln()) / 100.0)).abs() > 5e-11);
    let original = system.diffusion[0].material_sources().unwrap();
    let endpoint = original
        .resolve_inventories([VOLUME, VOLUME], [&left, &right], || false)
        .unwrap();
    for pair in endpoint.pairs() {
        assert!((pair.face.properties()[0].value_si() - coefficient).abs() < 1e-23);
        assert!(
            original
                .sources()
                .iter()
                .any(|source| source.selected_claim()
                    == pair.face.properties()[0].answer().receipt.selected)
        );
    }
    assert!(evaluations > 4);
    assert_eq!(bits(&system), saved);
}

#[test]
fn implicit_binary_root_refreshes_original_source_at_its_actual_endpoint() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let mut system = OsmoticNetwork::with_solution_diffusion(
        vec!["solute-0".into()],
        vec![cell(&model, &[100.0]), cell(&model, &[20.0])],
        vec![
            SolutionDiffusionLink::from_material_sources(
                0,
                1,
                0.01,
                0.001,
                material(&model, [1.0, 200.0]),
            )
            .unwrap(),
        ],
        physical(),
        1 << 22,
    )
    .unwrap()
    .with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls::default(), 1 << 22)
    .unwrap();
    let old = bits(&system);
    let mut solver = system.implicit_diffusion.take().unwrap();
    let mut evaluations = 0;
    let dt = 100.0;
    let transfer = solver
        .transfer(
            &system.cells[0].solution,
            &system.cells[1].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut || false,
        )
        .unwrap();
    let next = [100.0 * VOLUME - transfer, 20.0 * VOLUME + transfer];
    let (a, b) = (next[0] / VOLUME, next[1] / VOLUME);
    let mean = (a - b) / (a / b).ln();
    let coefficient = 1e-9 * (1.0 + mean / 100.0);
    let expected = dt * 0.01 * coefficient * (a - b) * (1.0 + VBAR * mean) / 0.001;
    assert!((transfer - expected).abs() <= 2.0 * (1e-20 + 1e-10 * 40.0 * VOLUME));
    let initial_mean = 80.0 / 5.0_f64.ln();
    assert!((coefficient - 1e-9 * (1.0 + initial_mean / 100.0)).abs() > 5e-11);
    let original = system.diffusion[0].material_sources().unwrap();
    let endpoint = original
        .resolve_inventories([VOLUME, VOLUME], [&next[..1], &next[1..]], || false)
        .unwrap();
    assert!((endpoint.pairs()[0].face.properties()[0].value_si() - coefficient).abs() < 1e-23);
    assert_eq!(
        endpoint.pairs()[0].face.properties()[0]
            .answer()
            .receipt
            .selected,
        original.sources()[0].selected_claim()
    );
    assert!(evaluations > 2);
    assert_eq!(bits(&system), old);
}

#[test]
fn g1_evolving_source_mobility_tracks_independent_unequal_volume_relaxation() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let sourced = material(&model, [1.0, 200.0]);
    let constant = sourced.reference_law().clone();
    let build = |link| {
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute-0".into()],
            vec![
                volume_cell(&model, VOLUME, 100.0),
                volume_cell(&model, 2.0 * VOLUME, 20.0),
            ],
            vec![link],
            physical(),
            1 << 22,
        )
        .unwrap()
    };
    let mut system =
        build(SolutionDiffusionLink::from_material_sources(0, 1, 0.01, 0.001, sourced).unwrap());
    let mut frozen = build(SolutionDiffusionLink::new(0, 1, 0.01, 0.001, constant).unwrap());
    // Independent scalar ODE uses conserved total amount and unequal volumes.
    // The source follows the same logarithmic-mean face as molecular friction.
    let rate = |left: f64| {
        let right = 0.5 * (140.0 - left);
        let mean = (left - right) / (left / right).ln();
        -0.01 / VOLUME * 1e-9 * (1.0 + mean / 100.0) * (left - right) * (1.0 + mean * VBAR) / 0.001
    };
    let mut exact = 100.0;
    let h = 20.0 / 4000.0;
    for _ in 0..4000 {
        let a = rate(exact);
        let b = rate(exact + 0.5 * h * a);
        let c = rate(exact + 0.5 * h * b);
        let d = rate(exact + h * c);
        exact += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let report = system.advance(20.0, numerical(1e-5), || false).unwrap();
    frozen.advance(20.0, numerical(1e-5), || false).unwrap();
    assert!((concentration(&system, 0, 0) - exact).abs() < 0.03);
    assert!((concentration(&system, 0, 0) - concentration(&frozen, 0, 0)).abs() > 0.1);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
    assert!(report.solute_balance_residual_mol[0].abs() <= physical().amount_balance_tolerance_mol);
    assert_eq!(report.water_balance_residual_m3, 0.0);
    let law = system.solution_diffusion_links()[0]
        .material_sources()
        .unwrap();
    let (left, right) = (
        &system.compartments()[0].solution,
        &system.compartments()[1].solution,
    );
    let terminal = law
        .resolve_inventories(
            [left.solvent_volume_m3(), right.solvent_volume_m3()],
            [left.solute_moles(), right.solute_moles()],
            || false,
        )
        .unwrap();
    let (a, b) = (concentration(&system, 0, 0), concentration(&system, 1, 0));
    let mean = (a - b) / (a / b).ln();
    assert!(
        (terminal.pairs()[0].face.properties()[0].value_si() - 1e-9 * (1.0 + mean / 100.0)).abs()
            < 1e-23
    );
    assert_eq!(
        terminal.pairs()[0].face.properties()[0]
            .answer()
            .receipt
            .selected,
        law.sources()[0].selected_claim()
    );
}

fn branched(range: [f64; 2]) -> OsmoticNetwork {
    let model = source(1.0, vec![0.0], vec![200.0]);
    OsmoticNetwork::with_solution_diffusion(
        vec!["solute-0".into()],
        vec![
            cell(&model, &[60.0]),
            cell(&model, &[60.0]),
            cell(&model, &[1.0]),
        ],
        vec![
            SolutionDiffusionLink::from_material_sources(
                0,
                1,
                0.001,
                0.001,
                material(&model, range),
            )
            .unwrap(),
            SolutionDiffusionLink::new(1, 2, 0.01, 0.001, law(&model)).unwrap(),
        ],
        physical(),
        1 << 22,
    )
    .unwrap()
}

#[test]
fn g4_later_edge_must_not_invalidate_an_earlier_links_source_support() {
    let mut system = branched([59.9, 100.0]);
    let original = bits(&system);
    // The sourced first link initially has zero flux and valid endpoints.
    // The later constant link withdraws from its endpoint; only a whole-sweep
    // source check can refuse that newly unsupported final state.
    assert!(matches!(
        system.advance(
            1.0,
            OsmoticStepControls {
                maximum_step_s: 1.0,
                ..numerical(1e-3)
            },
            || false
        ),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::Material(_)
        ))
    ));
    assert_eq!(bits(&system), original);
}

#[test]
fn g4_late_pinned_source_failure_preserves_the_last_complete_window() {
    let mut system = branched([50.0, 100.0]);
    let controls = OsmoticStepControls {
        maximum_step_s: 1.0,
        ..numerical(1e-4)
    };
    system.advance(2.0, controls, || false).unwrap();
    let original = bits(&system);
    assert_eq!(system.time_s(), 2.0);
    assert!(matches!(
        system.advance(50.0, controls, || false),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::Material(_)
        ))
    ));
    assert_eq!(bits(&system), original);
}

#[test]
fn g4_source_network_charges_retained_cards_and_refuses_initial_source_gaps() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let build = |budget, high| {
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute-0".into()],
            vec![cell(&model, &[high]), cell(&model, &[60.0])],
            vec![
                SolutionDiffusionLink::from_material_sources(
                    0,
                    1,
                    0.01,
                    0.001,
                    material(&model, [50.0, 100.0]),
                )
                .unwrap(),
            ],
            physical(),
            budget,
        )
    };
    let system = build(1 << 22, 80.0).unwrap();
    let source = system.solution_diffusion_links()[0]
        .material_sources()
        .unwrap();
    assert!(
        system.workspace_bytes()
            > source.sources()[0]
                .material()
                .estimated_retained_payload_bytes()
    );
    assert!(matches!(
        build(system.workspace_bytes() - 1, 80.0),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    // The logarithmic-mean face is admitted, but the endpoint c=120 is not.
    assert!(matches!(
        build(1 << 22, 120.0),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::Material(_)
        ))
    ));
}

#[test]
fn g1_source_bound_front_enters_pure_solvent_only_with_supported_zero_coordinates() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let build = |range| {
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute-0".into()],
            vec![cell(&model, &[100.0]), cell(&model, &[0.0])],
            vec![
                SolutionDiffusionLink::from_material_sources(
                    0,
                    1,
                    0.01,
                    0.001,
                    material(&model, range),
                )
                .unwrap(),
            ],
            physical(),
            1 << 22,
        )
    };
    // A material curve that starts at c=1 cannot acquire authority at c=0.
    assert!(matches!(
        build([1.0, 200.0]),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::Material(_)
        ))
    ));
    let mut sourced = build([0.0, 200.0]).unwrap();
    let mut constant = network(&model, &[&[100.0], &[0.0]], &[(0, 1, 0.01)]);
    let report = sourced.advance(10.0, numerical(1e-4), || false).unwrap();
    constant.advance(10.0, numerical(1e-4), || false).unwrap();
    assert!(concentration(&sourced, 1, 0) > concentration(&constant, 1, 0));
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
    assert!(report.solute_balance_residual_mol[0].abs() <= physical().amount_balance_tolerance_mol);
    let source = sourced.solution_diffusion_links()[0]
        .material_sources()
        .unwrap();
    let left = sourced.compartments()[0].solution();
    let right = sourced.compartments()[1].solution();
    let terminal = source
        .resolve_inventories(
            [left.solvent_volume_m3(), right.solvent_volume_m3()],
            [left.solute_moles(), right.solute_moles()],
            || false,
        )
        .unwrap();
    assert_eq!(
        terminal.pairs()[0].face.properties()[0]
            .answer()
            .receipt
            .selected,
        source.sources()[0].selected_claim()
    );
}
