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
    // Unlike an equal-volume pair, its arithmetic face concentration evolves.
    let rate = |left: f64| {
        let right = 0.5 * (140.0 - left);
        let mean = 0.5 * (left + right);
        -0.01 / VOLUME
            * 1e-9
            * (1.0 + mean / 100.0)
            * mean
            * (1.0 + mean * VBAR)
            * (left / right).ln()
            / 0.001
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
    let mean = 0.5 * (concentration(&system, 0, 0) + concentration(&system, 1, 0));
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
    // Arithmetic face c=90 is admitted, but the endpoint c=120 is not.
    assert!(matches!(
        build(1 << 22, 120.0),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::Material(_)
        ))
    ));
}
