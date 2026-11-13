//! G0/G1/G3: source-bound ductility at native shell thickness points.
//! Synthetic coefficients exercise mechanics and history, not measured alloys.
#![deny(unsafe_code)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::ductile::{DuctileDomain, JohnsonCookInitiation, JohnsonCookParameters};
use fs_material::ductile::source::{
    DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY, PLASTIC_STRAIN_RATE_AXIS,
    ResolvedDuctileLocus, STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::{
    LayeredPlaneStressSection, PlaneStressLayer, SectionBudget, ThicknessQuadrature,
};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::{SmallStrainLaw, Voigt};
use fs_plate::shell::{DRILLING_ALPHA, ShellMesh};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec};
use fs_solid::ductile::{
    DuctileAdvanceError, DuctileHistory, DuctileLawField, DuctileSolidError,
    DuctileStepConditions, DuctileTemperatures,
};
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;
use fs_solid::plastic_shell::{
    PlasticShellBudget, PlasticShellSections, SmallStrainShellProblem,
};
use fs_solid::IncrementSettings;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 83, kernel_id: 131, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&gate, &cx)
    })
}

fn close(actual: f64, expected: f64, allowance: f64) {
    assert!((actual - expected).abs() <= allowance,
        "actual {actual:e}, expected {expected:e}, allowance {allowance:e}");
}

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}

fn fixed_point() -> QueryPoint {
    QueryPoint::new().with_quantity("T", temperature(), 300.0).unwrap()
}

fn card(capacity: f64) -> MaterialCard {
    let fixed = ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 300.0);
    let provenance = Provenance {
        source: "synthetic ductile shell material".into(), license: "test-only".into(), artifact: None,
    };
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, 1000.0),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("yield_stress", Pressure::DIMS, 1.0),
        ("isotropic_hardening_modulus", Pressure::DIMS, 100.0),
        ("kinematic_hardening_modulus", Pressure::DIMS, 50.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: fixed.clone(), uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: provenance.clone(),
        }).unwrap();
    }
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE),
        value: PropertyValue::Curve {
            abscissa: STRESS_TRIAXIALITY_AXIS.into(), abscissa_dims: Dims::NONE,
            // epsilon_f = capacity * (1 - 0.2 * eta), so opposite bending
            // layers have different utilization despite equal plastic alpha.
            knots: vec![(-2.0, 1.4 * capacity), (2.0, 0.6 * capacity)], dims: Dims::NONE,
        },
        validity: fixed
            .with_quantity(STRESS_TRIAXIALITY_AXIS, QuantitySpec::dimensional(Dims::NONE), -2.0, 2.0)
            .with_quantity(PLASTIC_STRAIN_RATE_AXIS,
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])), 1e-10, 1.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::TabulatedOnly,
        observations: vec![], provenance,
    }).unwrap();
    MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic ductile sheet".into(), phase: "solid".into(),
        process: "synthetic coupon".into(), revision: 0,
    }, claims, vec![]).unwrap()
}

fn resolve(card: &MaterialCard) -> ResolvedJ2Plasticity {
    resolve_j2_plasticity(card, &fixed_point(), MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::MixedPlasticStrainModuli, 0.05).unwrap()
}

fn locus(card: &MaterialCard) -> ResolvedDuctileLocus {
    // The resolver selects the claim at one of its original knots (-2 and 2), never between them.
    let query = fixed_point()
        .with_quantity(STRESS_TRIAXIALITY_AXIS, QuantitySpec::dimensional(Dims::NONE), 2.0).unwrap()
        .with_quantity(PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])), 0.01).unwrap();
    resolve_ductile_locus(card, &query,
        &DuctileSourceRequest { temperature_axis: "T".into(), maximum_points: 2 },
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)).unwrap()
}

fn section(lower: ResolvedJ2Plasticity, upper: ResolvedJ2Plasticity)
    -> LayeredPlaneStressSection<ResolvedJ2Plasticity>
{
    LayeredPlaneStressSection::new(vec![
        PlaneStressLayer { bottom_z_m: -0.05, top_z_m: 0.0, subdivisions: 1, law: lower },
        PlaneStressLayer { bottom_z_m: 0.0, top_z_m: 0.05, subdivisions: 1, law: upper },
    ], ThicknessQuadrature::Gauss3,
        ReductionConfig { max_iterations: 48, stress_tolerance_pa: 1e-12, strain_limit: 0.05 },
        SectionBudget { maximum_layers: 2, maximum_quadrature_points: 6 }).unwrap()
}

fn mesh() -> ShellMesh {
    ShellMesh::new(vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]], vec![[0, 1, 2]]).unwrap()
}

fn bending(mesh: &ShellMesh, curvature: f64) -> Vec<(usize, f64)> {
    mesh.nodes.iter().enumerate().flat_map(|(node, x)| {
        [0.0, 0.0, 0.5 * curvature * x[0] * x[0], 0.0, -curvature * x[0], 0.0]
            .into_iter().enumerate().map(move |(c, value)| (6 * node + c, value))
    }).collect()
}

fn problem<'a>(mesh: &'a ShellMesh, section: &'a LayeredPlaneStressSection<ResolvedJ2Plasticity>,
    prescribed: &'a [(usize, f64)], loads: &'a [(usize, f64)])
    -> SmallStrainShellProblem<'a, ResolvedJ2Plasticity>
{
    SmallStrainShellProblem {
        mesh, sections: PlasticShellSections::Uniform(section), coordinate_length_m: 1.0,
        rotation_limit_rad: 0.05, drilling_alpha: DRILLING_ALPHA,
        prescribed, nodal_loads: loads, load_path: &[1.0], budget: PlasticShellBudget::default(),
        settings: IncrementSettings { max_newton_iterations: 40, ..IncrementSettings::default() },
    }
}

fn conditions(duration_s: f64) -> DuctileStepConditions<'static> {
    DuctileStepConditions { duration_s, temperatures_k: DuctileTemperatures::Uniform(300.0) }
}

fn triaxiality(stress: &Voigt) -> f64 {
    let q2 = 0.5 * ((stress[0] - stress[1]).powi(2) + (stress[1] - stress[2]).powi(2)
        + (stress[2] - stress[0]).powi(2)) + 3.0 * stress[3..].iter().map(|s| s * s).sum::<f64>();
    (stress[0] + stress[1] + stress[2]) / (3.0 * q2.sqrt())
}

#[test]
fn bending_keeps_opposite_layer_stresses_and_every_surface_point() {
    with_cx(|_, cx| {
        let card = card(0.02);
        let metal = resolve(&card);
        let material = section(metal.clone(), metal.clone());
        let locus = locus(&card);
        let mesh = mesh();
        let prescribed = bending(&mesh, 0.04);
        let p = problem(&mesh, &material, &prescribed, &[]);
        let history = DuctileHistory::virgin_shell(DuctileLawField::Uniform(&locus), &p, cx).unwrap();
        assert_eq!(history.points().len(), 18);
        let accepted = history.advance_shell(&p, None, conditions(1.0), cx).unwrap();
        let direct = p.solve(cx).unwrap();
        assert_eq!(accepted.mechanics, direct);
        let mut yielded = 0;
        for (surface, state) in accepted.mechanics.section_states.iter().enumerate() {
            close(state.generalized_strain[3], -0.04, 1e-15);
            for force in &accepted.mechanics.section_resultants[surface][..3] {
                close(*force, 0.0, 1e-12);
            }
            for (thickness, native) in state.points.iter().enumerate() {
                let plastic = native.history.plastic_strain;
                let elastic = core::array::from_fn(|i| native.strain[i] - plastic[i]);
                let stress = metal.law().elastic.stress(&elastic, &());
                let alpha = native.history.alpha;
                let usage = accepted.ductility.points()[6 * surface + thickness].utilization();
                close(accepted.ductility.points()[6 * surface + thickness].accumulated_plastic_strain(), alpha, 0.0);
                if alpha > 0.0 {
                    yielded += 1;
                    assert!(stress[0] * material.points()[thickness].z_m < 0.0);
                    close(usage, alpha / (0.02 * (1.0 - 0.2 * triaxiality(&stress))), 1e-13);
                } else {
                    assert_eq!(usage, 0.0);
                }
            }
            let bottom = &accepted.ductility.points()[6 * surface];
            let top = &accepted.ductility.points()[6 * surface + 5];
            close(bottom.accumulated_plastic_strain(), top.accumulated_plastic_strain(), 1e-15);
            assert!(bottom.utilization() > top.utilization());
        }
        assert_eq!(yielded, 6, "only the outermost pair yields at each surface point");
        assert!(history.points().iter().all(|p| p.utilization() == 0.0));
    });
}

#[test]
fn unload_reverse_and_restart_retain_irreversible_layer_history() {
    with_cx(|_, cx| {
        let card = card(0.02);
        let metal = resolve(&card);
        let material = section(metal.clone(), metal);
        let locus = locus(&card);
        let mesh = mesh();
        let prescribed = bending(&mesh, 0.04);
        let p = problem(&mesh, &material, &prescribed, &[]);
        let history = DuctileHistory::virgin_shell(DuctileLawField::Uniform(&locus), &p, cx).unwrap();
        let first = history.advance_shell(&p, None, conditions(1.0), cx).unwrap();
        let unloaded_target = bending(&mesh, 0.039);
        let unloaded = first.ductility.advance_shell(&problem(&mesh, &material, &unloaded_target, &[]),
            Some(&first.mechanics), conditions(1.0), cx).unwrap();
        assert_eq!(unloaded.ductility.points(), first.ductility.points());
        let reverse_target = bending(&mesh, -0.04);
        let reverse_problem = problem(&mesh, &material, &reverse_target, &[]);
        let restart = unloaded.ductility.clone();
        let reversed = unloaded.ductility.advance_shell(&reverse_problem, Some(&unloaded.mechanics), conditions(1.0), cx).unwrap();
        let replay = restart.advance_shell(&reverse_problem, Some(&unloaded.mechanics), conditions(1.0), cx).unwrap();
        assert_eq!(reversed.mechanics, replay.mechanics);
        assert_eq!(reversed.ductility.points(), replay.ductility.points());
        for (before, after) in unloaded.ductility.points().iter().zip(reversed.ductility.points()) {
            assert!(after.utilization() >= before.utilization());
            if before.utilization() > 0.0 {
                assert!(after.accumulated_plastic_strain() > before.accumulated_plastic_strain());
                assert!(after.utilization() > before.utilization());
            }
        }
        assert!(matches!(first.ductility.advance_shell(&p, None, conditions(1.0), cx),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))));
    });
}

#[test]
fn initiation_refuses_the_whole_shell_and_refinement_reuses_the_same_checkpoint() {
    with_cx(|gate, cx| {
        let card = card(0.001);
        let metal = resolve(&card);
        let material = section(metal.clone(), metal);
        let locus = locus(&card);
        let mesh = mesh();
        let first_target = bending(&mesh, 0.04);
        let p = problem(&mesh, &material, &first_target, &[]);
        let history = DuctileHistory::virgin_shell(DuctileLawField::Uniform(&locus), &p, cx).unwrap();
        let first = history.advance_shell(&p, None, conditions(1.0), cx).unwrap();
        let saved = first.mechanics.clone();
        let saved_history = first.ductility.points().to_vec();
        let excessive = bending(&mesh, 0.08);
        let refused = first.ductility.advance_shell(&problem(&mesh, &material, &excessive, &[]),
            Some(&first.mechanics), conditions(1.0), cx).unwrap_err();
        match refused {
            DuctileAdvanceError::Ductility(DuctileSolidError::Initiated { point, candidate }) => {
                assert_eq!(point % 6, 0, "tensile outer layer reaches its lower capacity first");
                assert!(candidate.state.initiated());
                assert!(candidate.initiation_fraction.unwrap() > 0.0);
                assert!(candidate.initiation_fraction.unwrap() < 1.0);
            },
            error => panic!("expected ductile initiation, got {error:?}"),
        }
        assert_eq!(first.mechanics, saved);
        assert_eq!(first.ductility.points(), saved_history);
        let refined_target = bending(&mesh, 0.042);
        let refined_problem = problem(&mesh, &material, &refined_target, &[]);
        let retry = first.ductility.advance_shell(&refined_problem, Some(&first.mechanics), conditions(1.0), cx).unwrap();
        assert!(retry.ductility.points()[0].utilization() > first.ductility.points()[0].utilization());
        assert!(retry.ductility.points().iter().all(|p| !p.initiated()));
        gate.request();
        assert!(matches!(first.ductility.advance_shell(&refined_problem, Some(&first.mechanics), conditions(1.0), cx),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Cancelled))));
        assert_eq!(first.mechanics, saved);
        assert_eq!(first.ductility.points(), saved_history);
    });
}

#[test]
fn heterogeneous_layers_bind_exact_source_cards_and_point_order() {
    with_cx(|_, cx| {
        let high_card = card(0.02);
        let low_card = card(0.001);
        let high = locus(&high_card);
        let low = locus(&low_card);
        let material = section(resolve(&high_card), resolve(&low_card));
        let mesh = mesh();
        let prescribed = bending(&mesh, 0.04);
        let p = problem(&mesh, &material, &prescribed, &[]);
        let wrong = DuctileHistory::virgin_shell(DuctileLawField::Uniform(&high), &p, cx).unwrap();
        assert!(matches!(wrong.advance_shell(&p, None, conditions(1.0), cx),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))));
        assert!(wrong.points().iter().all(|p| p.utilization() == 0.0));
        let point_laws: Vec<_> = (0..3).flat_map(|_| (0..6).map(|q| if q < 3 { high.clone() } else { low.clone() })).collect();
        let history = DuctileHistory::virgin_shell(DuctileLawField::PerPoint(&point_laws), &p, cx).unwrap();
        let accepted = history.advance_shell(&p, None, conditions(1.0), cx).unwrap();
        assert!(accepted.ductility.points()[5].utilization() > 10.0 * accepted.ductility.points()[0].utilization());
        let wrong_temperatures = [300.0; 17];
        let incomplete = DuctileStepConditions { duration_s: 1.0,
            temperatures_k: DuctileTemperatures::PerPoint(&wrong_temperatures) };
        assert!(matches!(history.advance_shell(&p, None, incomplete, cx),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))));
        let mut bounded = p;
        bounded.budget.maximum_material_points = 17;
        assert!(DuctileHistory::virgin_shell(DuctileLawField::Uniform(&high), &bounded, cx).is_err());
    });
}

#[test]
fn native_momentum_and_rate_share_one_duration_without_replaying_plasticity() {
    with_cx(|_, cx| {
        let source = card(0.05);
        let material = section(resolve(&source), resolve(&source));
        let mesh = mesh();
        let fixed: Vec<_> = (0..18).filter(|d| *d != 6 && *d != 13).map(|d| (d, 0.0)).collect();
        let loads = [(6, 0.024)];
        let p = problem(&mesh, &material, &fixed, &loads);
        let calibration = JohnsonCookInitiation::new(JohnsonCookParameters {
            d1: 0.02, d2: 0.0, d3: 0.0, d4: 0.02, d5: 0.0,
            reference_strain_rate_s_inv: 0.01,
            transition_temperature_k: 290.0, melting_temperature_k: 900.0,
        }, DuctileDomain {
            triaxiality: [-2.0, 2.0], strain_rate_s_inv: [1e-10, 1.0], temperature_k: [290.0, 310.0],
        }).unwrap();
        let history = DuctileHistory::virgin_shell(DuctileLawField::Uniform(&calibration), &p, cx).unwrap();
        let reference = p.prepare_reference_inertia(PlasticShellDensity::Uniform(2.0), cx).unwrap();
        let velocity = [[0.0; 6]; 3];
        let dt = 0.2;
        let direct = p.advance_backward_euler(None, &velocity, &reference, dt, cx).unwrap();
        let accepted = history.advance_shell_backward_euler(&p, None, &velocity, &reference, conditions(dt), cx).unwrap();
        assert_eq!(accepted.mechanics.solution, direct.solution);
        assert_eq!(accepted.mechanics.velocity, direct.velocity);
        assert_eq!(accepted.mechanics.inertia, direct.inertia);
        assert!(accepted.mechanics.inertia.physical.kinetic_after_j > 0.0);
        for (native, ductile) in accepted.mechanics.solution.section_states.iter().flat_map(|s| &s.points).zip(accepted.ductility.points()) {
            let alpha = native.history.alpha;
            assert!(alpha > 0.0);
            let fracture_strain = 0.02 * (1.0 + 0.02 * (alpha / (dt * 0.01)).ln());
            close(ductile.utilization(), alpha / fracture_strain, 1e-12);
        }
        let false_temperature = DuctileStepConditions { duration_s: dt,
            temperatures_k: DuctileTemperatures::Uniform(310.0) };
        // The numerical initiation domain allows 310 K; the actual frozen
        // mechanical source was queried at 300 K, so the adapter must refuse.
        assert!(matches!(history.advance_shell_backward_euler(&p, None, &velocity, &reference, false_temperature, cx),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))));
        assert!(history.points().iter().all(|p| p.utilization() == 0.0));
        assert_eq!(velocity, [[0.0; 6]; 3]);
    });
}
