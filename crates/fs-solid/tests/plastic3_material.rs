//! Material-card J2 data drive the existing P1 solid solver unchanged.
//!
//! G1: the same uniaxial coupon follows an independent bilinear stress law
//! after yield, hardening, elastic, or temperature data change.
//! G3: a plastic-strain hardening modulus and the equivalent measured
//! total-strain tangent produce the same loading and permanent deformation.
//! The cards are synthetic; these are coupling regressions, not calibrations.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::J2State;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY,
    POST_YIELD_TANGENT_MODULUS_PROPERTY, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, Pressure, Temperature};
use fs_solid::{
    IncrementRecord, IncrementSettings, SmallStrainLawField, SmallStrainTetProblem,
    SmallStrainTetSolution, TetAssemblyBudget,
};

#[derive(Clone, Copy)]
struct Parameters {
    young: f64,
    poisson: f64,
    yield_stress: f64,
    hardening: f64,
}

fn material_card(
    name: &str,
    cold: Parameters,
    hot: Parameters,
    basis: HardeningPropertyBasis,
) -> MaterialCard {
    let (hardening_name, source_hardening) = match basis {
        HardeningPropertyBasis::PlasticStrainModulus => (
            ISOTROPIC_HARDENING_MODULUS_PROPERTY,
            [cold.hardening, hot.hardening],
        ),
        HardeningPropertyBasis::UniaxialTangentModulus => (
            POST_YIELD_TANGENT_MODULUS_PROPERTY,
            [cold, hot].map(|p| p.young * p.hardening / (p.young + p.hardening)),
        ),
        HardeningPropertyBasis::MixedPlasticStrainModuli => {
            panic!("mixed coupon needs an explicitly added independent kinematic source")
        }
    };
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, [cold.young, hot.young]),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), [cold.poisson, hot.poisson]),
        (
            YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            [cold.yield_stress, hot.yield_stress],
        ),
        (hardening_name, Pressure::DIMS, source_hardening),
    ] {
        let (value, interpolation) = if values[0] == values[1] {
            (
                PropertyValue::Scalar { value: values[0], dims },
                InterpolationPolicy::ConstantWithinValidity,
            )
        } else {
            (
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (600.0, values[1])],
                    dims,
                },
                InterpolationPolicy::LinearInside,
            )
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value,
            validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic uniaxial J2 material regression".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "synthetic coupon".into(),
        revision: 0,
    }, claims, vec![]).unwrap()
}

fn resolve(
    card: &MaterialCard,
    temperature: f64,
    basis: HardeningPropertyBasis,
) -> ResolvedJ2Plasticity {
    let resolved = resolve_j2_plasticity(
        card,
        &QueryPoint::new().with("T", temperature).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        basis,
        0.05,
    ).unwrap();
    assert_eq!(resolved.material().card_identity(), card.content_hash());
    assert_eq!(
        resolved.material().query_point(),
        &[("T".to_owned(), temperature)],
    );
    assert_eq!(resolved.hardening_basis(), basis);
    resolved
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 5, kernel_id: 23, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

// One 2 x 1 x 1 m box, split into six positively oriented P1 tetrahedra.
// The exact solution is affine, so this small mesh tests actual equilibrium
// and free lateral motion without introducing a mesh-convergence assumption.
const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [2.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [2.0, 0.0, 1.0], [0.0, 1.0, 1.0], [2.0, 1.0, 1.0],
];
const TETS: [[usize; 4]; 6] = [
    [0, 1, 3, 7], [0, 1, 7, 5], [0, 2, 7, 3],
    [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 7, 6],
];
const PATH: [f64; 5] = [0.00025, 0.00075, 0.0015, 0.003, 0.006];

fn prescribed() -> Vec<(usize, f64)> {
    let mut prescribed = Vec::new();
    for (node, x) in NODES.iter().enumerate() {
        // Unit axial strain on the far face; the near face is a roller.
        prescribed.push((3 * node, x[0]));
        if x[1] == 0.0 {
            prescribed.push((3 * node + 1, 0.0));
        }
        if x[2] == 0.0 {
            prescribed.push((3 * node + 2, 0.0));
        }
    }
    prescribed
}

fn solve(resolved: &ResolvedJ2Plasticity, path: &[f64]) -> SmallStrainTetSolution<J2State> {
    let prescribed = prescribed();
    let problem = SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(resolved),
        prescribed_m: &prescribed,
        nodal_forces_n: &[],
        load_path: path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    with_cx(|cx| problem.solve(cx)).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.15e} vs {expected:.15e}, tolerance {tolerance:e}",
    );
}

fn axial_reaction(record: &IncrementRecord) -> f64 {
    let prescribed = prescribed();
    assert_eq!(record.reactions_n.len(), prescribed.len());
    prescribed.iter().zip(&record.reactions_n)
        .filter(|((dof, value), _)| dof % 3 == 0 && *value > 0.0)
        .map(|(_, reaction)| reaction)
        .sum()
}

// Uniaxial Hooke elasticity plus sigma = fy + H * alpha and
// epsilon = sigma / E + alpha gives this independent scalar closed form.
fn loading_response(p: Parameters, epsilon: f64) -> (f64, f64) {
    if p.young * epsilon <= p.yield_stress {
        (p.young * epsilon, 0.0)
    } else {
        let sigma = p.young * (p.yield_stress + p.hardening * epsilon)
            / (p.young + p.hardening);
        (sigma, epsilon - sigma / p.young)
    }
}

fn assert_fields(
    solution: &SmallStrainTetSolution<J2State>,
    p: Parameters,
    epsilon: f64,
    sigma: f64,
    alpha: f64,
) {
    let lateral = -p.poisson * sigma / p.young - alpha / 2.0;
    assert_eq!(solution.displacement_m.len(), NODES.len());
    for (x, u) in NODES.iter().zip(&solution.displacement_m) {
        for (actual, expected) in u.iter().zip([
            epsilon * x[0], lateral * x[1], lateral * x[2],
        ]) {
            close(*actual, expected, 1e-9);
        }
    }
    assert_eq!(solution.element_strain.len(), TETS.len());
    assert_eq!(solution.element_stress.len(), TETS.len());
    assert_eq!(solution.element_states.len(), TETS.len());
    for element in 0..TETS.len() {
        let expected_strain = [epsilon, lateral, lateral, 0.0, 0.0, 0.0];
        let expected_stress = [sigma, 0.0, 0.0, 0.0, 0.0, 0.0];
        let expected_plastic = [alpha, -alpha / 2.0, -alpha / 2.0, 0.0, 0.0, 0.0];
        for component in 0..6 {
            close(solution.element_strain[element][component], expected_strain[component], 1e-9);
            close(solution.element_stress[element][component], expected_stress[component], 2e-8);
            close(
                solution.element_states[element].plastic_strain[component],
                expected_plastic[component],
                1e-9,
            );
        }
        close(solution.element_states[element].alpha, alpha, 1e-9);
    }
}

fn assert_loading(resolved: &ResolvedJ2Plasticity, p: Parameters) -> SmallStrainTetSolution<J2State> {
    let solution = solve(resolved, &PATH);
    assert_eq!(solution.increments.len(), PATH.len());
    for (record, &epsilon) in solution.increments.iter().zip(&PATH) {
        // Area = 1 m2, so axial force equals the analytic stress numerically.
        close(axial_reaction(record), loading_response(p, epsilon).0, 2e-8);
    }
    assert_eq!(solution.increments[0].evolving_elements, 0);
    assert!(solution.increments.last().unwrap().evolving_elements > 0);
    let epsilon = PATH[PATH.len() - 1];
    let (sigma, alpha) = loading_response(p, epsilon);
    assert_fields(&solution, p, epsilon, sigma, alpha);
    solution
}

#[test]
fn data_only_yield_hardening_and_elastic_swaps_change_the_same_coupon() {
    let base = Parameters { young: 1000.0, poisson: 0.3, yield_stress: 1.0, hardening: 100.0 };
    let cases = [
        base,
        Parameters { yield_stress: 2.0, ..base },
        Parameters { hardening: 400.0, ..base },
        Parameters { hardening: 0.0, ..base },
        Parameters { young: 1500.0, poisson: 0.22, ..base },
    ];
    let basis = HardeningPropertyBasis::PlasticStrainModulus;
    let mut forces = Vec::new();
    for (index, &p) in cases.iter().enumerate() {
        let card = material_card(&format!("synthetic material {index}"), p, p, basis);
        let resolved = resolve(&card, 300.0, basis);
        close(resolved.hardening_modulus_pa(), p.hardening, 1e-12);
        let solution = assert_loading(&resolved, p);
        forces.push(axial_reaction(solution.increments.last().unwrap()));
    }
    assert!(forces[1] > forces[0]);
    assert!(forces[2] > forces[0]);
    assert!(forces[3] < forces[0]);
    assert!(forces[4] > forces[0]);
}

#[test]
fn sourced_mixed_hardening_changes_reverse_yield_on_the_original_tetrahedral_coupon() {
    let p = Parameters { young: 1000.0, poisson: 0.3, yield_stress: 1.0, hardening: 20.0 };
    let original = material_card("mixed-hardening coupon", p, p, HardeningPropertyBasis::PlasticStrainModulus);
    let mut claims = original.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS),
        value: PropertyValue::Scalar { value: 80.0, dims: Pressure::DIMS },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance { source: "authored cyclic Prager coupon".into(), license: "test-only".into(), artifact: None },
    }).unwrap();
    let source = MaterialCard::assemble(MaterialStateId {
        chemistry: "mixed-hardening coupon".into(), phase: "solid".into(),
        process: "synthetic cyclic coupon".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let mixed = resolve(&source, 300.0, HardeningPropertyBasis::MixedPlasticStrainModuli);
    let path = [0.0005, 0.006, 0.0038, -0.004, 0.002];
    let solution = solve(&mixed, &path);
    let mut plastic: f64 = 0.0;
    let mut accumulated: f64 = 0.0;
    for (&epsilon, record) in path.iter().zip(&solution.increments) {
        // Independent one-dimensional Prager/Hooke equations; the native
        // three-dimensional mesh must also solve its free lateral strains.
        let effective_trial = p.young * (epsilon - plastic) - 80.0 * plastic;
        let excess = effective_trial.abs() - (p.yield_stress + p.hardening * accumulated);
        if excess > 0.0 {
            let increment = excess / (p.young + p.hardening + 80.0);
            plastic += effective_trial.signum() * increment;
            accumulated += increment;
        }
        close(axial_reaction(record), p.young * (epsilon - plastic), 2e-8);
    }
    for state in &solution.element_states {
        close(state.plastic_strain[0], plastic, 1e-9);
        close(state.plastic_strain[1], -plastic / 2.0, 1e-9);
        close(state.alpha, accumulated, 1e-9);
    }
    let mono = Parameters { hardening: 100.0, ..p };
    let isotropic = resolve(&material_card("same monotonic tangent", mono, mono,
        HardeningPropertyBasis::PlasticStrainModulus), 300.0, HardeningPropertyBasis::PlasticStrainModulus);
    let other = solve(&isotropic, &path);
    close(axial_reaction(&solution.increments[1]), axial_reaction(&other.increments[1]), 2e-8);
    assert!((axial_reaction(&solution.increments[2]) - axial_reaction(&other.increments[2])).abs() > 1e-3);
}

#[test]
fn changing_only_temperature_requeries_all_plastic_and_elastic_properties() {
    let cold = Parameters { young: 1200.0, poisson: 0.28, yield_stress: 1.2, hardening: 120.0 };
    let hot = Parameters { young: 600.0, poisson: 0.34, yield_stress: 0.3, hardening: 20.0 };
    let midpoint = Parameters { young: 900.0, poisson: 0.31, yield_stress: 0.75, hardening: 70.0 };
    let basis = HardeningPropertyBasis::PlasticStrainModulus;
    let card = material_card("one temperature-dependent material", cold, hot, basis);
    let mut forces = Vec::new();
    let mut identities = Vec::new();
    for (temperature, p) in [(300.0, cold), (450.0, midpoint), (600.0, hot)] {
        let resolved = resolve(&card, temperature, basis);
        for (name, expected) in [
            (YOUNG_MODULUS_PROPERTY, p.young),
            (POISSON_RATIO_PROPERTY, p.poisson),
            (YIELD_STRESS_PROPERTY, p.yield_stress),
            (ISOTROPIC_HARDENING_MODULUS_PROPERTY, p.hardening),
        ] {
            close(resolved.material().property(name).unwrap().value_si(), expected, 1e-12);
        }
        identities.push(resolved.material().identity());
        let solution = assert_loading(&resolved, p);
        forces.push(axial_reaction(solution.increments.last().unwrap()));
    }
    assert!(forces.windows(2).all(|w| w[1] < w[0]), "{forces:?}");
    assert_ne!(identities[0], identities[1]);
    assert_ne!(identities[1], identities[2]);
}

#[test]
fn measured_total_strain_tangent_matches_plastic_hardening_through_unloading() {
    let p = Parameters { young: 1000.0, poisson: 0.3, yield_stress: 1.0, hardening: 100.0 };
    let h_basis = HardeningPropertyBasis::PlasticStrainModulus;
    let et_basis = HardeningPropertyBasis::UniaxialTangentModulus;
    let h_card = material_card("explicit plastic hardening", p, p, h_basis);
    let et_card = material_card("measured total-strain tangent", p, p, et_basis);
    let h = resolve(&h_card, 300.0, h_basis);
    let et = resolve(&et_card, 300.0, et_basis);
    close(h.hardening_modulus_pa(), p.hardening, 1e-12);
    close(et.hardening_modulus_pa(), p.hardening, 1e-12);
    let source_tangent = p.young * p.hardening / (p.young + p.hardening);
    close(et.material().property(POST_YIELD_TANGENT_MODULUS_PROPERTY).unwrap().value_si(),
        source_tangent, 1e-12);
    assert!(h.material().property(POST_YIELD_TANGENT_MODULUS_PROPERTY).is_none());
    assert!(et.material().property(ISOTROPIC_HARDENING_MODULUS_PROPERTY).is_none());

    let direct = assert_loading(&h, p);
    let measured = assert_loading(&et, p);
    for (a, b) in direct.increments.iter().zip(&measured.increments) {
        close(axial_reaction(a), axial_reaction(b), 2e-8);
    }

    // Returning to epsilon = alpha unloads elastically to zero force and
    // leaves the plastic axial strain and isochoric lateral contraction.
    let (_, alpha) = loading_response(p, PATH[PATH.len() - 1]);
    let mut unload_path = PATH.to_vec();
    unload_path.push(alpha);
    for resolved in [&h, &et] {
        let unloaded = solve(resolved, &unload_path);
        assert_eq!(unloaded.increments.len(), unload_path.len());
        close(axial_reaction(unloaded.increments.last().unwrap()), 0.0, 2e-8);
        assert_fields(&unloaded, p, alpha, 0.0, alpha);
    }
}
