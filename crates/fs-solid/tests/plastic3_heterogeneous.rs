//! Connected material regions exercise the native three-dimensional solve.
//! Synthetic cards isolate constitutive data changes; no material calibration
//! or mesh-convergence claim is made by these small affine fixtures.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::SmallStrainLaw;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, Pressure};
use fs_solid::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    TetAssemblyBudget,
};

#[derive(Clone, Copy)]
struct Parameters {
    young: f64,
    poisson: f64,
    yield_stress: f64,
    hardening: f64,
}

fn material(name: &str, p: Parameters) -> ResolvedJ2Plasticity {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, p.young),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), p.poisson),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, p.yield_stress),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, p.hardening),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic bonded-material regression".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "synthetic".into(),
        revision: 0,
    }, claims, vec![]).unwrap();
    let resolved = resolve_j2_plasticity(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::PlasticStrainModulus,
        0.05,
    ).unwrap();
    assert_eq!(resolved.material().card_identity(), card.content_hash());
    resolved
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 5, kernel_id: 23, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        f(&cx)
    })
}

// Two bonded unit cubes share the x=1 face and its actual mesh vertices.
// Each cube has six positively oriented tetrahedra; each region volume is 1 m3.
fn mesh() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let mut nodes = Vec::new();
    for z in 0..=1 {
        for y in 0..=1 {
            for x in 0..=2 {
                nodes.push([f64::from(x), f64::from(y), f64::from(z)]);
            }
        }
    }
    let mut tets = Vec::new();
    for x in 0..2 {
        let c = [x, x + 1, x + 3, x + 4, x + 6, x + 7, x + 9, x + 10];
        for tet in [
            [0, 1, 3, 7], [0, 1, 7, 5], [0, 2, 7, 3],
            [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 7, 6],
        ] {
            tets.push(tet.map(|i| c[i]));
        }
    }
    (nodes, tets)
}

fn field(left: &ResolvedJ2Plasticity, right: &ResolvedJ2Plasticity) -> Vec<ResolvedJ2Plasticity> {
    (0..12).map(|e| if e < 6 { left.clone() } else { right.clone() }).collect()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    laws: &'a [ResolvedJ2Plasticity],
    prescribed: &'a [(usize, f64)],
    path: &'a [f64],
) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        laws: SmallStrainLawField::PerElement(laws),
        prescribed_m: prescribed,
        nodal_forces_n: &[],
        load_path: path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:.15e} vs {expected:.15e}, tolerance {tolerance:e}");
}

#[test]
fn bonded_regions_follow_series_compliance_after_a_data_only_modulus_swap() {
    let (nodes, tets) = mesh();
    let left_p = Parameters { young: 1000.0, poisson: 0.0, yield_stress: 100.0, hardening: 0.0 };
    let left = material("left region", left_p);
    let mut prescribed = Vec::new();
    for (i, x) in nodes.iter().enumerate() {
        if x[0] == 0.0 { prescribed.push((3 * i, 0.0)); }
        if x[0] == 2.0 { prescribed.push((3 * i, 1.0)); }
        if x[1] == 0.0 { prescribed.push((3 * i + 1, 0.0)); }
        if x[2] == 0.0 { prescribed.push((3 * i + 2, 0.0)); }
    }
    let path = [0.001, 0.003];
    let mut reactions = Vec::new();
    for young in [1000.0, 4000.0] {
        let right = material("right region", Parameters { young, ..left_p });
        let laws = field(&left, &right);
        assert_ne!(laws[0].material().card_identity(), laws[6].material().card_identity());
        let p = problem(&nodes, &tets, &laws, &prescribed, &path);
        let result = with_cx(|cx| p.solve(cx)).unwrap();
        for (record, &extension) in result.increments.iter().zip(&path) {
            // Each series member has L=A=1. nu=0 gives an exact piecewise
            // affine displacement field with free lateral motion.
            let expected = extension / (1.0 / left_p.young + 1.0 / young);
            let force: f64 = prescribed.iter().zip(&record.reactions_n)
                .filter(|((dof, value), _)| dof % 3 == 0 && *value == 1.0)
                .map(|(_, force)| force).sum();
            close(force, expected, 1e-8);
        }
        let force = path[1] / (1.0 / left_p.young + 1.0 / young);
        for (x, u) in nodes.iter().zip(&result.displacement_m) {
            let expected = if x[0] <= 1.0 {
                force * x[0] / left_p.young
            } else {
                force / left_p.young + force * (x[0] - 1.0) / young
            };
            close(u[0], expected, 1e-10);
            close(u[1], 0.0, 1e-10);
            close(u[2], 0.0, 1e-10);
        }
        for stress in &result.element_stress {
            close(stress[0], force, 1e-8);
            for &value in &stress[1..] { close(value, 0.0, 1e-8); }
        }
        if young == left_p.young {
            let uniform = SmallStrainTetProblem {
                laws: SmallStrainLawField::Uniform(&left),
                ..p
            };
            assert_eq!(result, with_cx(|cx| uniform.solve(cx)).unwrap());
        }
        reactions.push(force);
    }
    assert!(reactions[1] > reactions[0]);
}

// Uniaxial STRAIN, with constrained lateral motion: q=2*mu*eps-3*mu*alpha,
// sigma_xx=K*eps+2*q/3. Monotonic peak strain gives alpha in closed form.
// The chosen short unloading leg stays inside the expanded yield surface.
fn constrained_response(p: Parameters, strain: f64, peak: f64) -> (f64, f64, f64) {
    let mu = p.young / (2.0 * (1.0 + p.poisson));
    let bulk = p.young / (3.0 * (1.0 - 2.0 * p.poisson));
    let alpha = ((2.0 * mu * peak - p.yield_stress) / (3.0 * mu + p.hardening)).max(0.0);
    let q = 2.0 * mu * strain - 3.0 * mu * alpha;
    assert!(q.abs() <= p.yield_stress + p.hardening * alpha + 1e-12);
    (bulk * strain + 2.0 * q / 3.0, bulk * strain - q / 3.0, alpha)
}

#[test]
fn material_data_move_the_yielding_region_and_preserve_unloading_history() {
    let (nodes, tets) = mesh();
    let left_p = Parameters { young: 1000.0, poisson: 0.25, yield_stress: 1.0, hardening: 100.0 };
    let right_p = Parameters { young: 1400.0, poisson: 0.3, yield_stress: 3.0, hardening: 200.0 };
    let left = material("soft region", left_p);
    let prescribed: Vec<_> = nodes.iter().enumerate()
        .flat_map(|(i, x)| [(3 * i, x[0]), (3 * i + 1, 0.0), (3 * i + 2, 0.0)])
        .collect();
    let path = [0.0005, 0.002, 0.004, 0.0035];
    for yield_stress in [3.0, 4.5] {
        let changed = Parameters { yield_stress, ..right_p };
        let right = material("replaceable region", changed);
        let laws = field(&left, &right);
        let result = with_cx(|cx| problem(&nodes, &tets, &laws, &prescribed, &path).solve(cx)).unwrap();
        let mut peak: f64 = 0.0;
        for (record, &strain) in result.increments.iter().zip(&path) {
            peak = peak.max(strain);
            let expected = constrained_response(left_p, strain, peak).0
                + constrained_response(changed, strain, peak).0;
            // Virtual work includes all imposed interface reactions; each
            // region contributes volume*sigma_xx for unit axial strain.
            let work: f64 = prescribed.iter().zip(&record.reactions_n)
                .map(|((_, displacement), force)| displacement * force).sum();
            close(work, expected, 1e-9);
        }
        assert_eq!(result.increments[0].evolving_elements, 0);
        assert_eq!(result.increments[1].evolving_elements, 6);
        assert_eq!(result.increments[2].evolving_elements, if yield_stress == 3.0 { 12 } else { 6 });
        assert_eq!(result.increments[3].evolving_elements, 0);
        for (e, state) in result.element_states.iter().enumerate() {
            let p = if e < 6 { left_p } else { changed };
            let (axial, lateral, alpha) = constrained_response(p, path[3], path[2]);
            close(state.alpha, alpha, 1e-12);
            for (actual, expected) in state.plastic_strain.iter().zip([alpha, -alpha / 2.0, -alpha / 2.0, 0.0, 0.0, 0.0]) {
                close(*actual, expected, 1e-12);
            }
            for (actual, expected) in result.element_stress[e].iter().zip([axial, lateral, lateral, 0.0, 0.0, 0.0]) {
                close(*actual, expected, 1e-9);
            }
        }
    }
}

#[test]
fn heterogeneous_force_and_tangent_preserve_tensor_shear_work() {
    let (nodes, tets) = mesh();
    let a = Parameters { young: 1000.0, poisson: 0.25, yield_stress: 1.0, hardening: 100.0 };
    let b = Parameters { young: 1400.0, poisson: 0.3, yield_stress: 3.0, hardening: 200.0 };
    let laws = field(&material("a", a), &material("b", b));
    let p = problem(&nodes, &tets, &laws, &[], &[1.0]);
    let virgin: Vec<_> = laws.iter().map(SmallStrainLaw::initial_state).collect();
    let gamma = 1e-5;
    let shear: Vec<_> = nodes.iter().flat_map(|x| [2.0 * gamma * x[1], 0.0, 0.0]).collect();
    let virtual_shear: Vec<_> = nodes.iter().flat_map(|x| [2.0 * x[1], 0.0, 0.0]).collect();
    let force = with_cx(|cx| p.internal_force(&shear, &virgin, cx)).unwrap();
    let work: f64 = force.iter().zip(&virtual_shear).map(|(f, d)| f * d).sum();
    let sum_mu = a.young / (2.0 * (1.0 + a.poisson)) + b.young / (2.0 * (1.0 + b.poisson));
    close(work, 4.0 * gamma * sum_mu, 1e-12);

    let states: Vec<_> = laws.iter().map(|law| {
        law.update_state(&[0.0015, -0.00075, -0.00075, 0.0, 0.0, 0.0], &law.initial_state())
    }).collect();
    assert!(states[0].alpha > 0.0);
    assert_eq!(states[6].alpha, 0.0);
    let u: Vec<_> = nodes.iter().flat_map(|x| [
        0.005 * x[0] + 0.002 * x[1] * x[2],
        -0.0007 * x[1],
        0.0015 * x[0] * x[2],
    ]).collect();
    let d: Vec<_> = (0..u.len()).map(|i| (i * 37 % 101) as f64 / 101.0 - 0.5).collect();
    let delta = 1e-8;
    let up: Vec<_> = u.iter().zip(&d).map(|(u, d)| u + delta * d).collect();
    let um: Vec<_> = u.iter().zip(&d).map(|(u, d)| u - delta * d).collect();
    let (k, fp, fm) = with_cx(|cx| (
        p.tangent_matrix(&u, &states, cx).unwrap(),
        p.internal_force(&up, &states, cx).unwrap(),
        p.internal_force(&um, &states, cx).unwrap(),
    ));
    let mut kd = vec![0.0; u.len()];
    k.spmv(&d, &mut kd);
    let error: f64 = kd.iter().zip(fp.iter().zip(&fm))
        .map(|(kd, (fp, fm))| (kd - (fp - fm) / (2.0 * delta)).powi(2)).sum();
    let scale: f64 = kd.iter().map(|v| v * v).sum();
    assert!(scale > 0.0);
    assert!((error / scale).sqrt() < 1e-6);
}

#[test]
fn wrong_law_counts_refuse_every_native_entry_point_before_indexing() {
    let (nodes, tets) = mesh();
    let law = material("one law", Parameters {
        young: 1000.0, poisson: 0.25, yield_stress: 1.0, hardening: 100.0,
    });
    let u = vec![0.0; 3 * nodes.len()];
    let states = vec![law.initial_state(); tets.len()];
    for count in [0, 11, 13] {
        let laws = vec![law.clone(); count];
        let p = problem(&nodes, &tets, &laws, &[], &[1.0]);
        let expected = SmallStrainTetError::LawCountMismatch { expected: 12, got: count };
        with_cx(|cx| {
            assert_eq!(p.solve(cx).unwrap_err(), expected);
            assert_eq!(p.internal_force(&u, &states, cx).unwrap_err(), expected);
            assert_eq!(p.tangent_matrix(&u, &states, cx).unwrap_err(), expected);
        });
    }
}
