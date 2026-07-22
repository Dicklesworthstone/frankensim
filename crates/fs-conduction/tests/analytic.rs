//! Analytic conduction solutions with DECLARED envelopes.
//!
//! Four families, and each states what kind of agreement it is claiming,
//! because they are not the same kind:
//!
//! | case | exact solution | why the envelope is what it is |
//! | --- | --- | --- |
//! | slab, Dirichlet–Dirichlet | linear in `x` | in the P₁ space, so agreement is ROUND-OFF-level and the envelope is `1e-9 K` |
//! | slab with a uniform source | quadratic | nodally reproduced by P₁ on this mesh, so again round-off at the nodes; the `L2` envelope is the interpolation error |
//! | slab, Dirichlet–Robin | linear in `x` | in the P₁ space; the envelope also pins the Robin heat rate against `k h ΔT/(k+h)` |
//! | annulus, radial log profile | `ln r` | NOT in the P₁ space: a DISCRETIZATION envelope, checked to shrink like `h²` under refinement. Its CONDUCTANCE envelope is looser and geometry-dominated: the annulus is meshed as a polygon, so the solved surface sits a chord sagitta inside the true cylinder |
//! | straight fin | 1-D fin equation | a MODEL comparison, not a discretization one: the envelope carries the fin model's own error, and the Biot number that bounds it is computed and printed |
//!
//! The fin row is the only one where "within envelope" includes a model
//! discrepancy the solver is not responsible for. That is stated here so
//! its number is never read as a discretization claim.

mod support;

use fs_conduction::assemble::{assemble_operator, full_residual};
use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{annulus_sector, box_grid, cylindrical_radius, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, ConductionSolution, InitialGuess, LinearConfig, Nonlinearity, SolveConfig,
    StopRule, element_heat_flux, solve,
};
use support::{l2_error, max_nodal_error, with_cx};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-conduction/analytic\",\"case\":\"{case}\",\
         \"verdict\":\"pass\",\"detail\":\"{}\"}}",
        support::json_escape(detail)
    );
}

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1e-11,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 12,
        },
        linear: LinearConfig {
            tolerance: 1e-13,
            max_iterations: 60_000,
            restart: 60,
        },
        initial: InitialGuess::DirichletMean,
    }
}

/// Nodal heat inflow at every vertex: the full-residual rows. On a
/// Dirichlet vertex this is the reaction — the heat entering the domain
/// through the prescribed row.
fn nodal_inflow(
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    source: &ScalarField,
    solution: &ConductionSolution,
) -> Vec<f64> {
    with_cx(|cx| {
        let system = assemble_operator(cx, mesh, boundary, material, source, &solution.temperature)
            .expect("assemble");
        full_residual(&system, &solution.temperature)
    })
}

// ------------------------------------------------------------ slab, D–D

#[test]
fn slab_dirichlet_dirichlet() {
    const K: f64 = 45.0;
    const T_HOT: f64 = 400.0;
    const T_COLD: f64 = 300.0;
    let (complex, positions) = box_grid([6, 3, 3], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(T_HOT).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(T_COLD).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            config(),
        )
        .expect("solve")
    });

    let exact = |p: [f64; 3]| T_HOT + (T_COLD - T_HOT) * p[0];
    let err = max_nodal_error(&mesh, &solution.temperature, &exact);
    assert!(
        err < 1e-9,
        "a linear profile lives in the P1 space; nodal error {err:e} must be round-off"
    );

    // Fourier's law through the slab: q_x = k ΔT / L, uniformly.
    let want_flux = K * (T_HOT - T_COLD);
    let flux = element_heat_flux(&mesh, &material, &solution.temperature).expect("flux");
    let worst_flux = flux
        .iter()
        .map(|q| (q[0] - want_flux).abs().max(q[1].abs()).max(q[2].abs()))
        .fold(0.0f64, f64::max);
    assert!(
        worst_flux < 1e-8 * want_flux,
        "recovered flux deviates by {worst_flux:e} from k ΔT/L = {want_flux:e}"
    );

    // The heat INTO the hot face must equal k A ΔT / L, and the heat out
    // of the cold face must match it: an independent check of the
    // Dirichlet reaction against the closed-form conductance.
    let inflow = nodal_inflow(&mesh, &boundary, &material, &source, &solution);
    let mut hot = 0.0f64;
    let mut cold = 0.0f64;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if on_box_face(p[0], 0.0) {
            hot += inflow[v];
        } else if on_box_face(p[0], 1.0) {
            cold += inflow[v];
        }
    }
    let want_q = K * (T_HOT - T_COLD);
    assert!(
        (hot - want_q).abs() < 1e-8 * want_q,
        "hot-face reaction {hot} != k A ΔT/L = {want_q}"
    );
    assert!(
        (cold + want_q).abs() < 1e-8 * want_q,
        "cold-face reaction {cold} != −k A ΔT/L"
    );
    verdict(
        "slab-dirichlet",
        &format!(
            "nodal_err={err:e} flux_err={worst_flux:e} Q_hot={hot} Q_analytic={want_q} \
             envelope=1e-9K/1e-8rel"
        ),
    );
}

// ----------------------------------------------------- slab with source

#[test]
fn slab_with_uniform_source() {
    const K: f64 = 20.0;
    const F: f64 = 4.0e4;
    const T_WALL: f64 = 300.0;
    let (complex, positions) = box_grid([8, 3, 3], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(F);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "walls",
            |f| on_box_face(f.centroid[0], 0.0) || on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(T_WALL).expect("bc"),
        )
        .expect("walls")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            config(),
        )
        .expect("solve")
    });

    // T(x) = T_wall + f x(1−x)/(2k); peak T_wall + f/(8k).
    let exact = |p: [f64; 3]| T_WALL + F * p[0] * (1.0 - p[0]) / (2.0 * K);
    let err = max_nodal_error(&mesh, &solution.temperature, &exact);
    assert!(
        err < 1e-8,
        "P1 reproduces a quadratic profile at the nodes on this mesh; got {err:e}"
    );
    let peak = solution
        .temperature
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let want_peak = F / (8.0 * K) + T_WALL;
    assert!((peak - want_peak).abs() < 1e-8);

    // The L2 envelope IS the interpolation error, because the nodal
    // values are exact: ‖T − I_h T‖ for this quadratic on h = 1/8. It is
    // stated RELATIVE to the peak temperature rise, because an absolute
    // kelvin envelope on an interpolation error is just a restatement of
    // the mesh size.
    let l2 = l2_error(&mesh, &solution.temperature, &exact);
    let rise = want_peak - T_WALL;
    assert!(
        l2 / rise < 2.0e-2,
        "L2 deviation {l2:e} K is {:.4} of the {rise} K rise, above the declared 2% \
         envelope for the h = 1/8 interpolation error",
        l2 / rise
    );

    // All the generated heat leaves through the two walls.
    let e = solution.report.energy;
    assert!(
        (e.source_w - F).abs() < 1e-6 * F,
        "the source integral {} should be f x volume = {F}",
        e.source_w
    );
    assert!(
        (e.dirichlet_in_w + F).abs() < 1e-6 * F,
        "steady state requires the walls to remove every watt generated"
    );
    verdict(
        "slab-with-source",
        &format!(
            "nodal_err={err:e} peak={peak} analytic_peak={want_peak} l2={l2:e} \
             l2_rel_rise={:.5} Q_gen={} Q_walls={} envelope=1e-8K nodal / 2% L2",
            l2 / rise,
            e.source_w,
            e.dirichlet_in_w
        ),
    );
}

// ------------------------------------------------------ slab, D–Robin

#[test]
fn slab_dirichlet_robin() {
    const K: f64 = 15.0;
    const H: f64 = 40.0;
    const T_HOT: f64 = 380.0;
    const T_INF: f64 = 295.0;
    let (complex, positions) = box_grid([6, 3, 3], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(T_HOT).expect("bc"),
        )
        .expect("hot")
        .region(
            "convective",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::robin(H, T_INF).expect("bc"),
        )
        .expect("convective")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            config(),
        )
        .expect("solve")
    });

    // −k C = h (T_hot + C − T_inf) ⇒ C = −h (T_hot − T_inf)/(k + h).
    let slope = -H * (T_HOT - T_INF) / (K + H);
    let exact = |p: [f64; 3]| T_HOT + slope * p[0];
    let err = max_nodal_error(&mesh, &solution.temperature, &exact);
    assert!(
        err < 1e-9,
        "the Dirichlet–Robin slab profile is linear, so P1 is exact; got {err:e}"
    );

    // The convective heat rate: k h ΔT/(k + h) per unit area.
    let want_q = K * H * (T_HOT - T_INF) / (K + H);
    let e = solution.report.energy;
    assert!(
        (e.robin_out_w - want_q).abs() < 1e-8 * want_q,
        "Robin heat rate {} != k h ΔT/(k+h) = {want_q}",
        e.robin_out_w
    );
    assert!(
        (e.dirichlet_in_w - want_q).abs() < 1e-8 * want_q,
        "the Dirichlet face must supply exactly what the convective face removes"
    );
    verdict(
        "slab-dirichlet-robin",
        &format!(
            "nodal_err={err:e} slope={slope} Q_robin={} Q_analytic={want_q} envelope=1e-9K",
            e.robin_out_w
        ),
    );
}

// ---------------------------------------------------- cylindrical shell

const R_IN: f64 = 0.5;
const R_OUT: f64 = 1.5;
const SWEEP: f64 = core::f64::consts::FRAC_PI_2;
const HEIGHT: f64 = 0.5;
const CYL_K: f64 = 12.0;
const T_IN: f64 = 420.0;
const T_OUT: f64 = 300.0;

fn cylinder_exact(p: [f64; 3]) -> f64 {
    let r = cylindrical_radius(p);
    T_IN + (T_OUT - T_IN) * fs_math::det::ln(r / R_IN) / fs_math::det::ln(R_OUT / R_IN)
}

fn run_cylinder(refine: usize) -> (f64, f64, f64) {
    let counts = [4 * refine, 6 * refine, 2 * refine];
    let (complex, positions) = annulus_sector(counts, R_IN, R_OUT, SWEEP, HEIGHT);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(CYL_K).expect("material");
    let source = ScalarField::Uniform(0.0);
    // Classify by the VERTEX radii, not the face centroid: the mesh is a
    // faceted approximation of the cylinder, so a chord's midpoint sits
    // inside the true radius while its vertices sit exactly on it.
    let radii: Vec<f64> = mesh
        .positions()
        .iter()
        .map(|&p| cylindrical_radius(p))
        .collect();
    let on_radius = |verts: [u32; 3], target: f64| {
        verts
            .iter()
            .all(|&v| (radii[v as usize] - target).abs() < 1e-9)
    };
    let radii_in = radii.clone();
    let radii_out = radii.clone();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "inner",
            |f| {
                f.vertices
                    .iter()
                    .all(|&v| (radii_in[v as usize] - R_IN).abs() < 1e-9)
            },
            ThermalBc::dirichlet(T_IN).expect("bc"),
        )
        .expect("inner")
        .region(
            "outer",
            |f| {
                f.vertices
                    .iter()
                    .all(|&v| (radii_out[v as usize] - R_OUT).abs() < 1e-9)
            },
            ThermalBc::dirichlet(T_OUT).expect("bc"),
        )
        .expect("outer")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    assert!(
        on_radius(mesh.boundary()[0].vertices, R_IN)
            || !on_radius(mesh.boundary()[0].vertices, R_IN),
        "radius classifier is total"
    );
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            config(),
        )
        .expect("solve")
    });
    let l2 = l2_error(&mesh, &solution.temperature, &cylinder_exact);
    let nodal = max_nodal_error(&mesh, &solution.temperature, &cylinder_exact);

    // Radial heat rate of the sector: k · height · sweep · ΔT / ln(r_o/r_i).
    let inflow = nodal_inflow(&mesh, &boundary, &material, &source, &solution);
    let mut q_in = 0.0f64;
    for (v, &r) in radii.iter().enumerate() {
        if (r - R_IN).abs() < 1e-9 {
            q_in += inflow[v];
        }
    }
    (l2, nodal, q_in)
}

#[test]
fn cylindrical_shell_radial_profile() {
    let want_q = CYL_K * HEIGHT * SWEEP * (T_IN - T_OUT) / fs_math::det::ln(R_OUT / R_IN);
    let (l2_coarse, _nodal_coarse, q_coarse) = run_cylinder(1);
    let (l2_fine, nodal_fine, q_fine) = run_cylinder(2);

    // ln r is NOT in the P1 space, so this is a DISCRETIZATION envelope:
    // it must shrink like h², i.e. roughly 4x per halving.
    let ratio = l2_coarse / l2_fine;
    assert!(
        ratio > 3.2,
        "L2 error ratio {ratio:.3} under a 2x refinement is not second order \
         ({l2_coarse:e} -> {l2_fine:e})"
    );
    let drop = T_IN - T_OUT;
    assert!(
        l2_fine / drop < 2.0e-3,
        "fine-grid L2 deviation {l2_fine:e} K is {:.5} of the {drop} K radial drop, \
         above the declared 0.2% envelope",
        l2_fine / drop
    );
    assert!(
        nodal_fine / drop < 1.0e-2,
        "fine-grid nodal deviation {nodal_fine:e} K is above the declared 1% envelope"
    );

    // The radial conductance is the classical log formula.
    let rel_coarse = (q_coarse - want_q).abs() / want_q;
    let rel_fine = (q_fine - want_q).abs() / want_q;
    // The 0.5% envelope on the CONDUCTANCE is dominated by GEOMETRY, not
    // by the PDE discretization: the annular boundary is meshed as a
    // polygon, so the solved domain's inner surface sits a chord sagitta
    // (≈ r Δθ²/8, here ≈ 0.2% of r) inside the true cylinder. Refinement
    // must shrink it, which is the assertion below.
    assert!(
        rel_fine < 5.0e-3,
        "fine-grid radial heat rate {q_fine} deviates {rel_fine:.4} from the analytic \
         {want_q} (envelope 0.5%, dominated by the polygonal boundary)"
    );
    assert!(
        rel_fine < rel_coarse,
        "refinement must improve the conductance: {rel_coarse:.5} -> {rel_fine:.5}"
    );
    verdict(
        "cylindrical-shell",
        &format!(
            "l2_coarse={l2_coarse:e} l2_fine={l2_fine:e} ratio={ratio:.3} \
             nodal_fine={nodal_fine:e} Q_fine={q_fine} Q_analytic={want_q} \
             rel={rel_fine:.5} envelopes=0.2%L2/1%nodal/0.5%Q(polygonal-boundary-dominated)"
        ),
    );
}

// ------------------------------------------------------------------ fin

#[test]
fn straight_fin_against_the_one_dimensional_model() {
    // Aluminium fin, forced-convection coefficient.
    const K: f64 = 200.0;
    const H: f64 = 25.0;
    const L: f64 = 0.05;
    const W: f64 = 0.02;
    const T: f64 = 0.002;
    const T_BASE: f64 = 350.0;
    const T_INF: f64 = 300.0;

    let (complex, positions) = box_grid([24, 6, 3], [L, W, T]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "base",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(T_BASE).expect("bc"),
        )
        .expect("base")
        .remainder("wetted", ThermalBc::robin(H, T_INF).expect("bc"))
        .expect("wetted")
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            config(),
        )
        .expect("solve")
    });

    // The 1-D fin model with a convective tip.
    let perimeter = 2.0 * (W + T);
    let area = W * T;
    let m = fs_math::det::sqrt(H * perimeter / (K * area));
    let hmk = H / (m * K);
    let ml = m * L;
    let cosh_ml = f64::midpoint(fs_math::det::exp(ml), fs_math::det::exp(-ml));
    let sinh_ml = (fs_math::det::exp(ml) - fs_math::det::exp(-ml)) / 2.0;
    let theta_b = T_BASE - T_INF;
    let denom = hmk.mul_add(sinh_ml, cosh_ml);
    let q_fin =
        fs_math::det::sqrt(H * perimeter * K * area) * theta_b * hmk.mul_add(cosh_ml, sinh_ml)
            / denom;
    // The Biot number that bounds the 1-D model's own error.
    let biot = H * (T / 2.0) / K;

    let q_solved = solution.report.energy.dirichlet_in_w;
    let rel = (q_solved - q_fin).abs() / q_fin;
    assert!(
        biot < 1.0e-3,
        "the 1-D fin model is only a fair comparison at small Biot; got {biot:e}"
    );
    assert!(
        rel < 2.0e-2,
        "3-D fin base heat rate {q_solved} W deviates {rel:.4} from the 1-D fin model \
         {q_fin} W (envelope 2%, which CARRIES the fin model's own error, not just \
         discretization)"
    );

    // Tip temperature from the same model.
    let tip_theta = theta_b / denom;
    let tip_numeric = mesh
        .positions()
        .iter()
        .zip(&solution.temperature)
        .filter(|(p, _)| on_box_face(p[0], L))
        .map(|(_, &t)| t)
        .fold(f64::NEG_INFINITY, f64::max);
    let tip_rel = ((tip_numeric - T_INF) - tip_theta).abs() / tip_theta;
    assert!(
        tip_rel < 2.0e-2,
        "tip excess temperature {} K deviates {tip_rel:.4} from the 1-D model {tip_theta} K",
        tip_numeric - T_INF
    );
    verdict(
        "straight-fin",
        &format!(
            "Bi={biot:e} mL={ml:.4} Q_3d={q_solved:.5}W Q_1d={q_fin:.5}W rel={rel:.5} \
             tip_3d={:.4}K tip_1d={tip_theta:.4}K envelope=2%(model+discretization)",
            tip_numeric - T_INF
        ),
    );
}
