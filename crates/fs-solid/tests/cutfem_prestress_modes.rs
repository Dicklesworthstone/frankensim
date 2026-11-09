//! G0/G1/G3: signed native prestress modes, physical mass and sparse whitening.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::mass::MaterialMass3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, PreparedPrestress3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_exec::CancelGate;
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solid::cutfront::dynamics3::CutElasticityRef3;
use fs_solid::cutfront::prestress3::{
    CutPrestressError3, CutPrestressLimits3, CutPrestressModeOptions3, CutPrestressModeReport3,
    CutPrestressModes3,
};
use fs_solver::LinearOp;

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn material() -> ElasticAcousticState {
    let diagonal = [12.0, 14.0, 16.0, 6.0, 7.0, 8.0];
    let tensor =
        std::array::from_fn(|i| std::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 }));
    ElasticAcousticState::try_new_mandel(tensor, 3.0, ContentHash([83; 32]), 0.1).unwrap()
}
fn control<T>(action: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    action(&mut control)
}
fn cartesian(clamped: bool) -> CutElasticity3 {
    control(|control| {
        CutElasticity3::build_with_material_field(
            HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
            [1; 3],
            &Top,
            ElasticMaterialField3::Uniform(&material()),
            &|p| clamped && p[0] == 0.0,
            Default::default(),
            control,
        )
        .unwrap()
    })
}
fn adaptive() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    control(|control| {
        AdaptiveElasticity3::build_with_material_field(
            HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
            &tree,
            &Top,
            ElasticMaterialField3::Uniform(&material()),
            &|p| p[0] == 0.0,
            Default::default(),
            control,
        )
        .unwrap()
    })
}
fn prestress(source: CutElasticityRef3<'_>, expansion: f64) -> PreparedPrestress3 {
    let field = |_, _| Ok([expansion, expansion, expansion, 0.0, 0.0, 0.0]);
    // q=0 is an equilibrium under its exact balancing dead force -b_free.
    // Positive prescribed expansion therefore gives compressive initial stress.
    match source {
        CutElasticityRef3::Cartesian(op) => {
            let free = op.prepare_eigenstrain(&field, 100_000, poll).unwrap();
            op.prepare_eigenstrain_prestress(
                &free,
                &vec![0.0; op.n()],
                0.05,
                100_000,
                100_000,
                poll,
            )
            .unwrap()
        }
        CutElasticityRef3::Adaptive(op) => {
            let free = op.prepare_eigenstrain(&field, 100_000, poll).unwrap();
            op.prepare_eigenstrain_prestress(
                &free,
                &vec![0.0; op.n()],
                0.05,
                100_000,
                100_000,
                poll,
            )
            .unwrap()
        }
    }
}
fn options(modes: usize) -> CutPrestressModeOptions3 {
    CutPrestressModeOptions3 {
        modes,
        max_ticks: 768,
        spectral_tolerance: 2e-8,
        residual_absolute_tolerance: 2e-8,
        ..Default::default()
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

// Independent bounded test oracle only. Production never builds these dense
// matrices. Native coordinate actions and a DENSE Cholesky reduction avoid
// sharing the adapter's scalar-pair assembly, permutation or LDL whitening.
fn dense_values(
    op: &dyn LinearOp,
    geometric: &PreparedPrestress3,
    mass: &MaterialMass3,
    occupancy: &[f64],
    free: &[usize],
) -> Vec<f64> {
    let n = free.len();
    assert!(n <= 256);
    let mut m = vec![0.0; n * n];
    let mut h = m.clone();
    let mut unit = vec![0.0; op.n()];
    let mut force = unit.clone();
    let mut g = unit.clone();
    for (column, &dof) in free.iter().enumerate() {
        unit.fill(0.0);
        unit[dof] = 1.0;
        mass.apply_into(occupancy, &unit, &mut force, poll).unwrap();
        for (row, &native) in free.iter().enumerate() {
            m[row * n + column] = force[native];
        }
        op.apply(&unit, &mut force);
        geometric.apply_into(&unit, &mut g, poll).unwrap();
        for (row, &native) in free.iter().enumerate() {
            h[row * n + column] = force[native] + g[native];
        }
    }
    let cholesky = fs_la::factor::cholesky(&m, n).unwrap();
    let mut s = vec![0.0; n * n];
    for column in 0..n {
        for row in (0..n).rev() {
            let mut value = if row == column { 1.0 } else { 0.0 };
            for k in row + 1..n {
                value -= cholesky.l(k, row) * s[k * n + column];
            }
            s[row * n + column] = value / cholesky.l(row, row);
        }
    }
    let mut hs = vec![0.0; n * n];
    let mut transformed = hs.clone();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                hs[i * n + j] += h[i * n + k] * s[k * n + j];
            }
        }
    }
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                transformed[i * n + j] += s[k * n + i] * hs[k * n + j];
            }
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            let value = f64::midpoint(transformed[i * n + j], transformed[j * n + i]);
            transformed[i * n + j] = value;
            transformed[j * n + i] = value;
        }
    }
    fs_la::eigen::householder_eigh(&transformed, n).0
}
fn check_physical(
    report: &CutPrestressModeReport3,
    op: &dyn LinearOp,
    geometric: &PreparedPrestress3,
    mass: &MaterialMass3,
    occupancy: &[f64],
    fixed: &[bool],
) {
    assert!(report.maximum_mass_orthogonality_error <= report.options.mass_orthogonality_tolerance);
    for (i, mode) in report.modes.iter().enumerate() {
        let mut inertia = vec![0.0; op.n()];
        let mut elastic = inertia.clone();
        let mut g = inertia.clone();
        mass.apply_into(occupancy, &mode.vector, &mut inertia, poll)
            .unwrap();
        op.apply(&mode.vector, &mut elastic);
        geometric.apply_into(&mode.vector, &mut g, poll).unwrap();
        close(dot(&mode.vector, &inertia), 1.0, 1e-10, 1e-10);
        for previous in &report.modes[..i] {
            assert!(dot(&previous.vector, &inertia).abs() < 1e-8);
        }
        let residual = elastic
            .iter()
            .zip(&g)
            .zip(&inertia)
            .fold(0.0_f64, |norm, ((k, g), m)| {
                norm.hypot(k + g - mode.eigenvalue * m)
            });
        close(residual, mode.native_residual_norm, 1e-5, 2e-12);
        assert!(
            mode.whitened_residual_norm
                <= report.options.residual_absolute_tolerance
                    + report.options.residual_relative_tolerance * mode.whitened_residual_scale
        );
        for (node, fixed) in fixed.iter().enumerate() {
            if *fixed {
                assert_eq!(&mode.vector[3 * node..3 * node + 3], &[0.0; 3]);
            }
        }
    }
}

#[test]
fn tensile_and_compressive_cut_pencils_match_independent_dense_generalized_modes() {
    let op = cartesian(true);
    let source = CutElasticityRef3::Cartesian(&op);
    let mass = op.prepare_material_mass(100_000, poll).unwrap();
    let occupancy = vec![0.7; op.cells()];
    let gate = CancelGate::new_clock_free();
    let mut spectra = Vec::new();
    for expansion in [0.03, 0.0, -0.03] {
        let g = prestress(source, expansion);
        let owner =
            CutPrestressModes3::new(source, &g, &mass, &occupancy, Default::default(), &gate)
                .unwrap();
        let report = owner.lowest_modes(options(3)).unwrap();
        let reference = dense_values(&op, &g, &mass, &occupancy, owner.free_dofs());
        for (actual, expected) in report.modes.iter().zip(&reference) {
            close(actual.eigenvalue, *expected, 2e-8, 2e-8);
        }
        check_physical(&report, &op, &g, &mass, &occupancy, op.fixed());
        spectra.push(
            report
                .modes
                .iter()
                .map(|m| m.eigenvalue)
                .collect::<Vec<_>>(),
        );
    }
    for i in 0..3 {
        assert!(spectra[0][i] < spectra[1][i]);
        assert!(spectra[1][i] < spectra[2][i]);
    }
}

#[test]
fn negative_and_zero_values_survive_without_frequency_conversion() {
    let op = cartesian(false);
    let source = CutElasticityRef3::Cartesian(&op);
    let mass = op.prepare_material_mass(100_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let gate = CancelGate::new_clock_free();
    for expansion in [0.02, 0.0] {
        let g = prestress(source, expansion);
        let owner =
            CutPrestressModes3::new(source, &g, &mass, &occupancy, Default::default(), &gate)
                .unwrap();
        let report = owner.lowest_modes(options(3)).unwrap();
        let reference = dense_values(&op, &g, &mass, &occupancy, owner.free_dofs());
        for (actual, expected) in report.modes.iter().zip(&reference) {
            close(actual.eigenvalue, *expected, 2e-8, 3e-8);
            if expansion > 0.0 {
                assert!(actual.eigenvalue < -1e-4);
            } else {
                assert!(actual.eigenvalue.abs() < 3e-8);
            }
        }
        check_physical(&report, &op, &g, &mass, &occupancy, op.fixed());
    }
}

#[test]
fn adaptive_condensation_and_nonuniform_physical_mass_survive_sparse_whitening() {
    let mut op = adaptive();
    assert!(op.physical_nodes().len() > op.nodes().len());
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.6 + 0.01 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let source = CutElasticityRef3::Adaptive(&op);
    let mass = op.prepare_material_mass(100_000, poll).unwrap();
    let occupancy: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.003 * i as f64).collect();
    let g = prestress(source, 0.003);
    let gate = CancelGate::new_clock_free();
    let owner =
        CutPrestressModes3::new(source, &g, &mass, &occupancy, Default::default(), &gate).unwrap();
    let report = owner.lowest_modes(options(2)).unwrap();
    let reference = dense_values(&op, &g, &mass, &occupancy, owner.free_dofs());
    for (actual, expected) in report.modes.iter().zip(&reference) {
        close(actual.eigenvalue, *expected, 2e-7, 2e-7);
    }
    check_physical(&report, &op, &g, &mass, &occupancy, op.fixed());
    let half: Vec<_> = occupancy.iter().map(|r| 0.5 * r).collect();
    let lighter =
        CutPrestressModes3::new(source, &g, &mass, &half, Default::default(), &gate).unwrap();
    let shifted = lighter.lowest_modes(options(2)).unwrap();
    for (light, full) in shifted.modes.iter().zip(&report.modes) {
        close(light.eigenvalue, 2.0 * full.eigenvalue, 2e-7, 2e-7);
    }
    assert_eq!(op.scales(), scales);
}

#[test]
fn stale_sources_mass_and_resource_refusals_do_not_publish_partial_modes() {
    let mut op = cartesian(true);
    let other = cartesian(true);
    let mass = op.prepare_material_mass(100_000, poll).unwrap();
    let foreign_mass = other.prepare_material_mass(100_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let source = CutElasticityRef3::Cartesian(&op);
    let g = prestress(source, 0.01);
    let foreign = prestress(CutElasticityRef3::Cartesian(&other), 0.01);
    let gate = CancelGate::new_clock_free();
    assert!(
        CutPrestressModes3::new(
            source,
            &foreign,
            &mass,
            &occupancy,
            Default::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        CutPrestressModes3::new(
            source,
            &g,
            &foreign_mass,
            &occupancy,
            Default::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        CutPrestressModes3::new(
            source,
            &g,
            &mass,
            &vec![0.0; op.cells()],
            Default::default(),
            &gate
        )
        .is_err()
    );
    let mut limits = CutPrestressLimits3::default();
    limits.max_mass_entries = 1;
    assert!(matches!(
        CutPrestressModes3::new(source, &g, &mass, &occupancy, limits, &gate),
        Err(CutPrestressError3::Budget(_))
    ));
    let mut limits = CutPrestressLimits3::default();
    limits.factor.max_symbolic_entries = 1;
    assert!(matches!(
        CutPrestressModes3::new(source, &g, &mass, &occupancy, limits, &gate),
        Err(CutPrestressError3::Factor(_))
    ));
    let limits = CutPrestressLimits3 {
        max_operator_applications: 1,
        ..Default::default()
    };
    let bounded = CutPrestressModes3::new(source, &g, &mass, &occupancy, limits, &gate).unwrap();
    assert!(matches!(
        bounded.lowest_modes(options(1)),
        Err(CutPrestressError3::Budget(_))
    ));
    let owner =
        CutPrestressModes3::new(source, &g, &mass, &occupancy, Default::default(), &gate).unwrap();
    let report = owner.lowest_modes(options(1)).unwrap();
    let repeated = owner.lowest_modes(options(1)).unwrap();
    assert_eq!(
        report.modes[0].eigenvalue.to_bits(),
        repeated.modes[0].eigenvalue.to_bits()
    );
    assert_eq!(report.modes[0].vector, repeated.modes[0].vector);
    gate.request();
    assert!(matches!(
        owner.lowest_modes(options(1)),
        Err(CutPrestressError3::Cancelled)
    ));
    drop(owner);
    drop(bounded);
    op.set_scales(&[0.5]).unwrap();
    let gate = CancelGate::new_clock_free();
    assert!(matches!(
        CutPrestressModes3::new(
            CutElasticityRef3::Cartesian(&op),
            &g,
            &mass,
            &occupancy,
            Default::default(),
            &gate
        ),
        Err(CutPrestressError3::Physics(_))
    ));
}
