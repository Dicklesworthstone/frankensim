use super::super::tests::{cartesian_field, reference, settings};
use super::*;

const SHAPE: [usize; 3] = [5, 6, 7];
const SPACING: [f64; 3] = [0.3, 0.7, 1.1];

fn config(averaging: DynamicAveraging) -> DynamicComponentsConfig {
    DynamicComponentsConfig {
        averaging,
        maximum_coefficient: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn fit(averaging: DynamicAveraging, components: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (rho, u, heat) = cartesian_field(SHAPE);
    let (mc, hc) = settings(averaging);
    let mut flow = PeriodicDynamicSgs::new_cartesian(SHAPE, SPACING, usize::MAX).unwrap();
    let mut thermal = PeriodicDynamicScalar::new(SHAPE, usize::MAX).unwrap();
    let mut group = PeriodicDynamicComponents::new(SHAPE, components.len(), usize::MAX).unwrap();
    flow.prepare_with_scalar_and_components(
        mc,
        hc,
        config(averaging),
        &rho,
        &u,
        &heat,
        &mut thermal,
        |component, cell| components[component][cell],
        &mut group,
        || false,
    )
    .unwrap();
    (
        flow.viscosities().to_vec(),
        thermal.diffusivities().to_vec(),
        group.diffusivities().to_vec(),
    )
}

fn close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
        assert!((a - b).abs() < tolerance, "cell {i}: {a} != {b}");
    }
}

#[test]
fn g1_one_component_matches_physical_reference_and_preserves_original_heat_and_flow_bits() {
    let (rho, u, heat) = cartesian_field(SHAPE);
    for averaging in [DynamicAveraging::LocalBox, DynamicAveraging::Domain] {
        let (nu, alpha, d) = fit(averaging, &[heat.clone()]);
        let (mc, hc) = settings(averaging);
        let mut flow = PeriodicDynamicSgs::new_cartesian(SHAPE, SPACING, usize::MAX).unwrap();
        let mut thermal = PeriodicDynamicScalar::new(SHAPE, usize::MAX).unwrap();
        flow.prepare_with_scalar(mc, hc, &rho, &u, &heat, &mut thermal, || false)
            .unwrap();
        assert_eq!(nu, flow.viscosities());
        assert_eq!(alpha, thermal.diffusivities());
        assert_eq!(d, alpha);
        let (_, expected) = reference(SHAPE, SPACING, &rho, &u, &heat, averaging);
        close(&d, &expected, 2e-13);
        // All components are same-unit affine copies, so the least-squares
        // numerator and denominator acquire the same sum of squared slopes.
        let collinear = vec![
            heat.clone(),
            heat.iter().map(|v| 2.0 - 3.0 * v).collect(),
            vec![7.0; heat.len()],
        ];
        close(&fit(averaging, &collinear).2, &d, 2e-13);
    }
}

#[test]
fn g3_same_unit_component_basis_and_species_order_do_not_change_common_diffusion() {
    let (_, _, phi) = cartesian_field(SHAPE);
    let first: Vec<_> = phi.iter().map(|v| 0.25 + 0.03 * v).collect();
    let second: Vec<_> = (0..phi.len())
        .map(|i| 0.35 + 0.02 * phi[(i + 17) % phi.len()])
        .collect();
    let third: Vec<_> = first
        .iter()
        .zip(&second)
        .map(|(a, b)| 1.0 - a - b)
        .collect();
    let fields = vec![first, second, third, vec![0.0; phi.len()]];
    let permuted = vec![
        fields[3].clone(),
        fields[2].clone(),
        fields[0].clone(),
        fields[1].clone(),
    ];
    let rotated = vec![
        fields[0]
            .iter()
            .zip(&fields[1])
            .map(|(a, b)| 0.6 * a + 0.8 * b)
            .collect(),
        fields[0]
            .iter()
            .zip(&fields[1])
            .map(|(a, b)| -0.8 * a + 0.6 * b)
            .collect(),
        fields[2].clone(),
        fields[3].clone(),
    ];
    let scaled: Vec<Vec<f64>> = fields
        .iter()
        .enumerate()
        .map(|(s, field)| field.iter().map(|v| -12.5 * v + s as f64).collect())
        .collect();
    for averaging in [DynamicAveraging::LocalBox, DynamicAveraging::Domain] {
        let reference = fit(averaging, &fields);
        for alternative in [&permuted, &rotated, &scaled] {
            let result = fit(averaging, alternative);
            assert_eq!(result.0, reference.0);
            assert_eq!(result.1, reference.1);
            close(&result.2, &reference.2, 3e-12);
        }
    }
}

#[test]
fn g3_trace_component_keeps_its_physical_weight_in_the_joint_objective() {
    let (_, _, phi) = cartesian_field(SHAPE);
    let other: Vec<_> = (0..phi.len()).map(|i| phi[(i + 17) % phi.len()]).collect();
    let averaging = DynamicAveraging::LocalBox;
    let reference = fit(averaging, &[phi.clone()]).2;
    let equal = fit(averaging, &[phi.clone(), other.clone()]).2;
    assert!(
        equal
            .iter()
            .zip(&reference)
            .any(|(a, b)| (a - b).abs() > 1e-5)
    );
    let trace: Vec<_> = other.iter().map(|v| 1e-6 * v).collect();
    let result = fit(averaging, &[phi.clone(), trace, vec![0.0; phi.len()]]).2;
    // A per-component normalization bug would reproduce the equal-weight
    // result above; the correct perturbation is of order 1e-12.
    close(&result, &reference, 2e-11);
}

#[test]
fn g0_independent_species_diffusion_reuses_buffers_and_refuses_incomplete_trials() {
    let n = 9;
    let cells = n * n * n;
    let rho = vec![1.0; cells];
    let heat = vec![300.0; cells];
    let mut u = vec![[0.0; 3]; cells];
    let mut y = vec![0.0; cells];
    for cell in 0..cells {
        let x = core::f64::consts::TAU * (cell % n) as f64 / n as f64;
        let z = core::f64::consts::TAU * (cell / n % n) as f64 / n as f64;
        u[cell][0] = 0.3 * z.sin();
        y[cell] = 0.5 + 0.1 * (x + z).sin();
    }
    let values = |s, cell| match s {
        0 => y[cell],
        1 => 1.0 - y[cell],
        _ => 0.0,
    };
    let averaging = DynamicAveraging::LocalBox;
    let (mc, hc) = settings(averaging);
    let cc = config(averaging);
    let mut flow = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    let mut thermal = PeriodicDynamicScalar::new([n; 3], usize::MAX).unwrap();
    let bytes = PeriodicDynamicComponents::required_bytes(cells).unwrap();
    assert_eq!(bytes, 80 * cells);
    assert!(matches!(
        PeriodicDynamicComponents::new([n; 3], 3, bytes - 1),
        Err(DynamicSgsError::WorkspaceBytes { .. })
    ));
    assert!(PeriodicDynamicComponents::new([n; 3], 0, bytes).is_err());
    assert!(PeriodicDynamicComponents::new([n; 3], usize::MAX, bytes).is_err());
    let mut group = PeriodicDynamicComponents::new([n; 3], 3, bytes).unwrap();
    assert_eq!(group.component_count(), 3);
    assert_eq!(group.workspace_bytes(), bytes);
    let allocation = group.diffusivities().as_ptr();
    let mut polls = 0;
    let reports = flow
        .prepare_with_scalar_and_components(
            mc,
            hc,
            cc,
            &rho,
            &u,
            &heat,
            &mut thermal,
            values,
            &mut group,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(reports.0.maximum_eddy_viscosity, 0.0);
    assert_eq!(reports.1.maximum_diffusivity, 0.0);
    assert!(reports.2.maximum_diffusivity > 1e-4);
    let expected = group.diffusivities().to_vec();
    for cancel_at in [1, polls - 1] {
        let mut seen = 0;
        assert_eq!(
            flow.prepare_with_scalar_and_components(
                mc,
                hc,
                cc,
                &rho,
                &u,
                &heat,
                &mut thermal,
                values,
                &mut group,
                || {
                    seen += 1;
                    seen == cancel_at
                },
            ),
            Err(DynamicSgsError::Cancelled)
        );
        let retry = flow
            .prepare_with_scalar_and_components(
                mc,
                hc,
                cc,
                &rho,
                &u,
                &heat,
                &mut thermal,
                values,
                &mut group,
                || false,
            )
            .unwrap();
        assert_eq!(retry, reports);
        assert_eq!(group.diffusivities(), expected);
        assert_eq!(group.diffusivities().as_ptr(), allocation);
        assert_eq!(group.workspace_bytes(), bytes);
    }
    assert!(matches!(
        flow.prepare_with_scalar_and_components(
            mc,
            hc,
            DynamicComponentsConfig {
                maximum_coefficient: 0.0,
                ..cc
            },
            &rho,
            &u,
            &heat,
            &mut thermal,
            values,
            &mut group,
            || false,
        ),
        Err(DynamicSgsError::Coefficient { .. })
    ));
    assert!(matches!(
        flow.prepare_with_scalar_and_components(
            mc,
            hc,
            cc,
            &rho,
            &u,
            &heat,
            &mut thermal,
            |s, cell| if s == 2 && cell + 1 == cells {
                f64::NAN
            } else {
                0.0
            },
            &mut group,
            || false,
        ),
        Err(DynamicSgsError::Cell { .. })
    ));
    let cleared = flow
        .prepare_with_scalar_and_components(
            mc,
            hc,
            cc,
            &rho,
            &u,
            &heat,
            &mut thermal,
            |s, _| if s == 0 { 1.0 } else { 0.0 },
            &mut group,
            || false,
        )
        .unwrap();
    assert_eq!(cleared.2.undetermined_cells, cells);
    assert!(group.diffusivities().iter().all(|v| *v == 0.0));
}
