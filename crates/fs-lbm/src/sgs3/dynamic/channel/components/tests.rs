use super::*;

const SHAPE: [usize; 3] = [5, 4, 6];
const SPACING: [f64; 3] = [0.3, 0.7, 1.1];

fn controls() -> (
    ChannelDynamicSgsConfig,
    ChannelDynamicScalarConfig,
    ChannelDynamicComponentsConfig,
) {
    (
        ChannelDynamicSgsConfig {
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        ChannelDynamicScalarConfig {
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        ChannelDynamicComponentsConfig {
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
    )
}

fn noise(i: usize, salt: u64) -> f64 {
    let mut x = (i as u64).wrapping_add(salt);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    ((x ^ (x >> 31)) >> 11) as f64 / ((1u64 << 53) as f64) - 0.5
}

fn fields(normal: usize) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>, [[f64; 3]; 2]) {
    let cells: usize = SHAPE.iter().product();
    let tangent = [(normal + 1) % 3, (normal + 2) % 3];
    let mut rho = Vec::new();
    let mut u = Vec::new();
    let mut phi = Vec::new();
    for cell in 0..cells {
        let p = coordinates(SHAPE, cell);
        let q = p[tangent[0]] + 11 * p[tangent[1]];
        let y = (p[normal] as f64 + 0.5) * SPACING[normal];
        let length = SHAPE[normal] as f64 * SPACING[normal];
        // Density varies within each plane but not along the normal. Thus
        // Favre filtering retains this exact quadratic normal profile.
        rho.push(1.2 + 0.2 * noise(q, 179));
        u.push(core::array::from_fn(|axis| {
            0.04 * noise(cell, 29 + axis as u64 * 79)
        }));
        phi.push(0.3 + 0.01 * noise(q, 431) * y * (length - y));
    }
    let mut walls = [[0.0; 3]; 2];
    walls[0][tangent[0]] = 0.013;
    walls[1][tangent[1]] = -0.017;
    (rho, u, phi, walls)
}

fn fit(normal: usize, components: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (rho, u, heat, walls) = fields(normal);
    let (mc, hc, cc) = controls();
    let mut flow = ChannelDynamicSgs::new_cartesian(SHAPE, SPACING, normal, usize::MAX).unwrap();
    let mut scalar = ChannelDynamicScalar::new(SHAPE, normal, usize::MAX).unwrap();
    let mut group =
        ChannelDynamicComponents::new(SHAPE, normal, components.len(), usize::MAX).unwrap();
    flow.prepare_with_scalar_and_components(
        mc,
        hc,
        cc,
        &rho,
        &u,
        walls,
        &heat,
        [ScalarWall::Value(0.3); 2],
        &mut scalar,
        |component, cell| components[component][cell],
        &mut group,
        || false,
    )
    .unwrap();
    (
        flow.viscosities().to_vec(),
        scalar.diffusivities().to_vec(),
        group.diffusivities().to_vec(),
    )
}

fn close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (cell, (&a, &b)) in a.iter().zip(b).enumerate() {
        assert!((a - b).abs() <= tolerance, "cell {cell}: {a} != {b}");
    }
}

#[test]
fn channel_components_reproduce_quadratic_wall_scalar_and_preserve_heat_flow() {
    for normal in 0..3 {
        let (rho, u, phi, walls) = fields(normal);
        let (mc, hc, _) = controls();
        let mut flow =
            ChannelDynamicSgs::new_cartesian(SHAPE, SPACING, normal, usize::MAX).unwrap();
        let mut scalar = ChannelDynamicScalar::new(SHAPE, normal, usize::MAX).unwrap();
        flow.prepare_with_scalar(
            mc,
            hc,
            &rho,
            &u,
            walls,
            &phi,
            [ScalarWall::Value(0.3); 2],
            &mut scalar,
            || false,
        )
        .unwrap();
        let (nu, heat, d) = fit(normal, &[phi.clone()]);
        assert_eq!(nu, flow.viscosities());
        assert_eq!(heat, scalar.diffusivities());
        // The declared-wall and interior stencils independently reproduce the
        // same quadratic's derivative, including a nonzero wall gradient.
        close(&d, &heat, 2e-12);
        assert!(d.iter().any(|v| *v > 0.0));
        let collinear = vec![
            phi.clone(),
            phi.iter().map(|v| 1.0 - 2.0 * v).collect(),
            vec![0.0; phi.len()],
            vec![7.0; phi.len()],
        ];
        close(&fit(normal, &collinear).2, &d, 2e-12);
    }
}

#[test]
fn channel_common_species_fit_preserves_basis_order_and_trace_weight() {
    let (_, _, first, _) = fields(2);
    let second: Vec<_> = (0..first.len())
        .map(|i| 0.25 + 0.03 * noise(i, 991))
        .collect();
    let third: Vec<_> = first
        .iter()
        .zip(&second)
        .map(|(a, b)| 1.0 - a - b)
        .collect();
    let components = vec![
        first.clone(),
        second.clone(),
        third.clone(),
        vec![0.0; first.len()],
    ];
    let baseline = fit(2, &components);
    assert!(baseline.2.iter().any(|d| *d > 0.0));
    let permuted = vec![
        components[3].clone(),
        third.clone(),
        first.clone(),
        second.clone(),
    ];
    close(&fit(2, &permuted).2, &baseline.2, 2e-12);
    let rotated = vec![
        first
            .iter()
            .zip(&second)
            .map(|(a, b)| 0.6 * a + 0.8 * b)
            .collect(),
        first
            .iter()
            .zip(&second)
            .map(|(a, b)| -0.8 * a + 0.6 * b)
            .collect(),
        third,
        vec![4.0; first.len()],
    ];
    let result = fit(2, &rotated);
    assert_eq!(result.0, baseline.0);
    assert_eq!(result.1, baseline.1);
    close(&result.2, &baseline.2, 2e-12);
    let dominant = fit(2, &[first.clone()]).2;
    let trace: Vec<_> = second.iter().map(|v| 1e-7 * v).collect();
    close(&fit(2, &[first, trace]).2, &dominant, 1e-11);
}

#[test]
fn interior_species_wall_gradient_allows_zero_flux_with_thermal_diffusion() {
    // j = -rho*D*grad(Y) - D_T*grad(log T). At rho*D=.02,
    // D_T=.004 and grad(log T)=.5, impermeability requires grad(Y)=-.1,
    // not zero. The scalar stencil must retain that composition gradient.
    for normal in 0..3 {
        let field: Vec<_> = (0..SHAPE.iter().product())
            .map(|cell| {
                let y = (coordinates(SHAPE, cell)[normal] as f64 + 0.5) * SPACING[normal];
                0.3 - 0.1 * y
            })
            .collect();
        let geometry = FilterGeometry::Channel(normal);
        for side in [0, SHAPE[normal] - 1] {
            let mut p = [1; 3];
            p[normal] = side;
            let cell = index(SHAPE, p);
            let stencil = geometry.stencil(SHAPE, SPACING.map(|h| 1.0 / h), cell, p, normal);
            let gradient = stencil.interior_scalar_derivative(&field).unwrap();
            assert!((gradient + 0.1).abs() < 2e-15);
            assert!((-0.02 * gradient - 0.004 * 0.5).abs() < 5e-17);
            let false_neumann = stencil
                .derivative(|i| field[i], [ScalarWall::ZeroGradient; 2])
                .unwrap();
            assert!((-0.02 * false_neumann - 0.004 * 0.5).abs() > 0.0009);
        }
    }
}

#[test]
fn channel_component_limits_cancellation_and_zero_retry_use_bounded_storage() {
    let (rho, u, phi, walls) = fields(2);
    let cells = phi.len();
    let bytes = ChannelDynamicComponents::required_bytes(cells).unwrap();
    assert_eq!(bytes, 80 * cells);
    assert!(ChannelDynamicComponents::new(SHAPE, 2, 3, bytes - 1).is_err());
    assert!(ChannelDynamicComponents::new([3, 3, 2], 2, 3, usize::MAX).is_err());
    assert!(ChannelDynamicComponents::new(SHAPE, 2, 0, usize::MAX).is_err());
    let mut flow = ChannelDynamicSgs::new_cartesian(SHAPE, SPACING, 2, usize::MAX).unwrap();
    let mut heat = ChannelDynamicScalar::new(SHAPE, 2, usize::MAX).unwrap();
    let mut group = ChannelDynamicComponents::new(SHAPE, 2, 2, bytes).unwrap();
    assert_eq!(group.workspace_bytes(), bytes);
    assert_eq!(group.component_count(), 2);
    let pointer = group.diffusivities().as_ptr();
    let (mc, hc, cc) = controls();
    let values = |s: usize, i: usize| if s == 0 { phi[i] } else { 1.0 - phi[i] };
    let run = |flow: &mut ChannelDynamicSgs,
               heat: &mut ChannelDynamicScalar,
               group: &mut ChannelDynamicComponents,
               config,
               cancel: &mut dyn FnMut() -> bool| {
        flow.prepare_with_scalar_and_components(
            mc,
            hc,
            config,
            &rho,
            &u,
            walls,
            &phi,
            [ScalarWall::Value(0.3); 2],
            heat,
            values,
            group,
            cancel,
        )
    };
    let report = run(&mut flow, &mut heat, &mut group, cc, &mut || false).unwrap();
    assert!(report.2.maximum_coefficient > 0.0);
    let accepted = group.diffusivities().to_vec();
    let mut limited = cc;
    limited.maximum_coefficient = 0.0;
    assert!(matches!(
        run(&mut flow, &mut heat, &mut group, limited, &mut || false),
        Err(DynamicSgsError::Coefficient { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        run(&mut flow, &mut heat, &mut group, cc, &mut || {
            polls += 1;
            polls == 15
        }),
        Err(DynamicSgsError::Cancelled)
    ));
    run(&mut flow, &mut heat, &mut group, cc, &mut || false).unwrap();
    assert_eq!(group.diffusivities(), accepted);
    assert_eq!(group.diffusivities().as_ptr(), pointer);
    flow.prepare_with_scalar_and_components(
        mc,
        hc,
        cc,
        &rho,
        &u,
        walls,
        &phi,
        [ScalarWall::Value(0.3); 2],
        &mut heat,
        |s, _| s as f64,
        &mut group,
        || false,
    )
    .unwrap();
    assert!(group.diffusivities().iter().all(|d| *d == 0.0));
    assert!(
        flow.prepare_with_scalar_and_components(
            mc,
            hc,
            cc,
            &rho,
            &u,
            walls,
            &phi,
            [ScalarWall::Value(0.3); 2],
            &mut heat,
            |_, _| f64::NAN,
            &mut group,
            || false
        )
        .is_err()
    );
}
