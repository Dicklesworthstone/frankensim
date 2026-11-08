use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::time::HeterogeneousTimeLimits;
use crate::highorder::simplex::LocalFn;
use fs_math::det;
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

const MEDIUM: Medium = Medium { density_kg_m3: 2.0, sound_speed_m_s: 2.0 };
const PLACEHOLDER: Medium = Medium { density_kg_m3: 1000.0, sound_speed_m_s: 1500.0 };
const POSITIONS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
];

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-12, "{actual:e} != {expected:e}");
}

fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(geometry, geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM), &mut || false).unwrap()
}

#[test]
fn coherent_surface_flux_retains_sub_ulp_pressure_differences() {
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&complex, &POSITIONS, 1, PLACEHOLDER).unwrap();
    let model = uniform(&geometry);
    let surface = model.surface_probe(&[ElementFace { element: 0, vertices: [0, 2, 3] }],
        2, SurfaceProbeLimits::default(), &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false).unwrap();
    let large = 2.0_f64.powi(53);
    let coefficients = [C64::new(large, large), C64::new(large - 2.0, large - 2.0),
        C64::new(large + 4.0, large), C64::new(large, large)];
    let power = surface.harmonic(&coefficients, 1.0 / core::f64::consts::TAU,
        &mut || false).unwrap();
    // On x=0: p=(large+4y)+i*large and outward velocity=-1+i.
    // Active intensity is -2y, so its unit-triangle integral is -1/3 W.
    close(power.active_power_w, -1.0 / 3.0);
    assert_eq!(surface.harmonic(&coefficients, 1.0 / core::f64::consts::TAU,
        &mut || false).unwrap(), power);
}

fn p2(geometry: &HelmholtzFe<'_>, field: impl Fn([f64; 3]) -> C64) -> Vec<C64> {
    let mut coefficients = vec![C64::ZERO; geometry.ndof()];
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        for (dof, shape) in geometry.space.element_dofs(element) {
            coefficients[dof] = match shape {
                LocalFn::Vertex(a) => field(geometry.positions()[tet[a] as usize]),
                LocalFn::Edge(a, b, 0) => {
                    let a = geometry.positions()[tet[a] as usize];
                    let b = geometry.positions()[tet[b] as usize];
                    let middle = core::array::from_fn(|axis| f64::midpoint(a[axis], b[axis]));
                    field(middle).scale(4.0) - (field(a) + field(b)).scale(2.0)
                }
                _ => panic!("P2 fixture requires original vertex and edge coordinates"),
            };
        }
    }
    coefficients
}

#[test]
fn interface_power_keeps_original_density_outward_sign_and_winding_invariance() {
    let positions = [
        [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0],
    ];
    let complex = TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]);
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let media = [MEDIUM, Medium { density_kg_m3: 8.0, sound_speed_m_s: 3.0 }];
    let queries = Cell::new(0);
    let mut source = |element: usize, _: [f64; 3]| {
        queries.set(queries.get() + 1);
        Ok::<_, AcousticFeError>(media[element])
    };
    let model = HeterogeneousHelmholtzFe::new(&geometry, 3, &mut source, &mut || false).unwrap();
    queries.set(0);
    let omega = 2.0;
    let velocity = 0.4;
    let pressure = [C64::ONE, C64::ONE, C64::ONE,
        C64::new(1.0, omega * media[0].density_kg_m3 * velocity),
        C64::new(1.0, -omega * media[1].density_kg_m3 * velocity)];
    for (element, vertices, sign) in [
        (0, [0, 1, 2], 1.0), (0, [2, 1, 0], 1.0), (1, [0, 1, 2], -1.0),
    ] {
        let face = ElementFace { element, vertices };
        let surface = model.surface_probe(&[face], 2, SurfaceProbeLimits::default(),
            &mut source, &mut || false).unwrap();
        close(surface.area_m2(), 0.5);
        assert_eq!(surface.sample_count(), 4);
        assert_eq!(surface.faces(), &[face]);
        let count = queries.get();
        let power = surface.harmonic(&pressure, omega / core::f64::consts::TAU,
            &mut || false).unwrap();
        close(power.active_power_w, sign * 0.25 * velocity);
        close(power.reactive_power_w, 0.0);
        assert_eq!(surface.harmonic(&pressure, omega / core::f64::consts::TAU,
            &mut || false).unwrap(), power);
        assert_eq!(queries.get(), count);
    }
    let faces = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
        .map(|vertices| ElementFace { element: 0, vertices });
    let closed = model.surface_probe(&faces, 2, SurfaceProbeLimits::default(),
        &mut source, &mut || false).unwrap();
    let power = closed.harmonic(&pressure, omega / core::f64::consts::TAU,
        &mut || false).unwrap();
    close(power.active_power_w, 0.0);
    close(power.reactive_power_w, -omega * media[0].density_kg_m3 * velocity * velocity / 12.0);
}

#[test]
fn high_order_surface_products_and_evanescent_reactive_power_are_integrated() {
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&complex, &POSITIONS, 2, PLACEHOLDER).unwrap();
    let model = uniform(&geometry);
    let face = ElementFace { element: 0, vertices: [0, 2, 3] }; // outward normal -x
    let surface = model.surface_probe(&[face], 3, SurfaceProbeLimits::default(),
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM), &mut || false).unwrap();
    assert_eq!(surface.quadrature_order(), 3);
    let omega = 2.0;
    let velocity = 0.7;
    let a = 2.0;
    let k = omega * MEDIUM.density_kg_m3 * velocity;
    let pressure = p2(&geometry, |p| C64::new(1.0, -k * p[0]).scale(1.0 + a * p[1]));
    let power = surface.harmonic(&pressure, omega / core::f64::consts::TAU,
        &mut || false).unwrap();
    // On the unit right triangle: integral(1+a*y)^2 = 1/2+a/3+a²/12.
    // Normal derivatives of functions with zero face trace are essential here.
    let squared_amplitude_area = 0.5 + a / 3.0 + a * a / 12.0;
    close(power.active_power_w, -0.5 * velocity * squared_amplitude_area);
    close(power.reactive_power_w, 0.0);

    let tangential_k = 1.4;
    let alpha = det::sqrt(tangential_k * tangential_k - 1.0);
    // Exact normal-decay jet of an evanescent field. At x=0, d_x p=-alpha p;
    // its normal intensity is purely reactive despite nonzero pressure/velocity.
    let evanescent = p2(&geometry, |p| {
        C64::new(1.0, -tangential_k * p[1]).scale(1.0 - alpha * p[0])
    });
    let power = surface.harmonic(&evanescent, omega / core::f64::consts::TAU,
        &mut || false).unwrap();
    close(power.active_power_w, 0.0);
    close(power.reactive_power_w,
        -alpha * (0.5 + tangential_k * tangential_k / 12.0) / (2.0 * omega * MEDIUM.density_kg_m3));
}

#[test]
fn transient_power_uses_pressure_potential_gauge_and_no_harmonic_half_factor() {
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&complex, &POSITIONS, 2, PLACEHOLDER).unwrap();
    let model = uniform(&geometry);
    let surface = model.surface_probe(&[ElementFace { element: 0, vertices: [0, 2, 3] }],
        3, SurfaceProbeLimits::default(), &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false).unwrap();
    let time = model.prepare_time(0.01, &[], HeterogeneousTimeLimits::default(),
        &mut || false).unwrap();
    let a = 2.0;
    let velocity = 0.7;
    let pressure = p2(&geometry, |p| C64::from_re(1.0 + a * p[1]))
        .into_iter().map(|v| v.re).collect::<Vec<_>>();
    for gauge in [0.0, 17.0] {
        let potential = p2(&geometry, |p| C64::from_re(
            gauge - MEDIUM.density_kg_m3 * velocity * p[0] * (1.0 + a * p[1])))
            .into_iter().map(|v| v.re).collect::<Vec<_>>();
        let state = time.state(&potential, &pressure).unwrap();
        close(surface.transient(&state, &mut || false).unwrap(),
            -velocity * (0.5 + a / 3.0 + a * a / 12.0));
        assert!(matches!(surface.transient(&state, &mut || true), Err(HeterogeneousError::Cancelled)));
    }
    let other = uniform(&geometry);
    let other_time = other.prepare_time(0.01, &[], HeterogeneousTimeLimits::default(),
        &mut || false).unwrap();
    assert!(surface.transient(&other_time.zero_state(), &mut || false).is_err());
}

#[test]
fn surface_admission_checks_the_whole_selection_before_queries_and_preserves_source_errors() {
    #[derive(Debug, PartialEq)]
    enum SourceError { Original, Interior(AcousticFeError) }
    impl From<AcousticFeError> for SourceError {
        fn from(error: AcousticFeError) -> Self { Self::Interior(error) }
    }
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&complex, &POSITIONS, 2, PLACEHOLDER).unwrap();
    let model = uniform(&geometry);
    let face = ElementFace { element: 0, vertices: [0, 2, 3] };
    let reversed = ElementFace { vertices: [3, 2, 0], ..face };
    let queries = Cell::new(0);
    let mut source = |_, _| { queries.set(queries.get() + 1); Ok::<_, AcousticFeError>(MEDIUM) };
    for faces in [vec![], vec![face, reversed],
        vec![face, ElementFace { element: 1, vertices: [0, 1, 2] }],
        vec![ElementFace { element: 0, vertices: [0, 0, 1] }]] {
        assert!(model.surface_probe(&faces, 3, SurfaceProbeLimits::default(),
            &mut source, &mut || false).is_err());
        assert_eq!(queries.get(), 0);
    }
    for (order, limits) in [
        (2, SurfaceProbeLimits::default()),
        (MAX_QUADRATURE_ORDER + 1, SurfaceProbeLimits::default()),
        (3, SurfaceProbeLimits { max_faces: 0, ..SurfaceProbeLimits::default() }),
        (3, SurfaceProbeLimits { max_quadrature_points: 8, ..SurfaceProbeLimits::default() }),
    ] {
        assert!(model.surface_probe(&[face], order, limits, &mut source, &mut || false).is_err());
        assert_eq!(queries.get(), 0);
    }
    assert!(matches!(model.surface_probe(&[face], 3, SurfaceProbeLimits::default(),
        &mut |_, _| Err::<Medium, _>(SourceError::Original), &mut || false),
        Err(SourceError::Original)));
    assert!(model.surface_probe(&[face], 3, SurfaceProbeLimits::default(),
        &mut source, &mut || queries.get() > 0).is_err());
    assert_eq!(queries.get(), 1);
    let surface = model.surface_probe(&[face], 3, SurfaceProbeLimits::default(),
        &mut source, &mut || false).unwrap();
    let pressure = vec![C64::ZERO; model.ndof()];
    assert!(matches!(surface.harmonic(&pressure, 1.0, &mut || true), Err(HeterogeneousError::Cancelled)));
    assert!(surface.harmonic(&pressure[..pressure.len() - 1], 1.0, &mut || false).is_err());
    assert!(surface.harmonic(&pressure, 0.0, &mut || false).is_err());
}
