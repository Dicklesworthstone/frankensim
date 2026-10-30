//! G1/G3: the rotary continuum against its independent Bessel boundary solution.
use super::*;

const RADIUS: f64 = 0.0003;
const DENSITY: f64 = 1.2;
const VISCOSITY: f64 = 1.8e-5;
fn resolution() -> CylinderFlowResolution {
    CylinderFlowResolution { minimum_wavenumber_radius: 1e-5,
        maximum_wavenumber_radius: 3000.0, intervals: 4096 }
}
fn flow(density: f64, viscosity: f64, radius: f64, resolution: CylinderFlowResolution)
    -> RotaryCylinderFlow<FluidTransportState>
{
    RotaryCylinderFlow::from_fluid(radius,
        FluidTransportState::try_new(density, 1000.0, 0.025, viscosity).unwrap(),
        resolution, &CancelGate::new_clock_free()).unwrap()
}
fn normalized(beta: f64, resolution: CylinderFlowResolution) -> C64 {
    let flow = flow(DENSITY, VISCOSITY, RADIUS, resolution);
    let omega = beta*beta*VISCOSITY/DENSITY/RADIUS/RADIUS;
    let z = flow.impedance_per_length(omega).unwrap();
    let scale = 2.0*core::f64::consts::PI*VISCOSITY*RADIUS*RADIUS;
    C64 { re: z.re/scale, im: z.im/scale }
}
fn distance(a: C64, b: C64) -> f64 { (a.re-b.re).hypot(a.im-b.im) }
fn near(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-14*a.abs().max(b.abs()), "{a:e} != {b:e}");
}

#[test]
fn continuum_converges_to_independent_modified_bessel_boundary_solution() {
    // z=exp(i*pi/4)*beta, target 2+z*K0(z)/K1(z). These constants were computed
    // from K_n(z)=integral_0^infinity exp(-z*cosh(t))*cosh(n*t) dt (DLMF10.32.9),
    // independent of the production J/Y spectrum. Composite Simpson at 32768
    // and 65536 panels, cut where Re(z)*cosh(t)=60, agreed in all shown digits.
    for (beta, re, im) in [
        (0.01, 2.000_078_405_334_562_6, 0.000_472_151_141_282_545_5),
        (0.1, 2.007_441_252_667_307, 0.024_408_656_227_612_784),
        (1.0, 2.400_235_576_852_202_3, 0.615_800_155_654_143_1),
        (10.0, 8.597_308_259_022_743, 7.047_962_398_373_279),
    ] {
        let expected = C64 { re, im };
        let actual = normalized(beta, resolution());
        let error = distance(actual, expected)/re.hypot(im);
        assert!(error < 0.002, "beta={beta}: {actual:?}, oracle={expected:?}, error={error:e}");
        assert!(actual.re > 2.0 && actual.im > 0.0);
    }
}

#[test]
fn each_spectral_resolution_axis_has_an_independent_ladder() {
    let target = C64 { re: 2.400_235_576_852_202_3, im: 0.615_800_155_654_143_1 };
    let high: Vec<_> = [30.0, 300.0, 3000.0].iter().map(|&cutoff|
        distance(normalized(1.0, CylinderFlowResolution {
            maximum_wavenumber_radius: cutoff, ..resolution() }), target)).collect();
    assert!(high[1] < 0.11*high[0] && high[2] < 0.11*high[1], "high cutoff: {high:?}");

    let low_target = C64 { re: 2.000_078_405_334_562_6, im: 0.000_472_151_141_282_545_5 };
    let low: Vec<_> = [0.5, 0.05, 0.005].iter().map(|&cutoff|
        distance(normalized(0.01, CylinderFlowResolution {
            minimum_wavenumber_radius: cutoff, ..resolution() }), low_target)).collect();
    assert!(low[1] < 0.05*low[0] && low[2] < 0.05*low[1], "low cutoff: {low:?}");

    // Fix both cutoffs: this isolates log-grid error from the omitted continuum.
    let fixed = CylinderFlowResolution { minimum_wavenumber_radius: 0.01,
        maximum_wavenumber_radius: 30.0, ..resolution() };
    let reference = normalized(1.0, fixed);
    let grid: Vec<_> = [32, 64, 128].iter().map(|&intervals|
        distance(normalized(1.0, CylinderFlowResolution { intervals, ..fixed }), reference)).collect();
    assert!(grid[1] < 0.3*grid[0] && grid[2] < 0.3*grid[1], "log grid: {grid:?}");
}

#[test]
fn material_transport_and_radius_scale_true_rotary_memory() {
    let res = CylinderFlowResolution { intervals: 7, ..resolution() };
    let a = flow(DENSITY, VISCOSITY, RADIUS, res);
    let rho = flow(4.0*DENSITY, VISCOSITY, RADIUS, res);
    let mu = flow(DENSITY, 4.0*VISCOSITY, RADIUS, res);
    let radius = flow(DENSITY, VISCOSITY, 2.0*RADIUS, res);
    near(a.steady_resistance_per_length_n_s(), 4.0*core::f64::consts::PI*VISCOSITY*RADIUS*RADIUS);
    near(rho.steady_resistance_per_length_n_s(), a.steady_resistance_per_length_n_s());
    near(mu.steady_resistance_per_length_n_s(), 4.0*a.steady_resistance_per_length_n_s());
    near(radius.steady_resistance_per_length_n_s(), 4.0*a.steady_resistance_per_length_n_s());
    assert_eq!(a.transport(), a.fluid());
    for (((a, rho), mu), radius) in a.branches().iter().zip(rho.branches()).zip(mu.branches()).zip(radius.branches()) {
        near(a.resistance_per_length_n_s/a.rate_per_s, a.inertia_per_length_kg_m);
        near(rho.resistance_per_length_n_s, a.resistance_per_length_n_s);
        near(rho.inertia_per_length_kg_m, 4.0*a.inertia_per_length_kg_m);
        near(mu.inertia_per_length_kg_m, a.inertia_per_length_kg_m);
        near(mu.rate_per_s, 4.0*a.rate_per_s);
        near(radius.inertia_per_length_kg_m, 16.0*a.inertia_per_length_kg_m);
        assert!(a.rate_per_s > 0.0 && a.resistance_per_length_n_s > 0.0);
    }
}

#[test]
fn invalid_geometry_resolution_frequency_and_cancellation_refuse() {
    let fluid = FluidTransportState::try_new(DENSITY, 1000.0, 0.025, VISCOSITY).unwrap();
    let gate = CancelGate::new_clock_free();
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE, f64::MAX] {
        assert!(RotaryCylinderFlow::from_fluid(radius, fluid, resolution(), &gate).is_err());
    }
    for intervals in [0, 4097] {
        assert!(RotaryCylinderFlow::from_fluid(RADIUS, fluid,
            CylinderFlowResolution { intervals, ..resolution() }, &gate).is_err());
    }
    let flow = flow(DENSITY, VISCOSITY, RADIUS, resolution());
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(flow.impedance_per_length(omega).is_err());
    }
    gate.request();
    assert!(matches!(RotaryCylinderFlow::from_fluid(RADIUS, fluid, resolution(), &gate),
        Err(AcousticRealizeError::Cancelled)));
}
