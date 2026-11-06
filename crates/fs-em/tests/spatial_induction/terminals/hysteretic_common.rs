//! Independent six-tetrahedron electrode/play reduction, shared with finite heat.
use super::data;
use fs_em::spatial::vector_potential::induction::hysteretic::SpatialHystereticProblem;
use fs_material::field_storage::hysteresis::vector::VectorHysteresisState;

#[derive(Clone, Copy, Default)]
pub struct Endpoint {
    pub h: f64,
    pub z: f64,
    pub b: f64,
    pub voltage: f64,
    pub incoming: f64,
}
pub fn electrodes(problem: &SpatialHystereticProblem) -> Vec<Vec<usize>> {
    [0.0, 1.0]
        .into_iter()
        .map(|x| {
            problem
                .exterior_faces()
                .into_iter()
                .filter(|&face| {
                    problem.mesh().faces[face]
                        .iter()
                        .all(|&v| problem.positions_m()[v as usize][0] == x)
                })
                .collect()
        })
        .collect()
}
pub fn at_field(old: Endpoint, sigma: f64, dt: f64, h: f64) -> Endpoint {
    let z = data::play(old.z, h);
    let b = data::MU * h + z;
    let rate = (b - old.b) / (2.0 * dt);
    // The actual field is A=a W, B=2a with int W_x=1/6 and
    // int |W|²=1/5. All eight vertices lie on the x electrodes.
    // 2H+sigma*a'/5=sigma*U/6; I=sigma*(U-a'/6).
    let voltage = 12.0 * h / sigma + 6.0 * rate / 5.0;
    let incoming = sigma * (voltage - rate / 6.0);
    Endpoint {
        h,
        z,
        b,
        voltage,
        incoming,
    }
}
pub fn current_endpoint(old: Endpoint, sigma: f64, dt: f64, incoming: f64) -> Endpoint {
    let h = data::bracket(-29.0, 29.0, |h| {
        at_field(old, sigma, dt, h).incoming - incoming
    });
    at_field(old, sigma, dt, h)
}
pub fn joule(old: Endpoint, next: Endpoint, sigma: f64, dt: f64) -> f64 {
    // Condensing the electrode voltage subtracts sigma/36 from sigma/5.
    dt * next.incoming.powi(2) / sigma + 31.0 * sigma * (next.b - old.b).powi(2) / (720.0 * dt)
}
pub fn assert_memory(states: &[VectorHysteresisState], expected: Endpoint) {
    assert_eq!(states.len(), 6);
    for (memory, direction) in states.iter().zip(data::DIRECTIONS) {
        for axis in 0..3 {
            assert!((memory.field_a_m()[axis] - expected.h * direction[axis]).abs() < 4e-7);
            assert!((memory.flux_density_t()[axis] - expected.b * direction[axis]).abs() < 3e-8);
            assert!(
                (memory.polarizations_t()[0][axis] - expected.z * direction[axis]).abs() < 3e-8
            );
        }
    }
}
