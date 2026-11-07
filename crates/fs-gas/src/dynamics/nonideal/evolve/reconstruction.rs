//! MUSCL traces admitted by the original real-gas constitutive source.

use super::*;

mod linear;

pub(super) type Trace = (Conserved, Primitive);
pub(super) type Traces = [Trace; 2];

pub(super) fn prepare(
    gas: &RealGas,
    states: &[Conserved],
    primitive: &[Primitive],
    mesh: &Mesh1d,
    boundaries: RealGasBoundaries,
    output: &mut [Traces],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, RealGasFlowError> {
    let n = states.len();
    let periodic = boundaries == RealGasBoundaries::Periodic;
    let width = |i: usize| {
        mesh.uniform_cell_width()
            .map_or_else(|| mesh.width(i), Length::value)
    };
    let values = |p: Primitive| [p.density, p.velocity, p.temperature];
    let mut limited = 0;
    for i in 0..n {
        poll(cancelled)?;
        let original = (states[i], primitive[i]);
        // End cells have no exterior center stencil. Retain the complete native
        // boundary relation; prescribed ghosts and physical walls are unchanged.
        if !periodic && (i == 0 || i + 1 == n) {
            output[i] = [original; 2];
            continue;
        }
        let left = (i + n - 1) % n;
        let right = (i + 1) % n;
        let increment = linear::increments(
            values(primitive[left]),
            values(primitive[i]),
            values(primitive[right]),
            [width(left), width(i), width(right)],
        );
        let pair = linear::admitted_pair(
            original,
            values(primitive[i]),
            increment,
            |[rho, velocity, temperature]| gas.state_at(rho, velocity, temperature, cancelled),
            |error| matches!(error, RealGasFlowError::Cancelled),
        )?;
        output[i] = pair.traces;
        limited += usize::from(pair.limited);
    }
    Ok(limited)
}
