//! Keep physical mouth flows until both radiation and the complete bore are loaded.
//!
//! Inverting A Z+B is only a coordinate choice. Its null vector can be fixed by
//! the bore equations even at an exact chimney resonance. The native graph
//! assembly retains the bore's two inward currents; this module adds each
//! mouth current and its original chimney/radiation constraint. No damping,
//! impedance perturbation or isolated admittance inverse is introduced.

use super::{CoupledFields, CoupledPort, OpenHoleInteraction, bad, chimney_matrices, finite};
use crate::DuctError;
use crate::network::{MAX_NETWORK_PORTS, WaveEdge, inverse_with_limit, wave_equations};
use fs_material::gas::GasState;
use fs_math::c64::C64;

pub(super) fn solve(
    nodes: usize,
    edges: &[WaveEdge],
    ports: &[CoupledPort],
    exterior: &GasState,
    omega: f64,
    interaction: OpenHoleInteraction,
    maximum_condition: f64,
) -> Result<CoupledFields, DuctError> {
    let q = ports.len();
    if q > MAX_NETWORK_PORTS {
        return Err(bad(
            "open-hole retained mouth-flow cap exceeded (32, including the bell)",
        ));
    }
    let (mut system, n, reference) =
        wave_equations(&vec![C64::ZERO; nodes * nodes], nodes, edges, q)?;
    let (pressure_map, inflow_map) = chimney_matrices(ports, exterior, omega, interaction)?;
    let first_mouth = n - q;
    for i in 0..q {
        let node = i + 1;
        let row = first_mouth + i;
        system[row * n + node] = C64::ONE;
        for j in 0..q {
            let column = first_mouth + j;
            // Unknown w = reference * U_mouth, like the native bore currents.
            // R*KCL receives R*U_chimney = (C Z+D)w; the additional row is
            // p_junction - (A Z+B)w/R = 0. Mouth pressure remains Z*U_mouth.
            system[node * n + column] = inflow_map[i * q + j];
            let z = pressure_map[i * q + j];
            system[row * n + column] = -C64::new(z.re / reference, z.im / reference);
        }
    }
    // Match the graph owner's pressure-unit row equilibration. The sole RHS
    // is R times the prescribed unit inlet flow, on the first KCL row.
    let mut drive = 0.0;
    for (row, values) in system.chunks_exact_mut(n).enumerate() {
        let scale = values.iter().map(|z| z.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale == 0.0 || values.iter().any(|&z| !finite(z)) {
            return Err(DuctError::Singular);
        }
        for value in values {
            *value = C64::new(value.re / scale, value.im / scale);
        }
        if row == 0 {
            drive = reference / scale;
        }
    }
    let (inverse, condition) = inverse_with_limit(&system, n, maximum_condition)?;
    let solution: Vec<C64> = (0..n).map(|i| inverse[i * n].scale(drive)).collect();
    if solution.iter().any(|&z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    let pressure = solution[..nodes].to_vec();
    let mouth: Vec<C64> = solution[first_mouth..]
        .iter()
        .map(|z| C64::new(z.re / reference, z.im / reference))
        .collect();
    let inflow: Vec<C64> = (0..q)
        .map(|i| (0..q).fold(C64::ZERO, |sum, j| sum + inflow_map[i * q + j] * mouth[j]))
        .collect();
    if mouth.iter().chain(&inflow).any(|&z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    Ok(CoupledFields {
        pressure,
        mouth,
        inflow,
        condition,
    })
}
