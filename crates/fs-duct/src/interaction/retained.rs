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

struct Loaded {
    inverse: Vec<C64>,
    n: usize,
    reference: f64,
    first_mouth: usize,
    inflow_map: Vec<C64>,
    row_scales: Vec<f64>,
    condition: f64,
}

#[allow(clippy::too_many_arguments)]
fn loaded(
    nodes: usize,
    edges: &[WaveEdge],
    ports: &[CoupledPort],
    exterior: &GasState,
    omega: f64,
    interaction: OpenHoleInteraction,
    maximum_condition: f64,
    inlet_admittance: f64,
) -> Result<Loaded, DuctError> {
    let q = ports.len();
    if q > MAX_NETWORK_PORTS {
        return Err(bad(
            "open-hole retained mouth-flow cap exceeded (32, including the bell)",
        ));
    }
    let mut boundary = vec![C64::ZERO; nodes * nodes];
    boundary[0] = C64::from_re(inlet_admittance);
    let (mut system, n, reference) = wave_equations(&boundary, nodes, edges, q)?;
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
    // Match the graph owner's pressure-unit row equilibration. Retain each
    // scale: prescribed inlet flow drives KCL, while an orifice pressure drop
    // drives its original chimney equation, p_junction-(A Z+B)U_mouth=dp.
    let mut row_scales = Vec::with_capacity(n);
    for values in system.chunks_exact_mut(n) {
        let scale = values.iter().map(|z| z.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale == 0.0 || values.iter().any(|&z| !finite(z)) {
            return Err(DuctError::Singular);
        }
        for value in values {
            *value = C64::new(value.re / scale, value.im / scale);
        }
        row_scales.push(scale);
    }
    let (inverse, condition) = inverse_with_limit(&system, n, maximum_condition)?;
    Ok(Loaded { inverse, n, reference, first_mouth, inflow_map, row_scales, condition })
}

pub(super) fn solve(
    nodes: usize,
    edges: &[WaveEdge],
    ports: &[CoupledPort],
    exterior: &GasState,
    omega: f64,
    interaction: OpenHoleInteraction,
    maximum_condition: f64,
) -> Result<CoupledFields, DuctError> {
    let Loaded { inverse, n, reference, first_mouth, inflow_map, row_scales, condition } =
        loaded(nodes, edges, ports, exterior, omega, interaction, maximum_condition, 0.0)?;
    let q = ports.len();
    let drive = reference / row_scales[0];
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

/// The complete coupled network driven by a unit pressure drop at each
/// nonlinear orifice. The inlet remains matched and every other hole/bell is
/// linear. Only the first `orifices` chimney currents are returned, with the
/// same U=-Y*dp convention as the reduced nodal path.
#[allow(clippy::too_many_arguments)]
pub(super) fn port_admittance(
    nodes: usize,
    edges: &[WaveEdge],
    ports: &[CoupledPort],
    exterior: &GasState,
    omega: f64,
    interaction: OpenHoleInteraction,
    maximum_condition: f64,
    inlet_reference: f64,
    orifices: usize,
) -> Result<Vec<C64>, DuctError> {
    let Loaded { inverse, n, reference, first_mouth, inflow_map, row_scales, .. } =
        loaded(nodes, edges, ports, exterior, omega, interaction, maximum_condition,
            1.0 / inlet_reference)?;
    let q = ports.len();
    let mut port = vec![C64::ZERO; orifices * orifices];
    for j in 0..orifices {
        let row = first_mouth + j;
        let drive = 1.0 / row_scales[row];
        for i in 0..orifices {
            let inflow = (0..q).fold(C64::ZERO, |sum, k| {
                let w = inverse[(first_mouth + k) * n + row].scale(drive);
                sum + inflow_map[i * q + k] * C64::new(w.re / reference, w.im / reference)
            });
            if !finite(inflow) {
                return Err(DuctError::Singular);
            }
            port[i * orifices + j] = -inflow;
        }
    }
    Ok(port)
}
