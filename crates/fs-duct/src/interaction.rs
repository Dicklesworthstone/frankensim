//! External radiation interaction between the open tone holes and the open end.
//!
//! The transfer-matrix method gives every open hole an independent radiation load. Lefebvre,
//! Scavone and Kergomard (Acta Acustica 99, 2013, arXiv 1207.5490, "TMMI") couple the openings
//! through one full radiation impedance matrix `Z`: the pressure at the mouth of opening `n` is
//! `Z_nn U_n + sum_m Z_nm U_m`, with the monopole mutual impedance of their Eq. 16 (Pritchard
//! 1960),
//!
//! `Z_nm = j k rho c exp(-j k d_nm) / (eps pi d_nm)`, `eps = 2` (half space) for every pair,
//!
//! written there for `exp(+j w t)`; this crate uses `exp(-i w t)`, so the conjugate is used.
//! Each hole keeps its own chimney (inner length, loss, wall) and its own self radiation as the
//! diagonal `Z_nn`; the bore between consecutive openings uses the original scattering
//! two-ports. A small nodal solve is retained away from component poles; singular isolated
//! admittances retain their pressure/flow equations until the complete load is applied.
//! With no coupling (`eps` -> infinity) the result is the
//! independent-hole response to rounding.
//!
//! Scope, stated plainly: the distances are the AXIAL distances between hole stations (the
//! azimuths of the holes are not known), which is an upper bound for holes on opposite sides of
//! the bore, but never less than the sum of the two mouth radii: two openings cannot overlap, so
//! their centres are at least that far apart, and two holes at ONE station (a doubled hole) are
//! coupled at that touching distance, the strongest coupling their geometry allows and the same
//! kind of upper bound (no spacing is invented; the single-station case used to refuse);
//! every open hole must be fully open (a vent fraction has no single mouth flow); the
//! bell is one more opening when the termination radiates. The feedback response and reflection
//! and the hole mouth flows use the coupled solution. The nonlinear orifice port reflection of
//! each hole keeps its independent chimney.

use fs_la::eigen_complex::lu_complex;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_math::det;

use crate::compose::{ComposeOptions, PreparedChain, PreparedOpenHole, TwoPort, duct_chain_with_wall};
use crate::{Duct, DuctError, LossModel, Segment, Termination, WallModel};

#[cfg(test)]
mod tests;

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Opt-in coupling of the open holes and the bell through their mutual radiation impedance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpenHoleInteraction {
    epsilon: f64,
}

impl OpenHoleInteraction {
    /// `epsilon` is the radiation-space factor of the mutual impedance: 2 (half space) as the
    /// TMMI paper prints it for every pair of openings, 4 for a complete space.
    ///
    /// # Errors
    /// A nonpositive or nonfinite factor refuses.
    pub fn new(epsilon: f64) -> Result<Self, DuctError> {
        if !(epsilon.is_finite() && epsilon > 0.0) {
            return Err(bad("open-hole interaction needs a positive finite radiation-space factor"));
        }
        Ok(Self { epsilon })
    }

    /// The radiation-space factor of the mutual impedance.
    #[must_use]
    pub const fn epsilon(self) -> f64 {
        self.epsilon
    }
}

/// Mutual radiation impedance of two monopole openings a distance `distance_m` apart, radiating
/// into the same exterior, under `exp(-i w t)`: `-i k rho c exp(+i k d) / (eps pi d)`.
pub(crate) fn mutual_radiation_impedance(
    exterior: &GasState,
    omega: f64,
    distance_m: f64,
    interaction: OpenHoleInteraction,
) -> Result<C64, DuctError> {
    if !(distance_m > 0.0 && distance_m.is_finite()) {
        return Err(bad("open-hole interaction needs distinct finite opening stations"));
    }
    let k = omega / exterior.sound_speed;
    let (s, c) = det::sin_cos(k * distance_m);
    let scale = k * exterior.density * exterior.sound_speed
        / (interaction.epsilon * core::f64::consts::PI * distance_m);
    let z = C64::new(scale * s, -scale * c);
    if finite(z) { Ok(z) } else { Err(DuctError::Singular) }
}

/// One radiating opening: its station, its chimney as `[p_bore, U_in] = M [p_mouth, U_mouth]`
/// (`U_mouth` leaves the mouth) and the self-radiation load that closes the chimney.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CoupledPort {
    pub position_m: f64,
    /// Radius of the mouth: two openings are never closer than the sum of their radii. Zero for the bore's open end.
    pub radius_m: f64,
    pub chimney: [C64; 4],
    pub self_radiation: C64,
}

impl CoupledPort {
    /// The open end of the bore: no chimney, the termination's load is its self radiation.
    pub(crate) fn open_end(position_m: f64, load: C64) -> Self {
        Self { position_m, radius_m: 0.0, chimney: [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE], self_radiation: load }
    }

    /// A fully open tone hole built exactly as `chimney_open_response_with_wall` builds it: the
    /// chimney of the given total `length` (height plus the hole law's lengths, as
    /// `PreparedOpenHole::chimney_length_m` carries it), the same loss, wall and flanged self radiation.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn open_hole(
        position_m: f64,
        state: &GasState,
        exterior: &GasState,
        hole_radius: f64,
        length: f64,
        omega: f64,
        loss: LossModel,
        wall: WallModel<'_>,
    ) -> Result<Self, DuctError> {
        let chimney_loss = if loss == LossModel::WideTube { LossModel::Bessel } else { loss };
        let duct = Duct { segments: vec![Segment::Cylinder { radius: hole_radius, length }] };
        let options = ComposeOptions::default();
        let chain = duct_chain_with_wall(&duct, state, omega, chimney_loss, wall, options, false, || false)?;
        let mut line = TwoPort::through(chain.reference_impedance)?;
        for block in &chain.blocks {
            line = line.star(*block, options.maximum_condition)?.0;
        }
        let (radiation, _) = crate::exterior::termination_impedance(
            Termination::FlangedOpen, state, exterior, hole_radius, omega)?;
        let self_radiation = radiation.ok_or_else(|| bad("an open hole needs a radiating mouth"))?;
        Ok(Self { position_m, radius_m: hole_radius, chimney: line.abcd()?, self_radiation })
    }
}

/// The distance two openings are coupled at: the axial distance between their stations, but never
/// less than the sum of their mouth radii (two mouths cannot overlap). Where the axial distance is
/// the larger, as for every pair of separate holes on a real bore, it is the original distance
/// bit for bit; two holes at one station are coupled at the touching distance.
fn pair_distance(a: &CoupledPort, b: &CoupledPort) -> f64 {
    (a.position_m - b.position_m).abs().max(a.radius_m + b.radius_m)
}

/// Everything the coupling gives, per opening in the order of the supplied ports.
struct Coupling {
    /// Flow into each chimney from its bore junction = `admittance` (q x q, row-major) times the junction pressures.
    admittance: Vec<C64>,
    /// Flow leaving each mouth = `mouth_admittance` times the junction pressures.
    mouth_admittance: Vec<C64>,
    /// Infinity-norm condition of the jointly loaded chimney matrix.
    condition: f64,
}

#[cfg(test)]
fn couple(
    ports: &[CoupledPort], exterior: &GasState, omega: f64, interaction: OpenHoleInteraction,
) -> Result<Coupling, DuctError> {
    couple_with_limit(ports, exterior, omega, interaction, ComposeOptions::default().maximum_condition)
}

/// Keep the complete mouth load before eliminating any chimney: with diagonal A, B, C, D,
/// `p = (A Z + B) U_mouth` and `U_bore = (C Z + D) U_mouth`. An isolated `A_i Z_ii + B_i`
/// can vanish even when `A Z + B` is invertible. Dividing by it first creates false poles
/// and makes the reconstructed mouth flow a cancellation of divergent terms near resonance.
fn couple_with_limit(
    ports: &[CoupledPort], exterior: &GasState, omega: f64, interaction: OpenHoleInteraction,
    maximum_condition: f64,
) -> Result<Coupling, DuctError> {
    let q = ports.len();
    let (mut m, mut inflow) = (vec![C64::ZERO; q * q], vec![C64::ZERO; q * q]);
    for i in 0..q {
        let [a, b, c, d] = ports[i].chimney;
        m[i * q + i] = a * ports[i].self_radiation + b;
        inflow[i * q + i] = c * ports[i].self_radiation + d;
        for j in i + 1..q {
            // The exterior is reciprocal; evaluate each mutual kernel once.
            let z = mutual_radiation_impedance(exterior, omega, pair_distance(&ports[i], &ports[j]), interaction)?;
            m[i * q + j] = a * z;
            m[j * q + i] = ports[j].chimney[0] * z;
            inflow[i * q + j] = c * z;
            inflow[j * q + i] = ports[j].chimney[2] * z;
        }
    }
    let (mouth_admittance, condition) = crate::network::inverse_with_limit(&m, q, maximum_condition)?;
    let mut admittance = vec![C64::ZERO; q * q];
    for i in 0..q {
        for j in 0..q {
            admittance[i * q + j] = (0..q).fold(C64::ZERO, |sum, k| {
                sum + inflow[i * q + k] * mouth_admittance[k * q + j]
            });
        }
    }
    if admittance.iter().any(|&z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    Ok(Coupling { admittance, mouth_admittance, condition })
}

/// The coupled solution for a unit input flow.
pub(crate) struct CoupledSolution {
    pub input_impedance: C64,
    /// Pressure at the bore end and flow leaving it (zero for a closed end).
    pub bell_pressure: C64,
    pub bell_flow: C64,
    /// Per open hole: flow leaving the mouth, and flow entering the chimney from the bore.
    pub hole_mouth_flows: Vec<C64>,
    pub hole_chimney_inflows: Vec<C64>,
    pub condition: f64,
}

/// Solve the bore between openings and the openings through [`couple_with_limit`].
/// A half-wave run has singular I+S, not necessarily a singular loaded instrument.
/// Such runs use the existing graph owner's retained scattering constraints.
/// `holes` and `ports` are the open holes in axial order; `bell` is the radiating open end, or
/// `None` for a closed end.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve(
    chain: &PreparedChain,
    holes: &[PreparedOpenHole],
    ports: &[CoupledPort],
    bell: Option<CoupledPort>,
    exterior: &GasState,
    omega: f64,
    interaction: OpenHoleInteraction,
    maximum_condition: f64,
) -> Result<CoupledSolution, DuctError> {
    let m = holes.len();
    if m == 0 || ports.len() != m {
        return Err(bad("open-hole interaction needs one coupled port per open hole"));
    }
    if holes.iter().any(|hole| hole.sigma != 1.0) {
        return Err(bad("open-hole interaction needs fully open holes; a vent fraction has no single mouth flow"));
    }
    let reference = chain.reference_impedance;
    let half = |hole: &PreparedOpenHole| {
        TwoPort::from_abcd([C64::ONE, hole.half_series_impedance, C64::ZERO, C64::ONE], reference)
    };
    let nodes = m + 2;
    let mut ytot = vec![C64::ZERO; nodes * nodes];
    let mut condition = 1.0_f64;
    let mut propagation_condition = 1.0_f64;
    let mut wave_edges = Vec::with_capacity(m + 1);
    let mut keep_waves = false;
    for k in 0..=m {
        let begin = if k == 0 { 0 } else { holes[k - 1].block_index + 1 };
        let end = if k == m { chain.blocks.len() } else { holes[k].block_index };
        if begin > end || end > chain.blocks.len() {
            return Err(bad("open-hole blocks must be ordered along the chain"));
        }
        let mut run = if k == 0 { TwoPort::through(reference)? } else { half(&holes[k - 1])? };
        let right = if k < m { Some(half(&holes[k])?) } else { None };
        for block in chain.blocks[begin..end].iter().copied().chain(right) {
            let (next, c) = run.star(block, maximum_condition)?;
            run = next;
            condition = condition.max(c);
            propagation_condition = propagation_condition.max(c);
        }
        wave_edges.push(crate::network::WaveEdge {
            ends: [Some(k), Some(k + 1)],
            two_port: run,
            condition: propagation_condition,
        });
        match run.admittance(maximum_condition) {
            Ok((y, c)) => {
                condition = condition.max(c);
                for (row, column, value) in [(k, k, y[0]), (k, k + 1, y[1]), (k + 1, k, y[2]), (k + 1, k + 1, y[3])] {
                    ytot[row * nodes + column] = ytot[row * nodes + column] + value;
                }
            }
            Err(DuctError::Singular) => keep_waves = true,
            Err(error) => return Err(error),
        }
    }
    let mut all: Vec<CoupledPort> = ports.to_vec();
    all.extend(bell);
    let q = all.len();
    let coupling = couple_with_limit(&all, exterior, omega, interaction, maximum_condition)?;
    condition = condition.max(coupling.condition);
    for i in 0..q {
        for j in 0..q {
            let at = (i + 1) * nodes + j + 1;
            ytot[at] = ytot[at] + coupling.admittance[i * q + j];
        }
    }
    let nodal = if keep_waves {
        Err(DuctError::Singular)
    } else {
        lu_complex(&ytot, nodes).map_err(|_| DuctError::Singular).and_then(|lu| {
            let mut p = vec![C64::ZERO; nodes];
            p[0] = C64::ONE;
            lu.solve(&mut p);
            if p.iter().all(|v| finite(*v)) { Ok(p) } else { Err(DuctError::Singular) }
        })
    };
    let p = match nodal {
        Ok(p) => p,
        Err(DuctError::Singular) => {
            // Only the aperture boundary is stamped here. The graph owner
            // retains every bore current, so no isolated I+S is inverted.
            let mut boundary_y = vec![C64::ZERO; nodes * nodes];
            for i in 0..q { for j in 0..q {
                boundary_y[(i + 1) * nodes + j + 1] = coupling.admittance[i * q + j];
            } }
            let (inverse, c) = crate::network::wave_pressure_inverse(
                &boundary_y, nodes, &wave_edges, maximum_condition)?;
            condition = propagation_condition.max(c).max(coupling.condition);
            (0..nodes).map(|i| inverse[i * nodes]).collect()
        }
        Err(error) => return Err(error),
    };
    let junction = &p[1..=q];
    let row = |matrix: &[C64], i: usize| (0..q).fold(C64::ZERO, |sum, j| sum + matrix[i * q + j] * junction[j]);
    let inflow: Vec<C64> = (0..q).map(|i| row(&coupling.admittance, i)).collect();
    let mouth: Vec<C64> = (0..q).map(|i| row(&coupling.mouth_admittance, i)).collect();
    // A passive one-port has a nonnegative input resistance.
    if p[0].re < -1.0e-9 * p[0].abs() {
        return Err(bad("coupled radiation is active at this frequency"));
    }
    Ok(CoupledSolution {
        input_impedance: p[0],
        bell_pressure: p[m + 1],
        bell_flow: if bell.is_some() { inflow[m] } else { C64::ZERO },
        hole_mouth_flows: mouth[..m].to_vec(),
        hole_chimney_inflows: inflow[..m].to_vec(),
        condition,
    })
}
