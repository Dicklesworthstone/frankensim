//! Molecular finite chimneys in the existing compact, independent-mouth model.
use super::*;
use crate::{HoleState, OpenHoleSite};

/// Physical branch geometry, in authored order. Radiation loads refer to the
/// fully open radius, before the existing vent-admittance weight is applied.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MolecularHoleSite {
    pub axial_position_m: f64,
    pub hole_radius_m: f64,
    pub bore_radius_m: f64,
    pub chimney_height_m: f64,
    pub open_chimney_length_m: f64,
    pub open_fraction: f64,
}

/// Linear branch response normalized to the same inlet volume flow as `duct`.
/// For a vent, open/closed flows carry their respective admittance weights;
/// mouth pressure does not. Thus mouth work is weighted once, not twice.
#[derive(Debug, Clone, Copy)]
pub struct MolecularHoleFlow {
    pub site: MolecularHoleSite,
    pub p_junction_over_u_in: C64,
    pub u_open_in_over_u_in: C64,
    pub u_closed_in_over_u_in: C64,
    pub u_mouth_over_u_in: C64,
    pub p_mouth_over_u_in: C64,
    /// Time-averaged branch input power divided by squared inlet flow [Pa s/m³].
    pub input_power_per_squared_inlet_flow: f64,
    /// Work delivered to the explicitly supplied mouth load, with the same units.
    pub mouth_power_per_squared_inlet_flow: f64,
    /// Molecular plus wall heat in the chimney; raw balance, never clipped.
    pub dissipated_power_per_squared_inlet_flow: f64,
}

#[derive(Debug, Clone)]
pub struct MolecularResponse {
    pub duct: ChainResult,
    /// Every authored tone hole, including sealed cavities with zero mouth flow.
    pub holes: Vec<MolecularHoleFlow>,
}

pub(super) fn opening(state: HoleState) -> Result<f64, DuctError> {
    match state {
        HoleState::Open => Ok(1.0),
        HoleState::Closed => Ok(0.0),
        HoleState::Vent(s) if s.is_finite() && (0.0..=1.0).contains(&s) => Ok(s),
        _ => Err(bad("molecular tone-hole opening must lie in [0,1]")),
    }
}

#[derive(Debug, Clone)]
pub(super) struct ToneHole {
    pub point: usize,
    pub site: MolecularHoleSite,
    series_length_m: f64,
}

// Projective input impedance: Z/R = n/d. One common branch amplitude X gives
// p_in=nX, R U_in=dX, p_mouth=mouth_pressure X, R U_mouth=mouth_flow X.
// No Z or Y is formed at a pole, and no growing exponential is evaluated.
struct Branch {
    n: C64,
    d: C64,
    mouth_flow: C64,
    mouth_pressure: C64,
}

impl Branch {
    fn prepare(
        radius: f64,
        length: f64,
        wave: &SegmentWave,
        reference: f64,
        load: Option<C64>,
    ) -> Result<Self, DuctError> {
        let phase = wave.wavenumber.scale(length);
        if !finite(phase) {
            return Err(DuctError::Singular);
        }
        let zc = wave
            .specific_impedance
            .scale(1.0 / (core::f64::consts::PI * radius * radius));
        let e = crate::exp_i(phase);
        let decay = det::exp(-2.0 * phase.im);
        let lost = -det::expm1(-2.0 * phase.im);
        let (s, c) = det::sin_cos(phase.re);
        // 1 +/- exp(2ikL), evaluated without cancellation at either a
        // low-frequency compliance or a quarter-wave pressure node.
        let plus = C64::new(lost + 2.0 * decay * c * c, 2.0 * decay * s * c);
        let minus = C64::new(lost + 2.0 * decay * s * s, -2.0 * decay * s * c);
        let ratio = C64::from_re(reference) / zc;
        let branch = match load {
            None => Self {
                n: plus,
                d: ratio * minus,
                mouth_pressure: e.scale(2.0),
                mouth_flow: C64::ZERO,
            },
            Some(z) => {
                let scale = z.abs().max(zc.abs()).max(reference);
                let load = z.scale(1.0 / scale);
                let characteristic = zc.scale(1.0 / scale);
                Self {
                    n: load * plus + characteristic * minus,
                    d: ratio * (load * minus + characteristic * plus),
                    mouth_pressure: (e * load).scale(2.0),
                    mouth_flow: e.scale(2.0 * (reference / scale)),
                }
            }
        };
        if [branch.n, branch.d, branch.mouth_pressure, branch.mouth_flow]
            .into_iter()
            .any(|z| !finite(z))
        {
            return Err(DuctError::Singular);
        }
        Ok(branch)
    }
}

pub(super) struct PreparedHole {
    block_index: usize,
    site: MolecularHoleSite,
    pressure: C64,
    open_flow: C64,
    closed_flow: C64,
    mouth_flow: C64,
    mouth_pressure: C64,
}

impl ToneHole {
    pub fn new(
        segment: &Segment,
        point: usize,
        x: f64,
        open_sites: &[OpenHoleSite],
    ) -> Result<Self, DuctError> {
        let Segment::ToneHole {
            hole_radius: b,
            chimney_height: h,
            bore_radius: a,
            state,
        } = *segment
        else {
            return Err(bad("molecular branch requires tone-hole geometry"));
        };
        let sigma = opening(state)?;
        let mut series_length_m = 0.0;
        if sigma > 0.0 {
            series_length_m = fs_phs::side_hole_series_length(b, a) * sigma;
            for site in open_sites {
                series_length_m += fs_phs::side_hole_mutual_length(
                    b * det::sqrt(sigma),
                    site.effective_radius_m,
                    a,
                    (site.axial_position_m - x).abs(),
                );
            }
        }
        Ok(Self {
            point,
            site: MolecularHoleSite {
                axial_position_m: x,
                hole_radius_m: b,
                bore_radius_m: a,
                chimney_height_m: h,
                open_chimney_length_m: h + fs_phs::side_hole_inner_length(b, a),
                open_fraction: sigma,
            },
            series_length_m,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        &self,
        wave: &SegmentWave,
        density: f64,
        omega: f64,
        mouth_load: Option<C64>,
        reference: f64,
        options: ComposeOptions,
        first_block: usize,
    ) -> Result<([TwoPort; 3], PreparedHole), DuctError> {
        let site = self.site;
        let sigma = site.open_fraction;
        let open = if sigma > 0.0 {
            Some(Branch::prepare(
                site.hole_radius_m,
                site.open_chimney_length_m,
                wave,
                reference,
                mouth_load,
            )?)
        } else {
            None
        };
        let closed = if sigma < 1.0 {
            Some(Branch::prepare(
                site.hole_radius_m,
                site.chimney_height_m,
                wave,
                reference,
                None,
            )?)
        } else {
            None
        };
        // Y = sigma Y_open + (1-sigma) Y_closed, in homogeneous coordinates.
        let (n, open_d, closed_d, mouth_flow, mouth_pressure) = match (open, closed) {
            (Some(o), Some(c)) => (
                o.n * c.n,
                (o.d * c.n).scale(sigma),
                (c.d * o.n).scale(1.0 - sigma),
                (o.mouth_flow * c.n).scale(sigma),
                o.mouth_pressure * c.n,
            ),
            (Some(o), None) => (o.n, o.d, C64::ZERO, o.mouth_flow, o.mouth_pressure),
            (None, Some(c)) => (c.n, C64::ZERO, c.d, C64::ZERO, C64::ZERO),
            (None, None) => unreachable!("admitted opening"),
        };
        let d = open_d + closed_d;
        let denominator = n.scale(2.0) + d;
        let condition = (2.0 * n.abs() + d.abs()) / denominator.abs();
        if !condition.is_finite() || condition > options.maximum_condition {
            return Err(DuctError::Singular);
        }
        let twice = C64::from_re(2.0) / denominator;
        let transmission = n * twice;
        let reflection = -d / denominator;
        let shunt = TwoPort::new(
            [reflection, transmission, transmission, reflection],
            reference,
        )?;
        let area = core::f64::consts::PI * site.bore_radius_m * site.bore_radius_m;
        let half = C64::new(0.0, -0.5 * omega * density * self.series_length_m / area);
        let series = TwoPort::from_abcd([C64::ONE, half, C64::ZERO, C64::ONE], reference)?;
        Ok((
            [series, shunt, series],
            PreparedHole {
                block_index: first_block + 1,
                site,
                pressure: transmission,
                open_flow: open_d * twice,
                closed_flow: closed_d * twice,
                mouth_flow: mouth_flow * twice,
                mouth_pressure: mouth_pressure * twice,
            },
        ))
    }
}

impl MolecularDuct {
    /// Open/vented mouths in the exact order required by `mouth_loads`.
    pub fn open_hole_sites(&self) -> Vec<MolecularHoleSite> {
        self.pieces
            .iter()
            .filter_map(|piece| match piece {
                Piece::ToneHole(hole) if hole.site.open_fraction > 0.0 => Some(hole.site),
                _ => None,
            })
            .collect()
    }

    /// Original stable blocks, with one explicit finite passive load per open
    /// or vented mouth. Loads are p/U in Pa s/m³ at the fully open radius b;
    /// the vent mixture applies its sigma weight exactly once. Zero denotes
    /// pressure release. Closed holes require no entry. A radiation, measured,
    /// tabulated or network owner supplies loads; no exterior gas is inferred.
    pub fn prepared_chain_with_hole_loads_controlled(
        &self,
        omega: f64,
        mouth_loads: &[C64],
        options: ComposeOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedChain, MolecularDuctError> {
        self.prepare_with_holes_controlled(omega, mouth_loads, options, cancelled)
            .map(|(chain, _)| chain)
    }

    /// Solved terminal and branch flows/work from one molecular scattering
    /// chain. Per-mouth scalar loads retain the independent-opening model;
    /// a coupled exterior radiation matrix requires an exposed-port network.
    pub fn response_with_holes(
        &self,
        omega: f64,
        terminal_load: Option<C64>,
        mouth_loads: &[C64],
        options: ComposeOptions,
    ) -> Result<MolecularResponse, MolecularDuctError> {
        self.response_with_holes_controlled(omega, terminal_load, mouth_loads, options, &mut || {
            false
        })
    }

    pub fn response_with_holes_controlled(
        &self,
        omega: f64,
        terminal_load: Option<C64>,
        mouth_loads: &[C64],
        options: ComposeOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularResponse, MolecularDuctError> {
        let (chain, prepared) =
            self.prepare_with_holes_controlled(omega, mouth_loads, options, cancelled)?;
        let duct =
            compose::compose_chain_controlled(&chain.blocks, terminal_load, options, cancelled)?;
        let incident = if prepared.is_empty() {
            Vec::new()
        } else {
            compose::incident_waves_controlled(&chain, terminal_load, options, cancelled)?
        };
        let mut holes = Vec::with_capacity(prepared.len());
        for (i, hole) in prepared.into_iter().enumerate() {
            poll(cancelled, i)?;
            let incoming = incident[hole.block_index];
            let sum = incoming[0] + incoming[1];
            let pressure = (hole.pressure * sum).scale(chain.reference_impedance);
            let open = hole.open_flow * sum;
            let closed = hole.closed_flow * sum;
            let mouth = hole.mouth_flow * sum;
            let mouth_pressure = (hole.mouth_pressure * sum).scale(chain.reference_impedance);
            let input_power = 0.5 * (pressure * (open + closed).conj()).re;
            let mouth_power = 0.5 * (mouth_pressure * mouth.conj()).re;
            if [pressure, open, closed, mouth, mouth_pressure]
                .into_iter()
                .any(|z| !finite(z))
                || !input_power.is_finite()
                || !mouth_power.is_finite()
            {
                return Err(DuctError::Singular.into());
            }
            holes.push(MolecularHoleFlow {
                site: hole.site,
                p_junction_over_u_in: pressure,
                u_open_in_over_u_in: open,
                u_closed_in_over_u_in: closed,
                u_mouth_over_u_in: mouth,
                p_mouth_over_u_in: mouth_pressure,
                input_power_per_squared_inlet_flow: input_power,
                mouth_power_per_squared_inlet_flow: mouth_power,
                dissipated_power_per_squared_inlet_flow: input_power - mouth_power,
            });
        }
        Ok(MolecularResponse { duct, holes })
    }
}
