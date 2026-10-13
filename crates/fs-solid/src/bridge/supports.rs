//! Frictionless point supports on the existing thermoviscoelastic beam.
//!
//! For the SAME discrete curvature basis w(x)=sum psi_k(x)*kappa_k,
//! a point reaction R_i adds moment psi_k(x_i)*R_i/h_k to cell k.
//! The implicit section owner supplies kappa_k=c_k*M_k+kappa_free,k.
//! Thus the contact flexibility is H_ij=sum psi_ki*c_k/h_k*psi_kj.
//! Its explicit Gram factor admits the existing L1 nonnegative quadratic solver:
//! min_{R>=0} 0.5 R^T H R + gap_free^T R. Its gradient is the physical gap.
//! No penalty stiffness, invented adhesion, finite-difference compliance or new
//! constitutive/optimization solver is introduced. This is discrete Signorini
//! contact at declared REFERENCE stations, not collision discovery or a proof
//! that a continuous obstacle is not penetrated between those stations.

use super::{Bridge, Config as BeamConfig, Support, finite, poll, positive};
use fs_solver::orthant::{GramOrthantConfig, GramOrthantError, solve_gram_orthant};
use std::collections::BTreeMap;

/// A horizontal, frictionless stop acting upward at a material reference station.
/// Axial motion is free; a fixed world-x obstacle needs a different kinematic map.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stop {
    pub id: u64,
    pub reference_x_m: f64,
    /// Minimum allowed REFERENCE-LINE vertical displacement w, in metres.
    /// For a contacted fiber at z, convert a surface elevation to elevation-z.
    /// Different endpoint bounds represent explicitly moved/removed supports.
    pub minimum_deflection_m: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// Gradient units are metres; columns are axial cells, variables are stops.
    pub contact: GramOrthantConfig,
    /// Fresh gap/KKT gate evaluated on the final, actually updated beam, m.
    pub gap_tolerance_m: f64,
    /// Fresh max |reaction * gap|, J. Not a time-integrated work allowance.
    pub complementarity_tolerance_j: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reaction {
    pub stop: Stop,
    pub upward_force_n: f64,
    pub gap_m: f64,
}

#[derive(Debug, Clone)]
pub struct Report {
    /// Existing physical beam report; includes support forces in its work.
    pub beam: super::Report,
    /// Stable stop-ID order, independent of input arrival order.
    pub reactions: Vec<Reaction>,
    pub contact_sweeps: usize,
    pub maximum_gap_error_m: f64,
    pub maximum_complementarity_j: f64,
    /// Endpoint reaction quadrature sum R_new*(w_new-w_old), J. For newly
    /// engaged or moved stops this need not vanish and is not exact path work.
    pub reaction_work_j: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Beam(super::Error),
    Contact(GramOrthantError),
    Invalid(&'static str),
    Residual { gap_m: f64, complementarity_j: f64 },
}
impl From<super::Error> for Error {
    fn from(error: super::Error) -> Self { Self::Beam(error) }
}
impl From<GramOrthantError> for Error {
    fn from(error: GramOrthantError) -> Self { Self::Contact(error) }
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "beam support refused: {self:?}")
    }
}
impl std::error::Error for Error {}

impl<'a> Bridge<'a> {
    /// Executor-scoped form with the same physical update and cancellation points.
    pub fn advance_supported_to(&self, cx: &fs_exec::Cx<'_>, target_s: f64,
        temperatures_k: &BTreeMap<(u64, u64), f64>, gravity_m_s2: f64,
        additional_downward_n_m: &[f64], stops: &[Stop], beam_config: BeamConfig,
        config: Config) -> Result<(Self, Report), Error> {
        self.advance_supported_controlled(target_s, temperatures_k, gravity_m_s2,
            additional_downward_n_m, stops, beam_config, config, &mut || cx.checkpoint().is_err())
    }

    /// Solve support reactions, then advance every actual material point ONCE.
    /// All constitutive probes use the original history. The free-beam prediction
    /// is not a material update and cannot wrongly refuse a supported solution
    /// solely because the unsupported shape would exceed its regime.
    ///
    /// Only the returned candidate commits. Missing temperatures, invalid stops,
    /// exhausted solve budgets, final strain/shape/energy/gap failures and
    /// cancellation return no changed beam. Config bounds at most 64 stops and
    /// 4096 cells; factor storage is stops*cells, setup work stops^2*cells, plus
    /// existing linear-time section/beam passes. No whole-field history log.
    ///
    /// Contact does not automatically change thermal area/conductance or bond
    /// strength. Built-in pin/clamp restraints remain bilateral. There is no
    /// large-deformation, impact, stability, continuous-obstacle or unique-force
    /// claim. Refine the spatial/time/contact grids independently.
    pub fn advance_supported_controlled(&self, target_s: f64,
        temperatures_k: &BTreeMap<(u64, u64), f64>, gravity_m_s2: f64,
        additional_downward_n_m: &[f64], stops: &[Stop], beam_config: BeamConfig,
        config: Config, cancelled: &mut impl FnMut() -> bool) -> Result<(Self, Report), Error> {
        poll(cancelled)?;
        positive(finite(target_s-self.time_s)?)?;
        positive(config.gap_tolerance_m)?;
        positive(config.complementarity_tolerance_j)?;
        let n = self.cells.len();
        let p = stops.len();
        let products = p.checked_mul(p).and_then(|v| v.checked_mul(n))
            .ok_or(Error::Invalid("contact setup overflow"))?;
        if p > config.contact.max_variables || p > 64 || n > config.contact.max_columns
            || n > 4096 || products > config.contact.max_setup_products {
            return Err(Error::Invalid("contact dimensions or setup exceed declared budget"));
        }
        // Delegate full policy validation to the numerical owner before setup.
        solve_gram_orthant(&[], &[], config.contact, &mut *cancelled)?;
        let layers = self.cells[0].section.layers().len();
        if layers == 0 || temperatures_k.len() != n.checked_mul(layers).ok_or(Error::Invalid("temperature count overflow"))? {
            return Err(Error::Invalid("born span and complete endpoint temperatures required"));
        }
        let mut ordered = BTreeMap::new();
        for &stop in stops {
            poll(cancelled)?;
            let x = stop.reference_x_m;
            if !x.is_finite() || x <= 0.0 || x > self.length_m
                || (self.support == Support::SimplySupported && x == self.length_m)
                || !stop.minimum_deflection_m.is_finite() || ordered.insert(stop.id, stop).is_some() {
                return Err(Error::Invalid("duplicate ID, nonfinite stop or station at a fixed/outside endpoint"));
            }
        }
        let ordered: Vec<_> = ordered.into_values().collect();
        let Prediction { free_curvature, roots } = predict(self, target_s, temperatures_k,
            gravity_m_s2, additional_downward_n_m, cancelled)?;
        let mut factor = Vec::with_capacity(p);
        let mut free_gaps = Vec::with_capacity(p);
        for stop in &ordered {
            poll(cancelled)?;
            let mut row = Vec::with_capacity(n);
            let mut gap = -stop.minimum_deflection_m;
            for (k, cell) in self.cells.iter().enumerate() {
                if k % 256 == 0 { poll(cancelled)?; }
                let psi = curvature_basis(self.support, self.length_m, cell.start_m,
                    cell.length_m(), stop.reference_x_m);
                row.push(finite(psi*roots[k])?);
                gap = finite(gap+psi*free_curvature[k])?;
            }
            factor.push(row);
            free_gaps.push(gap);
        }
        let solution = solve_gram_orthant(&factor, &free_gaps, config.contact, &mut *cancelled)?;
        let forces: Vec<_> = ordered.iter().zip(&solution.point)
            .map(|(stop, &force)| (stop.reference_x_m, force)).collect();
        let (next, beam) = self.advance_with_point_forces(target_s, temperatures_k,
            gravity_m_s2, additional_downward_n_m, &forces, beam_config, cancelled)?;
        let mut reactions = Vec::with_capacity(p);
        let mut maximum_gap_error_m = 0.0_f64;
        let mut maximum_complementarity_j = 0.0_f64;
        let mut reaction_work_j = 0.0;
        for (stop, force) in ordered.into_iter().zip(solution.point) {
            poll(cancelled)?;
            let (old, s_old) = self.locate(stop.reference_x_m)?;
            let (new, s_new) = next.locate(stop.reference_x_m)?;
            let gap = finite(new.height(s_new)-stop.minimum_deflection_m)?;
            maximum_gap_error_m = maximum_gap_error_m.max(if force > 0.0 { gap.abs() } else { (-gap).max(0.0) });
            maximum_complementarity_j = maximum_complementarity_j.max(finite(force*gap)?.abs());
            reaction_work_j = finite(reaction_work_j+force*(new.height(s_new)-old.height(s_old)))?;
            reactions.push(Reaction { stop, upward_force_n: force, gap_m: gap });
        }
        if maximum_gap_error_m > config.gap_tolerance_m
            || maximum_complementarity_j > config.complementarity_tolerance_j {
            return Err(Error::Residual { gap_m: maximum_gap_error_m, complementarity_j: maximum_complementarity_j });
        }
        poll(cancelled)?;
        Ok((next, Report { beam, reactions, contact_sweeps: solution.sweeps,
            maximum_gap_error_m, maximum_complementarity_j, reaction_work_j }))
    }
}

// Exact integral of a cell's constant curvature, satisfying the built-in
// support conditions. Units m2. This is the adjoint of the point-force moment
// transfer, not a fitted Green function or finite-difference beam probe.
pub(super) fn curvature_basis(support: Support, length: f64, start: f64, h: f64, x: f64) -> f64 {
    let s = (x-start).max(0.0);
    let primitive = if s < h { 0.5*s*s } else { h*(s-0.5*h) };
    match support {
        Support::CantileverLeft => primitive,
        Support::SimplySupported => primitive-x*h*((length-start-0.5*h)/length),
    }
}

// Shared by rigid-stop and deformable-neighbour contact. This is the SAME
// implicit section response, not a second beam or constitutive approximation.
pub(super) struct Prediction {
    pub free_curvature: Vec<f64>,
    pub roots: Vec<f64>,
}
pub(super) fn predict(beam: &Bridge<'_>, target_s: f64,
    temperatures_k: &BTreeMap<(u64, u64), f64>, gravity_m_s2: f64,
    additional_downward_n_m: &[f64], cancelled: &mut impl FnMut() -> bool,
) -> Result<Prediction, super::Error> {
    positive(finite(target_s-beam.time_s)?)?;
    let n = beam.cells.len();
    let layers = beam.cells[0].section.layers().len();
    if layers == 0 || Some(temperatures_k.len()) != n.checked_mul(layers) {
        return Err(super::Error::Invalid("born span and complete endpoint temperatures required"));
    }
    let load = beam.load_field(gravity_m_s2, additional_downward_n_m, &[], cancelled)?;
    let mut free_curvature = Vec::with_capacity(n);
    let mut roots = Vec::with_capacity(n);
    for (cell, &moment) in beam.cells.iter().zip(&load.moments) {
        poll(cancelled)?;
        let mut t = BTreeMap::new();
        for &layer in cell.section.layers().keys() {
            let value = *temperatures_k.get(&(cell.id, layer))
                .ok_or(super::Error::Invalid("missing layer temperature"))?;
            t.insert(layer, value);
        }
        let response = cell.section.bending_response(target_s, &t,
            || if cancelled() { Err(crate::deposition::Error::Cancelled) } else { Ok(()) })
            .map_err(|error| super::Error::Section { cell: cell.id, error })?;
        free_curvature.push(finite(response.zero_moment_curvature_per_m + response.compliance_per_n_m2*moment)?);
        let root = finite(fs_math::det::sqrt(response.compliance_per_n_m2/cell.length_m()))?;
        positive(root)?;
        roots.push(root);
    }
    Ok(Prediction { free_curvature, roots })
}

#[cfg(test)]
mod tests;
