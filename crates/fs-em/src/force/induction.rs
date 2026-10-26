//! Accepted induction histories use the same geometric Maxwell-stress owner.
//!
//! The history, not a detached recovered field, knows the actual impressed
//! current and every resolved eddy conductor. Those sources must be absent on
//! both sides of the contour; enclosed nonlinear materials and currents are
//! allowed. B/H are recovered from the accepted A without a static re-solve.

use super::{
    MaxwellForceLimits,
    contour::{self, ContourError},
};
use crate::{
    planar::{
        MagneticError,
        eddy::{
            EddyConductor, EddyError, PlanarEddy,
            nonlinear::{NonlinearEddyError, PlanarNonlinearEddy},
        },
        nonlinear::IsotropicMagneticMaterial,
        poll,
    },
    transfer::PlanarFrame,
};
use fs_ga::Wrench;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_rep_mesh::{Metric2, TriComplex2};

fn admit_recovery(
    mesh: &TriComplex2,
    limits: MaxwellForceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MagneticError> {
    poll(cancelled)?;
    if mesh.faces().len() > limits.max_cells || mesh.edges().len() > limits.max_edges {
        return Err(MagneticError::Invalid(
            "induction-force recovery exceeds mesh work bounds",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn wrench<'a>(
    mesh: &TriComplex2,
    conductors: &[EddyConductor],
    current: &[f64],
    mut material: impl FnMut(usize) -> Option<&'a ResolvedMaterialStatePoint>,
    b: &[[f64; 2]],
    h: &[[f64; 2]],
    enclosed_cells: &[usize],
    frame: &PlanarFrame,
    moment_origin_m: [f64; 3],
    limits: MaxwellForceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench, MagneticError> {
    let Metric2::Planar { thickness } = mesh.metric() else {
        return Err(MagneticError::Invalid(
            "induction force requires the original planar metric",
        ));
    };
    let mut conducting = vec![false; mesh.faces().len()];
    for source in conductors {
        poll(cancelled)?;
        conducting[source.triangle] = true;
    }
    let mut exterior = None;
    contour::integrate(
        mesh, thickness, enclosed_cells, frame, moment_origin_m, limits,
        |inner, outer, _| {
            let inside = material(inner).ok_or(ContourError::Invalid(
                "induction force contour requires an induced-linear exterior material",
            ))?;
            let outside = material(outer).ok_or(ContourError::Invalid(
                "induction force contour requires an induced-linear exterior material",
            ))?;
            let source = *exterior.get_or_insert(inside);
            if inside != source || outside != source
                || current[inner] != 0.0 || current[outer] != 0.0
                || conducting[inner] || conducting[outer]
            {
                return Err(ContourError::Invalid(
                    "induction contour needs one homogeneous nonconducting current-free material on both sides",
                ));
            }
            Ok((b[outer], h[outer]))
        },
        cancelled,
    ).map_err(|error| match error {
        ContourError::Invalid(reason) => MagneticError::Invalid(reason),
        ContourError::Cancelled => MagneticError::Cancelled,
    })
}

impl PlanarEddy<'_> {
    /// Integrate instantaneous Maxwell force from the accepted induction field.
    ///
    /// The original history supplies its actual current snapshot and conductor
    /// support. Both sides of every contour edge must be nonconducting, carry
    /// zero impressed current and share one complete linear exterior material.
    /// Resolved eddy currents may lie inside the contour. B/H come from the
    /// retained potential without resetting or re-equilibrating that history.
    ///
    /// The same closed-contour integration, extrusion depth, mechanical frame,
    /// moment convention and refinement limits as `PlanarField::maxwell_wrench`
    /// apply. The returned native wrench is an instantaneous observation: it
    /// neither debits field energy nor advances a moving body. Moving geometry
    /// needs its own work closure; do not add an overlapping Lorentz load.
    /// Failure or cancellation leaves the accepted field and time unchanged.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, EddyError> {
        admit_recovery(self.problem().mesh(), limits, cancelled)?;
        let field = self.field(cancelled)?;
        Ok(wrench(
            self.problem().mesh(),
            self.conductors(),
            self.impressed_current_density_a_m2(),
            |cell| Some(&self.problem().materials()[cell]),
            field.magnetic_flux_density_t(),
            field.magnetic_field_a_m(),
            enclosed_cells,
            frame,
            moment_origin_m,
            limits,
            cancelled,
        )?)
    }
}

impl PlanarNonlinearEddy<'_> {
    /// Integrate the actual nonlinear induction/circuit field into a native wrench.
    ///
    /// Enclosed bodies may use their original reversible H(B) laws, stranded
    /// impressed currents and resolved eddy currents. The integration contour
    /// must remain in one complete induced-linear, nonconducting exterior source
    /// with zero impressed current on both incident cells. This does not extend
    /// the linear exterior stress formula into nonlinear or conducting material.
    ///
    /// `NonlinearWindingCircuit::field()` exposes this SAME accepted history,
    /// so circuit voltage, source saturation and resistivity influence the
    /// computed load without fitted force coefficients or a static field reset.
    /// The linear induction owner shares all source admission and geometry.
    ///
    /// This is an instantaneous magnetoquasistatic net load in the fixed planar
    /// section, with no displacement-current/electric-stress, end-field,
    /// magnetostriction or exterior pressure-response claim. It does not debit
    /// magnetic energy, solve motion, or close moving-boundary mechanical work.
    /// Use the original frame/origin conventions and avoid duplicate body loads.
    pub fn maxwell_wrench(
        &self,
        enclosed_cells: &[usize],
        frame: &PlanarFrame,
        moment_origin_m: [f64; 3],
        limits: MaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench, NonlinearEddyError> {
        admit_recovery(self.problem().mesh(), limits, cancelled)?;
        let field = self.field(cancelled)?;
        Ok(wrench(
            self.problem().mesh(),
            self.conductors(),
            self.impressed_current_density_a_m2(),
            |cell| match &self.problem().materials()[cell] {
                IsotropicMagneticMaterial::Linear(source) => Some(source),
                IsotropicMagneticMaterial::Reversible(_) => None,
            },
            field.magnetic_flux_density_t(),
            field.magnetic_field_a_m(),
            enclosed_cells,
            frame,
            moment_origin_m,
            limits,
            cancelled,
        )?)
    }
}
