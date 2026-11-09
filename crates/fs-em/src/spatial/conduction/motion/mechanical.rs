//! Mechanical handoff of the current already accepted by the DC solve.

use super::{
    SpatialConduction, SpatialConductorMaterial, SpatialMovingConductionField,
    SpatialMovingConductionSnapshot,
};
use crate::spatial::vector_potential::induction::mechanical::{
    SpatialLorentzLoadError, SpatialLorentzLoadLimits, SpatialLorentzLoads,
    load::{CellLoad, integrate},
};

type Result<T> = core::result::Result<T, SpatialLorentzLoadError>;

impl<C: SpatialConductorMaterial> SpatialMovingConductionField<'_, C> {
    /// Transfer selected original-cell `J cross B` forces to native mechanics.
    ///
    /// The retained current and magnetic field own the force: extraction does
    /// not solve another electrical problem or reevaluate its force law. The
    /// original P1 mass integration reduces to `volume * force_density / 4`
    /// at each tetrahedron vertex for this cellwise-constant field. Shared
    /// vertices accumulate once, preserving force, torque and every P1 virtual
    /// velocity pairing. Wrenches reuse the canonical `fs_ga` frame/origin map.
    ///
    /// Physical work equals the selected `cell_mechanical_power_w` sum when
    /// the receiving P1 velocity has the prescribed mean in every selected
    /// cell. Use [`Self::lorentz_loads_for_velocity`] to check that condition;
    /// discontinuous cell velocities need not admit one continuous nodal trace.
    /// These are conduction-current volume forces only. Do not add overlapping
    /// induction/Maxwell loads, or interpret extraction as mechanical feedback.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzLoads> {
        self.recover_lorentz_loads(selected_cells, None, limits, cancelled)
    }

    /// Transfer the same forces only if the mechanical velocity reproduces
    /// the motion that generated their electrical current.
    ///
    /// `nodal_velocity_m_s` is in the original mesh coordinates/vertex order.
    /// Its four-node average must match each selected cell's prescribed
    /// velocity within `velocity_tolerance_m_s` in every component. This local
    /// check prevents opposite work defects cancelling in a global receipt.
    /// The nodal-power discrepancy is at most `velocity_tolerance_m_s` times
    /// the selected volume integral of the force density's componentwise
    /// absolute sum, plus integration roundoff.
    /// The tolerance is absolute, finite and nonnegative; zero admits exact
    /// means. Virtual velocities remain available through the returned loads.
    pub fn lorentz_loads_for_velocity(
        &self,
        selected_cells: &[usize],
        nodal_velocity_m_s: &[[f64; 3]],
        velocity_tolerance_m_s: f64,
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzLoads> {
        self.recover_lorentz_loads(
            selected_cells,
            Some((nodal_velocity_m_s, velocity_tolerance_m_s)),
            limits,
            cancelled,
        )
    }

    fn recover_lorentz_loads(
        &self,
        selected_cells: &[usize],
        velocity: Option<(&[[f64; 3]], f64)>,
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzLoads> {
        recover_lorentz_loads(
            self.problem(),
            self.velocity_m_s(),
            self.lorentz_force_density_n_m3(),
            selected_cells,
            velocity,
            limits,
            cancelled,
        )
    }
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    // Called only by the original nonlinear candidate, which retains the
    // operator/geometry paired with this immutable native moving snapshot.
    pub(in crate::spatial::conduction) fn moving_snapshot_lorentz_loads(
        &self,
        field: &SpatialMovingConductionSnapshot,
        selected_cells: &[usize],
        velocity: Option<(&[[f64; 3]], f64)>,
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzLoads> {
        recover_lorentz_loads(
            self,
            field.velocity_m_s(),
            field.lorentz_force_density_n_m3(),
            selected_cells,
            velocity,
            limits,
            cancelled,
        )
    }
}

// The borrowed fixed-source field and owned nonlinear candidate share the
// original cell integration, velocity-pairing check and canonical wrench map.
#[allow(clippy::too_many_arguments)]
fn recover_lorentz_loads<C: SpatialConductorMaterial>(
    problem: &SpatialConduction<C>,
    cell_velocity_m_s: &[[f64; 3]],
    force_density_n_m3: &[[f64; 3]],
    selected_cells: &[usize],
    velocity: Option<(&[[f64; 3]], f64)>,
    limits: SpatialLorentzLoadLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialLorentzLoads> {
    use SpatialLorentzLoadError::{Cancelled, Invalid};
    if cancelled() {
        return Err(Cancelled);
    }
    if selected_cells.is_empty()
        || selected_cells.len() > limits.max_cells
        || problem.positions_m().len() > limits.max_vertices
        || cell_velocity_m_s.len() != problem.mesh().tets.len()
        || force_density_n_m3.len() != problem.mesh().tets.len()
    {
        return Err(Invalid(
            "empty, incomplete or over-budget Lorentz load selection",
        ));
    }
    if let Some((nodal, tolerance)) = velocity {
        if nodal.len() != problem.positions_m().len() || !tolerance.is_finite() || tolerance < 0.0 {
            return Err(Invalid("invalid mechanical velocity trace or tolerance"));
        }
        for value in nodal {
            if cancelled() {
                return Err(Cancelled);
            }
            if value.iter().any(|component| !component.is_finite()) {
                return Err(Invalid("nonfinite mechanical velocity trace"));
            }
        }
    }
    let mut selected = Vec::with_capacity(selected_cells.len());
    for &cell in selected_cells {
        if cancelled() {
            return Err(Cancelled);
        }
        let vertices = *problem.mesh().tets.get(cell).ok_or(Invalid(
            "Lorentz selection requires original conducting cells",
        ))?;
        if let Some((nodal, tolerance)) = velocity {
            for axis in 0..3 {
                // Scale before summing: a representable mean must not overflow
                // merely because the unscaled node sum does.
                let mean = vertices
                    .iter()
                    .map(|&vertex| nodal[vertex as usize][axis] * 0.25)
                    .sum::<f64>();
                let residual = mean - cell_velocity_m_s[cell][axis];
                if !residual.is_finite() || residual.abs() > tolerance {
                    return Err(Invalid(
                        "mechanical velocity differs from prescribed cell motion",
                    ));
                }
            }
        }
        selected.push((
            cell,
            vertices,
            problem
                .cell_volume_m3(cell)
                .ok_or(Invalid("Lorentz cell is outside its original mesh"))?,
            [force_density_n_m3[cell]; 4],
        ));
    }
    selected.sort_unstable_by_key(|entry| entry.0);
    if selected.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(Invalid("Lorentz load selection contains a repeated cell"));
    }
    integrate(
        problem.positions_m(),
        selected.iter().map(|entry| CellLoad {
            vertices: entry.1,
            volume_m3: entry.2,
            force_density_n_m3: &entry.3,
        }),
        limits,
        cancelled,
    )
}
