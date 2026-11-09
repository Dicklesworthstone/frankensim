//! Held boundary inputs on the original open-grid thermal carrier.
//!
//! Setters change the next interval's forcing without advancing heat or
//! replacing its accepted history. Admission completes before any input is
//! changed. The ordinary step still owns resolved traces, energy accounting,
//! cancellation and joint publication with native momentum.

use super::conduction::boundary_cell;
use super::{
    ThermalOpenError3, ThermalOpenFace3, ThermalOpenField3, ThermalOpenLimits3, ThermalWall,
    check_boundary_temperature, invalid,
};
use crate::thermal::boundary::ThermalWallField;

fn admit_stateless(
    face: usize,
    wall: ThermalWall,
    limits: ThermalOpenLimits3,
) -> Result<(), ThermalOpenError3> {
    if matches!(wall, ThermalWall::FiniteSolid(_)) {
        return Err(invalid(
            "boundary input updates cannot create or replace finite solid storage",
        ));
    }
    wall.validate().map_err(ThermalOpenError3::Wall)?;
    if let ThermalWall::Temperature(temperature) = wall {
        check_boundary_temperature(face, temperature, limits)?;
    }
    Ok(())
}

impl ThermalOpenField3 {
    /// Replace the held law on an existing stateless external wall.
    ///
    /// `face` uses XMin/XMax/YMin/YMax/ZMin/ZMax order. Prescribed temperature,
    /// signed inward heat flux, ambient film data and insulation may change
    /// between accepted steps without resetting scalar populations, flux
    /// history or time. The law is uniform over the original face; existing
    /// native geometry and tangential wall motion remain unchanged.
    ///
    /// A finite wall cannot be replaced, and this setter cannot create finite
    /// storage. Use [`Self::set_finite_wall_heat`] to drive its retained cells.
    /// Film ambient temperature is source data; its actual fluid-side trace,
    /// rather than the ambient, must satisfy the carrier's temperature range.
    /// Validation failure leaves every input and physical field untouched.
    pub fn set_wall_boundary(
        &mut self,
        face: usize,
        wall: ThermalWall,
    ) -> Result<(), ThermalOpenError3> {
        if !matches!(
            self.boundaries.faces.get(face),
            Some(ThermalOpenFace3::Wall(_))
        ) {
            return Err(invalid("thermal target is not an external wall face"));
        }
        if self.walls[face]
            .iter()
            .any(|wall| matches!(wall, ThermalWall::FiniteSolid(_)))
        {
            return Err(invalid(
                "finite wall storage cannot be replaced by a boundary input",
            ));
        }
        admit_stateless(face, wall, self.limits)?;
        // The next begin_step copies these accepted inputs into its wall
        // scratch. No candidate reflection or accepted temperature is reused.
        self.walls[face].fill(wall);
        self.boundaries.faces[face] = ThermalOpenFace3::Wall(wall);
        Ok(())
    }

    /// Set a held spatial field on an existing stateless external wall.
    ///
    /// Spatial entries cover the full original x-fast plane: X faces use
    /// y+ny*z, Y faces x+nx*z, and Z faces x+nx*y. Each entry may independently
    /// prescribe temperature, inward heat flux, an ambient film or insulation.
    /// The original reflected pair resolves its surface temperature and actual
    /// heat exchange; no fluid cell temperature or wall heat is prescribed twice.
    ///
    /// Every value is admitted before any slot changes. A nonzero explicitly
    /// authored heat flux on an inactive solid-facing slot refuses. Temperatures
    /// and ambient film data in inactive slots still undergo ordinary admission.
    /// Finite storage cannot be created or replaced; its retained heat sources
    /// use set_finite_wall_heat instead.
    ///
    /// Spatial input reuses the already budgeted wall vectors and performs no
    /// allocation. wall_states returns the actual field; boundaries retains the
    /// last uniform declaration for face ownership. A later uniform setter
    /// replaces the entire field. Scalar/flux history, clocks and receipts stay
    /// unchanged until the original joint flow/heat step succeeds.
    pub fn set_wall_boundary_field(
        &mut self,
        face: usize,
        field: ThermalWallField<'_>,
    ) -> Result<(), ThermalOpenError3> {
        let ThermalWallField::Spatial(walls) = field else {
            return self.set_wall_boundary(face, field.at(0));
        };
        if !matches!(
            self.boundaries.faces.get(face),
            Some(ThermalOpenFace3::Wall(_))
        ) {
            return Err(invalid("thermal field target is not an external wall face"));
        }
        if walls.len() != self.walls[face].len() {
            return Err(invalid("thermal wall field shape mismatch"));
        }
        if self.walls[face]
            .iter()
            .any(|wall| matches!(wall, ThermalWall::FiniteSolid(_)))
        {
            return Err(invalid(
                "finite wall storage cannot be replaced by a boundary field",
            ));
        }
        for (link, &wall) in walls.iter().enumerate() {
            admit_stateless(face, wall, self.limits)?;
            if let ThermalWall::HeatFlux(heat) = wall
                && heat != 0.0
                && self.solid[boundary_cell(self.dims, face / 2, face % 2, link)]
            {
                return Err(invalid(
                    "spatial wall heat flux must be zero on inactive slots",
                ));
            }
        }
        self.walls[face].copy_from_slice(walls);
        Ok(())
    }

    /// Replace the held stateless law on unassigned immersed-solid links.
    ///
    /// Explicit immersed bodies and their conduction graphs retain ownership
    /// of every assigned link; this is only the fallback boundary condition.
    /// Prescribed temperatures, inward fluxes, ambient films and insulation
    /// retain the same units and trace checks as external walls. Finite storage
    /// must continue to use the original immersed-body installation APIs.
    pub fn set_internal_wall_boundary(
        &mut self,
        wall: ThermalWall,
    ) -> Result<(), ThermalOpenError3> {
        admit_stateless(6, wall, self.limits)?;
        self.boundaries.internal_wall = wall;
        Ok(())
    }

    /// Change external heat into each retained finite cell of one wall face.
    ///
    /// The field covers the complete original x-fast wall array. Each finite
    /// value is heat into THIS solid cell per lattice step, in the same scale
    /// as `FiniteThermalWall::heat_per_step`; cooling is negative. Slots next
    /// to solid grid cells are inactive and must receive zero. No inactive
    /// source is silently dropped or redistributed.
    ///
    /// Stored temperatures, heat capacities, contact resistances and prepared
    /// wall-conduction factors remain intact. Changing power adds no energy at
    /// setter time: the original finite-wall balance credits it once on the
    /// next successful step. Current values are exposed by `wall_states`;
    /// `boundaries` retains the original finite-wall declaration.
    ///
    /// Every entry is admitted before any source changes. A later refused or
    /// cancelled step retains these authored inputs and the accepted physical
    /// state, as it does for a changed reservoir temperature. This setter and
    /// subsequent stepping allocate no additional wall arrays.
    pub fn set_finite_wall_heat(
        &mut self,
        face: usize,
        heat_per_step: &[f64],
    ) -> Result<(), ThermalOpenError3> {
        if !matches!(
            self.boundaries.faces.get(face),
            Some(ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(_)))
        ) {
            return Err(invalid("thermal target is not a finite external wall"));
        }
        if heat_per_step.len() != self.walls[face].len() {
            return Err(invalid("finite wall source field shape mismatch"));
        }
        let axis = face / 2;
        let side = face % 2;
        for (link, (&heat, wall)) in heat_per_step.iter().zip(&self.walls[face]).enumerate() {
            if !heat.is_finite()
                || (self.solid[boundary_cell(self.dims, axis, side, link)] && heat != 0.0)
            {
                return Err(invalid(
                    "finite wall heat must be finite and zero on inactive slots",
                ));
            }
            if !matches!(wall, ThermalWall::FiniteSolid(_)) {
                return Err(invalid(
                    "finite wall source target has no retained solid cell",
                ));
            }
        }
        for (wall, &heat) in self.walls[face].iter_mut().zip(heat_per_step) {
            if let ThermalWall::FiniteSolid(solid) = wall {
                solid.heat_per_step = heat;
            }
        }
        Ok(())
    }
}
