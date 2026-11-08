//! Quasistatic radial heat transport in free spherical liquid inventories.
//!
//! Ordered material cells carry fixed mass. Their ORIGINAL caloric rho(h)
//! determines current volume, cumulative radius and radial resistance at every
//! native explicit substep. This is a spherical, low-Mach kinematic closure:
//! radial expansion inertia, viscous work, mixing, buoyant overturning,
//! capillarity and compatible solid strain are not solved. Distinct liquid
//! layers must remain ordered and unmixed under the caller's model assumptions.
//!
//! The original charts must describe the maintained pressure/composition
//! convention. This adapter does not qualify that pressure from metadata,
//! infer an EOS, or claim internal energy or pressure work. Native storage
//! remains m*h, with its existing energy audit and timestep owner.

use super::*;
use fs_thermal::enthalpy::{
    DepositionObserver,
    conduction::{CurrentConductionGeometry, CurrentConductiveBoundary, CurrentConductiveLink},
};
use std::collections::BTreeMap;

/// Current physical material shell, distinct from immutable reference geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalLiquidCellGeometry {
    pub cell_id: CellId,
    pub mass_kg: f64,
    pub density_kg_m3: f64,
    pub volume_m3: f64,
    pub inner_radius_m: f64,
    pub outer_radius_m: f64,
    pub centroid_radius_m: f64,
}

/// Accepted liquid geometry in original cell order. No new material is created.
#[derive(Clone, Debug, PartialEq)]
pub struct SphericalLiquidGeometry {
    pub cells: Vec<SphericalLiquidCellGeometry>,
    pub total_mass_kg: f64,
    pub total_volume_m3: f64,
    pub outer_radius_m: f64,
    pub outer_area_m2: f64,
}

/// Density and radial geometry on the SAME native explicit conduction timeline.
///
/// Cloning retains the original immutable sources and copies all mutable state.
/// Use this directly with DepositionSimulation::advance_observed. It deliberately
/// exposes no implicit-chain conversion: freezing these density-derived metrics
/// during Newton trials would omit their enthalpy dependence.
#[derive(Clone)]
pub struct SphericalLiquidDensityObserver<'s, 'b, 'a> {
    bound: &'s BoundPhaseLayeredThermalBody<'b, 'a>,
    conduction: ConductionObserver<'a>,
    geometry: SphericalLiquidGeometry,
}

impl<'b, 'a> BoundPhaseLayeredThermalBody<'b, 'a> {
    /// Bind free radial LIQUID volume to the original source-backed heat owner.
    ///
    /// Requires a complete sphere with centre r=0. Cylinders, hollow shells,
    /// solids and mushy states need another kinematic/mechanical closure.
    /// Every initial cell must satisfy mass/rho(h) = reference volume within
    /// geometric roundoff; enabling this mode cannot introduce an unexplained
    /// initial geometry jump. Later reference cells and inventories stay fixed.
    ///
    /// Each actual substep refreshes spherical areas and resistance coordinates
    /// from the same accepted enthalpies as conductivity. The outer film/bath,
    /// when bound, follows the entire moving surface. Without that bath an
    /// authored outer ambient htc applies per CURRENT area. Interior ambient
    /// exposure and additional contacts refuse; authored total cell watts remain
    /// additive. Boundary-temperature/source events retain their native clock.
    ///
    /// Every accepted endpoint must remain fully liquid. A source, phase,
    /// topology or cancellation refusal rolls back the complete thermal and
    /// geometric candidate. Refine time and space independently; the native
    /// frozen-coefficient CFL is not a general nonlinear accuracy bound.
    pub fn spherical_liquid_density_observer<'s>(
        &'s self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SphericalLiquidDensityObserver<'s, 'b, 'a>, DepositionError> {
        checkpoint()?;
        if !matches!(self.body.geometry, LayerStackGeometry::Spherical { inner_radius_m } if inner_radius_m == 0.0)
        {
            return Err(bad(
                "free liquid density geometry requires a complete sphere with centre r=0",
            ));
        }
        for (&id, link) in thermal.links() {
            checkpoint()?;
            if (owns_cell(self, link.a) || owns_cell(self, link.b)) && !owns_link(self, id) {
                return Err(bad("free liquid sphere cannot inherit an unowned contact"));
            }
        }
        let geometry = resolve_geometry(self, thermal, true, &mut checkpoint)?;
        let conduction = self
            .conduction_observer(thermal, &mut checkpoint)?
            .with_current_geometry(
                thermal,
                transport_geometry(self, &geometry, &mut checkpoint)?,
                &mut checkpoint,
            )?;
        checkpoint()?;
        Ok(SphericalLiquidDensityObserver {
            bound: self,
            conduction,
            geometry,
        })
    }
}

impl SphericalLiquidDensityObserver<'_, '_, '_> {
    /// Current accepted geometry; the timeline still retains its reference cells.
    pub fn current_geometry(&self) -> &SphericalLiquidGeometry {
        &self.geometry
    }

    fn refresh_geometry(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        let geometry = resolve_geometry(self.bound, thermal, false, &mut || checkpoint())?;
        checkpoint()?;
        self.geometry = geometry;
        Ok(())
    }
}

fn owns_cell(bound: &BoundPhaseLayeredThermalBody<'_, '_>, id: CellId) -> bool {
    id.0 >= bound.ids.first_cell && id.0 - bound.ids.first_cell < bound.body.cells.len() as u64
}
fn owns_link(bound: &BoundPhaseLayeredThermalBody<'_, '_>, id: u64) -> bool {
    id >= bound.ids.first_link && id - bound.ids.first_link < bound.body.links.len() as u64
}

fn resolve_geometry(
    bound: &BoundPhaseLayeredThermalBody<'_, '_>,
    thermal: &DepositionSimulation<'_>,
    initial: bool,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<SphericalLiquidGeometry, DepositionError> {
    bound.admit_timeline(thermal, checkpoint)?;
    let mut cells = Vec::with_capacity(bound.body.cells.len());
    let (mut mass, mut volume, mut inner) = (0.0, 0.0, 0.0);
    for (index, reference) in bound.body.cells.iter().enumerate() {
        checkpoint()?;
        let id = bound.cell_id(index);
        let cell = &thermal.cells()[&id];
        let chart = bound.body.layers[reference.layer]
            .conductivity
            .caloric_curve();
        let state = chart
            .values_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)?;
        if state.liquid_mass_fraction != 1.0 {
            return Err(bad(
                "free spherical density transport requires fully liquid cells",
            ));
        }
        let last = index + 1 == bound.body.cells.len();
        if cell.boundary.htc_w_m2_k != 0.0 && (!last || bound.reservoirs[1].is_some()) {
            return Err(bad(
                "free sphere ambient exchange requires its unoccupied outer face",
            ));
        }
        let cell_volume = finite(cell.mass_kg / state.bulk_density_kg_m3)?;
        if !positive(cell_volume)
            || initial
                && (cell_volume - reference.geometry.volume_m3).abs()
                    > 128.0 * f64::EPSILON * cell_volume.max(reference.geometry.volume_m3)
        {
            return Err(bad(
                "liquid mass and original density must agree with reference cell volume",
            ));
        }
        let next_volume = finite(volume + cell_volume)?;
        if next_volume <= volume {
            return Err(bad("liquid shell volume is below geometric resolution"));
        }
        let outer = (next_volume / (4.0 * core::f64::consts::PI / 3.0)).cbrt();
        if !positive(outer) || outer <= inner {
            return Err(bad("liquid shell radii are not representable and ordered"));
        }
        let centroid = bound
            .body
            .geometry
            .centroid_m([inner, outer])
            .map_err(bad)?;
        mass = finite(mass + cell.mass_kg)?;
        cells.push(SphericalLiquidCellGeometry {
            cell_id: id,
            mass_kg: cell.mass_kg,
            density_kg_m3: state.bulk_density_kg_m3,
            volume_m3: cell_volume,
            inner_radius_m: inner,
            outer_radius_m: outer,
            centroid_radius_m: centroid,
        });
        volume = next_volume;
        inner = outer;
    }
    let area = bound.body.geometry.face_area_m2(inner).map_err(bad)?;
    checkpoint()?;
    Ok(SphericalLiquidGeometry {
        cells,
        total_mass_kg: mass,
        total_volume_m3: volume,
        outer_radius_m: inner,
        outer_area_m2: area,
    })
}

fn transport_geometry(
    bound: &BoundPhaseLayeredThermalBody<'_, '_>,
    current: &SphericalLiquidGeometry,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<CurrentConductionGeometry, DepositionError> {
    let mut metrics = CurrentConductionGeometry::default();
    for (index, pair) in current.cells.windows(2).enumerate() {
        checkpoint()?;
        let shape = LayerGeometry::spherical(
            [pair[0].centroid_radius_m, pair[1].centroid_radius_m],
            pair[0].outer_radius_m,
        )?;
        metrics.links.insert(
            bound.ids.first_link + index as u64,
            CurrentConductiveLink {
                area_m2: shape.contact_area_m2(),
                half_distances_m: shape.equivalent_distances_m(),
            },
        );
    }
    if let Some(face) = bound.reservoirs[1] {
        checkpoint()?;
        let cell = &current.cells[face.index];
        let distance = bound
            .body
            .geometry
            .resistance_length_at_radius_m(
                [cell.centroid_radius_m, cell.outer_radius_m],
                cell.outer_radius_m,
            )
            .map_err(bad)?;
        metrics.reservoirs.insert(
            cell.cell_id,
            CurrentConductiveBoundary {
                area_m2: current.outer_area_m2,
                half_distance_m: distance,
            },
        );
    }
    checkpoint()?;
    Ok(metrics)
}

impl DepositionObserver for SphericalLiquidDensityObserver<'_, '_, '_> {
    type Error = DepositionError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.conduction.begin()
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        // The native heat owner calls this for every actual starting state,
        // before both bath channels and its conductive timestep restriction.
        let geometry = resolve_geometry(self.bound, thermal, false, &mut || checkpoint())?;
        let metrics = transport_geometry(self.bound, &geometry, &mut || checkpoint())?;
        let mut conduction =
            self.conduction
                .with_current_geometry(thermal, metrics, &mut || checkpoint())?;
        let values = conduction.conductances(thermal, checkpoint)?;
        checkpoint()?;
        self.conduction = conduction;
        self.geometry = geometry;
        Ok(values)
    }

    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let mut values = self.conduction.ambient_conductances(thermal, checkpoint)?;
        if self.bound.reservoirs[1].is_none() {
            let id = self.bound.cell_id(self.bound.body.cells.len() - 1);
            let g = finite(thermal.cells()[&id].boundary.htc_w_m2_k * self.geometry.outer_area_m2)?;
            if !nonnegative(g) {
                return Err(DepositionError::Arithmetic);
            }
            values.insert(id, g);
        }
        checkpoint()?;
        Ok(values)
    }

    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.conduction.reservoir_conductances(thermal, checkpoint)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        if let Event::Contact(link) = event {
            if (owns_cell(self.bound, link.a) || owns_cell(self.bound, link.b))
                && !owns_link(self.bound, link.id)
            {
                return Err(bad("free liquid sphere cannot acquire an unowned contact"));
            }
        }
        self.conduction.event(event, thermal)?;
        self.refresh_geometry(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.conduction
            .step(start_s, duration_s, old_temperatures_k, thermal, checkpoint)?;
        self.refresh_geometry(thermal, checkpoint)
    }
}

#[cfg(test)]
mod tests;
