//! Explicit exposed-cell heat to the SAME gas used by the native instrument.
use super::*;
use fs_couple::air_path::input::convection::{
    AidingVerticalPlateGasBoundary, AidingVerticalPlateSettings,
    ForcedGasBoundary, ForcedGasGeometry, GasConvectionBoundary,
    HorizontalPlateFacing, NaturalConvectionGeometry, NaturalGasBoundary,
};
use fs_material::gas_species::input::temperature::IsobaricGasLaw;

#[derive(Clone, Debug)]
pub(super) struct Route {
    pub cell: CellId,
    pub boundary: GasConvectionBoundary,
}

pub(super) fn parse(rows: &mut Rows<'_>, cells: &[Cell]) -> Result<Vec<Route>> {
    let count = rows.count("convection", MAX_CELLS)?;
    if count == 0 {
        return Err("gas-reservoir boundary requires at least one exposed convection cell".into());
    }
    let mut routes = Vec::<Route>::with_capacity(count);
    for _ in 0..count {
        let mut row = rows.row("convection-cell")?;
        let cell = CellId(row.parse()?);
        let geometry = row.word()?;
        let exposed_area_m2 = row.positive()?;
        let characteristic_length_m = row.positive()?;
        // Forced/natural rows end with speed/gravity. The aiding mixed row
        // supplies signed upward speed followed by gravity and density bound.
        let driving = row.parse()?;
        let forced = |geometry| GasConvectionBoundary::Forced(ForcedGasBoundary {
            geometry, exposed_area_m2, characteristic_length_m, speed_m_s: driving,
        });
        let natural = |geometry| GasConvectionBoundary::Natural(NaturalGasBoundary {
            geometry, exposed_area_m2, characteristic_length_m, gravity_m_s2: driving,
        });
        let boundary = match geometry {
            "cylinder-crossflow" => forced(ForcedGasGeometry::CylinderCrossflow),
            "laminar-flat-plate" => forced(ForcedGasGeometry::LaminarFlatPlate),
            "turbulent-flat-plate" => forced(ForcedGasGeometry::TurbulentFlatPlate),
            "sphere" => forced(ForcedGasGeometry::Sphere),
            "natural-vertical-plate" => natural(NaturalConvectionGeometry::VerticalPlate),
            "natural-horizontal-cylinder" => natural(NaturalConvectionGeometry::HorizontalCylinder),
            "natural-sphere" => natural(NaturalConvectionGeometry::Sphere),
            "natural-horizontal-plate-up" => natural(NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up)),
            "natural-horizontal-plate-down" => natural(NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down)),
            "aiding-vertical-plate" => GasConvectionBoundary::AidingVerticalPlate(
                AidingVerticalPlateGasBoundary {
                    exposed_area_m2, characteristic_length_m,
                    settings: AidingVerticalPlateSettings {
                        speed_up_m_s: driving,
                        gravity_m_s2: row.parse()?,
                        maximum_relative_density_departure: row.parse()?,
                    },
                },
            ),
            _ => return Err("convection-cell needs an explicit forced, natural or aiding vertical-plate geometry".into()),
        };
        row.end()?;
        boundary.validate()?;
        let original = cells.iter().find(|c| c.id == cell)
            .ok_or("convection-cell must name an original finite cell")?;
        if routes.iter().any(|r| r.cell == cell)
            || boundary.exposed_area_m2() > original.geometry.surface_area_m2
        {
            return Err("convection cells must be unique and exposed area cannot exceed the original cell surface".into());
        }
        routes.push(Route { cell, boundary });
    }
    Ok(routes)
}

// Resolve source inventory and geometry once, rather than searching all cells
// inside every nonlinear trial at the mechanical sample rate.
struct BoundRoute {
    cell: CellId,
    boundary: GasConvectionBoundary,
    geometry: CellGeometry,
    material_id: MaterialId,
    mass_kg: f64,
}
pub(super) struct Reservoir {
    law: IsobaricGasLaw,
    routes: Vec<BoundRoute>,
}
impl Reservoir {
    pub fn prepare(profile: &Profile, ready: &Ready, sources: &[CellSource],
        gate: &CancelGate) -> Result<Option<Self>> {
        if profile.convection.is_empty() { return Ok(None); }
        poll(gate)?;
        if ready.spec.gas != *ready.spec.gas_record.state() {
            return Err("gas cooling cannot replace the original mechanical/acoustic atmosphere".into());
        }
        // No constant-property approximation of a thermochemical snapshot.
        // The generic boundary also accepts a retained thermochemical law;
        // this frontend currently retains source models only for named/legacy gas.
        let law = IsobaricGasLaw::try_from_record(ready.spec.gas_record.clone())?;
        if sources.len() != profile.cells.len() {
            return Err("gas cooling requires the complete original cell source inventory".into());
        }
        let mut routes = Vec::with_capacity(profile.convection.len());
        for route in &profile.convection {
            let index = profile.cells.iter().position(|cell| cell.id == route.cell)
                .ok_or("convection lost its original receiving cell")?;
            let cell = &profile.cells[index];
            let mass_kg = sources[index].density * cell.geometry.volume_m3;
            if !mass_kg.is_finite() || mass_kg <= 0.0 {
                return Err("gas cooling needs the original positive representable cell mass".into());
            }
            for t in cell.temperatures {
                poll(gate)?;
                route.boundary.evaluate(&law, t)?;
            }
            routes.push(BoundRoute { cell: route.cell, boundary: route.boundary,
                geometry: cell.geometry, material_id: MaterialId(index as u64), mass_kg });
        }
        Ok(Some(Self { law, routes }))
    }

    pub fn sources(&self, prepared: &Prepared<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<BTreeMap<CellId, f64>>
    {
        let mut sources = BTreeMap::new();
        for route in &self.routes {
            prepared.checkpoint()?;
            let cell = thermal.cells().get(&route.cell)
                .ok_or("gas boundary lost its original thermal receiver")?;
            if cell.id != route.cell || cell.segment_id != route.cell.0
                || cell.body_id != route.cell.0 || cell.material_id != route.material_id
                || cell.mass_kg.to_bits() != route.mass_kg.to_bits()
                || cell.geometry != route.geometry
                || cell.boundary.htc_w_m2_k != 0.0
                || cell.boundary.reservoir.is_some()
                || cell.boundary.source_w != 0.0
            {
                return Err("gas boundary requires its original inventory/geometry and no duplicate boundary/source account".into());
            }
            // Pure trial evaluation: reservoir heat is external supply, not
            // material friction or viscosity. No reservoir is mutated on a trial.
            let sample = route.boundary.evaluate(&self.law, thermal.temperature_k(route.cell)?)?;
            sources.insert(route.cell, sample.heat_into_wall_w);
        }
        prepared.checkpoint()?;
        Ok(sources)
    }

    pub fn describe(&self, routes: &[Route]) {
        let gas = self.law.reference().state();
        eprintln!("contact_thermal_gas_reservoir=T_K={}; p_Pa={}; exposed_cells={}; original composition; correlation-specific property temperatures; prescribed infinite reservoir; no flow/finite-gas integration",
            gas.temperature, gas.pressure, routes.len());
        for route in routes {
            match route.boundary {
                GasConvectionBoundary::Forced(boundary) => {
                    eprintln!("contact_thermal_convection_cell={}; geometry={:?}; area_m2={}; length_m={}; free_stream_m_s={}",
                        route.cell.0, boundary.geometry, boundary.exposed_area_m2,
                        boundary.characteristic_length_m, boundary.speed_m_s);
                }
                GasConvectionBoundary::Natural(boundary) => {
                    eprintln!("contact_thermal_convection_cell={}; natural_geometry={:?}; area_m2={}; length_m={}; gravity_m_s2={}",
                        route.cell.0, boundary.geometry, boundary.exposed_area_m2,
                        boundary.characteristic_length_m, boundary.gravity_m_s2);
                }
                GasConvectionBoundary::AidingVerticalPlate(boundary) => {
                    eprintln!("contact_thermal_convection_cell={}; aiding_vertical_plate; area_m2={}; height_m={}; speed_up_m_s={}; gravity_m_s2={}; maximum_relative_density_departure={}",
                        route.cell.0, boundary.exposed_area_m2, boundary.characteristic_length_m,
                        boundary.settings.speed_up_m_s, boundary.settings.gravity_m_s2,
                        boundary.settings.maximum_relative_density_departure);
                }
            }
        }
    }
}
