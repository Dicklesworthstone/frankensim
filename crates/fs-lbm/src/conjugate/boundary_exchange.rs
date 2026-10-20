//! Face-local reservoir exchange through the existing heat solver's cell sinks.
//! Includes the solid half-cell resistance; the returned surface is NOT the cell
//! centre. A nonlinear owner can re-evaluate convection/radiation there.

use std::collections::BTreeSet;
use fs_exec::CancelGate;
use super::{CellSink, ChtError, EnergySolution, ThermalFace, ThermalSetup, VoxelDomain, poll};

/// Frozen Robin law q'' = h (T_surface - T_reservoir) on one exterior solid face.
/// Source-backed laws and any composition of reservoirs remain caller-owned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryExchange {
    /// Existing solid cell; neither its identity nor occupancy is changed.
    pub cell: usize,
    /// Exterior face in x-, x+, y-, y+, z-, z+ order.
    pub face: usize,
    /// Finite positive surface coefficient, W/(m² K).
    pub coefficient_w_m2_k: f64,
    /// Finite positive absolute reservoir temperature, K.
    pub reservoir_temperature_k: f64,
}

/// Observed exchange using the exact coefficient supplied to the thermal solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryExchangeReport {
    /// Original declared face and frozen surface law.
    pub boundary: BoundaryExchange,
    /// Actual Cartesian face area, m².
    pub area_m2: f64,
    /// Surface temperature after the conductive half-cell drop, K.
    pub surface_temperature_k: f64,
    /// Signed heat leaving the cell, W; negative heats it from the reservoir.
    pub outward_w: f64,
}

#[derive(Debug, Clone)]
struct Entry {
    boundary: BoundaryExchange,
    area: f64,
    conductance: f64,
    surface_fraction: f64,
}

/// One frozen boundary image for a steady or implicit-storage heat solve.
/// Its sinks already include the external heat exchange: never add it again.
#[derive(Debug, Clone)]
pub struct PreparedBoundaryExchange {
    setup: ThermalSetup,
    entries: Vec<Entry>,
    cells: usize,
}

impl PreparedBoundaryExchange {
    /// Lower external Robin faces to the original conservative cell assembly.
    ///
    /// G = A / (w/(2k) + 1/h). The original domain face must be Adiabatic,
    /// explicitly delegating its exchange to this adapter. A fixed/convective/
    /// flux face cannot acquire a second owner. Distinct physical cell sinks
    /// and contact faces remain intact. Only exterior solid faces are admitted;
    /// internal occlusion, fluid radiation and compact components need their
    /// own boundary models. The input setup and conductivity are never changed.
    ///
    /// # Errors
    /// Bad geometry/fields, duplicate or occupied faces, nonphysical or
    /// unrepresentable coefficients, and cancellation refuse before publication.
    pub fn new(
        domain: &VoxelDomain, setup: &ThermalSetup, conductivity: &[[f64; 3]],
        boundaries: &[BoundaryExchange], gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        let cells = domain.cell_count();
        if conductivity.len() != cells || boundaries.len() > cells.saturating_mul(6) {
            return Err(invalid("conductivity shape or exterior-face budget mismatch"));
        }
        let mut occupied = BTreeSet::new();
        let mut next = setup.clone();
        let mut entries = Vec::with_capacity(boundaries.len());
        for (index, &boundary) in boundaries.iter().enumerate() {
            if index.is_multiple_of(256) { poll(gate)?; }
            let BoundaryExchange { cell, face, coefficient_w_m2_k: h,
                reservoir_temperature_k: ambient } = boundary;
            if cell >= cells || face >= 6 || domain.is_fluid(cell)
                || domain.neighbor(cell, face).is_some()
            { return Err(invalid("exchange needs an exterior solid cell face")); }
            if !occupied.insert((cell, face)) || setup.faces[face] != ThermalFace::Adiabatic
                || setup.fixed_temperature.iter().any(|&(c, _)| c == cell)
            { return Err(invalid("boundary face already has an exchange or fixed-temperature owner")); }
            let at = domain.coords(cell);
            if setup.compact_components.iter().any(|p|
                (0..3).all(|a| at[a] >= p.lo[a] && at[a] < p.hi[a]))
            { return Err(invalid("compact-component cells do not own ordinary face transport")); }
            let k = conductivity[cell][face / 2];
            if [h, ambient, k].iter().any(|v| !v.is_finite() || *v <= 0.0) {
                return Err(invalid("coefficient, conductivity and reservoir must be finite positive SI values"));
            }
            let half = 0.5 * domain.widths(cell)[face / 2] / k;
            let film = 1.0 / h;
            let resistance = half + film;
            let area = domain.face_area(cell, face / 2);
            let conductance = area / resistance;
            let surface_fraction = film / resistance;
            if [half, film, resistance, area, conductance, surface_fraction]
                .iter().any(|v| !v.is_finite() || *v <= 0.0)
                || !(conductance * ambient).is_finite()
            { return Err(invalid("half-cell or surface exchange is unrepresentable")); }
            next.cell_sinks.push(CellSink {
                cell, conductance_w_k: conductance, temperature_k: ambient,
            });
            entries.push(Entry { boundary, area, conductance, surface_fraction });
        }
        poll(gate)?;
        Ok(Self { setup: next, entries, cells })
    }

    /// The original setup plus the resolved boundary sinks, for either the
    /// steady energy solve or `EnergyConfig::step_with_conductivity`.
    #[must_use]
    pub const fn setup(&self) -> &ThermalSetup { &self.setup }

    /// Recover the surface states and signed heat from the resulting cell field.
    /// `solution` must be from a solve using this prepared setup and conductivity.
    /// Its sink accounting already includes these watts; this is a disaggregation,
    /// not another heat source. The caller retains any nonlinear mismatch audit.
    pub fn report(&self, solution: &EnergySolution, gate: &CancelGate)
        -> Result<Vec<BoundaryExchangeReport>, ChtError>
    {
        poll(gate)?;
        if solution.temperature.len() != self.cells {
            return Err(invalid("solution temperature shape differs from the prepared domain"));
        }
        let mut reports = Vec::with_capacity(self.entries.len());
        for (index, entry) in self.entries.iter().enumerate() {
            if index.is_multiple_of(256) { poll(gate)?; }
            let t = solution.temperature[entry.boundary.cell];
            let ambient = entry.boundary.reservoir_temperature_k;
            let outward_w = entry.conductance * (t - ambient);
            let surface_temperature_k = entry.surface_fraction.mul_add(t - ambient, ambient);
            if !t.is_finite() || t <= 0.0 || !outward_w.is_finite()
                || !surface_temperature_k.is_finite() || surface_temperature_k <= 0.0
            { return Err(invalid("surface temperature or heat is not representable")); }
            reports.push(BoundaryExchangeReport { boundary: entry.boundary,
                area_m2: entry.area, surface_temperature_k, outward_w });
        }
        poll(gate)?;
        Ok(reports)
    }
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput { field: "thermal.boundary_exchange", reason: reason.into() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::{EnergyConfig, FlowField, FluidProperties, SolidMaterial, Voxel};

    fn fluid() -> FluidProperties { FluidProperties { density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0, conductivity_w_m_k: 1.0, kinematic_viscosity_m2_s: 1.0 } }
    fn boundary(h: f64) -> BoundaryExchange { BoundaryExchange { cell: 0, face: 1,
        coefficient_w_m2_k: h, reservoir_temperature_k: 300.0 } }

    #[test]
    fn g1_half_cell_drop_and_both_heat_signs_match_analytic_robin_solution() {
        let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let gate = CancelGate::new();
        for power in [-4.0, 4.0] {
            let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
            setup.power_w = vec![power];
            let before = setup.clone();
            let prepared = PreparedBoundaryExchange::new(&domain, &setup, &[[2.0; 3]],
                &[boundary(4.0)], &gate).unwrap();
            let solution = EnergyConfig::default().solve_with_conductivity(&domain, &fluid(),
                &[SolidMaterial::new("synthetic", 2.0)], &FlowField::quiescent(&domain),
                prepared.setup(), &[[2.0; 3]], &gate).unwrap();
            let report = prepared.report(&solution, &gate).unwrap();
            assert!((solution.temperature[0] - (300.0 + power * 0.5)).abs() < 1e-9);
            assert!((report[0].surface_temperature_k - (300.0 + power * 0.25)).abs() < 1e-9);
            assert!((report[0].outward_w - power).abs() < 1e-9);
            assert!(solution.report.balance.residual_w.abs() < 1e-9);
            assert_eq!(setup, before);
        }
    }

    #[test]
    fn g1_storage_and_boundary_loss_close_without_counting_storage_as_reservoir_heat() {
        let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let gate = CancelGate::new();
        let prepared = PreparedBoundaryExchange::new(&domain,
            &ThermalSetup::new([ThermalFace::Adiabatic; 6]), &[[2.0; 3]],
            &[boundary(4.0)], &gate).unwrap();
        let (solution, record) = EnergyConfig::default().step_with_conductivity(&domain,
            &fluid(), &[SolidMaterial::new("synthetic", 2.0)], &FlowField::quiescent(&domain),
            prepared.setup(), &[[2.0; 3]], &[4.0], &[310.0], 0.5, &gate).unwrap();
        let report = prepared.report(&solution, &gate).unwrap();
        assert!((solution.temperature[0] - 308.0).abs() < 1e-9);
        assert!((report[0].outward_w - 16.0).abs() < 1e-9);
        assert!((record.stored_energy_change_j + 8.0).abs() < 1e-9);
        assert!(record.closure_j.abs() < 1e-9);
    }

    #[test]
    fn g0_duplicate_occupied_internal_and_fluid_faces_refuse() {
        let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let gate = CancelGate::new();
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        let mut face = boundary(4.0); face.face = 0;
        assert!(PreparedBoundaryExchange::new(&domain, &setup, &[[2.0; 3]; 2],
            &[face, face], &gate).is_err());
        setup.faces[0] = ThermalFace::Temperature(300.0);
        assert!(PreparedBoundaryExchange::new(&domain, &setup, &[[2.0; 3]; 2], &[face], &gate).is_err());
        setup.faces[0] = ThermalFace::Adiabatic;
        assert!(PreparedBoundaryExchange::new(&domain, &setup, &[[2.0; 3]; 2], &[boundary(4.0)], &gate).is_err());
        let liquid = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
        assert!(PreparedBoundaryExchange::new(&liquid, &setup, &[[2.0; 3]], &[face], &gate).is_err());
        for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            face.coefficient_w_m2_k = h;
            assert!(PreparedBoundaryExchange::new(&domain, &setup, &[[2.0; 3]; 2], &[face], &gate).is_err());
        }
    }
}
