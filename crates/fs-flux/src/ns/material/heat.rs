//! Physical production integrated from retained original-stage powers.

use super::{MaterialFluxError, MaterialStokesSolution};

/// Navier-slip heat generated on one original boundary edge.
/// The thermal owner must assign any share to the fluid or wall explicitly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialBoundarySlipHeat {
    /// Original mesh boundary edge.
    pub edge: usize,
    /// Original connected fluid component.
    pub component: usize,
    /// Nonnegative generated heat per metre of out-of-plane span [J/m].
    /// Zero on non-slip and perfect-slip edges.
    pub energy_j_per_m: f64,
}

/// Physical molecular and Navier-slip production over accepted elapsed time.
///
/// Cell entries follow the original triangle order; boundary entries preserve
/// original edge/component identities. Positive SDIRK2 weights integrate the
/// SAME retained stage powers used by material flow and boundary reactions.
/// No source callback, constitutive fit or heating kernel is evaluated here.
///
/// This is generated heat, not automatic thermal deposition. Slip heat has no
/// inferred wall/fluid partition. Mechanical wall work, pressure work, kinetic
/// transport, numerical upwind/SIP damping and caloric evolution are excluded.
/// These integrals do not assert a full discrete energy balance or passivity.
/// Their temporal accuracy requires resolving the source and strain history;
/// a velocity-only adaptive tolerance does not certify heat accuracy.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialDissipatedHeat {
    duration_s: f64,
    cell_bulk_energy_j_per_m: Vec<f64>,
    boundary_slip_energy: Vec<MaterialBoundarySlipHeat>,
    bulk_energy_j_per_m: f64,
    slip_energy_j_per_m: f64,
}

impl MaterialDissipatedHeat {
    /// Elapsed physical interval represented by these integrals [s].
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Integral of original molecular power in each original cell [J/m].
    /// Divide by that original cell's area [m²] for integrated density [J/m³].
    #[must_use]
    pub fn cell_bulk_energy_j_per_m(&self) -> &[f64] {
        &self.cell_bulk_energy_j_per_m
    }

    /// Original boundary-edge order, including zero-power non-slip edges.
    #[must_use]
    pub fn boundary_slip_energy(&self) -> &[MaterialBoundarySlipHeat] {
        &self.boundary_slip_energy
    }

    /// Compensated sum of cell molecular production [J/m].
    #[must_use]
    pub const fn bulk_energy_j_per_m(&self) -> f64 {
        self.bulk_energy_j_per_m
    }

    /// Compensated sum of boundary slip production [J/m], unpartitioned.
    #[must_use]
    pub const fn slip_energy_j_per_m(&self) -> f64 {
        self.slip_energy_j_per_m
    }
}

// Internal to the original material/time owners: callers cannot attach a new
// mesh or invented interval to an already solved stage. Storage remains O(mesh)
// across any number of accepted adaptive steps.
pub(super) struct DissipatedHeatAccumulator {
    cells: Vec<Sum>,
    edges: Vec<(usize, usize, Sum)>,
}

impl DissipatedHeatAccumulator {
    pub(super) fn new(cells: usize) -> Self {
        Self {
            cells: vec![Sum::default(); cells],
            edges: Vec::new(),
        }
    }

    pub(super) fn add_stage(
        &mut self,
        areas_m2: &[f64],
        solution: &MaterialStokesSolution,
        weighted_time_s: f64,
    ) -> Result<(), MaterialFluxError> {
        let heating = solution.viscous_dissipation().cell_heating_w_m3();
        if areas_m2.len() != self.cells.len()
            || heating.len() != self.cells.len()
            || !weighted_time_s.is_finite()
            || weighted_time_s <= 0.0
        {
            return Err(invalid("incompatible original stage heat geometry or time"));
        }
        for ((sum, &area), &rate) in self.cells.iter_mut().zip(areas_m2).zip(heating) {
            if !area.is_finite() || area <= 0.0 {
                return Err(invalid("finite positive original heat cell area required"));
            }
            // This is the original cell-power conversion, not another strain
            // evaluation. The physical flow owner already checked rate*area.
            sum.add(energy(rate * area, weighted_time_s)?)?;
        }
        let powers = solution.boundary_slip_dissipation();
        let first = self.edges.is_empty();
        if !first && powers.len() != self.edges.len() {
            return Err(invalid("original heat boundary topology changed"));
        }
        for (index, item) in powers.enumerate() {
            let power = item?;
            if first {
                self.edges
                    .push((power.edge, power.component, Sum::default()));
            }
            let (edge, component, sum) = &mut self.edges[index];
            if (*edge, *component) != (power.edge, power.component) {
                return Err(invalid("original heat boundary identity changed"));
            }
            sum.add(energy(power.power_w_per_m, weighted_time_s)?)?;
        }
        Ok(())
    }

    pub(super) fn add_integral(
        &mut self,
        heat: &MaterialDissipatedHeat,
    ) -> Result<(), MaterialFluxError> {
        if heat.cell_bulk_energy_j_per_m.len() != self.cells.len() {
            return Err(invalid("original heat cell topology changed"));
        }
        for (sum, &value) in self.cells.iter_mut().zip(&heat.cell_bulk_energy_j_per_m) {
            sum.add(value)?;
        }
        let first = self.edges.is_empty();
        if !first && heat.boundary_slip_energy.len() != self.edges.len() {
            return Err(invalid("original heat boundary topology changed"));
        }
        for (index, heat) in heat.boundary_slip_energy.iter().enumerate() {
            if first {
                self.edges.push((heat.edge, heat.component, Sum::default()));
            }
            let (edge, component, sum) = &mut self.edges[index];
            if (*edge, *component) != (heat.edge, heat.component) {
                return Err(invalid("original heat boundary identity changed"));
            }
            sum.add(heat.energy_j_per_m)?;
        }
        Ok(())
    }

    pub(super) fn finish(
        self,
        duration_s: f64,
    ) -> Result<MaterialDissipatedHeat, MaterialFluxError> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(invalid("finite positive original heat interval required"));
        }
        let cell_bulk_energy_j_per_m: Vec<_> =
            self.cells.into_iter().map(|sum| sum.value).collect();
        let boundary_slip_energy: Vec<_> = self
            .edges
            .into_iter()
            .map(|(edge, component, sum)| MaterialBoundarySlipHeat {
                edge,
                component,
                energy_j_per_m: sum.value,
            })
            .collect();
        let mut bulk = Sum::default();
        let mut slip = Sum::default();
        for &value in &cell_bulk_energy_j_per_m {
            bulk.add(value)?;
        }
        for heat in &boundary_slip_energy {
            slip.add(heat.energy_j_per_m)?;
        }
        Ok(MaterialDissipatedHeat {
            duration_s,
            cell_bulk_energy_j_per_m,
            boundary_slip_energy,
            bulk_energy_j_per_m: bulk.value,
            slip_energy_j_per_m: slip.value,
        })
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    compensation: f64,
}

impl Sum {
    fn add(&mut self, value: f64) -> Result<(), MaterialFluxError> {
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("nonfinite or negative physical heat contribution"));
        }
        let corrected = value - self.compensation;
        let next = self.value + corrected;
        if !next.is_finite() || next < 0.0 {
            return Err(invalid("unrepresentable accumulated physical heat"));
        }
        self.compensation = (next - self.value) - corrected;
        self.value = next;
        Ok(())
    }
}

fn energy(power_w_per_m: f64, weighted_time_s: f64) -> Result<f64, MaterialFluxError> {
    let value = power_w_per_m * weighted_time_s;
    if !power_w_per_m.is_finite()
        || power_w_per_m < 0.0
        || !value.is_finite()
        || (power_w_per_m > 0.0 && value == 0.0)
    {
        return Err(invalid("unrepresentable original stage heat integral"));
    }
    Ok(value)
}

fn invalid(message: &'static str) -> MaterialFluxError {
    MaterialFluxError::Invalid(message)
}
