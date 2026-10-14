//! SI heat and momentum use the SAME retained material evaluation and units.

use super::{MaterialSlab, ThermalLbm, ThermalLbm3};
use crate::thermal::sources::{
    HeatBalance, HeatSourceLimits, ThermalSourceBalance, ThermalSourceError,
    ThermalSourceLimits, source_step_mapped,
};

struct Inputs<const D: usize> {
    heat: Vec<f64>,
    force: Vec<[f64; D]>,
    limits: ThermalSourceLimits,
    energy_unit: f64,
    prescribed_j: f64,
    tolerance_j: f64,
}

fn invalid(what: &'static str) -> ThermalSourceError { ThermalSourceError::Invalid(what) }
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }

impl<L> MaterialSlab<L> {
    fn heat_inputs<const D: usize>(
        &self, power: &[f64], force: &[[f64; D]], cells: usize,
        depth_m: f64, limits: HeatSourceLimits,
    ) -> Result<Inputs<D>, ThermalSourceError> {
        if cells > limits.max_cells || power.len() != cells || force.len() != cells {
            return Err(invalid("heat/force field shape or cell budget does not match the material slab"));
        }
        let [t_min, t_max] = limits.temperature_range_k;
        if !positive(t_min) || !positive(t_max) || t_min > t_max || !positive(depth_m)
            || !limits.balance_tolerance_j.is_finite() || limits.balance_tolerance_j < 0.0
        {
            return Err(invalid("heat limits require positive absolute temperatures/depth and a finite nonnegative energy tolerance"));
        }
        let units = self.force_units;
        let rho = units.transport.density_kg_m3();
        let cp = units.transport.specific_heat_j_kg_k();
        let volume = (units.cell_size_m * units.cell_size_m) * depth_m;
        let capacity = (rho * volume) * cp;
        let energy_unit = capacity * self.temperature_difference_k;
        if !positive(volume) || !positive(capacity) || !energy_unit.is_finite() || energy_unit == 0.0 {
            return Err(invalid("material cell volume or thermal energy unit is unrepresentable"));
        }
        // Keep the sign of DeltaT. Heating remains heating when the top wall is
        // warmer and the existing normalized temperature coordinate is reversed.
        let heat = power.iter().map(|&p| {
            let delta = (((p / rho) / cp) * units.timestep_s) / self.temperature_difference_k;
            if !p.is_finite() || !delta.is_finite() || (p != 0.0 && delta == 0.0) {
                Err(invalid("physical heating rate is not representable in lattice temperature units"))
            } else { Ok(delta) }
        }).collect::<Result<Vec<_>, _>>()?;
        let mut prescribed_j = 0.0;
        for &p in power {
            let energy = (p * volume) * units.timestep_s;
            if !energy.is_finite() || (p != 0.0 && energy == 0.0) {
                return Err(invalid("prescribed physical source energy is unrepresentable"));
            }
            prescribed_j += energy;
        }
        if !prescribed_j.is_finite() { return Err(invalid("summed physical source energy overflows")); }
        let mut range = [
            (t_min - self.top_temperature_k) / self.temperature_difference_k,
            (t_max - self.top_temperature_k) / self.temperature_difference_k,
        ];
        if !range.iter().all(|t| t.is_finite()) { return Err(invalid("physical temperature bounds are unrepresentable")); }
        if range[0] > range[1] { range.swap(0, 1); }
        let tolerance = limits.balance_tolerance_j / energy_unit.abs();
        if !tolerance.is_finite() || (limits.balance_tolerance_j != 0.0 && tolerance == 0.0) {
            return Err(invalid("physical energy tolerance is unrepresentable"));
        }
        Ok(Inputs {
            heat, force: units.field(force, cells)?,
            limits: ThermalSourceLimits { max_cells: limits.max_cells, temperature_range: range, balance_tolerance: tolerance },
            energy_unit, prescribed_j, tolerance_j: limits.balance_tolerance_j,
        })
    }
}

impl<const D: usize> Inputs<D> {
    // Runs BEFORE the common source step publishes its tentative state. An SI
    // conversion failure or SI balance failure therefore cannot partially step.
    fn balance(&self, raw: ThermalSourceBalance) -> Result<HeatBalance, ThermalSourceError> {
        let stored_change_j = raw.stored_change * self.energy_unit;
        let bottom_wall_j = raw.bottom_wall * self.energy_unit;
        let top_wall_j = raw.top_wall * self.energy_unit;
        let defect_j = stored_change_j - self.prescribed_j - bottom_wall_j - top_wall_j;
        if ![stored_change_j, bottom_wall_j, top_wall_j, defect_j].iter().all(|x| x.is_finite()) {
            return Err(invalid("physical thermal balance is unrepresentable"));
        }
        if defect_j.abs() > self.tolerance_j {
            return Err(invalid("physical thermal balance exceeds the requested joule tolerance"));
        }
        Ok(HeatBalance { source_j: self.prescribed_j, stored_change_j, bottom_wall_j, top_wall_j, defect_j })
    }
}

impl MaterialSlab<ThermalLbm> {
    /// Advance heat and flow under volumetric power [W/m³] and external force
    /// density [N/m³]. The supplied out-of-plane depth [m] determines the
    /// represented energy, NOT the 2-D temperature or velocity evolution.
    ///
    /// One entry per cell, including both wall rows; wall power must be zero.
    /// Both loads are prescribed afresh each step. `rho`, `cp`, `dx` and `dt`
    /// come from the slab's retained material snapshot, not caller constants.
    /// Positive power heats, negative power cools, regardless of wall ordering.
    ///
    /// All source, force, temperature and joule-balance checks precede publish.
    /// This is sensible heat in the frozen-density/heat-capacity Boussinesq
    /// regime, NOT total-energy, Joule/combustion source generation, latent heat
    /// or pressure-volume work. The caller owns those source laws and validity.
    pub fn step_with_power_density(
        &mut self, power_w_m3: &[f64], force_n_m3: &[[f64; 2]],
        depth_m: f64, limits: HeatSourceLimits,
    ) -> Result<HeatBalance, ThermalSourceError> {
        if self.lattice.t_bottom != 1.0 || self.lattice.t_top != 0.0 {
            return Err(invalid("material heating refuses independently changed normalized wall temperatures"));
        }
        let input = self.heat_inputs(power_w_m3, force_n_m3, self.lattice.gpop.len(), depth_m, limits)?;
        source_step_mapped::<ThermalLbm, 2, 5, _>(
            &mut self.lattice, &input.heat, &input.force, input.limits,
            |balance| input.balance(balance),
        )
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional counterpart of the 2-D physical heating step.
    /// Cell volume is the retained `dx³`; no extrusion-depth assumption.
    /// Uses identical material conversion, source coupling and energy admission.
    pub fn step_with_power_density(
        &mut self, power_w_m3: &[f64], force_n_m3: &[[f64; 3]], limits: HeatSourceLimits,
    ) -> Result<HeatBalance, ThermalSourceError> {
        if self.lattice.t_bottom != 1.0 || self.lattice.t_top != 0.0 {
            return Err(invalid("material heating refuses independently changed normalized wall temperatures"));
        }
        let input = self.heat_inputs(power_w_m3, force_n_m3, self.lattice.gpop.len(), self.force_units.cell_size_m, limits)?;
        source_step_mapped::<ThermalLbm3, 3, 7, _>(
            &mut self.lattice, &input.heat, &input.force, input.limits,
            |balance| input.balance(balance),
        )
    }
}
