//! One explicit baked load after all regional or continuous propagation.

use super::*;
use crate::driving_point::{LineDc, LineLadder, PASSIVITY_TOL};
use fs_duct::{
    TabulatedLoad,
    compose::{ComposeOptions, compose_chain_reflection},
};
use fs_math::c64::C64;

impl RegionalDuct {
    pub(super) fn loaded_reflection(
        &self,
        omega: f64,
        load: C64,
        wall: WallModel<'_>,
    ) -> Result<C64, DuctError> {
        self.admit_wall_override(wall)?;
        let zc = self.inlet_impedance();
        let exterior = self.exterior_gas();
        if let Some(field) = &self.continuous {
            let chain =
                field.chain_with_wall_exterior_controlled(omega, LossModel::Bessel,
                    Some(exterior), wall, &mut || false)?;
            return compose_chain_reflection(&chain.blocks, Some(load), ComposeOptions::default())?
                .referenced_to(zc);
        }
        let mut runs = self.regions.iter().rev();
        let terminal = runs.next().expect("admitted nonempty duct");
        let mut response =
            fs_duct::exterior::input_impedance_load(&terminal.duct, &terminal.gas,
                exterior,
                omega, LossModel::Bessel, load, Self::region_wall(terminal, wall))?;
        for region in runs {
            response = fs_duct::exterior::input_impedance_load(
                &region.duct,
                &region.gas,
                exterior,
                omega,
                LossModel::Bessel,
                response.impedance,
                Self::region_wall(region, wall),
            )?;
        }
        Ok(impedance_reflectance(&response, zc))
    }

    pub(super) fn tabulated_line(
        &self,
        table: &TabulatedLoad,
        rate: u32,
        n: usize,
        wall: WallModel<'_>,
        n_fft: Option<usize>,
    ) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
        let length = self
            .regions
            .iter()
            .map(|region| super::super::duct_axial_length(&region.duct))
            .sum();
        let radius = self
            .regions
            .last()
            .expect("admitted duct")
            .duct
            .segments
            .last()
            .expect("admitted outlet")
            .outlet_radius();
        let mouth = table
            .validate_mouth(length, radius, self.exterior_gas())
            .map_err(DrivingPointError::Duct)?;
        let band_n = n_fft.map_or(LineLadder::default().max_taps, |n| {
            n.clamp(256, 8192).next_power_of_two()
        });
        super::super::baked::validate_band(table, rate, band_n)?;
        // Keep the uniform path identical, including its arithmetic and refusals.
        if self.continuous.is_none()
            && !self.has_material_wall_regions()
            && self.regions.len() == 1
            && self.exterior_gas() == self.outlet_gas()
            && let WallModel::Rigid | WallModel::Pin(_) = wall
        {
            let pin = match wall {
                WallModel::Pin(pin) => Some(pin),
                _ => None,
            };
            let region = &self.regions[0];
            return super::super::baked::characteristic(
                &region.duct,
                &region.gas,
                table,
                rate,
                n,
                self.inlet_impedance(),
                pin,
                n_fft,
                None,
            );
        }
        let delay = 2.0 * self.one_way_time_s * f64::from(rate);
        let dc = self
            .loaded_reflection(1.0e-3, C64::new(mouth.dc_impedance, 0.0), wall)
            .map_err(DrivingPointError::Duct)?
            .re;
        characteristic_line_response_with_poles(
            delay,
            LineDc::Explicit(dc),
            rate,
            n,
            self.inlet_impedance(),
            n_fft,
            vec![bulk_delay_pole_state(delay)],
            |omega| {
                let r = self.loaded_reflection(omega, table.z_at(omega)?, wall)?;
                if !r.re.is_finite() || !r.im.is_finite() || r.abs() > 1.0 + PASSIVITY_TOL {
                    return Err(invalid(
                        "tabulated regional response is not finite and passive",
                    ));
                }
                Ok(r)
            },
        )
    }
}
