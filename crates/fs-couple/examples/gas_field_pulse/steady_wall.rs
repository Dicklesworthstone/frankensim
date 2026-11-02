//! Optional source-backed steady wall preparation for the existing DG example.
//! Parsing and pack ingress stay here; all heat, gas and acoustic laws keep
//! their original owners. No material number is supplied by this adapter.
use super::{Result, number, tag, word};
use fs_couple::material_source::thermal::{MaterialWallSource, ResolvedMaterialWallSource};
use fs_couple::material_source::{MAX_RECORD_AXES, MaterialRecord};
use fs_duct::{Duct, Segment};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::gas_field::AxialGasField;
use fs_material::gas_field::steady::wall::MaterialGasWall;
use fs_material::gas_field::steady::{
    AxialGasCell, SteadyAxialGas, SteadyGasSolution, SteadyGasSolve,
};
use fs_material::phase::SolidLiquidPhase;
use fs_material::thermal::conduction::{LayerConductionConfig, ReservoirLayerGeometry};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use std::{str::SplitAsciiWhitespace, sync::Arc};

#[derive(Debug)]
pub(super) struct Input {
    material: MaterialRecord,
    temperature_axis: String,
    phase: SolidLiquidPhase,
    property_names: [String; 3],
    source_span_k: [f64; 2],
    caloric_interpolation_tolerance_j_kg: f64,
    thickness_m: f64,
    film_resistance_m2_k_w: f64,
    reservoir_temperature_k: f64,
    mass_flow_kg_s: f64,
    gas_controls: SteadyGasSolve,
    wall_controls: LayerConductionConfig,
}

impl Input {
    /// Existing station-only inputs consume no extra token. The optional mode
    /// names every physical boundary and the independent numerical controls.
    pub(super) fn read_optional(words: &mut SplitAsciiWhitespace<'_>) -> Result<Option<Self>> {
        if words.clone().next() != Some("steady_wall") {
            return Ok(None);
        }
        tag(words, "steady_wall")?;
        tag(words, "wall_material")?;
        let pack = word(words, "wall material pack")?;
        let state = word(words, "wall material state")?;
        let axis_count: usize = number(words, "fixed material axis count")?;
        if axis_count > MAX_RECORD_AXES {
            return Err("wall material exceeds the existing fixed-axis budget".into());
        }
        let mut fields = vec![pack, state];
        for _ in 0..axis_count {
            fields.push(word(words, "fixed material AXIS=VALUE")?);
        }
        let material = MaterialRecord::parse(&fields)?;
        tag(words, "wall_source")?;
        let temperature_axis = word(words, "wall absolute-temperature axis")?.to_owned();
        if temperature_axis.is_empty()
            || temperature_axis.len() > 128
            || material
                .axes
                .iter()
                .any(|(axis, _)| axis == &temperature_axis)
        {
            return Err(
                "wall source owns one bounded temperature axis, absent from fixed coordinates"
                    .into(),
            );
        }
        let phase = match word(words, "wall pure phase")? {
            "solid" => SolidLiquidPhase::Solid,
            "liquid" => SolidLiquidPhase::Liquid,
            _ => return Err("wall source needs explicit solid or liquid phase".into()),
        };
        let property_names = [
            word(words, "wall density property")?.to_owned(),
            word(words, "wall specific heat capacity property")?.to_owned(),
            word(words, "wall thermal conductivity property")?.to_owned(),
        ];
        if property_names
            .iter()
            .any(|name| name.is_empty() || name.len() > 128)
        {
            return Err("wall property names must contain 1..=128 bytes".into());
        }
        let source_span_k = [
            number(words, "wall minimum kelvin")?,
            number(words, "wall maximum kelvin")?,
        ];
        let caloric_interpolation_tolerance_j_kg =
            number(words, "caloric interpolation tolerance J/kg")?;
        tag(words, "wall_geometry")?;
        let thickness_m = number(words, "wall thickness m")?;
        let film_resistance_m2_k_w = number(words, "areal gas-film resistance m2 K/W")?;
        tag(words, "wall_reservoir")?;
        let reservoir_temperature_k = number(words, "outer reservoir kelvin")?;
        tag(words, "wall_flow")?;
        let mass_flow_kg_s = number(words, "prescribed gas mass rate kg/s")?;
        tag(words, "gas_solve")?;
        let gas_controls = SteadyGasSolve {
            caloric: GasEnergySolve {
                minimum_temperature_k: number(words, "minimum gas kelvin")?,
                maximum_temperature_k: number(words, "maximum gas kelvin")?,
                energy_tolerance_j_kg: number(words, "gas caloric tolerance J/kg")?,
                maximum_iterations: number(words, "gas caloric iterations")?,
            },
            maximum_dt_s: number(words, "maximum thermal pseudo-step s")?,
            maximum_relaxations: number(words, "maximum thermal relaxations")?,
            maximum_substeps: number(words, "maximum thermal substeps")?,
            steady_power_tolerance_w: number(words, "steady power tolerance W")?,
            energy_tolerance_j: number(words, "thermal energy tolerance J")?,
            maximum_mean_mach: number(words, "maximum mean Mach")?,
        };
        tag(words, "wall_solve")?;
        let wall_controls = LayerConductionConfig {
            max_iterations: number(words, "wall solve iterations")?,
            absolute_heat_tolerance_w: number(words, "wall absolute heat tolerance W")?,
            relative_heat_tolerance: number(words, "wall relative heat tolerance")?,
        };
        let result = Self {
            material,
            temperature_axis,
            phase,
            property_names,
            source_span_k,
            caloric_interpolation_tolerance_j_kg,
            thickness_m,
            film_resistance_m2_k_w,
            reservoir_temperature_k,
            mass_flow_kg_s,
            gas_controls,
            wall_controls,
        };
        result.validate()?;
        Ok(Some(result))
    }

    fn validate(&self) -> Result<()> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        let nonnegative = |x: f64| x.is_finite() && x >= 0.0;
        let gas = self.gas_controls;
        if !self.source_span_k.into_iter().all(positive)
            || self.source_span_k[0] >= self.source_span_k[1]
            || !positive(self.caloric_interpolation_tolerance_j_kg)
            || !positive(self.thickness_m)
            || !nonnegative(self.film_resistance_m2_k_w)
            || !positive(self.reservoir_temperature_k)
            || !positive(self.mass_flow_kg_s)
            || !positive(gas.caloric.minimum_temperature_k)
            || !positive(gas.caloric.maximum_temperature_k)
            || gas.caloric.minimum_temperature_k >= gas.caloric.maximum_temperature_k
            || !positive(gas.caloric.energy_tolerance_j_kg)
            || !(1..=128).contains(&gas.caloric.maximum_iterations)
            || !positive(gas.maximum_dt_s)
            || !(1..=1_000_000).contains(&gas.maximum_relaxations)
            || !(1..=1_000_000).contains(&gas.maximum_substeps)
            || !positive(gas.steady_power_tolerance_w)
            || !nonnegative(gas.energy_tolerance_j)
            || !positive(gas.maximum_mean_mach)
            || gas.maximum_mean_mach > 0.1
            || !(1..=1024).contains(&self.wall_controls.max_iterations)
            || !nonnegative(self.wall_controls.absolute_heat_tolerance_w)
            || !nonnegative(self.wall_controls.relative_heat_tolerance)
        {
            return Err(
                "steady-wall data need finite physical scales and bounded solve controls".into(),
            );
        }
        Ok(())
    }

    /// Stations define the supplied inlet and each finite volume's initial
    /// temperature. No downstream composition or target temperature is fitted.
    pub(super) fn solve(
        &self,
        initial: &AxialGasField,
        duct: &Duct,
        mut load: impl FnMut(&str) -> std::result::Result<NormalizedMaterialCardPack, String>,
    ) -> Result<(SteadyGasSolution, String)> {
        self.validate()?;
        let (radius, length) = match duct.segments.as_slice() {
            [Segment::Cylinder { radius, length }] => (*radius, *length),
            _ => return Err("steady_wall currently requires one cylindrical bore".into()),
        };
        if initial.length_m() != length {
            return Err("steady_wall station endpoints must cover the exact bore length".into());
        }
        let states = initial.primitive_states();
        if states
            .iter()
            .any(|state| state.composition != states[0].composition)
        {
            return Err(
                "steady_wall has one inlet and no injections; station compositions must agree"
                    .into(),
            );
        }
        let source = ResolvedMaterialWallSource::new(
            MaterialWallSource {
                material: self.material.clone(),
                temperature_axis: self.temperature_axis.clone(),
                phase: self.phase,
                property_names: self.property_names.clone(),
                temperature_span_k: self.source_span_k,
                caloric_interpolation_tolerance_j_kg: self.caloric_interpolation_tolerance_j_kg,
            },
            Arc::new(load(&self.material.pack)?),
        )?;
        let layer = source.layer(|| Ok(()))?;
        let outer_radius = radius + self.thickness_m;
        let mut cells = Vec::with_capacity(states.len() - 1);
        let mut walls = Vec::with_capacity(states.len() - 1);
        for (index, pair) in initial.stations_m().windows(2).enumerate() {
            let dx = pair[1] - pair[0];
            let geometry = ReservoirLayerGeometry::cylindrical(radius, outer_radius, dx)?;
            let wall_area_m2 = geometry.surface_area_m2();
            cells.push(AxialGasCell {
                end_m: pair[1],
                volume_m3: core::f64::consts::PI * radius * radius * dx,
                wall_area_m2,
                wall_conductance_w_k: 0.0,
                wall_temperature_k: self.reservoir_temperature_k,
                initial_temperature_k: states[index + 1].temperature_k,
            });
            // Areal resistance keeps the physical film fixed under cell refinement.
            walls.push(MaterialGasWall::new(
                index,
                &layer,
                geometry,
                self.film_resistance_m2_k_w / wall_area_m2,
            )?);
        }
        let problem = SteadyAxialGas {
            source: initial.model().clone(),
            source_conditions: QueryPoint::new(),
            static_pressure_pa: initial.static_pressure_pa(),
            phase: initial.phase(),
            inlet: states[0].clone(),
            inlet_mass_flow_kg_s: self.mass_flow_kg_s,
            cells,
            injections: vec![],
        };
        let solution = problem.solve_with_material_walls(
            &walls,
            self.wall_controls,
            self.gas_controls,
            || false,
        )?;
        Ok((solution, self.material.provenance(source.pack())))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "steady_wall
        wall_material wall.fsmatpk first-heating 2 first-heating=1 source-pressure-known=0
        wall_source temperature solid density specific-heat-capacity thermal-conductivity 298.15 350 0.0001
        wall_geometry 0.003 0.01 wall_reservoir 300 wall_flow 0.0001
        gas_solve 298.15 350 0.0000001 64 0.01 2000 100000 0.00000001 0.00000001 0.05
        wall_solve 96 0.00000000001 0.000000000001";

    #[test]
    fn explicit_mode_preserves_next_acoustic_record_and_source_coordinates() {
        let text = format!("{INPUT} spatial_viscothermal 2 50 2000 0.01 4");
        let mut words = text.split_ascii_whitespace();
        let input = Input::read_optional(&mut words).unwrap().unwrap();
        assert_eq!(input.material.state, "first-heating");
        assert_eq!(input.material.axes.len(), 2);
        assert_eq!(input.temperature_axis, "temperature");
        assert_eq!(input.film_resistance_m2_k_w, 0.01);
        assert_eq!(words.next(), Some("spatial_viscothermal"));
        let mut original = "viscothermal 2 50 2000 0.01".split_ascii_whitespace();
        assert!(Input::read_optional(&mut original).unwrap().is_none());
        assert_eq!(original.next(), Some("viscothermal"));
    }

    #[test]
    fn invalid_sources_boundaries_and_work_are_refused_before_loading() {
        for text in [
            INPUT.replace("wall.fsmatpk", "../wall.fsmatpk"),
            INPUT.replace("first-heating=1", "temperature=300"),
            INPUT.replace(
                "thermal-conductivity 298.15 350",
                "thermal-conductivity 350 298.15",
            ),
            INPUT.replace("wall_geometry 0.003 0.01", "wall_geometry 0.003 -1"),
            INPUT.replace("wall_flow 0.0001", "wall_flow NaN"),
            INPUT.replace("wall_solve 96", "wall_solve 0"),
            INPUT.replace("2000 100000", "2000 1000001"),
        ] {
            assert!(
                Input::read_optional(&mut text.split_ascii_whitespace()).is_err(),
                "{text}"
            );
        }
    }
}
