//! Data-only fiber geometry and material selection for the existing renderer.
use super::*;
pub use fs_couple::fiber_ribbon::strain::FiberStrainLimits;
use fs_couple::fiber_ribbon::{FiberRibbonGeometry, FiberRibbonResolution, TensionedFiberRibbon};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, axial::YOUNG_VISCOSITY_PROPERTY,
};

#[derive(Clone, Debug)]
pub struct Input {
    pub pack: String,
    pub radius_m: f64,
    pub length_m: f64,
    pub fiber_count: usize,
    pub tension_n: f64,
    pub station_fraction: f64,
    pub axial_modes: usize,
    pub transverse_modes: usize,
    pub strain: Option<(FiberStrainLimits, NonlinearModalConfig)>,
}
impl Input {
    pub fn mode_count(&self, spec: &input::Spec) -> Result<usize> {
        let (_, cells) = spec
            .bow_ribbon
            .ok_or("bow_fibers requires bow_ribbon width/cells")?;
        if self.axial_modes == 0
            || self.transverse_modes == 0
            || self.fiber_count == 0
            || !(0.0 < self.station_fraction && self.station_fraction < 1.0)
        {
            return Err(
                "bow_fibers requires positive span modes/fiber count and interior fiber station"
                    .into(),
            );
        }
        if let Some((limits, config)) = self.strain {
            limits.validate()?;
            config.validate()?;
            if self.axial_modes > fs_nlmodal::longitudinal::MAX_LONGITUDINAL_MODES
                || self.transverse_modes > fs_nlmodal::longitudinal::MAX_TRANSVERSE_MODES
            {
                return Err(
                    "fiber strain requires an explicitly supported axial/normal resolution".into(),
                );
            }
            let string_config = spec
                .geometric_stretching
                .as_ref()
                .map(|(_, config)| config)
                .or_else(|| spec.longitudinal_waves.as_ref().map(|(_, config)| config));
            if string_config.is_some_and(|existing| *existing != config) {
                return Err(
                    "fiber and string strain must use the same joint solve configuration".into(),
                );
            }
            if spec
                .forces
                .as_ref()
                .is_some_and(|program| program.initialization == ForceInitialization::StaticPreload)
            {
                return Err("fiber strain requires retained initial state; nonlinear static preload is not implemented".into());
            }
        }
        // The current geometry freezes the contact along a uniformly translating
        // fiber. Do not silently model an accelerating/reversing bow as one.
        for track in spec.schedule.tracks() {
            let GestureValue::Bow {
                velocity_m_per_s: initial,
                ..
            } = track.initial
            else {
                return Err("fiber ribbon needs a complete bow drive".into());
            };
            for event in &track.events {
                if !matches!(event.value, GestureValue::Bow { velocity_m_per_s, .. }
                    if velocity_m_per_s.to_bits() == initial.to_bits())
                {
                    return Err("frozen-station fiber ribbon requires constant support speed; load may vary".into());
                }
            }
        }
        self.axial_modes
            .checked_add(self.transverse_modes)
            .and_then(|n| n.checked_mul(cells))
            .ok_or_else(|| "fiber mode count overflow".into())
    }

    pub fn compile(
        &self,
        spec: &input::Spec,
        catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<NormalizedPack>,
        gate: &CancelGate,
    ) -> Result<TensionedFiberRibbon> {
        self.mode_count(spec)?;
        poll(gate)?;
        let requirements = [
            Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?,
            Requirement::try_new("young_modulus", Pressure::DIMS, Admit::StrictlyPositive)?,
            Requirement::try_new(ETA, DynViscosity::DIMS, Admit::NonNegative)?,
            Requirement::try_new(
                YOUNG_VISCOSITY_PROPERTY,
                DynViscosity::DIMS,
                Admit::NonNegative,
            )?,
        ];
        let material = resolve_properties(&self.pack, catalog, loader, spec, &requirements)?;
        let section = ConcentricElasticSection::try_new(
            0.0,
            vec![ConcentricMaterialLayer {
                outer_radius_m: self.radius_m,
                material,
            }],
        )?;
        let (width_m, cells) = spec
            .bow_ribbon
            .ok_or("bow_fibers requires a physical ribbon")?;
        let ribbon = TensionedFiberRibbon::try_new(
            FiberRibbonGeometry {
                length_m: self.length_m,
                width_m,
                fiber_count: self.fiber_count,
                tension_n: self.tension_n,
            },
            FiberRibbonResolution {
                cells,
                axial_modes: self.axial_modes,
                transverse_modes: self.transverse_modes,
            },
            section
                .clone()
                .with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY)?,
            section.with_kelvin_voigt_bending(ETA)?,
            spec.rate,
            modal_budget(),
            gate,
        )?;
        match self.strain {
            Some((limits, config)) => Ok(ribbon.with_axial_transverse_strain(limits, config)?),
            None => Ok(ribbon),
        }
    }
}
