//! Original continuum permeability in the native liquid pressure/heat solver.
//!
//! A rigid, fully liquid-filled skeleton supplies only intrinsic geometry.
//! The existing resistance owner supplies the SAME current liquid rho/mu, Re,
//! pressure loss and power used by the thermal operating point. Gas slip and
//! gas relative saturation are never converted into a liquid constitutive law.
//! Absolute pressures admit the skeleton source; the liquid retains the native
//! low-compressibility, fixed-pressure transport approximation.

use fs_exec::Cx;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::porous::permeability::{ContinuumDarcyGeometry, DarcyGasBarrier};
use fs_qty::{Area, Length, MassFlowRate, Pressure, Temperature};

use super::super::LiquidHeatRegion;
use super::{
    LiquidFeedbackError, LiquidHeatPath, LiquidHydraulicRegion, LiquidResistanceSection,
    MaterialFlowResistance, PressureDrivenError, ResistanceError, ResistanceEvaluation, checkpoint,
};

mod sealed {
    pub trait Sealed {}
}

/// Native resistance owners accepted by the shared pressure/enthalpy paths.
/// Sealed to the existing geometry law and source-bound continuum owner.
pub trait LiquidResistanceModel: sealed::Sealed + Copy + core::fmt::Debug {
    /// Geometry-only native law; fluid properties are evaluated on every trial.
    fn geometry(&self) -> MaterialFlowResistance;
    /// Original porous source, when this section represents a porous skeleton.
    fn porous(&self) -> Option<PorousLiquidResistance<'_>> {
        None
    }
    /// Area shared with the thermal section.
    fn flow_area(&self) -> Area {
        self.geometry().flow_area()
    }
    /// Evaluate the native law and admit the actual bulk temperature/speed.
    /// Absolute source pressures are checked after the full trial is assembled.
    fn evaluate_material(
        &self,
        material: &ResolvedLiquidState,
        temperature_k: f64,
        mass: MassFlowRate,
    ) -> Result<ResistanceEvaluation, ResistanceError> {
        let loss = self.geometry().evaluate(material, mass)?;
        if let Some(porous) = self.porous() {
            porous.source.admit_temperature_speed(
                temperature_k,
                loss.volume_flow_m3_s / loss.law.flow_area().value(),
            )?;
        }
        Ok(loss)
    }
}
impl sealed::Sealed for MaterialFlowResistance {}
impl LiquidResistanceModel for MaterialFlowResistance {
    fn geometry(&self) -> MaterialFlowResistance {
        *self
    }
}

/// Original layered/tensor permeability with an explicit absolute outlet datum.
/// The datum is the downstream pressure of the WHOLE series path, not this
/// section. All porous sections of a path/network must use the same datum.
#[derive(Clone, Copy, Debug)]
pub struct PorousLiquidResistance<'a> {
    source: ContinuumDarcyGeometry<'a>,
    law: MaterialFlowResistance,
    outlet_pressure_pa: f64,
}
impl<'a> PorousLiquidResistance<'a> {
    /// Bind intrinsic permeability without sampling a fluid-dependent resistance.
    /// Reynolds length/band remain explicit geometry/model applicability data.
    ///
    /// # Errors
    /// Invalid geometry, pressure or Reynolds band, unrepresentable geometry, or
    /// a sealed barrier. A required positive-flow series path cannot cross a
    /// sealed layer; no permeability floor or infinite resistance is substituted.
    pub fn try_new(
        source: ContinuumDarcyGeometry<'a>,
        flow_area: Area,
        reynolds_length: Length,
        reynolds_range: (f64, f64),
        path_outlet_pressure: Pressure,
    ) -> Result<Self, ResistanceError> {
        let outlet_pressure_pa = path_outlet_pressure.value();
        if !outlet_pressure_pa.is_finite() || outlet_pressure_pa <= 0.0 {
            return Err(ResistanceError::Invalid(
                "porous liquid path requires a positive absolute outlet pressure",
            ));
        }
        let geometry = source
            .geometric_resistance_per_m()
            .ok_or(ResistanceError::Invalid(
                "sealed porous layer cannot carry positive series liquid flow",
            ))?;
        let law = MaterialFlowResistance::try_new(
            geometry / flow_area.value(),
            0.0,
            flow_area,
            reynolds_length,
            reynolds_range,
        )?;
        Ok(Self {
            source,
            law,
            outlet_pressure_pa,
        })
    }
    /// Borrow the unchanged original layer sources, frame and support.
    pub const fn source(&self) -> ContinuumDarcyGeometry<'a> {
        self.source
    }
    /// Absolute pressure at the common downstream end of the series path.
    pub const fn path_outlet_pressure_pa(&self) -> f64 {
        self.outlet_pressure_pa
    }
}
impl sealed::Sealed for PorousLiquidResistance<'_> {}
impl LiquidResistanceModel for PorousLiquidResistance<'_> {
    fn geometry(&self) -> MaterialFlowResistance {
        self.law
    }
    fn porous(&self) -> Option<PorousLiquidResistance<'_>> {
        Some(*self)
    }
}

/// Mixed ordinary channels and original porous materials in one native path.
#[derive(Clone, Copy, Debug)]
pub enum LiquidResistanceLaw<'a> {
    /// Existing authored/calibrated geometry.
    Geometry(MaterialFlowResistance),
    /// Original intrinsic permeability and applicability.
    Porous(PorousLiquidResistance<'a>),
}
impl sealed::Sealed for LiquidResistanceLaw<'_> {}
impl LiquidResistanceModel for LiquidResistanceLaw<'_> {
    fn geometry(&self) -> MaterialFlowResistance {
        match self {
            Self::Geometry(law) => *law,
            Self::Porous(law) => law.geometry(),
        }
    }
    fn porous(&self) -> Option<PorousLiquidResistance<'_>> {
        match self {
            Self::Geometry(_) => None,
            Self::Porous(law) => Some(*law),
        }
    }
}

/// Original skeleton and actual pressure faces of one accepted regional trial.
/// The hydraulic report alongside this receipt retains the liquid snapshot,
/// exact geometry, volume flow, Reynolds number and dissipated power.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousLiquidReceipt {
    /// Unchanged original layers, tensor frame, source claims and support.
    pub barrier: DarcyGasBarrier,
    /// Actual regional upstream absolute pressure, Pa.
    pub inlet_pressure_pa: f64,
    /// Actual regional downstream absolute pressure, Pa.
    pub outlet_pressure_pa: f64,
}

pub(super) fn path_outlet_pressure(
    sections: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
) -> Result<Option<f64>, PressureDrivenError> {
    let mut outlet = None;
    for section in sections {
        if let Some(source) = section.resistance.porous() {
            let value = source.path_outlet_pressure_pa();
            if outlet.is_some_and(|previous| previous != value) {
                return Err(PressureDrivenError::Invalid(
                    "porous sections require one absolute path outlet datum",
                ));
            }
            outlet = Some(value);
        }
    }
    Ok(outlet)
}

impl LiquidHeatPath<'_, '_> {
    // Called for EVERY assembled flow trial, including bracket endpoints and
    // solid-feedback candidates. Never certify only the final root.
    pub(super) fn admit_porous_trajectory(
        &self,
        cx: &Cx<'_>,
        sections: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        thermal: &[LiquidHeatRegion],
        hydraulic: &mut [LiquidHydraulicRegion],
    ) -> Result<(), PressureDrivenError> {
        let Some(mut downstream) = path_outlet_pressure(sections)? else {
            return Ok(());
        };
        if sections.len() != thermal.len() || sections.len() != hydraulic.len() {
            return Err(PressureDrivenError::Invalid(
                "incomplete porous hydraulic/thermal trajectory",
            ));
        }
        for ((section, heat), loss) in sections.iter().zip(thermal).zip(hydraulic).rev() {
            checkpoint(cx)?;
            let upstream = downstream + loss.loss.pressure_drop_pa;
            if !upstream.is_finite() || upstream <= downstream || downstream <= 0.0 {
                return Err(PressureDrivenError::Invalid(
                    "absolute porous pressure span is unrepresentable",
                ));
            }
            if let Some(porous) = section.resistance.porous() {
                let admit = |temperature_k, material: &ResolvedLiquidState| {
                    let speed = loss.loss.mass_flow_kg_per_s
                        / material.density_kg_m3()
                        / porous.law.flow_area().value();
                    porous
                        .source
                        .admit_flow(temperature_k, upstream, downstream, speed)
                        .map_err(|source| PressureDrivenError::Resistance {
                            region: section.region.to_owned(),
                            source: source.into(),
                        })
                };
                let inlet = heat.inlet.temperature_k();
                let outlet = heat.outlet.temperature_k();
                // Admit the actual endpoints as well as the regional bulk
                // snapshot. The native model does not resolve intra-region
                // pressure/temperature profiles or skeleton heat storage.
                admit(0.5 * inlet + 0.5 * outlet, &loss.material)?;
                for temperature in [inlet, outlet] {
                    let material = self
                        .law
                        .at_temperature(Temperature::new(temperature))
                        .map_err(LiquidFeedbackError::from)?;
                    admit(temperature, &material)?;
                }
                loss.porous = Some(PorousLiquidReceipt {
                    barrier: porous.source.barrier().clone(),
                    inlet_pressure_pa: upstream,
                    outlet_pressure_pa: downstream,
                });
            }
            downstream = upstream;
        }
        Ok(())
    }
}
