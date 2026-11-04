//! Gas-side convection at the original zero-storage liquid surface.
//!
//! A higher owner supplies its admitted correlation at EACH trial surface
//! temperature. This module owns no correlations, gas momentum or additional
//! heat integrator. The original liquid half-cell, radiation, source pinning,
//! finite-atmosphere balance and thermal step admission remain in force.

use super::*;
use fs_evidence::Evidence;

/// Gas-side heat-transfer coefficient [W/(m² K)], independent of a correlation crate.
pub type FreeSurfaceHeatTransferCoefficient = fs_qty::Qty<0, 1, -3, -1, 0, 0>;

/// Original gas states and actual geometry for one scalar surface-solve trial.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceConvectionQuery<'a> {
    /// Index into the original call's surface contacts.
    pub patch_index: usize,
    pub surface: FreeSurfaceInterfaceFace,
    /// Cartesian unit normal from liquid into gas. Reconstructed curved and
    /// height-graph surfaces are not admitted by this bounded interface.
    pub outward_normal: [f64; 3],
    pub surface_temperature_k: f64,
    /// Current prescribed reservoir or original finite-atmosphere parcel.
    pub bulk: &'a ResolvedThermochemicalGasState,
    /// The SAME gas source claims at (surface T + bulk gas T)/2 and current p.
    /// Use these held properties for the gas correlation; do not resolve a
    /// second ambient gas or replace a failed original claim with defaults.
    pub film: &'a ResolvedThermochemicalGasState,
}

/// A call-local, repeatable constitutive evaluation for a physical plate patch.
///
/// The provider must validate its geometry, physical characteristic length and
/// correlation regime on every query, including equilibrium. Its coefficient
/// is evaluated at the actual trial surface temperature, not the liquid cell T.
/// Returning `None` is allowed ONLY at exact surface/bulk-gas equilibrium (for
/// a natural-convection relation with no defined zero-Rayleigh coefficient).
/// Positive finite coefficients retain their exact original evidence.
///
/// Calls include rejected nonlinear trials, thermal endpoints and transactional
/// clones. A provider must not advance physical state or integrate heat. A
/// correlation refusal aborts the transaction; it is never a stagnant-film
/// fallback. The declared relative velocity/geometry are held for this split.
pub trait FreeSurfaceGasConvection {
    fn coefficient(
        &self,
        query: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>;
}

/// Exact gas-side law applied to this surface sample.
#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceGasHeatTransfer {
    /// Original source gas conductivity divided by the authored film thickness.
    StagnantFilm,
    /// Exact source-backed correlation evidence. None denotes ONLY exact
    /// surface/bulk equilibrium admitted by the provider, with zero gas heat.
    Convection(Option<Evidence<FreeSurfaceHeatTransferCoefficient>>),
}

pub(super) fn cartesian_normal(face: FreeSurfaceGasFace) -> Result<[f64; 3], Error> {
    let mut normal = [0.0; 3];
    let mut changed = false;
    for (axis, component) in normal.iter_mut().enumerate() {
        if face.liquid[axis] == face.gas[axis] {
            continue;
        }
        if changed || face.liquid[axis].abs_diff(face.gas[axis]) != 1 {
            return Err(Error::Invalid(
                "surface convection requires adjacent Cartesian cells",
            ));
        }
        changed = true;
        *component = if face.gas[axis] > face.liquid[axis] {
            1.0
        } else {
            -1.0
        };
    }
    if !changed {
        return Err(Error::Invalid(
            "surface convection requires distinct liquid and gas cells",
        ));
    }
    Ok(normal)
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Solve liquid/gas/radiation exchange with an admitted gas convection law.
    ///
    /// All surface contacts use the same provider; `patch_index` can select
    /// distinct explicitly declared physical plates. Plate length is a physical
    /// boundary-layer length, NEVER the cell width or a new entrance per face.
    /// The gas-film distance remains admitted geometry but contributes no
    /// second resistance: the supplied h replaces k_g/d_g in the shared root.
    /// The liquid half-cell resistance and original radiation remain present.
    ///
    /// Only Cartesian planar faces are admitted. This retains the original
    /// held-geometry thermal/native-flow split and its refinement requirement;
    /// it does not model evaporation, resolved gas flow or arbitrary curved
    /// deforming surfaces. Implicit trials, when selected, use this same law.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_surface_convection(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        surfaces: &[FreeSurfaceInterfaceContact<'_>],
        solve: FreeSurfaceInterfaceSolve,
        convection: &dyn FreeSurfaceGasConvection,
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        let applied = if with_viscous_heating {
            self.viscous_sources(gate)?
        } else {
            Vec::new()
        };
        let conduction = self.advance_with_thermal_boundaries(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            interfaces,
            &[],
            &[],
            surfaces,
            solve,
            Some(convection),
            None,
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }

    /// The same convective interface with one finite, well-mixed gas parcel.
    ///
    /// Its current original gas state feeds every trial; the parcel stores the
    /// exact opposite gas heat on each thermal substep. Its actual capacity and
    /// the eliminated gas conductance bound explicit substeps together with
    /// the liquid and finite walls. All states publish together after acceptance.
    /// `maximum_density_excursion`, when present, also applies this parcel's
    /// post-thermal EOS pressure through the existing native pressure path.
    /// The caller declares a single connected fixed-volume atmosphere; no gas
    /// compression work or gas momentum is inferred from liquid motion.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_surface_convection(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        gas: &mut IsochoricGasParcel,
        surfaces: &[FreeSurfaceInterfaceFace],
        gas_controls: GasEnergySolve,
        solve: FreeSurfaceInterfaceSolve,
        convection: &dyn FreeSurfaceGasConvection,
        maximum_density_excursion: Option<f64>,
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        self.step_with_finite_gas_surface_exchange_inner(
            faces,
            heat,
            controls,
            bodies,
            contacts,
            interfaces,
            gas,
            surfaces,
            gas_controls,
            solve,
            Some(convection),
            maximum_density_excursion,
            with_viscous_heating,
            gate,
        )
    }
}
