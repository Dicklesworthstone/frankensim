//! Finite atmospheric water coupled to the original swelling-metric diffusion.
//!
//! Each endpoint RH trial begins at the SAME accepted field. Inside that trial,
//! the original native solid resolves the source swelling, and the original
//! Kirchhoff transport uses its actual axial lengths and surface/contact area.
//! Only the original atmosphere owner debits vapor and reconstructs the gas.
//! No geometric, sorption, diffusion or gas model is fitted in this composition.

use super::*;
use crate::moisture::atmosphere::AtmosphereMoistureCandidate;
use crate::moisture::layered::atmosphere::{admit_exchange, surface_sorption, water_change};

/// Explicit physical admission and cumulative work for one proposed interval.
#[derive(Clone, Copy, Debug)]
pub struct MoistureSpatialAtmosphereControls {
    /// Original finite-water RH bracket, mass tolerances and root iteration cap.
    pub atmosphere: AtmosphereMoistureSolve,
    /// Native geometry's permitted departure from the separable axial map [m].
    /// Bending, shear, taper or transverse warping beyond this value refuse.
    pub maximum_geometry_residual_m: f64,
    /// Counts EVERY atmospheric trial, including initial/bracket and rejected
    /// trials. Each trial enters the original nonlinear diffusion/solid owner.
    pub maximum_transport_evaluations: usize,
    /// Counts EVERY native solid solve across ALL humidity and nonlinear
    /// geometry trials. The cap is checked before each solve, never reset when
    /// the atmosphere root discards a trial.
    pub maximum_native_solid_solves: usize,
}

/// Accepted spatial physics and actual cumulative work over the whole proposal.
#[derive(Clone, Debug)]
pub struct MoistureSpatialAtmosphereReport {
    /// Original accepted transport, axial metric, actual area and solid identity.
    /// Its native_solid_solves counts only the accepted humidity trial.
    pub spatial: MoistureSpatialDiffusionReport,
    /// Complete nonlinear transport evaluations, including rejected RH trials.
    pub transport_evaluations: usize,
    /// Complete native solid solves, including rejected RH/geometric trials.
    pub native_solid_solves: usize,
}

/// Original finite gas and layered water accepted with swelling-dependent flux.
pub type MoistureSpatialAtmosphereCandidate =
    AtmosphereMoistureCandidate<MoistureLayeredSlab, MoistureSpatialAtmosphereReport>;

impl MoistureLayeredSolid {
    /// Couple finite atmospheric depletion/enrichment to actual swelling metrics.
    ///
    /// Every RH trial uses the original accepted field and source-bound solid.
    /// Endpoint-dependent water transport and mechanics converge together before
    /// their water transfer enters the atmosphere balance. Finite surface and
    /// interface conductances use the candidate's current area; bulk transport
    /// uses its inverse axial stretch squared and the conserved dry skeleton.
    /// No equilibrium moisture is assigned throughout the specimen.
    ///
    /// The physical scope is the existing isothermal separable-prism reduction.
    /// Spatial isotropic diffusivity, reference moisture, principal swelling
    /// axes and their source domains remain explicitly authored. This API does
    /// not infer stress-assisted flux, hysteresis or a new source law.
    ///
    /// The gas FREE VOLUME stays fixed at the original atmosphere's volume.
    /// Specimen expansion is not subtracted from it. This describes a fixed
    /// reservoir-volume approximation; it does not model displacement of gas in
    /// a tightly enclosed chamber. Gas pressure is not a mechanical traction,
    /// and gas internal-energy change is not heat of sorption or boundary work.
    ///
    /// All inputs are borrowed. Water/carrier/total mass, source, geometry,
    /// mechanics, work-limit and cancellation checks precede the returned
    /// candidate; any refusal preserves all accepted physical owners.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_spatial_diffusion_finite_atmosphere(
        &self,
        atmosphere: &FiniteMoistureAtmosphere,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: MoistureSpatialAtmosphereControls,
        cx: &Cx<'_>,
    ) -> Result<(MoistureSpatialAtmosphereCandidate, MoistureSolidState), AtmosphereMoistureError>
    {
        poll(cx)?;
        if !positive(controls.maximum_geometry_residual_m)
            || controls.maximum_transport_evaluations == 0
            || controls.maximum_native_solid_solves == 0
        {
            return Err(AtmosphereMoistureError::Invalid(
                "positive spatial geometry tolerance and cumulative work budgets required",
            ));
        }
        let exposed = admit_exchange(
            moisture,
            atmosphere.gas().state().state().temperature,
            exposed,
            controls.atmosphere,
            cx,
        )?;
        let mut evaluations = 0usize;
        let mut work = SpatialSolveBudget {
            used: 0,
            maximum: controls.maximum_native_solid_solves,
        };
        let evaluate = |rh| -> Result<_, AtmosphereMoistureError> {
            poll(cx)?;
            if evaluations >= controls.maximum_transport_evaluations {
                return Err(AtmosphereMoistureError::Invalid(
                    "cumulative spatial atmosphere transport-evaluation budget exhausted",
                ));
            }
            evaluations += 1;
            let faces = exposed.map(|open| {
                if open {
                    MoistureLayeredFace::PrescribedRelativeHumidity(rh)
                } else {
                    MoistureLayeredFace::Sealed
                }
            });
            let candidate = self.advance_spatial_diffusion_budgeted(
                moisture,
                dt_s,
                &faces,
                transport,
                nonlinear,
                controls.maximum_geometry_residual_m,
                &mut work,
                cx,
            )?;
            let change = water_change(moisture, &candidate.moisture);
            let end_time = candidate.moisture.time_s();
            Ok((candidate, change, end_time))
        };
        let (atmosphere, candidate, exchange) = atmosphere.advance_exchange(
            moisture.time_s(),
            exposed == [false; 2],
            controls.atmosphere,
            cx,
            evaluate,
        )?;
        let sorption = surface_sorption(moisture, &candidate.transport.transport, cx)?;
        let MoistureLayeredSolidStep {
            moisture,
            transport: spatial,
            solid,
        } = candidate;
        let report = MoistureSpatialAtmosphereReport {
            spatial,
            transport_evaluations: evaluations,
            native_solid_solves: work.used,
        };
        let step = AtmosphereMoistureStep {
            transport: report,
            sorption,
            relative_humidity: exchange.relative_humidity,
            vapor_mass_kg: exchange.vapor_mass_kg,
            water_to_solid_kg: exchange.water_to_solids_kg,
            water_balance_residual_kg: exchange.water_balance_residual_kg,
            maximum_carrier_mass_residual_kg: exchange.maximum_carrier_mass_residual_kg,
            gas_internal_energy_change_j: exchange.gas_internal_energy_change_j,
            iterations: exchange.iterations,
        };
        poll(cx)?;
        Ok((
            AtmosphereMoistureCandidate::from_accepted(atmosphere, moisture, step),
            solid,
        ))
    }
}
