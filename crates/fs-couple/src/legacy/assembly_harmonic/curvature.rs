//! Analytic source-bound radiation-load curvature on the original shell assembly.

use super::{AssemblyMimoHarmonicError, PreparedAssemblyMimoHarmonicResponse, poll};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicFrequencyCurvatureControls;
use fs_plate::shell::preload::assembly::vibration::harmonic::ports::{
    AssemblyPortCurvatureBudget, AssemblyPortHarmonicCurvatureResponse,
};
use fs_vfit::mimo::load::PreparedMimoLoad;

impl PreparedAssemblyMimoHarmonicResponse<'_, '_, '_, '_, '_> {
    /// Evaluate motion, reactions and physical work curvature with the original
    /// admitted continuous rational load and the existing coupled factor.
    ///
    /// Both frequency derivatives come from the source's exact real/conjugate
    /// poles and all mutual entries. The structural exp(-i*w*t) chain rule is
    /// Z_w=-i Z_s and Z_ww=-Z_ss. Optional original component force derivatives
    /// retain their phase and all product-rule terms. The independent physical
    /// port recovery supplies actual signed-power extrema and linewidth queries.
    ///
    /// No fit, interpolation, pole matching, runtime-history reset, or promotion
    /// of a sampled passivity report occurs. Reusing this result for a finite
    /// parameter change or outside the original source band is not authorized by
    /// a local frequency derivative.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_frequency_curvature(
        &self,
        component_loads: &[&[C64]],
        component_load_derivatives: Option<&[&[C64]]>,
        component_load_second_derivatives: Option<&[&[C64]]>,
        controls: HarmonicFrequencyCurvatureControls,
        physical_budget: AssemblyPortCurvatureBudget,
        cx: &Cx<'_>,
    ) -> Result<AssemblyPortHarmonicCurvatureResponse, AssemblyMimoHarmonicError> {
        let (first, second, controls) =
            source_frequency_curvature(self.source, self.angular_frequency_rad_s, controls, cx)?;
        self.response
            .solve_with_frequency_curvature(
                component_loads,
                component_load_derivatives,
                component_load_second_derivatives,
                &first,
                &second,
                controls,
                physical_budget,
                cx,
            )
            .map_err(AssemblyMimoHarmonicError::Structure)
    }
}

/// Three source-owned complex matrices coexist during evaluation; the two
/// derivative matrices are converted in place. No extra matrix or sampled fit.
/// Work score 1024*p²*(expanded_order+1) reserves the scaled quotient chains.
fn source_frequency_curvature(
    source: &PreparedMimoLoad,
    omega: f64,
    controls: HarmonicFrequencyCurvatureControls,
    cx: &Cx<'_>,
) -> Result<(Vec<C64>, Vec<C64>, HarmonicFrequencyCurvatureControls), AssemblyMimoHarmonicError> {
    poll(cx)?;
    let entries = (source.ports() as u128) * (source.ports() as u128);
    let bytes = 48 * entries;
    let work = 1024 * entries * (source.source_model().order() as u128 + 1);
    if bytes > controls.maximum_additional_payload_bytes as u128
        || work > controls.maximum_work as u128
    {
        return Err(AssemblyMimoHarmonicError::Structure(
            "original MIMO curvature exceeds payload or work limits".into(),
        ));
    }
    source
        .continuous_operator()
        .admit_frequency(omega)
        .map_err(AssemblyMimoHarmonicError::Source)?;
    let mut jet = source
        .source_model()
        .eval_jet(C64::new(0.0, -omega))
        .map_err(|error| AssemblyMimoHarmonicError::Structure(error.to_string()))?;
    for value in &mut jet.laplace_derivative {
        *value = C64::new(value.im, -value.re);
    }
    for value in &mut jet.laplace_second_derivative {
        *value = -*value;
    }
    poll(cx)?;
    Ok((
        jet.laplace_derivative,
        jet.laplace_second_derivative,
        HarmonicFrequencyCurvatureControls {
            maximum_additional_payload_bytes: controls.maximum_additional_payload_bytes
                - bytes as usize,
            maximum_work: controls.maximum_work - work as u64,
        },
    ))
}

#[cfg(test)]
mod tests;
