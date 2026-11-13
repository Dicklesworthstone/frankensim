//! Full elastic stress at a fixed material point of a native Timoshenko frame.
//!
//! Accepted axial/bending strains and the canonical native transverse resultant
//! feed the original material owners. Circular Saint-Venant recovery supplies
//! both local shear components, including away from the bending plane; an
//! authored shear correction never becomes a point-stress factor. Local stress
//! recovery does not replace or add to native reduced mechanical storage.
//! This quasistatic cross-section model excludes end and cross-section inertia
//! corrections, plasticity, material viscosity and changing source conditions.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyCycle, FindleyCycleResult, FindleyError, FindleyLimits,
    StressClosureTolerance,
};
use fs_material::section::transverse_stress::CircularTransverseShearRecovery;
use fs_material::state_point::ResolvedMaterialStatePoint;

use super::elastic::{ElasticPoint, admit_scope};
use super::{Frame2d, FrameFatigueError, FrameFatiguePath, SolidError, poll};

/// Shear-sensitive fatigue-limit usage of a declared closed stress window.
/// The homogeneous circular/annular material must explicitly retain its
/// original isotropic parameter roles; unsupported sections refuse at admission.
#[derive(Debug, Clone)]
pub struct FrameFindleyGauge {
    pub element: usize,
    pub section: usize,
    pub layer: usize,
    /// Fixed material `[y,z]` in metres; +y is the native bending arm. With
    /// corotational geometry these axes follow the member, not the laboratory.
    pub cross_section_position_m: [f64; 2],
    pub temperature_axis: String,
    pub calibration: Arc<FindleyCalibration>,
    pub normals_in_material_frame: Vec<[f64; 3]>,
    pub limits: FindleyLimits,
    pub closure: StressClosureTolerance,
    /// Whole selected annulus engineering-shear bound, finite in (0, 0.05].
    pub maximum_engineering_shear: f64,
}

/// Actual accepted tensor at the gauge's permanent material point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameFindleySample {
    pub axial_mechanical_strain: f64,
    /// Tensor-Voigt `[xx,yy,zz,xy,yz,zx]` [Pa] in the co-rotating material frame.
    pub stress_pa: [f64; 6],
    /// Canonical native +y transverse resultant [N].
    pub transverse_force_n: f64,
    /// Maximum physical transverse engineering shear over the selected annulus.
    pub maximum_layer_engineering_shear: f64,
}

#[derive(Debug, Clone)]
pub(super) struct FindleyGaugeState {
    binding: FrameFindleyGauge,
    point: ElasticPoint,
    recovery: Arc<CircularTransverseShearRecovery>,
    cycle: FindleyCycle,
    sample: FrameFindleySample,
    index: usize,
}

pub(super) fn prepare(
    frame: &Frame2d,
    gauges: Vec<FrameFindleyGauge>,
    source_offset: usize,
    maximum_strain: f64,
    cx: &Cx<'_>,
) -> Result<Vec<FindleyGaugeState>, FrameFatigueError> {
    let mut prepared = Vec::with_capacity(gauges.len());
    let mut recoveries = BTreeMap::new();
    for (index, binding) in gauges.into_iter().enumerate() {
        poll(cx)?;
        if !binding.maximum_engineering_shear.is_finite()
            || binding.maximum_engineering_shear <= 0.0
            || binding.maximum_engineering_shear > 0.05
        {
            return Err(FrameFatigueError::Invalid(
                "declare a positive whole-annulus shear bound no greater than 0.05",
            ));
        }
        let point = ElasticPoint::new(
            frame,
            [binding.element, binding.section, binding.layer],
            binding.cross_section_position_m,
            maximum_strain,
        )?;
        let material = point.material.as_ref();
        let source_index = source_offset + index;
        admit_scope(
            material,
            &binding.temperature_axis,
            binding.calibration.scope(),
            &material.material().chemistry,
            source_index,
        )?;
        if let Some(source) = binding.calibration.source() {
            if binding.temperature_axis != source.request().temperature_axis
                || !source.matches_material_state(
                    material.card_identity(),
                    material.material(),
                    &material.query(),
                )
            {
                return Err(FrameFatigueError::SourceMismatch {
                    gauge: source_index,
                    what: "Findley source card, state or typed context differs from the original layer",
                });
            }
        }
        // One source projection per observed member, shared by its fixed gauges.
        let recovery = if let Some(recovery) = recoveries.get(&binding.element) {
            Arc::clone(recovery)
        } else {
            let source =
                frame
                    .material_shear_section(binding.element)
                    .ok_or(FrameFatigueError::Invalid(
                        "tensor fatigue requires an original material Timoshenko member",
                    ))?;
            let recovery = Arc::new(
                CircularTransverseShearRecovery::try_new(Arc::new(source.clone())).map_err(
                    |error| SolidError::MaterialRefused {
                        what: error.to_string(),
                    },
                )?,
            );
            recoveries.insert(binding.element, Arc::clone(&recovery));
            recovery
        };
        let cycle = FindleyCycle::new(
            binding.calibration.clone(),
            binding.normals_in_material_frame.clone(),
            binding.limits,
            binding.closure,
        )
        .map_err(|error| findley_error(index, error))?;
        let mut state = FindleyGaugeState {
            binding,
            point,
            recovery,
            cycle,
            index,
            sample: FrameFindleySample {
                axial_mechanical_strain: 0.0,
                stress_pa: [0.0; 6],
                transverse_force_n: 0.0,
                maximum_layer_engineering_shear: 0.0,
            },
        };
        state.observe(frame, cx)?;
        prepared.push(state);
    }
    Ok(prepared)
}

impl FindleyGaugeState {
    pub(super) fn observe(
        &mut self,
        frame: &Frame2d,
        cx: &Cx<'_>,
    ) -> Result<(), FrameFatigueError> {
        poll(cx)?;
        let normal = self.point.sample(frame, true)?;
        let force = frame.element_shear_response(self.binding.element)?.force_n;
        let shear = self
            .recovery
            .evaluate(
                self.binding.layer,
                self.binding.cross_section_position_m,
                [force, 0.0],
                self.binding.maximum_engineering_shear,
            )
            .map_err(|error| SolidError::MaterialRefused {
                what: error.to_string(),
            })?;
        let mut stress_pa = shear.stress_tensor_pa();
        stress_pa[0] = normal.normal_stress_pa;
        self.cycle
            .push_with_cancel(stress_pa, &mut || cx.checkpoint().is_err())
            .map_err(|error| findley_error(self.index, error))?;
        self.sample = FrameFindleySample {
            axial_mechanical_strain: normal.mechanical_strain,
            stress_pa,
            transverse_force_n: force,
            maximum_layer_engineering_shear: shear.maximum_layer_engineering_shear,
        };
        Ok(())
    }

    pub(super) fn admit_steps(&self, additional: u64) -> Result<(), FrameFatigueError> {
        if self.cycle.result().is_some() {
            return Err(findley_error(self.index, FindleyError::Finished));
        }
        let count = self
            .cycle
            .sample_count()
            .checked_add(additional)
            .ok_or(FrameFatigueError::Budget)?;
        if count > self.cycle.limits().maximum_samples
            || count
                .checked_mul(self.cycle.plane_count() as u64)
                .is_none_or(|n| n > self.cycle.limits().maximum_projection_work)
        {
            return Err(findley_error(self.index, FindleyError::Budget));
        }
        Ok(())
    }

    pub(super) fn checkpoint(&self) -> (u64, FrameFindleySample) {
        (self.cycle.sample_count(), self.sample)
    }
    pub(super) fn restore(&mut self, old: (u64, FrameFindleySample)) {
        self.cycle
            .discard_after(old.0)
            .expect("native observation only appends Findley samples");
        self.sample = old.1;
    }
}

pub(super) fn finish_bank(
    gauges: &mut [FindleyGaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<FindleyCycleResult>, FrameFatigueError> {
    let open: Vec<_> = gauges.iter().map(|g| g.cycle.result().is_none()).collect();
    let outcome = (|| {
        let mut results = Vec::with_capacity(gauges.len());
        for gauge in gauges.iter_mut() {
            results.push(
                gauge
                    .cycle
                    .finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| findley_error(gauge.index, error))?,
            );
        }
        poll(cx)?;
        Ok(results)
    })();
    if outcome.is_err() {
        for (gauge, was_open) in gauges.iter_mut().zip(open) {
            if was_open {
                gauge
                    .cycle
                    .discard_after(gauge.cycle.sample_count())
                    .expect("finalization retains all Findley samples");
            }
        }
    }
    outcome
}

impl FrameFatiguePath {
    pub fn findley_samples(&self) -> Vec<FrameFindleySample> {
        self.findley.iter().map(|g| g.sample).collect()
    }
    pub fn findley_results(&self) -> Vec<Option<FindleyCycleResult>> {
        self.findley
            .iter()
            .map(|g| g.cycle.result().copied())
            .collect()
    }
    pub fn findley_calibration(&self, gauge: usize) -> Option<&FindleyCalibration> {
        self.findley.get(gauge).map(|g| g.cycle.calibration())
    }
    pub fn findley_material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.findley.get(gauge).map(|g| g.point.material.as_ref())
    }

    /// Close all representative stress windows atomically; scalar histories
    /// remain open. Results are fatigue-limit usage, not Miner damage or life.
    pub fn close_findley_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FindleyCycleResult>, FrameFatigueError> {
        self.require_findley(cx)?;
        finish_bank(&mut self.findley, cx)
    }

    /// Start fresh windows at the latest accepted tensor. May explicitly discard
    /// a warm-up window; scalar histories and the global sample budget continue.
    pub fn restart_findley_cycle(&mut self, cx: &Cx<'_>) -> Result<(), FrameFatigueError> {
        self.require_findley(cx)?;
        let mut windows = Vec::with_capacity(self.findley.len());
        for gauge in &self.findley {
            poll(cx)?;
            let mut fresh = gauge
                .cycle
                .fresh()
                .map_err(|error| findley_error(gauge.index, error))?;
            fresh
                .push_with_cancel(gauge.sample.stress_pa, &mut || cx.checkpoint().is_err())
                .map_err(|error| findley_error(gauge.index, error))?;
            windows.push(fresh);
        }
        poll(cx)?;
        for (gauge, fresh) in self.findley.iter_mut().zip(windows) {
            gauge.cycle = fresh;
        }
        Ok(())
    }

    fn require_findley(&self, cx: &Cx<'_>) -> Result<(), FrameFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(FrameFatigueError::Finished);
        }
        if self.findley.is_empty() {
            return Err(FrameFatigueError::Invalid(
                "no Findley gauges were declared",
            ));
        }
        Ok(())
    }
}

fn findley_error(gauge: usize, error: FindleyError) -> FrameFatigueError {
    match error {
        FindleyError::Cancelled => FrameFatigueError::Cancelled,
        error => FrameFatigueError::Findley { gauge, error },
    }
}
