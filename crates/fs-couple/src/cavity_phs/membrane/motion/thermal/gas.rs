//! Original caloric candidates pulled through native membrane geometry.
use super::*;

pub(super) struct GasStep<R> {
    pub endpoint: MembraneGas,
    pub work_pressure_pa: f64,
    pub sampling_temperature_k: f64,
    pub sampling_pressure_pa: f64,
    pub sampling_state: Option<MembraneGasSample>,
    pub report: R,
}

/// Only the original source receipt selects its heat/volume law. Mechanics,
/// geometry, Newton actions and publication remain shared and monomorphized.
pub(super) trait HeatVolumeLaw: Copy {
    fn heat_volume(
        gas: &MembraneGas,
        volume: f64,
        heat: f64,
        controls: GasEnergySolve,
        retain_sample: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasStep<Self>, GasMembraneMotionError>;
    fn stored_energy_change_j(self) -> f64;
    fn work_by_gas_j(self) -> f64;
    fn pressure_final_volume_derivative_pa_m3(self) -> f64;
}

impl HeatVolumeLaw for GasHeatVolumeStep {
    fn heat_volume(
        gas: &MembraneGas,
        volume: f64,
        heat: f64,
        controls: GasEnergySolve,
        retain_sample: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasStep<Self>, GasMembraneMotionError> {
        let (endpoint, temperature, pressure, sample, report) = match gas {
            MembraneGas::Ideal(gas) => {
                let candidate = gas.midpoint_heat_volume(volume, heat, controls, cancelled)?;
                let temperature = candidate.midpoint().state().temperature;
                let pressure = candidate.midpoint().state().pressure;
                let sample =
                    retain_sample.then(|| MembraneGasSample::Ideal(candidate.midpoint().clone()));
                let (endpoint, report) = candidate.into_parts();
                (
                    MembraneGas::Ideal(endpoint),
                    temperature,
                    pressure,
                    sample,
                    report,
                )
            }
            MembraneGas::PengRobinson(gas) => {
                let candidate = gas
                    .clone()
                    .into_isochoric()
                    .midpoint_heat_volume(volume, heat, controls, cancelled)?;
                let temperature = candidate.midpoint().temperature_k();
                let pressure = candidate.midpoint().pressure_pa();
                let sample = retain_sample
                    .then(|| MembraneGasSample::PengRobinson(candidate.midpoint().clone()));
                let (endpoint, report) = candidate.into_parts();
                (
                    MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(endpoint)),
                    temperature,
                    pressure,
                    sample,
                    report,
                )
            }
        };
        Ok(GasStep {
            endpoint,
            work_pressure_pa: pressure,
            sampling_temperature_k: temperature,
            sampling_pressure_pa: pressure,
            sampling_state: sample,
            report,
        })
    }
    fn stored_energy_change_j(self) -> f64 {
        self.stored_energy_change_j
    }
    fn work_by_gas_j(self) -> f64 {
        self.work_by_gas_j
    }
    fn pressure_final_volume_derivative_pa_m3(self) -> f64 {
        self.pressure_final_volume_derivative_pa_m3
    }
}

impl HeatVolumeLaw for GasDiscreteGradientStep {
    fn heat_volume(
        gas: &MembraneGas,
        volume: f64,
        heat: f64,
        controls: GasEnergySolve,
        retain_sample: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasStep<Self>, GasMembraneMotionError> {
        let (endpoint, temperature, pressure, sample, report) = match gas {
            MembraneGas::Ideal(gas) => {
                let candidate =
                    gas.discrete_gradient_heat_volume(volume, heat, controls, cancelled)?;
                let temperature = candidate.thermodynamic_mean().state().temperature;
                let pressure = candidate.thermodynamic_mean().state().pressure;
                let sample = retain_sample
                    .then(|| MembraneGasSample::Ideal(candidate.thermodynamic_mean().clone()));
                let (endpoint, report) = candidate.into_parts();
                (
                    MembraneGas::Ideal(endpoint),
                    temperature,
                    pressure,
                    sample,
                    report,
                )
            }
            MembraneGas::PengRobinson(gas) => {
                let candidate = gas
                    .clone()
                    .into_isochoric()
                    .discrete_gradient_heat_volume(volume, heat, controls, cancelled)?;
                let temperature = candidate.thermodynamic_mean().temperature_k();
                let pressure = candidate.thermodynamic_mean().pressure_pa();
                let sample = retain_sample.then(|| {
                    MembraneGasSample::PengRobinson(candidate.thermodynamic_mean().clone())
                });
                let (endpoint, report) = candidate.into_parts();
                (
                    MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(endpoint)),
                    temperature,
                    pressure,
                    sample,
                    report,
                )
            }
        };
        Ok(GasStep {
            endpoint,
            // The PR material sampling state's EOS pressure is not the
            // discrete Gibbs effort. Only the source receipt drives momentum.
            work_pressure_pa: report.work_pressure_pa,
            sampling_temperature_k: temperature,
            sampling_pressure_pa: pressure,
            sampling_state: sample,
            report,
        })
    }
    fn stored_energy_change_j(self) -> f64 {
        self.stored_energy_change_j
    }
    fn work_by_gas_j(self) -> f64 {
        self.work_by_gas_j
    }
    fn pressure_final_volume_derivative_pa_m3(self) -> f64 {
        self.pressure_final_volume_derivative_pa_m3
    }
}

pub(super) struct Geometry {
    pub displacement: Vec<[f64; 3]>,
    pub internal: Vec<[f64; 3]>,
    pub strain_energy_j: f64,
    pub volume: f64,
    pub exterior_volume: Option<f64>,
}

pub(super) fn geometry(
    model: &GasMembrane<'_>,
    q: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Geometry, GasMembraneMotionError> {
    if q.len() != model.free_dofs.len() || q.iter().any(|x| !x.is_finite()) {
        return Err(GasMembraneMotionError::Invalid(
            "finite dimension-matched coordinates required",
        ));
    }
    let displacement = model.expand(q, false);
    model.admit_shape(&displacement)?;
    let membrane = model.patch.evaluate(&displacement, &mut *cancelled)?;
    let sweep =
        model
            .patch
            .pressure_work_j(&model.reference_m, &displacement, 1.0, &mut *cancelled)?;
    let admitted = |reference: f64, change: f64, interval: [f64; 2]| {
        let volume = reference + change;
        if !volume.is_finite()
            || volume < interval[0]
            || volume > interval[1]
            || (change != 0.0 && volume == reference)
        {
            Err(GasMembraneMotionError::Invalid(
                "thermal membrane volume outside its interval or unrepresentable",
            ))
        } else {
            Ok(volume)
        }
    };
    let volume = admitted(model.gas.volume_m3(), sweep, model.volume_interval_m3)?;
    let exterior_volume = model
        .sealed_exterior()
        .map(|gas| admitted(gas.gas.volume_m3(), -sweep, gas.volume_interval_m3))
        .transpose()?;
    Ok(Geometry {
        displacement,
        internal: membrane.internal_gradient_n,
        strain_energy_j: membrane.strain_energy_j,
        volume,
        exterior_volume,
    })
}

pub(super) struct Trial<R = GasHeatVolumeStep> {
    pub stage: Geometry,
    pub endpoint: Geometry,
    /// Exact path mean of grad(V), work-conjugate to the complete displacement.
    pub mean_volume_gradient: Vec<[f64; 3]>,
    /// Endpoint gradient for the native gas derivative and physical readback.
    pub endpoint_volume_gradient: Vec<[f64; 3]>,
    pub gas: GasStep<R>,
    pub exterior: Option<GasStep<R>>,
    pub exterior_work_pressure_pa: f64,
    pub net_pressure_pa: f64,
}

impl<R: HeatVolumeLaw> Trial<R> {
    #[cfg(test)]
    pub(super) fn new(
        model: &GasMembrane<'_>,
        initial: &GasMembraneThermalEvaluation,
        stage_q: &[f64],
        endpoint_q: &[f64],
        heat: f64,
        exterior_heat: Option<f64>,
        controls: GasMembraneThermalControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasMembraneMotionError> {
        Self::new_sampled(
            model,
            initial,
            stage_q,
            endpoint_q,
            heat,
            exterior_heat,
            controls,
            false,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn new_sampled(
        model: &GasMembrane<'_>,
        initial: &GasMembraneThermalEvaluation,
        stage_q: &[f64],
        endpoint_q: &[f64],
        heat: f64,
        exterior_heat: Option<f64>,
        controls: GasMembraneThermalControls,
        retain_sample: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasMembraneMotionError> {
        let stage = geometry(model, stage_q, &mut *cancelled)?;
        let endpoint = geometry(model, endpoint_q, &mut *cancelled)?;
        let mean_volume_gradient = model.patch.pressure_segment_load_n(
            &initial.displacement_m,
            &endpoint.displacement,
            1.0,
            &mut *cancelled,
        )?;
        let endpoint_volume_gradient =
            model
                .patch
                .pressure_load_n(&endpoint.displacement, 1.0, &mut *cancelled)?;
        let gas = R::heat_volume(
            &initial.gas,
            endpoint.volume,
            heat,
            controls.gas,
            retain_sample,
            &mut *cancelled,
        )?;
        let exterior = if let Some(old) = &initial.exterior_gas {
            Some(R::heat_volume(
                old,
                endpoint.exterior_volume.expect("sealed volume"),
                exterior_heat.expect("admitted exterior heat"),
                controls.exterior_gas.expect("admitted exterior controls"),
                false,
                &mut *cancelled,
            )?)
        } else {
            None
        };
        let exterior_work_pressure_pa = exterior
            .as_ref()
            .map_or(initial.exterior_pressure_pa, |gas| gas.work_pressure_pa);
        let net_pressure_pa = gas.work_pressure_pa - exterior_work_pressure_pa;
        Ok(Self {
            stage,
            endpoint,
            mean_volume_gradient,
            endpoint_volume_gradient,
            gas,
            exterior,
            exterior_work_pressure_pa,
            net_pressure_pa,
        })
    }

    pub(super) fn tangent(
        &self,
        model: &GasMembrane<'_>,
        initial_displacement: &[[f64; 3]],
        direction: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, GasMembraneMotionError> {
        if direction.len() != model.free_dofs.len() || direction.iter().any(|x| !x.is_finite()) {
            return Err(GasMembraneMotionError::Invalid(
                "finite thermal tangent direction required",
            ));
        }
        let full = model.expand(direction, true);
        let internal =
            model
                .patch
                .tangent_apply(&self.stage.displacement, &full, &mut *cancelled)?;
        let geometric = model.patch.pressure_segment_tangent_apply(
            initial_displacement,
            &self.endpoint.displacement,
            &full,
            self.net_pressure_pa,
            &mut *cancelled,
        )?;
        let mut dv = Dd::from_f64(0.0);
        for (&dof, &d) in model.free_dofs.iter().zip(direction) {
            dv = dv + Dd::from_f64(self.endpoint_volume_gradient[dof / 3][dof % 3] * d);
        }
        let derivative = self.gas.report.pressure_final_volume_derivative_pa_m3()
            + self.exterior.as_ref().map_or(0.0, |gas| {
                gas.report.pressure_final_volume_derivative_pa_m3()
            });
        // The source derivative holds interval heat fixed. q1=2*q_stage-q_old,
        // so BOTH endpoint-volume and mean-geometry derivatives receive a factor
        // of two. The pressure force uses the exact path mean of grad(V).
        // This source/geometry cross-gradient action is generally NONSYMMETRIC.
        let dp = derivative * (2.0 * dv.to_f64());
        let action: Vec<_> = model
            .free_dofs
            .iter()
            .map(|&i| {
                internal[i / 3][i % 3]
                    - 2.0 * geometric[i / 3][i % 3]
                    - dp * self.mean_volume_gradient[i / 3][i % 3]
            })
            .collect();
        if action.iter().any(|x| !x.is_finite()) {
            return Err(GasMembraneMotionError::Invalid(
                "thermal pressure tangent overflow",
            ));
        }
        Ok(action)
    }

    pub(super) fn into_evaluation(
        self,
        loads: &[[f64; 3]],
    ) -> Result<GasMembraneThermalEvaluation, GasMembraneMotionError> {
        let pressure = self
            .gas
            .endpoint
            .pressure_response()
            .map_err(GasMembraneError::from)?;
        let exterior_pressure_pa = match &self.exterior {
            Some(gas) => {
                gas.endpoint
                    .pressure_response()
                    .map_err(GasMembraneError::from)?
                    .pressure_pa
            }
            None => self.exterior_work_pressure_pa,
        };
        let mut residual_n = self.endpoint.internal;
        for ((row, gradient), loads) in residual_n
            .iter_mut()
            .zip(&self.endpoint_volume_gradient)
            .zip(loads)
        {
            for axis in 0..3 {
                row[axis] -=
                    loads[axis] + (pressure.pressure_pa - exterior_pressure_pa) * gradient[axis];
            }
        }
        if !finite_vectors(&residual_n, self.endpoint.displacement.len()) {
            return Err(GasMembraneMotionError::Invalid(
                "thermal endpoint force overflow",
            ));
        }
        Ok(GasMembraneThermalEvaluation {
            displacement_m: self.endpoint.displacement,
            gas: self.gas.endpoint,
            exterior_gas: self.exterior.map(|gas| gas.endpoint),
            pressure,
            exterior_pressure_pa,
            residual_n,
            strain_energy_j: self.endpoint.strain_energy_j,
        })
    }
}
