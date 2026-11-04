//! Original caloric midpoint candidates pulled through native membrane geometry.
use super::*;

pub(super) struct GasStep {
    pub endpoint: MembraneGas,
    pub midpoint_pressure_pa: f64,
    pub report: GasHeatVolumeStep,
}

fn heat_volume(
    gas: &MembraneGas,
    volume: f64,
    heat: f64,
    controls: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasStep, GasMembraneMotionError> {
    let (endpoint, midpoint_pressure_pa, report) = match gas {
        MembraneGas::Ideal(gas) => {
            let candidate = gas.midpoint_heat_volume(volume, heat, controls, cancelled)?;
            let pressure = candidate.midpoint().state().pressure;
            let (endpoint, report) = candidate.into_parts();
            (MembraneGas::Ideal(endpoint), pressure, report)
        }
        MembraneGas::PengRobinson(gas) => {
            let candidate = gas
                .clone()
                .into_isochoric()
                .midpoint_heat_volume(volume, heat, controls, cancelled)?;
            let pressure = candidate.midpoint().pressure_pa();
            let (endpoint, report) = candidate.into_parts();
            (
                MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(endpoint)),
                pressure,
                report,
            )
        }
    };
    Ok(GasStep {
        endpoint,
        midpoint_pressure_pa,
        report,
    })
}

pub(super) struct Geometry {
    pub displacement: Vec<[f64; 3]>,
    pub internal: Vec<[f64; 3]>,
    pub gradient: Vec<[f64; 3]>,
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
    let gradient = model
        .patch
        .pressure_load_n(&displacement, 1.0, &mut *cancelled)?;
    Ok(Geometry {
        displacement,
        internal: membrane.internal_gradient_n,
        gradient,
        strain_energy_j: membrane.strain_energy_j,
        volume,
        exterior_volume,
    })
}

pub(super) struct Trial {
    pub stage: Geometry,
    pub endpoint: Geometry,
    pub gas: GasStep,
    pub exterior: Option<GasStep>,
    pub exterior_midpoint_pa: f64,
    pub net_pressure_pa: f64,
}

impl Trial {
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
        let stage = geometry(model, stage_q, &mut *cancelled)?;
        let endpoint = geometry(model, endpoint_q, &mut *cancelled)?;
        let gas = heat_volume(
            &initial.gas,
            endpoint.volume,
            heat,
            controls.gas,
            &mut *cancelled,
        )?;
        let exterior = if let Some(old) = &initial.exterior_gas {
            Some(heat_volume(
                old,
                endpoint.exterior_volume.expect("sealed volume"),
                exterior_heat.expect("admitted exterior heat"),
                controls.exterior_gas.expect("admitted exterior controls"),
                &mut *cancelled,
            )?)
        } else {
            None
        };
        let exterior_midpoint_pa = exterior
            .as_ref()
            .map_or(initial.exterior_pressure_pa, |gas| gas.midpoint_pressure_pa);
        let net_pressure_pa = gas.midpoint_pressure_pa - exterior_midpoint_pa;
        Ok(Self {
            stage,
            endpoint,
            gas,
            exterior,
            exterior_midpoint_pa,
            net_pressure_pa,
        })
    }

    pub(super) fn tangent(
        &self,
        model: &GasMembrane<'_>,
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
        let geometric = model.patch.pressure_tangent_apply(
            &self.stage.displacement,
            &full,
            self.net_pressure_pa,
            &mut *cancelled,
        )?;
        let mut dv = Dd::from_f64(0.0);
        for (&dof, &d) in model.free_dofs.iter().zip(direction) {
            dv = dv + Dd::from_f64(self.endpoint.gradient[dof / 3][dof % 3] * d);
        }
        let derivative = self.gas.report.pressure_final_volume_derivative_pa_m3
            + self
                .exterior
                .as_ref()
                .map_or(0.0, |gas| gas.report.pressure_final_volume_derivative_pa_m3);
        // The native source derivative holds interval heat fixed. Its endpoint
        // is q1=2*q_stage-q_old, while pressure traction acts at q_stage.
        // This cross-gradient action is generally NONSYMMETRIC.
        let dp = derivative * (2.0 * dv.to_f64());
        let action: Vec<_> = model
            .free_dofs
            .iter()
            .map(|&i| {
                internal[i / 3][i % 3]
                    - geometric[i / 3][i % 3]
                    - dp * self.stage.gradient[i / 3][i % 3]
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
            None => self.exterior_midpoint_pa,
        };
        let mut residual_n = self.endpoint.internal;
        for ((row, gradient), loads) in residual_n
            .iter_mut()
            .zip(&self.endpoint.gradient)
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
