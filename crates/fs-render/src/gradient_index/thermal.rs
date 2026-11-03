//! Passive LTE source quadrature on the native guarded curved-ray stages.
//!
//! At fixed vacuum wavelength the reduced-radiance source is alpha B(T),
//! because the local n² in Kirchhoff emission cancels that in L/n². Integrate
//! alpha(s) B(T(s)) exp(-tau_abs(s)) along the same ray, using true absorption
//! only. Molecular free-flight sampling already accounts for scattering
//! survival. The camera applies its arriving Mueller intensity and n_camera².
//! See Zhao, Tan and Liu (2011), Eq. 21, <https://arxiv.org/abs/1110.5134>.
//!
//! Temperature is an immutable continuous physical reconstruction. Local
//! embedded estimates control quadrature and optics together; they do not
//! certify a global error or detect arbitrary unresolved spatial features.

use super::event::{
    BoundaryOracle, EventStepper, GradientBoundaryPropagation, GradientDeparture,
    GradientEventOptions, OpticalStepper, propagate_with_stepper,
};
use super::*;
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;

/// Original frozen thermodynamic temperature, in kelvin, at a world position.
/// A provider must refuse outside its supported domain; it must not clamp,
/// extrapolate, interpolate unrelated cell averages or infer T from index.
/// Smoothness and the selected optical maximum step must resolve its spatial
/// structure. Every nonzero-absorption query occurs at a guarded optical stage.
pub trait TemperatureField {
    fn temperature_k(&self, position_m: Point3) -> Result<f64, GradientIndexError>;
}

impl<F> TemperatureField for F
where
    F: Fn(Point3) -> Result<f64, GradientIndexError>,
{
    fn temperature_k(&self, position_m: Point3) -> Result<f64, GradientIndexError> {
        self(position_m)
    }
}

// Uniform optical coefficients stay prepared while only temperature is
// queried at the adaptive stages. Nonuniform optics retain their exact owner.
struct OpticalField<'a, F: GradientIndexField + ?Sized> {
    field: &'a F,
    uniform: Option<GradientIndexSample>,
}
impl<F: GradientIndexField + ?Sized> GradientIndexField for OpticalField<'_, F> {
    fn sample(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        match self.uniform {
            Some(sample) => Ok(sample),
            None => self.field.sample(p, wavelength_nm),
        }
    }
}

struct ThermalStepper<'a, F: GradientIndexField + ?Sized> {
    optical: OpticalStepper<'a, F>,
    field: OpticalField<'a, F>,
    temperature: &'a dyn TemperatureField,
    ray: Ray,
    wavelength_nm: f64,
    absolute_tolerance: f64,
    integral: f64,
    derivative: f64,
}

struct ThermalTrial {
    optical: Trial,
    increment: f64,
    error: f64,
    derivative: f64,
}

impl<'a, F: GradientIndexField + ?Sized> ThermalStepper<'a, F> {
    fn new(
        field: &'a F,
        temperature: &'a dyn TemperatureField,
        ray: Ray,
        axis: Vec3,
        wavelength_nm: f64,
        absolute_tolerance: f64,
    ) -> Result<Self, GradientIndexError> {
        if !absolute_tolerance.is_finite() || absolute_tolerance <= 0.0 {
            return Err(GradientIndexError::InvalidOptions);
        }
        let optical = OpticalStepper::new(field, ray, axis, wavelength_nm)?;
        let initial = optical.initial_sample();
        let derivative = source_derivative(temperature, ray.origin, wavelength_nm, initial, 0.0)?;
        if !derivative.is_finite() {
            return Err(GradientIndexError::NumericalRange);
        }
        Ok(Self {
            field: OpticalField {
                field,
                uniform: optical.is_uniform().then_some(initial),
            },
            optical,
            temperature,
            ray,
            wavelength_nm,
            absolute_tolerance,
            integral: 0.0,
            derivative,
        })
    }
}

impl<F: GradientIndexField + ?Sized> EventStepper for ThermalStepper<'_, F> {
    type Trial = ThermalTrial;

    fn state(&self) -> State {
        self.optical.state()
    }
    fn initial_sample(&self) -> GradientIndexSample {
        self.optical.initial_sample()
    }
    fn final_sample(&self) -> GradientIndexSample {
        self.optical.final_sample()
    }
    fn is_uniform(&self) -> bool {
        // A provably transparent homogeneous medium has no source even if its
        // declared temperature is unknown or inadmissible. Keep its exact path.
        self.optical.is_uniform() && self.initial_sample().absorption_per_m() == 0.0
    }
    fn trial_state(trial: &ThermalTrial) -> State {
        trial.optical.state
    }

    fn trial(
        &self,
        step: f64,
        new_distance: f64,
        guard: &mut impl FnMut(Point3, Vec3) -> bool,
        evaluations: &mut u32,
    ) -> Result<ThermalTrial, GradientIndexError> {
        if self.is_uniform() {
            return Ok(ThermalTrial {
                optical: self.optical.trial(step, new_distance, guard, evaluations)?,
                increment: 0.0,
                error: 0.0,
                derivative: 0.0,
            });
        }
        let state = self.state();
        let track_group = self.initial_sample().group_index().is_some();
        let mut q = [0.0; 7];
        q[0] = self.derivative;
        let mut stage_index = 1;
        let optical = trial_step_guarded_observed(
            &self.field,
            self.ray.origin,
            self.wavelength_nm,
            &state,
            rhs(&state, self.final_sample(), track_group)?,
            step,
            track_group,
            guard,
            evaluations,
            |point, stage, sample| {
                q[stage_index] = source_derivative(
                    self.temperature,
                    point,
                    self.wavelength_nm,
                    sample,
                    stage[11],
                )?;
                stage_index += 1;
                Ok(())
            },
        )?;
        let mut increment = 0.0;
        let mut error = 0.0;
        for i in 0..6 {
            increment += A[5][i] * q[i];
        }
        for i in 0..7 {
            error += ERROR[i] * q[i];
        }
        Ok(ThermalTrial {
            optical,
            increment: step * increment,
            error: step * error,
            derivative: q[6],
        })
    }

    fn normalized_error(
        &self,
        trial: &ThermalTrial,
        step: f64,
        options: GradientIndexOptions,
    ) -> Result<f64, GradientIndexError> {
        // Negative/overflowing passive quadrature is an unresolved trial, not
        // a physical source clamp. Refine with the same native step budget.
        if !trial.increment.is_finite()
            || trial.increment < 0.0
            || !trial.error.is_finite()
            || !trial.derivative.is_finite()
        {
            return Ok(f64::MAX);
        }
        let scale = self.absolute_tolerance + options.relative_tolerance * trial.increment;
        let source_error = trial.error.abs() / scale;
        Ok(self
            .optical
            .normalized_error(&trial.optical, step, options)?
            .max(source_error))
    }

    fn accept(&mut self, trial: ThermalTrial) -> Result<(), GradientIndexError> {
        let total = self.integral + trial.increment;
        if !total.is_finite() {
            return Err(GradientIndexError::NumericalRange);
        }
        self.optical.accept(trial.optical)?;
        self.integral = total;
        self.derivative = trial.derivative;
        Ok(())
    }
}

fn source_derivative(
    temperature: &dyn TemperatureField,
    position: Point3,
    wavelength_nm: f64,
    sample: GradientIndexSample,
    absorption_depth: f64,
) -> Result<f64, GradientIndexError> {
    let absorption = sample.absorption_per_m();
    // Kirchhoff emission is exactly zero before any temperature or Planck
    // query. Transparent original optics must not acquire a thermal refusal.
    if absorption == 0.0 {
        return Ok(0.0);
    }
    let kelvin = temperature.temperature_k(position)?;
    if !kelvin.is_finite() || kelvin <= 0.0 {
        return Err(GradientIndexError::InvalidSample);
    }
    let blackbody = planck_radiance_w_m2_sr_nm(kelvin, wavelength_nm)
        .map_err(|_| GradientIndexError::NumericalRange)?;
    if blackbody == 0.0 {
        return Ok(0.0);
    }
    let direct = absorption * blackbody * det::exp(-absorption_depth);
    if direct.is_finite() && direct > 0.0 {
        return Ok(direct);
    }
    // Beer survival can rescue an overflowing alpha*B product. Intermediate
    // RK depths need not be positive; an overflowing trial is refined above.
    Ok(det::exp(
        det::ln(absorption) + det::ln(blackbody) - absorption_depth,
    ))
}

/// Native event walk plus reduced source integral, before camera prefix/index.
#[allow(clippy::too_many_arguments)]
pub(crate) fn propagate_emission_to_boundary_with_control<F, O>(
    field: &F,
    temperature: &dyn TemperatureField,
    absolute_tolerance: f64,
    oracle: &mut O,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength_nm: f64,
    max_distance_m: f64,
    departure: Option<GradientDeparture>,
    options: GradientEventOptions,
    mut control: impl FnMut() -> bool,
) -> Result<(GradientBoundaryPropagation, f64), GradientIndexError>
where
    F: GradientIndexField + ?Sized,
    O: BoundaryOracle + ?Sized,
{
    let (propagation, stepper) = propagate_with_stepper(
        oracle,
        ray,
        transverse_axis,
        wavelength_nm,
        max_distance_m,
        departure,
        options,
        &mut control,
        |ray, axis| {
            ThermalStepper::new(
                field,
                temperature,
                ray,
                axis,
                wavelength_nm,
                absolute_tolerance,
            )
        },
    )?;
    Ok((propagation, stepper.integral))
}

/// Borrow an original source-resolved linear primitive gas field as temperature.
/// World positions project onto the caller's explicit unit axial direction;
/// radial uniformity is declared by this embedding, never inferred from optics.
/// Every query uses the material owner's `state_at` and retains its phase and
/// model-validity refusals. Cellwise-constant reconstructions are refused.
///
/// T is continuous and piecewise linear at original stations. Select/refine
/// the optical maximum step against those station spacings; the embedded local
/// error is not an accuracy guarantee across unresolved primitive spans.
#[cfg(feature = "material-state")]
pub struct AxialGasTemperatureField<'a> {
    field: &'a fs_material::gas_field::AxialGasField,
    origin: Point3,
    axis: Vec3,
}

#[cfg(feature = "material-state")]
impl<'a> AxialGasTemperatureField<'a> {
    pub fn try_new(
        field: &'a fs_material::gas_field::AxialGasField,
        origin: Point3,
        axis: Vec3,
    ) -> Result<Self, GradientIndexError> {
        if !finite_point(origin)
            || field.sampling() != fs_material::gas_field::AxialGasSampling::LinearPrimitives
        {
            return Err(GradientIndexError::InvalidSample);
        }
        let axis = unit(axis).ok_or(GradientIndexError::InvalidSample)?;
        Ok(Self {
            field,
            origin,
            axis,
        })
    }

    /// Original thermochemical owner and supplied primitive station map.
    pub const fn source(&self) -> &fs_material::gas_field::AxialGasField {
        self.field
    }
}

#[cfg(feature = "material-state")]
impl TemperatureField for AxialGasTemperatureField<'_> {
    fn temperature_k(&self, position_m: Point3) -> Result<f64, GradientIndexError> {
        let x = Vec3::new(
            position_m.x - self.origin.x,
            position_m.y - self.origin.y,
            position_m.z - self.origin.z,
        )
        .dot(self.axis);
        if !finite_point(position_m) || !x.is_finite() || x < 0.0 || x > self.field.length_m() {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        self.field
            .state_at(x)
            .map(|state| state.state().temperature)
            .map_err(|_| GradientIndexError::InvalidSample)
    }
}

#[cfg(test)]
mod tests;
