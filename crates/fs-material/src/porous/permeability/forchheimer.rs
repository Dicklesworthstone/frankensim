//! Calibrated superficial-velocity inertia in the original porous owner.
//!
//! The convention is -dp/dx = mu*q/k + rho*beta*|q|q, with beta in 1/m,
//! without an implicit factor of one half. beta is independent input, never
//! inferred from k, porosity, particle diameter or an unqualified Ergun law.
//! See <https://www.comsol.com/blogs/modeling-darcian-and-non-darcian-flow-in-porous-media>.
//!
//! At constant ideal-gas temperature/composition and viscosity, j = rho*q is
//! constant and integration gives
//! (p_left-p_right)*p_mean/(R*T) = mu*sum(L/k)*j + sum(beta*L)*|j|j.
//! This is a drag closure, not a compressible momentum-acceleration, choking,
//! heat-transfer or turbulence model. The explicit calibration bounds apply.

use super::{
    DarcyFlowConstraint, DarcyGasBarrier, GasState, IntrinsicPermeability, MAX_POROUS_LAYERS,
    MaterialError, PorousGasLayer, PorousGasSupport, dot, invalid, nonnegative, positive, unit,
};

mod source;

/// Explicit, immutable calibration of beta [1/m] in rho*beta*|q|q.
///
/// q is superficial velocity through the total area, not pore velocity. The
/// calibration binds the specimen normal (either flow sign), transverse
/// constraint, provenance and admitted absolute T/p/speed domain. It does not
/// claim universality across material directions, fluids or experimental
/// conditions absent from the supplied evidence. No gas slip or saturation
/// qualification is inferred. A zero coefficient retains the exact Darcy path.
#[derive(Clone, Debug, PartialEq)]
pub struct ForchheimerCalibration {
    beta_per_m: f64,
    normal: [f64; 3],
    flow_constraint: DarcyFlowConstraint,
    support: PorousGasSupport,
    provenance: String,
    source: Option<std::sync::Arc<source::ForchheimerSource>>,
}

impl ForchheimerCalibration {
    /// The normal and constraint must match the original barrier exactly;
    /// reversing the normal is permitted by the symmetric signed drag law.
    pub fn try_new(
        beta_per_m: f64,
        normal: [f64; 3],
        flow_constraint: DarcyFlowConstraint,
        support: PorousGasSupport,
        provenance: impl Into<String>,
    ) -> Result<Self, MaterialError> {
        support.validate()?;
        let provenance = provenance.into();
        if !nonnegative(beta_per_m)
            || !unit(normal)
            || provenance.trim().is_empty()
            || provenance.len() > 4096
        {
            return Err(invalid(
                "Forchheimer calibration needs nonnegative beta [1/m], a unit normal and bounded provenance",
            ));
        }
        Ok(Self {
            beta_per_m,
            normal,
            flow_constraint,
            support,
            provenance,
            source: None,
        })
    }

    pub const fn beta_per_m(&self) -> f64 {
        self.beta_per_m
    }
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }
    pub const fn flow_constraint(&self) -> DarcyFlowConstraint {
        self.flow_constraint
    }
    pub const fn support(&self) -> PorousGasSupport {
        self.support
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    fn admit_geometry(
        &self,
        tensor: &IntrinsicPermeability,
        normal: [f64; 3],
        constraint: DarcyFlowConstraint,
    ) -> Result<(), MaterialError> {
        if self.flow_constraint != constraint
            || normal != self.normal && normal != self.normal.map(|x| -x)
        {
            return Err(invalid(
                "Forchheimer calibration does not cover this normal/transverse constraint",
            ));
        }
        if self.beta_per_m > 0.0
            && constraint == DarcyFlowConstraint::NormalPressureGradient
            && !is_eigendirection(tensor, normal)?
        {
            return Err(invalid(
                "scalar Forchheimer normal-gradient flow requires a permeability eigen-direction; use explicitly calibrated confinement",
            ));
        }
        Ok(())
    }
}

// A normal pressure gradient produces K*n. The scalar |q_n|q_n law is valid
// only when that vector is normal. Scale its actual contributing components,
// not an unrelated large/small eigenvalue, so high contrast cannot hide
// transverse flow. The tolerance covers the admitted orthonormal-frame error.
fn is_eigendirection(
    tensor: &IntrinsicPermeability,
    normal: [f64; 3],
) -> Result<bool, MaterialError> {
    let projections = tensor.axes.map(|axis| dot(axis, normal));
    let mut velocity = [0.0_f64; 3];
    for i in 0..3 {
        velocity[i] = tensor.principal_m2[i] * projections[i];
    }
    let scale = velocity.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    if !nonnegative(scale) {
        return Err(invalid("Forchheimer tensor direction is unrepresentable"));
    }
    if scale == 0.0 {
        return Ok(true);
    }
    for component in &mut velocity {
        *component /= scale;
    }
    let axial = dot(velocity, projections);
    Ok((0..3).all(|i| (velocity[i] - axial * projections[i]).abs() <= 16e-12))
}

impl PorousGasLayer {
    /// Attach an explicit directional calibration without changing k, its
    /// source bundle, layer thickness or existing provenance.
    pub fn with_forchheimer_calibration(
        &self,
        calibration: ForchheimerCalibration,
    ) -> Result<Self, MaterialError> {
        calibration.admit_layer_source(self)?;
        if calibration.beta_per_m > 0.0
            && (self.klinkenberg_slip_pa != 0.0 || self.wet_state.is_some())
        {
            return Err(invalid(
                "Forchheimer calibration does not qualify gas slip or wet-pore inertia",
            ));
        }
        let mut layer = self.clone();
        layer.forchheimer = Some(calibration);
        Ok(layer)
    }

    pub const fn forchheimer_calibration(&self) -> Option<&ForchheimerCalibration> {
        self.forchheimer.as_ref()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct InertialGeometry {
    /// sum(beta_i*L_i), dimensionless; exactly zero is the Darcy limit.
    pub resistance: f64,
    pub fractions: [f64; MAX_POROUS_LAYERS],
}

pub(super) fn compile(
    layers: &[PorousGasLayer],
    normal: [f64; 3],
    constraint: DarcyFlowConstraint,
    support: &mut PorousGasSupport,
) -> Result<InertialGeometry, MaterialError> {
    let mut geometry = InertialGeometry {
        resistance: 0.0,
        fractions: [0.0; MAX_POROUS_LAYERS],
    };
    for (i, layer) in layers.iter().enumerate() {
        let Some(calibration) = &layer.forchheimer else {
            continue;
        };
        calibration.admit_geometry(&layer.permeability, normal, constraint)?;
        let domain = calibration.support;
        support.temperature_k = (
            support.temperature_k.0.max(domain.temperature_k.0),
            support.temperature_k.1.min(domain.temperature_k.1),
        );
        support.pressure_pa = (
            support.pressure_pa.0.max(domain.pressure_pa.0),
            support.pressure_pa.1.min(domain.pressure_pa.1),
        );
        support.maximum_superficial_speed_m_s = support
            .maximum_superficial_speed_m_s
            .min(domain.maximum_superficial_speed_m_s);
        let part = calibration.beta_per_m * layer.thickness_m;
        if !nonnegative(part)
            || calibration.beta_per_m > 0.0 && part == 0.0
            || !nonnegative(geometry.resistance + part)
        {
            return Err(invalid(
                "Forchheimer layer series resistance is unrepresentable",
            ));
        }
        geometry.fractions[i] = part;
        geometry.resistance += part;
    }
    if geometry.resistance > 0.0 {
        for fraction in &mut geometry.fractions[..layers.len()] {
            let nonzero = *fraction > 0.0;
            *fraction /= geometry.resistance;
            if nonzero && !positive(*fraction) {
                return Err(invalid("Forchheimer layer fraction is unrepresentable"));
            }
        }
    }
    Ok(geometry)
}

impl DarcyGasBarrier {
    /// Original sum(beta_i*L_i), dimensionless. There is no hidden coefficient
    /// or factor of one half. Zero remains the exact viscous Darcy limit.
    pub const fn inertial_geometric_resistance(&self) -> f64 {
        self.inertial_geometry.resistance
    }

    /// Attach one explicit calibration to each unchanged original layer.
    /// Select the desired transverse constraint before attaching calibrations;
    /// changing it afterwards requires new calibration evidence.
    pub fn with_forchheimer_calibrations(
        &self,
        calibrations: &[ForchheimerCalibration],
    ) -> Result<Self, MaterialError> {
        if calibrations.len() != self.layers.len() {
            return Err(invalid(
                "Forchheimer calibration must cover every original layer exactly once",
            ));
        }
        let layers = self
            .layers
            .iter()
            .zip(calibrations)
            .map(|(layer, calibration)| layer.with_forchheimer_calibration(calibration.clone()))
            .collect::<Result<Vec<_>, _>>()?;
        Self::compile(
            layers,
            self.normal,
            self.support,
            self.slip_gas_model.clone(),
            self.flow_constraint,
        )
    }

    /// Frozen-state coefficients (A, B) in dp = A*Q + B*Q*|Q|, using the
    /// supplied gas owner's actual mu and rho. A has units Pa s/m3, B Pa s2/m6.
    ///
    /// This is a local low-compressibility approximation, not the compressible
    /// mass-flow solve. The consumer must enforce its pressure-drop assumption
    /// and actual solved superficial-speed/face-pressure support. None means
    /// an exact seal or zero area. Gas slip refuses instead of being frozen.
    pub fn volume_resistance_coefficients(
        &self,
        gas: &GasState,
        area_m2: f64,
    ) -> Result<Option<(f64, f64)>, MaterialError> {
        self.support.admit(gas, gas.pressure, gas.pressure)?;
        if !nonnegative(area_m2) || self.maximum_slip_pa != 0.0 {
            return Err(invalid(
                "frozen porous resistance needs nonnegative area and no gas slip",
            ));
        }
        let Some(geometry) = self.geometric_resistance_per_m.filter(|_| area_m2 > 0.0) else {
            return Ok(None);
        };
        let linear = gas.dynamic_viscosity * geometry / area_m2;
        let quadratic = (gas.density * self.inertial_geometry.resistance / area_m2) / area_m2;
        if !positive(linear)
            || !nonnegative(quadratic)
            || self.inertial_geometry.resistance > 0.0 && quadratic == 0.0
        {
            return Err(invalid(
                "frozen porous volume resistance is unrepresentable",
            ));
        }
        Ok(Some((linear, quadratic)))
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct QuadraticFlow {
    pub flow: f64,
}

// Positive root of a*q+b*q²=c. Scaled dominant-term branches avoid a²,
// b*c, subtractive cancellation and an overflowing unopposed Darcy flow.
pub(super) fn positive_quadratic_flow(
    linear: f64,
    quadratic: f64,
    driving: f64,
) -> Result<QuadraticFlow, MaterialError> {
    if !positive(linear) || !nonnegative(quadratic) || !positive(driving) {
        return Err(invalid(
            "Forchheimer flow needs finite positive viscous resistance/driving and nonnegative inertia",
        ));
    }
    let flow = if quadratic == 0.0 {
        driving / linear
    } else {
        let sqrt_driving = fs_math::det::sqrt(driving);
        let sqrt_quadratic = fs_math::det::sqrt(quadratic);
        let cross = sqrt_driving * sqrt_quadratic;
        if !positive(cross) {
            return Err(invalid("Forchheimer quadratic scale is unrepresentable"));
        }
        if cross <= linear {
            let ratio = cross / linear;
            let denominator = 0.5 + fs_math::det::sqrt(0.25 + ratio * ratio);
            let unopposed = driving / linear;
            let flow = if unopposed.is_finite() {
                unopposed / denominator
            } else {
                (driving / denominator) / linear
            };
            flow
        } else {
            let ratio = linear / cross;
            let denominator = 0.5 * ratio + fs_math::det::sqrt(1.0 + 0.25 * ratio * ratio);
            let unopposed = sqrt_driving / sqrt_quadratic;
            let flow = if unopposed.is_finite() {
                unopposed / denominator
            } else {
                (sqrt_driving / denominator) / sqrt_quadratic
            };
            flow
        }
    };
    if !positive(flow) {
        return Err(invalid("Forchheimer flow is unrepresentable"));
    }
    Ok(QuadraticFlow { flow })
}

// area/(a+b*|j|), without overflowing the unopposed 1/a or the positive
// denominator sum. The quadratic solution, rather than an inadmissible Darcy
// trial, determines both the actual speed and the secant conductance.
pub(super) fn secant_volume_conductance(
    linear: f64,
    quadratic: f64,
    flux: f64,
    area: f64,
) -> Result<f64, MaterialError> {
    let inertial = quadratic * flux;
    if !positive(linear) || !nonnegative(inertial) || !positive(area) {
        return Err(invalid(
            "Forchheimer secant conductance coefficients are unrepresentable",
        ));
    }
    let scale = linear.max(inertial);
    let sum = linear / scale + inertial / scale;
    let ratio = area / scale;
    let conductance = if ratio.is_finite() {
        ratio / sum
    } else {
        (area / sum) / scale
    };
    if !positive(conductance) {
        return Err(invalid("Forchheimer secant conductance is unrepresentable"));
    }
    Ok(conductance)
}

pub(super) fn loss_fractions(
    linear: f64,
    quadratic: f64,
    flow: f64,
) -> Result<(f64, f64), MaterialError> {
    let inertial = quadratic * flow;
    if !positive(linear) || !nonnegative(inertial) {
        return Err(invalid(
            "Forchheimer pressure fractions are unrepresentable",
        ));
    }
    let scale = linear.max(inertial);
    let viscous = linear / scale;
    let inertial = inertial / scale;
    let sum = viscous + inertial;
    Ok((viscous / sum, inertial / sum))
}

#[cfg(test)]
mod tests;
