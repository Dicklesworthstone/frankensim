//! Point-junction momentum balance in the original owner's pressure chart.
//! The original hydraulic node has no fluid storage. These forces therefore
//! add no second fluid mass, pressure-work term or transported enthalpy owner.

use super::pipe::{
    PipeLoad, PipeLoadError, affine_mean_cube, affine_mean_square, check, cross, dot, number,
    positive, sub, vector,
};

type Result<T> = core::result::Result<T, PipeLoadError>;

/// Fluid reaction on one original zero-storage junction. Every incident pipe
/// supplies exactly one port, with its normal directed from the node into the
/// pipe. Its pressure traction and momentum flux oppose that pipe endpoint's
/// contributions exactly; do not apply this load twice when loading an assembly.
#[derive(Clone, Debug, PartialEq)]
pub struct JunctionLoad {
    pub node: usize,
    pub load: PipeLoad,
    /// Signed OUTWARD kinetic-energy transport over the original affine-flow
    /// interval: integral sum(rho*q_out^3/(2*A^2)) dt. In the default static
    /// chart, pressure work cancels through native volume continuity.
    ///
    /// This is an UNRESOLVED energy term in the existing point-junction model,
    /// not a nonnegative loss or an automatic thermal credit. Positive values
    /// require an energy supply or revised pressure/thermal closure; negative
    /// values require a declared destination. Zero net transport during reversal
    /// does not erase the nonzero Q^2 momentum exchange in the same interval.
    /// The explicit total-pressure owner instead balances this transport with
    /// per-port pressure work in its retained NodeStreamWork report.
    pub kinetic_energy_transport_j: f64,
}

#[derive(Clone, Copy)]
pub(super) struct PortInterval {
    pub outward_normal: [f64; 3],
    pub area_m2: f64,
    pub outflows_m3_s: [f64; 2],
}

#[derive(Clone, Copy)]
pub(super) struct JunctionState {
    pub node: usize,
    pub position_m: [f64; 3],
    pub about_m: [f64; 3],
    pub duration_s: f64,
    pub density_kg_m3: f64,
    pub pressure_pa: f64,
    pub exterior_pressure_pa: f64,
}

/// All rates, pressure and density are supplied by the original hydraulic
/// trial. Geometry admission and complete incident-port coverage live there.
pub(super) fn evaluate(
    state: JunctionState,
    ports: &[PortInterval],
    cancelled: impl FnMut() -> bool,
) -> Result<JunctionLoad> {
    evaluate_in_chart(state, ports, super::super::PressureChart::Static, cancelled)
}

pub(super) fn evaluate_in_chart(
    state: JunctionState,
    ports: &[PortInterval],
    chart: super::super::PressureChart,
    mut cancelled: impl FnMut() -> bool,
) -> Result<JunctionLoad> {
    check(&mut cancelled)?;
    if ports.is_empty()
        || !positive(state.duration_s)
        || !positive(state.density_kg_m3)
        || state
            .position_m
            .iter()
            .chain(&state.about_m)
            .any(|x| !x.is_finite())
    {
        return Err(PipeLoadError::Invalid(
            "finite complete junction and positive time/density required",
        ));
    }
    let pressure = number(state.pressure_pa - state.exterior_pressure_pa)?;
    let mut flow = [Sum::default(); 2];
    let mut magnitude = [Sum::default(); 2];
    let mut impulse = [Sum::default(); 3];
    let mut kinetic = Sum::default();
    for port in ports {
        check(&mut cancelled)?;
        if !positive(port.area_m2)
            || port
                .outward_normal
                .iter()
                .chain(&port.outflows_m3_s)
                .any(|x| !x.is_finite())
            || (dot(port.outward_normal, port.outward_normal) - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(PipeLoadError::Invalid(
                "positive area, unit outward normal and finite port flow required",
            ));
        }
        for endpoint in 0..2 {
            flow[endpoint].add(port.outflows_m3_s[endpoint])?;
            magnitude[endpoint].add(port.outflows_m3_s[endpoint].abs())?;
        }
        let square = affine_mean_square(port.outflows_m3_s)?;
        let cube = affine_mean_cube(port.outflows_m3_s)?;
        let dynamic_head = if chart == super::super::PressureChart::Total {
            number(0.5 * state.density_kg_m3 * (square / port.area_m2) / port.area_m2)?
        } else {
            0.0
        };
        let static_pressure = if chart == super::super::PressureChart::Total {
            number(number(state.pressure_pa - dynamic_head)? - state.exterior_pressure_pa)?
        } else {
            pressure
        };
        let pressure_force = number(static_pressure * port.area_m2)?;
        let momentum_force = number(state.density_kg_m3 * (square / port.area_m2))?;
        let normal_impulse = number(-state.duration_s * number(pressure_force + momentum_force)?)?;
        for axis in 0..3 {
            impulse[axis].add(number(normal_impulse * port.outward_normal[axis])?)?;
        }
        kinetic.add(number(
            0.5 * state.duration_s * state.density_kg_m3 * (cube / port.area_m2) / port.area_m2,
        )?)?;
    }
    // Endpoint continuity closes the whole affine trajectory, including a
    // reversal. A pressure datum with an unrepresented reservoir stream refuses.
    for endpoint in 0..2 {
        if flow[endpoint].value()?.abs() > 64.0 * f64::EPSILON * magnitude[endpoint].value()? {
            return Err(PipeLoadError::Invalid(
                "junction has unrepresented mass transfer or storage",
            ));
        }
    }
    let impulse = [
        impulse[0].value()?,
        impulse[1].value()?,
        impulse[2].value()?,
    ];
    let angular = vector(cross(
        vector(sub(state.position_m, state.about_m))?,
        impulse,
    ))?;
    let load = PipeLoad {
        about_m: state.about_m,
        duration_s: state.duration_s,
        force_n: vector(impulse.map(|x| x / state.duration_s))?,
        torque_nm: vector(angular.map(|x| x / state.duration_s))?,
        impulse_ns: impulse,
        angular_impulse_nms: angular,
    };
    check(&mut cancelled)?;
    Ok(JunctionLoad {
        node: state.node,
        load,
        kinetic_energy_transport_j: kinetic.value()?,
    })
}

/// Deterministic compensated sums retain small loads beside large components.
#[derive(Clone, Copy, Default)]
pub(super) struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    pub(super) fn add(&mut self, x: f64) -> Result<()> {
        let (sum, error) = fs_math::eft::two_sum(self.sum, number(x)?);
        self.sum = number(sum)?;
        self.correction = number(self.correction + error)?;
        Ok(())
    }
    pub(super) fn value(self) -> Result<f64> {
        number(self.sum + self.correction)
    }
}

pub(super) fn sum_loads<'a>(
    duration_s: f64,
    about_m: [f64; 3],
    loads: impl IntoIterator<Item = &'a PipeLoad>,
) -> Result<PipeLoad> {
    let mut impulse = [Sum::default(); 3];
    let mut angular = [Sum::default(); 3];
    for load in loads {
        if load.duration_s != duration_s || load.about_m != about_m {
            return Err(PipeLoadError::Invalid(
                "one interval and moment origin required for an assembly",
            ));
        }
        for axis in 0..3 {
            impulse[axis].add(load.impulse_ns[axis])?;
            angular[axis].add(load.angular_impulse_nms[axis])?;
        }
    }
    let impulse = [
        impulse[0].value()?,
        impulse[1].value()?,
        impulse[2].value()?,
    ];
    let angular = [
        angular[0].value()?,
        angular[1].value()?,
        angular[2].value()?,
    ];
    Ok(PipeLoad {
        about_m,
        duration_s,
        force_n: vector(impulse.map(|x| x / duration_s))?,
        torque_nm: vector(angular.map(|x| x / duration_s))?,
        impulse_ns: impulse,
        angular_impulse_nms: angular,
    })
}

#[cfg(test)]
mod tests;
