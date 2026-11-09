//! Physical axial-length sensitivities of the complete original DG resonance.
use super::*;
use fs_la::pole_response::PreparedPole;

/// Derivative with respect to uniform logarithmic axial stretch, theta=ln(L/L0).
/// Radii, boundary impedances and every local material law are held fixed;
/// spatial material data and ports move with their original axial stations.
/// These are local derivatives of one simple pole, not a finite-change tuning
/// rule, a mode-tracking result at crossings, or a forced-response peak shift.
#[derive(Clone, Debug, PartialEq)]
pub struct DgPoleLengthSensitivity {
    pub laplace_pole_per_s: C64,
    /// ds/dtheta [1/s]. Its real part changes growth, not just pitch.
    pub laplace_pole_derivative_per_s: C64,
    /// d(omega)/dtheta=i*ds/dtheta under exp(-i omega t) [rad/s].
    pub angular_frequency_derivative_rad_s: C64,
    /// One fixed energy coordinate is one; its derivative is zero. This gauge
    /// is shared by the reference and tangent, not a unit-energy normalization.
    pub normalization_pivot: usize,
    pub reference_energy_coordinates: Vec<C64>,
    pub energy_coordinate_derivative: Vec<C64>,
    /// Derivatives at corresponding MOVING DG nodes, including the derivative
    /// of the physical mass map. A fixed spatial observer needs a transport term.
    pub pressure_derivative_pa: Vec<C64>,
    pub volume_flow_derivative_m3_s: Vec<C64>,
    /// Pressure derivatives at the original, co-moving ports, in original order.
    pub port_pressure_derivatives: Vec<DgPolePort>,
    pub original_relative_residual: f64,
    /// Backward residual of the original bordered derivative equations.
    /// This is not a forward-error bound near a repeated or defective mode.
    pub derivative_relative_residual: f64,
}

impl DgPoleLengthSensitivity {
    /// Signed slope of amplitude decay -Re(s), per unit logarithmic stretch.
    pub fn amplitude_decay_derivative_per_s(&self) -> f64 {
        -self.laplace_pole_derivative_per_s.re
    }

    /// Slope of |Im(s)|/(2*pi), for the returned oscillatory branch.
    /// A real relaxation mode has no differentiated oscillation frequency.
    pub fn oscillation_frequency_derivative_hz(&self) -> Option<f64> {
        (self.laplace_pole_per_s.im != 0.0).then(|| {
            self.laplace_pole_per_s.im.signum() * self.laplace_pole_derivative_per_s.im
                / core::f64::consts::TAU
        })
    }
}

impl<M: DgMedium> PreparedDgPoles<'_, M> {
    /// Differentiate a simple original resonance under x -> exp(theta)*x.
    ///
    /// This changes axial lengths, not cross-sectional radii or material data.
    /// In the energy chart, the pressure/flow propagation and end-boundary
    /// operators scale as exp(-theta). Every LOCAL viscothermal, molecular,
    /// porous and compliant-wall memory operator stays fixed. Thus this method
    /// keeps their full dispersive effect on ds/dtheta instead of assuming
    /// that every lossy resonance simply scales as 1/L.
    ///
    /// The native mass and port maps scale as exp(-theta/2); physical field
    /// tangents therefore include that term in addition to the mode tangent.
    /// Spatial material profiles and side ports are pulled along with their
    /// original relative positions. This does not stretch radii, alter source
    /// composition/temperature, change a load, or model motion in time.
    ///
    /// Uses ONE shared bordered `PreparedPole` factorization and one tangent
    /// solve. The original mode and material band are independently rechecked;
    /// cached pole fields are not trusted. Explicit numerical payload/work
    /// budgets cover the entire query. No state or source history is changed.
    pub fn log_length_sensitivity(
        &self,
        pole: &DgResonancePole,
        options: &PoleResponseOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DgPoleLengthSensitivity, DgPoleError> {
        poll(&mut cancelled)?;
        let n = self.dimension();
        let ports = self.view.ports().len();
        let numerical_options = admit(n, ports, options)?;
        if pole.energy_coordinates.len() != n
            || pole.energy_coordinates.iter().any(|&z| !finite(z))
            || !finite(pole.laplace_pole_per_s)
        {
            return Err(invalid(
                "complete finite original DG pole required for length sensitivity",
            )
            .into());
        }
        if pole.laplace_pole_per_s.im != 0.0 {
            self.view
                .admit_frequency(pole.laplace_pole_per_s.im.abs())?;
        }
        let (prepared, residual_scale) = {
            let evaluated = self.evaluate(
                pole.laplace_pole_per_s,
                options.spectral_scale,
                &mut cancelled,
            )?;
            let prepared = PreparedPole::new(
                pole.laplace_pole_per_s,
                &pole.energy_coordinates,
                &evaluated,
                numerical_options,
                &mut cancelled,
            )?;
            (prepared, evaluated.residual_scale)
        };
        let reference = self.observe(
            NonlinearPole {
                value: prepared.value(),
                vector: prepared.right_vector().to_vec(),
                residual: prepared.pole_residual(),
                relative_step: prepared.relative_correction(),
                residual_scale,
                iterations: 0,
            },
            options.pole_residual_tolerance,
            &mut cancelled,
        )?;
        let action = length_action(
            self.view,
            self.scale,
            prepared.right_vector(),
            &mut cancelled,
        )?;
        let mut tangents = prepared.parameter_derivatives(&[&action], &mut cancelled)?;
        let tangent = tangents.pop().expect("one admitted length direction");
        let physical_tangent: Vec<_> = tangent
            .vector
            .iter()
            .zip(prepared.right_vector())
            .map(|(&dx, &x)| dx - x.scale(0.5))
            .collect();
        let real: Vec<_> = physical_tangent.iter().map(|x| x.re).collect();
        let imag: Vec<_> = physical_tangent.iter().map(|x| x.im).collect();
        let (rp, rq) = self.view.duct.from_energy_coords(&real);
        let (ip, iq) = self.view.duct.from_energy_coords(&imag);
        let pressure: Vec<_> = rp
            .into_iter()
            .zip(ip)
            .map(|(r, i)| C64::new(r, i))
            .collect();
        let flow: Vec<_> = rq
            .into_iter()
            .zip(iq)
            .map(|(r, i)| C64::new(r, i))
            .collect();
        let mut port_derivatives = Vec::with_capacity(ports);
        for (id, b) in self.view.ports() {
            poll(&mut cancelled)?;
            let mut sum = C64::ZERO;
            let mut correction = C64::ZERO;
            for (&b, &x) in b.iter().zip(&physical_tangent) {
                add(&mut sum, &mut correction, x.scale(b));
            }
            port_derivatives.push(DgPolePort {
                id: *id,
                pressure_pa: sum,
            });
        }
        if pressure.iter().chain(&flow).any(|&x| !finite(x))
            || port_derivatives.iter().any(|p| !finite(p.pressure_pa))
        {
            return Err(invalid("DG length field sensitivity is unrepresentable").into());
        }
        poll(&mut cancelled)?;
        Ok(DgPoleLengthSensitivity {
            laplace_pole_per_s: prepared.value(),
            laplace_pole_derivative_per_s: tangent.value,
            angular_frequency_derivative_rad_s: C64::new(-tangent.value.im, tangent.value.re),
            normalization_pivot: prepared.normalization_pivot(),
            reference_energy_coordinates: prepared.right_vector().to_vec(),
            energy_coordinate_derivative: tangent.vector,
            pressure_derivative_pa: pressure,
            volume_flow_derivative_m3_s: flow,
            port_pressure_derivatives: port_derivatives,
            original_relative_residual: reference.relative_residual,
            derivative_relative_residual: tangent.linear_residual,
        })
    }
}

fn admit(
    n: usize,
    ports: usize,
    options: &PoleResponseOptions,
) -> Result<PoleResponseOptions, DgPoleError> {
    if n == 0
        || n > MAX_NONLINEAR_POLE_DIMENSION
        || n > options.max_dimension
        || !(1..=256).contains(&options.max_vectors)
    {
        return Err(invalid(
            "complete DG length tangent exceeds the declared dimension/vector capacity",
        )
        .into());
    }
    let (n, p) = (n as u128, ports as u128);
    let bytes = 32 * n * n + 512 * n + 128 * p;
    let work = 64 * n * n + 32 * n * p;
    if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
        return Err(invalid("DG length tangent maps exceed their payload/work capacity").into());
    }
    let mut numerical = *options;
    numerical.max_dense_bytes -= bytes as usize;
    numerical.max_work -= work as u64;
    let d = n + 1;
    let k = options.max_vectors as u128;
    if 128 * d * d + 64 * d * (k + 4) > numerical.max_dense_bytes as u128
        || 32 * d * d * d + 128 * (k + 1) * d * d > numerical.max_work as u128
    {
        return Err(invalid("DG length tangent factor exceeds its payload/work capacity").into());
    }
    Ok(numerical)
}

fn add(sum: &mut C64, correction: &mut C64, value: C64) {
    let term = value - *correction;
    let next = *sum + term;
    *correction = (next - *sum) - term;
    *sum = next;
}

// (partial T / partial theta)*x = (L_propagation + L_boundary)*x / scale.
// Local memories do not cross the pressure/axial-flow blocks. Boundary loss
// is reconstructed from the original physical boundary law, never inferred
// by subtracting or symmetrizing material damping in the complete generator.
fn length_action<M: DgMedium>(
    view: DgContinuous<'_, M>,
    scale: f64,
    x: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, DgPoleError> {
    let duct = view.duct;
    let (n, m) = (view.dimension(), duct.field_dofs);
    let mut action = vec![C64::ZERO; n];
    for i in 0..2 * m {
        poll(cancelled)?;
        let mut correction = C64::ZERO;
        let range = if i < m { m..2 * m } else { 0..m };
        for j in range {
            let original = view.generator_per_s()[i * n + j];
            let scaled = original / scale;
            if original != 0.0 && scaled == 0.0 {
                return Err(invalid("DG length propagation derivative underflowed").into());
            }
            add(&mut action[i], &mut correction, x[j].scale(scaled));
        }
    }
    let local = duct.mesh.order + 1;
    for (bc, element, node, medium) in [
        (duct.inlet, 0, 0, &duct.gas),
        (
            duct.outlet,
            duct.elements.len() - 1,
            local - 1,
            &duct.outlet_gas,
        ),
    ] {
        poll(cancelled)?;
        let e = &duct.elements[element];
        let z = medium.density_kg_m3() * medium.sound_speed_m_s() / e.areas[node];
        if let Some((alpha, beta)) = bc.characteristic_weights(z) {
            for (offset, inverse_root, weight) in [
                (0, &e.mp_inv_sqrt, beta / z),
                (m, &e.mq_inv_sqrt, alpha * z),
            ] {
                let start = offset + element * local;
                for i in 0..local {
                    let mut correction = C64::ZERO;
                    for j in 0..local {
                        let original = -weight
                            * inverse_root[i * local + node]
                            * inverse_root[j * local + node];
                        let a = original / scale;
                        if !a.is_finite() || (original != 0.0 && a == 0.0) {
                            return Err(invalid(
                                "DG length boundary derivative is unrepresentable",
                            )
                            .into());
                        }
                        add(
                            &mut action[start + i],
                            &mut correction,
                            x[start + j].scale(a),
                        );
                    }
                }
            }
        }
    }
    if action.iter().any(|&z| !finite(z)) {
        return Err(invalid("DG length derivative action is unrepresentable").into());
    }
    Ok(action)
}
