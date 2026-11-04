//! Original port residues of an isolated simple mode, without fitting a response.
use super::*;
use fs_la::pole_response::PreparedPole;

/// Principal part of the original pressure/volume-flow transfer matrix.
///
/// For a simple pole `s0`, `Z(s) = residue / (s - s0) + regular(s)`.
/// This is a local Laurent coefficient, not a reduced time-domain model or a
/// guarantee that this pole dominates a measured response peak. Every original
/// acoustic and material-memory coordinate participates in the factorization.
/// Under `exp(-i omega t)`, the angular-frequency residue is `i*residue`.
#[derive(Clone, Debug, PartialEq)]
pub struct DgPoleResidue {
    pub laplace_pole_per_s: C64,
    /// Original boundary-then-side-port order for BOTH matrix axes.
    pub ports: Vec<PortId>,
    /// Row-major pressure-observer by inflow-driver residues [Pa/m3].
    /// Multiplication by `1/(s-s0)` gives pressure / volume flow [Pa s/m3].
    /// Zero or small entries are preserved; no observability threshold is used.
    pub pressure_per_volume_residue_pa_m3: Vec<C64>,
    /// Independent infinity-norm backward residual of the original DG mode.
    pub original_relative_residual: f64,
    /// Shared owner's coefficient-scaled residual of the supplied eigenpair.
    pub pole_relative_residual: f64,
    /// Relative Newton correction checked at the supplied pole, which is not moved.
    pub pole_relative_correction: f64,
    /// Original bordered-row backward residual for each original forcing port.
    /// These are not condition estimates or transfer forward-error bounds.
    pub forcing_relative_residuals: Vec<f64>,
}

impl<M: DgMedium> PreparedDgPoles<'_, M> {
    /// Resolve excitation and observation of a simple free mode at every port.
    ///
    /// One shared `PreparedPole` factorization supplies true left/right forcing
    /// participation for all original volume-flow columns. The original pencil
    /// and its forcing are scaled identically; pressure observations retain
    /// their physical scale. Thus `R_ij=(B_i^T x)(y^T B_j)/(y^T x)` with
    /// `L x=s0*x` and `L^T y=s0*y`, including nonnormal material damping.
    /// No right-mode pressure square or separate adjoint solver is substituted.
    ///
    /// `options` bounds the complete query, including evaluation matrices,
    /// physical maps and all port outputs. `max_vectors` must admit every port.
    /// Source caches, caller inputs and allocator overhead are outside the
    /// numerical payload. Work is O(n^3 + p*n^2 + p^2), for n original states
    /// and p ports. A foreign or edited mode must satisfy THIS original operator;
    /// its cached residuals and physical observations are never trusted.
    ///
    /// Oscillatory poles obey the original material-frequency admission. Real
    /// relaxation and exactly dark modes remain valid. Repeated/defective poles
    /// can refuse; this API does not regularize a singular simple-mode border.
    /// Cancellation or refusal publishes no partial result and changes no state.
    pub fn port_residue(
        &self,
        pole: &DgResonancePole,
        options: &PoleResponseOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DgPoleResidue, DgPoleError> {
        poll(&mut cancelled)?;
        let numerical_options = admit(self.dimension(), self.view.ports().len(), options)?;
        if pole.energy_coordinates.len() != self.dimension()
            || pole.energy_coordinates.iter().any(|&z| !finite(z))
            || !finite(pole.laplace_pole_per_s)
        {
            return Err(invalid("complete finite original DG pole coordinates required").into());
        }
        if pole.laplace_pole_per_s.im != 0.0 {
            self.view
                .admit_frequency(pole.laplace_pole_per_s.im.abs())?;
        }
        let (numerical, residual_scale) = {
            let evaluated = self.evaluate(
                pole.laplace_pole_per_s,
                options.spectral_scale,
                &mut cancelled,
            )?;
            let numerical = PreparedPole::new(
                pole.laplace_pole_per_s,
                &pole.energy_coordinates,
                &evaluated,
                numerical_options,
                &mut cancelled,
            )?;
            (numerical, evaluated.residual_scale)
        }; // Release the two evaluation matrices before mapping all ports.
        // Reuse the native original-row/physical observer. Its temporary state
        // has the shared pivot-one gauge; it is not published as a unit mode.
        let physical = self.observe(
            NonlinearPole {
                value: numerical.value(),
                vector: numerical.right_vector().to_vec(),
                residual: numerical.pole_residual(),
                relative_step: numerical.relative_correction(),
                residual_scale,
                iterations: 0,
            },
            options.pole_residual_tolerance,
            &mut cancelled,
        )?;
        let (residues, forcing_residuals) =
            observe_ports(&numerical, self.view.ports(), self.scale, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(DgPoleResidue {
            laplace_pole_per_s: numerical.value(),
            ports: self.view.ports().iter().map(|(id, _)| *id).collect(),
            pressure_per_volume_residue_pa_m3: residues,
            original_relative_residual: physical.relative_residual,
            pole_relative_residual: numerical.pole_residual(),
            pole_relative_correction: numerical.relative_correction(),
            forcing_relative_residuals: forcing_residuals,
        })
    }
}

fn admit(
    dimension: usize,
    ports: usize,
    options: &PoleResponseOptions,
) -> Result<PoleResponseOptions, DgPoleError> {
    if dimension == 0
        || dimension > MAX_NONLINEAR_POLE_DIMENSION
        || dimension > options.max_dimension
        || ports == 0
        || ports > options.max_vectors
        || !(1..=256).contains(&options.max_vectors)
        || !options.spectral_scale.is_finite()
        || options.spectral_scale <= 0.0
        || [
            options.pole_residual_tolerance,
            options.pole_step_tolerance,
            options.linear_residual_tolerance,
        ]
        .iter()
        .any(|tolerance| !tolerance.is_finite() || *tolerance <= 0.0 || *tolerance >= 1.0)
    {
        return Err(
            invalid("complete DG state and all original ports must fit response limits").into(),
        );
    }
    let (n, p) = (dimension as u128, ports as u128);
    // Original pencil/derivative, physical observer, scaled forcing columns,
    // complete port matrix and diagnostics coexist with the shared factor.
    // Shared PreparedPole separately charges every numerical batch vector.
    let bytes = 32 * n * n + 16 * n * p + 32 * p * p + 128 * (n + p + 1);
    let work = 64 * n * n + 64 * n * p + 16 * p * p;
    if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
        return Err(invalid("DG pole residue maps exceed their payload/work capacity").into());
    }
    let mut numerical = *options;
    numerical.max_dense_bytes -= bytes as usize;
    numerical.max_work -= work as u64;
    // Admit the existing owner's full envelope BEFORE allocating our matrices.
    // Its own admission remains authoritative when constructing the factor.
    let d = n + 1;
    let k = options.max_vectors as u128;
    if 128 * d * d + 64 * d * (k + 4) > numerical.max_dense_bytes as u128
        || 32 * d * d * d + 128 * (k + 1) * d * d > numerical.max_work as u128
    {
        return Err(invalid("complete DG pole response exceeds its payload/work capacity").into());
    }
    Ok(numerical)
}

fn observe_ports(
    numerical: &PreparedPole,
    ports: &[(PortId, Vec<f64>)],
    pencil_scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<f64>), DgPoleError> {
    let mut forces = Vec::with_capacity(ports.len());
    let mut observations = Vec::with_capacity(ports.len());
    for (_, b) in ports {
        poll(cancelled)?;
        let mut force = Vec::with_capacity(b.len());
        for &value in b {
            let scaled = value / pencil_scale;
            if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                return Err(
                    invalid("DG pole forcing scale cannot preserve the original port").into(),
                );
            }
            force.push(C64::from_re(scaled));
        }
        forces.push(force);
        let (mut observation, mut correction) = (C64::ZERO, C64::ZERO);
        for (&b, &x) in b.iter().zip(numerical.right_vector()) {
            let term = x.scale(b) - correction;
            let next = observation + term;
            correction = (next - observation) - term;
            observation = next;
        }
        if !finite(observation) {
            return Err(invalid("DG pole port observation is not representable").into());
        }
        observations.push(observation);
    }
    let refs: Vec<_> = forces.iter().map(Vec::as_slice).collect();
    let forced = numerical.forcing_residues(&refs, &mut *cancelled)?;
    let mut residues = Vec::with_capacity(ports.len() * ports.len());
    for observation in observations {
        poll(cancelled)?;
        for source in &forced {
            let residue = observation * source.coefficient;
            if !finite(residue) {
                return Err(invalid("DG pole port residue is not representable").into());
            }
            residues.push(residue);
        }
    }
    poll(cancelled)?;
    Ok((residues, forced.iter().map(|f| f.linear_residual).collect()))
}

#[cfg(test)]
mod tests;
