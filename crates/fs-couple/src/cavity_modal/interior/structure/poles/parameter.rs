//! Parameter tangents of the same complete coupled free mode.
//!
//! The caller supplies derivatives of original physical blocks per declared
//! parameter unit. One dB assembles BOTH opposite interface blocks. The shared
//! quadratic owner differentiates the full pencil, physical mass normalization
//! and a fixed original phase coordinate with one reusable bordered factor.
//! No constitutive derivative, mesh motion, mode crossing or finite change is
//! inferred from these local directions.

use super::{
    C64, Cx, Dd, PreparedStructuralInteriorModes, QuadraticModeError, StructuralInterior,
    StructuralInteriorError, StructuralInteriorPole, StructuralInteriorPoleError, poll,
};
use fs_modal::quadratic_modes::{PreparedQuadraticSensitivities, QuadraticParameterDirection};

/// Shared bordered-factor and native physical-mapping admission.
/// This adapter reserves its additional assembly/output resources from the
/// physical ceilings before passing the remainder to the shared owner.
pub use fs_modal::quadratic_modes::{
    PoleResponseOptions, QuadraticResponseOptions as PoleSensitivityOptions,
};

/// One original physical direction per caller-declared unit of theta.
/// All blocks are complete row-major real arrays. Structural blocks have
/// ns² entries, acoustic/wall blocks na², and coupling na*ns. Zero blocks must
/// be explicit. Matrix directions are signed, not new positive material laws.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorPoleDirection<'d> {
    /// dM_structure/dtheta, symmetric but not necessarily positive.
    pub structural_mass: &'d [f64],
    /// dC_structure/dtheta; all viscous and gyroscopic entries are retained.
    pub structural_damping: &'d [f64],
    /// dK_structure/dtheta, symmetric in the original displacement basis.
    pub structural_stiffness: &'d [f64],
    /// dM_acoustic/dtheta in the original pressure-potential basis.
    pub acoustic_mass: &'d [f64],
    /// dK_acoustic/dtheta in the original pressure-potential basis.
    pub acoustic_stiffness: &'d [f64],
    /// Derivative of the symmetric original relative-wall admittance matrix.
    pub wall_damping: &'d [f64],
    /// dB/dtheta, pressure rows by structural columns. Geometry directions
    /// must include all trace, normal and area-measure derivatives explicitly.
    pub coupling: &'d [f64],
}

/// Derivative of the physical mode and interface fields per unit theta.
/// The same original mass normalization and fixed phase apply to every batch.
#[derive(Clone, Debug, PartialEq)]
pub struct StructuralInteriorPoleDerivative {
    /// ds/dtheta, including the signed growth/decay and oscillation changes.
    /// This remains defined for a simple real pole, where d|Im(s)| may not be.
    pub laplace_pole_derivative_per_s: C64,
    /// Derivative of signed amplitude decay, -Re(ds/dtheta).
    pub amplitude_decay_rate_derivative_per_s: f64,
    /// Original displacement derivative, including differentiated modal mass.
    pub structural_displacement_derivative: Vec<C64>,
    /// `ds*u + s*du` in the original structural coordinates.
    pub structural_velocity_derivative: Vec<C64>,
    /// Derivative of every original pressure-potential coefficient.
    pub pressure_potential_derivative_pa_s: Vec<C64>,
    /// `ds*psi + s*dpsi`, retaining both oscillation and decay contributions.
    pub pressure_derivative_pa: Vec<C64>,
    /// `dB^T*p + B^T*dp` in original structural force coordinates.
    pub fluid_force_derivative: Vec<C64>,
    /// `dB*v + B*dv` in original pressure-trace coordinates.
    pub displaced_flow_derivative_m3_s: Vec<C64>,
    /// Shared bordered-system backward residual.
    pub linear_residual: f64,
    /// Independently checked differentiated original-polynomial row residual.
    pub polynomial_relative_residual: f64,
    /// Residual of `2 Re(phi^H M dphi) + phi^H dM phi = 0`.
    pub normalization_relative_residual: f64,
}

/// One physical mode and one shared factor reused for all parameter directions.
pub struct PreparedStructuralInteriorSensitivities<'p, 'm, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    sensitivities: PreparedQuadraticSensitivities<'p, 'm>,
    mode: StructuralInteriorPole,
    maximum_directions: usize,
}

impl<'m, 'a, 'f, 'c> PreparedStructuralInteriorModes<'m, 'a, 'f, 'c> {
    /// Refine the native mode ONCE, then prepare its shared parameter factor.
    /// The original PoleOptions continue to bound the mode refinement.
    ///
    /// In addition to native sensitivity resources, `options` reserves
    /// `24*k*n*n + 128*(k+1)*n + 512*(k+1)` physical bytes and work score
    /// `128*(k+1)*n*n + 256*(k+1)*n`, with k=maximum directions and n=ns+na.
    /// These cover full direction assembly, the retained physical baseline and
    /// each batch's original output recovery. Caller directions, the original
    /// model/source preparation and earlier returned batches are excluded.
    /// Shared numerical-factor limits remain separate and unchanged.
    ///
    /// Only a locally isolated simple mode is differentiated. Repeated roots,
    /// unresolved crossings and budget/cancellation failures return no prepared
    /// result. No numerical coordinate/rate scaling is differentiated.
    pub fn refine_with_sensitivities<'p>(
        &'p self,
        laplace_seed_per_s: C64,
        structural_displacement_seed: &[C64],
        pressure_potential_seed_pa_s: &[C64],
        options: PoleSensitivityOptions,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedStructuralInteriorSensitivities<'p, 'm, 'a, 'f, 'c>,
        StructuralInteriorPoleError,
    > {
        poll(cx)?;
        let n = self.model.structural_dofs + self.model.pressure_dofs();
        let maximum_directions = options.numerical.max_vectors;
        let native = reserve(n, options)?;
        let mut mode = self.refine_original(
            laplace_seed_per_s,
            structural_displacement_seed,
            pressure_potential_seed_pa_s,
            cx,
        )?;
        let sensitivities = self
            .modes
            .prepare_sensitivities(&mode, native, || cx.checkpoint().is_err())?;
        // The shared factor rechecks and canonically normalizes the genuine
        // native eigenpair. Observe ITS published shape, so all product rules
        // use exactly the same baseline as the parameter solver. No synthetic
        // mode or receipt is supplied to either numerical owner.
        let response = sensitivities.response();
        mode.laplace_frequency_per_s = response.laplace_pole_per_s();
        mode.displacement
            .copy_from_slice(response.right_displacement());
        mode.velocity.copy_from_slice(response.right_velocity());
        mode.relative_residual = response.physical_residual();
        mode.relative_correction = response.relative_correction();
        let mode = self.observe(mode, cx)?;
        poll(cx)?;
        Ok(PreparedStructuralInteriorSensitivities {
            model: self.model,
            sensitivities,
            mode,
            maximum_directions,
        })
    }
}

impl PreparedStructuralInteriorSensitivities<'_, '_, '_, '_, '_> {
    /// Physical baseline from the same rechecked mode as the derivative factor.
    #[must_use]
    pub fn mode(&self) -> &StructuralInteriorPole {
        &self.mode
    }

    /// Fixed phase component in concatenated [structural displacement, potential]
    /// order. This component stays real; a changing maximum is not differentiated.
    #[must_use]
    pub fn phase_reference_coordinate(&self) -> usize {
        self.sensitivities.phase_reference_coordinate()
    }

    /// Differentiate a complete nonempty batch using the SAME original factor.
    /// No eigensolve, factorization or finite-difference query occurs here.
    /// One dB contributes to both the equations and physical force/flow output.
    /// A refused direction cancels the whole batch without changing its baseline.
    pub fn parameter_derivatives(
        &self,
        directions: &[StructuralInteriorPoleDirection<'_>],
        cx: &Cx<'_>,
    ) -> Result<Vec<StructuralInteriorPoleDerivative>, StructuralInteriorPoleError> {
        poll(cx)?;
        if directions.is_empty() || directions.len() > self.maximum_directions {
            return Err(invalid(
                "nonempty original parameter batch within the declared direction cap required",
            ));
        }
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        let mut blocks = Vec::with_capacity(directions.len());
        for direction in directions {
            blocks.push(direction.assemble(ns, na, cx)?);
        }
        let native: Vec<_> = blocks
            .iter()
            .map(|b| QuadraticParameterDirection {
                mass: &b[0],
                damping: &b[1],
                stiffness: &b[2],
            })
            .collect();
        let derivatives = self
            .sensitivities
            .parameter_derivatives(&native, || cx.checkpoint().is_err())?;
        let mut results = Vec::with_capacity(directions.len());
        for (mut derivative, direction) in derivatives.into_iter().zip(directions) {
            poll(cx)?;
            let potential = derivative.displacement.split_off(ns);
            let pressure = derivative.velocity.split_off(ns);
            let (force, flow) =
                self.interface_derivative(direction.coupling, &derivative.velocity, &pressure, cx)?;
            results.push(StructuralInteriorPoleDerivative {
                laplace_pole_derivative_per_s: derivative.laplace_frequency_per_s,
                amplitude_decay_rate_derivative_per_s: -derivative.laplace_frequency_per_s.re,
                structural_displacement_derivative: derivative.displacement,
                structural_velocity_derivative: derivative.velocity,
                pressure_potential_derivative_pa_s: potential,
                pressure_derivative_pa: pressure,
                fluid_force_derivative: force,
                displaced_flow_derivative_m3_s: flow,
                linear_residual: derivative.linear_residual,
                polynomial_relative_residual: derivative.polynomial_relative_residual,
                normalization_relative_residual: derivative.normalization_relative_residual,
            });
        }
        poll(cx)?;
        Ok(results)
    }

    // Differentiate both reciprocal traces in one pass. Keep both product-rule
    // terms in the compensated accumulation before rounding their sum.
    fn interface_derivative(
        &self,
        db: &[f64],
        dv: &[C64],
        dp: &[C64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<C64>, Vec<C64>), StructuralInteriorPoleError> {
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        let mut force_re = vec![Dd::from_f64(0.0); ns];
        let mut force_im = vec![Dd::from_f64(0.0); ns];
        let mut flow = Vec::with_capacity(na);
        for a in 0..na {
            poll(cx)?;
            let p = self.mode.pressure_pa[a];
            let (mut re, mut im) = (Dd::from_f64(0.0), Dd::from_f64(0.0));
            for s in 0..ns {
                let b = Dd::from_f64(self.model.coupling[a * ns + s]);
                let db = Dd::from_f64(db[a * ns + s]);
                let v = self.mode.structural_velocity[s];
                force_re[s] = force_re[s] + b * Dd::from_f64(dp[a].re) + db * Dd::from_f64(p.re);
                force_im[s] = force_im[s] + b * Dd::from_f64(dp[a].im) + db * Dd::from_f64(p.im);
                re = re + b * Dd::from_f64(dv[s].re) + db * Dd::from_f64(v.re);
                im = im + b * Dd::from_f64(dv[s].im) + db * Dd::from_f64(v.im);
            }
            flow.push(C64::new(re.to_f64(), im.to_f64()));
        }
        let force: Vec<_> = force_re
            .into_iter()
            .zip(force_im)
            .map(|(re, im)| C64::new(re.to_f64(), im.to_f64()))
            .collect();
        if force
            .iter()
            .chain(&flow)
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid(
                "original interface force or flow derivative is unrepresentable",
            ));
        }
        Ok((force, flow))
    }
}

impl StructuralInteriorPoleDirection<'_> {
    fn assemble(
        &self,
        ns: usize,
        na: usize,
        cx: &Cx<'_>,
    ) -> Result<[Vec<f64>; 3], StructuralInteriorPoleError> {
        for (matrix, size, symmetric) in [
            (self.structural_mass, ns, true),
            (self.structural_damping, ns, false),
            (self.structural_stiffness, ns, true),
            (self.acoustic_mass, na, true),
            (self.acoustic_stiffness, na, true),
            (self.wall_damping, na, true),
        ] {
            if matrix.len() != size * size {
                return Err(invalid(
                    "parameter direction must cover every original block coordinate",
                ));
            }
            for i in 0..size {
                poll(cx)?;
                if (0..size).any(|j| {
                    !matrix[i * size + j].is_finite()
                        || symmetric && matrix[i * size + j] != matrix[j * size + i]
                }) {
                    return Err(invalid(
                        "finite directions with symmetric mass, stiffness and wall blocks required",
                    ));
                }
            }
        }
        if self.coupling.len() != na * ns {
            return Err(invalid(
                "one complete original dB pressure-by-structure direction required",
            ));
        }
        let n = ns + na;
        let mut blocks = [vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n]];
        for (order, offset, size, source) in [
            (0, 0, ns, self.structural_mass),
            (1, 0, ns, self.structural_damping),
            (2, 0, ns, self.structural_stiffness),
            (0, ns, na, self.acoustic_mass),
            (1, ns, na, self.wall_damping),
            (2, ns, na, self.acoustic_stiffness),
        ] {
            for i in 0..size {
                poll(cx)?;
                blocks[order][(offset + i) * n + offset..(offset + i) * n + offset + size]
                    .copy_from_slice(&source[i * size..(i + 1) * size]);
            }
        }
        for a in 0..na {
            poll(cx)?;
            for s in 0..ns {
                let value = self.coupling[a * ns + s];
                if !value.is_finite() {
                    return Err(invalid("finite original interface derivative required"));
                }
                blocks[1][(ns + a) * n + s] = value;
                blocks[1][s * n + ns + a] = -value;
            }
        }
        Ok(blocks)
    }
}

fn invalid(what: &'static str) -> StructuralInteriorPoleError {
    StructuralInteriorError::Invalid(what).into()
}

fn reserve(
    n: usize,
    mut options: PoleSensitivityOptions,
) -> Result<PoleSensitivityOptions, StructuralInteriorPoleError> {
    let (n, k) = (n as u128, options.numerical.max_vectors as u128);
    if k == 0 {
        return Err(invalid(
            "positive original parameter direction capacity required",
        ));
    }
    let bytes = 24 * k * n * n + 128 * (k + 1) * n + 512 * (k + 1);
    let work = 128 * (k + 1) * n * n + 256 * (k + 1) * n;
    for (quantity, requested, cap) in [
        (
            "coupled parameter assembly/output bytes",
            bytes,
            options.max_physical_bytes as u128,
        ),
        (
            "coupled parameter assembly/output work",
            work,
            options.max_physical_work as u128,
        ),
    ] {
        if requested > cap {
            return Err(QuadraticModeError::Budget {
                quantity,
                requested,
                cap,
            }
            .into());
        }
    }
    options.max_physical_bytes -= bytes as usize;
    options.max_physical_work -= work as u64;
    Ok(options)
}
