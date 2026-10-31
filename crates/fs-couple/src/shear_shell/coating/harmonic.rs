//! In-band forced response of source-backed films on an existing shell image.
//!
//! Peak phasors use exp(-i omega t), matching the original coupled harmonic
//! owner. For the equilibrium body basis this solves
//! `Z = diag(omega_mode²) - omega² I - i omega C_background`
//! `    + sum_f [(E'_f(omega) - E0_f) - i E''_f(omega)] H_f`.
//! Every complete film strain form H comes from the original MITC projector
//! with its equilibrium bubble transformation. Off-diagonal terms and signed
//! membrane/bending coupling are retained. The existing complex LU owner solves
//! Z q = F. No Prony fit, new eigensolve or extrapolation of the measured
//! constitutive law enters this fixed-basis response.
//!
//! Unselected wood/glue damping remains the fixed reduced operator selected
//! when the image was made. Explicit substrate memory adds its same causal
//! storage and loss at the drive frequency. The film law is evaluated at the common drive
//! frequency. This is an unloaded structural response: signed volume velocity
//! is exposed, but no cavity/radiation impedance, bridge/string feedback, sound
//! pressure or perceptual-quality claim is inferred from it.

use super::relaxation::ShellCoatingRelaxation;
use crate::acoustic_realize::AcousticRealizeError;
use crate::shear_shell::{
    MaterialShellImage, ResolvedShearShell, assembly::MaterialShellAssemblyImage, invalid, owner,
    poll,
};
use crate::thin_plate::CompactBody;
use fs_exec::CancelGate;
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;

/// Bounded dense work and a backward-error gate on the original equations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellCoatingHarmonicOptions {
    /// Admits n³+n²(2 films+interfaces+6)+n(2 ports+interfaces+3) scalar terms.
    /// No retained mode, film, interface or observation is silently discarded.
    pub maximum_work: usize,
    /// Infinity-norm backward error relative to `||Z|| ||q|| + ||F||`, in (0,1).
    pub solve_relative_tolerance: f64,
}

/// AC structural response in the unchanged mass-normalized body basis.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellCoatingHarmonicResponse {
    pub omega_rad_s: f64,
    /// Peak modal displacements [m sqrt(kg)], in original retained mode order.
    pub displacement_m_sqrt_kg: Vec<C64>,
    /// Peak modal velocities, exactly `-i omega q` [m sqrt(kg)/s].
    pub velocity_m_sqrt_kg_per_s: Vec<C64>,
    /// Work-conjugate physical port displacement [m], in original port order.
    pub force_port_displacement_m: Vec<C64>,
    /// Work-conjugate physical port velocity [m/s], in original port order.
    pub force_port_velocity_m_s: Vec<C64>,
    /// Signed exterior normal-volume velocity [m³/s], with the original skin.
    pub exterior_volume_velocity_m3_s: C64,
    /// Independently signed outward-from-cavity volume velocity [m³/s].
    pub cavity_volume_velocity_m3_s: C64,
    /// Mean supplied mechanical power, `Re(F* . v)/2` [W].
    pub external_power_w: f64,
    /// Fixed image wood/glue viscosity, excluding film and separate interfaces.
    pub background_dissipation_w: f64,
    /// Positive heat of explicitly selected substrate Maxwell histories [W].
    pub substrate_dissipation_w: f64,
    /// Original material-interface dashpots; elastic terms are already in K0.
    pub interface_dissipation_w: f64,
    /// `omega E''_f q* H_f q/2` [W], in projected source-film order.
    pub film_dissipation_w: Vec<f64>,
    /// Sum of physical mean losses minus external mean power [W].
    pub power_residual_w: f64,
    /// Backward residual of the original unscaled equations.
    pub solve_relative_residual: f64,
}

impl ShellCoatingHarmonicResponse {
    pub fn total_dissipation_w(&self) -> f64 {
        self.background_dissipation_w
            + self.substrate_dissipation_w
            + self.interface_dissipation_w
            + self.film_dissipation_w.iter().sum::<f64>()
    }
}

struct Image<'a> {
    modes: &'a [CompactBody],
    background: &'a [f64],
    force_ports: &'a [Vec<f64>],
    exterior: &'a [f64],
    cavity: &'a [f64],
    films: Option<&'a ShellCoatingRelaxation>,
    substrate: Option<&'a crate::shear_shell::constant_loss::relaxation::ShellSubstrateRelaxation>,
    interfaces: Vec<(&'a [f64], f64)>,
}

fn source_support(
    source: &ResolvedShearShell,
    omega: f64,
    gate: &CancelGate,
) -> Result<(), AcousticRealizeError> {
    poll(gate)?;
    if !omega.is_finite() || omega <= 0.0 {
        return Err(invalid(
            "coating harmonic frequency must be finite and positive",
        ));
    }
    if !source.coatings().is_empty() && source.coating_relaxation_options().is_none() {
        return Err(invalid(
            "coating harmonic response cannot mix frozen and hereditary films",
        ));
    }
    for band in source
        .regions()
        .iter()
        .map(|region| region.omega_band_rad_s())
        .chain(
            source
                .interfaces()
                .iter()
                .map(|link| link.omega_band_rad_s()),
        )
    {
        if omega < band.0 || omega > band.1 {
            return Err(invalid(
                "coating harmonic drive exceeds original substrate or interface support",
            ));
        }
    }
    Ok(())
}

impl MaterialShellImage {
    /// Evaluate each source film at this drive frequency and solve the original
    /// equilibrium-basis shell image. `force` contains one peak generalized
    /// force [N/sqrt(kg)] per mode. A physical force at port p is its unchanged
    /// `force_ports()[p]` row times that force's complex amplitude [N].
    ///
    /// A coated image must be reduced after `with_coating_relaxation`:
    /// equilibrium film stiffness occurs once and no frozen film dashpot is
    /// also present. An entirely uncoated image supplies the direct bare-body
    /// comparison. The selected Prony continuation is NOT used in this in-band
    /// calculation. Material histories and accepted states are never touched.
    pub fn coating_harmonic_response(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        options: ShellCoatingHarmonicOptions,
        gate: &CancelGate,
    ) -> Result<ShellCoatingHarmonicResponse, AcousticRealizeError> {
        source_support(self.source(), omega_rad_s, gate)?;
        let films = self.coating_relaxation();
        if films.map_or(0, |memory| memory.films().len()) != self.source().coatings().len() {
            return Err(invalid(
                "coating harmonic response requires every source film projection",
            ));
        }
        solve_image(
            Image {
                modes: self.radiators(),
                background: self.damping().matrix_per_s(),
                force_ports: self.force_ports(),
                exterior: self.volume_shapes_m2_per_sqrt_kg(),
                cavity: self.cavity_volume_shapes_m2_per_sqrt_kg(),
                films,
                substrate: self.substrate_relaxation(),
                interfaces: self
                    .source()
                    .interfaces()
                    .iter()
                    .zip(self.interface_columns())
                    .filter(|(source, _)| source.law().damping_n_s_m != 0.0)
                    .map(|(source, row)| (row.as_slice(), source.law().damping_n_s_m))
                    .collect(),
            },
            omega_rad_s,
            force,
            options,
            gate,
        )
    }
}

impl MaterialShellAssemblyImage {
    /// The same source-law response on one monolithic tied-shell body basis.
    /// Glue viscosity is already included in the image background; component
    /// interfaces are added once. Exterior/cavity rows preserve their original
    /// component orientations. No additional boundary or fluid model is made.
    pub fn coating_harmonic_response(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        options: ShellCoatingHarmonicOptions,
        gate: &CancelGate,
    ) -> Result<ShellCoatingHarmonicResponse, AcousticRealizeError> {
        poll(gate)?;
        for component in &self.source().components {
            source_support(&component.specimen, omega_rad_s, gate)?;
        }
        let films = self.coating_relaxation();
        if films.map_or(0, |memory| memory.films().len())
            != self
                .source()
                .components
                .iter()
                .map(|part| part.specimen.coatings().len())
                .sum::<usize>()
        {
            return Err(invalid(
                "coating harmonic response requires every source film projection",
            ));
        }
        solve_image(
            Image {
                modes: self.radiators(),
                background: self.damping().matrix_per_s(),
                force_ports: self.force_ports(),
                exterior: self.volume_shapes_m2_per_sqrt_kg(),
                cavity: self.cavity_volume_shapes_m2_per_sqrt_kg(),
                films,
                substrate: self.substrate_relaxation(),
                interfaces: self
                    .source()
                    .components
                    .iter()
                    .zip(self.interface_columns())
                    .flat_map(|(component, columns)| {
                        component.specimen.interfaces().iter().zip(columns)
                    })
                    .filter(|(source, _)| source.law().damping_n_s_m != 0.0)
                    .map(|(source, row)| (row.as_slice(), source.law().damping_n_s_m))
                    .collect(),
            },
            omega_rad_s,
            force,
            options,
            gate,
        )
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn projection(row: &[f64], state: &[C64]) -> C64 {
    row.iter()
        .zip(state)
        .fold(C64::ZERO, |sum, (&a, &q)| sum + q.scale(a))
}

fn quadratic(
    matrix: &[f64],
    state: &[C64],
    gate: &CancelGate,
) -> Result<f64, AcousticRealizeError> {
    let n = state.len();
    let mut sum = 0.0;
    let mut absolute = 0.0;
    for (i, &qi) in state.iter().enumerate() {
        poll(gate)?;
        for (j, &qj) in state.iter().enumerate() {
            let value = matrix[i * n + j] * (qi.re * qj.re + qi.im * qj.im);
            sum += value;
            absolute += value.abs();
        }
    }
    if !sum.is_finite() || !absolute.is_finite() || sum < -128.0 * f64::EPSILON * absolute {
        return Err(invalid(
            "coating harmonic quadratic loss is nonfinite or nonpassive",
        ));
    }
    Ok(sum.max(0.0))
}

fn lu_error(error: LuComplexError) -> AcousticRealizeError {
    match error {
        LuComplexError::Cancelled => AcousticRealizeError::Cancelled,
        LuComplexError::Singular { .. } => {
            invalid("coating harmonic dynamic stiffness is singular")
        }
    }
}

fn solve_image(
    image: Image<'_>,
    omega: f64,
    force: &[C64],
    options: ShellCoatingHarmonicOptions,
    gate: &CancelGate,
) -> Result<ShellCoatingHarmonicResponse, AcousticRealizeError> {
    poll(gate)?;
    let n = image.modes.len();
    let films = image.films.map_or(&[][..], ShellCoatingRelaxation::films);
    let work = n.checked_mul(n).and_then(|square| {
        let per_square = films
            .len()
            .checked_mul(2)?
            .checked_add(image.interfaces.len())?
            .checked_add(if image.substrate.is_some() { 8 } else { 6 })?
            .checked_add(n)?;
        let observations = image
            .force_ports
            .len()
            .checked_mul(2)?
            .checked_add(image.interfaces.len())?
            .checked_add(3)?;
        square
            .checked_mul(per_square)?
            .checked_add(n.checked_mul(observations)?)
    });
    if n == 0
        || n > crate::shear_shell::MAX_SHELL_MODES
        || image.films.is_some_and(|memory| memory.mode_count() != n)
        || image.substrate.is_some_and(|memory| memory.mode_count() != n
            || memory.loss_stiffness_per_s2().len() != n*n)
        || force.len() != n
        || force.iter().any(|&f| !finite(f))
        || !omega.is_finite()
        || omega <= 0.0
        || !options.solve_relative_tolerance.is_finite()
        || options.solve_relative_tolerance <= 0.0
        || options.solve_relative_tolerance >= 1.0
        || work.is_none_or(|w| w > options.maximum_work)
        || image.background.len() != n * n
        || image.exterior.len() != n
        || image.cavity.len() != n
        || image.force_ports.iter().any(|row| row.len() != n)
        || image
            .interfaces
            .iter()
            .any(|(row, damping)| row.len() != n || !damping.is_finite() || *damping < 0.0)
    {
        return Err(invalid(
            "coating harmonic response needs complete finite forces, bounded work and a solve tolerance",
        ));
    }
    if let Some(memory) = image.films {
        let band = memory.source_band_rad_s();
        if omega < band.0 || omega > band.1 {
            return Err(invalid(
                "coating harmonic drive exceeds the original measured film support",
            ));
        }
    }
    let mut operator: Vec<_> = image
        .background
        .iter()
        .map(|&c| C64::new(0.0, -omega * c))
        .collect();
    for (i, mode) in image.modes.iter().enumerate() {
        operator[i * n + i].re = mode.omega * mode.omega - omega * omega;
    }
    for &(row, damping) in &image.interfaces {
        for (i, &a) in row.iter().enumerate() {
            poll(gate)?;
            for (j, &b) in row.iter().enumerate() {
                operator[i * n + j].im -= omega * damping * a * b;
            }
        }
    }
    let substrate_loss = if let Some(memory) = image.substrate {
        let band = memory.source_band_rad_s();
        if omega < band.0 || omega > band.1 {
            return Err(invalid("harmonic substrate response exceeds original material/fit support"));
        }
        let (storage, loss) = memory.modulus_shapes(omega);
        let increment = C64::new(storage, -loss);
        for (value, &h) in operator.iter_mut().zip(memory.loss_stiffness_per_s2()) {
            *value = *value + increment.scale(h);
        }
        loss
    } else { 0.0 };
    let mut losses = Vec::with_capacity(films.len());
    for film in films {
        poll(gate)?;
        let source = film.source();
        let law = source.input().material.law();
        let band = source.omega_band_rad_s();
        let kernel = film.unit_stiffness_per_pa_s2();
        if omega < band.0
            || omega > band.1
            || source.reference_storage_pa() != law.equilibrium_modulus_pa()
            || kernel.len() != n * n
        {
            return Err(invalid(
                "coating harmonic film needs its exact equilibrium basis and source-supported frequency",
            ));
        }
        let (storage, loss) = law.modulus(omega).map_err(owner)?;
        let increment = C64::new(storage - source.reference_storage_pa(), -loss);
        for i in 0..n {
            poll(gate)?;
            for j in 0..n {
                operator[i * n + j] = operator[i * n + j] + increment.scale(kernel[i * n + j]);
            }
        }
        losses.push(loss);
    }
    if operator.iter().any(|&z| !finite(z)) {
        return Err(invalid("coating harmonic operator overflow"));
    }
    let scale = operator.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale == 0.0 {
        return Err(invalid("coating harmonic operator has no finite scale"));
    }
    let scaled: Vec<_> = operator
        .iter()
        .map(|z| C64::new(z.re / scale, z.im / scale))
        .collect();
    let factor = lu_complex_with_cancel(&scaled, n, || gate.is_requested()).map_err(lu_error)?;
    let mut q: Vec<_> = force
        .iter()
        .map(|f| C64::new(f.re / scale, f.im / scale))
        .collect();
    if q.iter().any(|&z| !finite(z)) {
        return Err(invalid("coating harmonic force scaling overflow"));
    }
    factor
        .solve_with_cancel(&mut q, || gate.is_requested())
        .map_err(lu_error)?;
    if q.iter().any(|&z| !finite(z)) {
        return Err(invalid("coating harmonic response is not representable"));
    }
    let q_norm = q.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let f_norm = force.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let mut residual = 0.0_f64;
    let mut z_norm = 0.0_f64;
    for i in 0..n {
        poll(gate)?;
        let mut row_residual = -force[i];
        let mut row_norm = 0.0;
        for (j, &value) in q.iter().enumerate() {
            row_residual = row_residual + operator[i * n + j] * value;
            row_norm += operator[i * n + j].abs();
        }
        if !finite(row_residual) || !row_norm.is_finite() {
            return Err(invalid("coating harmonic residual overflow"));
        }
        residual = residual.max(row_residual.abs());
        z_norm = z_norm.max(row_norm);
    }
    let denominator = z_norm * q_norm + f_norm;
    if !denominator.is_finite() {
        return Err(invalid("coating harmonic residual scale overflow"));
    }
    let relative = if denominator == 0.0 {
        0.0
    } else {
        residual / denominator
    };
    if relative > options.solve_relative_tolerance {
        return Err(owner(format!(
            "coating harmonic solve residual {relative:e} exceeds {}",
            options.solve_relative_tolerance
        )));
    }
    let v: Vec<_> = q
        .iter()
        .map(|q| C64::new(omega * q.im, -omega * q.re))
        .collect();
    let external_power_w = 0.5
        * force
            .iter()
            .zip(&v)
            .map(|(&f, &v)| (f.conj() * v).re)
            .sum::<f64>();
    let background_dissipation_w = 0.5 * quadratic(image.background, &v, gate)?;
    let substrate_dissipation_w = if let Some(memory) = image.substrate {
        0.5 * omega * substrate_loss * quadratic(memory.loss_stiffness_per_s2(), &q, gate)?
    } else { 0.0 };
    let mut interface_dissipation_w = 0.0;
    for &(row, damping) in &image.interfaces {
        poll(gate)?;
        interface_dissipation_w += 0.5 * damping * projection(row, &v).norm_sq();
    }
    let film_dissipation_w: Vec<_> = films
        .iter()
        .zip(losses)
        .map(|(film, loss)| {
            quadratic(film.unit_stiffness_per_pa_s2(), &q, gate)
                .map(|value| 0.5 * omega * loss * value)
        })
        .collect::<Result<_, _>>()?;
    let mut force_port_displacement_m = Vec::with_capacity(image.force_ports.len());
    let mut force_port_velocity_m_s = Vec::with_capacity(image.force_ports.len());
    for row in image.force_ports {
        poll(gate)?;
        force_port_displacement_m.push(projection(row, &q));
        force_port_velocity_m_s.push(projection(row, &v));
    }
    let exterior_volume_velocity_m3_s = projection(image.exterior, &v);
    let cavity_volume_velocity_m3_s = projection(image.cavity, &v);
    let power_residual_w =
        background_dissipation_w + substrate_dissipation_w + interface_dissipation_w + film_dissipation_w.iter().sum::<f64>()
            - external_power_w;
    if [
        external_power_w,
        background_dissipation_w,
        substrate_dissipation_w,
        interface_dissipation_w,
        power_residual_w,
    ]
    .iter()
    .chain(&film_dissipation_w)
    .any(|v| !v.is_finite())
        || force_port_displacement_m
            .iter()
            .chain(&force_port_velocity_m_s)
            .chain(&v)
            .any(|&z| !finite(z))
        || !finite(exterior_volume_velocity_m3_s)
        || !finite(cavity_volume_velocity_m3_s)
    {
        return Err(invalid("coating harmonic observation or power overflow"));
    }
    poll(gate)?;
    Ok(ShellCoatingHarmonicResponse {
        omega_rad_s: omega,
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
        force_port_displacement_m,
        force_port_velocity_m_s,
        exterior_volume_velocity_m3_s,
        cavity_volume_velocity_m3_s,
        external_power_w,
        background_dissipation_w,
        substrate_dissipation_w,
        interface_dissipation_w,
        film_dissipation_w,
        power_residual_w,
        solve_relative_residual: relative,
    })
}
