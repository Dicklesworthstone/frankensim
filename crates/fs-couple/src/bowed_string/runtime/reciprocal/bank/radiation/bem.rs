//! Exterior air reaction of the SAME material shell used by every string.
//!
//! The original global displacement modes are traced onto the explicitly owned
//! exterior triangles. The Helmholtz owner solves all normal-velocity fields in
//! one batch per frequency. Pressure is projected with the transpose of that
//! same trace and physical panel areas: Z_ij = sum_p phi_pi A_p p_pj. Thus
//! opposing faces, modal signs, mutual radiation and the complete gas survive;
//! no equivalent piston, copied damping ratio or fitted timbre is introduced.
//!
//! Scope: a closed, oriented midsurface and piecewise-constant acoustic panels,
//! with the existing near-triangle and seven-point far Burton--Miller quadrature.
//! This is a cold finite-band passive approximation of a BEM calculation, not a
//! continuum certificate or measured-instrument claim. The existing compact
//! microphone still observes loaded acceleration; this module does not turn it
//! into a finite-distance broadband receiver. Open holes need an explicit
//! interior/exterior boundary coupling and are refused, never silently capped.

use super::super::{ReciprocalStringBank, ReciprocalStringError};
use crate::render::plate::impact::radiation::{Model, fit};
use crate::shear_shell::ResolvedShearShell;
use fs_bem::helmholtz::{self, Formulation, HelmholtzError, Medium, PanelRule};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_plate::shell::pressure::ShellFace;
use std::collections::BTreeMap;

/// Independent geometry, frequency, work and causal-realization limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellBemOptions {
    /// Positive lower frequency [Hz].
    pub minimum_hz: f64,
    /// Upper frequency [Hz], covering every retained body reference mode.
    pub maximum_hz: f64,
    /// Odd 33..=257 frequencies: alternating training and withheld rows.
    pub frequency_samples: usize,
    /// All selected exterior triangles must fit; none are discarded.
    pub maximum_surface_panels: usize,
    /// Ceiling on frequencies * (panels^3 + panels^2*modes + panels*modes^2).
    pub maximum_sweep_terms: u64,
    /// Complete acoustic pole count, before the exact owner's stricter limits.
    pub maximum_poles: usize,
}
impl ShellBemOptions {
    /// Admission before geometry allocation or a Helmholtz factorization.
    pub fn validate(&self, rate: u32) -> Result<(), ShellBemError> {
        if !self.minimum_hz.is_finite()
            || self.minimum_hz <= 0.0
            || !self.maximum_hz.is_finite()
            || self.maximum_hz <= self.minimum_hz
            || !(33..=257).contains(&self.frequency_samples)
            || self.frequency_samples % 2 == 0
            || !(4..=2048).contains(&self.maximum_surface_panels)
            || !(1..=100_000_000_000).contains(&self.maximum_sweep_terms)
            || !(1..=256).contains(&self.maximum_poles)
            || rate == 0
            || 8.0 * core::f64::consts::TAU * self.maximum_hz
                >= 0.9 * core::f64::consts::PI * f64::from(rate)
        {
            return Err(ShellBemError::Invalid(
                "shell radiation needs a positive resolved band, odd 33..257 frequencies and bounded panel/work/pole allowances",
            ));
        }
        Ok(())
    }
}

/// Source binding and numerical diagnostics of the accepted physical load.
#[derive(Clone, Debug)]
pub struct ShellBemReport {
    /// Same signed material/geometry/modal basis pinned by the loaded bank.
    pub basis_identity: ContentHash,
    /// Complete atmosphere shared by strings, body, cavity and pressure output.
    pub gas: GasState,
    pub options: ShellBemOptions,
    pub body_modes: usize,
    pub surface_panels: usize,
    pub acoustic_poles: usize,
    pub sweep_terms: u64,
    /// Lowest panel-resolution diagnostic on the complete sampled band.
    pub minimum_panels_per_wavelength: f64,
    /// Largest original BEM probe-based condition LOWER bound, not an enclosure.
    pub maximum_condition_lower_bound: f64,
    /// Withheld whole-matrix complex fit error under the existing fitter's gate.
    pub peak_error: f64,
    pub rms_error: f64,
    /// Independently gated resistance error; weak loss cannot hide in inertia.
    pub resistance_peak_error: f64,
    pub resistance_rms_error: f64,
}

#[derive(Debug)]
pub enum ShellBemError {
    Invalid(&'static str),
    Geometry(fs_bem::BemError),
    Solve(HelmholtzError),
    Fit(String),
    Bank(ReciprocalStringError),
    Cancelled,
}
impl core::fmt::Display for ShellBemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material shell radiation: {self:?}")
    }
}
impl std::error::Error for ShellBemError {}
fn poll(gate: &CancelGate) -> Result<(), ShellBemError> {
    if gate.is_requested() {
        Err(ShellBemError::Cancelled)
    } else {
        Ok(())
    }
}

impl ReciprocalStringBank {
    /// Compute a closed shell's air load and attach it to the existing reciprocal
    /// mechanics. String/body/material changes rebuild the original modes and
    /// their radiation together. The same gas supplies c, rho and string loss.
    ///
    /// Body pressure reacts on the bridges and every string inside the existing
    /// exact propagator, friction solve, work balance and scheduled performance.
    /// No independent body oscillator, gain normalization or new stepper is used.
    pub fn with_material_shell_bem(
        self,
        options: ShellBemOptions,
        gate: &CancelGate,
    ) -> Result<(Self, ShellBemReport), ShellBemError> {
        let (model, report) = self.prepare_material_shell_bem(options, gate)?;
        let bank = self
            .with_body_radiation_load(model, options.maximum_poles, gate)
            .map_err(ShellBemError::Bank)?;
        poll(gate)?;
        Ok((bank, report))
    }

    /// Cold preparation without changing any mechanical state. Every original
    /// body mode participates, including those with zero net displaced volume.
    /// The result belongs to this exact material basis and gas; prefer the
    /// combined constructor above to preserve that binding automatically.
    ///
    /// A closed outward skin is required. Separate assembly parts may meet at
    /// identical coordinates; no tolerance-based welding or hole filling occurs.
    /// Exact edge closure is checked, not self-intersection or manifold quality.
    /// BEM is synchronous within one bounded frequency solve; cancellation is
    /// polled between frequencies, material parts and projection rows.
    pub fn prepare_material_shell_bem(
        &self,
        options: ShellBemOptions,
        gate: &CancelGate,
    ) -> Result<(Model, ShellBemReport), ShellBemError> {
        poll(gate)?;
        options.validate(self.source.options.sample_rate_hz)?;
        if self.system.samples_rendered() != 0
            || self.source.prepared_linear.is_some()
            || self.source.body_radiation.is_some()
        {
            return Err(ShellBemError::Invalid(
                "computed shell radiation requires an unstarted exact bank without a prior exterior load",
            ));
        }
        let (basis_identity, _) = self
            .source
            .body_radiation_basis(&self.system, gate)
            .map_err(ShellBemError::Bank)?;
        let modes = self
            .system
            .components()
            .get(self.source.strings.len())
            .ok_or(ShellBemError::Invalid(
                "missing original material body component",
            ))?
            .modes();
        if modes.is_empty()
            || modes.len() > crate::render::plate::impact::radiation::MAX_PORTS
            || modes.iter().any(|m| {
                let hz = m.angular_frequency_rad_s / core::f64::consts::TAU;
                hz < options.minimum_hz || hz > options.maximum_hz
            })
        {
            return Err(ShellBemError::Invalid(
                "shell radiation must cover every original body mode within the shared passive port capacity",
            ));
        }
        let mut parts = Vec::new();
        if let Some(image) = &self.source.material_shell_image {
            parts.push((image.source(), image.physical_modes()));
        } else if let Some(image) = &self.source.material_shell_assembly_image {
            for (part, rows) in image.source().components.iter().zip(image.physical_modes()) {
                parts.push((&part.specimen, rows.as_slice()));
            }
        } else {
            return Err(ShellBemError::Invalid(
                "computed closed-body radiation requires an original material shell basis",
            ));
        }
        let panels = parts
            .iter()
            .try_fold(0usize, |count, (shell, _)| {
                count.checked_add(
                    shell
                        .pressure_surfaces()
                        .map_or(shell.mesh().tris.len(), |s| s.exterior.len()),
                )
            })
            .ok_or(ShellBemError::Invalid(
                "shell exterior panel extent overflow",
            ))?;
        let p = panels as u64;
        let n = modes.len() as u64;
        let work = p
            .checked_mul(p)
            .and_then(|p2| p.checked_add(n).and_then(|sum| p2.checked_mul(sum)))
            .and_then(|v| {
                n.checked_mul(n)
                    .and_then(|n2| p.checked_mul(n2))
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| v.checked_mul(options.frequency_samples as u64))
            .ok_or(ShellBemError::Invalid("complete BEM sweep extent overflow"))?;
        if panels < 4
            || panels > options.maximum_surface_panels
            || work > options.maximum_sweep_terms
        {
            return Err(ShellBemError::Invalid(
                "complete shell exterior exceeds its panel or sweep-work allowance",
            ));
        }
        let (surface, fields) = surface_trace(&parts, modes.len(), gate)?;
        let gas = self.source.gas;
        gas.validate_coherence()
            .map_err(|_| ShellBemError::Invalid("incoherent shell exterior gas"))?;
        let medium = Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        };
        let low = core::f64::consts::TAU * options.minimum_hz;
        let high = core::f64::consts::TAU * options.maximum_hz;
        let mut omega: Vec<_> = (0..options.frequency_samples)
            .map(|i| low + (high - low) * i as f64 / (options.frequency_samples - 1) as f64)
            .collect();
        omega[0] = low;
        omega[options.frequency_samples - 1] = high;
        let mut impedance = Vec::with_capacity(omega.len());
        let mut ppw = f64::INFINITY;
        let mut condition = 0.0_f64;
        for &w in &omega {
            poll(gate)?;
            let (row, resolution, bound) = modal_impedance(&surface, &fields, w, medium, gate)?;
            impedance.push(row);
            ppw = ppw.min(resolution);
            condition = condition.max(bound);
        }
        let mut extent = 0.0_f64;
        for a in surface.centroids() {
            poll(gate)?;
            for b in surface.centroids() {
                extent = extent.max((a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]));
            }
        }
        let spacing = if extent > 0.0 {
            (core::f64::consts::TAU * gas.sound_speed / (6.0 * extent)).min(high - low)
        } else {
            high - low
        };
        let fitted = fit::fit_resolved(&omega, &impedance, modes.len(), spacing)
            .map_err(ShellBemError::Fit)?;
        poll(gate)?;
        if fitted.model.poles.len() > options.maximum_poles {
            return Err(ShellBemError::Invalid(
                "complete shell acoustic fit exceeds its pole allowance",
            ));
        }
        let report = ShellBemReport {
            basis_identity,
            gas,
            options,
            body_modes: modes.len(),
            surface_panels: panels,
            acoustic_poles: fitted.model.poles.len(),
            sweep_terms: work,
            minimum_panels_per_wavelength: ppw,
            maximum_condition_lower_bound: condition,
            peak_error: fitted.peak_error,
            rms_error: fitted.rms_error,
            resistance_peak_error: fitted.resistance_peak_error,
            resistance_rms_error: fitted.resistance_rms_error,
        };
        Ok((fitted.model, report))
    }
}

// BEM velocity and pressure use the SAME constant-panel normal trace. Rotations
// affect the elastic modes, but have no extra midsurface translation to add.
fn surface_trace(
    parts: &[(&ResolvedShearShell, &[Vec<f64>])],
    modes: usize,
    gate: &CancelGate,
) -> Result<(SpherePanels, Vec<Vec<C64>>), ShellBemError> {
    let mut triangles = Vec::new();
    let mut node_rows = Vec::new();
    for &(shell, rows) in parts {
        poll(gate)?;
        let mesh = shell.mesh();
        if rows.len() != modes || rows.iter().any(|r| r.len() != 6 * mesh.nodes.len()) {
            return Err(ShellBemError::Invalid(
                "exterior trace must retain the complete original shell basis",
            ));
        }
        let facets: Vec<_> = shell.pressure_surfaces().map_or_else(
            || {
                (0..mesh.tris.len())
                    .map(|i| (i, ShellFace::Positive))
                    .collect()
            },
            |s| s.exterior.iter().map(|f| (f.triangle, f.side)).collect(),
        );
        for (index, side) in facets {
            poll(gate)?;
            let mut nodes = mesh.tris[index];
            if side == ShellFace::Negative {
                nodes.swap(1, 2);
            }
            triangles.push(nodes.map(|n| mesh.nodes[n]));
            let mut velocities = Vec::with_capacity(modes);
            for row in rows {
                velocities.push(core::array::from_fn::<_, 3, _>(|axis| {
                    (row[6 * nodes[0] + axis] + row[6 * nodes[1] + axis] + row[6 * nodes[2] + axis])
                        / 3.0
                }));
            }
            node_rows.push(velocities);
        }
    }
    closed_outward(&triangles)?;
    let surface = SpherePanels::from_triangles(triangles).map_err(ShellBemError::Geometry)?;
    let mut fields = vec![Vec::with_capacity(surface.areas().len()); modes];
    for (normal, velocities) in surface.normals().iter().zip(node_rows) {
        poll(gate)?;
        for (field, velocity) in fields.iter_mut().zip(velocities) {
            let value = normal[0] * velocity[0] + normal[1] * velocity[1] + normal[2] * velocity[2];
            if !value.is_finite() {
                return Err(ShellBemError::Invalid(
                    "nonfinite material shell normal trace",
                ));
            }
            field.push(C64::new(value, 0.0));
        }
    }
    Ok((surface, fields))
}

fn modal_impedance(
    surface: &SpherePanels,
    fields: &[Vec<C64>],
    omega: f64,
    medium: Medium,
    gate: &CancelGate,
) -> Result<(Vec<C64>, f64, f64), ShellBemError> {
    let refs: Vec<_> = fields.iter().map(Vec::as_slice).collect();
    let solutions = helmholtz::solve_radiation_batch_with_rule(
        surface,
        omega / medium.sound_speed,
        medium,
        &refs,
        Formulation::BurtonMiller,
        PanelRule::Triangle {
            max_pairs: surface.areas().len().saturating_pow(2),
        },
    )
    .map_err(ShellBemError::Solve)?;
    let n = fields.len();
    let mut matrix = vec![C64::ZERO; n * n];
    let mut ppw = f64::INFINITY;
    let mut condition = 0.0_f64;
    for (j, solution) in solutions.iter().enumerate() {
        poll(gate)?;
        ppw = ppw.min(solution.panels_per_wavelength);
        condition = condition.max(solution.condition_lower_bound);
        for (i, field) in fields.iter().enumerate() {
            let mut value = C64::ZERO;
            for ((shape, &area), &pressure) in
                field.iter().zip(surface.areas()).zip(&solution.pressure)
            {
                value = value + pressure.scale(shape.re * area);
            }
            matrix[i * n + j] = value;
        }
    }
    // Keep the computed matrix untouched. The shared positive-real fit checks
    // against these same complex and resistance samples; it may refuse them.
    Ok((matrix, ppw, condition))
}

// Exact oriented edge closure, independently for each connected skin. This is
// deliberately stricter than mortar mechanics: an acoustic hole is not sealed
// merely because two mechanical nodes are constrained to move together.
fn closed_outward(triangles: &[[[f64; 3]; 3]]) -> Result<(), ShellBemError> {
    if triangles.len() < 4 {
        return Err(ShellBemError::Invalid(
            "closed exterior requires at least four triangles",
        ));
    }
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::<(usize, usize), Vec<(usize, bool)>>::new();
    let mut neighbors = vec![Vec::new(); triangles.len()];
    for (face, triangle) in triangles.iter().enumerate() {
        let mut ids = [0; 3];
        for (i, point) in triangle.iter().enumerate() {
            if point.iter().any(|x| !x.is_finite()) {
                return Err(ShellBemError::Invalid("nonfinite exterior vertex"));
            }
            let key = point.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
            let next = vertices.len();
            ids[i] = *vertices.entry(key).or_insert(next);
        }
        for i in 0..3 {
            let a = ids[i];
            let b = ids[(i + 1) % 3];
            if a == b {
                return Err(ShellBemError::Invalid("degenerate exterior edge"));
            }
            edges
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((face, a < b));
        }
    }
    for sides in edges.values() {
        if sides.len() != 2 || sides[0].1 == sides[1].1 {
            return Err(ShellBemError::Invalid(
                "BEM exterior requires a closed consistently oriented skin; open holes or mismatched seams need an explicit acoustic boundary",
            ));
        }
        neighbors[sides[0].0].push(sides[1].0);
        neighbors[sides[1].0].push(sides[0].0);
    }
    let mut seen = vec![false; triangles.len()];
    for first in 0..triangles.len() {
        if seen[first] {
            continue;
        }
        let origin = triangles[first][0];
        let mut stack = vec![first];
        let mut volume6 = 0.0;
        while let Some(face) = stack.pop() {
            if seen[face] {
                continue;
            }
            seen[face] = true;
            let [a, b, c] =
                triangles[face].map(|p| core::array::from_fn::<_, 3, _>(|i| p[i] - origin[i]));
            volume6 += a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            stack.extend(neighbors[face].iter().copied());
        }
        if !volume6.is_finite() || volume6 <= 0.0 {
            return Err(ShellBemError::Invalid(
                "every closed exterior skin must enclose positive outward-oriented volume",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::GasSpec;

    #[test]
    fn signed_zero_volume_fields_keep_mutual_air_reaction_and_the_selected_medium() {
        let gate = CancelGate::new_clock_free();
        let surface = SpherePanels::icosphere(0.1, 1).unwrap();
        // A translating sphere is a dipole: zero total displaced volume does
        // not mean zero exterior pressure, inertia, or radiation resistance.
        let field: Vec<_> = surface
            .normals()
            .iter()
            .map(|n| C64::from_re(n[0]))
            .collect();
        let volume: f64 = field
            .iter()
            .zip(surface.areas())
            .map(|(v, a)| v.re * a)
            .sum();
        assert!(volume.abs() < 1e-15);
        let opposite: Vec<_> = field.iter().map(|v| v.scale(-1.0)).collect();
        let fields = vec![field, opposite];
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let medium = Medium::from(&gas);
        let omega = medium.sound_speed / 0.1;
        let (z, _, _) = modal_impedance(&surface, &fields, omega, medium, &gate).unwrap();
        assert!(z[0].re > 0.0 && z[0].abs() > 1e-6);
        for (a, b) in [(z[0], z[3]), (z[0], z[1].scale(-1.0)), (z[1], z[2])] {
            assert!((a - b).abs() < 1e-12 * z[0].abs());
        }
        let (denser, _, _) = modal_impedance(
            &surface,
            &fields,
            omega,
            Medium {
                density: 2.0 * medium.density,
                ..medium
            },
            &gate,
        )
        .unwrap();
        for (a, b) in z.iter().zip(&denser) {
            assert!((*b - a.scale(2.0)).abs() < 1e-12 * a.abs());
        }
        let (faster, _, _) = modal_impedance(
            &surface,
            &fields,
            omega,
            Medium {
                sound_speed: 2.0 * medium.sound_speed,
                ..medium
            },
            &gate,
        )
        .unwrap();
        assert!((faster[0] - z[0]).abs() > 0.01 * z[0].abs());
        // This is a phase/impedance change, not a scalar gain adjustment.
        assert!((faster[0].re / faster[0].im - z[0].re / z[0].im).abs() > 0.01);
    }

    #[test]
    fn every_surface_edge_and_outward_component_is_required() {
        let surface = SpherePanels::icosphere(0.1, 0).unwrap();
        let triangles = surface.triangles().unwrap();
        closed_outward(triangles).unwrap();
        assert!(closed_outward(&triangles[1..]).is_err());
        let mut inside_out = triangles.to_vec();
        for tri in &mut inside_out {
            tri.swap(1, 2);
        }
        assert!(closed_outward(&inside_out).is_err());
        let mut bad_seam = triangles.to_vec();
        bad_seam[0][0][0] += 1e-12;
        assert!(closed_outward(&bad_seam).is_err());
        let mut duplicate = triangles.to_vec();
        duplicate.push(triangles[0]);
        assert!(closed_outward(&duplicate).is_err());
        assert!(closed_outward(&[]).is_err());
    }
}
