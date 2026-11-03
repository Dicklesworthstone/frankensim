//! Frequency-dependent material and boundary loading in one original shell solve.
//!
//! The shell's equilibrium basis is unchanged. Full film/substrate storage and
//! loss are evaluated at the drive frequency before the shared harmonic owner
//! adds -i omega B Z B^T and factors the joint operator. Boundary work remains
//! distinct from internal material heat. This does not infer a radiation law,
//! causality, broadband load passivity or modal truncation accuracy.

use super::{
    Image, MaterialOperator, ShellCoatingHarmonicOptions, ShellCoatingHarmonicResponse, finite,
    material_operator, observe_image, source_support,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::shear_shell::coating::strain::{
    PreparedShellCoatingStrain, ShellCoatingStrainBudget, ShellCoatingStrainReport,
};
use crate::shear_shell::{
    MaterialShellImage, assembly::MaterialShellAssemblyImage, invalid, owner, poll,
};
use fs_exec::{CancelGate, Cx};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicPortControls, HarmonicPortResponse, HarmonicResponseControls, HarmonicResponseError,
    PreparedHarmonicResponse,
};

/// Original work-conjugate row; input order is the impedance's port order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShellHarmonicPort {
    /// An original force port: opposing force [N], velocity [m/s].
    Force(usize),
    /// Selected exterior surface: opposing pressure [Pa], outward flow [m³/s].
    ExteriorVolume,
    /// Selected cavity surface: opposing pressure [Pa], outward-from-cavity
    /// volume flow [m³/s]. The structural reaction is minus this row times
    /// the opposing pressure; its sign is not an absolute cavity pressure.
    CavityVolume,
}

/// Separate source assembly, shared factor, physical load and strain budgets.
#[derive(Clone, Copy, Debug)]
pub struct ShellLoadedHarmonicOptions {
    /// Existing complete source material work allowance. Its solve tolerance
    /// also admits the source query; the shared solve uses `response` below.
    pub material: ShellCoatingHarmonicOptions,
    /// Additional numeric input payload: 40*n² + 8*n*p + 8*films bytes.
    /// Includes complex source assembly and real M/C/K/B inputs at their peak;
    /// excludes borrowed image storage and reference/vector bookkeeping.
    pub maximum_material_payload_bytes: usize,
    /// Complete shared harmonic factor/solve admission.
    pub response: HarmonicResponseControls,
    /// Complete physical impedance projection and pointwise PSD admission.
    pub ports: HarmonicPortControls,
    /// Original whole-film spatial/phase envelope plan and query allowance.
    pub strain: ShellCoatingStrainBudget,
}

/// One frequency and one unchanged full boundary impedance, reusable for loads.
pub struct PreparedLoadedShellCoatingResponse<'a> {
    image: Image<'a>,
    omega: f64,
    losses: Vec<f64>,
    substrate_loss: f64,
    response: PreparedHarmonicResponse,
    strain: PreparedShellCoatingStrain<'a>,
}

/// Source-admitted particular solution and separate physical power accounts.
#[derive(Debug)]
pub struct LoadedShellCoatingHarmonicResponse<'a> {
    /// Original shell observations and internal material heat. The embedded
    /// `power_residual_w` excludes boundary work; use the joint residual below.
    pub shell: ShellCoatingHarmonicResponse,
    /// Opposing load efforts, conjugate flows and net absorbed boundary work.
    /// Absorbed work is not automatically material heat or radiated sound.
    pub boundary: HarmonicPortResponse,
    /// Whole-film, whole-cycle bounds against every authored film envelope.
    pub coating_strain: ShellCoatingStrainReport<'a>,
    /// Internal material heat + absorbed boundary work - supplied work [W].
    pub power_balance_residual_w: f64,
}

fn harmonic_error(error: HarmonicResponseError) -> AcousticRealizeError {
    match error {
        HarmonicResponseError::Cancelled => AcousticRealizeError::Cancelled,
        error => owner(error),
    }
}

fn poll_both(gate: &CancelGate, cx: &Cx<'_>) -> Result<(), AcousticRealizeError> {
    poll(gate)?;
    cx.checkpoint().map_err(|_| AcousticRealizeError::Cancelled)
}

impl MaterialShellImage {
    /// Compose original force or signed pressure-volume ports with a full
    /// physical impedance evaluated by its owner at `omega_rad_s`.
    ///
    /// Peak phasors use exp(-i omega t). `impedance` is row-major p*p and
    /// maps each selected flow to its opposing effort; reactive, mutual and
    /// passive nonreciprocal terms are retained. Source films/substrate must
    /// admit this frequency. Pointwise load passivity is checked by fs-modal;
    /// the caller owns its physical units, frequency support and causality.
    ///
    /// `gate` controls original material/strain work; `cx` controls the shared
    /// solve. Both are checked before publication. Usually they share a gate.
    /// All returned motions must pass original whole-film strain admission.
    pub fn prepare_coating_harmonic_with_impedance(
        &self,
        omega_rad_s: f64,
        ports: &[ShellHarmonicPort],
        impedance: &[C64],
        options: ShellLoadedHarmonicOptions,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<PreparedLoadedShellCoatingResponse<'_>, AcousticRealizeError> {
        poll_both(gate, cx)?;
        source_support(self.source(), omega_rad_s, gate)?;
        if self.coating_relaxation().map_or(0, |m| m.films().len())
            != self.source().coatings().len()
        {
            return Err(invalid(
                "loaded shell requires every source film projection",
            ));
        }
        let image = Image::shell(self);
        let (response, losses, substrate_loss) =
            prepare(&image, omega_rad_s, ports, impedance, options, gate, cx)?;
        let strain = self.prepare_coating_strain(options.strain, gate)?;
        poll_both(gate, cx)?;
        Ok(PreparedLoadedShellCoatingResponse {
            image,
            omega: omega_rad_s,
            losses,
            substrate_loss,
            response,
            strain,
        })
    }
}

impl MaterialShellAssemblyImage {
    /// Same joint source/material/boundary solve on the original monolithic
    /// tied-shell basis, preserving every component's film and fluid orientation.
    /// See [`MaterialShellImage::prepare_coating_harmonic_with_impedance`].
    pub fn prepare_coating_harmonic_with_impedance(
        &self,
        omega_rad_s: f64,
        ports: &[ShellHarmonicPort],
        impedance: &[C64],
        options: ShellLoadedHarmonicOptions,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<PreparedLoadedShellCoatingResponse<'_>, AcousticRealizeError> {
        poll_both(gate, cx)?;
        for component in &self.source().components {
            source_support(&component.specimen, omega_rad_s, gate)?;
        }
        let expected: usize = self
            .source()
            .components
            .iter()
            .map(|component| component.specimen.coatings().len())
            .sum();
        if self.coating_relaxation().map_or(0, |m| m.films().len()) != expected {
            return Err(invalid(
                "loaded assembly requires every source film projection",
            ));
        }
        let image = Image::assembly(self);
        let (response, losses, substrate_loss) =
            prepare(&image, omega_rad_s, ports, impedance, options, gate, cx)?;
        let strain = self.prepare_coating_strain(options.strain, gate)?;
        poll_both(gate, cx)?;
        Ok(PreparedLoadedShellCoatingResponse {
            image,
            omega: omega_rad_s,
            losses,
            substrate_loss,
            response,
            strain,
        })
    }
}

#[allow(clippy::too_many_arguments)] // complete source and independent numerical owners
fn prepare(
    image: &Image<'_>,
    omega: f64,
    ports: &[ShellHarmonicPort],
    impedance: &[C64],
    options: ShellLoadedHarmonicOptions,
    gate: &CancelGate,
    cx: &Cx<'_>,
) -> Result<(PreparedHarmonicResponse, Vec<f64>, f64), AcousticRealizeError> {
    poll_both(gate, cx)?;
    let (n, p) = (image.modes.len(), ports.len());
    let films = image.films.map_or(0, |memory| memory.films().len());
    let square = n
        .checked_mul(n)
        .ok_or_else(|| invalid("loaded shell extent overflow"))?;
    let np = n
        .checked_mul(p)
        .ok_or_else(|| invalid("loaded shell port extent overflow"))?;
    let bytes = square
        .checked_mul(40)
        .and_then(|b| np.checked_mul(8).and_then(|x| b.checked_add(x)))
        .and_then(|b| films.checked_mul(8).and_then(|x| b.checked_add(x)))
        .ok_or_else(|| invalid("loaded shell source payload overflow"))?;
    if n == 0
        || n > options.response.maximum_coordinates
        || p == 0
        || p > options.ports.maximum_ports
        || p.checked_mul(p) != Some(impedance.len())
        || impedance.iter().any(|&z| !finite(z))
        || bytes > options.maximum_material_payload_bytes
    {
        return Err(invalid(
            "loaded shell requires complete ports and source payload allowance",
        ));
    }
    let mut force_map = vec![0.0; np];
    for (port, selection) in ports.iter().enumerate() {
        poll_both(gate, cx)?;
        let row = match *selection {
            ShellHarmonicPort::Force(index) => image
                .force_ports
                .get(index)
                .map(Vec::as_slice)
                .ok_or_else(|| invalid("unknown original shell force port"))?,
            ShellHarmonicPort::ExteriorVolume => image.exterior,
            ShellHarmonicPort::CavityVolume => image.cavity,
        };
        if row.len() != n {
            return Err(invalid("incomplete original shell load row"));
        }
        for i in 0..n {
            force_map[i * p + port] = row[i];
        }
    }
    let MaterialOperator {
        operator,
        losses,
        substrate_loss,
    } = material_operator(image, omega, options.material, false, gate)?;
    let mut mass = vec![0.0; square];
    for i in 0..n {
        mass[i * n + i] = 1.0;
    }
    let storage: Vec<_> = operator.iter().map(|z| z.re).collect();
    let mut damping = Vec::with_capacity(square);
    for z in &operator {
        let c = -z.im / omega;
        if !c.is_finite() || z.im != 0.0 && c == 0.0 {
            return Err(invalid(
                "loaded shell effective material loss is unrepresentable",
            ));
        }
        damping.push(c);
    }
    poll_both(gate, cx)?;
    let response = PreparedHarmonicResponse::try_with_port_impedance(
        &mass,
        &damping,
        &storage,
        n,
        &force_map,
        p,
        impedance,
        omega,
        options.response,
        options.ports,
        cx,
    )
    .map_err(harmonic_error)?;
    poll_both(gate, cx)?;
    Ok((response, losses, substrate_loss))
}

impl<'a> PreparedLoadedShellCoatingResponse<'a> {
    /// Actual drive frequency; changing it requires reevaluating original laws.
    pub fn omega_rad_s(&self) -> f64 {
        self.omega
    }

    /// Completed original strain queries, including sufficient-envelope refusals.
    pub fn completed_strain_evaluations(&self) -> usize {
        self.strain.completed_evaluations()
    }

    /// Solve another complete generalized peak-force phasor [N/sqrt(kg)] using
    /// the same factor. Material heat and boundary work use the resulting joint
    /// motion. A strain bound outside an authored film envelope refuses the
    /// response and consumes its successfully evaluated strain query; it is not
    /// evidence of actual material failure. Source histories remain unchanged.
    pub fn solve(
        &mut self,
        force: &[C64],
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<LoadedShellCoatingHarmonicResponse<'a>, AcousticRealizeError> {
        poll_both(gate, cx)?;
        let solved = self.response.solve(force, cx).map_err(harmonic_error)?;
        let boundary = solved
            .port_load
            .ok_or_else(|| invalid("loaded shell omitted boundary work"))?;
        let shell = observe_image(
            &self.image,
            self.omega,
            force,
            &self.losses,
            self.substrate_loss,
            solved.displacement,
            solved.velocity,
            solved.relative_residual,
            gate,
        )?;
        let power_balance_residual_w = shell.power_residual_w + boundary.average_absorbed_power_w;
        if !power_balance_residual_w.is_finite() {
            return Err(invalid("loaded shell power accounting overflow"));
        }
        poll_both(gate, cx)?;
        let coating_strain = self.strain.evaluate(&shell.displacement_m_sqrt_kg, gate)?;
        if !coating_strain.within_authored_envelope {
            return Err(invalid(
                "loaded shell motion exceeds the sufficient original film strain envelope",
            ));
        }
        poll_both(gate, cx)?;
        Ok(LoadedShellCoatingHarmonicResponse {
            shell,
            boundary,
            coating_strain,
            power_balance_residual_w,
        })
    }
}
