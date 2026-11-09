//! Moving translations and physical rotations in the original MITC coordinates.
//!
//! The native supported harmonic owner eliminates prescribed coordinates only
//! after assembling the complete dispersive material and reference inertia.
//! In particular, K_fc, C_fc AND M_fc drive the free response. Source admission
//! observes the complete displacement, including every moving support.

use super::{
    HarmonicMitcShell, PreparedHarmonicMitcOperator, Result, admit_displacement, cap, checkpoint,
    dense, invalid,
};
use fs_exec::Cx;
use fs_material::laminate::harmonic::PreparedHarmonicLaminate;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseControls, PreparedSupportedHarmonicResponse, SupportedHarmonicResponse,
};

/// One frequency's supported factor and its complete original material field.
/// Prescribed labels are original 6*node+component SUPPORT coordinates: the
/// original node frame governs both translations and axial rotations. These
/// labels are not global Cartesian components when that frame is oblique.
pub struct PreparedSupportedHarmonicMitc<'a, 'model> {
    source: &'a HarmonicMitcShell<'model>,
    omega: f64,
    materials: Vec<PreparedHarmonicLaminate<'a>>,
    original_support_dofs: Vec<usize>,
    response: PreparedSupportedHarmonicResponse,
}

fn admit_supported(
    source: &HarmonicMitcShell<'_>,
    omega: f64,
    original_support_dofs: &[usize],
    controls: HarmonicResponseControls,
    cx: &Cx<'_>,
) -> Result<usize> {
    let material_bytes = source.frequency_payload(omega, controls, false, cx)?;
    // Validate before allocating the map or evaluating material spectra.
    // Previously eliminated homogeneous constraints cannot become moving ones.
    let mut previous = None;
    for &label in original_support_dofs {
        checkpoint(cx)?;
        if previous.is_some_and(|old| old >= label)
            || !matches!(source.model.pencil.dof_map.get(label), Some(Some(_)))
        {
            return Err(invalid(
                "MITC moving supports require strictly increasing retained original labels",
            ));
        }
        previous = Some(label);
    }
    let n = source.model.pencil.free;
    // Native supported payload: 320*n*n + 512*n. Three input dense matrices
    // coexist with it; 32*n also covers both temporary and retained label maps.
    let bytes = n
        .checked_mul(n)
        .and_then(|v| v.checked_mul(344))
        .and_then(|v| n.checked_mul(544).and_then(|w| v.checked_add(w)))
        .and_then(|v| v.checked_add(material_bytes))
        .ok_or_else(|| invalid("MITC supported harmonic payload extent overflows"))?;
    cap(
        "MITC supported harmonic payload bytes",
        bytes as u128,
        controls.maximum_payload_bytes as u128,
    )?;
    Ok(material_bytes)
}

impl<'model> HarmonicMitcShell<'model> {
    /// Prepare harmonic motion driven by original forces AND moving supports.
    ///
    /// Labels are strictly increasing original 6*node+component support axes,
    /// and must still be retained by the original model. Existing eliminated
    /// coordinates stay homogeneous. For a support that will move, retain that
    /// coordinate in the original model and prescribe it here.
    ///
    /// The original mass, full tensor storage/loss, drilling and connections
    /// assemble before support elimination. Only the free principal block is
    /// factored; a singular bare full operator need not prevent this response.
    /// Frequency remains strictly positive and inside every material band.
    /// Empty and all-prescribed maps use the same native owner.
    pub fn prepare_supported_harmonic_response(
        &self,
        omega: f64,
        original_support_dofs: &[usize],
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedSupportedHarmonicMitc<'_, 'model>> {
        admit_supported(self, omega, original_support_dofs, controls, cx)?;
        self.prepare_harmonic_operator(omega, controls, cx)?
            .factor_supported(original_support_dofs, controls, cx)
    }
}

impl<'a, 'model> PreparedHarmonicMitcOperator<'a, 'model> {
    /// Reuse this exact prepared material pencil with prescribed support motion.
    /// No additional constitutive preparation or bare full factorization occurs.
    pub fn factor_supported(
        self,
        original_support_dofs: &[usize],
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedSupportedHarmonicMitc<'a, 'model>> {
        let material_bytes =
            admit_supported(self.source, self.omega, original_support_dofs, controls, cx)?;
        let n = self.source.model.pencil.free;
        let source_bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(24))
            .and_then(|v| n.checked_mul(32).and_then(|w| v.checked_add(w)))
            .and_then(|v| v.checked_add(material_bytes))
            .ok_or_else(|| invalid("MITC supported source payload extent overflows"))?;
        let mut support_dofs = Vec::with_capacity(original_support_dofs.len());
        for &label in original_support_dofs {
            checkpoint(cx)?;
            support_dofs.push(
                self.source.model.pencil.dof_map[label]
                    .ok_or_else(|| invalid("MITC moving support was eliminated"))?,
            );
        }
        let mass = dense(&self.source.model.pencil.m, n, cx)?;
        let storage = dense(&self.stiffness, n, cx)?;
        let viscosity = dense(&self.damping, n, cx)?;
        let Self {
            source,
            omega,
            materials,
            stiffness,
            damping,
        } = self;
        drop(stiffness);
        drop(damping);
        let response = PreparedSupportedHarmonicResponse::try_new(
            &mass,
            &viscosity,
            &storage,
            n,
            &support_dofs,
            omega,
            HarmonicResponseControls {
                maximum_payload_bytes: controls.maximum_payload_bytes - source_bytes,
                ..controls
            },
            cx,
        )?;
        checkpoint(cx)?;
        Ok(PreparedSupportedHarmonicMitc {
            source,
            omega,
            materials,
            original_support_dofs: original_support_dofs.to_vec(),
            response,
        })
    }
}

impl<'a, 'model> PreparedSupportedHarmonicMitc<'a, 'model> {
    pub const fn source(&self) -> &'a HarmonicMitcShell<'model> {
        self.source
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    /// Original 6*node+component support-axis labels, in drive/reaction order.
    pub fn original_support_dofs(&self) -> &[usize] {
        &self.original_support_dofs
    }
    /// Corresponding rows in the original model's retained-coordinate pencil.
    pub fn retained_support_dofs(&self) -> &[usize] {
        self.response.support_dofs()
    }
    /// Unprescribed rows in that same retained-coordinate pencil.
    pub fn free_dofs(&self) -> &[usize] {
        self.response.free_dofs()
    }

    /// Solve full retained-coordinate force/moment phasors and one displacement
    /// per prescribed original label, under exp(-i omega t). Returned motion
    /// retains ALL original pencil coordinates, including moving supports.
    ///
    /// Reactions are forces ON the structure, in original_support_dofs() order.
    /// Translation/force and axial rotation/moment share the original node axes.
    /// Applied forces at prescribed coordinates are subtracted exactly once.
    /// Signed support work, applied-force work and physical material heat remain
    /// independently computed by the shared supported owner. Fixed drilling and
    /// springs supply no heat. No actuator dynamics or startup history is added.
    ///
    /// Complete original layer-strain admission precedes publication. A failed
    /// drive or cancellation leaves the reusable factor and material unchanged.
    pub fn solve(
        &self,
        force: &[C64],
        support_displacement: &[C64],
        cx: &Cx<'_>,
    ) -> Result<SupportedHarmonicResponse> {
        let response = self.response.solve(force, support_displacement, cx)?;
        admit_displacement(self.source, &self.materials, &response.displacement, cx)?;
        Ok(response)
    }
}
