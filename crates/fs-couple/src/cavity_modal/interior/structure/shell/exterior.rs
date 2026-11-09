//! Interior pressure, original material-shell motion and exterior radiation.
//!
//! One native factor retains all free shell translations/rotations and all
//! original pressure-potential coordinates. The native exterior impedance adds
//! `+Z_ext v` to the structural equation, while interior pressure supplies
//! `-B^T p`. Thus fluid force on the structure is `B^T p - Z_ext v`.
//! Neither the dry shell nor the interior field is inverted independently.
//! The same exterior source columns supply coherent pressure, particle velocity
//! and outgoing flux after the jointly loaded velocity has been solved.
//!
//! Both sides use peak exp(-i omega t) phasors. Interior-outward normals must
//! agree with the original closed source's outward facet orientation; pressure
//! acts on opposite shell faces. The zero-offset source midsurface separates
//! the fluids, whose material properties may differ. This does not add an
//! aperture, open screen, thickness skin, material-memory time realization or
//! continuum error bound. The source owners retain their original band limits.

use super::{
    PreparedStructuralInterior, StructuralInteriorError, TensorShellInterior,
    TensorShellInteriorError, TensorShellInteriorResponse, cap, poll, relative,
};
use crate::tensor_shell_radiation::{
    PreparedTensorShellRadiation, TensorShellRadiationAcousticControls, TensorShellRadiationError,
};
use fs_bem::near_field::Options as ReceiverOptions;
use fs_bem::reciprocal_radiation::{
    CoherentRadiationPower, CoherentVelocityPortField, PreparedReceivers, VelocityPortSolution,
};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, VelocityTrace, VolumeVelocitySource,
};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicPortControls, PreparedHarmonicResponse};

/// Additional exterior source, physical load projection and observation limits.
/// Existing material and complete coupled-coordinate limits remain authoritative.
#[derive(Clone, Copy, Debug)]
pub struct TensorShellInteriorExteriorControls {
    pub acoustic: TensorShellRadiationAcousticControls,
    /// Native full physical-port projection and Hermitian PSD admission.
    pub ports: HarmonicPortControls,
    /// Additional temporary real map: 8*n*q bytes, with n complete coupled
    /// coordinates and q original shell coordinates. Native copies/projection
    /// are separately admitted by `ports`; allocator metadata is excluded.
    pub maximum_force_map_bytes: usize,
    /// Relative exterior port/surface and enclosing input/heat/flux defects,
    /// scaled by apparent work. The interior owner's tolerance independently
    /// applies to structural/interior balances including the exterior load.
    pub power_relative_tolerance: f64,
}

impl<'shell, 'a, 'f, 'c> TensorShellInterior<'shell, 'a, 'f, 'c> {
    /// Attach the closed exterior of this EXACT original material shell.
    /// All exterior panels remain, including fully supported scattering faces.
    /// Only the already declared interior interfaces move the interior fluid;
    /// an empty map is the explicit detached-cavity limit.
    ///
    /// No acoustic modal reduction, bare-system inverse, duplicate cavity
    /// spring or extra radiation debit is introduced. Repeated coherent force
    /// and inward-flow cases reuse one complete loaded factor and BEM source set.
    pub fn prepare_with_exterior<'p, 'geometry>(
        &'p self,
        exterior: &'p PreparedTensorShellRadiation<'geometry, 'shell>,
        controls: TensorShellInteriorExteriorControls,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedTensorShellInteriorExterior<'p, 'geometry, 'shell, 'a, 'f, 'c>,
        TensorShellInteriorError,
    > {
        poll(cx)?;
        if !core::ptr::eq(self.material.source(), exterior.source().shell()) {
            return Err(TensorShellInteriorError::Invalid(
                "interior and exterior must retain the same original material shell",
            ));
        }
        if !self.exterior_side_consistent {
            return Err(TensorShellInteriorError::Invalid(
                "cavity-outward normals must agree with the original exterior shell orientation",
            ));
        }
        let tolerance = controls.power_relative_tolerance;
        if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0 {
            return Err(TensorShellInteriorError::Invalid(
                "finite exterior power tolerance in (0,1) required",
            ));
        }
        let q = self.structural_dofs();
        let n = q + self.pressure_dofs();
        let (nn, qq) = (n as u128, q as u128);
        cap(
            "exterior original shell ports",
            qq,
            controls.ports.maximum_ports as u128,
        )?;
        cap(
            "exterior complete force-map bytes",
            8 * nn * qq,
            controls.maximum_force_map_bytes as u128,
        )?;
        // The same native dense-port bounds are checked before any BEM work.
        cap(
            "exterior full load projection",
            nn * qq * (nn + qq),
            controls.ports.maximum_projection_terms as u128,
        )?;
        cap(
            "exterior native load payload",
            64 * nn * qq + 256 * qq * qq + 128 * qq,
            controls.ports.maximum_payload_bytes as u128,
        )?;
        let omega = self.material.omega_rad_s();
        let radiation = exterior.prepare_acoustics(omega, controls.acoustic, cx)?;
        let mut force_map = vec![0.0; n * q];
        for j in 0..q {
            poll(cx)?;
            force_map[j * q + j] = 1.0;
        }
        let harmonic = PreparedHarmonicResponse::try_with_port_impedance(
            &self.model.mass,
            &self.model.damping,
            &self.model.stiffness,
            n,
            &force_map,
            q,
            radiation.impedance(),
            omega,
            self.model.controls,
            controls.ports,
            cx,
        )
        .map_err(StructuralInteriorError::from)?;
        poll(cx)?;
        Ok(PreparedTensorShellInteriorExterior {
            owner: self,
            exterior,
            radiation,
            coupled: PreparedStructuralInterior {
                model: &self.model,
                harmonic,
            },
            controls,
        })
    }
}

/// Full loaded spatial factor and the same source-bound exterior acoustic fields.
/// No frozen-coefficient time or frequency-derivative API is exposed.
pub struct PreparedTensorShellInteriorExterior<'p, 'geometry, 'shell, 'a, 'f, 'c> {
    owner: &'p TensorShellInterior<'shell, 'a, 'f, 'c>,
    exterior: &'p PreparedTensorShellRadiation<'geometry, 'shell>,
    radiation: VelocityPortSolution<'p, 'geometry>,
    coupled: PreparedStructuralInterior<'p, 'a, 'f, 'c>,
    controls: TensorShellInteriorExteriorControls,
}

impl<'p, 'geometry> PreparedTensorShellInteriorExterior<'p, 'geometry, '_, '_, '_, '_> {
    /// Original full shell-coordinate resisting impedance, with all mutual terms.
    pub fn impedance(&self) -> &[C64] {
        self.radiation.impedance()
    }

    /// Observe pressure/particle velocity from the same exterior source frequency
    /// and medium. Clearances and aggregate quadrature retain the native owner.
    pub fn prepare_receivers(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: ReceiverOptions,
        particle_velocity: bool,
        cx: &Cx<'_>,
    ) -> Result<PreparedReceivers<'p, 'geometry>, TensorShellInteriorError> {
        poll(cx)?;
        let medium = self.exterior.source().medium();
        // Bind the borrow to the original source cache, not this observation call.
        let boundary: &'p _ = self.exterior.boundary();
        Ok(boundary
            .receivers(points, minimum_clearance_m)
            .map_err(TensorShellRadiationError::from)?
            .prepare(
                self.owner.material.omega_rad_s() / medium.sound_speed,
                medium,
                options,
                particle_velocity,
                || cx.checkpoint().is_err(),
            )
            .map_err(TensorShellRadiationError::from)?)
    }

    /// Drive original distributed injection through the SAME complete loaded
    /// shell/interior factor and all exterior self/mutual reaction terms.
    /// Prepared shapes belong to the original heterogeneous interior; all
    /// amplitudes use exp(-i omega t). The real native projection adds no
    /// frequency, density, area normalization or phasor conjugation.
    /// Signed physical source work is included once in acoustic input, while
    /// original material heat, wall heat and independently observed exterior
    /// outgoing flux retain their existing separate admission.
    pub fn solve_with_sources(
        &self,
        structural_force: &[C64],
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cx: &Cx<'_>,
    ) -> Result<TensorShellInteriorExteriorResponse<'p, 'geometry>, TensorShellInteriorError> {
        let load = self.owner.model.project_sources(
            inward_velocity,
            boundary_sources,
            volume_sources,
            cx,
        )?;
        self.solve(structural_force, &load, cx)
    }

    /// Solve coherent external forces/torques and original inward flow loads.
    /// Intrinsic material heat, interior wall heat and exterior outgoing flux
    /// remain separate. Only accepted original strains and complete work balances
    /// publish a response; all state and source preparations are immutable.
    pub fn solve(
        &self,
        structural_force: &[C64],
        inward_flow_load: &[C64],
        cx: &Cx<'_>,
    ) -> Result<TensorShellInteriorExteriorResponse<'p, 'geometry>, TensorShellInteriorError> {
        poll(cx)?;
        let fields = self.coupled.solve(structural_force, inward_flow_load, cx)?;
        let coupled = self
            .owner
            .admit_response(fields, structural_force, inward_flow_load, cx)?;
        let exterior = self
            .radiation
            .coherent(
                &coupled.field.structural_velocity,
                self.controls.acoustic.field,
                || cx.checkpoint().is_err(),
            )
            .map_err(TensorShellRadiationError::from)?;
        let radiation = exterior
            .power(self.controls.acoustic.maximum_radiated_power_work, || {
                cx.checkpoint().is_err()
            })
            .map_err(TensorShellRadiationError::from)?;
        let port = coupled.field.structural_port_load.as_ref().ok_or(
            TensorShellInteriorError::Invalid("loaded spatial solve omitted its exterior work"),
        )?;
        let scale = coupled.apparent_power_scale_w + radiation.outgoing_power_w.abs();
        let port_work = relative(
            (port.average_absorbed_power_w - radiation.surface_power_w.re).abs(),
            scale,
        );
        let power_residual = coupled.field.average_structural_input_power_w
            + coupled.field.average_acoustic_input_power_w
            - coupled.field.average_structural_dissipation_w
            - coupled.field.average_wall_dissipation_w
            - radiation.outgoing_power_w;
        let total = relative(power_residual.abs(), scale);
        let tolerance = self.controls.power_relative_tolerance;
        if port_work > tolerance
            || total > tolerance
            || !radiation.relative_flux_defect.is_finite()
            || radiation.relative_flux_defect > tolerance
        {
            return Err(TensorShellRadiationError::PowerBalance {
                coupled: total,
                port_work,
                flux: radiation.relative_flux_defect,
                tolerance,
            }
            .into());
        }
        poll(cx)?;
        Ok(TensorShellInteriorExteriorResponse {
            coupled,
            exterior,
            radiation,
            power_balance_residual_w: power_residual,
            relative_power_balance: total,
            relative_port_work_defect: port_work,
        })
    }
}

/// One jointly loaded shell/interior field and its independently observable exterior.
pub struct TensorShellInteriorExteriorResponse<'p, 'geometry> {
    /// Original motion, pressure, material/wall heat and signed interior exchange.
    /// Its structural-port effort is the opposing exterior load, not extra heat.
    pub coupled: TensorShellInteriorResponse,
    /// Coherent original boundary field for native pressure/velocity receivers.
    pub exterior: CoherentVelocityPortField<'p, 'geometry>,
    /// Boundary pressure work and independently integrated outgoing flux.
    pub radiation: CoherentRadiationPower,
    /// Both source powers minus material heat, wall heat and outgoing flux.
    pub power_balance_residual_w: f64,
    pub relative_power_balance: f64,
    pub relative_port_work_defect: f64,
}
