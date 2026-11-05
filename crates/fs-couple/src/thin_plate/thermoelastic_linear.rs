//! Coupled thermoelastic memory on the original linear DKT modal basis.
//!
//! The elastic pencil and physical force/volume maps are harvested once.
//! Thermal factors are projected through those same eigenvectors, then the
//! existing complete held-force owner advances q/v and every thermal memory.
//! No diagonal thermal damping ratio or second eigensolve is substituted.

use super::thermoelastic_memory::{
    ThermoelasticMemoryOptions, admit_memory_band, project_thermal_relaxation,
};
use super::{
    AcousticRealizeError, AnisotropicPlateThermal, AssemblyOptions, CompactBody, EdgeSupport,
    PlateMesh, PlateProjection, ThinPlate, assemble, harvest_radiators_with_modes, map_plate,
    plane_stress_section, ss_omega11,
};
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use crate::render::schedule::force::coupled::{
    CoupledModalSystem, ModalCouplingConfig, PreparedModalBudget,
};
use fs_exec::CancelGate;
use fs_material::visco::slab_ted::SlabPronySpectrum;
use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;

/// A homogeneous plate's complete retained mechanical and thermal image.
///
/// Material parameters are caller-supplied isothermal elasticity and thermal
/// data at one frozen state. The original rectangular DKT mesh/window remains
/// the spatial approximation. The finite spectrum bounds only its normalized
/// constitutive response over the declared band; it is not a modal-truncation,
/// projection-roundoff or arbitrary-force time-domain error certificate.
#[derive(Clone)]
pub struct ThermoelasticPlateImage {
    radiators: Vec<CompactBody>,
    branches: Vec<RelaxationBranch>,
    force_shapes: Vec<f64>,
    volume_shapes: Vec<f64>,
    spectrum: SlabPronySpectrum,
}

impl ThermoelasticPlateImage {
    /// Retain all requested original in-window modes, up to 64, with the full
    /// thermal normal/shear channel and its constitutive memory.
    ///
    /// The original plate must be linear with no separately selected damping.
    /// Thermal memory is the sole internal loss here; additional contact,
    /// support, intrinsic or radiation mechanisms belong to their original
    /// explicit coupled owners. No exterior medium or pressure law is inferred.
    /// Cancellation or admission refusal publishes no image.
    ///
    /// # Errors
    /// Geometry/material/mode refusal, duplicate loss, exhausted projection or
    /// thermal spectrum allowance, or represented instantaneous linear stiffness
    /// outside the declared band.
    pub fn from_plate(
        plate: ThinPlate,
        thermal: AnisotropicPlateThermal,
        options: ThermoelasticMemoryOptions,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if plate.geometric_nonlinearity
            || !(1..=64).contains(&plate.n_modes)
            || plate.kelvin_voigt_bending.is_some()
        {
            return Err(invalid(
                "linear thermal memory needs 1..=64 original modes and no parallel scalar viscosity",
            ));
        }
        let section = plane_stress_section(plate)?;
        let slab = thermal.resolve_coupled(plate, &section)?;
        let spectrum = slab
            .lower_to_prony(
                options.band_rad_s,
                options.normalized_tolerance,
                options.max_thermal_terms,
            )
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let (nx, ny) = (5, 4);
        let mesh = PlateMesh::rectangle(plate.length_m, plate.width_m, nx, ny);
        let boundary = PlateMesh::rectangle_boundary(nx, ny);
        let model = assemble(
            &mesh,
            &section,
            &boundary,
            &[],
            &AssemblyOptions {
                pretension: plate.pretension_n_m,
                support: if plate.clamped {
                    EdgeSupport::Clamped
                } else {
                    EdgeSupport::SimplySupported
                },
            },
        )
        .map_err(map_plate)?;
        let omega11 = ss_omega11(&section, plate.length_m, plate.width_m);
        let window = (
            (0.25 * omega11).powi(2),
            (omega11 * (plate.n_modes as f64 + 2.0) * 4.0).powi(2),
        );
        let projection = PlateProjection::new(&mesh, plate.length_m / nx as f64);
        let mut modes = Vec::new();
        let radiators = harvest_radiators_with_modes(
            &model,
            &projection,
            window,
            plate.n_modes,
            0.0,
            None,
            |phi| {
                poll(gate)?;
                modes.push(phi.to_vec());
                Ok(())
            },
        )?;
        poll(gate)?;
        if radiators.len() != plate.n_modes {
            return Err(invalid(
                "plate modal window does not contain the complete requested thermal basis",
            ));
        }
        if radiators
            .iter()
            .any(|r| r.omega < options.band_rad_s.0 || r.omega > options.band_rad_s.1)
        {
            return Err(invalid(
                "retained plate mode lies outside the thermal spectrum band",
            ));
        }
        // CompactBody retains phi^T M phi, even when it differs slightly from
        // one. The exact same mass normalization applies to every physical map.
        let scales: Vec<_> = radiators
            .iter()
            .map(|r| 1.0 / det::sqrt(r.mass_kg))
            .collect();
        if scales.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(invalid(
                "thermal plate modal mass normalization is unrepresentable",
            ));
        }
        let columns: Vec<_> = modes.iter().map(Vec::as_slice).collect();
        let branches = project_thermal_relaxation(
            &mesh,
            &model.dof_map,
            &columns,
            &scales,
            slab,
            &spectrum,
            options.max_memory_states,
            options.max_projection_terms,
            gate,
        )?;
        let omegas: Vec<_> = radiators.iter().map(|r| r.omega).collect();
        admit_memory_band(&omegas, &branches, options.band_rad_s)?;
        let force_shapes: Vec<_> = radiators
            .iter()
            .zip(&scales)
            .map(|(r, scale)| r.drive_participation * scale)
            .collect();
        let volume_shapes: Vec<_> = radiators
            .iter()
            .zip(&scales)
            .map(|(r, scale)| r.area_m2 * scale)
            .collect();
        if force_shapes
            .iter()
            .chain(&volume_shapes)
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "thermal plate force or volume projection is unrepresentable",
            ));
        }
        poll(gate)?;
        Ok(Self {
            radiators,
            branches,
            force_shapes,
            volume_shapes,
            spectrum,
        })
    }

    /// Original geometric radiators; thermal loss is carried by the separate
    /// full memory image and is deliberately absent from their scalar zeta.
    #[must_use]
    pub fn radiators(&self) -> &[CompactBody] {
        &self.radiators
    }

    /// Positive original-basis arms with interleaved [q0,0,q1,0,...] projections.
    /// Initial memory for a fresh prepared system is fully relaxed at rest.
    #[must_use]
    pub fn relaxation_branches(&self) -> &[RelaxationBranch] {
        &self.branches
    }

    #[must_use]
    pub fn spectrum(&self) -> &SlabPronySpectrum {
        &self.spectrum
    }

    /// Generalized force per newton on the original unit-total-force footprint.
    /// The work-conjugate velocity is this same map applied to modal velocity.
    #[must_use]
    pub fn force_shapes_per_sqrt_kg(&self) -> &[f64] {
        &self.force_shapes
    }

    /// Outward volume flow per original mass-normalized modal velocity [m²/sqrt(kg)].
    /// This is the reciprocal pressure-force map for an explicit acoustic load.
    #[must_use]
    pub fn volume_shapes_m2_per_sqrt_kg(&self) -> &[f64] {
        &self.volume_shapes
    }

    /// Fill the existing runtime's generalized-force buffer without allocation.
    /// A refused input leaves the entire buffer unchanged.
    pub fn fill_force(
        &self,
        total_force_n: f64,
        output: &mut [f64],
    ) -> Result<(), AcousticRealizeError> {
        if output.len() != self.force_shapes.len()
            || !total_force_n.is_finite()
            || self
                .force_shapes
                .iter()
                .any(|b| !(b * total_force_n).is_finite())
        {
            return Err(invalid(
                "thermal plate drive requires a finite force and complete representable projection",
            ));
        }
        for (out, &shape) in output.iter_mut().zip(&self.force_shapes) {
            *out = shape * total_force_n;
        }
        Ok(())
    }

    /// Prepare the existing exact held-force system at zero mechanical and
    /// thermal state. All thermal cross terms and stiffness dispersion join its
    /// single augmented operator and original work/heat accounting.
    ///
    /// The original instantaneous-stiffness Nyquist guard and both setup
    /// budgets remain active. Thermal relaxation poles themselves are not
    /// oscillatory frequencies and need no artificial dt/tau restriction.
    /// A new preparation does not mutate this image or another live system.
    ///
    /// Pressure transfers are zero until an explicit acoustic owner is attached.
    /// Use the retained volume map for that reciprocal connection. Finite-band
    /// constitutive accuracy does not imply that arbitrary drive is bandlimited.
    pub fn prepare_coupled(
        &self,
        sample_rate_hz: u32,
        modal_budget: ModalAcousticTimeBudget,
        coupling: ModalCouplingConfig,
        preparation: PreparedModalBudget,
        gate: &CancelGate,
    ) -> Result<CoupledModalSystem, AcousticRealizeError> {
        poll(gate)?;
        let model = ModalAcousticTimeModel::try_new(
            sample_rate_hz,
            self.radiators
                .iter()
                .map(|r| ModalAcousticMode {
                    angular_frequency_rad_s: r.omega,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                })
                .collect(),
            modal_budget,
        )
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let system =
            CoupledModalSystem::new_prepared(vec![model], Vec::new(), coupling, preparation, gate)
                .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        system
            .with_material_relaxation(self.branches.clone(), None, gate)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))
    }
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plate() -> ThinPlate {
        ThinPlate {
            length_m: 0.20,
            width_m: 0.19,
            thickness_m: 0.002,
            density_kg_m3: 5000.0,
            e1_pa: 70e9,
            e2_pa: 53e9,
            nu12: 0.25,
            g12_pa: 20e9,
            material_angle_rad: 0.31,
            damping_ratio: 0.0,
            thermoelastic: None,
            kelvin_voigt_bending: None,
            n_modes: 4,
            geometric_nonlinearity: false,
            pretension_n_m: 10_000.0,
            clamped: false,
        }
    }

    fn thermal() -> AnisotropicPlateThermal {
        AnisotropicPlateThermal {
            temperature_k: 300.0,
            engineering_expansion_per_k: [2e-5, 7e-5, 1e-5],
            specific_heat_j_kg_k: 600.0,
            thickness_conductivity_w_m_k: 100.0,
        }
    }

    fn options() -> ThermoelasticMemoryOptions {
        ThermoelasticMemoryOptions {
            band_rad_s: (0.0, 30_000.0),
            normalized_tolerance: 1e-4,
            max_thermal_terms: 16,
            max_memory_states: 64,
            max_projection_terms: 1_000_000,
        }
    }

    fn config() -> ModalCouplingConfig {
        ModalCouplingConfig {
            max_modes: 8,
            max_connections: 8,
            max_setup_terms: 1_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1.0,
            maximum_abs_pressure_pa: 1e5,
            maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-14,
            energy_relative_tolerance: 1e-8,
        }
    }

    fn prepared() -> PreparedModalBudget {
        PreparedModalBudget {
            max_modes: 8,
            max_dense_bytes: 256 * 1024 * 1024,
            max_multiply_terms: 20_000_000_000,
        }
    }

    #[test]
    fn g1_full_linear_thermal_image_keeps_original_modes_ports_and_cross_terms() {
        let gate = CancelGate::new_clock_free();
        let original = super::super::certified_radiators(plate()).unwrap();
        let image =
            ThermoelasticPlateImage::from_plate(plate(), thermal(), options(), &gate).unwrap();
        assert_eq!(
            image.radiators().len(),
            4,
            "linear image must retain every requested original mode"
        );
        for (i, (a, b)) in image.radiators().iter().zip(&original).enumerate() {
            assert_eq!(a.omega.to_bits(), b.omega.to_bits());
            assert_eq!(a.mass_kg.to_bits(), b.mass_kg.to_bits());
            assert_eq!(a.zeta, 0.0);
            let expected_force = b.drive_participation / det::sqrt(b.mass_kg);
            let expected_volume = b.area_m2 / det::sqrt(b.mass_kg);
            assert!(
                (image.force_shapes_per_sqrt_kg()[i] - expected_force).abs()
                    <= 4.0 * f64::EPSILON * expected_force.abs().max(f64::MIN_POSITIVE)
            );
            assert!(
                (image.volume_shapes_m2_per_sqrt_kg()[i] - expected_volume).abs()
                    <= 4.0 * f64::EPSILON * expected_volume.abs().max(f64::MIN_POSITIVE)
            );
        }
        let mut off_diagonal = 0.0_f64;
        let mut diagonal = 0.0_f64;
        for i in 0..4 {
            for j in 0..4 {
                let entry: f64 = image
                    .relaxation_branches()
                    .iter()
                    .map(|b| b.stiffness * b.projection[2 * i] * b.projection[2 * j])
                    .sum();
                if i == j {
                    diagonal = diagonal.max(entry.abs());
                } else {
                    off_diagonal = off_diagonal.max(entry.abs());
                }
            }
        }
        assert!(diagonal > 0.0 && off_diagonal > 1e-6 * diagonal);
        let mut force = vec![0.0; 4];
        image.fill_force(0.25, &mut force).unwrap();
        let velocity = [0.2, -0.1, 0.4, -0.3];
        let mechanical_power: f64 = force.iter().zip(velocity).map(|(f, v)| f * v).sum();
        let footprint_velocity: f64 = image
            .force_shapes_per_sqrt_kg()
            .iter()
            .zip(velocity)
            .map(|(b, v)| b * v)
            .sum();
        assert!((mechanical_power - 0.25 * footprint_velocity).abs() < 1e-15);
    }

    #[test]
    fn g1_prepared_thermal_plate_closes_work_heat_and_retains_memory() {
        let gate = CancelGate::new_clock_free();
        let image =
            ThermoelasticPlateImage::from_plate(plate(), thermal(), options(), &gate).unwrap();
        let mut system = image
            .prepare_coupled(
                96_000,
                ModalAcousticTimeBudget::audible_reference(),
                config(),
                prepared(),
                &gate,
            )
            .unwrap();
        assert_eq!(
            system.material_memory().len(),
            image.relaxation_branches().len()
        );
        let mut force = vec![0.0; 4];
        let (mut work, mut heat) = (0.0, 0.0);
        for sample in 0..256 {
            image
                .fill_force(if sample < 24 { 0.1 } else { 0.0 }, &mut force)
                .unwrap();
            let step = system.step_under_gate(&force, &gate).unwrap();
            assert!(step.material_dissipation_j >= 0.0);
            assert_eq!(step.component_dissipation_j, 0.0);
            work += step.external_work_j;
            heat += step.material_dissipation_j;
        }
        let stored = system.total_energy_j().unwrap();
        assert!(heat > 0.0 && stored > 0.0);
        assert!((work - stored - heat).abs() < 2e-11 * work.abs().max(stored) + 1e-15);
        assert!(system.material_memory().iter().any(|z| *z != 0.0));
        let before = stored;
        for _ in 0..128 {
            system.step(&[0.0; 4]).unwrap();
        }
        assert!(system.total_energy_j().unwrap() < before);
    }

    #[test]
    fn g0_duplicate_loss_and_failed_force_projection_do_not_create_silent_fallbacks() {
        let gate = CancelGate::new_clock_free();
        let mut invalid_plate = plate();
        invalid_plate.damping_ratio = 0.01;
        assert!(
            ThermoelasticPlateImage::from_plate(invalid_plate, thermal(), options(), &gate)
                .is_err()
        );
        invalid_plate = plate();
        invalid_plate.geometric_nonlinearity = true;
        assert!(
            ThermoelasticPlateImage::from_plate(invalid_plate, thermal(), options(), &gate)
                .is_err()
        );
        let image =
            ThermoelasticPlateImage::from_plate(plate(), thermal(), options(), &gate).unwrap();
        let mut output = [7.0; 4];
        assert!(image.fill_force(f64::NAN, &mut output).is_err());
        assert_eq!(output, [7.0; 4]);
        let mut absent = thermal();
        absent.engineering_expansion_per_k = [0.0; 3];
        let zero = ThermoelasticPlateImage::from_plate(plate(), absent, options(), &gate).unwrap();
        assert!(zero.relaxation_branches().is_empty());
        assert_eq!(zero.radiators().len(), 4);
    }
}
