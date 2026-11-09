//! Finite-length, two-ended opening ducts in the original coupled problem.
use super::{
    AcousticMedium, C64, CylindricalCavityNeck, MAX_CAVITY_OPENINGS, ModalRadiationLoading,
    VibroError, VibroacousticModel, invalid, poll,
};
use fs_duct::compose::TwoPort;
use fs_duct::molecular::{MolecularDuct, MolecularDuctError, MolecularWallResolution};
use fs_exec::Cx;
use fs_material::{gas::GasState, relaxation::RelaxationState};

mod physics;
pub(super) use physics::DistributedOpeningFrequency;

/// Independent physical and numerical bounds for the finite opening ducts.
#[derive(Clone, Copy, Debug)]
pub struct DistributedOpeningControls {
    /// omega*radius/sqrt(R_specific*T) at every original frequency.
    /// A positive ceiling <= 1 retains the explicitly plane-pressure model
    /// below the circular transverse acoustic regime. This is a reference
    /// gas-domain bound, not an error certificate for the viscothermal field.
    /// Unlike the compact owner there is NO restriction on axial k*length.
    pub maximum_isothermal_ka: f64,
    /// Native Redheffer feedback bound when an optional entrance is composed.
    pub maximum_composition_condition: f64,
    /// Relative original scattering-face and independent boundary-work checks.
    /// Signed roundoff is retained; failed checks do not repair a duct or heat.
    pub boundary_relative_tolerance: f64,
}

#[derive(Clone, Debug)]
struct DuctRow {
    omega: f64,
    blocks: Vec<TwoPort>,
}

/// Original finite cylinders and one complete exterior radiation operator.
///
/// The two separate physical flows are positive from cavity INTO the duct and
/// from duct OUT into the exterior. Native scattering retains compressibility,
/// finite-length resonances, phase and attenuation without a growing ABCD
/// matrix or inversion at a pressure/flow node. The exterior sees only outer
/// flow; the cavity weak form sees only inner flow.
///
/// Cavity volume and modal norms must exclude the physical duct volume. The
/// optional inner-end correction remains restricted to a single uniform cavity
/// pressure mode and reference radius ka <= 0.25; resolved cavity entrance
/// inertia is not added again. All
/// outer end reaction belongs to the original joint exterior matrix.
///
/// Equilibrium cylinders use native Bessel shear/thermal response. Molecular
/// cylinders retain the complete original RelaxationState and the explicitly
/// selected native wall-diffusion resolution, with no independent equilibrium
/// heat-capacity shunt or added empirical attenuation. This remains frozen,
/// linear, quiescent plane-pressure acoustics, not a nonlinear jet, mean-flow
/// duct, evolving gas state or transverse-mode realization.
#[derive(Clone, Debug)]
pub struct DistributedCavityOpeningRadiation {
    structural_modes: usize,
    pressure_projection: Vec<Vec<f64>>,
    ducts: Vec<CylindricalCavityNeck>,
    gas: GasState,
    molecular: Option<RelaxationState>,
    molecular_wall: Option<MolecularWallResolution>,
    exterior: ModalRadiationLoading,
    controls: DistributedOpeningControls,
    rows: Vec<DuctRow>,
}

enum Material<'a> {
    Equilibrium(&'a GasState),
    Molecular(RelaxationState, Option<MolecularWallResolution>),
}
impl Material<'_> {
    fn gas(&self) -> &GasState {
        match self {
            Self::Equilibrium(gas) => gas,
            Self::Molecular(state, _) => state.gas(),
        }
    }
}

fn duct_error(error: fs_duct::DuctError) -> VibroError {
    match error {
        fs_duct::DuctError::Cancelled { .. } => VibroError::Cancelled,
        _ => invalid("native distributed opening duct refused"),
    }
}
fn molecular_error(error: MolecularDuctError) -> VibroError {
    match error {
        MolecularDuctError::Duct(error) => duct_error(error),
        _ => invalid("original molecular opening material refused"),
    }
}

impl DistributedCavityOpeningRadiation {
    /// Prepare finite equilibrium-gas Bessel cylinders at every exterior sample.
    ///
    /// The cavity uses the same gas and must have no molecular attachment.
    /// Use the molecular constructor when compression has finite-rate internal
    /// reservoirs; an equilibrium duct cannot silently replace that source.
    pub fn from_cylindrical_ducts(
        structural_modes: usize,
        pressure_projection: Vec<Vec<f64>>,
        ducts: &[CylindricalCavityNeck],
        gas: &GasState,
        exterior: ModalRadiationLoading,
        controls: DistributedOpeningControls,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        Self::prepare(
            structural_modes,
            pressure_projection,
            ducts,
            Material::Equilibrium(gas),
            exterior,
            controls,
            cx,
        )
    }

    /// Prepare native coupled molecular/shear/thermal propagation ONCE per duct.
    ///
    /// Some(wall) retains the original diffusion poles, frequency band and
    /// sampled wall-error admission; None explicitly selects adiabatic slip
    /// walls. The native prepared cylinders are reused over the complete grid.
    /// Attachment binds this SAME molecular state to the cavity compressibility,
    /// or requires exact equality with an existing cavity law.
    #[allow(clippy::too_many_arguments)]
    pub fn from_molecular_cylindrical_ducts(
        structural_modes: usize,
        pressure_projection: Vec<Vec<f64>>,
        ducts: &[CylindricalCavityNeck],
        state: RelaxationState,
        wall: Option<MolecularWallResolution>,
        exterior: ModalRadiationLoading,
        controls: DistributedOpeningControls,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        Self::prepare(
            structural_modes,
            pressure_projection,
            ducts,
            Material::Molecular(state, wall),
            exterior,
            controls,
            cx,
        )
    }

    fn prepare(
        structural_modes: usize,
        pressure_projection: Vec<Vec<f64>>,
        ducts: &[CylindricalCavityNeck],
        material: Material<'_>,
        exterior: ModalRadiationLoading,
        controls: DistributedOpeningControls,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        poll(cx)?;
        let count = ducts.len();
        if structural_modes == 0
            || count == 0
            || count > MAX_CAVITY_OPENINGS
            || structural_modes.checked_add(count) != Some(exterior.mode_count())
            || pressure_projection.is_empty()
            || pressure_projection.iter().any(|row| row.len() != count)
        {
            return Err(VibroError::ShapeMismatch {
                what: "distributed openings require original structure, cavity map and 1..=64 mouths",
            });
        }
        if pressure_projection.iter().flatten().any(|x| !x.is_finite())
            || !controls.maximum_isothermal_ka.is_finite()
            || controls.maximum_isothermal_ka <= 0.0
            || controls.maximum_isothermal_ka > 1.0
            || !controls.maximum_composition_condition.is_finite()
            || controls.maximum_composition_condition < 1.0
            || !controls.boundary_relative_tolerance.is_finite()
            || controls.boundary_relative_tolerance <= 0.0
            || controls.boundary_relative_tolerance >= 1.0
        {
            return Err(invalid(
                "finite opening map and physical/numerical duct controls required",
            ));
        }
        if pressure_projection.len() != 1 && ducts.iter().any(|duct| duct.inner_end.is_some()) {
            return Err(invalid(
                "compact inner-end correction is exclusive to one uniform cavity mode",
            ));
        }
        let gas = *material.gas();
        AcousticMedium::try_from_gas(&gas)?;
        let isothermal_speed = fs_math::det::sqrt(gas.specific_gas_constant * gas.temperature);
        if !isothermal_speed.is_finite() || isothermal_speed <= 0.0 {
            return Err(invalid(
                "original gas isothermal reference speed is unrepresentable",
            ));
        }
        let mut inner_lengths = Vec::with_capacity(count);
        let mut prepared = Vec::with_capacity(count);
        for duct in ducts {
            poll(cx)?;
            if !duct.radius_m.is_finite()
                || duct.radius_m <= 0.0
                || !duct.length_m.is_finite()
                || duct.length_m <= 0.0
            {
                return Err(invalid(
                    "distributed cavity duct needs finite positive radius and length",
                ));
            }
            inner_lengths.push(match duct.inner_end {
                Some(flange) => fs_phs::neck_end_correction(duct.radius_m, flange)
                    .map_err(|_| invalid("native inner duct entrance correction refused"))?,
                None => 0.0,
            });
            let owner = match &material {
                Material::Equilibrium(_) => None,
                Material::Molecular(state, wall) => Some(
                    MolecularDuct::new_controlled(
                        &fs_duct::Duct {
                            segments: vec![fs_duct::Segment::Cylinder {
                                radius: duct.radius_m,
                                length: duct.length_m,
                            }],
                        },
                        state.clone(),
                        *wall,
                        fs_duct::ConeSlices::Fixed(1),
                        1,
                        &mut || cx.checkpoint().is_err(),
                    )
                    .map_err(molecular_error)?,
                ),
            };
            prepared.push(owner);
        }
        let mut rows = Vec::new();
        rows.try_reserve_exact(exterior.angular_frequencies_rad_s().len())
            .map_err(|_| invalid("distributed opening frequency table allocation failed"))?;
        for omega in exterior.angular_frequencies_rad_s() {
            poll(cx)?;
            if let Material::Molecular(_, Some(wall)) = &material
                && (omega < wall.band_rad_s[0] || omega > wall.band_rad_s[1])
            {
                return Err(invalid(
                    "distributed opening frequency lies outside original wall-diffusion band",
                ));
            }
            let mut blocks = Vec::with_capacity(count);
            for (j, duct) in ducts.iter().enumerate() {
                poll(cx)?;
                let ka = omega * (duct.radius_m / isothermal_speed);
                if !ka.is_finite() || ka > controls.maximum_isothermal_ka {
                    return Err(invalid(
                        "distributed opening exceeds its declared transverse acoustic domain",
                    ));
                }
                if inner_lengths[j] != 0.0 && ka > 0.25 {
                    return Err(invalid(
                        "compact inner entrance correction requires isothermal radius ka <= 0.25",
                    ));
                }
                let area = core::f64::consts::PI * duct.radius_m * duct.radius_m;
                let mut block = if let Some(owner) = &prepared[j] {
                    let chain = owner
                        .prepared_chain_controlled(omega, &mut || cx.checkpoint().is_err())
                        .map_err(molecular_error)?;
                    if chain.blocks.len() != 1 {
                        return Err(invalid(
                            "original cylindrical molecular duct changed its one-block geometry",
                        ));
                    }
                    chain.blocks[0]
                } else {
                    let wave = fs_duct::segment_wave(
                        &gas,
                        duct.radius_m,
                        omega,
                        fs_duct::LossModel::Bessel,
                    )
                    .map_err(duct_error)?;
                    TwoPort::from_uniform_wave(
                        &wave,
                        area,
                        duct.length_m,
                        gas.characteristic_impedance / area,
                    )
                    .map_err(duct_error)?
                };
                if inner_lengths[j] != 0.0 {
                    let entrance = C64::new(0.0, -omega * (gas.density * inner_lengths[j] / area));
                    let local = TwoPort::from_abcd(
                        [C64::ONE, entrance, C64::ZERO, C64::ONE],
                        block.reference_impedance,
                    )
                    .map_err(duct_error)?;
                    block = local
                        .star(block, controls.maximum_composition_condition)
                        .map_err(duct_error)?
                        .0;
                }
                blocks.push(block);
            }
            rows.push(DuctRow { omega, blocks });
        }
        poll(cx)?;
        let (molecular, molecular_wall) = match material {
            Material::Equilibrium(_) => (None, None),
            Material::Molecular(state, wall) => (Some(state), wall),
        };
        Ok(Self {
            structural_modes,
            pressure_projection,
            ducts: ducts.to_vec(),
            gas,
            molecular,
            molecular_wall,
            exterior,
            controls,
            rows,
        })
    }

    pub fn opening_count(&self) -> usize {
        self.ducts.len()
    }
    pub fn cavity_mode_count(&self) -> usize {
        self.pressure_projection.len()
    }
    pub const fn structural_mode_count(&self) -> usize {
        self.structural_modes
    }
    pub fn gas(&self) -> &GasState {
        &self.gas
    }
    pub fn ducts(&self) -> &[CylindricalCavityNeck] {
        &self.ducts
    }
    pub fn pressure_projection(&self) -> &[Vec<f64>] {
        &self.pressure_projection
    }
    pub fn exterior_loading(&self) -> &ModalRadiationLoading {
        &self.exterior
    }
    pub fn molecular_state(&self) -> Option<&RelaxationState> {
        self.molecular.as_ref()
    }
    pub const fn molecular_wall_resolution(&self) -> Option<MolecularWallResolution> {
        self.molecular_wall
    }
    pub const fn controls(&self) -> DistributedOpeningControls {
        self.controls
    }

    /// Original native two-ports at one exact declared frequency.
    pub fn two_ports(&self, omega: f64) -> Result<&[TwoPort], VibroError> {
        Ok(self.at_frequency(omega)?.ducts)
    }

    pub(in crate::vibroacoustic) fn at_frequency(
        &self,
        omega: f64,
    ) -> Result<DistributedOpeningFrequency<'_>, VibroError> {
        let exterior = self.exterior.at_omega(omega)?;
        let row = self
            .rows
            .binary_search_by(|row| row.omega.total_cmp(&omega))
            .map_err(|_| invalid("distributed opening has no original frequency sample"))?;
        Ok(DistributedOpeningFrequency {
            structural_modes: self.structural_modes,
            pressure_projection: &self.pressure_projection,
            exterior,
            ducts: &self.rows[row].blocks,
            tolerance: self.controls.boundary_relative_tolerance,
        })
    }
}

impl VibroacousticModel {
    /// Attach finite ducts and one original joint exterior to a sealed basis.
    ///
    /// Both original flows remain unknowns even at transmission or pressure
    /// nodes. The first cavity mode must be uniform with zero frequency.
    /// The molecular constructor binds its full original law to the cavity;
    /// a different existing law or scalar cavity loss refuses. Equilibrium
    /// ducts cannot be attached to a molecular cavity. Later material
    /// attachment must preserve that exact source selection.
    pub fn with_distributed_cavity_openings(
        mut self,
        openings: DistributedCavityOpeningRadiation,
    ) -> Result<Self, VibroError> {
        if openings.structural_modes != self.omegas_s.len()
            || openings.cavity_mode_count() != self.omegas_a.len()
        {
            return Err(VibroError::ShapeMismatch {
                what: "distributed openings must retain the original structural and pressure bases",
            });
        }
        if self.z_modal.is_some() || self.openings.is_some() || self.distributed_openings.is_some()
        {
            return Err(invalid(
                "one complete compact or distributed exterior attachment is required",
            ));
        }
        if self.omegas_a[0] != 0.0 {
            return Err(invalid(
                "distributed openings require the rigid-wall cavity's uniform zero mode",
            ));
        }
        let medium = AcousticMedium::try_from_gas(&openings.gas)?;
        if medium.rho0 != self.medium_a.rho0 || medium.c0 != self.medium_a.c0 {
            return Err(invalid(
                "distributed duct and cavity basis must retain the same gas",
            ));
        }
        match openings.molecular_state() {
            Some(state) => {
                if self
                    .cavity_relaxation
                    .as_ref()
                    .is_some_and(|existing| existing != state)
                {
                    return Err(invalid(
                        "distributed duct and cavity must retain the same molecular reservoirs",
                    ));
                }
                self = self.with_cavity_relaxation(state.clone())?;
            }
            None if self.cavity_relaxation.is_some() => {
                return Err(invalid(
                    "equilibrium distributed ducts cannot replace molecular cavity compression",
                ));
            }
            None => {}
        }
        self.distributed_openings = Some(openings);
        Ok(self)
    }

    pub fn distributed_cavity_openings(&self) -> Option<&DistributedCavityOpeningRadiation> {
        self.distributed_openings.as_ref()
    }
}

#[cfg(test)]
mod tests;
