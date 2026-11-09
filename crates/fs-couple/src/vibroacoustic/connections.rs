//! Finite internal ducts coupling original sealed acoustic compartments.
use super::openings::{DuctMaterial, DuctRow, prepare_cylindrical_rows};
use super::{
    AcousticMedium, C64, CavityModes, CylindricalCavityNeck, DistributedOpeningControls,
    VibroError, VibroacousticModel,
};
use fs_duct::compose::TwoPort;
use fs_duct::molecular::MolecularWallResolution;
use fs_exec::Cx;
use fs_material::{gas::GasState, relaxation::RelaxationState};

mod basis;
pub use basis::CavityCompartmentBasis;
mod physics;
pub use physics::InternalCavityDuctResponse;
pub(super) use physics::flow_delta;

fn invalid(what: &'static str) -> VibroError {
    VibroError::BadParameter { what }
}
fn poll(cx: &Cx<'_>) -> Result<(), VibroError> {
    cx.checkpoint().map_err(|_| VibroError::Cancelled)
}

/// One straight finite duct between two DIFFERENT physical compartments.
///
/// Each face projection is the area mean of the ORIGINAL local pressure mode
/// against the unit-volume-flow profile, in that compartment's local mode
/// order. Geometry and modal norms exclude this duct's physical volume.
/// Both end flows remain unknowns. Positive axial flow is left to right;
/// the two flows generally differ because the finite duct stores compression.
#[derive(Debug, Clone)]
pub struct CavityDuctConnection {
    pub left_compartment: usize,
    pub right_compartment: usize,
    pub left_pressure_projection: Vec<f64>,
    pub right_pressure_projection: Vec<f64>,
    pub radius_m: f64,
    pub length_m: f64,
    /// Optional unresolved entrance inertia, ordered [left,right]. Each end
    /// requires ONE uniform mode in its own compartment and reference ka<=0.25.
    /// Resolved entrance fields use None, independently at each face.
    pub end_corrections: [Option<fs_phs::MouthFlange>; 2],
}

/// Original finite duct graph in the same complete structural/cavity solve.
///
/// Compatible with the original compact/distributed joint exterior owner or
/// a sealed exterior. Native scattering preserves internal duct resonances,
/// compressibility and material loss. All compartments and ducts share one
/// frozen original gas and molecular selection; different temperatures,
/// compositions, mean flows and nonlinear aperture exchange are not inferred.
#[derive(Debug, Clone)]
pub struct InternalCavityDucts {
    basis: CavityCompartmentBasis,
    connections: Vec<CavityDuctConnection>,
    gas: GasState,
    molecular: Option<RelaxationState>,
    molecular_wall: Option<MolecularWallResolution>,
    controls: DistributedOpeningControls,
    rows: Vec<DuctRow>,
    left_projection: Vec<Vec<f64>>,
    right_projection: Vec<Vec<f64>>,
}
impl InternalCavityDucts {
    pub fn from_cylindrical_ducts(
        basis: CavityCompartmentBasis,
        connections: Vec<CavityDuctConnection>,
        gas: &GasState,
        frequencies_rad_s: &[f64],
        controls: DistributedOpeningControls,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        Self::prepare(
            basis,
            connections,
            DuctMaterial::Equilibrium(gas),
            frequencies_rad_s,
            controls,
            cx,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn from_molecular_cylindrical_ducts(
        basis: CavityCompartmentBasis,
        connections: Vec<CavityDuctConnection>,
        state: RelaxationState,
        wall: Option<MolecularWallResolution>,
        frequencies_rad_s: &[f64],
        controls: DistributedOpeningControls,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        Self::prepare(
            basis,
            connections,
            DuctMaterial::Molecular(state, wall),
            frequencies_rad_s,
            controls,
            cx,
        )
    }
    fn prepare(
        basis: CavityCompartmentBasis,
        connections: Vec<CavityDuctConnection>,
        material: DuctMaterial<'_>,
        frequencies: &[f64],
        controls: DistributedOpeningControls,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        poll(cx)?;
        let count = connections.len();
        if count == 0 || count > super::openings::MAX_CAVITY_OPENINGS {
            return Err(invalid(
                "internal cavity graph requires 1..=64 finite ducts",
            ));
        }
        let gas = *material.gas();
        let medium = AcousticMedium::try_from_gas(&gas)?;
        if basis.modes.rho0 != medium.rho0 || basis.modes.c0 != medium.c0 {
            return Err(invalid(
                "internal ducts and every original pressure basis must share one gas",
            ));
        }
        let na = basis.modes.omegas.len();
        let mut left = vec![vec![0.0; count]; na];
        let mut right = vec![vec![0.0; count]; na];
        let mut cylinders = Vec::with_capacity(count);
        for (j, link) in connections.iter().enumerate() {
            poll(cx)?;
            if link.left_compartment == link.right_compartment
                || link.left_compartment >= basis.compartment_count()
                || link.right_compartment >= basis.compartment_count()
            {
                return Err(invalid(
                    "each internal duct must join two declared distinct compartments",
                ));
            }
            for (end, (compartment, projection, target)) in [
                (
                    link.left_compartment,
                    &link.left_pressure_projection,
                    &mut left,
                ),
                (
                    link.right_compartment,
                    &link.right_pressure_projection,
                    &mut right,
                ),
            ]
            .into_iter()
            .enumerate()
            {
                let map = &basis.mode_indices[compartment];
                if projection.len() != map.len() || projection.iter().any(|x| !x.is_finite()) {
                    return Err(invalid(
                        "duct face projection must retain its complete original compartment basis",
                    ));
                }
                let uniform =
                    basis.modes.interface[map[0]][basis.interface_ranges()[compartment].start];
                if (projection[0] - uniform).abs() > 64.0 * f64::EPSILON * uniform.abs() {
                    return Err(invalid(
                        "unit-volume-flow face must retain the original uniform pressure projection",
                    ));
                }
                if link.end_corrections[end].is_some() && map.len() != 1 {
                    return Err(invalid(
                        "unresolved duct entrance correction requires one uniform compartment mode",
                    ));
                }
                for (&q, &h) in map.iter().zip(projection) {
                    target[q][j] = h;
                }
            }
            cylinders.push(CylindricalCavityNeck {
                radius_m: link.radius_m,
                length_m: link.length_m,
                inner_end: link.end_corrections[0],
            });
        }
        // SAME native material/geometry owner as finite exterior opening ducts.
        let mut rows = prepare_cylindrical_rows(&cylinders, &material, frequencies, controls, cx)?;
        let isothermal = fs_math::det::sqrt(gas.specific_gas_constant * gas.temperature);
        for row in &mut rows {
            for (j, link) in connections.iter().enumerate() {
                poll(cx)?;
                if let Some(flange) = link.end_corrections[1] {
                    if row.omega * (link.radius_m / isothermal) > 0.25 {
                        return Err(invalid(
                            "compact right entrance correction requires isothermal radius ka <= 0.25",
                        ));
                    }
                    let delta = fs_phs::neck_end_correction(link.radius_m, flange)
                        .map_err(|_| invalid("native right duct entrance correction refused"))?;
                    let area = core::f64::consts::PI * link.radius_m * link.radius_m;
                    let impedance = C64::new(0.0, -row.omega * (gas.density * delta / area));
                    let entrance = TwoPort::from_abcd(
                        [C64::ONE, impedance, C64::ZERO, C64::ONE],
                        row.blocks[j].reference_impedance,
                    )
                    .map_err(|_| invalid("right duct entrance is unrepresentable"))?;
                    row.blocks[j] = row.blocks[j]
                        .star(entrance, controls.maximum_composition_condition)
                        .map_err(|_| {
                            invalid("right duct entrance feedback exceeds its original bound")
                        })?
                        .0;
                }
            }
        }
        let (molecular, molecular_wall) = match material {
            DuctMaterial::Equilibrium(_) => (None, None),
            DuctMaterial::Molecular(state, wall) => (Some(state), wall),
        };
        poll(cx)?;
        Ok(Self {
            basis,
            connections,
            gas,
            molecular,
            molecular_wall,
            controls,
            rows,
            left_projection: left,
            right_projection: right,
        })
    }
    pub fn basis(&self) -> &CavityCompartmentBasis {
        &self.basis
    }
    pub fn connections(&self) -> &[CavityDuctConnection] {
        &self.connections
    }
    pub fn gas(&self) -> &GasState {
        &self.gas
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
    pub fn angular_frequencies_rad_s(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.rows.iter().map(|row| row.omega)
    }
    pub fn two_ports(&self, omega: f64) -> Result<&[TwoPort], VibroError> {
        let at = self
            .rows
            .binary_search_by(|row| row.omega.total_cmp(&omega))
            .map_err(|_| {
                invalid("internal cavity ducts have no original sample at this frequency")
            })?;
        Ok(&self.rows[at].blocks)
    }
    pub(super) fn at_frequency(
        &self,
        omega: f64,
    ) -> Result<physics::InternalDuctFrequency<'_>, VibroError> {
        Ok(physics::InternalDuctFrequency {
            left_projection: &self.left_projection,
            right_projection: &self.right_projection,
            connections: &self.connections,
            compartments: self.basis.compartment_count(),
            ducts: self.two_ports(omega)?,
            tolerance: self.controls.boundary_relative_tolerance,
        })
    }
}

impl VibroacousticModel {
    /// Add internal finite ducts while retaining the ORIGINAL full modal solve
    /// and any complete joint exterior load. Basis geometry, mode normalization
    /// and coupling must match the declared compartment basis. The caller owns
    /// exact trace/basis identity: frequency/norm checks cannot identify a mode
    /// shape or prove that the supplied structural coupling uses that shape.
    /// No uniform
    /// compression coordinate may be removed by a subsequent truncation.
    pub fn with_internal_cavity_ducts(
        mut self,
        ducts: InternalCavityDucts,
    ) -> Result<Self, VibroError> {
        if self.internal_ducts.is_some() {
            return Err(invalid("one complete internal cavity graph is required"));
        }
        let basis = ducts.basis.modes();
        if self.omegas_a != basis.omegas
            || self.lambdas != basis.lambdas
            || self.eta_a != basis.loss_factor
            || self.medium_a.rho0 != basis.rho0
            || self.medium_a.c0 != basis.c0
        {
            return Err(invalid(
                "internal graph must retain the original combined pressure basis and gas",
            ));
        }
        if self
            .openings
            .as_ref()
            .is_some_and(|o| o.gas() != ducts.gas())
            || self
                .distributed_openings
                .as_ref()
                .is_some_and(|o| o.gas() != ducts.gas())
        {
            return Err(invalid(
                "internal and exterior opening ducts must retain the same original gas",
            ));
        }
        match ducts.molecular_state() {
            Some(state) => {
                if self
                    .cavity_relaxation
                    .as_ref()
                    .is_some_and(|original| original != state)
                {
                    return Err(invalid(
                        "internal ducts and cavity must retain the same molecular reservoirs",
                    ));
                }
                self = self.with_cavity_relaxation(state.clone())?;
            }
            None if self.cavity_relaxation.is_some() => {
                return Err(invalid(
                    "equilibrium internal ducts cannot replace molecular cavity compression",
                ));
            }
            None => {}
        }
        self.internal_ducts = Some(ducts);
        Ok(self)
    }
    pub fn internal_cavity_ducts(&self) -> Option<&InternalCavityDucts> {
        self.internal_ducts.as_ref()
    }
}

#[cfg(test)]
mod tests;
