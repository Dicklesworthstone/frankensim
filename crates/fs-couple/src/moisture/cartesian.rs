//! Conservative three-dimensional bound-water transport in a homogeneous solid.
//!
//! Uniform Cartesian reference cells retain one dry skeleton. Three explicitly
//! selected principal diffusivities are aligned with the actual grid axes; a
//! repeated source key is an explicit equal-diffusivity assignment, never an
//! inferred tensor. Internal faces use harmonic series resistance by default.
//! Opt-in Kirchhoff transport integrates each original principal D_i(M) across
//! the face, allowing a degenerate dry state to wet without an artificial
//! impermeable barrier. All six external faces use the original sorption and
//! surface-transfer owner, with the same nonlinear half-cell potential balance.
//! A moisture-independent source may explicitly assign zero to a principal
//! diffusivity: that direction is impermeable, including its external faces.
//!
//! Backward Euler resolves local source diffusivities at its accepted endpoint
//! using bounded Picard iteration and the native checked-CG solver. The compact
//! seven-point operator stores O(N) coefficients and has O(N) application cost.
//! This is isothermal, fixed-reference dilute bound moisture. Unlike-material
//! partitioning, rotated tensors, sorption heat and stress-assisted diffusion
//! are not supplied by this spatial owner. The original `atmosphere` and
//! `solution` controllers couple these same cells to finite external inventories.

use super::*;

/// Actual box geometry and an authored frame for its x/y/z principal directions.
#[derive(Clone, Copy, Debug)]
pub struct MoistureCartesianGeometry {
    pub lengths_m: [f64; 3],
    /// Equal cells per axis. Storage is x-fastest: x + nx*(y + ny*z).
    pub cells: [usize; 3],
    /// Nonzero identity of the grid frame; also used by the native solid mesh.
    pub frame: ContentHash,
    pub maximum_moisture: f64,
    pub maximum_cells: usize,
}

/// Accepted equal-volume spatial moisture field; cloning is a full checkpoint.
#[derive(Clone, Debug)]
pub struct MoistureCartesian {
    card: MaterialCard,
    reference: ResolvedMaterialStatePoint,
    transport: ResolvedMaterialStatePoint,
    properties: [String; 3],
    geometry: MoistureCartesianGeometry,
    dry_mass_per_cell_kg: f64,
    moisture: Vec<f64>,
    time_s: f64,
    identity: ContentHash,
    kirchhoff_diffusivity: bool,
}

/// Physical transfers and original source/solver evidence at one spatial endpoint.
#[derive(Clone, Debug)]
pub struct MoistureCartesianReport {
    pub interval_s: [f64; 2],
    /// Signed received water [kg], in -x,+x,-y,+y,-z,+z order.
    pub boundary_water_kg: [f64; 6],
    /// Actual interface content, (adjacent cell index, dry-basis moisture), in
    /// ascending cell order on each face. Sealed faces have an empty vector.
    pub surface_moisture: [Vec<(usize, f64)>; 6],
    pub reservoir_moisture: [Option<f64>; 6],
    pub surface_sorption: [Option<ResolvedMaterialStatePoint>; 6],
    pub surface_transfer: [Option<ResolvedMaterialStatePoint>; 6],
    pub stored_water_kg: [f64; 2],
    pub water_balance_residual_kg: f64,
    pub transport_start: Vec<ResolvedMaterialStatePoint>,
    pub transport_end: Vec<ResolvedMaterialStatePoint>,
    /// Last frozen-coefficient solve, independently checked by its native owner.
    /// The load and unknown use moisture relative to the initial first cell.
    pub linear: SolveReport,
    pub nonlinear_iterations: usize,
    /// Sum of native CG iterations across every Picard solve.
    pub linear_iterations: usize,
    /// True relative Euclidean residual of the endpoint-dependent equations,
    /// with the same shifted-load normalization as the linear solve.
    pub nonlinear_relative_residual: f64,
    /// Original prepared x/y/z source-chart identities for Kirchhoff transport.
    /// `None` denotes the ordinary harmonic-coefficient scheme.
    pub kirchhoff_sources: Option<[ContentHash; 3]>,
}

impl MoistureCartesian {
    /// Bind one material's source-selected diagonal diffusion law to real 3D cells.
    /// Property order is x/y/z; the caller asserts these are principal directions
    /// in `geometry.frame`. Scalar data never establish arbitrary tensor rotation.
    /// Density is evaluated only at the typed reference moisture, exactly as in
    /// `MoistureSlab`; subsequent uptake changes wet mass, never dry skeleton mass.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        card: MaterialCard,
        reference_point: QueryPoint,
        diffusivity_properties: [&str; 3],
        selection: MaterialPropertySelection,
        geometry: MoistureCartesianGeometry,
        initial_moisture: Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        Self::new_configured(
            card,
            reference_point,
            diffusivity_properties,
            selection,
            geometry,
            initial_moisture,
            false,
            cx,
        )
    }

    /// Bind source-integrated homogeneous diffusion in each principal direction.
    /// A nonnegative D_i(M) may vanish at the initial dry state. Every advance
    /// admits its complete moisture interval through the original material chart.
    /// Nonlinear accuracy and work limits remain explicit arguments to `advance`.
    #[allow(clippy::too_many_arguments)]
    pub fn new_kirchhoff(
        card: MaterialCard,
        reference_point: QueryPoint,
        diffusivity_properties: [&str; 3],
        selection: MaterialPropertySelection,
        geometry: MoistureCartesianGeometry,
        initial_moisture: Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        Self::new_configured(
            card,
            reference_point,
            diffusivity_properties,
            selection,
            geometry,
            initial_moisture,
            true,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_configured(
        card: MaterialCard,
        reference_point: QueryPoint,
        diffusivity_properties: [&str; 3],
        selection: MaterialPropertySelection,
        geometry: MoistureCartesianGeometry,
        initial_moisture: Vec<f64>,
        kirchhoff_diffusivity: bool,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let count = geometry
            .cells
            .iter()
            .try_fold(1usize, |n, &m| n.checked_mul(m))
            .ok_or(MoistureError::Invalid(
                "Cartesian moisture cell count overflow",
            ))?;
        if count == 0
            || count != initial_moisture.len()
            || count > geometry.maximum_cells
            || count > MAX_MOISTURE_CELLS
            || geometry.lengths_m.iter().any(|&x| !positive(x))
            || geometry.frame == ContentHash([0; 32])
            || !positive(geometry.maximum_moisture)
            || initial_moisture
                .iter()
                .any(|m| !m.is_finite() || *m < 0.0 || *m > geometry.maximum_moisture)
        {
            return Err(MoistureError::Invalid(
                "bounded Cartesian cells, positive lengths, identified axes and supported moisture required",
            ));
        }
        let (reference, transport, dry_density) = binding::bind_source(
            &card,
            &reference_point,
            &diffusivity_properties,
            selection,
            geometry.maximum_moisture,
            kirchhoff_diffusivity,
        )?;
        let spacing =
            std::array::from_fn::<_, 3, _>(|k| geometry.lengths_m[k] / geometry.cells[k] as f64);
        let dry_mass_per_cell_kg = dry_density * spacing[0] * spacing[1] * spacing[2];
        if spacing.iter().any(|&x| !positive(x))
            || !positive(dry_mass_per_cell_kg)
            || !positive(dry_mass_per_cell_kg * count as f64)
        {
            return Err(MoistureError::Invalid(
                "Cartesian dry inventory is unrepresentable",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.moisture-cartesian.v1");
        hash.update(reference.identity().as_bytes());
        hash.update(geometry.frame.as_bytes());
        for k in 0..3 {
            hash.update(&geometry.lengths_m[k].to_bits().to_le_bytes());
            hash.update(&(geometry.cells[k] as u64).to_le_bytes());
            hash.update(&(diffusivity_properties[k].len() as u64).to_le_bytes());
            hash.update(diffusivity_properties[k].as_bytes());
        }
        hash.update(&geometry.maximum_moisture.to_bits().to_le_bytes());
        hash.update(&dry_mass_per_cell_kg.to_bits().to_le_bytes());
        let field = Self {
            card,
            reference,
            transport,
            properties: diffusivity_properties.map(str::to_owned),
            geometry,
            dry_mass_per_cell_kg,
            moisture: initial_moisture,
            time_s: 0.0,
            identity: hash.finalize(),
            kirchhoff_diffusivity,
        };
        if !positive(field.wet_mass_kg()) || !field.water_mass_kg().is_finite() {
            return Err(MoistureError::Invalid(
                "Cartesian wet inventory is unrepresentable",
            ));
        }
        field.transport_states(&field.moisture, cx)?;
        poll(cx)?;
        Ok(field)
    }

    /// Select exact source-potential transport for this field and all existing
    /// reservoir/solid consumers. Geometry, material identity and accepted state
    /// stay the same; each report identifies the charts used by its interval.
    pub fn with_kirchhoff_diffusivity(mut self) -> Self {
        self.kirchhoff_diffusivity = true;
        self
    }

    pub const fn geometry(&self) -> MoistureCartesianGeometry {
        self.geometry
    }
    pub const fn definition_identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn card(&self) -> &MaterialCard {
        &self.card
    }
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }
    pub fn moisture(&self) -> &[f64] {
        &self.moisture
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn dry_mass_kg(&self) -> f64 {
        self.dry_mass_per_cell_kg * self.moisture.len() as f64
    }
    pub fn water_mass_kg(&self) -> f64 {
        self.moisture
            .iter()
            .map(|m| m * self.dry_mass_per_cell_kg)
            .sum()
    }
    pub fn wet_mass_kg(&self) -> f64 {
        self.dry_mass_kg() + self.water_mass_kg()
    }
    /// Native cell index, rejecting any coordinate outside the actual box.
    pub fn cell_index(&self, coordinate: [usize; 3]) -> Option<usize> {
        if coordinate
            .iter()
            .zip(self.geometry.cells)
            .any(|(&i, n)| i >= n)
        {
            return None;
        }
        Some(
            coordinate[0]
                + self.geometry.cells[0] * (coordinate[1] + self.geometry.cells[1] * coordinate[2]),
        )
    }

    pub(super) fn field(&self) -> binding::MoistureField<'_> {
        binding::MoistureField {
            card: &self.card,
            reference: &self.reference,
            moisture: &self.moisture,
            dry_mass_per_cell_kg: self.dry_mass_per_cell_kg,
            time_s: self.time_s,
            identity: self.identity,
        }
    }
    fn transport_states(
        &self,
        moisture: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<ResolvedMaterialStatePoint>, MoistureError> {
        moisture
            .iter()
            .map(|&m| {
                poll(cx)?;
                self.transport_state_at(m)
            })
            .collect()
    }

    fn transport_state_at(
        &self,
        moisture: f64,
    ) -> Result<ResolvedMaterialStatePoint, MoistureError> {
        let point = self.field().point_at(moisture)?;
        if self.kirchhoff_diffusivity {
            Ok(self.transport.at_query_point(&self.card, &point)?)
        } else {
            binding::resolve_transport_state(&self.card, &self.transport, &point)
        }
    }

    fn kirchhoff_profiles(
        &self,
        reservoir: [Option<f64>; 6],
        cx: &Cx<'_>,
    ) -> Result<[KirchhoffMoistureDiffusivity; 3], MoistureError> {
        let mut lower = f64::INFINITY;
        let mut upper: f64 = 0.0;
        for (i, m) in self
            .moisture
            .iter()
            .copied()
            .chain(reservoir.into_iter().flatten())
            .enumerate()
        {
            if i % 128 == 0 {
                poll(cx)?;
            }
            lower = lower.min(m);
            upper = upper.max(m);
        }
        let prepare = |axis: usize| {
            poll(cx)?;
            Ok::<_, MoistureError>(KirchhoffMoistureDiffusivity::from_binding(
                &self.card,
                &self.transport,
                &self.properties[axis],
                [lower, upper],
            )?)
        };
        let profiles = [prepare(0)?, prepare(1)?, prepare(2)?];
        poll(cx)?;
        Ok(profiles)
    }

    /// Propose a fully implicit interval on all three spatial directions.
    /// Faces are constant over this interval, in -x,+x,-y,+y,-z,+z order.
    /// The original maximum-content-change and water-balance gates apply;
    /// neither spatial nor temporal discretization error is inferred from CG.
    pub fn advance(
        &self,
        dt_s: f64,
        faces: &[MoistureFace; 6],
        config: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureCartesianReport), MoistureError> {
        poll(cx)?;
        let temperature = self.reference.query().axes()["temperature"];
        let boundaries = [
            faces[0].resolve(temperature)?,
            faces[1].resolve(temperature)?,
            faces[2].resolve(temperature)?,
            faces[3].resolve(temperature)?,
            faces[4].resolve(temperature)?,
            faces[5].resolve(temperature)?,
        ];
        self.advance_resolved(dt_s, boundaries, config, nonlinear, cx)
    }

    pub(super) fn advance_resolved(
        &self,
        dt_s: f64,
        boundaries: [ResolvedMoistureFace; 6],
        config: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureCartesianReport), MoistureError> {
        poll(cx)?;
        nonlinear.validate()?;
        let end = self.time_s + dt_s;
        if !positive(dt_s)
            || !end.is_finite()
            || end <= self.time_s
            || !positive(config.water_balance_tolerance_kg)
            || !positive(config.maximum_content_increment)
        {
            return Err(MoistureError::Invalid(
                "representable forward spatial interval and positive admission limits required",
            ));
        }
        let reservoir = std::array::from_fn(|k| boundaries[k].equilibrium);
        for m in reservoir.into_iter().flatten() {
            poll(cx)?;
            if m > self.geometry.maximum_moisture {
                return Err(MoistureError::Invalid(
                    "Cartesian surface exceeds bound-water applicability",
                ));
            }
            self.transport_state_at(m)?;
        }
        let profiles = if self.kirchhoff_diffusivity {
            Some(self.kirchhoff_profiles(reservoir, cx)?)
        } else {
            None
        };
        let start = self.transport_states(&self.moisture, cx)?;
        let mut system = self.system(dt_s, &boundaries, &start, profiles.as_ref(), cx)?;
        let mut iterate = self.moisture.clone();
        let mut linear_iterations = 0usize;
        for iteration in 1..=nonlinear.maximum_iterations {
            poll(cx)?;
            let solution = checked_cg(
                &system.operator,
                &system.operator,
                &system.rhs,
                config.linear,
                |_| poll(cx),
            )
            .map_err(|e| match e {
                CheckedCgError::Interrupted(error) => error,
                CheckedCgError::InvalidInput(reason) => MoistureError::Invalid(reason),
            })?;
            if !solution.report.converged_euclidean() {
                return Err(MoistureError::NotConverged(solution.report));
            }
            linear_iterations = linear_iterations.checked_add(solution.report.iters).ok_or(
                MoistureError::Invalid("Cartesian solver work count overflow"),
            )?;
            let candidate: Vec<_> = solution.x.iter().map(|m| m + self.moisture[0]).collect();
            let endpoint = self.transport_states(&candidate, cx)?;
            let endpoint_system =
                self.system(dt_s, &boundaries, &endpoint, profiles.as_ref(), cx)?;
            let relative = endpoint_system.relative_residual(&solution.x, cx)?;
            let transfer = self.transfers(&candidate, &boundaries, &endpoint_system);
            let water_after = candidate
                .iter()
                .map(|m| m * self.dry_mass_per_cell_kg)
                .sum::<f64>();
            let balance = water_after - self.water_mass_kg() - transfer.iter().sum::<f64>();
            if relative < nonlinear.relative_tolerance
                && balance.is_finite()
                && balance.abs() <= config.water_balance_tolerance_kg
            {
                let lower = self
                    .moisture
                    .iter()
                    .copied()
                    .chain(reservoir.into_iter().flatten())
                    .fold(f64::INFINITY, f64::min);
                let upper = self
                    .moisture
                    .iter()
                    .copied()
                    .chain(reservoir.into_iter().flatten())
                    .fold(0.0, f64::max);
                let roundoff = 128.0 * f64::EPSILON * upper.max(1.0);
                for (cell, (&old, &new)) in self.moisture.iter().zip(&candidate).enumerate() {
                    if !new.is_finite()
                        || new < 0.0
                        || new > self.geometry.maximum_moisture
                        || new < lower - roundoff
                        || new > upper + roundoff
                    {
                        return Err(MoistureError::Invalid(
                            "Cartesian candidate violates maximum principle or domain",
                        ));
                    }
                    let increment = (new - old).abs();
                    if increment > config.maximum_content_increment {
                        return Err(MoistureError::ContentIncrement {
                            cell,
                            increment,
                            maximum: config.maximum_content_increment,
                        });
                    }
                }
                if !positive(self.dry_mass_kg() + water_after) {
                    return Err(MoistureError::Invalid(
                        "Cartesian endpoint wet inventory is unrepresentable",
                    ));
                }
                let mut surface_moisture: [Vec<(usize, f64)>; 6] =
                    std::array::from_fn(|_| Vec::new());
                for face in &endpoint_system.faces {
                    poll(cx)?;
                    let equilibrium = reservoir[face.side].expect("active boundary");
                    let m = if face.fraction == 1.0 {
                        equilibrium
                    } else {
                        face.fraction
                            .mul_add(equilibrium - candidate[face.cell], candidate[face.cell])
                    };
                    self.transport_state_at(m)?;
                    surface_moisture[face.side].push((face.cell, m));
                }
                let report = MoistureCartesianReport {
                    interval_s: [self.time_s, end],
                    boundary_water_kg: transfer,
                    surface_moisture,
                    reservoir_moisture: reservoir,
                    surface_sorption: std::array::from_fn(|k| boundaries[k].sorption.clone()),
                    surface_transfer: std::array::from_fn(|k| boundaries[k].transfer.clone()),
                    stored_water_kg: [self.water_mass_kg(), water_after],
                    water_balance_residual_kg: balance,
                    transport_start: start,
                    transport_end: endpoint,
                    linear: solution.report,
                    nonlinear_iterations: iteration,
                    linear_iterations,
                    nonlinear_relative_residual: relative,
                    kirchhoff_sources: profiles
                        .as_ref()
                        .map(|p| std::array::from_fn(|k| p[k].identity())),
                };
                let mut next = self.clone();
                next.moisture = candidate;
                next.time_s = end;
                poll(cx)?;
                return Ok((next, report));
            }
            if iteration == nonlinear.maximum_iterations {
                return Err(MoistureError::SpatialNotConverged {
                    iterations: iteration,
                    linear_iterations,
                    relative_residual: relative,
                    water_balance_residual_kg: balance,
                });
            }
            if nonlinear.relaxation == 1.0 {
                iterate = candidate;
                system = endpoint_system;
            } else {
                for (old, candidate) in iterate.iter_mut().zip(&candidate) {
                    *old = nonlinear.relaxation.mul_add(candidate - *old, *old);
                }
                let states = self.transport_states(&iterate, cx)?;
                system = self.system(dt_s, &boundaries, &states, profiles.as_ref(), cx)?;
            }
        }
        unreachable!("validated positive spatial iteration budget")
    }

    fn system(
        &self,
        dt: f64,
        faces: &[ResolvedMoistureFace; 6],
        states: &[ResolvedMaterialStatePoint],
        profiles: Option<&[KirchhoffMoistureDiffusivity; 3]>,
        cx: &Cx<'_>,
    ) -> Result<CartesianSystem, MoistureError> {
        let n = self.moisture.len();
        let cells = self.geometry.cells;
        let strides = [1, cells[0], cells[0] * cells[1]];
        let spacing: [f64; 3] =
            std::array::from_fn(|k| self.geometry.lengths_m[k] / cells[k] as f64);
        let scale = spacing.map(|dx| dt / dx / dx);
        if scale.iter().any(|&a| !positive(a)) {
            return Err(MoistureError::Invalid(
                "Cartesian time/length scale is unrepresentable",
            ));
        }
        let mut system = CartesianSystem {
            operator: CartesianOperator {
                diagonal: vec![1.0; n],
                storage_boundary: vec![1.0; n],
                forward: vec![[0.0; 3]; n],
                strides,
            },
            // Solve relative to one common reference, preserving the scale of
            // small gradients independently of the uniform water inventory.
            rhs: self.moisture.iter().map(|m| m - self.moisture[0]).collect(),
            faces: Vec::new(),
        };
        for i in 0..n {
            poll(cx)?;
            for k in 0..3 {
                let coordinate = (i / strides[k]) % cells[k];
                let d = states[i]
                    .property(&self.properties[k])
                    .expect("pinned principal diffusivity")
                    .value_si();
                if coordinate + 1 < cells[k] {
                    let j = i + strides[k];
                    let a = if let Some(profiles) = profiles {
                        let mean = profiles[k].mean_diffusivity(
                            states[i]
                                .state_coordinate("moisture_content")
                                .expect("owned moisture coordinate"),
                            states[j]
                                .state_coordinate("moisture_content")
                                .expect("owned moisture coordinate"),
                        )?;
                        kirchhoff::scaled_diffusivity(scale[k], mean)?
                    } else {
                        let other = states[j]
                            .property(&self.properties[k])
                            .expect("pinned principal diffusivity")
                            .value_si();
                        interior_coupling(scale[k], d, other)?
                    };
                    system.operator.forward[i][k] = a;
                    system.operator.diagonal[i] += a;
                    system.operator.diagonal[j] += a;
                }
                for (side, active) in [
                    (2 * k, coordinate == 0),
                    (2 * k + 1, coordinate + 1 == cells[k]),
                ] {
                    if !active {
                        continue;
                    }
                    if let Some(equilibrium) = faces[side].equilibrium {
                        let (coefficient, fraction) = if let Some(profiles) = profiles {
                            kirchhoff::surface_coupling(
                                &profiles[k],
                                &faces[side],
                                dt,
                                spacing[k],
                                states[i]
                                    .state_coordinate("moisture_content")
                                    .expect("owned moisture coordinate"),
                                scale[k],
                                cx,
                            )?
                        } else {
                            faces[side].coupling(dt, spacing[k], d, scale[k])?
                        };
                        if !coefficient.is_finite()
                            || (d > 0.0
                                && faces[side].coefficient_m_s != Some(0.0)
                                && !positive(coefficient))
                        {
                            return Err(MoistureError::Invalid(
                                "Cartesian boundary coefficient overflow or underflow",
                            ));
                        }
                        system.operator.diagonal[i] += coefficient;
                        system.operator.storage_boundary[i] += coefficient;
                        system.rhs[i] += coefficient * (equilibrium - self.moisture[0]);
                        system.faces.push(CartesianFace {
                            side,
                            cell: i,
                            coefficient,
                            fraction,
                        });
                    }
                }
            }
        }
        if system.operator.diagonal.iter().any(|&a| !positive(a))
            || system.rhs.iter().any(|a| !a.is_finite())
        {
            return Err(MoistureError::Invalid(
                "Cartesian diffusion system overflow",
            ));
        }
        Ok(system)
    }
    fn transfers(
        &self,
        moisture: &[f64],
        faces: &[ResolvedMoistureFace; 6],
        system: &CartesianSystem,
    ) -> [f64; 6] {
        let mut water = [0.0; 6];
        for face in &system.faces {
            water[face.side] += self.dry_mass_per_cell_kg
                * face.coefficient
                * (faces[face.side].equilibrium.expect("active boundary") - moisture[face.cell]);
        }
        water
    }
}

struct CartesianFace {
    side: usize,
    cell: usize,
    coefficient: f64,
    fraction: f64,
}
struct CartesianSystem {
    operator: CartesianOperator,
    rhs: Vec<f64>,
    faces: Vec<CartesianFace>,
}
impl CartesianSystem {
    fn relative_residual(&self, x: &[f64], cx: &Cx<'_>) -> Result<f64, MoistureError> {
        let mut residual = vec![0.0; x.len()];
        LinearOp::apply(&self.operator, x, &mut residual);
        let mut rscale: f64 = 0.0;
        let mut bscale: f64 = 0.0;
        for (i, (r, &b)) in residual.iter_mut().zip(&self.rhs).enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            *r = b - *r;
            if !r.is_finite() {
                return Err(MoistureError::Invalid(
                    "nonfinite spatial diffusion residual",
                ));
            }
            rscale = rscale.max(r.abs());
            bscale = bscale.max(b.abs());
        }
        if rscale == 0.0 {
            return Ok(0.0);
        }
        if bscale == 0.0 {
            return Ok(f64::INFINITY);
        }
        let rn = fs_solver::norm2(&residual.iter().map(|r| r / rscale).collect::<Vec<_>>());
        let bn = fs_solver::norm2(&self.rhs.iter().map(|b| b / bscale).collect::<Vec<_>>());
        Ok((rscale / bscale) * (rn / bn))
    }
}
struct CartesianOperator {
    diagonal: Vec<f64>,
    storage_boundary: Vec<f64>,
    forward: Vec<[f64; 3]>,
    strides: [usize; 3],
}
impl LinearOp for CartesianOperator {
    fn n(&self) -> usize {
        self.diagonal.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for ((y, &a), &x) in y.iter_mut().zip(&self.storage_boundary).zip(x) {
            *y = a * x;
        }
        for (i, coefficients) in self.forward.iter().enumerate() {
            for (k, &a) in coefficients.iter().enumerate() {
                if a != 0.0 {
                    let j = i + self.strides[k];
                    let flux = a * (x[i] - x[j]);
                    y[i] += flux;
                    y[j] -= flux;
                }
            }
        }
    }
}
impl Precond for CartesianOperator {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((z, &r), &a) in z.iter_mut().zip(r).zip(&self.diagonal) {
            *z = r / a;
        }
    }
}
