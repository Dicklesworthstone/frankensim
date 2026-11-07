//! Two-sided radiation loading of the original tensor-Prony shell.
//!
//! The shell's actual P1 translations drive the pressure-jump screen equation.
//! Its transpose returns the acoustic reaction to the same physical coordinates,
//! including all translations, rotations, inertia and material relaxation. The
//! jump already accounts for both wet faces; a second surface or factor of two
//! would count the fluid twice. All phasors use exp(-i omega t), peak amplitude.
//!
//! This is the midsurface screen approximation, with planar acoustic refinement
//! independent of the original structural mesh. Physical thickness remains in
//! the shell mechanics. Both faces must contact the same homogeneous fluid;
//! cavity-backed faces and finite-thickness acoustic skins need different owners.
//! For the coupled screen formulation see Antil, Hardesty and Heinkenschloss,
//! *Shape Optimization of Shell Structure Acoustics*, sections 2.2–2.3:
//! <https://math.gmu.edu/~hantil/preprints/HAntil_MHeinkenschloss_SHardesty_2016.pdf>.

use super::{TensorShellRadiationError, admit, checkpoint, midpoint, relative_balance, reserved};
use fs_bem::helmholtz::{MAX_DENSE_PANELS, MAX_RADIATION_FIELDS_PER_BATCH, Medium};
use fs_bem::near_field::Options as ReceiverOptions;
use fs_bem::reciprocal_radiation::open_screen::{
    OpenOptions, OpenPanels, PreparedScreen, ScreenPortLimits, ScreenVelocityPortSolution,
};
use fs_exec::Cx;
use fs_math::{c64::C64, det};
use fs_modal::harmonic::{HarmonicPortControls, HarmonicResponse, HarmonicResponseControls};
use fs_plate::shell::viscosity::harmonic::{PreparedTensorPronyShell, TensorPronyShell};
use std::collections::BTreeMap;

/// Original shell-to-screen geometry and complete physical trace limits.
#[derive(Clone, Copy, Debug)]
pub struct OpenTensorShellRadiationControls {
    /// Conforming planar quartering levels, at most six. This is not edge grading.
    pub acoustic_refinements: u32,
    /// Complete refined triangle count, including entirely supported facets.
    pub maximum_panels: usize,
    /// All original free translations and rotations; no coordinate is truncated.
    pub maximum_coordinates: usize,
    /// Logical numeric payload bound: 24*q*P + 8*q*q + 640*P + 24*q + 256,
    /// with q original free coordinates and P final triangles. Allocator/map
    /// bookkeeping, borrowed shell and native BEM cache are separate.
    pub maximum_setup_bytes: usize,
    /// Conservative setup score: 512*P*(levels+1) + 64*q*P + q*q.
    pub maximum_setup_work: u64,
}

/// Native boundary and full-coordinate mechanical limits at one frequency.
#[derive(Clone, Copy, Debug)]
pub struct OpenTensorShellRadiationFrequencyControls {
    /// Original full-coordinate material preparation and solve admission.
    pub mechanical: HarmonicResponseControls,
    /// Full physical force/torque load projection and passive-work admission.
    pub ports: HarmonicPortControls,
    /// Complete self/mutual matrix and original unit-source pressure jumps.
    pub boundary: ScreenPortLimits,
    /// Additional numerical frequency payload, excluding the source geometry,
    /// static BEM cache and the independently bounded mechanical factor.
    /// Serial bound: 96*N*N + 64*q*(N+V) + 32*q*q + 256*P + 65536.
    /// With W>1 native workers, add 160*P*min(P,4*W) for triangle-pair
    /// results, (16*N+64)*A for retained angular coefficients/direction tables,
    /// and 128*(P+W+N) for row/index bookkeeping. Here N is active pressure
    /// unknowns, V vertices, P facets, q original coordinates and
    /// A=64*4^angular_refinement. The sum conservatively counts parallel
    /// assembly phases together; allocator overhead and thread stacks are excluded.
    pub maximum_boundary_payload_bytes: usize,
    /// Coherent pressure-jump assembly and antipodal far-field power work.
    pub maximum_radiated_power_work: u64,
    /// Relative input/heat/flux, port/jump and jump/flux defects, in (0,1).
    pub power_relative_tolerance: f64,
}

/// Immutable original shell, oriented open surface and facet-local P1 traces.
pub struct OpenTensorShellRadiation<'shell> {
    shell: &'shell TensorPronyShell,
    surface: OpenPanels,
    velocity_rows: Vec<Vec<[f64; 3]>>,
    identity_ports: Vec<f64>,
    medium: Medium,
}

impl<'shell> OpenTensorShellRadiation<'shell> {
    /// Bind both wet faces to one homogeneous fluid. Shared original edges are
    /// refined by topology, without welding, capping or changing the source mesh.
    /// Free-edge PRESSURE vanishes in the native screen space; physical velocity
    /// at a free edge is retained. Supported facets still scatter sound.
    pub fn new(
        shell: &'shell TensorPronyShell,
        medium: Medium,
        controls: OpenTensorShellRadiationControls,
        cx: &Cx<'_>,
    ) -> Result<Self, TensorShellRadiationError> {
        checkpoint(cx)?;
        if !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || controls.acoustic_refinements > 6
            || shell.model().free == 0
        {
            return Err(TensorShellRadiationError::Invalid(
                "open shell requires a positive medium, original free coordinates and at most six refinements",
            ));
        }
        let q = shell.model().free;
        let panels = shell
            .mesh()
            .tris
            .len()
            .checked_mul(4usize.pow(controls.acoustic_refinements))
            .ok_or(TensorShellRadiationError::Invalid(
                "open shell refinement overflow",
            ))?;
        admit(
            "open shell original coordinates",
            q as u128,
            controls
                .maximum_coordinates
                .min(MAX_RADIATION_FIELDS_PER_BATCH) as u128,
        )?;
        admit(
            "open shell panels",
            panels as u128,
            controls.maximum_panels.min(MAX_DENSE_PANELS) as u128,
        )?;
        let (n, p) = (q as u128, panels as u128);
        admit(
            "open shell setup bytes",
            24 * n * p + 8 * n * n + 640 * p + 24 * n + 256,
            controls.maximum_setup_bytes as u128,
        )?;
        admit(
            "open shell setup work",
            512 * p * (u128::from(controls.acoustic_refinements) + 1) + 64 * n * p + n * n,
            controls.maximum_setup_work as u128,
        )?;
        let (surface, traces) = refined_surface(shell, controls.acoustic_refinements, panels, cx)?;
        let mut velocity_rows = reserved(q, "open shell trace rows")?;
        for _ in 0..q {
            let mut row = reserved(panels, "open shell facet-corner trace")?;
            row.resize(panels, [0.0; 3]);
            velocity_rows.push(row);
        }
        for (panel, ((element, weights), triangle)) in
            traces.iter().zip(surface.triangles()).enumerate()
        {
            checkpoint(cx)?;
            let points = triangle.map(|node| surface.vertices()[node]);
            let a: [f64; 3] = core::array::from_fn(|d| points[1][d] - points[0][d]);
            let b: [f64; 3] = core::array::from_fn(|d| points[2][d] - points[0][d]);
            let cross = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            let length = det::sqrt(cross.iter().map(|v| v * v).sum());
            if !length.is_finite() || length <= 0.0 {
                return Err(TensorShellRadiationError::Invalid(
                    "unrepresentable open shell facet normal",
                ));
            }
            let normal = cross.map(|v| v / length);
            for (original, &node) in shell.mesh().tris[*element].iter().enumerate() {
                for (component, &direction) in normal.iter().enumerate() {
                    if let Some(free) = shell.model().dof_map[6 * node + component] {
                        for corner in 0..3 {
                            velocity_rows[free][panel][corner] +=
                                weights[corner][original] * direction;
                        }
                    }
                }
            }
        }
        let mut identity_ports = reserved(q * q, "open shell physical identity ports")?;
        identity_ports.resize(q * q, 0.0);
        for i in 0..q {
            identity_ports[i * q + i] = 1.0;
        }
        checkpoint(cx)?;
        Ok(Self {
            shell,
            surface,
            velocity_rows,
            identity_ports,
            medium,
        })
    }

    /// Immutable original material, geometry and support owner.
    pub const fn shell(&self) -> &'shell TensorPronyShell {
        self.shell
    }
    /// Refined original midsurface, with no artificial second wet face.
    pub const fn surface(&self) -> &OpenPanels {
        &self.surface
    }
    /// The same homogeneous fluid on both sides of the shell.
    pub const fn medium(&self) -> Medium {
        self.medium
    }
    /// [original free coordinate][refined facet][corner] normal velocity map.
    /// Rotation rows remain zero at the zero-offset midsurface, while their
    /// original structural coupling and inertia stay in the mechanical solve.
    pub fn normal_velocity_rows(&self) -> &[Vec<[f64; 3]>] {
        &self.velocity_rows
    }

    /// Prepare the original two-sided screen once over the entire source band.
    /// Edge, wavelength, singular quadrature, angular and work restrictions are
    /// those of the native screen solver, without an exterior closed-body test.
    pub fn prepare(
        &self,
        highest_omega_rad_s: f64,
        options: OpenOptions,
        cx: &Cx<'_>,
    ) -> Result<PreparedOpenTensorShellRadiation<'_, 'shell>, TensorShellRadiationError> {
        checkpoint(cx)?;
        if !highest_omega_rad_s.is_finite() || highest_omega_rad_s <= 0.0 {
            return Err(TensorShellRadiationError::Invalid(
                "positive finite open acoustic band required",
            ));
        }
        let screen = PreparedScreen::new(
            &self.surface,
            highest_omega_rad_s / self.medium.sound_speed,
            options,
            || cx.checkpoint().is_err(),
        )?;
        admit(
            "open shell complete source work",
            screen.velocity_port_work(self.velocity_rows.len())? as u128,
            options.maximum_work as u128,
        )?;
        checkpoint(cx)?;
        Ok(PreparedOpenTensorShellRadiation {
            source: self,
            screen,
            highest_omega_rad_s,
            options,
        })
    }
}

type FacetTrace = (usize, [[f64; 3]; 3]);

fn refined_surface(
    shell: &TensorPronyShell,
    levels: u32,
    panels: usize,
    cx: &Cx<'_>,
) -> Result<(OpenPanels, Vec<FacetTrace>), TensorShellRadiationError> {
    let mut vertices = reserved(3 * panels, "open shell refined vertices")?;
    vertices.extend_from_slice(&shell.mesh().nodes);
    let mut triangles = reserved(panels, "open shell refined triangles")?;
    triangles.extend_from_slice(&shell.mesh().tris);
    let mut traces = reserved(panels, "open shell original P1 traces")?;
    traces.extend(
        (0..triangles.len()).map(|e| (e, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])),
    );
    for _ in 0..levels {
        checkpoint(cx)?;
        let mut edges = BTreeMap::new();
        let mut next_triangles = reserved(4 * triangles.len(), "open shell child triangles")?;
        let mut next_traces = reserved(4 * traces.len(), "open shell child traces")?;
        for (triangle, (element, barycentric)) in triangles.iter().zip(&traces) {
            checkpoint(cx)?;
            let mut ids = [triangle[0], triangle[1], triangle[2], 0, 0, 0];
            for (edge, (a, b)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
                let (a, b) = (triangle[a], triangle[b]);
                ids[edge + 3] = *edges.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    let index = vertices.len();
                    vertices.push(midpoint(vertices[a], vertices[b]));
                    index
                });
            }
            let b = [
                barycentric[0],
                barycentric[1],
                barycentric[2],
                midpoint(barycentric[0], barycentric[1]),
                midpoint(barycentric[1], barycentric[2]),
                midpoint(barycentric[2], barycentric[0]),
            ];
            for child in [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]] {
                next_triangles.push(child.map(|i| ids[i]));
                next_traces.push((*element, child.map(|i| b[i])));
            }
        }
        triangles = next_triangles;
        traces = next_traces;
    }
    checkpoint(cx)?;
    let surface = OpenPanels::new(vertices, triangles)?;
    checkpoint(cx)?;
    Ok((surface, traces))
}

/// Shared geometry and singular quadrature for the original open shell.
pub struct PreparedOpenTensorShellRadiation<'geometry, 'shell> {
    source: &'geometry OpenTensorShellRadiation<'shell>,
    screen: PreparedScreen<'geometry>,
    highest_omega_rad_s: f64,
    options: OpenOptions,
}
impl<'geometry, 'shell> PreparedOpenTensorShellRadiation<'geometry, 'shell> {
    /// The exact source and physical trace retained by this preparation.
    pub const fn source(&self) -> &'geometry OpenTensorShellRadiation<'shell> {
        self.source
    }
    /// Native boundary owner, for independent source and work checks.
    pub const fn boundary(&self) -> &PreparedScreen<'geometry> {
        &self.screen
    }

    /// One acoustic factor for all original velocity columns, followed by one
    /// loaded native material factor. Repeated forces reuse both preparations.
    pub fn prepare_frequency(
        &self,
        omega_rad_s: f64,
        controls: OpenTensorShellRadiationFrequencyControls,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedOpenTensorShellRadiationFrequency<'_, 'geometry, 'shell>,
        TensorShellRadiationError,
    > {
        checkpoint(cx)?;
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || omega_rad_s > self.highest_omega_rad_s
            || !controls.power_relative_tolerance.is_finite()
            || controls.power_relative_tolerance <= 0.0
            || controls.power_relative_tolerance >= 1.0
        {
            return Err(TensorShellRadiationError::Invalid(
                "open shell frequency or power tolerance outside its admitted range",
            ));
        }
        let q = self.source.velocity_rows.len();
        let (n, p, v, ports) = (
            self.screen.unknowns() as u128,
            self.source.surface.triangles().len() as u128,
            self.source.surface.vertices().len() as u128,
            q as u128,
        );
        admit(
            "open shell harmonic coordinates",
            ports,
            controls.mechanical.maximum_coordinates as u128,
        )?;
        admit(
            "open shell harmonic ports",
            ports,
            controls.ports.maximum_ports as u128,
        )?;
        let angular = 64_u128 * (1_u128 << (2 * self.options.angular_refinement));
        let workers = self.options.workers as u128;
        let parallel_scratch = if workers > 1 {
            // The native screen retains up to four target rows per worker.
            // Each source/target pair is ten C64 values. Its parallel Gram
            // assembly also retains one N-entry C64 vector per direction,
            // beside the direction/weight and (weight, Vec) arrays. Both
            // phases are admitted before calling the original boundary owner.
            let block_rows = p.min(4 * workers);
            160 * p * block_rows + (16 * n + 64) * angular + 128 * (p + workers + n)
        } else {
            0
        };
        admit(
            "open shell frequency payload",
            96 * n * n
                + 64 * ports * (n + v)
                + 32 * ports * ports
                + 256 * p
                + 65536
                + parallel_scratch,
            controls.maximum_boundary_payload_bytes as u128,
        )?;
        admit(
            "open shell coherent power work",
            64 * ports * (v + ports) + 2048 * angular * (n + 16 * p),
            controls.maximum_radiated_power_work as u128,
        )?;
        let mut rows = reserved(q, "open shell complete source row references")?;
        rows.extend(self.source.velocity_rows.iter().map(Vec::as_slice));
        let radiation = self.screen.solve_velocity_ports(
            omega_rad_s / self.source.medium.sound_speed,
            self.source.medium,
            &rows,
            controls.boundary,
            || cx.checkpoint().is_err(),
        )?;
        let shell: &'shell TensorPronyShell = self.source.shell;
        let mechanical = shell.prepare_harmonic_response_with_port_impedance(
            omega_rad_s,
            &self.source.identity_ports,
            q,
            radiation.impedance(),
            controls.mechanical,
            controls.ports,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(PreparedOpenTensorShellRadiationFrequency {
            radiation,
            mechanical,
            controls,
        })
    }
}

/// Reusable loaded full-coordinate shell factor and original pressure-jump sources.
pub struct PreparedOpenTensorShellRadiationFrequency<'p, 'geometry, 'shell> {
    radiation: ScreenVelocityPortSolution<'p, 'geometry>,
    mechanical: PreparedTensorPronyShell<'shell>,
    controls: OpenTensorShellRadiationFrequencyControls,
}
impl<'p, 'geometry, 'shell> PreparedOpenTensorShellRadiationFrequency<'p, 'geometry, 'shell> {
    /// Loaded original material factor, reusable for stress and heat observation.
    pub const fn shell_preparation(&self) -> &PreparedTensorPronyShell<'shell> {
        &self.mechanical
    }
    /// Complete original physical unit-velocity fields and their work matrix.
    pub const fn radiation(&self) -> &ScreenVelocityPortSolution<'p, 'geometry> {
        &self.radiation
    }
    /// Full row-major physical impedance; zero direct rotational rows remain.
    pub fn impedance(&self) -> &[C64] {
        self.radiation.impedance()
    }

    /// Solve force/torque in the original coordinate chart and independently
    /// integrate the same coherent field over both radiating half-spaces.
    /// Input work = intrinsic material heat + radiated power; reaction is an
    /// internal transfer and is never charged a second time as material heat.
    pub fn solve(
        &self,
        force: &[C64],
        cx: &Cx<'_>,
    ) -> Result<OpenTensorShellRadiationResponse<'_, 'p, 'geometry>, TensorShellRadiationError>
    {
        checkpoint(cx)?;
        let motion = self.mechanical.solve(force, cx)?;
        let jump_power_w = self.radiation.input_power_w(&motion.velocity)?;
        let radiated_power_w = self
            .radiation
            .radiated_power_w(&motion.velocity, || cx.checkpoint().is_err())?;
        let port = motion
            .port_load
            .as_ref()
            .ok_or(TensorShellRadiationError::Invalid(
                "native open shell omitted original physical port work",
            ))?;
        let port_work = relative_balance(port.average_absorbed_power_w, jump_power_w, 0.0);
        let flux = relative_balance(jump_power_w, radiated_power_w, 0.0);
        let coupled = relative_balance(
            motion.average_input_power_w,
            motion.average_viscous_power_w,
            radiated_power_w,
        );
        let tolerance = self.controls.power_relative_tolerance;
        if !port_work.is_finite()
            || !flux.is_finite()
            || !coupled.is_finite()
            || port_work > tolerance
            || flux > tolerance
            || coupled > tolerance
        {
            return Err(TensorShellRadiationError::PowerBalance {
                coupled,
                port_work,
                flux,
                tolerance,
            });
        }
        checkpoint(cx)?;
        Ok(OpenTensorShellRadiationResponse {
            motion,
            radiation: &self.radiation,
            jump_power_w,
            radiated_power_w,
            relative_coupled_power_defect: coupled,
            relative_port_work_defect: port_work,
            relative_flux_defect: flux,
        })
    }
}

/// One accepted motion and its SAME acoustic field, retained by immutable borrow.
pub struct OpenTensorShellRadiationResponse<'r, 'p, 'geometry> {
    motion: HarmonicResponse,
    radiation: &'r ScreenVelocityPortSolution<'p, 'geometry>,
    jump_power_w: f64,
    radiated_power_w: f64,
    relative_coupled_power_defect: f64,
    relative_port_work_defect: f64,
    relative_flux_defect: f64,
}
impl OpenTensorShellRadiationResponse<'_, '_, '_> {
    /// Original displacement/rotation, velocity, source work and intrinsic heat.
    pub const fn motion(&self) -> &HarmonicResponse {
        &self.motion
    }
    /// Signed cycle-average work supplied to the two-sided pressure jump [W].
    pub const fn jump_power_w(&self) -> f64 {
        self.jump_power_w
    }
    /// Independently integrated coherent outgoing sound power [W].
    pub const fn radiated_power_w(&self) -> f64 {
        self.radiated_power_w
    }
    /// Relative defect in input = intrinsic heat + outgoing sound power.
    pub const fn relative_coupled_power_defect(&self) -> f64 {
        self.relative_coupled_power_defect
    }
    /// Relative defect between native mechanical port and pressure-jump work.
    pub const fn relative_port_work_defect(&self) -> f64 {
        self.relative_port_work_defect
    }
    /// Relative defect between pressure-jump work and integrated outgoing flux.
    pub const fn relative_flux_defect(&self) -> f64 {
        self.relative_flux_defect
    }

    /// Full nodal pressure difference in original acoustic-vertex order.
    pub fn pressure_jump(&self, cx: &Cx<'_>) -> Result<Vec<C64>, TensorShellRadiationError> {
        checkpoint(cx)?;
        Ok(self
            .radiation
            .pressure_jump(&self.motion.velocity, || cx.checkpoint().is_err())?)
    }
    /// Physical pressure on either side through native adaptive P1 integration.
    pub fn pressure_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: ReceiverOptions,
        cx: &Cx<'_>,
    ) -> Result<Vec<C64>, TensorShellRadiationError> {
        checkpoint(cx)?;
        Ok(self.radiation.pressure_at(
            &self.motion.velocity,
            points,
            minimum_clearance_m,
            options,
            || cx.checkpoint().is_err(),
        )?)
    }
    /// Coherent outgoing amplitude r exp(-ikr) p [Pa m], in a unit direction.
    pub fn far_field(
        &self,
        direction: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<C64, TensorShellRadiationError> {
        checkpoint(cx)?;
        Ok(self
            .radiation
            .far_field(&self.motion.velocity, direction, || {
                cx.checkpoint().is_err()
            })?)
    }
}
