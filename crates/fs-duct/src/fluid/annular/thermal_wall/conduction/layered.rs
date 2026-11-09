//! Original radial material layers with finite interlayer thermal contacts.
//!
//! Each layer retains BOTH Dirichlet face responses of the existing radial
//! wall owner. Their exact discrete thermal admittances assemble one native
//! tridiagonal chain of physical face temperatures. Perfect contacts share
//! a coordinate; a positive areal resistance retains a temperature jump.
//! No thermal transfer matrix or division by a small through-response occurs.
//!
//! All original radial coordinates are reconstructed coherently before the
//! shared positive stiffness/mass forms are evaluated. This is exact
//! elimination of the same P1 discretization at one frequency, not a fitted
//! lumped wall or a reduced caloric model. Layer conduction, interlayer contact
//! loss, bath exchange and reversible storage remain distinct.
//!
//! Layers are ordered by increasing radius, homogeneous within each selected
//! material at the common T0, and have perfect contact with gas at the terminal
//! face. A contact resistance is explicit boundary data at T0, with no inferred
//! material authority or heat storage. The original fixed-geometry, radial,
//! linear, single-phase assumptions and source domains remain applicable.

use super::{
    AnnularSection, DuctError, PreparedRadialThermalWall, PreparedThermalTransport,
    RadialThermalWallAdmission, RadialThermalWallError, RadialThermalWallField, RadialWallBath,
    RadialWallFace, ThermalTransportSample, bad, filled, finite, poll, positive, radial_forms,
    relative,
};
use fs_la::eigen_complex::solve_tridiagonal_complex;
use fs_math::{c64::C64, dd::Dd, det};

/// Bound the face-chain tile to at most 128 coordinates before perfect-contact merging.
pub const MAX_RADIAL_THERMAL_LAYERS: usize = 64;
/// Bound the combined original radial mesh, independently of per-layer limits.
pub const MAX_LAYERED_RADIAL_ELEMENTS: usize = 4096;

/// Actual geometry and original radial transport source of one material layer.
#[derive(Clone, Copy, Debug)]
pub struct RadialThermalLayer<'a> {
    pub section: &'a AnnularSection,
    pub source: &'a PreparedThermalTransport,
}

/// Physical coherent fields integrated over one original material layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialThermalLayerObservation {
    pub face_radii_m: [f64; 2],
    pub face_temperatures_k: [C64; 2],
    /// Positive INTO this layer at each face, per that face's physical area.
    pub inward_heat_flux_w_m2: [C64; 2],
    pub conduction_dissipation_w_m: f64,
    pub thermal_storage_j_m: f64,
    /// Original C_vol integral(T' dA) [J/m], a linear phasor rather than
    /// cycle-mean quadratic acoustic thermal storage.
    pub oscillatory_thermal_energy_j_m: C64,
}

/// One interface, ordered from smaller to larger radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialThermalContactObservation {
    pub radius_m: f64,
    pub resistance_m2_k_w: f64,
    pub inner_layer_temperature_k: C64,
    pub outer_layer_temperature_k: C64,
    /// Positive from the inner layer to the outer layer.
    pub outward_heat_flux_w_m2: C64,
    /// Positive availability loss P*|T_inner-T_outer|²/(2*T0*R).
    /// Exactly zero at a perfect contact; it stores no heat.
    pub contact_dissipation_w_m: f64,
}

/// Whole layered wall for a specified coherent gas-interface temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredRadialThermalWallObservation {
    pub interface_radius_m: f64,
    pub bath_radius_m: f64,
    pub interface_temperature_k: C64,
    pub bath_temperature_k: C64,
    pub fluid_to_wall_heat_flux_w_m2: C64,
    pub wall_to_bath_heat_flux_w_m2: C64,
    pub solid_conduction_w_m: f64,
    pub thermal_contact_w_m: f64,
    pub bath_exchange_w_m: f64,
    pub thermal_storage_j_m: f64,
    pub oscillatory_thermal_energy_j_m: C64,
}

#[derive(Debug)]
struct LayerResponse<'a> {
    basis: [PreparedRadialThermalWall<'a>; 2],
    // Native heat rate INTO each face, per axial metre, per face temperature.
    ports: [C64; 4],
    face_nodes: [usize; 2],
    temperatures: [C64; 2],
    inward_flux: [C64; 2],
    inward_flux_scales: [f64; 2],
    conduction_root: f64,
    storage_root: f64,
    energy: C64,
}
impl LayerResponse<'_> {
    fn node(&self, index: usize) -> C64 {
        self.basis[0].temperature_response[index] * self.temperatures[0]
            + self.basis[1].temperature_response[index] * self.temperatures[1]
    }
    fn section(&self) -> &AnnularSection {
        self.basis[0].section
    }
    fn sample(&self) -> &ThermalTransportSample {
        &self.basis[0].sample
    }
    fn weak_flux(&self, side: usize) -> (C64, f64) {
        let section = self.section();
        let argument =
            self.basis[0].omega / self.sample().diffusivity_m2_s() * section.gap_m * section.gap_m;
        let radius = if side == 0 {
            section.inner_m
        } else {
            section.outer_m
        };
        let terms =
            super::endpoint_reaction_terms(section, argument, side, |i| self.node(i).conj());
        let factor = self.sample().conductivity_w_m_k() / section.gap_m
            * (section.inner_m + 0.5 * section.gap_m)
            / radius;
        (
            (terms[0] + terms[1]).conj().scale(factor),
            (terms[0].abs() + terms[1].abs()) * factor,
        )
    }
}

#[derive(Clone, Copy, Debug)]
struct ContactResponse {
    radius: f64,
    resistance: f64,
    temperatures: [C64; 2],
    outward_flux: C64,
    loss_root: f64,
}

/// Retained original radial responses and one solved, coherent face chain.
/// Preparation is O(total original cells + layers). Local layer/contact and
/// whole-wall queries are O(1), with no allocation or source reevaluation.
#[derive(Debug)]
pub struct PreparedLayeredRadialThermalWall<'a> {
    layers: Vec<LayerResponse<'a>>,
    contacts: Vec<ContactResponse>,
    temperature: f64,
    omega: f64,
    face: RadialWallFace,
    bath: RadialWallBath,
    admittance: C64,
    bath_flux: C64,
    bath_temperature: C64,
    loss_roots: [f64; 3],
    storage_root: f64,
    energy: C64,
    relative_backward_residual: f64,
    relative_thermal_residual: f64,
}

impl<'a> PreparedLayeredRadialThermalWall<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        layers: &[RadialThermalLayer<'a>],
        contact_resistances_m2_k_w: &[f64],
        temperature_k: f64,
        omega: f64,
        face: RadialWallFace,
        bath: RadialWallBath,
        admission: RadialThermalWallAdmission,
    ) -> Result<Self, RadialThermalWallError> {
        Self::prepare_controlled(
            layers,
            contact_resistances_m2_k_w,
            temperature_k,
            omega,
            face,
            bath,
            admission,
            || false,
        )
    }

    /// At most 64 ordered contiguous layers and 4096 combined original cells.
    /// R=0 identifies a perfect contact; R>0 retains its actual temperature jump.
    /// Payload must cover 512*(total_cells+layers)+2048*layers+8192 bytes,
    /// including original geometry, retained bases, face-chain and radial solve
    /// scratch. Independently prepared source paths and allocator metadata are
    /// excluded. Polls surround each bounded original solve and occur every
    /// 64 radial coordinates; no partially prepared wall is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_controlled(
        layers: &[RadialThermalLayer<'a>],
        contact_resistances_m2_k_w: &[f64],
        temperature_k: f64,
        omega: f64,
        face: RadialWallFace,
        bath: RadialWallBath,
        admission: RadialThermalWallAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RadialThermalWallError> {
        poll(&mut cancelled, 0)?;
        let count = layers.len();
        if count == 0
            || count > MAX_RADIAL_THERMAL_LAYERS
            || contact_resistances_m2_k_w.len() != count - 1
            || !positive(temperature_k)
            || !positive(omega)
            || !positive(admission.maximum_relative_residual)
            || admission.maximum_relative_residual >= 1.0
            || contact_resistances_m2_k_w
                .iter()
                .any(|r| !r.is_finite() || *r < 0.0)
        {
            return Err(bad(
                "invalid radial thermal layer count, contacts, state or residual admission",
            )
            .into());
        }
        let cells = layers.iter().try_fold(0usize, |sum, layer| {
            sum.checked_add(layer.section.elements())
                .filter(|n| *n <= MAX_LAYERED_RADIAL_ELEMENTS)
                .ok_or_else(|| bad("layered radial wall exceeds its combined cell bound"))
        })?;
        if 512 * (cells + count) + 2048 * count + 8192 > admission.maximum_payload_bytes
            || layers
                .windows(2)
                .any(|pair| pair[0].section.outer_m != pair[1].section.inner_m)
        {
            return Err(bad(
                "radial thermal layers must be contiguous and fit the total payload budget",
            )
            .into());
        }
        let bath_conductance = bath.conductance()?;
        let mut responses = Vec::new();
        responses
            .try_reserve_exact(count)
            .map_err(|_| bad("radial layer allocation refused"))?;
        let mut chain_nodes = 0usize;
        let mut maximum_residual = 0.0_f64;
        for (index, layer) in layers.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let inner = layer.section.prepare_radial_thermal_wall_controlled(
                layer.source,
                temperature_k,
                omega,
                RadialWallFace::Inner,
                RadialWallBath::Isothermal,
                admission,
                &mut cancelled,
            )?;
            let outer = layer.section.prepare_radial_thermal_wall_controlled(
                layer.source,
                temperature_k,
                omega,
                RadialWallFace::Outer,
                RadialWallBath::Isothermal,
                admission,
                &mut cancelled,
            )?;
            let pi = core::f64::consts::TAU * layer.section.inner_m;
            let po = core::f64::consts::TAU * layer.section.outer_m;
            let ports = [
                inner.admittance.scale(pi),
                -outer.bath_flux.scale(pi),
                -inner.bath_flux.scale(po),
                outer.admittance.scale(po),
            ];
            // Keep both original cross reactions; do not divide by small
            // transmission or average away a reciprocity failure.
            maximum_residual = maximum_residual.max(relative(
                (ports[1] - ports[2]).abs(),
                ports.iter().map(|v| v.abs()).sum(),
            ));
            for basis in [&inner, &outer] {
                maximum_residual = maximum_residual
                    .max(basis.relative_reaction_defect)
                    .max(basis.relative_heat_balance_defect);
            }
            let first = if index == 0 {
                0
            } else if contact_resistances_m2_k_w[index - 1] == 0.0 {
                chain_nodes - 1
            } else {
                chain_nodes
            };
            let face_nodes = [first, first + 1];
            chain_nodes = first + 2;
            responses.push(LayerResponse {
                basis: [inner, outer],
                ports,
                face_nodes,
                temperatures: [C64::ZERO; 2],
                inward_flux: [C64::ZERO; 2],
                inward_flux_scales: [0.0; 2],
                conduction_root: 0.0,
                storage_root: 0.0,
                energy: C64::ZERO,
            });
        }
        let radii = [layers[0].section.inner_m, layers[count - 1].section.outer_m];
        let contact_side = face.side();
        let bath_side = 1 - contact_side;
        let contact_node = contact_side * (chain_nodes - 1);
        let bath_node = bath_side * (chain_nodes - 1);
        let mut diagonal = filled(chain_nodes, C64::ZERO)?;
        let mut lower = filled(chain_nodes - 1, C64::ZERO)?;
        let mut upper = filled(chain_nodes - 1, C64::ZERO)?;
        for response in &responses {
            let [a, b] = response.face_nodes;
            diagonal[a] = diagonal[a] + response.ports[0];
            diagonal[b] = diagonal[b] + response.ports[3];
            upper[a] = upper[a] + response.ports[1];
            lower[a] = lower[a] + response.ports[2];
        }
        for index in 0..count - 1 {
            poll(&mut cancelled, index)?;
            let resistance = contact_resistances_m2_k_w[index];
            if resistance > 0.0 {
                let a = responses[index].face_nodes[1];
                let b = responses[index + 1].face_nodes[0];
                let g = core::f64::consts::TAU * layers[index].section.outer_m / resistance;
                if !positive(g) {
                    return Err(bad("radial thermal contact conductance is unrepresentable").into());
                }
                diagonal[a] = diagonal[a] + C64::from_re(g);
                diagonal[b] = diagonal[b] + C64::from_re(g);
                lower[a] = C64::from_re(-g);
                upper[a] = C64::from_re(-g);
            }
        }
        diagonal[bath_node] = diagonal[bath_node]
            + C64::from_re(bath_conductance * core::f64::consts::TAU * radii[bath_side]);
        let (temperatures, backward) = solve_faces(
            diagonal,
            lower,
            upper,
            contact_node,
            bath_node,
            bath == RadialWallBath::Isothermal,
            &mut cancelled,
        )?;
        let (mut conduction, mut storage, mut energy_re, mut energy_im) =
            (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for (index, response) in responses.iter_mut().enumerate() {
            poll(&mut cancelled, index)?;
            response.temperatures = response.face_nodes.map(|node| temperatures[node]);
            let section = response.section();
            let sample = response.sample();
            let forms = radial_forms(section, |i| response.node(i), &mut cancelled)?;
            let kform =
                sample.conductivity_w_m_k() * section.area_m2 / section.gap_m / section.gap_m
                    * forms.stiffness;
            let cform = sample.volumetric_heat_capacity_j_m3_k() * section.area_m2 * forms.mass;
            let energy = forms
                .integral
                .scale(sample.volumetric_heat_capacity_j_m3_k() * section.area_m2);
            let reactions = [response.weak_flux(0), response.weak_flux(1)];
            let fluxes = reactions.map(|(flux, _)| flux);
            let flux_scales = reactions.map(|(_, scale)| scale);
            let inward = fluxes[0].scale(core::f64::consts::TAU * section.inner_m)
                + fluxes[1].scale(core::f64::consts::TAU * section.outer_m);
            let storage_rate = C64::new(0.0, omega) * energy;
            maximum_residual = maximum_residual.max(scaled_thermal_relative(
                (inward + storage_rate).abs(),
                inward.abs()
                    + storage_rate.abs()
                    + fluxes[0].abs() * core::f64::consts::TAU * section.inner_m
                    + fluxes[1].abs() * core::f64::consts::TAU * section.outer_m,
                response.ports.iter().map(|v| v.abs()).sum(),
            ));
            if ![kform, cform]
                .into_iter()
                .all(|v| v.is_finite() && v >= 0.0)
                || !finite(energy)
                || !fluxes.into_iter().all(finite)
            {
                return Err(
                    bad("coherent radial layer forms or heat flux are unrepresentable").into(),
                );
            }
            response.inward_flux = fluxes;
            response.inward_flux_scales = flux_scales;
            response.conduction_root = det::sqrt(kform / (2.0 * temperature_k));
            response.storage_root = det::sqrt(cform / (4.0 * temperature_k));
            response.energy = energy;
            conduction = conduction + Dd::from_f64(kform);
            storage = storage + Dd::from_f64(cform);
            energy_re = energy_re + Dd::from_f64(energy.re);
            energy_im = energy_im + Dd::from_f64(energy.im);
        }
        let mut contacts = Vec::new();
        contacts
            .try_reserve_exact(count - 1)
            .map_err(|_| bad("radial contact allocation refused"))?;
        let mut contact_form = Dd::ZERO;
        for index in 0..count - 1 {
            poll(&mut cancelled, index)?;
            let radius = layers[index].section.outer_m;
            let resistance = contact_resistances_m2_k_w[index];
            let temperatures = [
                responses[index].temperatures[1],
                responses[index + 1].temperatures[0],
            ];
            let left = -responses[index].inward_flux[1];
            let right = responses[index + 1].inward_flux[0];
            let jump = temperatures[0] - temperatures[1];
            let flux = if resistance == 0.0 {
                (left + right).scale(0.5)
            } else {
                jump.scale(1.0 / resistance)
            };
            maximum_residual = maximum_residual.max(scaled_thermal_relative(
                (left - flux).abs().max((right - flux).abs()),
                responses[index].inward_flux_scales[1]
                    + responses[index + 1].inward_flux_scales[0]
                    + flux.abs(),
                responses[index]
                    .ports
                    .iter()
                    .chain(&responses[index + 1].ports)
                    .map(|v| v.abs())
                    .sum::<f64>()
                    / (core::f64::consts::TAU * radius),
            ));
            let form = if resistance == 0.0 {
                0.0
            } else {
                core::f64::consts::TAU * radius / resistance * jump.norm_sq()
            };
            if !finite(flux) || !form.is_finite() || form < 0.0 {
                return Err(bad("radial contact flux or dissipation is unrepresentable").into());
            }
            contact_form = contact_form + Dd::from_f64(form);
            contacts.push(ContactResponse {
                radius,
                resistance,
                temperatures,
                outward_flux: flux,
                loss_root: det::sqrt(form / (2.0 * temperature_k)),
            });
        }
        let contact_layer = if contact_side == 0 { 0 } else { count - 1 };
        let bath_layer = if bath_side == 0 { 0 } else { count - 1 };
        let bath_temperature = temperatures[bath_node];
        let bath_flux = if bath == RadialWallBath::Isothermal {
            -responses[bath_layer].inward_flux[bath_side]
        } else {
            bath_temperature.scale(bath_conductance)
        };
        let contact_perimeter = core::f64::consts::TAU * radii[contact_side];
        let bath_perimeter = core::f64::consts::TAU * radii[bath_side];
        let bath_form = bath_conductance * bath_perimeter * bath_temperature.norm_sq();
        let forms = [conduction.to_f64(), contact_form.to_f64(), bath_form];
        let storage_form = storage.to_f64();
        let admittance = C64::new(
            forms.into_iter().sum::<f64>() / contact_perimeter,
            -omega * storage_form / contact_perimeter,
        );
        let reaction = responses[contact_layer].inward_flux[contact_side];
        maximum_residual = maximum_residual.max(relative(
            (admittance - reaction).abs(),
            admittance.abs() + reaction.abs(),
        ));
        let energy = C64::new(energy_re.to_f64(), energy_im.to_f64());
        let input = admittance.scale(contact_perimeter);
        let output = bath_flux.scale(bath_perimeter);
        let storage_rate = C64::new(0.0, omega) * energy;
        maximum_residual = maximum_residual.max(relative(
            (input - output + storage_rate).abs(),
            input.abs() + output.abs() + storage_rate.abs(),
        ));
        let loss_roots = forms.map(|v| det::sqrt(v / (2.0 * temperature_k)));
        let storage_root = det::sqrt(storage_form / (4.0 * temperature_k));
        if ![admittance, bath_flux, bath_temperature, energy]
            .into_iter()
            .all(finite)
            || admittance.re < 0.0
            || admittance.im >= 0.0
            || !loss_roots
                .into_iter()
                .chain([storage_root])
                .all(|v| v.is_finite() && v >= 0.0)
            || backward > admission.maximum_relative_residual
            || maximum_residual > admission.maximum_relative_residual
        {
            return Err(bad(
                "layered radial wall heat, work or algebraic residual exceeds admission",
            )
            .into());
        }
        poll(&mut cancelled, cells)?;
        Ok(Self {
            layers: responses,
            contacts,
            temperature: temperature_k,
            omega,
            face,
            bath,
            admittance,
            bath_flux,
            bath_temperature,
            loss_roots,
            storage_root,
            energy,
            relative_backward_residual: backward,
            relative_thermal_residual: maximum_residual,
        })
    }

    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer_section(&self, layer: usize) -> Option<&'a AnnularSection> {
        self.layers.get(layer).map(|r| r.basis[0].section)
    }
    pub fn layer_source(&self, layer: usize) -> Option<&'a PreparedThermalTransport> {
        self.layers.get(layer).map(|r| r.basis[0].source)
    }
    pub fn layer_sample(&self, layer: usize) -> Option<&ThermalTransportSample> {
        self.layers.get(layer).map(LayerResponse::sample)
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn gas_face(&self) -> RadialWallFace {
        self.face
    }
    pub const fn bath(&self) -> RadialWallBath {
        self.bath
    }
    pub fn interface_radius_m(&self) -> f64 {
        if self.face == RadialWallFace::Inner {
            self.layers[0].section().inner_m
        } else {
            self.layers[self.layers.len() - 1].section().outer_m
        }
    }
    pub fn bath_radius_m(&self) -> f64 {
        if self.face == RadialWallFace::Inner {
            self.layers[self.layers.len() - 1].section().outer_m
        } else {
            self.layers[0].section().inner_m
        }
    }
    pub const fn interface_admittance_w_m2_k(&self) -> C64 {
        self.admittance
    }
    pub const fn relative_backward_residual(&self) -> f64 {
        self.relative_backward_residual
    }
    /// Maximum original-basis reciprocity, weak flux, contact continuity,
    /// layer/global heat balance and independent positive-work defects.
    /// Contact continuity uses original endpoint-term backward-error scales,
    /// so a nearly insulated contact does not divide by its vanishing net flux.
    /// Layer/contact balances are normalized by their original unit-temperature
    /// admittances, with a smallest-normal denominator floor at extinguished tails.
    pub const fn relative_thermal_residual(&self) -> f64 {
        self.relative_thermal_residual
    }

    pub fn field_in_cell(
        &self,
        layer: usize,
        cell: usize,
        fraction: f64,
        interface_temperature_k: C64,
    ) -> Result<RadialThermalWallField, DuctError> {
        let response = self
            .layers
            .get(layer)
            .ok_or_else(|| bad("radial layer index is invalid"))?;
        let section = response.section();
        if cell >= section.elements()
            || !fraction.is_finite()
            || !(0.0..=1.0).contains(&fraction)
            || !finite(interface_temperature_k)
        {
            return Err(bad(
                "invalid layered radial cell, fraction or interface temperature",
            ));
        }
        let x = section.fractions[cell];
        let width = section.fractions[cell + 1] - x;
        let a = response.node(cell);
        let b = response.node(cell + 1);
        let temperature = (a.scale(1.0 - fraction) + b.scale(fraction)) * interface_temperature_k;
        let gradient = (b - a).scale(1.0 / (section.gap_m * width)) * interface_temperature_k;
        let sample = response.sample();
        let flux = gradient.scale(-sample.conductivity_w_m_k());
        let conduction = gradient
            .scale(response.basis[0].local_conduction_root)
            .norm_sq();
        let storage = temperature
            .scale(response.basis[0].local_storage_root)
            .norm_sq();
        if ![temperature, gradient, flux].into_iter().all(finite)
            || ![conduction, storage].into_iter().all(f64::is_finite)
        {
            return Err(bad("layered radial field is unrepresentable"));
        }
        Ok(RadialThermalWallField {
            radius_m: section.inner_m + section.gap_m * (x + fraction * width),
            temperature_k: temperature,
            radial_temperature_gradient_k_m: gradient,
            radial_heat_flux_w_m2: flux,
            conduction_dissipation_w_m3: conduction,
            thermal_storage_j_m3: storage,
        })
    }

    pub fn layer_observation(
        &self,
        layer: usize,
        interface_temperature_k: C64,
    ) -> Result<RadialThermalLayerObservation, DuctError> {
        let response = self
            .layers
            .get(layer)
            .ok_or_else(|| bad("radial layer index is invalid"))?;
        admit_temperature(interface_temperature_k)?;
        let temperatures = response.temperatures.map(|v| v * interface_temperature_k);
        let fluxes = response.inward_flux.map(|v| v * interface_temperature_k);
        let energy = response.energy * interface_temperature_k;
        let loss = interface_temperature_k
            .scale(response.conduction_root)
            .norm_sq();
        let storage = interface_temperature_k
            .scale(response.storage_root)
            .norm_sq();
        if !temperatures
            .into_iter()
            .chain(fluxes)
            .chain([energy])
            .all(finite)
            || !loss.is_finite()
            || !storage.is_finite()
        {
            return Err(bad("radial layer observation is unrepresentable"));
        }
        Ok(RadialThermalLayerObservation {
            face_radii_m: [response.section().inner_m, response.section().outer_m],
            face_temperatures_k: temperatures,
            inward_heat_flux_w_m2: fluxes,
            conduction_dissipation_w_m: loss,
            thermal_storage_j_m: storage,
            oscillatory_thermal_energy_j_m: energy,
        })
    }

    pub fn contact_observation(
        &self,
        contact: usize,
        interface_temperature_k: C64,
    ) -> Result<RadialThermalContactObservation, DuctError> {
        let response = self
            .contacts
            .get(contact)
            .ok_or_else(|| bad("radial contact index is invalid"))?;
        admit_temperature(interface_temperature_k)?;
        let temperatures = response.temperatures.map(|v| v * interface_temperature_k);
        let flux = response.outward_flux * interface_temperature_k;
        let loss = interface_temperature_k.scale(response.loss_root).norm_sq();
        if !temperatures.into_iter().chain([flux]).all(finite) || !loss.is_finite() {
            return Err(bad("radial contact observation is unrepresentable"));
        }
        Ok(RadialThermalContactObservation {
            radius_m: response.radius,
            resistance_m2_k_w: response.resistance,
            inner_layer_temperature_k: temperatures[0],
            outer_layer_temperature_k: temperatures[1],
            outward_heat_flux_w_m2: flux,
            contact_dissipation_w_m: loss,
        })
    }

    pub fn observation(
        &self,
        interface_temperature_k: C64,
    ) -> Result<LayeredRadialThermalWallObservation, DuctError> {
        admit_temperature(interface_temperature_k)?;
        let bath_temperature = self.bath_temperature * interface_temperature_k;
        let input = self.admittance * interface_temperature_k;
        let output = self.bath_flux * interface_temperature_k;
        let energy = self.energy * interface_temperature_k;
        let losses = self
            .loss_roots
            .map(|r| interface_temperature_k.scale(r).norm_sq());
        let storage = interface_temperature_k.scale(self.storage_root).norm_sq();
        if ![bath_temperature, input, output, energy]
            .into_iter()
            .all(finite)
            || !losses.into_iter().chain([storage]).all(f64::is_finite)
        {
            return Err(bad("layered wall observation is unrepresentable"));
        }
        Ok(LayeredRadialThermalWallObservation {
            interface_radius_m: self.interface_radius_m(),
            bath_radius_m: self.bath_radius_m(),
            interface_temperature_k,
            bath_temperature_k: bath_temperature,
            fluid_to_wall_heat_flux_w_m2: input,
            wall_to_bath_heat_flux_w_m2: output,
            solid_conduction_w_m: losses[0],
            thermal_contact_w_m: losses[1],
            bath_exchange_w_m: losses[2],
            thermal_storage_j_m: storage,
            oscillatory_thermal_energy_j_m: energy,
        })
    }
}

// Normalize dimensional heat equations before applying the same underflow
// safeguard as the unit-temperature radial/face solves. This only accounts
// for unavoidable subnormal rounding; no temperature or physical loss changes.
// The whole-wall heat/work check still uses its actual nonzero source response.
fn scaled_thermal_relative(error: f64, scale: f64, unit_response: f64) -> f64 {
    if !positive(unit_response) {
        return f64::INFINITY;
    }
    relative(
        error / unit_response,
        (scale / unit_response).max(f64::MIN_POSITIVE),
    )
}

fn admit_temperature(temperature: C64) -> Result<(), DuctError> {
    if finite(temperature) {
        Ok(())
    } else {
        Err(bad("layered wall interface temperature must be finite"))
    }
}

// Solve thermal ADMITTANCES directly; a vanishing cross-admittance simply
// decouples faces. The native pivoted tridiagonal solve handles the same
// complex phasors as the original radial owner.
fn solve_faces(
    mut diagonal: Vec<C64>,
    mut lower: Vec<C64>,
    mut upper: Vec<C64>,
    contact: usize,
    bath: usize,
    isothermal_bath: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64), DuctError> {
    let n = diagonal.len();
    let scale = 1.0
        / diagonal
            .iter()
            .chain(&lower)
            .chain(&upper)
            .map(|v| v.abs())
            .fold(1.0_f64, f64::max);
    if !positive(scale)
        || diagonal
            .iter()
            .chain(&lower)
            .chain(&upper)
            .any(|v| !finite(*v))
    {
        return Err(bad("layered wall face matrix is unrepresentable"));
    }
    for v in diagonal.iter_mut().chain(&mut lower).chain(&mut upper) {
        *v = v.scale(scale);
    }
    let mut drive = filled(n, C64::ZERO)?;
    let edge = if contact == 0 { 0 } else { n - 2 };
    let neighbor = if contact == 0 { 1 } else { n - 2 };
    drive[neighbor] = -if contact == 0 {
        lower[edge]
    } else {
        upper[edge]
    };
    drive[contact] = C64::ONE;
    diagonal[contact] = C64::ONE;
    lower[edge] = C64::ZERO;
    upper[edge] = C64::ZERO;
    if isothermal_bath {
        let edge = if bath == 0 { 0 } else { n - 2 };
        diagonal[bath] = C64::ONE;
        lower[edge] = C64::ZERO;
        upper[edge] = C64::ZERO;
        drive[bath] = C64::ZERO;
    }
    let mut values = filled(n, C64::ZERO)?;
    values.copy_from_slice(&drive);
    poll(cancelled, n)?;
    solve_tridiagonal_complex(&lower, &diagonal, &upper, &mut values, 1)
        .map_err(|_| DuctError::Singular)?;
    poll(cancelled, n)?;
    let mut residual = 0.0_f64;
    for i in 0..n {
        let mut error = diagonal[i] * values[i] - drive[i];
        let mut magnitude = diagonal[i].abs() * values[i].abs() + drive[i].abs();
        if i > 0 {
            error = error + lower[i - 1] * values[i - 1];
            magnitude += lower[i - 1].abs() * values[i - 1].abs();
        }
        if i + 1 < n {
            error = error + upper[i] * values[i + 1];
            magnitude += upper[i].abs() * values[i + 1].abs();
        }
        residual = residual.max(relative(error.abs(), magnitude.max(f64::MIN_POSITIVE)));
    }
    Ok((values, residual))
}
