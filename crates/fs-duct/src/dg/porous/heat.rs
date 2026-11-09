//! Independent positive porous source integrals on the native DG trajectory.
use super::*;
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedSourceLoss};

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps: 0 })
    } else {
        Ok(())
    }
}
fn integral_error(error: PreparedStepError) -> DuctError {
    match error {
        PreparedStepError::Cancelled => DuctError::Cancelled { completed_steps: 0 },
        PreparedStepError::Solver(reason) => DuctError::PorousHeatUnresolved { reason },
    }
}

/// Intrinsic irreversible porous losses integrated over the accepted interval.
/// The receiving thermal owner supplies its physical destinations. Reversible
/// memory energy, terminal work and numerical residual are not deposited here.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DgPorousHeating {
    pub viscous_j: f64,
    pub thermal_j: f64,
}
impl DgPorousHeating {
    pub fn total_j(&self) -> f64 {
        self.viscous_j + self.thermal_j
    }
}

/// Read-only, separately positive viscous and thermal heat caches for the full
/// original operator and held-input interval. Warm observation advances no
/// state and does not reconstruct heat from acoustic energy disappearance.
#[derive(Clone, Debug)]
pub struct DgPorousHeat {
    operator_identity: ContentHash,
    ports: Vec<PortId>,
    dimension: usize,
    viscous: PreparedSourceLoss,
    thermal: PreparedSourceLoss,
}

/// Original porous losses localized to one complete physical DG element.
/// The native mass transformations are element-local, so selecting ALL its
/// energy coordinates integrates its volume without another spatial rule.
/// Sources outside the element remain in the full propagation operator.
#[derive(Clone, Debug)]
pub struct DgPorousElementHeat {
    element: usize,
    interval_m: [f64; 2],
    area_m2: f64,
    heat: DgPorousHeat,
}
impl DgPorousElementHeat {
    /// Original left-to-right element index.
    pub const fn element(&self) -> usize {
        self.element
    }
    /// Physical near and far face positions along the original duct axis [m].
    pub const fn interval_m(&self) -> [f64; 2] {
        self.interval_m
    }
    /// Full specimen face area [m²]; porosity belongs to the material law.
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Both original mechanism integrals restricted to this element's volume.
    pub const fn heat(&self) -> &DgPorousHeat {
        &self.heat
    }
}

impl DgPorousHeat {
    pub const fn operator_identity(&self) -> ContentHash {
        self.operator_identity
    }
    pub fn duration(&self) -> f64 {
        self.viscous.duration()
    }
    /// Caller-owned warm scratch, including held efforts and source coordinates.
    pub fn scratch_len(&self) -> usize {
        self.ports.len() + 2 * self.dimension
    }

    /// Observe both original heat channels from the INITIAL complete state
    /// under the same held ports as the native step, without advancing it.
    pub fn integrated_heat_j(
        &self,
        state: &DgState,
        ports: &[PortInput],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgPorousHeating, DuctError> {
        let mut scratch = vec![0.0; self.scratch_len()];
        self.integrated_heat_into_j(state, ports, &mut scratch, cancelled)
    }

    /// The same positive observations with no warm allocation or exponential.
    /// Scratch may change on refusal; the original physical state never does.
    pub fn integrated_heat_into_j(
        &self,
        state: &DgState,
        ports: &[PortInput],
        scratch: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgPorousHeating, DuctError> {
        poll(cancelled)?;
        if state.operator_identity != self.operator_identity || scratch.len() != self.scratch_len()
        {
            return Err(invalid(
                "porous heat requires its complete original DG state and exact scratch",
            ));
        }
        let (efforts, coordinates) = scratch.split_at_mut(self.ports.len());
        efforts.fill(0.0);
        for port in ports {
            if !port.effort.is_finite() {
                return Err(invalid("porous heat held port effort must be finite"));
            }
            let slot = self
                .ports
                .iter()
                .position(|id| *id == port.id)
                .ok_or_else(|| invalid("porous heat names an undeclared native port"))?;
            efforts[slot] += port.effort;
            if !efforts[slot].is_finite() {
                return Err(invalid("porous heat summed port effort is unrepresentable"));
            }
        }
        let viscous_j = self
            .viscous
            .evaluate_into_controlled(&state.tilde, efforts, coordinates, &mut *cancelled)
            .map_err(integral_error)?;
        let thermal_j = self
            .thermal
            .evaluate_into_controlled(&state.tilde, efforts, coordinates, &mut *cancelled)
            .map_err(integral_error)?;
        if !(viscous_j + thermal_j).is_finite() {
            return Err(invalid("porous total irreversible heat is unrepresentable"));
        }
        Ok(DgPorousHeating {
            viscous_j,
            thermal_j,
        })
    }
}

fn row(rows: &mut Vec<f64>, n: usize) -> &mut [f64] {
    let start = rows.len();
    rows.resize(start + n, 0.0);
    &mut rows[start..]
}

impl<M: DgMedium> DgDuct<M> {
    /// Prepare independent viscous and thermal heat integrals on the complete
    /// acoustic-plus-memory trajectory. The original static Darcy drag and
    /// every retained arm contribute their own squared dissipative source.
    /// The unchanged native propagator supplies the same full operator and dt;
    /// attached structural memories and every declared port remain in the
    /// trajectory but contribute no porous heat rows.
    ///
    /// The complete state is bounded to 128 coordinates by the existing
    /// positive source-integral owner. Both caches and their preparation share
    /// the supplied total byte/work allowance. No spatial state, material arm
    /// or weak loss is dropped to fit the budget. An absent porous source gives
    /// exactly zero heat. Thermal deposition/temperature refresh is a separate
    /// receiving-owner transaction, not a mutation performed by this cache.
    pub fn prepare_porous_heat(
        &self,
        propagator: &DgPropagator,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgPorousHeat, DuctError> {
        self.prepare_porous_heat_on_dofs(propagator, 0..self.field_dofs, budget, cancelled)
    }

    /// Integrate viscous and thermal absorption in one physical DG element.
    ///
    /// The original static drag and memory source squares are selected by
    /// their element-local mass coordinates, BEFORE time integration. Their
    /// evolution uses the COMPLETE unchanged operator and held port inputs;
    /// interelement transport and every memory therefore remain coupled.
    /// Summing a complete element partition recovers the whole-specimen heat,
    /// including acoustic startup and ringdown, without assuming uniform heat
    /// or a cycle-mean tone. No endpoint power or total energy deficit is used.
    ///
    /// This optional cold preparation retains two native positive integral
    /// caches per requested element; ordinary whole-specimen queries allocate
    /// no additional storage. The supplied allowance covers this one element.
    /// A consumer preparing many elements must partition its total allowance.
    /// The same 128-state cap, frozen material and source-band limits apply.
    pub fn prepare_porous_element_heat(
        &self,
        propagator: &DgPropagator,
        element: usize,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgPorousElementHeat, DuctError> {
        poll(cancelled)?;
        let geometry = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("porous heat needs an original element index"))?;
        if self.gas.porous_diffusive_model().is_none() {
            return Err(invalid(
                "element heating requires an intrinsic porous material",
            ));
        }
        let width = self.mesh.order + 1;
        let start = element * width;
        let heat =
            self.prepare_porous_heat_on_dofs(propagator, start..start + width, budget, cancelled)?;
        poll(cancelled)?;
        Ok(DgPorousElementHeat {
            element,
            interval_m: [geometry.x_left, geometry.x_left + geometry.length],
            area_m2: geometry.areas[0],
            heat,
        })
    }

    fn prepare_porous_heat_on_dofs(
        &self,
        propagator: &DgPropagator,
        selected: core::ops::Range<usize>,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgPorousHeat, DuctError> {
        poll(cancelled)?;
        let n = self.total_dofs();
        let m = self.field_dofs;
        let source_width = selected.len();
        let ports = self.ports.len();
        if propagator.operator_identity != self.operator_identity
            || n > 128
            || 2 * n > budget.max_dimension
        {
            return Err(invalid(
                "porous heat needs its native propagator and complete bounded state",
            ));
        }
        let model = self.gas.porous_diffusive_model();
        let viscous_count = model
            .map_or(Some(0), |model| {
                model
                    .series_terms()
                    .len()
                    .checked_add(1)
                    .and_then(|v| v.checked_mul(source_width))
            })
            .ok_or_else(|| invalid("porous viscous source extent overflow"))?;
        let thermal_count = model
            .map_or(Some(0), |model| {
                model.shunt_terms().len().checked_mul(source_width)
            })
            .ok_or_else(|| invalid("porous thermal source extent overflow"))?;
        let bytes = viscous_count
            .checked_add(thermal_count)
            .and_then(|v| v.checked_add(ports))
            .and_then(|v| v.checked_mul(n))
            .and_then(|v| v.checked_mul(8));
        if viscous_count > n
            || thermal_count > n
            || bytes.is_none_or(|v| v > budget.max_dense_bytes)
        {
            return Err(invalid(
                "porous source rows exceed the complete byte allowance",
            ));
        }
        let channel_budget = DenseExponentialBudget {
            max_dimension: budget.max_dimension,
            max_dense_bytes: (budget.max_dense_bytes - bytes.unwrap()) / 2,
            max_multiply_terms: budget.max_multiply_terms / 2,
        };
        let mut viscous = Vec::with_capacity(viscous_count * n);
        let mut thermal = Vec::with_capacity(thermal_count * n);
        if let Some(model) = model {
            let drag_root = det::sqrt(model.series_drag_per_s());
            for i in selected.clone() {
                poll(cancelled)?;
                row(&mut viscous, n)[m + i] = drag_root;
            }
            for (branch, (offset, term, is_thermal)) in model
                .series_terms()
                .iter()
                .map(|t| (m, t, false))
                .chain(model.shunt_terms().iter().map(|t| (0, t, true)))
                .enumerate()
            {
                let root_gain = det::sqrt(term.gain_per_s);
                let root_rate = det::sqrt(term.rate_per_s);
                for i in selected.clone() {
                    poll(cancelled)?;
                    let source = if is_thermal {
                        row(&mut thermal, n)
                    } else {
                        row(&mut viscous, n)
                    };
                    source[offset + i] = root_gain;
                    source[2 * m + branch * m + i] = -root_rate;
                }
            }
        }
        let mut inputs = vec![0.0; n * ports];
        for (port, (_, column)) in self.ports.iter().enumerate() {
            poll(cancelled)?;
            for i in 0..n {
                inputs[i * ports + port] = column[i];
            }
        }
        let viscous = PreparedSourceLoss::new_controlled(
            n,
            ports,
            &self.l_tilde,
            &inputs,
            &viscous,
            propagator.dt,
            channel_budget,
            &mut *cancelled,
        )
        .map_err(integral_error)?;
        let thermal = PreparedSourceLoss::new_controlled(
            n,
            ports,
            &self.l_tilde,
            &inputs,
            &thermal,
            propagator.dt,
            channel_budget,
            &mut *cancelled,
        )
        .map_err(integral_error)?;
        Ok(DgPorousHeat {
            operator_identity: self.operator_identity,
            ports: self.ports.iter().map(|(id, _)| *id).collect(),
            dimension: n,
            viscous,
            thermal,
        })
    }
}
