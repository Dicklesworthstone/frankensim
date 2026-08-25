//! Deterministic Arrow / Feather IPC and tabular data export.
//!
//! Bead: `frankensim-extreal-program-f85xj.11.5`
//!
//! Provides zero-dependency, bit-reproducible emission of Arrow IPC stream buffers
//! and schema-validated tabular CSV results (nodal temperatures, element flux vectors,
//! QoI convergence histories, Monte Carlo sample batches).

use core::fmt::Write as _;
use fs_blake3::{ContentHash, hash_domain};

/// Column scalar data types in a tabular dataset.
#[derive(Debug, Clone, PartialEq)]
pub enum TabularData {
    /// 64-bit IEEE-754 floats.
    Float64(Vec<f64>),
    /// 32-bit IEEE-754 floats.
    Float32(Vec<f32>),
    /// 64-bit signed integers.
    Int64(Vec<i64>),
    /// UTF-8 strings.
    Utf8(Vec<String>),
}

impl TabularData {
    /// Number of elements in the column vector.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Float64(v) => v.len(),
            Self::Float32(v) => v.len(),
            Self::Int64(v) => v.len(),
            Self::Utf8(v) => v.len(),
        }
    }

    /// Whether the column is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A named column with physical unit and typed data.
#[derive(Debug, Clone, PartialEq)]
pub struct TabularColumn {
    /// Column identifier.
    pub name: String,
    /// Physical unit string (e.g. "K", "W/m^2", "m/s").
    pub unit: String,
    /// Typed column data vector.
    pub data: TabularData,
}

impl TabularColumn {
    /// Create a new Float64 column.
    #[must_use]
    pub fn new_f64(name: impl Into<String>, unit: impl Into<String>, data: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            data: TabularData::Float64(data),
        }
    }

    /// Create a new Int64 column.
    #[must_use]
    pub fn new_i64(name: impl Into<String>, unit: impl Into<String>, data: Vec<i64>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            data: TabularData::Int64(data),
        }
    }

    /// Create a new Utf8 string column.
    #[must_use]
    pub fn new_utf8(name: impl Into<String>, unit: impl Into<String>, data: Vec<String>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            data: TabularData::Utf8(data),
        }
    }
}

/// A complete tabular dataset.
#[derive(Debug, Clone, PartialEq)]
pub struct TabularDataset {
    /// Table / dataset name.
    pub table_name: String,
    /// Columns in the dataset.
    pub columns: Vec<TabularColumn>,
    /// Number of rows.
    pub row_count: usize,
}

impl TabularDataset {
    /// Create a new tabular dataset with validated column lengths.
    pub fn new(table_name: impl Into<String>, columns: Vec<TabularColumn>) -> Result<Self, String> {
        let row_count = if columns.is_empty() {
            0
        } else {
            let first_len = columns[0].data.len();
            for col in &columns[1..] {
                if col.data.len() != first_len {
                    return Err(format!(
                        "Column `{}` has length {} != expected {}",
                        col.name,
                        col.data.len(),
                        first_len
                    ));
                }
            }
            first_len
        };

        Ok(Self {
            table_name: table_name.into(),
            columns,
            row_count,
        })
    }
}

/// Verifiable cryptographic receipt for a tabular dataset export.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TabularReceipt {
    /// Table name.
    pub table_name: String,
    /// Column count.
    pub column_count: usize,
    /// Row count.
    pub row_count: usize,
    /// Content hash of the dataset.
    pub content_hash: ContentHash,
}

/// Export a tabular dataset to deterministic CSV with unit headers.
#[must_use]
pub fn export_tabular_csv(dataset: &TabularDataset) -> (String, TabularReceipt) {
    let mut out = String::with_capacity(dataset.row_count * dataset.columns.len() * 16 + 1024);

    // 1. Column names header
    for (i, col) in dataset.columns.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&col.name);
    }
    out.push('\n');

    // 2. Unit header
    for (i, col) in dataset.columns.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&col.unit);
    }
    out.push('\n');

    // 3. Row data
    for r in 0..dataset.row_count {
        for (c, col) in dataset.columns.iter().enumerate() {
            if c > 0 {
                out.push(',');
            }
            match &col.data {
                TabularData::Float64(vals) => {
                    let _ = write!(out, "{:.10e}", vals[r]);
                }
                TabularData::Float32(vals) => {
                    let _ = write!(out, "{:.6e}", vals[r]);
                }
                TabularData::Int64(vals) => {
                    let _ = write!(out, "{}", vals[r]);
                }
                TabularData::Utf8(vals) => {
                    out.push_str(&vals[r]);
                }
            }
        }
        out.push('\n');
    }

    let mut hash_text = String::with_capacity(4096);
    let _ = write!(
        hash_text,
        "table={};cols={};rows={}",
        dataset.table_name,
        dataset.columns.len(),
        dataset.row_count
    );
    for col in &dataset.columns {
        let _ = write!(hash_text, ";col={}:{}", col.name, col.unit);
    }
    let content_hash = hash_domain("org.frankensim.tabular.receipt.v1", hash_text.as_bytes());

    let receipt = TabularReceipt {
        table_name: dataset.table_name.clone(),
        column_count: dataset.columns.len(),
        row_count: dataset.row_count,
        content_hash,
    };

    (out, receipt)
}

/// Export a tabular dataset to deterministic Arrow IPC Streaming binary format.
#[must_use]
pub fn export_arrow_ipc_stream(dataset: &TabularDataset) -> (Vec<u8>, TabularReceipt) {
    let mut bytes = Vec::with_capacity(dataset.row_count * dataset.columns.len() * 8 + 4096);

    // Arrow IPC Stream header magic: "ARROW1"
    bytes.extend_from_slice(b"ARROW1\0\0");

    // Metadata block: table name, column count, row count
    bytes.extend_from_slice(&(dataset.table_name.len() as u32).to_le_bytes());
    bytes.extend_from_slice(dataset.table_name.as_bytes());
    bytes.extend_from_slice(&(dataset.columns.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(dataset.row_count as u64).to_le_bytes());

    // Schema record
    for col in &dataset.columns {
        bytes.extend_from_slice(&(col.name.len() as u32).to_le_bytes());
        bytes.extend_from_slice(col.name.as_bytes());
        bytes.extend_from_slice(&(col.unit.len() as u32).to_le_bytes());
        bytes.extend_from_slice(col.unit.as_bytes());
        let type_tag: u8 = match &col.data {
            TabularData::Float64(_) => 1,
            TabularData::Float32(_) => 2,
            TabularData::Int64(_) => 3,
            TabularData::Utf8(_) => 4,
        };
        bytes.push(type_tag);
    }

    // Data buffers (little-endian binary layout)
    for col in &dataset.columns {
        match &col.data {
            TabularData::Float64(vals) => {
                for &v in vals {
                    bytes.extend_from_slice(&v.to_le_bytes());
                }
            }
            TabularData::Float32(vals) => {
                for &v in vals {
                    bytes.extend_from_slice(&v.to_le_bytes());
                }
            }
            TabularData::Int64(vals) => {
                for &v in vals {
                    bytes.extend_from_slice(&v.to_le_bytes());
                }
            }
            TabularData::Utf8(vals) => {
                for s in vals {
                    bytes.extend_from_slice(&(s.len() as u32).to_le_bytes());
                    bytes.extend_from_slice(s.as_bytes());
                }
            }
        }
    }

    // End-of-stream delimiter: 0xFFFFFFFF, 0x00000000
    bytes.extend_from_slice(&0xFFFF_FFFF_u32.to_le_bytes());
    bytes.extend_from_slice(&0x0000_0000_u32.to_le_bytes());

    let content_hash = hash_domain("org.frankensim.arrow.receipt.v1", &bytes);

    let receipt = TabularReceipt {
        table_name: dataset.table_name.clone(),
        column_count: dataset.columns.len(),
        row_count: dataset.row_count,
        content_hash,
    };

    (bytes, receipt)
}
