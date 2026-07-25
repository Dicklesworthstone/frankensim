//! Synchronous raw-FrankenSQLite fixture access for integration tests.
//!
//! The public `Ledger` API is deliberately synchronous, while a few
//! conformance fixtures need to construct or corrupt an older database
//! through FrankenSQLite directly. The engine now exposes those raw storage
//! operations as futures. This adapter keeps the fixture call sites readable
//! without weakening production authority or pretending that these direct
//! writes are ledger operations.

// Each integration-test binary compiles this shared module independently and
// therefore uses a different subset of the raw fixture operations.
#![allow(dead_code)]

use std::marker::PhantomData;
use std::rc::Rc;

use fsqlite::{AsyncConnection, FrankenError, Row, SqliteValue};

/// Synchronous view of the raw engine connection used only by test fixtures.
pub(crate) struct SyncConnection {
    inner: AsyncConnection,
    _single_thread: PhantomData<Rc<()>>,
}

impl SyncConnection {
    pub(crate) fn open(path: impl Into<String>) -> Result<Self, FrankenError> {
        AsyncConnection::open_sync(path).map(|inner| Self {
            inner,
            _single_thread: PhantomData,
        })
    }

    pub(crate) fn query(&self, sql: &str) -> Result<Vec<Row>, FrankenError> {
        self.inner.query_sync(sql)
    }

    pub(crate) fn query_row(&self, sql: &str) -> Result<Row, FrankenError> {
        self.inner.query_row_sync(sql)
    }

    pub(crate) fn execute(&self, sql: &str) -> Result<usize, FrankenError> {
        self.inner.execute_sync(sql)
    }

    pub(crate) fn execute_with_params(
        &self,
        sql: &str,
        params: &[SqliteValue],
    ) -> Result<usize, FrankenError> {
        self.inner.execute_with_params_sync(sql, params)
    }

    pub(crate) fn execute_batch(&self, sql: &str) -> Result<(), FrankenError> {
        self.inner.execute_batch_sync(sql)
    }

    pub(crate) fn begin_transaction(&self) -> Result<(), FrankenError> {
        self.inner.begin_transaction_sync()
    }

    pub(crate) fn commit_transaction(&self) -> Result<(), FrankenError> {
        self.inner.commit_transaction_sync()
    }
}
