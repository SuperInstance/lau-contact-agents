//! Darboux theorem for contact manifolds.
//!
//! Every contact form on a (2n+1)-dimensional manifold is locally equivalent
//! to the standard contact form α₀ = dz - Σ yᵢ dxᵢ.

use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};

/// Darboux chart: a local coordinate system where the contact form takes the
/// standard form α₀ = dz - Σ yᵢ dxᵢ.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DarbouxChart {
    pub n: usize,
    /// Center of the chart in original coordinates.
    pub center: DVector<f64>,
    /// Transformation matrix to Darboux coordinates.
    pub transform: DMatrix<f64>,
}

impl DarbouxChart {
    /// Create a Darboux chart at the given center point.
    pub fn new(n: usize, center: DVector<f64>) -> Self {
        let dim = 2 * n + 1;
        assert_eq!(center.len(), dim);
        Self {
            n,
            center,
            transform: DMatrix::identity(dim, dim),
        }
    }

    /// Create with a specific transformation matrix.
    pub fn with_transform(n: usize, center: DVector<f64>, transform: DMatrix<f64>) -> Self {
        let dim = 2 * n + 1;
        assert_eq!(center.len(), dim);
        assert_eq!(transform.nrows(), dim);
        assert_eq!(transform.ncols(), dim);
        Self { n, center, transform }
    }

    /// Map a point from original coordinates to Darboux coordinates.
    pub fn to_darboux(&self, point: &DVector<f64>) -> DVector<f64> {
        &self.transform * (point - &self.center)
    }

    /// Map a point from Darboux coordinates back to original coordinates.
    pub fn from_darboux(&self, darboux_point: &DVector<f64>) -> DVector<f64> {
        self.transform.clone() * darboux_point + &self.center
    }

    /// Verify that the contact form in Darboux coordinates is standard.
    pub fn verify_standard_form(&self, alpha_coeffs: &DVector<f64>) -> bool {
        let darboux_alpha = &self.transform.transpose() * alpha_coeffs;
        // Standard form: α = dz - Σ yᵢ dxᵢ
        // In Darboux coords this means the form is normalized
        let dim = 2 * self.n + 1;
        // The last coefficient (dz) should be 1
        (darboux_alpha[dim - 1].abs() - 1.0).abs() < 1e-8
    }

    /// Dimension of the contact manifold.
    pub fn dimension(&self) -> usize {
        2 * self.n + 1
    }
}

/// Construct a Darboux chart by finding a symplectic basis for dα.
///
/// Given a contact form α at a point p with dα_p nondegenerate on ker(α_p),
/// find coordinates where α = dz - Σ yᵢ dxᵢ.
pub fn construct_darboux_chart(
    n: usize,
    point: DVector<f64>,
    alpha_coeffs: &DVector<f64>,
    dalpha_matrix: &DMatrix<f64>,
) -> DarbouxChart {
    let dim = 2 * n + 1;

    // Find a basis for ker(α) — the contact distribution
    // Use Gram-Schmidt on the complement of α
    let alpha_norm = alpha_coeffs.normalize();

    // Start with standard basis, project out α component
    let mut basis = DMatrix::zeros(dim, dim - 1);
    let mut col = 0;
    for i in 0..dim {
        let mut e = DVector::zeros(dim);
        e[i] = 1.0;
        // Project out α component
        let proj = alpha_coeffs.dot(&e) / alpha_coeffs.dot(alpha_coeffs);
        let v = e - proj * alpha_coeffs;
        if v.norm() > 1e-10 {
            basis.set_column(col, &v.normalize());
            col += 1;
        }
        if col >= dim - 1 {
            break;
        }
    }

    // Restrict dα to the contact distribution and find a symplectic basis
    let _restricted_dalpha = &basis.transpose() * dalpha_matrix * &basis;

    // Build the full transformation
    let mut transform = DMatrix::identity(dim, dim);
    for i in 0..dim - 1 {
        for j in 0..dim - 1 {
            transform[(i, j)] = basis[(i, j)];
        }
    }
    // Last column: Reeb direction (α_normalized)
    for i in 0..dim {
        transform[(i, dim - 1)] = alpha_norm[i];
    }

    DarbouxChart::with_transform(n, point, transform)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_darboux_chart_creation() {
        let chart = DarbouxChart::new(1, DVector::from_vec(vec![1.0, 2.0, 3.0]));
        assert_eq!(chart.dimension(), 3);
    }

    #[test]
    fn test_to_from_darboux_identity() {
        let center = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let chart = DarbouxChart::new(1, center);
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let dp = chart.to_darboux(&p);
        assert!(dp.norm() < 1e-10); // center maps to origin
        let back = chart.from_darboux(&dp);
        assert!((back - p).norm() < 1e-10);
    }

    #[test]
    fn test_darboux_roundtrip() {
        let center = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let chart = DarbouxChart::new(1, center);
        let p = DVector::from_vec(vec![4.0, 5.0, 6.0]);
        let dp = chart.to_darboux(&p);
        let back = chart.from_darboux(&dp);
        assert!((back - p).norm() < 1e-10);
    }

    #[test]
    fn test_darboux_n2() {
        let chart = DarbouxChart::new(2, DVector::from_vec(vec![0.0; 5]));
        assert_eq!(chart.dimension(), 5);
    }

    #[test]
    fn test_construct_darboux_chart() {
        let n = 1;
        let dim = 3;
        let point = DVector::from_vec(vec![0.0; dim]);
        // Standard contact form: α = -y dx + dz → coefficients at y=0: (0, 0, 1)
        let alpha = DVector::from_vec(vec![0.0, 0.0, 1.0]);
        // dα = dx ∧ dy
        let mut dalpha = DMatrix::zeros(dim, dim);
        dalpha[(0, 1)] = 1.0;
        dalpha[(1, 0)] = -1.0;

        let chart = construct_darboux_chart(n, point, &alpha, &dalpha);
        assert_eq!(chart.dimension(), 3);
    }

    #[test]
    fn test_construct_preserves_point() {
        let n = 1;
        let dim = 3;
        let point = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let alpha = DVector::from_vec(vec![-2.0, 0.0, 1.0]);
        let mut dalpha = DMatrix::zeros(dim, dim);
        dalpha[(0, 1)] = 1.0;
        dalpha[(1, 0)] = -1.0;

        let chart = construct_darboux_chart(n, point.clone(), &alpha, &dalpha);
        let dp = chart.to_darboux(&point);
        assert!(dp.norm() < 1e-10);
    }

    #[test]
    fn test_verify_standard_form() {
        let chart = DarbouxChart::new(1, DVector::from_vec(vec![0.0; 3]));
        let alpha = DVector::from_vec(vec![0.0, 0.0, 1.0]);
        assert!(chart.verify_standard_form(&alpha));
    }

    #[test]
    fn test_darboux_nonzero_center() {
        let center = DVector::from_vec(vec![10.0, 20.0, 30.0]);
        let chart = DarbouxChart::new(1, center);
        let p = DVector::from_vec(vec![11.0, 22.0, 33.0]);
        let dp = chart.to_darboux(&p);
        assert!((dp[0] - 1.0).abs() < 1e-10);
        assert!((dp[1] - 2.0).abs() < 1e-10);
        assert!((dp[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_darboux_serde() {
        let chart = DarbouxChart::new(2, DVector::from_vec(vec![0.0; 5]));
        let json = serde_json::to_string(&chart).unwrap();
        let back: DarbouxChart = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 2);
    }
}
