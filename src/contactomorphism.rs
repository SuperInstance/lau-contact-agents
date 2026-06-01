//! Contactomorphisms: diffeomorphisms preserving the contact structure.

use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};

/// A contactomorphism is a diffeomorphism φ: M → M such that φ*ξ = ξ
/// where ξ = ker(α) is the contact structure.
///
/// Equivalently, φ*α = f·α for some nowhere-zero function f.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Contactomorphism {
    pub n: usize,
    /// Jacobian matrix of the transformation at a point.
    pub jacobian: DMatrix<f64>,
    /// Whether it's strictly contact (φ*α = α, f = 1).
    pub strict: bool,
}

impl Contactomorphism {
    /// Identity contactomorphism.
    pub fn identity(n: usize) -> Self {
        let dim = 2 * n + 1;
        Self {
            n,
            jacobian: DMatrix::identity(dim, dim),
            strict: true,
        }
    }

    /// Create from a Jacobian matrix.
    pub fn from_jacobian(n: usize, jacobian: DMatrix<f64>, strict: bool) -> Self {
        let dim = 2 * n + 1;
        assert_eq!(jacobian.nrows(), dim);
        assert_eq!(jacobian.ncols(), dim);
        Self { n, jacobian, strict }
    }

    /// Compose two contactomorphisms.
    pub fn compose(&self, other: &Contactomorphism) -> Contactomorphism {
        assert_eq!(self.n, other.n);
        Contactomorphism {
            n: self.n,
            jacobian: &self.jacobian * &other.jacobian,
            strict: self.strict && other.strict,
        }
    }

    /// Invert a contactomorphism.
    pub fn inverse(&self) -> Option<Contactomorphism> {
        let inv = self.jacobian.clone().try_inverse()?;
        Some(Contactomorphism {
            n: self.n,
            jacobian: inv,
            strict: self.strict,
        })
    }

    /// Check if the map preserves the contact structure at a point.
    ///
    /// φ*α = f·α means the pullback of α is proportional to α.
    pub fn preserves_contact_structure(&self, alpha_coeffs: &DVector<f64>) -> bool {
        let pullback = &self.jacobian.transpose() * alpha_coeffs;
        let alpha_norm = alpha_coeffs.norm();
        let pullback_norm = pullback.norm();
        if alpha_norm < 1e-10 || pullback_norm < 1e-10 {
            return false;
        }
        let normalized_alpha = alpha_coeffs / alpha_norm;
        let normalized_pullback = &pullback / pullback_norm;
        // They should be parallel (up to sign)
        let cross = normalized_alpha.dot(&normalized_pullback);
        (cross.abs() - 1.0).abs() < 1e-8
    }

    /// Dimension of the ambient manifold.
    pub fn dimension(&self) -> usize {
        2 * self.n + 1
    }

    /// Generate the standard contactomorphism: a translation in z-direction.
    pub fn z_translation(n: usize, dz: f64) -> Self {
        let dim = 2 * n + 1;
        let jac = DMatrix::identity(dim, dim);
        // Translation doesn't change the Jacobian (it's the identity)
        let _ = dz;
        Self {
            n,
            jacobian: jac,
            strict: true,
        }
    }

    /// Rotation in the (x_i, y_i) plane — preserves contact structure.
    pub fn xy_rotation(n: usize, i: usize, theta: f64) -> Self {
        let dim = 2 * n + 1;
        let mut jac = DMatrix::identity(dim, dim);
        let c = theta.cos();
        let s = theta.sin();
        jac[(i, i)] = c;
        jac[(i, n + i)] = -s;
        jac[(n + i, i)] = s;
        jac[(n + i, n + i)] = c;
        Self {
            n,
            jacobian: jac,
            strict: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_preserves_contact() {
        let phi = Contactomorphism::identity(1);
        assert!(phi.strict);
        let alpha = DVector::from_vec(vec![-2.0, 0.0, 1.0]);
        assert!(phi.preserves_contact_structure(&alpha));
    }

    #[test]
    fn test_identity_compose() {
        let id = Contactomorphism::identity(2);
        let result = id.compose(&id);
        assert!(result.strict);
    }

    #[test]
    fn test_identity_inverse() {
        let id = Contactomorphism::identity(1);
        let inv = id.inverse().unwrap();
        assert!(inv.strict);
    }

    #[test]
    fn test_xy_rotation_preserves_contact() {
        // Rotation in (x,y) plane preserves the symplectic form dx∧dy,
        // hence preserves the contact structure.
        // At a point where y=0: α = (0, 0, 1). After rotation, still (0, 0, 1).
        let phi = Contactomorphism::xy_rotation(1, 0, std::f64::consts::PI / 4.0);
        let alpha = DVector::from_vec(vec![0.0, 0.0, 1.0]); // At y=0
        assert!(phi.preserves_contact_structure(&alpha));
    }

    #[test]
    fn test_xy_rotation_compose_identity() {
        let phi = Contactomorphism::xy_rotation(1, 0, 0.0);
        let id = Contactomorphism::identity(1);
        assert_eq!(phi.jacobian, id.jacobian);
    }

    #[test]
    fn test_compose_invertible() {
        let phi = Contactomorphism::xy_rotation(1, 0, 1.0);
        let inv = phi.inverse().unwrap();
        let product = phi.compose(&inv);
        let id = Contactomorphism::identity(1);
        for i in 0..3 {
            for j in 0..3 {
                assert!((product.jacobian[(i, j)] - id.jacobian[(i, j)]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_z_translation() {
        let phi = Contactomorphism::z_translation(1, 5.0);
        assert!(phi.strict);
        let alpha = DVector::from_vec(vec![-1.0, 0.0, 1.0]);
        assert!(phi.preserves_contact_structure(&alpha));
    }

    #[test]
    fn test_non_preserving_map() {
        let mut jac = DMatrix::identity(3, 3);
        jac[(0, 2)] = 1.0;  // Mix z into x
        let phi = Contactomorphism::from_jacobian(1, jac, false);
        let alpha = DVector::from_vec(vec![0.0, 0.0, 1.0]);
        // This should NOT preserve in general
        // but it depends on the specific alpha — let's use a tricky one
        let alpha2 = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        assert!(!phi.preserves_contact_structure(&alpha2));
    }

    #[test]
    fn test_dimension() {
        assert_eq!(Contactomorphism::identity(1).dimension(), 3);
        assert_eq!(Contactomorphism::identity(2).dimension(), 5);
        assert_eq!(Contactomorphism::identity(3).dimension(), 7);
    }

    #[test]
    fn test_serde() {
        let phi = Contactomorphism::xy_rotation(2, 0, 1.5);
        let json = serde_json::to_string(&phi).unwrap();
        let back: Contactomorphism = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 2);
    }
}
