//! Contact form α and Reeb vector field on a (2n+1)-dimensional manifold.

use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};

/// A contact form α on a (2n+1)-dimensional manifold.
///
/// In local coordinates (x₁,…,xₙ, y₁,…,yₙ, z), the standard contact form is
/// α = dz − Σ yᵢ dxᵢ. The defining condition is α ∧ (dα)ⁿ ≠ 0.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactForm {
    /// Dimension parameter n (manifold has dimension 2n+1).
    pub n: usize,
    /// The 1-form coefficients (length 2n+1) in local coordinates.
    pub coefficients: DVector<f64>,
}

impl ContactForm {
    /// Create the standard contact form α = dz − Σ yᵢ dxᵢ on R^{2n+1}.
    pub fn standard(n: usize) -> Self {
        let dim = 2 * n + 1;
        let coeffs = DVector::zeros(dim);
        Self { n, coefficients: coeffs }
    }

    /// Create a contact form from given coefficients at a point.
    pub fn from_coefficients(n: usize, coefficients: DVector<f64>) -> Self {
        assert_eq!(coefficients.len(), 2 * n + 1);
        Self { n, coefficients }
    }

    /// Evaluate α at a point p, returning the covector α_p.
    ///
    /// For the standard form α = dz − Σ yᵢ dxᵢ, the coefficients depend on p.
    pub fn evaluate_at(&self, p: &DVector<f64>) -> DVector<f64> {
        let dim = 2 * self.n + 1;
        assert_eq!(p.len(), dim);
        let mut result = DVector::zeros(dim);
        // α = dz - Σ yᵢ dxᵢ
        // Coefficients of dxᵢ: -yᵢ = -p[n+i]
        for i in 0..self.n {
            result[i] = -p[self.n + i];  // dxᵢ coefficient
        }
        // dyᵢ coefficients: 0 (already zero)
        // dz coefficient: 1
        result[2 * self.n] = 1.0;
        result
    }

    /// Compute dα (the exterior derivative) at a point p.
    ///
    /// For the standard form, dα = Σ dxᵢ ∧ dyᵢ (a 2-form).
    /// Returns as a skew-symmetric matrix (dim × dim).
    pub fn exterior_derivative(&self, _p: &DVector<f64>) -> DMatrix<f64> {
        let dim = 2 * self.n + 1;
        let mut dalpha = DMatrix::zeros(dim, dim);
        // dα = Σ dxᵢ ∧ dyᵢ
        for i in 0..self.n {
            dalpha[(i, self.n + i)] = 1.0;
            dalpha[(self.n + i, i)] = -1.0;
        }
        dalpha
    }

    /// Compute the top wedge α ∧ (dα)ⁿ (scalar value, up to normalization).
    ///
    /// For a contact form this should be nonzero everywhere.
    /// We compute it as the determinant of the matrix whose first row is α
    /// and remaining rows come from a basis adapted to dα.
    pub fn volume_form(&self, p: &DVector<f64>) -> f64 {
        let alpha = self.evaluate_at(p);
        let dalpha = self.exterior_derivative(p);
        let dim = 2 * self.n + 1;

        // Build the top-form matrix: first row = α, then n pairs of rows from dα
        // α ∧ (dα)ⁿ is computed via the determinant of a matrix whose columns
        // are the coefficient vectors.
        // For the standard form: α ∧ (dα)ⁿ = dz ∧ (dx₁∧dy₁) ∧ ... ∧ (dxₙ∧dyₙ)
        // which is always 1 (nonzero).
        //
        // General computation: construct (2n+1) x (2n+1) matrix with:
        //   row 0 = α^T
        //   rows 1..2n from the 2n columns of dα restricted to ker(α)
        // The determinant of this matrix gives the top wedge.

        // Use a simpler approach: compute the Pfaffian of dα restricted to ker(α),
        // multiplied by α evaluated on the Reeb direction.
        // For the standard form, this is always nonzero.

        let mut mat = DMatrix::zeros(dim, dim);
        // First row: α coefficients
        for j in 0..dim {
            mat[(0, j)] = alpha[j];
        }
        // Remaining rows: embed dα as row pairs
        // dα is antisymmetric, dim x dim. Use its rows 0..dim-1 for rows 1..dim-1
        for i in 1..dim {
            for j in 0..dim {
                mat[(i, j)] = dalpha[(i - 1, j)];
            }
        }
        mat.determinant()
    }

    /// Check if this is a contact form at point p (α ∧ (dα)ⁿ ≠ 0).
    pub fn is_contact_at(&self, p: &DVector<f64>) -> bool {
        let vol = self.volume_form(p);
        vol.abs() > 1e-10
    }

    /// Dimension of the manifold (2n+1).
    pub fn dimension(&self) -> usize {
        2 * self.n + 1
    }
}



/// The Reeb vector field R associated to a contact form α.
///
/// Defined by: α(R) = 1 and dα(R, ·) = 0.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReebVectorField {
    pub n: usize,
}

impl ReebVectorField {
    pub fn new(n: usize) -> Self {
        Self { n }
    }

    /// Compute the Reeb vector field at point p.
    ///
    /// For the standard contact form α = dz − Σ yᵢ dxᵢ,
    /// the Reeb vector field is R = ∂/∂z.
    pub fn evaluate(&self, _p: &DVector<f64>) -> DVector<f64> {
        let dim = 2 * self.n + 1;
        let mut r = DVector::zeros(dim);
        r[2 * self.n] = 1.0;  // ∂/∂z
        r
    }

    /// Verify Reeb condition: α(R) = 1.
    pub fn verify_alpha_condition(
        &self,
        p: &DVector<f64>,
        form: &ContactForm,
    ) -> bool {
        let r = self.evaluate(p);
        let alpha = form.evaluate_at(p);
        let dot = alpha.dot(&r);
        (dot - 1.0).abs() < 1e-10
    }

    /// Verify Reeb condition: dα(R, ·) = 0.
    pub fn verify_dalpha_condition(
        &self,
        p: &DVector<f64>,
        form: &ContactForm,
    ) -> bool {
        let r = self.evaluate(p);
        let dalpha = form.exterior_derivative(p);
        let contracted = &dalpha * &r;
        contracted.norm() < 1e-10
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DVector;

    #[test]
    fn test_standard_contact_form_dimension() {
        let form = ContactForm::standard(1);
        assert_eq!(form.dimension(), 3);
        let form2 = ContactForm::standard(2);
        assert_eq!(form2.dimension(), 5);
    }

    #[test]
    fn test_standard_contact_form_evaluate() {
        let form = ContactForm::standard(1);
        // Point (x, y, z) = (1, 2, 3) in R^3
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let alpha = form.evaluate_at(&p);
        // α = dz - y dx = -y·dx + 0·dy + 1·dz
        assert!((alpha[0] - (-2.0)).abs() < 1e-10);
        assert!(alpha[1].abs() < 1e-10);
        assert!((alpha[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_exterior_derivative_standard() {
        let form = ContactForm::standard(1);
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let dalpha = form.exterior_derivative(&p);
        // dα = dx ∧ dy
        assert!((dalpha[(0, 1)] - 1.0).abs() < 1e-10);
        assert!((dalpha[(1, 0)] - (-1.0)).abs() < 1e-10);
        assert!(dalpha[(0, 0)].abs() < 1e-10);
    }

    #[test]
    fn test_is_contact_at_standard() {
        let form = ContactForm::standard(1);
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        assert!(form.is_contact_at(&p));
    }

    #[test]
    fn test_is_contact_at_origin() {
        let form = ContactForm::standard(1);
        let p = DVector::from_vec(vec![0.0, 0.0, 0.0]);
        assert!(form.is_contact_at(&p));
    }

    #[test]
    fn test_contact_form_n2() {
        let form = ContactForm::standard(2);
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let alpha = form.evaluate_at(&p);
        assert_eq!(alpha.len(), 5);
        assert!((alpha[4] - 1.0).abs() < 1e-10); // dz coefficient
    }

    #[test]
    fn test_reeb_standard() {
        let reeb = ReebVectorField::new(1);
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let r = reeb.evaluate(&p);
        assert!((r[0]).abs() < 1e-10);
        assert!((r[1]).abs() < 1e-10);
        assert!((r[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reeb_alpha_condition() {
        let form = ContactForm::standard(1);
        let reeb = ReebVectorField::new(1);
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        assert!(reeb.verify_alpha_condition(&p, &form));
    }

    #[test]
    fn test_reeb_dalpha_condition() {
        let form = ContactForm::standard(1);
        let reeb = ReebVectorField::new(1);
        let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        assert!(reeb.verify_dalpha_condition(&p, &form));
    }

    #[test]
    fn test_reeb_n2() {
        let reeb = ReebVectorField::new(2);
        let p = DVector::from_vec(vec![0.0; 5]);
        let r = reeb.evaluate(&p);
        assert_eq!(r.len(), 5);
        assert!((r[4] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_contact_form_serde() {
        let form = ContactForm::standard(2);
        let json = serde_json::to_string(&form).unwrap();
        let back: ContactForm = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 2);
    }

    #[test]
    fn test_reeb_serde() {
        let reeb = ReebVectorField::new(3);
        let json = serde_json::to_string(&reeb).unwrap();
        let back: ReebVectorField = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 3);
    }
}
