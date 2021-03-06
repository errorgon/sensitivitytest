/*
 * Copyright (c) Numerical Method Inc.
 * http://www.numericalmethod.com/
 *
 * THIS SOFTWARE IS LICENSED, NOT SOLD.
 *
 * YOU MAY USE THIS SOFTWARE ONLY AS DESCRIBED IN THE LICENSE.
 * IF YOU ARE NOT AWARE OF AND/OR DO NOT AGREE TO THE TERMS OF THE LICENSE,
 * DO NOT USE THIS SOFTWARE.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITH NO WARRANTY WHATSOEVER,
 * EITHER EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION,
 * ANY WARRANTIES OF ACCURACY, ACCESSIBILITY, COMPLETENESS,
 * FITNESS FOR A PARTICULAR PURPOSE, MERCHANTABILITY, NON-INFRINGEMENT,
 * TITLE AND USEFULNESS.
 *
 * IN NO EVENT AND UNDER NO LEGAL THEORY,
 * WHETHER IN ACTION, CONTRACT, NEGLIGENCE, TORT, OR OTHERWISE,
 * SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIMS, DAMAGES OR OTHER LIABILITIES,
 * ARISING AS A RESULT OF USING OR OTHER DEALINGS IN THE SOFTWARE.
 */
package com.errorgon.sensitivitytest.math.matrix.doubles.linearsystem;

import com.errorgon.sensitivitytest.math.matrix.doubles.matrixtype.dense.triangle.UpperTriangularMatrix;
import com.errorgon.sensitivitytest.math.misc.SuanShuUtils;
import static com.errorgon.sensitivitytest.math.number.DoubleUtils.isZero;
import com.errorgon.sensitivitytest.math.vector.doubles.Vector;
import com.errorgon.sensitivitytest.math.vector.doubles.dense.DenseVector;

/**
 * Backward substitution solves a matrix equation in the form <i>Ux = b</i>
 * by an iterative process for an upper triangular matrix <i>U</i>.
 * The process is so called because for an upper triangular matrix, one first computes <i>x<sub>n</sub></i>,
 * then substitutes that backward into the next equation to solve for <i>x<sub>n-1</sub></i>,
 * and repeats until <i>x<sub>1</sub></i>.
 * Note that some diagonal entries in <i>U</i> can be 0s, provided that the system of equations is consistent.
 * For example,
 * \[
 * \begin{bmatrix}
 * 1 & 2 & 3\\
 * 0 & 0 & 5\\
 * 0 & 0 & 0
 * \end{bmatrix} \times
 * \begin{bmatrix}
 * 10\\
 * 0\\
 * 0
 * \end{bmatrix} =
 * \begin{bmatrix}
 * 10\\
 * 0\\
 * 0
 * \end{bmatrix}
 * \]
 *
 * @author Haksun Li
 * @see <a href="http://en.wikipedia.org/wiki/Triangular_matrix#Forward_and_Back_Substitution">Wikipedia: Forward and Back Substitution</a>
 */
public class BackwardSubstitution {

    /**
     * Solve <i>Ux = b</i>.
     *
     * @param U an upper triangular matrix, representing the system of linear equations (the homogeneous part)
     * @param b a vector
     * @return a solution <i>x</i> such that <i>Ux = b</i>
     * @throws LinearSystemSolver.NoSolution if there is no solution to the system
     */
    public Vector solve(UpperTriangularMatrix U, Vector b) {
        SuanShuUtils.assertArgument(U.nRows() == b.size(), "b must have the same length as U's dimension");

        int dim = b.size();
        DenseVector x = new DenseVector(dim);
        for (int i = dim; i >= 1; --i) {
            x.set(i, b.get(i));
            for (int j = i + 1; j <= dim; ++j) {
                x.set(i, x.get(i) - U.get(i, j) * x.get(j)); //x[i] -= U[i,j] * x[j];
            }

            if (!isZero(U.get(i, i), 0)) {
                x.set(i, x.get(i) / U.get(i, i));
            } else if (!isZero(x.get(i), 0)) {//U.get(i, i) == 0
                throw new LinearSystemSolver.NoSolution("no solution to this system of linear equations");
            }
        }

        return x;
    }
}
