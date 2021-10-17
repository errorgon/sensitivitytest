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
package com.errorgon.sensitivitytest.math.matrix.doubles.factorization.eigen;

import java.util.List;

/**
 * A spectrum is the set of eigenvalues of a matrix.
 *
 * @author Haksun Li
 * @see
 * <ul>
 * <li><a href="http://en.wikipedia.org/wiki/Eigenvalue,_eigenvector_and_eigenspace">Wikipedia: Spectrum and eigenvectors</a>
 * <li><a href="http://en.wikipedia.org/wiki/Eigenvalue,_eigenvector_and_eigenspace#Spectrum">Wikipedia: Spectrum</a>
 * </ul>
 */
public interface Spectrum {

    /**
     * Get all the eigenvalues.
     *
     * @return the eigenvalues
     */
    public List<Number> getEigenvalues();
}
